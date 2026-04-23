# -*- coding: utf-8 -*-
"""
module with functions for retrieving structures predicted by AlphaFold and ESMfold
"""
import asyncio
import json
import os
import platform
import re
import sys
import sysconfig
from logging import getLogger
from typing import Callable, List, Optional, Sequence, Tuple
from urllib import parse, request
from urllib.error import HTTPError, URLError

from .mol_viewer import (
    align_on_calpha,
    create_object_from_pdbid,
    create_object_from_pdbstr,
    display_model,
    get_pymol_version,
)
from .params import PluginParams

# asyncio.to_thread not available in Python3.8 and older
PY38 = sys.version_info.major == 3 and sys.version_info.minor <= 8

logger = getLogger(__name__)

# regex patterns for PDB and Uniprot accession codes
uniprot_pattern = re.compile(
    r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$",
    flags=re.IGNORECASE,
)
pdb_pattern = re.compile(r"[1-9][A-Z0-9]{3}$", flags=re.IGNORECASE)
# regex pattern to split FASTA title lines
fasta_pattern = re.compile(r"[^A-Z0-9_|.]", flags=re.IGNORECASE)


def retrieve_structures(
    query_input: str,
    server: str = "BioNeMo",
    align: bool = True,
    progress_callback: Optional[Callable] = None,
) -> None:
    """
    main dispatch function for structure retrieval, display and alignment
    """
    tasks = define_tasks(get_queries(query_input), server)

    results = asyncio.run(_retrieve_structures(tasks))

    obj_names = display_results(results)

    if obj_names:
        logger.info(f"{len(obj_names)} PyMOL object(s) created:")
        logger.info(f"{', '.join(obj_names)}")
        # if more than one PyMOL object has been created, align all objects to the first
        if len(obj_names) > 1 and align:
            align_on_calpha(obj_names)


async def _retrieve_structures(tasks: Sequence[Tuple]) -> Tuple:
    """
    asynchronous function for structure retrieval
    """
    return await asyncio.gather(*[func(data) for func, data in tasks])


def display_results(results: Tuple) -> List[str]:
    """
    display retrieved structures
    """
    # display functions depending on the type of structure
    display_funcs = {
        "model": display_model,
        "pdb_str": create_object_from_pdbstr,
        "pdb_id": create_object_from_pdbid,
    }

    obj_names = []
    for result in results:
        if result is None:  # we did not obtain a structure
            continue
        struct_type, name, struct_data = result
        obj_name = display_funcs[struct_type](struct_data, name)
        if obj_name is not None:
            obj_names.append(obj_name)

    return obj_names


def define_tasks(queries: Sequence[Tuple], server: str) -> Sequence[Tuple]:
    """
    create a list of tasks for asynchronous structure retrieval
    """
    # query functions
    query_funcs = {
        "sequence_ESMfold": query_esmfold,
        "sequence_BioNeMo": query_bionemo,
        "sequence_PyMolFold": query_pymolfold,
        "pdb": query_pdb,
        "uniprot": query_alphafold_db,
    }

    tasks = []
    for query_type, data in queries:
        if query_type == "sequence":
            query_type += f"_{server}"

        tasks.append((query_funcs[query_type], data))

    return tasks


def get_queries(query_input: str) -> Sequence[Tuple]:
    """
    interpret user input as different queries
    possible queries are pdb-id, uniprot-id, FASTA-formatted amino acid sequence
    everything else is interpreted as a (single) amino acid sequence
    """
    lines = [line.strip() for line in query_input.split("\n") if line.strip()]

    queries = []
    title = ""
    sequence = ""

    for line in lines:
        match_pdb = re.match(pdb_pattern, line)
        match_uniprot = re.match(uniprot_pattern, line)
        match_fasta = line.startswith(">")

        # delimiter, indicating a new sequence or new query
        if match_pdb or match_uniprot or match_fasta:
            if sequence:  # if we already assembled a sequence, store it and reset it
                queries.append(("sequence", (title, sequence)))
                sequence = ""
                title = ""
            if match_pdb:
                queries.append(("pdb", line))
            elif match_uniprot:
                queries.append(("uniprot", line))
            else:
                # split the FASTA title line at the first character that is non-alphanumeric
                # or in "_.|" and discard the leading ">"
                title = re.split(fasta_pattern, line, 2)[1]
                # deal with FASTA title from UniprotKB
                title = title.replace("|", "_")
                if title.startswith("tr_") or title.startswith("sp_"):
                    title = title[3:]

        else:  # everything else is interpreted as an amino acid sequence
            sequence += line

    # append the last sequence if applicable
    if sequence:
        queries.append(("sequence", (title, sequence)))

    return queries


def process_sequence(sequence: str, max_aa: int, title: str = "") -> str:
    """
    process and check a single sequence
    remove whitespace, concatenate lines, test for invalid characters
    """
    # remove "internal" whitespace characters
    sequence = "".join(re.split(r"\s", sequence))
    # ":" denotes a chain break, "/" is converted to ":"
    sequence = re.sub(r"/", ":", sequence)
    # consecutive colons are converted to a single ":"
    sequence = re.sub(r":+", ":", sequence)
    # colons at the beginning or the end of the sequence are removed
    sequence = re.sub(r"^:+", "", sequence)
    sequence = re.sub(r":+$", "", sequence)
    sequence = sequence.upper()

    diff = set(sequence) - PluginParams.AA_CODES
    if diff:
        logger.error(f"Sequence: {title}")
        logger.error(sequence)
        logger.error(f"contains invalid characters:\n{diff}")
        return ""

    num_aa = len(sequence)
    if num_aa > max_aa:
        logger.error(f"Sequence: {title}")
        logger.error(sequence)
        logger.error(f"is too long: {num_aa:d} > {max_aa:d} amino acids!")
        return ""

    return sequence


async def query_pdb(data: str) -> Optional[Tuple]:
    """
    get a PDB structure from the server
    """
    pdb_id = data.upper()
    logger.info(f"Fetching structure {pdb_id:s} from the Protein Data Bank (PDB)")

    if PluginParams.SERVER_KEY:
        response = await query_innocloud("pdb", pdb_id)
        if response is not None:
            pdb_string = response.decode("utf-8")
            if pdb_string.startswith(f"REMARK   1 {PluginParams.SERVER_KEY}"):
                return "pdb_str", pdb_id, pdb_string
            else:
                logger.error(f"Unexpected response from PDB-server:\n{pdb_string:s}")
        else:
            logger.error("Error accessing PDB-server.")
    else:
        return "pdb_id", pdb_id, pdb_id


def _write_message(name: str, sequence: str, server: str) -> None:
    logger.info(f"Predicting structure for sequence '{name}' using {server}")
    seq_length = len(sequence)
    seq_text = sequence if seq_length <= 30 else f"{sequence[:12]}...{sequence[-12:]}"
    logger.info(f"{seq_length} amino acids: {seq_text}")


async def query_esmfold(data: Tuple[str, str]) -> Optional[Tuple]:
    """
    predict protein structure using ESMFold
    """
    name = data[0]
    sequence = process_sequence(data[1], PluginParams.MAX_ESMFOLD, name)

    if not sequence:
        return

    if ":" in sequence:  # ESMFold cannot deal with chain breaks
        logger.error("Amino acid sequence contains chain breaks (':')")
        return

    if not name:
        name = "ESMfold_pred"

    _write_message(name, sequence, "ESMfold")

    if PluginParams.SERVER_KEY and "esmfold" in PluginParams.SERVER_LIST:
        response = await query_innocloud("esmfold", sequence)
    else:
        response = await make_request(
            PluginParams.URL_ESMFOLD,
            headers={"Content-Type": "application/x-www-form-urlencoded"},
            data=sequence.encode(),
        )

    if response is not None:
        pdb_string = response.decode("utf-8")
        if pdb_string.startswith("HEADER") or pdb_string.startswith(
            f"REMARK   1 {PluginParams.SERVER_KEY}"
        ):
            return "model", name, pdb_string
        else:
            logger.error(f"Unexpected response from ESMFold-server:\n{pdb_string:s}")
    else:
        logger.error("Error accessing ESMFold-server.")


async def query_bionemo(data: Tuple[str, str]) -> Optional[Tuple]:
    """
    predict protein structure using Nvidia's BioNeMo
    """
    name = data[0]
    sequence = process_sequence(data[1], PluginParams.MAX_BIONEMO, name)

    if not sequence:
        return

    if not name:
        name = "BioNeMo_pred"

    _write_message(name, sequence, "BioNeMo")

    if PluginParams.SERVER_KEY and "bionemo" in PluginParams.SERVER_LIST:
        response = await query_innocloud("bionemo", sequence)
        if response is not None:
            pdb_string = response.decode("utf-8")
            if pdb_string.startswith(f"REMARK   1 {PluginParams.SERVER_KEY}"):
                return "model", name, pdb_string
            else:
                logger.error(
                    f"Unexpected response from BioNeMo-server:\n{pdb_string:s}"
                )
        else:
            logger.error("Error accessing BioNeMo-server.")
    else:
        logger.error(
            "BioNeMo server is only accessible via the Innophore CavitOmiX cloud."
        )


async def query_pymolfold(
    data: Tuple[str, str], num_recycle: int = 3
) -> Optional[Tuple]:
    """
    predict protein structure using PyMolFold
    """
    name = data[0]
    sequence = process_sequence(data[1], PluginParams.MAX_ESMFOLD, name)

    if not sequence:
        return

    if not name:
        name = "PyMolFold_pred"

    _write_message(name, sequence, "PyMolFold")

    if PluginParams.SERVER_KEY and "pymolfold" in PluginParams.SERVER_LIST:
        response = await query_innocloud("pymolfold", sequence)
        if response is not None:
            pdb_string = response.decode("utf-8")
            if pdb_string.startswith(f"REMARK   1 {PluginParams.SERVER_KEY}"):
                return "model", name, pdb_string
            else:
                logger.error(
                    f"Unexpected response from PyMolFold-server:\n{pdb_string:s}"
                )
        else:
            logger.error("Error accessing PyMolFold-server.")
    else:
        headers = {
            "accept": "application/json",
            "content-type": "application/x-www-form-urlencoded",
        }
        params = {
            "sequence": "'" + sequence + "'",
            "num_recycles": int(num_recycle),
        }

        url = f"{PluginParams.URL_PYMOLFOLD}?{parse.urlencode(params)}"
        response = await make_request(url, headers=headers, data=b"")
        if response is not None:
            pdb_string = response.decode("utf-8").replace('"', "")
            if pdb_string.startswith("PARENT N/A\\n"):
                pdb_string = pdb_string.replace("PARENT N/A\\n", "")
                pdb_string = pdb_string.replace("\\n", "\n")
                return "model", name, pdb_string
            else:
                logger.error(
                    f"Unexpected response from PyMolFold-server:\n{pdb_string:s}"
                )
        else:
            logger.error("Error accessing PyMolFold-server.")


async def query_alphafold_db(data: str) -> Optional[Tuple]:
    """
    query the AlphaFold database @EBI, retrieve a predicted structure
    """
    uniprot_id = data.upper()
    logger.info(f"Fetching AlphaFold model for UniProt ID {uniprot_id:s}.")

    if PluginParams.SERVER_KEY:
        response = await query_innocloud("alphafold", uniprot_id)
    else:
        pdb_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
        response = await make_request(pdb_url)

    if not response:
        logger.error(f"No model found for accession number: {uniprot_id}")
    else:
        pdb_string = response.decode("utf-8")
        name = f"AF_{uniprot_id:s}"
        return "model", name, pdb_string


async def query_innocloud(source: str, data: str) -> Optional[bytes]:
    """
    general query to the Innophore server
    """
    headers = {"Content-Type": "application/x-www-form-urlencoded"}
    post_dict = {
        "server_key": PluginParams.SERVER_KEY,
        "plugin_key": PluginParams.PLUGIN_KEY,
        "data": data,
        "source": source,
    }
    response = await make_request(
        PluginParams.URL_INNOCLOUD,
        headers=headers,
        data=parse.urlencode(post_dict).encode(),
    )
    return response


async def make_request(
    url: str, headers: Optional[dict] = None, data: Optional[bytes] = None
) -> Optional[bytes]:
    """
    asynchronous wrapper of _make_request
    """
    if PY38:
        return await asyncio.get_running_loop().run_in_executor(
            None, _make_request, url, headers, data
        )
    else:
        return await asyncio.to_thread(_make_request, url, headers, data)


def _make_request(
    url: str, headers: Optional[dict] = None, data: Optional[bytes] = None
) -> Optional[bytes]:
    """
    send a general request to a server
    """
    req = request.Request(url, headers=headers or {}, data=data)
    try:
        with request.urlopen(req) as response:
            return response.read()
    except HTTPError as error:
        logger.error(f"HTTPError: {error.status}, {error.reason}")
    except URLError as error:
        logger.error(f"URLError: {error.reason}")
    except TimeoutError:
        logger.error("Request timed out")


def get_server_key() -> None:
    """
    obtain the server-key from the CavitOmiX-server
    """
    os_info = {
        "os.name": os.name,
        "sys.platform": sys.platform,
        "platform.system": platform.system(),
        "sysconfig.get_platform": sysconfig.get_platform(),
        "platform.machine": platform.machine(),
        "platform.architecture": platform.architecture(),
    }
    post_dict = {
        "server_key": "get",
        "python_info": sys.version,
        "plugin_info": PluginParams.__version__,
        "pymol_info": get_pymol_version(),
        "os_info": os_info,
    }
    response = _make_request(
        PluginParams.URL_INNOCLOUD,
        headers={"Content-Type": "application/x-www-form-urlencoded"},
        data=parse.urlencode(post_dict).encode(),
    )
    if response is not None:
        PluginParams.SERVER_KEY = response.decode("utf-8")
    else:
        logger.error("Error accessing CavitOmiX-server.")


def contact_inno_server() -> None:
    """
    contact the CavitOmiX-server with server-key and plugin-key
    """
    get_server_key()

    if PluginParams.SERVER_KEY:
        post_dict = {
            "source": "get",
        }
        response = _make_request(
            PluginParams.URL_INNOCLOUD,
            headers={"Content-Type": "application/x-www-form-urlencoded"},
            data=parse.urlencode(post_dict).encode(),
        )
        if response is not None:
            data = json.loads(response.decode("utf-8"))
            PluginParams.SERVER_LIST = data["server-list"]
        else:
            logger.error("Error accessing CavitOmiX-server.")
            PluginParams.SERVER_KEY = ""
            return

        post_dict = {
            "server_key": PluginParams.SERVER_KEY,
            "plugin_key": PluginParams.PLUGIN_KEY,
            "data": "info",
        }
        response = _make_request(
            PluginParams.URL_INNOCLOUD,
            headers={"Content-Type": "application/x-www-form-urlencoded"},
            data=parse.urlencode(post_dict).encode(),
        )
        if response is not None:
            logger.info(response.decode("utf-8"))
            logger.info("Available servers:")
            for server in PluginParams.SERVER_LIST:
                logger.info(server)
        else:
            logger.error("Error accessing CavitOmiX-server.")
            PluginParams.SERVER_KEY = ""
