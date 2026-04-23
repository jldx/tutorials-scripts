# -*- coding: utf-8 -*-
"""
module with pymol related functions
"""
from logging import getLogger
from typing import Any, List, Optional, Tuple

from pymol import cmd, util

from .cavfind import CavFind
from .util import get_unique_name

logger = getLogger(__name__)


def display_cavities(
    cavfind: CavFind,
    indices: List[int],
    annotations: List[str],
    show_env: bool = True,
    max_dist: float = 4.5,
) -> None:
    """
    display cavities stored in a CavFind object as non-bonded spheres and mesh
    """
    group_name = cavfind.name.upper()
    delete_group(group_name)

    # get current value of 'auto_zoom'
    auto_zoom = cmd.get("auto_zoom")
    cmd.set("auto_zoom", 0)

    obj_names = []
    for annotation in ["LIG"] + annotations:
        name = "{}_{}".format(cavfind.name, annotation).upper()
        obj_names.append(name)  # names of objects to be included in group

        pdbstr = "".join(
            [
                cav.get_pdbstr(annotation)
                for i, cav in enumerate(cavfind.cavities)
                if i in indices
            ]
        )

        cmd.read_pdbstr(pdbstr, name)  # create a new PyMOL object
        sele = f"({name})"
        cmd.unbond(sele, sele)  # remove all bonds between cavity points
        cmd.show_as("nb_spheres", sele)  # show cavities as non-bonded spheres
        cmd.show("mesh", sele)  # display mesh objects
        cmd.flag(24, sele, "clear")  # clear flags #24 and #25 to include cavity points
        cmd.flag(25, sele, "clear")  # in generating surface and mesh representations

        if annotation == "LIG":  # LigSite score
            cmd.spectrum(
                "b",
                palette="cyan_magenta",
                selection=sele,
                minimum=3.0,
                maximum=7.0,
            )

            if show_env:  # show amino acid side chains around cavities
                struct_sele = "".join(
                    [
                        f"{cavfind.obj_name}",
                        "" if cavfind.settings["keep_waters"] else " and not solvent",
                        "" if cavfind.settings["keep_hetatms"] else " and not het",
                        ""
                        if cavfind.settings["keep_hydrogens"]
                        else " and not hydrogen",
                    ]
                )
                env_name = f"{cavfind.name.upper()}_RES"
                env_selection = f"byres ({struct_sele}) near_to {max_dist} of {name}"
                cmd.select(env_name, env_selection)
                cmd.show("sticks", env_name)
                cmd.show("nb_spheres", env_name)
                cmd.color("white", env_name)  # just an example color
                util.cnc(env_name)  # color by atom type
                cmd.disable(env_name)
                obj_names.append(env_name)

        elif annotation == "HP":  # hydrophobicity
            cmd.spectrum("b", palette="rainbow", selection=sele)

        elif annotation == "CP":  # Coulomb potential
            cmd.spectrum(
                "b",
                palette="red_white_blue",
                selection=sele,
                minimum=-cavfind.settings["cp_limit"],
                maximum=cavfind.settings["cp_limit"],
            )

    cmd.group(group_name, " ".join(obj_names), action="add")  # create a new group
    cmd.group(group_name, action="open")

    cmd.set("auto_zoom", auto_zoom)  # reset 'auto_zoom'


def display_model(pdb_string: str, name: str) -> str:
    """
    display the predicted structure as cartoon colored according to pLDDT
    """
    obj_name = get_unique_name(name, cmd.get_names("public"))
    cmd.read_pdbstr(pdb_string, obj_name)

    plddt, factor = cal_plddt(pdb_string)
    logger.info(f"Mean Calpha-pLDDT of {obj_name}: {plddt:.2f}")

    if factor > 1.0:  # scaling of B-factor values
        cmd.alter(f"({obj_name})", f"b = b * {factor}")

    # AlphaFold coloring scheme for pLDDT
    cmd.set_color("high_lddt_c", [0.0, 0.3255, 0.8431])
    cmd.set_color("normal_lddt_c", [0.3412, 0.7922, 0.9765])
    cmd.set_color("medium_lddt_c", [1.0, 0.8588, 0.0706])
    cmd.set_color("low_lddt_c", [1.0, 0.49412, 0.2706])

    cmd.color("high_lddt_c", f"({obj_name}) and (b>90 or b=90)")
    cmd.color("normal_lddt_c", f"({obj_name}) and ((b<90 and b>70) or b=70)")
    cmd.color("medium_lddt_c", f"({obj_name}) and ((b<70 and b>50) or b=50)")
    cmd.color("low_lddt_c", f"({obj_name}) and (b<50)")

    return obj_name


def cal_plddt(pdb_string: str) -> Tuple[float, float]:
    """
    read b-factors of CA atoms from a PDB-formatted string and return
    the mean value
    used here to calculate the mean pLDDT of a predicted structure
    """

    plddts = [
        float(line[60:66])
        for line in pdb_string.split("\n")
        if line[:4] == "ATOM" and " CA " in line
    ]

    factor = 1.0
    if max(plddts) <= 1.0:
        factor = 100.0

    return sum(plddts) * factor / len(plddts), factor


def align_on_calpha(objects: List[str]) -> None:
    """
    superimpose a list of objects onto the first object
    """
    target = objects[0]
    logger.info(f"Align objects to {target} using CA positions.")

    cmd.zoom(target)
    for mobile in objects[1:]:
        cmd.align(f"{mobile} and name CA", f"{target} and name CA")


def delete_group(name) -> None:
    """
    delete a PyMOL group
    """
    groups = cmd.get_names_of_type("object:group")
    if name in groups:
        cmd.group(name, action="purge")  # delete all group objects


def get_coordinates(obj: str) -> Any:
    """
    get atom coordinates of a PyMOL object or selection as a numpy array
    """
    return cmd.get_coords(obj)


def get_object_list() -> List[str]:
    """
    get the list of current PyMOL objects
    """
    names = cmd.get_names("public")

    # identify cavity objects; they should not be used for cavity detection
    excluded_names = []
    for name in names:
        if name.endswith("_LIG"):
            excluded_names.append(name)
            root = name[:-4]
            excluded_names.extend([f"{root}_{tag}" for tag in ["CP", "HP", "RES"]])

    objects = []
    for name in names:
        # exclude cavity objects
        if name in excluded_names:
            continue

        obj_type = cmd.get_type(name)

        # only selections and molecular objects are added to the list
        if obj_type == "selection":
            objects.append(f"({name})")
        elif obj_type == "object:molecule":
            objects.append(name)

    return objects


def get_pdbstr(name: str) -> str:
    """
    get atoms of a PyMOL object or selection as PDB-string
    """
    return cmd.get_pdbstr(name)


def get_unique_object_name(name: str) -> str:
    """
    get a unique name for a PyMOL object
    """
    return get_unique_name(name, cmd.get_names("public"))


def create_object_from_pdbstr(pdbstr: str, name: str) -> str:
    """
    create a new PyMOL object from a PDB-formatted string
    """
    name = get_unique_object_name(name)
    cmd.read_pdbstr(pdbstr, name)
    return name


def create_object_from_pdbid(pdbid: str, name: str) -> Optional[str]:
    """
    create a new PyMOL object from a PDB-Id
    """
    name = get_unique_object_name(name)
    res = cmd.fetch(pdbid, name=name)
    return res if isinstance(res, str) else None


def get_pymol_version() -> Tuple:
    return cmd.get_version()


def is_valid_object(name: str) -> bool:
    """
    checks if an object or selection exists and contains atoms
    """
    return name in cmd.get_names("public") and cmd.count_atoms(name) > 0
