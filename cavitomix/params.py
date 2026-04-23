# -*- coding: utf-8 -*-


class PluginParams:
    """
    class to hold general constants and keys used in the plugin
    """

    __version__ = "1.0"
    __year__ = "2022"
    __copyright__ = f"(c) {__year__} Innophore GmbH, 8010 Graz, Austria"
    __developers__ = ["G. Steinkellner", "C.C. Gruber", "K. Gruber"]
    # constants
    MAX_BIONEMO = 400
    MAX_ESMFOLD = 400
    URL_ESMFOLD = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    MAX_PYMOLFOLD = 500
    URL_PYMOLFOLD = "http://region-8.seetacloud.com:17537/predict/"
    # parameters for accessing the CavitOmiX server
    URL_INNOCLOUD = "http://structures.cloud.innophore.com/"
    # URL_INNOCLOUD = "http://127.0.0.1:5000"
    PLUGIN_KEY = "2bd3b462-7d17-469e-8f27-ec8c8d95d2ab"
    SERVER_KEY = ""
    # allowed amino acid codes
    AA_CODES = set("ACDEFGHIKLMNPQRSTVWY:")

    # available servers
    SERVER_LIST = []

    # annotation names
    annotation_names = {"CP": "Coulomb potential", "HP": "Hydrophobicity"}
    annotation_keys = {value: key for key, value in annotation_names.items()}
