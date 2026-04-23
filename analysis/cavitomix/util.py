# -*- coding: utf-8 -*-
"""
module with utility functions
"""
from typing import Iterable


def get_unique_name(name: str, existing_names: Iterable[str]) -> str:
    """
    get a unique name
    if 'name' already exists f'_{i}' is appended to 'name' with increasing values of 'i'
    """
    i = 1
    tmp_name = name
    while tmp_name in existing_names:
        tmp_name = f"{name}_{i}"
        i += 1

    return tmp_name
