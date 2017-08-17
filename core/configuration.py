"""
Core function to be used for the Bioinformatics Core
"""

import os
import sys
import re
import subprocess
from yaml import Loader, Dumper, load, dump

def get_configuration(file):
    """
    Parse a YAML configuration file and get a dictionary containing
    configuration for the data analysis pipeline.

    Usage:
        get_configuraation(file='config.yaml')

    Input:
        * file: a YAML file containing key value configuration information

    Output:
    Returns a dictionary containing:
        samplename
        reference
        
    """

    with open(file, 'r') as stream:
        data = load(stream)

    return(data)

def is_validated(config):
    """
    Parse a configuration dictionary and validate the key-value pairs

    Usage:
        validate_configuration(config)

    Input:
        * config: configuration dictionary

    Output:
        Returns a boolean to determine whether the dictionary validates with
        all the correct keys.
    """

    # the dictionary should containing the following keys:
    #   * samplename
    #.  * reference
    keys = ['samplename', 'reference']

    for key, value in config.items():
        if key in keys:
            print(" ".join("Valid key: ", key))
        else:
            print("Invalid key: %s", key)
            invalid_keys.append(key)

    if (len(invalid_keys == 0)):
        return True
    else:
        return False
