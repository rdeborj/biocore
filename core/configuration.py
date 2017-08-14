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
