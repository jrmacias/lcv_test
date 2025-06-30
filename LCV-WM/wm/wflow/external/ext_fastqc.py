"""
Module that implements an external program to run FastQC
"""
from .ext_program import ExternalProgram


class FastQC(ExternalProgram):
    name = "FastQC"
    command_line = f'echo "{name}"'
