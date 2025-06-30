"""
Base class for external programs
"""
from ..program import Program


class ExternalProgram(Program):

    name = ""
    command_line = ""
