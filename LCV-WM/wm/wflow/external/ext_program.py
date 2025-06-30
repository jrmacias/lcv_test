"""
Base class for external programs
"""
from ..pipeline import Program


class ExternalProgram(Program):

    name = ""
    command_line = ""
