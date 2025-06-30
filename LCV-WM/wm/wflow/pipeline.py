"""
Base classes for pipelines
"""


class Program():
    """
    Base class for a running program, a step in the pipeline
    """

    name = ""
    command_line = ""


class Pipeline():
    """
    Base class for pipeline
    """
    name = ""
    imputs = []
    outputs = []
    parameters = []
    pipeline_steps = []

    def run(self):
        print(f"-> Running {self.name}")
        for step in self.pipeline_steps:
            print(f"--> Running step: {step.name}")

    def add_step(self, step: Program):
        self.pipeline_steps.append(step)
