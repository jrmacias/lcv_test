"""
Base class for pipelines
"""
from .program import Program


class Pipeline():
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
