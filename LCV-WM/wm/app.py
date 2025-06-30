"""
This is main module of the LCV-WM application
"""
from wflow.pipeline import Pipeline
from wflow.external import ext_fastqc


def main():
    """
    _summary_
    """
    pl = Pipeline()
    pl.name = "pipeline_illumina_virus"

    fastqc = ext_fastqc.FastQC()
    pl.add_step(fastqc)

    pl.run()


if __name__ == "__main__":
    main()
