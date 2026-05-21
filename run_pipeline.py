__doc__ = """ dummy run_pipeline.py script; apollo-reference doesn't have a CLI entry point to snakemake """

# Python Imports
from dataclasses import dataclass, field
from typing import Tuple
from rich.console import Console
from rich.markdown import Markdown

# Package Imports
from apollo_reference import __package_name__, __version__

class Pipeline:
    """ Dummy class to mimic the (working) Pipeline class from juno_library """
    # Required packages Imports
    # from juno_library import Pipeline
    pass

@dataclass
class ApolloReference(Pipeline):
    # see other juno_* and apollo_* repositories for examples how to augment the Pipeline class

    pipeline_name: str = __package_name__
    pipeline_version: str = __version__
    input_type: Tuple[str, ...] = ("fastq",)
    # ...

    def _add_args_to_parser(self) -> None:
        super()._add_args_to_parser()
        #self.parser.description = __description__
        # add custom argparse arguments here ...


def main() -> None:
    apollo_reference = ApolloReference()
    MARKDOWN_FILE = "README.md"
    console = Console()
    md = Markdown(open(MARKDOWN_FILE).read())
    console.print(md)
    import sys
    sys.exit()
    # "normal" activation of the pipeline
    apollo_reference.run()

if __name__ == "__main__":
    main()
