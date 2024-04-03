"""Command-line tool functionality for update database."""

from GeneEnrich.base_cli import AbstractCLI
from GeneEnrich.run_kegg.kegg_utilities import update_kegg_db
from GeneEnrich.run_go.go_utilities import update_go_db
from GeneEnrich.enricher import organism_mapper
import os


class CLI(AbstractCLI):
    """CLI implements AbstractCLI from the GeneEnrich package."""

    def __init__(self):
        self.name = 'update_db'
        self.args = None

    def get_name(self) -> str:
        return self.name

    def validate_args(self, args):
        """Validate parsed arguments."""

        assert args.species in ['human', 'mouse'], \
            "GeneEnrich only support human and mouse for analysis."

        self.args = args

        return args

    def run(self, args):
        """Run the main tool functionality on parsed arguments."""

        # Run the tool.
        main(args)


def update_db(args):
    """The full script for the command line tool to update database.
    Args:
        args: Inputs from the command line, already parsed using argparse.
    Note: Returns nothing, but writes output to a file(s) specified from
        command line.
    """

    organism = organism_mapper(args.species)

    if args.database == 'KEGG':
        update_kegg_db(organism=organism, force=False, version=args.version)
    elif args.database == 'GO':
        update_go_db(organism=organism, force=False, version='R')

    print(f'Database update finished. ')


def main(args):
    """Take command-line input, parse arguments, and run tests or tool."""

    # Run the tool.
    update_db(args)
