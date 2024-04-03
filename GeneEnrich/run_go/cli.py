"""Command-line tool functionality for run_go."""

from GeneEnrich.base_cli import AbstractCLI
from GeneEnrich.run_go.enrich_go import enrichgo
from GeneEnrich.preprocess import preprocess_deg
from GeneEnrich.enricher import organism_mapper
import os


class CLI(AbstractCLI):
    """CLI implements AbstractCLI from the GeneEnrich package."""

    def __init__(self):
        self.name = 'run_go'
        self.args = None

    def get_name(self) -> str:
        return self.name

    def validate_args(self, args):
        """Validate parsed arguments."""

        try:
            args.input_file = os.path.expanduser(args.input_file)
            os.makedirs(args.output_dir, exist_ok=True)
        except TypeError:
            raise ValueError("Problem with provided input and output paths.")

        assert args.species in ['human', 'mouse'], \
            "GeneEnrich only support human and mouse for analysis."

        args.updown = ['up', 'down'] if args.updown == 'both' else args.updown.split(',')

        self.args = args

        return args

    def run(self, args):
        """Run the main tool functionality on parsed arguments."""

        # Run the tool.
        main(args)


def run_go(args):
    """The full script for the command line tool to run go.
    Args:
        args: Inputs from the command line, already parsed using argparse.
    Note: Returns nothing, but writes output to a file(s) specified from
        command line.
    """

    organism = organism_mapper(args.species)

    for ud in args.updown:
        gene = preprocess_deg(deg_path=args.input_file,
                              organism=organism,
                              type=args.type,
                              updown=ud,
                              database='GO')
        res = enrichgo(
            gene,
            organism=organism,
            pvalueCutoff=float(args.pvalue),
            pAdjustMethod="fdr_bh",
            minGSSize=10,
            maxGSSize=500,
            qvalueCutoff=float(args.qvalue),
        )

        res.to_csv(f'{args.output_dir}/{args.prefix}_GO_Enrichment_{ud.upper()}_Result.xls', sep='\t', index=None)
        print(f'GO enrichment analysis for {ud}-regulated genes finished. ')


def main(args):
    """Take command-line input, parse arguments, and run tests or tool."""

    # Run the tool.
    run_go(args)
