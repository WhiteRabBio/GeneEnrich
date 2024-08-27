"""Command-line tool functionality for run_kegg."""

from GeneEnrich.base_cli import AbstractCLI
from GeneEnrich.run_kegg.enrich_kegg import enrichkegg
from GeneEnrich.preprocess import preprocess_deg, preprocess_protein
from GeneEnrich.enricher import organism_mapper
import os


class CLI(AbstractCLI):
    """CLI implements AbstractCLI from the GeneEnrich package."""

    def __init__(self):
        self.name = 'run_kegg'
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

        args.updown = ['up', 'down'] if args.updown == 'both' else args.updown.split(',')

        if args.type == 'uniprot':
            args.updown = ['up']

        self.args = args

        return args

    def run(self, args):
        """Run the main tool functionality on parsed arguments."""

        # Run the tool.
        main(args)


def run_kegg(args):
    """The full script for the command line tool to run kegg.
    Args:
        args: Inputs from the command line, already parsed using argparse.
    Note: Returns nothing, but writes output to a file(s) specified from
        command line.
    """

    organism = organism_mapper(args.species, database='KEGG')

    for ud in args.updown:
        if args.type != 'uniprot':
            gene = preprocess_deg(deg_path=args.input_file,
                                  organism=organism,
                                  type=args.type,
                                  updown=ud,
                                  database='KEGG')
            name = 'SYMBOL'
        else:
            gene = preprocess_protein(protein_path=args.input_file,
                                      organism=organism,
                                      database='KEGG')
            name = 'UNIPROT'

        res = enrichkegg(
            gene,
            organism=organism,
            pvalueCutoff=float(args.pvalue),
            pAdjustMethod="fdr_bh",
            minGSSize=10,
            maxGSSize=500,
            qvalueCutoff=float(args.qvalue),
        )

        res.to_csv(f'{args.output_dir}/{args.prefix}_KEGG_Enrichment_{ud.upper()}_Result.xls', sep='\t', index=None)
        print(f'KEGG enrichment analysis for {ud}-regulated {name} finished. ')


def main(args):
    """Take command-line input, parse arguments, and run tests or tool."""

    # Run the tool.
    run_kegg(args)
