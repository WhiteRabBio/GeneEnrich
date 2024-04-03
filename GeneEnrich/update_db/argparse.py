import argparse


def add_subparser_args(subparsers: argparse) -> argparse:
    """Add tool-specific arguments for update database.
    Args:
        subparsers: Parser object before addition of arguments specific to update database.
    Returns:
        parser: Parser object with additional parameters.
    """

    subparser = subparsers.add_parser("update_db",
                                      description="Update GO and KEGG database.",
                                      help="")

    subparser.add_argument("--database", nargs=None, type=str,
                           dest='database', default=None,
                           required=True,
                           choices=['KEGG', 'GO'],
                           help="Choosen from KEGG or GO database to update. ")
    subparser.add_argument("--species", nargs=None, type=str,
                           dest='species', default='human',
                           required=True,
                           help="Support human or mouse.")
    subparser.add_argument("--version", nargs=None, type=str,
                           dest='version', default=None,
                           required=True,
                           choices=['python', 'R'],
                           help="If python, will download latest database from website. " \
                                "If R, will extract database from existing database info.")

    return subparsers
