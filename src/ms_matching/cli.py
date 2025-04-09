# cli.py
import argparse

from ms_matching.match import do_match

def main():
    parser = argparse.ArgumentParser(description="ms-matching")
    subparsers = parser.add_subparsers(required=True)

    # Create the parser for the "gct_to_excel" command.
    parser_gct = subparsers.add_parser(
        "match", help="Convert a GCT file to Excel."
    )
    parser_gct.set_defaults(func=do_match)
    parser_gct.add_argument(
        "-o",
        "--outname",
        default="out.tsv",
        help="Output filename (default: %(default)s)",
    )
    parser_gct.add_argument("root_path", help="root_path")
    #

    # create the parser for something else
    parser_other = subparsers.add_parser("other", help="other")
    parser_other.set_defaults(func=lambda _: print("other"))
    #

    args = parser.parse_args()
    # Directly call the function associated with the chosen subcommand.
    args.func(args)


if __name__ == "__main__":
    main()

