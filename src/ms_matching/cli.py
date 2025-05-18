# cli.py
import argparse
#from ms_matching.match import do_match
from ms_matching import driver

def main():
    parser = argparse.ArgumentParser(description="ms-matching")
    subparsers = parser.add_subparsers(required=True)

    # Create the parser for the "gct_to_excel" command.
    parser_driver = subparsers.add_parser(
        "run", help="run"
    )




    # arguments
    # out
    parser_driver.add_argument(
        "-o",
        "--outname",
        default="out.tsv",
        help="Output filename (default: %(default)s)",
    )

    # config
    parser_driver.add_argument(
        "-s",
        "--config",
        default=None, # TODO set a default
        help="Input config yaml file (default: %(default)s)",
    )

    parser_driver.add_argument("-d", "--db-path", default="fragments.db", help="SQLite database path")
    parser_driver.add_argument("-f", "--fasta", required=True, help="Input FASTA file")


    #parser_gct.add_argument("root_path", help="root_path")
    #parser_gct.add_argument("root_path", help="root_path")


    parser_driver.set_defaults(func=lambda args: driver.run(
        config_path=args.config,
        db_path=args.db_path,
        fasta_path=args.fasta,
        outname=args.outname,
    ))
    # end of driver parser


    # create the parser for something else
    parser_other = subparsers.add_parser("other", help="other")
    parser_other.set_defaults(func=lambda _: print("other"))
    #

    args = parser.parse_args()
    # Directly call the function associated with the chosen subcommand.
    args.func(args)


if __name__ == "__main__":
    main()

