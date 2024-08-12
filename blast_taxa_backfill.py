#!/usr/bin/env python

"""blast_taxa_backfill.py: Fills in the taxonomic lineage of BLAST taxonomy IDs

This program accepts the tsv results of a BLASTn search. As long as the results are in tsv format and contain a taxaid columns, this program will backfill in their taxonomic lineage at the levels of Phylum, Class, Order, Family, and Genus.

The header argument should be formatted similar to this:
"6 qseqid sacc staxids sscinames bitscore pident qcovs"
The "6" indicates that the file is a tsv and "staxids" is the header for the taxonomy IDs off which this whole script works.
"""

# import the libraries
import argparse
from os import path
import sys
import pandas as pd
import numpy as np

try:
    from ete3 import NCBITaxa
except ImportError:
    raise Exception("This script requires the etetoolkit. Oh boy, looks like something broke")

__author__ = "Michael JV O'Mahoney"
__credits__ = ["Michael JV O'Mahoney", "Michael Trizna", "Matt Kweskin"]
__maintainer__ = "Michael JV O'Mahoney"
__status__ = "Prototype"


def get_args():
    # Create the parser
    parser = argparse.ArgumentParser(description="Fills out the higher taxonomic levels of BLAST matches",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Add the arguments
    parser.add_argument("INPUT", help='Input BLAST file')
    parser.add_argument("OUTPUT", help='Output BLAST file with full taxonomic lineages')
    parser.add_argument("header", type=str,
                        help='A string containing the headers you previously specified in the'
                             ' \'-outfmt\' argument of blastn')
    parser.add_argument("--overwrite",  action="store_true",
                        help='Automatically overwrite the output file and log file if they\'re '
                             'already present')

    # Parse the arguments
    return parser.parse_args()


# Some checks on the input/output/header/log files
def file_checks(args):
    if not path.exists(args.INPUT):
        raise Exception("Input file does not exist")

    # Store header argument and check if valid entry
    head = args.header
    if 'staxids' not in head:
        raise Exception("header argument is missing staxids. Please enter your "
                        "BLAST headers as they appear in the --outfmt argument of your blastn search")

    if path.exists(args.OUTPUT) and not args.overwrite:
        answer = input("Output file exists. " 
                       "(use --overwrite to Automatically overwrite the output file)\nOverwrite [y/n]? ")
        if answer.lower() != 'y':
            print("Exiting script.")
            sys.exit()

    outputdir = path.dirname(args.OUTPUT)
    if not path.exists(outputdir) and outputdir != "":
        raise Exception("Output directory does not exist")

    if head.startswith('6'):
        head = head.replace('6', '')
        head = head.split()
    elif head.startswith(('1', '2', '3', '4', '5', '7', '8', '9', '10','11', '12', '13', '14', '15', '16', '17', '18')):
        print('Input file is not in tsv format')
        sys.exit()
    else:
        head = head.split()

    return head


# I got the meat Jack!
def taxa_fill(original, modified, head):

    results = pd.read_csv(original, sep='\t', names=head)
    ncbi = NCBITaxa('/share/apps/bioinformatics/blast_taxa_backfill/taxa.sqlite')
    #ncbi.update_taxonomy_database()

    def grab_taxonomy(tax_id):
        tax_levels = ['phylum', 'class', 'order', 'family', 'genus', 'species']
        tax_dict = {'staxids': tax_id,
                    'phylum': np.nan,
                    'class': np.nan,
                    'order': np.nan,
                    'family': np.nan,
                    'genus': np.nan,
                    'species': np.nan}
        tax_id = str(int(tax_id))

        taxid_lineage = ncbi.get_lineage(tax_id)
        rank_dict = ncbi.get_rank(taxid_lineage)
        for k, v in rank_dict.items():
            if v in tax_levels:
                tax_dict[v] = ncbi.get_taxid_translator([k])[k]

        return tax_dict

    unique_taxids = results[pd.notnull(results['staxids'])]['staxids'].unique().tolist()

    tax_results = []
    for taxid in unique_taxids:
        tax_result = grab_taxonomy(taxid)
        tax_results.append(tax_result)

    tax_df = pd.DataFrame(tax_results)

    results_joined = results.merge(tax_df, on='staxids')

    results_joined.to_csv(modified, sep='\t', index=False)


# Let's make it work
def main():
    args = get_args()
    file_checks(args)
    head = file_checks(args)
    with open(args.INPUT) as original, open(args.OUTPUT, 'w') as modified:
        taxa_fill(original, modified, head)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()
