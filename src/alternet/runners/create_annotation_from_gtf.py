import argparse

import numpy as np
import pandas as pd
from gtfparse import read_gtf


FRAGMENT_TAGS = {'cds_start_NF', 'cds_end_NF', 'mRNA_start_NF', 'mRNA_end_NF'}


def read_gtf_transcripts(gtf_path):
    '''
    Reads the transcript lines of an Ensembl GTF into a biomart-like table.

    Parameters:
        gtf_path (str): Path to the Ensembl GTF (plain or .gz).

    Returns:
        pd.DataFrame: One row per transcript, with the column names of the old biomart export
            plus MANE, CCDS and fragment (`*_NF` tag) information.
    '''
    gtf = read_gtf(gtf_path, features={'transcript'}, result_type='pandas')
    gtf = gtf.replace('', np.nan)

    return gtf


def add_appris(annotation, appris_path):
    '''
    Adds APPRIS labels from an APPRIS data file (columns 'Transcript ID', 'APPRIS Annotation').

    'APPRIS Annotation' keeps the raw APPRIS label (PRINCIPAL:1, ALTERNATIVE:2, MINOR).
    'APPRIS annotation' uses the old biomart spelling (principal1, alternative2); MINOR is left empty,
    as in biomart.
    '''
    appris = pd.read_csv(appris_path, sep='\t', usecols=['Transcript ID', 'APPRIS Annotation'])
    appris = appris.drop_duplicates('Transcript ID')
    annotation = annotation.merge(appris, left_on='Transcript stable ID', right_on='Transcript ID', how='left')
    annotation = annotation.drop(columns='Transcript ID')

    label = annotation['APPRIS Annotation'].str.lower().str.replace(':', '', regex=False)
    annotation['APPRIS annotation'] = label.where(label.str.match(r'^(principal|alternative)', na=False))
    return annotation


def main():
    parser = argparse.ArgumentParser(description="Build a transcript annotation table from an Ensembl GTF (biomart replacement).")
    parser.add_argument("--gtf", type=str, required=True, help="Ensembl GTF, plain or .gz")
    parser.add_argument("--appris", type=str, required=False, help="APPRIS data file (e.g. data/appris_data.appris.txt)")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output TSV")
    args = parser.parse_args()

    annotation = read_gtf_transcripts(args.gtf)
    if args.appris is not None:
        annotation = add_appris(annotation, args.appris)

    annotation.to_csv(args.output, sep='\t', index=False)
    print(f"Wrote {len(annotation)} transcripts of {annotation['Gene stable ID'].nunique()} genes to {args.output}")


if __name__ == "__main__":
    main()
