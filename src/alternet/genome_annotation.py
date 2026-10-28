import numpy as np
from gtfparse import read_gtf


def read_annotation(gtf_file, features=None, usecols=None):
    """
    Reads a GENCODE/Ensembl GTF(.gz) file into a DataFrame using gtfparse.

    The attribute column is expanded into one column per attribute (gene_id, transcript_id,
    gene_name, ...). Attributes that occur several times per line (e.g. 'tag') are joined by ','.

    Parameters:
        gtf_file (str): Path to the .gtf or .gtf.gz file.
        features (set): Feature types to keep (e.g. {'exon'}). None keeps all lines.
        usecols (list): Columns to keep (e.g. ['seqname', 'start', 'gene_id']). None keeps all.

    Returns:
        pd.DataFrame: One row per GTF line, missing attributes as NaN.
    """
    gtf = read_gtf(gtf_file, features=features, usecols=usecols, result_type='pandas')
    return gtf.replace('', np.nan)
