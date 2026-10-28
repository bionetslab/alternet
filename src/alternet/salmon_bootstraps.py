import os.path as op

import numpy as np
import pandas as pd
from pytximport.importers import read_salmon
from pytximport.utils import convert_counts_to_tpm


def salmon_sample_name(quant_path):
    """
    Returns the sample name of a Salmon output: the name of the directory that holds quant.sf
    (nf-core/rnaseq layout salmon/<sample>/quant.sf).

    Parameters:
        quant_path (str): Salmon output directory or its quant.sf file.

    Returns:
        str: Sample name.
    """
    quant_dir = quant_path if op.isdir(quant_path) else op.dirname(quant_path)
    return op.basename(op.normpath(quant_dir))


def read_salmon_tpm(quant_dirs):
    """
    Reads the Salmon TPM point estimates of several samples.

    The version suffix is removed from the transcript IDs.

    Parameters:
        quant_dirs (dict): Sample name -> Salmon output directory (or its quant.sf file).

    Returns:
        pd.DataFrame: transcript_ids (index) x samples.
    """
    tpm = {}
    for sample, quant_dir in quant_dirs.items():
        quant = read_salmon(quant_dir)
        ids = pd.Index(quant['transcript_id']).str.split('.').str[0].rename('transcript_id')
        tpm[sample] = pd.Series(quant['abundance'], index=ids).groupby(level=0).sum()
    return pd.DataFrame(tpm)


def read_salmon_bootstraps(quant_dirs, transcripts):
    """
    Reads the Salmon bootstrap replicates (salmon quant --numBootstraps) of several samples as TPM.

    The bootstrap counts are converted to TPM with the effective lengths of quant.sf, over all
    transcripts, before they are restricted to `transcripts`. Bootstrap replicate b of the output
    holds the b-th bootstrap of every sample. All samples must have the same number of bootstraps.

    Parameters:
        quant_dirs (dict): Sample name -> Salmon output directory (or its quant.sf file), with aux_info/.
        transcripts (list): Transcript IDs without version to keep.

    Returns:
        list: One samples (index) x transcripts DataFrame of TPM per bootstrap replicate.
    """
    replicates = {}
    for sample, quant_dir in quant_dirs.items():
        quant = read_salmon(quant_dir, inferential_replicates=True)
        ids = pd.Index(quant['transcript_id']).str.split('.').str[0]
        counts = quant['inferential_replicates']['replicates']
        tpm = np.column_stack([convert_counts_to_tpm(counts[:, b], quant['length']) for b in range(counts.shape[1])])
        replicates[sample] = pd.DataFrame(tpm, index=ids).groupby(level=0).sum().reindex(transcripts)

    n_bootstraps = {sample: tpm.shape[1] for sample, tpm in replicates.items()}
    if len(set(n_bootstraps.values())) != 1:
        raise ValueError(f'Samples have different numbers of bootstraps: {n_bootstraps}')

    return [pd.DataFrame({sample: tpm[b] for sample, tpm in replicates.items()}).T
            for b in range(next(iter(n_bootstraps.values())))]
