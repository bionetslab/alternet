import numpy as np
import pandas as pd


def read_psi(psi_file):
    """
    Reads a Whippet .psi(.gz) file into a DataFrame.

    The 'Coord' column (e.g. 'chr1:179369975-179370078') is split into
    'chromosome', 'start' and 'end'.

    Parameters:
        psi_file (str): Path to the .psi or .psi.gz file.

    Returns:
        pd.DataFrame: One row per node, with the coordinate columns added.
    """
    psi = pd.read_csv(psi_file, sep='\t', na_values='NA')
    coords = psi['Coord'].str.extract(r'^(?P<chromosome>.+):(?P<start>\d+)-(?P<end>\d+)$')
    psi['chromosome'] = coords['chromosome']
    psi['start'] = coords['start'].astype(int)
    psi['end'] = coords['end'].astype(int)
    return psi


def classify_nodes(psi, low=0.05, high=0.95):
    """
    Labels each Whippet node by its Psi value in a new column 'inclusion'.

    Psi < low is 'excluded', Psi > high is 'included' and everything in between
    (bounds inclusive) is 'variable'. Nodes without a Psi value get NaN.

    Parameters:
        psi (pd.DataFrame): Whippet nodes from read_psi.
        low (float): Psi below this counts as not included.
        high (float): Psi above this counts as fully included.

    Returns:
        pd.DataFrame: psi with the 'inclusion' column added.
    """
    psi = psi.copy()
    psi['inclusion'] = np.select([psi['Psi'] < low, psi['Psi'] > high, psi['Psi'].notna()],
                                 ['excluded', 'included', 'variable'], default=None)
    return psi


def read_exon_nodes(exon_file):
    """
    Reads the Whippet exon-to-node mapping (graph.exons.tab(.gz) from the Whippet index).

    The 'Potential_Exon' column (e.g. 'chr1:155991390-155991510:+') is split into
    'chromosome', 'start', 'end' and 'strand'. 'Whippet_Nodes' (e.g. '1,2,3') becomes a
    tuple of node numbers in 'nodes'.

    Parameters:
        exon_file (str): Path to the graph.exons.tab or graph.exons.tab.gz file.

    Returns:
        pd.DataFrame: One row per potential exon, with the coordinate and 'nodes' columns added.
    """
    exons = pd.read_csv(exon_file, sep='\t', dtype={'Whippet_Nodes': str})
    coords = exons['Potential_Exon'].str.extract(
        r'^(?P<chromosome>.+):(?P<start>\d+)-(?P<end>\d+):(?P<strand>[+-])$')
    exons['chromosome'] = coords['chromosome']
    exons['start'] = coords['start'].astype(int)
    exons['end'] = coords['end'].astype(int)
    exons['strand'] = coords['strand']
    exons['nodes'] = exons['Whippet_Nodes'].str.split(',').apply(lambda n: tuple(int(x) for x in n))
    return exons


def read_transcript_tpm(tpm_file):
    """
    Reads a merged Salmon transcript TPM table (nf-core/rnaseq *_salmon_transcript_tpm.csv).

    The version suffix is removed from the transcript IDs in the 'Name' column.

    Parameters:
        tpm_file (str): Path to the comma-separated TPM table.

    Returns:
        pd.DataFrame: transcript_ids (index) x samples.
    """
    tpm = pd.read_csv(tpm_file, index_col='Name')
    tpm.index = tpm.index.str.split('.').str[0].rename('transcript_id')
    return tpm


def is_expressed(tpm, min_tpm=1.0, min_fraction=0.25):
    """
    Finds the rows with TPM >= min_tpm in at least min_fraction of the samples.

    Parameters:
        tpm (pd.DataFrame): features (index) x samples.
        min_tpm (float): Minimal TPM of a feature in a sample.
        min_fraction (float): Minimal fraction of the samples with TPM >= min_tpm.

    Returns:
        pd.Series: True for the features that pass, same index as tpm.
    """
    return (tpm >= min_tpm).mean(axis=1) >= min_fraction


def expressed_transcripts(tpm_files, min_tpm=1.0, min_fraction=0.25):
    """
    Finds the transcripts with TPM >= min_tpm in at least min_fraction of the samples, per study.

    Parameters:
        tpm_files (dict): Study name -> path of its transcript TPM table (see read_transcript_tpm).
        min_tpm (float): Minimal TPM of a transcript in a sample.
        min_fraction (float): Minimal fraction of the samples of a study with TPM >= min_tpm.

    Returns:
        pd.DataFrame: transcript_ids (index) x studies, True if the transcript passes in that
            study. Transcripts missing from a study count as not passing. Use .any(axis=1)
            for the union over studies.
    """
    expressed = {}
    for study, tpm_file in tpm_files.items():
        expressed[study] = is_expressed(read_transcript_tpm(tpm_file), min_tpm, min_fraction)
    return pd.DataFrame(expressed).fillna(False).astype(bool)


def variable_nodes(psi_files, min_fraction=0.25, min_reads=5, low=0.05, high=0.95):
    """
    Finds the Whippet nodes that are variable in at least min_fraction of the samples, per study.

    In each sample, nodes with at least min_reads reads (Total_Reads) are labelled with
    classify_nodes. Nodes that are missing from a sample or have fewer reads count as not
    variable in that sample.

    Parameters:
        psi_files (dict): Study name -> list of Whippet .psi(.gz) files, one per sample.
        min_fraction (float): Minimal fraction of the samples of a study with a variable node.
        min_reads (int): Minimal Total_Reads of a node to count as evidence.
        low (float): Psi below this counts as not included.
        high (float): Psi above this counts as fully included.

    Returns:
        pd.DataFrame: (Gene, Node) (index) x studies, True if the node passes in that study.
            Use .any(axis=1) for the union over studies.
    """
    variable = {}
    for study, files in psi_files.items():
        counts = None
        for psi_file in files:
            psi = read_psi(psi_file)
            nodes = classify_nodes(psi.loc[psi['Total_Reads'] >= min_reads,
                                           ['Gene', 'Node', 'Psi']], low, high)
            sample = (nodes.loc[nodes['inclusion'] == 'variable', ['Gene', 'Node']]
                      .drop_duplicates().value_counts())
            counts = sample if counts is None else counts.add(sample, fill_value=0)
        variable[study] = counts / len(files) >= min_fraction
    return pd.DataFrame(variable).fillna(False).astype(bool)


def transcript_nodes(gtf, exon_nodes):
    """
    Finds the Whippet nodes each transcript includes, using the Whippet exon-to-node mapping.

    Every GTF exon is matched to a potential exon of its gene with the same start and end.
    The nodes of a transcript are the union of the nodes of its exons.

    Parameters:
        gtf (pd.DataFrame): GTF annotation from read_annotation.
        exon_nodes (pd.DataFrame): Exon-to-node mapping from read_exon_nodes.

    Returns:
        pd.DataFrame: One row per transcript with columns 'gene_id', 'transcript_id',
            'nodes' (sorted tuple of node numbers) and 'n_unmapped_exons' (exons without a
            matching potential exon).
    """
    exons = gtf.loc[gtf['feature'] == 'exon', ['gene_id', 'transcript_id', 'start', 'end']]
    exons = exons.merge(exon_nodes[['Gene', 'start', 'end', 'nodes']], how='left',
                        left_on=['gene_id', 'start', 'end'], right_on=['Gene', 'start', 'end'])
    exons['unmapped'] = exons['nodes'].isna()
    transcripts = exons.groupby(['gene_id', 'transcript_id']).agg(
        nodes=('nodes', lambda n: tuple(sorted(set(x for e in n.dropna() for x in e)))),
        n_unmapped_exons=('unmapped', 'sum')).reset_index()
    return transcripts


def informative_nodes(variable, expressed, transcripts, regulator_genes):
    """
    Finds the variable nodes of regulator genes that split the expressed isoforms in two groups.

    A node is informative in a study if
    - it is variable (see variable_nodes),
    - its gene is in regulator_genes,
    - at least one expressed isoform of its gene contains it (see expressed_transcripts),
    - at least one other expressed isoform of its gene does not contain it.
    Expressed isoforms are the transcripts in transcripts that pass expressed_transcripts.

    Parameters:
        variable (pd.DataFrame): Output of variable_nodes.
        expressed (pd.DataFrame): Output of expressed_transcripts, with the same studies.
        transcripts (pd.DataFrame): Output of transcript_nodes.
        regulator_genes (set): Gene IDs without version (e.g. the TF and TF_SF genes of
            the regulator list).

    Returns:
        pd.DataFrame: One row per study and informative node with columns 'study', 'Gene',
            'Node', 'n_expressed' (expressed isoforms of the gene) and 'n_with_node'
            (expressed isoforms that contain the node).
    """
    tx = transcripts.loc[transcripts['gene_id'].str.split('.').str[0].isin(regulator_genes)]
    tx = tx.assign(tx=tx['transcript_id'].str.split('.').str[0])
    tx_nodes = tx[['gene_id', 'tx', 'nodes']].explode('nodes').dropna(subset=['nodes'])
    tx_nodes['nodes'] = tx_nodes['nodes'].astype(int)

    informative = []
    for study in variable.columns:
        is_expressed = tx['tx'].isin(expressed.index[expressed[study]])
        n_expressed = tx.loc[is_expressed].groupby('gene_id').size()
        hits = tx_nodes.loc[tx_nodes['tx'].isin(tx.loc[is_expressed, 'tx'])].merge(
            variable.index[variable[study]].to_frame(index=False),
            left_on=['gene_id', 'nodes'], right_on=['Gene', 'Node'])
        nodes = hits.groupby(['Gene', 'Node']).size().rename('n_with_node').reset_index()
        nodes['n_expressed'] = nodes['Gene'].map(n_expressed)
        nodes = nodes.loc[nodes['n_with_node'] < nodes['n_expressed']]
        informative.append(nodes.assign(study=study))
    columns = ['study', 'Gene', 'Node', 'n_expressed', 'n_with_node']
    return pd.concat(informative, ignore_index=True)[columns]


def node_psi_matrix(psi_files, nodes, min_reads=5):
    """
    Collects the Psi of the given Whippet nodes in every sample.

    Nodes with fewer than min_reads reads (Total_Reads) in a sample, or missing from it, get NaN.

    Parameters:
        psi_files (dict): Sample name -> Whippet .psi(.gz) file.
        nodes (pd.DataFrame): Nodes with columns 'Gene' and 'Node' (e.g. from informative_nodes).
        min_reads (int): Minimal Total_Reads of a node to use its Psi.

    Returns:
        pd.DataFrame: samples (index) x nodes, with (Gene, Node) MultiIndex columns.
    """
    keep = pd.MultiIndex.from_frame(nodes[['Gene', 'Node']].drop_duplicates())
    psi_values = {}
    for sample, psi_file in psi_files.items():
        psi = read_psi(psi_file)
        psi = (psi.loc[psi['Total_Reads'] >= min_reads].drop_duplicates(['Gene', 'Node'])
               .set_index(['Gene', 'Node'])['Psi'])
        psi_values[sample] = psi.reindex(keep)
    return pd.DataFrame(psi_values).T
