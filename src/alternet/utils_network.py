import pandas as pd 
import db as db
import os.path as op


def filter_aggregate(data, threshold_importance=0.3, threshold_frequency=5, importance_column='median_importance'):
    '''
    Filter out the top `threshold_importance` percent of data based on the importance column 
    and retain only edges that appear at least `threshold_frequency` times.

    Parameters:

        data : pd.DataFrame
            DataFrame containing the edges of the GRN inference, including importance values 
            and frequency of appearance after aggregation.

        threshold_importance : float
            Percentage (between 0 and 1) of top entries to retain based on the importance column.

        threshold_frequency : int
            Minimum number of times an edge must appear in the GRN inference to be retained.

        importance_column : str
            Name of the column containing importance values.

    Returns:
    
        pd.DataFrame
            Filtered DataFrame based on the specified importance and frequency thresholds.
    '''
    #only select edges that have been found in threshold_frequency number of times in grn inference
    freq_mask = data['frequency'] >= threshold_frequency
    
    # get value for threshold importance
    importance_threshold = data.loc[freq_mask, importance_column].quantile(1 - threshold_importance)
    
    return data.loc[freq_mask & (data[importance_column] >= importance_threshold)]

def filter_aggregated_networks(asaware_grn, canonical_grn, threshold_frequency=5, threshold_importance=0.3, importance_column='median_importance'):
    '''
    Filter out the top `threshold_importance` percent of data based on the importance column 
    and retain only edges that appear at least `threshold_frequency` times in both AS-aware 
    and canonical GRNs.

    Parameters:

        asaware_grn : pd.DataFrame
            DataFrame containing edges from the alternative splicing-aware GRN inference, 
            including importance and frequency values.

        canonical : pd.DataFrame
            DataFrame containing edges from the canonical GRN inference, including 
            importance and frequency values.

        threshold_importance : float
            Percentage (between 0 and 1) of top entries to retain based on the importance column.

        threshold_frequency : int
            Minimum number of times an edge must appear in the GRN inference to be retained.

        importance_column : str
            Name of the column containing importance values.

    Returns:

        tuple
            A tuple containing the filtered AS-aware and canonical GRN DataFrames.
    '''

    asaware_grn = filter_aggregate(asaware_grn, threshold_importance, threshold_frequency, importance_column)
    canonical_grn = filter_aggregate(canonical_grn, threshold_importance, threshold_frequency, importance_column)
    
    return asaware_grn, canonical_grn


def biomart_mapping(net, biomart, source_column='TF', target_column='target', as_aware=True):
    '''
    Map BioMart gene or transcript names to gene IDs in a gene regulatory network (GRN).

    Parameters:

        net : pd.DataFrame  
            DataFrame representing the gene regulatory network.

        biomart : pd.DataFrame  
            BioMart annotation DataFrame containing gene/transcript names and corresponding 
            gene/transcript IDs.

    Returns:

        pd.DataFrame  
            Updated GRN DataFrame with gene/transcript names mapped to their corresponding IDs.
    '''
    # Pre-select and de-duplicate biomart annotations
    biomart_source = biomart[["Transcript stable ID", "Gene stable ID", "Transcript name", "Gene name"]]
    biomart_target = biomart[["Gene stable ID", "Gene name"]].drop_duplicates("Gene stable ID")


    # Prepare masks 
    is_source_transcript = net[source_column].str.startswith('ENST')
    is_source_gene = net[source_column].str.startswith('ENSG')
    is_target_gene = net[target_column].str.startswith('ENSG')


    # Map source transcript -> names
    st_sub = net[is_source_transcript]
    st_sub = st_sub.merge(biomart_source, left_on=source_column, right_on="Transcript stable ID", how="left")

    # For st_sub where target is gene
    st_sub = st_sub.merge(biomart_target, left_on=target_column, right_on="Gene stable ID", how="left", suffixes=('_source', '_target'))


    # Map source gene -> names
    sg_sub = net[is_source_gene]
    sg_sub = sg_sub.merge(biomart_source.drop_duplicates("Gene stable ID"), left_on=source_column, right_on="Gene stable ID", how="left")

    # For sg_sub where target is gene
    sg_sub = sg_sub.merge(biomart_target, left_on=target_column, right_on="Gene stable ID", how="left", suffixes=('_source', '_target'))

    # Combine
    net = pd.concat([st_sub, sg_sub], ignore_index=True)

    # Rename columns uniformly
    net = net.rename(columns={
        "Transcript stable ID": "source_transcript",
        "Gene stable ID_source": "source_gene",
        "Gene stable ID_target": "target_gene",
        "Transcript name": "source_transcript_name",
        "Gene name_source": "source_gene_name",
        "Gene name_target": "target_gene_name"
    })

    return net



def add_edge_key(df, biomart, key='edge_key', type='AS', source_column='TF', target_column='target'):
    '''
    Add an edge key column to enable comparison between AS-aware network edges and canonical network edges.

    Parameters:

        df : pd.DataFrame
            DataFrame containing edges of the GRN inference.

        biomart : pd.DataFrame or str
            BioMart data (as a DataFrame or file path) used to map gene IDs to transcript IDs.

        key : str
            Name of the new edge key column to be added.

        type : str
            Specifies the network type: 'AS-Aware' for transcript-level edges or 'canonical' for gene-level edges.

    Returns:

        pd.DataFrame
            Input DataFrame with an added column containing the constructed edge key.

    '''
    df = biomart_mapping(df, biomart, source_column=source_column, target_column=target_column)

    # Pre-allocate columns (faster than assigning one-by-one)
    df = df.assign(
        type=type,
        **{key: df['source_gene'].astype(str) + '_' + df[target_column].astype(str)}
    )

    # Small sanity check
    if df[key].isnull().any():
        raise ValueError(f"Column '{key}' contains NaN values!")

    return df




def get_common_edges(gene_grn, transcript_grn, path=None, save=False):
    ''' 
    Compute the overlapping network between a gene regulatory network and a transcript regulatory network.

        Parameters:

            gene_grn : pd.DataFrame  
                gene-level regulatory network containing only gene nodes.

            transcript_grn : pd.DataFrame  
                isoform-level regulatory network containing genes and transcripts, with transcripts as transcription factors (TFs).

            path : str  
                File path to save the overlapping network.

        Returns:

            pd.DataFrame  
                Overlapping network containing all edges present in both the gene regulatory network and the transcript regulatory network.    
    ''' 
    
    # Use set intersection for faster membership checking
    gene_edges = set(gene_grn['edge_key'])
    transcript_edges = set(transcript_grn['edge_key'])

    common_edges = gene_edges.intersection(transcript_edges)

    # Only keep rows with common edge_keys
    overlap_gene_in_t = transcript_grn[transcript_grn['edge_key'].isin(common_edges)]
    overlap_gene_g = gene_grn[gene_grn['edge_key'].isin(common_edges)]

    # Merge results
    overlap = pd.concat([overlap_gene_in_t, overlap_gene_g], ignore_index=True)

    if save and path:
        overlap.to_csv(path, sep='\t', index=False)

    return overlap


def get_diff(gene_grn, transcript_grn, path=None, save=False):
    ''' 
    Compute two networks containing edges found exclusively in either the gene regulatory network 
    or the transcript regulatory network based on the edge key.

    Parameters:

        gene_grn : pd.DataFrame  
            gene-level regulatory network containing only gene nodes.

        transcript : pd.DataFrame  
            isoform-level regulatory network containing genes and transcripts, with transcripts as transcription factors (TFs).

        path : str  
            Path to save the resulting networks.

        save : bool  
            Flag indicating whether to save the resulting networks.

    Returns:

        tuple of pd.DataFrame  
            - DataFrame with edges found only in the canonical (gene-level) network.
            - DataFrame with edges found only in the AS-aware (isoform-level) network.
    '''
    # Convert edge_keys to sets for faster difference operation
    gene_edges = set(gene_grn['edge_key'])
    transcript_edges = set(transcript_grn['edge_key'])

    # Get the differences using set operations
    diff_gene_edges = gene_edges - transcript_edges
    diff_transcript_edges = transcript_edges - gene_edges

    # Filter rows based on the set differences
    diff_gene = gene_grn[gene_grn['edge_key'].isin(diff_gene_edges)]
    diff_transcript = transcript_grn[transcript_grn['edge_key'].isin(diff_transcript_edges)]

    # Save if required
    if save and path:
        diff_transcript.to_csv(f"{path}transcript_level_only.tsv", sep='\t', index=False, header=True)
        diff_gene.to_csv(f"{path}gene_level_only.tsv", sep='\t', index=False, header=True)

    return diff_gene, diff_transcript



def gene_isoform_mapping(data, gene_column='gene_id', transcript_column='transcript_id'):
    ''' 
    Map isoforms (transcripts) to their corresponding genes.

    Parameters:

        data : pd.DataFrame  
            DataFrame containing gene and transcript ID columns.

        gene_column : str  
            Name of the column with gene IDs.

        transcript_column : str  
            Name of the column with transcript IDs.

    Returns:

        pd.DataFrame  
            DataFrame representing the gene-to-isoform mapping.
    '''   
    # Group by gene_column and aggregate transcript IDs into lists
    gene_isoform_map = data.groupby(gene_column)[transcript_column].agg(['unique', 'count']).reset_index()

    # Convert to a dictionary format
    isoform_mapping = {
        gene: {'isoforms': isoforms.tolist(), 'nr_isoforms': count}
        for gene, isoforms, count in zip(gene_isoform_map[gene_column], gene_isoform_map['unique'], gene_isoform_map['count'])
    }

    return isoform_mapping



def get_isoform_distribution(tfs, grn, tf_columns = ['transcript_id', 'gene_id'], source_column = 'source', target_column='target'):
    
    '''
    Integrates transcription factor (TF) isoform information into a gene regulatory network (GRN) 
    and determines how many source nodes (transcripts) map to the same gene.

    Parameters:
    
        tfs : pd.DataFrame
            A DataFrame containing transcript and gene information for transcription factors (TFs).
        grn : pd.DataFrame
            The gene regulatory network (GRN), where the "source" column represents the transcript.
        tf_columns : list, optional (default=['transcript_id', 'gene_id'])
            Column names specifying transcript and gene identifiers in the TF DataFrame.

    Returns:
    
        pd.DataFrame
            A modified GRN DataFrame with an added 'nr_isoforms' column, which indicates 
            the number of isoforms associated with each source gene.

    '''
    grn_og = grn.copy(deep=True) 
    # Map gene to isoforms
    tfs_for_gene_mapping = tfs.loc[:, tf_columns]
    gene_mapping_tfs = gene_isoform_mapping(tfs_for_gene_mapping)
    gene_mapping_tfs = pd.DataFrame.from_dict(gene_mapping_tfs, orient='index').reset_index(names=['Gene'])
    gene_mapping_short = gene_mapping_tfs[['Gene', 'nr_isoforms']]

    # Separate the AS-aware and canonical networks
    as_aware_overlap = grn[grn[source_column].str.startswith('ENST')]
    can_overlap = grn[grn[source_column].str.startswith('ENSG')]
    
    # Ensure no data loss
    if len(as_aware_overlap) + len(can_overlap) != len(grn):
        raise ValueError('Invalid source ids in GRN!')

    # Get isoform distribution for AS-aware network
    if not as_aware_overlap.empty:
        as_aware_overlap = as_aware_overlap.merge(
            tfs_for_gene_mapping[['transcript_id', 'gene_id']], left_on=source_column, right_on='transcript_id', how='left'
        ).drop(columns=['transcript_id'])
        as_aware_overlap = as_aware_overlap.merge(
            gene_mapping_short, left_on='gene_id', right_on='Gene', how='left'
        ).drop(columns=['Gene'])
        
        if as_aware_overlap['nr_isoforms'].isnull().any():
            raise ValueError('Invalid Transcript Ids, cannot be mapped to source in AS-aware network.')

    # Get isoform distribution for canonical network
    if not can_overlap.empty:
        can_overlap['gene_id'] = can_overlap[source_column]
        can_overlap = can_overlap.merge(
            gene_mapping_short, left_on='gene_id', right_on='Gene', how='left'
        ).drop(columns=['Gene'])
        
        if can_overlap['nr_isoforms'].isnull().any():
            raise ValueError('Invalid Gene Ids, cannot be mapped to source in Canonical network.')

    # Combine both AS-aware and canonical network data
    if can_overlap.empty:
        grn = as_aware_overlap
    elif as_aware_overlap.empty:
        grn = can_overlap
    else:
        grn = pd.concat([as_aware_overlap, can_overlap], ignore_index=True)

    # Ensure no edges are lost in the process
    if len(grn) != len(grn_og):
        raise ValueError('Processing resulted in lost edges.')

    return grn




def isoform_categorization(transcript_tfs, gene_tfs, threshold_dominance=90, threshold_balanced=15):
    '''
    Categorizes transcript isoforms based on their contribution to total gene expression into 'dominant', 'balanced', 
    or 'non-dominant'. 

    Parameters
        transcript_tfs : pd.DataFrame
            A DataFrame containing transcript expression data with the following columns:
            - 'transcript_id': Unique identifier for each transcript isoform.
            - 'gene_id': Identifier for the gene associated with the transcript.
            - Expression columns: Multiple columns representing expression values for different samples.

        gene_tfs : pd.DataFrame
            A DataFrame containing gene expression data with the following columns:
            - 'gene_id': Identifier for the gene associated with the transcript.
            - Expression columns: Multiple columns representing expression values for different samples.

        threshold_dominance : int, optional (default=90)
            If an isoform contributes more than this percentage of the total gene expression, 
            it is classified as 'dominant'.

        threshold_balanced : int, optional (default=15)
            If the standard deviation of isoform expression percentages within a gene is below 
            this threshold, the isoforms are classified as 'balanced'.

    Returns:
    
    pd.DataFrame
        A DataFrame with the following additional columns:
        - 'median_expression_iso': Median expression value of the isoform.
        - 'median_expression_gene': Total median expression of the gene; sum of median_expression_iso of respectively mapped isoforms.
        - 'percentage': The contribution of each isoform to the total gene expression.
        - 'max_percentage': The highest isoform expression percentage for the gene.
        - 'min_percentage': The lowest isoform expression percentage for the gene.
        - 'std_percentage': The standard deviation of isoform expression percentages for the gene.
        - 'isoform_category': The classification of the isoform as:
            - 'dominant': If an isoform contributes more than `threshold_dominance%` of total gene expression.
            - 'balanced': If the standard deviation of expression percentages within the gene is below `threshold_balanced`.
            - 'semi-dominant': If an isoform has the highest percentage but does not exceed `threshold_dominance%`.
            - 'non-dominant': All other cases.
    '''

    expression_columns = transcript_tfs.columns[2:]
    df = transcript_tfs.copy(deep=True)
    df['sum_transcript_expression'] = df[expression_columns].sum(axis=1)

    isoform_expression = df.copy(deep=True).loc[:, ['transcript_id', 'gene_id', 'sum_transcript_expression']]


    expression_columns = gene_tfs.columns[1:]
    df = gene_tfs.copy(deep=True)
    df['sum_gene_expression'] = df[expression_columns].sum(axis=1)

    gene_expression = df.copy(deep=True).loc[:,['gene_id', 'sum_gene_expression']]
    
    isoform_categorization = pd.merge(isoform_expression, gene_expression, on='gene_id', suffixes=('_iso', '_gene'))
    isoform_categorization['percentage'] = isoform_categorization['sum_transcript_expression'] / isoform_categorization['sum_gene_expression'] * 100
    
    isoform_categorization['max_percentage'] = isoform_categorization.groupby('gene_id')['percentage'].transform('max')
    isoform_categorization['min_percentage'] = isoform_categorization.groupby('gene_id')['percentage'].transform('min')
    isoform_categorization['std_percentage'] = isoform_categorization.groupby('gene_id')['percentage'].transform('std')

    def classify_isoform(row):

        if row['percentage'] == row['max_percentage'] and row['percentage'] > threshold_dominance:
            return 'dominant'
        elif row['std_percentage'] < threshold_balanced:
            return 'balanced'
        else:
            return'non-dominant'

    isoform_categorization['isoform_category'] = isoform_categorization.apply(classify_isoform, axis=1)
    
    return isoform_categorization


def get_gene_cases(df):
    '''
    Categorize genes based on their isoform classifications.

    Categorization rules:
    - If at least one isoform is classified as 'dominant', the gene is labeled as 'dominant'.
    - If all isoforms are 'balanced', the gene is labeled as 'balanced'.
    - Otherwise, the gene is labeled as 'non-dominant'.

    Parameters:

        df : pd.DataFrame  
            DataFrame containing the following columns:
            - 'gene_id' (str): The identifier for the gene.
            - 'isoform_category' (str): The classification of each isoform.

    Returns:

        pd.DataFrame  
            DataFrame with two columns:
            - 'gene_id' (str): The gene identifier.
            - 'isoform_category' (str): The assigned category for the gene.
    '''
    def categorize_gene_cases(categories):
        
        if any(cat == 'dominant' for cat in categories):
            return 'dominant'
        elif all(cat == 'balanced' for cat in categories):
            
            return 'balanced'
        else: 
            return 'non-dominant'
        
    gene_cases = df.groupby('gene_id')['isoform_category'].apply(categorize_gene_cases).reset_index()
    gene_cases.rename(columns={'isoform_category' : 'gene_case'}, inplace=True)

    return gene_cases


def find_transcript(transcript_id, AS, common):
    '''
    Check in which of the given DataFrames the specified transcript is present.

    Parameters:

        transcript_id : str  
            The transcript ID to search for.

        AS : pd.DataFrame  
            DataFrame representing the AS-aware network.

        common : pd.DataFrame  
            DataFrame representing the common network.

    Returns:

        dict  
            Dictionary indicating the presence of the transcript in each DataFrame, 
            {'AS': True/False, 'common': True/False}.
    '''

    presence = {
        'AS': transcript_id in AS['source'].values,
        'common': transcript_id in common['source'].values
    }

    return presence
    
def find_gene(gene_id, AS, CAN, common, source_column = 'TF'):
    '''
    Check in which of the given DataFrames the specified gene is present.

    Parameters:

        gene_id : str  
            The gene ID to search for.

        AS : pd.DataFrame  
            DataFrame representing the AS-aware network.

        CAN : pd.DataFrame  
            DataFrame representing the canonical network.

        common : pd.DataFrame  
            DataFrame representing the common network.

    Returns:

        dict  
            Dictionary indicating the presence of the gene in each DataFrame, 
            {'AS': True/False, 'CAN': True/False, 'common': True/False}.
    '''

    presence = {
        'AS': gene_id in AS['source_gene'].values,
        'CAN' : gene_id in CAN[source_column].values,
        'common': gene_id in common[source_column].values
    }

    return presence
    
def get_isoforms(gene_id, gene_iso_mapping):
    '''
    Return all isoforms mapped to a given gene.

    Parameters:

        gene_id : str  
            Ensembl Gene ID.

        gene_iso_mapping : pd.DataFrame  
            DataFrame containing gene-to-isoform mappings.

    Returns:

        list  
            List of all isoforms mapped to the specified gene.
    '''
    isoforms = gene_iso_mapping[gene_id]['isoforms']
    return isoforms

def get_isoform_category(data, transcript_id):
    '''
    Retrieve the isoform category for a given transcript ID.

    Parameters:

        data : pd.DataFrame  
            DataFrame containing transcript-related information.

        transcript_id : str  
            Transcript ID to search for.

    Returns:

        str or None  
            The isoform category if found; otherwise, None.
    '''
    result = data.loc[data['transcript_id'] == transcript_id, 'isoform_category']
    return result.iloc[0] if not result.empty else None

def plausibility_filtering(config, isoform_unique, isoform_categories, tf_database):

    results_dir = op.join(config['results_dir'], config['tissue'])
    results_dir_grn = op.join(results_dir, 'grn')

    as_aware_prefix = f"{config['tissue']}_as-aware.network_"

    isoform_unique = isoform_unique.merge(isoform_categories.loc[:, ['transcript_id', 'isoform_category']], left_on='source', right_on='transcript_id', how='left').drop(columns=['transcript_id'])

    isoform_unique_annotated = db.merge_annotations_to_grn(isoform_unique, tf_database)

    file = op.join(results_dir_grn, as_aware_prefix + f"annotated_w_database.tsv")
    isoform_unique_annotated.to_csv(file, sep='\t', index=False)



    #filter
    protein_coding = isoform_unique_annotated[isoform_unique_annotated['Protein Coding'] == True]
    protein_coding = protein_coding.drop(columns=['mean_importance', 'source_transcript', 'target_gene', 'type', 'edge_key'])
    protein_coding = protein_coding[protein_coding['isoform_category'] != 'dominant']
    sort_pc = protein_coding.sort_values(by=['frequency', 'median_importance'], ascending=[False, False])

    file = op.join(results_dir_grn, as_aware_prefix + f"plausibility_filtered_iso_unique.tsv")
    sort_pc.to_csv(file, sep='\t', index=False)
    return sort_pc


    

