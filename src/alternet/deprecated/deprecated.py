
def filter_aggregated(data, top_percentile=0.3, threshold_frequency=5, importance_column='median_importance', frequency_column = 'frequency'):
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

    n_before = data.shape[0]
    #only select edges that have been found in threshold_frequency number of times in grn inference
    freq_mask = data[frequency_column] >= threshold_frequency
    
    # get value for threshold importance
    importance_threshold = data.loc[freq_mask, importance_column].quantile(1 - top_percentile)
    

    data =  data.loc[freq_mask & (data[importance_column] >= importance_threshold)]
    n_after = data.shape[0]

    filter_info = {'importance_freq': {'before': n_before, 'n_after': n_after}}
    return data.copy(), importance_threshold, filter_info




def map_transcript_to_gene(net, transcript_mapper, column = 'source', rename_column='source_transcript', output_column = 'source_gene'):
    net[output_column] = net[column].apply(lambda x: transcript_mapper[x])
    net = net.rename(columns={column:rename_column})
    return net


def create_ensg_to_geneid_mapping(biomart):
    transcript_mapper = dict(zip(biomart['Gene stable ID'], biomart['Gene name']))
    return transcript_mapper
    
def create_enst_to_tid_mapping(biomart):
    transcript_mapper = dict(zip(biomart['Transcript stable ID'], biomart['Transcript name']))
    return transcript_mapper


def add_transcript_names(edgelist, biomart, transcript_column = 'source_transcript'):
    transcript_name_mapper = create_enst_to_tid_mapping(biomart)
    edgelist[f'{transcript_column}_name'] = edgelist[transcript_column].apply(lambda x: transcript_name_mapper[x])
    return edgelist

def add_gene_names(edgelist, biomart, gene_column = 'target_gene'):
    gene_name_mapper = create_ensg_to_geneid_mapping(biomart)
    edgelist[f'{gene_column}_name'] = edgelist[gene_column].apply(lambda x: gene_name_mapper[x])
    return edgelist



def create_edge_key(net, source_column = 'source_gene', target_column = 'target'):
    net['edge_key'] = net[source_column]+'_' + net[target_column]
    return net



def get_common_edges(gene_grn, transcript_grn):
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

    return overlap


def get_diff(gene_grn, transcript_grn):
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

    return diff_gene, diff_transcript




def isoform_categorization(transcript_data, gene_data, tf_list, threshold_dominance=90, threshold_balanced=15):
    '''
    Categorizes transcript isoforms based on their contribution to total gene expression into 'dominant', 'balanced', 
    or 'non-dominant'. 

    Parameters
        transcript_data

        gene_data:

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

    gene_sums = pd.DataFrame(gene_data.sum())
    gene_sums.columns = ['gene_counts']

    transcript_sums = pd.DataFrame(transcript_data.sum())
    transcript_sums.columns = ['transcript_counts']

    transcript_sums = transcript_sums.merge(tf_list, left_index = True, right_on = 'Transcript stable ID')
    transcript_sums = transcript_sums.merge(gene_sums, left_on='Gene stable ID', right_index = True)
    transcript_sums['percentage'] = (transcript_sums['transcript_counts']/transcript_sums['gene_counts'])*100

    transcript_sums['max_percentage'] = transcript_sums.groupby('Gene stable ID')['percentage'].transform('max')
    transcript_sums['min_percentage'] = transcript_sums.groupby('Gene stable ID')['percentage'].transform('min')
    transcript_sums['std_percentage'] = transcript_sums.groupby('Gene stable ID')['percentage'].transform('std')

    def classify_isoform(row, threshold_balanced =15 , threshold_dominance = 80):

        if row['percentage'] == 100:
            return 'single'
        elif row['percentage'] == row['max_percentage'] and row['percentage'] > threshold_dominance:
            return 'dominant'
        elif row['std_percentage'] < threshold_balanced:
            return 'balanced'
        else:
            return'non-dominant'

    transcript_sums['isoform_category'] = transcript_sums.apply(classify_isoform, axis=1)

    transcript_sums = transcript_sums.loc[:, ['Gene stable ID', 'Transcript stable ID', 'percentage', 'isoform_category']]


    return transcript_sums



    
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




def filter_importance(data, absolute_treshold=0.3, importance_column='median_importance'):
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

    n_before = data.shape[0]

    data =  data[ (data[importance_column] >= absolute_treshold)]
    n_after = data.shape[0]

    filter_info = {'n_before': n_before, 'n_after_cutoff': n_after}
    return data.copy(), filter_info

def frequency_filter(data,  threshold_frequency=5,frequency_column = 'frequency'):
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

    n_before = data.shape[0]
    #only select edges that have been found in threshold_frequency number of times in grn inference
    data = data[data[frequency_column] >= threshold_frequency]
    
    n_after = data.shape[0]

    filter_info = {'before': n_before, 'n_after_frequency': n_after}
    return data.copy(), filter_info

   

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

        if any(cat == 'single' for cat in categories):
            return 'single'
        elif any(cat == 'dominant' for cat in categories):
            return 'dominant'
        elif all(cat == 'balanced' for cat in categories):
            
            return 'balanced'
        else: 
            return 'non-dominant'
        
    gene_cases = df.groupby('Gene stable ID')['isoform_category'].apply(categorize_gene_cases).reset_index()
    gene_cases.rename(columns={'isoform_category' : 'gene_category'}, inplace=True)
    return gene_cases



def plausibility_filtering_gene_unique(gene_unique, gene_categories):

    n_before = gene_unique.shape[0]
    gene_unique = gene_unique[~gene_unique.source_gene.isin(gene_categories[gene_categories.gene_category.isin(['single'])]['Gene stable ID'])]
    n_after_single = gene_unique.shape[0]
    gene_unique = gene_unique[~gene_unique.source_gene.isin(gene_categories[gene_categories.gene_category.isin(['dominant'])]['Gene stable ID'])]
    n_after_dominant = gene_unique.shape[0]
    filter_info = {'n_before': n_before, 'n_after_equivalence': n_after_single,  'n_after_dominance': n_after_dominant}
    gene_unique = gene_unique.sort_values('median_importance', ascending=False)
    return gene_unique, filter_info



def create_common_edge_dataframe(common_edges):
    common_edges_ge = common_edges[pd.isna(common_edges.source_transcript)].copy()
    common_edges_t = common_edges[~pd.isna(common_edges.source_transcript)].copy()
    merged_edges = common_edges_t.merge(common_edges_ge, on='edge_key', suffixes=['_te', '_ge'])
    merged_edges = merged_edges.rename(columns={'source_transcript_te': 'source_transcript', 'source_gene_te': 'source_gene', 'target_te': 'target_gene'})
    merged_edges = merged_edges.loc[:, ['source_transcript', 'source_gene', 'target_gene', 'edge_key', 'frequency_te', 'mean_importance_te', 'median_importance_te',  'frequency_ge', 'mean_importance_ge', 'median_importance_ge']]
    return merged_edges


def frequency_filtering_common_edges_dominant(consistent, threshold_frequency = 10):
    """
    Plausbility filtering for edges which are found in both networks
    As the isoforms are either dominant or there is a single isoform it is required
        - the are at the same high frequency in both networks
        - fold change of importances is about 1.0

    """
    n_before = consistent.shape[0]
    consistent = consistent[(consistent['frequency_te']>=threshold_frequency )& (consistent['frequency_ge']>=threshold_frequency)]
    n_after_frequency = consistent.shape[0]
    
    filter_info = {
        'before': n_before, 
        'n_after_frequency': n_after_frequency, 
        }

    return consistent.copy(), filter_info


def plausibility_filtering_common_edges_dominant(consistent, threshold_upper=1.5, threshold_lower=0.5):
    """
    Plausbility filtering for edges which are found in both networks
    As the isoforms are either dominant or there is a single isoform it is required
        - the are at the same high frequency in both networks
        - fold change of importances is about 1.0

    """
    n_before = consistent.shape[0]
    consistent['fc'] = consistent['mean_importance_te']/consistent['mean_importance_ge']
    consistent['mean_importance_te_ge'] = consistent[['median_importance_te', 'median_importance_ge']].mean(axis=1)

    consistent = consistent[(consistent['fc']<threshold_upper) & (consistent['fc']>threshold_lower)].sort_values('mean_importance_te', ascending=False)
    n_after_importance = consistent.shape[0]

    consistent = consistent.reset_index()
    consistent = consistent.drop(columns = ['index'])
    consistent = consistent.sort_values('mean_importance_te_ge',ascending=False)
    
    filter_info = {
        'before': n_before, 
        'n_after_importance': n_after_importance
        }

    return consistent.copy(), filter_info


def split_by_isoform_category(merged_edges, gene_categories):
    merged_edges = merged_edges.merge(gene_categories, left_on='source_gene', right_on='Gene stable ID')
    merged_edges = merged_edges.loc[:, ['source_transcript', 'source_gene', 'target_gene', 'edge_key', 'frequency_te', 'mean_importance_te', 'median_importance_te',  'frequency_ge', 'mean_importance_ge', 'median_importance_ge', 'gene_category']]

    consistent = merged_edges[merged_edges.gene_category.isin(['single', 'dominant'])].copy()
    ambigous = merged_edges[merged_edges.gene_category.isin(['balanced', 'non-dominant'])].copy()
    return consistent, ambigous


def find_likely_isoform_specific(ambigous, lf_threshold = 2, frequency_threshold = 10):
    ambigous['fc']  = ambigous['mean_importance_te']/ambigous['mean_importance_ge']
    
    n_ambigous_before = ambigous.shape[0]
    likely_isoform_specific = ambigous[(ambigous['fc']>lf_threshold)  ]
    n_isoform_likely_before_frequency_filtering = likely_isoform_specific.shape[0]
    likely_isoform_specific = likely_isoform_specific[(likely_isoform_specific['frequency_te']>=frequency_threshold)].copy()
    n_isoform_likely_after_frequency_filtering = likely_isoform_specific.shape[0]
    remove_keys = likely_isoform_specific.edge_key
    likely_isoform_specific = likely_isoform_specific.sort_values('median_importance_te', ascending=False)
    likely_isoform_specific = likely_isoform_specific.reset_index()
    likely_isoform_specific = likely_isoform_specific.drop(columns = ['index'])
    ambigous = ambigous[~ambigous.edge_key.isin(remove_keys)].copy()

    filter_info = {
        'n_ambigous_before': n_ambigous_before,  
        'n_isoform_likely_before_frequency_filtering': n_isoform_likely_before_frequency_filtering, 
        'n_after_frequency': n_isoform_likely_after_frequency_filtering,
        'n_ambigous_after': ambigous.shape[0]}

    return likely_isoform_specific, ambigous, filter_info


def find_likely_gene_specific(ambigous, frequency_threshold = 10, lf_threshold = 0.5):
    ambigous['fc']  = ambigous['mean_importance_te']/ambigous['mean_importance_ge']

    n_ambigous_before = ambigous.shape[0]

    ambigous['gene_specific'] = ambigous.groupby('edge_key')['fc'].transform('max') < lf_threshold
    
    likely_gene_specific = ambigous[ambigous.gene_specific]
    n_gene_likely_before_frequency_filtering = likely_gene_specific.shape[0]
    remove_keys_edge  = likely_gene_specific.edge_key
    likely_gene_specific = likely_gene_specific[(likely_gene_specific['frequency_ge']>=frequency_threshold)].copy()
    n_gene_likely_after_frequency_filtering = likely_gene_specific.shape[0]

    likely_gene_specific.sort_values('median_importance_ge')
    likely_gene_specific = likely_gene_specific.reset_index()
    likely_gene_specific = likely_gene_specific.drop(columns = ['index'])

    n_gene_likely = likely_gene_specific.shape[0]

    ambigous = ambigous[~ambigous.edge_key.isin(remove_keys_edge)].copy()
    
    filter_info = {
        'n_ambigous_before': n_ambigous_before,  
        'n_gene_likely_before_frequency_filtering': n_gene_likely_before_frequency_filtering, 
        'n_after_frequency': n_gene_likely_after_frequency_filtering,
        'n_ambigous_after': ambigous.shape[0]}


    return likely_gene_specific, ambigous, filter_info
    