import pandas as pd 
from collections import defaultdict
import numpy as np

GENE_STABLE_ID = 'Gene stable ID'
TRANSCRIPT_STABLE_ID = 'Transcript stable ID'

def map_tf_ids(tf_list, biomart):
    ''' 
    Reads a transcription factor list from file and merges it with gene information from a BioMart DataFrame.

    Parameters:
        tf_list: (pd.DataFrame): dataFrame with a single column containing the gene name
        biomart (pd.DataFrame): DataFrame with gene and transcript information from BioMart.

    Returns:
        pd.DataFrame: Transcription factors with corresponding gene and transcript information from BioMart.
   '''
    tf_list.columns  = ['TF'] # rename column to TF 
    tf_list = tf_list.merge(biomart, left_on = 'TF', right_on = 'Gene name') # merge tf list with biomart: join 'TF' with 'Gene name'
    tf_list = tf_list.loc[:, ['TF', GENE_STABLE_ID, TRANSCRIPT_STABLE_ID]].drop_duplicates() #remove individual versions
    return tf_list


def map_sf_ids(sf_list_raw: pd.DataFrame, biomart: pd.DataFrame) -> pd.DataFrame:
    """Map SF gene names to Ensembl gene and transcript IDs.

    Args:
        sf_list_raw: DataFrame with a single column of SF gene names.
        biomart: BioMart reference DataFrame.

    Returns:
        DataFrame with columns ``SF``, ``Gene stable ID``, ``Transcript stable ID``.
    """
    sf_list_raw = sf_list_raw.copy()
    sf_list_raw.columns = ["SF"]
    sf_list = sf_list_raw.merge(biomart, left_on="SF", right_on="Gene name")
    sf_list = sf_list.loc[
        :, ["SF", "Gene stable ID", "Transcript stable ID"]
    ].drop_duplicates()
    return sf_list



def combine_tf_sf_lists(
    tf_list: pd.DataFrame, sf_list: pd.DataFrame
) -> pd.DataFrame:
    """Combine TF and SF lists, marking overlapping genes as TF_SF.

    Args:
        tf_list: TF DataFrame with columns ``TF``, ``Gene stable ID``,
            ``Transcript stable ID``.
        sf_list: SF DataFrame with columns ``SF``, ``Gene stable ID``,
            ``Transcript stable ID``.

    Returns:
        Combined DataFrame with columns ``Regulator_name``, ``Gene stable ID``,
        ``Transcript stable ID``, ``Regulator_type``.
    """
    tf_list = tf_list.copy()
    sf_list = sf_list.copy()

    tf_list["Regulator_type"] = "TF"
    sf_list["Regulator_type"] = "SF"
    tf_list = tf_list.rename(columns={"TF": "Regulator_name"})
    sf_list = sf_list.rename(columns={"SF": "Regulator_name"})

    tf_genes = set(tf_list["Gene stable ID"])
    sf_genes = set(sf_list["Gene stable ID"])
    overlap_genes = tf_genes & sf_genes

    combined = pd.concat([tf_list, sf_list], ignore_index=True)
    combined.loc[
        combined["Gene stable ID"].isin(overlap_genes), "Regulator_type"
    ] = "TF_SF"
    combined = combined.drop_duplicates(subset=["Transcript stable ID"], keep="first")
    return combined



def create_transcript_mapping(biomart):
    transcript_mapper = dict(zip(biomart['Transcript stable ID'], biomart['Gene stable ID']))
    return transcript_mapper
    





def create_filtered_gene_to_transcripts_mapping(biomart, gene_list, transcript_list):
    """
    Creates a dictionary mapping 'Gene stable ID' to a list of associated 
    'Transcript stable ID's, filtering the pairs against a provided 
    list of valid genes and a list of valid transcripts.

    Args:
        biomart (pd.DataFrame): A pandas DataFrame with 'Gene stable ID' 
                                and 'Transcript stable ID' columns.
        gene_list (list or set): A collection of gene IDs to include.
        transcript_list (list or set): A collection of transcript IDs to include.

    Returns:
        dict: A dictionary where keys are filtered gene IDs (str) and values are 
              lists of filtered transcript IDs (list of str).
    """
    # Convert lists to sets for O(1) average time complexity lookups
    valid_genes = set(gene_list)
    valid_transcripts = set(transcript_list)
    
    gene_to_transcripts = {}
    
    for gene_id, transcript_id in zip(biomart['Gene stable ID'], biomart['Transcript stable ID']):
        # Apply the exclusion criteria: Both gene and transcript must be in their respective lists
        if gene_id in valid_genes and transcript_id in valid_transcripts:
            if gene_id not in gene_to_transcripts:
                gene_to_transcripts[gene_id] = []
            gene_to_transcripts[gene_id].append(transcript_id)
            
    return gene_to_transcripts
    

def create_transcipt_annotation_database(tf_list, appris_df, digger):
    '''
    Creates an annotation database for transcription factor (TF) isoforms by integrating data from 
    APPRIS and DIGGER sources.

    Parameters:

        tf_list : pd.DataFrame
            A DataFrame containing transcript-level information for transcription factors, 
            including a 'Transcript stable ID' column. May include an 'index' column that will be dropped.

        appris_path : str
            Path to file containing APPRIS annotations for transcripts.

        digger_path : str
            Path to file containing DIGGER annotations including domain mappings to transcripts.

    Returns
   
        pd.DataFrame
            A merged DataFrame that includes original TF transcript information along with additional 
            functional annotations from APPRIS and DIGGER, aggregated at the transcript level.
    '''



    # preprocessing

    # merge appris
    tf_database = tf_list.merge(appris_df, left_on='Transcript stable ID', right_on='Transcript ID', how='left')
    # drop unrelevant columns
    tf_database = tf_database.drop(columns=['Ensembl Gene ID', 'Transcript ID'])

    #preprocess digger
    digger = digger.drop(columns=['CDS start', 'CDS end', 'Pfam start','Pfam end', 'Genomic coding start', 'Genomic coding end', 'Strand', 'Chromosome/scaffold name', 'Strand'])
    digger = digger.groupby('Transcript stable ID', as_index=False).agg(lambda x: list(x.dropna()))

    #merge digger
    tf_database = tf_database.merge(digger, on='Transcript stable ID', how='left')

    return tf_database




## ANNOTATION FOR ISOFORM EXLUSIVE EDGES

def get_unique_items(transcript_items, related_df, column_name):
    '''
    Find items that are present only in specific transcript and not in related isoforms.

    Parameters:

        transcript_items : List
            List of characteristics from the transcript (e.g., exon IDs, Pfam IDs).

        related_df : pd.DataFrame
            DataFrame of related transcripts.

        column_name : str
            The name of the column in related_df containing the comparable lists.

    Returns:
    
        list
            Items that are unique to the input transcript.
    '''
    if not isinstance(transcript_items, list):
        return []
    
    tr_set = set(transcript_items)

    if column_name in related_df.columns:
        other_items = related_df[column_name].explode()
        other_items = set(other_items.dropna()) # avoid nans
        
        return list(tr_set - other_items)

    return list(tr_set)


def get_missing_items(transcript_items, related_df, column_name):
    '''
    Find items that are present only in specific transcript and not in related isoforms.

    Parameters:

        transcript_items : List
            List of characteristics from the transcript (e.g., exon IDs, Pfam IDs).

        related_df : pd.DataFrame
            DataFrame of related transcripts.

        column_name : str
            The name of the column in related_df containing the comparable lists.

    Returns:
    
        list
            Items that are unique to the input transcript.
    '''
    if not isinstance(transcript_items, list):
        return []
    
    tr_set = set(transcript_items)

    if column_name in related_df.columns:
        other_items = related_df[column_name].explode()
        other_items = set(other_items.dropna()) # avoid nans
        return list(other_items - tr_set)

    return list(tr_set)



def compare_values(transcript_data, related_transcripts):
    '''
    Compares 'Exon stable ID' and 'Pfam ID' from a transcript dictionary and a DataFrame
    of related isoforms, and returns a dictionary with unique values for each of these columns.

    Parameters:

        transcript_dict : dict
            Dictionary containing 'Exon stable ID' and 'Pfam ID' for a single transcript.

        isoforms_df : pd.DataFrame
            DataFrame containing related isoforms with 'Exon stable ID' and 'Pfam ID' columns.

    Returns:
    
        dict
            A dictionary with keys 'unique_exons' and 'unique_pfam', containing lists of unique values.

    
    '''


    for col, unique_key in [('Exon stable ID', 'unique Exon stable ID'),
                            ('Pfam ID', 'unique Pfam ID'),
                            ('Exon stable ID', 'missing Exon stable ID'),
                            ('Pfam ID', 'missing Pfam ID')]:
        transcript_data[unique_key] = get_unique_items(transcript_data.get(col), related_transcripts, col)

    for col, unique_key in [('Exon stable ID', 'missing Exon stable ID'),
                            ('Pfam ID', 'missing Pfam ID')]:
        transcript_data[unique_key] = get_missing_items(transcript_data.get(col), related_transcripts, col)
    return transcript_data


def check_annotations(transcript_id, annotation_database):

    '''
    Get the available annotations of the transcript and compile a unique list of exon IDs and Pfam IDs.

    Parameters:

        transcript_id : str
            Ensembl ID of the transcript.

        annotation_database : pd.DataFrame
            DataFrame containing APPRIS, and DIGGER annotation data.

    Returns:
    
        dict
            A dictionary containing all annotations of the transcript, including unique exon IDs and Pfam IDs.
  
    '''

    transcript = annotation_database[annotation_database['Transcript stable ID'] == transcript_id]
    
    if transcript.empty:
        # return same structure with None values if transcript not found
        return pd.Series({
            'Protein Coding': False,
            'Transcript type': None,
            'APPRIS Annotation': None,
            'Exon stable ID': None,
            'unique Exon stable ID': None,
            'Pfam ID': None,
            'unique Pfam ID': None
        })
    
    t_row = transcript.iloc[0]

    if t_row['Transcript type'] != 'protein_coding':
        # no further investigation because not plausible
        return pd.Series({
            'Protein Coding': False,
            'Transcript type': t_row.get('Transcript type'),
            'APPRIS Annotation': t_row.get('APPRIS Annotation'),
            'Trifid Score': t_row.get('Trifid Score'),
            'Exon stable ID': t_row.get('Exon stable ID'),
            'unique Exon stable ID': None,
            'Pfam ID': t_row.get('Pfam ID'),
            'unique Pfam ID': None
        })
    
    gene_id = t_row['Gene stable ID']
    # get other transcripts from the same gene
    related_transcripts = annotation_database[(annotation_database['Gene stable ID'] == gene_id) & (annotation_database['Transcript stable ID'] != transcript_id)]

    #Check Appris Annotation
    if pd.notna(t_row['APPRIS Annotation']) and ('PRINCIPAL' not in t_row['APPRIS Annotation']):
        # if the current transcript is not the principal isoform, get the annotaiton for the principal isoform to compare
        # against, but if it does not exist, use all related transcripts.
        principal_isos = related_transcripts[related_transcripts['APPRIS Annotation'].str.contains('PRINCIPAL', na=False)]
        comparison_df = principal_isos if not principal_isos.empty else related_transcripts
    else:
        # if the current transcript is annotated as the principal isoform, compare against all other transcripts
        comparison_df = related_transcripts

    transcript_data = {
        'Protein Coding': True,
        'Transcript type': t_row.get('Transcript type'),
        'APPRIS Annotation': t_row.get('APPRIS Annotation'),
        'Trifid Score': t_row.get('Trifid Score'),
        'Exon stable ID': t_row.get('Exon stable ID'),
        'Pfam ID': t_row.get('Pfam ID'),
    }

    # for compare the transcript data against the chosen set to find out what makes the transcript unique.
    unique_items = compare_values(transcript_data, comparison_df)
    return pd.Series(unique_items)



def build_transcript_annotation_table_for_unique_tfs(unique_tfs, annotation_database):
    '''
    Build a transcript-level annotation table for a set of unique transcription factor (TF) isoforms.

    Parameters:

        unique_tfs : list
            List of Ensembl transcript IDs corresponding to unique TF isoforms.

        annotation_database : pd.DataFrame
            DataFrame containing APPRIS, DIGGER, and Ensembl annotation data.

    Returns:
    
        pd.DataFrame
            DataFrame indexed by transcript ID, containing annotation information for each TF isoform.

    
    '''
    # Precompute annotations for all transcripts
    annotations = []
    for tid in unique_tfs:

        annot = check_annotations(tid, annotation_database)
        annot['Transcript stable ID'] = tid
        annotations.append(annot)

    annotation_df = pd.DataFrame(annotations).set_index('Transcript stable ID')
    return annotation_df




def annotate_isoform_exclusive_edges(grn, annotation_database, transcript_column = 'source_transcript', suffixes = None):
    '''
    Merge transcript-level annotations into a gene regulatory network (GRN) based on source transcript IDs.

    Parameters:

        grn : pd.DataFrame
            DataFrame representing the gene regulatory network, containing a 'source' column with transcript IDs.

        annotation_database : pd.DataFrame
            DataFrame containing APPRIS, and DIGGER annotation data.

    Returns:
    
        pd.DataFrame
            GRN DataFrame with additional columns from the annotation database merged on the 'source' transcript ID.

    '''

    unique_transcripts = grn[transcript_column].unique()
    annot_df = build_transcript_annotation_table_for_unique_tfs(unique_transcripts, annotation_database)
    if suffixes is not None:
        grn_annot = grn.merge(annot_df, how='left', left_on=transcript_column, right_index=True, suffixes =suffixes)
    else:
        grn_annot = grn.merge(annot_df, how='left', left_on=transcript_column, right_index=True)
        

    return grn_annot


#### ANNOTATE CONSITSTENT EDGES 

def get_annotation(transcript_id, annotation_database):

    '''
    Get the available annotations of the transcript and compile a unique list of exon IDs and Pfam IDs.

    Parameters:

        transcript_id : str
            Ensembl ID of the transcript.

        annotation_database : pd.DataFrame
            DataFrame containing APPRIS, and DIGGER annotation data.

    Returns:
    
        dict
            A dictionary containing all annotations of the transcript, including unique exon IDs and Pfam IDs.
  
    '''

    transcript = annotation_database[annotation_database['Transcript stable ID'] == transcript_id]
    
    if transcript.empty:
        # return same structure with None values if transcript not found
        return pd.Series({
            'Protein Coding': False,
            'Transcript type': None,
            'APPRIS Annotation': None,
            'Exon stable ID': None,
            'unique Exon stable ID': None,
            'Pfam ID': None,
            'unique Pfam ID': None
        })
    
    t_row = transcript.iloc[0]


    transcript_data = {
        'Protein Coding': True,
        'Transcript type': t_row.get('Transcript type'),
        'APPRIS Annotation': t_row.get('APPRIS Annotation'),
        'Trifid Score': t_row.get('Trifid Score'),
        'Exon stable ID': t_row.get('Exon stable ID'),
        'Pfam ID': t_row.get('Pfam ID'),
    }

    return transcript_data



def get_transcript_annotation_table_for_unique_tfs(unique_tfs, annotation_database):
    '''
    Build a transcript-level annotation table for a set of unique transcription factor (TF) isoforms.

    Parameters:

        unique_tfs : list
            List of Ensembl transcript IDs corresponding to unique TF isoforms.

        annotation_database : pd.DataFrame
            DataFrame containing APPRIS, DIGGER, and Ensembl annotation data.

    Returns:
    
        pd.DataFrame
            DataFrame indexed by transcript ID, containing annotation information for each TF isoform.

    
    '''
    # Precompute annotations for all transcripts
    annotations = []
    for tid in unique_tfs:

        annot = get_annotation(tid, annotation_database)
        annot['Transcript stable ID'] = tid
        annotations.append(annot)

    annotation_df = pd.DataFrame(annotations).set_index('Transcript stable ID')
    return annotation_df


def annotate_consistent_edges(grn, annotation_database, transcript_column = 'source_transcript'):
    '''
    Merge transcript-level annotations into a gene regulatory network (GRN) based on source transcript IDs.

    Parameters:

        grn : pd.DataFrame
            DataFrame representing the gene regulatory network, containing a 'source' column with transcript IDs.

        annotation_database : pd.DataFrame
            DataFrame containing APPRIS, and DIGGER annotation data.

    Returns:
    
        pd.DataFrame
            GRN DataFrame with additional columns from the annotation database merged on the 'source' transcript ID.

    '''

    unique_transcripts = grn[transcript_column].unique()
    annot_df = get_transcript_annotation_table_for_unique_tfs(unique_transcripts, annotation_database)
    grn_annot = grn.merge(annot_df, how='left', left_on=transcript_column, right_index=True)

    return grn_annot





def get_common_annotation_dataframe(unique_genes, annotation_database):
    '''
    Build a transcript-level annotation table for a set of unique transcription factor (TF) isoforms.

    Parameters:

        unique_tfs : list
            List of Ensembl transcript IDs corresponding to unique TF isoforms.

        annotation_database : pd.DataFrame
            DataFrame containing APPRIS, DIGGER, and Ensembl annotation data.

    Returns:
    
        pd.DataFrame
            DataFrame indexed by transcript ID, containing annotation information for each TF isoform.

    
    '''
    # Precompute annotations for all transcripts
    annotations = []
    for gid in unique_genes.keys():
        annot = get_common_annotations(gid, unique_genes[gid], annotation_database)
        annot['Gene stable ID'] = gid
        annotations.append(annot)
    annotation_df = pd.DataFrame(annotations).set_index('Gene stable ID')
    return annotation_df



def annotate_gene_exclusive_edges(grn, annotation_database, gene_transcript_mapping, gene_column = 'source_gene', transcript_column = 'source_transcript'):
    '''
    Merge transcript-level annotations into a gene regulatory network (GRN) based on source transcript IDs.

    Parameters:

        grn : pd.DataFrame
            DataFrame representing the gene regulatory network, containing a 'source' column with transcript IDs.

        annotation_database : pd.DataFrame
            DataFrame containing APPRIS, and DIGGER annotation data.

    Returns:
    
        pd.DataFrame
            GRN DataFrame with additional columns from the annotation database merged on the 'source' transcript ID.

    '''

    #gene_transcript_mapping = grn.groupby(gene_column)[transcript_column].unique().apply(list).to_dict()
    annot_df = get_common_annotation_dataframe(gene_transcript_mapping, annotation_database)
    grn_annot = grn.merge(annot_df, how='left', left_on=gene_column, right_index=True)
    return grn_annot



def get_intersection(relevant_items):

    relevant_items = [set(tl)  for tl in relevant_items]

    if len(relevant_items)>1:
        common_elements = relevant_items[0].intersection(*relevant_items[1:])
    elif len(relevant_items)>0:
        common_elements = relevant_items[0]
    else:
        common_elements = []
    return list(common_elements)



def get_intersection_wrapper(transcript_data, transcript_info):
    '''
    Compares 'Exon stable ID' and 'Pfam ID' from a transcript dictionary and a DataFrame
    of related isoforms, and returns a dictionary with unique values for each of these columns.

    Parameters:

        transcript_dict : dict
            Dictionary containing 'Exon stable ID' and 'Pfam ID' for a single transcript.

        isoforms_df : pd.DataFrame
            DataFrame containing related isoforms with 'Exon stable ID' and 'Pfam ID' columns.

    Returns:
    
        dict
            A dictionary with keys 'unique_exons' and 'unique_pfam', containing lists of unique values.
    '''


    for col, unique_key in [('Exon stable ID', 'common Exon stable ID'),
                            ('Pfam ID', 'common Pfam ID')]:
        transcript_data[unique_key] = get_intersection(transcript_info.get(col).dropna())

    return transcript_data


def get_common_annotations(gene_id, transcript_ids, annotation_database):

    '''
    Get the available annotations of the transcript and compile a unique list of exon IDs and Pfam IDs.

    Parameters:

        transcript_id : str
            Ensembl ID of the transcript.

        annotation_database : pd.DataFrame
            DataFrame containing APPRIS, and DIGGER annotation data.

    Returns:
    
        dict
            A dictionary containing all annotations of the transcript, including unique exon IDs and Pfam IDs.
  
    '''

    gene = annotation_database[annotation_database['Gene stable ID'] == gene_id]

    if gene.empty:
        # return same structure with None values if transcript not found
        return pd.Series({
            'Protein Coding': False,
            'Transcript type': None,
            'APPRIS Annotation': None,
            'Exon stable ID': None,
            'common Exon stable ID': None,
            'Pfam ID': None,
            'common Pfam ID': None
        })
    

    
    relevant_transcripts = gene[gene['Transcript stable ID'].isin(transcript_ids)]

    if relevant_transcripts.shape[0]==0:
        return pd.Series({
            'Protein Coding': False,
            'Transcript type': None,
            'APPRIS Annotation': None,
            'Exon stable ID': None,
            'common Exon stable ID': None,
            'Pfam ID': None,
            'common Pfam ID': None
        })

    g_row = relevant_transcripts.iloc[0]

    transcript_data = {
        'Protein Coding': True,
        'Transcript type': g_row.get('Transcript type'),
        'APPRIS Annotation': g_row.get('APPRIS Annotation'),
        'Trifid Score': g_row.get('Trifid Score'),
        'Exon stable ID': g_row.get('Exon stable ID'),
        'Pfam ID': g_row.get('Pfam ID'),
    }


    # for compare the transcript data against the chosen set to find out what makes the transcript unique.
    unique_items = get_intersection_wrapper(transcript_data, relevant_transcripts)
    return pd.Series(unique_items)







def compute_isoform_gene_correlations(transcript_data_cp_scaled, gene_data_cp_scaled, gene_to_transcript_mapping):

    correlation_collector = []

    for g in gene_to_transcript_mapping:
        if len(gene_to_transcript_mapping[g]) == 1:
            continue
        subi = pd.concat([transcript_data_cp_scaled.loc[:, gene_to_transcript_mapping[g]].T, gene_data_cp_scaled.loc[:, [g]].T], ignore_index=True)
        corre = np.corrcoef(subi)[0][1:]
        for c in range(len(corre)):
            correlation_collector.append([g, gene_to_transcript_mapping[g][c], corre[c]])
    correlation_collector = pd.DataFrame(correlation_collector)
    correlation_collector.columns = ['gene','transcript', 'correlation']

    return correlation_collector