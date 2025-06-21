import pandas as pd 
from collections import defaultdict
import numpy as np



def create_transcipt_annotation_database(tf_list, appris_path, digger_path):
    '''
    Creates an annotation database for transcription factor (TF) isoforms by integrating data from 
    APPRIS and DIGGER sources.

    Parameters:

        tf_list : pd.DataFrame
            A DataFrame containing transcript-level information for transcription factors, 
            including a 'Transcript stable ID' column. May include an 'index' column that will be dropped.

        appris_path : str
            Path to a tab-separated file containing APPRIS annotations for transcripts.

        digger_path : str
            Path to a CSV file containing DIGGER annotations including domain mappings to transcripts.

    Returns
   
        pd.DataFrame
            A merged DataFrame that includes original TF transcript information along with additional 
            functional annotations from APPRIS and DIGGER, aggregated at the transcript level.
    '''
    
    appris_df = pd.read_csv(appris_path, sep='\t')
    digger = pd.read_csv(digger_path)


    # preprocessing
    tf_database_og = tf_list.drop(columns=['index'])

    # merge appris
    tf_database = tf_database_og.merge(appris_df, left_on='Transcript stable ID', right_on='Transcript ID', how='left')
    # drop unrelevant columns
    tf_database = tf_database.drop(columns=['Ensembl Gene ID', 'Transcript ID'])

    #preprocess digger
    digger = digger.drop(columns=['CDS start', 'CDS end', 'Pfam start','Pfam end', 'Genomic coding start', 'Genomic coding end', 'Strand', 'Chromosome/scaffold name', 'Strand'])
    digger = digger.groupby('Transcript stable ID', as_index=False).agg(lambda x: list(x.dropna()))

    #merge digger
    tf_database = tf_database.merge(digger, on='Transcript stable ID', how='left')

    return tf_database


def get_unique_items(transcript_items, related_df, column_name):
    '''
    Find items that are present only in the transcript and not in related isoforms.

    Parameters:

        transcript_items : list
            List of items from the transcript (e.g., exon IDs, Pfam IDs).

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
                            ('Pfam ID', 'unique Pfam ID')]:
        transcript_data[unique_key] = get_unique_items(transcript_data.get(col), related_transcripts, col)
    return transcript_data


def check_annotations(transcript_id, annotation_database):

    '''
    Get the available annotations of the transcript and compile a unique list of exon IDs and Pfam IDs.

    Parameters:

        transcript_id : str
            Ensembl ID of the transcript.

        annotation_database : pd.DataFrame
            DataFrame containing APPRIS, DIGGER, and Ensembl annotation data.

    Returns:
    
        dict
            A dictionary containing all annotations of the transcript, including unique exon and Pfam IDs.
  
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
    related_transcripts = annotation_database[
        (annotation_database['Gene stable ID'] == gene_id) &
        (annotation_database['Transcript stable ID'] != transcript_id)
    ]

    #Check Appris Annotation
    if pd.notna(t_row['APPRIS Annotation']) and 'PRINCIPAL' not in t_row['APPRIS Annotation']:
        principal_isos = related_transcripts[
            related_transcripts['APPRIS Annotation'].str.contains('PRINCIPAL', na=False)
        ]
        comparison_df = principal_isos if not principal_isos.empty else related_transcripts
    else:
        comparison_df = related_transcripts

    transcript_data = {
        'Protein Coding': True,
        'Transcript type': t_row.get('Transcript type'),
        'APPRIS Annotation': t_row.get('APPRIS Annotation'),
        'Trifid Score': t_row.get('Trifid Score'),
        'Exon stable ID': t_row.get('Exon stable ID'),
        'Pfam ID': t_row.get('Pfam ID'),
    }

    return pd.Series(compare_values(transcript_data, comparison_df))



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


def merge_annotations_to_grn(grn, annotation_database):
    '''
    Merge transcript-level annotations into a gene regulatory network (GRN) based on source transcript IDs.

    Parameters:

        grn : pd.DataFrame
            DataFrame representing the gene regulatory network, containing a 'source' column with transcript IDs.

        annotation_database : pd.DataFrame
            DataFrame containing APPRIS, DIGGER, and Ensembl annotation data.

    Returns:
    
        pd.DataFrame
            GRN DataFrame with additional columns from the annotation database merged on the 'source' transcript ID.

    '''

    unique_transcripts = grn['source'].unique()
    annot_df = build_transcript_annotation_table_for_unique_tfs(unique_transcripts, annotation_database)
    grn_annot = grn.merge(annot_df, how='left', left_on='source', right_index=True)

    return grn_annot