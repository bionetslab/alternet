import pandas as pd
from alternet.edge_categorization import *

class Alternet:
    __slots__ = [
        'canonical', 
        'as_source', 
        'as_full', 
        'transcript_data', 
        'transcript_col',
        'gene_col',
        'gene_data', 
        'usage_df', 
        'reliability_df', 
        'sample_cols', 
        'set_a', 
        'set_b', 
        'set_c', 
        'set_d',
        'set_d_full',
        'MIN_FREQUENCY',
        'IMPORTANCE_PERCENTILE',
        'regulator_list',
        'gene2tx',
        'tx2gene'
    ]

    def __init__(
        self, 
        canonical: pd.DataFrame, 
        as_source: pd.DataFrame, 
        as_full: pd.DataFrame, 
        transcript_data: pd.DataFrame,
        regulator_list: pd.DataFrame,
        gene_col: str ='gene_id', 
        transcript_col: str='transcript_id'
        
    ):
        self.canonical = canonical
        self.as_source = as_source
        self.as_full = as_full
        self.transcript_data = transcript_data
        self.regulator_list = regulator_list
        self.transcript_col = transcript_col
        self.gene_col = gene_col

        self.MIN_FREQUENCY = 10
        self.IMPORTANCE_PERCENTILE = 0.7

        self.gene_data, self.sample_cols = self._compute_gene_counts()
        self.gene2tx = self._gene2tx()
        self.tx2gene = self._tx2gene()
        self.usage_df, self.reliability_df = self._calculate_transcript_usage()
        self.canonical, self. as_source, self.as_full = self._filter_edges()
        self.set_a = self._compute_set_a()
        self.set_d_full = self._compute_set_d_full()
        self.set_d = self._compute_set_d()
        self.set_b = self._compute_set_b()
        self.set_c = self._compute_set_c()







    def _compute_gene_counts(self):
        sample_cols = [c for c in self.transcript_data.columns if c not in [self.transcript_col, self.gene_col]]
        gene_tpm = self.transcript_data.groupby('gene_id')[sample_cols].sum().reset_index()
        return gene_tpm, sample_cols

        
    def _gene2tx(self):
        gene_to_transcripts = self.transcript_data.groupby(self.gene_col)[self.transcript_col].apply(list).to_dict()
        return gene_to_transcripts

    def _tx2gene(self):
        tx2gene =  dict(zip(self.transcript_data[self.transcript_col], self.transcript_data[self.gene_col]))
        return tx2gene


    def _filter_edges(self):
        canonical_grn_filtered = filter_edges(self.canonical, self.MIN_FREQUENCY, self.IMPORTANCE_PERCENTILE)
        as_source_grn_filtered = filter_edges(self.as_source, self.MIN_FREQUENCY, self.IMPORTANCE_PERCENTILE)
        fully_as_grn_filtered = filter_edges(self.as_full, self.MIN_FREQUENCY, self.IMPORTANCE_PERCENTILE)
        return canonical_grn_filtered, as_source_grn_filtered, fully_as_grn_filtered


    def _calculate_transcript_usage(self, min_gene_tpm=1.0, epsilon=1e-8):
        """
        Calculates transcript usage (fraction of gene expression) and reliability masks.
        
        Returns:
            usage_df: DataFrame of transcript usage (0.0 to 1.0)
            reliability_df: Binary mask (1.0 if gene expression >= min_gene_tpm)
            metadata: Dictionary containing gene-to-transcript mapping and mean expressions
        """
        genes = self.transcript_data[self.gene_col].values

        tpm_numeric = self.transcript_data[self.sample_cols].astype(float)
        gene_numeric_total  = self.transcript_data[[self.gene_col]+self.sample_cols].groupby(self.gene_col).transform('sum').astype(float)

        usage_values = tpm_numeric / (gene_numeric_total + epsilon)
        usage_df = pd.concat([self.transcript_data[[self.transcript_col, self.gene_col]], usage_values], axis=1)

        reliability_df = (gene_numeric_total >= min_gene_tpm).astype(float)
        reliability_df.index = self.transcript_data[self.transcript_col].values

        return usage_df, reliability_df


    def _compute_set_a(self):
        print('Computing set A')
        net1_tf = self.canonical[self.canonical['reg_type'].isin(['TF', 'TF_SF'])].copy()
        net2_tf = self.as_source[self.as_source['reg_type'].isin(['TF', 'TF_SF'])].copy()
        set_a = canonical_vs_source_as(net1_tf,net2_tf)
        return set_a

    def _compute_set_d_full(self):
        print('Set D diambiguation')
        # Add information about TF likelyhood to TFs
        net3_tfsf = self.as_full[self.as_full['reg_type'] == 'TF_SF'].copy()
        t_temps = self.transcript_data.drop(columns = {self.gene_col}).set_index(self.transcript_col)
        usage_temp = self.usage_df.drop(columns = {self.gene_col}).set_index(self.transcript_col)
        set_d_full = tf_sf_disambigouation_fully_as_aware(net3_tfsf, self.regulator_list, t_temps , usage_temp, self.reliability_df)

        return set_d_full

    def _compute_set_d(self):
        print('Computing set D')
        tfsf_joint = self.set_d_full[self.set_d_full['tfsf_category'] == 'tfsf_joint'].copy()
        tfsf_ambiguous = self.set_d_full[self.set_d_full['tfsf_category'] == 'tfsf_ambiguous'].copy()
        set_d = pd.concat([tfsf_joint, tfsf_ambiguous], ignore_index=True)
        # Sort by median importance
        set_d = set_d.sort_values('median_importance', ascending=False)
        return set_d

    def _compute_set_b(self):
        print('Computeing set B')
        # Net2: ALL TF-like edges (TF + TF_SF) - they all act as TFs in gene-level target network
        net2_for_t2 = self.as_source[self.as_source['reg_type'].isin(['TF', 'TF_SF'])].copy()
        # Net3: TF only (TF_SF is handled via Set D)
        net3_tf_only = self.as_full[self.as_full['reg_type'] == 'TF'].copy()
        # Add the ones from set d which are tf edges
        tfsf_tf_like = self.set_d_full[self.set_d_full['tfsf_category'] == 'tfsf_tf_like'].copy()
        
        net3_tfsf = self.as_full[self.as_full['reg_type'] == 'TF_SF'].copy()
        set_b = as_source_vs_as_full(net2_for_t2, net3_tf_only, tfsf_tf_like, net3_tfsf, self.tx2gene)
        print(set_b)
        set_b_unpacked = unpack_set_b(net3_tf_only, tfsf_tf_like, net3_tfsf, set_b)
        print(set_b_unpacked)
        set_b_unpacked = annotate_set_b(set_b,set_b_unpacked)
        print(set_b_unpacked)
        return set_b_unpacked

    def _compute_set_c(self):
        net3_tfsf = self.as_full[self.as_full['reg_type'] == 'TF_SF'].copy()
        tfsf_sf_like = self.set_d_full[self.set_d_full['tfsf_category'] == 'tfsf_sf_like'].copy()
        net3_sf = self.as_full[self.as_full['reg_type'] == 'SF'].copy()
        udf = self.usage_df.drop(columns = {'gene_id'}).set_index("transcript_id")
        t_temps = self.transcript_data.drop(columns = {self.gene_col}).set_index(self.transcript_col)

        set_c = compute_set_c(net3_sf, t_temps, self.gene2tx, udf, self.reliability_df, tfsf_sf_like, self.sample_cols, net3_tfsf)
        set_c_unpacked = unpack_set_c(net3_sf, tfsf_sf_like, net3_tfsf, set_c)
        set_c = annotate_set_c(set_c, set_c_unpacked)
        return set_c
        