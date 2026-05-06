import os, sys, pandas as pd

from ..rsrc import bioinfo_tabs, prog_tabs, dir_tabs
from . import fa_ops


# BASIC FUNCTION #
# USAGE:    load the cluster table file to a pandas dataframe.
# INPUT:    clust_tab_file: the path of the cluster table file.
#           tab_type: the type of the cluster table file, default is 'orthof_g'.
# OUTPUT:   a pandas dataframe of the cluster table file.
def load_clust_tab_file(clust_tab_file, tab_type='orthof_g'):
    ## set the column names, data types and columns to use according to the tab type.
    if tab_type == 'orthof_g':
        clust_tab_df = pd.read_csv(clust_tab_file, sep='\t', dtype=str)
    elif tab_type == 'orthof_c':
        clust_tab_df = pd.read_csv(clust_tab_file, sep='\t', dtype={0: str})
    else:
        raise ValueError('The tab type is not supported.')

    ## return the cluster table dataframe.
    return clust_tab_df













#def load_standard_cluster_file(standard_cluster_file):
    #data_type_dict = {
        #'Genome':str, 'Orthogroup':str, 'Gene_id':str
    #}
    #col_list = [0,1,2]
    #standard_cluster_file_df = pd.read_csv(standard_cluster_file, sep='\t', header=0, dtype=data_type_dict, usecols=col_list)
    #return standard_cluster_file_df

#def load_mmseqs_cluster_file(mmseqs_cluster_file):
    #name_list = ['cluster_id', 'seq_id']
    #data_type_dict = {
        #'cluster_id':str, 'seq_id':str
    #}
    #col_list = [0,1]
    #mmseqs_cluster_file_df = pd.read_csv(mmseqs_cluster_file, sep='\t', header=None, names=name_list, dtype=data_type_dict, usecols=col_list)
    #mmseqs_cluster_file_groupby = mmseqs_cluster_file_df.groupby('cluster_id')
    #return mmseqs_cluster_file_groupby

#def load_orthomcl_cluster_file(orthomcl_cluster_file):
    #orthomcl_cluster_file_df = pd.read_csv(orthomcl_cluster_file, sep=' ', header=None, dtype=str)
    #return orthomcl_cluster_file_df

#def load_orthofinder_orthogroups_tsv_file(orthogroups_tsv_file):
    #orthogroups_tsv_df = pd.read_csv(orthogroups_tsv_file, sep='\t', header=0, dtype=str)
    #return orthogroups_tsv_df

#def load_orthofinder_genecount_tsv_file(genecount_tsv_file):
    #genecount_tsv_df = pd.read_csv(genecount_tsv_file, sep='\t', header=0, dtype=str)
    #return genecount_tsv_df

#def modify_mmseqs_tsv_df_2_wide(mmseqs_cluster_file_df, genome_2_seqid_tsv_df):
    #modified_mmseqs_cluster_file_df = mmseqs_cluster_file_df.merge(genome_2_seqid_tsv_df, on='seq_id', how='left')
    #modified_wide_mmseqs_cluster_file_df = modified_mmseqs_cluster_file_df.pivot_table(index='cluster_id', columns='genome_id', values='seq_id', aggfunc=lambda seq_id_sers: ', '.join(seq_id_sers))
    #modified_wide_mmseqs_cluster_file_df = modified_wide_mmseqs_cluster_file_df.reset_index()
    #return modified_wide_mmseqs_cluster_file_df

#def parse_orthogroups_tsv_df_2_long(orthogroups_tsv_df):
    #melted_long_orthogroups_tsv_df = orthogroups_tsv_df.melt(id_vars='Orthogroup', var_name='Genome', value_name='Gene_ids')
    #splitted_melted_long_orthogroups_tsv_df = melted_long_orthogroups_tsv_df.assign(Gene_ids=melted_long_orthogroups_tsv_df['Gene_ids'].str.split(', ')).explode('Gene_ids')
    #splitted_melted_long_orthogroups_tsv_df = splitted_melted_long_orthogroups_tsv_df.rename(columns={'Orthogroup':'Orthogroup', 'Genome':'Genome', 'Gene_ids':'Gene_id'})
    #splitted_melted_long_orthogroups_tsv_df = splitted_melted_long_orthogroups_tsv_df[['Genome', 'Orthogroup', 'Gene_id']]
    #return splitted_melted_long_orthogroups_tsv_df

#def generate_genome_2_seqid_tsv_df(genome_2_seq_fasta_file_list):
    # parse genome to seq fasta file list to dict
    #genome_2_seq_fasta_file_dict = {pair.split(',')[0]: pair.split(',')[1] for pair in genome_2_seq_fasta_file_list}
    # initialize dataframe
    #genome_2_seq_fasta_file_content_id_df_df = pd.DataFrame(columns=['genome_id', 'fasta_file', 'fasta_content', 'id_df'])
    #genome_2_seqid_tsv_df = pd.DataFrame(columns=['genome_id', 'seq_id'])
    # parse genome to seq fasta file dict to dataframe
    #genome_2_seq_fasta_file_content_id_df_df['genome_id'] = list(genome_2_seq_fasta_file_dict.keys())
    #genome_2_seq_fasta_file_content_id_df_df['fasta_file'] = list(genome_2_seq_fasta_file_dict.values())
    # load seq fasta file content, parse all seq id of each genome to a single dataframe in corresponding row
    #genome_2_seq_fasta_file_content_id_df_df['fasta_content'] = genome_2_seq_fasta_file_content_id_df_df['fasta_file'].apply(fasta_operator.load_fasta_file_content)
    #genome_2_seq_fasta_file_content_id_df_df['id_df'] = genome_2_seq_fasta_file_content_id_df_df.apply(lambda row: fasta_operator.load_fasta_file_header(row['seq_fasta_content'], mode='id'), axis=1)
    # generate genome to seq id tsv df
    #for _, row in genome_2_seq_fasta_file_content_id_df_df.iterrows():
        #genome_id = row['genome_id']
        #id_df = row['id_df']
        #id_df['genome_id'] = genome_id
        #genome_2_seqid_tsv_df = pd.concat([genome_2_seqid_tsv_df, id_df], ignore_index=True)
    #return genome_2_seqid_tsv_df
