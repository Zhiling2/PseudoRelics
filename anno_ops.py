#——————————————————————————————————————————————————————————————————————————————————#
# PROJ:     ToolKits                                                               #
# VER:      1.0                                                                    #
# TYPE:     Module                                                                 #
# DESC:     A module for annotation operations.                                    #
# AUTHOR:                                                                  #
# AFFIL:    Wuhan Botanical Garden, Chinese Academy Of Sciences                    #
# E-MAIL:                       #
# DATE:     2024-6-14                                                              #
#——————————————————————————————————————————————————————————————————————————————————#

import os, pandas as pd, portion as P, operator
from intervaltree import IntervalTree
from functools import reduce

# BASIC FUNCTION #
# USAGE:    load a gff tab file to a pandas dataframe.
# INPUT:    gff_tab_file: the path of the gff tab file.
# OUTPUT:   gff_tab_df: a pandas dataframe created from the gff tab file.
def load_gff_tab_file(gff_tab_file):
    col_name_list = ['region', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    data_type_dict = {
        'region':str, 'source':str, 'feature':str, 
        'start':int, 'end':int, 'score':str, 
        'strand':str, 'frame':str, 'attribute':str
    }
    col_list = [0, 1, 2, 3, 4, 5, 6, 7, 8]
    gff_tab_df = pd.read_csv(gff_tab_file, sep='\t', header=None, comment='#', names=col_name_list, dtype=data_type_dict, usecols=col_list, encoding_errors='replace', encoding='utf-8')
    return gff_tab_df

# BASIC FUNCTION #
# USAGE:    load a bed tab file to a pandas dataframe.
# INPUT:    bed_tab_file: the path of the bed tab file.
# OUTPUT:   bed_tab_df: a pandas dataframe created from the bed tab file.
def load_bed_tab_file(bed_tab_file):
    col_name_list = [
        'chrom', 'chromStart', 'chromEnd', 
        'name', 'score', 'strand', 
        'thickStart', 'thickEnd', 'itemRgb', 
        'blockCount', 'blockSizes', 'blockStarts'
    ]
    data_type_dict = {
        'chrom':str, 'chromStart':int, 'chromEnd':int, 
        'name':str, 'score':int, 'strand':str, 
        'thickStart':int, 'thickEnd':int, 'itemRgb':str, 
        'blockCount':int, 'blockSizes':str, 'blockStarts':str
    }
    col_list = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    bed_tab_df = pd.read_csv(bed_tab_file, sep='\t', header=None, names=col_name_list, dtype=data_type_dict, usecols=col_list)
    return bed_tab_df

# BASIC FUNCTION #
# USAGE:    load a gls align tab file to a pandas dataframe.
# INPUT:    gls_align_tab_file: the path of the gls align tab file.
# OUTPUT:   gls_align_tab_df: a pandas dataframe created from the gls align tab file.
def load_pr_briefing_tsv_file(pr_briefing_tsv_file):
    col_name_list = [
        'id', 'type', 
        'query', 'qstart', 'qend', 'qlen', 
        'region', 'strand', 'start', 'end', 
        'idt', 'cov', 'raws', 'frag', 'cds', 
        'miss', 'ins', 'del', 'preterm', 
        'missinit', 'missterm', 
        'fsh', 'intron', 'lintron', 'ner', 
        'polyA', 'dirrep'
    ]
    data_type_dict = {
        'id':str, 'type':str, 
        'query':str, 'qstart':int, 'qend':int, 'qlen':int, 
        'region':str, 'strand':str, 'start':int, 'end':int, 
        'idt':float, 'cov':float, 'raws':int, 'frag':str, 'cds':str, 
        'miss':int, 'ins':int, 'del':int, 'preterm':int, 
        'missinit':bool, 'missterm':bool, 
        'fsh':int, 'intron':int, 'lintron':int, 'ner':int, 
        'polyA':bool, 'dirrep':bool
    }
    col_list = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26]
    pr_briefing_tsv_df = pd.read_csv(pr_briefing_tsv_file, sep='\t', header=0, names=col_name_list, dtype=data_type_dict, usecols=col_list)
    return pr_briefing_tsv_df

# DATA TYPE DESCRIPTION: <feat_loc_tab_df>
# The <feat_loc_tab_df> expects to be a gff-like pandas dataframe with the following columns:
#   region:   the region id (necessary);
#   strand:   the strand of the feature (necessary);
#   start:    the start location of the feature (necessary; 1-based; like gff format the start value 
#             smaller than the end value);
#   end:      the end location of the feature (necessary; 1-based; like gff format the start value 
#             smaller than the end value);
#   query:    the query id of the feature (optional);
#   score:    the score of the feature, can be any custom value (optional);
#   evalue:   the evalue of the feature (optional).

# BASIC FUNCTION #
# USAGE:    create a IntervalTree object for the feat_loc_tab_df.
# INPUT:    feat_loc_tab_df: a pandas dataframe with specific columns.
# OUTPUT:   feat_loc_itvtree: an IntervalTree object created from the feat_loc_tab_df.
def create_feat_loc_itvtree(feat_loc_tab_df):
    ## initialize an empty IntervalTree object
    feat_loc_itvtree = IntervalTree()
    ## iterate the rows of the feat_loc_tab_df and add current row as an Interval object to the IntervalTree object
    ## along with the row as data.
    for _, row in feat_loc_tab_df.iterrows():
        feat_loc_itvtree[row['start']:row['end']] = row
    return feat_loc_itvtree

# BASIC FUNCTION #
# USAGE:    create a dictionary with key as (region, strand) and value as an IntervalSet object, 
#           witch is the union of all the intervals of the features in the same region and strand.
# INPUT:    feat_loc_tab_df: a pandas dataframe with specific columns.
# OUTPUT:   feat_loc_itvset_dict: a dictionary with key as (region, strand) and value as an IntervalSet object.
def create_feat_loc_itvset_dict(feat_loc_tab_df):
    ## initialize an empty dictionary.
    feat_loc_itvset_dict = {}
    ## group_df the feat_loc_tab_df by region and strand.
    feat_loc_tsv_groupby = feat_loc_tab_df.groupby(['region', 'strand'])
    ## iterate the groupby object and create an IntervalSet object for each group_df.
    for (region, strand), group_df in feat_loc_tsv_groupby:
        ### add a column 'itv' to the group_df, the value of the column is an Interval object.
        group_df = add_itv_obj_2_feat_loc_tsv_df(group_df)
        ### merge the intervals of the group_df to a single IntervalSet object.
        group_feat_loc_itvset = merge_itv_4_feat_loc_tsv_df(group_df)
        ### add the IntervalSet object to the dictionary with key as (region, strand).
        feat_loc_itvset_dict[(region, strand)] = group_feat_loc_itvset
    ## return the dictionary with key as (region, strand) and value as an IntervalSet object.
    return feat_loc_itvset_dict

# BASIC FUNCTION #
# USAGE:    add a column 'itv' to the feat_loc_tab_df, the value of the column is an Interval object 
#           created from the start and end columns.
# INPUT:    feat_loc_tab_df: a pandas dataframe with specific columns.
# OUTPUT:   feat_loc_tab_df: the input feat_loc_tab_df with an additional column 'itv'.
def add_itv_obj_2_feat_loc_tsv_df(feat_loc_tab_df):
    feat_loc_tab_df['itv'] = feat_loc_tab_df.apply(lambda row: P.closed(row['start'], row['end']), axis=1)
    return feat_loc_tab_df

# BASIC FUNCTION #
# USAGE: merge the intervals of the feat_loc_tab_df to a single IntervalSet object.
# INPUT: feat_loc_tab_df: a pandas dataframe with specific columns and an additional column 'itv'.
# OUTPUT: feat_loc_itvset: an IntervalSet object created from the 'itv' column of the feat_loc_tab_df.
def merge_itv_4_feat_loc_tsv_df(feat_loc_tab_df):
    ## create a list of Interval objects from the 'itv' column of the feat_loc_tab_df
    itv_list = list(feat_loc_tab_df['itv'])
    ## use reduce function and | operator to merge the Interval objects in the list to an IntervalSet object
    feat_loc_itvset = reduce(operator.or_, itv_list)
    ## return the IntervalSet object created from the 'itv' column of the feat_loc_tab_df
    return feat_loc_itvset

# BASIC FUNCTION #
# USAGE:    slice the gff tab dataframe by region and output the sliced dataframes to gff files.
#           The function is mainly used for data preprocessing of pseudopipe software.
# INPUT:    gff_tab_df: a pandas dataframe created from a gff tab file.
#           output_dir: the directory to output the sliced gff files.
# OUTPUT:   None, output the sliced gff files to the output_dir.
def slice_output_gff_tab_df(gff_tab_df, output_dir):
    ## group the gff tab dataframe by region.
    gff_tab_groupby = gff_tab_df.groupby('region')
    ## iterate the groupby object and output the group_df to a gff file.
    for region, group_df in gff_tab_groupby:
        group_df.to_csv(os.path.join(output_dir, f'{region}.gff'), sep='\t', index=False)
    return None

# BASIC FUNCTION #
# USAGE:    extract the feature location information from the gff tab dataframe.
#           The function is mainly used for data preprocessing of pseudopipe software.
# INPUT:    gff_tab_df: a pandas dataframe created from a gff tab file.
#           strand: the strand of the feature to extract.
#           feat_kwd_list: a list of keywords to extract the dataframe by the 'feature' column.
#           reg_kwd_list: a list of keywords to extract the dataframe by the 'region' column.
#           id_kwd_list: a list of keywords to extract the dataframe by the 'parsed_id' column.
#           id_attr_kwd_list: a list of keywords to extract the feature id from the attribute column.
#           ign_pfix: a boolean value to indicate whether to ignore the prefix of the feature id.
#           ign_sfix: a boolean value to indicate whether to ignore the suffix of the feature id.
# OUTPUT:   gff_feat_loc_groupby: a groupby object of the feature location information extracted from the gff tab dataframe.
def extr_gff_tab_feat_loc(gff_tab_df, strand=None, feat_kwd_list=None, reg_kwd_list=None, id_kwd_list=None, id_attr_kwd_list=None, ign_pfix=False, ign_sfix=False):
    ## filter the gff tab dataframe by the keyword list in the specific columns.
    gff_tab_df = filt_gff_tab_df_by_kwd(gff_tab_df, strand, feat_kwd_list, reg_kwd_list, id_kwd_list, id_attr_kwd_list, ign_pfix, ign_sfix)
    ## create a feature location tab dataframe from the gff tab dataframe.
    gff_feat_loc_df = gff_tab_df[['region', 'source', 'start', 'end']].copy()
    ## group the feature location tab dataframe by region.
    gff_feat_loc_groupby = gff_feat_loc_df.groupby('region')
    ## return the groupby object of the feature location information extracted from the gff tab dataframe.
    return gff_feat_loc_groupby

# BASIC FUNCTION #
# USAGE:    format the gff tab dataframe to a bed tab dataframe.
# INPUT:    gff_tab_df: a pandas dataframe created from a gff tab file.
#           id_attr_kwd_list: a list of keywords to extract the feature id from the attribute column.
#           ign_pfix: a boolean value to indicate whether to ignore the prefix of the feature id.
#           ign_sfix: a boolean value to indicate whether to ignore the suffix of the feature id.
# OUTPUT:   bed_tab_df: a pandas dataframe formatted from the gff tab dataframe.
def format_gff_2_bed_tab_df(gff_tab_df, id_attr_kwd_list=None, ign_pfix=False, ign_sfix=False):
    ## parse the feature id from the attribute column of the gff tab dataframe.
    gff_tab_df = parse_feat_id_4_gff_tab(gff_tab_df, id_attr_kwd_list, ign_pfix, ign_sfix)
    ## create a bed tab dataframe from the gff tab dataframe.
    bed_tab_df = gff_tab_df[['region', 'start', 'end', 'id', 'score', 'strand']].copy()
    ## convert the start value to 0-based(for gff is 1-based).
    bed_tab_df['start'] = bed_tab_df['start'] - 1
    ## return the bed tab dataframe.
    return bed_tab_df

# BASIC FUNCTION #
# USAGE:    format the gff tab dataframe to a feature location tab dataframe.
# INPUT:    gff_tab_df: a pandas dataframe created from a gff tab file.
# OUTPUT:   feat_loc_tab_df: a pandas dataframe formatted from the gff tab dataframe.
def format_gff_2_feat_loc_tsv_df(gff_tab_df):
    ## create a feature location tab dataframe from the gff tab dataframe.
    feat_loc_tab_df = gff_tab_df[['region', 'feature', 'start', 'end']].copy()
    ## return the feature location tab dataframe.
    return feat_loc_tab_df

# BASIC FUNCTION #
# USAGE:    format the gff tab dataframe to a gls align tab dataframe.
# INPUT:    gff_tab_df: a pandas dataframe created from a gff tab file.
# OUTPUT:   gls_align_tab_df: a pandas dataframe formatted from the gff tab dataframe.
def format_gff_2_gls_align_tab_df(gff_tab_df):
    ## create a gls align tab dataframe from the gff tab dataframe.
    gls_align_tab_df = gff_tab_df[['id', 'region', 'strand', 'start', 'end']].copy()
    ## rename the columns of the gls align tab dataframe.
    gls_align_tab_df = gls_align_tab_df.rename(
        columns={
            'id': 'query', 
            'region': 'subject', 
            'strand': 'sstrand', 
            'start': 'sstart', 
            'end': 'send'}
            )
    ## add additional columns to the gls align tab dataframe.
    gls_align_tab_df = gls_align_tab_df.assign(qstart=None, qend=None)
    ## return the gls align tab dataframe.
    return gls_align_tab_df

# BASIC FUNCTION #
# USAGE:    filter the gff tab dataframe by the keyword list in the specific columns.
# INPUT:    gff_tab_df: a pandas dataframe created from a gff tab file.
#           strand: the strand of the feature to filter.
#           feat_kwd_list: a list of keywords to filter the dataframe by the 'feature' column.
#           reg_kwd_list: a list of keywords to filter the dataframe by the 'region' column.
#           id_kwd_list: a list of keywords to filter the dataframe by the 'parsed_id' column.
#           id_attr_kwd_list: a list of keywords to extract the feature id from the attribute column.
#           ign_pfix: a boolean value to indicate whether to ignore the prefix of the feature id.
#           ign_sfix: a boolean value to indicate whether to ignore the suffix of the feature id.
# OUTPUT:   gff_tab_df: the input gff tab dataframe filtered by the keyword list in the specific columns.
def filt_gff_tab_df_by_kwd(gff_tab_df, strand=None, feat_kwd_list=None, reg_kwd_list=None, id_kwd_list=None, id_attr_kwd_list=None, ign_pfix=False, ign_sfix=False):
    ## parse the feature id from the attribute column of the gff tab dataframe.
    gff_tab_df = parse_feat_id_4_gff_tab(gff_tab_df, id_attr_kwd_list, ign_pfix, ign_sfix)

    ## filter the gff tab dataframe by strand column.
    if strand == '+' or strand == '-':
        gff_tab_df = gff_tab_df[gff_tab_df['strand'] == strand]

    ## filter the gff tab dataframe by the keyword list in the specific columns.
    gff_tab_df = filt_df_by_kwd(gff_tab_df, 'feature', feat_kwd_list)
    gff_tab_df = filt_df_by_kwd(gff_tab_df, 'region', reg_kwd_list)
    gff_tab_df = filt_df_by_kwd(gff_tab_df, 'parsed_id', id_kwd_list)

    ## drop the 'id' column.
    gff_tab_df.drop('id', axis=1)
    ## return the gff tab dataframe filtered by the keyword list in the specific columns.
    return gff_tab_df

# BASIC FUNCTION #
# USAGE:    filter the pr briefing tsv dataframe by the keyword list in the specific columns.
# INPUT:    pr_briefing_tsv_df: a pandas dataframe created from a pr briefing tsv file.
#           strand: the strand of the feature to filter.
#           type_kwd_list: a list of keywords to filter the dataframe by the 'type' column.
#           reg_kwd_list: a list of keywords to filter the dataframe by the 'region' column.
#           id_kwd_list: a list of keywords to filter the dataframe by the 'id' column.
# OUTPUT:   pr_briefing_tsv_df: the input pr briefing tsv dataframe filtered by the keyword list in the specific columns.
def filt_pr_briefing_tsv_df_by_kwd(pr_briefing_tsv_df, strand=None, type_kwd_list=None, reg_kwd_list=None, id_kwd_list=None):
    ## filter the pr briefing tsv dataframe by strand column.
    if strand == '+' or strand == '-':
        pr_briefing_tsv_df = pr_briefing_tsv_df[pr_briefing_tsv_df['strand'] == strand]

    ## filter the pr briefing tsv dataframe by the keyword list in the specific columns.
    pr_briefing_tsv_df = filt_df_by_kwd(pr_briefing_tsv_df, 'type', type_kwd_list)
    pr_briefing_tsv_df = filt_df_by_kwd(pr_briefing_tsv_df, 'region', reg_kwd_list)
    pr_briefing_tsv_df = filt_df_by_kwd(pr_briefing_tsv_df, 'id', id_kwd_list)

    ## return the pr briefing tsv dataframe filtered by the keyword list in the specific columns.
    return pr_briefing_tsv_df

# BASIC FUNCTION #
# USAGE:    filter the dataframe by the keyword list in the specific column.
# INPUT:    df: a pandas dataframe.
#           column: the column name to filter.
#           kwd_list: a list of keywords to filter the dataframe.
# OUTPUT:   df: the input dataframe filtered by the keyword list in the specific column.
def filt_df_by_kwd(df, column, kwd_list):
    ## if kwd_list is provided, filter the dataframe by the keyword list in the specific column.
    if kwd_list:
        ### create a pattern to filter the dataframe by the keyword list in the specific column.
        kwd_ptn_str = '|'.join(kwd_list)
        ### filter the dataframe by the keyword list in the specific column.
        df = df[df[column].str.match(f'{kwd_ptn_str}')].copy()
    return df

# BASIC FUNCTION #
# USAGE:    filter the pr briefing tsv dataframe by the quality of the feature.
# INPUT:    pr_briefing_tsv_df: a pandas dataframe created from a pr briefing tsv file.
#           ident_thold: the threshold of the identity of the feature.
#           cov_thold: the threshold of the coverage of the feature.
#           raws_thold: the threshold of the raw score of the feature.
# OUTPUT:   pr_briefing_tsv_df: the input pr briefing tsv dataframe filtered by the quality of the feature.
def filt_pr_briefing_tsv_df_by_qlt(pr_briefing_tsv_df, ident_thold=None, cov_thold=None, raws_thold=None):
    ## filter the pr briefing tsv dataframe by the quality of the feature.
    if ident_thold:
        pr_briefing_tsv_df = pr_briefing_tsv_df[pr_briefing_tsv_df['idt'] >= ident_thold]
    if cov_thold:
        pr_briefing_tsv_df = pr_briefing_tsv_df[pr_briefing_tsv_df['cov'] >= cov_thold]
    if raws_thold:
        pr_briefing_tsv_df = pr_briefing_tsv_df[pr_briefing_tsv_df['raws'] >= raws_thold]
    ## return the pr briefing tsv dataframe filtered by the quality of the feature.
    return pr_briefing_tsv_df

# BASIC FUNCTION #
# USAGE:    parse the feature id from the attribute column of the gff tab dataframe.
# INPUT:    gff_tab_df: a pandas dataframe created from a gff tab file.
#           id_attr_kwd_list: a list of keywords to extract the feature id from the attribute column.
#           ign_pfix: a boolean value to indicate whether to ignore the prefix of the feature id.
#           ign_sfix: a boolean value to indicate whether to ignore the suffix of the feature id.
# OUTPUT:   gff_tab_df: the input gff tab dataframe with an additional column 'id'.
def parse_feat_id_4_gff_tab(gff_tab_df, id_attr_kwd_list=None, ign_pfix=False, ign_sfix=False):
    ## parse the feature id from the attribute column of the gff tab dataframe
    ### if id_attr_kwd_list is provided, extract the feature id from the attribute column based on the keywords.
    if id_attr_kwd_list:
        #### create a pattern to extract the feature id from the attribute column based on all keywords.
        id_attr_kwd_ptn_str = '|'.join([f'[^;=]*{kwd}[^;=]*=([^;=]+);*' for kwd in id_attr_kwd_list])
        #### extract all possible feature ids with all keywords each row, store in a dataframe.
        extr_id_df = gff_tab_df['attribute'].str.extractall(id_attr_kwd_ptn_str).groupby(level=0).first()
        #### the first column of non-null value in each row is considered as the feature id.
        extr_id_df['id'] = extr_id_df.bfill(axis=1).iloc[:, 0]
        #### join the extracted feature id to the gff tab dataframe.
        gff_tab_df = gff_tab_df.join(extr_id_df[['id']])
    ### if id_attr_kwd_list is not provided, extract the feature id from the attribute column based on the default pattern.
    else:
        #### extract the feature id from the attribute column using split function.
        gff_tab_df['id'] = gff_tab_df['attribute'].str.split(';').str[0].str.split('=').str[1]

    ## ignore the prefix and suffix of the feature id if necessary
    if ign_pfix:
        ### ignore the prefix of the feature id
        gff_tab_df['id'] = gff_tab_df['id'].str.replace('^.*:', '', regex=True)
    if ign_sfix:
        ### ignore the suffix of the feature id
        gff_tab_df['id'] = gff_tab_df['id'].str.replace('\..*$', '', regex=True)
    
    ## return the gff tab dataframe with an additional column 'id'
    return gff_tab_df
