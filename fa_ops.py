#——————————————————————————————————————————————————————————————————————————————————#
# PROJ:     ToolKits                                                               #
# VER:      1.0                                                                    #
# TYPE:     Module                                                                 #
# DESC:     A module for fasta file operations                                     #
# AUTHOR:                                                                  #
# AFFIL:    Wuhan Botanical Garden, Chinese Academy Of Sciences                    #
# E-MAIL:                       #
# DATE:     2024-6-14                                                              #
#——————————————————————————————————————————————————————————————————————————————————#

import re, pandas as pd
import sys
sys.path.append(r'/public/home/WT_lius/MySoftware/PseudoRelics/pkgs/rsrc')
import bioinfo_tabs

# RESOURCE CLASS #
# USAGE: define the fa_ptns class, which contains all the regular expression patterns for fasta file processing.
class fa_ptns:
    mult_line_break_ptn = re.compile('\n{2,}')
    line_break_ptn = re.compile(r'\n')
    fa_id_ptn = re.compile(r'>(\S+)')
    fa_header_ptn = re.compile(r'>(\S+)(.*)')
    fa_item_ptn = re.compile(r'>(\S+)(.*)\n([^>]+)')
    fa_seq_end_ptn = re.compile(r'\*|\.')
    fa_seq_slice_ptn = re.compile(r'(.{60})')
    fa_id_ptn = re.compile(r'([^\s.]+)(?:\.\d+)?(?:[^\s\.]*)')
    #fa_id_ptn = re.compile(r'([^\s.)]+\.\d+)(?:.*)')

# BASIC FUNCTION #
# USAGE:    load fasta file content from a fasta file.
# INPUT:    fa_file: the path of the fasta file.
# OUTPUT:   fa_file_cont: the content of the fasta file, stored in a string.
def load_fa_file(fa_file):
    with open(fa_file) as file:
        fa_file_cont = file.read()
    return fa_file_cont

# BASIC FUNCTION #
# USAGE:    merge the content of multiple fasta files.
# INPUT:    fa_file_cont_list: a list of fasta file content.
# OUTPUT:   merged_fa_file_cont: the merged content of multiple fasta files, stored in a string.
def merge_fa_file_cont(fa_file_cont_list):
    merged_fa_file_cont = '\n'.join(fa_file_cont_list) + '\n'
    merged_fa_file_cont = re.sub(fa_ptns.mult_line_break_ptn, '\n', merged_fa_file_cont)
    return  merged_fa_file_cont

# BASIC FUNCTION #
# USAGE:    load the header of the fasta file content.
# INPUT:    fa_file_cont: the content of the fasta file, stored in a string.
#           ret_type: the type of the return value, 'line' or 'id'.
# OUTPUT:   fa_file_header_df: the header of the fasta file content, stored in a pandas dataframe.
def load_fa_file_header(fa_file_cont, ret_type='line'):
    if ret_type == 'line':
        fa_file_header_list = ['>' + match.group(1) + match.group(2) for match in re.finditer(fa_ptns.fa_item_ptn, fa_file_cont)]
    if ret_type == 'id':
        fa_file_header_list = [match.group(1) for match in re.finditer(fa_ptns.fa_item_ptn, fa_file_cont)]
    fa_file_header_df = pd.DataFrame(fa_file_header_list, columns=['header'])
    print('fasta item count: ' + str(fa_file_header_df.shape[0]))
    return fa_file_header_df

# BASIC FUNCTION #
# USAGE:    load the id map file.
# INPUT:    id_map_file: the path of the id map file.
# OUTPUT:   id_map_df: the id map file content, stored in a pandas dataframe.
def load_fa_id_map_file(id_map_file):
    id_map_df = pd.read_csv(id_map_file, sep='\t', header=None, names=['org_id', 'targ_id'])
    return id_map_df

# BASIC FUNCTION #
# USAGE:    format the content of the fasta file into a dictionary.
# INPUT:    fa_file_cont: the content of the fasta file, stored in a string;
#           rmv_line_break: whether to remove the line break in the sequence;
#           rmv_seq_end: whether to remove the end symbol of the sequence.
# OUTPUT: fa_item_dict: the content of the fasta file, stored in a dictionary.
def format_fa_file_cont_2_dict(fa_file_cont, rmv_line_break=True, rmv_seq_end=True, case='nochg'):
    fa_item_dict = {}
    rmv_line_break_func = (lambda x: re.sub(fa_ptns.line_break_ptn, '', x)) if rmv_line_break else (lambda x: x)
    rmv_seq_end_func = (lambda x: re.sub(fa_ptns.fa_seq_end_ptn, '', x)) if rmv_seq_end else (lambda x: x)
    if case == 'upper':
        conv_seq_func = (lambda x: x.upper())
    elif case == 'lower':
        conv_seq_func = (lambda x: x.lower())
    elif case == 'nochg':
        conv_seq_func = (lambda x: x)
    else:
        raise ValueError('sequence case not supported.')
    for match in fa_ptns.fa_item_ptn.finditer(fa_file_cont):
        fa_item_id = match.group(1)
        header = '>' + match.group(1) + match.group(2)
        seq = match.group(3)
        seq = rmv_line_break_func(seq)
        seq = rmv_seq_end_func(seq)
        seq = conv_seq_func(seq)
        fa_item_dict[fa_item_id] = [header, seq]
    return fa_item_dict

# BASIC FUNCTION #
# USAGE:    add an end symbol to the sequence of each fasta item sequence.
# INPUT:    fa_item_dict: the content of the fasta file, stored in a dictionary.
# OUTPUT:   fa_item_dict: the content of the fasta file, stored in a dictionary.
def add_end_4_fa_seq(fa_item_dict):
    fa_item_dict = {fa_item_id: [fa_item[0], fa_item[1] + '*'] for fa_item_id, fa_item in fa_item_dict.items()}
    return fa_item_dict

# BASIC FUNCTION #
# USAGE:    convert the id of the fasta item dictionary.
# INPUT:    fa_item_dict: the content of the fasta file, stored in a dictionary;
#           id_map_df: the id map file content, stored in a pandas dataframe.
# OUTPUT:   fa_item_dict: the content of the fasta file, stored in a dictionary.
def conv_fa_id(fa_item_dict, id_map_df):
    fa_item_df = format_fa_item_dict_2_df(fa_item_dict)
    fa_item_df = fa_item_df.merge(id_map_df, left_on='id', right_on='org_id', how='left')
    fa_item_df = fa_item_df.drop(columns=['id', 'org_id'])

    fa_item_df['header'] = fa_item_df.apply(lambda row: row['header'].replace(row['header'].split(' ')[0], '>' + row['targ_id']), axis=1)
    fa_item_df = fa_item_df.rename(columns={'targ_id': 'id'})
    fa_item_dict = fa_item_df.set_index('id')[['header', 'seq']].T.to_dict('list')

    return fa_item_dict

# BASIC FUNCTION #
# USAGE:    format the content of the fasta file stored in a dictionary into a pandas dataframe.
# INPUT:    fa_item_dict: the content of the fasta file, stored in a dictionary.
# OUTPUT:   fa_item_df: the content of the fasta file, stored in a pandas dataframe.
def format_fa_item_dict_2_df(fa_item_dict):
    fa_item_df = pd.DataFrame.from_dict(fa_item_dict, orient='index', columns=['header', 'seq'])
    fa_item_df = fa_item_df.reset_index().rename(columns={'index': 'id'})
    return fa_item_df

# BASIC FUNCTION #
# USAGE:    extract sub sequence from the fasta item dictionary by the location.
# INPUT:    fa_item_dict: the content of the fasta file, stored in a dictionary;
#           fa_item_id: the id of the fasta item;
#           loc_start: the start coordinate of the sub sequence;
#           loc_end: the end coordinate of the sub sequence;
#           strand: the strand of the sequence, '+' or '-';
#           seq_type: the type of the sequence, 'DNA' or 'RNA'.
# OUTPUT:   fa_item_seq: the sub sequence.
def extr_fa_seq_by_loc(fa_item_dict, fa_item_id, loc_start, loc_end, strand='+', seq_type='DNA'):
    fa_item = fa_item_dict[fa_item_id]
    fa_item_seq = fa_item[1]
    fa_item_seq = extr_sub_seq_by_loc(fa_item_seq, loc_start, loc_end, strand, seq_type)
    return fa_item_seq

# BASIC FUNCTION #
# USAGE:    extract the sub sequence by the location.
# INPUT:    seq: the sequence;
#           loc_start: the start coordinate of the sub sequence;
#           loc_end: the end coordinate of the sub sequence;
#           strand: the strand of the sequence, '+' or '-';
#           seq_type: the type of the sequence, 'DNA' or 'RNA'.
# OUTPUT:   seq: the sub sequence.
def extr_sub_seq_by_loc(seq, loc_start, loc_end, strand='+', seq_type='DNA'):
    if strand == '+':
        seq = extr_fwd_sub_seq_by_coord(seq, loc_start, loc_end)
    if strand == '-':
        seq = extr_rev_compl_sub_seq_by_coord(seq, loc_start, loc_end, seq_type)
    return seq

# BASIC FUNCTION #
# USAGE:    extract reverse complementary sub sequence by pure coordinate number, note that the loc_start 
#           and loc_end are expected to be 1-based, and will be converted to 0-based in the following code.
# INPUT:    seq: the sequence;
#           loc_start: the start coordinate of the sub sequence;
#           loc_end: the end coordinate of the sub sequence;
#           seq_type: the type of the sequence, 'DNA' or 'RNA'.
# OUTPUT:   seq: the reverse complementary sub sequence.
def extr_rev_compl_sub_seq_by_coord(seq, loc_start, loc_end, seq_type='DNA'):
    ## extract the forward sub sequence by the coordinate number.
    seq = seq[loc_start-1:loc_end]
    ## convert the sequence to reverse complementary sequence.
    seq = conv_nt_seq_2_rev_compl(seq, seq_type)
    return seq

# BASIC FUNCTION #
# USAGE:    extract forward sub sequence by pure coordinate number, note that the loc_start and loc_end 
#           are expected to be 1-based, and will be converted to 0-based in the following code.
# INPUT:    seq: the sequence;
#           loc_start: the start coordinate of the sub sequence;
#           loc_end: the end coordinate of the sub sequence.
# OUTPUT:   seq: the forward sub sequence.
def extr_fwd_sub_seq_by_coord(seq, loc_start, loc_end):
    ## extract the forward sub sequence by the coordinate number.
    seq = seq[loc_start-1:loc_end]
    return seq

# BASIC FUNCTION #
# USAGE:    convert nt sequence to reverse complementary sequence.
# INPUT:    nt_seq: the nucleotide sequence;
#           seq_type: the type of the sequence, 'DNA' or 'RNA'.
# OUTPUT:   nt_seq: the reverse complementary sequence of the input sequence
#                   (the sequence is reversed and then converted to complementary sequence)
def conv_nt_seq_2_rev_compl(nt_seq, seq_type='DNA'):
    ## reverse the sequence.
    nt_seq = nt_seq[::-1]
    ## convert the nt sequence to forward complementary sequence.
    nt_seq = conv_nt_seq_2_fwd_compl(nt_seq, seq_type)
    return nt_seq

# BASIC FUNCTION #
# USAGE:    convert nt sequence to forward complementary sequence, but not reverse.
# INPUT:    nt_seq: the nucleotide sequence;
#           seq_type: the type of the sequence, 'DNA' or 'RNA'.
# OUTPUT:   nt_seq: the forward complementary sequence of the input sequence.
def conv_nt_seq_2_fwd_compl(nt_seq, seq_type='DNA'):
    ## initialize the translation table
    if seq_type == 'DNA':
        nt_trans_tab = str.maketrans(bioinfo_tabs.genetic_code_tabs.dna_bp_tab)
    elif seq_type == 'RNA':
        nt_trans_tab = str.maketrans(bioinfo_tabs.genetic_code_tabs.rna_bp_tab)
    else:
        raise ValueError('sequence type not supported.')
    ## convert the nt sequence to complementary sequence, 
    ## only base converted, sequence not reversed.
    nt_seq = nt_seq.translate(nt_trans_tab)
    return nt_seq

# BASIC FUNCTION #
# USAGE:    extract the sequence from the fasta item dictionary by the id.
# INPUT:    fa_item_dict: the content of the fasta file, stored in a dictionary;
#           fa_item_id: the id of the fasta item.
# OUTPUT:   fa_item_seq: the sequence of the fasta item.
def extr_fa_seq_by_id(fa_item_dict, fa_item_id):
    fa_item = fa_item_dict[fa_item_id]
    fa_item_seq = fa_item[1]
    return fa_item_seq

# BASIC FUNCTION #
# USAGE:    filter the fasta item dictionary by the keyword.
# INPUT:    fa_item_dict: the content of the fasta file, stored in a dictionary; 
#           kwd_list: the list of the keyword; 
#           mode: the filter mode, 'kp' for keep, 'rm' for remove.
# OUTPUT:   fa_item_dict: the content of the fasta file, stored in a dictionary.
def filt_fa_item_by_kwd(fa_item_dict, kwd_list=[], mode='kp'):
    if len(kwd_list) == 0:
        return fa_item_dict

    keyword_ptn = '|'.join(kwd_list)
    if mode == 'kp':
        fa_item_dict = {fa_item_id:fa_item for fa_item_id, fa_item in fa_item_dict.items() if re.search(keyword_ptn, fa_item_id, re.IGNORECASE)}
    elif mode == 'rm':
        fa_item_dict = {fa_item_id:fa_item for fa_item_id, fa_item in fa_item_dict.items() if not re.search(keyword_ptn, fa_item_id, re.IGNORECASE)}
    else:
        raise ValueError('filter mode not supported.')

    return fa_item_dict

# BASIC FUNCTION #
# USAGE:    extract the fasta item by ids to a new dictionary.
# INPUT:    fa_item_dict: the content of the fasta file, stored in a dictionary; 
#           id_list: the list of the id;
#           mode: the extract mode, 'kp' for keep, 'rm' for remove.
# OUTPUT:   extr_fa_item_dict: the extracted fasta item dictionary.
def extr_fa_item_by_id(fa_item_dict, id_list=[], mode='kp'):
    if len(id_list) == 0:
        extr_fa_item_dict = fa_item_dict

    if mode == 'kp':
        extr_fa_item_dict = {id: fa_item_dict[id] for id in id_list if id in fa_item_dict}
    elif mode == 'rm':
        for id in id_list:
            fa_item_dict.pop(id, None)
        extr_fa_item_dict = fa_item_dict
    else:
        raise ValueError('extract mode not supported.')

    return extr_fa_item_dict

# BASIC FUNCTION #
# USAGE:    eliminate the redundant fasta item by the sequence.
# INPUT:    fa_item_dict: the content of the fasta file, stored in a dictionary;
#           adj_id_header: whether to adjust the id and header of the fasta item.
# OUTPUT:   fa_item_dict: the content of the fasta file, stored in a dictionary.
def elim_redun_fa_item_by_seq(fa_item_dict, adj_id_header=False):
    ## format the fasta item dictionary to a pandas dataframe.
    fa_item_df = format_fa_item_dict_2_df(fa_item_dict)
    ## delete the original fasta item dictionary to save memory.
    del fa_item_dict

    ## determine the groupby and get sequence functions.
    group_df_func = (lambda x: x.groupby('seq'))
    get_seq_func = 'first'
    ## determine the get id and get header functions.
    if adj_id_header:
        ### get id function: the id of the fasta item will be the merged ids.
        get_id_func = (lambda x: '|'.join(fa_item_df.loc[x.index, 'id'].astype(str)))
        ### get header function: The header of the fasta item will be composed of the merged ids as the sequence id 
        ### and no sequence info.
        get_header_func = (lambda x: '>' + '|'.join(fa_item_df.loc[x.index, 'id'].astype(str)))
    else:
        ### get id function: the id of the fasta item will be the first sequence's id.
        get_id_func = 'first'
        ### get header function: The header of the fasta item will be composed of the first sequence's id as the 
        ### sequence id and the merged ids as the sequence info.
        get_header_func = (lambda x: '>' + fa_item_df.loc[x.index[0], 'id'] + ' ' + '|'.join(fa_item_df.loc[x.index, 'id'].astype(str)))

    ## eliminate the redundant fasta item using the determined groupby and aggregation functions.
    fa_item_df = group_df_func(fa_item_df).agg({
        'id': get_id_func, 
        'header': get_header_func, 
        'seq': get_seq_func, 
    }).reset_index(drop=True)

    ## convert the pandas dataframe back to the fasta item dictionary.
    fa_item_dict = fa_item_df.set_index('id')[['header', 'seq']].T.to_dict('list')
    ## return the fasta item dictionary.
    return fa_item_dict

# BASIC FUNCTION #
# USAGE:    eliminate the redundant fasta item by the parent id.
# INPUT:    fa_item_dict: the content of the fasta file, stored in a dictionary;
#           adj_id_header: whether to adjust the id and header of the fasta item.
# OUTPUT:   fa_item_dict: the content of the fasta file, stored in a dictionary.
def elim_redun_fa_item_by_id(fa_item_dict, adj_id_header=False):
    ## format the fasta item dictionary to a pandas dataframe.
    fa_item_df = format_fa_item_dict_2_df(fa_item_dict)
    ## delete the original fasta item dictionary to save memory.
    del fa_item_dict
    ## extract the parent id and sequence length of the fasta item.
    fa_item_df['parent_id'] = fa_item_df['id'].str.extract(fa_ptns.fa_id_ptn, expand=False)
    fa_item_df['seq_len'] = fa_item_df['seq'].str.len()
    ## group the dataframe by the parent id.
    fa_item_groupby = fa_item_df.groupby('parent_id')

    ## determine the get id and header functions.
    if adj_id_header:
        get_id_header_func = get_adj_id_header
    else:
        get_id_header_func = get_raw_id_header

    ## initialize a dataframe to store the result.
    proc_fa_item_df = pd.DataFrame(columns=['id', 'header', 'seq'])

    ## iterate the groupby object.
    for parent_id, group_fa_item_df in fa_item_groupby:
        sotd_group_fa_item_df = group_fa_item_df.reset_index(drop=True).sort_values(by='seq_len', ascending=False)
        lgst_seq = sotd_group_fa_item_df['seq'].iloc[0]
        lgst_seq_id = sotd_group_fa_item_df['id'].iloc[0]
        merged_id = '|'.join(group_fa_item_df['id'].astype(str))
        id, header = get_id_header_func(parent_id, lgst_seq_id, merged_id)
        ## append the result to the dataframe.
        new_row = pd.DataFrame([{'id': id, 'header': header, 'seq': lgst_seq}])
        proc_fa_item_df = pd.concat([proc_fa_item_df, new_row], ignore_index=True)
            
    ## convert the pandas dataframe back to the fasta item dictionary.
    proc_fa_item_dict = proc_fa_item_df.set_index('id')[['header', 'seq']].T.to_dict('list')

    ## return the fasta item dictionary.
    return proc_fa_item_dict

# SUB FUNCTION #
# USAGE:    get the adjusted id and header of the fasta item.
# INPUT:    parent_id: the parent id of the fasta item; 
#           lgst_seq_id: the id of the largest sequence in the group; 
#           merged_ids: the merged ids of the fasta item.
# OUTPUT:   A tuple of the new id and header of the fasta item.
def get_adj_id_header(parent_id, lgst_seq_id, merged_id):
    return parent_id, f'>{parent_id} {merged_id}'

# SUB FUNCTION #
# USAGE:    get the raw id and header of the fasta item.
# INPUT:    parent_id: the parent id of the fasta item; 
#           lgst_seq_id: the id of the largest sequence in the group; 
#           merged_ids: the merged ids of the fasta item.
# OUTPUT:   A tuple of the new id and header of the fasta item.
def get_raw_id_header(parent_id, lgst_seq_id, merged_id):
    return lgst_seq_id, f'>{lgst_seq_id} {merged_id}'

# BASIC FUNCTION #
# USAGE:    calculate the length of the fasta item sequence and store the length in a dictionary.
# INPUT:    fa_item_dict: the content of the fasta file, stored in a dictionary.
# OUTPUT:   fa_item_len_dict: the length of the fasta item sequence, stored in a dictionary.
def calc_fa_item_seq_len(fa_item_dict):
    fa_item_len_dict = {}
    for fa_item_id, fa_item in fa_item_dict.items():
        fa_item_seq = fa_item[1]
        fa_item_seq_len = len(fa_item_seq)
        fa_item_len_dict[fa_item_id] = fa_item_seq_len
    return fa_item_len_dict

# SUB FUNCTION #
# USAGE:    format the fasta file content for output.
# INPUT:    fa_item_match: the match object of the fasta item.
# OUTPUT:   fmtd_fa_item_match: the formatted fasta item match.
def _format_fa_file_cont_4_output(fa_item_match):
    match_header = '>' + fa_item_match.group(1) + fa_item_match.group(2)
    match_seq = fa_item_match.group(3)
    match_header = re.sub('\t', ' ', match_header)
    match_seq = re.sub(fa_ptns.fa_seq_slice_ptn, '\\1\n', match_seq)
    match_seq = '\n' + match_seq + '\n'
    fmtd_fa_item_match = match_header + match_seq
    return fmtd_fa_item_match

# BASIC FUNCTION #
# USAGE:    output the fasta item dictionary to a fasta file.
# INPUT:    fa_item_dict: the content of the fasta file, stored in a dictionary;
#           output_file: the path of the output file;
#           write_mode: the mode of the write operation.
# OUTPUT:   None, the fasta item dictionary will be written to the output file.
def output_fa_item_dict_2_file(fa_item_dict, output_file, write_mode='w'):
    fa_file_cont = '\n'.join([fa_item[0] + '\n' + fa_item[1] for fa_item in fa_item_dict.values()])
    fa_file_cont = re.sub(fa_ptns.fa_item_ptn, _format_fa_file_cont_4_output, fa_file_cont)
    fa_file_cont = re.sub(fa_ptns.mult_line_break_ptn, '\n', fa_file_cont)
    with open(output_file, write_mode) as file:
        file.write(fa_file_cont)
    return None

# BASIC FUNCTION #
# USAGE:    format the sequence to the fasta file content.
# INPUT:    seq: the sequence;
#           item_id: the id of the sequence.
# OUTPUT:   seq_fa_file_cont: the content of the sequence, stored in a string.
def format_seq_2_fa_file_cont(seq, item_id='seq'):
    seq = re.sub(fa_ptns.fa_seq_slice_ptn, '\\1\n', seq)
    seq_fa_file_cont = f'>{item_id}\n{seq}\n'
    return seq_fa_file_cont

# BASIC CLASS #
# USAGE:     a class for creating fasta database and access fasta sequence by location(fai index) 
#            or access fasta location by sub sequence(fm index).
# INPUT:     see the __init__ method.
# OUTPUT:    see the extraction methods and the location methods depending on the index type.
class fa_database:
    # INIT METHOD #
    # USAGE:     store the input parameters, initialize the public intermediate data set. 
    # INPUT:     see the input parameters.
    # OUTPUT:    none.
    # UPDATE:    all the input parameters; 
    #            all the public intermediate data set.
    # CLEAR:     none.
    def __init__(self, fa_file_i, idx_type_i, wdir_i=None):
        self.fa_file = fa_file_i
        self.idx_type = idx_type_i
        self.wdir = wdir_i

        self.fai_idx_df = None

    # BASIC METHOD #
    # USAGE:    build index for the fasta file in fai format.
    # INPUT:    none.
    # OUTPUT:   none.
    # UPDATE:   the fai index dataframe.
    # CLEAR:    none.
    def build_fai_idx(self):
        ## initialize the fasta index dataframe.
        ## including 5 standard columns: NAME, LENGTH, OFFSET, LINEBASES, LINEWIDTH, 
        ## and 2 additional columns: LINEBASESCHG, LINEWIDTHCHG for checking the 
        ## consistency of the sequence line length.
        self.fai_idx_df = pd.DataFrame(
            columns=[
                'NAME', 'LENGTH', 'OFFSET', 'LINEBASES', 'LINEWIDTH', 
                'LINEBASESCHG', 'LINEWIDTHCHG'
                ]
            )

        ## open the fasta file in read mode.
        with open(self.fa_file, 'r') as file:
            ## iterate each line of the fasta file.
            for line in file:
                ### if the line is a header line.
                if line.startswith('>'):
                    #### accquire the fasta item id.
                    fa_item_id = line.strip().split(' ')[0][1:]
                    #### accquire the offset of the fasta item 
                    #### which is the first character of the fasta item sequence(not the header).
                    offset = file.tell()
                    #### initialize a index item in the fai index dataframe for the fasta item.
                    self.fa_idx_df = pd.concat(
                        [
                            self.fa_idx_df, 
                            pd.DataFrame([{
                                'NAME': fa_item_id, 
                                'LENGTH': 0, 
                                'OFFSET': offset, 
                                'LINEBASES': 0, 
                                'LINEWIDTH': 0, 
                                'LINEBASESCHG': 0, 
                                'LINEWIDTHCHG': 0
                            }])
                        ], 
                        ignore_index=True
                        )
                ### if the line is a sequence line.
                else:
                    #### acquire the sub sequence length and full line length of the sequence line.
                    currt_linebases = len(line.strip())
                    currt_linewidth = len(line)

                    #### check the consistency of the current sequence line length.
                    if self.fa_idx_df.iloc[-1]['LENEBASES'] != currt_linebases:
                        self.fa_idx_df.iloc[-1]['LINEBASESCHG'] += 1
                    if self.fa_idx_df.iloc[-1]['LINEWIDTH'] != currt_linewidth:
                        self.fa_idx_df.iloc[-1]['LINEWIDTHCHG'] += 1

                    #### update the length of the fasta item sequence.
                    self.fa_idx_df.loc[self.fa_idx_df.shape[0]-1, 'LENGTH'] += currt_linebases
                    self.fa_idx_df.loc[self.fa_idx_df.shape[0]-1, 'LINEBASES'] = currt_linebases
                    self.fa_idx_df.loc[self.fa_idx_df.shape[0]-1, 'LINEWIDTH'] = currt_linewidth

        ## check the consistency of the sequence line length.
        if (self.fa_idx_df['LINEBASESCHG'] > 2).any():
            raise ValueError('inconsistent sequence line sub sequence length, index building failed.')
        if (self.fa_idx_df['LINEWIDTHCHG'] > 2).any():
            raise ValueError('inconsistent sequence line full length, index building failed.')

        ## drop the additional columns.
        self.fa_idx_df = self.fa_idx_df.drop(columns=['LINEBASESCHG', 'LINEWIDTHCHG'])
        ## set the index of the fai index dataframe to the NAME column.
        self.fa_idx_df = self.fa_idx_df.set_index('NAME')

        ## save the fai index dataframe to the index file.
        if self.wdir is not None:
            self.fa_idx_df.to_csv(self.wdir + 'loc_index.fai', sep='\t', header=False, index=False)

        return None

    # BASIC METHOD #
    # USAGE:    extract fasta sequence by item id & location using the fai index.
    # INPUT:    fa_item_id: the id of the fasta item; 
    #           loc_start: the start coordinate of the sub sequence, 1-based;
    #           loc_end: the end coordinate of the sub sequence, 1-based;
    #           strand: the strand of the sequence, '+' or '-'; 
    #           seq_type: the type of the sequence, 'DNA', 'RNA', or 'PEP'.
    # OUTPUT:   seq: the extracted sequence.
    # UPDATE:   the fai index dataframe.
    # CLEAR:    none.
    def extr_fa_seq_by_fai(self, fa_item_id, loc_start=1, loc_end=-1, strand='+', seq_type='DNA'):
        ## check the fai index dataframe.
        if self.fai_idx_df is None:
            raise ValueError('fai index not exist, please build the index first.')

        ## acquire the index values of the fasta item.
        item_length, item_offset, item_linebases, item_linewidth = self.fai_idx_df.loc[fa_item_id, ['LENGTH', 'OFFSET', 'LINEBASES', 'LINEWIDTH']]

        ## check and correct the input location values.
        if loc_start < 1 or loc_start > item_length:
            loc_start = 1
        if loc_end < 1 or loc_end > item_length:
            loc_end = item_length
        if loc_start > loc_end:
            loc_start, loc_end = loc_end, loc_start

        ## calculate the start and end offset of the sub sequence.
        start_offset = item_offset + (loc_start-1)//item_linebases*item_linewidth + (loc_start-1)%item_linebases
        end_offset = item_offset + (loc_end-1)//item_linebases*item_linewidth + (loc_end-1)%item_linebases
        read_offset = end_offset - start_offset + 1

        ## open the fasta file in read mode.
        with open(self.fa_file, 'r') as file:
            ## move the file pointer to the start offset of the sub sequence.
            file.seek(start_offset)
            ## read the sub sequence.
            seq = file.read(read_offset)

        seq = seq.strip()
        ## convert the sequence to reverse complementary sequence if the strand is '-'
        if strand == '-' and seq_type != 'PEP':
            seq = conv_nt_seq_2_rev_compl(seq, seq_type)

        return seq
