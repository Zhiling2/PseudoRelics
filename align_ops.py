#——————————————————————————————————————————————————————————————————————————————————#
# PROJ:     ToolKits                                                               #
# VER:      1.0                                                                    #
# TYPE:     Module                                                                 #
# DESC:     A module for alignment operations                                      #
# AUTHOR:   Zhilin Guan, Zihe Wang,Shuo Liu                                                           #
# AFFIL:    Wuhan Botanical Garden, Chinese Academy Of Sciences                    #
# E-MAIL:   zhilinggg@163.com                    #
# DATE:     2025-10-18                                                              #
#——————————————————————————————————————————————————————————————————————————————————#

import re, os, pandas as pd, parasail, subprocess, tempfile, shutil, random, numpy as np, portion as P
import sys
sys.path.append(r'/public/home/WT_lius/MySoftware/PseudoRelics/pkgs/rsrc')
import bioinfo_tabs, prog_tabs, dir_tabs
sys.path.append(r'/public/home/WT_lius/MySoftware/PseudoRelics/pkgs/mods')
import fa_ops

# RESOURCE CLASS #
# USAGE:    define the align_tab_types class, which contains the align tab types for different alignment programs.
class align_tab_types:
    m6_like = ['m6', 'cust_hs']
    exone_sugar = ['exone_sugar']

# RESOURCE CLASS #
# USAGE:    define the align_rst_ptns class, which contains the regular expression patterns for parsing alignment 
#           results.
class align_rst_ptns:
    ## exonerate protein2genome alignment custom result patterns.
    exone_cust_align_rst_ptn = re.compile(r'( +\d+ : (?:(?!C4 Alignment:).)+ : +\d+)\n\n(cust_tab(?:\t[^\t\n]+)+)', re.DOTALL)
    exone_cust_tab_ptn = re.compile(r'^(cust_tab(\t[^\t\n]+)+)', re.MULTILINE)
    exone_visu_str_ptn = re.compile(r'( +\d+ : ((?!C4 Alignment:).)+ : +\d+)', re.DOTALL)
    ## exonerate protein2genome alignment visual string prefix and suffix patterns.
    exone_visu_str_pfix_ptn = re.compile(r' +\d+ : ')
    exone_visu_str_sfix_ptn = re.compile(r' : +\d+')
    ## exonerate protein2genome alignment visual string splitcodon patterns.
    exone_visu_str_splitcodon_ptn = re.compile(r'\{|\}')
    ## exonerate protein2genome alignment visual peptide string termcodon patterns.
    exone_visu_pep_str_aa_ptn = re.compile(r'[A-Z][a-z]{2}')
    exone_visu_pep_str_splitaa_ptn = re.compile(r'\{[A-Za-z]+\}')
    exone_visu_pep_str_allaa_ptn = re.compile(r'([A-Z][a-z]{2})|\{([A-Za-z]+)\}')
    exone_visu_pep_str_termcodon_ptn = re.compile(r'\*\*\*')
    ## exonerate protein2genome alignment visual query peptide string intron, target deletion, 
    ## non-equivalenced region, frame shift, insertion and non-amino acid patterns.
    exone_visu_query_pep_str_intron_ptn = re.compile(r' *>>>> Target Intron \d+ >>>> *')
    exone_visu_query_pep_str_tdel_ptn = re.compile(r' *>>>> Target Deletion \d+ >>>> *')
    exone_visu_query_pep_str_ner_ptn = re.compile(r' *>>>> Non-equivalenced Region \d+ >>>> *')
    exone_visu_query_pep_str_frsh_ptn = re.compile(r'(?<!<)-+(?!>)')
    exone_visu_query_pep_str_ins_ptn = re.compile(r'(?:<->)+')
    exone_visu_query_pep_str_nonaa_ptn = re.compile(r' *>>>> Target Intron \d+ >>>> *| *>>>> Target Deletion \d+ >>>> *| *>>>> Non-equivalenced Region \d+ >>>> *|(?<!<)-+(?!>)|\{|\}')
    ## exonerate protein2genome alignment visual alignment link string intron, target deletion, 
    ## non-equivalenced region and frame shift patterns.
    exone_visu_alignlink_str_intron_ptn = re.compile(r'(?<!T:)  *(\d+) bp *')
    exone_visu_alignlink_str_tdel_ptn = re.compile(r'(?<!Q:)  *(\d+) aa *')
    exone_visu_alignlink_str_ner_ptn = re.compile(r' *Q: (\d+) aa / T: (\d+) bp *')
    exone_visu_alignlink_str_frsh_ptn = re.compile(r'#+')
    ## exonerate protein2genome alignment visual target DNA2pep string intron, target deletion, 
    ## non-equivalenced region, frame shift, deletion, non-amino acid(& preterm), complemented peptide, 
    ## truncated peptide, sub peptide, start amino acid, end amino acid and cds intrupted patterns.
    exone_visu_target_dna2pep_str_intron_ptn = re.compile(r'[\+-]* +[\+-]*')
    exone_visu_target_dna2pep_str_tdel_ptn = re.compile(r' *TD *')
    exone_visu_target_dna2pep_str_ner_ptn = re.compile(r' *NER *')
    exone_visu_target_dna2pep_str_frsh_ptn = re.compile(r'#+')
    exone_visu_target_dna2pep_str_del_ptn = re.compile(r'(?:---)+')
    exone_visu_target_dna2pep_str_nonaa_ptn = re.compile(r'[\+-]* +[\+-]*| *TD *| *NER *|#+|\{|\}|(?:---)+')
    exone_visu_target_dna2pep_str_nonaa_preterm_ptn = re.compile(r'[\+-]* +[\+-]*| *TD *| *NER *|#+|\{|\}|(?:---)+|\*\*\*(?!\Z)')
    exone_visu_target_dna2pep_str_complpep_ptn = re.compile(r'^Met(?:[A-Z][a-z]{2}|---|\{[A-Za-z]+\}|[\+-]* +[\+-]*| *TD *)*\*\*\*$')
    exone_visu_target_dna2pep_str_truncpep_ptn = re.compile(r'^Met(?:[A-Z][a-z]{2}|---|\{[A-Za-z]+\}|[\+-]* +[\+-]*| *TD *)*\*\*\*')
    exone_visu_target_dna2pep_str_pep_ptn = re.compile(r'Met(?:[A-Z][a-z]{2}|---|\{[A-Za-z]+\}|[\+-]* +[\+-]*| *TD *)*\*\*\*')
    exone_visu_target_dna2pep_str_start_aa_ptn = re.compile(r'^[A-Z][a-z]{2}|^\*\*\*')
    exone_visu_target_dna2pep_str_end_aa_ptn = re.compile(r'[A-Z][a-z]{2}$|\*\*\*$')
    exone_visu_target_dna2pep_str_cds_itrup_ptn = re.compile(r'\*\*\*|#+| *NER *')
    ## exonerate protein2genome alignment visual target genome string intron, target deletion, 
    ## non-equivalenced region and deletion patterns.
    exone_visu_target_geno_str_intron_tdel_ner_ptn = re.compile(r'[a-z]*\.+[a-z]*')
    exone_visu_target_geno_str_intron_del_ptn = re.compile(r'(?:---)+')

# BASIC FUNCTION #
# USAGE:    load the alignment table file to a pandas dataframe.
# INPUT:    align_tab_file: the path to the alignment table file; 
#           tab_type: the type of the alignment table file, default is 'm6'.
# OUTPUT:   a pandas dataframe of the alignment table.
def load_align_tab_file(align_tab_file, tab_type='m6'):
    ## set the column names, data types and columns to use according to the tab type.
    if tab_type == 'm6':
        col_name_list = ['query', 'subject', 'identity', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
        data_type_dict = {
            'query':str, 'subject':str, 'identity':float, 
            'length':int, 'mismatch':int, 'gapopen':int, 
            'qstart':int, 'qend':int, 'sstart':int, 
            'send':int, 'evalue':float, 'bitscore':float
        }
        col_list = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    elif tab_type == 'cust_hs':
        col_name_list = ['query', 'subject', 'identity', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qcovs']
        data_type_dict = {
            'query':str, 'subject':str, 'identity':float, 
            'length':int, 'mismatch':int, 'gapopen':int, 
            'qstart':int, 'qend':int, 'sstart':int, 
            'send':int, 'evalue':float, 'bitscore':float, 
            'qcovs':float
        }
        col_list = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    elif tab_type == 'exone_sugar':
        col_name_list = ['query', 'qstart', 'qend', 'qstrand', 'target', 'tstart', 'tend', 'tstrand', 'score']
        data_type_dict = {
            'query':str, 'qstart':int, 'qend':int, 
            'qstrand':str, 'target':str, 'tstart':int, 
            'tend':int, 'tstrand':str, 'score':float
        }
        col_list = [1, 2, 3, 4, 5, 6, 7, 8, 9]

    ## read the alignment table file to a pandas dataframe.
    if tab_type in align_tab_types.m6_like:
        align_tab_df = pd.read_csv(align_tab_file, sep='\t', header=None, names=col_name_list, dtype=data_type_dict, usecols=col_list)
    elif tab_type in align_tab_types.exone_sugar:
        align_tab_df = pd.read_csv(align_tab_file, sep=' ', header=None, names=col_name_list, dtype=data_type_dict, usecols=col_list)

    ## return the alignment table dataframe.
    return align_tab_df

# BASIC FUNCTION #
# USAGE:    add strand information for subject or query to the alignment table dataframe.
# INPUT:    align_tab_df: the alignment table dataframe, expect to be m6-like format;
#           seq_mode: the sequence mode to indicate subject or query including: 'subject', 'query', 'both', default is 'subject'.
#               'subject' mode expect to have columns: sstart, send;
#               'query' mode expect to have columns: qstart, qend;
#               'both' mode expect to have columns: sstart, send, qstart, qend.
# OUTPUT:   the alignment table dataframe with strand information added.
def add_strand_4_align_tab_df(align_tab_df, seq_mode='subject'):
    ## seq modes that will be executed are stored in seq mode list, and later iterate to use.
    seq_mode_list = ['subject', 'query'] if seq_mode == 'both' else [seq_mode]

    ## iterate the seq mode list to add strand information for each seq mode.
    for seq_mode in seq_mode_list:
        ### select the strand column name, start column name and end column name according to the seq mode.
        if seq_mode == 'subject':
            strand_col = 'sstrand'
            start_col = 'sstart'
            end_col = 'send'
        elif seq_mode == 'query':
            strand_col = 'qstrand'
            start_col = 'qstart'
            end_col = 'qend'
        else:
            raise Exception('Unsupported seq mode')

        ### add strand information to the alignment table dataframe.
        align_tab_df[strand_col] = np.where(align_tab_df[start_col] <= align_tab_df[end_col], '+', '-')

    ## return the alignment table dataframe with strand information added.
    return align_tab_df

# BASIC FUNCTION #
# USAGE:    add closed interval object from portion module for subject or query to the alignment result tab.
# INPUT:    align_tab_df: the alignment table dataframe, expect to be m6-like format and using 1-based coordinates;
#           seq_mode: the sequence mode to indicate subject or query including: 'subject', 'query', 'both', default is 'subject'.
#               'subject' mode expect to have columns: sstart, send, sstrand;
#               'query' mode expect to have columns: qstart, qend, qstrand;
#               'both' mode expect to have columns: sstart, send, sstrand, qstart, qend, qstrand.
# OUTPUT:   the alignment table dataframe with closed interval object added.
def add_itv_obj_4_align_tab_df(align_tab_df, seq_mode='subject'):
    ## add closed interval object for subject or query according to the seq mode.
    if seq_mode == 'subject':
        ### add closed interval object for subject location.
        align_tab_df['subject_itv'] = align_tab_df.apply(
            lambda row: P.closed(row['sstart'], row['send']) if row['sstrand'] in ['+', '.'] else P.closed(row['send'], row['sstart']), 
            axis=1
        )
    elif seq_mode == 'query':
        ### add closed interval object for query location.
        align_tab_df['query_itv'] = align_tab_df.apply(
            lambda row: P.closed(row['qstart'], row['qend']) if row['qstrand'] in ['+', '.'] else P.closed(row['qend'], row['qstart']), 
            axis=1
        )
    elif seq_mode == 'both':
        ### add closed interval object for both subject and query location.
        align_tab_df[['subject_itv', 'query_itv']] = align_tab_df.apply(
            lambda row: pd.Series({
                'subject_itv': P.closed(row['sstart'], row['send']) if row['sstrand'] in ['+', '.'] else P.closed(row['send'], row['sstart']),
                'query_itv': P.closed(row['qstart'], row['qend']) if row['qstrand'] in ['+', '.'] else P.closed(row['qend'], row['qstart'])
            }), 
            axis=1
        )
    else:
        raise Exception('Unsupported mode')

    ## return the alignment table dataframe with closed interval object added.
    return align_tab_df

# BASIC FUNCTION #
# USAGE:    convert the alignment table dataframe from exonerate to standard m6-like format.
# INPUT:    exone_align_tab_df: the alignment table dataframe from exonerate, expect to be exonerate alignment tab format.
# OUTPUT:   the alignment table dataframe in standard m6-like format.
def format_exone_align_tab_df_2_m6like(exone_align_tab_df):
    ## rename columns for exonerate alignment tab dataframe to standard m6-like format.
    exone_align_tab_df = exone_align_tab_df.rename(columns={'target':'subject', 'tstart':'sstart', 'tend':'send', 'tstrand':'sstrand'})
    ## convert coordinates from 0-based to 1-based.
    exone_align_tab_df.loc[exone_align_tab_df['qstrand'].isin(['+', '.']), 'qstart'] += 1
    exone_align_tab_df.loc[exone_align_tab_df['qstrand'] == '-', 'qend'] += 1
    exone_align_tab_df.loc[exone_align_tab_df['sstrand'].isin(['+', '.']), 'sstart'] += 1
    exone_align_tab_df.loc[exone_align_tab_df['sstrand'] == '-', 'send'] += 1

    ## return the alignment table dataframe in standard m6-like format.
    return exone_align_tab_df

# BASIC FUNCTION #
# USAGE:    filter the alignment table dataframe by quality thresholds.
# INPUT:    align_tab_df: the alignment table dataframe, expect to be m6-like format;
#           id_thold: the identity threshold, default is None;
#           len_thold: the length threshold, default is None;
#           e_thold: the evalue threshold, default is None;
#           bits_thold: the bitscore threshold, default is None;
#           qcovs_thold: the query coverage threshold, default is None.
# OUTPUT:   the alignment table dataframe filtered by quality thresholds.
def filt_align_tab_df_by_qlty(align_tab_df, id_thold=None, len_thold=None, e_thold=None, bits_thold=None, qcovs_thold=None):
    ## filter the alignment table dataframe by quality thresholds.
    if id_thold:
        align_tab_df = align_tab_df[align_tab_df['identity'] >= id_thold].copy()
    if len_thold:
        align_tab_df = align_tab_df[align_tab_df['length'] >= len_thold].copy()
    if e_thold:
        align_tab_df = align_tab_df[align_tab_df['evalue'] <= e_thold].copy()
    if bits_thold:
        align_tab_df = align_tab_df[align_tab_df['bitscore'] >= bits_thold].copy()
    if qcovs_thold and ('qcovs' in align_tab_df.columns):
        align_tab_df = align_tab_df[align_tab_df['qcovs'] >= qcovs_thold].copy()

    ## return the alignment table dataframe filtered by quality thresholds.
    return align_tab_df

# BASIC FUNCTION #
# USAGE:    perform a semi-global seq alignment for two DNA sequences using parasail module.
# INPUT:    query_dna_seq: the query DNA sequence;
#           ref_dna_seq: the reference DNA sequence;
#           gap_open: the gap open penalty, default is 5;
#           gap_extend: the gap extend penalty, default is 1.
# OUTPUT:   the alignment result object.
def sg_align_4_dna(query_dna_seq, ref_dna_seq, gap_open=5, gap_extend=1):
    ## perform a semi-global seq alignment for two DNA sequences using parasail module.
    align_rst = parasail.sg_trace(query_dna_seq, ref_dna_seq, gap_open, gap_extend, parasail.blosum62)

    ## return the alignment result object.
    return align_rst

# BASIC FUNCTION #
# USAGE:    calculate the identity and coverage of the semi-global seq alignment.
# INPUT:    query_tbk: the query traceback of the alignment result;
#           ref_tbk: the reference traceback of the alignment result;
#           ref_len: the length of the reference sequence.
# OUTPUT:   a tuple of identity and coverage.
def calc_sg_align_qlty(query_tbk, ref_tbk, ref_len):
    ## count the matches and mismatches in the alignment.
    mch_mismch_count = sum(c1 == c2 or (c1 != '-' and c2 != '-') for c1, c2 in zip(query_tbk, ref_tbk))
    mch_count = sum(c1 == c2 for c1, c2 in zip(query_tbk, ref_tbk) if c1 != '-' and c2 != '-')

    ## calculate the identity and coverage of the alignment.
    ident = mch_count / mch_mismch_count
    cov = mch_mismch_count / ref_len

    ## return the identity and coverage of the alignment.
    return ident, cov

# BASIC FUNCTION #
# USAGE:    perform exonerate <protein2genome:bestfit> or <protein2genome:local> alignment for query peptide sequence 
#           against target DNA sequence.
# INPUT:    query_pep_seq: the query peptide sequence;
#           target_dna_seq: the target DNA sequence;
#           exone_mode: the exonerate alignment mode, including 'bf' for bestfit and 'lc' for local, default is 'bf';
#           exone_prog: the exonerate program path, default is None, and will use the default path in prog_tabs;
#           wdir: the working directory for exonerate to store temporary files, default is '/dev/shm/exonerate_tmp';
#           del_wdir: whether to delete the working directory after the alignment, default is False.
# OUTPUT:   the exonerate alignment standard output in custom format, stored as a string.
def run_exone_p2g(query_pep_seq, target_dna_seq, exone_mode='bf', exone_prog=None, wdir=None, del_wdir=False):
    ## set the default values for exonerate program and working directory.
    if not exone_prog:
        exone_prog = prog_tabs.align_prog_tabs.loc_align_prog_tab['exonerate']
    if not wdir:
        wdir = dir_tabs.prog_temp_wdir_tab['exonerate']

    ## create the working directory if not exists.
    if not os.path.exists(wdir):
        os.makedirs(wdir)

    ## convert the query peptide sequence and target DNA sequence to fasta file content format.
    query_fa_file_content = fa_ops.format_seq_2_fa_file_cont(query_pep_seq, 'query')
    target_fa_file_content = fa_ops.format_seq_2_fa_file_cont(target_dna_seq, 'target')
    ## delete the input sequences to save memory.
    del query_pep_seq, target_dna_seq

    ## write the query peptide sequence and target DNA sequence to temporary fasta files in working directory, 
    ## and perform exonerate alignment with selected mode on them.
    ### open temporary fasta files for writing.
    with tempfile.NamedTemporaryFile(dir=wdir, mode='w', delete=True) as temp_query_fa_file_obj, tempfile.NamedTemporaryFile(dir=wdir, mode='w', delete=True) as temp_target_fa_file_obj:
        ### write the query peptide sequence and target DNA sequence to the temporary fasta files.
        temp_query_fa_file_obj.write(query_fa_file_content)
        temp_query_fa_file_obj.flush()
        temp_target_fa_file_obj.write(target_fa_file_content)
        temp_target_fa_file_obj.flush()
        ### delete the query_fa_file_content and target_fa_file_content to save memory.
        del query_fa_file_content, target_fa_file_content

        ### acquire the temporary fasta file paths.
        temp_query_fa_file = temp_query_fa_file_obj.name
        temp_target_fa_file = temp_target_fa_file_obj.name
        ### format the command line list for exonerate according to the selected mode.
        if exone_mode == 'bf':
            exone_cmd_list = [exone_prog, 
                              '--exhaustive', 
                              '--model', 'protein2genome:bestfit', 
                              '--query', temp_query_fa_file, 
                              '--target', temp_target_fa_file, 
                              '--querytype', 'protein', 
                              '--targettype', 'dna', 
                              '--showsugar', 'no', 
                              '--showvulgar', 'no', 
                              '--showalignment', 'yes', 
                              '--showquerygff', 'no',
                              '--ryo', r'cust_tab\t%ql\t%qal\t%qab\t%qae\t%tab\t%tae\t%s\t%et\t%ei\t%V\n'
                              ]
        elif exone_mode == 'lc':
            exone_cmd_list = [exone_prog, 
                              '--model', 'protein2genome', 
                              '--query', temp_query_fa_file, 
                              '--target', temp_target_fa_file, 
                              '--querytype', 'protein', 
                              '--targettype', 'dna', 
                              '--showsugar', 'no', 
                              '--showvulgar', 'no', 
                              '--showalignment', 'yes', 
                              '--showquerygff', 'no', 
                              '--ryo', r'cust_tab\t%ql\t%qal\t%qab\t%qae\t%tab\t%tae\t%s\t%et\t%ei\t%V\n'
                              ]
        ### perform exonerate alignment on the temporary fasta files.
        exone_subprocess = subprocess.Popen(exone_cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        ### acquire the exonerate alignment standard output in custom format.
        exone_raw_stdout, error = exone_subprocess.communicate()
        exone_cust_stdout = exone_raw_stdout.decode()
        ### delete the non-decoded output to save memory.
        del exone_raw_stdout
    ## delete the working directory if required.
    if del_wdir:
        shutil.rmtree(wdir)

    ## return the exonerate alignment standard output in custom format.
    return exone_cust_stdout

# BASIC FUNCTION #
# USAGE:    parse the exonerate protein2genome alignment standard output to a tuple containing the formatted 
#           visual string and custom table(which are also stored as tuples) as a parsed result. The parsed 
#           result can be either the best alignment result or the merged alignment result depending on the parse mode.
# INPUT:    exone_cust_stdout: the exonerate protein2genome standard output, stored as a string.
#           parse_mode: the parse mode, including 'best' and 'merge', default is 'merge', see INFO for details.
# OUTPUT:   a tuple of two tuples:
#               the formatted visual string tuple with four strings; 
#               the formatted custom table tuple with parsed alignment information.
# INFO:     I assume that the exonerate run in protein2genome:bestfit or protein2genome:local mode. In some cases of 
#           the local mode, there may be multiple alignment results due to two reasons: (1) large segment deletions that 
#           separate the query and target sequence into multiple non-overlapping and collinear alignments; (2) tandem 
#           repeats or satellites in the target sequence that cause multiple overlapping and non-collinear alignments. 
#           So I provide two modes to parse the exonerate standard output: 'best' and 'merge'. The 'best' mode will find 
#           the best alignment result from the standard output using the raw score as parsed result; The 'merge' mode will
#           merge all eligible alignment results to the parsed result by checking if current neighboring alignments are 
#           non-overlapping and collinear to last one. if so, merge it to the parsed result and consider the gap between 
#           them as a non-equivalenced region; if not, skip the current alignment result and continue.
def parse_exone_p2g_cust_stdout(exone_cust_stdout, parse_mode='merge'):
    ## initialize a tuple to store the parsed result.
    parsed_rst_tup = ()

    ## find all the custom alignment results in the exonerate standard output.
    align_rst_list = re.findall(align_rst_ptns.exone_cust_align_rst_ptn, exone_cust_stdout)
    ## format the visual string and custom table for each custom alignment result, and store them in a new list 
    ## of tuples. each tuple contains two tuples: visual string tuple and custom table tuple(just like parsed result).
    fmtd_align_rst_list = [(format_exone_p2g_visu_str(align_rst[0]), format_exone_p2g_cust_tab(align_rst[1])) for align_rst in align_rst_list]
    ## filter out all alignment results on reverse strand of the target sequence.
    fmtd_align_rst_list = [align_rst for align_rst in fmtd_align_rst_list if align_rst[1][4] <= align_rst[1][5]]

    ## return None if no alignment result is found.
    if len(fmtd_align_rst_list) == 0:
        return None

    ## if use 'best' mode, find the best alignment result using the raw score as parsed result.
    if parse_mode == 'best':
        ### find the best alignment result using the raw score.
        parsed_rst_tup = max(fmtd_align_rst_list, key=lambda x: x[1][6])
        ### return the parsed result tuple.
        return parsed_rst_tup

    ## if use 'merge' mode, merge all eligible alignment results to parsed result.
    if parse_mode == 'merge':
        ### sort the alignment results list by the query alignment start and raw score.
        fmtd_align_rst_list.sort(key=lambda x: (x[1][2], x[1][6]))
        ### iterate the alignment results list to acquire the parsed result.
        for i, align_rst in enumerate(fmtd_align_rst_list):
            #### if the parsed result is empty, append the current alignment result to it, and continue.
            if parsed_rst_tup == ():
                parsed_rst_tup = align_rst
                continue

            #### calculate the step length between the end of the previous parsed result and the start of the 
            #### current alignment result.
            target_step = align_rst[1][4] - parsed_rst_tup[1][5]
            query_step = align_rst[1][2] - parsed_rst_tup[1][3]
            #### if the step length is less than 0, skip the current alignment, for they are overlapping.
            if target_step < 0 or query_step < 0:
                continue
            #### merge the current alignment result to the previous parsed result.
            upd_visu_str_tup = _merge_exone_p2g_visu_str(i, query_step, target_step, parsed_rst_tup[0], align_rst[0])
            upd_cust_tab_tup = _merge_exone_p2g_cust_tab(query_step, target_step, parsed_rst_tup[1], align_rst[1])
            #### update the parsed result tuple.
            parsed_rst_tup = (upd_visu_str_tup, upd_cust_tab_tup)

        ### return the parsed result tuple.
        return parsed_rst_tup

# BASIC FUNCTION #
# USAGE:    format the exonerate protein2genome alignment result's visual string to a tuple containing four strings.
# INPUT:    exone_visu_str: the exonerate visual string, stored as a string.
# OUTPUT:   a tuple of four formatted strings: 
#               query peptide sequence, 
#               alignment link, 
#               target DNA2pep sequence,
#               target genome sequence.
def format_exone_p2g_visu_str(exone_visu_str):
    ## initialize four lists for: 
    ## query peptide sequence, alignment link, target DNA2pep sequence and target genome sequence.
    query_pep_str_list = []
    alignlink_str_list = []
    target_dna2pep_str_list = []
    target_geno_str_list = []

    ## acquire the length of the line prefix and suffix.
    line_pfix_len = len(re.search(align_rst_ptns.exone_visu_str_pfix_ptn, exone_visu_str).group(0))
    line_sfix_len = len(re.search(align_rst_ptns.exone_visu_str_sfix_ptn, exone_visu_str).group(0))

    ## split the exone visual str to blocks, each block contains four lines: 
    ## query peptide, alignment link, target DNA2pep and target genome.
    exone_visu_str_block_list = exone_visu_str.split('\n\n')
    ## iterate the blocks to acquire the formatted exone visual str.
    for block in exone_visu_str_block_list:
        ### split the block to lines.
        block_line_list = block.split('\n')
        ### acquire the exone visual str for each block, and append them to the lists.
        query_pep_str_list.append(block_line_list[0][line_pfix_len:-line_sfix_len])
        alignlink_str_list.append(block_line_list[1][line_pfix_len:])
        target_dna2pep_str_list.append(block_line_list[2][line_pfix_len:])
        target_geno_str_list.append(block_line_list[3][line_pfix_len:-line_sfix_len])

    ## concatenate the all substrings to the formatted strings.
    query_pep_str = ''.join(query_pep_str_list)
    alignlink_str = ''.join(alignlink_str_list)
    target_dna2pep_str = ''.join(target_dna2pep_str_list)
    target_geno_str = ''.join(target_geno_str_list)

    ## pack the formatted exone visual str to a tuple and return.
    exone_visu_str_tup = (query_pep_str, alignlink_str, target_dna2pep_str, target_geno_str)
    return exone_visu_str_tup

# BASIC FUNCTION #
# USAGE:    format the exonerate protein2genome alignment result's custom table to a tuple.
# INPUT:    exone_cust_tab: the exonerate custom table, stored as a string.
# OUTPUT:   a tuple of parsed alignment custom table.
def format_exone_p2g_cust_tab(exone_cust_tab):
    ## split the exone custom table to a list.
    exone_cust_tab_list = exone_cust_tab.split('\t')[1:]
    ## convert all elements to int except the last one
    exone_cust_tab_list[:-1] = map(int, exone_cust_tab_list[:-1])
    ## convert the list to a tuple and return
    return tuple(exone_cust_tab_list)

# SUB FUNCTION #
# USAGE:    merge the previous and current exonerate visual strings using visual link strings representing specific regions.
# INPUT:    iter: the current iteration number;
#           query_step: the step length between two neighboring alignments on the query sequence;
#           target_step: the step length between two neighboring alignments on the target sequence;
#           prev_visu_str_tup: the previous exonerate visual string tuple;
#           curr_visu_str_tup: the current exonerate visual string tuple.
# OUTPUT:   a tuple of merged visual strings.
def _merge_exone_p2g_visu_str(iter, query_step, target_step, prev_visu_str_tup, curr_visu_str_tup):
    ## if query_step or target_step is less than 0, two neighboring alignments are overlapping and will not be merged, 
    ## return the previous visual strings tuple.
    if query_step < 0 or target_step < 0:
        return prev_visu_str_tup

    ## determine contents and fillchars for visual link string according to the query_step and target_step.
    ### if query_step and target_step are both 0, the two neighboring alignments will be merged directly.
    if query_step == 0 and target_step == 0:
        query_pep_link_str = ''
        alignlink_link_str = ''
        target_dna2pep_link_str = ''
        target_geno_link_str = ''
    ### if query_step is 0 and target_step is greater than 0, the two neighboring alignments will be merged with 
    ### target intron link string representing the intron region on the target sequence.
    elif query_step == 0 and target_step > 0:
        query_pep_link_str = f'  >>>> Target Intron {iter} >>>>  '
        alignlink_link_str = f' {target_step} bp '
        target_dna2pep_link_str = ''
        target_geno_link_str = ''
    ### if query_step is greater than 0 and target_step is 0, the two neighboring alignments will be merged with 
    ### target deletion link string representing the deletion region on the target sequence.
    elif query_step > 0 and target_step == 0:
        query_pep_link_str = f'  >>>> Target Deletion {iter} >>>>  '
        alignlink_link_str = f' {query_step} aa '
        target_dna2pep_link_str = ' TD '
        target_geno_link_str = ''
    ### if query_step and target_step are both greater than 0, the two neighboring alignments will be merged with
    ### non-equivalenced region link string representing the non-equivalenced region between the two alignments.
    else:
        query_pep_link_str = f'  >>>> Non-equivalenced Region {iter} >>>>  '
        alignlink_link_str = f' Q: {query_step} aa / T: {target_step} bp '
        target_dna2pep_link_str = ' NER '
        target_geno_link_str = ''

    ## format the visual link strings.
    ### acquire the maximum length of the visual link strings.
    link_str_max_len = max(len(query_pep_link_str), len(alignlink_link_str), len(target_dna2pep_link_str), len(target_geno_link_str))
    ### center the visual link strings to the maximum length.
    query_pep_link_str = query_pep_link_str.center(link_str_max_len)
    alignlink_link_str = alignlink_link_str.center(link_str_max_len)
    target_dna2pep_link_str = target_dna2pep_link_str.center(link_str_max_len)
    target_geno_link_str = target_geno_link_str.center(link_str_max_len, '.')

    ## merge the current visual strings to the previous visual strings.
    merged_query_pep_str = prev_visu_str_tup[0] + query_pep_link_str + curr_visu_str_tup[0]
    merged_alignlink_str = prev_visu_str_tup[1] + alignlink_link_str + curr_visu_str_tup[1]
    merged_target_dna2pep_str = prev_visu_str_tup[2] + target_dna2pep_link_str + curr_visu_str_tup[2]
    merged_target_geno_str = prev_visu_str_tup[3] + target_geno_link_str + curr_visu_str_tup[3]
    merged_visu_str_tup = (merged_query_pep_str, merged_alignlink_str, merged_target_dna2pep_str, merged_target_geno_str)
    ## return the merged visual strings tuple.
    return merged_visu_str_tup

# SUB FUNCTION #
# USAGE:    merge the previous and current exonerate custom tables.
# INPUT:    query_step: the step length between two neighboring alignments on the query sequence;
#           target_step: the step length between two neighboring alignments on the target sequence;
#           prev_cust_tab_tup: the previous exonerate custom table tuple;
#           curr_cust_tab_tup: the current exonerate custom table tuple.
# OUTPUT:   a tuple of merged custom tables.
def _merge_exone_p2g_cust_tab(query_step, target_step, prev_cust_tab_tup, curr_cust_tab_tup):
    ## if query_step or target_step is less than 0, two neighboring alignments are overlapping and will not be merged, 
    ## return the previous custom table tuple.
    if query_step < 0 or target_step < 0:
        return prev_cust_tab_tup

    ## determine the vulgar link string according to the query_step and target_step.
    if query_step == 0 and target_step == 0:
        vulgar_link_str = ' '
    elif query_step == 0 and target_step > 0:
        vulgar_link_str = f' I {query_step} {target_step} '
    elif query_step > 0 and target_step == 0:
        vulgar_link_str = f' G {query_step} {target_step} '
    else:
        vulgar_link_str = f' N {query_step} {target_step} '

    ## merge the current custom table tuple to the previous custom table tuple.
    ### acquire updated data for merged custom table tuple.
    query_len = curr_cust_tab_tup[0]
    query_align_len = prev_cust_tab_tup[1] + curr_cust_tab_tup[1]
    query_align_begin = prev_cust_tab_tup[2]
    query_align_end = curr_cust_tab_tup[3]
    target_align_begin = prev_cust_tab_tup[4]
    target_align_end = curr_cust_tab_tup[5]
    raw_score = prev_cust_tab_tup[6] + curr_cust_tab_tup[6]
    total_count = prev_cust_tab_tup[7] + curr_cust_tab_tup[7]
    ident_count = prev_cust_tab_tup[8] + curr_cust_tab_tup[8]
    vulgar_str = prev_cust_tab_tup[9] + vulgar_link_str + curr_cust_tab_tup[9]
    ### create the merged custom table tuple.
    merged_cust_tab_tup = (
        query_len, query_align_len, query_align_begin, query_align_end, target_align_begin, target_align_end, 
        raw_score, total_count, ident_count, vulgar_str
        )
    ## return the merged custom table tuple.
    return merged_cust_tab_tup

# BASIC FUNCTION #
# USAGE:    analyze the exonerate protein2genome alignment result for detailed information.
# INPUT:    exone_parsed_rst_tup: the exonerate custom result tuple, including formatted visual string tuple and custom table tuple;
# OUTPUT:   a dictionary of alignment detailed information.
def anal_exone_p2g_align_rst(exone_parsed_rst_tup):
    ## initialize all detailed information variables for analysis.
    ### identity, coverage, fragment, CDS
    ident, cov, frag, cds = 0, 0, None, None
    ### match, missense, insertion, deletion, premature termination, missing initiation, missing termination
    mch_count, missen_count, ins_count, del_count, prem_term_count = 0, 0, 0, 0, 0
    miss_init, miss_term = False, False
    ### frameshift, intron, long intron, non-equivalenced region
    frsh_count, intron_count, lintron_count, ner_count = 0, 0, 0, 0
    ### intermediate variables: status of 5' and 3' alignments or CDS prime
    align_5_prime, align_3_prime, cds_5_prime, cds_3_prime = False, False, False, None

    ## parse the input exonerate parsed result tuple, and acquire all necessary information to analyze the 
    ## alignment result.
    ### parse the exonerate visual string tuple and custom table tuple from exonerate parsed result tuple.
    (exone_visu_str_tup, exone_cust_tab_tup) = exone_parsed_rst_tup
    ### acquire alignment information from the exonerate custom table tuple.
    (
        query_len, _, 
        query_align_begin, query_align_end, _, _, 
        _, align_total_count, align_ident_count, _
    ) = exone_cust_tab_tup

    ## perform the analysis.
    ### analysis step 1: analyze exonerate visual string tuple, and parse analysis results.
    visu_str_detail_tup = anal_exone_p2g_visu_str(exone_visu_str_tup)
    (
        mch_count, missen_count, ins_count, del_count, prem_term_count, 
        frsh_count, intron_count, lintron_count, ner_count, 
        tstart_aa, tend_aa, first_cds_itrup_struct, first_cds_itrup_struct_at_end
    ) = visu_str_detail_tup
    ### analysis step 2: calculate the identity and coverage.
    ident = round((align_ident_count / align_total_count) * 100, 1)
    cov = round((mch_count / query_len) * 100, 1)
    ### analysis step 3: analyze the fragment, CDS, missing initiation and termination status.
    #### analyze the 5' and 3' alignments or CDS prime status.
    align_5_prime = bool(query_align_begin == 0 and tstart_aa)
    align_3_prime = bool(query_align_end == query_len and tend_aa)
    cds_5_prime = align_5_prime and tstart_aa == 'Met'
    if first_cds_itrup_struct is None:
        cds_3_prime = 'O'
    elif '*' in first_cds_itrup_struct:
        cds_3_prime = 'C' if first_cds_itrup_struct_at_end and query_align_end == query_len else 'T'
    else:
        cds_3_prime = 'O'
    #### analyze the fragment status.
    frag_status_dict = {(True, True): 'C', (True, False): '3M', 
                        (False, True): '5M', (False, False): 'BM'}
    frag = frag_status_dict[(align_5_prime, align_3_prime)]
    #### analyze the CDS status.
    cds_status_dict = {(True, 'C'): 'I', (True, 'T'): '3T', (True, 'O'): '3O', 
                       (False, 'C'): '5M', (False, 'T'): '5M3T', (False, 'O'): '5M3O'}
    cds = cds_status_dict[(cds_5_prime, cds_3_prime)]
    #### analyze the missing initiation and termination status.
    miss_init = not cds_5_prime
    miss_term = not cds_3_prime == 'C'

    ## pack the alignment detailed information to a tuple.
    align_detail_tup = (
        ident, cov, frag, cds, 
        mch_count, missen_count, ins_count, del_count, prem_term_count, miss_init, miss_term, 
        frsh_count, intron_count, lintron_count, ner_count
    )

    ## return the alignment detailed information tuple.
    return align_detail_tup

# BASIC FUNCTION #
# USAGE:    analyze the exonerate protein2genome alignment result's visual string tuple for detailed information.
# INPUT:    exone_visu_str_tup: the exonerate visual string tuple, including query peptide, alignment link,
#                                 target DNA2pep and target genome strings.
# OUTPUT:   a tuple of alignment visual string detailed information.
def anal_exone_p2g_visu_str(exone_visu_str_tup):
    ## initialize all detailed information variables for analysis.
    ### match, missense, insertion, deletion, premature termination
    mch_count, missen_count, ins_count, del_count, prem_term_count = 0, 0, 0, 0, 0
    ### frameshift, intron, long intron, non-equivalenced region
    frsh_count, intron_count, lintron_count, ner_count = 0, 0, 0, 0
    ### amino acid of target dna2pep string start and end, first CDS interruption structure, 
    ### and whether the first CDS interruption structure is at the end of the target dna2pep string.
    tstart_aa, tend_aa, first_cds_itrup_struct, first_cds_itrup_struct_at_end = None, None, None, False
    ### two intermediate variables: single deletion, target deletion
    sgl_del_count, tdel_count = 0, 0

    ## parse and preprocess the input exonerate visual string tuple to acquire first 3 visual strings, 
    ## the query and target amino acid lists.
    ### acquire the query peptide, alignment link, target DNA2pep strings from the visual string tuple.
    query_pep_str = exone_visu_str_tup[0]
    alignlink_str = exone_visu_str_tup[1]
    target_dna2pep_str = exone_visu_str_tup[2]
    ### acquire the amino acid strings and lists for both query pep and target DNA2pep strings.
    query_aa_str = extr_aa_4_exone_visu_pep_str(query_pep_str, pep_type='query')
    target_aa_str = extr_aa_4_exone_visu_pep_str(target_dna2pep_str, pep_type='target')
    query_aa_list = [query_aa_str[i:i+3] for i in range(0, len(query_aa_str), 3)]
    target_aa_list = [target_aa_str[i:i+3] for i in range(0, len(target_aa_str), 3)]

    ## analysis step 1: iterate the query and target amino acid lists to count the match, missense, insertion, 
    ## single deletion and premature termination.
    for q_aa, t_aa in zip(query_aa_list, target_aa_list):
        if q_aa == t_aa:
            mch_count += 1
        elif t_aa == '---':
            sgl_del_count += 1
        elif q_aa == '<->':
            ins_count += 1
        elif t_aa == '***':
            mch_count += 1
            prem_term_count += 1
        else:
            mch_count += 1
            missen_count += 1

    ## analysis step 2: count all frameshift, intron, long intron, target deletion and non-equivalenced 
    ## region in the visual strings.
    ### count the frameshift from the target DNA2pep string.
    frsh_count = len(re.findall(align_rst_ptns.exone_visu_target_dna2pep_str_frsh_ptn, target_dna2pep_str))
    ### count the intron and long intron from the alignment link string.
    intron_mch_list = re.findall(align_rst_ptns.exone_visu_alignlink_str_intron_ptn, alignlink_str)
    intron_len_list = [int(mch) for mch in intron_mch_list]
    intron_count = len(intron_len_list)
    lintron_count = len([len for len in intron_len_list if len >= 2000])
    ### count the target deletion from the alignment link string.
    target_del_mch_list = re.findall(align_rst_ptns.exone_visu_alignlink_str_tdel_ptn, alignlink_str)
    target_del_len_list = [int(mch) for mch in target_del_mch_list]
    tdel_count = sum(target_del_len_list)
    del_count = tdel_count + sgl_del_count
    ### count the non-equivalenced region from the query peptide string.
    ner_count = len(re.findall(align_rst_ptns.exone_visu_query_pep_str_ner_ptn, query_pep_str))

    ## analysis step 3: try to find the start and end amino acid, the first CDS interruption structure in the 
    ## target DNA2pep string.
    tstart_aa_mch = re.search(align_rst_ptns.exone_visu_target_dna2pep_str_start_aa_ptn, target_dna2pep_str)
    tend_aa_mch = re.search(align_rst_ptns.exone_visu_target_dna2pep_str_end_aa_ptn, target_dna2pep_str)
    tstart_aa = tstart_aa_mch.group(0) if tstart_aa_mch else None
    tend_aa = tend_aa_mch.group(0) if tend_aa_mch else None
    first_cds_itrup_struct_mch = re.search(align_rst_ptns.exone_visu_target_dna2pep_str_cds_itrup_ptn, target_dna2pep_str)
    first_cds_itrup_struct = first_cds_itrup_struct_mch.group(0) if first_cds_itrup_struct_mch else None
    first_cds_itrup_struct_at_end = first_cds_itrup_struct_mch.end() == len(target_dna2pep_str) if first_cds_itrup_struct_mch else None

    ## pack the visual string detailed information to a tuple.
    visu_str_detail_tup = (
        mch_count, missen_count, ins_count, del_count, prem_term_count, 
        frsh_count, intron_count, lintron_count, ner_count, 
        tstart_aa, tend_aa, first_cds_itrup_struct, first_cds_itrup_struct_at_end
    )

    ## return the visual string detailed information tuple.
    return visu_str_detail_tup

# BASIC FUNCTION #
# USAGE:    extract the amino acid sequence from the exonerate protein2genome alignment result's visual peptide string.
# INPUT:    exone_visu_pep_str: the exonerate visual peptide string, stored as a string;
#           pep_type: the type of the peptide, 'query' or 'target', default is 'query';
#           drop_term: whether to drop the terminal codon, default is False;
#           short_aa: whether to convert the amino acid to short form, default is False.
# OUTPUT:   the amino acid sequence extracted from the visual peptide string.
def extr_aa_4_exone_visu_pep_str(exone_visu_pep_str, pep_type='query', drop_term=False, short_aa=False):
    ## determine the regular expression patterns for removing non-amino acid substrings(intron, target deletion, 
    ## non-equivalenced region, split codon and frameshift)
    if pep_type == 'query':
        nonaa_ptn = align_rst_ptns.exone_visu_query_pep_str_nonaa_ptn
    elif pep_type == 'target':
        nonaa_ptn = align_rst_ptns.exone_visu_target_dna2pep_str_nonaa_ptn

    ## remove all non-amino acid substrings from the visual peptide string.
    aa_str = re.sub(nonaa_ptn, '', exone_visu_pep_str)

    ## drop the terminal codon if necessary.
    if drop_term:
        termcodon_ptn = align_rst_ptns.exone_visu_pep_str_termcodon_ptn
        aa_str = re.sub(termcodon_ptn, '', aa_str)

    ## convert the amino acid to short form if necessary.
    if short_aa:
        aa_str = conv_aa_seq_2_short(aa_str)

    ## return the extracted amino acid sequence.
    return aa_str

# BASIC FUNCTION #
# USAGE:    convert the amino acid sequence to short form.
# INPUT:    aa_seq: the amino acid sequence in full form, stored as a string.
# OUTPUT:   the amino acid sequence in short form.
def conv_aa_seq_2_short(aa_seq):
    ## check if the input amino acid sequence is empty, if so, return an empty string.
    if aa_seq == '':
        return ''

    ## convert the amino acid sequence to short form.
    assert len(aa_seq) % 3 == 0
    f2s_aa_conv_tab = bioinfo_tabs.genetic_code_tabs.f2s_aa_conv_tab
    aa_list = [aa_seq[i:i+3] for i in range(0, len(aa_seq), 3)]
    short_aa_seq = ''.join([f2s_aa_conv_tab.get(aa, 'X') for aa in aa_list])
    ## return the amino acid sequence in short form.
    return short_aa_seq

# BASIC FUNCTION #
# USAGE:    extract the sub visual strings tuple from the exonerate protein2genome alignment result's visual strings tuple.
# INPUT:    exone_visu_str_tup: the exonerate visual strings tuple, including query peptide, alignment link, 
#                                 target DNA2pep and target genome strings;
#           ref_str:    the reference string for extracting the sub visual strings, default is an empty string;
#           mch_ptn:    the regular expression pattern for the matched region, default is '.*';
#           mch_op:     the match operation, 'ps' for preserving the matched region, 'rm' for removing the matched region, 
#                       default is 'ps'.
# OUTPUT:   a tuple of extracted sub visual strings.
def extr_sub_visu_str_tup(visu_str_tup, ref_str='', mch_ptn=r'.*', mch_op='ps'):
    ## acquire all matched regions in the reference string.
    mch_loc_list = [(mch.start(), mch.end()) for mch in re.finditer(mch_ptn, ref_str)]

    ## acquire the extracted location list according to the match operation.
    ### initialize the extracted location list.
    extr_loc_list = []
    ### determine the extracted location list according to the match operation.
    #### if the match operation is 'ps', the extracted location list is the matched location list.
    if mch_op == 'ps':
        extr_loc_list = mch_loc_list
    #### if the match operation is 'rm', the extracted location list is the non-matched location list.
    elif mch_op == 'rm':
        ##### initialize the last end position variable.
        last_end = 0
        ##### iterate the matched location list to acquire the non-matched location list.
        for start, end in mch_loc_list:
            extr_loc_list.append((last_end, start))
            last_end = end
        extr_loc_list.append((last_end, len(ref_str)))
    #### if the match operation is invalid, raise a ValueError.
    else:
        raise ValueError('Invalid match operation!')

    ## extract the sub visual strings tuple according to the extracted location list.
    ### initialize the sub visual strings list.
    sub_visu_str_list = []
    ### iterate the visual strings tuple to extract the sub visual strings.
    for visu_str in visu_str_tup:
        sub_visu_str = ''.join([visu_str[start:end] for start, end in extr_loc_list])
        sub_visu_str_list.append(sub_visu_str)
    ### pack the sub visual strings list to a tuple.
    sub_visu_str_tup = tuple(sub_visu_str_list)

    ## return the extracted sub visual strings tuple.
    return sub_visu_str_tup

# BASIC FUNCTION #
# USAGE:    find the longest common subsequence between two sequences, which is usually direct repeat of genomic sequences.
# INPUT:    seq_1: the first sequence;
#           seq_2: the second sequence.
# OUTPUT:   the longest common subsequence between the two sequences.
def find_lgst_comm_subseq(seq_1, seq_2):
    ## check if the two input sequences are empty, if so, return an empty string.
    if seq_1 == '' or seq_2 == '':
        return ''

    ## perform a local alignment between the two sequences using parasail module.
    align_result = parasail.sw_trace(seq_1, seq_2, 5, 2, parasail.blosum62)
    ## acquire the query and reference traceback strings from the alignment result.
    query_tbk = align_result.traceback.query
    ref_tbk = align_result.traceback.ref

    ## iterate the query and reference traceback strings to turn 
    ## the mismatched nucleotides to '-' in the query traceback strings.
    ## the transformed query traceback strings will be used to find the longest common subsequence.
    trans_query_traceback = ''
    for query_nt, ref_nt in zip(query_tbk, ref_tbk):
        if query_nt != ref_nt:
            query_nt = '-'
        trans_query_traceback += query_nt

    ## split the transformed query traceback strings to substrings by '-',
    ## and find the longest substring as the longest common subsequence.
    comm_subseq_list = trans_query_traceback.split('-')
    lgst_comm_subseq = max(comm_subseq_list, key=len)

    ## return the longest common subsequence between the two sequences.
    return lgst_comm_subseq

# BASIC FUNCTION #
# USAGE:    calculate the p-value of the longest common subsequence length between two sequences by Monte Carlo simulation.
# INPUT:    seq1: the first sequence;
#           seq2: the second sequence;
#           iters: the number of Monte Carlo iters, default is 1000.
# OUTPUT:   the p-value of the longest common subsequence length between the two sequences.
def monte_carlo_test_4_lcs(seq1, seq2, iters=1000):
    ## acquire the length of the original longest common subsequence between the two sequences.
    org_lcs_len = len(find_lgst_comm_subseq(seq1, seq2))
    ## perform Monte Carlo simulation to acquire the p-value of the longest common subsequence length.
    random_lcs_len_list = []
    for _ in range(iters):
        random_seq1 = list(seq1)
        random.shuffle(random_seq1)
        random_lcs_len_list.append(len(find_lgst_comm_subseq(''.join(random_seq1), seq2)))
    ## calculate the p-value of the longest common subsequence length.
    p_value = np.mean([l > org_lcs_len for l in random_lcs_len_list])
    ## return the p-value of the longest common subsequence length between the two sequences.
    return p_value

# BASIC FUNCTION #
# USAGE:    merge HSPs of same query that share collinearity on same region and strand to a gene-like structure.
# INPUT:    align_tab_df: a dataframe expect to have columns: query, subject, sstrand, qstart, qend, sstart, send;
#           sbj_gap_thold: the threshold of the gap between two HSPs on subject, default is 5000;
#           mp_pool: the multiprocessing pool object, default is None.
# OUTPUT:   a dataframe of merged gene-like structure intervals with following columns: query, subject, sstrand, 
#           qstart, qend, sstart, send.
def merge_align_hsps_2_gls(align_tab_df, sbj_gap_thold=5000, mp_pool=None):
    ## For two HSPs, it is necessary to have the relative positional information on both the query and subject 
    ## (so called collinearity) in order to determine if they should be merged. For that, we add some auxiliary 
    ## columns to the dataframe.

    ## int tag values start from 0 to infinity, rows with same tag value will be merged.
    align_tab_df['merge_tag'] = None
    ## send of the previous neighboring HSP for calculating subject_step.
    align_tab_df['prev_send'] = None
    ## interval median of current and previous HSPs on query for collinearity check between these two subject 
    ## neighboring HSPs.
    align_tab_df['query_loc_median'] = (align_tab_df['qstart'] + align_tab_df['qend']) / 2
    align_tab_df['prev_query_loc_median'] = None
    ## the distance between previous HSP's end and current HSP's start on subject.
    align_tab_df['subject_step'] = None
    ## the distance between previous and current HSP's interval median on query.
    align_tab_df['query_loc_median_step'] = None

    ## group the dataframe by query, subject and sstrand.
    align_tab_groupby = align_tab_df.groupby(['query', 'subject', 'sstrand'])
    del align_tab_df

    ## convert the groupby object to a iterator of group tuples, each tuple contains the group's key, dataframe 
    ## and the subject gap threshold.
    align_tab_group_iter = ((group_key, group_df, sbj_gap_thold) for group_key, group_df in align_tab_groupby)

    ## main functional procedure, merge the HSPs for each group.
    ### if mp_pool is None, we will use single process to merge the HSPs as default setting.
    if mp_pool == None:
        #### initialize a list to store the merged group dataframes.
        group_gls_align_tab_df_list = []
        #### iterate the list of group tuples to run the _merge_hsps_4_group function for each group.
        for group_align_tab_tup in align_tab_group_iter:
            group_gls_align_tab_df_list.append(_merge_hsps_4_group(group_align_tab_tup))
        #### concatenate the merged group dataframes to a output dataframe.
        gls_align_tab_df = pd.concat(group_gls_align_tab_df_list)
    ### if mp_pool is not None, and it is a Pool object, we will use multiple processes to merge the HSPs.
    elif mp_pool != None:
        #### Use the pool's imap function to apply _merge_hsps_4_group to each group, and get the merged group 
        #### dataframes to a list.
        group_gls_align_tab_df_list = mp_pool.imap(_merge_hsps_4_group, align_tab_group_iter)
        #### Concatenate the merged group dataframes to a output dataframe.
        gls_align_tab_df = pd.concat(group_gls_align_tab_df_list)
    ### if mp_pool is not None, but it is not a Pool object, we will raise an exception and exit the program.
    else:
        raise Exception('Unsupported multiprocessing pool object')

    ## drop the merge_tag column and reset the index.
    gls_align_tab_df.drop(columns=['merge_tag'], inplace=True)
    gls_align_tab_df.reset_index(drop=True, inplace=True)

    ## return the dataframe of merged gene-like structure intervals.
    return gls_align_tab_df

# SUB FUNCTION #
# USAGE:    merge the HSPs of the same query that share collinearity on the same region and strand to a gene-like structure.
# INPUT:    group_align_tab_tup: a tuple of group key, dataframe and subject gap threshold;
# OUTPUT:   the dataframe of merged gene-like structure intervals.
def _merge_hsps_4_group(group_align_tab_tup):
    ## unpack the group's key and dataframe.
    group_key, group_align_tab_df, sbj_gap_thold = group_align_tab_tup
    ## acquire group's query, subject, sstrand information.
    query = group_key[0]
    subject = group_key[1]
    sstrand = group_key[2]

    ## perform three tasks according to the sstrand value:
    ## 1.sort group_align_tab_df by subject position, HSP closer to the 5' end of the subject strand will be placed first;
    ## 2.determine the subject_step_ccl_func which will be used to calculate the subject_step;
    ## 3.determine the agg function's parameters to acquire the sstart and send position of the merged 
    ## GLSs(gene-like structure), which will be used in the groupby aggregation later.
    if sstrand == '+':
        group_align_tab_df = group_align_tab_df.sort_values(by=['sstart', 'send'])  #task no.1
        subject_step_ccl_func = lambda sstart_sers, prev_send_sers: sstart_sers - prev_send_sers    #task no.2
        sstart_acq_key = 'min'  #task no.3
        send_acq_key = 'max'    #task no.3
    elif sstrand == '-':
        group_align_tab_df = group_align_tab_df.sort_values(by=['sstart', 'send'], ascending=False) #task no.1
        subject_step_ccl_func = lambda sstart_sers, prev_send_sers: prev_send_sers - sstart_sers    #task no.2
        sstart_acq_key = 'max'  #task no.3
        send_acq_key = 'min'    #task no.3
    ## reset the index of group_align_tab_df.
    group_align_tab_df.reset_index(drop=True, inplace=True)

    ## calculate the distance between two neighboring HSPs on subject and query.
    ### assign values to auxiliary columns.
    group_align_tab_df['prev_send'] = group_align_tab_df['send'].shift(1)
    group_align_tab_df['prev_query_loc_median'] = group_align_tab_df['query_loc_median'].shift(1)
    ### calculate the subject_step.
    group_align_tab_df['subject_step'] = subject_step_ccl_func(group_align_tab_df['sstart'], group_align_tab_df['prev_send']).fillna(0).astype(int)
    ### calculate the query_loc_median_step.
    group_align_tab_df['query_loc_median_step'] = (group_align_tab_df['query_loc_median'] - group_align_tab_df['prev_query_loc_median']).fillna(0).astype(int)

    ## calculate the merge_tag value.
    ### initialize the merge_tag value for the first row.
    group_align_tab_df.loc[group_align_tab_df.index[0], 'merge_tag'] = 0
    ### calculate the merge_tag value from the second to last row, here we use vector operation to avoid the loop.
    #### determine whether the current HSP should be merged with the previous HSP, boolean value will be assigned first, 
    #### True if no merging(or splitting) between the current and previous HSP;
    #### False if merging(or no splitting) between the current and previous HSP.
    group_align_tab_df.loc[1:, 'merge_tag'] = ((group_align_tab_df['query_loc_median_step'] < 0) | (group_align_tab_df['subject_step'] > sbj_gap_thold)).iloc[1:]
    #### convert the boolean value to integer merge_tag value using cumsum, cumsum will accumulate: 
    #### the True value to 1, which means updating the merge_tag value due to no merging(or splitting);
    #### the False value to 0, which means continuing the previous merge_tag value due to merging(or no splitting).
    group_align_tab_df['merge_tag'] = group_align_tab_df['merge_tag'].astype(int).cumsum()

    ## merge the HSPs with the same merge_tag value to GLS(gene-like structure), and format the GLSs to the output dataframe.
    ### group the dataframe by merge_tag and aggregate information to acquire the merged GLSs each group.
    group_gls_align_tab_df = group_align_tab_df.groupby('merge_tag').agg({
        'qstart': 'min',
        'qend': 'max',
        'sstart': sstart_acq_key,
        'send': send_acq_key,
    }).reset_index()
    ### add the common query, subject, sstrand information to the output dataframe.
    group_gls_align_tab_df['query'] = query
    group_gls_align_tab_df['subject'] = subject
    group_gls_align_tab_df['sstrand'] = sstrand

    ## return the dataframe of merged gene-like structure intervals.
    return group_gls_align_tab_df
