#——————————————————————————————————————————————————————————————————————————————————#
# PROJ:     PseudoRelics                                                           #
# VER:      1.0                                                                    #
# TYPE:     Work Flow                                                              #
# DESC:     The main workflow for PseudoRelics                                     #
# AUTHOR:                                                                  #
# AFFIL:    Wuhan Botanical Garden, Chinese Academy Of Sciences                    #
# E-MAIL:                       #
# DATE:     2024-6-18                                                              #
#——————————————————————————————————————————————————————————————————————————————————#

import os, sys, pandas as pd, multiprocessing as mp
from intervaltree import Interval, IntervalTree
from datetime import datetime
import sys
sys.path.append(r'/public/home/WT_lius/MySoftware/PseudoRelics/pkgs/rsrc')
import prog_tabs, dir_tabs
sys.path.append(r'/public/home/WT_lius/MySoftware/PseudoRelics/pkgs/mods')
import align_ops, anno_ops, fa_ops, motif_ops
import pickle

# WORKFLOW CLASS #
# USAGE:     manage the entire pseudogene annotation pipeline.
# INPUT:     see the __init__ method.
# OUTPUT:    see the format_output_rst method.
class workflow_pseudo_relics:

    # CLASS VARIABLES #
    ## initialize the exonerate program path and working directory.
    exone_prog = None
    exone_wdir = None

    # INIT METHOD #
    # USAGE:     store the input parameters from command line, initialize the public intermediate data set 
    #            and multiprocessing pool.
    # INPUT:     see the input parameters.
    # OUTPUT:    none.
    # UPDATE:    all the input parameters; 
    #            all the public intermediate data set; 
    #            the multiprocessing pool.
    # CLEAR:     none.
    def __init__(self, 
                 proc_num_i=1, 
                 ppl_mode_i='IGR', 
                 align_tab_file_i=None, gff_tab_file_i=None, geno_fa_file_i=None, pep_fa_file_i=None, 
                 output_dir_i=None, 
                 output_gff_i=False, output_bed_i=False, 
                 output_raw_fa_i=False, output_cds_fa_i=False, output_pep_fa_i=False, 
                 exone_prog_i=None, exone_wdir_i=None, 
                 align_tab_type_i='m6', sbj_gap_thold_i=5000, 
                 gff_avd_feat_kwd_list_i=[], gff_id_kwd_list_i=None, gff_id_ign_sfix_i=False, gff_id_ign_pfix_i=False, 
                 gls_ups_len_i=1000, gls_downs_len_i=1000, gls_ups_ext_len_i=100, gls_downs_ext_len_i=100, 
                 idt_thold_i=60, cov_thold_i=40, qlt_mode_i='idt_cov', 
                 cds_pep_preterm_op_i='dp'
                 ):

        ## default gff feature keyword list for different pipeline mode.
        mode_2_feat_kwd_list_dict = {'IGR': ['gene'], 'GR': ['gene']}

        ## input parameters from command line.
        ### process number
        self.proc_num = proc_num_i
        ### pipeline mode
        self.ppl_mode = ppl_mode_i
        ### input files
        self.align_tab_file = align_tab_file_i
        self.gff_tab_file = gff_tab_file_i
        self.geno_fa_file = geno_fa_file_i
        self.pep_fa_file = pep_fa_file_i
        ### output directory and output files
        self.output_dir = output_dir_i
        self.output_gff = output_gff_i
        self.output_bed = output_bed_i
        self.output_raw_fa = output_raw_fa_i
        self.output_cds_fa = output_cds_fa_i
        self.output_pep_fa = output_pep_fa_i
        ### exonerate program path and working directory
        workflow_pseudo_relics.exone_prog = exone_prog_i if exone_prog_i else prog_tabs.align_prog_tabs.loc_align_prog_tab['exonerate']
        workflow_pseudo_relics.exone_wdir = exone_wdir_i if exone_wdir_i else dir_tabs.prog_temp_wdir_tab['exonerate']
        ### alignment tab file type and subject gap threshold
        self.align_tab_type = align_tab_type_i
        self.sbj_gap_thold = sbj_gap_thold_i
        ### gff feature keyword list and id keyword list and ignore flags
        self.gff_avd_feat_kwd_list = mode_2_feat_kwd_list_dict.get(self.ppl_mode) if not gff_avd_feat_kwd_list_i else gff_avd_feat_kwd_list_i
        self.gff_id_kwd_list = gff_id_kwd_list_i
        self.gff_id_ign_sfix = gff_id_ign_sfix_i
        self.gff_id_ign_pfix = gff_id_ign_pfix_i
        ### GLS upstream, downstream and extension length
        self.gls_ups_len = gls_ups_len_i
        self.gls_downs_len = gls_downs_len_i
        self.gls_ups_ext_len = gls_ups_ext_len_i
        self.gls_downs_ext_len = gls_downs_ext_len_i
        ### identity and coverage threshold, quality mode
        self.idt_thold = idt_thold_i
        self.cov_thold = cov_thold_i
        self.qlt_mode = qlt_mode_i
        ### cds and pep premature termination operation flag
        self.cds_pep_preterm_op = cds_pep_preterm_op_i

        ## variables of intermediate data set.
        self.gls_align_tab_df = None
        self.gls_align_tab_groupby = None
        self.gene_loc_itvset_dict = None
        self.pep_fa_item_dict = None
        self.geno_fa_item_dict = None
        self.geno_fa_item_len_dict = None
        self.gls_info_df = pd.DataFrame(
            columns=[
                'query', 'id', 'type', 
                'region', 'strand', 'start', 'end', 
                'qstart', 'qend', 'qlen', 
                'idt', 'cov', 'raws', 'frag', 'cds', 
                'miss', 'ins', 'del', 'preterm', 
                'missinit', 'missterm', 
                'frsh', 'intron', 'lintron', 'nonequal', 
                'polyA', 'dirrep', 
                'visu_str_tup'
                ])
        self.gls_info_groupby = None

        ## processing pool from multiprocessing module.
        self.mp_pool = None

        return None

    # MAIN METHOD #
    # USAGE:     call necessary functions and initiate the entire pseudogene annotation pipeline.
    # INPUT:     none.
    # OUTPUT:    none.
    # UPDATE:    none.
    # CLEAR:     none.
    def perform_flow(self):
        ## acquire and print the initiation time.
        print(f'PseudoRelics pipeline initiated. | Time: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}', flush=True)

        ## initiate the multiprocessing module and manage the multiprocessing pool.
        self.mana_mp()
        ## start the pipeline according to the pipeline mode.
        ### IGR mode: discover pseudogenes in the intergenic region.
        if self.ppl_mode == 'IGR':
            self.flow_mode_igr()
        ### GR mode: discover pseudogenes in the gene region.
        elif self.ppl_mode == 'GR':
            self.flow_mode_gr()
        ### mode not supported, exit the program.
        else:
            print(f'Pipeline mode not supported, exit the program. | Time: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}', flush=True)
            sys.exit(1)

        ## acquire and print the finish time.
        print(f'PseudoRelics pipeline finished. | Time: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}', flush=True)
        return None

    # PIPELINE METHOD #
    # USAGE:     organize the pipeline method for IGR mode.
    # INPUT:     none.
    # OUTPUT:    none.
    # UPDATE:    none.
    # CLEAR:     none.
    def flow_mode_igr(self):
        self.load_prep_align_tab_file()
        self.load_prep_gff_tab_file()
        self.elim_gls_ovlp_with_gene()
        self.load_prep_seq()
        self.extr_anal_gls_seq()
        self.filt_gls_by_qual()
        self.elim_redun_gls()
        self.assign_id_4_gls()
        self.format_output_rst()
        return None

    # PIPELINE METHOD #
    # USAGE:     organize the pipeline method for GR mode.
    # INPUT:     none.
    # OUTPUT:    none.
    # UPDATE:    none.
    # CLEAR:     none.
    def flow_mode_gr(self):
        self.load_prep_gff_tab_file()
        self.load_prep_align_tab_file()
        self.load_prep_seq()
        self.extr_anal_gls_seq()
        self.assign_id_4_gls()
        self.format_output_rst()
        return None

    # MANAGER METHOD #
    # USAGE:     initialize the multiprocessing module and manage the multiprocessing pool.
    # INPUT:     none.
    # OUTPUT:    none.
    # UPDATE:    self.mp_pool.
    # CLEAR:     none.
    def mana_mp(self):
        ## initialize the multiprocessing using forkserver method.
        mp.set_start_method('forkserver')
        ## create a multiprocessing pool object, stored in self.mp_pool.
        self.mp_pool = mp.Pool(self.proc_num)
        return None

    # TASK METHOD #
    # USAGE:     load alignment tab file, convert it to m6-like format, merge HSPs to gene-like 
    #            structures, add subject interval column, final result is a dataframe containing 
    #            gene-like structure info stored in self.gls_align_tab_df.
    # INPUT:     none.
    # OUTPUT:    none.
    # UPDATE:    self.gls_align_tab_df.
    # CLEAR:     none.
    def load_prep_align_tab_file(self):
        ## load alignment tab file.
        align_tab_df = align_ops.load_align_tab_file(self.align_tab_file, self.align_tab_type)

        ## convert the alignment tab df to standard m6-like format, add sstrand column for m6-like 
        ## format, and by the way check if the alignment tab type is supported.
        ### add sstrand column for m6-like alignment tab.
        if self.align_tab_type in align_ops.align_tab_types.m6_like:
            align_tab_df = align_ops.add_strand_4_align_tab_df(align_tab_df, 'subject')
        ### format the custom_hs alignment result to m6-like format.
        elif self.align_tab_type in align_ops.align_tab_types.exone_sugar:
            align_tab_df = align_ops.format_exone_align_tab_df_2_m6like(align_tab_df)
        ### format the exonerate sugar alignment result to m6-like format.
        else:
            print(f'Alignment tab type not supported, exit the program. | Time: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}', flush=True)
            sys.exit(1)

        ## merge HSPs of same query that share collinearity on same region and strand to a gene-like 
        ## structure.
        self.gls_align_tab_df = align_ops.merge_align_hsps_2_gls(align_tab_df, self.sbj_gap_thold, self.mp_pool)

        ## delete the original alignment tab df to save memory.
        del align_tab_df

        ## add subject interval column(called 'subject_itv') for gls align tab df.
        self.gls_align_tab_df = align_ops.add_itv_obj_4_align_tab_df(self.gls_align_tab_df)

        print(f'Alignment table file loaded and preprocessed, gene-like structure acquired. | Time: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}', flush=True)
        return None

    # TASK METHOD #
    # USAGE:     load gff file, filter out gene features, create gene location intervalset dict, final 
    #            result is a dictionary containing gene location intervalset info stored in self.gene_loc_itvset_dict.
    # INPUT:     none.
    # OUTPUT:    none.
    # UPDATE:    self.gene_loc_itvset_dict.
    # CLEAR:     none.
    def load_prep_gff_tab_file(self):
        ## load gff file.
        gff_tab_df = anno_ops.load_gff_tab_file(self.gff_tab_file)

        ## preprocess the gff df according to the pipeline mode.
        ### IGR mode: filter out cds features.
        if self.ppl_mode == 'IGR':
            #### filter out cds features from the gff df.
            gff_tab_df = anno_ops.filt_gff_tab_df_by_kwd(gff_tab_df, feat_kwd_list=self.gff_avd_feat_kwd_list)
            #### create gene location intervalset dict, the key is (region, strand) and the value is an 
            #### IntervalSet object.
            self.gene_loc_itvset_dict = anno_ops.create_feat_loc_itvset_dict(gff_tab_df)
            print(f'Gff table file loaded and preprocessed, gene location intervalset dict acquired. | Time: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}', flush=True)

        ### GR mode: filter out gene features.
        elif self.ppl_mode == 'GR':
            #### filter out gene features from the gff df.
            gff_tab_df = anno_ops.filt_gff_tab_df_by_kwd(gff_tab_df, feat_kwd_list=self.gff_avd_feat_kwd_list)
            #### parse the gff df to acquire the gene id and convert the gff df to m6-like format gls 
            #### alignment tab df.
            id_parsed_gff_df = anno_ops.parse_feat_id_4_gff_tab(gff_tab_df, self.gff_id_kwd_list, self.gff_id_ign_sfix, self.gff_id_ign_pfix)
            self.gls_align_tab_df = anno_ops.format_gff_2_gls_align_tab_df(id_parsed_gff_df)
            print(f'Gff table file loaded and preprocessed, gene-like structure acquired. | Time: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}', flush=True)

        return None

    # TASK METHOD #
    # USAGE:     eliminate the merged GLS(gene-like structure) that overlap with gene features, the 
    #            remaining GLSs are cadidate pseudogenes, which still stored in self.gls_align_tab_df.
    # INPUT:     none.
    # OUTPUT:    none.
    # UPDATE:    self.gls_align_tab_df.
    # CLEAR:     self.gls_align_tab_groupby; 
    #            self.gene_loc_itvset_dict.
    def elim_gls_ovlp_with_gene(self):
        ## group the merged gls alignment tab df by subject and sstrand.
        self.gls_align_tab_groupby = self.gls_align_tab_df.groupby(['subject', 'sstrand'])
        ## clear the gls alignment tab df to save memory.
        self.gls_align_tab_df = None

        ## create a generator object to generate task data for parallel processing.
        ## each task data is a tuple containing (group_gls_align_tab_df, group_gene_loc_itvset).
        mp_task_data_genor = self.gen_task_data_4_gls_gene_itvset()

        ## eliminate the merged gls that overlap with gene features, the elimination is done in parallel, 
        ## each group is processed by a process.
        group_gls_align_tab_df_genor = self.mp_pool.imap(self._elim_gls_ovlp_with_itvset, mp_task_data_genor)

        ## collect the results of parallel processing.
        self.gls_align_tab_df = pd.concat(group_gls_align_tab_df_genor)

        ## clear variables to save memory.
        self.gls_align_tab_groupby = None
        self.gene_loc_itvset_dict = None

        print(f'GLSs overlaping with gene eliminated, candidate GLSs acquired. | Time: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}', flush=True)
        return None

    # TASK METHOD #
    # USAGE:     load genome and pep fasta file, preprocess them to fasta item dict stored in 
    #            self.geno_fa_item_dict and self.pep_fa_item_dict which will be used in later steps.
    # INPUT:     none.
    # OUTPUT:    none.
    # UPDATE:    self.geno_fa_item_dict; 
    #            self.pep_fa_item_dict; 
    #            self.geno_fa_item_len_dict.
    # CLEAR:     none.
    def load_prep_seq(self):
        ## first load entire fasta file content as string.
        pep_fa_cont = fa_ops.load_fa_file(self.pep_fa_file)
        geno_fa_cont = fa_ops.load_fa_file(self.geno_fa_file)
        ## convert the fasta file content to fasta item dict.
        self.pep_fa_item_dict = fa_ops.format_fa_file_cont_2_dict(pep_fa_cont, rmv_seq_end=False)
        self.geno_fa_item_dict = fa_ops.format_fa_file_cont_2_dict(geno_fa_cont)
        ## acquire the length of each genome sequence and store in a dict.
        self.geno_fa_item_len_dict = fa_ops.calc_fa_item_seq_len(self.geno_fa_item_dict)

        print(f'Fasta files loaded and preprocessed, fasta item dict acquired. | Time: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}', flush=True)
        return None

    # TASK METHOD #
    # USAGE:     quantify candidate GLSs' sequence features, including: location, sequence details, alignment details, 
    #            mutation details and pseudogene type. the final result is a dataframe containing candidate pseudogene 
    #            info stored in self.gls_info_df.
    # INPUT:     none.
    # OUTPUT:    none.
    # UPDATE:    self.gls_info_df.
    # CLEAR:     self.gls_align_tab_groupby.
    def extr_anal_gls_seq(self):
        ## group gls alignment tab df by query, subject and sstrand, 
        ## for each group, the subject and query sequence will be shared.
        self.gls_align_tab_groupby = self.gls_align_tab_df.groupby(['query', 'subject', 'sstrand'])

        ## clear the gls alignment tab df to save memory.
        self.gls_align_tab_df = None

        ## create a working directory to store the exonerate temporary files prior to parallel processing, 
        ## mainly to avoid the file conflict when multiple processes create the same file.
        os.makedirs(workflow_pseudo_relics.exone_wdir, exist_ok=True)

        ## create a generator object to generate task data for parallel processing.
        ## each task data is a tuple containing (row, query_seq, gls_ext_seq, gls_ups_seq, gls_downs_seq, gls_ext_sitv).
        mp_task_data_genor = self.gen_task_data_4_gls_seq()

        ## analyze and quantify the candidate GLSs sequence features, 
        ## the analysis is done in parallel.
        row_gls_info_df_genor = self.mp_pool.imap(self._anal_gls_feat, mp_task_data_genor)

        ## collect the results of parallel processing.
        self.gls_info_df = pd.concat(row_gls_info_df_genor)

        ## clear variables to save memory.
        self.gls_align_tab_groupby = None

        self.pep_fa_item_dict = None
        self.geno_fa_item_len_dict = None

        print(f'Candidate GLSs sequence realigned to parent gene, detailed features quantified. | Time: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}', flush=True)
        return None

    # TASK METHOD #
    # USAGE:     filter the candidate GLSs by quality threshold, the remaining are high quality GLSs.
    # INPUT:     none.
    # OUTPUT:    none.
    # UPDATE:    self.gls_info_df.
    # CLEAR:     none.
    def filt_gls_by_qual(self):
        ## filter the candidate GLSs by identity and coverage threshold.
        self.gls_info_df = self.gls_info_df[(self.gls_info_df['idt'] > self.idt_thold) & (self.gls_info_df['cov'] > self.cov_thold)]

        print(f'GLSs filtered by quality, high quality GLSs set acquired. | Time: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}', flush=True)
        return None

    # TASK METHOD #
    # USAGE:     eliminate the redundancy of candidate GLSs based on selected quality mode, the remaining are 
    #            non-overlapping and high quality GLSs.
    # INPUT:     none.
    # OUTPUT:    none.
    # UPDATE:    self.gls_info_df.
    # CLEAR:     self.gls_info_groupby.
    def elim_redun_gls(self):
        ## select the quality mode to create the quality column, 
        ## the quality column is used to eliminate the redundancy of GLSs.
        if self.qlt_mode == 'idt_cov':
            self.gls_info_df['quality'] = self.gls_info_df['idt'] + self.gls_info_df['cov']
        elif self.qlt_mode == 'raws':
            self.gls_info_df['quality'] = self.gls_info_df['raws']
        elif self.qlt_mode == 'complex':
            self.gls_info_df['quality'] = self.gls_info_df['idt']*self.gls_info_df['cov']*self.gls_info_df['raws']
        else:
            print(f'Quality mode not supported, exit the program. | Time: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}', flush=True)
            sys.exit(1)

        ## reset the index of the GLS info df.
        self.gls_info_df = self.gls_info_df.reset_index(drop=True)
        ## group the GLS info df by region and strand.
        self.gls_info_groupby = self.gls_info_df.groupby(['region', 'strand'])
        ## clear the GLS info df to save memory.
        self.gls_info_df = None

        ## create a generator object to generate task data for parallel processing.
        ## each task data is a dataframe containing the GLS info of a group.
        mp_task_data_genor = self.gen_task_data_4_redun_gls()

        ## eliminate the redundancy of candidate GLSs, 
        ## the elimination is done in parallel, each group is processed by a process.
        group_gls_info_df_genor = self.mp_pool.imap(self._elim_redun_gls_4_strand, mp_task_data_genor)

        ## collect the results of parallel processing.
        self.gls_info_df = pd.concat(group_gls_info_df_genor)

        ## clear variables to save memory.
        self.gls_info_groupby = None

        print(f'GLSs redundancy eliminated by {self.qlt_mode}, final non-redundant GLSs set acquired. | Time: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}', flush=True)
        return None

    # TASK METHOD #
    # USAGE:     assign the GLS id based on the query and region for gls_info_df.
    # INPUT:     none.
    # OUTPUT:    none.
    # UPDATE:    self.gls_info_df.
    # CLEAR:     none.
    def assign_id_4_gls(self):
        ## group the GLS info df by query, and assign number suffix to each group.
        self.gls_info_df['suffix'] = self.gls_info_df.groupby('query').cumcount() + 1
        ## assign the GLS id based on the query and suffix.
        self.gls_info_df['id'] = self.gls_info_df['query'] + '_pg_' + self.gls_info_df['suffix'].astype(str)
        ## drop the suffix column to save memory.
        self.gls_info_df = self.gls_info_df.drop(columns=['suffix'])

        print(f'GLSs id assigned. | Time: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}', flush=True)
        return None

    # TASK METHOD #
    # USAGE:     format the final data and output to the output directory.
    # INPUT:     none.
    # OUTPUT:    none.
    # UPDATE:    none.
    # CLEAR:     none.
    def format_output_rst(self):
        ## create the output directory if not exists.
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir, exist_ok=True)

        ## output the final result to the output files.
        ### output default files.
        self.output_briefing_tsv_file()
        self.output_detail_align_file()
        ### output the gff file if required.
        if self.output_gff:
            self.output_briefing_gff_file()
        ### output the bed file if required.
        if self.output_bed:
            self.output_briefing_bed_file()
        ### output the sequence files if required.
        if self.output_raw_fa or self.output_cds_fa or self.output_pep_fa:
            self.output_seq_fa_files()

        print(f'Data formatted and output, result files stored in {self.output_dir}. | Time: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}', flush=True)
        return None

    # TASK DATA GENERATOR METHOD #
    # USAGE:     generate the task data for parallel processing, return a generator object, used in the 
    #            task method <elim_gls_ovlp_with_gene>.
    # INPUT:     none.
    # OUTPUT:    a generator object to generate tuple containing (group_gls_align_tab_df, group_gene_loc_itvset).
    # UPDATE:    none.
    # CLEAR:     none.
    def gen_task_data_4_gls_gene_itvset(self):
        for key, group_gls_align_tab_df in self.gls_align_tab_groupby:
            yield (group_gls_align_tab_df, self.gene_loc_itvset_dict[key])

    # TASK DATA GENERATOR METHOD #
    # USAGE:     generate the task data for parallel processing, return a generator object, used in the 
    #            task method <extr_anal_gls_seq>.
    # INPUT:     none.
    # OUTPUT:    a generator object to generate tuple containing (row, query_seq, gls_ext_seq, gls_ups_seq, 
    #            gls_downs_seq, gls_ext_sitv).
    # UPDATE:    none.
    # CLEAR:     none.
    def gen_task_data_4_gls_seq(self):
        ## iterate the groupby object to aquire the group key and group df.
        for key, group_gls_align_tab_df in self.gls_align_tab_groupby:
            ### acquire the group query, subject and sstrand.
            query_id = key[0]
            subject_id = key[1]
            sstrand = key[2]

            ### extract the subject and query sequence from fasta item dict.
            query_seq = fa_ops.extr_fa_seq_by_id(self.pep_fa_item_dict, query_id)
            subject_seq = fa_ops.extr_fa_seq_by_id(self.geno_fa_item_dict, subject_id)
            ### acquire the subject sequence length.
            subject_len = self.geno_fa_item_len_dict[subject_id]

            ### perform three tasks according to the sstrand value:
            ### 1. if sstrand is '-', convert the subject sequence to the forward complementary sequence.
            ### 2. determine the actual start and end position of the subject interval.
            ### 3. determine the function to extract the GLS's original, upstream and downstream sequence.
            if sstrand == '+':
                extr_sub_seq_func = fa_ops.extr_fwd_sub_seq_by_coord
                gls_org_sitv_start_col, gls_org_sitv_end_col = 'sstart', 'send'
                get_gls_ext_sitv_func = lambda gls_org_sitv: (max(1, gls_org_sitv[0]-self.gls_ups_ext_len), min(subject_len, gls_org_sitv[1]+self.gls_downs_ext_len))
                get_gls_ups_sitv_func = lambda gls_org_sitv: (gls_org_sitv[0]-self.gls_ups_len, gls_org_sitv[0]-1)
                get_gls_downs_sitv_func = lambda gls_org_sitv: (gls_org_sitv[1]+1, gls_org_sitv[1]+self.gls_downs_len)
            elif sstrand == '-':
                extr_sub_seq_func = fa_ops.extr_rev_compl_sub_seq_by_coord
                gls_org_sitv_start_col, gls_org_sitv_end_col = 'send', 'sstart'
                get_gls_ext_sitv_func = lambda gls_org_sitv: (max(1, gls_org_sitv[0]-self.gls_downs_ext_len), min(subject_len, gls_org_sitv[1]+self.gls_ups_ext_len))
                get_gls_ups_sitv_func = lambda gls_org_sitv: (gls_org_sitv[1]+1, gls_org_sitv[1]+self.gls_ups_len)
                get_gls_downs_sitv_func = lambda gls_org_sitv: (gls_org_sitv[0]-self.gls_downs_len, gls_org_sitv[0]-1)

            ### iterate the group df to acquire and generate the task data for each GLS.
            for _, row in group_gls_align_tab_df.iterrows():
                ### acquire the GLS's original, extended, upstream and downstream intervals.
                ### note that in m6-like format when strand is '-', the start position is the larger one, which means 
                ### the start position is the interval end position and the end position is the interval start position.
                gls_org_sitv = row[gls_org_sitv_start_col], row[gls_org_sitv_end_col]
                gls_ext_sitv = get_gls_ext_sitv_func(gls_org_sitv)
                gls_ups_sitv = get_gls_ups_sitv_func(gls_org_sitv)
                gls_downs_sitv = get_gls_downs_sitv_func(gls_org_sitv)
                ### extract the GLS's extended, upstream and downstream sequence.
                gls_ext_seq = extr_sub_seq_func(subject_seq, *gls_ext_sitv)
                gls_ups_seq = extr_sub_seq_func(subject_seq, *gls_ups_sitv)
                gls_downs_seq = extr_sub_seq_func(subject_seq, *gls_downs_sitv)
                ### generate the task data tuple and return it, use yield to save memory.
                yield (row, query_seq, gls_ext_seq, gls_ups_seq, gls_downs_seq, gls_ext_sitv)

    # TASK DATA GENERATOR METHOD #
    # USAGE:     generate task data for parallel processing, return a generator object, used in the
    #            task method <elim_redun_gls>.
    # INPUT:     none.
    # OUTPUT:    a generator object to generate dataframe containing the GLS info of a group.
    # UPDATE:    none.
    # CLEAR:     none.
    def gen_task_data_4_redun_gls(self):
        for _, group_gls_info_df in self.gls_info_groupby:
            yield group_gls_info_df

    # SUB TASK METHOD #
    # USAGE:     static method for eliminating the GLSs overlaping with gene/cds features. used in 
    #            the task method <elim_gls_ovlp_with_gene>.
    # INPUT:     a tuple containing (group_gls_align_tab_df, group_gene_loc_itvset).
    # OUTPUT:    a dataframe containing the GLSs that do not overlap with gene features.
    # UPDATE:    none.
    # CLEAR:     none.
    @staticmethod
    def _elim_gls_ovlp_with_itvset(task_data_tup):
        ## unpack the task data tuple.
        group_gls_align_tab_df, group_gene_loc_itvset = task_data_tup
        ## delete task data tuple to save memory.
        del task_data_tup
        ## eliminate the overlaping GLSs with gene features by comparing the gene location intervalset.
        group_gls_align_tab_df = group_gls_align_tab_df[group_gls_align_tab_df.apply(lambda row: not row['subject_itv'] & group_gene_loc_itvset, axis=1)]
        return group_gls_align_tab_df

    # SUB TASK METHOD #
    # USAGE:     static method for analyzing the GLS's sequence features. used in the task method 
    #            <extr_anal_gls_seq>.
    # INPUT:     a tuple containing (row, query_seq, gls_ext_seq, gls_ups_seq, gls_downs_seq, gls_ext_sitv).
    # OUTPUT:    a one-row dataframe containing the current GLS's sequence features.
    # UPDATE:    none.
    # CLEAR:     none.
    @staticmethod
    def _anal_gls_feat(task_data_tup):
        ## unpack the task data tuple.
        (
            row, 
            query_seq, gls_ext_seq, gls_ups_seq, gls_downs_seq, 
            gls_ext_sitv, 
        ) = task_data_tup
        ## realign the query pep to the GLS dna sequence to characterize the GLS's sequence features.
        ### perform the exonerate alignment using both the protein2genome model 'bf' and 'lc'.
        exone_bf_cust_stdout = align_ops.run_exone_p2g(query_seq, gls_ext_seq, 'bf', workflow_pseudo_relics.exone_prog, workflow_pseudo_relics.exone_wdir)
        exone_lc_cust_stdout = align_ops.run_exone_p2g(query_seq, gls_ext_seq, 'lc', workflow_pseudo_relics.exone_prog, workflow_pseudo_relics.exone_wdir)
        ### try to parse the exonerate alignment stdout to acquire the alignment result.
        exone_bf_cust_rst_tup = align_ops.parse_exone_p2g_cust_stdout(exone_bf_cust_stdout)
        exone_lc_cust_rst_tup = align_ops.parse_exone_p2g_cust_stdout(exone_lc_cust_stdout)
        ### if both exonerate alignment result is None, return an empty dataframe.        
        if exone_bf_cust_rst_tup is None and exone_lc_cust_rst_tup is None:
            return pd.DataFrame()
        ### otherwise, select the best exonerate alignment result based on the raw score.
        best_exone_cust_rst_tup = max(
            (exone_bf_cust_rst_tup, exone_lc_cust_rst_tup), 
            key=lambda x: x[1][6] if x is not None else 0
        )
        ### analyze the exonerate alignment result to acquire the detailed alignment information.
        align_detail_tup = align_ops.anal_exone_p2g_align_rst(best_exone_cust_rst_tup)
        ### parse to acquire the alignment detail information.
        (
            idt, cov, frag, cds, 
            _, missen_count, ins_count, del_count, prem_term_count, miss_init, miss_term, 
            frsh_count, intron_count, lintron_count, ner_count
        ) = align_detail_tup
        (best_exone_visu_str_tup, best_exone_cust_tab_tup) = best_exone_cust_rst_tup
        (
            query_len, _, 
            query_align_begin, query_align_end, target_align_begin, target_align_end, 
            align_raw_score, _, _, _
        ) = best_exone_cust_tab_tup

        ## detect the downstream polyA signal in the GLS's downstream sequence.
        hmm_polya_model_obj = motif_ops.hmm_polynt_model(dna_seq_i=gls_downs_seq)
        polynt_max_len = hmm_polya_model_obj.viterbi_detect_polynt()[1]
        if polynt_max_len >= 12:
            polya_exist = True
        else:
            polya_exist = False
        ## detect the up&downstream direct repeat in the GLS's upstream and downstream sequence.
        dirt_repe_subseq = align_ops.find_lgst_comm_subseq(gls_ups_seq, gls_downs_seq)
        if len(dirt_repe_subseq) >= 9:
            dirt_repe_exist = True
        else:
            dirt_repe_exist = False

        ## determine the type of the GLS based on the sequence features.
        if intron_count == 0 and cov < 50:
            gls_type = 'FRAG'
        elif intron_count == 0 and (polya_exist or dirt_repe_exist):
            gls_type = 'PSSD'
        else:
            gls_type = 'DUP'

        ## rectify the start and end position of the GLS's alignment using the exonerate alignment result.
        if row['sstrand'] == '+':
            rect_start = gls_ext_sitv[0] + target_align_begin
            rect_end = gls_ext_sitv[0] + target_align_end - 1
        elif row['sstrand'] == '-':
            rect_start = gls_ext_sitv[1] - target_align_end + 1
            rect_end = gls_ext_sitv[1] - target_align_begin

        ## format information of the GLS to a one-row dataframe, 
        ## which will be used to concatenate to the final result dataframe.
        row_gls_info_df = pd.DataFrame([{
            'query': row['query'], 'id': None, 'type': gls_type, 
            'region': row['subject'], 'strand': row['sstrand'], 'start': rect_start, 'end': rect_end, 
            'qstart': query_align_begin + 1, 'qend': query_align_end, 'qlen': query_len, 
            'idt': idt, 'cov': cov, 'raws': align_raw_score, 'frag': frag, 'cds': cds, 
            'miss': missen_count, 'ins': ins_count, 'del': del_count, 'preterm': prem_term_count, 
            'missinit': miss_init, 'missterm': miss_term, 
            'fsh': frsh_count, 'intron': intron_count,'lintron': lintron_count, 'ner': ner_count, 
            'polyA': polya_exist, 'dirrep': dirt_repe_exist, 
            'visu_str_tup': best_exone_visu_str_tup
        }])

        return row_gls_info_df

    # SUB TASK METHOD #
    # USAGE:     static method for eliminating the redundancy of GLSs. used in the task method 
    #            <elim_redun_gls>.
    # INPUT:     a dataframe containing the GLS info of a group.
    # OUTPUT:    a dataframe containing the non-redundant GLSs of the group.
    # UPDATE:    none.
    # CLEAR:     none.
    @staticmethod
    def _elim_redun_gls_4_strand(group_gls_info_df):

        ## sort the GLS info df by quality in descending order.
        group_gls_info_df = group_gls_info_df.sort_values(by='quality', ascending=False)
        ## reset the index of the GLS info df.
        group_gls_info_df = group_gls_info_df.reset_index(drop=True)
        ## create an IntervalTree object to store the GLS info df.
        group_gls_itvtree = IntervalTree(Interval(row['start'], row['end'], index) for index, row in group_gls_info_df.iterrows())

        ## initialize the dropped GLS index set.
        dropped_gls_index_set = set()
        ## iterate the GLS info df to eliminate the redundancy of GLSs.
        for index, row in group_gls_info_df.iterrows():
            ### if the GLS has been dropped, skip it.
            if index in dropped_gls_index_set:
                continue
            ### find the GLSs that overlap with the current GLS using the IntervalTree object.
            ovlp_gls_set = group_gls_itvtree.overlap(row['start'], row['end'])
            ### if the GLS set is not empty, eliminate the redundancy of GLSs.
            ### note to avoid the self-overlap, the current GLS is not in the ovlp_gls_set.
            if ovlp_gls_set:
                dropped_gls_index_set.update([ovlp_gls.data for ovlp_gls in ovlp_gls_set if ovlp_gls.data != index])

        ## drop the redundancy GLSs from the GLS info df using the dropped GLS index set.
        group_gls_info_df = group_gls_info_df.drop(list(dropped_gls_index_set))

        return group_gls_info_df

    # SUB TASK METHOD #
    # USAGE:     method to output the final result: briefing.tsv. used in the task method 
    #            <format_output_rst>.
    # INPUT:     none.
    # OUTPUT:    none.
    # UPDATE:    none.
    # CLEAR:     none.
    def output_briefing_tsv_file(self):
        ## initialize the output file path.
        output_file = os.path.join(self.output_dir, 'briefing.tsv')
        ## output the final result to the output file.
        self.gls_info_df.to_csv(
            output_file, 
            sep='\t', 
            index=False, 
            columns=[
                'id', 'type', 
                'query', 'qstart', 'qend', 'qlen', 
                'region', 'strand', 'start', 'end', 
                'idt', 'cov', 'raws', 'frag', 'cds', 
                'miss', 'ins', 'del', 'preterm', 
                'missinit', 'missterm', 
                'fsh', 'intron', 'lintron', 'ner', 
                'polyA', 'dirrep'
                ]
            )

        return None

    # SUB TASK METHOD #
    # USAGE:     method to output the final result: briefing.gff. used in the task method 
    #            <format_output_rst>.
    # INPUT:     none.
    # OUTPUT:    none.
    # UPDATE:    none.
    # CLEAR:     none.
    def output_briefing_gff_file(self):
        ## add columns to store the gff info.
        self.gls_info_df['source'] = 'PseudoRelics'
        self.gls_info_df['feature'] = 'pseudogene'
        self.gls_info_df['frame'] = '.'
        self.gls_info_df['attribute'] = 'ID=' + self.gls_info_df['id'] + ';locus_type=pseudogene' + ';Query=' + self.gls_info_df['query'] + ';Type=' + self.gls_info_df['type'] + ';IDT=' + self.gls_info_df['idt'].astype(str) + ';COV=' + self.gls_info_df['cov'].astype(str) + ';RAWS=' + self.gls_info_df['raws'].astype(str) + ';FRAG=' + self.gls_info_df['frag'].astype(str) + ';CDS=' + self.gls_info_df['cds'].astype(str)
        ## initialize the output file path.
        output_file = os.path.join(self.output_dir, 'briefing.gff')
        ## output the final result to the output file.
        self.gls_info_df.to_csv(
            output_file, 
            sep='\t', 
            index=False, 
            header=False, 
            columns=[
                'region', 'source', 'feature', 
                'start', 'end', 'raws', 
                'strand', 'frame', 'attribute'
            ]
        )
        ## drop the gff info columns to save memory.
        self.gls_info_df.drop(['source', 'feature', 'frame', 'attribute'], axis=1, inplace=True)

        return None

    # SUB TASK METHOD #
    # USAGE:     method to output the final result: briefing.bed. used in the task method 
    #            <format_output_rst>.
    # INPUT:     none.
    # OUTPUT:    none.
    # UPDATE:    none.
    # CLEAR:     none.
    def output_briefing_bed_file(self):
        ## add a column to store the 0-based start position of the GLS.
        self.gls_info_df['0_start'] = self.gls_info_df['start'] - 1
        ## initialize the output file path.
        output_file = os.path.join(self.output_dir, 'briefing.bed')
        ## output the final result to the output file.
        self.gls_info_df.to_csv(
            output_file, 
            sep='\t', 
            index=False, 
            header=False, 
            columns=[
                'region', '0_start', 'end', 
                'id', 'raws', 'strand', 
            ]
        )
        ## drop the 0-based start position column to save memory.
        self.gls_info_df.drop('0_start', axis=1, inplace=True)

        return None

    # SUB TASK METHOD #
    # USAGE:     method to output the final result: detail.align. used in the task method 
    #            <format_output_rst>.
    # INPUT:     none.
    # OUTPUT:    none.
    # UPDATE:    none.
    # CLEAR:     none.
    def output_detail_align_file(self):
        ## initialize the output file path and create the file.
        output_file = os.path.join(self.output_dir, 'detail.align')
        with open(output_file, 'w') as file:
            pass
        ## append the detailed alignment info to the file.
        with open(output_file, 'a') as file:
            for _, row in self.gls_info_df.iterrows():
                ### acquire the detailed information of the GLS of each row.
                id = row['id']
                query = row['query']
                qstart = row['qstart']
                qend = row['qend']
                qlen = row['qlen']
                region = row['region']
                start = row['start']
                end = row['end']
                strand = row['strand']
                idt = row['idt']
                cov = row['cov']
                raws = row['raws']
                frag = row['frag']
                cds = row['cds']
                pseudo_type = row['type']
                visu_str_tup = row['visu_str_tup']
                align_str = '\n'.join(visu_str_tup)
                ### format the detailed information to a row string.
                row_gls_detail_info = f'>{id}\t{query}|{qstart}->{qend}/{qlen}\t{region}|{start}->{end}|{strand}\tidt={idt}|cov={cov}|raws={raws}|frag={frag}|cds={cds}|type={pseudo_type}\n{align_str}\n'
                ### write the row string to the file.
                file.write(row_gls_detail_info)

        return None

    # SUB TASK METHOD #
    # USAGE:     method to output the final result: pseudo_raw.fa, pseudo_cds.fa, pseudo_pep.fa. used in
    #            task method <format_output_rst>.
    # INPUT:     none.
    # OUTPUT:    none.
    # UPDATE:    none.
    # CLEAR:     none.
    def output_seq_fa_files(self):
        ## initialize the output file path and dictinary to store the fasta item.
        if self.output_raw_fa:
            raw_fa_file = os.path.join(self.output_dir, 'pseudo_raw.fa')
            raw_fa_item_dict = {}
        if self.output_cds_fa:
            cds_fa_file = os.path.join(self.output_dir, 'pseudo_cds.fa')
            cds_fa_item_dict = {}
        if self.output_pep_fa:
            pep_fa_file = os.path.join(self.output_dir, 'pseudo_pep.fa')
            pep_fa_item_dict = {}

        ## determine the sub visual string tuple match pattern based on the cds_pep_preterm_op.
        if self.cds_pep_preterm_op == 'dp':
            sub_tup_mch_ptn = align_ops.align_rst_ptns.exone_visu_target_dna2pep_str_nonaa_preterm_ptn
        elif self.cds_pep_preterm_op == 'kp':
            sub_tup_mch_ptn = align_ops.align_rst_ptns.exone_visu_target_dna2pep_str_nonaa_ptn

        ## iterate the GLS info df to acquire the GLS's raw, cds and pep sequence.
        for _, row in self.gls_info_df.iterrows():
            ### acquire the GLS's id, query, and create the fasta item header.
            gls_id = row['id']
            query = row['query']
            fa_item_header = f'>{gls_id} {query}'

            ### extract the GLS's raw sequence if required.
            if self.output_raw_fa:
                #### extract the GLS's raw sequence from the genome fasta item dict.
                raw_seq = fa_ops.extr_fa_seq_by_loc(self.geno_fa_item_dict, row['region'], row['start'], row['end'], row['strand'])
                #### store the raw sequence to the raw fasta item dict.
                raw_fa_item_dict[gls_id] = [fa_item_header, raw_seq]

            ### acquire the GLS's visu string tuple and extract the sub visu string tuple witch drop all 
            ### matched sub-string that is not in the cds sequence.
            if self.output_cds_fa or self.output_pep_fa:
                visu_str_tup = row['visu_str_tup']
                sub_visu_str_tup = align_ops.extr_sub_visu_str_tup(visu_str_tup, visu_str_tup[2], sub_tup_mch_ptn, 'rm')

            ### extract the GLS's cds sequence if required.
            if self.output_cds_fa:
                #### extract the GLS's cds sequence from the sub visu string tuple.
                cds_seq = sub_visu_str_tup[3]
                #### store the cds sequence to the cds fasta item dict.
                cds_fa_item_dict[gls_id] = [fa_item_header, cds_seq]

            ### extract the GLS's pep sequence if required.
            if self.output_pep_fa:
                #### extract the GLS's pep sequence from the sub visu string tuple.
                pep_seq = sub_visu_str_tup[2]
                #### convert the pep sequence to short form.
                pep_seq = align_ops.conv_aa_seq_2_short(pep_seq)
                #### store the pep sequence to the pep fasta item dict.
                pep_fa_item_dict[gls_id] = [fa_item_header, pep_seq]

        ## output the raw, cds and pep fasta item dict to the corresponding fasta files.
        if self.output_raw_fa:
            fa_ops.output_fa_item_dict_2_file(raw_fa_item_dict, raw_fa_file)
        if self.output_cds_fa:
            fa_ops.output_fa_item_dict_2_file(cds_fa_item_dict, cds_fa_file)
        if self.output_pep_fa:
            fa_ops.output_fa_item_dict_2_file(pep_fa_item_dict, pep_fa_file)

        return None
