#——————————————————————————————————————————————————————————————————————————————————#
# PROJ:     PseudoRelics                                                           #
# VER:      1.0                                                                    #
# TYPE:     Command Line Interface                                                 #
# DESC:     The command line interface for PseudoRelics pipeline                   #
# AUTHOR:                                                                  #
# AFFIL:    Wuhan Botanical Garden, Chinese Academy Of Sciences                    #
# E-MAIL:                       #
# DATE:     2024-6-18                                                              #
#——————————————————————————————————————————————————————————————————————————————————#

import sys, argparse, textwrap
from ..wkfls import wkfl_pseudo_relics

# Command Line Interface Class #
# Usage:    recieve command line arguments and perform the pipeline workflow
# Input:    none
# Output:   none
class scpt_main:
    # INIT METHOD #
    # Usage:    initialize the command line interface, recieve command line arguments and perform the pipeline workflow.
    # Input:    none
    # Output:   none
    # Update:   none
    # Clear:    none
    def __init__(self):

        ## set script information
        scpt_info = textwrap.dedent('''
            Script description: Discover pseudogenes in the intergenic region from alignment result and gff file.
            Usage:
                python PseudoRelics.py -ppl_mode IGR -align alignment.m6 -gff genome.gff -geno genome.fasta -pep pep.fasta -out output_dir
            ''')

        ## set command line arguments
        self.parser = argparse.ArgumentParser(description=scpt_info)
        self.parser.add_argument('-proc_num', dest='proc_num', type=int, default=1, help='number of processes, default is 1')
        self.parser.add_argument('-ppl_mode', dest='ppl_mode', default='IGR', choices=['IGR', 'GR'], help='pipeline mode, IGR for intergenic region, GR for gene region')
        self.parser.add_argument('-align', dest='align', required=True, help='alignment result file')
        self.parser.add_argument('-gff', dest='gff', required=True, help='gff file at least containing gene feature')
        self.parser.add_argument('-geno', dest='geno', required=True, help='genome fasta file')
        self.parser.add_argument('-pep', dest='pep', required=True, help='pep fasta file')
        self.parser.add_argument('-out', dest='out', required=True, help='output directory')
        self.parser.add_argument('-out_gff', dest='out_gff', action='store_true', default=False, help='gff file output flag, default is False')
        self.parser.add_argument('-out_bed', dest='out_bed', action='store_true', default=False, help='bed file output flag, default is False')
        self.parser.add_argument('-out_raw_fa', dest='out_raw_fa', action='store_true', default=False, help='raw fasta file output flag, default is False')
        self.parser.add_argument('-out_cds_fa', dest='out_cds_fa', action='store_true', default=False, help='cds fasta file output flag, default is False')
        self.parser.add_argument('-out_pep_fa', dest='out_pep_fa', action='store_true', default=False, help='pep fasta file output flag, default is False')
        self.parser.add_argument('-exone', dest='exone', help='exonerate program path')
        self.parser.add_argument('-wdir', dest='wdir', help='working directory')
        self.parser.add_argument('-align_type', dest='align_type', default='m6', choices=['m6', 'cust_hs', 'exone_sugar'], help='alignment result file type, default is m6')
        self.parser.add_argument('-gap', dest='gap', type=int, default=5000, help='subject gap threshold for merging HSPs, default is 5000')
        self.parser.add_argument('-gff_feat', dest='gff_feat', nargs='+', help='gff feature keyword list')
        self.parser.add_argument('-gff_id', dest='gff_id', nargs='+', default=['ID'], help='gff feature id keyword list, default is ID')
        self.parser.add_argument('-gff_ign_sfix', dest='gff_ign_sfix', action='store_true', default=False, help='gff feature id suffix ignore flag, default is False')
        self.parser.add_argument('-gff_ign_pfix', dest='gff_ign_pfix', action='store_true', default=False, help='gff feature id prefix ignore flag, default is False')
        self.parser.add_argument('-ups', dest='ups', type=int, default=1000, help='length for GLS upstream sequence extraction, default is 1000')
        self.parser.add_argument('-dws', dest='dws', type=int, default=1000, help='length for GLS downstream sequence extraction, default is 1000')
        self.parser.add_argument('-ups_ext', dest='ups_ext', type=int, default=100, help='upstream extension length for GLS sequence extraction, default is 100')
        self.parser.add_argument('-dws_ext', dest='dws_ext', type=int, default=100, help='downstream extension length for GLS sequence extraction, default is 100')
        self.parser.add_argument('-idt', dest='idt', type=float, default=60, help='minimum identity for pseudogene')
        self.parser.add_argument('-cov', dest='cov', type=float, default=40, help='minimum coverage for pseudogene')
        self.parser.add_argument('-qlt_mode', dest='qlt_mode', default='idt_cov', choices=['idt_cov', 'raws', 'complex'], help='quality mode for pseudogene redundancy elimination, default is ident_cov')
        self.parser.add_argument('-pterm', dest='pterm', default='dp', choices=['dp', 'kp'], help='premature termination option, dp for drop, kp for keep, default is dp')

        ## check command line arguments
        try:
            self.args = self.parser.parse_args()
        except (argparse.ArgumentTypeError, argparse.ArgumentError):
            self.parser.print_help()
            sys.exit(1)

        ## acquire command line arguments
        self.proc_num = self.args.proc_num
        self.ppl_mode = self.args.ppl_mode
        self.align_tab_file = self.args.align
        self.gff_tab_file = self.args.gff
        self.geno_fa_file = self.args.geno
        self.pep_fa_file = self.args.pep
        self.output_dir = self.args.out
        self.output_gff = self.args.out_gff
        self.output_bed = self.args.out_bed
        self.output_raw_fa = self.args.out_raw_fa
        self.output_cds_fa = self.args.out_cds_fa
        self.output_pep_fa = self.args.out_pep_fa
        self.exone_prog = self.args.exone
        self.exone_wdir = self.args.wdir
        self.align_tab_type = self.args.align_type
        self.sbj_gap_thold = self.args.gap
        self.gff_feat_list = self.args.gff_feat
        self.gff_id_kwd_list = self.args.gff_id
        self.gff_ign_sfix = self.args.gff_ign_sfix
        self.gff_ign_pfix = self.args.gff_ign_pfix
        self.gls_ups_len = self.args.ups
        self.gls_downs_len = self.args.dws
        self.gls_ups_ext_len = self.args.ups_ext
        self.gls_downs_ext_len = self.args.dws_ext
        self.idt = self.args.idt
        self.cov = self.args.cov
        self.qlt_mode = self.args.qlt_mode
        self.pterm = self.args.pterm

        ## perform the pipeline workflow
        self.workflow = wkfl_pseudo_relics.workflow_pseudo_relics(self.proc_num, self.ppl_mode, self.align_tab_file, self.gff_tab_file, self.geno_fa_file, self.pep_fa_file, self.output_dir, self.output_gff, self.output_bed, self.output_raw_fa, self.output_cds_fa, self.output_pep_fa, self.exone_prog, self.exone_wdir, self.align_tab_type, self.sbj_gap_thold, self.gff_feat_list, self.gff_id_kwd_list, self.gff_ign_sfix, self.gff_ign_pfix, self.gls_ups_len, self.gls_downs_len, self.gls_ups_ext_len, self.gls_downs_ext_len, self.idt, self.cov, self.qlt_mode, self.pterm).perform_flow()

        return None
