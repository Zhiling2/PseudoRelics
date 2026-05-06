#——————————————————————————————————————————————————————————————————————————————————#
# PROJ:     ToolKits                                                               #
# VER:      1.0                                                                    #
# TYPE:     Resource                                                               #
# DESC:     A resource file to store the bioinformatics tables.                    #
# AUTHOR:                                                                  #
# AFFIL:    Wuhan Botanical Garden, Chinese Academy Of Sciences                    #
# E-MAIL:                       #
# DATE:     2024-6-14                                                              #
#——————————————————————————————————————————————————————————————————————————————————#

# RESOURCES CLASS #
class genetic_code_tabs:
        
        std_code_tab = {
            'TTT':'F', 'TCT':'S', 'TAT':'Y', 'TGT':'C', 
            'TTC':'F', 'TCC':'S', 'TAC':'Y', 'TGC':'C', 
            'TTA':'L', 'TCA':'S', 'TAA':'*', 'TGA':'*', 
            'TTG':'L', 'TCG':'S', 'TAG':'*', 'TGG':'W', 

            'CTT':'L', 'CCT':'P', 'CAT':'H', 'CGT':'R', 
            'CTC':'L', 'CCC':'P', 'CAC':'H', 'CGC':'R', 
            'CTA':'L', 'CCA':'P', 'CAA':'Q', 'CGA':'R', 
            'CTG':'L', 'CCG':'P', 'CAG':'Q', 'CGG':'R', 

            'ATT':'I', 'ACT':'T', 'AAT':'N', 'AGT':'S', 
            'ATC':'I', 'ACC':'T', 'AAC':'N', 'AGC':'S', 
            'ATA':'I', 'ACA':'T', 'AAA':'K', 'AGA':'R', 
            'ATG':'M', 'ACG':'T', 'AAG':'K', 'AGG':'R', 

            'GTT':'V', 'GCT':'A', 'GAT':'D', 'GGT':'G', 
            'GTC':'V', 'GCC':'A', 'GAC':'D', 'GGC':'G', 
            'GTA':'V', 'GCA':'A', 'GAA':'E', 'GGA':'G', 
            'GTG':'V', 'GCG':'A', 'GAG':'E', 'GGG':'G', 

            'GCR':'A', 'CGR':'R', 'GGR':'G', 'CTR':'L', 'CCR':'P', 'TCR':'S', 'ACR':'T', 'GTR':'V', 
            'GCY':'A', 'CGY':'R', 'GGY':'G', 'CTY':'L', 'CCY':'P', 'TCY':'S', 'ACY':'T', 'GTY':'V', 
            'GCM':'A', 'CGM':'R', 'GGM':'G', 'CTM':'L', 'CCM':'P', 'TCM':'S', 'ACM':'T', 'GTM':'V', 
            'GCK':'A', 'CGK':'R', 'GGK':'G', 'CTK':'L', 'CCK':'P', 'TCK':'S', 'ACK':'T', 'GTK':'V', 
            'GCS':'A', 'CGS':'R', 'GGS':'G', 'CTS':'L', 'CCS':'P', 'TCS':'S', 'ACS':'T', 'GTS':'V', 
            'GCW':'A', 'CGW':'R', 'GGW':'G', 'CTW':'L', 'CCW':'P', 'TCW':'S', 'ACW':'T', 'GTW':'V', 
            'GCH':'A', 'CGH':'R', 'GGH':'G', 'CTH':'L', 'CCH':'P', 'TCH':'S', 'ACH':'T', 'GTH':'V', 
            'GCB':'A', 'CGB':'R', 'GGB':'G', 'CTB':'L', 'CCB':'P', 'TCB':'S', 'ACB':'T', 'GTB':'V', 
            'GCV':'A', 'CGV':'R', 'GGV':'G', 'CTV':'L', 'CCV':'P', 'TCV':'S', 'ACV':'T', 'GTV':'V', 
            'GCD':'A', 'CGD':'R', 'GGD':'G', 'CTD':'L', 'CCD':'P', 'TCD':'S', 'ACD':'T', 'GTD':'V', 
            'GCN':'A', 'CGN':'R', 'GGN':'G', 'CTN':'L', 'CCN':'P', 'TCN':'S', 'ACN':'T', 'GTN':'V', 

            'ATY':'I', 
            'ATW':'I', 
            'ATM':'I', 
            'ATH':'I', 

            'AGR':'R', 
            'MGG':'R', 
            'MGA':'R', 
            'AAY':'N', 
            'GAY':'D', 
            'TGY':'C', 
            'GAR':'E', 
            'CAR':'Q', 
            'CAY':'H', 
            'TTR':'L', 
            'YTA':'L', 
            'YTG':'L', 
            'AAR':'K', 
            'TTY':'F', 
            'AGY':'S', 
            'TAY':'Y', 
            'TAR':'*', 
            'TRA':'*' 
        }

        dna_bp_tab = {
            'A':'T', 'T':'A', 'G':'C', 'C':'G', 
            'R':'Y', 'Y':'R', 'M':'K', 'K':'M', 
            'H':'D', 'D':'H', 'B':'V', 'V':'B', 
            'S':'S', 'W':'W', 'N':'N'
        }

        rna_bp_tab = {
            'A':'U', 'U':'A', 'G':'C', 'C':'G', 
            'R':'Y', 'Y':'R', 'M':'K', 'K':'M', 
            'H':'D', 'D':'H', 'B':'V', 'V':'B', 
            'S':'S', 'W':'W', 'N':'N'
        }

        dna2rna_bp_tab = {
            'A':'U', 'T':'A', 'G':'C', 'C':'G', 
            'R':'Y', 'Y':'R', 'M':'K', 'K':'M', 
            'H':'D', 'D':'H', 'B':'V', 'V':'B', 
            'S':'S', 'W':'W', 'N':'N'
        }

        rna2dna_bp_tab = {
            'A':'T', 'U':'A', 'G':'C', 'C':'G', 
            'R':'Y', 'Y':'R', 'M':'K', 'K':'M', 
            'H':'D', 'D':'H', 'B':'V', 'V':'B', 
            'S':'S', 'W':'W', 'N':'N'
        }

        dna2cdna2rna_bp_tab = {
            'A':'A', 'T':'U', 'G':'G', 'C':'C', 
            'R':'R', 'Y':'Y', 'M':'M', 'K':'K', 
            'H':'H', 'D':'D', 'B':'B', 'V':'V', 
            'S':'S', 'W':'W', 'N':'N'
        }

        f2s_aa_conv_tab = {
            'Ala':'A', 'Cys':'C', 'Asp':'D', 'Glu':'E', 
            'Phe':'F', 'Gly':'G', 'His':'H', 'Lys':'K', 
            'Ile':'I', 'Leu':'L', 'Met':'M', 'Asn':'N', 
            'Pro':'P', 'Gln':'Q', 'Arg':'R', 'Ser':'S', 
            'Thr':'T', 'Val':'V', 'Tyr':'Y', 'Trp':'W', 
            '***':'*', 'Unk':'X', '---':'-', '<->':'-'
        }

        s2f_aa_conv_tab = {
            'A':'Ala', 'C':'Cys', 'D':'Asp', 'E':'Glu', 
            'F':'Phe', 'G':'Gly', 'H':'His', 'K':'Lys', 
            'I':'Ile', 'L':'Leu', 'M':'Met', 'N':'Asn', 
            'P':'Pro', 'Q':'Gln', 'R':'Arg', 'S':'Ser', 
            'T':'Thr', 'V':'Val', 'Y':'Tyr', 'W':'Trp', 
            '*':'***', 'X':'Unk', '-':'---'
        }
