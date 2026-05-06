#——————————————————————————————————————————————————————————————————————————————————#
# PROJ:     ToolKits                                                               #
# VER:      1.0                                                                    #
# TYPE:     Module                                                                 #
# DESC:     A module for motif operations.                                         #
# AUTHOR:                                                                 #
# AFFIL:    Wuhan Botanical Garden, Chinese Academy Of Sciences                    #
# E-MAIL:                       #
# DATE:     2024-6-14                                                              #
#——————————————————————————————————————————————————————————————————————————————————#

import numpy as np

# MODEL CLASS #
# USAGE:     create a hidden markov model for poly-nucleotide detection.
# INPUT:     nt_i: the poly-nucleotide to be detected;
#            dna_seq_i: the input dna sequence.
# OUTPUT:    see the methods.
class hmm_polynt_model:

    # INIT METHOD #
    # USAGE:     acquire the input parameters and initialize the model variables.
    # INPUT:     nt_i: the poly-nucleotide to be detected; 
    #            dna_seq_i: the input dna sequence.
    # OUTPUT:    none.
    # UPDATE:    self.nt; 
    #            self.dna_seq; 
    #            self.hdn_states_tup; 
    #            self.obser_states_tup; 
    #            self.init_prob_ary;
    #            self.trans_prob_mtx; 
    #            self.emis_prob_mtx;
    #            self.obser_ary.
    # CLEAR:     none.
    def __init__(self, nt_i='A', dna_seq_i=None):
        ## store the input parameters.
        self.nt = nt_i
        self.dna_seq = dna_seq_i

        ## initial hmm model variables.
        ### hidden states and observation states.
        self.hdn_states_tup = (f'{self.nt}', 'N')
        self.obser_states_tup = (f'{self.nt}', 'X')
        ### initial probabilities.
        self.init_prob_ary = np.array([0.01, 0.99])
        ### transition probabilities
        self.trans_prob_mtx = np.array([
            [0.99, 0.01],
            [0.01, 0.99]
        ])
        ### emission probabilities.
        self.emis_prob_mtx = np.array([
            [0.99, 0.01],
            [0.25, 0.75]
        ])
        ### format dna sequence to observation array.
        self.obser_ary = self._dna_seq_2_obser_ary(self.nt, self.dna_seq)

    # SUB METHOD #
    # USAGE:     format the dna sequence to observation array.
    # INPUT:     nt: the poly-nucleotide to be detected;
    #            dna_seq: the input dna sequence.
    # OUTPUT:    obser_ary: the observation array.
    # UPDATE:    none.
    def _dna_seq_2_obser_ary(self, nt, dna_seq):
        ## initialize the observation sequence map.
        obser_seq_map = {f'{nt}': 0, 'X': 1}
        ## format the dna sequence to observation sequence.
        obser_seq = ''.join(f'{nt}' if base == f'{nt}' else 'X' for base in dna_seq)
        ## format the observation sequence to observation array.
        obser_ary = np.array([obser_seq_map[obser_state] for obser_state in obser_seq])
        ## return the observation array.
        return obser_ary

    # TASK METHOD #
    # USAGE:     detect the poly-nucleotide using viterbi algorithm.
    # INPUT:     none.
    # OUTPUT:    a tuple of four elements:
    #            1. polynt_exist: a boolean value to indicate the existence of the poly-nucleotide;
    #            2. polynt_max_len: the maximum length of the poly-nucleotide;
    #            3. best_hdn_state_prob: the probability of the best hidden state sequence;
    #            4. best_hdn_state_seq: the best hidden state sequence.
    # UPDATE:    none.
    def viterbi_detect_polynt(self):
        ## check the observation array, return default values if empty.
        if len(self.obser_ary) == 0:
            return (False, 0, 0, '')

        ## initialize the viterbi probability matrix and path probability matrix.
        vtb_prob_mtx = np.zeros((len(self.obser_ary), len(self.init_prob_ary)))
        path_prob_mtx = np.zeros((len(self.obser_ary), len(self.init_prob_ary)), dtype=int)
        ## calculate the viterbi probability matrix and path probability matrix.
        ### viterbi probability matrix is the probability of the best hidden state sequence ending 
        ### with the hidden state at the current time step.
        vtb_prob_mtx[0, :] = self.init_prob_ary * self.emis_prob_mtx[:, self.obser_ary[0]]
        ### path probability matrix is the hidden state index of the previous time step that leads to 
        ### the hidden state at the current time step with the maximum probability.
        for time_step in range(1, len(self.obser_ary)):
            for hdn_state in range(len(self.init_prob_ary)):
                prob = vtb_prob_mtx[time_step-1] * self.trans_prob_mtx[:, hdn_state] * self.emis_prob_mtx[hdn_state, self.obser_ary[time_step]]
                vtb_prob_mtx[time_step, hdn_state] = np.max(prob)
                path_prob_mtx[time_step, hdn_state] = np.argmax(prob)

        ## backtracking to get the best hidden state sequence.
        ### initialize the best hidden state sequence array.
        best_path_ary = np.zeros(len(self.obser_ary), dtype=int)
        ### get the best hidden state sequence.
        best_path_ary[-1] = np.argmax(vtb_prob_mtx[-1])
        for time_step in reversed(range(1, len(self.obser_ary))):
            best_path_ary[time_step-1] = path_prob_mtx[time_step, best_path_ary[time_step]]

        ## get the detection results.
        best_hdn_state_prob = np.max(vtb_prob_mtx[-1])
        best_hdn_state_seq = ''.join(self.hdn_states_tup[i] for i in best_path_ary)
        polynt_exist = True if self.nt in best_hdn_state_seq else False
        polynt_max_len = max(len(seq) for seq in best_hdn_state_seq.split('N'))

        ## return the detection results.
        return (polynt_exist, polynt_max_len, best_hdn_state_prob, best_hdn_state_seq)

    # TASK METHOD #
    # USAGE:     calculate the forward probability of the given observation sequence under the model.
    # INPUT:     none.
    # OUTPUT:    a float value of the forward probability.
    # UPDATE:    none.
    def forward_calc_prob(self):
        ## initialize the forward probability matrix.
        fwd_prob_mtx = np.zeros((len(self.obser_ary), len(self.init_prob_ary)))
        ## calculate the forward probability matrix.
        fwd_prob_mtx[0, :] = self.init_prob_ary * self.emis_prob_mtx[:, self.obser_ary[0]]

        ## calculate the forward probability.
        for time_step in range(1, len(self.obser_ary)):
            for hdn_state in range(len(self.init_prob_ary)):
                fwd_prob_mtx[time_step, hdn_state] = np.sum(fwd_prob_mtx[time_step-1] * self.trans_prob_mtx[:, hdn_state] * self.emis_prob_mtx[hdn_state, self.obser_ary[time_step]])

        ## return the forward probability.
        return np.sum(fwd_prob_mtx[-1, :])
