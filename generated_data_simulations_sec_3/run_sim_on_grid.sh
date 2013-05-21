#!/bin/bash

#$ -j y
#$ -N seq_rct_matching_power_simulations
############ The final task number must be MANUALLY changed based on how big the beta_mat permutations is
#$ -t 1-6
#$ -q all.q

echo "starting R for seq_rct_matching power simulation #$SGE_TASK_ID"
R --no-save --args iter_num=$SGE_TASK_ID < all_match_sims.R

