Sequential Matching Simulations
==============================

This repository contains information to replicate figures and tables in the paper "Matching on-the-fly in Sequential Experiments for Higher Power and Effciency."

Section 3
---------

This section contained simulations for three scenarios: non-linear, linear, and zero effects for many different pair matching methods and competitors. 

To replicate, run all_match_sims.R in an R environment. It's faster using grid computing (see the bash file "run_sim_on_grid.sh"). We have included our raw ourput in csv and xlsx formats for generating tables 3 and 4. 

After simulations are finished, you can aggregate the power rows into a file "all_powers.csv" and run "plots_for_paper.R" to generate figure 1.

Section 4
---------

This section contained simulations to illustrate sequential matching's effectiveness in historic RCT data.

### Section 4.1 ###

Tables 5 / 6 can be replicated by running "sec_4_1_analyze_kapcha_data.R" and 
un/commenting lines 372/378.

### Section 4.2 ###

Tables 7 / 8 can be replicated by running "sec_4_2_analyze_cceb_data.R" and 
un/commenting lines 282/288. The raw sas file is semi-private. The authors can furnish it via email.
