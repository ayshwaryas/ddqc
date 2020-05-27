#!/bin/bash
cd ~/method_comparison
bash_scripts/scheduler.sh mc+mc_plot ebi
bash_scripts/scheduler.sh mc+mc_plot mca
bash_scripts/scheduler.sh mc+mc_plot tm
bash_scripts/scheduler.sh mc+mc_plot ts30
bash_scripts/scheduler.sh mc+mc_plot other
bash_scripts/scheduler.sh mc+mc_plot ebi_tm
bash_scripts/scheduler.sh mc+mc_plot other_10X
bash_scripts/scheduler.sh mc+mc_plot PanglaoDB