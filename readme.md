Required languages:
	R-3.5
	Python3.6 with umap-learn installed

Required R packages:
	Seurat
	ggplot2
	dplyr
	ggridges


***Only files from bash_scripts/ directory should be run. Before running them, cd to project root directory***


Files in job_scripts/ are for scheduling jobs:
	initializer.py - python file which generates and submits the job
	mc_init.sh - requests resources and starts mc R script
	mc_plot_init.sh - requests resources and starts fc_plots R script


Files in scripts/ directory are R scripts for analysis:
	fc_plots.R - filtered cell plots functions
	mc_functions - method comparison functions
	mc_plot.R - filtered cells plot initializer
	mc.R - method comparison initializer
	readers.R - dataset reading functions
	settings.R - file with settings


Files in bash_scripts/ directory:
	greperr.sh: searches if there are any errors in job logs.
		Usage: grepper.sh -j $jid, where $jid is the job id

	shceduler.sh: loads Python and UGER and runs initializer.py 
		Usage: shceduler.sh $script $dataset **args
			script: the script name, available options: mc (method comparison), mc_plot (filtered cells plots), mc+mc_plot(both comparison and plots, plots will be submitted with -hold_jid)
			
			dataset: name of the dataset, available options: ebi, ebi_tm, mca, tm, ts24, ts30, other
			**args - any additional args passed to qsub (example -hold_jid)


Saving Files
	If you want to save seurat objects as RDS files after script completion you should add the following line into "files_to_save.txt" (all lowercase):
	project-tissue-res-method-param


R script settings (scripts/settings.R):

	data.dir - absolute path of data location
	
	output.dir - absolute path of output location
	
	source.dir - is used by fc_plots, should be the same as output.dir	

	cells.filter - filtering parameter that indicates minimum number of cells in which gene should be present (3 by default)

	genes.filter - filtering parameter that sets the minimum gene cutoff (100 by default)

	save.res.1 - TRUE/FALSE. Indicates whether all files of resolution 1 and method none should be saved (FALSE by default)

Logs:
	all job logs are located in the logs/ directory. Name structure is %script_name%.o%job_id%.%task_id%

	file submissions.txt contains information about all submitted jobs in the following format: submission_time|job.id|tasks|job_name 


Output directory structure:
	Location of output dir is set in R script settings. It has the following structure:

	output/robjs/script_name/*.rds - location where all RDS objects are saved

	output/script_name/tissue/res-method-param/!filtered.csv - list of filtered cells (if file is missing no cells were filtered)

	output/script_name/tissue/res-method-param/!markers.csv - list of DE genes

	output/script_name/tissue/res-method-param/!clusters.csv - clusters summary

	output/script_name/tissue/res-method-param/!cells.csv - cells summary

	output/script_name/tissue/res-method-param/*.pdf and output/script_name/tissue/filtered_cells_plots/*.pdf - all method specific plots

	output/script_name/tissue/score_summary.csv - average cell type scores for all methods

	output/script_name/stats_summary.csv - table of QC statistics for each cell.







