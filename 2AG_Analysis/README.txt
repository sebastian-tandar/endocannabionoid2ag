##########################
# S. T. Tandar; 2024 #
# Analysis package for 2AG signaling 
#  > Semi-mechanistic model depicting 2AG accumulation and release from N2A cells + 
#    uptake by and signal generation in HEK293T.GRABeCB2.0 cells
##########################

# DIRECTORIES #
# 1_RawData: (Raw data was not included in this analysis package. Preprocessed data are available on the next directory.)
# 2_CreatedData: Datasets, figures, and model fit files generated during the analysis. Files produced by one script and used in another was ordered in sequence by the first letter of the file name (i.e. A, B, and so on).
# 3_Scripts: R scripts used to perform data preprocessing and analysis. Scripts was ordered in the sequence that is relevant for the analysis by the first letter of the file name (i.e. A, B...) and is connected to the letters used to name the output files under '2_CreatedData' (script starting with 'A_~' produces an output with filename starting with 'A_~'). Scripts starting with 'FS' are concluding files used to: 1) one of the final figures, and 2) compile information needed for model simulation in the shinyApp (see next point). 
# Signaling2AG: This folder contains the final model and simulation package used to run the shinyApp. The application can be used to further explore the model (modify conditions, modify system parameters, etc...)
###############

# RUNNING THE SCRIPTS #
# To run the script, the following prerequisites are needed:
#    > Installation of R and relevant packages ('dplyr', 'reshape2', 'rlist', 'readxl', 'pso', 'nlmixr2', 'rxode2', 'ggplot2', 'cowplot', 'scales', 'shiny')
#    > Directory adjustment - Path to the main '2AG_Analysis' directory need to be mentioned on L.2 of all scripts as 'mainwd' (except for the shinyApp under 'Signaling2AG' folder; see next point). For example, if the folder is put under 'C:\Users\GentleReader\Documents\2AG_Analysis', then you need to change:

				mainwd <- "path_to_folder/"

to

				mainwd <- "C:\\Users\\Reader\\Documents\\2AG_Analysis\\"


		>> If you copy-paste the path directly from windows, it will give something like 'C:\Users\Reader\Documents\2AG_Analysis'. Using this as is will give an error.
			>>> Change the diagonal slash '\' to EITHER '\\' or '/'
			>>> Don't forget to add the extra slash(es) at the end


# ShinyApp
#	> The shiny app can be opened by running either 'ui.R' or 'server.R'. 
#	> Path adjustment for the ShinyApp needs to be performed on L10 of 'server.R' script.