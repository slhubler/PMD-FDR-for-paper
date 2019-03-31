###############################################################################
# main.R                                                                      #
#                                                                             #
# Project 020 - Create GitHub module for paper                                #
#                                                                             #
# Description: Analyze data and produce figures for paper                     #
#                                                                             #
# Started: 2019-02-28                                                         #
#                                                                             #
# Version: (1.1)                                                              #
#                                                                             #
###############################################################################

# NEW REQUIREMENTS
#   - All analysis but be within the same directory and self-contained
#   - Main should only contain what is necessary to run analysis
#   - ... including all images used in paper

###############################################################################
####### Load libraries etc.
###############################################################################

#source("016 - functions.R")
#source("016 - Actual_Distributions.R")
source("functions.R")

#source(file.path("..","Project 004 - myRUnit",  "RCRUnit.R"))

###############################################################################
####### 2018-03-21 Create full analysis structure for PMD paper
###############################################################################

###############################################################################
# Reprocess data
###############################################################################

# 
# process_parameters = list(load_collection = "skip",
#                           data            = "load",
#                           clean_data      = "skip",
#                           grouping        = "make",
#                           densities       = "make",
#                           fdr_approx      = "make",
#                           i_fdr           = "make",
#                           ci_results      = "skip",
#                           one_percent_FDR = "skip",
#                           save_collection = "save")
# 
# 
# # 
# # apply_process_to_all_collections(list_of_data_collections = list_of_data_collections,
# #                                  process_parameters = process_parameters,
# #                                  over_write = FALSE )
# 
# #result <- apply_process_to_one_collection(data_collection = data_collection_oral_737_NS_combined, process_parameters = process_parameters, over_write = TRUE)
# create_full_data_collection_all(list_of_data_collections = list_of_data_collections, over_write = TRUE)

###############################################################################
####### Finalize plots
###############################################################################

plots_for_paper <- Plots_for_Paper$new()
plots_for_paper$create_data_collections(list_of_data_collections=list_of_data_collections)
#plots_for_paper$create_plots_for_paper(include_main = TRUE , finalize=FALSE)
plots_for_paper$create_plots_for_paper(include_main = FALSE, finalize=TRUE)

