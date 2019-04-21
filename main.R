###############################################################################
# main.R                                                                      #
#                                                                             #
# Project 019 - Convert 017 to Galaxy Module                                  #
#                                                                             #
# Description: Combine projects 017 and 016 to form a Galaxy module           #
#                                                                             #
# Started: 2019-02-12                                                         #
#                                                                             #
# Version: (1.0)                                                              #
#                                                                             #
###############################################################################

###############################################################################
####### Load libraries etc.
###############################################################################
source("functions.R")
#source("data.R") # Loads list_of_collections

###############################################################################
# Display plots
###############################################################################
plots_for_paper <- Plots_for_Paper$new()
plots_for_paper$create_plots_for_paper(include_main = TRUE , finalize=FALSE)
#plots_for_paper$create_plots_for_paper(include_main = FALSE, finalize=TRUE)

###############################################################################
# Goals
###############################################################################
# One or more Galaxy modules