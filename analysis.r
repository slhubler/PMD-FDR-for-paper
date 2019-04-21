###############################################################################
# analysis.R                                                                  #
#                                                                             #
# Project 019 - Convert 017 to Galaxy Module                                  #
#                                                                             #
# Description: Separate analyses based on main                                #
#                                                                             #
# Started: 2019-02-12                                                         #
#                                                                             #
# Version: (1.0)                                                              #
#                                                                             #
###############################################################################

###############################################################################
# xxxxxx
###############################################################################

data_processor <- Data_Processor$new(p_info = Data_Object_Info_737_two_step$new())
data_processor$alpha$ensure()
head(data_processor$alpha$df)
