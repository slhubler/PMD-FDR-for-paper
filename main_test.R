###############################################################################
# main_test.R                                                                 #
#                                                                             #
# Project 019 - Convert 017 to Galaxy Module                                  #
#                                                                             #
# Description: Main code for testing project entities                         #
#                                                                             #
# Started: 2019-02-15                                                         #
#                                                                             #
# Version: (1.0)                                                              #
#                                                                             #
###############################################################################

library(RUnit)

source(file.path("..","Project 004 - myRUnit",  "RCRUnit.R"))
source(file.path("..","Project 005 - ArgParser","ArgParser.R"))

source("functions.R")
source("test.R")

############################ TESTING
testCase <- Test_Analyze_Collection$new()
cat(testCase$completeTest()$report())
#testCase <- Test_List_of_Collections$new()
#cat(testCase$completeTest()$report())
#testCase <- Test_Collection_Info$new()
#cat(testCase$completeTest()$report())
testCase <- Test_Data_Object$new()
cat(testCase$completeTest()$report())
#testCase <- Test_Data_Collection$new()
#cat(testCase$completeTest()$report())
testCase <- Test_Data_Object_Info$new()
cat(testCase$completeTest()$report())
testCase <- Test_Data_Processor$new()
cat(testCase$completeTest()$report())
