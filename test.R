###############################################################################
# test.R                                                                      #
#                                                                             #
# Project 019 - Convert 017 to Galaxy Module                                  #
#                                                                             #
# Description: Test cases for this project                                    #
#                                                                             #
# Started: 2019-02-14                                                         #
#                                                                             #
# Version: (1.0)                                                              #
#                                                                             #
###############################################################################

###############################################################################
# Tests
###############################################################################
FAKE_VALUE = 42

#####
# Mock_Data_Object
#####
Mock_Data_Object <- setRefClass("Mock_Data_Object", 
                                contains = "Data_Object",
                                fields = list(m_has_been_loaded = "logical",
                                              m_ensure_has_run  = "logical") )
Mock_Data_Object$methods(
  initialize = function(){
    callSuper()
    set_loaded(FALSE)
    set_ensure_has_run(FALSE)
  },
  m_load_data = function(){
    set_loaded(TRUE)
  },
  ensure = function(){
    callSuper()
    set_ensure_has_run(TRUE)
  },
  set_loaded = function(new_value){
    m_has_been_loaded <<- new_value
  },
  set_ensure_has_run = function(new_value){
    m_ensure_has_run <<- new_value
  }
)

#####
# Test_Data_Object_Info
#####
Test_Data_Object_Info <- setRefClass("Test_Data_Object_Info", 
                                     contains = "RCTestCase",
                                     fields = list() )
Test_Data_Object_Info$methods(
  test_load_data = function(){
    DATA_FILE_NAME = "test.dat"
    DATA_PATH_NAME = "."
    EXPERIMENT_NAME = "Test"
    DESIGNATION = "analysis1"
    DATA_TEXT = "test2"
    
    info <- Data_Object_Info$new()

    info$data_file_name  = DATA_FILE_NAME
    info$data_path_name  = DATA_PATH_NAME
    info$experiment_name = EXPERIMENT_NAME
    info$designation     = DESIGNATION

    checkEquals( msg     = "info$data_file_name should be stored",
                 current =  info$data_file_name,
                 target  =  DATA_FILE_NAME)
    checkEquals( msg     = "info$data_path_name should be stored",
                 current =  info$data_path_name,
                 target  =  DATA_PATH_NAME)
    checkEquals( msg     = "info$experiment_name should be stored",
                 current =  info$experiment_name,
                 target  =  EXPERIMENT_NAME)
    checkEquals( msg     = "info$designation is should be stored",
                 current =  info$designation,
                 target  =  DESIGNATION)


  },
  test_load_only_validates_info = function(){
    silent=TRUE # refers to display the (expected) error messages - these need to be reviewed with silent set to FALSE
    
    GOOD_DATA_FILE_NAME = "test.dat"
    BAD_DATA_FILE_NAME = "junkjunk"
    DATA_PATH_NAME = "."
    EXPERIMENT_NAME = "Test"
    DESIGNATION = "analysis1"
    DATA_TEXT = "test2"
    
    info <- Data_Object_Info$new()
    
    # BUGBUG: Valdidating the data sounds nice but doesn't work when there are multiple files of interest - find a better way to manage this
    
    # 
    # msg = "Data_Object_Info$load_data() without information should lead to an error with a (useful) message about lack of info"
    # checkException(msg=msg,
    #                info$load_data(), 
    #                silent = silent)
    
    # BUGBUG: Again, nice idea but doesn't work well with other file names
    
    # msg = "Data_Object_Info$load_data() with a bad file name (one that can't be opened) should lead to a useful message about a missing file"
    # info <- Data_Object_Info$new()
    # info$data_file_name  = BAD_DATA_FILE_NAME
    # info$data_path_name  = DATA_PATH_NAME
    # info$experiment_name = EXPERIMENT_NAME
    # info$designation     = DESIGNATION
    # checkException(msg=msg,
    #                info$load_data(), 
    #                silent = silent)
    
    msg = "Data_Object_Info$load_data() with a good file name should simply continue (no data is loaded by the info object itself)"
    make_test_dat = function(){
      df <- data.frame(list(junk=DATA_TEXT), stringsAsFactors = FALSE)
      write.table(df, file=GOOD_DATA_FILE_NAME, sep = "\t", row.names = FALSE, col.names=TRUE)
    }
    make_test_dat()
    info$data_file_name  = GOOD_DATA_FILE_NAME
    info$load_data()
    
    target  = sprintf("%s_%s", info$experiment_name, info$designation)
    current = info$collection_name()
    msg = sprintf("Data_Object_Info$collection_name() should return a combination of experiment and designation; expected '%s' but got '%s'", target, current)
    checkEquals(msg=msg, target=target, current=current)
    
    # data_collection$load_data()
    # df_data <- data_collection$df_data
    # checkEquals("Data collection can load 'test' data.frame - it has a field named 'junk'",
    #             target  = "test2",
    #             current = df_data$junk)
    
  }
)
#####
# Test_Data_Processor
#####
Test_Data_Processor <- setRefClass("Test_Data_Processor", 
                                    contains = "RCTestCase",
                                    fields = list() )
Test_Data_Processor$methods(
  test_load_data = function(){
    processor <- Data_Processor$new(p_info = Data_Object_Info_737_two_step$new())
    
    processor$raw_data$load_data()
    df <- processor$raw_data$df
    
    # Validate raw data ...
    
    # Number of rows
    target  = 122589
    current = nrow(df)
    msg = sprintf("Raw data file should be loaded; expected %d rows but data.frame has %d", target, current)
    checkEquals(msg=msg, target=target, current=current)
  },
  test_load_i_fdr = function(){
    processor <- Data_Processor$new(p_info = Data_Object_Info_737_two_step$new())
    
    # This should cause a chain reaction of loads
    processor$i_fdr$load_data()
    
    raw_data    <- processor$raw_data   $df
    data_groups <- processor$data_groups$df
    densities   <- processor$densities  $df
    alpha       <- processor$alpha      $df
    i_fdr       <- processor$i_fdr      $df
    
    # Validate data ...
    
    # Number of rows
    check_n_rows = function(name_of_df=NULL, target=NULL){
      current = nrow(processor[[name_of_df]]$df)
      msg = sprintf("%s should be generated; expected %d rows but data.frame has %d", name_of_df, target, current)
      checkEquals(msg=msg, target=target, current=current)
    }
    
    TOTAL_ROWS_DATA      <- 122589
    TOTAL_ROWS_DENSITIES <- 512
    TOTAL_ROWS_ALPHA     <- 10
    check_n_rows("raw_data",    TOTAL_ROWS_DATA)      # Original PSM data
    check_n_rows("data_groups", TOTAL_ROWS_DATA)      # Each row describes all the groups/data for each PSM
    check_n_rows("densities",   TOTAL_ROWS_DENSITIES) # As deep as a "density" function
    check_n_rows("alpha",       TOTAL_ROWS_ALPHA)     # Number of groups analyzed
    check_n_rows("i_fdr",       TOTAL_ROWS_DATA)      # A different fdr for each individual PSM
    
  },
  test_load_1_percent_when_file_exists = function(){
    DATA_TEXT           = "stuff in file"
    GOOD_DATA_FILE_NAME = "junk.txt"
    make_test_dat = function(){
      df <- data.frame(list(junk=DATA_TEXT), stringsAsFactors = FALSE)
      write.table(df, file=GOOD_DATA_FILE_NAME, sep = "\t", row.names = FALSE, col.names=TRUE)
    }
    
    make_test_dat()
    info <- Data_Object_Info$new()
    info$data_path_name               <- "."
    info$data_file_name_1_percent_FDR <- GOOD_DATA_FILE_NAME
    data_processor <- Data_Processor$new(p_info = info)
    
    data_processor$raw_1_percent$load_data()
    checkEquals(msg = "Data_Object_Raw_1_Percent$exists() should be TRUE when data can be loaded",
                target  = TRUE,
                current = data_processor$raw_1_percent$exists())
    checkEquals(msg = "Data_Object_Raw_1_Percent$load_data() should have loaded something",
                target  = DATA_TEXT,
                current = data_processor$raw_1_percent$df[1,1])
  },
  test_load_1_percent_when_file_name_is_empty = function(){
    info <- Data_Object_Info$new()
    info$data_path_name               <- "."
    data_processor <- Data_Processor$new(p_info = info)
    
    # If we don't bother to set the 1 percent file name, we don't have one.
    checkEquals(msg = "If no 1 percent filename is set, info$file_path_1_percent_FDR() should be empty",
                target  = "",
                current = info$file_path_1_percent_FDR())
    checkEquals(msg = "If no 1 percent filename is set, Data_Object_Raw_1_Percent$exists() should be False",
                target  = FALSE,
                current = data_processor$raw_1_percent$exists())
    
    # Also allow for setting 1 percent file name to "" rather than leaving it blank
    checkEquals(msg = "If 1 percent filename is set to empty string, info$file_path_1_percent_FDR() should be empty",
                target  = "",
                current = info$file_path_1_percent_FDR())
    checkEquals(msg = "If 1 percent filename is set to empty string, Data_Object_Raw_1_Percent$exists() should be False",
                target  = FALSE,
                current = data_processor$raw_1_percent$exists())
  }
  #,
  # test_load_1_percent_when_file_does_not_exist = function(){
  #   
  # }
)
# Create Data_Groups

#####
# Test_Data_Object
#####
Test_Data_Object <- setRefClass("Test_Data_Object", 
                                contains = "RCTestCase",
                                fields = list(data_object="Data_Object") )
Test_Data_Object$methods(
  test_Mock_Data_Object_sets_has_loaded_on_load = function(){
    data_object <<- Mock_Data_Object$new()
    
    checkEquals(msg = "Default value of has_been_loaded should be FALSE",
                target =  FALSE,
                data_object$m_has_been_loaded)
    data_object$load_data()
    checkEquals(msg = "has_been_loaded should be TRUE after loading",
                target =  TRUE,
                data_object$m_has_been_loaded)
  },
  test_data_load_fails_in_base_class = function(){
    data_object <<- Data_Object$new()
    
    msg = "Running 'load_data' in base class should produce (useful) error message - this is an abstract method"
    checkException(msg=msg,
                   expr = data_object$load_data(),
                   silent = TRUE)
  },
  test_ensure_runs_load_iff_dirty = function(){
    data_object <<- Mock_Data_Object$new()
    
    msg = "ensure() should load data when data is dirty"
    data_object$set_loaded(FALSE)
    data_object$set_dirty(TRUE)
    data_object$ensure()
    checkEquals(msg=msg,
                target =  TRUE,
                data_object$m_has_been_loaded)
    
    msg = "ensure() should NOT load data when data is NOT dirty"
    data_object$set_loaded(FALSE)
    data_object$set_dirty(FALSE)
    data_object$ensure()
    checkEquals(msg=msg,
                target = FALSE,
                data_object$m_has_been_loaded)

  },
  test_is_dirty = function(){
    data_object <<- Mock_Data_Object$new() # BUGBUG: Error in MyRUnit - test independence is not maintained!!! (setup() isn't run every time?)
    
    msg = "is_dirty defaults to TRUE"
    checkEquals(msg=msg,
                target  = TRUE,
                current = data_object$m_is_dirty)
  },
  test_load_causes_ensure_in_parent = function(){
    child  <- Mock_Data_Object$new()
    parent <- Mock_Data_Object$new()
    
    child$append_parent(parent)
    parent$append_child(child)
    
    msg = "Mock_Data_Object$m_ensure_has_run defaults to FALSE"
    checkEquals(msg=msg,
                target  = FALSE,
                parent$m_ensure_has_run)
    
    msg = "Running load_data() in child should cause ensure() to run on parent"
    child$load_data()
    checkEquals(msg=msg,
                target  = TRUE,
                parent$m_ensure_has_run)
  },
  test_changing_dirty_makes_children_dirty = function(){
    child  <- Mock_Data_Object$new()
    parent <- Mock_Data_Object$new()
    
    child$append_parent(parent)
    parent$append_child(child)
    
    msg = "Changing state of is_dirty from TRUE to FALSE in the parent sets is_dirty to TRUE in child"
    child$m_is_dirty <- FALSE
    parent$set_dirty(FALSE)
    checkEquals(msg=msg,
                target  = TRUE,
                child$m_is_dirty)
    
    msg = "Changing state of is_dirty from FALSE to TRUE in the parent sets is_dirty to TRUE in child"
    child$m_is_dirty <- FALSE
    parent$set_dirty(TRUE)
    checkEquals(msg=msg,
                target  = TRUE,
                child$m_is_dirty)
    
    msg = "Calling set_dirty without changing state in parent leaves is_dirty alone in child"
    child$m_is_dirty <- FALSE
    parent$set_dirty(TRUE)
    checkEquals(msg=msg,
                target  = FALSE,
                child$m_is_dirty)
  }
    
)
#####
# Test_Analyze_Collection
#####
Test_Analyze_Collection <- setRefClass("Test_Analyze_Collection", 
                                    contains = "RCTestCase",
                                    fields = list() )
Test_Analyze_Collection$methods(
  test_no_op = function(){
    
  }
)
# Load data_groups
