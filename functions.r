###############################################################################
# analysis.R                                                                  #
#                                                                             #
# Project 020 - Create GitHub module for paper                                #
#                                                                             #
# Support functions for main.r                                                #
#                                                                             #
# Started: 2019-02-28                                                         #
#                                                                             #
# Version: (1.1)                                                              #
#                                                                             #
###############################################################################
library(stringr)

###############################################################################
# Global parameters
###############################################################################
MIN_GOOD_PEPTIDE_LENGTH <- 11

###############################################################################
# Generic Functions
###############################################################################
rename_column <- function(df=NULL, name_before=NULL, name_after=NULL, suppressWarnings=FALSE){
  if (is.null(df)){
    stop("Dataframe (df) does not exist - unable to rename column")
  }
  if (name_before %in% colnames(df)){
    df[,name_after]  <- df[,name_before]
    df[,name_before] <- NULL
  } else if (!suppressWarnings){
    warning(sprintf("'%s' is not a field in the data frame and so has not been renamed", name_before))
  }
  return(df)
}
rename_columns <- function(df=NULL, names_before=NULL, names_after=NULL){
  for (i in safe_iterator(length(names_before))){
    df <- rename_column(df, names_before[i], names_after[i])
  }
  return(df)
}
write_standard_df <- function(file_path=NULL, data=NULL){
  # Support function, for "clarity"
  create_blank_df_with_colnames <- function(colnames){
    df <- read.table(text="", col.names = colnames)
    return(df)
  }
  
  # Main code
  cols <- create_blank_df_with_colnames(colnames(data))
  cols[1,]=colnames(data)
print(file_path)
  write.table(x = cols, file = file_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE) # Was attempting to save real column names.  Didn't work.
  write.table(x = data, file = file_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
}
function_to_tolerance <- function(x=NULL, tolerance=NULL, FUN=NULL, ...) { return(FUN(x/tolerance, ...) * tolerance) }
round_to_tolerance    <- function(x=NULL, tolerance=NULL, ...){ return(function_to_tolerance(x=x, tolerance=tolerance, FUN=round, ...)) }
floor_to_tolerance    <- function(x=NULL, tolerance=NULL, ...){ return(function_to_tolerance(x=x, tolerance=tolerance, FUN=floor, ...)) }
safe_median <- function(x) median(x, na.rm=TRUE)
safe_mean   <- function(x) mean(  x, na.rm=TRUE)
normalize_density <- function(d){
  # Normalizes y-values in density function
  # so that the integral under the curve is 1
  # (uses rectangles to approximate area)
  delta_x               <- diff(range(d$x)) / length(d$x)
  unnormalized_integral <- delta_x * sum(d$y)
  new_d   <- d
  new_d$y <- with(new_d, y )
  
  return(new_d)
}
verify_data_frame <- function(df = NULL){
  if (is.null(df)) {
    stop("data.frame variable is not set")
  } else if (! is.data.frame(df)) {
    stop("Object is not a data.frame but one was expected")
  } else if (nrow(df) == 0){
    stop("data.frame is empty")
  }
}

if_null <- function(cond=NULL, null_result=NULL, not_null_result=NULL){
  return(switch(1+is.null(cond), 
                not_null_result, 
                null_result))
}
rainbow_with_fixed_intensity <- function(n=NULL, goal_intensity_0_1=NULL, alpha=NULL){
  goal_intensity <- 255*goal_intensity_0_1
  hex_colors <- rainbow(n)
  rgb_colors <- col2rgb(hex_colors)
  df_colors <- data.frame(t(rgb_colors))
  df_colors$intensity <- with(df_colors, 0.2989*red + 0.5870*green + 0.1140*blue)
  
  df_colors$white_black <- with(df_colors, ifelse(intensity < goal_intensity, 255, 0))
  df_colors$mix_level   <- with(df_colors, (white_black - goal_intensity) / (white_black - intensity  ) )
  df_colors$new_red     <- with(df_colors, mix_level*red   + (1-mix_level)*white_black)
  df_colors$new_green   <- with(df_colors, mix_level*green + (1-mix_level)*white_black)
  df_colors$new_blue    <- with(df_colors, mix_level*blue  + (1-mix_level)*white_black)
  names_pref_new <- c("new_red", "new_green", "new_blue")
  names_no_pref  <- c("red", "green", "blue")
  df_colors <- df_colors[,names_pref_new]
  df_colors <- rename_columns(df_colors, names_before = names_pref_new, names_after = names_no_pref)
  rgb_colors <-as.matrix(df_colors/255 )
  
  return(rgb(rgb_colors, alpha=alpha))
}
safe_iterator <- function(n_steps = NULL){
  if (n_steps < 1){
    result = numeric(0)
  } else {
    result = 1:n_steps
  }
  return(result)
}
col2hex <- function(cols=NULL, col_alpha=255){
  if (all(col_alpha<=1)){
    col_alpha <- round(col_alpha*255)
  }
  col_matrix <- t(col2rgb(cols))
  results <- rgb(col_matrix, alpha=col_alpha, maxColorValue = 255)
  return(results)
}
credible_interval <- function(x=NULL, N=NULL, precision=0.001, alpha=0.05){
  # Approximates "highest posterior density interval"
  # Uses exact binomial but with a finite list of potential values (1/precision)
  
  p <- seq(from=0, to=1, by=precision)
  d <- dbinom(x = x, size = N, prob = p)
  d <- d / sum(d)
  df <- data.frame(p=p, d=d)
  df <- df[order(-df$d),]
  df$cumsum <- cumsum(df$d)
  max_idx <- sum(df$cumsum < (1-alpha)) + 1
  max_idx <- min(max_idx, nrow(df))
  
  lower <- min(df$p[1:max_idx])
  upper <- max(df$p[1:max_idx])
  
  return(c(lower,upper))
}
safe_is_true <- function(x=NULL){
  if (is.logical(x)){
    if (length(x) == 1){
      return(x)
    }
  } 
  
  return(FALSE)
}

###############################################################################
# Load data Descriptions
###############################################################################
# load_data_collection <- function(data_collection = NULL, process_parameters = NULL ){
#   
#   data_collection$data_original <- read.table(data_collection$file_name_dataset)
#   file_path_data
#   
#   return(data_collection)
# }

###############################################################################
# Data Handling
###############################################################################
load_standard_df <- function(file_path=NULL){
  if (file.exists(file_path)){
    data <- read.table(file = file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, blank.lines.skip = TRUE)#, check.names = FALSE)
  } else {
    stop(sprintf("File path does not exist: '%s'", file_path))
  }
  return(data)
}
# load_data <- function(data_collection=NULL){
#   data_collection$data_original <- load_standard_df(file_path = data_collection$file_path_data)
#   return(data_collection)
# }
standardize_fields <- function(data=NULL){
  data_new <- data
  manage_decoy_column <- function(data=NULL){
    is_column_of <- function(colname=NULL, data=NULL){
      return(colname %in% colnames(data))
    }
    data_new <- data
    if (is_column_of("Decoy", data_new)){
      data_new <- rename_column(data_new, "Decoy" , "decoy")
    } else {
      data_new$decoy <- 0
    }
    
    return(data_new)
  }
  
  #data_new <- rename_column(data_new, "Variable.Modifications"   , "ptm_list")
  data_new <- rename_column(data_new, "Confidence...."           , "confidence")
  data_new <- rename_column(data_new, "Precursor.m.z.Error..ppm.", "precursor_ppm")
  #data_new <- rename_column(data_new, "Isotope.Number"           , "isotope_number")
  #data_new <- rename_column(data_new, "m.z"                      , "m_z")
  #data_new <- rename_column(data_new, "Measured.Charge"          , "charge")
  data_new <- rename_column(data_new, "Spectrum.File"            , "spectrum_file")
  data_new <- rename_column(data_new, "Protein.s."               , "proteins")
  data_new <- rename_column(data_new, "Spectrum.Title"           , "spectrum_title")
  data_new <- manage_decoy_column(data_new)
  
  data_new$value          <- data_new$precursor_ppm
  data_new$peptide_length <- str_length(data_new$Sequence)
  #data_new$charge_value   <- with(data_new, as.numeric(substr(charge, start=1, stop=str_length(charge)-1)))
  #data_new$measured_mass  <- with(data_new, m_z*charge_value)
  data_new$spectrum_index <- NA
  data_new$spectrum_index[order(data_new$spectrum_title, na.last = TRUE)] <- 1:nrow(data_new)
  
  return(data_new)
}

###############################################################################
# Processing Data Collections
###############################################################################
process_data_collection <- function(process_parameters = NULL, data_collection = NULL, verify=TRUE, function_grouping = group_by_confidence2,  function_make_ci=NULL, nReps=NULL, ...){
  
  #### Support code for process_data_collection()  
  
  report_invalid_approach <- function(...){
    stop("Cannot 'make' original data")
  }
  load_data_file <- function(file_type = NULL, file_path = NULL){
    
    if        (file_type == "data.frame"){
      new_data <- load_standard_df(file_path = file_path)
    } else if (file_type == "RDS"){
      new_data <- readRDS(file = file_path)
    } else {
      print(sprintf("Failed to READ file -- invalid file type: '%s'", file_type))
    }
    new_data <- load_standard_df(file_path = file_path)
    
    return(new_data)
  }
  load_data_collection <- function(data_collection = NULL, process_parameters = NULL ){
    process_command <- process_parameters$load_collection
    if (process_command == "load"){
      file_path <- data_collection$file_path_collection
      data_collection <- readRDS(file = file_path)
    } else if (process_command == "skip") {
      # Do nothing
    }
    
    return(data_collection)
  }
  load_new_data_into_collection <- function(data_collection   = NULL, varname_file_path  = NULL, varname_new_data = NULL, file_type         = NULL){
    # Support functions (for load_new_data_into_collection)
    
    # Main code (for load_new_data_into_collection)
    file_path <- data_collection[[varname_file_path]]
    new_data  <- load_data_file(file_type = file_type, file_path = file_path)
    data_collection[[varname_new_data]] <- new_data
    
    return(data_collection)
  }
  make_new_data_for_collection <- function( data_collection   = NULL, varname_input_data = NULL, varname_new_data = NULL, varname_file_path = NULL, function_make = NULL, file_type = NULL, ...){
    
    # Support functions (for make_new_data_for_collection)
    save_data <- function(file_type = NULL, file_path = NULL, new_data  = NULL) {
      if        (file_type == "data.frame"){
        write_standard_df(file_path = file_path, data = new_data)
      } else if (file_type == "RDS"){
        saveRDS(new_data, file = file_path)
      } else if (file_type == "skip"){
        # Do nothing
      } else {
        print(sprintf("Failed to SAVE file -- invalid file type: '%s'", file_type))
      }
    }
    
    # Main code (for make_new_data_for_collection)
    if (varname_input_data == "data_collection"){
      old_data <- data_collection
    } else {
      old_data <- data_collection[[varname_input_data]]
    }
    new_data <- function_make(old_data, ...)
    data_collection[[varname_new_data]] <- new_data
    save_data(file_type = file_type,
              file_path = data_collection[[varname_file_path]],
              new_data  = new_data)
    return(data_collection)
  } 
  make_densities <- function(data_groups=NULL){
    
    # Support functions for make_densities()
    set_values_of_interest <- function(data_groups=NULL){
      field_decoy_group = "group_decoy_confidence"
      field_value       = "value_norm"
      new_data_groups <- data_groups
      new_data_groups$value_of_interest <- new_data_groups[,field_value]
      new_data_groups$group_of_interest <- new_data_groups[,field_decoy_group]
      #groups <- sort(unique(new_data_groups$group_of_interest))
      
      return(new_data_groups)
    }
    get_ylim <- function(data_groups=NULL){
      ylim <- range(data_groups$value_of_interest)
      
      return(ylim)
    }
    make_hit_density <- function(data_subset=NULL, ylim=NULL){
      uniformalize_density <- function(d){
        # Reorganizes y-values of density function so that 
        # function is monotone increasing to mode
        # and monotone decreasing afterwards
        idx_mode   <- which.max(d$y)
        idx_lower <- 1:(idx_mode-1)
        idx_upper <- idx_mode:length(d$y)
        
        values_lower <- d$y[idx_lower]
        values_upper <- d$y[idx_upper]
        
        new_d   <- d
        new_d$y <- c(sort(values_lower, decreasing = FALSE), 
                     sort(values_upper, decreasing = TRUE))
        
        return(new_d)
      }
      
      MIN_PEPTIDE_LENGTH = 11
      d <- with(subset(data_subset,  
                       (peptide_length >= MIN_PEPTIDE_LENGTH) & 
                         (used_to_find_middle == FALSE)), 
                density(value_of_interest, 
                        from = ylim[1], 
                        to   = ylim[2]))
      d <- normalize_density(   d)
      d <- uniformalize_density(d)
      
      return(d)
    }
    make_true_hit_density <- function(data_groups=NULL){
      d_true  <- make_hit_density(data_subset = subset(data_groups, (group_training_class == "good_testing") ),
                                  ylim        = get_ylim(data_groups))
      return(d_true)
    }
    make_false_hit_density <- function(data_groups=NULL){
      d_false <- make_hit_density(data_subset = subset(data_groups, (group_training_class == "bad_long") ),
                                  ylim        = get_ylim(data_groups))
      
      return(d_false)
    }
    add_v_densities <- function(data_groups=NULL, densities=NULL){
      groups <- sort(unique(data_groups$group_of_interest))
      
      new_densities <- densities
      
      for (local_group in groups){
        d_v <- make_hit_density(data_subset = subset(data_groups, (group_of_interest == local_group)),
                                ylim        = get_ylim(data_groups))
        new_densities[local_group] <- d_v$y
      }
      
      return(new_densities)
    }
    
    # Main section of make_densities()
    
    data_groups <- set_values_of_interest(data_groups=data_groups)
    d_true  <- make_true_hit_density( data_groups)
    d_false <- make_false_hit_density(data_groups)
    
    densities <- data.frame(x=d_true$x, 
                            t=d_true$y, 
                            f=d_false$y)
    densities <- add_v_densities(data_groups=data_groups, densities=densities)
    
    return(densities)
  }
  make_alpha_from_densities <- function(densities=NULL){
    max_of_density = apply(X = densities, MARGIN = 2, FUN = max)
    df_alpha <- data.frame(stringsAsFactors = FALSE,
                           list(v = max_of_density,
                                group_of_interest = names(max_of_density)))
    df_alpha <- subset(df_alpha, group_of_interest != "x")
    t <- with(subset(df_alpha, group_of_interest=="t"), v)
    f <- with(subset(df_alpha, group_of_interest=="f"), v)
    df_alpha$alpha <- with(df_alpha, (t-v)/(t-f))
    
    alpha <- df_alpha[,c("group_of_interest", "alpha")]
    alpha <- subset(alpha, (group_of_interest != "t") & (group_of_interest != "f"))
    
    return(alpha)
  }
  make_i_fdr <- function(data_collection = NULL){
    add_FDR_to_data_groups <- function(data_groups=NULL, densities=NULL, alpha=NULL, field_value=NULL, field_decoy_group=NULL, set_decoy_to_1=FALSE){
      
      # Support functions for add_FDR_to_data_groups()
      get_group_fdr <- function(group_stats = NULL, data_groups = NULL, densities=NULL){
        group_fdr <- apply(X = densities, MARGIN = 2, FUN = max)
        df_group_fdr <- data.frame(group_fdr)
        df_group_fdr <- rename_column(df_group_fdr, "group_fdr", "v")
        df_group_fdr$group_of_interest <- names(group_fdr)
        t <- df_group_fdr[df_group_fdr$group_of_interest == "t", "v"]
        f <- df_group_fdr[df_group_fdr$group_of_interest == "f", "v"]
        df_group_fdr <- subset(df_group_fdr, !(group_of_interest %in% c("x", "t", "f")))
        df_group_fdr$group_fdr <-(df_group_fdr$v - t) / (f - t)
        
        return(df_group_fdr)
      }
      
      get_mode <- function(x){
        d <- density(x)
        return(d$x[which.max(d$y)])
      }
      
      # Main code for add_FDR_to_data_groups()
      
      # Set up analysis
      data_new <- data_groups
      data_new$value_of_interest <- data_new[,field_value]
      data_new$group_of_interest <- data_new[,field_decoy_group]
      
      data_subset <- subset(data_new, peptide_length >= 11)
      
      # Identify mean confidence per group
      
      group_confidence <- aggregate(confidence~group_of_interest, data=data_subset, FUN=mean)
      group_confidence <- rename_column(group_confidence, "confidence", "group_confidence")
      
      #group_fdr   <- get_group_fdr(data_groups = data_subset, densities=densities)
      group_stats <- merge(alpha, group_confidence)
      group_stats <- subset(group_stats, group_of_interest != "decoy")
      
      x=c(0,group_stats$group_confidence)
      y=c(1,group_stats$alpha)
      FUN_interp <- approxfun(x=x,y=y)
      
      data_new$interpolated_groupwise_FDR <- FUN_interp(data_new$confidence)
      if (set_decoy_to_1){
        data_new$interpolated_groupwise_FDR[data_new$decoy == 1] <- 1
      }
      
      return(data_new)
    }
    
    data_groups = data_collection$data_groups
    densities   = data_collection$densities
    alpha       = data_collection$alpha
    
    d_true  <- densities[,c("x", "t")]
    d_false <- densities[,c("x", "f")]

    i_fdr <- add_FDR_to_data_groups(data_groups       = data_groups, 
                                    densities         = densities,
                                    alpha             = alpha,
                                    field_value       ="value_norm", 
                                    field_decoy_group = "group_decoy_confidence")
    # Derive local t
    interp_t <- splinefun(x=d_true$x,  y=d_true$t) #approxfun(x=d_true$x, y=d_true$y)
    
    # Derive local f
    interp_f <- splinefun(x=d_false$x, y=d_false$f) #approxfun(x=d_true$x, y=d_true$y)
    
    # Derive local FDR
    i_fdr$t     <- interp_t(i_fdr$value_of_interest)
    i_fdr$f     <- interp_f(i_fdr$value_of_interest)
    i_fdr$alpha <- i_fdr$interpolated_groupwise_FDR
    i_fdr$i_fdr <- with(i_fdr, (alpha*f) / (alpha*f + (1-alpha)*t)) 
    
    return(i_fdr)
  }
  make_data_groups <- function(data_original=NULL){
    # Functions supporting make_data_groups()
    
    add_grouped_variable <- function(data_groups = data_groups, field_name_to_group = NULL, vec.length.out = NULL, vec.tolerance = NULL, value_format = NULL){
      
      # Support functions for add_grouped_variable()
      find_interval_vec <- function(x=NULL, length.out = NULL, tolerance = NULL){
        q <- quantile(x = x, probs = seq(from=0, to=1, length.out = length.out), na.rm=TRUE)
        q <- round_to_tolerance(q, tolerance = tolerance)
        return(q)
      }
      get_group_data_frame <- function(vec=NULL, value_format = NULL){
        n <- length(vec)
        a <- vec[-n]
        b <- vec[-1]
        
        lower      <- ifelse(a == b           , "eq", NA)
        lower      <- ifelse(is.na(lower     ), "ge", lower)
        upper      <- ifelse(a == b           , "eq", NA)
        upper[n-1] <- ifelse(is.na(upper[n-1]), "le", "eq")
        upper      <- ifelse(is.na(upper     ), "lt", upper)
        group <- data.frame(list(idx=1:(n-1), a=a, b=b, lower=lower, upper=upper))
        
        name_format <- sprintf("%%%s_%%%s_%%s_%%s", value_format, value_format)
        group$new_var <- with(group, sprintf(name_format, a, b, lower, upper))
        
        return(group)
      }
      merge_group_with_data <- function(data_groups = NULL, group = NULL, vec = NULL, field_name_to_group = NULL){
        field_name_new <- sprintf("group_%s", field_name_to_group)
        group_idx      <- findInterval(x = data_groups[,field_name_to_group], 
                                       vec = vec, 
                                       all.inside=TRUE)
        
        data_groups$new_var <- group$new_var[group_idx]
        data_groups         <- rename_column(data_groups, "new_var", field_name_new)
      }
      # Body of add_grouped_variable()
      
      vec    <- find_interval_vec(x          = data_groups[[field_name_to_group]], 
                                  length.out = vec.length.out, 
                                  tolerance  = vec.tolerance )
      group  <- get_group_data_frame(vec          = vec, 
                                     value_format = value_format)
      df_new <- merge_group_with_data(data_groups         = data_groups, 
                                      group               = group, 
                                      vec                 = vec,
                                      field_name_to_group = field_name_to_group)
      df_new <- add_group_decoy(df_new, field_name_to_group)
      
      return(df_new)
    }
    add_already_grouped_variable <- function(field_name_to_group = NULL, data_groups = NULL ){
      old_name <- field_name_to_group
      new_name <- sprintf("group_%s", old_name)
      df_new <- data_groups
      df_new[[new_name]] <- data_groups[[old_name]]
      
      df_new <- add_group_decoy(data_groups = df_new, field_name_to_group = field_name_to_group)
      
      return(df_new)
    }
    add_value_norm <- function(data_groups = NULL){
      
      df_new            <- data_groups
      df_new$value_norm <- with(df_new, value - median_of_group_index)
      
      return(df_new)
    }
    add_protein_group <-function(data_groups = NULL){
      data_new <- data_groups
      df_group_def <- data.frame(stringsAsFactors = FALSE,
                                 list(pattern    = c(""     , "pfu_"      , "cRAP"),
                                      group_name = c("human", "pyrococcus", "trash")))
      for (i in 1:nrow(df_group_def)){
        idx <- grepl(pattern = df_group_def$pattern[i],
                     x       = data_new$proteins)
        data_new$group_proteins[idx] <- df_group_def$group_name[i]
      }
      
      data_new <- add_group_decoy(data_groups = data_new, field_name_to_group = "proteins")
      return(data_new)
    }
    add_group_decoy <- function(data_groups=NULL, field_name_to_group=NULL){
      field_name_decoy <- sprintf("group_decoy_%s", field_name_to_group)
      field_name_group <- sprintf("group_%s",       field_name_to_group)
      
      data_groups[[field_name_decoy]] <- with(data_groups, ifelse(decoy, "decoy", data_groups[[field_name_group]]))
      
      return(data_groups)
    }
    add_group_training_class <- function(data_groups = NULL){
      df_new <- data_groups
      
      lowest_confidence_group <- min(data_groups$group_confidence)
      
      is_long_enough   <- with(df_new, (peptide_length >= MIN_GOOD_PEPTIDE_LENGTH)    )
      is_good          <- with(df_new, (decoy == 0) & (confidence == 100)             )
      is_bad           <- with(df_new, (decoy == 1) )
      #is_used_to_train <- with(df_new, used_to_find_middle) # BUGBUG: circular definition
      
      idx_good         <- which(is_good         ) # & is_long_enough)
      n_good           <- length(idx_good)
      idx_testing      <- idx_good[c(TRUE,FALSE)] # Selects every other item
      idx_training     <- setdiff(idx_good, idx_testing)
      
      #is_good_short    <- with(df_new,  is_good      & !is_long_enough                )
      #is_good_long     <- with(df_new,  is_good      &  is_long_enough                )
      is_bad_short     <- with(df_new,  is_bad       & !is_long_enough                )
      is_bad_long      <- with(df_new,  is_bad       &  is_long_enough                )
      #is_good_training <- with(df_new,  is_good_long & (used_to_find_middle == TRUE ) )
      #is_good_testing  <- with(df_new,  is_good_long & (used_to_find_middle == FALSE) )
      
      df_new$group_training_class                   <- "other_short"   # Default
      df_new$group_training_class[is_long_enough  ] <- "other_long"    # Default (if long enough)
      df_new$group_training_class[idx_training    ] <- "good_training" # Length does not matter (anymore)
      df_new$group_training_class[idx_testing     ] <- "good_testing"  # Ditto
      #df_new$group_training_class[is_good_short   ] <- "good_short"
      df_new$group_training_class[is_bad_long     ] <- "bad_long"      # ...except for "bad"
      df_new$group_training_class[is_bad_short    ] <- "bad_short"
      
      df_new <- add_used_to_find_middle( data_groups = df_new ) # Guarantees consistency between duplicated definitions
      
      return(df_new)
    }
    add_used_to_find_middle <- function(data_groups = NULL){
      df_new    <- data_groups
      idx_used  <- which(data_groups$group_training_class == "good_training")
      
      df_new$used_to_find_middle           <- FALSE
      df_new$used_to_find_middle[idx_used] <- TRUE
      
      return(df_new)
    }
    add_group_spectrum_index <- function(data_groups = NULL){
      
      # Supporting functions for add_group_spectrum_index()
      
      get_breaks_all <- function(df_new){
        # Supporting function(s) for get_breaks_all()
        
        get_cut_points <- function(data_subset){
          
          # Supporting function(s) for get_cut_points()
          
          cut_values <- function(data=NULL, minimum_segment_length=NULL){
            # using cpt.mean -- Appears to have a memory leak
            #results_cpt <- cpt.mean(data=data, method="PELT", minimum_segment_length=minimum_segment_length)
            #results <- results_cpt@cpts
            
            # Just look at the end
            #results <- c(length(data))
            
            # regularly spaced, slightly larger than minimum_segment_length
            n_points <- length(data)
            n_regions <- floor(n_points / minimum_segment_length)
            n_regions <- ifelse(n_regions == 0, 1, n_regions)
            results <- round(seq(1, n_points, length.out = n_regions + 1))
            results <- results[-1]
            return(results)
          }
          remove_last <- function(x){
            return(x[-length(x)] )
          }
          
          # Main code of for get_cut_points()
          max_idx = max(data_subset$spectrum_index)
          data_sub_sub <- subset(data_subset, group_training_class == "good_training") #(confidence==100) & (decoy==0))
          minimum_segment_length = 50
          
          values <- data_sub_sub$value
          n_values <- length(values)
          local_to_global_idx <- data_sub_sub$spectrum_index
          if (n_values <= minimum_segment_length){
            result <- c()
          } else {
            local_idx <- cut_values(data=values, minimum_segment_length=minimum_segment_length)
            result <- local_to_global_idx[local_idx]
            result <- remove_last(result)
          }
          result <- c(result, max_idx)
          return(result)
        }
        remove_last <- function(vec) {
          return(vec[-length(vec)])
        }
        
        # Main code of get_breaks_all()
        
        breaks <- 1
        
        files <- unique(df_new$spectrum_file)
        
        for (local_file in files){
          data_subset <- subset(df_new, (spectrum_file==local_file))
          if (nrow(data_subset) > 0){
            breaks <- c(breaks, get_cut_points(data_subset))
          }
        }
        breaks <- sort(unique(breaks))
        breaks <- remove_last(breaks)
        breaks <- c(breaks, max(df_new$spectrum_index + 1))
        
        return(breaks)
      }
      
      # Main code of add_group_spectrum_index()
      
      field_name_to_group <- "spectrum_index"
      
      df_new <- data_groups[order(data_groups[[field_name_to_group]]),]
      breaks <- get_breaks_all(df_new)
      
      df_new$group_spectrum_index <- cut(x = df_new[[field_name_to_group]], breaks = breaks, right = FALSE, dig.lab = 6)
      df_new <- add_group_decoy(data_groups = df_new, field_name_to_group = field_name_to_group)
      
      return(df_new)
    }
    add_median_of_group_index <-function(data_groups = NULL){
      field_median <- "median_of_group_index"
      data_good <- subset(data_groups, used_to_find_middle )
      med <- aggregate(value~group_spectrum_index, data=data_good, FUN=safe_median)
      med <- rename_column(med, "value", field_median)
      
      data_groups[[field_median]] <- NULL
      df_new <- merge(data_groups, med)
      
      return(df_new)
    }
    
    # Main code of make_data_groups()
    verify_data_frame(data_original)
    print(dim(data_original))
    data_groups <- standardize_fields(data_original)
    
    data_groups <- add_grouped_variable(field_name_to_group = "confidence", 
                                        data_groups         = data_groups, 
                                        vec.length.out      = 14, 
                                        vec.tolerance       = 1, 
                                        value_format        = "03d")
    
    data_groups <- add_grouped_variable(field_name_to_group = "precursor_ppm", 
                                        data_groups         = data_groups, 
                                        vec.length.out      = 21, 
                                        vec.tolerance       = 0.1, 
                                        value_format        = "+05.1f")
    
    data_groups <- add_grouped_variable(field_name_to_group = "peptide_length", 
                                        data_groups         = data_groups, 
                                        vec.length.out      = 11, 
                                        vec.tolerance       = 1, 
                                        value_format        = "02d")
    
    # data_groups <- add_grouped_variable(field_name_to_group = "m_z", 
    #                                     data_groups         = data_groups, 
    #                                     vec.length.out      = 11, 
    #                                     vec.tolerance       = 10, 
    #                                     value_format        = "04.0f")
    # 
    # data_groups <- add_grouped_variable(field_name_to_group = "measured_mass", 
    #                                     data_groups         = data_groups, 
    #                                     vec.length.out      = 11, 
    #                                     vec.tolerance       = 1, 
    #                                     value_format        = "04.0f")
    # 
    # data_groups <- add_already_grouped_variable(field_name_to_group = "isotope_number",
    #                                             data_groups         = data_groups )
    # 
    # data_groups <- add_already_grouped_variable(field_name_to_group = "charge",
    #                                             data_groups         = data_groups )
    # 
    data_groups <- add_already_grouped_variable(field_name_to_group = "spectrum_file",
                                                data_groups         = data_groups )
    data_groups <- add_protein_group(data_groups = data_groups)
    data_groups <- add_group_training_class(  data_groups = data_groups)
    data_groups <- add_group_spectrum_index(  data_groups = data_groups)
    data_groups <- add_median_of_group_index( data_groups = data_groups)
    data_groups <- add_value_norm(            data_groups = data_groups)
    
    # fields_of_interest <- c("confidence", "precursor_ppm", "m_z", "peptide_length", "isotope_number", "charge", "spectrum_file", "measured_mass", "spectrum_index", "proteins")
    # fields_of_interest <- c("value", 
    #                         "decoy",
    #                         "spectrum_title",
    #                         "median_of_group_index",
    #                         "value_norm",
    #                         "used_to_find_middle",
    #                         "group_training_class",
    #                         fields_of_interest, 
    #                         sprintf("group_%s"      , fields_of_interest),
    #                         sprintf("group_decoy_%s", fields_of_interest))
    
    fields_of_interest <- c("confidence", "precursor_ppm", "peptide_length", "spectrum_file", "spectrum_index", "proteins")
    fields_of_interest <- c("value",
                            "decoy",
                            "spectrum_title",
                            "median_of_group_index",
                            "value_norm",
                            "used_to_find_middle",
                            "group_training_class",
                            fields_of_interest,
                            sprintf("group_%s"      , fields_of_interest),
                            sprintf("group_decoy_%s", fields_of_interest))
    
  data_groups <- data_groups[,fields_of_interest]
    
    return(data_groups)
  }
  add_1_percent_to_data_groups <- function(data_collection=NULL){
    
    # Support function(s) for get_1percent_data()
    verify_requirements_for_1_percent <- function(data_collection=NULL){
      verified_data_groups <- !is.null(data_collection$data_groups)
      verified_1percent    <- !is.null(data_collection$file_name_dataset_1_percent)
      return(verified_data_groups & verified_1percent)
    }
    
    # Main code for get_1percent_data()
    
    # Is there a file to load?
    data_new <- data_collection$data_groups
    
    if (verify_requirements_for_1_percent(data_collection)){
      # Load 1 percent file
      file_type = "data.frame"
      file_path = data_collection$file_path_dataset_1_percent
      df_1_percent <- load_data_file(file_type = file_type, file_path = file_path)
      
      # Get relevant fields
      df_1_percent$is_in_1percent <- TRUE
      df_1_percent                <- rename_column(df_1_percent, "Spectrum.Title", "spectrum_title")
      df_1_percent                <- df_1_percent[,c("spectrum_title", "is_in_1percent")]

      # Merge with data_groups
      data_new <- merge(data_new, df_1_percent, all.x=TRUE)
      data_new$is_in_1percent[is.na(data_new$is_in_1percent)] <- FALSE
    }
    
    # Save results
    return(data_new)
    
  }
  
  verify_original_data <- function(data_collection = NULL){
    
    d <- dim(data_collection$data_original)
    print(sprintf("size of original data: %d x %d", d[1], d[2]))
  }
  verify_densities <- function(data_collection=NULL){
    densities <- data_collection$densities
    groups    <- colnames(densities)
    groups    <- setdiff(groups, c("x", "t", "f"))
    main <- sprintf("Distribution of t, v, and f\n(Dataset: %s)",
                    data_collection$dataset_name)
    ylab="Density"
    xlab="PMD"
    plot(t~x, data=densities, type="l", lwd=5, main=main, xlab=xlab, ylab=ylab)
    for (local_group in groups){
      densities$local_group <- densities[[local_group]]
      lines(local_group~x, data=densities)
    }
    
  }
  verify_alpha_from_densities <- function(data_collection = NULL){
    alpha <- data_collection$alpha
    
    print(data_collection$alpha)
    #plot_PMA_distributions_with_peak_heights(data_collection=data_collection)
    #plot_FDR_vs_confidence(data_collection = data_collection)
  }
  verify_i_fdr <- function(data_collection = NULL){
    i_fdr <- data_collection$i_fdr
    
    print(head(i_fdr))
    main = sprintf("Comparison of Confidence and i-FDR\n(Dataset: %s)", data_collection$dataset_name)
    xlab="Individual PMD-FDR"
    ylab="Confidence Score (from Peptide Shaker)"
    ylim=c(0, 100)
    pch="."
    plot( confidence~i_fdr, data=i_fdr, pch=pch, xlab=xlab, ylab=ylab, ylim=ylim, main=main)
    points(confidence~alpha, data=i_fdr, col="purple")
    abline(a=100, b=-100)
    
    print(head(i_fdr[order(i_fdr$i_fdr),]))
  }
  verify_1_percent_TD_FDR <- function(data_collection=NULL){
    print("Table for data_collection$data_groups$is_in_1percent")
    print(table(data_collection$data_groups$is_in_1percent))
  }
  verify_data_groups <- function(data_collection=NULL){
    data_groups <- data_collection$data_groups
    
    fields_of_interest <- c("confidence", "precursor_ppm", "m_z", "peptide_length", "isotope_number", "charge", "spectrum_file", "measured_mass", "spectrum_index")
    for (field_name in fields_of_interest){
      group_name <- sprintf("group_%s"      , field_name)
      decoy_name <- sprintf("group_decoy_%s", field_name)
      
      print(table(data_groups[[group_name]]))
      print(table(data_groups[[decoy_name]]))
    }
    plot(value~spectrum_index, data=data_groups, pch=".")
    #plot_group_decoy_spectrum_index(data_collection)
    
    #plot_good_and_decoy_value(     data_collection)
    #plot_good_and_decoy_value_norm(data_collection)
    
    with(data_groups, print(table(decoy                     )))
  }
  save_data_collection <- function(data_collection = NULL, process_parameters = NULL ){
    if (process_parameters$save_collection == "save"){
      file_path <- data_collection$file_path_collection
      saveRDS(data_collection, file = file_path)
    }
  }
  
  process_data <- function(data_collection=NULL, process_parameters=NULL, varname_input_data=NULL, varname_process=NULL, varname_file_path=NULL, varname_new_data=NULL, function_make=NULL, function_verify=NULL, file_type=NULL, ...){
    
    # Main code (process_data)
    
    process_type <- process_parameters[[varname_process]]
    
    if        (is.null(process_type) ){
      stop (sprintf("Process '%s' was not defined in the process parameters", varname_process))
    }
    
    if        (process_type == "load"){
      data_collection <- load_new_data_into_collection(data_collection   = data_collection,
                                                       varname_file_path = varname_file_path,
                                                       varname_new_data  = varname_new_data,
                                                       file_type         = file_type)
    } else if (process_type == "make"){
      data_collection <- make_new_data_for_collection(data_collection     = data_collection,
                                                      varname_input_data  = varname_input_data,
                                                      varname_new_data    = varname_new_data,
                                                      varname_file_path   = varname_file_path,
                                                      function_make       = function_make,
                                                      file_type           = file_type,
                                                      ... )
    } else if (process_type == "skip"){
      # Do nothing
      data_collection <- data_collection
    } else {
      stop(sprintf("'%s' is not an accepted setting for process_parameters$%s", process_type, varname_process))
    }
    
    if (is.null(data_collection[[varname_new_data]])){
      print(sprintf("skipping verification of %15s - it has not been created", varname_new_data ))
    } else if (verify){
      function_verify(data_collection=data_collection)
    }
    
    return(data_collection)
  }
  
  # Main code for process_data_collection()
  
  err_msg <- ""
 
  time_to_process = system.time({
    # Collection
    # data_collection <- process_data(data_collection    = data_collection,
    #                                 process_parameters = process_parameters,
    #                                 varname_process    = "load_collection",
    #                                 varname_new_data   = "data_original",
    #                                 file_type          = "RDS",
    #                                 varname_file_path  = "file_path_collection",
    #                                 varname_input_data = "",
    #                                 function_make      = report_invalid_approach,
    #                                 function_verify    = verify_original_data)
    if (process_parameters$load_collection == "load"){
      data_collection <- readRDS(file=data_collection$file_path_collection)
    }
    
    if ( ! safe_is_true(data_collection$file_names_made)){
      if (safe_is_true(data_collection$dataset_experiment_name != "")){
        err_msg <- (sprintf("File names for data collection '%s' has not been made", data_collection$dataset_experiment_name))
      } else {
        err_msg <- (sprintf("Data collection not initialized - missing experiment name as well as file names"))
      }
      print(err_msg) # The timer eats error messages
      stop(err_msg)
    }
    
    # BUGBUG: Need to verify structure of "process_parameters" before using it
    # BUGBUG: Need more graceful failure when "make_file_names" has not been run
    
    # Original Data
    data_collection <- process_data(data_collection    = data_collection,
                                    process_parameters = process_parameters,
                                    varname_process    = "data",
                                    varname_new_data   = "data_original",
                                    file_type          = "data.frame",
                                    varname_file_path  = "file_path_data",
                                    varname_input_data = "",
                                    function_make      = report_invalid_approach,
                                    function_verify    = verify_original_data)
    
    # # Clean data
    # data_collection <- process_data(data_collection     = data_collection,
    #                                 process_parameters  = process_parameters,
    #                                 varname_process     = "clean_data",
    #                                 varname_new_data    = "clean_data",
    #                                 file_type           = "data.frame",
    #                                 varname_file_path   = "file_path_clean_data",
    #                                 varname_input_data  = "data_original",
    #                                 function_make       = make_clean_data,
    #                                 function_verify     = verify_clean_data,
    #                                 function_grouping   = function_grouping,
    #                                 ...)

    # Create data_groups
    data_collection <- process_data(data_collection     = data_collection,
                                    process_parameters  = process_parameters,
                                    varname_process     = "grouping",
                                    varname_new_data    = "data_groups",
                                    file_type           = "data.frame",
                                    varname_file_path   = "file_path_groups",
                                    varname_input_data  = "data_original",
                                    function_make       = make_data_groups,
                                    function_verify     = verify_data_groups )

    # Add 1% TD-FDR to data_groups
    data_collection <- process_data(data_collection     = data_collection,
                                    process_parameters  = process_parameters,
                                    varname_process     = "one_percent_FDR",
                                    varname_new_data    = "data_groups",
                                    file_type           = "data.frame",
                                    varname_file_path   = "file_path_groups",
                                    varname_input_data  = "data_collection",
                                    function_make       = add_1_percent_to_data_groups,
                                    function_verify     = verify_1_percent_TD_FDR )

    # Create densities
    data_collection <- process_data(data_collection     = data_collection,
                                    process_parameters  = process_parameters,
                                    varname_process     = "densities",
                                    varname_new_data    = "densities",
                                    file_type           = "data.frame",
                                    varname_file_path   = "file_path_densities",
                                    varname_input_data  = "data_groups",
                                    function_make       = make_densities,
                                    function_verify     = verify_densities )

    # Approximate alpha
    data_collection <- process_data(data_collection     = data_collection,
                                    process_parameters  = process_parameters,
                                    varname_process     = "fdr_approx",
                                    varname_new_data    = "alpha",
                                    file_type           = "skip",
                                    varname_file_path   = "file_path_alpha",
                                    varname_input_data  = "densities",
                                    function_make       = make_alpha_from_densities,
                                    function_verify     = verify_alpha_from_densities )
    # 
    # data_collection <- process_data(data_collection     = data_collection,
    #                                 process_parameters  = process_parameters,
    #                                 varname_process     = "fdr_approx",
    #                                 varname_new_data    = "alpha",
    #                                 file_type           = "skip",
    #                                 varname_file_path   = "file_path_alpha",
    #                                 varname_input_data  = "clean_data",
    #                                 function_make       = make_alpha_approx,
    #                                 function_verify     = verify_alpha )
    
    # Local FDR (i-FDR, "individual")
    data_collection <- process_data(data_collection     = data_collection,
                                    process_parameters  = process_parameters,
                                    varname_process     = "i_fdr",
                                    varname_new_data    = "i_fdr",
                                    file_type           = "data.frame",
                                    varname_file_path   = "file_path_i_fdr",
                                    varname_input_data  = "data_collection",
                                    function_make       = make_i_fdr,
                                    function_verify     = verify_i_fdr )
    # Calculate confidence intervals
    
    # data_collection <- process_data(data_collection    = data_collection,
    #                                 process_parameters = process_parameters,
    #                                 varname_process    = "ci_results",
    #                                 varname_new_data   = "ci_results",
    #                                 file_type          = "RDS",
    #                                 varname_file_path  = "file_path_CI",
    #                                 varname_input_data = "clean_data",
    #                                 function_make      = function_make_ci,
    #                                 function_verify    = verify_confidence_intervals,
    #                                 nReps              = nReps )
    
    save_data_collection(data_collection    = data_collection,
                         process_parameters = process_parameters )
    
    return(data_collection)
    
  })
  
  print(time_to_process)
  
  return(data_collection)
}
apply_process_to_one_collection <- function(data_collection = NULL, process_parameters=NULL, over_write=NULL, return_data_collection=TRUE){
  #data_collection <- make_file_names(data_collection = data_collection)
  
  file_path <- data_collection$file_path_collection
  create_file <- over_write | ! file.exists(file_path)
  print(sprintf("File %s: %s",
                ifelse(file.exists(file_path), "exists", "does not exist"),
                file_path))
  if (create_file){
    print(sprintf("    Loading and running %s", data_collection$dataset_designation))
    data_collection <- process_data_collection( data_collection    = data_collection, 
                                                process_parameters = process_parameters, 
                                                function_grouping  = group_by_confidence3,  
                                                function_make_ci   = get_ci_results, 
                                                nReps              = 10 )    }
  return(data_collection)
}

create_data_collection_from_scratch <- function(data_collection=NULL, save_results=TRUE){
  save_collection <- ifelse(save_results, "save", "skip")
  process_parameters = list(load_collection = "skip",
                            data            = "load",
                            clean_data      = "skip",
                            grouping        = "make",
                            densities       = "make",
                            fdr_approx      = "make",
                            i_fdr           = "make",
                            ci_results      = "skip",
                            one_percent_FDR = "make",
                            save_collection = save_collection)
  result <- apply_process_to_one_collection(data_collection = data_collection, process_parameters = process_parameters, over_write = TRUE)
  return(result)
}
#data_collection = data_collection_a

###############################################################################
# Plotting Classes
###############################################################################
###############################################################################
######## Class: Plot_Image
###############################################################################
Plot_Image = setRefClass("Plot_Image",
                         fields = list(data_collection    = "list",
                                       plot_title         = "character",
                                       include_text       = "logical",
                                       include_main       = "logical", 
                                       x.intersp          = "numeric",
                                       y.intersp          = "numeric",
                                       scale              = "numeric",
                                       main               = "character",
                                       is_image_container = "logical"))
Plot_Image$methods(
  initialize = function(p_data_collection = list(), 
                        p_include_main = TRUE, 
                        p_include_text = TRUE,
                        p_is_image_container = FALSE){
    include_main    <<- p_include_main
    include_text    <<- p_include_text
    data_collection <<- p_data_collection
    is_image_container <<- p_is_image_container
  },
  plot_image = function(){
    plot(main="Define plot_image() for subclass") # Abstract function
  },
  get_n = function(){
    stop("Need to define function get_n() for subclass") #Abstract function
  },
  create_standard_main = function(){
    needs_main <- function(){
      return(include_text & include_main & !is_image_container)
    }
    if (needs_main()){
      main <<- sprintf("%s\n(Dataset: %s; n=%s)", plot_title, data_collection$dataset_name, format(get_n(), big.mark = ","))
    }
  },
  plot_image_in_window = function(p_scale=NULL, window_height=NULL, window_width=NULL){
    scale <<- p_scale
    SIZE_AXIS      <- 2.5 * scale # in the units used by mar
    SIZE_MAIN      <- 2.5 * scale
    SIZE_NO_MARGIN <- 0.1 * scale
    FONT_SIZE      <- 8   * scale
    WINDOW_WIDTH   <- window_width  * scale
    WINDOW_HEIGHT  <- window_height * scale
    X_INTERSP      <- 0.5 * scale + 0.4 # manages legend text spacing
    Y_INTERSP      <- 0.5 * scale + 0.4 # manages
    
    if (include_main){
      mar = c(SIZE_AXIS, SIZE_AXIS, SIZE_MAIN     , SIZE_NO_MARGIN)
    } else {
      mar = c(SIZE_AXIS, SIZE_AXIS, SIZE_NO_MARGIN, SIZE_NO_MARGIN)
    }
    mgp = c(SIZE_AXIS/2, SIZE_AXIS/4, 0) # Margin line (mex units) for axis title, axis labels, axis lines
    ps  = FONT_SIZE
    x.intersp <<- X_INTERSP
    y.intersp <<- Y_INTERSP
    
    windows(width = WINDOW_WIDTH, height=WINDOW_HEIGHT)
    
    old_par  <- par(mar=mar, ps=ps, mgp=mgp)
    create_standard_main()
    
    plot_image()
    if (!is_image_container){
      axis(side=1, labels=include_text, tcl=-0.5, lwd=scale)
      axis(side=2, labels=include_text, tcl=-0.5, lwd=scale)
      box(lwd=scale)
    }
    par(old_par)
  },
  plot_image_in_small_window = function(p_scale=1){
    plot_image_in_window(p_scale=p_scale, window_height=2, window_width=3.25)
  },
  plot_image_in_large_window = function(p_scale=1, window_height=NULL){
    plot_image_in_window(p_scale=p_scale, window_height=window_height, window_width=7)
  }
)
###############################################################################
######## Class: Legend_Object
###############################################################################
Legend_Object = setRefClass("Legend_Object",
                            contains = "Plot_Image",
                            fields = list(user_params = "list",
                                          scale       = "numeric"))
Legend_Object$methods(
  initialize = function(p_user_params = NULL, p_scale = NULL){
    if (is.null(p_user_params)){
      user_params <<- list()
    } else {
      user_params <<- p_user_params
    }
    if (is.null(p_scale)){
      stop("Legend_Object must have a valid scale")
    } else {
      scale <<- p_scale
    }
    user_params$lwd       <<- if_null(user_params$lwd      , NULL, user_params$lwd        * scale) # Because we allow NULL, scale must be inside parens
    user_params$seg.len   <<- if_null(user_params$seg.len  ,    3, user_params$seg.len  ) * scale
    user_params$box.lwd   <<- if_null(user_params$box.lwd  ,    1, user_params$box.lwd  ) * scale
    user_params$x.intersp <<- if_null(user_params$x.intersp,  0.6, user_params$x.intersp) * scale
    user_params$y.intersp <<- if_null(user_params$y.intersp,  0.4, user_params$y.intersp) * scale + 0.2
  },
  show = function(){
    first_legend = legend(x         = "topleft", 
                          title     = "", 
                          legend    = user_params$leg, 
                          col       = user_params$col, 
                          lty       = user_params$lty, 
                          lwd       = user_params$lwd, 
                          seg.len   = user_params$seg.len, 
                          box.lwd   = user_params$box.lwd, 
                          x.intersp = user_params$x.intersp, 
                          y.intersp = user_params$y.intersp)
    new_x = first_legend$rect$left 
    new_y = first_legend$rect$top + first_legend$rect$h * ifelse(scale==1, 0.07, 0.03 - (scale * 0.02)) #switch(scale, 0.01, -0.01, -0.03, -0.05)# (0.07 - 0.09 * ((scale-1)^2))#(0.15 - 0.08*scale)#.07 * (2 - scale)
    legend(x=new_x, y=new_y, title = user_params$title, legend = "", cex=1.15, bty="n")
    
  }
)
###############################################################################
######## Class: Plot_Compare_PMD_and_Norm_Density
###############################################################################
Plot_Compare_PMD_and_Norm_Density = setRefClass("Plot_Compare_PMD_and_Norm_Density",
                                                contains = "Plot_Image",
                                                fields = list(show_norm      = "logical",
                                                              display_n_psms = "logical"))
Plot_Compare_PMD_and_Norm_Density$methods(
  initialize = function(p_show_norm=TRUE, p_display_n_psms=TRUE, ...){
    show_norm       <<- p_show_norm
    display_n_psms  <<- p_display_n_psms
    plot_title      <<- "True Hit and False Hit Distributions"
    
    callSuper(...)
  },
  plot_image = function(){
    
    # Support functions for plot_compare_PMD_and_norm_density()
    
    get_densities <- function(data_subset = NULL, var_value = NULL){
      data_subset$value_of_interest <- data_subset[,var_value]
      from <- min(data_subset$value_of_interest)
      to   <- max(data_subset$value_of_interest)
      xlim = range(data_subset$value_of_interest)
      data_true  <- subset(data_subset, (decoy==0) & (confidence==100))
      data_false <- subset(data_subset, (decoy==1))       
      d_true  <- with(data_true , density(value_of_interest, from = from, to = to))
      d_false <- with(data_false, density(value_of_interest, from = from, to = to))
      d_true  <- normalize_density(d_true)
      d_false <- normalize_density(d_false)
      
      densities <- list(d_true=d_true, d_false=d_false, var_value = var_value, n_true = nrow(data_true), n_false = nrow(data_false))
      
      return(densities)
    }
    get_xlim <- function(densities_a = NULL, densities_b = NULL, show_norm=NULL){
      xlim   <- range(c(      densities_a$d_true$x, densities_a$d_false$y))
      if (show_norm){
        xlim <- range(c(xlim, densities_b$d_true$x, densities_b$d_false$y))
      }
      return(xlim)
    }
    get_ylim <- function(densities_a = NULL, densities_b = NULL, show_norm=NULL){
      ylim   <- range(c(      densities_a$d_true$y, densities_a$d_false$y))
      if (show_norm){
        ylim <- range(c(ylim, densities_b$d_true$y, densities_b$d_false$y))
      }
      return(ylim)
    }
    plot_distributions <- function(densities = NULL, var_value= NULL, dataset_name = NULL, ...){
      leg = list()
      leg$leg = c("Good", "Bad")
      if (display_n_psms){
        leg$leg = sprintf("%s (%d PSMs)", 
                          leg$leg,
                          c(densities$n_true, densities$n_false))
        
      }
      leg$col = c("black", "red")
      leg$lwd = c(3      ,     3)
      leg$lty = c(1      ,     2)
      leg$title = "Hit Category"
      xlab = ifelse(var_value == "value",
                    "PMD (ppm)",
                    "PMD - normalized (ppm)")
      ylab = "Density"
      if (!include_text){
        xlab = ""
        ylab = ""
      }
      plot( densities$d_true , col=leg$col[1], lwd=leg$lwd[1] * scale, lty=leg$lty[1], xaxt = "n", yaxt = "n", main=main, xlab = xlab, ylab=ylab, ...)
      lines(densities$d_false, col=leg$col[2], lwd=leg$lwd[2] * scale, lty=leg$lty[2])
      abline(v=0, h=0, col="gray", lwd=1*scale)
      if (include_text){
        legend_object <- Legend_Object$new(leg, scale)
        legend_object$show()
        #legend("topleft", legend=leg.leg, col=leg.col, lwd=leg.lwd, lty=leg.lty, x.intersp = x.intersp, y.intersp = y.intersp)
      }
    }
    
    # Main code block for plot_compare_PMD_and_norm_density
    data_groups <- data_collection$data_groups
    
    data_subset_a <- subset(data_groups  , used_to_find_middle == FALSE)
    data_subset_b <- subset(data_subset_a, peptide_length > 11)
    
    densities_a <- get_densities(data_subset = data_subset_a, var_value = "value")
    densities_b <- get_densities(data_subset = data_subset_b, var_value = "value_norm")
    
    xlim=get_xlim(densities_a, densities_b, show_norm = show_norm)
    ylim=get_ylim(densities_a, densities_b, show_norm = show_norm)
    
    dataset_name <- data_collection$dataset_name
    plot_distributions(densities=densities_a, var_value = "value"     , dataset_name = dataset_name, xlim=xlim, ylim=ylim)
    if (show_norm){
      plot_distributions(densities_b, var_value = "value_norm", dataset_name = dataset_name, xlim=xlim, ylim=ylim)
    }
  },
  get_n = function(){
    data_subset_a <- subset(data_collection$data_groups  , used_to_find_middle == FALSE)
    data_subset_b <- subset(data_subset_a, peptide_length > 11)
    
    if (show_norm){
      data_subset <- data_subset_a
    } else {
      data_subset <- data_subset_b
    }
    
    data_true  <- subset(data_subset, (decoy==0) & (confidence==100))
    data_false <- subset(data_subset, (decoy==1))       
    
    return(nrow(data_true) + nrow(data_false))
  }
)
###############################################################################
######## Class: Plot_Multiple_Images
###############################################################################
Plot_Multiple_Images = setRefClass("Plot_Multiple_Images",
                                   contains = "Plot_Image",
                                   fields = list(n_images_wide = "numeric",
                                                 n_images_tall = "numeric",
                                                 image_list    = "list"))
Plot_Multiple_Images$methods(
  initialize = function(p_n_images_wide=1, p_n_images_tall=2, p_image_list=NULL, ...){
    n_images_wide  <<- p_n_images_wide
    n_images_tall  <<- p_n_images_tall
    image_list     <<- p_image_list
    #plot_title      <<- "True Hit and False Hit Distributions"
    
    callSuper(p_is_image_container=TRUE, ...)
  },
  plot_image = function(){
    # Support functions
    apply_mtext <- function(letter=NULL){
      line=1.3*scale
      mtext(letter, side=1, line=line, adj=0)
    }
    # main code
    old_par <- par(mfrow=c(n_images_tall, n_images_wide))
    i=0
    n_images <- length(image_list)
    
    for (i in 1:n_images){
      image <- image_list[[i]]
      image$create_standard_main()
      image$scale <- scale
      image$plot_image()
      axis(side=1, labels=include_text, tcl=-0.5, lwd=scale)
      axis(side=2, labels=include_text, tcl=-0.5, lwd=scale)
      box(lwd=scale)
      apply_mtext(letter=sprintf("(%s)", letters[i]))
      
    }
    par(old_par)
    
  }
)
###############################################################################
######## Class: Plot_Time_Invariance_Alt
###############################################################################
Plot_Time_Invariance_Alt = setRefClass("Plot_Time_Invariance_Alt",
                                       contains = "Plot_Image",
                                       fields = list(show_norm      = "logical",
                                                     display_n_psms = "logical",
                                                     training_class = "character",
                                                     ylim           = "numeric",
                                                     field_of_interest = "character"))
Plot_Time_Invariance_Alt$methods(
  initialize = function(p_ylim=NULL, p_training_class=NULL, p_field_of_interest="value_norm", ...){
    get_subset_title <- function(training_class=NULL){
      if        (training_class == "bad_long"){
        subset_title="bad only"
      } else if (training_class == "good_testing"){
        subset_title="good-testing only"
      } else if (training_class == "good_training"){
        subset_title="good-training only"
      } else if (training_class == "other"){
        subset_title="other only"
      } else {
        stop("Unexpected training_class in plot_time_invariance")
      }
      return(subset_title)
    }
    
    ylim <<- p_ylim
    training_class <<- p_training_class
    field_of_interest <<- p_field_of_interest
    subset_title <- get_subset_title(training_class=training_class)
    backup_title <- sprintf("Middle 25%% PMD for spectra sorted by index%s", 
                            ifelse(is.null(subset_title),
                                   "",
                                   sprintf(" - %s", subset_title)))
    #plot_title <<- get_main(main_title=main, backup_title=backup_title, data_collection = data_collection)
    plot_title <<- backup_title
    
    callSuper(...)
  },
  plot_image = function(){
    # Support functions for plot_time_invariance()
    
    # Main code of plot_time_invariance()
    data_subset = get_data_subset()
    plot_group_spectrum_index_from_subset_boxes(data_subset = data_subset)
    abline(h=0, col="blue", lwd=scale)
  },
  get_n = function(){
    return(nrow(get_data_subset()))
  },
  get_data_subset = function(){
    return(subset(data_collection$data_groups, (group_training_class==training_class)))
  },
  plot_group_spectrum_index_from_subset_boxes = function(data_subset = NULL){
    n_plot_groups <- 100
    
    field_name_text <- ifelse(field_of_interest=="value", "PMD", "Translated PMD")
    new_subset                   <- data_subset
    new_subset$value_of_interest <- new_subset[,field_of_interest]
    new_subset                   <- new_subset[order(new_subset$spectrum_index),]
    
    idxs <- round_to_tolerance(seq(from=1, to=nrow(new_subset), length.out = n_plot_groups+1), 1)
    idxs_left  <- idxs[-(n_plot_groups+1)]
    idxs_right <- idxs[-1] - 1
    idxs_right[n_plot_groups] <- idxs_right[n_plot_groups] + 1
    
    new_subset$plot_group <- NA
    for (i in 1:n_plot_groups){
      new_subset$plot_group[idxs_left[i]:idxs_right[i]] <- i 
    }
    xleft   <- aggregate(spectrum_index   ~plot_group, data=new_subset, FUN=min)
    xright  <- aggregate(spectrum_index   ~plot_group, data=new_subset, FUN=max)
    ybottom <- aggregate(value_of_interest~plot_group, data=new_subset, FUN=function(x){quantile(x, probs = 0.5 - (0.25/2))})
    ytop    <- aggregate(value_of_interest~plot_group, data=new_subset, FUN=function(x){quantile(x, probs = 0.5 + (0.25/2))})
    boxes <- merge(            rename_column(xleft  , "spectrum_index"   , "xleft"),
                               merge(      rename_column(xright , "spectrum_index"   , "xright"),
                                           merge(rename_column(ybottom, "value_of_interest", "ybottom"),
                                                 rename_column(ytop   , "value_of_interest", "ytop"))))
    
    xlab <- "Spectrum Index"
    ylab <- sprintf("%s (ppm)", field_name_text )
    if (is.null(ylim)){
      ylim <<- range(new_subset$value_of_interest)
    }
    if (!include_text){
      xlab=""
      ylab=""
    }
    plot(value_of_interest~spectrum_index, data=new_subset, type="n", ylim=ylim, xlab = xlab, ylab=ylab, main=main, xaxt="n", yaxt="n")
    with(boxes, rect(xleft = xleft, ybottom = ybottom, xright = xright, ytop = ytop, lwd=scale))
    #points(median_of_group_index~spectrum_index, data=data_subset, cex=.5, pch=15)
    axis(1, labels=include_text, lwd=scale)
    axis(2, labels=include_text, lwd=scale)
    box(lwd=scale) #box around plot area
  }
  
)
###############################################################################
######## Class: Plot_Time_Invariance_Alt_Before_and_After
###############################################################################
Plot_Time_Invariance_Alt_Before_and_After = setRefClass("Plot_Time_Invariance_Alt_Before_and_After",
                                                        contains = "Plot_Multiple_Images",
                                                        fields = list())
Plot_Time_Invariance_Alt_Before_and_After$methods(
  initialize = function(p_data_collection = data_collection, 
                        p_include_text=TRUE, 
                        p_include_main=FALSE,
                        p_ylim = c(-4,4), ...){
    plot_object1 <- Plot_Time_Invariance_Alt$new(p_data_collection = p_data_collection, 
                                                 p_include_text=p_include_text, 
                                                 p_include_main=p_include_main,
                                                 p_training_class = "good_testing",
                                                 p_field_of_interest = "value",
                                                 p_ylim = p_ylim)
    
    plot_object2 <- Plot_Time_Invariance_Alt$new(p_data_collection = p_data_collection, 
                                                 p_include_text=p_include_text, 
                                                 p_include_main=p_include_main,
                                                 p_training_class = "good_testing",
                                                 p_field_of_interest = "value_norm",
                                                 p_ylim = p_ylim)
    
    callSuper(p_n_images_wide=1, 
              p_n_images_tall=2, 
              p_include_text=p_include_text,
              p_include_main=p_include_main,
              p_image_list = list(plot_object1, plot_object2), ...)
  }
)

###############################################################################
######## Class: Plot_Density_PMD_and_Norm_Decoy_by_AA_Length
###############################################################################
Plot_Density_PMD_and_Norm_Decoy_by_AA_Length = setRefClass("Plot_Density_PMD_and_Norm_Decoy_by_AA_Length",
                                                           contains = "Plot_Image",
                                                           fields = list(show_norm = "logical"))
Plot_Density_PMD_and_Norm_Decoy_by_AA_Length$methods(
  initialize = function(p_show_norm=FALSE, ...){
    plot_title <<- "The Decoy Bump: PMD Distribution of Decoy matches by peptide length"
    show_norm  <<- p_show_norm
    callSuper(...)
  },
  get_n = function(){
    data_subset <- subset(data_collection$data_groups, (decoy == 1))
    return(nrow(data_subset))
  },
  plot_image = function(){
    
    # Support functions for plot_density_PMD_and_norm_decoy_by_aa_length()
    
    add_group_peptide_length_special <- function(){
      data_groups <- data_collection$data_groups
      data_groups <- subset(data_groups, used_to_find_middle == FALSE)
      
      df_group_definition <- data.frame(stringsAsFactors = FALSE,
                                        list(group_peptide_length_special = c("06-08", "09-10", "11-12", "13-15", "16-20", "21-50"),
                                             min                          = c(  6    ,   9    ,  11    ,  13    ,  16    ,  21    ),
                                             max                          = c(     8 ,     10 ,     12 ,     15 ,     20 ,     50 ) ))
      group_peptide_length_special     <- data.frame(list(peptide_length = 6:50))
      group_peptide_length_special$min <- with(group_peptide_length_special, sapply(peptide_length, FUN = function(i) max(df_group_definition$min[df_group_definition$min <= i])))
      group_peptide_length_special     <- merge(group_peptide_length_special, df_group_definition)
      
      data_groups$group_peptide_length_special <- NULL
      new_data_groups <- (merge(data_groups, 
                                group_peptide_length_special[,c("peptide_length", 
                                                                "group_peptide_length_special")]))
      return(new_data_groups)
    }
    get_densities <- function(data_subset = NULL, field_value = NULL, field_group=NULL){
      get_density_from_subset <- function(data_subset=NULL, xlim=NULL){
        
        d_group            <- with(data_subset , density(value_of_interest, from = xlim[1], to = xlim[2]))
        d_group            <- normalize_density(d_group)
        
        return(d_group)
      }
      
      data_temp                   <- data_subset
      data_temp$value_of_interest <- data_temp[[field_value]]
      data_temp$group_of_interest <- data_temp[[field_group]]
      
      xlim = range(data_temp$value_of_interest)
      
      groups      <- sort(unique(data_temp$group_of_interest))
      n_groups    <- length(groups)
      
      d_group <- get_density_from_subset( data_subset=data_temp, xlim = xlim )
      densities <- list("All decoys" = d_group)
      for (i in 1:n_groups){
        group <- groups[i]
        
        d_group <- get_density_from_subset( data_subset=subset(data_temp, (group_of_interest == group)), 
                                            xlim = xlim )
        densities[[group]] <- d_group
      }
      
      return(densities)
    }
    get_limits <- function(densities_a = NULL, densities_b = NULL){
      xlim = c()
      ylim = c(0)
      for (single_density in densities_a){
        xlim=range(c(xlim, single_density$x))
        ylim=range(c(ylim, single_density$y))
      }
      for (single_density in densities_b){
        xlim=range(c(xlim, single_density$x))
        ylim=range(c(ylim, single_density$y))
      }
      
      return(list(xlim=xlim, ylim=ylim))
    }
    plot_distributions <- function(data_groups = NULL, xlim=NULL, ylim=NULL, densities = NULL, field_group= NULL, field_value = "value", xlab_modifier = "", var_value= NULL, include_peak_dots=TRUE, dataset_name = NULL, ...){
      data_groups$group_of_interest <- data_groups[[field_group]]
      data_groups$value_of_interest <- data_groups[[field_value]]
      
      # Main body of plot_decoy_distribution_by_field_of_interest()
      FIXED_LWD=3
      
      groups <- sort(unique(data_groups$group_of_interest))
      n      <- length(groups)
      
      df_leg <- data.frame(stringsAsFactors = FALSE,
                           list(leg = groups,
                                col = rainbow_with_fixed_intensity(n = n, goal_intensity_0_1 = 0.4),
                                lty = rep(1:6, length.out=n),
                                lwd = rep(FIXED_LWD , n)) )
      
      d <- densities[["All decoys"]]
      
      xlab = sprintf("Precursor Mass Discrepancy%s (ppm)", xlab_modifier)
      ylab = "Density"
      
      if (!include_text){
        xlab=""
        ylab=""
      }
      plot(d, lwd=FIXED_LWD * scale, main=main, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, xaxt="n", yaxt="n")
      
      ave_peak <- max(d$y)
      max_peak <- 0
      
      for (local_group in groups){
        data_subset <- subset(data_groups, group_of_interest == local_group)
        data_info   <- subset(df_leg     , leg               == local_group)
        col <- data_info$col[1]
        lty <- data_info$lty[1]
        lwd <- data_info$lwd[1]
        if (nrow(data_subset) > 100){
          d <- densities[[local_group]]  #density(data_subset[[field_value]])
          lines(d, col=col, lty=lty, lwd=lwd * scale)
          peak <- max(d$y)
          max_peak <- max(max_peak, peak)
        }
      }
      abline(v=0, h=0, lwd=scale)
      leg <- list(title = "Peptide length (aa)", 
                  leg = c("All decoys"     , df_leg$leg),
                  col = c(col2hex("black") , df_leg$col),
                  lty = c(1                , df_leg$lty),
                  lwd = c(FIXED_LWD        , df_leg$lwd)
      )
      if (include_text){
        legend_object = Legend_Object$new(leg, scale)
        legend_object$show()
        #first_legend = legend(x="topleft", title = "", legend = leg$leg, col = leg$col, lty = leg$lty, lwd = leg$lwd, seg.len=leg$seg.len, box.lwd=leg$box.lwd, x.intersp = leg$x.intersp, y.intersp = leg$y.intersp)
        #new_x = first_legend$rect$left 
        #new_y = first_legend$rect$top + first_legend$rect$h * .07 * (2 - scale)
        #legend(x=new_x, y=new_y, title = leg$title, legend = "", cex=1.15, bty="n")
      }
      
      box(lwd=scale) #box around plot area
      
    }
    
    # Main body for plot_density_PMD_and_norm_decoy_by_aa_length()
    
    data_mod <- add_group_peptide_length_special()
    data_mod <- subset(data_mod, decoy==1)
    
    densities_a <- get_densities(data_subset = data_mod, field_value = "value"     , field_group = "group_peptide_length_special")
    densities_b <- get_densities(data_subset = data_mod, field_value = "value_norm", field_group = "group_peptide_length_special")
    
    dataset_name <- data_collection$dataset_name
    
    limits <- get_limits(densities_a, densities_b)
    xlim   <- limits$xlim
    ylim   <- limits$ylim
    
    if (show_norm){
      plot_distributions(data_groups = data_mod, densities=densities_b, field_value = "value_norm", xlab_modifier = " - normalized", field_group = "group_peptide_length_special", dataset_name=dataset_name, xlim=xlim, ylim=ylim)
    } else {
      plot_distributions(data_groups = data_mod, densities=densities_a, field_value = "value"     , xlab_modifier = ""             , field_group = "group_peptide_length_special", dataset_name=dataset_name, xlim=xlim, ylim=ylim)
    }
  }
  
)

###############################################################################
######## Class: Plot_Bad_CI
###############################################################################
Plot_Bad_CI = setRefClass("Plot_Bad_CI",
                          contains = "Plot_Image",
                          fields = list(breaks = "numeric",
                                        ylim   = "numeric"))
Plot_Bad_CI$methods(
  initialize = function(p_breaks=20, p_ylim=NULL, ...){
    if (is.null(p_ylim)){
      ylim <<- numeric(0)
    } else {
      ylim <<- p_ylim
    }
    breaks <<- p_breaks
    plot_title <<- "Credible Intervals for proportion within range - bad"
    callSuper(...)
  },
  get_n = function(){
    return(nrow(subset(data_collection$data_groups, (decoy == 1))))
  },
  plot_image = function(){
    data_groups <- data_collection$data_groups
    data_decoy <- subset(data_groups, data_groups$group_training_class == "bad_long")
    data_decoy$region <- cut(x = data_decoy$value, breaks = breaks)
    table(data_decoy$region)
    regions <- unique(data_decoy$region)
    
    N = nrow(data_decoy)
    find_lower_ci_bound <- function(x){
      ci <- credible_interval(length(x), N, precision = 0.001, alpha=0.05)
      return(ci[1])
    }
    find_upper_ci_bound <- function(x){
      ci <- credible_interval(length(x), N, precision = 0.001, alpha=0.05)
      return(ci[2])
    }
    xleft   <- aggregate(value~region, data=data_decoy, FUN=min)
    xright  <- aggregate(value~region, data=data_decoy, FUN=max)
    ytop    <- aggregate(value~region, data=data_decoy, FUN=find_upper_ci_bound)
    ybottom <- aggregate(value~region, data=data_decoy, FUN=find_lower_ci_bound)
    
    xleft   <- rename_column(xleft  , "value", "xleft"  )
    xright  <- rename_column(xright , "value", "xright" )
    ytop    <- rename_column(ytop   , "value", "ytop"   )
    ybottom <- rename_column(ybottom, "value", "ybottom")
    
    boxes <- merge(merge(xleft, xright), merge(ytop, ybottom))
    
    
    xlab <- "Precursor Mass Discrepancy (ppm)"
    ylab <- "Proportion of PSMs\nin subgroup"
    xlim=range(data_decoy$value)
    get_ylim(boxes=boxes)
    if (!include_text){
      xlab=""
      ylab=""
    }
    
    plot(x=c(-10,10), y=c(0,1), type="n", ylim=ylim, xlim=xlim, xlab=xlab, ylab=ylab, main=main, xaxt="n", yaxt="n")
    
    with(boxes, rect(xleft=xleft, xright=xright, ytop=ytop, ybottom=ybottom, lwd=scale))
    
    abline(h=1/breaks, col="blue", lwd=scale)
  },
  get_ylim = function(boxes=NULL){
    is_valid_range <- function(r=NULL){
      return(length(r) == 2)
    }
    if (! is_valid_range(ylim)){
      ylim <<- range(c(0,boxes$ytop, boxes$ybottom))
    }
  }
  
)
###############################################################################
######## Class: Plot_Compare_iFDR_Confidence_1_Percent_TD_FDR
###############################################################################
Plot_Compare_iFDR_Confidence_1_Percent_TD_FDR = setRefClass("Plot_Compare_iFDR_Confidence_1_Percent_TD_FDR",
                                                            contains = "Plot_Image",
                                                            fields = list())
Plot_Compare_iFDR_Confidence_1_Percent_TD_FDR$methods(
  initialize = function( ...){
    plot_title <<- "Precursor Mass Discrepance i-FDR for 1% Target-Decoy FDR PSMs"
    callSuper(...)
  },
  get_n = function(){
    if (one_percent_calculation_exists()){
      data_subset <- subset(data_collection$i_fdr, is_in_1percent==TRUE)
      n <- nrow(data_subset)
    } else {
      n <- 0
    }
    
    return (n)
  },
  plot_image = function(){
    if (one_percent_calculation_exists()){
      i_fdr        <- get_modified_fdr()
      report_good_discrepancies(i_fdr)
      data_TD_good <- get_data_TD_good(i_fdr)
      mean_results <- get_mean_results(data_TD_good)
      boxes        <- mean_results
      boxes        <- rename_columns(df = boxes, 
                                     names_before = c("min_conf", "max_conf", "lower"  , "upper"),
                                     names_after  = c("xleft"   , "xright"  , "ybottom", "ytop" ))
      xlim <- range(boxes[,c("xleft", "xright")])
      ylim <- range(boxes[,c("ybottom", "ytop")])
      
      #head(mean_results)
      
      xlab = "Confidence Score (Peptide Shaker)"
      ylab = "Mean PMD i-FDR"
      
      if (!include_text){
        xlab=""
        ylab=""
      }
      
      plot(mean_i_fdr~mean_conf, data=mean_results, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=main, xaxt="n", yaxt="n", cex=scale, lwd=scale)
      with(boxes, rect(xleft = xleft, ybottom = ybottom, xright = xright, ytop = ytop, lwd=scale))
      abline(b=-1, a=100, lwd=4*scale, col="dark gray")
      abline(h=0, v=100, lwd=1*scale)
      
    } else {
      stop(sprintf("Dataset '%s' does not include 1%% FDR data", data_collection$dataset_name))
    }
  },
  get_mean_results = function(data_TD_good = NULL){
    mean_i_fdr <- aggregate(i_fdr~conf_group, data=data_TD_good, FUN=mean)
    mean_i_fdr <- rename_column(mean_i_fdr, "i_fdr", "mean_i_fdr")
    sd_i_fdr <- aggregate(i_fdr~conf_group, data=data_TD_good, FUN=sd)
    sd_i_fdr <- rename_column(sd_i_fdr, "i_fdr", "sd_i_fdr")
    n_i_fdr <- aggregate(i_fdr~conf_group, data=data_TD_good, FUN=length)
    n_i_fdr <- rename_column(n_i_fdr, "i_fdr", "n")
    mean_conf <- aggregate(confidence~conf_group, data=data_TD_good, FUN=mean)
    mean_conf <- rename_column(mean_conf, "confidence", "mean_conf")
    min_conf <- aggregate(confidence~conf_group, data=data_TD_good, FUN=min)
    min_conf <- rename_column(min_conf, "confidence", "min_conf")
    max_conf <- aggregate(confidence~conf_group, data=data_TD_good, FUN=max)
    max_conf <- rename_column(max_conf, "confidence", "max_conf")
    
    mean_results <-                     mean_i_fdr
    mean_results <- merge(mean_results, sd_i_fdr)
    mean_results <- merge(mean_results, n_i_fdr)
    mean_results <- merge(mean_results, mean_conf)
    mean_results <- merge(mean_results, min_conf)
    mean_results <- merge(mean_results, max_conf)
    
    mean_results$se    <- with(mean_results, sd_i_fdr / sqrt(n - 1))
    mean_results$lower <- with(mean_results, mean_i_fdr - 2*se)
    mean_results$upper <- with(mean_results, mean_i_fdr + 2*se)
    return(mean_results)
  },
  get_data_TD_good = function(i_fdr=NULL){
    data_TD_good <- subset(i_fdr, TD_good==TRUE)
    data_TD_good <- data_TD_good[order(data_TD_good$confidence),]
    n <- nrow(data_TD_good)
    data_TD_good$conf_group <- cut(1:n, breaks=floor(n/100))
    data_TD_good$i_fdr <- 100 * data_TD_good$i_fdr
    return(data_TD_good)
  },
  get_modified_fdr = function(){
    i_fdr <- data_collection$i_fdr
    i_fdr$PMD_good  <- i_fdr$i_fdr < 0.01
    i_fdr$TD_good   <- i_fdr$is_in_1percent == TRUE
    i_fdr$conf_good <- i_fdr$confidence == 100
    return(i_fdr)
  },
  one_percent_calculation_exists = function(){
    return( "is_in_1percent" %in% colnames(data_collection$i_fdr))
  },
  report_good_discrepancies = function(i_fdr=NULL){
    with(subset(i_fdr,                                        (decoy == 0)), print(table(TD_good, PMD_good)))
    with(subset(i_fdr, (confidence==100)                    & (decoy == 0)), print(table(TD_good, PMD_good)))
    with(subset(i_fdr, (confidence>= 99) & (confidence<100) & (decoy == 0)), print(table(TD_good, PMD_good)))
    with(subset(i_fdr, (confidence>= 99) & (confidence<100) & (decoy == 0)), print(table(TD_good, PMD_good)))
    with(subset(i_fdr, (confidence>= 90) & (confidence< 99) & (decoy == 0)), print(table(TD_good, PMD_good)))
  }
  
)

###############################################################################
######## Class: Plot_Selective_Loss
###############################################################################
Plot_Selective_Loss = setRefClass("Plot_Selective_Loss",
                                  contains = "Plot_Image",
                                  fields = list())
Plot_Selective_Loss$methods(
  initialize = function( ...){
    plot_title <<- "PMD-FDR Selectively removes Bad Hits"
    callSuper(...)
  },
  get_n = function(){
    data_subset <- data_collection$i_fdr
    return(nrow(data_subset))
  },
  plot_image = function(){
    # Support functions for plot_selective_loss()
    
    samples_lost_by_threshold <- function(updated_i_fdr=NULL, score_threshold=NULL){
      data_subset <- subset(updated_i_fdr, confidence >= score_threshold)
      tbl <- with(updated_i_fdr, 
                  table(confidence >= score_threshold, 
                        new_confidence < score_threshold, 
                        group_decoy_proteins))
      df <- data.frame(tbl)
      df_n <- aggregate(Freq~group_decoy_proteins+Var1, data=df, FUN=sum)
      df_n <- rename_column(df_n, name_before = "Freq", "n")
      df <- merge(df, df_n)
      df$rate_of_loss <- with(df, Freq/n)
      df <- subset(df, (Var1==TRUE) & (Var2==TRUE))
      df <- df[,c("group_decoy_proteins", "rate_of_loss", "n", "Freq")]
      if (nrow(df) > 0){
        df$score_threshold <- score_threshold
      }
      return(df)
    }
    get_loss_record <- function(updated_i_fdr=NULL, score_thresholds=NULL){
      df=data.frame()
      for (score_threshold in score_thresholds){
        df_new_loss <- samples_lost_by_threshold(updated_i_fdr, score_threshold)
        df <- rbind(df, df_new_loss)
      }
      return(df)
    }
    
    # Main code for plot_selective_loss()
    
    updated_i_fdr                <- data_collection$i_fdr
    updated_i_fdr$new_confidence <- with(updated_i_fdr, 100 * (1-i_fdr)) #ifelse((1-i_fdr) < (confidence / 100), (1-i_fdr), (confidence/100)))
    loss_record <- get_loss_record(updated_i_fdr=updated_i_fdr, score_thresholds = 1:100)
    xlim <- with(loss_record, range(score_threshold))
    ylim <- c(0,1)
    xlab <- "Fixed Confidence threshold (PeptideShaker score)"
    ylab <- "Rate of PSM disqualification from PMD-FDR"
    lwd  <- 4
    plot(x=xlim, y=ylim, type="n", main=main, xlab=xlab, ylab=ylab)
    
    groups <- sort(unique(loss_record$group_decoy_proteins))
    n_g    <- length(groups)
    
    cols <- rainbow_with_fixed_intensity(n=n_g, goal_intensity_0_1 = 0.5, alpha = 1)
    ltys <- rep(1:6, length.out=n_g)
    
    leg     <- list(leg=groups, col=cols, lty=ltys, lwd=lwd, title="Species/Category")
    
    for (i in 1:n_g){
      lines(rate_of_loss~score_threshold, data=subset(loss_record, group_decoy_proteins==leg$leg[i]), col=leg$col[i], lwd=leg$lwd * scale, lty=leg$lty[i])
    }
    abline(h=0, v=100, lwd=scale)
    abline(h=c(0.1, 0.8), col="gray", lwd=scale)
    
    #leg = list(leg=group, col=col, lty=lty, lwd=lwd)
    #with(leg, legend(x = "topleft", legend = group, col = col, lty = lty, lwd = lwd, seg.len = seg.len))
    legend_object <- Legend_Object$new(leg, scale)
    legend_object$show()
  }
  
)
###############################################################################
######## Class: Plot_Dataset_Description
###############################################################################
Plot_Dataset_Description = setRefClass("Plot_Dataset_Description",
                                       contains = "Plot_Multiple_Images",
                                       fields = list(ylim_time_invariance = "numeric"))
Plot_Dataset_Description$methods(
  initialize = function(p_data_collection = NULL, 
                        p_include_text=TRUE, 
                        p_include_main=FALSE,
                        p_ylim_time_invariance = c(-4,4), ...){
    plot_object_r1_c1 <- Plot_Time_Invariance_Alt$new(p_data_collection = p_data_collection, 
                                                      p_include_text=p_include_text, 
                                                      p_include_main=p_include_main,
                                                      p_training_class = "good_testing",
                                                      p_field_of_interest = "value",
                                                      p_ylim = p_ylim_time_invariance)
    
    plot_object_r1_c2 <- Plot_Time_Invariance_Alt$new(p_data_collection = p_data_collection, 
                                                      p_include_text=p_include_text, 
                                                      p_include_main=p_include_main,
                                                      p_training_class = "good_testing",
                                                      p_field_of_interest = "value_norm",
                                                      p_ylim = p_ylim_time_invariance)
    plot_object_r2_c1 <- Plot_Density_PMD_by_Score$new(p_data_collection = p_data_collection, 
                                                       p_show_norm=FALSE, 
                                                       p_include_text=p_include_text, 
                                                       p_include_main=p_include_main)
    
    plot_object_r2_c2 <- Plot_Density_PMD_and_Norm_Decoy_by_AA_Length$new(p_data_collection = p_data_collection, 
                                                                          p_show_norm=FALSE,
                                                                          p_include_text=p_include_text, 
                                                                          p_include_main=p_include_main)
    
    plot_object_r3_c1 <- Plot_Density_PMD_by_Score$new(p_data_collection = p_data_collection, 
                                                       p_show_norm=TRUE, 
                                                       p_include_text=p_include_text, 
                                                       p_include_main=p_include_main)
    plot_object_r3_c2 <- Plot_Density_PMD_and_Norm_Decoy_by_AA_Length$new(p_data_collection = p_data_collection, 
                                                                          p_show_norm=TRUE,
                                                                          p_include_text=p_include_text, 
                                                                          p_include_main=p_include_main)
    callSuper(p_n_images_wide=2, 
              p_n_images_tall=3, 
              p_include_text=p_include_text,
              p_include_main=p_include_main,
              p_image_list = list(plot_object_r1_c1, plot_object_r1_c2,
                                  plot_object_r2_c1, plot_object_r2_c2,
                                  plot_object_r3_c1, plot_object_r3_c2), ...)
    
  }
)
###############################################################################
######## Class: Plot_Density_PMD_by_Score
###############################################################################
Plot_Density_PMD_by_Score = setRefClass("Plot_Density_PMD_by_Score",
                                        contains = "Plot_Image",
                                        fields = list(show_norm = "logical"))
Plot_Density_PMD_by_Score$methods(
  initialize = function(p_show_norm=FALSE, ...){
    show_norm <<- p_show_norm
    plot_title <<- "PMD distribution, by Confidence ranges"
    callSuper(...)
    
  },
  get_n = function(){
    return(nrow(data_collection$data_groups))
    #data_subset <- data_collection$i_fdr
    #return(nrow(data_subset))
  },
  get_modified_data_groups = function(var_value = NULL){
    # Note: Filters out used_to_find_middle
    # Note: Creates "value_of_interest" field
    # Note: Remakes "group_decoy_confidence" field
    data_new                   <- data_collection$data_groups
    data_new                   <- subset(data_new, !used_to_find_middle )
    data_new$value_of_interest <- data_new[, var_value]
    
    cutoff_points <- c(100, 100, 95, 80, 50, 0, 0)
    n <- length(cutoff_points)
    uppers <- cutoff_points[-n]
    lowers <- cutoff_points[-1]
    
    for (i in 1:(n-1)){
      upper <- uppers[i]
      lower <- lowers[i]
      
      
      if (lower==upper){
        idx <- with(data_new, which(                        (confidence == upper) & (decoy == 0)))
        cat_name <- sprintf("%d", upper)
      } else {
        idx <- with(data_new, which((confidence >= lower) & (confidence <  upper) & (decoy == 0)))
        cat_name <- sprintf("%02d - %2d", lower, upper)
      }
      data_new$group_decoy_confidence[idx] <- cat_name
    }
    
    return(data_new)
  },
  plot_image = function(){
    
    # Support functions for plot_density_PMD_by_score()
    
    get_densities <- function(data_subset = NULL, var_value = NULL){
      
      # Support functions for get_densities()
      
      # New version
      
      # Main body of get_densities()
      
      data_subset <- get_modified_data_groups(var_value=var_value)
      #data_subset$value_of_interest <- data_subset[,var_value]
      from <- min(data_subset$value_of_interest)
      to   <- max(data_subset$value_of_interest)
      xlim = range(data_subset$value_of_interest)     
      
      groups   <- sort(unique(data_subset$group_decoy_confidence), decreasing = TRUE)
      n_groups <- length(groups)
      
      densities <- list(var_value = var_value, groups=groups)
      
      for (i in 1:n_groups){
        group <- groups[i]
        
        data_group_single  <- subset(data_subset, (group_decoy_confidence == group))
        d_group            <- with(data_group_single , density(value_of_interest, from = from, to = to))
        d_group            <- normalize_density(d_group)
        
        densities[[group]] <- d_group
      }
      
      return(densities)
      
    }
    get_xlim <- function(densities_a = NULL, densities_b = NULL){
      groups <- densities_a$groups
      
      xlim <- 0
      for (group in groups){
        xlim <- range(xlim, densities_a[[group]]$x, densities_b[[group]]$x)
      }
      
      return(xlim)
      
    }
    get_ylim <- function(densities_a = NULL, densities_b = NULL){
      groups <- densities_a$groups
      
      ylim <- 0
      for (group in groups){
        ylim <- range(ylim, densities_a[[group]]$y, densities_b[[group]]$y)
      }
      
      return(ylim)
      
    }
    plot_distributions <- function(densities = NULL, var_value= NULL,include_peak_dots=TRUE, xlab_modifier="", xlim=NULL, ylim=NULL, ...){
      data_groups <- get_modified_data_groups(var_value=var_value)
      groups      <- sort(unique(data_groups$group_decoy_confidence))
      n_groups    <- length(groups)
      
      groups_std   <- setdiff(groups, c("100", "decoy", "0") )
      groups_std   <- sort(groups_std, decreasing = TRUE)
      groups_std   <- c(groups_std, "0")
      n_std        <- length(groups_std)
      cols <- rainbow_with_fixed_intensity(n = n_std, goal_intensity_0_1 = 0.5, alpha=0.5)
      
      leg <- list(group = c("100"             , groups_std   , "decoy"                           ),
                  leg   = c("100"             , groups_std   , "All Decoys"                      ),
                  col   = c(col2hex("black")  , cols         , col2hex("purple", col_alpha = 0.5)), 
                  lwd   = c(4                 , rep(2, n_std), 4                                 ), 
                  title = "Confidence Score")
      
      xlab = sprintf("Precursor Mass Discrepancy%s (ppm)",
                     xlab_modifier)
      ylab = "Density"
      if (!include_text){
        xlab=""
        ylab=""
      }
      
      
      plot( x=xlim, y=ylim, col=leg$col[1], lwd=leg$lwd[1] * scale, main=main, xlab=xlab, ylab=ylab, xaxt="n", yaxt="n", cex=scale, type="n")#, lty=leg.lty[1], ...)
      
      include_peak_dots = FALSE # BUGBUG: Disabling this for now.  Need to move this to class parameter
      
      for (i in 1:length(leg$group)){
        group <- leg$group[i]
        d     <- densities[[group]]
        lines(d, col=leg$col[i], lwd=leg$lwd[i] * scale)
        if (include_peak_dots){
          x=d$x[which.max(d$y)]
          y=max(d$y)
          points(x=c(x,x), y=c(0,y), pch=19, col=leg$col[i], cex=scale)
        }
      }
      
      abline(v=0, lwd=scale)
      
      if (include_text){
        legend_object = Legend_Object$new(leg, scale)
        legend_object$show()
      }
      
    }
    
    # Main body for plot_density_PMD_by_score()
    
    data_groups <- data_collection$data_groups
    
    data_subset_a <- subset(data_groups  , used_to_find_middle == FALSE)
    data_subset_b <- subset(data_subset_a, peptide_length > 11)
    
    densities_a <- get_densities(data_subset = data_subset_a, var_value = "value")        
    densities_b <- get_densities(data_subset = data_subset_b, var_value = "value_norm")
    
    xlim=get_xlim(densities_a, densities_b)
    ylim=get_ylim(densities_a, densities_b)
    
    dataset_name <- data_collection$dataset_name
    if (show_norm){
      plot_distributions(densities=densities_b, var_value = "value_norm", xlab_modifier = " - normalized", xlim=xlim, ylim=ylim)
    } else {
      plot_distributions(densities=densities_a, var_value = "value"     , xlab_modifier = ""             , xlim=xlim, ylim=ylim)
    }
  }
)
###############################################################################
######## Class: Plots_for_Paper
###############################################################################
Plots_for_Paper <- setRefClass("Plots_for_Paper", fields =list(list_of_data_collections = "list",
                                                               data_collection_a = "list",
                                                               data_collection_b = "list",
                                                               data_collection_c = "list",
                                                               data_collection_d = "list",
                                                               include_text      = "logical",
                                                               include_main      = "logical", 
                                                               mai               = "numeric"))
Plots_for_Paper$methods(
  report_PMD_confidence_comparison_from_subset = function(data_subset=NULL, group_name=NULL){
    print(group_name)
    print(sprintf("    Number of PSMs: %d", nrow(data_subset)))
    mean_confidence <- mean(data_subset$confidence)
    print(sprintf("    Mean Confidence Score: %3.1f", mean_confidence))
    print(sprintf("    PeptideShaker g-FDR: %3.1f", 100-mean_confidence))
    mean_PMD_FDR = mean(data_subset$i_fdr)
    print(sprintf("    PMD g-FDR: %3.1f", 100*mean_PMD_FDR))
    #col <- col2hex("black", 0.2)
    #plot(data_subset$i_fdr, pch=".", cex=2, col=col)
    #abline(h=0)
  },
  report_comparison_of_Confidence_and_PMD = function (i_fdr = NULL, min_conf=NULL, max_conf=NULL, include_max=FALSE){
    
    if (is.null(max_conf)) {
      data_subset <- subset(i_fdr, confidence == min_conf)
      group_name <- sprintf("Group %d", min_conf)
    } else if (include_max){
      data_subset <- subset(i_fdr, (confidence >= min_conf) & (confidence <= max_conf))
      group_name <- sprintf("Group %d through %d", min_conf, max_conf)
    } else {
      data_subset <- subset(i_fdr, (confidence >= min_conf) & (confidence < max_conf))
      group_name <- sprintf("Group %d to %d", min_conf, max_conf)
    }
    
    report_PMD_confidence_comparison_from_subset(data_subset=data_subset, group_name=group_name)
  },
  create_stats_for_grouping_figure = function(data_collection=NULL){
    
    aug_i_fdr                      <- data_collection$i_fdr
    aug_i_fdr$group_good_bad_other <- gsub("_.*", "", aug_i_fdr$group_training_class) 
    aug_i_fdr$group_null           <- "all"
    table(aug_i_fdr$group_training_class)
    table(aug_i_fdr$group_good_bad_other)
    table(aug_i_fdr$group_null)
    
    create_agg_fdr_stats <- function(i_fdr=NULL, grouping_var_name = NULL){
      formula_fdr <- as.formula(sprintf("%s~%s", "i_fdr", grouping_var_name))
      formula_len <- as.formula(sprintf("%s~%s", "peptide_length", grouping_var_name))
      agg_fdr <- aggregate(formula=formula_fdr, data=i_fdr, FUN=mean)
      agg_n   <- aggregate(formula=formula_fdr, data=i_fdr, FUN=length)
      agg_len <- aggregate(formula=formula_len, data=i_fdr, FUN=mean)
      agg_fdr <- rename_columns(df = agg_fdr, 
                                names_before = c(grouping_var_name, "i_fdr"), 
                                names_after  = c("group"          , "fdr"))
      agg_n   <- rename_columns(df = agg_n, 
                                names_before = c(grouping_var_name, "i_fdr"), 
                                names_after  = c("group"          , "n"))
      agg_len <- rename_columns(df = agg_len, 
                                names_before = c(grouping_var_name), 
                                names_after  = c("group"          ))
      agg <- merge(agg_fdr, agg_n)
      agg <- merge(agg    , agg_len)
      
      return(agg)
    }
    
    agg_detail  <- create_agg_fdr_stats(i_fdr = aug_i_fdr, grouping_var_name = "group_training_class")
    agg_grouped <- create_agg_fdr_stats(i_fdr = aug_i_fdr, grouping_var_name = "group_good_bad_other")
    agg_all     <- create_agg_fdr_stats(i_fdr = aug_i_fdr, grouping_var_name = "group_null")
    
    agg <- rbind(agg_detail, agg_grouped)
    agg <- rbind(agg, agg_all)
    
    agg$fdr <- ifelse(agg$fdr < 1, agg$fdr, 1)
    
    linear_combo <- function(x=NULL, a0=NULL, a1=NULL){
      result <- (a0 * (1-x) + a1 * x)
      return(result)
    }
    
    agg$r <- linear_combo(agg$fdr, a0=197, a1= 47)
    agg$g <- linear_combo(agg$fdr, a0= 90, a1= 85)
    agg$b <- linear_combo(agg$fdr, a0= 17, a1=151)
    
    return(agg)
  },
  report_ranges_of_comparisons = function(data_collection=NULL){
    i_fdr <- data_collection$i_fdr
    print(sprintf("PMD and Confidence comparison for -- %s", data_collection$dataset_name))
    report_comparison_of_Confidence_and_PMD(i_fdr = i_fdr, min_conf=100, max_conf=NULL, include_max=TRUE)
    report_comparison_of_Confidence_and_PMD(i_fdr = i_fdr, min_conf= 99, max_conf=100 , include_max=FALSE)
    report_comparison_of_Confidence_and_PMD(i_fdr = i_fdr, min_conf= 90, max_conf= 99 , include_max=FALSE)
    report_comparison_of_Confidence_and_PMD(i_fdr = i_fdr, min_conf=  0, max_conf=100 , include_max=TRUE)
  },
  print_figure_2_data = function(){
    print(create_stats_for_grouping_figure(data_collection = data_collection_a))
  },
  print_table_4_data = function(){
    report_ranges_of_comparisons(data_collection=data_collection_a)
    report_ranges_of_comparisons(data_collection=data_collection_c)
  },
  plot_figure_3 = function(p_scale=NULL, p_include_main=NULL){
    plot_object <- Plot_Compare_PMD_and_Norm_Density$new(p_data_collection = plots_for_paper$data_collection_a,
                                                         p_show_norm       = FALSE,
                                                         p_include_text    = TRUE,
                                                         p_include_main    = p_include_main,
                                                         p_display_n_psms  = FALSE)
    plot_object$plot_image_in_small_window(p_scale=p_scale)
  },
  plot_figure_4 = function(p_scale=NULL, p_include_main=NULL){
    plot_object <- Plot_Time_Invariance_Alt_Before_and_After$new(p_data_collection = plots_for_paper$data_collection_a, 
                                                                 p_include_text=TRUE, 
                                                                 p_include_main=p_include_main,
                                                                 p_ylim = c(-4,4))
    plot_object$plot_image_in_large_window(window_height=4, p_scale=p_scale)
    
  },
  plot_figure_5 = function(p_scale=NULL, p_include_main=NULL){
    plot_object <- Plot_Density_PMD_and_Norm_Decoy_by_AA_Length$new(p_data_collection = plots_for_paper$data_collection_a, 
                                                                    p_include_text=TRUE, 
                                                                    p_include_main=p_include_main)
    plot_object$plot_image_in_large_window(window_height=4, p_scale=p_scale)
  },
  plot_figure_6 = function(p_scale=NULL, p_include_main=NULL){
    plot_object <- Plot_Bad_CI$new(p_data_collection = plots_for_paper$data_collection_a, 
                                   p_include_text=TRUE, 
                                   p_include_main=p_include_main)
    plot_object$plot_image_in_small_window(p_scale=p_scale)
  },
  plot_figure_7 = function(p_scale=NULL, p_include_main=NULL){
    plot_object <- Plot_Compare_iFDR_Confidence_1_Percent_TD_FDR$new(p_data_collection = plots_for_paper$data_collection_a, 
                                                                     p_include_text=TRUE, 
                                                                     p_include_main=p_include_main)
    plot_object$plot_image_in_large_window(window_height=4, p_scale=p_scale)
  },
  plot_figure_8 = function(p_scale=NULL, p_include_main=NULL){
    plot_object <- Plot_Selective_Loss$new(p_data_collection = plots_for_paper$data_collection_c, 
                                           p_include_text=TRUE, 
                                           p_include_main=p_include_main)
    plot_object$plot_image_in_large_window(window_height=4, p_scale=p_scale)
  },
  plot_figure_A = function(p_scale=NULL, p_include_main=NULL){
    plot_object <- Plot_Dataset_Description$new(p_data_collection=plots_for_paper$data_collection_a, 
                                                p_include_text=TRUE,
                                                p_include_main=p_include_main,
                                                p_ylim_time_invariance=c(-4,4) )
    plot_object$plot_image_in_large_window(window_height=4, p_scale=p_scale)
  },
  plot_figure_B = function(p_scale=NULL, p_include_main=NULL){
    plot_object <- Plot_Dataset_Description$new(p_data_collection=plots_for_paper$data_collection_b, 
                                                p_include_text=TRUE,
                                                p_include_main=p_include_main,
                                                p_ylim_time_invariance=c(-4,4) )
    plot_object$plot_image_in_large_window(window_height=4, p_scale=p_scale)
  },
  plot_figure_C = function(p_scale=NULL, p_include_main=NULL){
    plot_object <- Plot_Dataset_Description$new(p_data_collection=plots_for_paper$data_collection_c, 
                                                p_include_text=TRUE,
                                                p_include_main=p_include_main,
                                                p_ylim_time_invariance=c(-4,4) )
    plot_object$plot_image_in_large_window(window_height=4, p_scale=p_scale)
  },
  plot_figure_D = function(p_scale=NULL, p_include_main=NULL){
    plot_object <- Plot_Dataset_Description$new(p_data_collection=plots_for_paper$data_collection_d, 
                                                p_include_text=TRUE,
                                                p_include_main=p_include_main,
                                                p_ylim_time_invariance=c(-4,4) )
    plot_object$plot_image_in_large_window(window_height=4, p_scale=p_scale)
  },
  
  
  create_data_collections = function(list_of_data_collections=NULL){
    get_big_list_of_data_collections <- function(){
      make_file_names <- function(data_collection=NULL){
        dataset_experiment_name  = data_collection$dataset_experiment_name
        dataset_designation = data_collection$dataset_designation
        file_name_dataset   = data_collection$file_name_dataset
        dir_dataset         = data_collection$dir_dataset
        dir_results         = data_collection$dir_results
        
        DIR_GRIFFIN     = "."
        DIR_PROJECTS    = ".."
        
        dataset_name                 <- sprintf("%s_%s", dataset_experiment_name, dataset_designation)#  "Sectioning_01"
        #file_name_clean_data         <- sprintf("clean_%s"         , dataset_name)
        #file_name_CI                 <- sprintf("CI_results_%s.rds", dataset_name)
        file_name_collection         <- sprintf("collection_%s.rds", dataset_name)
        file_name_groups             <- sprintf("groupings_%s.rds", dataset_name)
        file_name_densities          <- sprintf("densities_%s.rds", dataset_name)
        file_name_i_fdr              <- sprintf("i_fdr_%s.rds", dataset_name)
        #file_name_plots              <- sprintf("plots_%s.pdf", dataset_name)
        data_collection$dataset_name         <- dataset_name
        data_collection$file_path_data       <- file.path(dir_dataset, file_name_dataset)
        #data_collection$file_path_clean_data <- file.path(DIR_PROJECTS, dir_results, file_name_clean_data)
        #data_collection$file_path_CI         <- file.path(DIR_PROJECTS, dir_results, file_name_CI        )
        data_collection$file_path_collection <- file.path(dir_results, file_name_collection)
        data_collection$file_path_groups     <- file.path(dir_results, file_name_groups)
        data_collection$file_path_densities  <- file.path(DIR_PROJECTS, dir_results, file_name_densities)
        data_collection$file_path_i_fdr      <- file.path(DIR_PROJECTS, dir_results, file_name_i_fdr)
        #data_collection$file_path_plots      <- file.path(DIR_PROJECTS, dir_results, file_name_plots)
        
        if (! is.null(data_collection$file_name_dataset_1_percent)){
          file_name_1_percent <- data_collection$file_name_dataset_1_percent
          data_collection$file_path_dataset_1_percent <- file.path(dir_dataset, file_name_1_percent)
        }
        
        data_collection$file_names_made <- TRUE
        
        return(data_collection)
      }
      
      data_collection_group_all = list()
      data_collection_group_all$dir_results = file.path(".", "Results")
      data_collection_group_all$dir_dataset = file.path(".", "Data")
      
      # Data collections: Mouse Mutations
      data_collection_group_mutations = data_collection_group_all
      data_collection_group_mutations$dataset_experiment_name  = "Mouse Mutations"
      data_collection_mouse_combined_05 = data_collection_group_mutations
      data_collection_mouse_combined_05$dataset_designation = "combined_05" ######################################### This one
      data_collection_mouse_combined_05$file_name_dataset   = "Combined_DB_Mouse_5PTM.tabular"
      
      # Data collections: Pyrococcus
      data_collection_group_pyrococcus = data_collection_group_all
      data_collection_group_pyrococcus$dataset_experiment_name  = "Pyrococcus"
      data_collection_pyrococcus_tr = data_collection_group_pyrococcus  ######################################### This one
      data_collection_pyrococcus_tr$dataset_designation = "tr"
      data_collection_pyrococcus_tr$file_name_dataset   = "Pfu_traditional_Extended_PSM_Report.tabular"
      
      # Data collections: Oral_737_NS
      data_collection_group_oral_737_NS = data_collection_group_all
      data_collection_group_oral_737_NS$dataset_experiment_name  = "Oral_737_NS"
      data_collection_oral_737_NS_combined = data_collection_group_oral_737_NS
      data_collection_oral_737_NS_two_step = data_collection_group_oral_737_NS
      data_collection_oral_737_NS_combined$dataset_designation = "combined"  ######################################### This one
      data_collection_oral_737_NS_two_step$dataset_designation = "two_step"  ######################################### This one
      data_collection_oral_737_NS_combined$file_name_dataset           = "737_NS_Peptide_Shaker_Extended_PSM_Report_CombinedDB.tabular"
      data_collection_oral_737_NS_two_step$file_name_dataset           = "737_NS_Peptide_Shaker_Extended_PSM_Report_Multi_Stage_Two_Step.tabular.tabular"
      data_collection_oral_737_NS_combined$file_name_dataset_1_percent = "737_NS_Peptide_Shaker_PSM_Report_CombinedDB.tabular"
      data_collection_oral_737_NS_two_step$file_name_dataset_1_percent = "737_NS_Peptide_Shaker_PSM_Report_Multi_Stage_Two_Step.tabular"
      
      data_collection_mouse_combined_05    <- make_file_names(data_collection_mouse_combined_05   )
      data_collection_pyrococcus_tr        <- make_file_names(data_collection_pyrococcus_tr       )
      data_collection_oral_737_NS_combined <- make_file_names(data_collection_oral_737_NS_combined)
      data_collection_oral_737_NS_two_step <- make_file_names(data_collection_oral_737_NS_two_step)
      
      # List of all datasets
      big_list <- list(data_collection_mouse_combined_05,
                       data_collection_pyrococcus_tr,
                       data_collection_oral_737_NS_combined,
                       data_collection_oral_737_NS_two_step
      )
      
      return(big_list)
    }
    create_data_collection <- function(big_list = NULL, experiment_name = NULL, designation = NULL){
      
      get_collection = function(big_list=NULL, experiment_name=NULL, designation=NULL){
        
        for (data_collection in big_list){
          if ((data_collection$dataset_designation     == designation) &
              (data_collection$dataset_experiment_name == experiment_name)){
            
            return(data_collection)
          }
        }
        stop(sprintf("Data collection not found: %s_%s",
                     experiment_name,
                     designation))
      }
      process_data <- function(data_collection=NULL){
        data_original <- data_collection$data_original
        verify_data_frame(data_original)
        #data_collection$data_groups <- make_data_groups(data_original)
        return(data_collection)
      }
      
      data_collection <- get_collection(big_list, 
                                        experiment_name = experiment_name, 
                                        designation     = designation)
      #data_collection <- make_file_names(data_collection)
      #data_collection <- load_data(      data_collection)
      data_collection <- create_data_collection_from_scratch(data_collection = data_collection, save_results = FALSE)
      
      data_collection <- process_data(data_collection = data_collection)
      return(data_collection)
    }
    big_list <- get_big_list_of_data_collections()
    
    data_collection_a <<- create_data_collection(big_list, 
                                                 experiment_name = "Oral_737_NS", 
                                                 designation     = "two_step")
    data_collection_b <<- create_data_collection(big_list,
                                                 experiment_name = "Oral_737_NS",
                                                 designation     = "combined")
    data_collection_c <<- create_data_collection(big_list,
                                                 experiment_name = "Pyrococcus",
                                                 designation     = "tr")
    data_collection_d <<- create_data_collection(big_list,
                                                 experiment_name = "Mouse Mutations",
                                                 designation     = "combined_05")
    
    
  },
  create_plots_for_paper = function(include_main=TRUE, finalize=TRUE){
    plot_figure_3(p_scale=ifelse(finalize, 4, 1), p_include_main = include_main)
    plot_figure_4(p_scale=ifelse(finalize, 2, 1), p_include_main = include_main)
    plot_figure_5(p_scale=ifelse(finalize, 2, 1), p_include_main = include_main)
    plot_figure_6(p_scale=ifelse(finalize, 4, 1), p_include_main = include_main)
    plot_figure_7(p_scale=ifelse(finalize, 2, 1), p_include_main = include_main)
    plot_figure_8(p_scale=ifelse(finalize, 2, 1), p_include_main = include_main)
    plot_figure_A(p_scale=ifelse(finalize, 2, 1), p_include_main = include_main)
    plot_figure_B(p_scale=ifelse(finalize, 2, 1), p_include_main = include_main)
    plot_figure_C(p_scale=ifelse(finalize, 2, 1), p_include_main = include_main)
    plot_figure_D(p_scale=ifelse(finalize, 2, 1), p_include_main = include_main)
    print_figure_2_data()
    print_table_4_data()
  }
)
