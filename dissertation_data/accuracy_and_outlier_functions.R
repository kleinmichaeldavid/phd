


slow_resp_cutoff <- function(trials_per_block, num_blocks, pval = 0.05, prob_slower = 0.5){
  
  #########################
  ##
  ## Assuming the probability of a participant responding more slowly than the 
  ## stimulus follows a binomial distribution, this function returns the number of 
  ## slow responses corresponding to the inputted p-value. 
  ##
  ## Inputs:
  ## trials_per_block - number of trials in a block
  ## num_blocks - number of blocks of interest
  ## pval - desired alpha level
  ## prob_slower - probability of a slow response if participants follow instructions
  ## 
  #########################
  
  p_1_block <- pbinom(c(0:trials_per_block), trials_per_block, prob_slower, lower.tail = FALSE)
  p_x_blocks <- 1-(1-p_1_block)^num_blocks
  cut_off <- (which(p_x_blocks < pval)[1]-1) / trials_per_block
  print(paste("Remove participants with a proportion of slow trials greater than",cut_off,"in one block."))
  cut_off
}

find_slow_subjects <- function(df, subject_col, rt_col, duration_col, condition_cols, which_conditions, cutoff){
  
  ##########################
  ##
  ## Returns list of subjects who responded slowly more often than
  ## the inputted cutoff.
  ##
  ## Inputs:
  ## df - dataframe containing the data
  ## subject_col - name of column containing subject labels
  ## rt_col - name of column containing RT
  ## duration_col - name of column containing stimulus durations
  ## condition_cols - name of column containing condition labels
  ## which_conditions - the conditions that this function should work on
  ## cutoff - cutoff number of trials
  ##
  ##########################
  
  df$slow_resp <- df[,rt_col] > df[,duration_col]
  summary <- df %>% group_by_at(c(subject_col,condition_cols)) %>%
    summarize(slow_prop = mean(slow_resp)) %>% as.data.frame()
  for(i in 1:length(condition_cols)){
    summary <- summary[summary[,condition_cols[i]] %in% which_conditions[[i]],]
  }
  slow_responders <- unique(summary$subject[summary$slow_prop > cutoff])
  slow_responders
}





check_accuracies <- function(df, variable = 'acc', subject = 'subjects', cut_off = 0.7){
  
  ##############################################################
  ## 
  ## Returns a list containing an array of the low accuracy
  ## participants as well as a dataframe summarizing the
  ## accuracy of all participants. Assumes we're interested
  ## in overall accuracy rather than in any given condition.
  ##
  ## Inputs:
  ## df - dataframe containing data of interest
  ## variable - name of column containing accuracy
  ## subject - name of column containing subject labels
  ## cut_off - accuracy cutoff to be included in the array of 
  ##            low accuracy participants
  ## 
  ###############################################################
  
  
  ## should add in conditions
  
  library(rlang)
  
  # get acc for each participant
  accuracy_df <- group_by(df,!!sym(subject)) %>% summarize(acc_mean = mean(!!sym(variable), na.rm = TRUE)) %>% as.data.frame()
  
  # create list of participants who fail to reach the cut-off
  low_acc_subjects <- accuracy_df[,subject][which(accuracy_df$acc_mean < cut_off)]
  
  # return
  list(low_acc_subjects, accuracy_df)
}


outlier_trials <- function(df, variables, separators, cutoff=2.5){
  
  ###############################################
  ##
  ## Given a dataframe containing experimental data,
  ## returns that dataframe with additional columns
  ## indicating whether a given trial as a
  ## high/low outlier.
  ##
  ## Inputs:
  ## df - dataframe containing the data
  ## variables - a list of variables to check for outliers
  ## separators - a list of variables to group the outlier check by
  ##  (i.e. find outliers only WITHIN levels of these variables)
  ## cutoff - number of s.d. from the mean
  ##              to be considered an outlier
  ## 
  ################################################

  len <- length(variables)
  means_and_sds <- list() #actually no reason to store as list instead of overwriting
  for(var in 1:len){
    means_and_sds[[var]] <- df %>% group_by_at(separators) %>%
      summarize(mean = mean(!!sym(variables[var])), sd = sd(!!sym(variables[var]))) %>%
      as.data.frame() %>%
      mutate(lower_cutoff = mean - cutoff*sd) %>%
      mutate(upper_cutoff = mean + cutoff*sd)
    
    trial_condition_args <- c(df[,separators],sep="")
    summary_condition_args <- c(means_and_sds[[var]][,separators],sep="")
    trial_summary_match <- match(do.call(paste,trial_condition_args),do.call(paste,summary_condition_args))
    df[,paste("low_outlier",var,sep="")] <- df[,variables[var]] < means_and_sds[[var]]$lower_cutoff[trial_summary_match]
    df[,paste("high_outlier",var,sep="")] <- df[,variables[var]] > means_and_sds[[var]]$upper_cutoff[trial_summary_match]
  }
  df
}


outlier_subjects <- function(df, variables, true_separator, other_separators, cutoff=2.5, only_conditions = NULL){
  
  ##########################################
  ##
  ## Returns list of subjects who were subject-level outliers
  ## 
  ## Inputs:
  ## df - dataframe containing the data
  ## variables - a list of variables to check for outliers (i.e. probably "RT")
  ## true_separator - the column in df representing what is being tested
  ##      for whether it is an outlier (i.e. probably "subject")
  ## other_separators - columns columns indicating the various experimental
  ##      conditions that must be checked
  ## cutoff - number of s.d. from the mean
  ##      to be considered an outlier
  ## 
  ######################################
  
  len <- length(variables)
  condition_means <- list() #actually no reason to store as list instead of overwriting
  outlier_subjects <- list()
  for(var in 1:len){
    condition_means[[var]] <- df %>% group_by_at(c(true_separator,other_separators)) %>%
      summarize(mean = mean(!!sym(variables[var]))) %>%
      as.data.frame()
    
    if(length(other_separators)>0){
      condition_means[[var]] <- condition_means[[var]] %>%
        unite('condition',other_separators) %>% spread(condition, mean)
    } #else {
    #condition_means[[var]] <- condition_means[[var]] %>% spread(true_separator, mean)
    #}
    
    only_conditions <- c(true_separator, only_conditions)
    if (!is.null(only_conditions)) condition_means[[var]] <- condition_means[[var]][,only_conditions]
    
    subject_mean_df <- data.frame(subject = condition_means[[var]][,true_separator],
                                  mean = rowMeans(data.frame(condition_means[[var]][,-1])),
                                  outlier = rep(FALSE, nrow(condition_means[[var]])))
    subject_sd <- sd(subject_mean_df$mean)
    subject_mean <- mean(subject_mean_df$mean)
    
    subject_mean_df$outlier <- subject_mean_df$mean < subject_mean - cutoff * subject_sd | subject_mean_df$mean > subject_mean + cutoff * subject_sd
    outlier_count <- sum(subject_mean_df$outlier)
    new_outliers <- outlier_count - 0
    while(new_outliers > 0){
      subject_sd <- sd(subject_mean_df$mean[!subject_mean_df$outlier])
      subject_mean <- mean(subject_mean_df$mean[!subject_mean_df$outlier])
      subject_mean_df$outlier <- subject_mean_df$mean < subject_mean - cutoff * subject_sd | subject_mean_df$mean > subject_mean + cutoff * subject_sd
      new_outliers <- sum(subject_mean_df$outlier) - outlier_count
      outlier_count <- sum(subject_mean_df$outlier)
    }
    
    outlier_subjects[[var]] <- subject_mean_df$subject[subject_mean_df$outlier]
  }
  
  outlier_subjects
}

summarizer <- function(df, variable, rows, columns, within_err = TRUE){
  
  ########################################
  ##
  ## Given trial-level data, returns a dataframe of 
  ## means for each subject in each condition.
  ##
  ## Inputs: 
  ## df - dataframe of trial-level data
  ## variable - variable to summarize
  ## rows - column representing the variable to 
  ##    group within (probably "subject")
  ## columns - names of columns that contain labels for
  ##    the various conditions
  ## within_err - boolean representing whether between-subjects
  ##    variability should be removed (e.g. TRUE if you want to
  ##    create a plot including within-subjects confidence intervals)
  ## 
  #########################################
  
  conditions_summary <- df %>% group_by_at(c(rows,columns)) %>%
    summarize(mean = mean(!!sym(variable))) %>%
    as.data.frame() %>%
    unite(condition, columns, remove=TRUE) %>%
    spread(condition, mean)
  
  for(i in ncol(conditions_summary):1){
    if(mean(is.na(conditions_summary[,i]))==1){
      conditions_summary <- conditions_summary[,-i]
    }
  }
  
  if(within_err==TRUE){
    # remove between subjects error
    row_means <- rowMeans(conditions_summary[,-1])
    overall_mean <- mean(row_means)
    adjustment <- replicate(ncol(conditions_summary)-1, row_means) - overall_mean
    conditions_summary[,-1] <- conditions_summary[,-1] - adjustment
  }
  
  conditions_summary
}


bar_plotter <- function(df,
                        within_ivs = c(),
                        between_ivs = c(),
                        within_types = c(),
                        between_types = c(),
                        within_order = NULL,
                        between_order = NULL,
                        xgroups = NULL,
                        xlab = NULL, ylab = NULL,
                        yrange = NULL,
                        leg_title = NULL,
                        leg_labels = NULL,
                        leg_pos = c(0.2,0.2),
                        fill_colours = NULL){
  
  #############
  ##
  ## Returns a bar plot that uses the formatting
  ## I'd like to use for my dissertation.
  ##
  ## Inputs:
  ## df - dataframe containing summary of group means for 
  ##    each condition (e.g. as from the 'summarizer' function)
  ## within_ivs - names of columns containing within subjects IVs
  ## between_ivs - names of columns containing between subjects IVs
  ## within_types - how you want each variable represented in the plot
  ##      (e.g. 'x' = as the x variable, 'fill' = as the fill colour)
  ## between_types - see 'within_types'
  ## within_order - desired order of the within subjects variables
  ## between_order - desired order of the between subjects variables
  ## ...others presumable self explanatory
  ##
  #############
  
  ## preprocessing
  ivs <- c(within_ivs, between_ivs)
  iv_class <- c(rep("within",length(within_ivs)),rep("between",length(between_ivs)))
  iv_type = c(within_types, between_types)
  
  between_vars <- ivs[which(iv_class %in% "between")]
  between_indices <- which(names(df) %in% between_vars)
  
  num_within_levels <- ncol(df[,-c(1,between_indices)])
  num_between_levels <- 1
  for(var in between_vars){
    num_between_levels <- num_between_levels * length(unique(df[,var]))
  }
  
  within_count <- length(within_ivs)
  between_count <- length(between_ivs)
  iv_count <- within_count + between_count
  plot_df <- data.frame(matrix(NA, nrow = num_within_levels * num_between_levels, ncol = iv_count+2))
  names(plot_df) <- c(ivs,"mean","se")
  
  ## fill in the within var names
  for(i in 1:within_count){
    #within_index <- sum(iv_class[1:i] %in% "within")
    plot_df[,i] <- df[,-c(1,between_indices)] %>% names() %>% strsplit("_") %>%
      lapply("[",i) %>% unlist() %>% factor()
    
    if(is.null(within_order)){
      ord <- levels(plot_df[,i])
    } else {
      ord <- within_order[[i]]
    }
    
    plot_df[,i] <- factor(plot_df[,i], levels <- ord)
  }
  
  ## fill in the between var names
  if(between_count>0){
    levels_remaining <- num_between_levels
    for(i in (1+within_count):(iv_count)){
      
      levels <- df[,ivs[i]] %>% unique()
      levels_remaining <- levels_remaining / length(levels)
      plot_df[,i] <- as.factor(rep(levels,times=levels_remaining,each=num_within_levels*(i-within_count)))
      
      if(is.null(between_order)){
        ord <- levels(plot_df[,i])
      } else {
        ord <- between_order[i-within_count]
      }
      
      plot_df[,i] <- factor(plot_df[,i], levels <- ord[[1]])
      
    }
  }
  
  if(between_count>0){
    
    ## dump NA cols
    plot_df <- plot_df[,!is.na(names(plot_df))]
    
    
    ## fill in the means and SEs
    merge_df <- plot_df %>% unite(new_col, between_vars, remove = TRUE) %>% unique()
    between_trialtypes <- merge_df$new_col
    working_df <- df %>% unite(new_col, between_vars, remove = FALSE)
    
    for(b_tt in between_trialtypes){
      cur_df <- working_df[working_df$new_col %in% b_tt,]
      plot_df[which(merge_df$new_col %in% b_tt),]$mean <- apply(cur_df[-c(1:(2+between_count))], 2, mean)
      plot_df[which(merge_df$new_col %in% b_tt),]$se <- apply(cur_df[-c(1:(2+between_count))], 2, sd) / sqrt(nrow(cur_df))
    }
  } else if(between_count ==0){
    plot_df[,iv_count+1] <- apply(df[,-1], 2, mean)
    plot_df[,iv_count+2] <- apply(df[,-1], 2, sd)/ sqrt(nrow(df))
    
    ## dump NA cols
    plot_df <- plot_df[,!is.na(names(plot_df))]
  }
  
  ## plotting part
  limits <- aes(ymax = plot_df$mean + plot_df$se,
                ymin = plot_df$mean - plot_df$se)
  
  if ("fill" %in% iv_type){
    filler = ivs[match("fill",iv_type)]
  } else {
    filler = NULL
  }
  
  if ("x" %in% iv_type) xer = ivs[match("x",iv_type)] else xer = NULL
  if ("alpha" %in% iv_type) alpher = ivs[match("alpha",iv_type)] else alpher = NULL
  
  p <- ggplot(data=plot_df,
              aes_string(y="mean",
                         x=xer,
                         fill=filler,
                         alpha = alpher))
  
  p <- p + geom_bar(stat = "identity",
                    position = position_dodge(0.9), colour="black",
                    width = 0.9, lwd = 0.7) +
    geom_errorbar(limits, position = position_dodge(0.9),
                  width = 0.1, lwd = 1) + theme_classic() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text = element_text(size=16, colour = "black"),
          axis.text.x=element_text(colour="black", size = 16),
          axis.text.y=element_text(colour="black", size = 16),
          axis.title=element_text(colour="black", size = 16, face = "bold"),
          legend.text=element_text(size=16),
          legend.justification=c(1,1),
          legend.position=leg_pos,
          axis.line=element_line(size = 1.2))
  
  # plot axes
  if(!is.null(xlab)) p <- p + labs(x=xlab)
  if(!is.null(ylab)) p <- p + labs(y=ylab)
  if(!is.null(yrange)) p <- p + coord_cartesian(ylim=yrange)
  
  # plot legend
  if(is.null(leg_title)) leg_title <- ""
  p <- p + guides(fill=guide_legend(title=leg_title))
  
  # x groups
  if(!is.null(xgroups)) p <- p + scale_x_discrete(labels = xgroups)
  
  if(!is.null(fill_colours)){
    if(!is.null(leg_labels)){
      p <- p + scale_fill_manual(values=fill_colours, labels=leg_labels)
    } else {
      p <- p + scale_fill_manual(values=fill_colours)
    }
  } else {
    if(!is.null(leg_labels)) p <- p + scale_fill_discrete(labels=leg_labels)
  }
  
  ## x labels
  ## change order of x axis
  
  p
}

match_distributions <- function(df, match_within, match_across,
                                match_variable, estimate_variable,
                                proximity_cut_off){
  
  #########################################
  ##
  ## Given a dataframe containing experimental data,
  ## returns the dataframe, matches up sets of rows
  ## based on their similarity to one another on one
  ## variable, and returns the dataframe, but with
  ## only similar items.
  ##
  ## Inputs:
  ## df - dataframe containing data
  ## match_within - name of column containing variable to be
  ##      matched within (i.e. a given pair/triplet must both have
  ##      the same value within this column) - probably "subject"
  ## match_across - name of column containing conditions
  ##      to be matched across (i.e. a given pair/triplet must
  ##      cover all possible values of this column)
  ## match_variable - the IV being matched
  ## estimate_variable - the DV
  ## proximity_cut_off - required similarity for trials to be
  ##      considered a match
  ## 
  ########################################
  
  #set up df
  final_names <- c(match_within,"condition","count","actual_mean",
                   "estimate_mean","actual_sd","estimate_sd","r")
  matched_df <- data.frame(matrix(rep(NA,7+length(match_within)),nrow=1))
  names(matched_df) <- final_names
  matched_df <- matched_df[-1,]
  
  
  groups <- unique(df[,match_across])
  n_groups <- length(groups)
  
  # merge the "match_within" variables into single column that can be cycled through
  merge_df <- df %>% unite(new_col, match_within, remove = TRUE)
  
  subsets <- unique(merge_df$new_col)
  matched_data <- setNames(vector("list", length(subsets)),subsets)
  for(i in 1:length(subsets)){
    
    currdata <- merge_df[merge_df$new_col == subsets[i],]
    
    matches <- list()
    for(n in 1:n_groups){
      matches[[n]] <- currdata[currdata[,match_across] == groups[n],][,match_variable]
    }
    
    if(n_groups == 3){
      
      lens <- list(length(matches[[1]]),length(matches[[2]]),length(matches[[3]]))
      arrays <- list(array(rep(matches[[1]],each=1,times=lens[[2]]*lens[[3]]),
                           dim = unlist(lens)),
                     array(rep(matches[[2]],each=lens[[1]],times=lens[[3]]),
                           dim = unlist(lens)),
                     array(rep(matches[[3]],each=lens[[1]]*lens[[2]],times=1),
                           dim = unlist(lens)))
      ratios <- (abs(arrays[[1]]-arrays[[2]])+
                   abs(arrays[[2]]-arrays[[3]])+
                   abs(arrays[[3]]-arrays[[1]]))/
        (arrays[[1]]+arrays[[2]]+arrays[[3]])
      
      ## now get the smallest triplets from the ratio array
      matches_df <- data.frame(a = as.numeric(), b = as.numeric(), c = as.numeric())
      while(min(ratios) < proximity_cut_off){
        ## find indices of the smallest values in the ratio matrix
        minima <- which(ratios == min(ratios), arr.ind = TRUE) 
        minimum <- minima[sample(1:nrow(minima),1),] # pick randomly if more than one
        matches_df <- rbind(matches_df,
                            data.frame(a = minimum[1], b = minimum[2], c = minimum[3]))
        ## remove the trials corresponding the the triplet from the array
        ratios[minimum[1],,] <- rep(Inf, length(ratios[minimum[1],,]))
        ratios[,minimum[2],] <- rep(Inf, length(ratios[,minimum[2],]))
        ratios[,,minimum[3]] <- rep(Inf, length(ratios[,,minimum[3]]))
      }
      
    } else if (n_groups == 2){

      ratios <- abs(outer(matches[[1]],matches[[2]],"-")) /
        (outer(matches[[1]],matches[[2]],"+")/2)
      
      ## now get the smallest pairs from the ratio array
      matches_df <- data.frame(a = as.numeric(), b = as.numeric())
      while(min(ratios) < proximity_cut_off){
        ## find indices of the smallest values in the ratio matrix
        minima <- which(ratios == min(ratios), arr.ind = TRUE) 
        minimum <- minima[sample(1:nrow(minima),1),] # pick randomly if more than one
        matches_df <- rbind(matches_df,
                            data.frame(a = minimum[1], b = minimum[2]))
        ## remove the trials corresponding the the triplet from the array
        ratios[minimum[1],] <- rep(Inf, length(ratios[minimum[1],]))
        ratios[,minimum[2]] <- rep(Inf, length(ratios[,minimum[2]]))
      }
    } else if (n_groups == 4){
      lens <- list(length(matches[[1]]),length(matches[[2]]),length(matches[[3]]),length(matches[[4]]))
      arrays <- list(array(rep(matches[[1]],each=1,times=lens[[2]]*lens[[3]]*lens[[4]]),
                           dim = unlist(lens)),
                     array(rep(matches[[2]],each=lens[[1]],times=lens[[3]]*lens[[4]]),
                           dim = unlist(lens)),
                     array(rep(matches[[3]],each=lens[[1]]*lens[[2]],times=lens[[4]]),
                           dim = unlist(lens)),
                     array(rep(matches[[4]],each=lens[[1]]*lens[[2]]*lens[[3]],times=1),
                           dim = unlist(lens)))
      
      ratios <- (abs(arrays[[1]]-arrays[[2]]) + abs(arrays[[1]]-arrays[[3]]) + abs(arrays[[1]]-arrays[[4]]) +
                   abs(arrays[[2]]-arrays[[3]]) + abs(arrays[[2]]-arrays[[4]]) + abs(arrays[[3]]-arrays[[4]])) /
        (arrays[[1]]+arrays[[2]]+arrays[[3]]+arrays[[4]])
      
      ## now get the smallest triplets from the ratio array
      matches_df <- data.frame(a = as.numeric(), b = as.numeric(), c = as.numeric(), d = as.numeric())
      while(min(ratios) < proximity_cut_off){
        ## find indices of the smallest values in the ratio matrix
        minima <- which(ratios == min(ratios), arr.ind = TRUE)
        minimum <- minima[sample(1:nrow(minima),1),] # pick randomly if more than one
        matches_df <- rbind(matches_df, data.frame(a = minimum[1], b = minimum[2], c = minimum[3], d = minimum[4]))
        ## remove the trials corresponding the the triplet from the array
        ratios[minimum[1],,,] <- rep(Inf, length(ratios[minimum[1],,,]))
        ratios[,minimum[2],,] <- rep(Inf, length(ratios[,minimum[2],,]))
        ratios[,,minimum[3],] <- rep(Inf, length(ratios[,,minimum[3],]))
        ratios[,,,minimum[4]] <- rep(Inf, length(ratios[,,,minimum[4]]))
        }
    }
    
    counts <- nrow(matches_df)
    actuals <- list()
    estimates <- list()
    for(n in 1:n_groups){
      actuals[[n]] <- currdata[currdata[,match_across] == groups[n],][,match_variable][matches_df[,n]]
      estimates[[n]] <- currdata[currdata[,match_across] == groups[n],][,estimate_variable][matches_df[,n]]
    }
    #subset_data <- list(count = counts, actual = actuals, estimated = estimates)
    #matched_data[[i]] <- subset_data
    
    ## now analyse
    actual_means <- actuals %>% lapply(mean) %>% unlist()
    estimated_means <- estimates %>% lapply(mean) %>% unlist()
    actual_sds <- actuals %>% lapply(sd) %>% unlist()
    estimated_sds <- estimates %>% lapply(sd) %>% unlist()
    rs <- c()
    for(j in 1:length(actuals)){
      rs[j] <- cor(actuals[[j]],estimates[[j]])
    }
    
    df_add <- strsplit(subsets[i],"_")[[1]] %>% rep(each=n_groups) %>%
      matrix(nrow=n_groups,dimnames=list(NULL,match_within)) %>%
      as.data.frame() %>% cbind(data.frame(condition=groups,
                                           count=counts,
                                           actual_mean=actual_means,
                                           estimate_mean=estimated_means,
                                           actual_sd=actual_sds,
                                           estimate_sd=estimated_sds,
                                           r=rs))
    matched_df <- rbind(matched_df,df_add)
    
  }
  #matched_data
  matched_df
}

trimmer <- function(df, column, cutoff, remove_if = "<", extend_to){ 
  
  ############################################
  ##
  ## Trims data from dataframe when value of some variable
  ## doesn't meet the given cutoff.
  ##
  ## Inputs:
  ## df - dataframe containing experimental data
  ## column - name of column to use as cutoff variable
  ## cutoff - the cutoff
  ## remove_if - function to use to compare the column and cutoff
  ## extend_to - name of column to extend the cutoff to,
  ##      if a row is removed, then all other rows with
  ##      the same value in this column are also removed
  ##      --currently must include exactly 1 col name--
  ##
  #############################################
  
  missed_tf <- sapply(paste(df[,column],remove_if,cutoff), function(x) eval(parse(text=x)))
  to_remove <- unique(df[missed_tf,extend_to])
  df[!(df[,extend_to] %in% to_remove),]
}

get_actuals <- function(df, by, levels, columns){
  
  #############################################
  ##
  ## For cleaning up time estimation data in which
  ## participants estimated multiple types of intervals.
  ## Given a df of data, returns the df with a new
  ## column represented the to-be-estimated duration
  ## for each trial.
  ##
  ## Inputs:
  ## df - dataframe containing the data
  ## by - name of column containing the different
  ##    conditions in which participants estimated
  ##    different types of durations
  ## levels - the names of the different conditions
  ## columns - the names of the columns corresponding
  ##    to the durations in each condition (! must be
  ##    in the proper order given the order inputted
  ##    to the 'levels' variable !)
  ## 
  #############################################
  
  for(i in 1:length(levels)){
    trial_type <- df[,by] == levels[i]
    df$actual[trial_type] <- df[,columns[i]][trial_type]
  }
  
  df
}

ttest <- function(in1,in2,paired = TRUE){
  
  model <- t.test(in1, in2, paired = paired)

  df <- model$parameter
  t <- sprintf('%.2f',abs(model$statistic))
  p <- sprintf('%.3f',model$p.value) %>% substr(2,5)
  se <- sprintf('%.1f',sd(in1 - in2) / sqrt(length(in1)))
  sd1 <- sd(in1)
  sd2 <- sd(in2)
  r <- cor(in1, in2)
  drm <- sprintf('%.2f', abs(model$estimate) / sqrt(sd1^2 + sd2^2 - 2*r*sd1*sd2) * sqrt(2*(1-r)))

  print(paste('t(',df,') = ',t,', SE = ', se, ', p = ', p, ', drm = ', drm, sep = ''))
  
  c(df, t, se, p, drm)
  
}

