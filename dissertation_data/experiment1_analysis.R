
##################
## Experiment 1 ##
##################

## load libraries
library(dplyr)
library(tidyr)
library(ggplot2)

functions_url <- 'https://raw.githubusercontent.com/kleinmichaeldavid/phd/master/dissertation_data/accuracy_and_outlier_functions.R'
source(functions_url)

## open the data
data_url <- 'https://raw.githubusercontent.com/kleinmichaeldavid/phd/master/dissertation_data/experiment1_anon.txt'
raw_exp1 <- read.table(data_url,
                       fileEncoding = "UCS2", header = TRUE,
                       stringsAsFactors = FALSE, sep = "\t")

column_url <- 'https://raw.githubusercontent.com/kleinmichaeldavid/phd/master/dissertation_data/experiment0_colnames.txt'
col_names <- read.table(column_url,header=TRUE) #info on cols to keep

## pare down the data to the neccessities
raw_exp1 <- raw_exp1[,names(raw_exp1) %in% col_names$old]
names(raw_exp1) <- col_names$new

## remove practice trials and extract info from version column
raw_exp1 <- raw_exp1[!(raw_exp1$block %in% 1:2),]
exts <- c("PerceptionOfTime_DIFFAT2_1st","PerceptionOfTime_DIFFBT2_2nd",
          "PerceptionOfTime_DIFFCT2_1st","PerceptionOfTime_DIFFDT2_2nd")
raw_exp1$est_type[raw_exp1$version %in% exts] <- "External"
raw_exp1$est_type[!(raw_exp1$version %in% exts)] <- "Internal"

raw_exp1$acc <-  raw_exp1 %>% with(resp == cresp)
raw_exp1$difficulty[raw_exp1$subtrahend > 10] <- 'high'
raw_exp1$difficulty[raw_exp1$subtrahend < 10] <- 'low'

raw_exp1 <- raw_exp1[raw_exp1$subject != 677696150,] # this subject didn't complete the experiment


## accuracy
acc_info <- check_accuracies(raw_exp1, 'acc', 'subject', cut_off = 0.7)
low_acc_subs <- acc_info[[1]]
acc_summary <- acc_info[[2]]
raw_exp1 <- raw_exp1[!(raw_exp1$subject %in% low_acc_subs),] # remove inaccurate subjects
accurate_exp1 <- raw_exp1[raw_exp1$acc,] # remove inaccurate trials

####### OUTLIERS #######

# trial level
accurate_exp1 <- outlier_trials(accurate_exp1,
                                variables=c("rt"),
                                separators=c("subject","timing_required","est_type"))
inliers_exp1 <- accurate_exp1[!(accurate_exp1$low_outlier1 %in% TRUE) &
                                !(accurate_exp1$high_outlier1 %in% TRUE),]

# subject level
subject_outliers <- outlier_subjects(inliers_exp1, variables = c("rt"),
                                     true_separator = "subject",
                                     other_separators = c("timing_required","est_type"),
                                     cutoff = 2.5)[[1]]

clean_exp1 <- inliers_exp1[!(inliers_exp1$subject %in% c(subject_outliers)),] 


####### RT ANALYSIS #######

# start without difficulty
rt_exp1 <- summarizer(clean_exp1,
                      variable = "rt",
                      rows = "subject",
                      columns = c("timing_required","est_type"),
                      within_err = TRUE)
# the plot of RTs
bar_plotter(rt_exp1,
            within_ivs = c("timing_required", "est_type"),
            within_types = c("fill","x"),
            within_order = list(c("NoTime","Time"),c("External","Internal")),
            xlab = "Session Type",
            ylab = "RT to Equation Verification Task (ms)",
            yrange=c(0,3000),
            leg_labels=c("No Time","Time"),
            leg_title="Timing",
            #fill_colours=c("#39d2dd", "#5bef20", "#dd3862"),
            leg_pos=c(0.4,0.98))


# again with difficulty in there
rt_exp1_wdiff <- summarizer(clean_exp1,
                            variable = "rt",
                            rows = "subject",
                            columns = c("timing_required","est_type","difficulty"),
                            within_err = TRUE)

bar_plotter(rt_exp1_wdiff, ###### change order of entre variables instead of levels within variables?
            within_ivs = c("timing_required", "est_type","difficulty"),
            within_types = c("alpha","x","fill"),
            within_order = list(c("NoTime","Time"),c("External","Internal"),c("low","high")),
            xlab = "Session Type",
            ylab = "RT to Equation Verification Task (ms)",
            yrange=c(0,5000),
            leg_labels=c("Low","High"),
            leg_title="Equation Difficulty",
            #fill_colours=c("#39d2dd", "#5bef20", "#dd3862"),
            leg_pos=c(0.9,1))

####### ACCURACY ANALYSIS #######

acc_exp1 <- summarizer(raw_exp1,
                       variable = "acc",
                       rows = "subject",
                       columns = c("timing_required","est_type","difficulty"),
                       within_err = TRUE)

bar_plotter(acc_exp1, ###### change order of entre variables instead of levels within variables?
            within_ivs = c("timing_required", "est_type","difficulty"),
            within_types = c("alpha","x","fill"),
            within_order = list(c("NoTime","Time"),c("External","Internal"),c("low","high")),
            xlab = "Session Type",
            ylab = "Accuracy to Equation Verification Task",
            yrange=c(0,1),
            leg_labels=c("Low","High"),
            leg_title="Equation Difficulty",
            #fill_colours=c("#39d2dd", "#5bef20", "#dd3862"),
            leg_pos=c(0.5,0.4))


####### TIME ANALYSIS #######

## clean up
internal_trials <- clean_exp1$est_type == "Internal"
external_trials <- clean_exp1$est_type == "External"
clean_exp1$actual[internal_trials] <- clean_exp1$rt[internal_trials]
clean_exp1$actual[external_trials] <- clean_exp1$stim_duration[external_trials]
clean_exp1$estimate <- clean_exp1$estimate %>% as.numeric()

## matching procedure
matches_exp1 <- match_distributions(clean_exp1[clean_exp1$timing_required != "NoTime",],
                                    match_within = c("subject"),
                                    match_across = c("est_type"),
                                    match_variable = "actual",
                                    estimate_variable = "estimate",
                                    proximity_cut_off = 0.1)

## trim subjects without enough matched pairs
enough_matches_exp1 <- trimmer(matches_exp1, column = "count",
                               cutoff = 10, extend_to = "subject")
## plot actual means
actual_means_exp1 <- summarizer(enough_matches_exp1,
                                variable = "actual_mean",
                                rows = c("subject"),
                                columns = c("condition"),
                                within_err = TRUE)

bar_plotter(actual_means_exp1,within_ivs = c("sess_type"),
            within_types = c("x"),
            within_order = list(c("External","Internal")),
            xlab = "Session Type", ylab = "Mean Actual Duration (ms)",
            yrange = c(0,3000),
            fill_colours = NULL)

## plot estimated means
estimated_means_exp1 <- summarizer(enough_matches_exp1,
                                   variable = "estimate_mean",
                                   rows = c("subject"),
                                   columns = c("condition"),
                                   within_err = TRUE)

bar_plotter(estimated_means_exp1,within_ivs = c("sess_type"),
            within_types = c("x"),
            within_order = list(c("External","Internal")),
            xlab = "Session Type", ylab = "Mean Estimated Duration (ms)",
            yrange = c(0,3000),
            fill_colours = NULL) ## maybe put a line across here for actual mean

## plot correlations
r_exp1 <- summarizer(enough_matches_exp1,
                     variable = "r",
                     rows = c("subject"),
                     columns = c("condition"),
                     within_err = TRUE)

bar_plotter(r_exp1,within_ivs = c("gtedfededegt"),## !!!
            within_types = c("x"),
            within_order = list(c("External","Internal")),
            xlab = "Session Type", ylab = "correlation",
            yrange = c(0,1),
            fill_colours = NULL)




