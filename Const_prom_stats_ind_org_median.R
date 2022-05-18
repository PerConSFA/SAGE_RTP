###############################################################################
# Code written by Joshua Elmore (joshua.elmore@pnnl.gov for SAGE: Serine      #
# recombinase-assisted genome engineering manuscript                          #
# Note, this version of the promoter analysis worksheet is for processing     #
# prenormalized, non-log scale values.                                        #
###############################################################################
library(ggplot2)
library(plyr)
library(Rmisc)
library(dbplyr)
#library(devtools)
library(devtools, quietly=TRUE)
# some of the packages below may be unneccessary
library(growthcurver)
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(scales)
library(svglite)
library(Cairo)
library(gdtools)
library(ggnewscale)
library(viridis)
library(ggridges)
library(tidyverse)


##* modified on 7-17-21 to replace all instances of summarySE with custom 
## function summarySE_median that also adds the median value of the data
source("summarySE_median.R")

# source_gist(
# "https://gist.github.com/eliocamp/eabafab2825779b88905954d84c82b32")
# package that enables multiple variable chart shading. 
# See website for reference:
# https://eliocamp.github.io/codigo-r/2018/09/multiple-color-and-fill-scales-with-ggplot2/
# github_install("eliocamp/ggnewscale")
# install_github("easyGgplot2", "kassambara")

# Plotting references:
# https://ggplot2.tidyverse.org/reference/geom_boxplot.html
# https://www.r-graph-gallery.com/240-custom-layout-background-ggplot2.html
# http://r-statistics.co/Top50-Ggplot2-Visualizations-MasterList-R-Code.html
# https://www.r-graph-gallery.com/boxplot.html
# https://www.r-graph-gallery.com/70-boxplot-with-categories-on-multiple-lines.html
# 
# 
################################################################################
# Data import, sorting and processing                                          #
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## # #
#
# data frame is like a spreadsheet
# Columns (aka: variables) are vectors
# Rows (aka: observations) are lists and must contain an equal number of columns
# see file Pflu_promoter_full_dataframe_NGS_cutoff5_3MAD_for_R.csv for an 
# example of what the csv needs to look like before input into R.


################################################################################
## Import the data from a file in the R project working directory rename the   #
## actual promoter data file to 'd' dataframe.                                 #
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##


Const_promoter_stats_indiv_org_median <- function(
 Input_file_label = "Pflu",
  Output_file_label = "Pflu",
  num_cond = 4,
  num_bc_rep = 5,
  num_exp_rep = 4,
  UTR_cut = 3,
  Cond_cut = 4
){
  
  
 
  ##########################################################
  ##  Set function input values for testing modifications  #       
  ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## #
  #Input_file_label = "Pflu"
  #Output_file_label = "Pflu"
  #num_cond = 4
  #num_bc_rep = 5
  #num_exp_rep = 4
  #UTR_cut = 3
  #Cond_cut = 4
  
  
  ## Number of conditions tested
  num_conditions = num_cond
  ## Number of barcodes used per promoter
  num_barcode_replicates = num_bc_rep
  ## Number of experimental replicates per condition
  num_experimental_replicates = num_exp_rep
  ## Cutoff for considering promoter UTR insensitive
  UTR_cutoff = UTR_cut
  ## Cutoff for considering promoter condition insensitive
  Cond_cutoff = Cond_cut
  
  
  file_name <- 
    paste("./R_input_data/",
          Input_file_label,
          "_promoter_full_dataframe_NGS_cutoff5_3MAD_min3reps_for_R.csv",
          sep="")
  
  d <- read.table(file_name, header=TRUE, sep = ",", stringsAsFactors = FALSE)
  
  
  
  file_name2 <- 
    "./R_input_data/Promoter_Barcode_lookup_table.csv"
  BC_lookUp <- read.table(file_name2, header=TRUE, sep = ",", 
                          stringsAsFactors = FALSE)
  
  
  
  
  
################################################################################
## calculate stats for various data subsets prior to barcode outlier excision, #
## remove NA values, and sort data in some cases                               #
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## # #
  
  ## Omit promoters with NA values from the original 'd' dataframe
  d <-na.omit(d)
  
  ## remove barcodes with insufficient datapoints for analysis
  ## first generate stats for each barcode. The cutoff we will use is: less than 
  ## 3/4 of the total number of possible samples for a given barcode. 
  ## Max number of samples = (# of conditions) * (# of replicate experiments) 
  ## E.g. for P. putida this is: (3 conditions) * (4 replicates).
  ## So our cutoff here will require the number of samples needs to be >8.
  d_barcode_filter <- summarySE_median(d, measurevar = "Promoter_Strength", 
                                groupvars = "Barcode",na.rm=TRUE)
  
  NxE_cutoff = (((num_conditions * num_experimental_replicates) * 0.75) - 1)
  
  d_barcode_filter <- d_barcode_filter %>% filter(
    d_barcode_filter$N > (NxE_cutoff)
  )
## Generate a dataframe with the names added to the the barcode filter dataframe
  
  d_barcode_filter <- merge(d_barcode_filter,BC_lookUp,by ="Barcode")
  
  
  ## Generates a dataframe of barcodes from d_stats that remain after the 
  ## upstream promoter filtering step. 
  ## Also renamed the N (total counts of data points for a barcode) 
  ## to N_barcode_samples so that this data is added to d_barcode_stats
  filtered_barcodes <- subset(d_barcode_filter, select =c("Barcode","N"))
  colnames(filtered_barcodes)<- c("Barcode","N_barcode_samples")
  
################################################################################
## Removes promoters from dataframe d that failed to pass the N > 75% of
## max barcode data point threshold. Also adds the #N for the promoter as
## the variable N_total.
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## # #
  d <- merge(filtered_barcodes, d, by="Barcode")
  
  
  ## Generate statistics for the promoter strength values with samples 
  ## grouped by the name variable
  d_stats <- summarySE_median(d, measurevar="Promoter_Strength", 
                       groupvars="Name",na.rm=TRUE)
  
  ## Filters out any promoters in which the total number of samples containing 
  ## data is less than 50% of the max number of samples. Any less than this will 
  ## provide too sparse of a data set to perform the neccessary stats.
  ## Max number of samples = (# of conditions) * (# of replicate experiments) *
  ## (# of barcodes).
  ## E.g. for P. putida this is: (3 conditions) * (4 replicates) * (5 barcodes).
  ## So our cutoff will be that the number of samples needs to be >29.
  d_stats <- d_stats %>% 
    filter(d_stats$N > (
      ((num_conditions * num_experimental_replicates * num_barcode_replicates)
       / 2 ) - 1
    )
    )
  
  
  ## Generates a dataframe of promoters from d_stats that remain after the 
  ## upstream promoter filtering step. 
  ## Also renamed the N (total counts of data points for a promoter) to 
  ## N_promoter_samples so that this data is added to d_barcode_stats
  filtered_promoters <- subset(d_stats, select =c("Name","N"))
  colnames(filtered_promoters)<- c("Name","N_promoter_samples")
  
  ## Generate statistics for the promoter strength values with samples 
  ## grouped by the barcode variable
  d_barcode_stats <- summarySE_median(d, measurevar = "Promoter_Strength", 
                               groupvars = "Barcode",na.rm=TRUE)
  
  ## Generate a dataframe with the names added to the the barcode stats dataframe
  d_barcode_stats <- merge(d_barcode_stats,BC_lookUp,by ="Barcode")
  
  ## Removes promoters or barcodes lacking a value from each datafram
  d_stats <- na.omit(d_stats)
  d_barcode_stats <- na.omit(d_barcode_stats)
  
###############################################################################
# Removes promoters from d_barcode_stats that failed to pass the N > 50% of
# max promoter data point threshold. Also adds the #N for the promoter as
# the variable N_promoter_samples.
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## # #
  
## Removes promoters in d_barcode_stats that are not in the filtered_promoters 
## dataframe 
  d_barcode_stats <- merge(filtered_promoters, d_barcode_stats, by="Name")
  
################################################################################
## Removes promoters from dataframe d that failed to pass the N > 50% of
## max promoter data point threshold. Also adds the #N for the promoter as
## the variable N_total.
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## # #
  d_N_filter <- merge(filtered_promoters, d, by="Name")
  
################################################################################
## Removes barcode outliers using the 1.5x IQR threshold method to determine 
## which barcodes are potential outliers. 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
  
## Defining two functions IQRout_top and IQRout_bot to calculate a value for the
## 1.5x IQR value. The top and bottom functions return the high and low values
## beyond which any sample values are considered outliers
  IQRout_top <- function(x){
    Q3<-quantile(x,0.75)
    Q1<-quantile(x,0.25)
    IQR<-(Q3-Q1)
    top_val<- (Q3+(1.5*IQR))
    top_val
  }
  IQRout_bot <- function(x){
    Q3<-quantile(x,0.75)
    Q1<-quantile(x,0.25)
    IQR<-(Q3-Q1)
    bot_val<- (Q1-(1.5*IQR))
    bot_val
  }
  
  
## Determine the 1.5x IQR values for each promoter in the d_barcode_stats
## dataframe. Generates a new variable in the data frame for each. The ave()
## function allows us to iteratively apply the function across subsets of
## the datasat (as described in the function), rather than across the entire
## data set.
  d_barcode_stats$top_quantile <- 
    ave(d_barcode_stats$Promoter_Strength, d_barcode_stats$Name, FUN=IQRout_top) 
  d_barcode_stats$bottom_quantile <- 
    ave(d_barcode_stats$Promoter_Strength, d_barcode_stats$Name, FUN=IQRout_bot) 
  
## Assign a string value, 'No' or 'Yes', to new variable 'bc_outlier' indicating
## whether the barcode is an outlier as determined using the 1.5x IQR method
  d_barcode_stats <- d_barcode_stats %>%
    mutate(
      bc_outlier = case_when(
        d_barcode_stats$Promoter_Strength > d_barcode_stats$top_quantile |
          d_barcode_stats$Promoter_Strength <= d_barcode_stats$bottom_quantile 
        ~ "Yes",
        T ~ "No"
      )
    )
  
## Generates a dataframe that filters out barcodes from d_barcode_stats 
## that are considered outliers using the 1.x5 IQR method
  d_barcode_stats_no_outliers <- 
    d_barcode_stats %>% filter(d_barcode_stats$bc_outlier == "No")
  
## Generate a list of barcodes and their outlier status. Will use to filter out
## outlier barcodes from other dataframes
  barcode_outlier_status <- 
    subset(d_barcode_stats, select =c("Barcode","bc_outlier"))
  
################################################################################
## Labeling barcodes as outliers or not in dataframe d_N_filter, then generate 
## a subset of d_N_filter that only includes non-outlier barcodes
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
  
## Generating new dataframe that adds the bc_outlier variable and its values to
## N-filtered dataframe d_N_filter
  d_filter <- merge(barcode_outlier_status, d_N_filter, by="Barcode")
  
## Generates a copy of d_filter that retains outlier barcode data for future use
## while we are going to remove that data from d_filter below
  
  d_filter_w_outlier_bc <- d_filter 
  
## Generates a dataframe that filters out barcodes from d_barcode_stats 
## that are considered outliers using the 1.x5 IQR method
  d_filter <- 
    d_filter %>% filter(d_filter$bc_outlier == "No")
  
## These 4 lines of code generate stats that lump together all of the data
## that has the same promoter name (d_filter_stats), or lump together the data
## that has the same name, but does not lump together data that have distinct
## condition (d_filter_con_stats) or barcode values (d_filter_bar_stats, 
## d_filter_bc_stats). Note, Replicate variable refers to the barcode variant.
## Replicate number can range from 1 to 5, and each value represents one of the
## 5 barcode variants in a given promoter.
  d_filter_bc_stats <- summarySE_median(d_filter, measurevar="Promoter_Strength", 
                                 groupvars = c("Replicate","Name"), na.rm=TRUE)
  d_filter_con_stats <- summarySE_median(d_filter, measurevar="Promoter_Strength",
                                  groupvars = c("Condition","Name"), na.rm=TRUE)
  d_filter_stats <- summarySE_median(d_filter, measurevar="Promoter_Strength", 
                              groupvars = "Name", na.rm=TRUE)
  
## These few lines of codes generate stats as above with d_filter, but with
## d_filter_w_outlier_bc.
  d_filter_out_bc_stats <- summarySE_median(d_filter_w_outlier_bc, 
                                     measurevar="Promoter_Strength", 
                                     groupvars = c("Replicate","Name"), 
                                     na.rm=TRUE)
  d_filter_out_con_stats <- summarySE_median(d_filter_w_outlier_bc, 
                                      measurevar="Promoter_Strength",
                                      groupvars = c("Condition","Name"), 
                                      na.rm=TRUE)
  d_filter_out_stats <- summarySE_median(d_filter_w_outlier_bc, 
                                  measurevar="Promoter_Strength", 
                                  groupvars = "Name", 
                                  na.rm=TRUE)
  
  
## Removes NA values from data. 
  d_filter_bc_stats <-na.omit(d_filter_bc_stats)
  d_filter_con_stats <- na.omit(d_filter_con_stats)
  d_filter_stats <- na.omit(d_filter_stats)
  d_filter_out_bc_stats <-na.omit(d_filter_out_bc_stats)
  d_filter_out_con_stats <- na.omit(d_filter_out_con_stats)
  d_filter_out_stats <- na.omit(d_filter_out_stats)
  ####d_filter_bar_stats <- na.omit(d_barcode_stats)
  
##* code added on  7/18/21 to determine whether a barcode is problematic -
##* in other words, if it were the sole barcode would its inclusion have
##* substantially changed the reported promoter strength that we determined
##* by using multiple barcodes per promoter. The cutoff we will use is the 
##* UTR_cutoff value - this way we can also do a parameter sweep if desired.

##* d_sufficient is a dataframe containing promoters that have at least 2 
##* barcode's worth of data, and thus can be used for the prob. barcode analysis
  d_sufficient <- summarySE(d_filter_bc_stats, 
                            measurevar = "Promoter_Strength",
                            groupvars = "Name",
                            na.rm = TRUE)
  d_sufficient <- subset(d_sufficient, select = c("Name",
                                                  "N"))
  d_sufficient <- d_sufficient %>% filter(d_sufficient$N > 2)
  
##* generate overall means for comparison with individual barcode means
  d_problematic_means <- subset(d_filter_stats, select = c("Name", 
                                                           "Promoter_Strength"))
  d_problematic_means <- merge(d_problematic_means,d_sufficient, by="Name")
  names(d_problematic_means)[names(d_problematic_means) == "Promoter_Strength"] 
      <- "Overall_PS"
  names(d_problematic_means)[names(d_problematic_means) == "N"] <- "N_barcodes"
  ##* generate dataframe combining o
  d_problematic <- merge(d_problematic_means,d_filter_out_bc_stats,by="Name")
  
  d_problematic <- d_problematic %>%
    mutate(
      problematic_bc = case_when(
        (d_problematic$Promoter_Strength / d_problematic$Overall_PS) 
        > (UTR_cutoff/2) |
          (d_problematic$Promoter_Strength / d_problematic$Overall_PS) 
        < (1 / (UTR_cutoff/2)
        )
        
        ~ 'Yes',
        T ~ 'No'
      )
    )
    
  

 ## Sort d_stats columns by promoter strength values. Reminder, d_stats includes 
 ## values for barcodes that are outliers. Unlike d_filter_w_outlier_bc, the 
 ## stats are not performed
  d_sort_all <- d_stats[order(d_stats$Promoter_Strength,decreasing=TRUE),]
  d_sort_all$Name <- factor(d_sort_all$Name, levels = 
                              d_sort_all$Name[order(d_sort_all$Promoter_Strength,
                                                    decreasing=TRUE)])
  
  
 ## Takes barcode values and generates stats for the promoter_strength for each 
 ## 'Promoter' at the 'Name' after first generating stats within barcodes.
 ##
 ## This is distinct from similar data structues in which all replicates 
 ## for all barcodes are lumped together in the statistical analyses.
  d_barcode_comp <- summarySE_median(d_barcode_stats_no_outliers, 
                              measurevar="Promoter_Strength", 
                              groupvars = "Name", na.rm=TRUE)
  d_barcode_comp <- na.omit(d_barcode_comp)
  
 ## rename N (number of barcodes for a promoter) to N_promoter_barcodes
  names(d_barcode_comp)[names(d_barcode_comp) == "N"] <- "N_promoter_barcodes"
  
 ## Sort barcode outlier filtered promoter data by promoter strength 
 ## for use in barcharts with Name on the X-axis
  d_barcode_comp <- 
    d_barcode_comp[order(d_barcode_comp$Promoter_Strength,decreasing=TRUE),]
  d_barcode_comp$Name <- 
    factor(d_barcode_comp$Name, 
           levels = d_barcode_comp$Name[order(d_barcode_comp$Promoter_Strength,
                                              decreasing=TRUE)]) 
  
  
  
  
  
  
 ###############################################################################
 # Calculate promoter condition and 5'UTR independence statuses                 
 # ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## # 
  
  ####################################################################
  ## Promoters are considered to be independent of 5'UTR sequence if #
  ## the ratio of Max barcode RTA / Min barcode RTA for a given      #
  ## promoter is less than what is defined below as 'UTR_cutoff'     #
  ## Promoters are considered to provide consistent expression       #
  ## conditions if the ratio of Max condition RTA / Min condition    #
  ## RTA is less than what is defined as 'Cond_cutoff'               #
  # ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ###
  
################################################################################
##  screening for 5'UTR sensitivity
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
  UTR_data <- summarySE_median(d_filter_bc_stats, measurevar = "Promoter_Strength", 
                        groupvars = c("Name"), na.rm = TRUE)
  UTR_data <- subset(UTR_data, select =c("Name","Promoter_Strength"))
  UTR_data$max_RTA <- by(d_filter_bc_stats$Promoter_Strength, 
                         d_filter_bc_stats$Name, max)
  UTR_data$min_RTA <- by(d_filter_bc_stats$Promoter_Strength, 
                         d_filter_bc_stats$Name, min)
  UTR_data$UTR_max_min_ratio <- UTR_data$max_RTA / UTR_data$min_RTA
  
## Cutoff fold change between maximum RTA and minimum RTA barcode RTAs that for
## considering a promoter to have an insensitivie 5'UTR
  
  
  UTR_data <- UTR_data %>%
    mutate(
      UTR_sensitive = case_when(
        UTR_data$UTR_max_min_ratio >= UTR_cutoff
        ~ 'Yes',
        T ~ 'No'
      )
    )
  ## Subset UTR_data 
  UTR_data <- subset(UTR_data, select =c("Name",
                                         "UTR_sensitive",
                                         "UTR_max_min_ratio"))
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
  
  
  
################################################################################
##  screening for condition independent expression
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
  Cond_data <- summarySE_median(d_filter_con_stats, measurevar = "Promoter_Strength", 
                         groupvars = c("Name"), na.rm = TRUE)
  Cond_data <- subset(Cond_data, select =c("Name","Promoter_Strength"))
  Cond_data$max_RTA <- by(d_filter_con_stats$Promoter_Strength, 
                          d_filter_con_stats$Name, max)
  Cond_data$min_RTA <- by(d_filter_con_stats$Promoter_Strength, 
                          d_filter_con_stats$Name, min)
  Cond_data$Cond_max_min_ratio <- Cond_data$max_RTA / Cond_data$min_RTA
## Cutoff fold change between maximum RTA and minimum RTA barcode RTAs that for
## considering a promoter to perform independent of condition
  
  
  
  Cond_data <- Cond_data %>%
    mutate(
      Cond_sensitive = case_when(
        Cond_data$Cond_max_min_ratio >= Cond_cutoff
        ~ 'Yes',
        T ~ 'No'
      )
    )
  ## Subset Cond_data
  Cond_data <- subset(Cond_data, select =c("Name",
                                           "Cond_sensitive",
                                           "Cond_max_min_ratio"))
  ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
  
  
  ## Takes barcode values and generates stats for the promoter_strength for each 
  ## 'Promoter' at the 'Name' after first generating stats within barcodes.
  ##
  ## This is distinct from similar data structues in which all replicates 
  ## for all barcodes are lumped together in the statistical analyses.
  d_filter_bc_comp <- summarySE_median(d_filter_bc_stats, 
                                measurevar="Promoter_Strength", 
                                groupvars = "Name", na.rm=TRUE)
  d_filter_bc_comp <- na.omit(d_filter_bc_comp)
  d_filter_bc_comp <- merge(d_filter_bc_comp,Cond_data, by="Name")
  d_filter_bc_comp <- merge(d_filter_bc_comp,UTR_data, by="Name")
  
## The following code adds the consistency variable, which will be used to 
## label promoters that we consider robust for consistent expression regardless 
## of condition and independent of downstream sequence.
  d_filter_bc_comp <- d_filter_bc_comp %>%
    mutate(
      consistent = case_when(
        Cond_sensitive == "No" & UTR_sensitive == "No" ~ "(1) yes",
        Cond_sensitive == "No" & UTR_sensitive == "Yes" ~ "(2) no - UTR",
        Cond_sensitive == "Yes" & UTR_sensitive == "Yes" ~ "(2) no - UTR",
        Cond_sensitive == "Yes" & UTR_sensitive == "No" ~ "(3) no - Cond"
      )
    )
  
## The keepers dataframe contains the names of the promoters and the consistency
## classification generated above.
  keepers <- subset(d_filter_bc_comp, select = c("Name", "consistent"))
  
  ## merge keepers with various data frames to add the consistency label to
  ## different data frames
  d_filter <- merge(d_filter,keepers, by="Name")
  d_filter_bc_stats <- merge(d_filter_bc_stats, keepers, by="Name")
  d_filter_con_stats <- merge(d_filter_con_stats, keepers, by="Name")
  d_filter_stats <- merge(d_filter_stats,keepers, by="Name")
  
  d_filter_w_outlier_bc <- merge(d_filter_w_outlier_bc, keepers, by="Name")
  ##potentially remove
  d_filter_out_bc_stats <- merge(d_filter_out_bc_stats, keepers, by="Name")
  d_filter_out_con_stats <- merge(d_filter_out_con_stats, keepers, by="Name")
  d_filter_out_stats <- merge(d_filter_out_stats,keepers, by="Name")
  
################################################################################
## CREATE A DATAFRAME WITH ROBUST LOW NOISE, 5'UTR INDEPENDENT PROMOTERS      ##
## This code generates a data frame that is a subset of d_barcode_comp that   ##
## has excluded promoters that are classified as                              ##
## non_noisy_and_condition_independent                                        ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
  consistent_promoters_bc <- d_filter_bc_comp %>% filter(
    d_filter_bc_comp$consistent == "(1) yes")
  consistent_promoters <- d_filter_stats %>% filter(
    d_filter_stats$consistent == "(1) yes")
  
  ## Sort consistent_promoters columns by promoter strength values. 
  consistent_promoters <- 
    consistent_promoters[order(consistent_promoters$Promoter_Strength,
                               decreasing=TRUE),]
  consistent_promoters$Name <- 
    factor(consistent_promoters$Name, 
           levels = 
             consistent_promoters$Name[order(
               consistent_promoters$Promoter_Strength,
               decreasing=TRUE)]) 
  
  ## Sort consistent_promoters_bc columns by promoter strength values. 
  consistent_promoters_bc <- 
    consistent_promoters_bc[order(consistent_promoters_bc$Promoter_Strength,
                                  decreasing=TRUE),]
  consistent_promoters_bc$Name <- 
    factor(consistent_promoters_bc$Name, 
           levels = 
             consistent_promoters_bc$Name[order(
               consistent_promoters_bc$Promoter_Strength,
               decreasing=TRUE)]) 
  
################################################################################
## Define promoters of interest to highlight in charts                        ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
  
## The "notable" variable created below is one which can be filled with specific
## promoters that one may want to highlight in charts for quick reference.
## By default it highlights a set of JEx promoters and the tac promoter.
  
  consistent_promoters <-  consistent_promoters %>%
    mutate(
      notable = case_when(
        Name == "JEa1" ~ "JEa1",
        Name == "tac" ~ "tac",
        Name == "JEa3" ~ "JEa3",
        Name == "JEa2" ~ "JEa2",
        Name == "JEc1" ~ "JEc1",
        Name == "JEc3" ~ "JEc3",
        Name == "JEc2" ~ "JEc2",
        Name == "JEb1" ~ "JEb1",
        Name == "JEb2" ~ "JEb2",
        T ~ "else"
      )
    )
  
  ## Dataframe containing the "noise type" for each promoter
consistency_vector <- subset(d_filter_bc_comp, select = c("Name","consistent"))
  
  
################################################################################
## Data preparation for generation of large data-dense graphics of promoter 
## data with either barcode variants shown or condition data shown 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
  
## Data frame that collects the average promoter strength (when collapsing all
## variables. Use for adding Overall PS to other dataframes, and this value
## is useful for arranging the order of samples on the x-axis in charts
d_overall_PS <- subset(d_filter_stats, select = c("Name", "Promoter_Strength"))
  colnames(d_overall_PS) <- c("Name","Overall_PS")
  
  ## Merge outlier status and overall PS variables with each of the three 
  ## following data frames. Use this for ordering the data prior to plotting.
  d_filter <- merge(d_overall_PS, d_filter, by="Name")
  d_filter_bc_stats <- merge(d_overall_PS, d_filter_bc_stats, by="Name")
  d_filter_con_stats <- merge(d_overall_PS, d_filter_con_stats, by="Name")
  
  ##potentially remove
d_filter_w_outlier_bc <- merge(d_overall_PS, d_filter_w_outlier_bc, by="Name")
d_filter_out_bc_stats <- merge(d_overall_PS, d_filter_out_bc_stats, by="Name")
d_filter_out_con_stats <- merge(d_overall_PS, d_filter_out_con_stats, by="Name")
  
## In the original code the original code the following occured in this section: 
## 'outlier_status' dataframe was merged with the dataframes below to allow
## the x-axis promoter name labels to be collored according to whether they
## had expression that was generally consistent across conditions.
## The correspodning variable is the Cond_sensitive variable which is in the
## UTR_data dataframe. We'll merge this here.
  d_filter <- merge(Cond_data, d_filter, by= "Name")
  d_filter_bc_stats <- merge(Cond_data, d_filter_bc_stats, by= "Name")
  d_filter_con_stats <- merge(Cond_data, d_filter_con_stats, by= "Name")
  
## The following three functions arrange the three dataframes such that 
## promoters are ordered according to 'consistent' and within the 'consistent'
## they are arranged according to overall PS
  d_filter_chart <-d_filter %>% 
    arrange(consistent, desc(Overall_PS)) %>% 
    mutate(Name=factor(Name, unique(Name)))
  
  d_filter_bc_stats <-d_filter_bc_stats %>% 
    arrange(consistent, desc(Overall_PS)) %>% 
    mutate(Name=factor(Name, unique(Name)))
  
  d_filter_con_stats <-d_filter_con_stats %>% 
    arrange(consistent, desc(Overall_PS)) %>% 
    mutate(Name=factor(Name, unique(Name)))
  
  ##potentially remove
  d_filter_out_chart <-d_filter_w_outlier_bc %>% 
    arrange(consistent, desc(Overall_PS)) %>% 
    mutate(Name=factor(Name, unique(Name)))
  
  d_filter_out_bc_stats <-d_filter_out_bc_stats %>% 
    arrange(consistent, desc(Overall_PS)) %>% 
    mutate(Name=factor(Name, unique(Name)))
  
  d_filter_out_con_stats <-d_filter_out_con_stats %>% 
    arrange(consistent, desc(Overall_PS)) %>% 
    mutate(Name=factor(Name, unique(Name)))
  
## data frame that collects the noise type values and can be used to generate 
## a vector with the colors of promtoers on the x-axis. Factored according to 
## the overall promoter strenths. The following stats summary function is a 
## lazy way to collapse all of the data points into one data point per promoter. 
  d_color_DF <- summarySE_median(d_filter_bc_stats, measurevar = "Overall_PS", 
                          groupvars = c("Name","Cond_sensitive"))
  
  ## vector with the color status if we want to color promoters according to 
  ## condition outlier status. If never_outlier, colored Orange. If condition 
  ## sensitive color orange,  otherwise color black.
  Cond_color <- ifelse(d_color_DF$Cond_sensitive == "No", "Black", "Orange")
  
  ## Add a new variable for plotting error bars for error values lower 
  ## than 1 that would potentially mess up the axes on log scale plots.
  d_filter_bc_stats <- d_filter_bc_stats %>%
    mutate(
      error_bar_ymin_se = case_when(
        Promoter_Strength-se > 1 ~ Promoter_Strength-se,
        Promoter_Strength-se <= 1 ~ 1.00001
      )
    )
  
  d_filter_con_stats <- d_filter_con_stats %>%
    mutate(
      error_bar_ymin_se = case_when(
        Promoter_Strength-se > 1 ~ Promoter_Strength-se,
        Promoter_Strength-se <= 1 ~ 1.00001
      )
    )
  d_filter_bc_stats <- d_filter_bc_stats %>%
    mutate(
      error_bar_ymin_sd = case_when(
        Promoter_Strength-sd > 1 ~ Promoter_Strength-sd,
        Promoter_Strength-sd <= 1 ~ 1.00001
      )
    )
  
  d_filter_con_stats <- d_filter_con_stats %>%
    mutate(
      error_bar_ymin_sd = case_when(
        Promoter_Strength-sd > 1 ~ Promoter_Strength-sd,
        Promoter_Strength-sd <= 1 ~ 1.00001
      )
    )
  
################################################################################
## Supplemental charts with data for all promoters, with promoters sorted      #
## by 5' UTR sensitivity (left set of samples are 5' UTR insensitive), and     #
## then sorted by promoter strength. Condition independence is indicated by    #
## the color of the promoter label on the bottom. If they are considered       #
## condition independent they are colored black, and if not they are colored   #
## orange.                                                                     #
## Some of the charts will display data prior to outlier barcode exclusion and #
## and others will display data post filtering. Another differentiating factor #
## is that some figures will show individual data points, and others standard  #
## error or standard deviation as error bars.                                  #
## Finally, mean promoter strength for each condition (mean of all barcodes)   #
## or for each barcode (mean of all condition samples) will be indicated as a  #
## primary data point.                                                         #
################################################################################
  
################################################################################
## START charts using data where outlier barcodes have been removed           ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

  
  ##### chart with individual data points and mean values - 
  ## separates promoter data by barcodes ####
  all_bc_by_barcode_points <- 
    ggplot(d_filter_chart, aes(x=Name, 
                               y=Promoter_Strength, 
                               color=as.factor(Replicate),
                               fill=as.factor(Replicate) )
    ) +
    scale_color_viridis(discrete = TRUE, 
                        option = "plasma"
    )+
    geom_point(position=position_dodge(width=0.5), 
               size= 0.4
    )+
    new_scale_color()+
    scale_fill_viridis(discrete = TRUE, 
                       option = "plasma")+
    geom_point(data=d_filter_bc_stats, 
               size=3,
               position=position_dodge(width=0.5),
               shape=21,
               aes(x=Name, 
                   y=Promoter_Strength, 
                   fill=as.factor(Replicate)), 
               color="black"
    )+
    theme_bw()+ 
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color="gray90",size=0.2),
      strip.background = element_rect(colour = "black", size=1),
      legend.position = "none",
      axis.text.x = element_text(angle = 90, color = Cond_color) 
    ) +
    scale_y_log10(  
      breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000,10000000),
      labels=trans_format('log10',math_format(10^.x)),
      limits=c(1,10000000), 
      expand = c(0,0)
    )  
  
  ##### chart with individual data points and mean values - 
  ## separates promoter data by conditions ####
  all_bc_by_condition_points <- 
    ggplot(d_filter_chart, aes(x=Name, 
                               y=Promoter_Strength, 
                               color=as.factor(Condition),
                               fill=as.factor(Condition) )
    ) +
    scale_color_viridis(discrete = TRUE, 
                        option = "mako",
                        begin = 0.1,
                        end= 0.8
    )+
    geom_point(position=position_dodge(width=0.5), 
               size= 0.4
    )+
    new_scale_color()+
    scale_fill_viridis(discrete = TRUE, 
                       option = "mako",
                       begin = 0.1,
                       end= 0.8
    )+
    geom_point(data=d_filter_con_stats, 
               size=3,
               position=position_dodge(width=0.5),
               shape=21,
               aes(x=Name, 
                   y=Promoter_Strength, 
                   fill=as.factor(Condition)), 
               color="black"
    )+
    theme_bw()+ 
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color="gray90",size=0.2),
      strip.background = element_rect(colour = "black", size=1),
      legend.position = "none",
      axis.text.x = element_text(angle = 90, color = Cond_color) 
    ) +
    scale_y_log10(  
      breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000,10000000),
      labels=trans_format('log10',math_format(10^.x)),
      limits=c(1,10000000), 
      expand = c(0,0)
    )  
  
  ### chart with mean values and error bars - 
  ## separates promoter data by barcodes ####
  all_bc_by_barcode_error_bars <- 
    ggplot(d_filter_bc_stats, aes(x=Name, 
                                  y=Promoter_Strength, 
                                  color=as.factor(Replicate), 
                                  fill=as.factor(Replicate))
    ) +
    geom_errorbar(size=0.75, 
                  position=position_dodge(width=0.5),
                  aes(ymin=error_bar_ymin_se, 
                      ymax=Promoter_Strength+se, 
                      color=as.factor(Replicate)),
                  width=0.1
    )+
    geom_point(size=3,
               position=position_dodge(width=0.5),
               shape=21,
               aes(x=Name, 
                   y=Promoter_Strength, 
                   fill=as.factor(Replicate)), 
               color="black"
    )+
    scale_color_viridis(discrete = TRUE, 
                        option = "plasma"
    )+
    scale_fill_viridis(discrete = TRUE, 
                       option = "plasma"
    )+
    theme_bw()+ 
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color="gray90",size=0.2),
      strip.background = element_rect(colour = "black", size=1),
      legend.position = "none",
      axis.text.x = element_text(angle = 90, color = Cond_color) 
    ) +
    scale_y_log10(  
      breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000,10000000),
      labels=trans_format('log10',math_format(10^.x)),
      limits=c(1,10000000), 
      expand = c(0,0)
    )  
  
  ##### chart with mean values and error bars - 
  ## separates promoter data by conditions ####
  all_bc_by_condition_error_bars <- 
    ggplot(d_filter_con_stats, aes(x=Name, 
                                   y=Promoter_Strength, 
                                   color=as.factor(Condition), 
                                   fill=as.factor(Condition))
    ) +
    geom_errorbar(size=0.75, 
                  position=position_dodge(width=0.5),
                  aes(ymin=error_bar_ymin_se, 
                      ymax=Promoter_Strength+se, 
                      color=as.factor(Condition)),
                  width=0.1
    )+
    geom_point(size=3,
               position=position_dodge(width=0.5),
               shape=21,
               aes(x=Name, 
                   y=Promoter_Strength, 
                   fill=as.factor(Condition)), 
               color="black"
    )+
    scale_color_viridis(discrete = TRUE, 
                        option = "mako",
                        begin = 0.1,
                        end = 0.8
    )+
    scale_fill_viridis(discrete = TRUE, 
                       option = "mako",
                       begin = 0.1,
                       end = 0.8
    )+
    theme_bw()+ 
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color="gray90",size=0.2),
      strip.background = element_rect(colour = "black", size=1),
      legend.position = "none",
      axis.text.x = element_text(angle = 90, color = Cond_color) 
    ) +
    scale_y_log10(  
      breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000,10000000),
      labels=trans_format('log10',math_format(10^.x)),
      limits=c(1,10000000), 
      expand = c(0,0)
    )  
  
######## print charts #######
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## END charts using data where outlier barcodes have been removed             ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
  
################################################################################
## START charts using data where outlier barcodes have been included          ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
  
##### chart with individual data points and mean values - 
## separates promoter data by barcodes ####
  all_bc_by_barcode_points_outliers <- 
    ggplot(d_filter_out_chart, aes(x=Name, 
                                   y=Promoter_Strength, 
                                   color=as.factor(Replicate),
                                   fill=as.factor(Replicate) )
    ) +
    scale_color_viridis(discrete = TRUE, 
                        option = "plasma"
    )+
    geom_point(position=position_dodge(width=0.5), 
               size= 0.4
    )+
    new_scale_color()+
    scale_fill_viridis(discrete = TRUE, 
                       option = "plasma")+
    geom_point(data=d_filter_out_bc_stats, 
               size=3,
               position=position_dodge(width=0.5),
               shape=21,
               aes(x=Name, 
                   y=Promoter_Strength, 
                   fill=as.factor(Replicate)), 
               color="black"
    )+
    theme_bw()+ 
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color="gray90",size=0.2),
      strip.background = element_rect(colour = "black", size=1),
      legend.position = "none",
      axis.text.x = element_text(angle = 90, 
                                 #color = Cond_color
      ) 
    ) +
    scale_y_log10(  
      breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000,10000000),
      labels=trans_format('log10',math_format(10^.x)),
      limits=c(1,10000000), 
      expand = c(0,0)
    )  
  
  ##### chart with individual data points and mean values - 
  ## separates promoter data by condition ####
  all_bc_by_condition_points_outliers <- 
    ggplot(d_filter_out_chart, aes(x=Name, 
                                   y=Promoter_Strength, 
                                   color=as.factor(Condition),
                                   fill=as.factor(Condition) )
    ) +
    scale_color_viridis(discrete = TRUE, 
                        option = "mako",
                        begin = 0.1,
                        end= 0.8
    )+
    geom_point(position=position_dodge(width=0.5), 
               size= 0.4
    )+
    new_scale_color()+
    scale_fill_viridis(discrete = TRUE, 
                       option = "mako",
                       begin = 0.1,
                       end= 0.8
    )+
    geom_point(data=d_filter_out_con_stats, 
               size=3,
               position=position_dodge(width=0.5),
               shape=21,
               aes(x=Name, 
                   y=Promoter_Strength, 
                   fill=as.factor(Condition)), 
               color="black"
    )+
    theme_bw()+ 
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color="gray90",size=0.2),
      strip.background = element_rect(colour = "black", size=1),
      legend.position = "none",
      axis.text.x = element_text(angle = 90, 
                                 #color = Cond_color
      ) 
    ) +
    scale_y_log10(  
      breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000,10000000),
      labels=trans_format('log10',math_format(10^.x)),
      limits=c(1,10000000), 
      expand = c(0,0)
    )  
  
  
  dotplot_file_name <- 
    paste("./",
          Output_file_label,
          "_graphics/",
          Output_file_label,
          "_Supplemental scatter_median.pdf",
          sep="")
  
  pdf(dotplot_file_name, width = 80, height = 7 )
  print(all_bc_by_barcode_points)
  print(all_bc_by_barcode_points_outliers)
  print(all_bc_by_condition_points)
  print(all_bc_by_condition_points_outliers)
  print(all_bc_by_barcode_error_bars)
  print(all_bc_by_condition_error_bars)
  dev.off()
  
  
################################################################################
##           Begin capture of pdf for charts                                   #
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ### ## ## ## ##
  
  
################################################################################
##  Generate bar charts with condition independent, low noise promoters        ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
  
  Promoter_overview <- ggplot(consistent_promoters,
                              aes(x = Name, y = Promoter_Strength)) +
    ## the identity value tells it to use the number's we're inputting rather
    ## than counting some sort of other values to generate the data being plotted.
    ## plots the data as a bar chart. Here the fill is colored accoording to the
    ## name and the data comes directly from the dataframe (e.g. no calculations
    ## done to data for the charting)
    geom_bar(stat = "identity", width = 0.8, aes(fill = Name)) +
    
    # Repeating gradient color fill according to the x-axis location (i.e. Name)
    scale_fill_viridis(discrete = TRUE,
                       option = "plasma",
                       begin = 0.22) +
    #### uses the ggnewscale package function new_scale() to essentially reset 
    #### the aesthetic/scale indicated, fill here, so it can be reapplied in a 
    #### fresh chart overlay
    #### https://eliocamp.github.io/codigo-r/2018/09/multiple-color-and-fill-scales-with-ggplot2/
    
    new_scale("fill") +
    
    # generates error bars with the standard deviation (sd) in the promoter strength as the error bars
    geom_errorbar(
      aes(ymin = Promoter_Strength - sd, ymax = Promoter_Strength + sd),
      size = 0.3,
      width = .1
    ) +
    
    # Add the label promoter to the x-axis
    xlab("Promoters") +
    theme(
      axis.line = element_line(size = 0.5, colour = "black"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(size = 0.1, colour = "black"),
      #panel.grid.minor.y = element_line(size = 0.1, colour = "black"),
      panel.background = element_rect(fill = "white"),
      panel.border = element_rect(
        size = 0.5,
        colour = "black",
        fill = NA
      ),
      # removes the legend entirely
      legend.position = "none"
    ) +
    
    # Edit parameters of the y-axis
    scale_y_log10(
      breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000, 10000000),
      labels = trans_format('log10', math_format(10 ^ .x)),
      limits = c(1, 10000000),
      expand = c(0, 0)
    ) +
    
    ggtitle("consistent promoters, standard deivation bars")
  
  ################################################################################
  Promoter_overview_bc <- ggplot(consistent_promoters_bc,
                                 aes(x = Name, y = Promoter_Strength)
  ) +
    ## the identity value tells it to use the number's we're inputting rather
    ## than counting some sort of other values to generate the data being plotted.
    ## plots the data as a bar chart. Here the fill is colored accoording to the
    ## name and the data comes directly from the dataframe (e.g. no calculations
    ## done to data for the charting)
    geom_bar(stat = "identity", width = 0.8, aes(fill = Name)) +
    
    # Repeating gradient color fill according to the x-axis location (i.e. Name)
    scale_fill_viridis(discrete = TRUE,
                       option = "plasma",
                       begin = 0.22) +
    #### uses the ggnewscale package function new_scale() to essentially 
    #### reset the aesthetic/scale indicated, fill here, so it can be reapplied 
    #### in a fresh chart overlay
    #### https://eliocamp.github.io/codigo-r/2018/09/multiple-color-and-fill-scales-with-ggplot2/
    
    new_scale("fill") +
    
    #### recoloring the chart with noise as the variable determining shading.
    
    geom_bar(
      stat = "identity",
      aes(fill = consistent, alpha = consistent),
      data = consistent_promoters_bc,
      width = 0.87
    ) +
    
    scale_fill_viridis(discrete = TRUE, option = "inferno") +
    
    scale_alpha_manual("test", values = c(0, 1, 1, 1)) +
    
    # generates error bars with the standard deviation (sd) in the promoter strength as the error bars
    geom_errorbar(
      aes(ymin = Promoter_Strength - sd, ymax = Promoter_Strength + sd),
      size = 0.3,
      width = .1
    ) +
    
    # Add the label promoter to the x-axis
    xlab("Promoters") +
    theme(
      axis.line = element_line(size = 0.5, colour = "black"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(size = 0.1, colour = "black"),
      #panel.grid.minor.y = element_line(size = 0.1, colour = "black"),
      panel.background = element_rect(fill = "white"),
      panel.border = element_rect(
        size = 0.5,
        colour = "black",
        fill = NA
      ),
      # removes the legend entirely
      legend.position = "none"
    ) +
    
    # Edit parameters of the y-axis
    scale_y_log10(
      breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000, 10000000),
      labels = trans_format('log10', math_format(10 ^ .x)),
      limits = c(1, 10000000),
      expand = c(0, 0)
    ) +
    
    ggtitle("consistent promoters_bc, standard deviation bars")
  
  
  ## 'pdf(' Activates capture of a pdf file from the graph described below.
  ## Need to cap the graph with dev.off() to end capture.
  bar_and_chart_file_name <- 
    paste("./",
          Output_file_label,
          "_graphics/",
          Output_file_label,
          "_bar-and-cond-charts_median.pdf",
          sep="")
  pdf(
    bar_and_chart_file_name,
    width = 6,
    height = 3
  )
  
  ## Print plots
  print(Promoter_overview)
  print(Promoter_overview_bc)
  # dev.off() turns off the pdf recording.
  dev.off()
  
  ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
  ## Generate a scatter plot comparing performance of individual promoters across
  ## conditions. Dots are mean Promoter_Strength for a given condition. Dots are
  ## colored according to condition. A regression line will be generated for each
  ## condition to highlight both the general trend in promoter strengths across
  ## conditions. The steepness of the regression line slope will also serve as an
  ## indicator of the dynamic range of promoter performance in a given condition.
  ## Generally, a shallow slope (sub 1) means there is less difference in promoter
  ## strength from strongest to weakest promoter for the x-axis condition than for
  ## the y-axis condition. The converse is true for a >1 slope. If all slopes are
  ## close to 1, it will suggest that the relative performance of the promoters
  ## in the library are pretty consitutive and insensitive to condition. If the
  ## slopes are close to 1 (suggesting most promoters are generally condition
  ## insensitive, an alternative method for calling a promoter condition 
  ## insensitive could be some fractional difference between promoter performance
  ## and predicted promoter performance (somehow based upon slope?)
  ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
  
  
  
  ## Generate an array containing the distinct conditions used in the experiment.
  ## This will be used to filter a data frame into subsets of the data containing
  ## only the data for that condition
  condition_array <- unique(subset(d_filter_con_stats, 
                                   select = c("Condition")
  )
  )
  ## Generate a subset of the d_filter_con_stats that pulls out the first
  ## condition in the array. Doing this allows us to set it as the X-axis in 
  ## promoter comparison charts
  d_cond_1 <- subset(d_filter_con_stats, 
                     Condition == condition_array$Condition[1], 
                     select = c("Name","Condition","Promoter_Strength","se"))
  
  ## These generate a distinct label for the promoteer strength and standard error
  ## for the condition excised in the previous line of code. Together these three
  ## lines allow us to keep the data for this condition as a distinct variable
  ## when we recombine all of the condition data into a single reformatted
  ## dataframe below
  condition_label = paste(condition_array$Condition[1], "_PS", sep="")
  
  condition_se_label = paste(condition_array$Condition[1], "_se", sep="")
  
  colnames(d_cond_1)<-c("Name","Condition", condition_label, condition_se_label)
  
  ## Generate a subset of the d_filter_con_stats that pulls out the all condition
  ## condition in the array that is not in the first condition in the array. 
  d_cond_others <- subset(d_filter_con_stats, 
                          Condition != condition_array$Condition[1], 
                          select=c("Name","Condition","Promoter_Strength","se"))
  
  ## These generate a distinct label for the promoteer strength and standard error
  ## for the condition excised in the previous line of code. Together these three
  ## lines allow us to keep the data for this condition as a distinct variable
  ## when we recombine all of the condition data into a single reformatted
  ## dataframe below  
  condition_label2 = paste("not_",
                           condition_array$Condition[1],
                           "_PS",
                           sep="")
  
  condition_se_label2 = paste("not_",
                              condition_array$Condition[1],
                              "_se",
                              sep="")
  
  colnames(d_cond_others)<-c("Name",
                             "Condition", 
                             condition_label2, 
                             condition_se_label2)
  
  ## merge the dataframes generated above. Merging these two allows us to use the
  ## single dataframe for charting promoter strengths for one condition on x-axis
  ## and promoter strengths for the other conditions on the y-axis with different
  ## colored markers for each condition
  d_cond_comp <- merge(d_cond_1, d_cond_others, by = "Name")
  
  ################################################################################
  ## Plot chart conmparing the mean promoteer strengths for each promoter in a 
  ## given condition against its promoter strength under different conditionos
  ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
  
  ## the aes() call here uses a reference to a position in the array, as this data
  ## should be in the same position for each dataset that gets run through the
  ## function. Howeever, the columns name will be distinct each time, as it is 
  ## dependent on the conditions being tested. So we cannot simply reference the 
  ## name of the coluumn.
  ## Color and fill will be determined by the Condition variable.
  Cond_compare_chart <- ggplot(d_cond_comp,aes(x=d_cond_comp[,3], 
                                               y=d_cond_comp[,6]), 
                               color=as.factor(Condition),
                               fill=as.factor(Condition),
  ) +
    
    
    
    # add linear regression line. Color according to Condition.y variable.
    geom_smooth(method="lm", se=F, aes(color = Condition.y)) + 
    # use the geom_point() to add points to graph
    geom_point(aes(color = Condition.y, fill = Condition.y), shape = 21,
               size = 2, stroke=.4) +
    # determines the color of the fill for the datapoints
    scale_fill_viridis(discrete = TRUE, 
                       option = "plasma",
                       begin = 0.1,
                       end= 0.8,
                       alpha = 0.3
    )+
    # determines the color of both the outlines of the data points and the
    # regression line.
    scale_color_viridis(discrete = TRUE, 
                        option = "plasma",
                        begin = 0.1,
                        end= 0.8,
    )+
    
    ## The two additions below (new_scale_color and scale_color_viridis) are likely
    ## no longer needed, and will probably need to / should be deprecated.
    new_scale_color()+
    scale_color_viridis(discrete = TRUE, 
                        option = "plasma",
                        begin = 0.1,
                        end= 0.8,
    )+
    
    # Edit parameters of the y-axis
    scale_y_log10(  
      breaks = c(10, 100, 1000, 10000, 100000, 1000000,10000000),
      labels=trans_format('log10',math_format(10^.x)),
      limits=c(10,10000000), 
      expand = c(0,0)
    ) +
    ylab("other condition promoter activity") +
    # Edit parameters of the x-axis
    scale_x_log10(  
      breaks = c(10, 100, 1000, 10000, 100000, 1000000,10000000),
      labels=trans_format('log10',math_format(10^.x)),
      limits=c(10,10000000), 
      expand = c(0,0)
    ) +
    xlab(paste(condition_array$Condition[1]," condition promoter activity")) +
    # Defines graphical elements of the plot.  
    theme(
      
      panel.grid.major.y=element_line(color="gray20", size = 0.1),
      panel.grid.minor.y=element_blank(),
      panel.grid.major.x=element_line(color="gray20", size = 0.1),
      panel.grid.minor.x=element_blank(),
      panel.background = element_rect(fill = "white"),
      panel.border = element_rect(fill=NA, color = "black" ),
      axis.ticks = element_line(size=0.25)
      
    )
  
  # force plot to be square and not reformat shape as a function of data values
  # (note: if the scales for x and y-axis values were different, 
  # we may have to fiddle with this.)
  Cond_compare_chart <- Cond_compare_chart + coord_fixed()
  
  # Print the chart as a pdf
  bar_and_chart_file_name <- 
    paste("./",
          Output_file_label,
          "_graphics/",
          Output_file_label,
          "_cond_compare_chart_median.pdf",
          sep="")
  pdf(
    bar_and_chart_file_name,
    width = 5,
    height = 4
  )
  
  print(Cond_compare_chart)
  dev.off()
  # end generation of and printing of chart.
  
  
  #############################################################################
  ## Output data / stats for generating a table containing information on    ##
  ## promoter stats. Basically, number of barcodes in dataset, etc.          ##
  ## Generating a dataframe to hold the data.
  ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
  column_names <- c("Organism",
                    "num_bc_pass",
                    "num_prom_in_bc_pass",
                    "num_prom_suff_data",
                    "num_bc_pass_outlier",
                    "num_UTR_insensitive",
                    "num_Cond_insensitive",
                    "num_consistent",
                    "UTR_cutoff",
                    "Cond_cutoff",
                    "expression_range",
                    "Problematic_percent",
                    "Number_of_promoters_w_prob_bc"
  )
  df_output <-data.frame(matrix(ncol = 13, nrow = 1))
  colnames(df_output) <- column_names
  
  ## assigned Organism as input_file_label. Will replace this based on a lookup
  ## table
  df_output$Organism <- Input_file_label
  
  ## count number of barcodes that pass the minimal amount of data filter
  ## Many of these are removed due to lack of RNA reads, as we did not include 0
  ## reads for analysis, as imputing an arbitrary non-zero value for analysis 
  ## added a lot of noisy values for weak promoters
  df_output$num_bc_pass <- nrow(filtered_barcodes)
  ## count number of promoters that pass the minimal amount of data filter
  df_output$num_prom_in_bc_pass <- n_distinct (d_barcode_filter$Name)
  ## count number of promoters with sufficient data to continue analysis
  df_output$num_prom_suff_data <- nrow(filtered_promoters)
  ## number of barcodes that are not considered outliers
  df_output$num_bc_pass_outlier <- nrow(d_barcode_stats_no_outliers)
  ## number of promoters that are considered insensitive to UTR sequences
  df_output$num_UTR_insensitive <- sum(
    with(d_filter_bc_comp, 
         d_filter_bc_comp$UTR_sensitive == "No")
  )
  ## number of promoters that are considered insensitive to Condition
  df_output$num_Cond_insensitive <- sum(
    with(d_filter_bc_comp, 
         d_filter_bc_comp$Cond_sensitive == "No")
  )
  ## number of promoters that are considered both UTR insensitive and 
  ## Condition insensitive
  df_output$num_consistent <- nrow(consistent_promoters)
  ## input the user defined or default (if used) UTR and Cond cutoff values
  df_output$UTR_cutoff <- UTR_cutoff
  df_output$Cond_cutoff <- Cond_cutoff
  ## export range of expression values (mean) for consistent promoters
  df_output$expression_range <- (
    max(consistent_promoters$Promoter_Strength) / 
      min(consistent_promoters$Promoter_Strength)
  )
  ##* % total number of barcodes (in the assessment) that are considered
  ##* "problematic" barcodes: those whose expression data would have given a value
  ##* greater than (UTR_cutoff/2) fold different than the mean promoter expression.
  ##* Problematic_percent includes all non-outlier barcodes when calculating median
  ##* expression, and then evaluates all barcodes including outlier barcodes when
  ## performing assessment, as they would have been included in a 'single barcode' 
  ## based method.
  
  ##* calculate % of promoters that are problematic
  bc_no_filter_prob <- sum(
    with(d_problematic, 
         d_problematic$problematic_bc == "Yes")
  )
  percent_problematic_nf <- 
    round(bc_no_filter_prob / nrow(d_problematic) * 100, digits = 2)
  
  df_output$Problematic_percent <- paste(percent_problematic_nf,"%",sep="")
  
  ## determine number of promoters that have a problematic barcode.
  ## Filtering out barcodes that are not problematic
  prob_promoter <- d_problematic %>% filter(d_problematic$problematic_bc == "Yes")
  
  ## collapsing the dataset to only include a single row per "Name" value.
  prob_promoter <- unique(subset(prob_promoter, 
                                 select = c("Name"))
  )
  
  ## Determine number of problematic promoters by counting number of rows in the
  ## dataset generated above.
  num_prob_promoter <- (nrow(prob_promoter) / 1)
  
  df_output$Number_of_promoters_w_prob_bc <- num_prob_promoter

  return(df_output)
}

################################################################################
##  Promoter data needed for cross-species comparisons
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
##  
##  Candidate dataframe as starting point -> d_filter_out_chart. Can collapse
##  this data as needed to look at mean values for barcodes, conditions, etc.
##  while maintaining all of the data points for big charts.
##
##    1. Organism
##        - add as a column to whatever dataframe is generated. Use function
##          input variable input File lable as the value. If more detailed name
##          is needed, it will be added elsewhere using a lookup table
##    2. Promoter name
##        - Use the "Name" column from the most appropriate source data frame(s)
##          for the data below.
##    3. Barcode replicate #
##        - Use the "Replicate" column from the most appropriate source data
##          frame(s) for the data below.
##    4. Barcode strength
##        - Use the "Promoter Strength " column from the source data frame(s) 
##    5. Barcode outlier status
##        - Use the "bc_outlier" column from the source data frame(s) 
##    6. Promoter and barcode presence in data set (is it filtered out prior to 
##       barcode outlier analysis?)
##        - This may be most readily obtained from Rob's scripts if we want to
##          include RNA = 0 samples
##    7. Consistent variable value
##        - Use the "consistent" column from the source data frame(s) 
##    8. Condition-based max_min_ratio
##        - Take the data that is contained for the "Cond_data" dataframe and 
##          transfer it to a different dataframe prior to subsetting the the 
##          Cond_data to only containing Name and Cond_sensitive
##    9. Barcode-based max_min_ratio
##        - Take the data that is contained for the "UTR_data" dataframe and 
##          transfer it to a different dataframe prior to subsetting the the 
##          UTR_data to only containing Name and UTR_sensitive
##    10. Pseudomonad status
##        - Compare condition and barcode independent status among all three 
##          Pseudomonads. Look just needs to be done post-merging of the output
##          from this script.
##    11. Proteobacteria status
##        - See above but include Rpal.
##    12. Promoter present in organism sample? (two versions? See below)
##      a. Are promoter DNA barcodes present at >5 counts in >75% of samples? In 100% of samples?
##        - Need data from Rob's script for this
##      b. Did the promoter have sufficient data to be analyzed in a given organism?
##        - These may require a distinct output from the other data included in
##          this list. We would need the DNA counts from Rob's script, and then
##          for those that do pass the DNA count filter, whether they have 
##          sufficient samples with non-zero RNA reads for a given barcode to be
##          used. Then determine if there are at least 3 barcodes for a given
##          promoter that passes all of these filters. This may be best output
##          from a modified version of Rob's scripts.
##    13. If promoter is present, was its expression high enough to be assessed?
##        - This is basically lumped in with the broader question in 12b. See
##          above.
##    14. All data points for each organism (post barcode outlier filtering)
##        - Add the following to d_filter_out_chart and then output:
##            a. Organism label (see #1 above)
##            b. UTR_sensitive
##            c. Promoter_Strength from d_filter_bc_stats -> PS_bc_no_outlier_bc
##            d. Cond_sensitive
##            c. Promoter_Strength from d_filter_con_stats -> 
##               PS_con_no_outlier_bc




















