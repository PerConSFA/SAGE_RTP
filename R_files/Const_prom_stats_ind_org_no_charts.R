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
library(devtools)
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


## for ANOVA testing
## library(broom)
## library(AICcmodavg)

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


Const_prom_stats_ind_org_no_charts <- function(Input_file_label = "Pflu",
                                           Output_file_label = "Pflu",
                                           num_cond = 3,
                                           num_bc_rep = 5,
                                           num_exp_rep = 4,
                                           UTR_cut = 3,
                                           Cond_cut = 4
                                           ){

  #######################################
  ## Test values to validate function  ##
  ## ## ## ## ## ## ## ## ## ## ## ## ###
  ## Input_file_label = "RHA1"
  ## Output_file_label = "Rjostii"
  ## ## ## ## ## ## ## ## ## ## ## ## ###

  #######################################
  ##  Set function variable values  #       
  ## ## ## ## ## ## ## ## ## ## ## ## ###
  Input_file_label = "Pflu"
  Output_file_label = "Pflu"
  num_cond = 4
  num_bc_rep = 5
  num_exp_rep = 4
  UTR_cut = 3.5
  Cond_cut = 4
  
   
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
d_barcode_filter <- summarySE(d, measurevar = "Promoter_Strength", 
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
d_stats <- summarySE(d, measurevar="Promoter_Strength", 
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
## Also renamed the N (total counts of data points for a promoter) to N_promoter_samples 
## so that this data is added to d_barcode_stats
filtered_promoters <- subset(d_stats, select =c("Name","N"))
colnames(filtered_promoters)<- c("Name","N_promoter_samples")

## Generate statistics for the promoter strength values with samples 
## grouped by the barcode variable
d_barcode_stats <- summarySE(d, measurevar = "Promoter_Strength", 
                             groupvars = "Barcode",na.rm=TRUE)

## Generate a dataframe with the names added to the the barcode stats dataframe
d_barcode_stats <- merge(d_barcode_stats,BC_lookUp,by ="Barcode")

## Removes promoters or barcodes lacking a value from each datafram
d_stats <- na.omit(d_stats)
d_barcode_stats <- na.omit(d_barcode_stats)

################################################################################
## Removes promoters from d_barcode_stats that failed to pass the N > 50% of
## max promoter data point threshold. Also adds the #N for the promoter as
## the variable N_promoter_samples.
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
## New section added on 6-17-21 that aims to remove barcode outliers using the 
## 1.5x IQR threshold method to determine which barcodes are potential outliers. 
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
d_filter_bc_stats <- summarySE(d_filter, measurevar="Promoter_Strength", 
                               groupvars = c("Replicate","Name"), na.rm=TRUE)
d_filter_con_stats <- summarySE(d_filter, measurevar="Promoter_Strength",
                                groupvars = c("Condition","Name"), na.rm=TRUE)
d_filter_stats <- summarySE(d_filter, measurevar="Promoter_Strength", 
                            groupvars = "Name", na.rm=TRUE)

## These few lines of codes generate stats as above with d_filter, but with
## d_filter_w_outlier_bc.
d_filter_out_bc_stats <- summarySE(d_filter_w_outlier_bc, 
                                   measurevar="Promoter_Strength", 
                                   groupvars = c("Replicate","Name"), 
                                   na.rm=TRUE)
d_filter_out_con_stats <- summarySE(d_filter_w_outlier_bc, 
                                    measurevar="Promoter_Strength",
                                    groupvars = c("Condition","Name"), 
                                    na.rm=TRUE)
d_filter_out_stats <- summarySE(d_filter_w_outlier_bc, 
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

##* test code added on  7/18/21 to determine whether a barcode is problematic -
##* in other words, if it were the sole barcode would its inclusion have
##* substantially changed the reported promoter strength that we determined
##* by using multiple barcodes per promoter. The cutoff we will use is the 
##* UTR_cutoff value - this way we can also do a parameter sweep if desired.

##* d_sufficient is a dataframe containing promoters that have at least 3 
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
names(d_problematic_means)[names(d_problematic_means) == "Promoter_Strength"] <- "Overall_PS"
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
d_barcode_comp <- summarySE(d_barcode_stats_no_outliers, 
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






################################################################################
# Calculate promoter condition and 5'UTR independence statuses                 #
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## # #

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
UTR_data <- summarySE(d_filter_bc_stats, measurevar = "Promoter_Strength", 
                      groupvars = c("Name"), na.rm = TRUE)
UTR_data <- subset(UTR_data, select =c("Name","Promoter_Strength"))
UTR_data$max_RTA <- by(d_filter_bc_stats$Promoter_Strength, 
                       d_filter_bc_stats$Name, max)
UTR_data$min_RTA <- by(d_filter_bc_stats$Promoter_Strength, 
                       d_filter_bc_stats$Name, min)
UTR_data$max_min_ratio <- UTR_data$max_RTA / UTR_data$min_RTA

## Cutoff fold change between maximum RTA and minimum RTA barcode RTAs that for
## considering a promoter to have an insensitivie 5'UTR


UTR_data <- UTR_data %>%
  mutate(
    UTR_sensitive = case_when(
      UTR_data$max_min_ratio >= UTR_cutoff
      ~ 'Yes',
      T ~ 'No'
    )
  )
## Subset UTR_data 
UTR_data <- subset(UTR_data, select =c("Name","UTR_sensitive"))
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##



################################################################################
##  screening for condition independent expression
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
Cond_data <- summarySE(d_filter_con_stats, measurevar = "Promoter_Strength", 
                       groupvars = c("Name"), na.rm = TRUE)
Cond_data <- subset(Cond_data, select =c("Name","Promoter_Strength"))
Cond_data$max_RTA <- by(d_filter_con_stats$Promoter_Strength, 
                        d_filter_con_stats$Name, max)
Cond_data$min_RTA <- by(d_filter_con_stats$Promoter_Strength, 
                        d_filter_con_stats$Name, min)
Cond_data$max_min_ratio <- Cond_data$max_RTA / Cond_data$min_RTA
## Cutoff fold change between maximum RTA and minimum RTA barcode RTAs that for
## considering a promoter to perform independent of condition



Cond_data <- Cond_data %>%
  mutate(
    Cond_sensitive = case_when(
      Cond_data$max_min_ratio >= Cond_cutoff
      ~ 'Yes',
      T ~ 'No'
    )
  )
## Subset Cond_data
Cond_data <- subset(Cond_data, select =c("Name","Cond_sensitive"))
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##


## Takes barcode values and generates stats for the promoter_strength for each 
## 'Promoter' at the 'Name' after first generating stats within barcodes.
##
## This is distinct from similar data structues in which all replicates 
## for all barcodes are lumped together in the statistical analyses.
d_filter_bc_comp <- summarySE(d_filter_bc_stats, 
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
d_color_DF <- summarySE(d_filter_bc_stats, measurevar = "Overall_PS", 
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
##* Problematic_percent includes all non-outlier barcodes when calculating mean
##* expression, and then evaluates all barcodes )including outlier barcodes when
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
