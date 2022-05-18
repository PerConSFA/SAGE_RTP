###############################################################################
# Code written by Joshua Elmore (joshua.elmore@pnnl.gov) for SAGE: Serine     #
# recombinase-assisted genome engineering manuscript                          #
# Note, this version of the promoter analysis worksheet is for processing     #
# non-log scale values.                                                       #
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
## added on 1-12-2022 to plot lin reg eq on chart
library(ggpmisc)



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


generate_org_compare_charts_1 <- function(source_data_frame, 
                                          organismX = "X_test", 
                                          organismY ="Y_test",
                                          array_of_organisms){  
################################################################################
## Plot chart conmparing the mean promoteer strengths for each promoter in a 
## given organism against its promoter strength under different organism
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

## the aes() call here uses a reference to a position in the array, as this data
## should be in the same position for each dataset that gets run through the
## function. Howeever, the columns name will be distinct each time, as it is 
## dependent on the organisms being tested. So we cannot simply reference the 
## name of the column.
## Color and fill will be determined by the organism variable.

  d_org_comp = source_data_frame
  
  #should remove duplicates
  
  d_org_comp <- unique(d_org_comp)
  
  
  orgX = organismX
  orgY = organismY
  ##this bit is used for the regression line printing
  my.formula <- y ~ x
  ## end bit
  
  ## Use 3 and 6 to compare Overall_PS rather than Overall_PS_BC
  org_compare_chart <- ggplot(d_org_comp,aes(x=d_org_comp[,4], 
                                             y=d_org_comp[,7]), 
                              color=as.factor(organism),
                              fill=as.factor(organism),
  ) +
    
    
    # add linear regression line. Color according to organism.y variable.
    geom_smooth(method="lm", se=F, aes(color = organism.y),
                ##this bit is used for the regression line printing
                formula = my.formula
                ) + 
    
    # use the geom_point() to add points to graph
    geom_point(aes(color = organism.y, fill = organism.y), shape = 21,
               size = 2, stroke=.4) +
    
    ##prints the regression line
    stat_poly_eq(formula = my.formula, 
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                 parse = TRUE) +    
    
    # determines the color of the fill for the datapoints
    scale_fill_viridis(discrete = TRUE, 
                       option = "plasma",
                       begin = 0.5,
                       end= 0.5,
                       alpha = 0.3
    )+
    # determines the color of both the outlines of the data points and the
    # regression line.
    scale_color_viridis(discrete = TRUE, 
                        option = "plasma",
                        begin = 0.1,
                        end= 0.8,
    )+
    
 
    
    # Edit parameters of the y-axis
    scale_y_log10(  
      breaks = c(100, 1000, 10000, 100000, 1000000,10000000),
      labels=trans_format('log10',math_format(10^.x)),
      limits=c(100,10000000), 
      expand = c(0,0)
    ) +
    ylab(paste(orgY," promoter activity")) +
    # Edit parameters of the x-axis
    scale_x_log10(  
      breaks = c(100, 1000, 10000, 100000, 1000000,10000000),
      labels=trans_format('log10',math_format(10^.x)),
      limits=c(100,10000000), 
      expand = c(0,0)
    ) +
    xlab(paste(orgX," promoter activity")) +
    
    # Defines graphical elements of the plot.  
    theme(
      
      panel.grid.major.y=element_line(color="gray20", size = 0.1),
      panel.grid.minor.y=element_blank(),
      panel.grid.major.x=element_line(color="gray20", size = 0.1),
      panel.grid.minor.x=element_blank(),
      panel.background = element_rect(fill = "white"),
      panel.border = element_rect(fill=NA, color = "black" ),
      axis.ticks = element_line(size=0.25),
      ## removes the legend, which is just a distraction here.
      legend.position = "none"
      
    )
  
  # force plot to be square and not reformat shape as a function of data values
  # note: if the scales for x and y-axis values were different, 
  # we may have to fiddle with this.
  org_compare_chart <- org_compare_chart + coord_fixed()
  
 
  
  # Print the chart as a pdf
  chart_file_name <- 
    paste("./Multi_org_graphics/",orgX, " vs ", orgY, " comp chart.pdf",
          sep="")
  pdf(
    chart_file_name,
    width = 5,
    height = 4
  )
  
  #log10 conversion prior to statistical analysis (linear regression)
  
  
  
  #linear regression stats summary
  csv_file_name <-
    paste("./Multi_org_graphics/",orgX, " vs ", orgY, " comp chart_stats.csv",
          sep="")
  
  ## Use 3 and 6 to compare Overall_PS rather than Overall_PS_BC
  
  x_test = log(d_org_comp[,4],10)
  y_test = log(d_org_comp[,7],10)
  
  
  sink(csv_file_name)
  print(summary(lm(y_test ~ x_test, data=d_org_comp)))
  sink()  # returns output to the console
  
  ## write summary of input data

  
   csv_file_name2 <-
    paste("./Multi_org_graphics/",orgX, " vs ", orgY, " comp chart_values.csv",
          sep="")
  
  write_csv(d_org_comp,
            csv_file_name2)
  
  print(org_compare_chart)
  dev.off()
  # end generation of and printing of chart.
  return(org_compare_chart)
}