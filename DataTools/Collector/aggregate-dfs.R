###########################################################################
###                                                                     ###
###          DATA SCRIPT TO COMBINE ALL DATA INTO ONE DATAFRAME         ###
###                                                                     ###
###########################################################################

# libraries we are using
rm(list = ls())
library(ggplot2)
library(cowplot)
library(dplyr)
library(PupillometryR)

# run setup script
setwd("C:/Users/jgh9094/Desktop/Research/Projects/SelectionDiagnostics/Diagnosing-Island-Structures/")

p_theme <- theme(
  text = element_text(size = 28),
  plot.title = element_text( face = "bold", size = 22),
  panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  legend.title=element_text(size=18),
  legend.text=element_text(size=14),
  axis.title = element_text(size=18),
  axis.text = element_text(size=16),
  legend.position="bottom",
  panel.background = element_rect(fill = "#f1f2f5",
                                  colour = "white",
                                  size = 0.5, linetype = "solid")
)


# default variables

REPLICATES = 100
MODEL = c('EA','IS','NMIS')
SCHEME = c('TRUNCATION','TOURNAMENT','LEXICASE')
DIAGNOSTIC = c('EXPLOITATION_RATE', 'ORDERED_EXPLOITATION', 'CONTRADICTORY_OBJECTIVES', 'MULTIPATH_EXPLORATION')
DIMENSIONALITY = 100
cb_palette <- c('#D81B60','#1E88E5','#FFC107')
SHAPE = c(15,16,17)
TSIZE = 20
GENERATIONS = 50000

# data related
DATA_DIR = './ISLAND-STRUCTURE-DATA-1-20-22/'

# go through each diagnostic and collect over time data for cross comparison (cc)
over_time_df = data.frame()
for(model in MODEL)
{
  print(model)
  for(scheme in SCHEME)
  {
    dir = paste(DATA_DIR,model,'/over-time-',scheme, '.csv', sep = "", collapse = NULL)
    print(paste('DIRECTORY:',dir, sep = "", collapse = NULL))
    over_time_df = rbind(over_time_df, read.csv(dir, header = TRUE, stringsAsFactors = FALSE))
  }
}


# go through each diagnostic and collect best over time for cross comparison (cc)
best_df = data.frame()
for(model in MODEL)
{
  print(model)
  for(scheme in SCHEME)
  {
    dir = paste(DATA_DIR,model,'/best-',scheme, '.csv', sep = "", collapse = NULL)
    print(paste('DIRECTORY:',dir, sep = "", collapse = NULL))
    best_df = rbind(best_df, read.csv(dir, header = TRUE, stringsAsFactors = FALSE))
  }
}

# get generation a satisfactory solution is found for cross comparison (cc)
sati_sol_df = data.frame()
for(model in MODEL)
{
  print(model)
  for(scheme in SCHEME)
  {
    dir = paste(DATA_DIR,model,'/ssf-',scheme, '.csv', sep = "", collapse = NULL)
    print(paste('DIRECTORY:',dir, sep = "", collapse = NULL))
    sati_sol_df = rbind(sati_sol_df, read.csv(dir, header = TRUE, stringsAsFactors = FALSE))
  }
}

colnames(over_time_df)[colnames(over_time_df) == "Selection.Scheme"] = 'Selection\nScheme'
over_time_df$Structure <- factor(over_time_df$Structure, levels = MODEL)
over_time_df$sel_pre = over_time_df$sel_pre * -1.0

colnames(best_df)[colnames(best_df) == "SEL"] = 'Selection\nScheme'
colnames(best_df)[colnames(best_df) == "MOD"] = 'Structure'
best_df$Structure <- factor(best_df$Structure, levels = MODEL)

colnames(sati_sol_df)[colnames(sati_sol_df) == "SEL"] = 'Selection\nScheme'
colnames(sati_sol_df)[colnames(sati_sol_df) == "GEN"] = 'Generations'
colnames(sati_sol_df)[colnames(sati_sol_df) == "MOD"] = 'Structure'
sati_sol_df$Structure <- factor(sati_sol_df$Structure, levels = MODEL)
