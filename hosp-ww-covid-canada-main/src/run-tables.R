## Generate tables for LaTeX document ##
library(tidyverse)
library(patchwork)
library(kableExtra)

source('utils.R')

# Retrieve start and end dates
dates = readRDS(file = '../out/var-analysis.rds')[['df.dates']] |>
  rename(City = "city",
         Province = "prov",
         Variant = "variant",
         Start = "start",
         End = "end")

tbl.dates = knitr::kable(dates, 
                         format = 'latex',
                         vline = '',
                         linesep = ''
                         )
writeLines(tbl.dates, '../ms/tables/dates-tbl.tex')

# Summary table
df.summary = readRDS('../out/mcmcobj.rds')[["postsum"]]

tbl.summary = create_sum_table(df.summary)|>
  knitr::kable(format = 'latex', 
               align = c('l','c','c'),
               linesep = '', 
               vline= '',
               escape = FALSE)
writeLines(tbl.summary, '../ms/tables/summary-tbl.tex')
