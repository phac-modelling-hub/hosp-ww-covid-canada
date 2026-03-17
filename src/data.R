# ---- Retrieve and prep data for model ---- 

library(tidyr)
library(dplyr)
library(ggplot2) ; theme_set(theme_bw())
library(lubridate)
library(patchwork)
library(CCT) # remotes::install_github("TheZetner/cct")

source("utils.R")
source("var.R")
source("plot.R")

# Wastewater surveillance of respiratory infections
dat.ww   = readRDS('../data/nml-database.RDS') |>
  prep.ww(lag = 0)


# Hospital weekly admissions. 
dat.hosp = readRDS('../data/ari_hosp.rds') |>
  prep.hosp(lag = 0, src = 'DAD')


# SARSCOV2 variants in Canada
dat.gen  = readRDS('../data/allprop.rds') |>
  prep.var()

# Merged ww and hosp dataframe
df.merged = inner_join(dat.ww, dat.hosp, 
                       by = join_by(city, date, virus, prov))

# Save merged data frame (for sensitivity analysis)
saveRDS(df.merged, file = '../out/df.merged.rds')

# ---- Prepare model data ----

# Run variant analysis to obtain date windows for model
run_var_analysis(variant.start = "B.1.1.7")

# Retrieve variant date windows
var.dates = readRDS("../out/var-analysis.rds")[["df.dates"]]

# Filter merged df with variant date windows for each city/prov
df.model = create.mdl.obj(df = df.merged,
                          dates = var.dates)

# Retrieve population info (for ms)
pop = calc.pop()

# Calculate counts information (for ms)
counts = summarise.counts(dat.gen)
saveRDS(counts, file = "../out/var_counts.rds")

# Retrieve sublineage data (for ms)
sublineage = retrieve.sublineage()
write.csv(sublineage,
          file = "../ms/tables/Supplementary_Data_1.csv",
          row.names = FALSE)
