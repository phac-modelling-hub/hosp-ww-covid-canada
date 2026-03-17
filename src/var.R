# ---- Functions for variant analysis and inclusion into model ---- 

# Calculating peak proportions of dom variants
#' @param dat Dataframe of variant counts
#' @param c City name
calc.prop <- function(dat, c){
  message(paste0("Calculating variant proportions for ", c))
  df = dat |> 
    filter(city == c) |>
    group_by(date, city, prov) |>
    summarise(
      dom.variant = unique(dom.variant),
      dom.prop = var.prop[variant == dom.variant],
      ndom.prop = sum(var.prop, na.rm = TRUE) - dom.prop,
      .groups = 'drop') |> 
    group_by(dom.variant, city, prov) |>
    summarise(
      maxdomprp = max(dom.prop), # max proportion of dominant variant 
      date.maxprp = first(date[dom.prop ==
                                 maxdomprp]), # date of max proportion
      ndom.prop = unique(ndom.prop[date == 
                                    date.maxprp]),# proportion of other variants
                                                  # on day of max dominant prop
      pct.diff = 
        (maxdomprp - ndom.prop)/((maxdomprp + 
                                    ndom.prop)/2), # percentage difference
                                                   # of dominant prop
      .groups = 'drop'
    )
  
  return(df)
  
}

#' @description
#' Execution of variant analysis
#' 
#' @param variant.start Starting dominant variant in dataframe
#' @param p Province
#' @param var.rm Variants to remove from dataframe
#'
#' @return
run_var_analysis <- function(variant.start, 
                             p = c("BC", "AB", "SK", "MB", 
                                    "ON", "NS", "NL"),
                             end.date = '2024-03-31', # for visualization purposes, 
                                                      # variant data beyond DAD end date will not be included
                             var.rm = c("Unassigned", "None")) {
  
  if(0){  #DEBUG
    variant.start = "B.1.1.7"
    var.rm = c("Unassigned", "None")
  }
  
  # Remove the variants that were not assigned a lineage.
  dat = readRDS('../out/var.rds') |>
    filter(!variant %in% var.rm,
           prov %in% p,
           date <= end.date,
           !dom.variant %in% var.rm) |>
    mutate(dom.variant = factor(dom.variant, levels = unique(dom.variant)))
  
  var.analysis(dat, variant.start)
}


#' @description
#' Analysis and processing of variant data for inclusion in model.
#' 
#' @param d Variant dataframe
#' @param startvar Starting dominant variant
#'
#' @return
#' 
var.analysis <- function(d, startvar){
  message("Analyzing variant data...")
  cities = unique(d$city)
  
  # Filter data based on starting variant
  d.filtered = filter.var(d, startvar) |>
    filter(var.prop >= 0.02)
  
  # Calculate peak proportion for each variant
  df.dom = lapply(cities, calc.prop, dat = d.filtered) |>
    bind_rows()
  
  # filter ts with only dom.variants
  df = right_join(d.filtered, df.dom, by = c("city", "prov",
                                          "variant" = "dom.variant")) |>
    filter(variant == dom.variant)
  
  # Save dom variant df (for sensitivity analysis)
  saveRDS(df, file = '../out/domvars.rds')
                  
  # calc start and end dates of dominant variants 
  dates = get_var_dates(df)

  # get dataframe of variants used in model
  df.model = dates |>
    full_join(df, by = c("city", "prov", "variant")) |>
    group_by(city, prov, variant) |>
    filter(date >= start & date <= end) |>
    ungroup()
  
  # ---- plots ----
  
  message("Generating plots...")

  g.periods = plot_domvariant_periods(dates)
  
  pdf(file = paste0("../figs/dom-variant-periods-", today(), ".pdf"), 
      height = 9, width = 9)
  plot(g.periods[[2]])
  dev.off()
  ggsave(paste0("../ms/figs/dom-variant-periods-", today(), ".png"),
         g.periods[[2]],
         height = 9, width = 9)
  pdf(file = paste0("../figs/dom-variant-periods-grid-", today(), ".pdf"), 
      height = 9, width = 22)
  plot(g.periods[[1]])
  dev.off() 
  ggsave(paste0("../ms/figs/dom-variant-periods-grid-", today(), ".png"),
         g.periods[[1]],
         height = 9, width = 12)
  
  out = list(ts.dom = df.model,
             df.dates = dates)
  
  fname = 'out/var-analysis.rds'
  
  saveRDS(out, file = paste0('../', fname))
  message(paste("Analysis saved in", fname))
  message("Complete.")
}

# Calculate date range of variants for each location
#' @param df Dataframe
#' @param sensitivity Logical flag for function use in sensitivity analysis
#' @param var.prop Minimum variant prop. threshold (for sensitivity analysis)
#' @param var.window Minimum date window (for sensitivity analysis)
#' @param var.pctabove Minimum weeks above threshold (for sensitivity analysis)
get_var_dates <- function(df,
                          sensitivity = FALSE,
                          var.prop = NULL,
                          var.window = NULL,
                          var.pctabove = NULL){
  
  if(!sensitivity){
    date.prm = fetch_var_params()
  }
  
  else if(sensitivity){
    date.prm = tibble(
      var.prop = var.prop,
      var.window = var.window,
      var.pctabove = var.pctabove
    )
  }
  
  dates = df |>
    group_by(city, prov, variant) |>
    complete(date = seq(min(date), max(date), by = "week")) |> # avoid skipping weeks
                                                              # in row calculations
    replace_na(list(var.prop = 0)) |> # arbitrary 0 for dom calc
                                # variant for that week isn't actually 0, 
                                # but it's not a dom variant that week so should
                                # not be captured in week range
    summarize(
      start_end_dates = list(find_valid_date_range(
        data = pick(everything()), 
        column = "var.prop", 
        threshold = date.prm$var.prop, 
        min_weeks = date.prm$var.window,
        pct.above = date.prm$var.pctabove
      )),
      .groups = 'drop'
    ) |> 
    unnest_wider(start_end_dates) |>
    drop_na(start) |>
    arrange(start, prov) |>
    mutate(variant = factor(variant, levels = unique(variant)))
  
  if(sensitivity){
    dates = dates |>
      mutate(
        threshold = date.prm$var.prop,
        min_weeks = date.prm$var.window,
        pct.above = date.prm$var.pctabove
      )
  }
  
  return(dates)
}

#' @description
#' Finding valid date range for a specified variant
#' @param data Dataframe of variant
#' @param column Column name used for analysis
#' @param threshold Minimum proportion
#' @param min_weeks Minimum weeks window
#' @param pct.above Percentage of weeks in window above threshold
find_valid_date_range <- function(data, column, threshold,
                                  min_weeks = 8, pct.above = 0.9) {
  
  # Ensure the data is sorted by date
  data = data[order(data$date), ]
  n = nrow(data)
  
  # If the dataset is too small, return NULL immediately
  if (n < min_weeks) return(NULL)
  
  # Store the best range found
  best_range = NULL
  best_length = 0
  
  # Iterate over all possible subranges
  for (start_row in 1:(n - min_weeks + 1)) {
    for (end_row in (start_row + min_weeks - 1):n) {
      subset_data = data[start_row:end_row, ]
      
      # Check if the subset meets the required percentage above the threshold
      if (sum(subset_data[[column]] > threshold) >= nrow(subset_data) * pct.above) {
        
        # If the first row is below the threshold, find the first valid row
        if(subset_data[[column]][1] < threshold){
          first_valid_index = min(which(subset_data[[column]] >= threshold))
          subset_data = subset_data[first_valid_index:nrow(subset_data),]
        }
        
        # If the last row is below the threshold, find the last valid row
        if (subset_data[[column]][nrow(subset_data)] < threshold) {
          last_valid_index = max(which(subset_data[[column]] >= threshold))
          subset_data = subset_data[1:last_valid_index, ]  # Trim the subset
        }
        
        # Check if trimmed subset meets minimum required rows
        if(nrow(subset_data) < min_weeks * pct.above){
          next
        }
        
        # If this range is longer than the best found so far, update it
        if ((nrow(subset_data)) > best_length) {
          best_range = list(
            start = first(subset_data$date),
            end = last(subset_data$date)
          )
          best_length = nrow(subset_data)
        }
      }
    }
  }
  
  # Return the best range found, or NULL if none found
  return(best_range)
}
