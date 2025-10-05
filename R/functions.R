# R/functions.R
# Functions for last criteria project

# === UTILITY FUNCTIONS ===

#' Function for deriving number of specs per case
#' @param df Character vector of diagnosis text
#' @param string Regex pattern to match specimen numbers
#' @return Numeric vector of maximum specimen numbers
spector <- function(df, string) {
  df %>%
    str_extract_all(pattern = string) %>%
    str_extract_all(pattern = "[1-9]") %>%
    map(max) %>%
    as.double()
}

#' Function to make normal case numbers
#' @param x Character vector of case IDs to normalize
#' @return Character vector of normalized case numbers in format "TYPE-YY-NNNNNNN"
caserize <- function(x) {
  casetype <- str_extract(x, "[:alpha:]+")
  yr_case <- str_extract(x, "[:digit:]{9}$")
  yr <- str_extract(yr_case, "^[:digit:]{2}")
  case <- str_extract(yr_case, "[:digit:]{7}$")
  case_number <- as.character(glue("{casetype}-{yr}-{case}"))
  return(case_number)
}

#' Function to get pathologist from the cerner database search results
#' @param path Character vector of pathologist names to clean
#' @return Character vector of cleaned pathologist surnames
pathologer <- function(path) {
  pathologist <- str_replace_all(path, "[:punct:]", "")
  pathologist <- str_replace_all(
    pathologist, 
    "(JD|MD|PhD|DO|MSIV|MB ChB|MBChB|DPhil)", 
    ""
  ) %>% 
    str_trim(side = "both")
  
  begin <- str_locate_all(pathologist, "[:blank:]") %>% 
    map(rbind) %>% 
    map_dbl(max)
  
  end <- str_length(pathologist)
  pathologist <- substr(pathologist, begin, end)
  str_trim(pathologist)
}

#' Filter diagnosis strings by minimum length
#' @param x Character vector of diagnosis strings
#' @return Character vector with strings < 5 characters removed
dxs_map <- function(x) {
  x[str_length(x) >= 5]
}



# === DATA PROCESSING FUNCTIONS ===

#' Create diagnosis-level data from combined case data
#' @param combined_data Data frame with combined case information
#' @return Data frame with one row per individual diagnosis
create_diagnosis_level_data <- function(combined_data) {
  combined_data %>%
    mutate(dxs_tot = map(dxs, dxs_map)) %>%
    unnest(dxs_tot) %>%
    mutate(
      dxs = str_extract(dxs_tot, "(?<=:[[:space:][:cntrl:]])[[:cntrl:][:print:]]+")
    ) %>%
    mutate(
      dxs = ifelse(
        is.na(dxs),
        str_extract(dxs_tot, "(?<=[[:cntrl:]])[[:cntrl:][:print:]]+"),
        dxs
      )
    ) %>%
    filter(
      !str_detect(dxs_tot, "\\[Final Diagnosis\\]"),
      !str_detect(dxs_tot, "DEEPERS PENDING"),
      !str_detect(dxs_tot, "^[[:space:]]+$"),
      !str_detect(dxs_tot, "NAME\\["),
      !is.na(dxs_tot)
    ) %>%
    group_by(case_num) %>%
    mutate(dx_num = row_number()) %>%
    ungroup() %>%
    mutate(rownum = row_number())
}

#' Create corpus data for text analysis
#' @param diagnosis_level_data Data frame with diagnosis-level data
#' @param year_filter Integer year to filter for
#' @return Data frame formatted for corpus creation with text and docid_field columns
create_corpus_data <- function(diagnosis_level_data, year_filter) {
  diagnosis_level_data %>%
    mutate(docid_field = row_number()) %>%
    ungroup() %>%
    filter(year == year_filter) %>%
    filter(!is.na(dxs)) %>%
    filter(dxs != "") %>%
    filter(!is.na(year)) %>%
    select(text = dxs, docid_field)
}

#' Perform hierarchical clustering on diagnosis text
#' @param diagnosis_level_data Data frame with diagnosis-level data
#' @param year_filter Integer year to filter for
#' @return hclust object with hierarchical clustering results
perform_clustering <- function(diagnosis_level_data, year_filter) {
  library(quanteda)
  library(quanteda.textstats)
  
  # Create corpus data
  corpus_df <- diagnosis_level_data %>%
    mutate(docid_field = row_number()) %>%
    ungroup() %>%
    filter(year == year_filter) %>%
    filter(!is.na(dxs)) %>%
    select(text = dxs, docid_field)
  
  # Create corpus
  corp <- corpus(corpus_df)
  
  # Tokenize
  tokens_obj <- tokens(corp, remove_punct = TRUE)
  
  # Create dfm
  dfm_obj <- dfm(tokens_obj) %>%
    dfm_remove(stopwords("english")) %>%
    dfm_wordstem()
  
  # Calculate distances
  dfm_matrix <- as.matrix(dfm_obj)
  simil <- dist(dfm_matrix, method = "euclidean")
  
  # Hierarchical clustering
  cluster_result <- hclust(simil)
  cluster_result$labels <- docnames(dfm_obj)
  
  cluster_result
}

#' Create dendrogram plot from clustering results
#' @param clustering_result hclust object from perform_clustering
#' @param year Integer year for plot title
#' @param output_file Character path to save plot
#' @return ggplot object
create_cluster_plot <- function(clustering_result, year, output_file) {
  library(ggdendro)
  
  dhc <- as.dendrogram(clustering_result)
  ddata <- dendro_data(dhc, type = "rectangle")
  
  p <- ggplot(segment(ddata)) + 
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
    coord_flip() + 
    scale_y_reverse(expand = c(0.2, 0)) +
    theme_dendro() +
    theme(axis.text.x = element_text()) +
    labs(title = paste0("Relationship between diagnoses rendered in ", year, 
                        " in terms of character proximity")) +
    annotate("text", x = -100, y = 5, label = "Euclidean Distance")
  
  ggsave(output_file, p, height = 8, width = 12)
  p
}

# === ANALYSIS FUNCTIONS ===

#' Analyze HPV term usage in diagnoses
#' @param diagnosis_level_data Data frame with diagnosis-level data
#' @return Data frame with HPV indicator column added
analyze_hpv <- function(diagnosis_level_data) {
  diagnosis_level_data %>%
    filter(!is.na(dxs)) %>%
    mutate(
      hpv = str_detect(str_to_upper(dxs), "HPV|HUMAN PAPILLOMA[[:blank:]]*VIRUS|CONDYL|PAPILLOMA")
    )
}

#' Analyze P16 staining usage in reports
#' @param combined_data Data frame with combined case data
#' @return Data frame with P16 indicator column added
analyze_p16 <- function(combined_data) {
  combined_data %>%
    filter(!is.na(report)) %>%
    mutate(p16 = str_detect(str_to_upper(report), "P16"))
}

#' Analyze CIN/dysplasia classification in diagnoses
#' @param diagnosis_level_data Data frame with diagnosis-level data
#' @return Data frame with cin and dysplasia classification columns added
analyze_cin <- function(diagnosis_level_data) {
  diagnosis_level_data %>%
    filter(!is.na(dxs)) %>%
    mutate(dxs = str_replace_all(dxs, "[[:space:]]{5,10}", " ")) %>%
    mutate(
      cin = case_when(
        str_detect(str_to_upper(dxs), "CIN[[:blank:]]*II-III|CIN[[:blank:]]*2-3") ~ "CIN II-III",
        str_detect(str_to_upper(dxs), "CIN[[:blank:]]*III|CIN[[:blank:]]*3") ~ "CIN III",
        str_detect(str_to_upper(dxs), "CIN[[:blank:]]*II(?!I)|CIN[[:blank:]]*2") ~ "CIN II",
        str_detect(str_to_upper(dxs), "CIN[[:blank:]]*I(?!I)|CIN[[:blank:]]*1") ~ "CIN I",
        str_detect(str_to_lower(dxs), "(?<!no[[:space:]]{0,1}|negative for[[:space:]]{0,1})high[[:space:]]{0,1}-{0,1}grade (squamous ){0,1}dysplasia(?! no)") ~ "High Grade",
        TRUE ~ "No CIN"
      )
    ) %>%
    mutate(
      dysplasia = case_when(
        cin == "No CIN" ~ "No CIN",
        cin == "CIN I" ~ "Low Grade",
        TRUE ~ "High Grade"
      )
    )
}

#' Combine CIN analysis with P16 data
#' @param cin_analysis Data frame with CIN analysis results
#' @param p16_analysis Data frame with P16 analysis results
#' @return Data frame with both CIN and P16 information merged
combine_p16_dysplasia <- function(cin_analysis, p16_analysis) {
  cin_analysis %>%
    left_join(
      p16_analysis %>% select(case_num, year, p16),
      by = c("case_num", "year")
    ) %>%
    mutate(p16_ondx = p16)
}

# === STATISTICAL HELPER FUNCTIONS ===

#' Calculate statistics by pre/post period
#' @param data Data frame to analyze
#' @param year_col Column name containing year variable (unquoted)
#' @param metric_col Column name containing metric variable (unquoted)
#' @param cutoff_year Integer year dividing pre and post periods (default 2012)
#' @return Data frame with period statistics
calculate_period_stats <- function(data, year_col, metric_col, cutoff_year = 2012) {
  data %>%
    mutate(
      period = case_when(
        {{year_col}} < cutoff_year ~ "pre",
        {{year_col}} > cutoff_year ~ "post",
        TRUE ~ "washout"
      )
    ) %>%
    filter(period != "washout") %>%
    group_by(period) %>%
    summarise(
      total = n(),
      metric_cases = sum({{metric_col}}, na.rm = TRUE),
      metric_rate = metric_cases / total,
      .groups = "drop"
    )
}

#' Calculate percent drop between two rates
#' @param pre_rate Numeric value for pre-period rate
#' @param post_rate Numeric value for post-period rate
#' @return Numeric percent drop (positive = decrease, negative = increase)
calculate_percent_drop <- function(pre_rate, post_rate) {
  ((pre_rate - post_rate) / pre_rate) * 100
}

#' Calculate fold change between two rates
#' @param post_rate Numeric value for post-period rate
#' @param pre_rate Numeric value for pre-period rate
#' @return Numeric fold change
calculate_fold_change <- function(post_rate, pre_rate) {
  post_rate / pre_rate
}

#' Calculate percentage point change between two rates
#' @param post_rate Numeric value for post-period rate
#' @param pre_rate Numeric value for pre-period rate
#' @return Numeric percentage point change
calculate_percentage_point_change <- function(post_rate, pre_rate) {
  (post_rate - pre_rate) * 100
}

#' Perform chi-square test comparing pre/post periods
#' @param data Data frame to analyze
#' @param year_col Column name containing year variable (unquoted)
#' @param metric_col_name Character string name of metric column
#' @param cutoff_year Integer year dividing pre and post periods (default 2012)
#' @return List containing test results, contingency table, and p-value
perform_chi_square_test <- function(data, year_col, metric_col_name, cutoff_year = 2012) {
  test_data <- data %>%
    mutate(
      period = case_when(
        {{year_col}} < cutoff_year ~ "pre",
        {{year_col}} > cutoff_year ~ "post",
        TRUE ~ "washout"
      )
    ) %>%
    filter(period != "washout", !is.na(.data[[metric_col_name]]))
  
  test_result <- chisq.test(table(test_data$period, test_data[[metric_col_name]]))
  
  chi_table <- test_data %>%
    count(period, .data[[metric_col_name]]) %>%
    pivot_wider(names_from = period, values_from = n, values_fill = 0)
  
  list(
    test = test_result,
    table = chi_table,
    p_value = test_result$p.value
  )
}

# === PLOTTING FUNCTIONS ===



#' Plot HPV diagnosis distribution by year
#' @param hpv_analysis Data frame with HPV analysis results
#' @param output_file Character path to save plot
#' @return ggplot object
plot_hpv_by_year <- function(hpv_analysis, output_file) {
  hpv_data <- hpv_analysis %>%
    filter(!is.na(hpv)) %>%
    count(year, hpv) %>%
    complete(year, hpv, fill = list(n = 0))
  
  p <- ggplot(hpv_data) +
    scale_y_continuous(expand = c(0, 0)) +
    geom_bar(aes(x = year, y = n, fill = hpv), stat = "identity", position = "fill") +
    scale_x_continuous(breaks = seq(2006, 2017)) +
    theme_bw() +
    labs(
      title = "HPV Diagnosis Distribution by Year",
      x = "Year",
      y = "Count",
      fill = "HPV"
    )
  
  ggsave(output_file, p, height = 12, width = 8)
  p
}

#' Plot P16 staining distribution by year
#' @param p16_analysis Data frame with P16 analysis results
#' @param output_file Character path to save plot
#' @return ggplot object
plot_p16_by_year <- function(p16_analysis, output_file) {
  p16_data <- p16_analysis %>%
    count(year, p16) %>%
    complete(year, p16, fill = list(n = 0))
  
  p <- ggplot(p16_data) +
    scale_y_continuous(expand = c(0, 0)) +
    geom_bar(aes(x = year, y = n, fill = p16), stat = "identity", position = "fill") +
    scale_x_continuous(breaks = seq(2006, 2017)) +
    theme_bw() +
    labs(
      title = "P16 Distribution by Year",
      x = "Year",
      y = "Proportion",
      fill = "P16"
    )
  
  ggsave(output_file, p, height = 8, width = 10)
  p
}

#' Plot CIN/dysplasia distribution by year
#' @param cin_analysis Data frame with CIN analysis results
#' @param output_file Character path to save plot
#' @return ggplot object
plot_cin_by_year <- function(cin_analysis, output_file) {
  cin_data <- cin_analysis %>%
    group_by(year) %>%
    count(dysplasia) %>%
    mutate(frac = n/sum(n)) %>%
    ungroup()
  
  p <- ggplot(cin_data, aes(year, y = frac, fill = dysplasia)) +
    scale_x_continuous(breaks = seq(2006, 2017), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(panel.grid = element_blank(), panel.border = element_blank(),
          axis.line = element_blank()) +
    geom_bar(color = "white", stat = "identity", position = 'fill', width = 1) +
    labs(x = "Year", y = "Proportion") +
    theme(legend.title = element_blank())
  
  save_plot(output_file, p, base_height = 12)
  p
}

#' Plot P16 and dysplasia panel figure
#' @param p16_dysplasia_combined Data frame with combined P16 and dysplasia data
#' @param output_file Character path to save plot
#' @return ggplot object
plot_p16_dysplasia_panel <- function(p16_dysplasia_combined, output_file) {
  p <- p16_dysplasia_combined %>%
    filter(year > 2009) %>%
    ggplot() +
    geom_bar(aes(x = year, fill = dysplasia), position = "fill", color = "white", width = 0.9) +
    facet_grid(rows = vars(p16_ondx), 
               labeller = labeller(p16_ondx = c("FALSE" = "FALSE", "TRUE" = "TRUE"))) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1, 0.25)) +
    scale_x_continuous(breaks = seq(2010, 2017, 2)) +
    scale_fill_manual(
      values = c("High Grade" = "#E57373", "Low Grade" = "#66BB6A", "No CIN" = "#64B5F6"),
      labels = c("High Grade" = "High grade squamous\nintraepithelial lesion (HSIL)",
                 "Low Grade" = "Low grade squamous\nintraepithelial lesion (LSIL)",
                 "No CIN" = "No squamous\nintraepithelial lesion")
    ) +
    labs(
      x = "Year (2000's)",
      y = "Fraction of Total Cases",
      fill = NULL
    ) +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
      strip.background = element_rect(fill = "grey80", color = NA),
      strip.text = element_text(size = 14, face = "bold"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.5, "lines"),
      axis.text = element_text(size = 13),
      axis.title = element_text(size = 14, face = "bold")
    )
  
  save_plot(output_file, p, base_height = 12)
  p
}


#' Create pathologist streamgraph visualization
#' @param pathologist_volume_by_year Data frame with year, pathologist, and case counts
#' @return List with html file path, png file path, and streamgraph widget
create_pathologist_streamgraph <- function(pathologist_volume_by_year) {
  stream_data <- pathologist_volume_by_year %>%
    filter(year >= 2006, year <= 2016) %>%
    group_by(pathologist) %>%
    filter(sum(n) >= 50) %>%
    ungroup() %>%
    complete(year, pathologist, fill = list(n = 0)) %>%
    arrange(year, pathologist) %>%
    rename(key = pathologist, value = n, date = year)
  
  sg <- stream_data %>%
    streamgraph(key = "key", value = "value", date = "date", offset = "zero") %>%
    sg_legend(show = FALSE) %>%
    sg_axis_x(1, "year", "%Y") %>%
    sg_colors("category20")
  
  html_file <- "outputs/plots/pathologist_streamgraph.html"
  htmlwidgets::saveWidget(sg, file = html_file, selfcontained = TRUE)
  
  png_file <- "outputs/plots/figure_2.png"
  webshot2::webshot(
    html_file, 
    png_file,
    delay = 1,
    vwidth = 1200,
    vheight = 800,
    cliprect = c(0, 50, 1200, 750)
  )
  
  list(html = html_file, png = png_file, widget = sg)
}

# === EXPORT FUNCTIONS ===

#' Export HPV analysis summary to CSV
#' @param hpv_analysis Data frame with HPV analysis results
#' @param output_path Character path for output CSV file
#' @return Character path to saved file
export_hpv_summary <- function(hpv_analysis, output_path) {
  hpv_summary <- hpv_analysis %>%
    group_by(year, hpv) %>%
    summarise(count = n(), .groups = "drop") %>%
    pivot_wider(names_from = hpv, values_from = count, values_fill = 0)
  
  write_csv(hpv_summary, output_path)
  output_path
}

#' Export CIN analysis summary to CSV
#' @param cin_analysis Data frame with CIN analysis results
#' @param output_path Character path for output CSV file
#' @return Character path to saved file
export_cin_summary <- function(cin_analysis, output_path) {
  cin_summary <- cin_analysis %>%
    count(year, dysplasia) %>%
    pivot_wider(names_from = dysplasia, values_from = n, values_fill = 0)
  
  write_csv(cin_summary, output_path)
  output_path
}

#' Export statistical test results summary
#' @param hpv_test List with HPV chi-square test results
#' @param cin_test List with CIN chi-square test results
#' @param output_path Character path for output CSV file
#' @return Character path to saved file
export_statistical_test_summary <- function(hpv_test, cin_test, output_path) {
  require(readr)
  test_summary <- tibble(
    test_name = c("HPV by Era", "CIN by Era")) %>% 
    mutate(p_value = c(hpv_test$test$p.value, cin_test$test$p.value)
  ) %>%
    mutate(significant = p_value < 0.05)
  
  write_csv(test_summary, output_path)
  output_path
}

# === REPORTING FUNCTIONS ===

#' Create abstract statistics report
#' @param hpv_analysis Data frame with HPV analysis results
#' @param p16_analysis Data frame with P16 analysis results
#' @param cin_analysis Data frame with CIN analysis results
#' @param output_txt Character path for output text file
#' @param output_csv Character path for output CSV file
#' @return List with txt_file path, csv_file path, and summary stats tibble
create_abstract_statistics_report <- function(hpv_analysis, p16_analysis, cin_analysis, 
                                              output_txt, output_csv) {
  # Total diagnoses (excluding washout year 2012)
  total_diagnoses <- hpv_analysis %>% filter(year != 2012) %>% nrow()
  
  # HPV term usage with helper functions
  hpv_stats <- calculate_period_stats(hpv_analysis, year, hpv)
  hpv_change_pct <- calculate_percent_drop(
    hpv_stats %>% filter(period == "pre") %>% pull(metric_rate),
    hpv_stats %>% filter(period == "post") %>% pull(metric_rate)
  )
  
  # P16 staining with helper functions
  p16_stats <- calculate_period_stats(p16_analysis, year, p16)
  p16_fold_change <- calculate_fold_change(
    p16_stats %>% filter(period == "post") %>% pull(metric_rate),
    p16_stats %>% filter(period == "pre") %>% pull(metric_rate)
  )
  
  # CIN diagnosis changes (excluding 2012)
  cin_stats <- cin_analysis %>%
    mutate(period = case_when(
      year < 2012 ~ "pre",
      year > 2012 ~ "post",
      TRUE ~ "washout"
    )) %>%
    filter(period != "washout") %>%
    count(period, dysplasia) %>%
    group_by(period) %>%
    mutate(
      total = sum(n),
      rate = n / total
    ) %>%
    ungroup()
  
  cin_wide <- cin_stats %>%
    select(period, dysplasia, rate) %>%
    pivot_wider(names_from = dysplasia, values_from = rate)
  
  low_grade_change <- calculate_percentage_point_change(
    cin_wide %>% filter(period == "post") %>% pull(`Low Grade`),
    cin_wide %>% filter(period == "pre") %>% pull(`Low Grade`)
  )
  
  no_cin_change <- calculate_percentage_point_change(
    cin_wide %>% filter(period == "post") %>% pull(`No CIN`),
    cin_wide %>% filter(period == "pre") %>% pull(`No CIN`)
  )
  
  high_grade_change <- calculate_percentage_point_change(
    cin_wide %>% filter(period == "post") %>% pull(`High Grade`),
    cin_wide %>% filter(period == "pre") %>% pull(`High Grade`)
  )
  
  # Statistical tests (excluding 2012)
  cin_test_data <- cin_analysis %>%
    mutate(period = case_when(
      year < 2012 ~ "pre",
      year > 2012 ~ "post",
      TRUE ~ "washout"
    )) %>%
    filter(period != "washout")
  
  chisq_low <- cin_test_data %>%
    mutate(is_low_grade = dysplasia == "Low Grade") %>%
    {chisq.test(table(.$period, .$is_low_grade))}
  
  chisq_no_cin <- cin_test_data %>%
    mutate(is_no_cin = dysplasia == "No CIN") %>%
    {chisq.test(table(.$period, .$is_no_cin))}
  
  chisq_high <- cin_test_data %>%
    mutate(is_high_grade = dysplasia == "High Grade") %>%
    {chisq.test(table(.$period, .$is_high_grade))}
  
  # Create report text
  report <- str_c(
    "ABSTRACT STATISTICS REPORT\n",
    "Generated: ", Sys.time(), "\n",
    strrep("=", 80), "\n\n",
    
    "STUDY DESIGN\n",
    strrep("-", 40), "\n",
    "Note: Year 2012 excluded as washout period\n",
    "Pre-LAST period: 2006-2011\n",
    "Post-LAST period: 2013-2016\n\n",
    
    "TOTAL DIAGNOSES\n",
    strrep("-", 40), "\n",
    "Total cervical biopsy diagnoses analyzed: ", 
    format(total_diagnoses, big.mark = ","), "\n",
    "Study period: 11 years (2006-2016, excluding 2012)\n",
    "Pre-LAST period (2006-2011): ", 
    format(hpv_stats %>% filter(period == "pre") %>% pull(total), big.mark = ","), 
    " diagnoses\n",
    "Post-LAST period (2013-2016): ", 
    format(hpv_stats %>% filter(period == "post") %>% pull(total), big.mark = ","), 
    " diagnoses\n\n",
    
    "HPV TERM USAGE\n",
    strrep("-", 40), "\n",
    "Pre-LAST HPV term rate: ", 
    round((hpv_stats %>% filter(period == "pre") %>% pull(metric_rate)) * 100, 2), "%\n",
    "Post-LAST HPV term rate: ", 
    round((hpv_stats %>% filter(period == "post") %>% pull(metric_rate)) * 100, 2), "%\n",
    "Percent drop in HPV term usage: ", round(hpv_change_pct, 1), "%\n\n",
    
    "P16 STAINING\n",
    strrep("-", 40), "\n",
    "Pre-LAST p16 rate: ", 
    round((p16_stats %>% filter(period == "pre") %>% pull(metric_rate)) * 100, 2), "%\n",
    "Post-LAST p16 rate: ", 
    round((p16_stats %>% filter(period == "post") %>% pull(metric_rate)) * 100, 2), "%\n",
    "Fold increase in p16 use: ", round(p16_fold_change, 1), "x\n\n",
    
    "DIAGNOSIS RATE CHANGES (Post - Pre)\n",
    strrep("-", 40), "\n",
    "Low-grade diagnoses: ", round(low_grade_change, 1), 
    " percentage points (p ", format.pval(chisq_low$p.value, digits = 3, eps = 0.001), ")\n",
    "Non-CIN diagnoses: ", round(no_cin_change, 1), 
    " percentage points (p ", format.pval(chisq_no_cin$p.value, digits = 3, eps = 0.001), ")\n",
    "High-grade diagnoses: ", round(high_grade_change, 1), 
    " percentage points (p ", format.pval(chisq_high$p.value, digits = 3, eps = 0.001), ")\n\n",
    
    "DIAGNOSIS DISTRIBUTION BY PERIOD\n",
    strrep("-", 40), "\n"
  )
  
  # Add detailed breakdown table
  detailed_table <- cin_stats %>%
    mutate(
      rate_pct = str_c(round(rate * 100, 2), "%"),
      count = str_c(n, " / ", total)
    ) %>%
    select(period, dysplasia, count, rate_pct)
  
  # Write report
  write_lines(report, output_txt)
  write_tsv(detailed_table, output_txt, append = TRUE)
  
  # Create CSV with key statistics
  summary_stats <- tibble(
    metric = c("Total Diagnoses (excl 2012)", "HPV Term Drop (%)", "P16 Fold Increase", 
               "Low Grade Change (pp)", "No CIN Change (pp)", "High Grade Change (pp)",
               "Low Grade p-value", "No CIN p-value", "High Grade p-value"),
    value = c(total_diagnoses, round(hpv_change_pct, 1), round(p16_fold_change, 1),
              round(low_grade_change, 1), round(no_cin_change, 1), round(high_grade_change, 1),
              chisq_low$p.value, chisq_no_cin$p.value, chisq_high$p.value)
  )
  
  write_csv(summary_stats, output_csv)
  
  list(txt_file = output_txt, csv_file = output_csv, stats = summary_stats)
}

#' Create results section report
#' @param hpv_analysis Data frame with HPV analysis results
#' @param p16_analysis Data frame with P16 analysis results
#' @param cin_analysis Data frame with CIN analysis results
#' @param output_txt Character path for output text file
#' @param output_csv Character path for output CSV file
#' @return List with txt_file path, csv_file path, and summary stats tibble
create_results_section_report <- function(hpv_analysis, p16_analysis, cin_analysis,
                                          output_txt, output_csv) {
  # PERIOD COUNTS
  pre_2012_count <- hpv_analysis %>% filter(year <= 2012) %>% nrow()
  post_2012_count <- hpv_analysis %>% filter(year > 2012) %>% nrow()
  final_diagnoses <- nrow(hpv_analysis)
  
  # HPV TERM USAGE
  hpv_stats <- calculate_period_stats(hpv_analysis, year, hpv)
  hpv_pre_count <- hpv_stats$metric_cases[hpv_stats$period == "pre"]
  hpv_pre_rate <- hpv_stats$metric_rate[hpv_stats$period == "pre"]
  hpv_post_count <- hpv_stats$metric_cases[hpv_stats$period == "post"]
  hpv_post_rate <- hpv_stats$metric_rate[hpv_stats$period == "post"]
  hpv_drop <- calculate_percent_drop(hpv_pre_rate, hpv_post_rate)
  
  hpv_test <- perform_chi_square_test(hpv_analysis, year, "hpv")
  
  # P16 STAINING
  p16_stats <- calculate_period_stats(p16_analysis, year, p16)
  p16_pre_count <- p16_stats$metric_cases[p16_stats$period == "pre"]
  p16_pre_rate <- p16_stats$metric_rate[p16_stats$period == "pre"]
  p16_post_count <- p16_stats$metric_cases[p16_stats$period == "post"]
  p16_post_rate <- p16_stats$metric_rate[p16_stats$period == "post"]
  p16_fold <- calculate_fold_change(p16_post_rate, p16_pre_rate)
  
  p16_test <- perform_chi_square_test(p16_analysis, year, "p16")
  
  # CIN DIAGNOSIS CHANGES
  cin_stats <- cin_analysis %>%
    mutate(period = ifelse(year <= 2012, "pre", "post")) %>%
    group_by(period, dysplasia) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(period) %>%
    mutate(
      total = sum(n),
      rate = n / total
    )
  
  cin_wide <- cin_stats %>%
    select(period, dysplasia, rate) %>%
    pivot_wider(names_from = dysplasia, values_from = rate)
  
  low_grade_change <- calculate_percentage_point_change(
    cin_wide$`Low Grade`[cin_wide$period == "post"],
    cin_wide$`Low Grade`[cin_wide$period == "pre"]
  )
  
  no_cin_change <- calculate_percentage_point_change(
    cin_wide$`No CIN`[cin_wide$period == "post"],
    cin_wide$`No CIN`[cin_wide$period == "pre"]
  )
  
  high_grade_change <- calculate_percentage_point_change(
    cin_wide$`High Grade`[cin_wide$period == "post"],
    cin_wide$`High Grade`[cin_wide$period == "pre"]
  )
  
  # Statistical tests for CIN changes
  cin_test_data <- cin_analysis %>%
    mutate(period = ifelse(year <= 2012, "pre", "post"))
  
  chisq_low <- chisq.test(table(cin_test_data$period, cin_test_data$dysplasia == "Low Grade"))
  chisq_no_cin <- chisq.test(table(cin_test_data$period, cin_test_data$dysplasia == "No CIN"))
  chisq_high <- chisq.test(table(cin_test_data$period, cin_test_data$dysplasia == "High Grade"))
  
  # CREATE REPORT
  report <- paste0(
    "RESULTS SECTION REPORT\n",
    "Generated: ", Sys.time(), "\n",
    strrep("=", 80), "\n\n",
    
    "TEMPORAL DISTRIBUTION\n",
    strrep("-", 80), "\n",
    "Total individual diagnoses analyzed: ", format(final_diagnoses, big.mark = ","), "\n",
    "Study period: 11 years (2006-2016)\n",
    "Individual specimen diagnoses before 2012: ", format(pre_2012_count, big.mark = ","), "\n",
    "Individual specimen diagnoses after 2012: ", format(post_2012_count, big.mark = ","), "\n\n",
    
    "HPV TERM USAGE\n",
    strrep("-", 80), "\n",
    "Before 2012:\n",
    "  - HPV term diagnoses: ", format(hpv_pre_count, big.mark = ","), " (", 
    round(hpv_pre_rate * 100, 1), "%)\n",
    "After 2012:\n",
    "  - HPV term diagnoses: ", format(hpv_post_count, big.mark = ","), " (", 
    round(hpv_post_rate * 100, 1), "%)\n",
    "Change: ", round(hpv_drop, 0), "% drop from baseline (P ", 
    format.pval(hpv_test$p.value, digits = 3, eps = 0.001), ")\n\n",
    
    "P16 IMMUNOHISTOCHEMICAL STAINING\n",
    strrep("-", 80), "\n",
    "Before 2012:\n",
    "  - P16 positive diagnoses: ", format(p16_pre_count, big.mark = ","), " (", 
    round(p16_pre_rate * 100, 1), "%)\n",
    "After 2012:\n",
    "  - P16 positive diagnoses: ", format(p16_post_count, big.mark = ","), " (", 
    round(p16_post_rate * 100, 1), "%)\n",
    "Change: ", round(p16_fold, 1), "-fold increase (P ", 
    format.pval(p16_test$p.value, digits = 3, eps = 0.001), ")\n\n",
    
    "DIAGNOSTIC CATEGORY CHANGES (After LAST Publication)\n",
    strrep("-", 80), "\n",
    "Low-grade diagnoses:\n",
    "  - Change: ", round(low_grade_change, 1), " percentage points\n",
    "  - P-value: ", format.pval(chisq_low$p.value, digits = 3, eps = 0.001), "\n",
    "Non-CIN diagnoses:\n",
    "  - Change: ", round(no_cin_change, 1), " percentage points\n",
    "  - P-value: ", format.pval(chisq_no_cin$p.value, digits = 3, eps = 0.001), "\n",
    "High-grade diagnoses:\n",
    "  - Change: ", round(high_grade_change, 1), " percentage points\n",
    "  - P-value: ", format.pval(chisq_high$p.value, digits = 3, eps = 0.001), "\n\n",
    
    "DETAILED BREAKDOWN BY PERIOD\n",
    strrep("-", 80), "\n"
  )
  
  # Detailed table
  detailed_table <- cin_stats %>%
    mutate(
      rate_pct = paste0(round(rate * 100, 2), "%"),
      count = paste0(format(n, big.mark = ","), " / ", format(total, big.mark = ","))
    ) %>%
    select(Period = period, Dysplasia = dysplasia, Count = count, Rate = rate_pct)
  
  # Write text report
  writeLines(report, output_txt)
  write.table(detailed_table, output_txt, append = TRUE, quote = FALSE, 
              row.names = FALSE, col.names = TRUE, sep = "\t")
  
  # Create detailed CSV
  results_summary <- tibble(
    section = c(rep("Temporal", 3), 
                rep("HPV Terms", 4), rep("P16 Staining", 4),
                rep("CIN Changes", 6)),
    metric = c("Total Diagnoses", "Pre-2012 Diagnoses", "Post-2012 Diagnoses",
               "Pre-2012 HPV Count", "Pre-2012 HPV Rate (%)", 
               "Post-2012 HPV Count", "Post-2012 HPV Rate (%)",
               "Pre-2012 P16 Count", "Pre-2012 P16 Rate (%)",
               "Post-2012 P16 Count", "Post-2012 P16 Rate (%)",
               "Low Grade Change (pp)", "Low Grade P-value",
               "No CIN Change (pp)", "No CIN P-value",
               "High Grade Change (pp)", "High Grade P-value"),
    value = c(final_diagnoses, pre_2012_count, post_2012_count,
              hpv_pre_count, round(hpv_pre_rate * 100, 1),
              hpv_post_count, round(hpv_post_rate * 100, 1),
              p16_pre_count, round(p16_pre_rate * 100, 1),
              p16_post_count, round(p16_post_rate * 100, 1),
              round(low_grade_change, 1), chisq_low$p.value,
              round(no_cin_change, 1), chisq_no_cin$p.value,
              round(high_grade_change, 1), chisq_high$p.value)
  )
  
  write_csv(results_summary, output_csv)
  
  list(txt_file = output_txt, csv_file = output_csv, stats = results_summary)
}