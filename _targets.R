# _targets.R
library(targets)
library(tarchetypes)

# Source functions
source("R/functions.R")

# Set target options
tar_option_set(
  packages = c(
    "readxl", "tibble", "dplyr", "stringr", "anytime", "purrr", 
    "tidyr", "ggplot2", "viridis", "testthat", "glue", "scales",
    "cowplot", "RColorBrewer", "readr", "streamgraph", "htmlwidgets", 
    "webshot2", "quanteda", "quanteda.textstats", "ggdendro"
  )
)

# Create output directories
if (!dir.exists("outputs")) dir.create("outputs")
if (!dir.exists("outputs/plots")) dir.create("outputs/plots")
if (!dir.exists("outputs/data")) dir.create("outputs/data")
if (!dir.exists("outputs/reports")) dir.create("outputs/reports")

# Define targets
list(
  tar_target(
    functions_file,
    {
      source("R/functions.R")
      "functions loaded"
    }
  ),
  # Data file discovery targets
  tar_target(
    cerner_files,
    {
      files <- dir(path = "data/cerner") %>%
        data_frame() %>%
        filter(. != "real") %>%
        unlist() %>%
        as.character()
      files
    }
  ),
  
  tar_target(
    ties_files,
    {
      files <- dir(path = "data/ties") %>%
        data_frame() %>%
        filter(. != "old") %>%
        unlist() %>%
        as.character()
      files
    }
  ),
  
  # Raw data loading targets
  tar_target(
    cerner_raw,
    {
      cerner_files %>%
        map(function(x) {
          read_xlsx(
            file.path("data/cerner", x),
            col_types = c("text", "skip", "skip", "skip", "skip", "skip", 
                          "text", "text", "text", "skip", "skip", "skip", 
                          "skip", "numeric")
          )
        }) %>%
        reduce(rbind)
    }
  ),
  
  tar_target(
    ties_raw,
    {
      map_df(ties_files, function(x) {
        read_xlsx(file.path("data/ties", x))
      }) %>%
        select(
          record = `Record ID`,
          year = `Event Year`,
          date = `Event Date`,
          report = 8
        )
    }
  ),
  
  # Cleaned data targets
  tar_target(
    cerner_clean,
    {
      cer <- cerner_raw %>%
        filter(
          !str_detect(str_to_upper(`Accession#`), "PS"),
          !str_detect(str_to_upper(`Accession#`), "SB"),
          !str_detect(str_to_upper(`Accession#`), "OS")
        ) %>%
        select(-`Verified Date`, -Note, -Addendum) %>%
        distinct() %>%
        filter(!str_detect(`Accession#`, "HC-"))
      
      # Handle duplicates by selecting rows with minimum NA count
      cer$na_count <- apply(cer, 1, function(x) sum(is.na(x)))
      cer <- cer %>% 
        group_by(`Accession#`) %>%
        filter(na_count == min(na_count)) %>%
        filter(row_number() == 1) %>%
        ungroup()
      
      cer
    }
  ),
  
  tar_target(
    ties_clean,
    {
      ties <- ties_raw %>%
        mutate(year = as.numeric(year)) %>%
        # Make diagnosis column
        mutate(
          diagnosis = str_extract(
            report, 
            "\\[Final Diagnosis\\](.|[[:cntrl:]])+?(Microscopic Description|Gross Description)"
          ),
          repo_len = str_length(diagnosis),
          report_len = str_length(report)
        ) %>%
        # Filter outside cases
        filter(
          !str_detect(str_to_upper(diagnosis), "SLIDES"),
          !str_detect(str_to_upper(record), "PS"),
          !str_detect(str_to_upper(record), "SB"),
          !str_detect(str_to_upper(record), "OS"),
          !str_detect(str_to_upper(report), "OUTSIDE SURGICAL")
        ) %>%
        distinct()
      
      # Identify report types
      ties <- ties %>%
        mutate(
          pre_07 = ifelse(
            str_detect(diagnosis, "(FINAL DIAGNOSIS|MICROSCOPIC DIAGNOSIS):[[:cntrl:]]"), 
            1, 0
          ),
          is_07 = ifelse(
            str_detect(diagnosis, "FINAL DIAGNOSIS:[^[[:cntrl:]]]"), 
            1, 0
          ),
          post_07 = ifelse(
            str_detect(diagnosis, " \\.{0,1}Final Diagnosis") | year == 2016, 
            1, 0
          )
        ) %>%
        mutate(
          pre_07 = ifelse(year > 2008, 0, pre_07),
          is_07 = ifelse(year %in% c(2007, 2008) & (pre_07 + is_07 + post_07 < 1), 1, is_07)
        ) %>%
        mutate(
          reporttype = case_when(
            pre_07 == 1 ~ 'pre_07',
            is_07 == 1 ~ 'is_07',
            post_07 == 1 ~ 'post_07'
          )
        )
      
      ties
    }
  ),
  
  # Parse specimens and diagnoses
  tar_target(
    ties_parsed,
    {
      ties <- ties_clean %>%
        # Additional filtering
        filter(
          !str_detect(str_to_upper(diagnosis), "LYMPH"),
          !str_detect(str_to_upper(diagnosis), "OVAR"),
          !str_detect(str_to_upper(diagnosis), "UTERUS AND CERVIX"),
          !str_detect(str_to_upper(diagnosis), "HYSTERECTOMY"),
          record != 'null',
          !str_detect(record, "^NA")
        ) %>%
        # Get number of specimens and split diagnoses
        mutate(
          specs = case_when(
            reporttype == "pre_07" ~ spector(diagnosis, "\\n[0-9]\\.[[:space:]](?!o'clock)"),
            reporttype == "is_07" | is.na(reporttype) ~ spector(diagnosis, "[[:space:]][0-9][\\.]?[[:space:]](?!o'clock)"),
            reporttype == "post_07" ~ spector(diagnosis, "\\n.{1,5}[\\.]?[1-9][\\.|[[:space:]]](?!o'clock)")
          ),
          dxs = case_when(
            reporttype == "pre_07" ~ str_split(diagnosis, "(\\n[0-9]\\.[[:space:]](?!o'clock))|Gross Description"),
            reporttype == "is_07" | is.na(reporttype) ~ str_split(diagnosis, "([[:space:]][0-9][\\.]?[[:space:]](?!o'clock))|Gross Description"),
            reporttype == "post_07" ~ str_split(diagnosis, "(\\n.{1,5}[\\.]?[1-9][\\.|[[:space:]]](?!o'clock))|The case material was reviewed and the report verified by")
          ),
          specs_dx = (dxs %>% map_int(length)) - 2,
          dxs = map(dxs, as.character)
        ) %>%
        # Make normalized case numbers
        mutate(
          case_num = case_when(
            reporttype %in% c("pre_07", "is_07") | is.na(reporttype) ~ record,
            reporttype == "post_07" ~ caserize(record)
          )
        ) %>%
        filter(case_num != 'NA-NA-NA') %>%
        # Handle duplicates
        group_by(case_num) %>%
        filter(repo_len == max(repo_len)) %>%
        filter(report_len == max(report_len)) %>%
        filter(row_number() == 1) %>%
        ungroup()
      
      ties
    }
  ),
  
  # Combine TIES and Cerner data
  tar_target(
    combined_data,
    {
      tot <- ties_parsed %>%
        left_join(cerner_clean, by = c("case_num" = "Accession#")) %>%
        mutate(pathologist = pathologer(Pathologist)) %>%
        mutate(
          pathologist = case_when(
            str_detect(pathologist, "Carolin") ~ "Reyes",
            str_detect(pathologist, "^E$") ~ "Schwartz",
            TRUE ~ pathologist
          )
        ) %>%
        mutate(dxs_len = repo_len / specs) %>%
        filter(!is.na(pathologist)) %>%
        # Clean diagnosis text
        mutate(
          diagnosis_mut = str_replace(
            diagnosis,
            regex("^\\[Final Diagnosis\\][[:cntrl:][:print:]]+Final Diagnosis[:]?", ignore_case = TRUE),
            ""
          )
        ) %>%
        mutate(
          diagnosis_mut = str_replace(
            diagnosis_mut,
            regex("^\\[Final Diagnosis\\][[:cntrl:][:print:]]+MICROSCOPIC DIAGNOSIS[:]?", ignore_case = TRUE),
            ""
          )
        ) %>%
        mutate(
          diagnosis_mut = str_replace(
            diagnosis_mut,
            regex("The case material was[[:cntrl:][:print:]]+|\\[Gross Description+"),
            ""
          )
        ) %>%
        mutate(
          diagnosis_mut = str_replace(diagnosis_mut, "The diagnosis has been revised as follows:", ""),
          diagnosis_mut = str_replace(diagnosis_mut, "_", "")
        ) %>%
        mutate(
          dxs = case_when(
            reporttype == "pre_07" ~ str_split(diagnosis_mut, "(\\n[0-9]\\.[[:space:]](?!o'clock))|Gross Description"),
            reporttype == "is_07" ~ str_split(diagnosis_mut, "([[:space:]][0-9][\\.]?[[:space:]](?!o'clock))|Gross Description"),
            reporttype == "post_07" ~ str_split(diagnosis_mut, "(\\n.{1,5}[\\.]?[1-9][\\.|[[:space:]]](?!o'clock))|The case material was reviewed and the report verified by")
          )
        ) %>%
        mutate(month_year = format(as.POSIXct(date), format = "%Y-%m")) %>%
        filter(!str_detect(diagnosis_mut, "RESIDENT"))
      
      # Create pathologist deidentification
      pathologist_df <- tot %>% 
        select(pathologist) %>% 
        distinct() %>%
        mutate(path_deid = paste("pathologist-", row.names(.), sep = ""))
      
      tot <- left_join(tot, pathologist_df)
      tot
    }
  ),
  
  # Expand to individual diagnoses
  tar_target(
    diagnosis_level_data,
    create_diagnosis_level_data(combined_data)
  ),
  
  # Text corpus preparation for clustering
  tar_target(
    corpus_data_2006,
    create_corpus_data(diagnosis_level_data, 2006)
  ),
  
  tar_target(
    corpus_data_2015,
    create_corpus_data(diagnosis_level_data, 2015)
  ),
  
  
  tar_target(
    clustering_2006,
    perform_clustering(diagnosis_level_data, 2006)
  ),

  tar_target(
    cluster_plot_2006,
    create_cluster_plot(clustering_2006, 2006, "outputs/plots/figure_3a.png")
  ),
    
  tar_target(
    clustering_2015,
    perform_clustering(diagnosis_level_data, 2015)
  ),
  
  tar_target(
    cluster_plot_2015,
    create_cluster_plot(clustering_2015, 2015, "outputs/plots/figure_3b.png")
  ),
  
  tar_target(
    hpv_analysis,
    analyze_hpv(diagnosis_level_data)
  ),
  
  tar_target(
    p16_analysis,
    analyze_p16(combined_data)
  ),
  
  # CIN classification
  tar_target(
    cin_analysis,
    analyze_cin(diagnosis_level_data)
  ),
  

  
  tar_target(
    pathologist_volume_by_year,
    {
      vol_stats <- combined_data %>%
        count(year, pathologist) %>%
        arrange(year, desc(n))
      
      vol_stats
    }
  ),
  
  # Statistical Testing Targets
  tar_target(
    hpv_chi_square_test,
    perform_chi_square_test(hpv_analysis, year, "hpv")
  ),
  
  tar_target(
    cin_chi_square_test,
    perform_chi_square_test(cin_analysis, year, "dysplasia")
  ),
  
  tar_target(
    pathologist_volume_plot,
    plot_pathologist_volume(pathologist_volume_by_year, "outputs/plots/pathologist_volume_comparison.png")
  ),
  
  tar_target(
    hpv_trends_plot,
    plot_hpv_trends(hpv_analysis, "outputs/plots/hpv_trends.png")
  ),
  
  tar_target(
    hpv_yr_plot_1,
    plot_hpv_by_year(hpv_analysis, "outputs/plots/figure_4.png")
  ),
  
  tar_target(
    p16_yr_plot,
    plot_p16_by_year(p16_analysis, "outputs/plots/figure_5.png")
  ),
  
  tar_target(
    cin_year_per_color_plot,
    plot_cin_by_year(cin_analysis, "outputs/plots/figure_6.png")
  ),
  
  tar_target(
    p16_dysplasia_combined,
    combine_p16_dysplasia(cin_analysis, p16_analysis)
  ),
  
  tar_target(
    p16_prop_panel_plot,
    plot_p16_dysplasia_panel(p16_dysplasia_combined, "outputs/plots/figure_7.png")
  ),
  
  tar_target(
    cin_distribution_plot,
    plot_cin_distribution(cin_analysis, "outputs/plots/cin_distribution.png")
  ),
  
  tar_target(
    pathologist_streamgraph,
    create_pathologist_streamgraph(pathologist_volume_by_year)
  ),
  

  
  tar_target(
    export_hpv_analysis,
    export_hpv_summary(hpv_analysis, "outputs/data/hpv_analysis.csv"),
    format = "file"
  ),
  
  tar_target(
    export_cin_analysis,
    export_cin_summary(cin_analysis, "outputs/data/cin_analysis.csv"),
    format = "file"
  ),
  
  tar_target(
    abstract_statistics_report,
    create_abstract_statistics_report(
      hpv_analysis, 
      p16_analysis, 
      cin_analysis,
      "outputs/reports/abstract_statistics.txt",
      "outputs/reports/abstract_statistics_summary.csv"
    )
  ),
  tar_target(
    results_section_report,
    create_results_section_report(
      hpv_analysis,
      p16_analysis,
      cin_analysis,
      "outputs/reports/results_section.txt",
      "outputs/reports/results_section_summary.csv"
    )
  ),
  
  tar_target(
    export_statistical_tests,
    export_statistical_test_summary(
      hpv_chi_square_test,
      cin_chi_square_test,
      "outputs/data/statistical_test_results.csv"
    ),
    format = "file"
  )
)