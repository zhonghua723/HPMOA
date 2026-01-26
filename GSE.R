
library(readxl)
GEO_190639 <- read_xlsx('~/GSE190639_.xlsx',sheet = 2,col_names = F) |>
  as.data.frame()


library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(purrr)
library(tibble)


all_sample_ids <- GEO_190639 |>
  dplyr::slice(6) |>
  dplyr::select(-`...1`) |>
  unlist(use.names = FALSE)

group_info_strings <- GEO_190639 |>
  dplyr::slice(1) |>
  dplyr::select(-`...1`) |>
  unlist(use.names = FALSE)

sample_to_group_mapping <- tibble::tibble(
  sample_id = all_sample_ids,
  group_raw = group_info_strings
) |>
  dplyr::mutate(group = stringr::str_remove(group_raw, "group: "))

healthy_sample_ids <- sample_to_group_mapping |>
  dplyr::filter(group == "healthy") |>
  dplyr::pull(sample_id)

if (length(healthy_sample_ids) == 0) {
  stop("No samples found in the 'healthy' group.")
}

birth_weight_all_samples_raw <- GEO_190639 |>
  dplyr::slice(4) |>
  dplyr::select(-`...1`) |>
  unlist(use.names = FALSE)

birth_weight_data <- tibble::tibble(
  sample_id = all_sample_ids,
  birth_weight_raw_str = birth_weight_all_samples_raw
) |>
  dplyr::filter(sample_id %in% healthy_sample_ids) |>
  dplyr::mutate(
    birth_weight = readr::parse_number(stringr::str_extract(birth_weight_raw_str, "(?<=birth weight \\(grams\\): ).*"))
  ) |>
  dplyr::select(sample_id, birth_weight) |>
  dplyr::filter(!is.na(birth_weight))

gene_identifiers_col <- GEO_190639 |>
  dplyr::slice(7:dplyr::n()) |>
  dplyr::pull(`...1`)

expression_values_all_samples <- GEO_190639 |>
  dplyr::slice(7:dplyr::n()) |>
  dplyr::select(-`...1`)

colnames(expression_values_all_samples) <- all_sample_ids

expression_values_healthy_samples <- expression_values_all_samples |>
  dplyr::select(dplyr::all_of(healthy_sample_ids))

gene_expression_long <- dplyr::bind_cols(
  tibble::tibble(gene_raw = gene_identifiers_col),
  expression_values_healthy_samples
) |>
  dplyr::mutate(gene = stringr::str_remove(gene_raw, "-mRNA$")) |>
  dplyr::select(gene, dplyr::all_of(healthy_sample_ids)) |>
  tidyr::pivot_longer(
    cols = -gene,
    names_to = "sample_id",
    values_to = "expression_value_char"
  ) |>
  dplyr::mutate(
    expression_value = readr::parse_number(expression_value_char)
  ) |>
  dplyr::select(gene, sample_id, expression_value) |>
  dplyr::filter(!is.na(expression_value))

merged_data <- gene_expression_long |>
  dplyr::inner_join(birth_weight_data, by = "sample_id")

correlation_results_GEO190639 <- merged_data |>
  dplyr::group_by(gene) |>
  dplyr::summarise(
    n_pairs = dplyr::n(),
    correlation_test_result = if (dplyr::n() >= 3) {
      list(stats::cor.test(expression_value, birth_weight, method = "pearson"))
    } else {
      list(NULL)
    },
    .groups = "drop"
  ) |>
  dplyr::filter(!purrr::map_lgl(correlation_test_result, is.null)) |>
  dplyr::mutate(
    correlation_coefficient = purrr::map_dbl(correlation_test_result, ~ .x$estimate),
    p_value = purrr::map_dbl(correlation_test_result, ~ .x$p.value)
  ) |>
  dplyr::select(gene, n_pairs, correlation_coefficient, p_value) |>
  dplyr::arrange(p_value)

write_csv(correlation_results_GEO190639,'GEO190639_birthweight.csv')



all_sample_ids <- GEO_190639 |>
  dplyr::slice(6) |>
  dplyr::select(-`...1`) |>
  unlist(use.names = FALSE)

group_info_strings <- GEO_190639 |>
  dplyr::slice(1) |>
  dplyr::select(-`...1`) |>
  unlist(use.names = FALSE)

sample_to_group_mapping <- tibble::tibble(
  sample_id = all_sample_ids,
  group_raw = group_info_strings
) |>
  dplyr::mutate(group = stringr::str_remove(group_raw, "group: "))

healthy_sample_ids <- sample_to_group_mapping |>
  dplyr::filter(group == "healthy") |>
  dplyr::pull(sample_id)

if (length(healthy_sample_ids) == 0) {
  stop("No samples found in the 'healthy' group.")
}

GA_all_samples_raw <- GEO_190639 |>
  dplyr::slice(2) |>
  dplyr::select(-`...1`) |>
  unlist(use.names = FALSE)

GA_data <- tibble::tibble(
  sample_id = all_sample_ids,
  GA_raw_str = GA_all_samples_raw
) |>
  dplyr::filter(sample_id %in% healthy_sample_ids) |>
  dplyr::mutate(
    GA = readr::parse_number(stringr::str_extract(GA_raw_str, "(?<=gestational age \\(weeks\\): ).*"))
  ) |>
  dplyr::select(sample_id, GA) |>
  dplyr::filter(!is.na(GA))

gene_identifiers_col <- GEO_190639 |>
  dplyr::slice(7:dplyr::n()) |>
  dplyr::pull(`...1`)

expression_values_all_samples <- GEO_190639 |>
  dplyr::slice(7:dplyr::n()) |>
  dplyr::select(-`...1`)

colnames(expression_values_all_samples) <- all_sample_ids

expression_values_healthy_samples <- expression_values_all_samples |>
  dplyr::select(dplyr::all_of(healthy_sample_ids))

gene_expression_long <- dplyr::bind_cols(
  tibble::tibble(gene_raw = gene_identifiers_col),
  expression_values_healthy_samples
) |>
  dplyr::mutate(gene = stringr::str_remove(gene_raw, "-mRNA$")) |>
  dplyr::select(gene, dplyr::all_of(healthy_sample_ids)) |>
  tidyr::pivot_longer(
    cols = -gene,
    names_to = "sample_id",
    values_to = "expression_value_char"
  ) |>
  dplyr::mutate(
    expression_value = readr::parse_number(expression_value_char)
  ) |>
  dplyr::select(gene, sample_id, expression_value) |>
  dplyr::filter(!is.na(expression_value))

merged_data <- gene_expression_long |>
  dplyr::inner_join(GA_data, by = "sample_id")

correlation_results_GEO190639_GA <- merged_data |>
  dplyr::group_by(gene) |>
  dplyr::summarise(
    n_pairs = dplyr::n(),
    correlation_test_result = if (dplyr::n() >= 3) {
      list(stats::cor.test(expression_value, GA, method = "pearson"))
    } else {
      list(NULL)
    },
    .groups = "drop"
  ) |>
  dplyr::filter(!purrr::map_lgl(correlation_test_result, is.null)) |>
  dplyr::mutate(
    correlation_coefficient = purrr::map_dbl(correlation_test_result, ~ .x$estimate),
    p_value = purrr::map_dbl(correlation_test_result, ~ .x$p.value)
  ) |>
  dplyr::select(gene, n_pairs, correlation_coefficient, p_value) |>
  dplyr::arrange(p_value)



GSE190971_data <- openxlsx::read.xlsx('~/GSE190971.xlsx',sheet = 2)
data_tbl <- GSE190971_data |>
  tibble::as_tibble()

all_sample_ids <- colnames(data_tbl)[-1]

disease_state_info_row <- data_tbl |>
  dplyr::slice(1) |>
  dplyr::select(dplyr::all_of(all_sample_ids))

sample_disease_states <- disease_state_info_row |>
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "sample_id",
    values_to = "disease_state_raw"
  ) |>
  dplyr::mutate(
    disease_state = stringr::str_remove(disease_state_raw, "disease state \\(normal or pe\\): ")
  )

np_sample_ids <- sample_disease_states |>
  dplyr::filter(disease_state %in% c("NP")) |>
  dplyr::pull(sample_id)

if (length(np_sample_ids) == 0) {
  stop("No samples found with disease state 'NP'.")
}

birth_weight_info_row_np <- data_tbl |>
  dplyr::slice(16) |>
  dplyr::select(dplyr::all_of(np_sample_ids))

birth_weight_long_np <- birth_weight_info_row_np |>
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "sample_id",
    values_to = "birth_weight_raw"
  ) |>
  dplyr::mutate(
    birth_weight = readr::parse_number(stringr::str_remove(birth_weight_raw, "birth weight: "))
  ) |>
  dplyr::select(sample_id, birth_weight) |>
  dplyr::filter(!is.na(birth_weight))

gene_expression_data_np <- data_tbl |>
  dplyr::slice(18:dplyr::n()) |>
  dplyr::rename(gene = `!Sample_title`) |>
  dplyr::select(gene, dplyr::all_of(np_sample_ids))

gene_expression_long_np <- gene_expression_data_np |>
  tidyr::pivot_longer(
    cols = -gene,
    names_to = "sample_id",
    values_to = "expression_value_char"
  ) |>
  dplyr::mutate(
    expression_value = readr::parse_number(expression_value_char)
  ) |>
  dplyr::select(gene, sample_id, expression_value) |>
  dplyr::filter(!is.na(expression_value))

merged_data_np <- gene_expression_long_np |>
  dplyr::inner_join(birth_weight_long_np, by = "sample_id")

correlation_results_GSE190971 <- merged_data_np |>
  dplyr::group_by(gene) |>
  dplyr::summarise(
    n_pairs = dplyr::n(),
    correlation_test_result = {
      if (dplyr::n() >= 3 &&
          stats::sd(expression_value) > 0 &&
          stats::sd(birth_weight) > 0) {
        list(stats::cor.test(expression_value, birth_weight, method = "spearman"))
      } else {
        list(NULL)
      }
    },
    .groups = "drop"
  ) |>
  dplyr::filter(!purrr::map_lgl(correlation_test_result, is.null)) |>
  dplyr::mutate(
    correlation_coefficient = purrr::map_dbl(correlation_test_result, ~ .x$estimate),
    p_value = purrr::map_dbl(correlation_test_result, ~ .x$p.value)
  ) |>
  dplyr::select(gene, n_pairs, correlation_coefficient, p_value) |>
  dplyr::arrange(p_value)


GA_info_row_np <- data_tbl |>
  dplyr::slice(14) |>
  dplyr::select(dplyr::all_of(np_sample_ids))

GA_long_np <- GA_info_row_np |>
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "sample_id",
    values_to = "GA_raw"
  ) |>
  dplyr::mutate(
    ga_value = stringr::str_remove(GA_raw, "gestational age at delivery: ")
  ) |>
  tidyr::separate_wider_delim(
    ga_value,
    delim = "+",
    names = c("weeks", "days"),
    too_few = "align_start"
  ) |>
  dplyr::mutate(
    weeks = readr::parse_double(weeks),
    days = readr::parse_double(days) |> tidyr::replace_na(0),
    GA_days = weeks * 7 + days
  ) |>
  dplyr::select(sample_id, GA_days)

gene_expression_data_np <- data_tbl |>
  dplyr::slice(18:dplyr::n()) |>
  dplyr::rename(gene = `!Sample_title`) |>
  dplyr::select(gene, dplyr::all_of(np_sample_ids))

gene_expression_long_np <- gene_expression_data_np |>
  tidyr::pivot_longer(
    cols = -gene,
    names_to = "sample_id",
    values_to = "expression_value_char"
  ) |>
  dplyr::mutate(
    expression_value = readr::parse_number(expression_value_char)
  ) |>
  dplyr::select(gene, sample_id, expression_value) |>
  dplyr::filter(!is.na(expression_value))

merged_data_np <- gene_expression_long_np |>
  dplyr::inner_join(GA_long_np, by = "sample_id")

correlation_results_GSE190971_GA <- merged_data_np |>
  dplyr::group_by(gene) |>
  dplyr::summarise(
    n_pairs = dplyr::n(),
    correlation_test_result = {
      if (dplyr::n() >= 3 &&
          stats::sd(expression_value) > 0 &&
          stats::sd(GA_days) > 0) {
        list(stats::cor.test(expression_value, GA_days, method = "pearson"))
      } else {
        list(NULL)
      }
    },
    .groups = "drop"
  ) |>
  dplyr::filter(!purrr::map_lgl(correlation_test_result, is.null)) |>
  dplyr::mutate(
    correlation_coefficient = purrr::map_dbl(correlation_test_result, ~ .x$estimate),
    p_value = purrr::map_dbl(correlation_test_result, ~ .x$p.value)
  ) |>
  dplyr::select(gene, n_pairs, correlation_coefficient, p_value) |>
  dplyr::arrange(p_value)




GSE98224 <- openxlsx::read.xlsx('~/GSE98224.xlsx',sheet = 2)
data_tbl <- GSE98224 |>
  tibble::as_tibble()

all_sample_ids <- colnames(data_tbl)[-1]

diagnosis_info_row <- data_tbl |>
  dplyr::slice(2)

sample_diagnoses <- diagnosis_info_row |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(all_sample_ids),
    names_to = "sample_id",
    values_to = "diagnosis_raw"
  ) |>
  dplyr::mutate(
    diagnosis = stringr::str_remove(diagnosis_raw, "diagnosis: ")
  )

non_pe_sample_ids <- sample_diagnoses |>
  dplyr::filter(diagnosis == "non-PE") |>
  dplyr::pull(sample_id)

if (length(non_pe_sample_ids) == 0) {
  stop("No samples found with diagnosis 'non-PE'.")
}

z_score_info_row_non_pe <- data_tbl |>
  dplyr::slice(22) |>
  dplyr::select(dplyr::all_of(non_pe_sample_ids))

z_score_long_non_pe <- z_score_info_row_non_pe |>
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "sample_id",
    values_to = "z_score_raw"
  ) |>
  dplyr::mutate(
    z_score = readr::parse_number(stringr::str_remove(z_score_raw, "newborn weight z-score: "))
  ) |>
  dplyr::select(sample_id, z_score) |>
  dplyr::filter(!is.na(z_score))

gene_expression_data_non_pe <- data_tbl |>
  dplyr::slice(30:dplyr::n()) |>
  dplyr::rename(gene = `!Sample_title`) |>
  dplyr::select(gene, dplyr::all_of(non_pe_sample_ids))

gene_expression_long_non_pe <- gene_expression_data_non_pe |>
  tidyr::pivot_longer(
    cols = -gene,
    names_to = "sample_id",
    values_to = "expression_value_char"
  ) |>
  dplyr::mutate(
    expression_value = readr::parse_number(expression_value_char)
  ) |>
  dplyr::select(gene, sample_id, expression_value) |>
  dplyr::filter(!is.na(expression_value))

merged_data_non_pe <- gene_expression_long_non_pe |>
  dplyr::inner_join(z_score_long_non_pe, by = "sample_id")


correlation_results_GSE98224 <- merged_data_non_pe |>
  dplyr::group_by(gene) |>
  dplyr::summarise(
    n_pairs = dplyr::n(),
    correlation_test_result = {
      if (dplyr::n() >= 3 &&
          stats::sd(expression_value) > 0 &&
          stats::sd(z_score) > 0) {
        list(stats::cor.test(expression_value, z_score, method = "pearson"))
      } else {
        list(NULL)
      }
    },
    .groups = "drop"
  ) |>
  dplyr::filter(!purrr::map_lgl(correlation_test_result, is.null)) |>
  dplyr::mutate(
    correlation_coefficient = purrr::map_dbl(correlation_test_result, ~ .x$estimate),
    p_value = purrr::map_dbl(correlation_test_result, ~ .x$p.value)
  ) |>
  dplyr::select(gene, n_pairs, correlation_coefficient, p_value) |>
  dplyr::arrange(p_value)


GSE98224_gene_map <- openxlsx::read.xlsx('~/GPL6244-17930.xlsx')
GSE98224_gene_map2 <- GSE98224_gene_map |>
  dplyr::select(ID,Gene) |>
  na.omit() |>
  mutate(ID = as.character(ID))

correlation_results_GSE98224_genename <- left_join(correlation_results_GSE98224,GSE98224_gene_map2, by = c('gene'='ID'),keep = T) |>
  mutate(Gene = trimws(Gene))



GA_week_info_row_non_pe <- data_tbl |>
  dplyr::slice(19) |>
  dplyr::select(dplyr::all_of(non_pe_sample_ids))
GA_day_info_row_non_pe <- data_tbl |>
  dplyr::slice(20) |>
  dplyr::select(dplyr::all_of(non_pe_sample_ids))

ga_weeks_long <- GA_week_info_row_non_pe |>
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "sample_id",
    values_to = "week_raw"
  ) |>
  dplyr::mutate(weeks = readr::parse_number(week_raw)) |>
  dplyr::select(sample_id, weeks)

ga_days_long <- GA_day_info_row_non_pe |>
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "sample_id",
    values_to = "day_raw"
  ) |>
  dplyr::mutate(days = readr::parse_number(day_raw)) |>
  dplyr::select(sample_id, days)

GA_long_non_pe <- ga_weeks_long |>
  dplyr::left_join(ga_days_long, by = "sample_id") |>
  dplyr::mutate(
    GA_days = (weeks * 7) + days
  ) |>
  dplyr::select(sample_id, GA_days)

gene_expression_data_non_pe <- data_tbl |>
  dplyr::slice(30:dplyr::n()) |>
  dplyr::rename(gene = `!Sample_title`) |>
  dplyr::select(gene, dplyr::all_of(non_pe_sample_ids))

gene_expression_long_non_pe <- gene_expression_data_non_pe |>
  tidyr::pivot_longer(
    cols = -gene,
    names_to = "sample_id",
    values_to = "expression_value_char"
  ) |>
  dplyr::mutate(
    expression_value = readr::parse_number(expression_value_char)
  ) |>
  dplyr::select(gene, sample_id, expression_value) |>
  dplyr::filter(!is.na(expression_value))

merged_data_non_pe <- gene_expression_long_non_pe |>
  dplyr::inner_join(GA_long_non_pe, by = "sample_id")


correlation_results_GSE98224_GA <- merged_data_non_pe |>
  dplyr::group_by(gene) |>
  dplyr::summarise(
    n_pairs = dplyr::n(),
    correlation_test_result = {
      if (dplyr::n() >= 3 &&
          stats::sd(expression_value) > 0 &&
          stats::sd(GA_days) > 0) {
        list(stats::cor.test(expression_value, GA_days, method = "pearson"))
      } else {
        list(NULL)
      }
    },
    .groups = "drop"
  ) |>
  dplyr::filter(!purrr::map_lgl(correlation_test_result, is.null)) |>
  dplyr::mutate(
    correlation_coefficient = purrr::map_dbl(correlation_test_result, ~ .x$estimate),
    p_value = purrr::map_dbl(correlation_test_result, ~ .x$p.value)
  ) |>
  dplyr::select(gene, n_pairs, correlation_coefficient, p_value) |>
  dplyr::arrange(p_value)


GSE98224_gene_map <- read.xlsx('~/GPL6244-17930.xlsx')
GSE98224_gene_map2 <- GSE98224_gene_map |>
  dplyr::select(ID,Gene) |>
  na.omit() |>
  mutate(ID = as.character(ID))

correlation_results_GSE98224_GA_genename <- left_join(correlation_results_GSE98224_GA,GSE98224_gene_map2, by = c('gene'='ID'),keep = T) |>
  mutate(Gene = trimws(Gene))
write_csv(correlation_results_GSE98224_GA_genename,'GSE98224_GA.csv')


library(openxlsx)
GSE100415_df <- read.xlsx('~/GSE100415.xlsx',sheet = 2)

library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(purrr)
library(tibble)


data_tbl <- GSE100415_df |>
  tibble::as_tibble()

sample_ids <- colnames(data_tbl)[-1]

z_score_info_row <- data_tbl |>
  dplyr::slice(15) |>
  dplyr::select(dplyr::all_of(sample_ids))

z_score_long <- z_score_info_row |>
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "sample_id",
    values_to = "z_score_raw"
  ) |>
  dplyr::mutate(
    z_score = readr::parse_number(stringr::str_remove(z_score_raw, "newborn weight z-score: "))
  ) |>
  dplyr::select(sample_id, z_score) |>
  dplyr::filter(!is.na(z_score))

gene_expression_data <- data_tbl |>
  dplyr::slice(26:dplyr::n()) |>
  dplyr::rename(gene_symbol = `!Sample_title`) |>
  dplyr::select(gene_symbol, dplyr::all_of(sample_ids))

gene_expression_long <- gene_expression_data |>
  tidyr::pivot_longer(
    cols = -gene_symbol,
    names_to = "sample_id",
    values_to = "expression_value_char"
  ) |>
  dplyr::mutate(
    expression_value = readr::parse_number(expression_value_char)
  ) |>
  dplyr::select(gene_symbol, sample_id, expression_value) |>
  dplyr::filter(!is.na(expression_value))

merged_data <- gene_expression_long |>
  dplyr::inner_join(z_score_long, by = "sample_id")

correlation_results_GSE100415 <- merged_data |>
  dplyr::group_by(gene_symbol) |>
  dplyr::summarise(
    n_pairs = dplyr::n(),
    correlation_test_result = {
      if (dplyr::n() >= 3 &&
          stats::sd(expression_value) > 0 &&
          stats::sd(z_score) > 0) {
        list(stats::cor.test(expression_value, z_score, method = "pearson"))
      } else {
        list(NULL)
      }
    },
    .groups = "drop"
  ) |>
  dplyr::filter(!purrr::map_lgl(correlation_test_result, is.null)) |>
  dplyr::mutate(
    correlation_coefficient = purrr::map_dbl(correlation_test_result, ~ .x$estimate),
    p_value = purrr::map_dbl(correlation_test_result, ~ .x$p.value)
  ) |>
  dplyr::select(gene_symbol, n_pairs, correlation_coefficient, p_value) |>
  dplyr::arrange(p_value)


write_csv(correlation_results_GSE100415,'GSE100415_birthweight.csv')




GA_weeks_info_row <- data_tbl |>
  dplyr::slice(10) |>
  dplyr::select(dplyr::all_of(sample_ids))

GA_days_info_row <- data_tbl |>
  dplyr::slice(11) |>
  dplyr::select(dplyr::all_of(sample_ids))

ga_weeks_long <- GA_weeks_info_row |>
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "sample_id",
    values_to = "week_raw"
  ) |>
  dplyr::mutate(weeks = readr::parse_number(week_raw)) |>
  dplyr::select(sample_id, weeks)

ga_days_long <- GA_days_info_row |>
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "sample_id",
    values_to = "day_raw"
  ) |>
  dplyr::mutate(days = readr::parse_number(day_raw)) |>
  dplyr::select(sample_id, days)

GA_info <- ga_weeks_long |>
  dplyr::left_join(ga_days_long, by = "sample_id") |>
  dplyr::mutate(
    GA_days = (weeks * 7) + days
  ) |>
  dplyr::select(sample_id, GA_days)

gene_expression_data <- data_tbl |>
  dplyr::slice(26:dplyr::n()) |>
  dplyr::rename(gene_symbol = `!Sample_title`) |>
  dplyr::select(gene_symbol, dplyr::all_of(sample_ids))

gene_expression_long <- gene_expression_data |>
  tidyr::pivot_longer(
    cols = -gene_symbol,
    names_to = "sample_id",
    values_to = "expression_value_char"
  ) |>
  dplyr::mutate(
    expression_value = readr::parse_number(expression_value_char)
  ) |>
  dplyr::select(gene_symbol, sample_id, expression_value) |>
  dplyr::filter(!is.na(expression_value))

merged_data <- gene_expression_long |>
  dplyr::inner_join(GA_info, by = "sample_id")

correlation_results_GSE100415_GA <- merged_data |>
  dplyr::group_by(gene_symbol) |>
  dplyr::summarise(
    n_pairs = dplyr::n(),
    correlation_test_result = {
      if (dplyr::n() >= 3 &&
          stats::sd(expression_value) > 0 &&
          stats::sd(GA_days) > 0) {
        list(stats::cor.test(expression_value, GA_days, method = "pearson"))
      } else {
        list(NULL)
      }
    },
    .groups = "drop"
  ) |>
  dplyr::filter(!purrr::map_lgl(correlation_test_result, is.null)) |>
  dplyr::mutate(
    correlation_coefficient = purrr::map_dbl(correlation_test_result, ~ .x$estimate),
    p_value = purrr::map_dbl(correlation_test_result, ~ .x$p.value)
  ) |>
  dplyr::select(gene_symbol, n_pairs, correlation_coefficient, p_value) |>
  dplyr::arrange(p_value)


GSE75010_df <- read.xlsx('GSE75010.xlsx')
data_tbl <- GSE75010_df |>
  tibble::as_tibble()

all_sample_ids <- colnames(data_tbl)[-1]

diagnosis_info_row <- data_tbl |>
  dplyr::slice(2)

sample_diagnoses <- diagnosis_info_row |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(all_sample_ids),
    names_to = "sample_id",
    values_to = "diagnosis_raw"
  ) |>
  dplyr::mutate(
    diagnosis = stringr::str_remove(diagnosis_raw, "diagnosis: ")
  )

non_pe_sample_ids <- sample_diagnoses |>
  dplyr::filter(diagnosis %in% c("non-PE")) |>
  dplyr::pull(sample_id)

if (length(non_pe_sample_ids) == 0) {
  stop("No samples found with diagnosis 'non-PE'.")
}

z_score_info_row_non_pe <- data_tbl |>
  dplyr::slice(23) |>
  dplyr::select(dplyr::all_of(non_pe_sample_ids))

z_score_long_non_pe <- z_score_info_row_non_pe |>
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "sample_id",
    values_to = "z_score_raw"
  ) |>
  dplyr::mutate(
    z_score = readr::parse_number(stringr::str_remove(z_score_raw, "newborn weight z-score: "))
  ) |>
  dplyr::select(sample_id, z_score) |>
  dplyr::filter(!is.na(z_score))

gene_expression_data_non_pe <- data_tbl |>
  dplyr::slice(31:dplyr::n()) |>
  dplyr::rename(gene = `!Sample_title`) |>
  dplyr::select(gene, dplyr::all_of(non_pe_sample_ids))

gene_expression_long_non_pe <- gene_expression_data_non_pe |>
  tidyr::pivot_longer(
    cols = -gene,
    names_to = "sample_id",
    values_to = "expression_value_char"
  ) |>
  dplyr::mutate(
    expression_value = readr::parse_number(expression_value_char)
  ) |>
  dplyr::select(gene, sample_id, expression_value) |>
  dplyr::filter(!is.na(expression_value))

merged_data_non_pe <- gene_expression_long_non_pe |>
  dplyr::inner_join(z_score_long_non_pe, by = "sample_id")

correlation_results_GSE75010 <- merged_data_non_pe |>
  dplyr::group_by(gene) |>
  dplyr::summarise(
    n_pairs = dplyr::n(),
    correlation_test_result = {
      if (dplyr::n() >= 3 &&
          stats::sd(expression_value) > 0 &&
          stats::sd(z_score) > 0) {
        list(stats::cor.test(expression_value, z_score, method = "pearson"))
      } else {
        list(NULL)
      }
    },
    .groups = "drop"
  ) |>
  dplyr::filter(!purrr::map_lgl(correlation_test_result, is.null)) |>
  dplyr::mutate(
    correlation_coefficient = purrr::map_dbl(correlation_test_result, ~ .x$estimate),
    p_value = purrr::map_dbl(correlation_test_result, ~ .x$p.value)
  ) |>
  dplyr::select(gene, n_pairs, correlation_coefficient, p_value) |>
  dplyr::arrange(p_value)

GSE75010_genemap <- read.xlsx('GPL6244-17930.xlsx')
GSE75010_genemap <- GSE75010_genemap |> dplyr::select(ID,gene_assignment)
GSE75010_genemap_cleaned <- GSE75010_genemap |>
  dplyr::mutate(
    gene_symbol = stringr::str_split(gene_assignment, "//") |>
      purrr::map_chr(~ stringr::str_trim(.x[2]))
  )

GSE75010_genemap_cleaned2 <- GSE75010_genemap_cleaned |>
  dplyr::select(ID, gene_symbol) |>
  na.omit() |>
  mutate(ID = as.character(ID))

correlation_results_GSE75010_withgenename <- merge(correlation_results_GSE75010,
                                                  GSE75010_genemap_cleaned2,
                                                  by.x = 'gene',
                                                  by.y = 'ID',
                                                  all.x = T)
correlation_results_GSE75010_withgenename_noNA <- correlation_results_GSE75010_withgenename |>
  na.omit()



GA_weeks_info_row_non_pe <- data_tbl |>
  dplyr::slice(20) |>
  dplyr::select(dplyr::all_of(non_pe_sample_ids))

GA_days_info_row_non_pe <- data_tbl |>
  dplyr::slice(21) |>
  dplyr::select(dplyr::all_of(non_pe_sample_ids))

ga_weeks_long <- GA_weeks_info_row_non_pe |>
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "sample_id",
    values_to = "week_raw"
  ) |>
  dplyr::mutate(weeks = readr::parse_number(week_raw)) |>
  dplyr::select(sample_id, weeks)

ga_days_long <- GA_days_info_row_non_pe |>
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "sample_id",
    values_to = "day_raw"
  ) |>
  dplyr::mutate(days = readr::parse_number(day_raw)) |>
  dplyr::select(sample_id, days)

GA_info_non_pe <- ga_weeks_long |>
  dplyr::left_join(ga_days_long, by = "sample_id") |>
  dplyr::mutate(
    GA_days = (weeks * 7) + days
  ) |>
  dplyr::select(sample_id, GA_days)

gene_expression_data_non_pe <- data_tbl |>
  dplyr::slice(31:dplyr::n()) |>
  dplyr::rename(gene = `!Sample_title`) |>
  dplyr::select(gene, dplyr::all_of(non_pe_sample_ids))

gene_expression_long_non_pe <- gene_expression_data_non_pe |>
  tidyr::pivot_longer(
    cols = -gene,
    names_to = "sample_id",
    values_to = "expression_value_char"
  ) |>
  dplyr::mutate(
    expression_value = readr::parse_number(expression_value_char)
  ) |>
  dplyr::select(gene, sample_id, expression_value) |>
  dplyr::filter(!is.na(expression_value))

merged_data_non_pe <- gene_expression_long_non_pe |>
  dplyr::inner_join(GA_info_non_pe, by = "sample_id")

correlation_results_GSE75010_GA <- merged_data_non_pe |>
  dplyr::group_by(gene) |>
  dplyr::summarise(
    n_pairs = dplyr::n(),
    correlation_test_result = {
      if (dplyr::n() >= 3 &&
          stats::sd(expression_value) > 0 &&
          stats::sd(GA_days) > 0) {
        list(stats::cor.test(expression_value, GA_days, method = "pearson"))
      } else {
        list(NULL)
      }
    },
    .groups = "drop"
  ) |>
  dplyr::filter(!purrr::map_lgl(correlation_test_result, is.null)) |>
  dplyr::mutate(
    correlation_coefficient = purrr::map_dbl(correlation_test_result, ~ .x$estimate),
    p_value = purrr::map_dbl(correlation_test_result, ~ .x$p.value)
  ) |>
  dplyr::select(gene, n_pairs, correlation_coefficient, p_value) |>
  dplyr::arrange(p_value)

GSE75010_genemap <- read.xlsx('GPL6244-17930.xlsx')
GSE75010_genemap <- GSE75010_genemap |> dplyr::select(ID,gene_assignment)
GSE75010_genemap_cleaned <- GSE75010_genemap |>
  dplyr::mutate(
    gene_symbol = stringr::str_split(gene_assignment, "//") |>
      purrr::map_chr(~ stringr::str_trim(.x[2]))
  )

GSE75010_genemap_cleaned2 <- GSE75010_genemap_cleaned |>
  dplyr::select(ID, gene_symbol) |>
  na.omit() |>
  mutate(ID = as.character(ID))

correlation_results_GSE75010_GA_withgenename <- merge(correlation_results_GSE75010_GA,
                                                   GSE75010_genemap_cleaned2,
                                                   by.x = 'gene',
                                                   by.y = 'ID',
                                                   all.x = T)
correlation_results_GSE75010_GA_withgenename_noNA <- correlation_results_GSE75010_GA_withgenename |>
  na.omit()


GSE30032_df <- read.xlsx('~/GSE30032.xlsx',sheet = 2)
data_tbl <- GSE30032_df |>
  tibble::as_tibble()

all_sample_ids <- colnames(data_tbl)[-1]
placenta_sample_ids <- all_sample_ids |>
  stringr::str_subset("non-smoker")  |>
  stringr::str_subset("\\.placenta$")

if (length(placenta_sample_ids) == 0) {
  stop("No placenta samples were found.")
}

newborn_weight_info_row <- data_tbl |>
  dplyr::slice(10) |>
  dplyr::select(dplyr::all_of(placenta_sample_ids))

newborn_weight_long <- newborn_weight_info_row |>
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "sample_id",
    values_to = "weight_raw"
  ) |>
  dplyr::mutate(
    newborn_weight = readr::parse_number(stringr::str_remove(weight_raw, "newborn weight \\(g\\): "))
  ) |>
  dplyr::select(sample_id, newborn_weight) |>
  dplyr::filter(!is.na(newborn_weight))

gene_expression_data <- data_tbl |>
  dplyr::slice(17:dplyr::n()) |>
  dplyr::rename(gene = `!Sample_title`) |>
  dplyr::select(gene, dplyr::all_of(placenta_sample_ids))

gene_expression_long <- gene_expression_data |>
  tidyr::pivot_longer(
    cols = -gene,
    names_to = "sample_id",
    values_to = "expression_value_char"
  ) |>
  dplyr::mutate(
    expression_value = readr::parse_number(expression_value_char)
  ) |>
  dplyr::select(gene, sample_id, expression_value) |>
  dplyr::filter(!is.na(expression_value))

merged_data <- gene_expression_long |>
  dplyr::inner_join(newborn_weight_long, by = "sample_id")

correlation_results_GSE30032 <- merged_data |>
  dplyr::group_by(gene) |>
  dplyr::summarise(
    n_pairs = dplyr::n(),
    correlation_test_result = {
      if (dplyr::n() >= 3 &&
          stats::sd(expression_value) > 0 &&
          stats::sd(newborn_weight) > 0) {
        list(stats::cor.test(expression_value, newborn_weight, method = "pearson"))
      } else {
        list(NULL)
      }
    },
    .groups = "drop"
  ) |>
  dplyr::filter(!purrr::map_lgl(correlation_test_result, is.null)) |>
  dplyr::mutate(
    correlation_coefficient = purrr::map_dbl(correlation_test_result, ~ .x$estimate),
    p_value = purrr::map_dbl(correlation_test_result, ~ .x$p.value)
  ) |>
  dplyr::select(gene, n_pairs, correlation_coefficient, p_value) |>
  dplyr::arrange(p_value)


GSE30032_genemap <- read.xlsx('~/GSE30032.xlsx')
GSE30032_genemap2 <- GSE30032_genemap |>
  dplyr::select(Probe_Id, Symbol)

correlation_results_GSE30032_withgenename <- merge(correlation_results_GSE30032,
                                                   GSE30032_genemap2,
                                                   by.x = 'gene',
                                                   by.y = 'Probe_Id',all.x = T)

GSE190639_birthweight <- read_csv('GEO190639_birthweight.csv') %>% mutate(GSE = 'GSE190639') |>
  dplyr::select(correlation_coefficient, p_value, gene,GSE)


GSE98224_birthweight <- read_csv('GSE982224_birthweight.csv') %>% mutate(GSE = 'GSE98224') |>
  dplyr::select(correlation_coefficient, p_value, Gene,GSE) |>
  dplyr::rename(gene = Gene)

GSE100415_birthweight <- read_csv('GSE100415_birthweight.csv') %>% mutate(GSE = 'GSE100415') |>
  dplyr::select(correlation_coefficient, p_value, gene_symbol,GSE) |>
  dplyr::rename(gene = gene_symbol)

GSE75010_birthweight <- read_csv('GSE75010_birthweight.csv') %>% mutate(GSE = 'GSE75010') |>
  dplyr::select(correlation_coefficient, p_value, gene_symbol,GSE) |>
  dplyr::rename(gene = gene_symbol)

GSE30032_birthweight <- read_csv('GSE30032_birthweight.csv') %>% mutate(GSE = 'GSE30032') |>
  dplyr::select(correlation_coefficient, p_value, Symbol,GSE) |>
  dplyr::rename(gene = Symbol)

manhat_plot_data <- bind_rows(GSE190639_birthweight,GSE98224_birthweight,GSE100415_birthweight,GSE75010_birthweight,GSE30032_birthweight)
p_m <- manhat_plot_data |>
  dplyr::mutate(
    color_group = dplyr::case_when(
      p_value > 0.05 ~ "#F4F3EE",
      p_value <= 0.05 & correlation_coefficient > 0 ~ "#EF767B",
      p_value <= 0.05 & correlation_coefficient < 0 ~ "#43A3EF"
    )
  ) |>
  ggplot2::ggplot(ggplot2::aes(x = GSE, y = correlation_coefficient, color = color_group)) +
  ggplot2::geom_jitter(width = 0.25, alpha = 0.6, size = 1.5) +
  ggplot2::scale_color_identity(
    name = "Legend",
    guide = "legend",
    labels = c(
      red = "Positive (p <= 0.05)",
      blue = "Negative (p <= 0.05)",
      grey = "No (p > 0.05)"
    )
  ) +
  ggplot2::labs(
    x = "GSE",
    y = "Correlation Coefficient"
  ) +
  ggthemes::theme_base() +
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    legend.position = "top"
  )


GSE190639_GA <- read_csv('GEO190639_GA.csv') %>% mutate(GSE = 'GSE190639') |>
  dplyr::select(correlation_coefficient, p_value, gene,GSE)


GSE98224_GA <- read_csv('GSE98224_GA.csv') %>% mutate(GSE = 'GSE982224') |>
  dplyr::select(correlation_coefficient, p_value, Gene,GSE) |>
  dplyr::rename(gene = Gene)

GSE100415_GA <- read_csv('GSE100415_GA.csv') %>% mutate(GSE = 'GSE100415') |>
  dplyr::select(correlation_coefficient, p_value, gene_symbol,GSE) |>
  dplyr::rename(gene = gene_symbol)

GSE75010_GA <- read_csv('GSE75010_GA.csv') %>% mutate(GSE = 'GSE75010') |>
  dplyr::select(correlation_coefficient, p_value, gene_symbol,GSE) |>
  dplyr::rename(gene = gene_symbol)

manhat_plot_data <- bind_rows(GSE190639_GA,GSE98224_GA,GSE100415_GA,GSE75010_GA)
p_m <- manhat_plot_data |>
  dplyr::mutate(
    color_group = dplyr::case_when(
      p_value > 0.05 ~ "#F4F3EE",
      p_value <= 0.05 & correlation_coefficient > 0 ~ "#EF767B",
      p_value <= 0.05 & correlation_coefficient < 0 ~ "#43A3EF"
    )
  ) |>
  ggplot2::ggplot(ggplot2::aes(x = GSE, y = correlation_coefficient, color = color_group)) +
  ggplot2::geom_jitter(width = 0.25, alpha = 0.6, size = 1.5) +
  ggplot2::scale_color_identity(
    name = "Legend",
    guide = "legend",
    labels = c(
      red = "Positive (p <= 0.05)",
      blue = "Negative (p <= 0.05)",
      grey = "No (p > 0.05)"
    )
  ) +
  ggplot2::labs(
    x = "GSE",
    y = "Correlation Coefficient"
  ) +
  ggthemes::theme_base() +
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    legend.position = "top"
  )
