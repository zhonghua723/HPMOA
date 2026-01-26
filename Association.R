library(tidyverse)
library(readxl)
library(openxlsx)
uniprot_protein_map <- read.xlsx('~/idmapping.xlsx')

uniprot_protein_map2 <- uniprot_protein_map |>
  dplyr::select(From, `Subcellular.location.[CC]`)
clean_subcellular <- function(text) {
  if (is.na(text)) return(NA)

  if (grepl(" Note=", text)) {
    text <- sub(" Note=.*$", "", text)
  }

  text <- gsub("SUBCELLULAR LOCATION:\\s*", "", text)

  text <- gsub("\\[[^\\]]+\\]: ?", "", text, perl=TRUE)


  entries <- strsplit(text, "\\. ")[[1]]

  locations <- c()
  for (entry in entries) {
    if (entry == "") next

    entry <- sub(" \\{.*$", "", entry)

    parts <- strsplit(entry, "; ")[[1]]

    for (part in parts) {
      main_loc <- strsplit(part, ", ")[[1]]

      locations <- c(locations, trimws(main_loc))
    }
  }

  unique_locations <- unique(locations)
  unique_locations <- gsub("\\.$", "", unique_locations)

  return(paste(unique_locations, collapse = "; "))
}


uniprot_protein_map2$clean_subcellular <- sapply(uniprot_protein_map2$`Subcellular.location.[CC]`, clean_subcellular)
uniprot_protein_map3 <- uniprot_protein_map2 |> dplyr::select(From,clean_subcellular)
subcellular_locations <- unlist(strsplit(as.character(uniprot_protein_map3$clean_subcellular), "; "))
subcellular_locations <- subcellular_locations[!is.na(subcellular_locations)]
location_counts <- table(subcellular_locations)
location_counts <- sort(location_counts, decreasing = TRUE)
classify_subcellular_location <- function(location) {
  if(is.na(location) || location == "") return(NA)

  classification_rules <- list(
    "Secreted" = c("Secreted", "secretory vesicle", "secretory vesicle membrane", "secretory vesicle lumen"),

    "Extracellular matrix" = c("Extracellular matrix", "extracellular matrix", "basement membrane", "matrix", "interphotoreceptor matrix"),

    "Extracellular space" = c("Extracellular", "extracellular space", "Extracellular side", "Extracellular vesicle",
                             "Synaptic cleft", "extracellular exosome"),

    "Cell membrane" = c("Cell membrane", "cell surface", "Apical cell membrane", "Basolateral cell membrane",
                        "Lateral cell membrane", "Basal cell membrane", "Apicolateral cell membrane",
                        "Target cell membrane", "membrane protein", "Membrane raft", "Cell surface",
                        "sarcolemma", "caveola", "Single-pass", "Multi-pass", "GPI-anchor", "Lipid-anchor"),

    "Golgi apparatus" = c("Golgi apparatus", "Golgi stack", "Golgi apparatus membrane", "Golgi apparatus lumen",
                         "Golgi outpost", "trans-Golgi network", "cis-Golgi network", "trans-Golgi network membrane",
                         "cis-Golgi network membrane", "Golgi stack membrane"),

    "Melanosome" = c("Melanosome", "Melanosome membrane"),

    "Cytoplasm" = c("Cytoplasm", "Cytoplasmic vesicle", "Cytoplasmic granule", "Cytoplasmic side",
                   "Perikaryon", "Cell tip", "perinuclear region", "Stress granule", "P-body",
                   "Cytoplasmic ribonucleoprotein granule", "cell cortex", "bleb", "uropodium","Cytosol", "cytosol"),

    "Endoplasmic reticulum" = c("Endoplasmic reticulum", "ER ", "Sarcoplasmic reticulum",
                               "Endoplasmic reticulum membrane", "Endoplasmic reticulum lumen",
                               "Rough endoplasmic reticulum", "Rough endoplasmic reticulum membrane",
                               "Rough endoplasmic reticulum lumen", "Smooth endoplasmic reticulum membrane",
                               "Sarcoplasmic reticulum membrane", "Sarcoplasmic reticulum lumen"),

    "Nucleus" = c("Nucleus", "Nuclear", "Chromosome", "nucleolus", "nucleoplasm", "Nucleus speckle",
                 "Nucleus matrix", "Nucleus lamina", "Nucleus inner membrane", "Nucleus outer membrane",
                 "Nucleus membrane", "Nucleus envelope", "nuclear pore complex", "PML body",
                 "Cajal body", "telomere", "centromere", "nuclear body", "nucleolus fibrillar center",
                 "Nucleus intermembrane space", "gem"),

    "Lysosome" = c("Lysosome", "Lysosome membrane", "Lysosome lumen", "Autolysosome", "Autolysosome membrane"),

    "Peroxisome" = c("Peroxisome", "Peroxisome membrane", "Peroxisome matrix", "^Peroxisome$"),

    "Cytoskeleton" = c("Cytoskeleton", "cytoskeleton", "Cell projection", "Midbody", "Cleavage furrow",
                      "Dynein", "spindle", "lamellipodium", "focal adhesion", "dendrite", "axon",
                      "ruffle", "ruffle membrane", "filopodium", "stress fiber", "myofibril", "sarcomere",
                      "Z line", "I band", "M line", "cilium", "stereocilium", "microvillus", "flagellum",
                      "spindle pole", "Midbody ring", "cilium axoneme", "neuron projection", "dendritic spine",
                      "filopodium membrane", "filopodium tip", "stereocilium membrane", "Dynein axonemal particle"),

    "Mitochondrion" = c("Mitochondrion", "Mitochondrial", "Mitochondrion inner membrane", "Mitochondrion outer membrane",
                       "Mitochondrion matrix", "Mitochondrion membrane", "Mitochondrion intermembrane space",
                       "mitochondrion nucleoid"),

    "Centrosome" = c("Centrosome", "centrosome", "microtubule organizing center", "centriole", "centriolar satellite",
                    "cilium basal body", "spindle pole body"),

    "Endosome" = c("Endosome", "Endosome membrane", "Early endosome", "Early endosome membrane",
                  "Late endosome", "Late endosome membrane", "Recycling endosome", "Recycling endosome membrane",
                  "Endosome lumen", "multivesicular body", "multivesicular body membrane")
  )

  for(category in names(classification_rules)) {
    patterns <- classification_rules[[category]]
    for(pattern in patterns) {
      if(grepl(pattern, location, ignore.case = TRUE)) {
        return(category)
      }
    }
  }

  if(grepl("membrane", location, ignore.case = TRUE)) {
    return("Cell membrane")
  }
  if(grepl("vesicle", location, ignore.case = TRUE)) {
    return("Cytoplasm")
  }
  if(grepl("junction", location, ignore.case = TRUE)) {
    return("Cell membrane")
  }
  if(grepl("phagosome", location, ignore.case = TRUE) ||
     grepl("autophagosome", location, ignore.case = TRUE)) {
    return("Cytoplasm")
  }

  return("Other")
}

classify_multiple_locations_all <- function(locations) {
  if(is.na(locations) || locations == "") return(NA)

  location_list <- unlist(strsplit(locations, ";"))
  location_list <- trimws(location_list)

  classifications <- sapply(location_list, classify_subcellular_location)

  unique_classifications <- unique(classifications)

  if(length(unique_classifications) == 0) {
    return("Other")
  } else {
    return(paste(unique_classifications, collapse = "; "))
  }
}

test_location <- uniprot_protein_map3$clean_subcellular[11]
result_all <- classify_multiple_locations_all(test_location)
print(result_all)

uniprot_protein_map3$clean_subcellular_classfied <- lapply(uniprot_protein_map3$clean_subcellular,
                                                           classify_multiple_locations_all)
subcellular_df <- as.data.frame(table(unlist(strsplit(as.character(uniprot_protein_map3$clean_subcellular_classfied), "; "))))
protein_subcelluar <- cbind(uniprot_protein_map2 |> dplyr::select(From,`Subcellular.location.[CC]`),
                            uniprot_protein_map3 |> dplyr::select(From,clean_subcellular_classfied))
protein_subcelluar$clean_subcellular_classfied <- unlist(protein_subcelluar$clean_subcellular_classfied)


rm(list = ls())
gc()
uniprot_protein_map <- read.xlsx('~/idmapping.xlsx')
df1 <- uniprot_protein_map |>
  dplyr::select(From,`Gene.Ontology.(biological.process)`)
df2 <- df1 |>
  dplyr::mutate(
    go_ids = stringr::str_extract_all(`Gene.Ontology.(biological.process)`, "GO:\\d{7}")
  )

library(ontologyIndex)
library(igraph)
library(dplyr)
library(stringr)
library(readr)
library(tidyr)

load("~/GO_find_Ancestors.RData")

find_GO_father <- function(x){
  x1 <- unlist(x)
  results <- process_go_terms(x1, bp_data, mf_data, cc_data, ontology)
  results_bp <- results[["biological_process"]]
  unique_ancestors <- results_bp$Ancestors_1_Name |>
    stringr::str_split(pattern = ";") |>
    base::unlist() |>
    stringr::str_trim() |>
    base::unique()
  ancestors <- stringr::str_c(unique_ancestors, collapse = ";")
  return(ancestors)
}
df3 <- df2 |>
  mutate(GO_ancestor = purrr::map(go_ids,find_GO_father))

terms_to_other <- c(
  "viral process" = "other",
  "rhythmic process" = "other",
  "pigmentation" = "other",
  "locomotion" = "other"
)

df4 <- df3 |>
  dplyr::mutate(
    GO_ancestor = stringr::str_replace_all(GO_ancestor, terms_to_other)
  )
df4$GO_ancestor <- unlist(df4$GO_ancestor)
df4 <- df4 |>
  dplyr::mutate(
    go_ids2 = purrr::map_chr(go_ids, ~ stringr::str_c(.x, collapse = ";"))
  ) |> dplyr::select(From,`Gene.Ontology.(biological.process)`,go_ids2,GO_ancestor)

BP_df <- as.data.frame(table(unlist(strsplit(as.character(df4$GO_ancestor), ";"))))

BP_df_selected <- BP_df |>
  dplyr::filter(Var1 != '')
p_bp <- BP_df_selected |>
  dplyr::mutate(
    Var1 = dplyr::case_when(
      Var1 == "biological process involved in interspecies interaction between organisms" ~ "interspecies interaction\nbetween organisms",
      Var1 == "biological process involved in intraspecies interaction between organisms" ~ "intraspecies interaction\nbetween organisms",
      Var1 == "multicellular organismal process" ~ "multicellular process",
      TRUE ~ Var1
    ),
    Var1 = forcats::fct_reorder(Var1, Freq),
    Var1 = forcats::fct_relevel(Var1, "other")
  ) |>
  ggplot2::ggplot(ggplot2::aes(x = Freq, y = Var1, fill = Freq)) +
  ggplot2::geom_col() +
  ggplot2::geom_text(ggplot2::aes(label = Freq), hjust = 0, nudge_x = 50, size = 4.5) +
  ggplot2::scale_fill_gradient(low = "#f7ddd8", high = "#fbb1a2") +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.01, 0.1))) +
  ggplot2::labs(x = NULL, y = NULL) +
  ggplot2::theme(legend.position = "none") +
  ggthemes::theme_base() +
  ggplot2::theme(axis.text = ggplot2::element_text(size = 16))

p_bp

rm(list = ls())
gc()
uniprot_protein_map <- read.xlsx('~/idmapping.xlsx')
df1 <- uniprot_protein_map |>
  dplyr::select(From,`Gene.Ontology.(molecular.function)`)
df2 <- df1 |>
  dplyr::mutate(
    go_ids = stringr::str_extract_all(`Gene.Ontology.(molecular.function)`, "GO:\\d{7}")
  )

library(ontologyIndex)
library(igraph)
library(dplyr)
library(stringr)
library(readr)
library(tidyr)

load("~/GO_find_Ancestors.RData")

find_GO_father_MF <- function(x){
  x1 <- unlist(x)
  results <- process_go_terms(x1, bp_data, mf_data, cc_data, ontology)
  results_mf <- results[["molecular_function"]]
  unique_ancestors <- results_mf$Ancestors_1_Name |>
    stringr::str_split(pattern = ";") |>
    base::unlist() |>
    stringr::str_trim() |>
    base::unique()
  ancestors <- stringr::str_c(unique_ancestors, collapse = ";")
  return(ancestors)
}
df3 <- df2 |>
  mutate(GO_ancestor = purrr::map(go_ids,find_GO_father_MF))

terms_to_other <- c(
  "protein-containing complex stabilizing activity" = "other",
  "molecular template activity" = "other",
  "fusogenic activity" = "other",
  "protein-containing complex destabilizing activity" = "other",
  'molecular tag activity'= 'other',
  'cargo receptor activity' = 'other',
  'general transcription initiation factor activity' = 'other'
)

df4 <- df3 |>
  dplyr::mutate(
    GO_ancestor = stringr::str_replace_all(GO_ancestor, terms_to_other)
  )
df4$GO_ancestor <- unlist(df4$GO_ancestor)
df4 <- df4 |>
  dplyr::mutate(
    go_ids2 = purrr::map_chr(go_ids, ~ stringr::str_c(.x, collapse = ";"))
  ) |> dplyr::select(From,`Gene.Ontology.(molecular.function)`,go_ids2,GO_ancestor)

MF_df <- as.data.frame(table(unlist(strsplit(as.character(df4$GO_ancestor), ";"))))

MF_df_selected <- MF_df |>
  dplyr::filter(Var1 != '')
p_mf <- MF_df_selected |>
  dplyr::mutate(
    Var1 = dplyr::case_when(
      Var1 == "molecular adaptor activity" ~ "adaptor activity",
      Var1 == "molecular function regulator activity" ~ "regulator activity",
      Var1 == "molecular carrier activity" ~ "carrier activity",
      Var1 == "molecular sequestering activity" ~ "sequestering activity",
      Var1 == "molecular transducer activity" ~ "transducer activity",
      Var1 == "structural molecule activity" ~ "structural activity",
      Var1 == "molecular sequestering activity" ~ "sequestering activity",
      Var1 == "molecular sequestering activity" ~ "sequestering activity",
      TRUE ~ Var1
    ),
    Var1 = forcats::fct_reorder(Var1, Freq),
    Var1 = forcats::fct_relevel(Var1, "other")
  ) |>
  ggplot2::ggplot(ggplot2::aes(x = Freq, y = Var1, fill = Freq)) +
  ggplot2::geom_col() +
  ggplot2::geom_text(ggplot2::aes(label = Freq), hjust = 0, nudge_x = 50, size = 4.5) +
  ggplot2::scale_fill_gradient(low = "#ebf5e9", high = "#9dd193") +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.01, 0.1))) +
  ggplot2::labs(x = NULL, y = NULL) +
  ggplot2::theme(legend.position = "none") +
  ggthemes::theme_base() +
  ggplot2::theme(axis.text = ggplot2::element_text(size = 16))
p_mf

library(tidyverse)
library(gghalves)
library(circlize)
library(writexl)
library(gridExtra)
proteome_dat <- read_csv('cleaned_data.csv')

phenotype_dat <- read_csv('phenotype.csv')

phenotype_dat_birthweight <- phenotype_dat |> dplyr::select(number,birth_weight)

birthweight_cor_dat <- merge(phenotype_dat_birthweight,
                             proteome_dat,
                             by.x = 'number',
                             by.y = 'number_identifier') |>
  na.omit()


phenotype_heatmap_df <- phenotype_dat |>
  dplyr::filter(number %in% birthweight_cor_dat$number) |>
  group_by(number) |>
  slice_head(n = 1) |>
  ungroup() |>
  dplyr::select(number, Age,BMI,Psmoke,Drink,Edu,Income,GA_age)


age_color <- '#f93c56'
bmi_color <- '#7067d0'
psmoke_color <- c("0" = '#91CAE8', "1" = '#F48892')
drink_color <- c("0" = '#a3c785', "1" = '#BDB9B8')
edu_color <- c("1" = '#D8B365', "2" = '#5BB5AC', "3" = '#DE526C')
income_color <- c("1" = '#ECCBAE', "2" = '#0B775E', "3" = '#D69C4E')

df <- data.frame(
  factors = as.character(phenotype_heatmap_df$number),
  x = 1,
  y = 1,
  phenotype_heatmap_df,
  stringsAsFactors = FALSE
) %>%
  arrange(Age) %>%
  mutate(factors = factor(factors, levels = factors),
         Psmoke = as.character(Psmoke),
         Drink = as.character(Drink),
         Edu = as.character(Edu),
         Income = as.character(Income))


circos.clear()
circos.par(
  "track.height" = 0.1,
  start.degree = 90,
  clock.wise = TRUE,
  gap.after = c(rep(0, nrow(df) - 1), 90),
  circle.margin = c(0.05, 0.05, 0.05, 0.05),
  cell.padding = c(0, 0, 0, 0)
)
circos.initialize(factors = df$factors, x = df$x, xlim = c(0.5, 1.5))

plot_circos_track <- function(temp_value, color, ylab, ylim_mult = c(0.8, 1.1), track_height = 0.12) {
  circos.track(
    factors = df$factors,
    y = temp_value,
    ylim = range(temp_value, na.rm = TRUE) * ylim_mult,
    bg.border = "#f2f2f2",
    track.height = track_height,
    panel.fun = function(x, y) {
      name = get.cell.meta.data("sector.index")
      i = get.cell.meta.data("sector.numeric.index")
      xlim = get.cell.meta.data("xlim")
      ylim = get.cell.meta.data("ylim")

      circos.yaxis(
        side = "left",
        at = c(ceiling(0.8*min(temp_value, na.rm = TRUE)),
               round((min(temp_value, na.rm = TRUE) + max(temp_value, na.rm = TRUE)) / 2, 0),
               round(max(temp_value, na.rm = TRUE), 0)),
        sector.index = get.all.sector.index()[1],
        labels.cex = 0.6,
        labels.niceFacing = FALSE
      )

      circos.lines(
        x = mean(xlim, na.rm = TRUE),
        y = temp_value[i],
        type = "h",
        col = color,
        lwd = 0.8
      )

      circos.points(
        x = mean(xlim),
        y = temp_value[i],
        pch = 16,
        cex = 0.4,
        col = color
      )
    }
  )
}

plot_categorical_track <- function(temp_var, colors, track_height = 0.08) {
  temp_var[is.na(temp_var)] <- "Unknown"
  mapped_colors <- colors[temp_var]

  circos.track(
    factors = df$factors,
    y = df$y,
    ylim = c(0, 1),
    bg.border = "#f2f2f2",
    track.height = track_height,
    panel.fun = function(x, y) {
      i = get.cell.meta.data("sector.numeric.index")
      xlim = get.cell.meta.data("xlim")
      ylim = get.cell.meta.data("ylim")

      circos.rect(
        xleft = xlim[1],
        ybottom = ylim[1],
        xright = xlim[2],
        ytop = ylim[2],
        col = mapped_colors[i],
        border = "#f2f2f2",
        lwd = 0.3
      )
    }
  )
}


plot_circos_track(df$Age, age_color, "Age", track_height = 0.15)
plot_circos_track(df$BMI, bmi_color, "BMI", track_height = 0.15)

plot_categorical_track(df$Psmoke, psmoke_color,track_height = 0.1)
plot_categorical_track(df$Edu, edu_color,track_height = 0.1)
plot_categorical_track(df$Income, income_color,track_height = 0.1)



plot_box_dot <- function(data, color, binwidth) {
  ggplot(data.frame(class = "class", value = data), aes(x = class, y = value)) +
    geom_boxplot(outlier.shape = NA, width = 0.2) +
    geom_dotplot(binaxis = "y",
                 color = color,
                 fill = color,
                 dotsize = 0.5,
                 binwidth = binwidth,
                 stackdir = "center",
                 stackratio = 0.5,
                 method = "histodot") +
    ggthemes::theme_base() +
    labs(x = "", y = "") +
    scale_x_discrete(expand = expansion(mult = c(0, 0))) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
}

plot_bar <- function(data, colors, levels_order) {
  ggplot(data.frame(class = "class", value = data), aes(x = class)) +
    geom_bar(
      aes(fill = factor(value, levels = levels_order)),
      color = "grey90",
      position = "stack",
      show.legend = FALSE,
      width = 0.8
    ) +
    scale_fill_manual(values = colors) +
    ggthemes::theme_base() +
    labs(x = "", y = "") +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
}


plot_age <- plot_box_dot(df$Age, age_color, binwidth = 0.3)
plot_age
plot_bmi <- plot_box_dot(df$BMI, bmi_color, binwidth = 0.3)
plot_bmi

plot_birth_weight <- plot_box_dot(df$birth_weight, birth_weight_color, binwidth = 50)
plot_birth_weight

plot_birth_length <- plot_box_dot(df$birth_length, birth_length_color, 0.3)
plot_birth_length

plot_psmoke <- plot_bar(df$Psmoke, psmoke_color, c("0", "1"))
plot_psmoke

plot_edu <- plot_bar(df$Edu, edu_color, c("1", "2", "3", "4"))
plot_edu
plot_income <- plot_bar(df$Income, income_color, c("1", "2", "3"))
plot_income

time_data <- readxl::read_xlsx('~/date.xlsx')
proteome_dat <- read_csv('cleaned_data.csv')
time_data_filtered <- time_data |>
  dplyr::filter(number %in% proteome_dat$number_identifier) |>
  dplyr::filter(lubridate::year(day_start) >= 2015)

x_axis_breaks <- seq(
  from = lubridate::as_datetime("2015-01-01"),
  to = lubridate::as_datetime("2017-07-01"),
  by = "6 months"
)

p_t <- time_data_filtered |>
  dplyr::filter(lubridate::year(day_start) >= 2015) |>
  dplyr::mutate(number = forcats::fct_reorder(factor(number), number, .desc = TRUE)) |>
  ggplot2::ggplot(ggplot2::aes(y = number)) +
  ggplot2::geom_segment(
    ggplot2::aes(x = day_start, xend = day_birth),
    color = "gray90",
    linewidth = 0.2,
    alpha = 0.9
  ) +
  ggplot2::geom_point(ggplot2::aes(x = day_start), color = "#91CAE8", size = 1) +
  ggplot2::geom_point(ggplot2::aes(x = day_birth), color = "#F48892", size = 1) +
  ggplot2::ylab(NULL) +
  xlab(NULL)+
  ggplot2::scale_x_datetime(
    breaks = x_axis_breaks,
    date_labels = "%Y-%m-%d",
    limits = c(min(x_axis_breaks), max(x_axis_breaks))
  ) +
  ggthemes::theme_base() +
  ggplot2::theme(
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank()
  )
p_t


time_data_filtered2 <- time_data_filtered %>%
  filter(lubridate::year(day_start) >= 2015) %>%
  mutate(number = fct_reorder(factor(number), number, .desc = TRUE))
p_h <- ggplot(time_data_filtered2, aes(x = 1, y = as.numeric(number), fill = total_days)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#F48892") +
  labs(x = "",
       y = "Number",
       fill = "Total Days") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank()) +
  coord_fixed(ratio = 0.1)
p_h


rm(list = ls())
gc()
library(readr)
proteome_dat <- read_csv('cleaned_data.csv')
phenotype_dat <- read_csv('phenotype.csv')
uniprot_chromosome <- read.delim('idmapping_uniprot_chromosome.tsv')
uniprot_chromosome2 <- uniprot_chromosome |>
  dplyr::select(From,Length,Proteomes)

uniprot_chromosome2$Chromosome <- gsub(".*: (.*)", "\\1", uniprot_chromosome2$Proteomes)

uniprot_chromosome2 <- uniprot_chromosome2 |> na.omit()
uniprot_ids <- uniprot_chromosome2$From

phenotype_dat_birthweight <- phenotype_dat |> dplyr::select(number,birth_weight)

birthweight_cor_dat <- merge(phenotype_dat_birthweight,
                             proteome_dat |> dplyr::select(number_identifier,all_of(uniprot_ids)),
                             by.x = 'number',
                             by.y = 'number_identifier') |>
  na.omit()

protein_data <- birthweight_cor_dat |>
  dplyr::select(number,uniprot_ids)

library(ggplot2)
library(dplyr)
library(tidyr)
library(gtools)

protein_data_long <- protein_data %>%
  pivot_longer(cols = -number,
               names_to = "protein_id",
               values_to = "expression")

protein_chrom_data <- protein_data_long %>%
  left_join(uniprot_chromosome2 %>% select(From, Chromosome),
            by = c("protein_id" = "From"))

protein_chrom_data$Chr_label <- gsub("Chromosome (\\d+|X|Y|MT|Unplaced)", "Chr \\1", protein_chrom_data$Chromosome)

chrom_order <- mixedsort(unique(protein_chrom_data$Chr_label))
protein_chrom_data$Chr_label <- factor(protein_chrom_data$Chr_label,levels = chrom_order)

protein_chrom_data$log_expression <- log10(protein_chrom_data$expression + 1)
df_chr <- as.data.frame(table(uniprot_chromosome2$Chromosome))

p_chr <- ggplot(protein_chrom_data, aes(x = Chr_label, y = log_expression)) +
  geom_boxplot(fill = "#FFE7E7", alpha = 0.7, outlier.shape = NA) +
  ggthemes::theme_clean() +
  labs(title = "Protein Expression by Chromosome (log10 scale)",
       x = "Chromosome",
       y = "Expression Level (log10)") +
  ggthemes::theme_base() +

  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

q1 <- quantile(protein_chrom_data$log_expression, 0.25)
q3 <- quantile(protein_chrom_data$log_expression, 0.75)
iqr <- q3 - q1
lower_bound <- q1 - 1.5 * iqr
upper_bound <- q3 + 1.5 * iqr

p_chr <- p_chr + coord_cartesian(ylim = c(lower_bound, upper_bound))
p_chr

protein_count_by_chr <- protein_chrom_data %>%
  distinct(protein_id, Chr_label) %>%
  count(Chr_label, name = "protein_count")

p_chr_combined <- p_chr +
  geom_line(data = protein_count_by_chr,
            aes(x = Chr_label, y = scales::rescale(protein_count, to = c(lower_bound, upper_bound)),
                group = 1),
            color = "#ebb9b9", size = 1) +
  geom_point(data = protein_count_by_chr,
             aes(x = Chr_label, y = scales::rescale(protein_count, to = c(lower_bound, upper_bound))),
             color = "#ebb9b9", size = 2) +
  scale_y_continuous(
    name = "Expression (log10)",
    sec.axis = sec_axis(~ scales::rescale(., from = c(lower_bound, upper_bound), to = c(min(protein_count_by_chr$protein_count), max(protein_count_by_chr$protein_count))),
                        name = "number")
  )

p_chr_combined

library(tidyverse)
hpa_placenta_detected_df <- read.delim('HPA-placenta-detected.tsv')
hpa_placenta_detected_df2 <- hpa_placenta_detected_df |>
  dplyr::filter(Uniprot != '')


hpa_placenta_elevated_df <- read.delim('HPA-placenta-elevated.tsv')
hpa_placenta_elevated_df2 <- hpa_placenta_elevated_df |>
  dplyr::filter(Uniprot != '')

raw_protein_data <- openxlsx::read.xlsx('~/Placenta_4Dproteome.xlsx',sheet = 2)
raw_gene <- raw_protein_data$Gene.name

df1 <- data.frame(
  current_detected = 7108,
  HPA_detected = nrow(hpa_placenta_detected_df2)
)

proteome_dat <- read_csv('cleaned_data.csv')
self_protein <- names(proteome_dat)[2:ncol(proteome_dat)]
length(intersect(self_protein,hpa_placenta_elevated_df2$Uniprot))

df2 <- data.frame(
  current_placenta_elevated = length(intersect(self_protein,hpa_placenta_elevated_df2$Uniprot)),
  HPA_placenta_elevated = nrow(hpa_placenta_elevated_df2)
)

detected_percentage <- round(df1$current_detected / df1$HPA_detected * 100, 2)
elevated_percentage <- round(df2$current_placenta_elevated / df2$HPA_placenta_elevated * 100, 2)

detected_plot <- ggplot() +
  geom_bar(data = data.frame(x = "HPA Database", y = df1$HPA_detected),
           aes(x = x, y = y), stat = "identity", fill = "#7FB3D5", width = 0.6) +
  geom_bar(data = data.frame(x = "Our Study", y = df1$current_detected),
           aes(x = x, y = y), stat = "identity", fill = "#EC7063", width = 0.6) +
  geom_text(data = data.frame(x = "Our Study", y = df1$current_detected + 200,
                              label = paste0(df1$current_detected, " (", detected_percentage, "%)")),
            aes(x = x, y = y, label = label), size = 5) +
  geom_text(data = data.frame(x = "HPA Database", y = df1$HPA_detected + 200,
                              label = as.character(df1$HPA_detected)),
            aes(x = x, y = y, label = label), size = 5) +
  labs(title = "Comparison of Detected Proteins in Placenta",
       subtitle = paste0("Our study detected ", detected_percentage, "% of proteins in HPA database"),
       x = "", y = "Number of Proteins") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )

elevated_plot <- ggplot() +
  geom_bar(data = data.frame(x = "HPA Database", y = df2$HPA_placenta_elevated),
           aes(x = x, y = y), stat = "identity", fill = "#7FB3D5", width = 0.6) +
  geom_bar(data = data.frame(x = "Our Study", y = df2$current_placenta_elevated),
           aes(x = x, y = y), stat = "identity", fill = "#EC7063", width = 0.6) +
  geom_text(data = data.frame(x = "Our Study", y = df2$current_placenta_elevated + 20,
                              label = paste0(df2$current_placenta_elevated, " (", elevated_percentage, "%)")),
            aes(x = x, y = y, label = label), size = 5) +
  geom_text(data = data.frame(x = "HPA Database", y = df2$HPA_placenta_elevated + 20,
                              label = as.character(df2$HPA_placenta_elevated)),
            aes(x = x, y = y, label = label), size = 5) +
  labs(title = "Comparison of Placenta-Elevated Proteins",
       subtitle = paste0("Our study detected ", elevated_percentage, "% of placenta-elevated proteins in HPA database"),
       x = "", y = "Number of Proteins") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )

library(patchwork)
combined_plot <- detected_plot / elevated_plot
combined_plot

stacked_data <- data.frame(
  Source = rep(c("HPA Database", "Our Study"), each = 2),
  Type = rep(c("Elevated", "Other"), 2),
  Count = c(
    df2$HPA_placenta_elevated,
    df1$HPA_detected - df2$HPA_placenta_elevated,
    df2$current_placenta_elevated,
    df1$current_detected - df2$current_placenta_elevated
  )
)
total_heights <- aggregate(Count ~ Source, stacked_data, sum)

integrated_plot <- ggplot(stacked_data, aes(x = Source, y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  scale_fill_manual(values = c("Elevated" = "#EC7063", "Other" = "#7FB3D5")) +
  geom_text(data = data.frame(
    Source = c("HPA Database", "HPA Database", "Our Study", "Our Study"),
    Type = c("Other", "Elevated", "Other", "Elevated"),
    y = c(
      df1$HPA_detected + 500,
      df2$HPA_placenta_elevated / 2,
      df1$current_detected + 500,
      df2$current_placenta_elevated / 2
    ),
    label = c(
      as.character(df1$HPA_detected),
      paste0(df2$HPA_placenta_elevated, " (elevated)"),
      paste0(df1$current_detected, " (", detected_percentage, "%)"),
      paste0(df2$current_placenta_elevated, " (elevated)")
    )
  ), aes(x = Source, y = y, label = label), size = 4, vjust = 0.5) +
  labs(       x = "", y = "Number of Proteins") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "bottom",
    legend.title = element_blank()
  )

print(integrated_plot)


library(tidyverse)
hpa_singlecell_blood_immune <- read.delim('blood & immune cell-specific proteome.tsv')
hpa_singlecell_blood_immune2 <- hpa_singlecell_blood_immune |>
  dplyr::select(Gene,Uniprot,RNA.blood.cell.distribution, RNA.single.cell.type.specific.nTPM)

hpa_singlecell_blood_immune2 <- hpa_singlecell_blood_immune2 |>
  mutate(Hofbauer_Tcells_info = sapply(RNA.single.cell.type.specific.nTPM, function(x) {
    has_hofbauer <- grepl("Hofbauer cells", x)
    has_tcells <- grepl("T-cells", x)

    result <- c()
    if (has_hofbauer) result <- c(result, "Hofbauer cells")
    if (has_tcells) result <- c(result, "T-cells")

    if (length(result) > 0) {
      paste(result, collapse=";")
    } else {
      NA_character_
    }
  })) |>
  filter(!is.na(Hofbauer_Tcells_info)) |>
  dplyr::select(Gene,Uniprot,Hofbauer_Tcells_info)

head(hpa_singlecell_blood_immune2)

hpa_singlecell_blood_immune_long <- hpa_singlecell_blood_immune2 |>
  tidyr::separate_rows(Hofbauer_Tcells_info, sep=";")

head(hpa_singlecell_blood_immune_long)


proteome_dat <- read_csv('cleaned_data.csv')
self_protein <- names(proteome_dat)[2:ncol(proteome_dat)]

hpa_singlecell_blood_immune_self <- hpa_singlecell_blood_immune_long |>
  dplyr::filter(Uniprot %in% self_protein)



hpa_singlecell_trophoblast <- read.delim('trophoblast cell-specific proteome.tsv')
hpa_singlecell_trophoblast2 <- hpa_singlecell_trophoblast |>
  dplyr::select(Gene,Uniprot,RNA.trophoblast.cell.distribution, RNA.single.cell.type.specific.nTPM)

hpa_singlecell_trophoblast2 <- hpa_singlecell_trophoblast2 |>
  mutate(Trophoblast_info = sapply(RNA.single.cell.type.specific.nTPM, function(x) {
    has_cytotrophoblasts <- grepl("Cytotrophoblasts", x)
    has_syncytiotrophoblasts <- grepl("Syncytiotrophoblasts", x)
    has_extravillous_trophoblasts <- grepl("Extravillous trophoblasts", x)

    result <- c()
    if (has_cytotrophoblasts) result <- c(result, "Cytotrophoblasts")
    if (has_syncytiotrophoblasts) result <- c(result, "Syncytiotrophoblasts")
    if (has_extravillous_trophoblasts) result <- c(result, "Extravillous trophoblasts")

    if (length(result) > 0) {
      paste(result, collapse=";")
    } else {
      NA_character_
    }
  })) |>
  filter(!is.na(Trophoblast_info)) |>
  dplyr::select(Gene,Uniprot,Trophoblast_info)

hpa_singlecell_trophoblast_long <- hpa_singlecell_trophoblast2 |>
  tidyr::separate_rows(Trophoblast_info, sep=";")

hpa_singlecell_trophoblast_self <- hpa_singlecell_trophoblast_long |>
  dplyr::filter(Uniprot %in% self_protein)


hpa_singlecell_endothelial <- read.delim('endothelial cell-specific proteome.tsv')
hpa_singlecell_endothelial2 <- hpa_singlecell_endothelial |>
  dplyr::select(Gene,Uniprot,RNA.single.cell.type.specific.nTPM)

hpa_singlecell_endothelial2 <- hpa_singlecell_endothelial2 |>
  mutate(Endothelial_info = sapply(RNA.single.cell.type.specific.nTPM, function(x) {
    has_endothelial <- grepl("Endothelial cells", x)

    result <- c()
    if (has_endothelial) result <- c(result, "Endothelial cells")
    if (length(result) > 0) {
      paste(result, collapse=";")
    } else {
      NA_character_
    }
  })) |>
  filter(!is.na(Endothelial_info)) |>
  dplyr::select(Gene,Uniprot,Endothelial_info)

hpa_singlecell_endothelial_long <- hpa_singlecell_endothelial2 |>
  tidyr::separate_rows(Endothelial_info, sep=";")

hpa_singlecell_endothelial_self <- hpa_singlecell_endothelial_long |>
  dplyr::filter(Uniprot %in% self_protein)
table(hpa_singlecell_endothelial_self$Endothelial_info)

hpa_singlecell_muscle <- read.delim('Muscle cell-specific proteome.tsv')
hpa_singlecell_muscle2 <- hpa_singlecell_muscle |>
  dplyr::select(Gene,Uniprot,RNA.single.cell.type.specific.nTPM)

hpa_singlecell_muscle2 <- hpa_singlecell_muscle2 |>
  mutate(smooth_muscle_info = sapply(RNA.single.cell.type.specific.nTPM, function(x) {
    has_smooth_muscle <- grepl("Smooth muscle cells", x)

    result <- c()
    if (has_smooth_muscle) result <- c(result, "Smooth muscle cells")
    if (length(result) > 0) {
      paste(result, collapse=";")
    } else {
      NA_character_
    }
  })) |>
  filter(!is.na(smooth_muscle_info)) |>
  dplyr::select(Gene,Uniprot,smooth_muscle_info)

hpa_singlecell_muscle_long <- hpa_singlecell_muscle2 |>
  tidyr::separate_rows(smooth_muscle_info, sep=";")

hpa_singlecell_muscle_self <- hpa_singlecell_muscle_long |>
  dplyr::filter(Uniprot %in% self_protein)
table(hpa_singlecell_muscle_self$smooth_muscle_info)



hpa_singlecell_mesenchymal <- read.delim('mesenchymal cell-specific proteome.tsv')
hpa_singlecell_mesenchymal2 <- hpa_singlecell_mesenchymal |>
  dplyr::select(Gene,Uniprot,RNA.single.cell.type.specific.nTPM)

hpa_singlecell_mesenchymal2 <- hpa_singlecell_mesenchymal2 |>
  mutate(mesenchymal_info = sapply(RNA.single.cell.type.specific.nTPM, function(x) {
    has_fibroblasts <- grepl("Fibroblasts", x)

    result <- c()
    if (has_fibroblasts) result <- c(result, "Fibroblasts")
    if (length(result) > 0) {
      paste(result, collapse=";")
    } else {
      NA_character_
    }
  })) |>
  filter(!is.na(mesenchymal_info)) |>
  dplyr::select(Gene,Uniprot,mesenchymal_info)

hpa_singlecell_mesenchymal_long <- hpa_singlecell_mesenchymal2 |>
  tidyr::separate_rows(mesenchymal_info, sep=";")

hpa_singlecell_mesenchymal_self <- hpa_singlecell_mesenchymal_long |>
  dplyr::filter(Uniprot %in% self_protein)
table(hpa_singlecell_mesenchymal_self$mesenchymal_info)


hpa_singlecell_mesenchymal <- read.delim('mesenchymal cell-specific proteome.tsv')
hpa_singlecell_mesenchymal2 <- hpa_singlecell_mesenchymal |>
  dplyr::select(Gene,Uniprot,RNA.single.cell.type.specific.nTPM)

hpa_singlecell_mesenchymal2 <- hpa_singlecell_mesenchymal2 |>
  mutate(mesenchymal_info = sapply(RNA.single.cell.type.specific.nTPM, function(x) {
    has_fibroblasts <- grepl("Fibroblasts", x)

    result <- c()
    if (has_fibroblasts) result <- c(result, "Fibroblasts")
    if (length(result) > 0) {
      paste(result, collapse=";")
    } else {
      NA_character_
    }
  })) |>
  filter(!is.na(mesenchymal_info)) |>
  dplyr::select(Gene,Uniprot,mesenchymal_info)

hpa_singlecell_mesenchymal_long <- hpa_singlecell_mesenchymal2 |>
  tidyr::separate_rows(mesenchymal_info, sep=";")

hpa_singlecell_mesenchymal_self <- hpa_singlecell_mesenchymal_long |>
  dplyr::filter(Uniprot %in% self_protein)
table(hpa_singlecell_mesenchymal_self$mesenchymal_info)

single_cell_data1 <- read.xlsx('endothelial.xlsx')
names(single_cell_data1)[3] <- 'single_cell_location'
single_cell_data2 <- read.xlsx('Fibroblasts.xlsx')
names(single_cell_data2)[3] <- 'single_cell_location'
single_cell_data3 <- read.xlsx('Smooth muscle cells.xlsx')
names(single_cell_data3)[3] <- 'single_cell_location'
single_cell_data4 <- read.xlsx('trophoblast.xlsx')
names(single_cell_data4)[3] <- 'single_cell_location'
single_cell_data5 <- read.xlsx('immune_cells.xlsx')
names(single_cell_data5)[3] <- 'single_cell_location'
single_cell_data <- bind_rows(single_cell_data1,single_cell_data2,single_cell_data3,single_cell_data4,single_cell_data5)

wgcna_results <- read.xlsx('WGANA-module.xlsx')
table(wgcna_results$moduleColor)
wgcna_results_royalblue <- wgcna_results |> dplyr::filter(moduleColor == 'royalblue') |>
  dplyr::left_join(single_cell_data,by = c('geneSymbol' = 'Uniprot'),keep = T)

wgcna_results_green <- wgcna_results |> dplyr::filter(moduleColor == 'green') |>
  dplyr::left_join(single_cell_data,by = c('geneSymbol' = 'Uniprot'),keep = T)

wgcna_results_midnightblue <- wgcna_results |> dplyr::filter(moduleColor == 'midnightblue') |>
  dplyr::left_join(single_cell_data,by = c('geneSymbol' = 'Uniprot'),keep = T)

visual_order <- c(
  "T-cells", "Hofbauer cells", "Cytotrophoblasts", "Syncytiotrophoblasts",
  "Extravillous trophoblasts", "Endothelial cells", "Fibroblasts",
  "Smooth muscle cells", "Not Annotated"
)

custom_colors <- c(
  "Cytotrophoblasts" = "#ddc2c0",
  "Endothelial cells" = "#139ba2",
  "Extravillous trophoblasts" = "#F2A1A7",
  "Fibroblasts" = "#93a7b2",
  "Hofbauer cells" = "#c2b6e7",
  "Smooth muscle cells" = "#7f5e46",
  "Syncytiotrophoblasts" = "#FBDDDD",
  "T-cells" = "#83477d",
  "Not Annotated" = "#f2f2f2"
)

factor_level_order <- rev(visual_order)

plot_data <- wgcna_results |>
  dplyr::filter(moduleColor %in% c("royalblue", "green", "midnightblue")) |>
  dplyr::left_join(single_cell_data, by = c("geneSymbol" = "Uniprot")) |>
  dplyr::mutate(
    single_cell_location = forcats::fct_explicit_na(single_cell_location, na_level = "Not Annotated")
  ) |>
  dplyr::count(moduleColor, single_cell_location, name = "gene_count") |>
  dplyr::mutate(
    single_cell_location = forcats::fct_relevel(single_cell_location, factor_level_order)
  )
p_single <- ggplot2::ggplot(plot_data, ggplot2::aes(x = moduleColor, y = gene_count, fill = single_cell_location)) +
  ggplot2::geom_col(position = "fill") +
  ggplot2::scale_y_continuous(labels = scales::percent) +
  ggplot2::scale_fill_manual(values = custom_colors) +
  ggplot2::coord_flip() +
  ggplot2::labs(
    x = "Module Color",
    y = "Proportion of Genes",
    fill = "Single-Cell Location",
    title = "Proportional Distribution of Genes Across Single-Cell Locations"
  ) +
  ggplot2::theme_bw()


subcellular_df <- read.xlsx('cellcular.xlsx')
wgcna_results <- read.xlsx('WGANA-module.xlsx')
wgcna_results_royalblue_subcellular <- wgcna_results |> dplyr::filter(moduleColor == 'royalblue') |>
  dplyr::left_join(subcellular_df,by = c('geneSymbol' = 'From'),keep = T) |>
  dplyr::select(geneSymbol,subcellular_classfied)

wgcna_results_green_subcellular <- wgcna_results |> dplyr::filter(moduleColor == 'green') |>
  dplyr::left_join(subcellular_df,by = c('geneSymbol' = 'From'),keep = T)|>
  dplyr::select(geneSymbol,subcellular_classfied)

wgcna_results_midnightblue_subcellular <- wgcna_results |> dplyr::filter(moduleColor == 'midnightblue') |>
  dplyr::left_join(subcellular_df,by = c('geneSymbol' = 'From'),keep = T)|>
  dplyr::select(geneSymbol,subcellular_classfied)



visual_order_subcellular <- c(
  "Nucleus", "Endoplasmic reticulum", "Golgi apparatus", "Mitochondrion",
  "Peroxisome", "Lysosome", "Endosome", "Melanosome", "Cytoplasm",
  "Cytoskeleton", "Centrosome", "Cell membrane", "Secreted",
  "Extracellular matrix", "Extracellular space", "Not Annotated"
)

factor_level_order_subcellular <- rev(visual_order_subcellular)

custom_colors_subcellular <- c(
  "Nucleus" = "#9bbbe1",
  "Endoplasmic reticulum" = "#e8e8e8",
  "Golgi apparatus" = "#c3c3c4",
  "Mitochondrion" = "#a5a5a7",
  "Peroxisome" = "#878789",
  "Lysosome" = "#69696c",
  "Endosome" = "#4b4b4e",
  "Melanosome" = "#2d2d31",
  "Cytoplasm" = "#b7b7eb",
  "Cytoskeleton" = "#c9c9f0",
  "Centrosome" = "#ededfa",
  "Cell membrane" = "#eab883",
  "Secreted" = "#f09ba0",
  "Extracellular matrix" = "#f4b4b8",
  "Extracellular space" = "#fbe6e7",
  "Not Annotated" = "#f2f2f2"
)

plot_data_subcellular <- wgcna_results |>
  dplyr::filter(moduleColor %in% c("royalblue", "green", "midnightblue")) |>
  dplyr::left_join(subcellular_df, by = c("geneSymbol" = "From")) |>
  dplyr::mutate(
    subcellular_classfied = tidyr::replace_na(subcellular_classfied, "Not Annotated")
  ) |>
  tidyr::separate_rows(subcellular_classfied, sep = ";") |>
  dplyr::mutate(
    subcellular_classfied = stringr::str_trim(subcellular_classfied)
  ) |>
  dplyr::count(moduleColor, subcellular_classfied, name = "location_count") |>
  dplyr::filter(subcellular_classfied != "Other") |>
  dplyr::mutate(
    subcellular_classfied = forcats::fct_relevel(subcellular_classfied, factor_level_order_subcellular)
  )

p_subcellcular <- ggplot2::ggplot(plot_data_subcellular, ggplot2::aes(x = moduleColor, y = location_count, fill = subcellular_classfied)) +
  ggplot2::geom_col(position = "fill") +
  ggplot2::scale_y_continuous(labels = scales::percent) +
  ggplot2::scale_fill_manual(values = custom_colors_subcellular) +
  ggplot2::coord_flip() +
  ggplot2::labs(
    x = "Module Color",
    y = "Proportion of Occurrences",
    fill = "Subcellular Location",
    title = "Proportional Subcellular Location Distribution by Module"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "none")


royalblue_bp_df <- read.delim('royalblue_BP_level2.txt')
royalblue_bp_df_sig <- royalblue_bp_df |> dplyr::filter(FDR<= 0.05) |> mutate(module = 'M1') |> dplyr::select(module,Term,FDR)
midnightblue_bp_df <- read.delim('midnightblue_BP_level2.txt')
midnightblue_bp_df_sig <- midnightblue_bp_df |> dplyr::filter(FDR<= 0.05) |> mutate(module = 'M2') |> dplyr::select(module,Term,FDR)
green_bp_df <- read.delim('green_BP_level2.txt')
green_bp_df_sig <- green_bp_df |> dplyr::filter(FDR<= 0.05) |> mutate(module = 'M3') |> dplyr::select(module,Term,FDR)

GO_BP_df <- bind_rows(royalblue_bp_df_sig, midnightblue_bp_df_sig, green_bp_df_sig) |>
  dplyr::mutate(
    Term = stringr::str_remove(Term, "GO:\\d+~")
  )

royalblue_mf_df <- read.delim('royalblue_MF_level2.txt')
royalblue_mf_df_sig <- royalblue_mf_df |> dplyr::filter(FDR<= 0.05) |> mutate(module = 'M1') |> dplyr::select(module,Term,FDR)
midnightblue_mf_df <- read.delim('midnightblue_MF_level2.txt')
midnightblue_mf_df_sig <- midnightblue_mf_df |> dplyr::filter(FDR<= 0.05)|> mutate(module = 'M2') |> dplyr::select(module,Term,FDR)
green_mf_df <- read.delim('green_MF_level2.txt')
green_mf_df_sig <- green_mf_df |> dplyr::filter(FDR<= 0.05)|> mutate(module = 'M3') |> dplyr::select(module,Term,FDR)

GO_MF_df <- bind_rows(royalblue_mf_df_sig, midnightblue_mf_df_sig, green_mf_df_sig) |>
  dplyr::mutate(
    Term = stringr::str_remove(Term, "GO:\\d+~")
  )

GO_df <- bind_rows(GO_BP_df,GO_MF_df)

go_term_bp_mf <- bind_rows(go_terms_bp,go_terms_mf)

GO_df_withCategory <- merge(GO_df,go_term_bp_mf,by.x = 'Term',by.y = 'GOTerm')

plot_data <- GO_df_withCategory |>
  dplyr::mutate(neg_log10_fdr = -log10(FDR)) |>
  dplyr::mutate(Category = factor(Category, levels = class_order)) |>
  dplyr::arrange(Category, FDR) |>
  dplyr::mutate(Term = forcats::fct_inorder(Term))

p_go <- plot_data |>
  ggplot2::ggplot(ggplot2::aes(x = module, y = Term)) +
  ggplot2::geom_tile(ggplot2::aes(fill = neg_log10_fdr), color = "white", size = 0.1) +
  ggplot2::facet_grid(
    Category ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  ggplot2::scale_fill_gradient(low = "#ffd8c6", high = "#ff4f00") +
  ggplot2::labs(
    x = NULL,
    y = "",
    fill = "-log10(q-value)"
  ) +
  ggthemes::theme_base() +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(size = 10),
    strip.placement = "outside",
    strip.background = ggplot2::element_blank(),
    strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 0, size = 11)
  )
p_go
topptx(p_go,'GO_analysis_of_3_modules_2.pptx',width = 6.5,height = 20)



phenotype_dat_birthweight <- phenotype_dat |> dplyr::select(number,birth_weight)

birthweight_cor_dat <- merge(phenotype_dat_birthweight,
                             proteome_dat,
                             by.x = 'number',
                             by.y = 'number_identifier') |>
  na.omit()

protein_columns <- colnames(birthweight_cor_dat)[3:ncol(birthweight_cor_dat)]

correlation_results <- purrr::map_dfr(protein_columns, ~{
  test_result <- stats::cor.test(birthweight_cor_dat[['birth_weight']], birthweight_cor_dat[[.x]], method = "pearson")
  tibble::tibble(
    protein = .x,
    correlation = test_result$estimate,
    p_value = test_result$p.value
  )
})

length(which(correlation_results$p_value<=0.05))
uniprot_ids <- correlation_results |>
  dplyr::distinct(protein) |>
  dplyr::pull(protein)

gene_mapping <- AnnotationDbi::select(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = uniprot_ids,
  columns = "SYMBOL",
  keytype = "UNIPROT"
) |>
  dplyr::as_tibble() |>
  dplyr::rename(protein = UNIPROT, gene_name = SYMBOL) |>
  dplyr::distinct(protein, .keep_all = TRUE)

correlation_results_with_gene_names <- correlation_results |>
  dplyr::left_join(gene_mapping, by = "protein")

correlation_results_with_gene_names_sig <- correlation_results_with_gene_names |> filter(p_value <= 0.05)
correlation_results_with_gene_names |> filter(p_value <= 0.05) |>
  write_csv('test_birthweight.csv')

reactome_enrichment_df <- read_xlsx('~/birth_weight_enrichment_results.xlsx')
uniprot_pathway_mapping <- reactome_enrichment_df |>
  dplyr::select(pathway, members_input_overlap) |>
  tidyr::separate_longer_delim(cols = members_input_overlap, delim = ";") |>
  dplyr::mutate(members_input_overlap = stringr::str_trim(members_input_overlap)) |>
  dplyr::filter(members_input_overlap != "") |>
  dplyr::rename(uniprot_id = members_input_overlap)

correlation_results_with_gene_names_sig_withpathway = merge(correlation_results_with_gene_names_sig,
                                                            uniprot_pathway_mapping,
                                                            by.x = 'protein',
                                                            by.y = 'uniprot_id',all.x = T)


phenotype_dat <- read_csv('phenotype.csv')
phenotype_dat_birthweight <- phenotype_dat |> dplyr::select(number,birth_weight,birth_length,Weight.Z, Length.Z,GA_age) |>
  mutate(Ponderal_Index = (birth_weight/(birth_length ^3))*100,
         BMI = birth_weight/(birth_length/100)^2/1000)

metabolite_data <- openxlsx::read.xlsx('~/placenta_metabolite.xlsx')
metabolite_data_info <- openxlsx::read.xlsx('~/placenta_metabolite.xlsx',sheet = 'For_Figure')
metabolite_selected <- openxlsx::read.xlsx('~/placenta_metabolite.xlsx',sheet = 3)
metabolite_data_selected <- metabolite_data |> select(number,metabolite_selected$metabo_id)


phenotype_metabolite <- merge(phenotype_dat_birthweight,
                             metabolite_data_selected,
                             by.x = 'number',
                             by.y = 'number')
library(tidyverse)
library(broom)


phenotype_vars <- c("birth_weight", "birth_length", "Weight.Z", "Length.Z", "GA_age", "Ponderal_Index", "BMI")
metabolite_vars <- colnames(phenotype_metabolite) |>
  stringr::str_subset("^Metabo_")

variable_pairs <- tidyr::expand_grid(
  phenotype = phenotype_vars,
  metabolite = metabolite_vars
)

correlation_table <- variable_pairs |>
  purrr::pmap_dfr(\(phenotype, metabolite) {
    test_result <- stats::cor.test(
      phenotype_metabolite[[phenotype]],
      phenotype_metabolite[[metabolite]],
      method = "pearson"
    )

    broom::tidy(test_result) |>
      dplyr::mutate(phenotype = phenotype, metabolite = metabolite)
  }) |>
  dplyr::select(
    phenotype,
    metabolite,
    correlation_coefficient = estimate,
    p_value = p.value
  )

correlation_table_weight <- correlation_table |>
  dplyr::filter(phenotype == 'Weight.Z') |>
  dplyr::left_join(
    metabolite_selected |> dplyr::select(metabo_id, id_kegg, Standardized.name),
    by = c("metabolite" = "metabo_id")
  )
plot_data <- read.xlsx('weight_metabolite_relationship.xlsx',sheet = 'deduplicated')
plot_data_birthweight <- plot_data |> dplyr::filter(p_value< 0.05) |> pull(Standardized.name)
plot_data_birthweight_metabo_id <- plot_data |> dplyr::filter(p_value< 0.05) |> pull(metabolite)
p <- plot_data |>
  dplyr::mutate(
    significance_group = dplyr::case_when(
      p_value < 0.05 & correlation_coefficient < 0 ~ "Negative Significant",
      p_value < 0.05 & correlation_coefficient > 0 ~ "Positive Significant",
      TRUE ~ "Not Significant"
    )
  ) |>
  ggplot2::ggplot(ggplot2::aes(x = correlation_coefficient, y = -log10(p_value), color = significance_group)) +
  ggplot2::geom_point(alpha = 0.7, size = 2) +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey30") +
  ggplot2::scale_color_manual(
    values = c(
      "Positive Significant" = "#d62728",
      "Negative Significant" = "#31599f",
      "Not Significant" = "grey"
    ),
    breaks = c("Positive Significant", "Negative Significant", "Not Significant")
  ) +
  ggplot2::labs(
    x = "Correlation Coefficient",
    y = bquote(-log[10]~'(P-value)'),
    title = "Volcano Plot",
    color = "Significance"
  ) +
  theme_bw()
plot_data |>
  dplyr::mutate(
    significance_group = dplyr::case_when(
      p_value < 0.05 & correlation_coefficient < 0 ~ "Negative Significant",
      p_value < 0.05 & correlation_coefficient > 0 ~ "Positive Significant",
      TRUE ~ "Not Significant"
    )
  ) |> dplyr::count(significance_group)

topptx(p,'weight_metabolite_volcano.pptx',height = 4,width = 5)

correlation_table_BMI <- correlation_table |> dplyr::filter(phenotype == 'BMI') |>
  dplyr::left_join(
    metabolite_selected |> dplyr::select(metabo_id, id_kegg, Standardized.name),
    by = c("metabolite" = "metabo_id")
  )
plot_data <- read.xlsx('BMI_metabolite_relationship.xlsx',sheet = 'deduplicated')
plot_data_BMI <- plot_data |> dplyr::filter(p_value< 0.05) |> pull(Standardized.name)
plot_data_BMI_metabo_id <- plot_data |> dplyr::filter(p_value< 0.05) |> pull(metabolite)
p <- plot_data |>
  dplyr::mutate(
    significance_group = dplyr::case_when(
      p_value < 0.05 & correlation_coefficient < 0 ~ "Negative Significant",
      p_value < 0.05 & correlation_coefficient > 0 ~ "Positive Significant",
      TRUE ~ "Not Significant"
    )
  ) |>
  ggplot2::ggplot(ggplot2::aes(x = correlation_coefficient, y = -log10(p_value), color = significance_group)) +
  ggplot2::geom_point(alpha = 0.7, size = 2) +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey30") +
  ggplot2::scale_color_manual(
    values = c(
      "Positive Significant" = "#d62728",
      "Negative Significant" = "#31599f",
      "Not Significant" = "grey"
    ),
    breaks = c("Positive Significant", "Negative Significant", "Not Significant")
  ) +
  ggplot2::labs(
    x = "Correlation Coefficient",
    y = bquote(-log[10]~'(P-value)'),
    title = "Volcano Plot",
    color = "Significance"
  ) +
  theme_bw()

plot_data |>
  dplyr::mutate(
    significance_group = dplyr::case_when(
      p_value < 0.05 & correlation_coefficient < 0 ~ "Negative Significant",
      p_value < 0.05 & correlation_coefficient > 0 ~ "Positive Significant",
      TRUE ~ "Not Significant"
    )
  ) |> dplyr::count(significance_group)

topptx(p,'BMI_metabolite_volcano.pptx',height = 4,width = 5)

correlation_table_pi <- correlation_table |> dplyr::filter(phenotype == 'Ponderal_Index') |>
  dplyr::left_join(
    metabolite_selected |> dplyr::select(metabo_id, id_kegg, Standardized.name),
    by = c("metabolite" = "metabo_id")
  )
openxlsx::write.xlsx(correlation_table_pi,'PI_metabolite_relationship.xlsx')

plot_data <- read.xlsx('PI_metabolite_relationship.xlsx',sheet = 'deduplicated')
plot_data_PI <- plot_data |> dplyr::filter(p_value< 0.05) |> pull(Standardized.name)
plot_data_PI_metabo_id <- plot_data |> dplyr::filter(p_value< 0.05) |> pull(metabolite)
p <- plot_data |>
  dplyr::mutate(
    significance_group = dplyr::case_when(
      p_value < 0.05 & correlation_coefficient < 0 ~ "Negative Significant",
      p_value < 0.05 & correlation_coefficient > 0 ~ "Positive Significant",
      TRUE ~ "Not Significant"
    )
  ) |>
  ggplot2::ggplot(ggplot2::aes(x = correlation_coefficient, y = -log10(p_value), color = significance_group)) +
  ggplot2::geom_point(alpha = 0.7, size = 2) +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey30") +
  ggplot2::scale_color_manual(
    values = c(
      "Positive Significant" = "#d62728",
      "Negative Significant" = "#31599f",
      "Not Significant" = "grey"
    ),
    breaks = c("Positive Significant", "Negative Significant", "Not Significant")
  ) +
  ggplot2::labs(
    x = "Correlation Coefficient",
    y = bquote(-log[10]~'(P-value)'),
    title = "Volcano Plot",
    color = "Significance"
  ) +
  theme_bw()

plot_data |>
  dplyr::mutate(
    significance_group = dplyr::case_when(
      p_value < 0.05 & correlation_coefficient < 0 ~ "Negative Significant",
      p_value < 0.05 & correlation_coefficient > 0 ~ "Positive Significant",
      TRUE ~ "Not Significant"
    )
  ) |> dplyr::count(significance_group)

topptx(p,'PI_metabolite_volcano.pptx',height = 4,width = 5)


correlation_table_length <- correlation_table |>
  dplyr::filter(phenotype == 'Length.Z') |>
  dplyr::left_join(
    metabolite_selected |> dplyr::select(metabo_id, id_kegg, Standardized.name),
    by = c("metabolite" = "metabo_id")
  )
plot_data <- read.xlsx('length_metabolite_relationship.xlsx',sheet = 'deduplicated')
plot_data_birthlength <- plot_data |> dplyr::filter(p_value< 0.05) |> pull(Standardized.name)
plot_data_birthlength_metabo_id <- plot_data |> dplyr::filter(p_value< 0.05) |> pull(metabolite)
p <- plot_data |>
  dplyr::mutate(
    significance_group = dplyr::case_when(
      p_value < 0.05 & correlation_coefficient < 0 ~ "Negative Significant",
      p_value < 0.05 & correlation_coefficient > 0 ~ "Positive Significant",
      TRUE ~ "Not Significant"
    )
  ) |>
  ggplot2::ggplot(ggplot2::aes(x = correlation_coefficient, y = -log10(p_value), color = significance_group)) +
  ggplot2::geom_point(alpha = 0.7, size = 2) +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey30") +
  ggplot2::scale_color_manual(
    values = c(
      "Positive Significant" = "#d62728",
      "Negative Significant" = "#31599f",
      "Not Significant" = "grey"
    ),
    breaks = c("Positive Significant", "Negative Significant", "Not Significant")
  ) +
  ggplot2::labs(
    x = "Correlation Coefficient",
    y = bquote(-log[10]~'(P-value)'),
    title = "Volcano Plot",
    color = "Significance"
  ) +
  theme_bw()
plot_data |>
  dplyr::mutate(
    significance_group = dplyr::case_when(
      p_value < 0.05 & correlation_coefficient < 0 ~ "Negative Significant",
      p_value < 0.05 & correlation_coefficient > 0 ~ "Positive Significant",
      TRUE ~ "Not Significant"
    )
  ) |> dplyr::count(significance_group)

topptx(p,'length_metabolite_volcano.pptx',height = 4,width = 5)


correlation_table_GA <- correlation_table |>
  dplyr::filter(phenotype == 'GA_age') |>
  dplyr::left_join(
    metabolite_selected |> dplyr::select(metabo_id, id_kegg, Standardized.name),
    by = c("metabolite" = "metabo_id")
  )
plot_data <- read.xlsx('GA_metabolite_relationship.xlsx',sheet = 'deduplicated')
plot_data_GA <- plot_data |> dplyr::filter(p_value< 0.05) |> pull(Standardized.name)
plot_data_GA_metabo_id <- plot_data |> dplyr::filter(p_value< 0.05) |> pull(metabolite)
p <- plot_data |>
  dplyr::mutate(
    significance_group = dplyr::case_when(
      p_value < 0.05 & correlation_coefficient < 0 ~ "Negative Significant",
      p_value < 0.05 & correlation_coefficient > 0 ~ "Positive Significant",
      TRUE ~ "Not Significant"
    )
  ) |>
  ggplot2::ggplot(ggplot2::aes(x = correlation_coefficient, y = -log10(p_value), color = significance_group)) +
  ggplot2::geom_point(alpha = 0.7, size = 2) +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey30") +
  ggplot2::scale_color_manual(
    values = c(
      "Positive Significant" = "#d62728",
      "Negative Significant" = "#31599f",
      "Not Significant" = "grey"
    ),
    breaks = c("Positive Significant", "Negative Significant", "Not Significant")
  ) +
  ggplot2::labs(
    x = "Correlation Coefficient",
    y = bquote(-log[10]~'(P-value)'),
    title = "Volcano Plot",
    color = "Significance"
  ) +
  theme_bw()
plot_data |>
  dplyr::mutate(
    significance_group = dplyr::case_when(
      p_value < 0.05 & correlation_coefficient < 0 ~ "Negative Significant",
      p_value < 0.05 & correlation_coefficient > 0 ~ "Positive Significant",
      TRUE ~ "Not Significant"
    )
  ) |> dplyr::count(significance_group)

topptx(p,'GA_metabolite_volcano.pptx',height = 4,width = 5)

upset_data_list <- list(plot_data_birthweight,plot_data_birthlength,
                        plot_data_GA,plot_data_PI,plot_data_BMI
)

max_length <- upset_data_list |>
  purrr::map_int(length) |>
  max()

result_df <- upset_data_list |>
  purrr::map(~`length<-`(.x, max_length)) |>
  dplyr::bind_cols() |>
  dplyr::rename_with(~ paste0("Sample", seq_along(.x)))
write.xlsx(result_df,'metabolite_upset_plot.xlsx')


metabo_id_list <- unique(c(plot_data_birthweight_metabo_id,plot_data_birthlength_metabo_id,plot_data_GA_metabo_id,plot_data_PI_metabo_id,plot_data_BMI_metabo_id))
sig_metabolite_data <- phenotype_metabolite |> dplyr::select(number,metabo_id_list)
write_xlsx(sig_metabolite_data,'significant_metabolites_with_phenotype.xlsx')



reactome_pathway_level1_df <- read_csv('Reactome_pathwayID_to_level1ID_mapping.csv')
birthweight_reactome <- read.delim('birthweight_reactome.tab') |>
  dplyr::mutate(neg_log10_qval = -log10(`q.value`)) |> dplyr::select(pathway,external_id,neg_log10_qval)
birthweight_reactome$outcome <- 'birthweight'
birthlength_reactome <- read.delim('birthlength_reactome.tab') |>
  dplyr::mutate(neg_log10_qval = -log10(`q.value`)) |> dplyr::select(pathway,external_id,neg_log10_qval)
birthlength_reactome$outcome <- 'birthlength'
GA_reactome <- read.delim('GA_reactome.tab')|>
  dplyr::mutate(neg_log10_qval = -log10(`q.value`)) |> dplyr::select(pathway,external_id,neg_log10_qval)
GA_reactome$outcome <- 'GA'
BMI_reactome <- read.delim('BMI_reactome.tab')|>
  dplyr::mutate(neg_log10_qval = -log10(`q.value`)) |> dplyr::select(pathway,external_id,neg_log10_qval)
BMI_reactome$outcome <- 'BMI'
PI_reactome <- read.delim('PI_reactome.tab')|>
  dplyr::mutate(neg_log10_qval = -log10(`q.value`)) |> dplyr::select(pathway,external_id,neg_log10_qval)
PI_reactome$outcome <- 'PI'

reactome_df <- bind_rows(birthweight_reactome,birthlength_reactome, GA_reactome, PI_reactome, BMI_reactome)
reactome_df_with_catrgory <- merge(reactome_df,reactome_pathway_level1_df,by.x = 'external_id',by.y = 'term',all.x = T)
exclude_list <- c('Cell Cycle','Gene expression (Transcription)','Transport of small molecules','Disease','DNA Replication')
reactome_df_with_catrgory_filtered <- reactome_df_with_catrgory |> dplyr::filter(!pathway_name %in% exclude_list)

reactome_df_with_catrgory_filtered2 <- read_csv('metabolite_reactome_pathways_to_filter.csv')
reactome_class_order <- c('Metabolism','Metabolism of RNA','Metabolism of proteins','Developmental Biology','Cell Cycle','Signal Transduction','responses to stimuli','Hemostasis')

plot_data <- reactome_df_with_catrgory_filtered2 |>
  dplyr::rename(Category = pathway_name,Term = pathway) |>
  dplyr::mutate(neg_log10_fdr = neg_log10_qval,
                outcome = factor(outcome, levels = c("birthweight", "birthlength", "GA", "PI", "BMI"))) |>
  dplyr::mutate(Category = factor(Category, levels = reactome_class_order)) |>
  dplyr::arrange(Category, neg_log10_fdr) |>
  dplyr::mutate(Term = forcats::fct_inorder(Term))

p_reactome <- plot_data |>
  ggplot2::ggplot(ggplot2::aes(x = outcome, y = Term)) +
  ggplot2::geom_tile(ggplot2::aes(fill = neg_log10_fdr), color = "white", size = 0.1) +
  ggplot2::scale_fill_gradient(low = "#f8e6eb", high = "#c22147") +
  ggplot2::facet_grid(
    Category ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  ggplot2::labs(
    x = NULL,
    y = NULL,
    fill = "-log10(FDR)"
  ) +
  ggthemes::theme_base() +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(size = 10),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    strip.placement = "outside",
    strip.background = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank(),
    strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 0, size = 11)
  )
p_reactome

library(dplyr)
library(igraph)
library(ggraph)
library(ggthemes)
library(writexl)
library(eoffice)
library(readxl)
proteome_dat <- read_csv('cleaned_data.csv')
module_info <- read_xlsx('WGANA-module.xlsx')
names(proteome_dat)[1] <- 'number'
proteome_dat <- proteome_dat |> dplyr::select(number,module_info$geneSymbol)

metabo_dat <- read_xlsx('sig_metabolites.xlsx')
prote_meta_dat <- merge(proteome_dat,metabo_dat,by = 'number',all.x = T)

protein_cols <- names(proteome_dat) |>
  dplyr::setdiff("number")

metabolite_cols <- names(metabo_dat) |>
  dplyr::setdiff("number")

correlation_results <- tidyr::expand_grid(
  protein = protein_cols,
  metabolite = metabolite_cols
) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    test_output = list(stats::cor.test(prote_meta_dat[[protein]], prote_meta_dat[[metabolite]], method = "pearson")),
    Correlation = test_output$estimate,
    PValue = test_output$p.value
  ) |>
  dplyr::ungroup() |>
  dplyr::select(
    Protein = protein,
    Metabolite = metabolite,
    Correlation,
    PValue
  ) |>
  dplyr::arrange(PValue)

correlation_results_sig <- correlation_results |>
  dplyr::group_by(Protein) |>
  dplyr::mutate(FDR = p.adjust(PValue,method = 'fdr')) |>
  ungroup() |> dplyr::filter(FDR<0.01)

correlation_results_sig_gene <- merge(correlation_results_sig,module_info |> dplyr::select(geneSymbol,Gene),by.x = 'Protein',by.y = 'geneSymbol',all.x = T)
correlation_results_sig_gene <- read_xlsx('correlation_results_sig_gene_e2.xlsx')

royalblue_string_ppi <- read.delim('royalblue_string_interactions_short.tsv')
graph_full <- royalblue_string_ppi |>
  dplyr::select(from = X.node1, to = node2, weight = combined_score) |>
  dplyr::mutate(class = 'Protein') |>
  igraph::graph_from_data_frame(directed = FALSE) |>
  igraph::simplify(
    remove.multiple = TRUE,
    remove.loops = TRUE,
    edge.attr.comb = 'first'
  )

degrees_full <- igraph::degree(graph_full, mode = "all")
vertices_to_keep <- V(graph_full)[degrees_full > 1]
graph_full <- igraph::induced_subgraph(graph_full, vertices_to_keep)
node_table <- igraph::as_data_frame(graph_full, what = "vertices") |>
  tibble::as_tibble()

protein_metabolite_edges <- correlation_results_sig_gene |>
  dplyr::transmute(
    from = Gene,
    to = Metabolite,
    weight = abs(Correlation)/max(abs(Correlation))
  ) |> dplyr::mutate(class = 'Metabolite') |>
  dplyr::filter(from %in% node_table$name)

graph_protein_metabolite <- protein_metabolite_edges |>
  igraph::graph_from_data_frame(directed = FALSE)
graph_combined <- graph_full + graph_protein_metabolite
graph_full <- graph_combined |>
  igraph::simplify(
    remove.multiple = TRUE,
    remove.loops = TRUE,
    edge.attr.comb = "first"
  )

module_info <- read_xlsx('WGANA-module.xlsx')
midnightblue_string_ppi <- read.delim('midnightblue_string_interactions_short.tsv')
graph_full <- midnightblue_string_ppi |>
  dplyr::select(from = X.node1, to = node2, weight = combined_score) |>
  dplyr::mutate(class = 'Protein') |>
  igraph::graph_from_data_frame(directed = FALSE) |>
  igraph::simplify(
    remove.multiple = TRUE,
    remove.loops = TRUE,
    edge.attr.comb = 'first'
  )

degrees_full <- igraph::degree(graph_full, mode = "all")
vertices_to_keep <- V(graph_full)[degrees_full > 1]
graph_full <- igraph::induced_subgraph(graph_full, vertices_to_keep)
node_table <- igraph::as_data_frame(graph_full, what = "vertices") |>
  tibble::as_tibble()

protein_metabolite_edges <- correlation_results_sig_gene |>
  dplyr::transmute(
    from = Gene,
    to = Metabolite,
    weight = abs(Correlation)/max(abs(Correlation))
  ) |> dplyr::mutate(class = 'Metabolite') |>
  dplyr::filter(from %in% node_table$name)

graph_protein_metabolite <- protein_metabolite_edges |>
  igraph::graph_from_data_frame(directed = FALSE)
graph_combined <- graph_full + graph_protein_metabolite
graph_full <- graph_combined |>
  igraph::simplify(
    remove.multiple = TRUE,
    remove.loops = TRUE,
    edge.attr.comb = "first"
  )



module_info <- read_xlsx('WGANA-module.xlsx')
green_string_ppi <- read.delim('green_string_interactions_short.tsv')

graph_full <- green_string_ppi |>
  dplyr::select(from = X.node1, to = node2, weight = combined_score) |>
  dplyr::mutate(class = 'Protein') |>
  igraph::graph_from_data_frame(directed = FALSE) |>
  igraph::simplify(
    remove.multiple = TRUE,
    remove.loops = TRUE,
    edge.attr.comb = 'first'
  )

degrees_full <- igraph::degree(graph_full, mode = "all")
vertices_to_keep <- V(graph_full)[degrees_full > 1]
graph_full <- igraph::induced_subgraph(graph_full, vertices_to_keep)
node_table <- igraph::as_data_frame(graph_full, what = "vertices") |>
  tibble::as_tibble()

protein_metabolite_edges <- correlation_results_sig_gene |>
  dplyr::transmute(
    from = Gene,
    to = Metabolite,
    weight = abs(Correlation)/max(abs(Correlation))
  ) |> dplyr::mutate(class = 'Metabolite') |>
  dplyr::filter(from %in% node_table$name)

graph_protein_metabolite <- protein_metabolite_edges |>
  igraph::graph_from_data_frame(directed = FALSE)
graph_combined <- graph_full + graph_protein_metabolite
graph_full <- graph_combined |>
  igraph::simplify(
    remove.multiple = TRUE,
    remove.loops = TRUE,
    edge.attr.comb = "first"
  )

degrees_full <- igraph::degree(graph_full, mode = "all")
vertices_to_keep <- V(graph_full)[degrees_full > 1]
graph_obj <- igraph::induced_subgraph(graph_full, vertices_to_keep)
library(tidygraph)
graph_obj <- as_tbl_graph(graph_obj)
graph_obj <- graph_obj %>%
  activate(edges) %>%
  mutate(weight = coalesce(weight_1, weight_2))

graph_obj <- graph_obj %>%
  activate(nodes) %>%
  mutate(degree = centrality_degree(mode = 'all'))

graph_obj <- graph_obj |>
  tidygraph::activate(nodes) |>
  dplyr::mutate(
    class = dplyr::if_else(name %in% node_table$name, "Protein", "Metabolite"),
    degree = tidygraph::centrality_degree(mode = 'all'),
    community = as.factor(tidygraph::group_louvain(weights = .E()$weight))
  )

community_obj <- cluster_louvain(graph_obj, weights = E(graph_obj)$weight)
V(graph_obj)$community <- as.factor(community_obj$membership)

edge_data <- graph_obj %>%
  activate(edges) %>%
  as_tibble()

node_data <- graph_obj %>%
  activate(nodes) %>%
  as_tibble()


layout_fr <- ggraph::create_layout(graph_obj, layout = 'fr', weights = igraph::E(graph_obj)$weight)

set.seed(123)
a <- ggraph::ggraph(layout_fr) +
  ggraph::geom_edge_link(ggplot2::aes(edge_width = weight, alpha = weight), color = "grey") +
  ggraph::scale_edge_width_continuous(range = c(0.3, 1)) +
  ggraph::geom_node_point(ggplot2::aes(size = degree, color = community)) +
  ggraph::geom_node_point(ggplot2::aes(size = degree), shape = 1, color = "black") +
  ggplot2::scale_size_continuous(range = c(2, 10)) +
  ggraph::geom_node_text(ggplot2::aes(label = name, color = community), repel = TRUE, max.overlaps = 15) +
  ggplot2::guides(
    edge_width = "none",
    edge_alpha = "none",
    color = "none",
    size = ggplot2::guide_legend(title = "Degree")
  ) +
  ggraph::theme_graph(base_family = "sans")
a

class_colors <- c("Protein" = "#fb9f93", "Metabolite" = "#b4cbe3")
b <- ggraph::ggraph(layout_fr) +
  ggraph::geom_edge_link(ggplot2::aes(edge_width = weight, alpha = weight), color = "grey") +
  ggraph::scale_edge_width_continuous(range = c(0.3, 1)) +
  ggraph::geom_node_point(ggplot2::aes(size = degree, color = class)) +
  ggraph::geom_node_point(ggplot2::aes(size = degree), shape = 1, color = "black") +
  ggplot2::scale_size_continuous(range = c(2, 10)) +
  ggplot2::scale_color_manual(values = class_colors, name = "Class") +
  ggplot2::guides(
    edge_width = "none",
    edge_alpha = "none",
    color = "none",
    size = ggplot2::guide_legend(title = "Degree")
  ) +
  ggraph::theme_graph(base_family = "sans")+
  theme(legend.position = "none")
b

node_community_table <- dplyr::tibble(
  node = V(graph_obj)$name,
  community = V(graph_obj)$community,
  degree = V(graph_obj)$degree
) |>
  dplyr::arrange(community, dplyr::desc(degree))
node_community_table2 <- merge(node_community_table,module_info,by.x = 'node',by.y = 'Gene',all.x = T)
node_community_table3 <- node_community_table2 |> dplyr::select(node,community,degree,geneSymbol)
