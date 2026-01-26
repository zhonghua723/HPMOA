library(tidyverse)
library(openxlsx)
library(ggplot2)
library(dplyr)
library(eoffice)
library(ggthemes)
library(pheatmap)
dat2 <- read_csv('cleaned_data') |> as.data.frame()


protein_abundance <- dat2 |>
  tidyr::pivot_longer(
    cols = -number_identifier,
    names_to = "Protein",
    values_to = "Intensity"
  ) |>
  dplyr::mutate(
    Log10_Intensity = log10(Intensity + 1)
  ) |>
  dplyr::group_by(Protein) |>
  dplyr::summarise(
    Average_log10_Intensity = mean(Log10_Intensity, na.rm = TRUE)
  ) |>
  dplyr::arrange(dplyr::desc(Average_log10_Intensity)) |>
  dplyr::mutate(Rank = dplyr::row_number())

p1 <- ggplot(protein_abundance, aes(x = Rank, y = Average_log10_Intensity)) +
  geom_point(col = '#F48892') +
  ggthemes::theme_base() +
  labs(x = "Abundance Rank", y = "MS signal [Log10]")


library(dplyr)
rm(list = ls())

dat1 <- read.xlsx('~/Hela_proteinGroups.xlsx')
dat2 <- dat1 %>% select(Protein.IDs,LFQ.intensity.Hela_1:LFQ.intensity.Hela_9) %>%
  gather(key = "Sample", value = "Intensity", -Protein.IDs) %>%
  mutate(Intensity_log10 = log10(as.numeric(Intensity)))

dat2_cv <- dat2 %>%
  group_by(Protein.IDs) %>%
  summarise(
    Mean_Intensity = mean(Intensity_log10, na.rm = TRUE),
    SD_Intensity = sd(Intensity_log10, na.rm = TRUE)
  ) %>%
  filter(Mean_Intensity > 0) %>%
  mutate(CV = (SD_Intensity / Mean_Intensity) * 100)

dat2_cv <- dat2_cv %>%
  mutate(color = ifelse(CV < 15, "#91CAE8", "#F48892"))

p2 <- ggplot(dat2_cv, aes(x = Mean_Intensity, y = CV, color = color)) +
  geom_point() +
  theme_base() +
  labs(x = "MS signal [Log10]", y = "Coefficients of variation (CV)") +
  scale_color_identity()
p2


mean_cv <- mean(dat2_cv$CV, na.rm = TRUE)
print(mean_cv)

rm(list = ls())
dat1 <- read.xlsx('~/Hela_proteinGroups.xlsx')
dat2 <- dat1 %>% select(LFQ.intensity.Hela_1:LFQ.intensity.Hela_9)
names(dat2) <- paste0("QA",1:14)
head(dat2)
library(pheatmap)

dat2 <- dat2[ , sapply(dat2, is.numeric)]

cor_matrix <- cor(dat2, use = "pairwise.complete.obs")


axis_order <- paste0("QA", 1:14)
cor_matrix_ordered <- cor_matrix[axis_order, axis_order]

p_mat <- purrr::map_dfc(as.data.frame(cor_matrix_ordered), ~ ifelse(abs(.x) >= 0.9, 0, 1)) |>
  as.matrix()
rownames(p_mat) <- rownames(cor_matrix_ordered)


color_palette <- colorRampPalette(c("#91CAE8", "white", "#F48892"))(200)

corrplot::corrplot(
  cor_matrix_ordered,
  method = "circle",
  type = "lower",
  order = "original",
  col = color_palette,
  p.mat = p_mat,
  insig = "blank",
  tl.col = "black",
  tl.srt = 45,
  tl.cex = 1.4,
  diag = T
)

corrplot::corrplot(
  cor_matrix_ordered,
  add = TRUE,
  method = "number",
  type = "upper",
  order = "original",
  p.mat = p_mat,
  insig = "blank",
  number.cex = 1.3,
  number.digits = 2,
  col = color_palette,
  diag = T,
  tl.pos = "n",
  cl.pos = "n"
)

