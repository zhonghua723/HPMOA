library(tidyverse)
library(readxl)
library(openxlsx)
library(WGCNA)
proteome_dat <- read_csv('cleaned_data.csv')
phenotype_dat <- read_csv('phenotype.csv')

phenotype_dat_birthweight <- phenotype_dat |> dplyr::select(number,birth_weight,birth_length,Weight.Z, Length.Z,GA_age) |>
  mutate(Ponderal_Index = (birth_weight/(birth_length ^3))*100,
         BMI = birth_weight/(birth_length/100)^2/1000)

birthweight_cor_dat <- merge(phenotype_dat_birthweight,
                             proteome_dat,
                             by.x = 'number',
                             by.y = 'number_identifier') |>
  na.omit() |>
  group_by(number) |>
  slice_head(n =1) |>
  ungroup()

neonatal_feature_df <- birthweight_cor_dat |> dplyr::select(number,birth_weight,birth_length, GA_age, Ponderal_Index, BMI)
p_box <- neonatal_feature_df |>
  ggplot2::ggplot(ggplot2::aes(x = "birth_length", y = Ponderal_Index)) +
  ggplot2::geom_jitter(width = 0.1, height = 0, alpha = 0.8, color = "#237B9F") +
  ggplot2::geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.5, fill = "#237B9F") +

  ggplot2::labs(
    x = NULL,
    y = NULL
  ) +
  ggthemes::theme_base() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank()
  )
p_box


cor_data <- neonatal_feature_df |>
  dplyr::select(birth_weight, birth_length, GA_age, Ponderal_Index, BMI)

cor_matrix <- cor(cor_data)

library(Hmisc)
res <- rcorr(as.matrix(cor_data))
res$r
res$P

custom_colors <- colorRampPalette(c("#5ba8ff", "white", "#c22147"))(100)

a <- corrplot::corrplot.mixed(
  cor_matrix,
  lower = "circle",
  upper = "number",
  tl.col = "black",
  tl.srt = 45,
  number.cex = 0.8,
  lower.col = custom_colors,
  upper.col = custom_colors
)

M <- birthweight_cor_dat |> dplyr::select(A0AV96:Q96T59) |> as.data.frame()
rownames(M) <- birthweight_cor_dat$number
datExpr0 <- as.data.frame(M)
colnames(datExpr0)
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust(dist(datExpr0), method = "average");
sizeGrWindow(12,9)
par(cex = 0.45);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers",
     sub="", xlab="",
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
abline(h = 100000000, col = "red");

clust = cutreeStatic(sampleTree, cutHeight = 100000000, minSize = 10)
table(clust)
keepSamples = (clust == 1)
datExpr = datExpr0[keepSamples, ]
nGenes =ncol(datExpr)
nSamples = nrow(datExpr)

library(dplyr)
traitData <- birthweight_cor_dat[,2:8] |> as.data.frame()
rownames(traitData) <- birthweight_cor_dat$number
colnames(traitData)
library(stringr)
tumorSamples = rownames(datExpr);
traitRows = match(tumorSamples, rownames(traitData));
datTraits = traitData[traitRows,];
rownames(datTraits) = rownames(traitData)[traitRows];
collectGarbage();

sampleTree2 = hclust(dist(datExpr), method = "average")
datTraitsColor <- numbers2colors(datTraits, signed = FALSE);
sizeGrWindow(12,9)
plotDendroAndColors(sampleTree2, datTraitsColor,
                    groupLabels = names(datTraits),
                    colorHeight = 0.2,
                    colorHeightBase = 0.2,
                    colorHeightMax = 0.4,
                    rowWidths = NULL,
                    dendroLabels = NULL,
                    addGuide = FALSE, guideAll = FALSE,
                    guideCount = 50, guideHang = 0.2,
                    addTextGuide = FALSE,
                    cex.colorLabels = 0.8,
                    cex.dendroLabels = 0.7,
                    cex.rowText = 0.8,
                    marAll = c(1, 5, 3, 1), saveMar = TRUE,
                    main = "Sample dendrogram and trait heatmap")
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
library(patchwork)
plot_data <- sft$fitIndices |>
  tibble::as_tibble() |>
  dplyr::mutate(signed_R2 = -sign(slope) * SFT.R.sq)

p1 <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Power, y = signed_R2)) +
  ggplot2::geom_text(ggplot2::aes(label = Power), color = "red") +
  ggplot2::geom_hline(yintercept = 0.90, color = "red", linetype = "dashed") +
  ggplot2::labs(
    title = "Scale independence",
    x = "Soft Threshold (power)",
    y = "Scale Free Topology Model Fit, signed R^2"
  ) +
  theme_bw()

p2 <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Power, y = mean.k.)) +
  ggplot2::geom_text(ggplot2::aes(label = Power), color = "red") +
  ggplot2::labs(
    title = "Mean connectivity",
    x = "Soft Threshold (power)",
    y = "Mean Connectivity"
  ) +
  theme_bw()

p <- p1 + p2
p
ADJ1=abs(cor(datExpr,use="p"))^4
k=as.vector(apply(ADJ1,2,sum, na.rm=T))
p_hist <- ggplot2::ggplot(data.frame(k = k), ggplot2::aes(x = k)) +
  ggplot2::geom_histogram(bins = 50, fill = "skyblue", color = "black") +
  ggplot2::labs(
    title = "Histogram of Connectivity",
    x = "Connectivity (k)",
    y = "Frequency"
  ) +
  ggthemes::theme_clean()

scale_free_data <- hist(k, breaks = 20, plot = FALSE) |>
  (\(h) tibble::tibble(k = h$mids, p_k = h$density))() |>
  dplyr::filter(p_k > 0 & k > 0) |>
  dplyr::mutate(
    log10_k = log10(k),
    log10_p_k = log10(p_k)
  )

fit <- lm(log10_p_k ~ log10_k, data = scale_free_data)
r_squared <- summary(fit)$r.squared
plot_title <- paste("Check Scale Free Topology, R^2 =", round(r_squared, 2))

p_scale_free <- ggplot2::ggplot(scale_free_data, ggplot2::aes(x = log10_k, y = log10_p_k)) +
  ggplot2::geom_point(color = "black") +
  ggplot2::geom_smooth(method = "lm", se = FALSE, color = "#bf3e37") +
  ggplot2::labs(
    title = plot_title,
    x = "log10(k)",
    y = "log10(p(k))"
  ) +
  theme_bw()

p_scale_free

softPower = 4;
adjacency = adjacency(datExpr, power = softPower)
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average");
sizeGrWindow(12,12)
plot(geneTree, xlab="", sub="", main = "Protein clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
minModuleSize = 30;
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro =
                              FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
sizeGrWindow(8,12)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average");
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.25
abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight =
                            MEDissThres, verbose = 3)
mergedColors = merge$colors;
mergedMEs = merge$newMEs;
sizeGrWindow(12, 9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
moduleColors = mergedColors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
unique_module_color <- unique(mergedColors)

nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


sizeGrWindow(10,8)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(9, 8, 1, 1));
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

Surtime = as.data.frame(datTraits$Weight.Z);
names(Surtime) = "count"
modNames = substring(names(MEs), 3)
geneTraitSignificance = as.data.frame(cor(datExpr, Surtime, use =
                                            "p"));

GSPvalue =
  as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance),
                                 nSamples));
names(geneTraitSignificance) = paste("GS.", names(Surtime), sep="");
names(GSPvalue) = paste("p.GS.", names(Surtime), sep="");

geneInfo0 = data.frame(geneSymbol = rownames(geneTraitSignificance),
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)

table(geneInfo0$moduleColor)
