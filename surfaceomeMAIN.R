# This script compiles all major work done for Stehn et al. 2026
# For further questions, comments, or critiques please email stehn007@umn.edu
# Christopher Stehn
# August 18, 2025

wd <- '~/Desktop/Larg Lab/MPNST Surfaceome/GITHUB UPLOAD/' # CHANGE FOR INDIVIDUAL USER
setwd(wd)

source('./surfaceomeFUNC.R')
library(tidyverse)
library(ComplexHeatmap)
library(org.Hs.eg.db)
library(ggrepel)
library(clusterProfiler)
library(VennDiagram)
library(matrixStats)
library(limma)
library(ggrepel)
select <- dplyr::select # saves time by defining conflicted functions ahead of time
filter <- dplyr::filter
mutate <- dplyr::mutate


## ----- define variables in later use -----
# data files
MPNSTCLfile <- './larga002_wangx933_19253_20230726_EwIM_hcdIT_MUltiConsensus_wQuant_Proteins.txt'
WUpdxFile <- './larga002_wangx933_19513_20240221_EwIM_hcdIT_MultiConsensus_wQuant_Proteins.txt'
MNpdxFile <- './larga002_wangx933_19494_20240212_EnF_hcdIT_MultiConsensus_wQuant_Proteins.txt'
ic50dbFile <- './MPNST_IC50db.xlsx'
QuantibriteFile <- './PE_stats_Quantibrite_20-May-2025.txt'
allLinesPEFile <- './PE_stats_20-May-2025.txt'
CSPA <- './CSPA_human.xlsx'
MPNST.RNA <- './tpmEngineeredMat.txt'
S462TY.pairwiseCor <- './S462TY-SUZ12_pairwiseCor.txt'
# plot names
# plots were created as SVGs for easy editing
MPNST007.v.S462TY.VENN <- './MPNST007.v.S462TY_PRC2inactiveUP.png'
MPNST007volcanoPlot <- './Fig1A.svg'
S462TYvolcanoPlot <- './Fig1B.svg'
RNA.v.protALLPlot <- './Fig2A.svg'
RNA.v.protDecPlot <- './Fig2B.svg'
RNA.v.protFCPlot <- './Fig2C.svg'
protSampleCorPlot <- './Fig2D.svg'
allSampleHeatmapPlot <- './Fig4.svg'
PDXspecFltHeatmapPlot <- './PDXspecFltHeatmap.png'
PreFltTop50HeatmapPlot <- './FigS1.svg'
ic50Plot <- './Fig5B.svg'
ic50.v.QuantPlot <- './Fig5C.svg'
QuantBeadValidPlot <- './FigS3.svg'



## ----- Define consensus surfaceome -----
# commands below are specific to ProteomeDiscoverer output files
# read in files, change column names to match cell lines directly
MPNST <- read_tsv(MPNSTCLfile) %>%
  filter(Master == "IsMasterProtein")
colnames(MPNST)[grepl('^Abundance', colnames(MPNST))] <- c('Ctrl','JH2-002','JH2-079','JH2-103','sNF96.2','ST88-14','STS26T','MPNST007','MP3')
MPNST %>%
  filter(is.na(Ctrl) | `Found in Sample: F1: Sample, A` == 'Peak Found') %>%
  filter(Ctrl < 1000000 | is.na(Ctrl)) -> MPNST
MPNST %>%
  mutate(foundIn = rowSums(!is.na(select(., matches('26T|JH2|sNF|ST88|MPNST')))),
         medIntensity = rowMedians(as.matrix(select(., matches('26T|JH2|sNF|ST88|MPNST'))), na.rm = TRUE)) %>%
  filter(foundIn > 2 & !is.na(`Gene Symbol`)) -> MPNST
WU_PDX <- read_tsv(WUpdxFile) %>%
  filter(Master == "IsMasterProtein")
colnames(WU_PDX)[grepl('^Abundance', colnames(WU_PDX))] <- c('WU_225re', 'WU_386', 'WU_536')
WU_PDX %>%
  mutate(foundIn = rowSums(!is.na(select(., matches('WU')))),
         totInt = rowSums(select(., matches('WU')), na.rm = TRUE)) %>%
  filter(foundIn > 1) -> WU_PDX
MN_PDX <- read_tsv(MNpdxFile) %>%
  filter(Master == "IsMasterProtein")
colnames(MN_PDX)[grepl('^Abundance', colnames(MN_PDX))] <- c('MN_02_1', 'MN_02_2', 'MN_02_3')
# for MN_02, additional step is creating an average (or grouped) abundance column by taking the average intensity of found samples
# since samples foundIn < 2 will be discarded anyway, it is okay to do this en masse
MN_PDX %>%
  mutate(foundIn = rowSums(!is.na(select(., matches('MN_02')))),
         totInt = rowSums(select(., matches('MN_02')), na.rm = TRUE),
         MN_02 = rowMeans(select(., matches('MN_02')), na.rm = TRUE)) %>%
  arrange(-totInt) %>%
  filter(foundIn > 1) %>%
  distinct(`Gene Symbol`, .keep_all = TRUE) -> MN_PDX

PDX <- select(WU_PDX, `Gene Symbol`, matches('WU')) %>%
  left_join(select(MN_PDX, `Gene Symbol`, MN_02), by = 'Gene Symbol', relationship = "many-to-many") %>%
  mutate(foundIn = rowSums(!is.na(select(., matches('WU|MN')))),
         medIntensity = rowMedians(as.matrix(select(., matches('WU|MN'))), na.rm = TRUE),
         totInt = rowSums(select(., matches('WU|MN')))) %>%
  arrange(-totInt) %>%
  filter(foundIn > 2 & !is.na(`Gene Symbol`)) %>%
  distinct(`Gene Symbol`, .keep_all = TRUE) %>%
  select(-totInt)
MPNSTint <- MPNST %>%
  select(`Gene Symbol`, matches('26T|JH2|sNF|ST88|MPNST'), foundIn, medIntensity) %>%
  arrange(-medIntensity)

# need to median normalize numbers given that all of these proteins have
# plot medians of log10 before and after median to see why
PDX %>%
  pivot_longer(cols = matches('WU|MN'), names_to = 'sample', values_to = 'int') %>%
  ggplot(aes(x = sample, y = log10(int), color = sample)) + geom_boxplot()
MPNST %>%
  pivot_longer(cols = matches('26T|JH2|sNF|ST88|MPNST'), names_to = 'sample', values_to = 'int') %>%
  ggplot(aes(x = sample, y = log10(int), color = sample)) + geom_boxplot()
PDXmeds <- c(median(PDX$WU_225re, na.rm = T), median(PDX$WU_386, na.rm = T), median(PDX$WU_536, na.rm = T), median(PDX$MN_02, na.rm = T))
MPNSTmeds <- c(median(MPNSTint$`JH2-002`, na.rm = T), median(MPNSTint$`JH2-079`, na.rm = T), median(MPNSTint$`JH2-103`, na.rm = T), median(MPNSTint$sNF96.2, na.rm = T), median(MPNSTint$`ST88-14`, na.rm = T), median(MPNSTint$STS26T, na.rm = T), median(MPNSTint$MPNST007, na.rm = T))

# use limma normalizebetweenarrays to scale values for quantile before log10 norm
MPNSTint[,2:8] <- log10(limma::normalizeBetweenArrays(MPNSTint[,2:8], method = 'quantile'))
MPNSTint %>%
  pivot_longer(cols = matches('26T|JH2|sNF|ST88|MPNST'), names_to = 'sample', values_to = 'int') %>%
  ggplot(aes(x = sample, y = int, color = sample)) + geom_boxplot()

PDX[,2:5] <- log10(limma::normalizeBetweenArrays(PDX[,2:5], method = 'quantile'))
PDX %>%
  pivot_longer(cols = matches('WU|MN'), names_to = 'sample', values_to = 'int') %>%
  ggplot(aes(x = sample, y = int, color = sample)) + geom_boxplot()


MPNSTint %>%
  select(-medIntensity) %>%
  mutate(meanInt = rowMeans(select(., matches('26T|JH2|sNF|ST88|MPNST')), na.rm = TRUE)) %>%
  arrange(-foundIn, -meanInt) -> MPNSTsort

PDX %>%
  select(-medIntensity) %>%
  mutate(meanInt = rowMeans(select(., matches('WU|MN')), na.rm = TRUE)) %>%
  arrange(-foundIn, -meanInt) -> PDXsort

# we have plenty of proteins that are found in every sample, so subset to just those
MPNSTsort %>% filter(foundIn == max(foundIn)) -> MPNSTenrich
PDXsort %>% filter(foundIn == max(foundIn)) -> PDXenrich

MPNSTenrich %>% slice(1:50) %>% pull(`Gene Symbol`) -> MPNSTtop50
PDXenrich %>% slice(1:50) %>% pull(`Gene Symbol`) -> PDXtop50
top50genes <- c(MPNSTtop50, PDXtop50)
top50genes <- top50genes[!duplicated(top50genes)]

MPNSTenrich %>%
  filter(`Gene Symbol` %in% top50genes) %>%
  select(1:8) -> MPNSTtop50
PDXenrich %>%
  filter(`Gene Symbol` %in% top50genes) %>%
  select(1:5) -> PDXtop50

top50total <- full_join(MPNSTtop50, PDXtop50, by = 'Gene Symbol') %>%
  filter(!duplicated(`Gene Symbol`)) %>%
  column_to_rownames(var = 'Gene Symbol') %>%
  as.matrix()
# plot
range <- range(top50total, na.rm = TRUE)
col_fun <- circlize::colorRamp2(breaks = seq(range[1], range[2], length.out = 11),
                                colors = rev(RColorBrewer::brewer.pal(11, 'RdBu')))
col_split <- c(rep(1, 7), rep(2, 4))
column_ha <- HeatmapAnnotation(Model = c(rep('MPNST', 7), rep('PDX', 4)),
                               col = list(Model = c('MPNST' = 'black', 'PDX' = 'grey')))

png(PreFltTop50HeatmapPlot, width = 1.5625*1200, height = 1200)
ComplexHeatmap::Heatmap(top50total, cluster_rows = FALSE, cluster_columns = FALSE,
                        col = col_fun, row_names_gp = gpar(fontsize = 9),
                        heatmap_legend_param = list(title = expression('log'[10]*' iBAQ')),
                        top_annotation = column_ha, column_split = col_split, column_title = NULL,
                        column_names_gp = gpar(fontsize = 18, fontface = 'italic'), 
                        column_names_rot = 60)
graphics.off()

# now that we see that, let's narrow things down to those that actually cell surface proteins
# then, select for the most highly expressed
cspa <- readxl::read_xlsx(CSPA) %>%
  filter(`CSPA category` != '3 - unspecific')
MPNSTenrich %>%
  filter(`Gene Symbol` %in% cspa$`ENTREZ gene symbol`) -> MPNSTflt1
PDXenrich %>%
  filter(`Gene Symbol` %in% cspa$`ENTREZ gene symbol`) -> PDXflt1

totFlt1 <- inner_join(MPNSTflt1, PDXflt1, by = 'Gene Symbol') %>%
  select(-c('foundIn.x', 'foundIn.y', 'meanInt.x', 'meanInt.y'))
totFlt1mat <- totFlt1 %>% filter(!duplicated(`Gene Symbol`)) %>%
  column_to_rownames(var = 'Gene Symbol') %>% as.matrix()

range <- range(totFlt1mat, na.rm = TRUE)
col_fun <- circlize::colorRamp2(breaks = seq(range[1], range[2], length.out = 11),
                                colors = rev(RColorBrewer::brewer.pal(11, 'RdBu')))
col_split <- c(rep(1, 7), rep(2, 4))
column_ha <- HeatmapAnnotation(Model = c(rep('MPNST', 7), rep('PDX', 4)),
                               col = list(Model = c('MPNST' = 'black', 'PDX' = 'grey')))

png(allSampleHeatmapPlot, width = 1.5625*1200, height = 1200)
ComplexHeatmap::Heatmap(totFlt1mat, cluster_rows = FALSE, cluster_columns = FALSE,
                        col = col_fun, row_names_gp = gpar(fontsize = 9),
                        heatmap_legend_param = list(title = expression('log'[10]*' iBAQ')),
                        top_annotation = column_ha, column_split = col_split, column_title = NULL,
                        column_names_gp = gpar(fontsize = 18, fontface = 'italic'), 
                        column_names_rot = 60)
graphics.off()

sum(MPNSTflt1$`Gene Symbol` %in% PDXflt1$`Gene Symbol`)

PTK7mat <- t(as.matrix(totFlt1mat['PTK7',]))
rownames(PTK7mat) <- 'PTK7'
range <- range(PTK7mat, na.rm = TRUE)
col_fun <- circlize::colorRamp2(breaks = seq(range[1], range[2], length.out = 11),
                                colors = rev(RColorBrewer::brewer.pal(11, 'RdBu')))
png('./fltAllModelsPTK7.png', width = 500, height = 150)
ComplexHeatmap::Heatmap(PTK7mat, cluster_columns = FALSE,
                        col = col_fun, row_names_gp = gpar(fontsize = 9),
                        top_annotation = column_ha, column_split = col_split,
                        heatmap_legend_param = list(title = expression('log'[10]*' LFQ')), 
                        column_title = NULL,
                        column_names_gp = gpar(fontsize = 9, fontface = 'italic'), 
                        column_names_rot = 60)
graphics.off()

# show differences in how many proteins were selected for
myCol <- brewer.pal(3, 'Pastel1')
darkCol <- brewer.pal(3, 'Set1')
venn.diagram(x = list(MPNSTenrich$`Gene Symbol`, PDXenrich$`Gene Symbol`), 
             category.names = c('MPNST', 'PDX'), 
             filename = './MPNST_PDX_surfaceomeOverlapNEW.png', output = FALSE,
             fill = myCol[1:2], fontfamily = 'sans',
             cat.fontfamily = 'sans')

venn.diagram(x = list(MPNSTflt1$`Gene Symbol`, PDXflt1$`Gene Symbol`), 
             category.names = c('MPNST', 'PDX'), 
             filename = './MPNST_PDX_surfaceomeOverlapPOSTFLT.png', output = FALSE,
             fill = myCol[1:2], fontfamily = 'sans',
             cat.fontfamily = 'sans')

PDXflt2 <- PDXflt1 %>%
  filter(`Gene Symbol` %in% MPNSTflt1$`Gene Symbol` == 0)
PDXflt2 %>%
  select(-c('foundIn', 'meanInt')) %>%
  column_to_rownames(var = 'Gene Symbol') %>%
  as.matrix() -> PDXflt2mat

range <- range(PDXflt2mat, na.rm = TRUE)
col_fun <- circlize::colorRamp2(breaks = seq(range[1], range[2], length.out = 11),
                                colors = rev(RColorBrewer::brewer.pal(11, 'RdBu')))
# plot PDX-specific heatmap
png(PDXspecFltHeatmapPlot, width = 1.5625*1200, height = 1200)
ComplexHeatmap::Heatmap(PDXflt2mat, cluster_rows = FALSE, cluster_columns = FALSE,
                        col = col_fun, row_names_gp = gpar(fontsize = 9),
                        heatmap_legend_param = list(title = expression('log'[10]*' LFQ')),
                        column_title = NULL,
                        column_names_gp = gpar(fontsize = 18, fontface = 'italic'), 
                        column_names_rot = 60)
graphics.off()

# compare RNA and protein data in expression
# read in TPM, use regex to remove lines starting with N[0-9] (engineered HSCs)
tpm <- read.table(MPNST.RNA, sep = '\t', row.names = 1) %>%
  rownames_to_column(var = 'id') %>%
  select(-matches('^N[0-9]+.*'))
# change keys of RNA data to match protein data
keys <- tpm$id
annot <- AnnotationDbi::select(org.Hs.eg.db, keys = keys, columns = 'SYMBOL',
                               keytype = 'ENSEMBL') %>%
  filter(!duplicated(ENSEMBL))
tpm <- tpm %>% left_join(annot, by = c('id' = 'ENSEMBL'))
# rearrange columns
tpm <- tpm %>% dplyr::select(14, 2:13)

# now calculate correlations between protein data for samples and the associated RNA (for samples in which it's available)
# remove duplicates from MPNST, first sorting to keep highest expressed
MPNST <- MPNST %>%
  arrange(-medIntensity) %>%
  filter(!duplicated(`Gene Symbol`))
ST88.14 <- MPNST %>%
  dplyr::select(`Gene Symbol`, `ST88-14`) %>%
  left_join(dplyr::select(tpm, SYMBOL, ST88.14.S1), by = c('Gene Symbol' = 'SYMBOL')) %>%
  mutate(id = paste0(`Gene Symbol`, '_ST88-14'))
colnames(ST88.14) <- c('gene', 'protein', 'RNA', 'id')
MP3 <- MPNST %>%
  dplyr::select(`Gene Symbol`, `MP3`) %>%
  left_join(dplyr::select(tpm, SYMBOL, MP3.S1), by = c('Gene Symbol' = 'SYMBOL')) %>%
  mutate(id = paste0(`Gene Symbol`, '_MP3'))
colnames(MP3) <- c('gene', 'protein', 'RNA', 'id')
STS26T <- MPNST %>%
  dplyr::select(`Gene Symbol`, `STS26T`) %>%
  left_join(dplyr::select(tpm, SYMBOL, X26T.S1), by = c('Gene Symbol' = 'SYMBOL')) %>%
  mutate(id = paste0(`Gene Symbol`, '_STS26T'))
colnames(STS26T) <- c('gene', 'protein', 'RNA', 'id')
sNF96.2 <- MPNST %>%
  dplyr::select(`Gene Symbol`, `sNF96.2`) %>%
  left_join(dplyr::select(tpm, SYMBOL, SNF.S1), by = c('Gene Symbol' = 'SYMBOL')) %>%
  mutate(id = paste0(`Gene Symbol`, '_sNF96.2'))
colnames(sNF96.2) <- c('gene', 'protein', 'RNA', 'id')

# create data frame of all four lines
cor <- rbind(MP3, sNF96.2, ST88.14, STS26T) %>%
  filter(!is.na(protein) & !is.na(RNA)) %>%
  arrange(-c(RNA)) %>%
  filter(!duplicated(id)) %>%
  mutate(logRNA = log10(RNA),
         logprot = log10(protein))

RNAcorPlot <- cor %>%
  filter(RNA != 0) %>%
  ggplot(aes(x = logRNA, y = logprot)) + geom_point() +
  xlab(expression('RNA (log'[10]*' TPM)')) + 
  ylab(expression('Protein (log'[10]*' LFQ)')) +
  geom_smooth(method = 'lm', linetype = 'longdash', se = FALSE) +
  geom_text(aes(x = -2, y = 9, label = '\u03c1 = 0.02')) +
  theme_minimal()
ggsave(RNAcorPlot, file = RNA.v.protALLPlot, width = 6, height = 4)

# create deciles and sort
cor %>%
  arrange(-logprot) %>%
  mutate(n = row_number(),
         percentile = 100*n*(1/8288),
         decile = floor(percentile/10) + 1) -> sortedcor
# with the way I calculated deciles, the last row gets put into the "11th decile", change manually
sortedcor$decile[8288] <- 10
# plot decile of protein abundance against paired RNA expression
sortedcorplot <- sortedcor %>%
  filter(RNA != 0) %>%
  mutate(decile = as.factor(decile)) %>%
  ggplot(aes(x = decile, y = logRNA, group = decile, color = decile)) +
  geom_boxplot() + ylab('RNA expression (log(TPM))') +
  xlab('Surface protein abundance (Decile)') +
  theme_bw() + scale_color_brewer(type = 'div', palette = 'RdBu', 
                                  direction = -1) +
  guides(color = 'none')
ggsave(sortedcorplot, file = RNA.v.protDecPlot, width = 6, height = 4)

sortedcor %>%
  filter(RNA != 0 & percentile < 10) -> sortedtopdecilecor

sortedcor %>%
  filter(RNA != 0) %>%
  mutate(decile = floor(percentile)/10) -> sortedcor

rho <- cor.test(x = sortedtopdecilecor$logRNA,
                y = sortedtopdecilecor$logprot)

## ----- PRC2 COMPARISONS -----
# comparisons were undertaken using ProteomeDiscoverer, so script is primarily for plotting results

# sample pairwise correlations were derived from outside software and recorded in .txt, so read in and plot
# keep colnames but change them to be nicer for plotting
mat <- read.table(S462TY.pairwiseCor, header = TRUE)
rownames(mat) <- c(paste0('S462TY DMSO', 1:3), paste0('S462TY DOX', 1:3))
colnames(mat) <- c(paste0('S462TY DMSO', 1:3), paste0('S462TY DOX', 1:3))


# calculate correlation between RNA and protein fold changes in S462TY
fcRNA <- read_tsv('./S462TY_addback_Dox.v.DMSO.txt')
fcProt <- read_tsv('./larga002_wangx933_19345_20230929_LFQ_DMSO_DOX- ANOVA_wImputation_Proteins.txt') %>%
  mutate(lfc = log2(`Abundance Ratio DOX  DMSO`))
fcALL <- select(fcProt, `Gene Symbol`, lfc) %>%
  left_join(select(fcRNA, genename, log2FoldChange), 
            by = c('Gene Symbol' = 'genename')) %>%
  filter(!is.na(lfc) & !is.na(log2FoldChange))
colnames(fcALL) <- c('gene', 'prot', 'RNA')


fcPlot <- fcALL %>%
  filter(prot > -6 & prot < 6) %>%
  ggplot(aes(x = RNA, y = prot)) + geom_point() +
  geom_smooth(method = 'lm', se = FALSE, linetype = 'longdash') +
  theme_minimal() + xlab('RNA Fold Change') +
  ylab('Surface Protein Fold Change') +
  geom_text(aes(x = -3, y = 3, label = '\u03c1 = 0.25'))
ggsave(fcPlot, file = RNA.v.protFCPlot, width = 6, height = 4)

# get overlap of MPNST007 and S462TY PRC2 on/off experiments
fcProtMPNST007 <- read_tsv('./larga002_wangx933_19661_20240722_FAIMS_hcdIT_LFQ_PD30_proteinExport.xlsx - Proteins.tsv') %>%
  mutate(lfc = log2(`Abundance Ratio: (EPZ) / (DMSO)`),
         pval = `Abundance Ratio Adj. P-Value: (EPZ) / (DMSO)`,
         result = case_when(pval < 0.05 & lfc > 1 ~ 'up',
                            pval < 0.05 & lfc < -1 ~ 'down',
                            .default = 'NA'),
         name = `Gene Symbol`) %>%
  separate_rows(name, sep = '; ')
# S462TY: multiply lfc by -1 as directionality of experiments is opposite
# want positive values to correspond to higher abundance when PRC2 is INACTIVE
fcProtS462TY <- read_tsv('./larga002_wangx933_19345_20230929_LFQ_DMSO_DOX- ANOVA_wImputation_Proteins.txt') %>%
  mutate(lfc = -1*log2(`Abundance Ratio DOX  DMSO`),
         pval = `Abundance Ratio Adj P-Value DOX  DMSO`,
         result = case_when(pval < 0.05 & lfc > 1 ~ 'up',
                            pval < 0.05 & lfc < -1 ~ 'down',
                            .default = 'NA'),
         name = `Gene Symbol`) %>%
  separate_rows(name, sep = '; ')

sigUpMPNST007 <- fcProtMPNST007 %>% filter(pval < 0.05 & lfc > 0)
sigUpS462TY <- fcProtS462TY %>% filter(pval < 0.05 & lfc > 0 & Contaminant == 'FALSE')
sigUpNamesMPNST007 <- na.omit(sigUpMPNST007$`Gene Symbol`[!duplicated(sigUpMPNST007$`Gene Symbol`)])
sigUpNamesS462TY <- na.omit(sigUpS462TY$`Gene Symbol`[!duplicated(sigUpS462TY$`Gene Symbol`)])

# plot MPNST007
# red means higher expression when PRC2 is active, blue is when PRC2 is inactive
# use ggrepel to label most differentially abundant proteins (|lfc| > 4)
# label plot with cell line
fcProtMPNST007graph <- fcProtMPNST007 %>%
  ggplot(aes(x = lfc, y = -log10(pval), color = result)) + geom_point() +
  xlab(expression('log'[2]*' FC')) + ylab(expression('-log'[10]*' p-value')) +
  geom_vline(xintercept = c(-1, 1), linetype = 'longdash') + 
  geom_hline(yintercept = -log10(0.05), linetype = 'longdash') +
  theme_bw() + ggtitle('MPNST007') + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c('red', 'black', 'blue')) + guides(color = 'none') +
  geom_label_repel(data = filter(fcProtMPNST007, abs(lfc) > 3.5 & pval < 0.05), 
                   inherit.aes = TRUE, aes(label = name),
                   max.overlaps = 20)
ggsave(fcProtMPNST007graph, file = MPNST007volcanoPlot, width = 6, height = 4)  

# plot S462TY, same format as MPNST007
fcProtS462TYgraph <- fcProtS462TY %>%
  ggplot(aes(x = lfc, y = -log10(pval), color = result)) + geom_point() +
  xlab(expression('log'[2]*' FC')) + ylab(expression('-log'[10]*' p-value')) +
  geom_vline(xintercept = c(-1, 1), linetype = 'longdash') + 
  geom_hline(yintercept = -log10(0.05), linetype = 'longdash') +
  theme_bw() + ggtitle('S462TY') + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c('red', 'black', 'blue')) + guides(color = 'none') +
  geom_label_repel(data = filter(fcProtS462TY, abs(lfc) > 3.5 & pval < 0.05), 
                   inherit.aes = TRUE, aes(label = name),
                   max.overlaps = 20)
ggsave(fcProtS462TYgraph, file = S462TYvolcanoPlot, width = 6, height = 4)  

myCol <- RColorBrewer::brewer.pal(3, 'Pastel1')
# Package VennDiagram uses futile.logger to create logs for all actions by default
# Change threshold to ERROR to ensure log files made only when error occurs
futile.logger::flog.threshold(ERROR, name = 'VennDiagramLogger')
venn.diagram(x = list(sigUpNamesMPNST007, sigUpNamesS462TY), 
             category.names = c('MPNST007', 'S462TY'), 
             filename = MPNST007.v.S462TY.VENN, output = FALSE,
             fill = c(myCol[1], myCol[2]), fontfamily = 'sans')
sigUpNamesS462TY[sigUpNamesS462TY %in% sigUpNamesMPNST007]


## ----- Flow cytometry data -----
quantBeads <- read_tsv(QuantibriteFile) %>%
  dplyr::slice(1) %>%
  as.data.frame() %>%
  column_to_rownames(var = 'Sample:')
# get info from CytoFlex column names (tier and stat)
colnames(quantBeads) <- gsub('Quant Beads/Tier ([0-9]) \\| (.*) \\(.*', 'Tier_\\1_\\2', colnames(quantBeads))
quantBeads <- t(quantBeads) %>%
  as.data.frame() %>%
  rownames_to_column(var = 'metric')
# make df with all our info
# nPE is based on mean PE molecules per bead by quantibrite lot number
# QuantiBRITE PE Quantitation Kit (Lot #96597)
beadsDf <- data.frame(nPE = c(638, 6930, 31008, 67601),
                      level = c('Low', 'Medium-Low', 'Medium-High', 'High'),
                      median = quantBeads$`01-PE quantification beads-A3.fcs`[c(2,1,3,4)],
                      geomean = quantBeads$`01-PE quantification beads-A3.fcs`[c(5:8)]) %>%
  mutate(logMedian = log10(median),
         logGeoMean = log10(geomean),
         logPE = log10(nPE))

beadsFit <- lm(logGeoMean ~ logPE, data = beadsDf)
beadsR2 <- summary(beadsFit)$r.squared
beadsCoef <- summary(beadsFit)$coefficients[1:2,1]

beadsDf %>%
  ggplot(aes(x = logPE, y = logMedian)) + geom_point() +
  geom_smooth(method = 'lm', se = FALSE, color = 'black') +
  theme_bw() + xlab(expression('log'[10]*' PE molecules/bead')) +
  ylab(expression('log'[10]*' median fluorescence')) +
  annotate('text', label = paste0('y = ', round(beadsCoef[2], 3), 'x + ', round(beadsCoef[1], 3)),
           x = 3.5, y = 5.5) +
  annotate('text', label = paste0('R\u00b2 = ', round(beadsR2, 4)),
           x = 3.5, y = 5.35)

abPEr <- data.frame(target = c('PTK7', 'EGFR', 'HER2', 'MET', 'Isotype'),
                    ratio = c(1.445, 1.029, 1.139, 1.430, 1.576))

# read in all data
allLines <- read_tsv(allLinesPEFile) %>%
  dplyr::slice(1:48)
colnames(allLines) <- c('sample', 'median', 'geomean')
wellToSamp <- data.frame(well = c(paste0('B', 1:11), paste0('C', 1:7),
                                  paste0('D', 1:7), paste0('E', 1:7),
                                  paste0('F', 1:7), paste0('G', 1:4),
                                  paste0('H', 1:4)),
                         line = c(rep('S462TY', 7), rep('SK-OV-3', 4),
                                  rep('ST88-14', 7), rep('STS26T', 7),
                                  rep('sNF96.2', 7), rep('JH-2-002', 7),
                                  rep('A549', 4), rep('RD', 4)))
allLines <- allLines %>%
  filter(!grepl('quantification', sample)) %>%
  mutate(logMed = log10(median),
         logMean = log10(geomean),
         well = gsub('01-.*-(.*).fcs', '\\1', sample),
         target = gsub('01-(.*)-.*', '\\1', sample),
         nPE = 10^((logMed - 0.921)/1.022)) %>%
  left_join(wellToSamp, by = 'well') %>%
  left_join(abPEr, by = 'target') %>%
  mutate(nTarget = nPE/ratio)
# given STS26T having bimodal PTK7 PE dist, manually adding these
ptk7addendum <- data.frame(sample = c('01-PTK7-D3-negative.fcs', '01-PTK7-D3-positive.fcs'),
                           median = c(16439, 215754), geomean = c(16315, 212164)) %>%
  mutate(logMed = log10(median),
         logMean = log10(geomean),
         well = 'D3', target = 'PTK7',
         nPE = 10^((logMed - 0.921)/1.022),
         line = paste0('STS26T-', c('negative', 'positive')),
         ratio = 1.445, nTarget = nPE/ratio)
allLines <- bind_rows(allLines, ptk7addendum)

# read in IC50 data, append, and plot
ic50 <- readxl::read_xlsx(ic50dbFile)

allLines <- allLines %>%
  left_join(ic50, by = c('line', 'target'))
# remember though, since we have distinct STS26T populations, we need to manually add those IC50s
allLines$IC50[48:49] <- allLines$IC50[allLines$line == 'STS26T' & allLines$target == 'PTK7']

surfFit <- lm(IC50 ~ nTarget, data = subset(allLines, line %in% c('S462TY', 'ST88-14', 'STS26T', 'sNF96.2', 'JH-2-002', 'STS26T-positive') &
                                              target %in% c('EGFR', 'HER2', 'PTK7', 'MET')))
a <- subset(allLines, line %in% c('S462TY', 'ST88-14', 'STS26T', 'sNF96.2', 'JH-2-002', 'STS26T-positive') &
              target %in% c('EGFR', 'HER2', 'PTK7', 'MET') & !is.na(IC50))
cor.test(a$IC50, a$nTarget, method = 'spearman') -> b

antigenCorPlot <- allLines %>%
  filter(line %in% c('S462TY', 'ST88-14', 'STS26T', 'sNF96.2', 'JH-2-002', 'STS26T-positive') &
           target %in% c('EGFR', 'HER2', 'PTK7', 'MET')) %>%
  ggplot(aes(x = nTarget, y = IC50, color = target)) + geom_point(size = 2) + 
  scale_x_log10() + scale_y_log10() +
  geom_smooth(method = 'lm', aes(color = NULL), color = 'black', 
              linetype = 'longdash', alpha = 0.1) + 
  theme_bw() + labs(x = 'Antigens per Cell', y = expression('IC'[50]*' (nM)'),
                    color = 'Target') +
  annotate('text', x = 20000, y = 300, label = paste0('\u03c1 = ', round(b$estimate, 3))) +
  annotate('text', x = 20000, y = 175, label = paste0('p < 0.05')) +
  scale_color_manual(values = c('#00D600', '#FF7F00','#FF4848','#00D6D6'))
ggsave(antigenCorPlot,file = ic50.v.QuantPlot)

# replot the original IC50 boxplots
ic50boxplot <- ic50 %>%
  ggplot(aes(x = target, y = IC50, fill = target)) + geom_boxplot() +
  scale_y_log10() + theme_bw() + labs(x = 'Target', y = expression('IC'[50]*' (nM)')) +
  scale_fill_manual(values = c('pink', '#00D600', '#FF7F00', '#FF4848', '#00D6D6', 'mediumpurple')) +
  guides(fill = 'none') + geom_point()
ggsave(ic50boxplot, file = ic50Plot)
