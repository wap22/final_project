# work up motif in R
# Pair up equivalently marked regions
Bothon <- read.delim("both.motifs", header=FALSE)
BothonSyn <- read.delim("bothsyn.motifs", header=FALSE)
names(Bothon)[4]<-"motif"
names(Bothon)[6]<-"L_score"
names(BothonSyn)[4]<-"motif"
names(BothonSyn)[6]<-"S_score"
Bothon$motif <- ifelse(Bothon$motif == 'RAR:RXR(NR),DR5/ES-RAR-ChIP-Seq(GSE56893)/Homer' & Bothon$V3 == 'RGGTCADNNAGAGGTCAV', 'RAR:RXR(NR),DR5/ES-RAR-ChIP-Seq(GSE56893)/Homer+2', as.character(Bothon$motif))
BothonSyn$motif <- ifelse(BothonSyn$motif == 'RAR:RXR(NR),DR5/ES-RAR-ChIP-Seq(GSE56893)/Homer' & BothonSyn$V3 == 'RGGTCADNNAGAGGTCAV', 'RAR:RXR(NR),DR5/ES-RAR-ChIP-Seq(GSE56893)/Homer+2', as.character(BothonSyn$motif))
Bothon <- Bothon[order(Bothon$motif, Bothon$V1),]
BothonSyn <- BothonSyn[order(BothonSyn$motif, BothonSyn$V1),]
Bothon$S_score <- BothonSyn$S_score
# Pair up marked and unmarked scores
Diffon <- read.delim("diff.motifs", header=FALSE)
DiffonSyn <- read.delim("diffsyn.motifs", header=FALSE)
names(Diffon)[4]<-"motif"
names(Diffon)[6]<-"L_score"
names(DiffonSyn)[4]<-"motif"
names(DiffonSyn)[6]<-"S_score"
DiffonSyn$motif <- ifelse(DiffonSyn$motif == 'RAR:RXR(NR),DR5/ES-RAR-ChIP-Seq(GSE56893)/Homer' & DiffonSyn$V3 == 'RGGTCADNNAGAGGTCAV', 'RAR:RXR(NR),DR5/ES-RAR-ChIP-Seq(GSE56893)/Homer+2', as.character(DiffonSyn$motif))
Diffon$motif <- ifelse(Diffon$motif == 'RAR:RXR(NR),DR5/ES-RAR-ChIP-Seq(GSE56893)/Homer' & Diffon$V3 == 'RGGTCADNNAGAGGTCAV', 'RAR:RXR(NR),DR5/ES-RAR-ChIP-Seq(GSE56893)/Homer+2', as.character(Diffon$motif))
DiffonSyn <- DiffonSyn[order(DiffonSyn$motif, DiffonSyn$V1),]
Diffon <- Diffon[order(Diffon$motif, Diffon$V1),]
Diffon$S_score <- DiffonSyn$S_score
library(dplyr)
library(kSamples)
# Run paired t test and save p value
t_BL <- Bothon %>% group_by(motif) %>% do(w = t.test(.$L_score, .$S_score, paired=TRUE)) %>% summarise(motif, T.BL = w$p.value)
t_L <- Diffon %>% group_by(motif) %>% do(w = t.test(.$L_score, .$S_score, paired=TRUE)) %>% summarise(motif, T.L = w$p.value)
t_total <- merge(t_BL, t_L, by = "motif")
# multiple test correction adjustment
t_total$T.BL.adj <- p.adjust(t_total$T.BL, "fdr")
t_total$T.L.adj <- p.adjust(t_total$T.L, "fdr")
# Print plots as pdf
X11()
plot(-log10(t_total[order(t_total$T.L.adj),]$T.L.adj), type = "h", xlim = c(0,400), ylim = c(0,25), ylab = "Motif Enrichment", axes = FALSE, main = "Motif enrichment in differentially marked enhancers")
axis(1, cex.axis = 1.5)
axis(2, cex.axis = 1.5)
dev.copy(pdf, "motif_full.pdf", width=7, height=7)
dev.off()
plot(-log10(t_total[order(t_total$T.L.adj),]$T.L.adj), type = "h", xlim = c(0,20), ylim = c(17,25), ylab = "Motif Enrichment", axes = FALSE, main = "Motif enrichment in differentially marked enhancers")
axis(1, cex.axis = 1.5)
axis(2, cex.axis = 1.5, at = c(17,19,21,23,25))
dev.copy(pdf, "motif_top.pdf", width=7, height=7)
dev.off()
# save output
write.table(t_total, "output.txt", sep="\t", row.names=FALSE, quote=FALSE)
