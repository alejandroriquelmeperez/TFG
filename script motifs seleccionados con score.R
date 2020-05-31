#librerias
library(universalmotif)
library(heatmaply)
library(Biostrings)
#Cálculo raw
setwd("/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_raw/homerResults/")
ruta_raw <- "/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_raw/homerResults/"
files_raw <- list.files(".")
  
scores_raw <- c() #vector que tendra los ficheros y las scores
file_raw <- c() #vector con nombre del fichero
score_raw <- c() #vector con las scores por fichero

for(i in 1:length(files_raw)){
  file_raw <- append(file_raw, files_raw[i])
  ruta_file_raw <- paste0(ruta_raw, files_raw[i]) 
  homer_raw <- read_homer(ruta_file_raw)
  valor_raw <- -(log10(homer_raw@pval))*log10(length(files_raw))
  score_raw <- append(score_raw, valor_raw)
}

scores_raw$file_raw <- file_raw
scores_raw$score_raw <- score_raw
scores_raw <- as.data.frame(scores_raw)
scores_raw <- head(scores_raw[order(-score_raw),],8)

#Cálculo 50
setwd("/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_50/homerResults/")
ruta_50 <- "/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_50/homerResults/"
files_50 <- list.files(".")

scores_50 <- c() #vector que tendra los ficheros y las scores
file_50 <- c() #vector con nombre del fichero
score_50 <- c() #vector con las scores por fichero

for(i in 1:length(files_50)){
  file_50 <- append(file_50, files_50[i])
  ruta_file_50 <- paste0(ruta_50, files_50[i]) 
  homer_50 <- read_homer(ruta_file_50)
  valor_50 <- -(log10(homer_50@pval))*log10(length(files_50))
  score_50 <- append(score_50, valor_50)
}

scores_50$file_50 <- file_50
scores_50$score_50 <- score_50
scores_50 <- as.data.frame(scores_50)
scores_50 <- head(scores_50[order(-score_50),],8)

#Cálculo 100
setwd("/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_100/homerResults/")
ruta_100 <- "/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_100/homerResults/"
files_100 <- list.files(".")

scores_100 <- c() #vector que tendra los ficheros y las scores
file_100 <- c() #vector con nombre del fichero
score_100 <- c() #vector con las scores por fichero

for(i in 1:length(files_100)){
  file_100 <- append(file_100, files_100[i])
  ruta_file_100 <- paste0(ruta_100, files_100[i]) 
  homer_100 <- read_homer(ruta_file_100)
  valor_100 <- -(log10(homer_100@pval))*log10(length(files_100))
  score_100 <- append(score_100, valor_100)
}

scores_100$file_100 <- file_100
scores_100$score_100 <- score_100
scores_100 <- as.data.frame(scores_100)
scores_100 <- head(scores_100[order(-score_100),],8)

#Cálculo 200
setwd("/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_200/homerResults/")
ruta_200 <- "/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_200/homerResults/"
files_200 <- list.files(".")

scores_200 <- c() #vector que tendra los ficheros y las scores
file_200 <- c() #vector con nombre del fichero
score_200 <- c() #vector con las scores por fichero

for(i in 1:length(files_200)){
  file_200 <- append(file_200, files_200[i])
  ruta_file_200 <- paste0(ruta_200, files_200[i]) 
  homer_200 <- read_homer(ruta_file_200)
  valor_200 <- -(log10(homer_200@pval))*log10(length(files_200))
  score_200 <- append(score_200, valor_200)
}

scores_200$file_200 <- file_200
scores_200$score_200 <- score_200
scores_200 <- as.data.frame(scores_200)
scores_200 <- head(scores_200[order(-score_200),],8)

#motifs seleccionados con mejor score

motif1_raw <- read_homer("/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_raw/homerResults/motif6.motif")
motif1_raw@name <- "motif1_raw"
motif2_raw <- read_homer("/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_raw/homerResults/motif2.motif")
motif2_raw@name <- "motif2_raw"
motif3_raw <- read_homer("/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_raw/homerResults/motif15.motif")
motif3_raw@name <- "motif3_raw"
motif4_raw <- read_homer("/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_raw/homerResults/motif1.motif")
motif4_raw@name <- "motif4_raw"
motif5_raw <- read_homer("/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_raw/homerResults/motif11.motif")
motif5_raw@name <- "motif5_raw"
motif6_raw <- read_homer("/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_raw/homerResults/motif14.motif")
motif6_raw@name <- "motif6_raw"
motif7_raw <- read_homer("/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_raw/homerResults/motif8.motif")
motif7_raw@name <- "motif7_raw"
motif8_raw <- read_homer("/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_raw/homerResults/motif5.motif")
motif8_raw@name <- "motif8_raw"

motif1_50 <- read_homer("/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_50/homerResults/motif8.motif")
motif1_50@name <- "motif1_50"
motif2_50 <- read_homer("/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_50/homerResults/motif6.motif")
motif2_50@name <- "motif2_50"
motif3_50 <- read_homer("/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_50/homerResults/motif7.motif")
motif3_50@name <- "motif3_50"
motif4_50 <- read_homer("/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_50/homerResults/motif3.motif")
motif4_50@name <- "motif4_50"
motif5_50 <- read_homer("/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_50/homerResults/motif4.motif")
motif5_50@name <- "motif5_50"
motif6_50 <- read_homer("/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_50/homerResults/motif5.motif")
motif6_50@name <- "motif6_50"
motif7_50 <- read_homer("/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_50/homerResults/motif1.motif")
motif7_50@name <- "motif7_50"
motif8_50 <- read_homer("/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_50/homerResults/motif2.motif")
motif8_50@name <- "motif8_50"

motif1_100 <- read_homer("/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_100/homerResults/motif8.motif")
motif1_100@name <- "motif1_100"
motif2_100 <- read_homer("/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_100/homerResults/motif2.motif")
motif2_100@name <- "motif2_100"
motif3_100 <- read_homer("/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_100/homerResults/motif6.motif")
motif3_100@name <- "motif3_100"
motif4_100 <- read_homer("/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_100/homerResults/motif3.motif")
motif4_100@name <- "motif4_100"
motif5_100 <- read_homer("/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_100/homerResults/motif4.motif")
motif5_100@name <- "motif5_100"
motif6_100 <- read_homer("/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_100/homerResults/motif5.motif")
motif6_100@name <- "motif6_100"
motif7_100 <- read_homer("/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_100/homerResults/motif7.motif")
motif7_100@name <- "motif7_100"
motif8_100 <- read_homer("/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_100/homerResults/motif1.motif")
motif8_100@name <- "motif8_100"

motif1_200 <- read_homer("/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_200/homerResults/motif28.motif")
motif1_200@name <- "motif1_200"
motif2_200 <- read_homer("/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_200/homerResults/motif29.motif")
motif2_200@name <- "motif2_200"
motif3_200 <- read_homer("/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_200/homerResults/motif26.motif")
motif3_200@name <- "motif3_200"
motif4_200 <- read_homer("/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_200/homerResults/motif30.motif")
motif4_200@name <- "motif4_200"
motif5_200 <- read_homer("/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_200/homerResults/motif5.motif")
motif5_200@name <- "motif5_200"
motif6_200 <- read_homer("/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_200/homerResults/motif3.motif")
motif6_200@name <- "motif6_200"
motif7_200 <- read_homer("/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_200/homerResults/motif32.motif")
motif7_200@name <- "motif7_200"
motif8_200 <- read_homer("/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_200/homerResults/motif14.motif")
motif8_200@name <- "motif8_200"

#Stage 4: Comparacion de los motifs
motifs_compared <- compare_motifs(list(motif1_raw, motif1_50, motif1_100, motif1_200, 
                                    motif2_raw, motif2_50, motif2_100, motif2_200, 
                                    motif3_raw, motif3_50, motif3_100, motif3_200,
                                    motif4_raw, motif4_50, motif4_100, motif4_200,
                                    motif5_raw, motif5_50, motif5_100, motif5_200, 
                                    motif6_raw, motif6_50, motif6_100, motif6_200,
                                    motif7_raw, motif7_50, motif7_100, motif7_200,
                                    motif8_raw, motif8_50, motif8_100, motif8_200),
                                  use.type = "ICM", method = "PCC", tryRC = T)
motifs_compared <- as.data.frame(motifs_compared)
#Stage 4: Representación ggheatmap
ggheatmap(motifs_compared, xaxis_font_size = 6)
motif_cluster<- view_motifs(list(motif2_raw, motif7_50, motif8_100, motif6_200))
plot(motif_cluster)
#Stage 5: escaneo de secuencias del cluster seleccionado
seq_raw <- read.table("/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_raw/seq_raw_solo_secuencia")
seq_raw <- as.matrix.Vector(seq_raw)
seq_raw <- DNAStringSet(seq_raw)
seq_50 <- read.table("/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_50/seq_50_solo_secuencia")
seq_50 <- as.matrix.Vector(seq_50)
seq_50 <- DNAStringSet(seq_50)
seq_100 <- read.table("/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_100/seq_100_solo_secuencia")
seq_100 <- as.matrix.Vector(seq_100)
seq_100 <- DNAStringSet(seq_100)
seq_200 <- read.table("/home/alejandro/TFG/RELA/1H/Homer/RELA_1H_200/seq_200_solo_secuencia")
seq_200 <- as.matrix.Vector(seq_200)
seq_200 <- DNAStringSet(seq_200)

seq_raw_scanned <- scan_sequences(motif2_raw, seq_raw, threshold = 0.80, threshold.type = "logodds", RC= T)
seq_50_scanned <- scan_sequences(motif7_50, seq_50, threshold = 0.80, threshold.type = "logodds", RC= T)
seq_100_scanned <- scan_sequences(motif8_100, seq_100, threshold = 0.80, threshold.type = "logodds", RC= T)
seq_200_scanned <- scan_sequences(motif6_200, seq_200, threshold = 0.80, threshold.type = "logodds", RC= T)

#Stage 5: crear motifs nuevos refinados
bkg <- readRDS("/home/alejandro/Homer/bkg.rds")
motif_raw_refined <- create_motif(seq_raw_scanned@listData[["match"]], alphabet = "DNA", type = "PCC", name = "motif_raw_refined", pseudocount = 1/3, bkg)
motif_raw_refined_enriched <- enrich_motifs(motif_raw_refined, seq_raw, qval.method = "BH", threshold = 0.80, threshold.type = "logodds", RC = T)
motif_raw_refined@name <- "motif_raw_refined"
motif_raw_refined@pval <- motif_raw_refined_enriched@listData[["Pval"]]
motif_raw_refined@qval <- motif_raw_refined_enriched@listData[["Qval"]]
motif_raw_refined@eval <- motif_raw_refined_enriched@listData[["Eval"]]
score_raw_refined <- -(log10(motif_raw_refined@pval))*log10(motif_raw_refined_enriched@listData$target.hits)
motif_raw_refined@extrainfo <- paste0("target:", motif_raw_refined_enriched@listData$target.hits, '/', motif_raw_refined_enriched@listData$target.seq.hits, '/', motif_raw_refined_enriched@listData$target.seq.count, ' bkg:', 
                      motif_raw_refined_enriched@listData$bkg.hits, '/', motif_raw_refined_enriched@listData$bkg.seq.hits, '/', motif_raw_refined_enriched@listData$bkg.seq.count, " score: ", score_raw_refined)

motif_50_refined <- create_motif(seq_50_scanned@listData[["match"]], alphabet = "DNA", type = "PCC", name = "motif_50_refined", pseudocount = 1/3, bkg)
motif_50_refined_enriched <- enrich_motifs(motif_50_refined, seq_50, qval.method = "BH", threshold = 0.80, threshold.type = "logodds", RC = T)
motif_50_refined@name <- "motif_50_refined"
motif_50_refined@pval <- motif_50_refined_enriched@listData[["Pval"]]
motif_50_refined@qval <- motif_50_refined_enriched@listData[["Qval"]]
motif_50_refined@eval <- motif_50_refined_enriched@listData[["Eval"]]
score_50_refined <- -(log10(motif_50_refined@pval))*log10(motif_50_refined_enriched@listData$target.hits)
motif_50_refined@extrainfo <- paste0("target:", motif_50_refined_enriched@listData$target.hits, '/', motif_50_refined_enriched@listData$target.seq.hits, '/', motif_50_refined_enriched@listData$target.seq.count, ' bkg:', 
                                      motif_50_refined_enriched@listData$bkg.hits, '/', motif_50_refined_enriched@listData$bkg.seq.hits, '/', motif_50_refined_enriched@listData$bkg.seq.count, " score: ", score_50_refined)

motif_100_refined <- create_motif(seq_100_scanned@listData[["match"]], alphabet = "DNA", type = "PCC", name = "motif_100_refined", pseudocount = 1/3, bkg)
motif_100_refined_enriched <- enrich_motifs(motif_100_refined, seq_100, qval.method = "BH", threshold = 0.80, threshold.type = "logodds", RC = T)
motif_100_refined@name <- "motif_100_refined"
motif_100_refined@pval <- motif_100_refined_enriched@listData[["Pval"]]
motif_100_refined@qval <- motif_100_refined_enriched@listData[["Qval"]]
motif_100_refined@eval <- motif_100_refined_enriched@listData[["Eval"]]
score_100_refined <- -(log10(motif_100_refined@pval))*log10(motif_100_refined_enriched@listData$target.hits)
motif_100_refined@extrainfo <- paste0("target:", motif_100_refined_enriched@listData$target.hits, '/', motif_100_refined_enriched@listData$target.seq.hits, '/', motif_100_refined_enriched@listData$target.seq.count, ' bkg:', 
                                      motif_100_refined_enriched@listData$bkg.hits, '/', motif_100_refined_enriched@listData$bkg.seq.hits, '/', motif_100_refined_enriched@listData$bkg.seq.count, " score: ", score_100_refined)

motif_200_refined <- create_motif(seq_200_scanned@listData[["match"]], alphabet = "DNA", type = "PCC", name = "motif_200_refined", pseudocount = 1/3, bkg)
motif_200_refined_enriched <- enrich_motifs(motif_200_refined, seq_200, qval.method = "BH", threshold = 0.80, threshold.type = "logodds", RC = T)
motif_200_refined@name <- "motif_200_refined"
motif_200_refined@pval <- motif_200_refined_enriched@listData[["Pval"]]
motif_200_refined@qval <- motif_200_refined_enriched@listData[["Qval"]]
motif_200_refined@eval <- motif_200_refined_enriched@listData[["Eval"]]
score_200_refined <- -(log10(motif_200_refined@pval))*log10(motif_200_refined_enriched@listData$target.hits)
motif_200_refined@extrainfo <- paste0("target:", motif_200_refined_enriched@listData$target.hits, '/', motif_200_refined_enriched@listData$target.seq.hits, '/', motif_200_refined_enriched@listData$target.seq.count, ' bkg:', 
                                      motif_200_refined_enriched@listData$bkg.hits, '/', motif_200_refined_enriched@listData$bkg.seq.hits, '/', motif_200_refined_enriched@listData$bkg.seq.count, " score: ", score_200_refined)
write_motifs(motif_raw_refined, file = "/home/alejandro/TFG/RELA/1H/Universalmotif/motif_raw_refined")
write_motifs(motif_50_refined, file = "/home/alejandro/TFG/RELA/1H/Universalmotif/motif_50_refined")
write_motifs(motif_100_refined, file = "/home/alejandro/TFG/RELA/1H/Universalmotif/motif_100_refined")
write_motifs(motif_200_refined, file = "/home/alejandro/TFG/RELA/1H/Universalmotif/motif_200_refined")
#Stage 6: juntar todos los motifs y eliminar los extremos con poca informacion
motif_merged <- merge_motifs(list(motif_raw_refined, motif_50_refined, motif_100_refined, motif_200_refined), method = "PCC", use.type = "PPM", tryRC = F)
motif_merged@name <- "motif_merged"
motif_merged <- trim_motifs(motif_merged, min.ic = 0.25)

#Stage 6: refinar motif definitivo
motif_definitivo <- scan_sequences(motif_merged, seq_raw, threshold = 0.80, threshold.type = "logodds", RC= T)
bkg <- readRDS("/home/alejandro/Homer/bkg.rds")
motif_definitivo <- create_motif(motif_definitivo@listData[["match"]], alphabet = "DNA", type = "PCC", name = "motif_raw_refined", pseudocount = 1/3, bkg)
motif_definitivo_enriched <- enrich_motifs(motif_definitivo, seq_raw, qval.method = "BH", threshold = 0.80, threshold.type = "logodds", RC = T)
motif_definitivo@name <- "motif_merged"
motif_definitivo@pval <- motif_definitivo_enriched@listData[["Pval"]]
motif_definitivo@qval <- motif_definitivo_enriched@listData[["Qval"]]
motif_definitivo@eval <- motif_definitivo_enriched@listData[["Eval"]]
score_definitivo <- -(log10(motif_definitivo@pval))*log10(motif_definitivo_enriched@listData$target.hits)
motif_definitivo@extrainfo <- paste0("target:", motif_definitivo_enriched@listData$target.hits, '/', motif_definitivo_enriched@listData$target.seq.hits, '/', motif_definitivo_enriched@listData$target.seq.count, ' bkg:', 
                                      motif_definitivo_enriched@listData$bkg.hits, '/', motif_definitivo_enriched@listData$bkg.seq.hits, '/', motif_definitivo_enriched@listData$bkg.seq.count, " score: ", score_definitivo)
#representar alineamiento:
logo_motifs_refined <- view_motifs(list(motif_raw_refined, motif_50_refined, motif_100_refined, motif_200_refined, motif_definitivo), use.type = "PPM")
plot(logo_motifs_refined)
logo_motif_definitivo <- view_motifs(motif_definitivo, use.type = "PPM")
plot (logo_motif_definitivo)
write_motifs(motif_definitivo, file = "/home/alejandro/TFG/RELA/1H/Universalmotif/motif_definitivo")
write_jaspar(motif_definitivo, file = "/home/alejandro/TFG/RELA/1H/Universalmotif/motif_definitivo_jaspar")
write_jaspar(motif_raw_refined, file = "/home/alejandro/TFG/RELA/1H/Universalmotif/motif_raw_refined_jaspar")
write_jaspar(motif2_raw, file = "/home/alejandro/TFG/RELA/1H/Universalmotif/motif2_raw_jaspar")
