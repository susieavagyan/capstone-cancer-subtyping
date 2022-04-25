# Loading requirements
library(argparser, quietly=TRUE)
library(oposSOM)
library(pals)
library(fdrtool)

#modified functions of the package 
pipeline.summarySheetsSamples.2 <- function (env) 
{
  if (ncol(env$indata) >= 1000) 
    return()
  dir.create(env$output.paths["Summary Sheets Samples"], showWarnings = FALSE)
  ylim.max <- 0
  for (m in 1:ncol(env$indata)) {
    h <- hist(env$p.g.m[, m], bre = 20, plot = FALSE)
    y.max <- max(h$density)
    if (y.max > ylim.max) {
      ylim.max <- y.max
    }
  }
  util.info("Writing:", file.path(env$output.paths["Summary Sheets Samples"], 
                                  "*.pdf"))
  for (m in 1:ncol(env$indata)) {
    basename <- paste(make.names(make.unique(colnames(env$indata))[m]), 
                      ".pdf", sep = "")
    pdf(file.path(env$output.paths["Summary Sheets Samples"], 
                  basename), 29.7/2.54, 21/2.54, useDingbats = FALSE)
    layout(matrix(c(1, 2, 0, 1, 3, 4, 5, 5, 6, 7, 7, 8), 
                  3, 4), widths = c(1, 1, 2, 2), heights = c(2, 1, 
                                                             1))
    par(mar = c(0, 0, 0, 0))
    plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", 
         xlim = c(0, 1), ylim = c(0, 1))
    text(0.1, 0.94, colnames(env$indata)[m], cex = 3, adj = 0)
    text(0.1, 0.8, "Global Summary", cex = 1.8, adj = 0)
    text(0.1, 0.7, paste("%DE =", round(env$perc.DE.m[colnames(env$indata)[m]], 
                                        2)), adj = 0)
    all.fdr.genes <- which(env$fdr.g.m[, m] < 0.2)
    plus.fdr.genes <- which(env$indata[all.fdr.genes, m] > 
                              0)
    minus.fdr.genes <- which(env$indata[all.fdr.genes, m] <= 
                               0)
    text(0.1, 0.65, paste("# genes with fdr < 0.2  =", length(all.fdr.genes), 
                          " (", length(plus.fdr.genes), "+ /", length(minus.fdr.genes), 
                          " -)"), adj = 0)
    all.fdr.genes <- which(env$fdr.g.m[, m] < 0.1)
    plus.fdr.genes <- which(env$indata[all.fdr.genes, m] > 
                              0)
    minus.fdr.genes <- which(env$indata[all.fdr.genes, m] <= 
                               0)
    text(0.1, 0.6, paste("# genes with fdr < 0.1  =", length(all.fdr.genes), 
                         " (", length(plus.fdr.genes), "+ /", length(minus.fdr.genes), 
                         " -)"), adj = 0)
    all.fdr.genes <- which(env$fdr.g.m[, m] < 0.05)
    plus.fdr.genes <- which(env$indata[all.fdr.genes, m] > 
                              0)
    minus.fdr.genes <- which(env$indata[all.fdr.genes, m] <= 
                               0)
    text(0.1, 0.55, paste("# genes with fdr < 0.05  =", 
                          length(all.fdr.genes), " (", length(plus.fdr.genes), 
                          "+ /", length(minus.fdr.genes), " -)"), adj = 0)
    all.fdr.genes <- which(env$fdr.g.m[, m] < 0.01)
    plus.fdr.genes <- which(env$indata[all.fdr.genes, m] > 
                              0)
    minus.fdr.genes <- which(env$indata[all.fdr.genes, m] <= 
                               0)
    text(0.1, 0.5, paste("# genes with fdr < 0.01 =", length(all.fdr.genes), 
                         " (", length(plus.fdr.genes), "+ /", length(minus.fdr.genes), 
                         " -)"), adj = 0)
    text(0.1, 0.35, paste("<FC> =", round(mean(env$indata[, 
                                                          m]), 2)), adj = 0)
    text(0.1, 0.3, paste("<p-value> =", round(10^mean(log10(env$p.g.m[, 
                                                                      m])), 2)), adj = 0)
    text(0.1, 0.25, paste("<fdr> =", round(mean(env$fdr.g.m[, 
                                                            m]), 2)), adj = 0)
    par(mar = c(2, 3, 3, 1))
    image(matrix(env$metadata[, m], env$preferences$dim.1stLvlSom, 
                 env$preferences$dim.1stLvlSom), axes = FALSE, col = env$color.palette.portraits(1000), 
          main = "Portrait", cex.main = 1.5)
    axis(1, seq(0, 1, length.out = env$preferences$dim.1stLvlSom/10 + 
                  1), c(1, seq(10, env$preferences$dim.1stLvlSom, 
                               length.out = env$preferences$dim.1stLvlSom/10)), 
         cex.axis = 1)
    axis(2, seq(0, 1, length.out = env$preferences$dim.1stLvlSom/10 + 
                  1), c(1, seq(10, env$preferences$dim.1stLvlSom, 
                               length.out = env$preferences$dim.1stLvlSom/10)), 
         cex.axis = 1, las = 1)
    box()
    n.genes <- min(100, nrow(env$indata))
    n.map <- matrix(0, env$preferences$dim.1stLvlSom, env$preferences$dim.1stLvlSom)
    set.genes <- names(sort(env$p.g.m[, m]))[1:n.genes]
    gs.nodes <- env$som.result$feature.BMU[set.genes]
    n.map[as.numeric(names(table(gs.nodes)))] <- table(gs.nodes)
    n.map[which(n.map == 0)] <- NA
    n.map <- matrix(n.map, env$preferences$dim.1stLvlSom)
    lim <- c(1, env$preferences$dim.1stLvlSom) + env$preferences$dim.1stLvlSom * 
      0.01 * c(-1, 1)
    colr <- env$color.palette.heatmaps(1000)[(na.omit(as.vector(n.map)) - 
                                                min(n.map, na.rm = TRUE))/max(1, (max(n.map, na.rm = TRUE) - 
                                                                                    min(n.map, na.rm = TRUE))) * 999 + 1]
    par(mar = c(2, 3, 3, 1))
    plot(which(!is.na(n.map), arr.ind = TRUE), xlim = lim, 
         ylim = lim, pch = 16, axes = FALSE, xlab = "", ylab = "", 
         xaxs = "i", yaxs = "i", col = colr, main = "Top 100 DE genes", 
         cex.main = 1.5, cex = 0.5 + na.omit(as.vector(n.map))/max(n.map, 
                                                                   na.rm = TRUE) * 2.8)
    axis(1, c(1, seq(10, env$preferences$dim.1stLvlSom, 
                     length.out = env$preferences$dim.1stLvlSom/10)), 
         c(1, seq(10, env$preferences$dim.1stLvlSom, length.out = env$preferences$dim.1stLvlSom/10)), 
         cex.axis = 1)
    axis(2, c(1, seq(10, env$preferences$dim.1stLvlSom, 
                     length.out = env$preferences$dim.1stLvlSom/10)), 
         c(1, seq(10, env$preferences$dim.1stLvlSom, length.out = env$preferences$dim.1stLvlSom/10)), 
         cex.axis = 1, las = 1)
    box()
    par(mar = c(3, 3, 3, 1))
    plot(env$indata[set.genes, m], -log10(env$p.g.m[set.genes, 
                                                    m]), xlim = c(-1, 1) * max(abs(env$indata[set.genes, 
                                                                                              m])), xlab = "", ylab = "", pch = 16, col = "gray30", 
         las = 1)
    mtext("log FC", 1, line = 2, cex = 0.6)
    mtext("-log10(p)", 2, line = 2, cex = 0.6)
    abline(v = 0, lty = 3, col = "gray80", lwd = 1.5)
    abline(h = -log10(0.05), lty = 3, col = "gray80", lwd = 1.5)
    n.genes <- 20
    de.genes <- names(sort(env$p.g.m[, m]))
    de.genes <- de.genes[which(env$indata[de.genes, m] > 
                                 0)][1:n.genes]
    de.genes <- de.genes[!(is.na(de.genes))]
    de.genes.labels <- env$gene.info$names[de.genes]
    de.genes.labels[which(de.genes.labels == "")] <- de.genes[which(de.genes.labels == 
                                                                      "")]
    par(mar = c(0, 0, 0, 0))
    x.coords <- c(0, 0.06, 0.2, 0.28, 0.36, 0.44, 0.52)
    y.coords <- seq(0.75, 0.4, length.out = n.genes)
    plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", 
         xlim = c(0, 1), ylim = c(0, 1))
    text(0, 0.88, "Differentially expressed genes", cex = 1.8, 
         adj = 0)
    text(x.coords, rep(c(0.82, 0.8), 4)[1:7], c("Rank", 
                                                "ID", "log(FC)", "p-value", "fdr", "Metagene", "Description"), 
         cex = 1, adj = 0)
    text(x.coords[1], 0.77, "Overexpressed", cex = 0.8, 
         adj = 0, font = 3)
    text(x.coords[1], y.coords, c(1:n.genes), adj = 0)
    text(x.coords[2], y.coords, de.genes.labels, cex = 0.6, 
         adj = 0)
    rect(x.coords[3] - 0.02, y.coords[1] + 0.01, 1, 0, border = "white", 
         col = "white")
    text(x.coords[3], y.coords, round(env$indata[de.genes, 
                                                 m], 2), cex = 0.6, adj = 0)
    text(x.coords[4], y.coords, format(env$p.g.m[de.genes, 
                                                 m], digits = 1), cex = 0.6, adj = 0)
    text(x.coords[5], y.coords, format(env$fdr.g.m[de.genes, 
                                                   m], digits = 1), cex = 0.6, adj = 0)
    text(x.coords[6], y.coords, env$gene.info$coordinates[de.genes], 
         cex = 0.6, adj = 0)
    text(x.coords[7], y.coords, env$gene.info$descriptions[de.genes], 
         cex = 0.6, adj = 0)
    de.genes <- names(sort(env$p.g.m[, m]))
    de.genes <- de.genes[which(env$indata[de.genes, m] < 
                                 0)][1:n.genes]
    de.genes <- de.genes[!(is.na(de.genes))]
    de.genes.labels <- env$gene.info$names[de.genes]
    de.genes.labels[which(de.genes.labels == "")] <- de.genes[which(de.genes.labels == 
                                                                      "")]
    y.coords <- seq(0.35, 0.02, length.out = n.genes)
    text(x.coords[1], 0.37, "Underexpressed", cex = 0.8, 
         adj = 0, font = 3)
    text(x.coords[1], y.coords, c(1:n.genes), adj = 0)
    text(x.coords[2], y.coords, de.genes.labels, cex = 0.6, 
         adj = 0)
    rect(x.coords[3] - 0.02, y.coords[1] + 0.01, 1, 0, border = "white", 
         col = "white")
    text(x.coords[3], y.coords, round(env$indata[de.genes, 
                                                 m], 2), cex = 0.6, adj = 0)
    text(x.coords[4], y.coords, format(env$p.g.m[de.genes, 
                                                 m], digits = 1), cex = 0.6, adj = 0)
    text(x.coords[5], y.coords, format(env$fdr.g.m[de.genes, 
                                                   m], digits = 1), cex = 0.6, adj = 0)
    text(x.coords[6], y.coords, env$gene.info$coordinates[de.genes], 
         cex = 0.6, adj = 0)
    text(x.coords[7], y.coords, env$gene.info$descriptions[de.genes], 
         cex = 0.6, adj = 0)
    par(mar = c(3, 6, 2, 6))
    hist(env$p.g.m[, m], bre = 20, freq = FALSE, xlab = "p-value", 
         ylab = "", main = "p-values", ylim = c(0, ylim.max), 
         las = 1, cex.main = 1.5, cex.lab = 1, cex.axis = 1)
    box()
    mtext("Density", side = 2, line = 3, cex = 1)
    mtext("FDR", side = 4, line = 3, cex = 1)
    mtext(paste("%DE =", round(env$perc.DE.m[colnames(env$indata)[m]], 
                               2)), line = -1.2, cex = 0.5)
    abline(h = env$n.0.m[colnames(env$indata)[m]], col = "gray", 
           lwd = 2)
    par(new = TRUE)
    plot(0, type = "n", xlim = c(0, 1), ylim = c(0, 1), 
         xlab = "", ylab = "", axes = FALSE)
    axis(4, seq(0, 1, 0.2), seq(0, 1, 0.2), las = 1, cex.axis = 1)
    o <- order(env$p.g.m[, m])
    lines(env$p.g.m[o, m], env$fdr.g.m[o, m], lty = 3, lwd = 3)
    legend("topright", c("p", expression(env$eta[0]), "fdr"), 
           col = c("black", "gray", "black"), lty = c(1, 1, 
                                                      3), lwd = c(1, 1, 2), cex = 0.7)
    if (env$preferences$activated.modules$geneset.analysis) {
      n.sets <- 20
      top.gs.score <- sort(env$spot.list.samples[[m]]$GSZ.score, 
                           decreasing = TRUE)[1:n.sets]
      top.gs.score <- top.gs.score[!(is.na(top.gs.score))]
      top.gs.p <- env$spot.list.samples[[m]]$GSZ.p.value[names(top.gs.score)]
      par(mar = c(0, 0, 0, 0))
      x.coords <- c(0, 0.1, 0.18, 0.3, 0.39, 0.47)
      y.coords <- seq(0.75, 0.4, length.out = n.sets)
      plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", 
           xlim = c(0, 1), ylim = c(0, 1))
      text(0, 0.88, "Differentially expressed gene sets", 
           cex = 1.8, adj = 0)
      text(x.coords, 0.82, c("Rank", "GSZ", "p-value", 
                             "#all", "Geneset", ""), cex = 1, adj = 0)
      text(x.coords[1], 0.77, "Overexpressed", cex = 0.8, 
           adj = 0, font = 3)
      text(x.coords[1], y.coords, c(1:n.genes), adj = 0)
      text(x.coords[2], y.coords, round(top.gs.score, 
                                        2), cex = 0.6, adj = 0)
      text(x.coords[3], y.coords, format(top.gs.p, digits = 1), 
           cex = 0.6, adj = 0)
      text(x.coords[4], y.coords, sapply(env$gs.def.list[names(top.gs.score)], 
                                         function(x) {
                                           length(x$Genes)
                                         }), cex = 0.6, adj = 0)
      text(x.coords[5], y.coords, sapply(env$gs.def.list, 
                                         function(x) {
                                           x$Type
                                         })[names(top.gs.score)], cex = 0.6, adj = 0)
      text(x.coords[6], y.coords, names(top.gs.score), 
           cex = 0.6, adj = 0)
      top.gs.score <- sort(env$spot.list.samples[[m]]$GSZ.score, 
                           decreasing = FALSE)[1:n.sets]
      top.gs.score <- top.gs.score[!(is.na(top.gs.score))]
      top.gs.p <- env$spot.list.samples[[m]]$GSZ.p.value[names(top.gs.score)]
      y.coords <- seq(0.35, 0.02, length.out = n.sets)
      text(x.coords[1], 0.37, "Underexpressed", cex = 0.8, 
           adj = 0, font = 3)
      text(x.coords[1], y.coords, c(1:n.genes), adj = 0)
      text(x.coords[2], y.coords, round(top.gs.score, 
                                        2), cex = 0.6, adj = 0)
      text(x.coords[3], y.coords, format(top.gs.p, digits = 1), 
           cex = 0.6, adj = 0)
      text(x.coords[4], y.coords, sapply(env$gs.def.list[names(top.gs.score)], 
                                         function(x) {
                                           length(x$Genes)
                                         }), cex = 0.6, adj = 0)
      text(x.coords[5], y.coords, sapply(env$gs.def.list, 
                                         function(x) {
                                           x$Type
                                         })[names(top.gs.score)], cex = 0.6, adj = 0)
      text(x.coords[6], y.coords, names(top.gs.score), 
           cex = 0.6, adj = 0)
      p <- env$spot.list.samples[[m]]$GSZ.p.value
      fdrtool.result <- suppressWarnings(fdrtool(p, statistic = "pvalue", 
                                                 verbose = FALSE, plot = FALSE))
      fdr.spot.list.samples <- fdrtool.result$lfdr
      Fdr.spot.list.samples <- fdrtool.result$qval
      n.0.spot.list.samples <- fdrtool.result$param[1, 
                                                    "eta0"]
      perc.DE.spot.list.samples <- 1 - n.0.spot.list.samples
      par(mar = c(3, 6, 2, 6))
      hist(p, bre = 20, freq = FALSE, xlab = "p-value", 
           ylab = "", main = "p-values", las = 1, cex.main = 1.5, 
           cex.lab = 1, cex.axis = 1)
      box()
      mtext("Density", side = 2, line = 3, cex = 1)
      mtext("FDR", side = 4, line = 3, cex = 1)
      mtext(paste("%DE =", round(perc.DE.spot.list.samples, 
                                 2)), line = -1.2, cex = 0.5)
      abline(h = n.0.spot.list.samples, col = "gray", 
             lwd = 2)
      par(new = TRUE)
      plot(0, type = "n", xlim = c(0, 1), ylim = c(0, 
                                                   1), xlab = "", ylab = "", axes = FALSE)
      axis(4, seq(0, 1, 0.2), seq(0, 1, 0.2), las = 1, 
           cex.axis = 1)
      o = order(p)
      lines(p[o], Fdr.spot.list.samples[o], lty = 2, lwd = 2)
      lines(p[o], fdr.spot.list.samples[o], lty = 3, lwd = 3)
      legend("topright", c("p", expression(env$eta[0]), 
                           "Fdr", "fdr"), col = c("black", "gray", "black", 
                                                  "black"), lty = c(1, 1, 2, 3), lwd = c(1, 1, 
                                                                                         1, 2), cex = 0.7)
    }
    dev.off()
  }
}

opossom.run.2 <- function (env) 
{
  util.info("Started:", env$preferences$started)
  util.info("Name:", env$preferences$dataset.name)
  env <- pipeline.checkInputParameters(env)
  if (!env$passedInputChecking) {
    return()
  }
  if (env$preferences$activated.modules$primary.analysis) {
    env$preferences$system.info <- Sys.info()
    env$preferences$session.info <- sessionInfo()
    env$preferences$started <- format(Sys.time(), "%a %d %b %Y %X")
  }
  if (env$preferences$activated.modules$reporting) {
    dir.create(paste(env$files.name, "- Results"), showWarnings = FALSE)
    dir.create(paste(env$files.name, "- Results/CSV Sheets"), 
               showWarnings = FALSE)
    if (env$preferences$activated.modules$primary.analysis) {
      pipeline.qualityCheck(env)
    }
  }
  if (env$preferences$activated.modules$primary.analysis || 
      env$preferences$activated.modules$geneset.analysis) {
    util.info("Loading gene annotation data.")
    env <- pipeline.prepareAnnotation(env)
  }
  if (env$preferences$activated.modules$primary.analysis) {
    util.info("Processing SOM. This may take several time until next notification.")
    env <- pipeline.prepareIndata(env)
    env <- pipeline.generateSOM(env)
    filename <- paste(env$files.name, "pre.RData")
    util.info("Saving environment image:", filename)
    save(env, file = filename)
    util.info("Processing Differential Expression Statistics")
    env <- pipeline.diffExpressionStatistics(env)
    util.info("Detecting Spots")
    env <- pipeline.detectSpotsSamples(env)
    env <- pipeline.detectSpotsModules(env)
    env <- pipeline.patAssignment(env)
    env <- pipeline.groupAssignment(env)
  }
  if (env$preferences$activated.modules$geneset.analysis) {
    util.info("Calculating Geneset Enrichment")
    env <- pipeline.genesetStatisticSamples(env)
    env <- pipeline.genesetStatisticModules(env)
  }
  if (env$preferences$activated.modules$psf.analysis) {
    util.info("Calculating Pathway Signal Flow (PSF)")
    env <- pipeline.PSFcalculation(env)
  }
  if (env$preferences$activated.modules$primary.analysis || 
      env$preferences$activated.modules$geneset.analysis) {
    filename <- paste(env$files.name, ".RData", sep = "")
    util.info("Saving environment image:", filename)
    save(env, file = filename)
    if (file.exists(paste(env$files.name, "pre.RData")) && 
        file.exists(filename)) {
      file.remove(paste(env$files.name, "pre.RData"))
    }
  }
  if (env$preferences$activated.modules$reporting) {
    util.info("Plotting Supporting Information")
    pipeline.supportingMaps(env)
    pipeline.entropyProfiles(env)
    pipeline.topologyProfiles(env)
    if (length(env$chromosome.list) > 0) {
      util.info("Plotting Chromosome Expression Reports")
      pipeline.chromosomeExpressionReports(env)
    }
    if (ncol(env$indata) < 1000) {
      util.info("Plotting Sample Portraits")
      pipeline.sampleExpressionPortraits(env)
    }
    if (env$preferences$activated.modules$sample.similarity.analysis && 
        ncol(env$indata) > 2) {
      util.info("Plotting Sample Similarity Analysis")
      dir.create(file.path(paste(env$files.name, "- Results"), 
                           "Sample Similarity Analysis"), showWarnings = FALSE)
      pipeline.sampleSimilarityAnalysisED(env)
      pipeline.sampleSimilarityAnalysisCor(env)
      pipeline.sampleSimilarityAnalysisICA(env)
      pipeline.sampleSimilarityAnalysisSOM(env)
    }
    if (env$preferences$activated.modules$geneset.analysis) {
      dir.create(paste(env$files.name, "- Results/Geneset Analysis"), 
                 showWarnings = FALSE)
      util.info("Plotting Geneset Enrichment Heatmaps")
      pipeline.genesetOverviews(env)
      util.info("Plotting Geneset Profiles and Maps")
      pipeline.genesetProfilesAndMaps(env)
      util.info("Calculating Cancer Hallmark Enrichment")
      pipeline.cancerHallmarks(env)
    }
    if (env$preferences$activated.modules$psf.analysis) {
      util.info("Plotting PSF results")
      pipeline.PSFoutput(env)
    }
    util.info("Writing Gene Lists")
    pipeline.geneLists(env)
    util.info("Plotting Summary Sheets (Samples)")
    pipeline.summarySheetsSamples.2(env)
    util.info("Plotting Summary Sheets (Modules & PATs)")
    pipeline.summarySheetsModules(env)
    pipeline.summarySheetsPATs(env)
    if (env$preferences$activated.modules$group.analysis && 
        length(unique(env$group.labels)) >= 2) {
      util.info("Processing Group-centered Analyses")
      pipeline.groupAnalysis(env)
    }
    if (env$preferences$activated.modules$difference.analysis) {
      util.info("Processing Difference Analyses")
      pipeline.differenceAnalyses(env)
    }
    util.info("Generating HTML Report")
    pipeline.htmlSampleSummary(env)
    pipeline.htmlModuleSummary(env)
    pipeline.htmlGenesetAnalysis(env)
    pipeline.htmlPsfAnalysis(env)
    pipeline.htmlSummary(env)
  }
  util.info("Finished:", format(Sys.time(), "%a %b %d %X"))
}



# Parsing command line arguments
p <- arg_parser("Functional Enrichment Analysis")

p <- add_argument(p, "--type", help="Cancer type", default = "colorectal",
                  "mut_proc_colorectal_coadread_tcga_pan_can_atlas_2018_latest.csv")
p <- add_argument(p, "--altfile", help="Path to input data file of genomic alterations", default = file.path("/media/disk_2/data_science", 
                  "clustering_test/data/datasets",
                  "mut_proc_colorectal_coadread_tcga_pan_can_atlas_2018_latest.csv"))
p <- add_argument(p, "--clustfile", help="Path to input data file of cluster assignments", default = file.path("/media/disk_2/data_science", 
                  "clustering_test/data/outputs/TCGA/colorectal",
                  "colorectal_coadread_tcga_pan_can_atlas_2018_cluster_assignment.txt"))
p <- add_argument(p, "--outfolder", help="Name of the output folder, saved in the same directory",
                  default = paste0(format(Sys.Date(), format="%b%d"), "_SOM_Analysis_"))
p <- add_argument(p, "--reporting", help = "Generate a report", default = TRUE)
p <- add_argument(p, "--primary.analysis", help = "Perform primary analysis", default = TRUE)                                  
p <- add_argument(p, "--sample.similarity.analysis", help = "Perform sample similarity analysis", default = FALSE)                                              
p <- add_argument(p, "--geneset.analysis", help = "Perform geneset analysis", default = TRUE)                               
p <- add_argument(p, "--geneset.analysis.exact", help = "Perform exact geneset analysis", default = FALSE)             
p <- add_argument(p, "--group.analysis", help = "Perform group analysis", default = FALSE)                                         
p <- add_argument(p, "--difference.analysis", help = "Perform difference analysis", default = FALSE)                                            

                  
args <- parse_args(p)


# Preprocessing data
altData <- read.csv(args$altfile)
altData <- altData[rowSums(altData[-1]) > 1,]
clustData <- read.csv(args$clustfile, sep = "\t")
myData <- merge(clustData, altData )
indata <- t(as.matrix(myData[, -(1:2)]))
metadata <- myData[, (1:2)]
colnames(indata) <- rownames(metadata) <- metadata$Gene_Expession_Barcode
ordering <- order(metadata$cluster)
metadata$Gene_Expession_Barcode <- NULL
my_genes <- rownames(indata)
merged_indata <- t(do.call(cbind,(by(indata, my_genes, colSums))))

env <- opossom.new(list(dataset.name = paste0(args$outfolder, args$type),
                        
                        dim.1stLvlSom = "auto",
                        dim.2ndLvlSom = "auto",
                        
                        database.biomart = "ENSEMBL_MART_ENSEMBL",
                        database.dataset = "hsapiens_gene_ensembl",
                        database.id.type = "hgnc_symbol",
                        
                        activated.modules = list( "reporting" = args$reporting,
                                                  "primary.analysis" = args$primary.analysis, 
                                                  "sample.similarity.analysis" = args$sample.similarity.analysis,
                                                  "geneset.analysis" = args$geneset.analysis, 
                                                  "geneset.analysis.exact" = args$geneset.analysis.exact,
                                                  "group.analysis" = args$group.analysis,
                                                  "difference.analysis" = args$difference.analysis),
                        
                        
                        feature.centralization = T, ##switch off if needed
                        sample.quantile.normalization = F,
                        
                        pairwise.comparison.list = list() ) )


env$indata <- merged_indata[, ordering]
env$group.labels <- paste0("C" ,metadata$cluster[ordering])
env$group.colors <- glasbey(length(unique(env$group.labels)))[ match( env$group.labels,unique(env$group.labels))
]

# execute
# opossom.run()
opossom.run.2(env)

## extracting geneset propagated data for BP

selected_gs <- c()
for (i in rownames(tb)) {
 if (env$gs.def.list[[i]]$Type == 'BP') {
   selected_gs <- c(i, selected_gs)
 }
}

tb_selected <- tb[selected_gs, ]
colnames(tb_selected) <- str_replace(colnames(tb_selected), "V", "")
tb_selected_ordered <- tb_selected[ , order(as.numeric(colnames(tb_selected)))]
colnames(tb_selected_ordered) <- altData$patient_id
tb_final <- rownames_to_column(as.data.frame(t(tb_selected_ordered)), "patient_id")
write.csv(tb_final, "/media/disk_2/data_science/clustering_test/data/datasets/propagated_mut_colorectal_coadread_2018.csv")
write.csv()