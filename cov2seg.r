
main <- function(fn, output_dir) {
  # Load necessary libraries
  library(facets)
  library(dplyr)
  library(outliers)
  library(DNAcopy)
  library(GenomicRanges)
  library(future.apply)
  library(Rsamtools)
  

  options(digits = 4)
  options(scipen = 999)
  set.seed(123)

  fn_base_name <- basename(fn)

  site_density <- 12
  min_makers <- 40
  
  # Process SNP data
  cnv_data <- process_snp_data(fn)
  probe_data <- na.omit(cnv_data$probe)
  snp_data <- na.omit(cnv_data$snp)
  
  # Save processed data to files
  write.table(probe_data, file = paste(output_dir, paste(fn_base_name, "probes", sep="."), sep="/"), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  write.table(snp_data, file = paste(output_dir, paste(fn_base_name, "snps", sep="."), sep="/"), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  

  # CNV and LOH calling
  segments <- segment_data_parallel(probe_data, fn, "DNAcopy")
  
  snp_data$imb <- abs(0.5 - snp_data$baf)
  lohs <- segment_data_parallel(snp_data, fn, "LOH")
  
  # Save raw segment and LOH data
  write.table(segments, file=paste(output_dir, paste(fn_base_name, ".sd1.sdundo.raw.seg", sep=""), sep="/"), quote=FALSE, sep="\t", row.names=FALSE)
  write.table(lohs, file=paste(output_dir, paste(fn_base_name, ".sd1.sdundo.raw.lohs", sep=""), sep="/"), quote=FALSE, sep="\t", row.names=FALSE)
  
  # Filter segments and LOH data
  segments_filtered <- segments %>%
  filter(Markers >= min_makers, Size > 100000, (Size / Markers / 1000) < site_density) %>%
  bind_rows(segments %>% filter(Chromosome == "Y", Markers > 100))

  lohs_processed <- process_lohs(lohs, min_makers)
  write.table(lohs_processed, file=paste(output_dir, paste(fn_base_name, ".sd1.sdundo.loh", sep=""), sep="/"), quote=FALSE, sep="\t", row.names=FALSE)

  # Adjust 'Start' values in 'segments_filtered' to avoid overlaps
  overlap_indices <- match(segments_filtered$Start, segments_filtered$End)
  if (length(na.omit(overlap_indices)) > 0) {
    segments_filtered[na.omit(overlap_indices) + 1, "Start"] <- segments_filtered[na.omit(overlap_indices) + 1, "Start"] + 2
  }


  # Combine 'segments_filtered' and 'lohs' if 'lohs' is not empty
  if (nrow(lohs_processed) > 0) {
    segments_loh <- combine_cn_log(segments_filtered, lohs_processed)
    res <- data.frame(
      seqnames = seqnames(segments_loh),
      starts = start(segments_loh),
      ends = end(segments_loh),
      Markers = mcols(segments_loh)$extraInfo.Markers,
      Value = mcols(segments_loh)$extraInfo.Value,
      Size = mcols(segments_loh)$extraInfo.Size
    )  
  } else {
    res <- segments_filtered
  }

  write.table(res, file=paste(output_dir, paste(fn_base_name, ".sd1.sdundo.seg", sep=""), sep="/"), quote=FALSE, sep="\t", row.names=FALSE)
}

process_snp_data <- function(fn) {
  high_cut <- 2  
  
  normalize <- function(z) {
    return ((z - median(z)) / (max(z) - min(z)))
  }

  x <- readSnpMatrix(filename = fn, perl.pileup = TRUE)
  
  cxy <- x[x$Chromosome %in% c('X', 'Y'), ]
  cn <- x[!x$Chromosome %in% c('X', 'Y'), ]

  # Process autosomes
  min_ND <- min(10, as.integer(quantile(cn$NOR.DP, 0.3) * 0.9 + 0.5))
  xx <- preProcSample(cn, gbuild = "hg38", ndepth = min_ND, het.thresh = 0.22, snp.nbhd = 50)
  xx2 <- subset(xx$pmat, keep == 1)

  # Process sex chromosomes
  xmin_ND <- as.integer(quantile(cxy$NOR.DP, 0.3) * 0.9 + 0.5)
  xx3 <- preProcSample(cxy, gbuild = "hg38", ndepth = xmin_ND, het.thresh = 0.22, snp.nbhd = 50)
  xx4 <- subset(xx3$pmat, keep == 1)
  xx4$chrom[xx4$chrom == 23] <- "X"

  # Process Y chromosomes specifically
  if ("Y" %in% cxy$Chromosome) {
    datay <- cxy %>%
      filter(Chromosome == "Y") %>%
      mutate(
        vafN = (NOR.DP - NOR.RD) / NOR.DP,
        vafT = (TUM.DP - TUM.RD) / TUM.DP,
        het = 0,
        keep = 1
      ) %>%
      select(chrom = Chromosome, maploc = Position, rCountN = NOR.DP, rCountT = TUM.DP, vafN, vafT, het, keep)
    xx4 <- bind_rows(xx4, datay)
  }

  # Remove "Y" if fewer than 100 after addition
  if ("Y" %in% xx4$chrom) {
    count_y <- sum(xx4$chrom == "Y")
    if (count_y < 100) {
      xx4 <- filter(xx4, chrom != "Y")
    }
  }

  combined_xx <- rbind(xx2, xx4) %>% mutate(chrom = as.character(chrom))
  
  # Join with original data and calculate metrics
  ns <- sum(xx2$rCountN) / 1e6
  ts <- sum(xx2$rCountT) / 1e6
  combined_xx_joined <- inner_join(x, combined_xx, by = c("Chromosome" = "chrom", "Position" = "maploc"))
  combined_xx_joined$lrr <- atan(log2((combined_xx_joined$TUM.DP / ts) / (combined_xx_joined$NOR.DP / ns))) * 2 / pi
  combined_xx_joined$baf <- 1 - combined_xx_joined$TUM.RD / combined_xx_joined$TUM.DP

  # Final filtering and selection based on 'keep' flag
  tns <- sum(x$NOR.DP) / 1e6
  x$n_ratio <- x$NOR.DP / tns
  x$nsd <- scores(x$n_ratio, type = "t")
  
  x_filtered <- x[abs(x$nsd) < high_cut, -grep("n_ratio|nsd", names(x))]
  x2 <- inner_join(combined_xx_joined, x_filtered, by = c("Chromosome", "Position"))
  x2 <- x2[x2$keep == 1, c("Chromosome", "Position", "lrr", "baf")]
  
  # Extract specific columns for 'probe' and 'snp'
  probe_columns <- c("Chromosome", "Position", "lrr")
  snp_columns <- c("Chromosome", "Position", "baf")
  at_probe <- x2[probe_columns]
  at_snp <- x2[snp_columns]
  
  return(list(probe = at_probe, snp = at_snp))
}

segment_data <- function(data, sampleID, mode="DNAcopy") {
    segs <- NULL
    chroms <- unique(data$Chromosome)
    for (c in 1:length(chroms)) {
        print(c)
        tdata <- data[data$Chromosome == chroms[c],]
        if (nrow(tdata) > 0) {
            if (mode == "DNAcopy") {
                cnaObject <- segment(smooth.CNA(CNA(genomdat=tdata$lrr, chrom = tdata$Chromosome, maploc = tdata$Position, data.type='logratio', sampleid=sampleID)), undo.splits='sdundo', undo.SD=1, verbose=0)
            } else if (mode == "LOH") {
                cnaObject <- segment(CNA(genomdat = tdata$imb, chrom = tdata$Chromosome, maploc = tdata$Position, data.type="binary", sampleid=sampleID), undo.splits='sdundo', undo.SD=1.5, verbose=0)
            } else {
                stop("Invalid mode specified. Use 'DNAcopy' or 'LOH'.")
            }
            segs <- rbind(segs, cnaObject$output)
        }
    }
    segs <- segs[,2:6]
    colnames(segs) <- c('Chromosome', 'Start', 'End', 'Markers', 'Value')
    segs$Size <- segs$End - segs$Start
    return(segs)
}


segment_data_parallel <- function(data, sampleID, mode="DNAcopy") {  
  plan(multisession,  workers = 8)  # Adjust based on your available resources
  
  chroms <- unique(data$Chromosome)
  process_chromosome <- function(chrom) {
    tdata <- data[data$Chromosome == chrom, ]
    if (nrow(tdata) > 0) {
      if (mode == "DNAcopy") {
        cnaObject <- segment(smooth.CNA(CNA(genomdat=tdata$lrr, chrom = tdata$Chromosome, maploc = tdata$Position, data.type='logratio', sampleid=sampleID)), undo.splits='sdundo', undo.SD=1, verbose=0)
      } else if (mode == "LOH") {
        cnaObject <- segment(CNA(genomdat = tdata$imb, chrom = tdata$Chromosome, maploc = tdata$Position, data.type="binary", sampleid=sampleID), undo.splits='sdundo', undo.SD=1.5, verbose=0)
      } else {
        stop("Invalid mode specified. Use 'DNAcopy' or 'LOH'.")
      }
      return(cnaObject$output)
    } else {
      return(NULL)
    }
  }
  
  # Apply the function to each chromosome in parallel
  segs_list <- future_lapply(chroms, process_chromosome)
  segs <- do.call(rbind, segs_list)
  
  # Finalize the result
  segs <- segs[,2:6]
  colnames(segs) <- c('Chromosome', 'Start', 'End', 'Markers', 'Value')
  segs$Size <- segs$End - segs$Start
  
  return(segs)
}


combine_cn_log <- function(seg, LOH) {
        gseg<-makeGRangesFromDataFrame(seg, keep.extra.columns=TRUE, ignore.strand=TRUE)
        gloh<-makeGRangesFromDataFrame(LOH, keep.extra.columns=TRUE, ignore.strand=TRUE)

        g_non_loh<-setdiff(gseg, gloh, ignore.strand=TRUE)
        overlapping_indexes <- findOverlaps(g_non_loh, gseg)
        overlapping_granges <- gseg[subjectHits(overlapping_indexes)]
        extra_info <- mcols(overlapping_granges)
        non_LOH_seg <- GRanges(
                seqnames = seqnames(g_non_loh),
                ranges = ranges(g_non_loh),
                extraInfo = extra_info
        )

        g_loh<-intersect(gseg, gloh, ignore.strand=TRUE)
        overlapping_indexes_loh <- findOverlaps(g_loh, gseg)
        overlapping_granges_loh <- gseg[subjectHits(overlapping_indexes_loh)]
        extra_info <- mcols(overlapping_granges_loh)
        LOH_res <- GRanges(
                seqnames = seqnames(g_loh),
                ranges = ranges(g_loh),
                extraInfo = extra_info
                )

        wgs_cn_LOH <- sort(c(non_LOH_seg, LOH_res))
        return(wgs_cn_LOH)
}


process_lohs <- function(lohs, min_makers) {
  # Filter 'lohs' for small and large categories based on conditions
  lohs_small <- lohs[lohs$Value > 0.485 & lohs$Markers >= min_makers & lohs$Size > 1000000, ]
  lohs_large <- lohs[lohs$Value > 0.47 & lohs$Markers >= min_makers & lohs$Size > 4000000, ]

  # Combine 'lohs_small' and 'lohs_large', remove duplicates, and order
  lohs_combined <- rbind(lohs_small, lohs_large)
  lohs_combined <- unique(lohs_combined)
  lohs_combined <- lohs_combined[order(lohs_combined$Chromosome, lohs_combined$Start), ]

  return(lohs_combined)
}


args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) {
  fn <- args[1]
  output_dir <- args[2]
  main(fn, output_dir)
} else {
  stop("Insufficient arguments provided. Expected filename and output directory.")
}

