#!/usr/bin/R

#######################################
# Quality Control of single-cell data #
#######################################

# This script contains the QC steps for single-cell analysis and gets data
# ready for downstream protocols; clustering, scDGEA, etc.

library(optparse)

optlist <- list(
  make_option(
    opt_str = c("-y", "--yaml"), type = "character",
    help = "Configuration file: Instructions in YAML format."
  ),
  make_option(
    opt_str = c("-l", "--log"), default = TRUE,
    help = "Log file: Create a log file rather to stdout."
  ),
  make_option(
    opt_str = c("-v", "--verbose"), default = TRUE,
    help = "Verbose: Show progress."
  )
)

# Getting arguments from command line and setting their values to their respective variable names.
optparse <- OptionParser(option_list = optlist)
opt <- parse_args(optparse)
if(interactive()){ # Example/manually
  opt$yaml = "/home/ciro/amica/scripts/SiEs12_qc.yaml"
}
config = yaml::read_yaml(opt$yaml)
options(warn = -1)

if(suppressMessages(require(crayon))){
  cyan = crayon::cyan; red_bold = crayon::red$bold
}else{ cyan = red_bold = c }

if(opt$verbose){
  cat(cyan("\n************ Vijay Lab - LJI\n"))
  cat(cyan("-------------------------------------\n"))
  cat(red_bold("------------ Quality Control analysis\n"))
}

#### Directories structure ####-------------------------------------------------
output_dir <- paste0(sub("\\/$", "", config$output_dir), "/", config$project_name, "/")
if(!grepl("scratch|beegfs", getwd())){
  cat("No scratch folder involved; careful about temp files...\n")
  dir.create(output_dir, recursive = TRUE); setwd(output_dir)
}
cat("Working in:", getwd(), "\n")
#### ########### ######### ####-------------------------------------------------
### Log file ####---------------------------------------------------------------
log_file <- paste0(output_dir, '_qc.log')
register_log <- !interactive() && !grepl("scratch|beegfs", getwd()) && opt$log
# register_log <- (!interactive() && grepl("scratch|beegfs", getwd())) && opt$log
if(register_log){
  if(opt$verbose) cat("Check log file:", log_file, "\n")
  if(!file.exists(log_file)) file.create(log_file)
  out.file <- file(log_file, open = 'wt')
  sink(out.file) ; sink(out.file, type = 'message')
}
cat('Date and time:\n') ; st.time <- timestamp();
### ### #### ####---------------------------------------------------------------

#### Functions ####-------------------------------------------------------------
cat("Loading functions\n")
resources = c(
  "/home/ciro/scripts/handy_functions/devel/utilities.R",
  # dircheck show_parameters newlines
  "/home/ciro/scripts/handy_functions/devel/file_reading.R",
  # readfile file.find
  "/home/ciro/scripts/clustering/R/utilities.R", # get_source_data add_pcts filters_qc_limits
  "/home/ciro/scripts/handy_functions/devel/filters.R",
  # filters_complex
  "/home/ciro/scripts/clustering/R/plotting.R", # qc_scatters qc_violin simple_violin
  "/home/ciro/scripts/handy_functions/devel/plots.R"
  # make_grid
)
for(i in resources){ source(i) }
#### ######### ####-------------------------------------------------------------

#### Loading dependencies ####--------------------------------------------------
packages <- c("ggplot2", "cowplot")
loaded <- lapply(X = packages, FUN = function(x){
  if(opt$verbose) cat("*", x, "\n")
  suppressMessages(require(package = x, quietly = TRUE, character.only = TRUE))
})
theme_set(theme_cowplot())
#### ####### ############ ####--------------------------------------------------

#### Presenting Parameters ####-------------------------------------------------
str(opt); str(config)
if(opt$verbose) cat(cyan("----------- Digesting parameters -----------------\n"))
qc_filters <- filters_qc_limits(config$filtering)
if(is.null(config$ident_vars)) config$ident_vars = "orig~"
parameters2set = c("filtering$nSamples_expressed" = "0.001")
for(i in 1:length(parameters2set)){
  param = paste0("config$", names(parameters2set)[i])
  command = paste0(param, " = ", parameters2set[[i]])
  if(is.null(eval(parse(text = param)))){
    if(opt$verbose) cat(" ", command, "\n"); eval(parse(text = command))
  }
}
#### ########## ########## ####-------------------------------------------------

#### Getting data ####----------------------------------------------------------
if(opt$verbose) cat(cyan('----------- Loading data -------------------------\n'))
Sys.time()
# Find 10X directory (barcodes, gene_names, counts), Seurat object, CSV file, TXT, RDS
expr_data <- get_source_data(
  xpath = config$input_expression,
  pj = config$project_name,
  metadata = config$metadata[1],
  merge_counts = grepl("meco|merge", config$project_name) || !dir.exists(config$input_expression),
  verbose = opt$verbose
)
Sys.time()
config$input_expression <- expr_data$source
addmetadataf <- if(length(config$metadata[-1]) > 0) config$metadata[-1] else "no_file"
meta_data <- remove.factors(expr_data$mdata)
expr_data <- expr_data$edata

addmetadataf <- unlist(strsplit(path.expand(addmetadataf), "~"))
tvar <- file.exists(addmetadataf)
if(any(tvar)){
  addmetadataf <- addmetadataf[tvar]
  if(opt$verbose) cat("Adding metadata from:", addmetadataf, sep = "\n")
  for(addmetadataf_i in addmetadataf){
    addannot <- remove.factors(readfile(addmetadataf_i))
    print(dim(meta_data))
    # partial matching problem addressed in /home/ciro/covid19/scripts/partial_matching.R
    meta_data <- joindf(x = meta_data, y = addannot)
  }
  print(dim(meta_data))
}; str(meta_data)
meta_data = ident_tag(
  mdata = meta_data,
  tag = config$ident_vars[1],
  pattern = config$ident_vars[2],
  exclude = "^RNA",
  verbose = opt$verbose
)
config$ident_vars[1] = gsub("~", ".", config$ident_vars[1])

fname <- paste0('metadata.rdata'); save(meta_data, file = fname)
#### ####### #### ####----------------------------------------------------------

if(opt$verbose) cat('----------------------- Initial filters -------------------------\n')
str(meta_data)
filters_init = if(!is.null(config$filtering$subset)){
  lapply(names(config$filtering$subset), function(x) c(x, config$filtering$subset[[x]]) )
}
filtereddata <- filters_complex( # if you won't use all the data, you need to filter it first
  mdata = meta_data,
  filters = filters_init,
  verbose = opt$verbose
)
meta_data <- filtereddata[[1]]; rm(filtereddata)
# # this pre-filter would change the metrics
# # for every different subsampling of the metadata
# expr_data = expr_data[, rownames(meta_data)]

### Minimal filter: cells and genes ###
if(opt$verbose){
  cat("\nCalculating number of expressed features and expressing samples.\n")
  cat("Samples/cells:", ncol(expr_data), "\n"); cat("Features:", nrow(expr_data), "\n")
}
# Slow and will probably be re-calculated
tvar <- which(names(qc_filters[[1]]) == 'nFeature_RNA')
nFeature_min = ifelse(isTRUE(qc_filters$low[tvar] > 0), qc_filters$low[tvar], 0)
nFeature_min = min(c(nFeature_min, 200))
if(opt$verbose)
  cat('Minimal filtering: keeping cells with >=', nFeature_min, 'expressed genes.\n')
nFeatures_expressd <- Matrix::colSums(expr_data > 0); summary(nFeatures_expressd)
# filter cells first; and so features reflect the final set
expr_data <- expr_data[, nFeatures_expressd >= nFeature_min]

nSamples_expressed <- Matrix::rowSums(expr_data > 0); summary(nSamples_expressed)
nSamples_min = ceiling(config$filtering$nSamples_expressed[[1]] * ncol(expr_data))
nSamples_min = min(c(nSamples_min, 30))
tvar <- 'Minimal filtering: keeping genes with >='
if(opt$verbose) cat(tvar, nSamples_min, 'samples/cells expressing them.\n')
expr_data <- expr_data[nSamples_expressed >= nSamples_min, ]
if(opt$verbose){
  cat("Samples/cells:", ncol(expr_data), "\n");
  cat("Features:", nrow(expr_data), "\n")
}
### ####### ####### ##### ### ##### ##

cellsf <- show_found(colnames(expr_data), rownames(meta_data), v = opt$verbose)
meta_data <- meta_data[cellsf, , drop = FALSE]
expr_data <- expr_data[, cellsf]

# # This changes the clustering names and makes it harder to subset
# if("orig.ident" %in% colnames(meta_data)){
#   # if the meta data comes from a Seurat object, check if it has cluster columns
#   # previous project name and its clusters are renamed
#   if(opt$verbose) cat("Preserve previous analyses clusters\n")
#   for(orig in grep('^res|.*_snn_res', colnames(meta_data), value = TRUE)){
#     meta_data[, orig] <- paste0(meta_data[, "orig.ident"], "_", meta_data[, orig])
#   }
# }
# This line changes column names' case... it probably shouldn't!
# sub('^res.|.*_snn_res.', config$ident_vars[1], )
colnames(meta_data) <- make.names(names = colnames(meta_data), unique = TRUE)
tvar <- grep(pattern = "orig.ident|seurat_clusters", x = colnames(meta_data), value = TRUE, invert = TRUE)
if(opt$verbose) cat("Keeping", length(tvar), "of", ncol(meta_data), "columns\n")
meta_data <- meta_data[, tvar, drop = FALSE]
#### ####### ######### ####-----------------------------------------------------

#### Adding variables of interest ####------------------------------------------
if(opt$verbose) cat("\n\nAdding QC metrics from features (percentages of counts)\n")
meta_data <- add_pcts(
  mdata = meta_data,
  edata = expr_data,
  feature_pattern = list(
    percent.mt = c('^mt-', '^m-', '^hg19_mt-', '^mm10_mt-'),
    percent.ribo = paste0(c('^rps', '^rpl'), "[0-9]"),
    percent.hs = paste0("^hsp", letters[1:5])
  ),
  v = TRUE
)
meta_data$nFeature_RNA <- Matrix::colSums(expr_data[, rownames(meta_data)] > 0)
meta_data$nCount_RNA <- Matrix::colSums(expr_data[, rownames(meta_data)])
meta_data$Data = casefold(x = config$project_name, upper = TRUE)
#### ###### ######### ## ######## ####------------------------------------------
if(opt$verbose){
  cat("Presenting data:\nAnnotation:\n"); str(meta_data)
  cat("Matrix:"); str(expr_data)
  cat("Aprox. expression range"); print(summary(as.matrix(expr_data[, sample(1:ncol(expr_data), 4)])))
}; fname <- paste0('metadata_prefilter.rdata'); save(meta_data, file = fname)

#### Filtering variables ####---------------------------------------------------
if(opt$verbose) cat('\n\nDistribution of filtering criteria\n')
qc_filters = lapply(qc_filters, function(i) i[names(i) %in% colnames(meta_data)] )
meta_data[, names(qc_filters[[1]])] <- apply(meta_data[, names(qc_filters[[1]]), drop = FALSE], 2, function(x){
  y <- x; y[is.na(y)] <- Inf; return(y); #any(is.na(x))
})
summ <- sapply(meta_data[, names(qc_filters[[1]]), drop = FALSE], function(x){
  c(rev(summary(x)),
    quantile(x[!is.na(x)], prob = c(.001, 0.1, 0.9, 0.992)))
})
tvar <- t(data.frame(row.names = names(qc_filters[[1]]), qc_filters$low, qc_filters$high))
tvar <- ifelse(is.infinite(tvar), summ[c('Min.', 'Max.'), ], tvar)
print(rbind(summ, tvar))
mytab = rbind(summ, tvar)
rs = c(
  "qc_filters.low", "qc_filters.high", "Min.", "0.1%", "10%", "1st Qu.",
  "Median", "Mean", "3rd Qu.", "90%", "99.2%", "Max."
)
rs = rs[rs %in% rownames(mytab)]
if(length(rs) > (nrow(mytab) / 2)) mytab = mytab[rs, ]
write.csv(mytab, file = "filters_summary.csv")
summ <- rbind(summ[c("0.1%", "99.2%"), ], tvar) # change infinites for boundaries and merge quantiles
head(summ)

tvar <- sapply(1:length(names(qc_filters[[1]])), function(x){ # checking for metric unique values
  length(table(meta_data[, names(qc_filters[[1]])[x]])) < 3
})
if(opt$verbose) cat("Removing parameters with constant value accrosss cells:", sum(tvar), "\n")
qc_filters$low <- qc_filters$low[!tvar]; qc_filters$high <- qc_filters$high[!tvar]
names(qc_filters$low) <- names(qc_filters$high) <- names(qc_filters[[1]]) <- names(qc_filters[[1]])[!tvar]
# Checking for automatic/default filtering
tvar <- list(names(qc_filters[[1]]), c("qc_filters$low", "qc_filters$high"))
default_filters <- matrix(nrow = length(tvar[[1]]), ncol = length(tvar[[2]]), dimnames = tvar)
default_filters[, "qc_filters$low"] <- qc_filters$low[c(names(qc_filters[[1]]))]
default_filters[, "qc_filters$high"] <- qc_filters$high[c(names(qc_filters[[1]]))]
# Boundaries when indicated (-1)
qc_filters$low <- ifelse(qc_filters$low == -1, summ[1, names(qc_filters[[1]])], qc_filters$low)
qc_filters$high <- ifelse(qc_filters$high == -1, summ[2, names(qc_filters[[1]])], qc_filters$high)
# Defaults when indicated (-2)
qc_filters$low <- ifelse(qc_filters$low == -2, default_filters[names(qc_filters[[1]]), 1], qc_filters$low)
qc_filters$high <- ifelse(qc_filters$high == -2, default_filters[names(qc_filters[[1]]), 2], qc_filters$high)
filtersdf <- data.frame(qc_filters$low, qc_filters$high, qc_filters$apply)
summ[3:4, ] <- t(filtersdf[, -3])
fname <- paste0('filters_table.csv'); if(!file.exists(fname)) write.csv(filtersdf, file = fname)
print(filtersdf)

if(opt$verbose) cat("\n\nGroups found in data:\n")
tvar <- paste0("^RNA|", config$ident_vars[1])
orignames <- grep(pattern = tvar, x = colnames(meta_data), value = TRUE)
tvar <- lapply(meta_data[, orignames, drop = FALSE], table, useNA = 'always')
print(tvar); orignames <- orignames[sapply(tvar, length) > 1]

### Classifying cell quality ###------------------------------------------------
if(opt$verbose){
  cat("\n\nClassification of cell quality\n")
  cat("--- Zhang Lab; Guo et al 2018\n")
  # they removed ribosomal reads before mapping to the reference
}
class_params <- c('nFeature_RNA', 'nCount_RNA', 'percent.mt')
tvar <- t(sapply(class_params, function(x){
  if(x == "percent.mt") return(c(-Inf, 10))
  y <- (3 * mad(meta_data[, x], constant = 1, low = TRUE))
  c(median(meta_data[, x]) - y, median(meta_data[, x]) + y)
}))# average count > 1
meta_data$qc_tag_guo <- "Good"
for(i in 1:nrow(tvar)){
  meta_data$qc_tag_guo[meta_data[, names(tvar[i, 2])] > tvar[i, 2]] <- "Bad"
  meta_data$qc_tag_guo[meta_data[, names(tvar[i, 1])] < tvar[i, 1]] <- "Bad"
}; table(meta_data$qc_tag_guo)
meta_data$orig.qc_tag <- "ALL"

## Classification — Greg's
# filtered_reads became final_STAR_counts
class_params <- c('Total_genes', 'final_STAR_counts', 'uniquely_mapped_reads_perc',
'bias_5to3_prim', 'percent.mt')
plate_name = "orig.Plate"
if(all(class_params %in% colnames(meta_data))){
  if(opt$verbose) cat("--- Vijay Lab: Smart-Seq2\n")
  meta_data$orig.qc_tag <- "Manual"
  tvar <- meta_data$Total_genes < 200 & meta_data$final_STAR_counts >= 50000
  meta_data$orig.qc_tag[tvar] <- "1.Bad"
  reseq_def <- sapply(unique(meta_data$orig.Plate), function(x){
    y <- meta_data[meta_data[, plate_name] == x, ]
    tvar <- y$Total_genes < 200 & y$final_STAR_counts < 50000 &
      y$uniquely_mapped_reads_perc >= 60
    (sum(tvar) / nrow(y) * 100)
  })
  meta_data$reseq_def <- reseq_def[meta_data[, plate_name]]
  meta_data$orig.qc_tag[meta_data$reseq_def >= 20] <- "Reseq_plate"
  meta_data$orig.qc_tag[meta_data$reseq_def < 20] <- "2.Bad"
  tvar <- meta_data$Total_genes >= 200 & meta_data$final_STAR_counts >= 50000 &
    meta_data$bias_5to3_prim <= 2
  tvar <- tvar & meta_data$percent.mt <= 20 &
    meta_data$uniquely_mapped_reads_perc >= 60
  meta_data$orig.qc_tag[tvar] <- "Good"
  meta_data = meta_data[, !colnames(meta_data) %in% c('reseq_def')]
  print(table(meta_data$orig.qc_tag))
  table(meta_data$orig.qc_tag, meta_data$qc_tag_guo)
}; fname <- paste0('metadata_prefilter.rdata'); save(meta_data, file = fname)
### ########### #### ####### ###------------------------------------------------

# #### Qlucore file ####----------------------------------------------------------
# if(opt$verbose) cat("\n\n----- Qlucore input -------\n")
# passed <- if(is.null(meta_data$orig.qc_tag)){
#   rownames(meta_data)
# }else{
#   if(opt$verbose) cat("Using only the 'good' ones.\n")
#   rownames(meta_data[meta_data$orig.qc_tag == "Good", ])
# }
# tvar <- gsub("raw\\.", "TPM\\.", config$input_expression)
# tvar <- if(!file.exists(tvar)) gsub("_raw", "_TPM", config$input_expression) else tvar
# suffixtag = "_counts"
# expr_data_t <- if(file.exists(tvar) && grepl("raw|TPM", tvar, ignore.case = TRUE)){
#   suffixtag = "_TPM" # DATA DUPLICATION
#   expr_data_t <- readfile(tvar)
# }else{ expr_data }
# source("/home/ciro/scripts/handy_functions/R/qlucore_file.R")
# for(thresh in c(0, 5)){
#   gmeans <- sort(Matrix::rowMeans(expr_data_t[, passed]), decreasing = TRUE)
#   genes <- rownames(expr_data_t[names(gmeans[gmeans > thresh]), ])
#   if(opt$verbose) cat(length(genes), "of", nrow(expr_data_t), 'features\n')
#   mysamples <- sample_grp(annot = meta_data, cname = 'Data', maxln = '5000')
#   if(opt$verbose) cat(length(mysamples), "of", ncol(expr_data_t), 'samples\n')
#   suffix <- paste0("qlucore", if(is.numeric(thresh)){
#     paste0("_gt", thresh, "mean_")
#   }else{ paste0("_", thresh, "_")) }
#   suffix <- paste0(suffix, length(genes), "genes", nrow(meta_data), "samples")
#   suffix <- paste0(suffix, suffixtag, ".txt")
#   nameout <- paste0(suffix)
#   if(opt$verbose) cat("Qlucore file:", nameout, "\n")
#   if(!file.exists(nameout)){
#     qf <- qlucore_format(
#       mat = expr_data_t, metadata = meta_data, rnames = genes, cnames = , v = TRUE)
#     write.table(qf, file = nameout, sep = "\t", quote = FALSE, row.names = FALSE)
#   }else{
#     if(opt$verbose) cat("Exists\n")
#   }
# }; if(opt$verbose) cat(list.files(path = "./", pattern = "qlucore"), sep = "\n")
# #### ####### #### ####----------------------------------------------------------

if(opt$verbose)
  cat('----------------------- Plots and final filters -----------------\n')
if(opt$verbose) cat("Scatter filters\n")
qc_scatters(dat = meta_data, thresholds = summ[, qc_filters$apply], prefix = "0_qc")

if(!exists('meta_databk')) meta_databk <- meta_data
vlns <- lapply(names(qc_filters[[1]]), function(x){
  qc_violin(
    dat = meta_data, xax = "orig.qc_tag", yax = x,
    lb_filt = qc_filters$low, hb_filt = qc_filters$high, filtby = qc_filters$apply)
})
is_simple <- length(table(meta_data$orig.qc_tag)) < 3 & length(vlns) < 4
mygrid <- make_grid(length(names(qc_filters[[1]])), ncol = if(is_simple) 3)
pdf(paste0('0_qc_all_qc_tag.pdf'), width = 5.5 * mygrid[2], height = 5 * mygrid[1])
print(cowplot::plot_grid(plotlist = vlns, ncol = mygrid[2]))
graphics.off()

for(prefix in c("1_qc_all", "2_filtered")){
  if(opt$verbose) cat("=======================", prefix, "\n")
  vlns <- lapply(names(qc_filters[[1]]), function(x){
    qc_violin(dat = meta_data, yax = x,
      lb_filt = qc_filters$low, hb_filt = qc_filters$high, filtby = qc_filters$apply)
  })
  is_simple <- length(table(meta_data$orig.qc_tag)) < 3 & length(vlns) < 4
  mygrid <- make_grid(length(names(qc_filters[[1]])), ncol = if(is_simple) 3)
  tvar <- ifelse(is_simple, 7, 5) * mygrid[1]
  pdf(paste0(prefix, '_data.pdf'), width = 5.5 * mygrid[2], height = tvar)
  print(cowplot::plot_grid(plotlist = vlns, ncol = mygrid[2]))
  graphics.off()
  orignamesp <- orignames[!orignames %in% names(qc_filters[[1]])]
  if(length(orignamesp) > 0){
    for(origs in orignamesp){
      if(opt$verbose) cat(paste0(origs, ", "))
      tvar <- if(origs != "origlib") sub(config$ident_vars[1], "", origs) else "lib"
      fname <- paste0(prefix, '_data_confounder_', tvar, '.pdf')
      qc_violin_list <- lapply(names(qc_filters[[1]]), function(x){
        qc_violin(
          dat = meta_data, xax = origs, yax = x,
          lb_filt = qc_filters$low,
          hb_filt = qc_filters$high,
          filtby = qc_filters$apply
        )
      }); sfact <- if(length(table(meta_data[, origs])) < 10) 5.5 else 3
      tvar <- ifelse(is_simple, 7, sfact) * mygrid[1]
      pdf(fname, width = sfact * mygrid[2], height = tvar)
      if(length(table(meta_data[, origs])) < 10){
        print(cowplot::plot_grid(plotlist = qc_violin_list, ncol = mygrid[2]))
      }else{
        tvar <- lapply(qc_violin_list, print)
      }; graphics.off()
    }; if(opt$verbose) cat("\n")
  }
  if(prefix == "1_qc_all"){
    if(opt$verbose) cat('Filtering cells\n')
    print(filtersdf[qc_filters$apply, 1:2])
    tvar <- names(qc_filters[[1]])[qc_filters$apply]
    all_filters <- paste(
      paste(tvar, ">=", qc_filters$low[qc_filters$apply], collapse = " & "),
      paste(tvar, "<=", qc_filters$high[qc_filters$apply], collapse = " & "),
    sep = " & ")
    if(opt$verbose) cat("Criteria:", all_filters, "\n")
    expr = grep(pattern = "^expr|addit", x = names(config$filtering), value = TRUE)
    if(length(config$filtering[[expr]]) > 0){
      tvar <- unlist(config$filtering[[expr]]);
      cat("Specific filters:", tvar, sep = "\n");
      tvar <- paste(paste0("(", tvar, ")", collapse = " & "))
      all_filters <- paste("(", tvar, ") & ", all_filters)
    }
    # !!!!!!!!!!!! choose wisely !!!!!!!!!!!!
    sset <- paste0("meta_data <- subset(x = meta_data, subset = ", all_filters, ")")
    tvar <- c(qc_filters$high[qc_filters$apply], qc_filters$low[qc_filters$apply])
    if(!all(is.infinite(tvar)) || length(config$filtering[[expr]]) > 0){
      if(opt$verbose) cat('- Before:', nrow(meta_databk), '\n')
      eval(expr = parse(text = sset))
      if(opt$verbose) cat('- After:', nrow(meta_data), '\n')
    }; if(nrow(meta_databk) == nrow(meta_data)) break
    if(opt$verbose) cat("------------- Per quality plots\n")
    if(is.null(meta_databk$orig.qc_tag)){
      meta_databk$orig.qc_tag <- "Good"
      meta_databk[!rownames(meta_databk) %in% rownames(meta_data), ]$orig.qc_tag <- "bad"
    }
    tvar <- unique(meta_databk$orig.qc_tag); if(length(tvar) < 2) next
    for(qc_tag in tvar){
      if(opt$verbose) cat("By quality\n");
      meta_data_ss <- meta_databk[meta_databk$orig.qc_tag == qc_tag, ]
      orignamesp <- if(plate_name %in% colnames(meta_data_ss)) plate_name else "origlib"
      if(length(orignamesp) > 0){
        for(origs in orignamesp){
          if(opt$verbose) cat(".")
          tvar <- if(origs != "origlib") sub(config$ident_vars[1], "", origs) else "lib"
          fname <- paste0(prefix, '_data_confounders_', qc_tag, tvar, '.pdf')
          qc_violin_list <- lapply(names(qc_filters[[1]]), function(x){
            qc_violin(
              dat = meta_data_ss, xax = origs, yax = x,
              lb_filt = qc_filters$low,
              hb_filt = qc_filters$high,
              filtby = qc_filters$apply
            )
          }); sfact <- if(length(table(meta_data_ss[, origs])) < 10) 5.5 else 3
          pdf(fname, width = sfact * mygrid[2], height = ifelse(is_simple, 7, sfact) * mygrid[1])
          if(length(table(meta_data_ss[, origs])) < 10){
            print(cowplot::plot_grid(plotlist = qc_violin_list, ncol = mygrid[2]))
          }else{
            tvar <- lapply(qc_violin_list, print)
          }; graphics.off()
        }; if(opt$verbose) cat("\n")
      }
    }
  } # Finish filtering
}
qc_scatters(dat = meta_data, thresholds = summ[, qc_filters$apply], prefix = "3_qc")
# ############# ####------------------------------------------------------------

save(meta_data, file = "metadata_filtered.rdata")
# write.table(
#   rownames(expr_data), file = "feature_names.txt",
#   quote = FALSE, row.names = FALSE, col.names = FALSE);
#### ######### ######### ####---------------------------------------------------

if(opt$verbose){
  cat('\n\n*******************************************************************\n')
  cat('Starting time:\n'); cat(st.time, '\n')
  cat('Finishing time:\n'); timestamp()
  cat('*******************************************************************\n')
  cat('SESSION INFO:\n'); print(sessionInfo()); cat("\n")
  cat('Pipeline finished successfully\n')
}
if(register_log){
  sink(type = 'message')
  sink()
}
