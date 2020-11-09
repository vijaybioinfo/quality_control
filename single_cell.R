#!/usr/bin/R

#######################################
# Quality Control of single-cell data #
#######################################

# This script contains the QC steps for single-cell analysis and gets data
# ready for downstream protocols; clustering, scDGEA, etc.

library(optparse)
library(yaml)

optlist <- list(
  make_option(
    opt_str = c("-y", "--yaml"), type = "character",
    help = "Configuration file: Instructions in YAML format."
  ),
  make_option(
    opt_str = c("-v", "--verbose"), default = TRUE,
    help = "Verbose: Show progress."
  )
)

# Getting arguments from command line and setting their values to their respective variable names.
optparse <- OptionParser(option_list = optlist)
defargs <- parse_args(optparse)
if(interactive()){ # Example/manually
  defargs$yaml = "/home/ciro/scripts/clustering/qc.yaml"
}
config = read_yaml(defargs$yaml)
str(config)
options(warnings = FALSE)
options(warning = FALSE)

#### Functions ####-------------------------------------------------------------
cat("Loading functions\n")
resources = c(
  "https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/devel/utilities.R", # dircheck show_parameters newlines
  "https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/devel/file_reading.R", # readfile
  "https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/R/clustering_utilities.R", # get_source_data add_pcts
  "https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/devel/filters.R", # filters_complex
  "https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/R/clustering_plotting.R", # feature_scatter qc_violin
  "https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/devel/plots.R", # make_grid
  # Soon to be changed # ---
  "https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/devel/plotting_dismantle.R" # plot_corr simple_violin
)
for(i in resources){ source(i) }
qc_scatters <- function(
  dat,
  thresholds,
  prefix = "0_qc"
){
  apfilters <- unique(c(colnames(thresholds), "final_STAR_counts"))
  apfilters <- gtools::combinations(
    length(apfilters), r = 2,
    v = apfilters,
    set = TRUE, repeats.allowed = FALSE
  )
  for(i in 1:nrow(apfilters)){
    axes <- c(apfilters[i, ], "uniquely_mapped_reads_perc") # var. of interes
    axes <- unique(c(axes, "percent.mt")) # var. of interes
    axes <- axes[axes %in% colnames(dat)]
    axes <- axes[axes %in% colnames(thresholds)]
    if(length(axes) < 2) next
    cat(show_commas(axes))
    p <- try(feature_scatter(
      object = dat,
      variables = axes,
      thresh = thresholds[, axes],
      verbose = !TRUE
    ), silent = TRUE)
    if(class(p)[1] != "gg"){ cat(" - failed\n"); next }; cat(" - plotting\n")
    fname <- paste0(c(prefix, axes, 'all_data.pdf'), collapse = "_")
    pdf(fname, width = 10, height = 10)
    print(try(p, silent = TRUE))
    dev.off()
  }
}
#### ######### ####-------------------------------------------------------------
output_dir <- dircheck(paste0(dircheck(config$output_dir), config$project_name, "/"))
setwd(output_dir)

#### Directories structure ####-------------------------------------------------
# defargs$output_dir <- paste0(defargs$output_dir, "qc/")
#### ########### ######### ####-------------------------------------------------
### Log file ####---------------------------------------------------------------
zetinfo <- ""
log_file <- paste0(output_dir, 'a1_qc.log')
register_log <- !interactive() && TRUE
if(register_log){
  if(!file.exists(log_file)) file.create(log_file)
  out.file <- file(log_file, open = 'wt')
  sink(out.file) ; sink(out.file, type = 'message')
}
cat('Date and time:\n') ; st.time <- timestamp();
### ### #### ####---------------------------------------------------------------

#### Loading dependencies ####--------------------------------------------------
packages <- c("ggplot2", "cowplot")
loaded <- lapply(X = packages, FUN = function(x){
  suppressMessages(require(package = x, quietly = TRUE, character.only = TRUE))
})
theme_set(theme_cowplot())
#### ####### ############ ####--------------------------------------------------

#### Presenting Parameters ####-------------------------------------------------
str(config)
## Digesting some first
cat('\n\n----------------------- Digesting parameters --------------------\n')
mytab <- if(file.exists(config$filtering$file)){ # QC filters first
  readfile(config$filtering$file, stringsAsFactors = FALSE)
}else{
  tvar <- names(config$filtering)[!names(config$filtering) %in% c("file", "subset")]
  data.frame(sapply(config$filtering[names(config$filtering) %in% tvar], function(x) as.numeric(unlist(x)[1:3]) ))
}
mytab[is.na(mytab)] <- -1; # gets quantiles
ivariables <- colnames(mytab)
low.filt.cells <- unlist(mytab[1, ]); high.filt.cells <- unlist(mytab[2, ])
apply_filter <- unlist(mytab[3, ])
# apply all of them by default; except when some of them are specified
apply_filter[is.na(apply_filter)] <- ifelse(all(is.na(apply_filter)), 1, 0)
apply_filter <- apply_filter > 0
fix_boundary <- function(x, y, replacement){ x <- x[y]; x[is.na(x)] <- replacement; names(x) <- y; x }
low.filt.cells <- fix_boundary(x = low.filt.cells, y = ivariables, -Inf)
high.filt.cells <- fix_boundary(x = high.filt.cells, y = ivariables, Inf)
#### ########## ########## ####-------------------------------------------------

#### Getting data ####----------------------------------------------------------
cat('----------------------- Loading data -----------------------\n')
cat('======================= Expression matrix\n')
# Find whatever you give, 10X directory (barcodes, gene_names, counts), Seurat object, CSV file, TXT, RDS
is_vector_counts <- grepl("meco", config$project_name) | !dir.exists(config$input_expression)
Sys.time()
mycellsdata <- get_source_data(config$input_expression, pj = config$project_name, merge_counts = is_vector_counts, v = TRUE)
Sys.time()
config$input_expression <- mycellsdata[[2]]
if(length(mycellsdata) > 2) annottab <- mycellsdata[[3]]
mycellsdata <- if(length(mycellsdata) > 1 && is.null(dim(mycellsdata[[1]]))){
  cat("Keeping only expression data to store in another assay\n")
  mycellsdata[[grep("expression", names(mycellsdata), ignore.case = TRUE)]]
}else{ mycellsdata[[1]] }
if(is.data.frame(mycellsdata)) mycellsdata <- Matrix(as.matrix(mycellsdata), sparse = TRUE)
str(mycellsdata, max.level = 2); gc()
cat('======================= Annotation/metadata\n')
# This file must be a table with the following structure
# Library, Condition_1, Condition_2, etc, or it can be the annotation/meta.data
addmetadataf <- if(length(config$metadata[-1]) > 0) config$metadata[-1] else "no_file"
grpsamples <- if(file.exists(config$metadata)){
  readfile(config$metadata[1], row.names = 1, stringsAsFactors = FALSE, check.names = FALSE, header = TRUE)
}else{
  data.frame(result = 'none', stringsAsFactors = FALSE)
}
if(any(grepl("cell_id", colnames(grpsamples)))) rownames(grpsamples) <- grpsamples$cell_id
tvar <- apply(grpsamples, 2, function(x){ # find which column contains barcodes
  if(any(x %in% colnames(mycellsdata))) TRUE else FALSE # then check if they're in the rownames
}); tvar['rowname_barcode'] <- any(rownames(grpsamples) %in% colnames(mycellsdata))
if(isTRUE(tvar['rowname_barcode'])) cat("Sample names in rows\n")
if(any(tvar)){ # check if any of the columns in the table has the cell names (barcodes)
  cat('- provided from file\n');
  annottab <- grpsamples; rm(grpsamples) # this could mean that it's the annotation itself!
  # if the barcodes are not in the rownames, then find the column with them
  if(!tvar['rowname_barcode']) rownames(annottab) <- annottab[, head(names(tvar[tvar]), 1)]
}else if(exists('annottab')){
  cat('- given in expression input\n');
  print(str(annottab))
}
if(!exists('annottab')){
  cat('- building table...\n');
  annottab <- data.frame(row.names = colnames(mycellsdata), stringsAsFactors = FALSE, check.names = FALSE)
  if(grepl('-', colnames(mycellsdata)[1])){ # You do need to find aggregation_csv to know the order of the libraries
    # You do need to find aggregation_csv to know the order of the libraries
    aggname <- if(!dir.exists(config$input_expression)) basename(config$input_expression) else 'aggregation*.csv'
    input_tab <- if(!dir.exists(config$input_expression)) dirname(config$input_expression) else config$input_expression
    aggregated <- findfile(name = aggname, path = input_tab, read = TRUE, v = TRUE, stringsAsFactors = FALSE)
    cell_gem <- sub('.*\\-', '', colnames(mycellsdata), perl = TRUE) # getting cells GEM
    cell_lib <- rownames(aggregated); names(cell_lib) <- aggregated[, 1] # and vector named as the library name
    cat('**\nAssigning names per library:\n', show_commas(cell_lib), '\n', show_commas(names(cell_lib)), '\n***\n')

    grpsamples$its_row_names <- rownames(grpsamples)
    libname <- names(which.max(sapply(grpsamples, function(x) sum(x %in% names(cell_lib)) )))[1]
    cat("Libraries column:", libname, "\n")
    rownames(grpsamples) <- grpsamples[, libname] # taking the library names as rownames for ordering

    if(!all(names(cell_lib) %in% rownames(grpsamples))){
      cat("Missing:", show_commas(names(cell_lib)[!names(cell_lib) %in% rownames(grpsamples)]), '\n')
      cat("Given:", show_commas(grpsamples[!rownames(grpsamples) %in% names(cell_lib), libname]), '\n')
    }
    grpsamples_selected <- grpsamples <- grpsamples[names(cell_lib), , drop = FALSE] # ordering accorgin to GEM from aggregation file
    rownames(grpsamples) <- NULL; colnames(grpsamples) <- gsub("its_row_names", "Library", colnames(grpsamples))
    grpsamples_selected <- grpsamples_selected[, !colnames(grpsamples_selected) %in% "its_row_names", drop = FALSE]
    cell_lib <- rownames(grpsamples_selected) # witching to GEM as names in the named vector
    names(cell_lib) <- rownames(grpsamples_selected) <- rownames(aggregated)
    group_names <- grpsamples_selected

    print(data.frame(GEM = rownames(aggregated), Library = cell_lib), row.names = F)
    cat('\nAssigning group names\n')
    print(group_names)
    track_names <- which(sapply(group_names, is.character) & sapply(group_names, function(x) length(table(x)) ) > 1)
    cat("Traking variables:", show_commas(names(track_names)), "\n")
    colnames(group_names)[track_names] <- paste0('orig.', colnames(group_names)[track_names])
    colnames(group_names) <- make.names(casefold(colnames(group_names)))
    annottab <- remove.factors(cbind(annottab, origlib = cell_lib[cell_gem], group_names[cell_gem, , drop = FALSE]))
    write.csv(grpsamples, file = paste0(zetinfo, 'groups.csv'), row.names = FALSE)
  }else if(length(nrow(grpsamples))){
    cat('No pattern in cell names:', show_commas(colnames(mycellsdata)), '\n')
    if(grpsamples[1, 1] != "void")
      cat('Library names:', show_commas(grpsamples[, ifelse(is.null(ncol(grpsamples)), 1, 2)]), '\n')
    annottab$origlib <- rep('L1', ncol(mycellsdata))
  }
};

annottab <- remove.factors(annottab)
addmetadataf <- unlist(strsplit(path.expand(addmetadataf), "~"))
tvar <- file.exists(addmetadataf)
str(addmetadataf)
if(any(tvar)){
  addmetadataf <- addmetadataf[tvar]
  cat("Adding metadata from:", addmetadataf, sep = "\n")
  for(addmetadataf_i in addmetadataf){
    addannot <- remove.factors(readfile(addmetadataf_i))
    print(dim(annottab))
    # partial matching problem addressed in /home/ciro/covid19/scripts/partial_matching.R
    annottab <- joindf(x = annottab, y = addannot)
  }
  print(dim(annottab))
}; fname <- paste0(zetinfo, 'metadata.rdata'); if(!file.exists(fname)) save(annottab, file = fname)
#### ####### #### ####----------------------------------------------------------

#### Initial Filtering ####-----------------------------------------------------
filtereddata <- filters_complex(
  mdata = annottab,
  filters = lapply(names(config$filtering$subset), function(x) c(x, config$filtering$subset[[x]]) ),
  v = TRUE
)
annottab <- filtereddata[[1]]; rm(filtereddata)

### Minimal filter: cells and genes ###
mincs = ceiling(0.001 * ncol(mycellsdata))
tvar <- which(ivariables == 'nFeature_RNA')
mings = ifelse(isTRUE(low.filt.cells[tvar] > 0), low.filt.cells[tvar], 200)
cat("Calculating number of expressed features and expressing samples.\n")
# Slow and will probably be re-calculated
efeatures <- Matrix::colSums(mycellsdata > 0); esamples <- Matrix::rowSums(mycellsdata > 0)
# head(reshape2::melt(sort(efeatures[efeatures > 0])))
# head(reshape2::melt(sort(esamples[esamples > 0])))
summary(efeatures); summary(esamples)
# mycellsdata <- mycellsdata[esamples >= mincs, efeatures >= mings]
cat('Minimal filtering: ~0.1% of the data min.cells >=', mincs, '& min.genes (not applied) >=', mings,'\n')
cat("Keeping genes:", sum(esamples >= mincs), "/", length(esamples) , "\n")
mycellsdata <- mycellsdata[esamples >= mincs, ]
### ####### ####### ##### ### ##### ##

cellsf <- show_found(colnames(mycellsdata), rownames(annottab), element = 'samples/cells', v = TRUE)
annottab <- annottab[cellsf, , drop = FALSE]
mycellsdata <- mycellsdata[, cellsf]

if("orig.ident" %in% colnames(annottab)){ # if the meta data comes from a Seurat object, check if it has cluster columns
  cat("Preserve previous analyses clusters\n")
  for(orig in grepatrn(colnames(annottab), accept = '^res|.*_snn_res')){
    # pastes previous project name and its clusters
    annottab[, orig] <- paste0(annottab[, "orig.ident"], "_", annottab[, orig])
  }
}; colnames(annottab) <- casefold(make.names(sub('^res|.*_snn_res', 'orig', colnames(annottab)), unique = TRUE))
tvar <- grep(pattern = "orig.ident|seurat_clusters", x = colnames(annottab), value = TRUE, invert = TRUE)
cat("Keeping", length(tvar), "of", ncol(annottab), "columns\n")
annottab <- remove.factors(annottab[, tvar, drop = FALSE])
#### ####### ######### ####-----------------------------------------------------

#### Adding variables of interest ####------------------------------------------
cat("\n\nAdding QC metrics from features (percentages of counts)\n")
annottab <- add_pcts(
  mdata = annottab,
  edata = mycellsdata,
  feature_pattern = list(
    percent.mt = c('^mt-', '^m-', '^hg19_mt-', '^mm10_mt-'),
    percent.ribo = paste0(c('^rps', '^rpl'), "[0-9]"),
    percent.hs = paste0("^hsp", letters[1:5])
  ),
  v = TRUE
)
annottab$nFeature_RNA <- Matrix::colSums(mycellsdata[, rownames(annottab)] > 0)
annottab$nCount_RNA <- Matrix::colSums(mycellsdata[, rownames(annottab)])
annottab$Data = casefold(x = config$project_name, upper = TRUE)
#### ###### ######### ## ######## ####------------------------------------------
if(1){
  cat("Presenting data:\nAnnotation:\n"); str(annottab[, headtail(1:ncol(annottab), 5)])
  cat("Matrix:"); str(mycellsdata)
  cat("Aprox. expression range"); print(summary(as.matrix(mycellsdata[, sample(1:ncol(mycellsdata), 4)])))
}; fname <- paste0(zetinfo, 'metadata_prefilter.rdata'); save(annottab, file = fname)

#### Filtering variables ####---------------------------------------------------
cat('\n\nDistribution of filtering criteria\n')
annottab[, ivariables] <- apply(annottab[, ivariables, drop = FALSE], 2, function(x){
  y <- x; y[is.na(y)] <- Inf; return(y); #any(is.na(x))
})
summ <- sapply(annottab[, ivariables, drop = FALSE], function(x){
  c(rev(summary(x)),
    quantile(x[!is.na(x)], prob = c(.001, 0.992)))
})
# print(summ[c('Min.', 'Max.'), ]) # adding selected thresholds
tvar <- t(data.frame(row.names = ivariables, low.filt.cells, high.filt.cells))
tvar <- ifelse(is.infinite(tvar), summ[c('Min.', 'Max.'), ], tvar)
print(rbind(summ, tvar))
write.csv(rbind(summ, tvar), file = "filters_summary.csv")
summ <- rbind(summ[7:8, ], tvar) # change infinites for boundaries and merge quantiles
head(summ)

tvar <- sapply(1:length(ivariables), function(x){ # checking for metric unique values
  length(table(annottab[, ivariables[x]])) < 3
}); cat("Removing parameters with constant value accrosss cells:", sum(tvar), "\n")
low.filt.cells <- low.filt.cells[!tvar]; high.filt.cells <- high.filt.cells[!tvar]
names(low.filt.cells) <- names(high.filt.cells) <- ivariables <- ivariables[!tvar]
# Checking for automatic/default filtering
tvar <- list(ivariables, c("low.filt.cells", "high.filt.cells"))
default_filters <- matrix(nrow = length(tvar[[1]]), ncol = length(tvar[[2]]), dimnames = tvar)
default_filters[, "low.filt.cells"] <- low.filt.cells[c(ivariables)]
default_filters[, "high.filt.cells"] <- high.filt.cells[c(ivariables)]
# Boundaries when indicated (-1)
low.filt.cells <- ifelse(low.filt.cells == -1, summ[1, ivariables], low.filt.cells)
high.filt.cells <- ifelse(high.filt.cells == -1, summ[2, ivariables], high.filt.cells)
# Defaults when indicated (-2)
low.filt.cells <- ifelse(low.filt.cells == -2, default_filters[ivariables, 1], low.filt.cells)
high.filt.cells <- ifelse(high.filt.cells == -2, default_filters[ivariables, 2], high.filt.cells)
filtersdf <- data.frame(low.filt.cells, high.filt.cells, apply_filter)
summ[3:4, ] <- t(filtersdf[, -3])
fname <- paste0(zetinfo, 'filters_table.csv'); if(!file.exists(fname)) write.csv(filtersdf, file = fname)
print(filtersdf)

cat("\n\nGroups found in data:\n")
orignames <- grep(pattern = "orig", x = colnames(annottab), value = TRUE)
tvar <- lapply(annottab[, orignames], table, useNA = 'always')
print(tvar); orignames <- orignames[sapply(tvar, length) > 1]

### Classifying cell quality ###------------------------------------------------
cat("\n\nClassification of cell quality\n")
cat("--- Zhang Lab; Guo et al 2018\n") # they removed ribosomal reads before mapping to the reference
class_params <- c('nFeature_RNA', 'nCount_RNA', 'percent.mt')
tvar <- t(sapply(class_params, function(x){
  if(x == "percent.mt") return(c(-Inf, 10))
  y <- (3 * mad(annottab[, x], constant = 1, low = TRUE))
  c(median(annottab[, x]) - y, median(annottab[, x]) + y)
}))# average count > 1
annottab$orig.qc_tag_guo <- "Good"
for(i in 1:nrow(tvar)){
  annottab$orig.qc_tag_guo[annottab[, names(tvar[i, 2])] > tvar[i, 2]] <- "Bad"
  annottab$orig.qc_tag_guo[annottab[, names(tvar[i, 1])] < tvar[i, 1]] <- "Bad"
}; table(annottab$orig.qc_tag_guo)
annottab$orig.qc_tag <- "ALL"

## Classification — Greg's
# filtered_reads became final_STAR_counts
class_params <- c('Total_genes', 'final_STAR_counts', 'uniquely_mapped_reads_perc', 'bias_5to3_prim', 'percent.mt')
plate_name = "orig.Plate"
if(all(class_params %in% colnames(annottab))){
  cat("--- Vijay Lab: Smart-Seq2\n")
  annottab$orig.qc_tag <- "Manual"
  annottab$orig.qc_tag[annottab$Total_genes < 200 & annottab$final_STAR_counts >= 50000] <- "1.Bad"
  reseq_def <- sapply(unique(annottab$orig.Plate), function(x){
    y <- annottab[annottab[, plate_name] == x, ]
    tvar <- y$Total_genes < 200 & y$final_STAR_counts < 50000 & y$uniquely_mapped_reads_perc >= 60
    (sum(tvar) / nrow(y) * 100)
  })
  annottab$reseq_def <- reseq_def[annottab[, plate_name]]
  annottab$orig.qc_tag[annottab$reseq_def >= 20] <- "Reseq_plate"
  annottab$orig.qc_tag[annottab$reseq_def < 20] <- "2.Bad"
  tvar <- annottab$Total_genes >= 200 & annottab$final_STAR_counts >= 50000 & annottab$bias_5to3_prim <= 2
  tvar <- tvar & annottab$percent.mt <= 20 & annottab$uniquely_mapped_reads_perc >= 60
  annottab$orig.qc_tag[tvar] <- "Good"
  annottab = annottab[, !colnames(annottab) %in% c('reseq_def')]
  print(table(annottab$orig.qc_tag))
  table(annottab$orig.qc_tag, annottab$orig.qc_tag_guo)
}; fname <- paste0(zetinfo, 'metadata_prefilter.rdata'); save(annottab, file = fname)
### ########### #### ####### ###------------------------------------------------

# #### Qlucore file ####----------------------------------------------------------
# cat("\n\n----- Qlucore input -------\n")
# passed <- if(is.null(annottab$orig.qc_tag)){
#   rownames(annottab)
# }else{
#   cat("Using only the 'good' ones.\n")
#   rownames(annottab[annottab$orig.qc_tag == "Good", ])
# }
# tvar <- gsub("raw\\.", "TPM\\.", config$input_expression)
# tvar <- if(!file.exists(tvar)) gsub("_raw", "_TPM", config$input_expression) else tvar
# suffixtag = "_counts"
# mycellsdata_t <- if(file.exists(tvar) && grepl("raw|TPM", tvar, ignore.case = TRUE)){ # DATA DUPLICATION
#   suffixtag = "_TPM"
#   mycellsdata_t <- readfile(tvar)
# }else{ mycellsdata }
# source("https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/R/qlucore_file.R")
# for(thresh in c(0, 5)){
#   gmeans <- sort(Matrix::rowMeans(mycellsdata_t[, passed]), decreasing = TRUE)
#   genes <- rownames(mycellsdata_t[names(gmeans[gmeans > thresh]), ])
#   cat(length(genes), "of", nrow(mycellsdata_t), 'features\n')
#   mysamples <- sample_grp(annot = annottab, cname = 'Data', maxln = '5000')
#   cat(length(mysamples), "of", ncol(mycellsdata_t), 'samples\n')
#   suffix <- paste0(zetinfo, "qlucore", if(is.numeric(thresh)) paste0("_gt", thresh, "mean_") else paste0("_", thresh, "_"))
#   suffix <- paste0(suffix, length(genes), "genes", nrow(annottab), "samples")
#   suffix <- paste0(suffix, suffixtag, ".txt")
#   nameout <- paste0(suffix)
#   cat("Qlucore file:", nameout, "\n")
#   if(!file.exists(nameout)){
#     qf <- qlucore_format(mat = mycellsdata_t, metadata = annottab, rnames = genes, cnames = , v = TRUE)
#     write.table(qf, file = nameout, sep = "\t", quote = FALSE, row.names = FALSE)
#   }else{
#     cat("Exists\n")
#   }
# }; cat(list.files(path = zetinfo, pattern = "qlucore"), sep = "\n")
# #### ####### #### ####----------------------------------------------------------

# visualisation ####------------------------------------------------------------
cat("Scatter filters\n")
qc_scatters(dat = annottab, thresholds = summ[, apply_filter], prefix = "0_qc")

if(!exists('annottabbk')) annottabbk <- annottab
vlns <- lapply(ivariables, function(x){
  qc_violin(mydat = annottab, xax = "orig.qc_tag", yax = x, lb_filt = low.filt.cells, hb_filt = high.filt.cells, filtby = apply_filter)
})
is_simple <- length(table(annottab$orig.qc_tag)) < 3 & length(vlns) < 4
mygrid <- make_grid(length(ivariables), ncol = if(is_simple) 3)
pdf(paste0('0_qc_all_qc_tag.pdf'), width = 5.5 * mygrid[2], height = 5 * mygrid[1])
print(cowplot::plot_grid(plotlist = vlns, ncol = mygrid[2]))
graphics.off()

for(prefix in c("1_qc", "2_filtered")){
  cat("=======================", prefix, "\n")
  vlns <- lapply(ivariables, function(x){
    qc_violin(mydat = annottab, yax = x, lb_filt = low.filt.cells, hb_filt = high.filt.cells, filtby = apply_filter)
  })
  is_simple <- length(table(annottab$orig.qc_tag)) < 3 & length(vlns) < 4
  mygrid <- make_grid(length(ivariables), ncol = if(is_simple) 3)
  pdf(paste0(prefix, '_all_data.pdf'), width = 5.5 * mygrid[2], height = 5 * mygrid[1])
  print(cowplot::plot_grid(plotlist = vlns, ncol = mygrid[2]))
  graphics.off()
  orignamesp <- orignames[!orignames %in% ivariables]
  if(length(orignamesp) > 0){
    for(origs in orignamesp){
      cat(paste0(origs, ", "))
      fname <- paste0(prefix, '_all_data_confounder', sub("orig", "", origs), '.pdf')
      void <- lapply(ivariables, function(x){
        qc_violin(
          mydat = annottab, xax = origs, yax = x,
          lb_filt = low.filt.cells,
          hb_filt = high.filt.cells,
          filtby = apply_filter
        )
      }); sfact <- if(length(table(annottab[, origs])) < 10) 5.5 else 3
      pdf(fname, width = sfact * mygrid[2], height = sfact * mygrid[1])
      if(length(table(annottab[, origs])) < 10){
        print(cowplot::plot_grid(plotlist = void, ncol = mygrid[2]))
      }else{
        tvar <- lapply(void, print)
      }; graphics.off()
    }; cat("\n")
  }
  if(prefix == "1_qc"){
    cat('Filtering cells\n')
    print(filtersdf[apply_filter, 1:2])
    tvar <- paste(
      paste(ivariables[apply_filter], ">=", low.filt.cells[apply_filter], collapse = " & "),
      paste(ivariables[apply_filter], "<=", high.filt.cells[apply_filter], collapse = " & "),
    sep = " & ")
    cat("Criteria:", tvar, "\n")
    sset <- paste0("annottab <- subset(x = annottab, subset = ", tvar, ")") # << choose wisely
    if(!all(is.infinite(c(high.filt.cells[apply_filter], low.filt.cells[apply_filter])))){
      cat('- Before:', nrow(annottabbk), '\n')
      eval(expr = parse(text = sset))
      cat('- After:', nrow(annottab), '\n')
    }; if(nrow(annottabbk) == nrow(annottab)) break
    cat("------------- Per quality plots\n")
    if(is.null(annottabbk$orig.qc_tag)){
      annottabbk$orig.qc_tag <- "Good"
      annottabbk[!rownames(annottabbk) %in% rownames(annottab), ]$orig.qc_tag <- "bad"
    }
    tvar <- unique(annottabbk$orig.qc_tag); tvar <- tvar[tvar != "Good"]; if(length(tvar) == 0) next
    for(qc_tag in tvar){
      annottab_ss <- annottabbk[annottabbk$orig.qc_tag == qc_tag, ]
      orignamesp <- if(plate_name %in% colnames(annottab_ss)) plate_name else "origlib"
      if(length(orignamesp) > 0){
        for(origs in orignamesp){
          cat(".")
          fname <- paste0(prefix, '_all_data_confounders_', qc_tag, sub("orig", "", origs), '.pdf')
          void <- lapply(ivariables, function(x){
            qc_violin(
              mydat = annottab_ss, xax = origs, yax = x,
              lb_filt = low.filt.cells,
              hb_filt = high.filt.cells,
              filtby = apply_filter
            )
          }); sfact <- if(length(table(annottab_ss[, origs])) < 10) 5.5 else 3
          pdf(fname, width = sfact * mygrid[2], height = sfact * mygrid[1])
          if(length(table(annottab_ss[, origs])) < 10){
            print(cowplot::plot_grid(plotlist = void, ncol = mygrid[2]))
          }else{
            tvar <- lapply(void, print)
          }; graphics.off()
        }; cat("\n")
      }
    }
  } # Finish filtering
}
qc_scatters(dat = annottab, thresholds = summ[, apply_filter], prefix = "3_qc")
# ############# ####------------------------------------------------------------

save(annottab, file = "metadata_filtered.rdata")
write.table(rownames(mycellsdata), file = "feature_names.txt", quote = FALSE, row.names = FALSE, col.names = FALSE);
#### ######### ######### ####---------------------------------------------------

cat('\nDONE\nStarting time:\n', st.time, '\n')
cat('Finishing time:\n') ; timestamp()
cat('\n\n*******************************************************************\n')
cat('SESSION INFO:\n') ; print(sessionInfo())
cat('*******************************************************************\n\n')
if(register_log){
  sink(type = 'message')
  sink()
}
