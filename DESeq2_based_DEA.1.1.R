cat('\n\n')
############    --   Differential Expression Analysis   --    ############
############    -----------  Based on DESeq2   -----------    ############
cat('############    --   Differential Expression Analysis   --    ############\n')
cat('############    -----------  Based on DESeq2   -----------    ############\n')

# By: Vicente Fajardo

# Version: 1
# Version updates:
#   First version.
# Subversion: 0
# Subversion updates:
#   First subversion.

### -------------------------- Description -------------------------- ###

# Given a dataset with two or more condition groups, this program will apply Differential (Expression) Analysis (DE/DEA) among all different groups (a formula for this comparisons must be given) and will create a summary of the multiple pairwise comparisons. Along the way, the program will also output the results for each specific comparison. Appropriate QCs are also output.
# NOTICE: Please check the QC results and the model fit output look good previous to make any conclusions out of the DEA results. Check any further documentation on their interpretation in this Google docs file: https://docs.google.com/document/d/1G7udcsKj_4048sYEkS7ObapbPMzgLcy0vlFhC6omQp0/edit
# If you have no access to the documentation, please just ask for it and permision should be granted soonly upon your request.

cat('\n\n')
### -------------------------- Dependencies ------------------------- ###
cat('### -------------------------- Dependencies ------------------------- ###\n')
library(optparse)
library(DESeq2)
library(stringr)
library(NMF)
library(ggplot2)
library(data.table)
source('/home/vfajardo/scripts/functions/R_handy_functions.0.4.R')
cat('Dependencies imported!\n')

cat('\n\n')
### --------------------------- Arguments --------------------------- ###
cat('### --------------------------- Arguments --------------------------- ###\n')
option.list <- list(
  make_option(opt_str="--ReportsPath", type="character", default=NULL, dest="reports.path", help="Absolute path to directory to save DE genes table.\n"),
  make_option(opt_str="--CountData", type="character", default=NULL, dest="counts.data.file", help="Absolute path to file saving counts data. Their column names (different than the first one, which'll turn into the row names) must match the row names of the metadata file.\n"),
  make_option(opt_str="--MetaData", type="character", default=NULL, dest="meta.data.file", help="Absolute path to file storing metadata. Their row names (defined in the first column of the file regardless of the name it has) must match the column names of the counts data file.\n"),
  make_option(opt_str="--ImpSamples", type="character", default=NULL, dest="samples.file", help="Absolute path to csv file (single column) listing the samples of interest to consider during the analysis.\n"),
  make_option(opt_str="--Design", type="character", default=NULL, dest="design.formula", help="Character, design formula to create the DESeq2 object. See DESeq2 documentation for details, specifically for the function DESeqDataSetFromMatrix.\n"),
  make_option(opt_str="--Contrast", type="character", default=NULL, dest="contrast.element", help="Character, element in the design formula to use for contrast for pairwise comparisons. If NULL, the program will take the last element in the formula assuming it's a single one (i.e., no ':' in the last element).\n"),
  make_option(opt_str="--LFCShrink", type="logical", default=TRUE, dest="lfc.shrink", help="Logical, indicates whether LFC shrinkage should be applied per pairwise comparison.\n"),
  # ---> Significance thresholds.
  make_option(opt_str="--PThold", type="numeric", default=0.05, dest="p.thold", help="Numeric, p-adjusted value trheshold to call signifcant DEGs.\n"),
  make_option(opt_str="--LFCThold", type="numeric", default=1, dest="lfc.thold", help="Numeric, LFC value trheshold to call signifcant DEGs (used for both, upregulated genes (> threshold) and downregulated genes (< minus threshold)).\n"),
  # ---> Not required.
  make_option(opt_str="--TPMData", type="character", default=NULL, dest="tpms.data.file", help="Absolute path to file saving normalized counts data. Their column names (different than the first one, which'll turn into the row names) must match the row names of the metadata file. If NULL, this option will be disregarded.\n"),
  make_option(opt_str="--FeatAnns", type="character", default=NULL, dest="feat.anns.file", help="Absolute path to file saving feature annotations. Relevant annotations (up to this program version) are 'gene_type' and 'gene_name'; therefore, any other annotations in the file will be disregarded. If NULL, this option will be disregarded.\n")
)
# Getting arguments from command line and setting their values to their respective variable names.
opt.parser = OptionParser(option_list=option.list);
opt = parse_args(opt.parser);
# Moving options to their own variables
reports.path <- opt$reports.path
counts.data.file <- opt$counts.data.file
meta.data.file <- opt$meta.data.file
samples.file <- opt$samples.file
design.formula <- opt$design.formula
contrast.element <- opt$contrast.element
lfc.shrink <- opt$lfc.shrink
p.thold <- opt$p.thold
lfc.thold <- opt$lfc.thold
tpms.data.file <- opt$tpms.data.file
feat.anns.file <- opt$feat.anns.file
# Temp.
# reports.path <- '/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/tmp/deseq2'
# meta.data.file <- '/mnt/BioAdHoc/Groups/vd-vijay/Cristian/RNA-seq/projects/10X_SMARTSeq2/metadata_all.csv'
# counts.data.file <- '/mnt/BioAdHoc/Groups/vd-vijay/Cristian/RNA-seq/projects/10X_SMARTSeq2/4.Output/counts/raw_counts.csv'
# samples.file <- '/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/tmp/deseq2/input/stim_samples.csv'
# design.formula <- '~ Cell.Type'
# contrast.element <- NULL
# lfc.shrink <- TRUE
# tpms.data.file <- '/mnt/BioAdHoc/Groups/vd-vijay/Cristian/RNA-seq/projects/10X_SMARTSeq2/4.Output/counts/TPM_counts.csv'
# feat.anns.file <- '/mnt/BioAdHoc/Groups/vd-ay/RNASeq_Workflow/Reference/GRCh37_annotation.csv'
# ---> Show argument values.
cat('Reports path: ', reports.path, '\n')
cat('Counts data file: ', counts.data.file, '\n')
cat('Metadata file: ', meta.data.file, '\n')
cat('Samples of interest file (if any): ', samples.file, '\n')
cat('Design formula: ', design.formula, '\n')
cat('Contrast element (none if NULL): ', contrast.element, '\n')
cat('Should LFC shrinkage be applied? ', lfc.shrink, '\n')
cat('TPMs data file: ', tpms.data.file, '\n')
cat('Features annotations file: ', feat.anns.file, '\n')
cat('\n\n')

# Defaults ------------------------------------------------------------->
# Others --------------------------------------------------------------->
# Column names for the DEA results.
col.replace <- c('base.mean', 'lfc', 'lfc.se', 'test.stat', 'p.val', 'p.adj')
names(col.replace) <- c('baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj')

### --------------------------- Functions --------------------------- ###

cat('\n\n')
### --------------------------- Load data --------------------------- ###
cat('### --------------------------- Load data --------------------------- ###\n')
# ---> Counts data.
counts.data <- read.csv(file=counts.data.file, stringsAsFactors=FALSE, row.names=1, check.names=FALSE)
# ---> Metadata.
meta.data <- read.csv(file=meta.data.file, row.names=1)
# ---> Samples.
# Samples to keep (if specified).
if(!is.null(samples.file)){
  samples.to.keep <- read.csv(file=samples.file, stringsAsFactors=FALSE)[, 1]
  if(!(all(samples.to.keep %in% colnames(counts.data)) & all(samples.to.keep %in% rownames(meta.data)))) stop('Samples to take into account were indicated in samples file (path below), but they were not found as appropriately defined in either metadata, counts data or both.\n\tSamples file: ', samples.file, '\n\n')
  counts.data <- counts.data[, samples.to.keep]
  meta.data <- meta.data[samples.to.keep, ]
}else{
  samples.to.keep <- colnames(counts.data)
}
# ---> Check.
# @ Metadata samples match those from counts data.
if(all(colnames(counts.data) %in% rownames(meta.data)) & all(rownames(meta.data) %in% colnames(counts.data))){
  # Check one matches the other in the same order. If not, order accordingly.
  if(!all(colnames(counts.data) == rownames(meta.data))){
    new.idxs <- match(x=colnames(counts.data), table=rownames(meta.data))
    meta.data <- meta.data[new.idxs, ]
    cat('Metadata and counts data samples were modified to match accordingly.\n')
  }else{
    cat('Metadata and counts data samples matched accordingly.\n')
  }
}else{
  stop('Samples defined for counts data and metadata do not match. For both files, first column, regardless its name, will be deinfed as the row names for each table. Check this makes sense.\n')
}
# @ Design elements are part of the metadata columns.
design.elements <- str_replace(string=design.formula, pattern='~\\s*', replacement='')
design.elements <- str_split(string=design.elements, pattern='\\s*\\+\\s*')
design.elements <- unique(unlist(lapply(X=design.elements, FUN=function(tmp.string) str_split(string=tmp.string, pattern='\\:'))))
if(!all(design.elements %in% colnames(meta.data))) stop('Not all elements defining the design formula input are part of the samples\' metadata.\n')
# ---> TPMs data.
if(!is.null(tpms.data.file)){
  tpms.data <- read.csv(file=tpms.data.file, stringsAsFactors=FALSE, row.names=1, check.names=FALSE)
  if(all(samples.to.keep %in% colnames(tpms.data))) tpms.data <- tpms.data[, samples.to.keep] else stop('TPMs file was indicated (path below), but all important samples were not found as appropriately defined.\n\tTPMs file: ', tpms.data.file, '\n\n')
}else{
  tpms.data <- NULL
}
# ---> Feature annotations.
if(!is.null(feat.anns.file)){
  feat.anns <- read.csv(file=feat.anns.file, stringsAsFactors=FALSE, row.names=1)
  feat.anns$gene.id <- rownames(feat.anns)
  feat.anns <- as.data.table(feat.anns)
  # Keep only relevant info.
  if(all(c('gene_name', 'gene_type') %in% colnames(feat.anns))){
    feat.anns <- feat.anns[, c('gene.id', 'gene_name', 'gene_type')]
  }else{
    feat.anns <- NULL
  }
}else{
  feat.anns <- NULL
}
cat('All data appropriately loaded!\n')

cat('\n\n')
### ------------------------- Predefinitions ------------------------ ###
cat('### ------------------------- Predefinitions ------------------------ ###')
# ---> DESeq2 object
deseq.obj <- DESeqDataSetFromMatrix(countData=counts.data, colData=meta.data, design=as.formula(design.formula))
# ---> Normalization  by DESeq2 normalization.
# Estimate size factors per library.
deseq.obj <- estimateSizeFactors(deseq.obj)
cat('DESeq2 object created succesfully!\n')
# ---> Contrast element.
# Determine if necessary.
if(is.null(contrast.element)){
  contrast.element <- str_replace(string=design.formula, pattern='~\\s*', replacement='')
  contrast.element <- str_split(string=contrast.element, pattern='\\s*\\+\\s*')
  contrast.element <- contrast.element[[length(contrast.element)]]
  cat(paste0('Contrast element determined here: ', contrast.element), '\n')
}else{
  if(!contrast.element %in% design.elements) stop('Contrast element input not part of design formula.\n')
}

cat('\n\n')
### ------------------------- Main program -------------------------- ###
cat('### ------------------------- Main program -------------------------- ###\n')

cat('\n\n')
### ------------------------ Quality control ------------------------ ###
cat('### ------------------------ Quality control ------------------------ ###\n')

qc.path <- paste0(reports.path, '/QC')
create.dir(qc.path, 'QC')

# ---> Log-transformed, normaized data by VST.
vst.mat <- DESeq2::vst(object=deseq.obj, blind=TRUE)

# ---> Colors definition.
elements.cols <- lapply(X=design.elements, FUN=function(tmp.element){
  element.vals <- unique(as.character(meta.data[, tmp.element]))
  to.output <- sample(x=favorite.colors, size=length(element.vals), replace=FALSE)
  names(to.output) <- element.vals
  return(to.output)
})
names(elements.cols) <- design.elements

# ---> Heatmap.
# Get correlations.
corr.mat <- cor(assay(vst.mat))
# Prepare MetaData
tmp.meta.data <- meta.data[, design.elements]
if(is.null(dim(tmp.meta.data))) tmp.meta.data <- data.frame(tmp.meta.data); colnames(tmp.meta.data) <- design.elements
# Heatmap and output.
tmp.file.name <- paste0(qc.path, '/CorrelationsHeatmap.pdf')
aheatmap(x=corr.mat, scale='none', Rowv=TRUE, Colv=TRUE, annCol=tmp.meta.data, annColors=elements.cols, filename=tmp.file.name, width=10)
# Heatmap with all metadata elements.
tmp.file.name <- paste0(qc.path, '/CorrelationsHeatmapAllColumns.pdf')
aheatmap(x=corr.mat, scale='none', Rowv=TRUE, Colv=TRUE, annCol=meta.data, annColors=elements.cols, filename=tmp.file.name, width=10)

# ---> PCA
tmp.file.name <- paste0(qc.path, '/PCAExplorations.pdf')
pdf(file=tmp.file.name)
for(tmp.element in design.elements){
  tmp.ggplot <- DESeq2::plotPCA(object=vst.mat, intgroup=tmp.element)
  tmp.ggplot <- tmp.ggplot + scale_color_manual(values=elements.cols[[tmp.element]])
  print(tmp.ggplot + theme_minimal())
}
dev.off()

cat(paste0('Hierarchical clustering and PCA applied considering the elements that are part of the design formula. Please check results make sense.\n'))

cat('\n\n')
### --------------------- Differential analysis --------------------- ###
cat('### --------------------- Differential analysis --------------------- ###\n')

model.path <- paste0(reports.path, '/model_fitting')
create.dir(model.path, 'Model fitting')

# ---> DEA
deseq.obj <- DESeq(object=deseq.obj)

# ---> Dispersion estimates.
tmp.file.name <- paste0(model.path, '/DispersionEstimatesModel.pdf')
pdf(file=tmp.file.name)
print(plotDispEsts(object=deseq.obj))
dev.off()

cat('NB Model fitted to data (check dispersion estimates fit to the model appropriately) and DEA applied.\n')

cat('\n\n')
### ----------------- Multiple pairwise comparisons ----------------- ###
cat('### ----------------- Multiple pairwise comparisons ----------------- ###\n')

dea.path <- paste0(reports.path, '/differential_analysis')
create.dir(dea.path, 'DEA')

# ---> Multiple pairwise comparisons.

# Get values for the contrast element.
contrast.vals <- unique(as.character(meta.data[, contrast.element]))

# Stats of normalized values across groups.
tpm.stats <- lapply(X=contrast.vals, FUN=function(tmp.val){
  tmp.samples <- rownames(meta.data)[meta.data[, contrast.element]==tmp.val]
  tmp.output <- rowMeans(tpms.data[, tmp.samples])
  tmp.output <- data.table(tmp.output)
  colnames(tmp.output) <- paste0('mean.tpm.', tmp.val)
  return(tmp.output)
})
tpm.stats <- Reduce(x=tpm.stats, f=cbind)
tpm.stats <- cbind(data.table(gene.id=rownames(tpms.data)), tpm.stats)

# ---> DEA per level to compare.

all.dea.results <- lapply(X=contrast.vals, FUN=function(lvl.to.comp){

  lvl.dea.path <- paste0(dea.path, '/', lvl.to.comp)
  create.dir(lvl.dea.path, paste0('DEA for level to compare: ', lvl.to.comp))

  # ---> Individual comparisons.
  ind.dea.path <- paste0(lvl.dea.path, '/individual_comparisons')
  create.dir(ind.dea.path, 'Individual comparisons')

  # Against the rest of the values for multiple pairwise comparisons.
  base.lvls <- setdiff(x=contrast.vals, y=lvl.to.comp)
  dea.results <- lapply(X=base.lvls, FUN=function(base.lvl){
    # Directory for specific comparison.
    base.dea.path <- paste0(ind.dea.path, '/', base.lvl)
    create.dir(base.dea.path, paste0('Specific comparison ', lvl.to.comp, ' vs ', base.lvl))
    # Create contrast.
    contrast.vec <- c(contrast.element, lvl.to.comp, base.lvl)
    # Retrieve results.
    to.output <- results(object=deseq.obj, contrast=contrast.vec, alpha=0.05)
    # Means vs LFCs (MA plot).
    tmp.file.name <- paste0(base.dea.path, '/NormalizedCountsMean_vs_LFCs_NoShrinkage.pdf')
    pdf(file=tmp.file.name, width=10)
    print(plotMA(object=to.output))
    dev.off()
    # LFC shrinkage (if required)
    if(lfc.shrink){
      # Shrinkage
      to.output <- lfcShrink(dds=deseq.obj, contrast=contrast.vec, res=to.output)
      # Means vs LFCs (MA plot).
      tmp.file.name <- paste0(base.dea.path, '/NormalizedCountsMean_vs_LFCs_Shrinkage.pdf')
      pdf(file=tmp.file.name, width=10)
      print(plotMA(object=to.output))
      dev.off()
    }
    # Results summary.
    colnames(to.output) <- col.replace[colnames(to.output)]
    gene.ids <- data.table(gene.id=row.names(to.output))
    to.output <- as.data.table(to.output)
    to.output <- cbind(gene.ids, to.output)
    to.output[, is.significant:=ifelse(test=p.adj<0.05, yes='significant', no='non-significant')]
    to.output[, lfc.type:=ifelse(test=lfc>0, yes='up', no='down')]
    res.summ <- to.output[, .(is.significant, lfc.type)][, .(genes=.N), by=.(is.significant, lfc.type)]
    setorderv(x=res.summ, cols=c('is.significant', 'lfc.type'), order=1)
    tmp.file.name <- paste0(base.dea.path, '/ResultsSummary.csv')
    fwrite(file=tmp.file.name, x=res.summ, na=NA, quote=FALSE)
    # Volcano plot.
    tmp.cols <- c('firebrick', 'black')
    names(tmp.cols) <- c('significant', 'non-significant')
    tmp.ggplot <- ggplot(to.output[!is.na(is.significant)], aes(x=lfc, y=-log10(p.adj), col=is.significant)) + geom_point(alpha=0.8, size=0.7) + scale_color_manual(values=tmp.cols) + labs(x='Log(2) Fold Change', y='Significance (-log10(p-adjusted value))') + theme_minimal()
    tmp.file.name <- paste0(base.dea.path, '/VolcanoPlot.pdf')
    pdf(file=tmp.file.name)
    print(tmp.ggplot)
    dev.off()
    # Volcano plot with top 20 genes.
    setorderv(x=to.output, cols='p.adj', order=1, na.last=TRUE)
    top.genes <- to.output[1:20]
    if(!is.null(feat.anns)){
      top.genes <- merge(x=top.genes, y=feat.anns, by='gene.id', all.x=TRUE, all.y=FALSE)
    }else{
      top.genes$gene_name <- top.genes$gene.id
    }
    tmp.ggplot <- tmp.ggplot + geom_text(data=top.genes, aes(label=gene_name), color='black', show.legend=FALSE, check_overlap=TRUE, size=2.5)
    tmp.file.name <- paste0(base.dea.path, '/VolcanoPlotWithLabels.pdf')
    pdf(file=tmp.file.name)
    print(tmp.ggplot)
    dev.off()
    # Remove non-significant genes and non-important columns.
    to.output[, lfc.type:=NULL]
    to.output[, is.significant:=NULL]
    to.output <- to.output[p.adj<p.thold & (lfc>lfc.thold | lfc<(-lfc.thold))]
    # Bind important stats and metadata.
    tmp.data <- to.output
    tmp.data <- as.data.frame(tmp.data)
    # Stats
    if(!is.null(tpm.stats)){
      tmp.data <- merge(x=tmp.data, y=tpm.stats, by='gene.id', all.x=TRUE, all.y=FALSE)
    }
    # Metadata
    if(!is.null(feat.anns)){
      tmp.data <- merge(x=feat.anns, y=tmp.data, by='gene.id', all.x=FALSE, all.y=TRUE)
    }
    # Output final table saving to a file and out of the function...
    tmp.file.name <- paste0(base.dea.path, '/DEGsInfo.csv')
    fwrite(x=tmp.data, file=tmp.file.name, quote=FALSE, na= NA)
    # ... disregarding redundant info.
    to.output <- to.output[, c('gene.id', 'lfc', 'lfc.se', 'test.stat', 'p.val', 'p.adj')]
    colnames(to.output) <- c('gene.id', paste(colnames(to.output)[-1], base.lvl, sep='.'))
    return(to.output)
  })

  # ---> Merge pairwise comparisons.
  merged.dea.path <- paste0(lvl.dea.path, '/merged_results')
  create.dir(merged.dea.path, 'Merged results')

  # Merge results.
  dea.results <- Reduce(x=dea.results, f=function(x, y){
    merge(x=x, y=y, by='gene.id', all=FALSE)
  })

  # Calculate minimum LFC and maximum p-value across comparisons.
  # Minimum lFC.
  lfc.cols <- colnames(dea.results)[str_detect(string=colnames(dea.results), pattern='lfc\\.[^se\\.]')]
  dea.results$lfc..minimum <- rowMins(as.matrix(dea.results[, ..lfc.cols]))
  # Maximum p-value
  pval.cols <- colnames(dea.results)[str_detect(string=colnames(dea.results), pattern='p.adj')]
  dea.results$p.adj..maximum <- rowMax(as.matrix(dea.results[, ..pval.cols]))
  # Order columns.
  cols.order <- sort(colnames(dea.results))
  dea.results <- dea.results[, ..cols.order]

  # Bind important stats and metadata.
  # Stats
  if(!is.null(tpm.stats)){
    dea.results <- merge(x=dea.results, y=tpm.stats, by='gene.id', all.x=TRUE, all.y=FALSE)
  }
  # Metadata
  if(!is.null(feat.anns)){
    dea.results <- merge(x=feat.anns, y=dea.results, by='gene.id', all.x=FALSE, all.y=TRUE)
  }

  # Output
  tmp.file.name <- paste0(merged.dea.path, '/MergedResults_DEGsInfo.csv')
  fwrite(file=tmp.file.name, x=dea.results, quote=FALSE, na=NA)

  return(dea.results)
})
names(all.dea.results) <- contrast.vals

# ---> Merge DEA results for all levels.
all.dea.path <- paste0(dea.path, '/1_multiple_comparisons')
create.dir(all.dea.path, 'Multiple comparisons')

# Merge results.
all.dea.results <- rbindlist(l=all.dea.results, use.names=TRUE, fill=TRUE, idcol=contrast.element)

# Order column names.
cols.to.order <- colnames(all.dea.results)
cols.to.exclude <- c(contrast.element, 'gene.id', 'gene_name', 'gene_type', cols.to.order[str_detect(string=cols.to.order, pattern='mean.tpm')])
cols.to.order <- setdiff(x=cols.to.order, y=cols.to.exclude)
final.cols.order <- c(cols.to.exclude, sort(cols.to.order))
all.dea.results <- all.dea.results[, ..final.cols.order]

# Output final result.
tmp.file.name <- paste0(all.dea.path, '/MergedResults_DEGsInfo_MultipleComparisons.csv')
fwrite(file=tmp.file.name, x=all.dea.results, quote=FALSE, na=NA)

cat('Multiple pairwise comparisons carried out. Table summary output as: ', tmp.file.name, '\n')

sessionInfo()

cat('\n\nAll finished!\n\n')
