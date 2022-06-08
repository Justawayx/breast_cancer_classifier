############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############
############    --   Data to post in DICE's website    --    ############
############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############


# ---> About the script.
# Version: 0
# Subversion: 1
# Prepare data for DEA for feature engineering.


### -------------------------- Description -------------------------- ###


### --------------------------- Libraries --------------------------- ###
library(data.table)
library(stringr)
library(gtools)


### --------------------------- Functions --------------------------- ###


### ----------------------- General Arguments ----------------------- ###
# ---> Path definitions.
gen.path <- '/mnt/hpcscratch/vfajardo/BISB/BENG-203/course_project'
data.path <- paste0(gen.path, '/data')
reports.path <- data.path
# ---> File definitions.
# Healthy control samples.
hlty.raw.file <- paste0(data.path, '/pnas_normal_readcounts.txt')
hlty.tpm.file <- paste0(data.path, '/pnas_normal_tpm.txt')
# Breast cancer samples.
cancer.meta.file <- paste0(data.path, '/pnas_patient_info.csv')
cancer.raw.file <- paste0(data.path, '/pnas_readcounts_96_nodup.txt')
cancer.tpm.file <- paste0(data.path, '/pnas_tpm_96_nodup.txt')
# ---> Check directories and files.
essential.files <- c(
  hlty.raw.file, hlty.tpm.file,
  cancer.raw.file, cancer.tpm.file
)
essential.files <- essential.files[!unlist(lapply(X=essential.files, FUN=file.exists))]
if(length(essential.files) > 0) stop(paste0('Next essential files are not appropriately defined -make sure you\'ve got their names right-:\n', paste0(essential.files, collapse='\n'), '\n'))


### ------------------------- Data Loading -------------------------- ###

# Metadata.

# Healthy control samples.
hlty.raw.data <- fread(file=hlty.raw.file)
colnames(hlty.raw.data)[1] <- 'geneid'
hlty.tpm.data <- fread(file=hlty.tpm.file)
colnames(hlty.tpm.data)[1] <- 'geneid'
# Breast cancer samples.
cancer.meta <- fread(file=cancer.meta.file)
cancer.raw.data <- fread(file=cancer.raw.file)
colnames(cancer.raw.data) <- c('geneid', cancer.meta[, paste('BC', poiseid, fu, sep='.')])
cancer.tpm.data <- fread(file=cancer.tpm.file)
colnames(cancer.tpm.data) <- c('geneid', cancer.meta[, paste('BC', poiseid, fu, sep='.')])


### ---------------------- Data preprocessing ----------------------- ###
### ------------------------- Main program -------------------------- ###

# ----> Merge data from both types.
# Raw counts.
raw.data <- merge(x=hlty.raw.data, y=cancer.raw.data, by='geneid')
hlty.raw.data[!geneid %in% raw.data[, geneid], .N]
cancer.raw.data[!geneid %in% raw.data[, geneid], .N]
raw.data <- raw.data[str_starts(string=geneid, pattern='ENSG')]
# TPM counts.
tpm.data <- merge(x=hlty.tpm.data, y=cancer.tpm.data, by='geneid')
hlty.tpm.data[!geneid %in% tpm.data[, geneid], .N]
cancer.tpm.data[!geneid %in% tpm.data[, geneid], .N]
tpm.data <- tpm.data[str_starts(string=geneid, pattern='ENSG')]
# Output merged data.
tmp.file.name <- paste0(reports.path, '/MergedRawData.csv')
fwrite(file=tmp.file.name, x=raw.data)
tmp.file.name <- paste0(reports.path, '/MergedTPMData.csv')
fwrite(file=tmp.file.name, x=tpm.data)


# ---> Metadata for all samples.
# Get variables.
sample.id <- c(colnames(hlty.tpm.data)[2:ncol(hlty.tpm.data)], cancer.meta[, paste('BC', poiseid, fu, sep='.')])
disease.status <- c(
  rep(x='control', times=ncol(hlty.tpm.data)-1),
  rep(x='cancer', times=nrow(cancer.meta))
)
donor.id <- c(colnames(hlty.tpm.data)[2:ncol(hlty.tpm.data)], cancer.meta[, poiseid])
all.meta <- data.table(
  sample.id, disease.status, donor.id
)
# Output whole metadata.
tmp.file.name <- paste0(reports.path, '/MergedMetadata.csv')
fwrite(file=tmp.file.name, x=all.meta)


# ---> Cancer samples data for comparison between recurrence and non-recurrence.
# Raw counts.
cancer.raw.data <- cancer.raw.data[str_starts(string=geneid, pattern='ENSG')]
tmp.file.name <- paste0(reports.path, '/CancerRawData.csv')
fwrite(file=tmp.file.name, x=cancer.raw.data)
# TPM counts.
cancer.tpm.data <- cancer.tpm.data[str_starts(string=geneid, pattern='ENSG')]
tmp.file.name <- paste0(reports.path, '/CancerTPMData.csv')
fwrite(file=tmp.file.name, x=cancer.raw.data)


# ---> Metadata for cancer samples.
# Get variables.
sample.id <- cancer.meta[, paste('BC', poiseid, fu, sep='.')]
recurr.status <- cancer.meta[, recurStatus]
disease.status <- c(
  rep(x='cancer', times=nrow(cancer.meta))
)
donor.id <- c(cancer.meta[, poiseid])
cancer.meta <- data.table(
  sample.id, recurr.status, disease.status, donor.id
)
# Output whole metadata.
tmp.file.name <- paste0(reports.path, '/CancerMetadata.csv')
fwrite(file=tmp.file.name, x=cancer.meta)
