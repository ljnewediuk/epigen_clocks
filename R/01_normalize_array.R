

library(tidyverse)
library(sesame)
library(minfi)

# Install sesame and minfi developmental versions
# Must first install SeSaMe and minfi packages from Bioconductor
# BiocManager::install("minfi")
# BiocManager::install("sesame")
#
# install.packages("SeSaMeRpackages/HorvathMammalMethylChip40manifest_0.2.2.tar.gz", repos=NULL, type="source")
# install.packages("SeSaMeRpackages/HorvathMammalMethylChip40anno.test.unknown_0.2.2.tar.gz", repos=NULL, type="source")


sample_sheet_file_name <- 'input/PB_array1_sample_sheet_rmv_nonrun.csv'

manifest_file_name_sesame <- 'input/HorvathMammal40.CanonicalManifest.3.2019.sesame.csv'

# NORMALIZE USING MINFI

# Sample sheet needs a column 'Basename' pointing to the basename of a two-colour
# IDAT file (i.e., either _Red.idat or _Grn.idat). Need to add the col Basename
# with the file path for each sample in the corresponding row. Then use 'read.metharray.exp'
# to find the corresponding files using the sample sheet.
sample_sheet <- read.csv(sample_sheet_file_name) %>%
  mutate(Basename = paste0('input/iscans/', chip.ID, '/', chip.ID, '_', stripe))

# Get sample IDs and chip ID
sample_IDs <- sample_sheet %>%
  mutate(chip.ID.loc = paste(chip.ID, stripe, sep = '_')) %>%
  select(Sample_Name, chip.ID.loc)

# Creates an RGChannelSet object containing the raw red green channel information from the idat files
RGset <- minfi::read.metharray.exp(base=NULL, targets=sample_sheet, recursive=TRUE) 

# Annotates the RGset object with probe coordinates. This line is currently just a place holder as the annotation is empty, but needs to be run anyway.
RGset@annotation = c(array='HorvathMammalMethylChip40', annotation="test.unknown")

# Calling getBeta on the RGset object will return a data frame of raw beta values for each CG site
# I.e., get the ß values or DNA methylation values at each CG site
# Are the R##C## separate probe sites for each sample?
raw_betas_minfi <- as_tibble(minfi::getBeta(RGset), rownames="CGid")

# Calling preprocessNoob on RGset will return a MethylSet object, which contains normalized beta values along with a few other things
Mset = minfi::preprocessNoob(RGset)

# Calling getBeta on Mset will return a data frame of normalized beta values for each CG Site
normalized_betas_minfi = as_tibble(minfi::getBeta(Mset), rownames="CGid")

# Save normalized betas
# saveRDS(normalized_betas_minfi, 'output/nbetas_PB_array1.rds')

# Load normalized betas
n_betas <- readRDS('output/nbetas_PB_array1.rds')

# Subsetting probes is based on a file containing which probes map to which genomes, 
# created by Mike Thompson. There are currently two versions of this file, one 
# including only probes that map uniquely to each genome, and one that also includes 
# the probes that may map to multiple locations in each genome. These other probes 
# might still provide meaningful signal, but it is potentially an average of multiple 
# loci from the genome. As we don't yet fully understand the behaviors of those probes
# , below I am using the file containing uniquely mapping probes, but you should 
# consider which file makes sense for your own particular analysis.
# 
# Read probe annotation file
probe_mapping_file <- read_csv("input/HorvathMammalChip_SpeciesGenomeCoordsUnique_v1.0.csv.gz", col_types = cols(.default = "c"))


# Obtain a list of probes and coordinates for your species of interest
species_probes <- probe_mapping_file %>%
  # select just the coordinates for your species and the CG id columns
  select(one_of("probeID", 'UrsusMaritimus')) %>% 
  # filter out probes with NA in the coordinates column because those don't map to your species
  filter(!is.na(!!as.name('UrsusMaritimus')))

# Filter probes from normalized betas
n_betas_filtered <- n_betas %>%
  filter(CGid %in% species_probes$probeID)

# Transpose matrix after removing CGid
n_betas_t <- n_betas %>%
  select(! CGid) %>%
  # column_to_rownames(var = CGid)
  t() %>%
  as.data.frame() %>%
  # janitor::row_to_names(1) %>%
  # Make rownames (chip IDs) into column
  rownames_to_column(var = 'chip.ID.loc')

# Rename rest of columns to Cpg sites
colnames(n_betas_t)[2:ncol(n_betas_t)] <- n_betas$CGid

# Add sample info (ALSO NEED TO ADD ACTUAL BEAR AGE, YEAR, SEX, PLUS ANAGE INFO)
n_betas_info <- n_betas_t %>%
  left_join(sample_IDs) %>%
  relocate(Sample_Name, .after = chip.ID.loc)

# Save transposed betas
# USE IN CLOCK SCRIPT TO ESTIMATE AGE
saveRDS(n_betas_info, 'output/beta_info_PB_array1.rds')



