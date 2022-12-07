
library(tidyverse)
library(glmnet)

# Load sample specs (i.e., sample type)
samp_specs <- read.csv('input/batch1_samples.csv') %>%
  select(sampleId, Spec, age, sex) %>%
  dplyr::rename('SampleID' = sampleId, 'Age' = age, 'Sex' = sex)

# Load capture information (to get ages and how they were estimated)
samp_metrics <- read.csv('input/bear_capture_info.csv') %>%
  mutate(BC_yr = paste(BearCode, substr(Date, 1, 4), sep = '_')) %>%
  select(BC_yr, AQual) %>% 
  distinct()

# This dataset includes species characters (needed for clocks) and CpGs 
meth_betas <- readRDS('output/beta_info_PB_array1.rds') %>%
  dplyr::rename('SampleID' = Sample_Name) %>%
  left_join(samp_specs) %>%
  relocate(c(Age, Sex, Spec), .after = SampleID)

# Get matrix of betas 
# ****** NEED TO CONVERT TO MATRIX WITHOUT EVER ADDING EXTRA COLS
meth_betas_m <- meth_betas %>%
  # select(! chip.ID.loc:Sample_Name) %>%
  # as.matrix() %>%
  # Make chip positions rownames
  column_to_rownames('chip.ID.loc') %>%
  # # Remove extra cols
  select(! SampleID:Spec) %>%
  # # Convert to matrix; problems with NA introduced by coersion, and makeX
  # # function
  as.matrix()

# Randomize order of samples
# order_betas <- sample(1:nrow(meth_betas), nrow(meth_betas), replace = FALSE)
# meth_betas <- meth_betas[order_betas, ]
# meth_betas_m <- meth_betas_m[order_betas, ]

# Add column for predictions
age_preds <- meth_betas %>%
  mutate(AgePredict = 0) %>%
  select(SampleID:Spec, AgePredict)

# Sample cutoffs
sample_cuts <- seq(10, 70, 10)

for (i in sample_cuts) {
  ifelse(i == 70, 
         betasLoop <- meth_betas_m[-c(71:79), ], 
         betasLoop <- meth_betas_m[-c((i - 9):i), ])

  ifelse(i == 70, 
         ageLoop <- as.numeric(age_preds[-c(71:79), ]$Age), 
         ageLoop <-  as.numeric(age_preds[-c((i - 9):i), ]$Age))
  
  # Cross validation using x = betasLoop as matrix of predictors (Cpgs) and
  # y = agesLoop as vector of responses (chronological ages)
  cvfit <- cv.glmnet(betasLoop, ageLoop, nfolds = 10, alpha = .5)
 
  # Make predictions
  if (i %in% seq(10, 60, 10)) {
    age_preds$AgePredict[c((i - 9):i)] <- predict(cvfit, newx = meth_betas_m[c((i - 9):i), ], type = "response", s = "lambda.min")
  }
  if (i == 70) {
    age_preds[c(71:79),]$AgePredict <- predict(cvfit, newx = meth_betas_m[c(71:79),], type = "link", s = "lambda.min")
  }
}






