
library(tidyverse)
library(glmnet)

# Sample age information
samp_age <- read.csv('batches/batch1_samples.csv') %>%
  select(sampleId, age, sex, Spec) %>%
  dplyr::rename('SampleID' = sampleId, 'Age' = age, 'Sex' = sex)

# This dataset includes species characters (needed for clocks) and CpGs 
info <- readRDS('output/beta_info_PB_array1.rds') %>%
  dplyr::rename('SampleID' = Sample_Name) %>%
  mutate(GestationTimeInYears = 65/365,
         averagedMaturity.yrs = 1734/365,
         maxAgeCaesar = 43.8) %>%
  left_join(samp_age) %>%
  relocate(c(GestationTimeInYears:maxAgeCaesar, Age, Sex, Spec), .after = SampleID)

# S3.1 to S3.3 tables in SupplmentaryData of Lu et al.
glmnet.csv <- c('input/clock1.csv', 'input/clock2.csv', 'input/clock3.csv')

# Variables
beta.name <- c('beta_clock1','beta_clock2','beta_clock3')
y.name <- c('Y.pred1','Y.pred2','Y.pred3') 
age.name <- c('DNAmAgeClock1','DNAmAgeClock2','DNAmAgeClock3')

# Clock 2 function
F2_antitrans_clock2 <- function(y, y.maxAge, y.gestation, const = 1) {
  x0 <- const * exp(- exp(- 1 * y))
  x1 <- x0 * (y.maxAge + y.gestation)
  x <- x1 - y.gestation
  x
}

#clock3
F1_logli <- function(age1, m1, m2 = m1, c1 = 1) {
  ifelse(age1 >= m1, (age1 - m1) / m2 , c1 * log((age1 - m1) / m2 / c1 + 1)) 
}

#RelativeAdultAge
F2_revtrsf_clock3 <- function(y.pred, m1, m2 = m1, c1 = 1) {
  ifelse(y.pred < 0, (exp(y.pred / c1) - 1) * m2 * c1 + m1, y.pred * m2 + m1 ) 
}


# The `loglifn` function shows how to calculate m1 for the transformation 
# It is the `a_Logli` in the function
F3_loglifn <- function(dat1, b1 = 1, max_tage = 4, c1 = 5, c2 = 0.38, c0  = 0) {
  n <- nrow(dat1)
  age1 <- (dat1$maxAgeCaesar+dat1$GestationTimeInYears) / (dat1$averagedMaturity.yrs + dat1$GestationTimeInYears)
  a1 <- age1 / (1 + max_tage)
  dat1$a1_Logli <- a1 #x/m1 in manuscript
  a2 = (dat1$GestationTimeInYears + c0) / (dat1$averagedMaturity.yrs) 
  #m=5*(G/ASM)^0.38 from regression analysis/formula(7)
  dat1$a_Logli = a_Logli = c1*a2^c2
  x <- dat1$Age + dat1$GestationTimeInYears
  t2 = dat1$averagedMaturity.yrs * b1 + dat1$GestationTimeInYears 
  x2 = x/t2 #### log(x/t2)
  y = F1_logli(x2, a_Logli, a_Logli)
  dat1$LogliAge <- y
  return(dat1)
}

#
info$Intercept <- 1

MYMAX=1.3

info$HighmaxAgeCaesar = MYMAX * info$maxAgeCaesar 
info$HighmaxAgeCaesar[info$SpeciesLatinName == 'Ursus maritimus'] <- info$maxAgeCaesar[info$SpeciesLatinName == 'Ursus maritimus']

# predict RelativeAge

#
for(k in 1:3){
  glmnet <- read.csv(glmnet.csv[k])
  glmnet$beta <- glmnet[, beta.name[k]] 
  glmnet$var[1] <- ifelse(glmnet$var[1] == "(Intercept)", 'Intercept', glmnet$var[1]) 
  temp <- info %>%
    select(as.character(glmnet$var)) %>%
    as.matrix()
  # temp <- as.matrix(subset(info, select=as.character(glmnet$var))) 
  info[,y.name[k]]=as.numeric(as.matrix(subset(info,select=as.character(glmnet$var))) %*% glmnet$beta)
}

#(1) Clock 1
info[,age.name[1]] <- exp(info[, y.name[k]]) - 2
#(2) Clock 2
info$DNAmRelativeAge <- exp(- exp(- 1 * info[, y.name[2]]))
info[, age.name[2]] <- F2_antitrans_clock2(info[, y.name[2]], info$HighmaxAgeCaesar, info$GestationTimeInYears, const = 1) #(3) Clock 3
info <- F3_loglifn(info)#to compute m estimate for tuning point in the log-linear transformation 
info$m1 = info$a_Logli
# This option kind of makes it more accurate
# info$m1 = info$LogliAge
# Here 
info$DNAmRelativeAdultAge <- F2_revtrsf_clock3(info[, y.name[3]], info$m1)
info[, age.name[3]]<- info$DNAmRelativeAdultAge * (info$averagedMaturity.yrs + info$GestationTimeInYears) - info$GestationTimeInYears
#
#final output
output <- subset(info,select=c('SampleID','Age', 'Sex', 'Spec', 'DNAmRelativeAge','DNAmRelativeAdultAge',age.name))

# Calculate age error
age_error <- output %>%
  mutate(C1 = DNAmAgeClock1 - Age,
         C2 = DNAmAgeClock2 - Age,
         C3 = DNAmAgeClock3 - Age) %>%
  mutate(ID = substr(SampleID, 1, 6),
         Year = as.numeric(substr(SampleID, 8, 11))) %>%
  select(ID, Year, Age, Sex, Spec, C1:C3) %>%
  pivot_longer(C1:C3, names_to = 'Clock', values_to = 'Error') %>%
  mutate(abs_Error = abs(Error))

ggplot(age_error, aes(x = Year, y = Error, colour = Clock)) +
  geom_point() 

ggplot(age_error, aes(x = Age, y = Error, colour = Clock)) +
  geom_point() 

ggplot(age_error, aes(x = Sex, y = Error)) + geom_boxplot()

ggplot(age_error, aes(x = Spec, y = Error)) + geom_boxplot()

