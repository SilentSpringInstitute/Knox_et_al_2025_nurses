#Application of a non-targeted biomonitoring method to characterize occupational chemical exposures of women nurses relative to office workers 
#Kristin E. Knox, Dimitri Abrahamsson, Jessica Trowbridge, June-Soo Park, Miaomiao Wang, Erin Carrera, Lisa Hartmayer, Rachel Morello-Frosch, R. A. Rudel

# Additional processing/analysis of MS-DIAL data files, ESI-
# AUTHOR: K Knox
# WRITTEN IN R VERSION: R version 4.4.0 (2024-04-24)


library(tidyverse)
library(readxl)
library(truncnorm)
library(magrittr)
library(sva)
library(FactoMineR)
library(factoextra)
library(broom)
library(rstatix)
library(janitor)
library(patchwork)


# Set the working directory
workingdir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(workingdir)
setwd('..') 


options(scipen=999)

##########

#map between id and batch
map <-read_xlsx("data/sample_batch_map.xlsx", sheet = "clean") %>%
  unique() %>% #only 105; we are missing a batch for W131
  rename(participant_id=sample)

# MS Dial files (2/22/24 run, converted to .csv)
negative <- read_csv("data/Feb2024_files/MSDial_resultsNurse_Neg_20242231353.csv")


# Just keep the variables we need
negative_small <- negative %>%
  select(-`Metabolite name`, -`MS/MS assigned`, -`Reference RT`, -`Reference m/z`, -Formula, -Ontology, 
         -INCHIKEY, -SMILES, -`Annotation tag (VS1.0)`, -`RT matched`, -`m/z matched`, -`MS/MS matched`, 
         -Comment, -`Manually modified for quantification`, -`Manually modified for annotation`, 
         -`Isotope tracking parent ID`, -`Isotope tracking weight number`, -`Total score`, -`RT similarity`,
         -`Dot product`, -`Reverse dot product`, -`Fragment presence %`, -`MS/MS spectrum`)

# take these out because these are just water blanks using different tips (gel tips) to test the background 
negative_small <- negative_small %>%
  select(-`DIWB-Gel1`, -`DIWB-Gel2`, -`DIWB-Gel3`, -`DIWB-Norm`)

# NoINJ=no-injections - injected before each run to equilibrate the system - disregard
negative_small <- negative_small %>%
  select(-`NoINJ-r001`, -`NoINJ-r002`, -`NoINJ-r003`, -`NoINJ-r004`, -`NoINJ-r005`)

# also drop the pooled samples
negative_small <- negative_small %>%
  select(-`OW19-r1`, -`OW19-r2`, -pool_N19_1, -pool_N19_2, -`poolN19-r1`, -`poolN19-r2`, -`pool_N9-r1`, -`pool_N9-r2`, -pool_OW19_1, -pool_OW19_2, -`poolOW19-r1`, -`poolOW19-r2`)
#Nurses: 133:136, 139:142, 144:146, 150:160, 162, 165:168, 174:182, 184:187, 189:193, 195:196, 198:204, 207:215 (63 people)
#OW: 101:132, 137:138, 147, 163:164, 170:173, 183, 194 (43 people)
#Blanks: DIWB1_  negB1, negB2, negB3, negB4, negB5
#        DIWB2_  negB1, negB2, negb3, negB4, negB5
# Also have MB1, MB1_B1, MB2, MB2_B2 and QC1_2, QC1-1, QC1-2, QC1-2-r1, QC2-1, QC2-2, QC2-2-r2 

# Note: Data has field blanks (they look like people from the IDs)
# Field blanks are: 0147, N150, O163, N165,O194, N196
# name the field blanks as such
#"W147", "N150", "W163", "N165","W194", "N196"
# We will designate the office worker field blanks FB1, FB2, FB3, and 
# the nurse field blanks FB4, FB5, FB6
negative_small <- negative_small %>%
  rename(`FB1-r1` = `W147-r1`,
         `FB1-r2` = `W147-r2`,
         `FB2-r1` = `W163-r1`,
         `FB2-r2` = `W163-r2`,
         `FB3-r1` = `W194-r1`,
         `FB3-r2` = `W194-r2`,
         `FB4-r1` = `N150-r1`,
         `FB4-r2` = `N150-r2`,
         `FB5-r1` = `N165-r1`,
         `FB5-r2` = `N165-r2`,
         `FB6-r1` = `N196-r1`,
         `FB6-r2` = `N196-r2`)


# Now pivot longer so that we can run through the calcs on the two runs
negative_long <- negative_small %>%
  pivot_longer(cols=c("N133-r1":"N146-r2", "N151-r1":"N162-r2", "N166-r1":"N195-r2", "N198-r1":"N215-r2", "W101-r1":"W138-r2", "W164-r1":"W183-r2"), names_to=c("participant_id", "run"), names_sep="-", values_to="value")

# this gives us one row per person-run
# now pivot wider again so that we get one row per person, with a column for each run
negative_working <- negative_long %>%
  pivot_wider(names_from="run", values_from="value")

#merge in the batch numbers (as set up now we won't have the batch numbers for the field blanks)
negative_working <- negative_working %>%
  left_join(map, by="participant_id")

#MISSING batch number for W131 - that person is in batch 1
negative_working <- negative_working %>%
  mutate(batch=case_when(participant_id=="W131" ~ 1,
                         TRUE ~ batch))


#Consolidate the blank data
negative_working <- negative_working %>%
  mutate(office_field_blank_mean = (`FB1-r1`+`FB1-r2`+`FB2-r1`+`FB2-r2`+`FB3-r1`+`FB3-r2`)/6,
         nurse_field_blank_mean = (`FB4-r1`+`FB4-r2`+`FB5-r1`+`FB5-r2`+`FB6-r1`+`FB6-r2`)/6) %>%
  mutate(field_blank_mean = (office_field_blank_mean + nurse_field_blank_mean)/2) %>%
  mutate(DIW_blank_mean = (DIWB1_negB1+DIWB1_negB2+DIWB1_negB3+DIWB1_negB4+DIWB1_negB5+DIWB2_negB1+DIWB2_negB2+DIWB2_negb3+DIWB2_negB4+DIWB2_negB5)/10) %>%
  mutate(combined_blank_mean = (6/11)*field_blank_mean + (5/11)*DIW_blank_mean) %>% # 6 field blank samples and 5 DIWB samples - confirm this is how we want to weight
  mutate(office_field_blank_max = pmax(`FB1-r1`, `FB1-r2`, `FB2-r1`, `FB2-r2`, `FB3-r1`, `FB3-r2`),
         nurse_field_blank_max = pmax(`FB4-r1`, `FB4-r2`, `FB5-r1`, `FB5-r2`, `FB6-r1`, `FB6-r2`)) %>%
  mutate(field_blank_max = pmax(office_field_blank_max, nurse_field_blank_max)) %>%
  mutate(DIW_blank_max = pmax(DIWB1_negB1, DIWB1_negB2, DIWB1_negB3, DIWB1_negB4, DIWB1_negB5, DIWB2_negB1, DIWB2_negB2, DIWB2_negb3, DIWB2_negB4, DIWB2_negB5)) %>%
  mutate(max_blank = pmax(field_blank_max, DIW_blank_max))


#pull out the blank data we will need later 
negative_blanks <- negative_working %>%
  select(`Alignment ID`, max_blank, office_field_blank_mean, nurse_field_blank_mean, DIW_blank_mean) %>%
  distinct() %>%
  clean_names()
  

# standard deviation of 2 numbers simplifies to |a-b|/sqrt(2)
# relative standard deviation is sd/mean
# create an average of the two runs for each person-feature and relative standard deviation between the two runs
# if r1=0 and r2=0 then the relative sd will be NaN; replace these with zeros
negative_working <- negative_working %>%
  mutate(mean_intensity=(r1+r2)/2, 
         sd=(abs(r1-r2))/sqrt(2)) %>%
  mutate(relative_sd_intensity=case_when(
    mean_intensity==0 ~ 0,
    TRUE ~ sd/mean_intensity)) %>%
  select(-sd)


#drop the variables we don't need anymore and clean up the variable names
negative_working <- negative_working %>%
  select(`Alignment ID`, `Average Rt(min)`, `Average Mz`, participant_id, batch, office_field_blank_mean:relative_sd_intensity) %>%
  clean_names()
  
# Calculate feature max because we will need it to filter out noisy features
negative_features <- negative_working %>%
  group_by(alignment_id) %>%
  summarize(feature_max=max(mean_intensity))


# Remove "noisy" features
negative_features <- negative_features %>%
  left_join(negative_blanks, by="alignment_id") %>%
  filter(feature_max>=3000) %>% 
  filter((feature_max)>=(3*max_blank)) 
# Original negative data has 11,831 features
# After we remove the noisy features, we are left with 9828 features


# make a list of features we are keeping
keep_features<- negative_features$alignment_id

# let's look at the blanks
negative_blanks <- negative_blanks %>%
  filter(alignment_id %in% keep_features) 

negative_blanks <- negative_blanks %>%
  mutate(flag1 = case_when(
    nurse_field_blank_mean/office_field_blank_mean >2 ~ 1,
    TRUE ~ 0),
    flag2 = case_when(
      nurse_field_blank_mean/office_field_blank_mean <.5 ~ 1,
      TRUE ~ 0),
    flag3=case_when(
      diw_blank_mean>nurse_field_blank_mean & diw_blank_mean>office_field_blank_mean ~ 1,
      TRUE ~0
  ))

table(negative_blanks$flag1)
table(negative_blanks$flag2)
table(negative_blanks$flag3)


# Now we need to go back to the person-feature level data
# But we only want to do this for our features that made it past our filtering step above
negative_for_analysis <- negative_working %>%
  filter(alignment_id %in% keep_features) %>%
  select(alignment_id, average_rt_min, average_mz, participant_id, batch, max_blank, mean_intensity, relative_sd_intensity) %>% 
  mutate(worker=case_when(
    str_detect(participant_id, "N") ~ "nurse",
    TRUE ~ "office"))

# Calculate our two measures of detection frequency (only det1 is used in the current analysis)
negative_for_analysis <- negative_for_analysis %>%
  mutate(det1=case_when(
    mean_intensity>3000 ~ 1,
    TRUE ~ 0),
    det2 = case_when(
      mean_intensity>(3*max_blank) ~ 1,
      TRUE ~ 0)
  )

# Now calculate the detection frequency for each feature
# Also calculate other summary stats: the standard deviation (across participants), mean, 95th percentile 
# for each feature. 
negative_features2 <- negative_for_analysis %>%
  group_by(alignment_id) %>%
  summarize(feature_df1=mean(det1),
            feature_df2=mean(det2),
            feature_mean=mean(mean_intensity), 
            feature_sd=sd(mean_intensity),
            feature_95p=quantile(mean_intensity, prob=.95),
            feature_mean_rel_sd=mean(relative_sd_intensity)) 

# Flag those features with poor reproducibility -  greater than 50% average relative standard deviation between runs 1 and 2
negative_features2 <- negative_features2 %>%
  mutate(poor_repro_flag = case_when(
    feature_mean_rel_sd > .5 ~ 1,
    TRUE ~ 0
  ))
 
table(negative_features2$poor_repro_flag) 

#How many of our features have a df >= 65%?
negative_features2 <- negative_features2 %>%
  mutate(det1_65_flag=case_when(
    feature_df1>=.65 ~ 1,
    TRUE ~ 0))

table(negative_features2$det1_65_flag) 


# Also compute summary stats separately for nurses and office workers
# Compare the detection frequencies between nurses and office workers across all the features
negative_features3 <- negative_for_analysis %>%
  select(alignment_id, participant_id, mean_intensity, worker, det1, det2) %>%
  group_by(alignment_id, worker) %>%
  summarize(fw_df1=mean(det1),
            fw_df2=mean(det2),
            fw_mean=mean(mean_intensity),
            fw_sd=sd(mean_intensity),
            fw_95p=quantile(mean_intensity, prob=.95)) %>%
  pivot_wider(names_from="worker", values_from = c("fw_df1", "fw_df2", "fw_mean", "fw_sd", "fw_95p"))



# Replace the 0s with NaNs
negative_for_imputation <- negative_for_analysis %>%
  mutate(mean_intensity = case_when(mean_intensity > 0 ~ mean_intensity,
                                    TRUE ~ NaN))

#Imputation
# Next step is to do the imputation, only for those features with df >= 65%

# First, take ln of the average intensity
# Then, merge in the df for each feature
negative_for_imputation <- negative_for_imputation %>%
  mutate(ln_intensity=log(mean_intensity)) %>%
  left_join(negative_features2, by="alignment_id")

#Now filter to only select those features with df >= 65%
negative_for_imputation_65 <- negative_for_imputation %>%
  filter(det1_65_flag==1)


# Function to do the imputation
fillNaN_with_unifrand <- function(df) {
  lower <- 0
  upper <- min(df, na.rm=TRUE)
  a <- as.matrix(df)
  m <- is.nan(a)
  mu <- min(df, na.rm=TRUE)
  sigma <- sd(df, na.rm=TRUE)
  a[m] <- rtruncnorm(sum(m), a = lower, b = upper, mean = mu, sd = sigma)
  return(a[,1])
}
#rtruncnorm(n, a=-Inf, b=Inf, mean = 0, sd = 1)
# generates random deviates from the truncated normal distribution
# n is the number of observations
# a is vector of lower bounds
# b is vector of upper bounds
# mean is vector of means
# sd is vector of standard deviations
#The numerical arguments other than n are recycled to the length of the result


# do the imputation 3 times and then take the average
# need to lock this down (because it changes every time)
# going to output a .rds file, then comment out the imputation code, then read in the .rds file
#negative_imputed_65 <- negative_for_imputation_65 %>%
 # group_by(alignment_id) %>%
#  mutate(ln_imputed_value1=fillNaN_with_unifrand(ln_intensity),
#         ln_imputed_value2=fillNaN_with_unifrand(ln_intensity),
#         ln_imputed_value3=fillNaN_with_unifrand(ln_intensity)) %>%
#  mutate(ln_imputed_value=(ln_imputed_value1+ln_imputed_value2+ln_imputed_value3)/3) %>%
#  ungroup()
#write_rds(negative_imputed_65, "data/imputed data/negative_imputed_65.rds")

negative_imputed_65<- read_rds("data/imputed data/negative_imputed_65.rds")


# Next step is batch correction using ComBat() 
#https://academic.oup.com/nargab/article/2/3/lqaa078/5909519
#https://rdrr.io/bioc/sva/man/ComBat.html

# note: to install sva package need the following chunk of code:
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("sva")

negative_batch_corrected <- negative_imputed_65 %>%
  select(alignment_id, participant_id, batch, ln_imputed_value, worker)

n_feat=length(unique(negative_batch_corrected$alignment_id))

#first, set up the matrix
neg_mat1 <- negative_batch_corrected %>%
  pivot_wider(names_from=alignment_id, names_prefix="feature_", values_from=ln_imputed_value)

#make the batch variable
batch<-neg_mat1%>%
  use_series(batch)

#pull out the participant details - we will need to merge them back in later
negative_participant_variables <- neg_mat1 %>%
  select(participant_id, batch, worker)

#specify the model with adjustment variables
mod=model.matrix(~as.factor(worker), data=neg_mat1)

neg_mat2 <- neg_mat1 %>%
  column_to_rownames(var="participant_id") %>%
  select(3:(n_feat+2)) %>%
  t() %>%
  as.data.frame()%>%
  ComBat(batch=batch, mod=mod, par.prior=TRUE, prior.plots=T) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "participant_id") %>%
  left_join(negative_participant_variables, by="participant_id")


# now do pca and make scatter plots so we can be sure the batch correction worked
# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/#pca-data-format

neg_pca_batch_corrected <- neg_mat2 %>%
  select(-participant_id, -batch, -worker) 

neg_pca_batch_corrected <- PCA(neg_pca_batch_corrected, scale.unit=TRUE, ncp=5, graph=FALSE)

get_eigenvalue(neg_pca_batch_corrected)
fviz_eig(neg_pca_batch_corrected, addlabels = TRUE, ylim = c(0, 50))

# get results for individuals
ind <- get_pca_ind(neg_pca_batch_corrected)

# color by worker
worker_corrected <-
  fviz_pca_ind(neg_pca_batch_corrected,
             geom.ind = "point", 
             col.ind = neg_mat2$worker, # color by worker
            # palette = c("#00AFBB", "#E7B800", "#FC4E07"),
            # addEllipses = TRUE, # Concentration ellipses
             legend.title = "") +
  labs(title="")

# color by batch
batch_corrected <-
  fviz_pca_ind(neg_pca_batch_corrected,
             geom.ind = "point", 
             col.ind = neg_mat2$batch, # color by worker
             # palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             # addEllipses = TRUE, # Concentration ellipses
             legend.title = "batch") +
  labs(title="") #can also do x= and y= inside labs()


#now do the PCA for the non-batch corrected data
neg_pca_non_batch_corrected <- neg_mat1 %>%
  select(4:(n_feat+3))

neg_pca_non_batch_corrected <- PCA(neg_pca_non_batch_corrected, scale.unit=TRUE, ncp=5, graph=FALSE)

fviz_eig(neg_pca_non_batch_corrected, addlabels = TRUE, ylim = c(0, 50))

worker_orig <-
fviz_pca_ind(neg_pca_non_batch_corrected,
             geom.ind = "point", 
             col.ind = neg_mat2$worker, # color by worker
             # palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             # addEllipses = TRUE, # Concentration ellipses
             legend.title = "") +
  labs(title="")

batch_orig <-
fviz_pca_ind(neg_pca_non_batch_corrected,
             geom.ind = "point", 
             col.ind = neg_mat2$batch, # color by worker
             # palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             # addEllipses = TRUE, # Concentration ellipses
             legend.title = "batch") +
  labs(title="")

# make the plot for the supplemental
(worker_orig | worker_corrected ) / (batch_orig | batch_corrected) + plot_layout(guides = "collect") + plot_annotation(title="              before batch correction                        after batch correction")
#ggsave(filename='figures/final paper figures/figureS1.svg', width=7, height=8.666, units="in", dpi=300)


# Now we are finally ready to compare nurses and office workers

# Gate 1 - Compare Mean Abundances
# Using the batch-corrected data (only features with df >= 65%)
neg_batch_corrected_long <- neg_mat2 %>%
  pivot_longer(cols = starts_with("feature_"), names_to = "feature_ID", values_to = "ln_intensity") %>%
  mutate(nurse=case_when(
    worker=="nurse" ~ 1,
    TRUE ~ 0
  ))

# run the glm models but correct the p-values
#https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust
# Benjamini & Hochberg (1995)
glm_model_neg <- neg_batch_corrected_long %>%
  nest_by(feature_ID) %>%
  mutate(model=list(glm(ln_intensity ~ nurse, data=data))) %>%
  ungroup() %>%
  transmute(feature_ID, modelcoef=map(model, tidy)) %>%
  unnest(modelcoef) %>%
  mutate(p.value.adj=p.adjust(p.value, method="BH"))


#Make a volcano plot
# https://en.wikipedia.org/wiki/Volcano_plot_(statistics)
# exp(coef) is the fold change
# But we want to plot negative logarithm (usually base 10) of the p-value on 
# y-axis - this results in data points with low p values (highly significant)
# appearing towards the top of the plot; x-axis is log 2 fold change
volcano_neg <- glm_model_neg %>%
  filter(term=="nurse") %>%
  mutate(fold=exp(estimate)) %>%
  mutate(neg_log10_padj=-log10(p.value.adj),
         log2fold=log2(fold)) %>%
  mutate(pflag=case_when(
    p.value<0.05 ~"unadj p-value sig",
    TRUE ~ "not sig"))
#Export data for final paper figures
#write_rds(volcano_neg, "data/volcano_neg.rds")


#Pull out the features that make it past the gate (use 5% significance, adjusted p-value
# and only select those features where the coefficient on nurse is positive)
gate1 <- glm_model_neg %>%
  filter(term=="nurse") %>%
  mutate(gate1a_flag = case_when(
    p.value.adj<0.05 & estimate >0 ~ "g1a_pass",
    TRUE ~ "g1a_fail")) %>%
  select(feature_ID, gate1a_flag)


# Gate 2 - Compare the 95th percentile between nurses and office workers
# Look at the fold change between the nurses and office workers (so can use 
# any feature that has a detected 95th percentile for nurses and for office workers) 
# flag those with fold change >=3
gate2 <- negative_features3 %>%
  ungroup() %>%
  filter(fw_95p_nurse>=3000) %>%
  filter(fw_95p_office>=3000) %>% #8151 features eligible
  mutate(fold_change_95=fw_95p_nurse/fw_95p_office) %>%
  mutate(gate2_flag = case_when(
    fold_change_95>=3 ~ "g2_pass",
    TRUE ~ "g2_fail"))

gate2 <- gate2  %>%
  mutate(feature_ID=paste0("feature_", as.numeric(alignment_id))) %>%
  select(feature_ID, gate2_flag)


#Gate 3 - Include all features except those that are 100% detected
# because fisher.test will not work if all detects, or all non-detects (but we don't have any that are all non-detects)
align_ids_for_det_analysis <- negative_features2 %>%
  filter(feature_df1!=1) %>%
  select(alignment_id)
#8599

df_compare <- negative_for_analysis %>%
  filter(alignment_id %in% align_ids_for_det_analysis$alignment_id) %>%
  select(alignment_id, participant_id, worker, det1)

num_feat_all <- length(unique(df_compare$alignment_id))

df_compare <- df_compare %>%
  nest_by(alignment_id) 

gate3_id<-NULL
gate3_padj<-NULL 

#Want a one-sided test that nurses are detected more frequently 
#because of how worker coded, we specify alternative="less"
for(i in 1:num_feat_all) {
  gate3_id[i]<-df_compare[[1]][[i]]
  df<-df_compare[[2]][[i]]
  fishtest<- fisher.test(df$worker, df$det1, alternative="less") 
  p<- fishtest$p.value
  gate3_padj[i]<-p.adjust(p, method="BH")
}  

gate3 <- as.data.frame(cbind(gate3_id, gate3_padj))
colnames(gate3) <- c("alignment_id", "adjusted_p")

gate3 <- gate3 %>%
  mutate(gate3_flag = case_when(
    adjusted_p<0.05 ~ "g3_pass",
    TRUE ~ "g3_fail"
  ))
table(gate3$gate3_flag)
#558 features make it through the gate

gate3 <- gate3 %>%
  mutate(feature_ID=paste0("feature_", as.numeric(alignment_id))) %>%
  select(feature_ID, gate3_flag)

#Combine all the gates
gates <-  gate1 %>%
  full_join(gate2, by="feature_ID") %>%
  full_join(gate3, by="feature_ID")

# pull out the features that pass at least one gate
gates_passers <- gates %>%
  filter(gate1a_flag=="g1a_pass" | gate2_flag=="g2_pass" | gate3_flag=="g3_pass") %>%
  select(feature_ID)


###################################################################

# Tentative Identification of Prioritized Features

###################################################################

# Bring in the MSMS data from the pooled samples
# First, this is the DDA ("autofrag") data from December 2023
autofrag_neg <- read_csv("data/Feb2024_files/AutoFrag_NOW_MSDial_Area_NEG_20242231315.csv") %>%
  clean_names()

# We just want to use those lines where ms_ms_matched==TRUE (n=121)
#Step 1 is to merge based on mass. First, round mass in the autofrag data
autofrag_neg_small <- autofrag_neg %>%
  filter(ms_ms_matched==TRUE) %>%
  mutate(mass_round=round(average_mz))

#pull out feature data 
neg_features_pooled_matching <- negative_for_analysis %>%
  select(alignment_id, average_mz, average_rt_min) %>%
  distinct() %>%
  mutate(feature_ID=paste0("feature_", as.numeric(alignment_id))) %>%
  filter(feature_ID %in% gates_passers$feature_ID) %>%
  rename(alignment_id_ms1run=alignment_id,
         average_mz_ms1run=average_mz, 
         average_rt_ms1run=average_rt_min) %>%
  mutate(mass_round=round(average_mz_ms1run))

neg_features_pooled_matching_autofrag <- neg_features_pooled_matching %>%
  inner_join(autofrag_neg_small, by="mass_round") 

# this formula calculates the mass difference between the detected mass and the theoretical mass in parts per million (ppm)
neg_features_pooled_matching_autofrag <- neg_features_pooled_matching_autofrag %>%
  mutate(mass_diff=(abs(average_mz_ms1run-average_mz))/average_mz*(10**6))

# for a QTOF, we want to set a threshold for 10 ppm
# anything over 10 ppm is not a match
neg_features_pooled_matching_autofrag <- neg_features_pooled_matching_autofrag %>%
  mutate(mass_flag=case_when(
    mass_diff>10 ~ 1,
    TRUE ~ 0
  ))

# so we want the mass_flag=0 - these are the matches
# and differentiate those that match within 5 versus up to 10 ppm
neg_features_pooled_matching_autofrag <- neg_features_pooled_matching_autofrag %>%
  filter(mass_flag==0) %>%
  mutate(mass_match = case_when(
    mass_diff<=5 ~ "Z",
    TRUE ~ "z"
  )) 
# this still has some duplicated feature IDs (in part because the MSMS has multiple lines with same metabolite name but different post_curation_result) 

# Step 2 - also match on retention time - the retention times have to be within 0.5 min of each other
# to be a weak match "r" and within 0.1 min of each other to be a strong match "R"
# we've decided that if the retention time difference is more than 0.5 then it's not a match
neg_features_pooled_matching_autofrag <- neg_features_pooled_matching_autofrag %>%
  mutate(rt_diff=abs(average_rt_ms1run-average_rt_min)) 

neg_features_pooled_matching_autofrag <- neg_features_pooled_matching_autofrag %>%
  filter(rt_diff<=0.5) %>%
  mutate(rt_match = case_when(
    rt_diff<=0.1 ~ "R",
    TRUE ~ "r"
  ))

# multiple matches each way
# we are going to need to look at all of them
# for example, feature 4907, 4910, 4919 matched to 8-HETE twice - 
# it's because there are multiple alignment ids in the autofrag data that matched
# to start, will need all the feature IDs in the MS1 data that had at least one match, 
# and we do not want double-counting of the feature IDs here
pooled_matchers_autofrag <- neg_features_pooled_matching_autofrag %>%
  select(feature_ID, average_mz_ms1run, average_rt_ms1run, alignment_id, average_rt_min, average_mz, metabolite_name, inchikey, dot_product, reverse_dot_product, mass_match, rt_match) %>%
  rename(MS2_id=alignment_id,
         MS2_average_rt_min=average_rt_min, 
         MS2_average_mz=average_mz) %>%
  mutate(MS2_type="autofrag")


# Now, do this all again using the 2019 targeted MS2 data
msms2019_neg <- read_csv("data/2019 MS2 data/NOW_NEG_Area_0_20245241334.csv") %>%
  clean_names()

# We just want to use those lines where ms_ms_matched==TRUE (n=25)
#Step 1 is to merge based on mass. First, round mass in the MS2 data
msms2019_neg_small <- msms2019_neg %>%
  filter(ms_ms_matched==TRUE) %>%
  mutate(mass_round=round(average_mz))

# merge with our nurse/office worker data
neg_features_pooled_matching_2019 <- neg_features_pooled_matching %>%
  inner_join(msms2019_neg_small, by="mass_round") 

# this formula calculates the mass difference between the detected mass and the theoretical mass in parts per million (ppm)
neg_features_pooled_matching_2019 <- neg_features_pooled_matching_2019 %>%
  mutate(mass_diff=(abs(average_mz_ms1run-average_mz))/average_mz*(10**6))

# for a QTOF, we want to set a threshold for 10 ppm
# anything over 10 ppm is not a match
neg_features_pooled_matching_2019 <- neg_features_pooled_matching_2019 %>%
  mutate(mass_flag=case_when(
    mass_diff>10 ~ 1,
    TRUE ~ 0
  ))

# so we want the mass_flag=0 - these are the matches
# also want to flag a mass match within 10 versus within 5
neg_features_pooled_matching_2019 <- neg_features_pooled_matching_2019 %>%
  filter(mass_flag==0) %>%
  mutate(mass_match = case_when(
    mass_diff<=5 ~ "Z",
    TRUE ~ "z"
  ))

# Step 2 - also match on retention time - the retention times have to be within 0.5 min of each other
# to be a weak match "r" and within 0.1 min of each other to be a strong match "R"
# we've decided that if the retention time difference is more than 0.5 then it's not a match
neg_features_pooled_matching_2019 <- neg_features_pooled_matching_2019 %>%
  mutate(rt_diff=abs(average_rt_ms1run-average_rt_min)) 

neg_features_pooled_matching_2019 <- neg_features_pooled_matching_2019 %>%
  filter(rt_diff<=0.5) %>%
  mutate(rt_match = case_when(
    rt_diff<=0.1 ~ "R",
    TRUE ~ "r"
  ))

pooled_matchers_2019 <- neg_features_pooled_matching_2019 %>%
  select(feature_ID, average_mz_ms1run, average_rt_ms1run, alignment_id, average_rt_min, average_mz, metabolite_name, inchikey, dot_product, reverse_dot_product, mass_match, rt_match) %>%
  rename(MS2_id=alignment_id,
         MS2_average_rt_min=average_rt_min, 
         MS2_average_mz=average_mz) %>%
  mutate(MS2_type="2019")

pooled_matchers <- rbind(pooled_matchers_autofrag, pooled_matchers_2019) %>%
  distinct()


###############################################################
###############################################################

# Aside - match again, but this time against the MONA database
# First, with the DDA ("autofrag") data
msmsMONA_neg <- read_csv("data/autofrag_mona_MS2_data/Area_2_2024109615_NEG.csv") %>%
  clean_names()

# We just want to use those lines where ms_ms_matched==TRUE (n=69)
#Step 1 is to merge based on mass. First, round mass in the MS2 data
msmsMONA_neg_small <- msmsMONA_neg %>%
  filter(ms_ms_matched==TRUE) %>%
  mutate(mass_round=round(average_mz))

# merge with our nurse/office worker data
neg_features_pooled_matching_MONA <- neg_features_pooled_matching %>%
  inner_join(msmsMONA_neg_small, by="mass_round") 

# this formula calculates the mass difference between the detected mass and the theoretical mass in parts per million (ppm)
neg_features_pooled_matching_MONA <- neg_features_pooled_matching_MONA %>%
  mutate(mass_diff=(abs(average_mz_ms1run-average_mz))/average_mz*(10**6))

# for a QTOF, we want to set a threshold for 10 ppm
# anything over 10 ppm is not a match
neg_features_pooled_matching_MONA <- neg_features_pooled_matching_MONA %>%
  mutate(mass_flag=case_when(
    mass_diff>10 ~ 1,
    TRUE ~ 0
  ))

# so we want the mass_flag=0 - these are the matches
# also want to flag a mass match within 10 versus within 5
neg_features_pooled_matching_MONA <- neg_features_pooled_matching_MONA %>%
  filter(mass_flag==0) %>%
  mutate(mass_match = case_when(
    mass_diff<=5 ~ "Z",
    TRUE ~ "z"
  ))

# Step 2 - also match on retention time - the retention times have to be within 0.5 min of each other
# to be a weak match "r" and within 0.1 min of each other to be a strong match "R"
# we've decided that if the retention time difference is more than 0.5 then it's not a match
neg_features_pooled_matching_MONA <- neg_features_pooled_matching_MONA %>%
  mutate(rt_diff=abs(average_rt_ms1run-average_rt_min)) 

neg_features_pooled_matching_MONA <- neg_features_pooled_matching_MONA %>%
  filter(rt_diff<=0.5) %>%
  mutate(rt_match = case_when(
    rt_diff<=0.1 ~ "R",
    TRUE ~ "r"
  ))

pooled_matchers_MONA <- neg_features_pooled_matching_MONA %>%
  select(feature_ID, average_mz_ms1run, average_rt_ms1run, alignment_id, average_rt_min, average_mz, metabolite_name, inchikey, dot_product, reverse_dot_product, mass_match, rt_match) %>%
  rename(MS2_id=alignment_id,
         MS2_average_rt_min=average_rt_min, 
         MS2_average_mz=average_mz) %>%
  mutate(MS2_type="MONA")
# Above manually curated in MS-Dial; no new matches

#############################################
# Do again - 2019 targeted MSMS data with MONA
targeted_MONA_neg <- read_csv("data/2019 MS2 data/MONA/Area_3_202412201440_NEG.csv") %>%
  clean_names()

# We just want to use those lines where ms_ms_matched==TRUE (n=12)
#Step 1 is to merge based on mass. First, round mass in the MS2 data
target_MONA_neg_small <- targeted_MONA_neg %>%
  filter(ms_ms_matched==TRUE) %>%
  mutate(mass_round=round(average_mz))

# merge with our nurse/office worker data
neg_features_pooled_matching_MONA2019 <- neg_features_pooled_matching %>%
  inner_join(target_MONA_neg_small, by="mass_round") 

# this formula calculates the mass difference between the detected mass and the theoretical mass in parts per million (ppm)
neg_features_pooled_matching_MONA2019 <- neg_features_pooled_matching_MONA2019 %>%
  mutate(mass_diff=(abs(average_mz_ms1run-average_mz))/average_mz*(10**6))

# for a QTOF, we want to set a threshold for 10 ppm
# anything over 10 ppm is not a match
neg_features_pooled_matching_MONA2019 <- neg_features_pooled_matching_MONA2019 %>%
  mutate(mass_flag=case_when(
    mass_diff>10 ~ 1,
    TRUE ~ 0
  ))

# so we want the mass_flag=0 - these are the matches
# also want to flag a mass match within 10 versus within 5
neg_features_pooled_matching_MONA2019 <- neg_features_pooled_matching_MONA2019 %>%
  filter(mass_flag==0) %>%
  mutate(mass_match = case_when(
    mass_diff<=5 ~ "Z",
    TRUE ~ "z"
  ))

# Step 2 - also match on run time - the retention times have to be within 0.5 min of each other
# to be a weak match "r" and within 0.1 min of each other to be a strong match "R"
# we've decided that if the retention time difference is more than 0.5 then it's not a match
neg_features_pooled_matching_MONA2019 <- neg_features_pooled_matching_MONA2019 %>%
  mutate(rt_diff=abs(average_rt_ms1run-average_rt_min)) 

neg_features_pooled_matching_MONA2019 <- neg_features_pooled_matching_MONA2019 %>%
  filter(rt_diff<=0.5) %>%
  mutate(rt_match = case_when(
    rt_diff<=0.1 ~ "R",
    TRUE ~ "r"
  ))

pooled_matchers_MONA2019 <- neg_features_pooled_matching_MONA2019 %>%
  select(feature_ID, average_mz_ms1run, average_rt_ms1run, alignment_id, average_rt_min, average_mz, metabolite_name, inchikey, dot_product, reverse_dot_product, mass_match, rt_match) %>%
  rename(MS2_id=alignment_id,
         MS2_average_rt_min=average_rt_min, 
         MS2_average_mz=average_mz) %>%
  mutate(MS2_type="MONA")
# Above manually curated in MS-Dial; no new matches



###############################################################
###############################################################

# Match against DA's standards data, negative file
#Load in the standards
#Need "neg" sheet
standards_neg <- read_xlsx("data/Dimitri MSMS Standards Data/ENTACT_phthalates_MSMS_preliminary.xlsx", sheet="neg") %>%
  clean_names() %>%
  select(group_id,compound_id, m_h, rt_neg) %>%
  rename(compound=compound_id,
         m_h_standard=m_h, 
         rt_standard=rt_neg)


#no separate phthalates sheet for the negative files
#(although there are phthalates included in the negative sheet, there's just not a separate sheet of them)






standards_neg <- standards_neg %>%
  mutate(mass_round=round(m_h_standard))

# merge with our MS1 nurse/office worker data on the gatepassers
neg_features_against_standards <- neg_features_pooled_matching %>%
  inner_join(standards_neg, by="mass_round") 

# this formula calculates the mass difference between the detected mass and the theoretical mass in parts per million (ppm)
neg_features_against_standards <- neg_features_against_standards %>%
  mutate(mass_diff=(abs(average_mz_ms1run-m_h_standard))/m_h_standard*(10**6))

# for a QTOF, we want to set a threshold for 10 ppm
# anything over 10 ppm is not a match
neg_features_against_standards <- neg_features_against_standards %>%
  mutate(mass_flag_check=case_when(
    mass_diff>10 ~ 1,
    TRUE ~ 0
  ))

# so we want the mass_flag_check=0 - these are the matches
neg_features_against_standards <- neg_features_against_standards %>%
  filter(mass_flag_check==0) 

# Step 2 - also match on retention time - the retention times have to be within 0.5 min of each other
neg_features_against_standards <- neg_features_against_standards %>%
  mutate(rt_diff=abs(average_rt_ms1run-rt_standard)) 

neg_features_against_standards <- neg_features_against_standards %>%
  filter(rt_diff<=0.5) 

# we want to export this
neg_features_against_standards <- neg_features_against_standards %>%
  rename(ms1_id=alignment_id_ms1run,
         ms1_mz=average_mz_ms1run,
         ms1_rt=average_rt_ms1run,
         standard_group=group_id, 
         standard_name=compound,
         standard_mz=m_h_standard,
         standard_rt=rt_standard) %>%
  select(ms1_id, ms1_mz, ms1_rt, standard_group, standard_name, standard_mz, standard_rt, mass_round)
#write_csv(neg_features_against_standards, "data/negative_ms1_features_matched_to_standards.csv")


# Now let's take the small set of standards that matched to the MS1 gatepassers and also try to match them to the pooled samples
neg_standards_matched_to_gatepassers <- neg_features_against_standards %>%
  select(standard_group, standard_name, standard_mz, standard_rt, mass_round) %>%
  distinct() #7 negative DA standards that matched to a MS1 gatepasser

# get the pooled MS2 data ready - 2 files, autofrag and 2019; want ms_ms_asigned==TRUE
pooled_autofrag <- autofrag_neg %>%
  filter(ms_ms_assigned==TRUE) %>%
  select(alignment_id, average_rt_min, average_mz) %>%
  mutate(mass_round=round(average_mz))

pooled_2019 <- msms2019_neg %>%
  filter(ms_ms_assigned==TRUE) %>%
  select(alignment_id, average_rt_min, average_mz) %>%
  mutate(mass_round=round(average_mz))

#Now matched the two sets of MSMS data on the pooled samples to the standards
#First the autofrag
pooled_autofrag_matched_to_standards <- pooled_autofrag %>%
  inner_join(neg_standards_matched_to_gatepassers, by="mass_round") 

# this formula calculates the mass difference between the detected mass and the theoretical mass in parts per million (ppm)
pooled_autofrag_matched_to_standards <- pooled_autofrag_matched_to_standards %>%
  mutate(mass_diff=(abs(average_mz-standard_mz))/standard_mz*(10**6))

# for a QTOF, we want to set a threshold for 10 ppm
# anything over 10 ppm is not a match
pooled_autofrag_matched_to_standards <- pooled_autofrag_matched_to_standards %>%
  mutate(mass_flag_check=case_when(
    mass_diff>10 ~ 1,
    TRUE ~ 0
  ))

# so we want the mass_flag_check=0 - these are the matches
pooled_autofrag_matched_to_standards <- pooled_autofrag_matched_to_standards %>%
  filter(mass_flag_check==0) 

# Step 2 - also match on retention time - the retention times have to be within 0.5 min of each other
pooled_autofrag_matched_to_standards <- pooled_autofrag_matched_to_standards %>%
  mutate(rt_diff=abs(average_rt_min-standard_rt)) 

pooled_autofrag_matched_to_standards <- pooled_autofrag_matched_to_standards %>%
  filter(rt_diff<=0.5) 
# This is what we want to use as our basis for checking spectra in MS-DIAL


#Now need to do this again with the targeted MSMS data on the pooled samples
pooled_2019_matched_to_standards <- pooled_2019 %>%
  inner_join(neg_standards_matched_to_gatepassers, by="mass_round") 

# this formula calculates the mass difference between the detected mass and the theoretical mass in parts per million (ppm)
pooled_2019_matched_to_standards <- pooled_2019_matched_to_standards %>%
  mutate(mass_diff=(abs(average_mz-standard_mz))/standard_mz*(10**6))

# for a QTOF, we want to set a threshold for 10 ppm
# anything over 10 ppm is not a match
pooled_2019_matched_to_standards <- pooled_2019_matched_to_standards %>%
  mutate(mass_flag_check=case_when(
    mass_diff>10 ~ 1,
    TRUE ~ 0
  ))

# so we want the mass_flag_check=0 - these are the matches
pooled_2019_matched_to_standards <- pooled_2019_matched_to_standards %>%
  filter(mass_flag_check==0) 

# Step 2 - also match on retention time - the retention times have to be within 0.5 min of each other
pooled_2019_matched_to_standards <- pooled_2019_matched_to_standards %>%
  mutate(rt_diff=abs(average_rt_min-standard_rt)) 

pooled_2019_matched_to_standards <- pooled_2019_matched_to_standards %>%
  filter(rt_diff<=0.5) 
# This is what we want to use as our basis for checking spectra in MS-DIAL


#################

# to make office worker/nurse plots for any MS1 feature use code below:
temp_plot_neg  <- negative_for_analysis %>%
  mutate(feature_ID=paste0("feature_", as.numeric(alignment_id))) %>%
  select(alignment_id, feature_ID, participant_id, mean_intensity, worker) %>%
  mutate(detect=case_when(
    mean_intensity<3000 ~ "non-detect",
    TRUE ~ "detect"
  ))
#save this df so we can use for final paper figures
#write_rds(temp_plot_neg, "data/data_to_plot_neg.rds")


###############################################################
###############################################################

#Let's make the plots we will need for the manual curation
pooled_matchers_to_plot <- pooled_matchers %>%
  select(feature_ID) %>%
  distinct()

msms_matchers_plot <- negative_for_analysis %>%
  mutate(feature_ID=paste0("feature_", as.numeric(alignment_id))) %>%
  filter(feature_ID %in% pooled_matchers_to_plot$feature_ID) %>%
  select(alignment_id, feature_ID, participant_id, mean_intensity, worker) %>%
  mutate(detect=case_when(
    mean_intensity<3000 ~ "non-detect",
    TRUE ~ "detect"
  ))


# Code below plots all of the potentially-matched features, which were used for manual curation
plotnames <- msms_matchers_plot %>%
  pull(feature_ID) %>%
  unique()

feature_numbers <- msms_matchers_plot %>%
  pull(feature_ID) %>%
  unique()


num_features_plot <- length(plotnames)

for (i in 1:num_features_plot) {
  name <- plotnames[i]
  featnumber <- feature_numbers[i] 
  plotdata <- msms_matchers_plot %>% filter(feature_ID==featnumber)
  ggplot(plotdata, aes(x=worker, y=mean_intensity, color=worker, shape=detect)) +
    geom_point(position=position_jitter(h=NULL, w=0.075)) +
    scale_shape_manual(values=c(16,1)) +
    scale_y_continuous(trans="log10") +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) +
    theme(panel.border=element_rect(size=0.5, linetype="solid", colour="black", fill="transparent")) +
    theme(plot.title = element_text(hjust=0.5)) +
    ylab("mean intensity") +
    xlab("") + 
    theme(legend.position = "none")+
    ggtitle(name)
  ggsave(filename=paste0("figures/negative/",name,".jpg"),width=4, height=6, units="in", dpi=300)
}


#################################################

# Pull together all the feature summary data and output as a .csv to include in the supplemental
summary_of_features_NEG <- negative_features2 %>%
  select(-det1_65_flag, -feature_mean_rel_sd, -feature_df2) %>%
  left_join(negative_features3, by="alignment_id") %>%
  mutate(feature_ID=paste0("feature_", as.numeric(alignment_id))) %>%
  left_join(gates, by="feature_ID") %>%
  select(-feature_ID, -fw_df2_nurse, -fw_df2_office) %>%
  mutate(gate1=case_when(
    gate1a_flag=="g1a_fail" ~ "fail",
    gate1a_flag=="g1a_pass" ~ "pass",
    TRUE ~ "not eligible"
  ),
  gate2=case_when(
    gate2_flag=="g2_fail" ~ "fail",
    gate2_flag=="g2_pass" ~ "pass",
    TRUE ~ "not eligible"
  ),
  gate3=case_when(
    gate3_flag=="g3_fail" ~ "fail",
    gate3_flag=="g3_pass" ~ "pass",
    TRUE ~ "not eligible"
  )) %>%
  select(alignment_id, feature_df1, feature_mean, feature_sd, feature_95p, fw_df1_nurse, fw_mean_nurse,  fw_sd_nurse, fw_95p_nurse, fw_df1_office, fw_mean_office, fw_sd_office, fw_95p_office, gate1, gate2, gate3, poor_repro_flag)
#write_csv(summary_of_features_NEG, "data/summary_of_features_NEG.csv")
