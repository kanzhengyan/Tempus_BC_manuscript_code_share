#!/usr/bin/env Rscript
#
# Zhengyan 'George' Kan
# June 16, 2023
# 
# table of cohort summary statistics
#
source('Code/TEMPUS_BC.R');
library(dplyr);

datadir = 'Input';
outdir = 'Output.Manuscript';

cvfile <- paste(datadir, 'Tempus_CDKi_feasibility_covariates.20240515.csv', sep='/');
outfile <- paste(outdir, 'tempus_ms_cohort_summary_statistics.txt', sep='/');
outfile2 <- paste(outdir, 'tempus_ms_patient_sample_annot.txt', sep='/');

rdata = tempus_bc();
sdata = rdata$sdata;

# patient-level annotation (RWE "covariate" calculations)
pannot <- read.delim(cvfile, sep=',', header=T)

#-------------------------------------------
# (1) CDK4/6i & ET breakdown
#
df <- NULL;
vars <- c('received_palbociclib', 'received_abemaciclib', 'received_ribociclib',
          'letrozole_during_cdk', 'fulvestrant_during_cdk', 'anastrozole_during_cdk',
          'tamoxifen_during_cdk' , 'exemestane_during_cdk' , 'leuprolide_during_cdk' , 'goserelin_during_cdk');
pdata <- distinct(data.frame(patient_id=pannot$patient_id, pannot[,vars]));
foo <- merge(sdata[, c('sample_id', 'patient_id', 'prepost')], pdata, by='patient_id', all.x=T);
for (i in 1:length(vars)){
  var <- vars[i];
  tmp <- table(foo$prepost, foo$patient_id, foo[,var])['Pre',,]
  num1 <- table(tmp[,'1'] > 0)['TRUE'];
  tmp <- table(foo$prepost, foo$patient_id, foo[,var])['Post',,]
  num2 <- table(tmp[,'1'] > 0)['TRUE'];
  tmp <- table(foo$patient_id, foo[,var])
  num3 <- table(tmp[,'1'] > 0)['TRUE'];
  df <- rbind(df, data.frame(variable=var, pre_count=num1, post_count=num2, count=num3,
                             pre_pct=num1/200, post_pct=num2/227, pct=num3/400));
}

# (2) Prior exposure
# - prior exposure in the adjuvant/neoadjuvant setting
# - prior exposure in the metastatic setting
#
df2 <- NULL;
vars <- c('exposure_adjuvant', 'exposure_adjuvant_chemo', 'exposure_adjuvant_endocrine',
          'exposure_met_before_cdki_chemo', 'exposure_met_before_cdki_endocrine', 'exposure_met_before_cdki');
for (i in 1:length(vars)){
  var <- vars[i];
  tmp <- table(sdata$prepost, sdata$patient_id, sdata[,var])['Pre',,]
  num1 <- table(tmp[,'1'] > 0)['TRUE'];
  tmp <- table(sdata$prepost, sdata$patient_id, sdata[,var])['Post',,]
  num2 <- table(tmp[,'1'] > 0)['TRUE'];
  tmp <- table(sdata$patient_id, sdata[,var])
  num3 <- table(tmp[,'1'] > 0)['TRUE'];
  df2 <- rbind(df2, data.frame(variable=var, pre_count=num1, post_count=num2, count=num3,
                               pre_pct=num1/200, post_pct=num2/227, pct=num3/400));
}

# (3) biopsy tissue site distribution
#
tmp <- rbind(table(sdata$prepost, sdata$tissue_site_simplified),
             table(sdata$tissue_site_simplified))
rownames(tmp)[3] <- 'All';
df3 <- cbind(t(tmp), t(tmp/c(227, 200, 427)))
colnames(df3) <- c('post_count', 'pre_count', 'count',
                   'post_pct', 'pre_pct', 'pct');
df3 <- as.data.frame(df3);
df3$variable <- rownames(df3)
df3 <- df3[,c(7,2,1,3,5,4,6)];

# (4) patient age distribution
#
tmp <- distinct(data.frame(patient_id=sdata$patient_id, subcohort=sdata$subcohorts,
                           prepost=sdata$prepost));
pdata <- merge(pannot[,c('patient_id', 'age_at_cdk')], tmp, by='patient_id');
pdata$age_at_cdk <- round(pdata$age_at_cdk, 1)

m1 <- median(pdata[pdata$prepost == 'Pre',]$age_at_cdk, na.rm=T);
sd1 <- sd(pdata[pdata$prepost == 'Pre',]$age_at_cdk, na.rm=T);
m2 <- median(pdata[pdata$prepost == 'Post',]$age_at_cdk, na.rm=T);
sd2 <- sd(pdata[pdata$prepost == 'Post',]$age_at_cdk, na.rm=T);
m <- median(pdata$age_at_cdk, na.rm=T);
sd <- sd(pdata$age_at_cdk, na.rm=T);
#tmp <- as.data.frame(t(c(m1, m2, m, sd1, sd2, sd)));
#names(tmp) <- c('pre_median', 'post_median', 'median',
#                'pre_sd', 'post_sd', 'sd');
#df4 <- data.frame(variable='age_at_cdk', tmp);

res1 <- quantile(pdata[pdata$prepost == 'Pre',]$age_at_cdk, na.rm=T);
iqr1 <- paste(res1['25%'], res1['75%'], sep=',');
res2 <- quantile(pdata[pdata$prepost == 'Post',]$age_at_cdk, na.rm=T);
iqr2 <- paste(res2['25%'], res2['75%'], sep=',');
res <- quantile(pdata$age_at_cdk, na.rm=T);
iqr <- paste(res['25%'], res['75%'], sep=',');

tmp <- as.data.frame(t(c(m1, m2, m, iqr1, iqr2, iqr)));
names(tmp) <- c('pre_median', 'post_median', 'median',
                'pre_IQR', 'post_IQR', 'IQR');
df4 <- data.frame(variable='age_at_cdk', tmp);

# (5) time biopsy to start-of-cdk4/6i
#
m1 <- median(sdata[sdata$prepost == 'Pre',]$time_biopsy_to_cdki_start, na.rm=T);
sd1 <- sd(sdata[sdata$prepost == 'Pre',]$time_biopsy_to_cdki_start, na.rm=T);
m2 <- median(sdata[sdata$prepost == 'Post',]$time_biopsy_to_cdki_end, na.rm=T);
sd2 <- sd(sdata[sdata$prepost == 'Post',]$time_biopsy_to_cdki_end, na.rm=T);
#tmp <- as.data.frame(t(c(m1, m2, sd1, sd2)));
#names(tmp) <- c('pre_median', 'post_median', 'pre_sd', 'post_sd');
#df5 <- data.frame(variable='time_biopsy_to_cdki', tmp);

# interquartile range (25% - 75%)
res1 <- quantile(sdata[sdata$prepost == 'Pre',]$time_biopsy_to_cdki_start, na.rm=T);
iqr1 <- paste(res1['25%'], ',', res1['75%'], sep='');
res2 <- quantile(sdata[sdata$prepost == 'Post',]$time_biopsy_to_cdki_end, na.rm=T);
iqr2 <- paste(res2['25%'], ',', res2['75%'], sep='');
tmp <- as.data.frame(t(c(m1, m2, iqr1, iqr2)));
names(tmp) <- c('pre_median', 'post_median', 'pre_IQR', 'post_IQR');
df5 <- data.frame(variable='time_biopsy_to_cdki', tmp);

#------------------------------
# combine output
#
foo <- rbind(df, df2, df3);
foo$pre_descr <- paste(foo$pre_count, ' (', sprintf("%1.1f%%", foo$pre_pct*100), ')', sep='');
foo$post_descr <- paste(foo$post_count, ' (', sprintf("%1.1f%%", foo$post_pct*100), ')', sep='');
foo$descr <- paste(foo$count, ' (', sprintf("%1.1f%%", foo$pct*100), ')', sep='');

foo2 <- df4;
foo2$pre_descr <- paste(foo2$pre_median, ' (', foo2$pre_IQR, ')', sep='');
foo2$post_descr <- paste(foo2$post_median, ' (', foo2$post_IQR, ')', sep='');
foo2$descr <- paste(foo2$median, ' (', foo2$IQR, ')', sep='');

foo3 <- df5;
foo3$pre_descr <- paste(foo3$pre_median, ' (', foo3$pre_IQR, ')', sep='');
foo3$post_descr <- paste(foo3$post_median, ' (', foo3$post_IQR, ')', sep='');
foo3$descr <- '';

# print descriptions for table
#
out <- rbind(foo[,c('variable', 'pre_descr', 'post_descr', 'descr')],
             foo2[,c('variable', 'pre_descr', 'post_descr', 'descr')],
             foo3[,c('variable', 'pre_descr', 'post_descr', 'descr')]);
write.table(out, file=outfile, quote=F, row.names=F, sep='\t');


#------------------------------------------
# (6) sample/patient annotation table
# Supplementary Table 1
#
vars <- c('sample_id', 'patient_id', 'prepost', 'organ', 'pfs', 'pfs_event',
          'dnaseq_molecular_assay',	'rnaseq_molecular_assay',
          'exposure_adjuvant', 'exposure_adjuvant_chemo', 'exposure_adjuvant_endocrine',
          'exposure_met_before_cdki_chemo', 'exposure_met_before_cdki_endocrine', 'exposure_met_before_cdki',
          'time_biopsy_to_cdki_start', 'time_biopsy_to_cdki_end');
sdata <- sdata[,vars];

vars <- c('age_at_cdk', 'received_palbociclib', 'received_abemaciclib', 'received_ribociclib',
          'letrozole_during_cdk', 'fulvestrant_during_cdk', 'anastrozole_during_cdk',
          'tamoxifen_during_cdk' , 'exemestane_during_cdk' , 'leuprolide_during_cdk' , 'goserelin_during_cdk');
pdata <- distinct(data.frame(patient_id=pannot$patient_id, pannot[,vars]));

foo <- merge(sdata, pdata, by='patient_id');
write.table(foo, file=outfile2, quote=F, row.names=F, sep='\t');


