#!/usr/bin/env Rscript
#
# Zhengyan 'George' Kan
# May 10, 2023
#
# script to create boxplots to show association
# for genomic alterations (discrete) vs. molecular features (continuous)
#
# tempus_ms_assoc_genomic_vs_molecular_feature.R
#
datadir = 'Input';
outdir = 'Output.Manuscript';

library(reshape2)
library(ggpubr);
source('Code/TEMPUS_BC.R');

outfile <- paste(outdir, "tempus_ms_assoc_genomic_vs_molecular_feature.txt", sep='/')
plotfile <- paste(outdir, 'tempus_ms_assoc_genomic_vs_molecular_feature.pdf', sep='/');

#--------------------------------
# (1) get molecular features
#
molfeatfile <- paste(datadir, 'Tempus_molecular_features_extended_v3.txt', sep='/');
mdata <- read.table(molfeatfile, header=T, sep='\t');

# revise data
mdata$pfs <- mdata$pfs1;
mdata$pfs_censored <- mdata$pfs1_censored;
mdata$MUT.BRCA.Pathogenic <- mdata$brca_somatic_germline;
colnames(mdata)[grepl('h.HALLMARK', colnames(mdata))] <- gsub('h.HALLMARK_', '', colnames(mdata)[grepl('h.HALLMARK_', colnames(mdata))])
colnames(mdata)[grepl('c2_cgp.', colnames(mdata))] <- gsub('c2_cgp.', '', colnames(mdata)[grepl('c2_cgp.', colnames(mdata))])

# binary: genomic alteration, PAM50 subtype
# continuous: gene expression, gene signature, PAM50 scores, NMF factors
#
features_d = c('MUT.TP53', 'MUT.RB1', 'MUT.PIK3CA', 'MUT.AKT1', 'MUT.ESR1',
                    'MUT.BRCA.Pathogenic', 'MUT.PTEN', 'MUT.CDKN2A',
                    'AMP.CCND1', 'AMP.FGFR1', 'AMP.MDM2');
features_c = c('proliferative_index', 'cyt_score', 'hrd_index', 'tumor_purity_tempus', 'mutsig_13_sigma',
                        'pam50_cor_basal', 'pam50_cor_her2', 'pam50_cor_luma', 'pam50_cor_lumb', 'pam50_cor_normal',
                        'Paloma3_F1_MYC_E2F_activation', 'Paloma3_F10_estrogen_response', 'Paloma3_F8_EMT', 'Paloma3_F7_IFNG_response', 
                        'MYC_TARGETS_V1', 'ESTROGEN_RESPONSE_EARLY', 'E2F_TARGETS', 'SMID_BREAST_CANCER_ERBB2_UP',
                        'PI3K_AKT_MTOR_SIGNALING', 'CHOLESTEROL_HOMEOSTASIS',
                        'EXPR.CCNE1', 'EXPR.ESR1', 'EXPR.PGR');

#----------------------------------------------
# (2) calculate association statistics
#

#----------------------------------------------
# (2) create boxplots to show association
# mutations vs. expression features
#

# graphics output file
pdf(file=plotfile, onefile=T, height=9, width=10, pointsize=12);

genomic_muts = c('MUT.ESR1', 'MUT.RB1', 'MUT.TP53');
features = c('proliferative_index', 'ESTROGEN_RESPONSE_EARLY', 'E2F_TARGETS', 'SMID_BREAST_CANCER_ERBB2_UP', 
             'EXPR.CCNE1', 'EXPR.PGR', 'EXPR.ESR1', 'ESTROGEN_RESPONSE_EARLY',
             'ENRS_CERES_CDK4', 'ENRS_CERES_CDK2', 'ENRS_CERES_CDK6', 
             'ENRS_CERES_ESR1', 'ENRS_PALBOCICLIB');
colors_alter = c('deepskyblue3', 'firebrick3')

mm <- melt(mdata[,c('sample_id', genomic_muts)], id.vars=1)
mm$value <- ifelse(mm$value == 1, 'MUT', 'WT')

plist <- list();
comps <- list(c('WT', 'MUT'));
for (i in 1:length(features)){
  yvar <- features[i];
  foo <- merge(mm, mdata[,c('sample_id', yvar)], by='sample_id')
  plist[[i]] <- ggboxplot(foo, 'value', yvar, xlab='', ylab=yvar, title='', outlier.shape=NA, remove=NA,
          fill='value', palette=colors_alter, facet.by='variable') + 
    stat_compare_means(size=4, label='p.format', comp=comps) + 
    theme(text = element_text(size = 12, face='bold')) + 
    geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=2);
}

ggarrange(plotlist=plist, nrow=2, ncol=2, legend=F);

#----------------------------------------------
# (3) create boxplots to show association
# Pre/Post vs. expression features
#
features = c('EXPR.CCNE1', 'Paloma3_F1_MYC_E2F_activation', 'Paloma3_F10_estrogen_response')
colors_alter = c('deepskyblue3', 'firebrick3')

mm <- melt(mdata[,c('sample_id', features)], id.vars=1)

plist <- list();
comps <- list(c('Pre', 'Post'));

foo <- merge(mm, mdata[,c('sample_id', 'pre_post')], by='sample_id')
plist[[1]] <- ggboxplot(foo, 'pre_post', 'value', xlab='', ylab='', title='', outlier.shape=NA, remove=NA,
                          fill='pre_post', palette='Set1', facet.by='variable', ncol=4, scales='free') + 
    stat_compare_means(size=4, label='p.format', comp=comps) + 
    theme(text = element_text(size = 12, face='bold')) + 
    geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=2);

ggarrange(plotlist=plist, nrow=2, ncol=1, legend=F);

dev.off();

