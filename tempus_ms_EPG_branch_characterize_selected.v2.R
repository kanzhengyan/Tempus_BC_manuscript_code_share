#!/usr/bin/env Rscript
#
# Zhengyan 'George' Kan
# June 29, 2023
#
# for EPG 20-node branches
# association analyses - KM plots, boxplot, barplot
#
datadir = 'Input';
outdir = 'Output.Manuscript';

library('ggpubr');
library(survival);
library(survminer);
library(ggplot2);
library(reshape2);
library(RColorBrewer);
source('Code/TEMPUS_BC.R');

# input/output files
plotfile <- paste(outdir, 'tempus_ms_branch_characterize_selected.v2.pdf', sep='/');

# load EPG branch & node data
infile <- paste(datadir, 'tempus_epg_node_branch.txt', sep='/');
epgdata <- read.table(infile, sep="\t", header=T)

#-------------------------------------------------------
# (1) Load updated features, add missing features
#
molfeatfile <- paste(datadir, 'Tempus_molecular_features_extended_v3.txt', sep='/');
mdata <- read.table(molfeatfile, header=T, sep='\t');
mdata <- merge(mdata, epgdata, by='sample_id')

# revise molecular feature data
mdata$pfs <- mdata$pfs1;
mdata$pfs_event <- mdata$pfs1_censored;
mdata$tissue_origin[grepl('Lymph node', mdata$tissue_origin)] = 'Lymph nodes';
mdata$organ <- ifelse(mdata$tissue_origin %in% c('Lymph nodes', 'Liver', 'Breast'), mdata$tissue_origin, 'Others');

# remove B3 (too few data points)
mdata <- mdata[mdata$icluster %in% c('IC1', 'IC2', 'IC3', 'IC4'),];
colnames(mdata) <- gsub('^h.HALLMARK', 'HALLMARK', colnames(mdata));

# specifiy color scheme
colors_branch = brewer.pal(4, 'Set1')[c(4, 3,2)]

#-----------------------------------
# (2) KM plot w/ log-rank test
#

# Graphics output file
pdf(file=plotfile, onefile=T, height=10, width=14, pointsize=12);

# % IC, Pre vs. Post
print('Pre:')
print(table(mdata[mdata$pre_post == 'Pre',]$branch));
print(prop.table(table(mdata[mdata$pre_post == 'Pre',]$branch)));
print ('Post:')
print(table(mdata[mdata$pre_post == 'Post',]$branch));
print(prop.table(table(mdata[mdata$pre_post == 'Post',]$branch)));
print ('All:')
print(table(mdata$branch))
print(prop.table(table(mdata$branch)));

# Survival analysis - only baseline features
# HR 95% CI
#
foo <- mdata[mdata$pre_post == 'Pre',];
foo$branch <- factor(foo$branch, levels = c('A', 'B', 'C'));
model <- coxph(Surv(pfs, pfs_event) ~ branch, data = foo);
HR_CI_label1 <- paste('(', round(exp(confint(model))['branchB', '2.5 %' ], 2), '-', 
                      round(exp(confint(model))['branchB', '97.5 %' ], 2), ')', sep='')
HR_CI_label2 <- paste('(', round(exp(confint(model))['branchC', '2.5 %' ], 2), '-', 
                      round(exp(confint(model))['branchC', '97.5 %' ], 2), ')', sep='')

# different way of calculating HR & p-value
label1 <- paste("B: HR=", format(coef(summary(model))["branchB", "exp(coef)"], digits = 4), 
                ", p=", format(coef(summary(model))["branchB", "Pr(>|z|)"], digits = 4, sci = TRUE), sep = "");
label2 <- paste("C: HR=", format(coef(summary(model))["branchC", "exp(coef)"], digits = 4), 
                ", p=", format(coef(summary(model))["branchC", "Pr(>|z|)"], digits = 4, sci = TRUE), sep = "");

# KM plot
fit <- survfit(Surv(pfs, pfs_event) ~ branch, data = foo);
plot <- ggsurvplot(fit, data=foo, legend.title="branch",  risk.table=T,
                   fontsize=6, tables.theme=theme_classic2(base_size=20),
                   legend.labs=levels(as.factor(foo$branch)),
                   legend="right", pval=FALSE, palette=colors_branch, size=2);
plot$plot  <- plot$plot +  
  ggplot2::annotate("text", x = 1500, y = 0.9, label = label1, size = 7) + # x and y coordinates of the text
  ggplot2::annotate("text", x = 1500, y = 0.8, label = label2, size = 7)   # x and y coordinates of the text

plot$plot <- ggpar(plot$plot, 
                   font.main = c(20, "bold"),
                   font.x = c(20, "bold"),
                   font.y = c(20, "bold"),
                   font.caption = c(20, "bold"), 
                   font.legend = c(20, "bold"), 
                   font.tickslab = c(20, "bold"))

plot;

#------------------------------------------------------
# (3) branch vs. continuous variables (boxplots)
#
vars = c('proliferative_index', 'hrd_index', 'cyt_score', 
         'estimate_tumor_purity', 'estimate_stromalscore', 'estimate_immunescore',
         'mutsig_3_sigma', 'mutsig_2_sigma', 'mutsig_13_sigma',
         'pam50_cor_luma', 'pam50_cor_lumb', 'pam50_cor_her2', 'pam50_cor_basal',
         'EXPR.ESR1', 'EXPR.PGR', 'EXPR.CCNE1', 'EXPR.CCNE2', 'EXPR.CDK6', 'EXPR.CCND1',
         'HALLMARK_ESTROGEN_RESPONSE_EARLY', 'HALLMARK_E2F_TARGETS', 'HALLMARK_MYC_TARGETS_V1',
         'HALLMARK_INFLAMMATORY_RESPONSE', 
         'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',
         'Paloma3_F1_MYC_E2F_activation', 'Paloma3_F10_estrogen_response', 'Paloma3_F8_EMT');
titles = c('Proliferative Index', 'HRD Index', 'CYT Score',
           'Tumor Purity (ESTIMATE)', 'Stromal Score (ESTIMATE)', 'Immune Score (ESTIMATE)',
           'S3 (HRD)', 'S2 (APOBEC)', 'S13 (APOBEC)', 
           'LumA', 'LumB', 'Her2E', 'Basal',
           'ESR1 (Log2TPM)', 'PGR (Log2TPM)', 'CCNE1 (Log2TPM)', 'CCNE2 (Log2TPM)', 'CDK6 (Log2TPM)', 'CCND1 (Log2TPM)',
           'Estrogen Response (Early)', 'E2F Targets', 'MYC Targets', 'Inflammatory Response', 
           'EPITHELIAL_MESENCHYMAL_TRANSITION', 'MYC/E2F Activation (F1)', 'Estrogen Response (F10)', 'EMT (F8)');
plist = alist();
comps = list(c('A', 'B'), c('B', 'C'), c('A', 'C'));
orders = c('A', 'B', 'C');
for (i in 1:length(vars)){
  var = vars[i];
  ymin = min(mdata[,var], na.rm=T);
  ymax = max(mdata[,var], na.rm=T);
  if (i > 7) { ratio = 1.6; } else { ratio = 1.4;}
  p<- ggboxplot(mdata, 'branch', var, xlab='', ylab='', title=titles[i], outlier.shape=NA, remove=NA,
                #                fill='branch', palette=colors_branch, ylim=c(ymin, ymax*ratio),
                fill='branch', palette=colors_branch, order=orders) + 
    stat_compare_means(size=4, label='p.format', comp=comps) + 
    theme(text = element_text(size = 12, face='bold')) + 
    geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);
  plist[[i]] = p;
}
print(ggarrange(plotlist=plist, nrow=2, ncol=4, legend=F))

#------------------------------------------------------------
# (4) branch vs. categorical variables (mosaic plots)
#
annot_colors = list(
  time = c('Post'='dimgrey', 'Pre'='white'),
  mut = c('mut'='firebrick2', 'wt'='lightskyblue'),
  amp = c('amp'='orange3', 'wt'='lightskyblue'),
  del = c('wt'='lightskyblue', 'del'='darkblue'),
  subtype = c('Basal'='purple', 'Her2'='green', 'LumA'='yellow', 'LumB'='orange', 'Normal'='pink'),
  brca_status = c('0'='lightskyblue', '1'='firebrick2'),
  icluster = brewer.pal(4, 'Set1'),
  mutsig_status = rainbow(6),
  histology = rainbow(3),
  organ = rainbow(4)
)

# pre/post vs. branch
var = 'branch';
mdata$time = mdata$pre_post;
res <- chisq.test(table(mdata$time, mdata[,var]));
title = paste(var, ': p=', round(res$p.value, 5), sep='');
par(mfrow=c(2,3), cex.main=1.5, cex.sub=1.5, cex.lab=1.5, font=2);
mosaicplot(table(mdata$time, mdata[,var])[c(2,1),], main=title, col=colors_branch, las=1, cex.axis=2, border='grey');

# icluster vs. branch
res <- chisq.test(table(mdata$icluster, mdata[,var]));
title = paste(var, ': p=', round(res$p.value, 5), sep='');
mosaicplot(table(mdata$icluster, mdata[,var]), main=title, col=colors_branch, las=1, cex.axis=2, border='grey');

# branch vs. categorical variables
vars = c('time', 'histology', 'organ', 
         'subtype_pam50', 'mutsig_sigma_status', 'icluster');
colnames = c('time', 'histology', 'organ',
             'subtype', 'mutsig_status', 'icluster');
for (i in 1:length(vars)){
  var = vars[i];
  colname = colnames[i];
  res <- chisq.test(table(mdata[,var], mdata[,'branch']));
  title = paste(var, ': p=', round(res$p.value, 5), sep='');
  par(cex.main=1.5, cex.sub=1.5, cex.lab=1.5, font=2);
  if (colname %in% c('mut', 'amp')) { 
    y = ifelse(mdata[,var] == 1, colname, 'wt');
  } else {
    y = mdata[,var];
  }
  mosaicplot(table(mdata[,'branch'], y), main=title, col=annot_colors[[colname]], ylab='', las=1, cex.axis=1.7, border='grey');
}

#-------------------------------------
# branch vs. genomic features

# all samples
par(mfrow=c(2,5), cex.main=1.5, cex.sub=1.5, cex.lab=1.5, font=2);

vars = c('brca_somatic_germline',
         'MUT.TP53', 'MUT.ESR1', 'MUT.RB1', 'AMP.MYC', 'MUT.PIK3CA', 'MUT.PTEN', 'MUT.GATA3',
         'AMP.FGFR1', 'AMP.CCND1', 'AMP.MDM2', 'DEL.PTEN');
colnames = c('brca_status', 'mut', 'mut', 'mut', 'amp', 'mut', 'mut', 'mut',
             'amp', 'amp', 'amp', 'del');
for (i in 1:length(vars)){
  var = vars[i];
  colname = colnames[i];
  res <- chisq.test(table(mdata[,var], mdata[,'branch']));
  #title = paste(var, ': p=', round(res$p.value, 5), sep='');
  title = var;
  par(cex.main=1.5, cex.sub=1.5, cex.lab=1.5, font=2);
  if (colname %in% c('mut', 'amp')) { 
    y = ifelse(mdata[,var] == 1, colname, 'wt');
  } else {
    y = mdata[,var];
  }
  mosaicplot(table(mdata[,'branch'], y), main=title, col=annot_colors[[colname]], 
             ylab='', las=1, cex.axis=1.7, border='grey');
}

# Post only
mdata2 <- mdata[mdata$pre_post == 'Post',];
for (i in 1:length(vars)){
  var = vars[i];
  colname = colnames[i];
  res <- chisq.test(table(mdata2[,var], mdata2[,'branch']));
  title = paste(var, ': p=', round(res$p.value, 5), sep='');
  #title = var;
  par(cex.main=1.5, cex.sub=1.5, cex.lab=1.5, font=2);
  if (colname %in% c('mut', 'amp')) { 
    y = ifelse(mdata2[,var] == 1, colname, 'wt');
  } else {
    y = mdata2[,var];
  }
  mosaicplot(table(mdata2[,'branch'], y), main=title, col=annot_colors[[colname]], 
             ylab='', las=1, cex.axis=1.7, border='grey');
}

#------------------------------------------------------
# (5) branch vs. continuous variables (boxplots)
# multi-panel
#
vars_exp = c('EXPR.ESR1', 'EXPR.PGR', 'EXPR.CCNE1');
varlabels_exp = c('ESR1', 'PGR', 'CCNE1' );
vars_gsva = c('HALLMARK_ESTROGEN_RESPONSE_EARLY', 'HALLMARK_E2F_TARGETS', 'HALLMARK_MYC_TARGETS_V1');
varlabels_gsva = c('Estrogen Response', 'E2F Targets', 'MYC Targets');
vars_paloma3 = c('Paloma3_F1_MYC_E2F_activation', 'Paloma3_F10_estrogen_response');
varlabels_paloma3 = c('MYC/E2F Activation (F1)', 'Estrogen Response (F10)');

list_vars <- list(vars_exp, vars_gsva, vars_paloma3);
list_varlabels <- list(varlabels_exp, varlabels_gsva, varlabels_paloma3);

plist = alist();
comps = list(c('A', 'B'), c('B', 'C'), c('A', 'C'));
orders = c('A', 'B', 'C');
for (i in 1:3){
  
  vars <- list_vars[[i]];
  varlabels <- list_varlabels[[i]];
  if (i == 1) {ylabel = 'Log2(TPM)';} else if(i == 2) { ylabel = 'GSVA';} else { ylabel = 'Factor weight'}
  
  # append data into one
  for (j in 1:length(vars)){
    var = vars[j];
    df <- data.frame(mdata, value=mdata[,var], varlabel=varlabels[j]);
    if (j == 1) { foo <- df; } else { foo <- rbind(foo, df)}
  }
  foo$varlabel = factor(foo$varlabel, levels=varlabels)
  
  plist[[i]]<- ggboxplot(foo, 'branch', 'value', xlab='', ylab=ylabel, title='', outlier.shape=NA, remove=NA,
                         facet.by='varlabel', scales='free_y', ncol=length(vars),
                         fill='branch', palette=colors_branch, order=orders) + 
    stat_compare_means(size=5, label='p.signif', comp=comps) + 
    theme(text = element_text(size = 12, face='bold')) + 
    scale_y_continuous(expand = expansion(mult = c(.1, .15))) +
    geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);
}
print(ggarrange(plotlist=plist[1:2], nrow=2, ncol=2, legend=F))
print(ggarrange(plotlist=plist[3], nrow=2, ncol=2, legend=F, widths=c(2,3)))


#-------------------------------------
# (6) branch vs. genomic features
# grouped stacked barplots
#

# Mutation

vars = c('MUT.TP53', 'MUT.ESR1', 'MUT.RB1', 'MUT.GATA3');
for (i in 1:length(vars)){
  tmp <- melt(prop.table(table(mdata[,'branch'], mdata[,vars[i]]), 1))
  tmp$Var3 <- vars[i]
  if (i == 1){
    foo <- tmp;
  } else {
    foo <- rbind(foo, tmp)
  }
}
foo$Var2 <- ifelse(foo$Var2 == 0, 'N', 'Y')
foo$Var3 <- gsub('MUT.', '', foo$Var3);
colnames(foo) <- c('Branch', 'Mutation', 'Proportion', 'Gene')
foo$Gene <- factor(foo$Gene, levels = c('TP53', 'ESR1', 'RB1', 'GATA3'));

p1 <- ggplot(foo, aes(x = Branch, y = Proportion, fill = Mutation)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values=c('lightskyblue', 'firebrick2')) +
  facet_grid(~ Gene) + ylab('% Samples') +
  theme(text = element_text(size = 15, face='bold')) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))

# Amplification

vars = c('AMP.FGFR1', 'AMP.CCND1', 'AMP.MYC');
for (i in 1:length(vars)){
  tmp <- melt(prop.table(table(mdata[,'branch'], mdata[,vars[i]]), 1))
  tmp$Var3 <- vars[i]
  if (i == 1){
    foo <- tmp;
  } else {
    foo <- rbind(foo, tmp)
  }
}
foo$Var2 <- ifelse(foo$Var2 == 0, 'N', 'Y')
foo$Var3 <- gsub('AMP.', '', foo$Var3);
colnames(foo) <- c('Branch', 'Amp', 'Proportion', 'Gene')
foo$Gene <- factor(foo$Gene, levels = c('CCND1', 'MYC', 'FGFR1'));

p2 <- ggplot(foo, aes(x = Branch, y = Proportion, fill = Amp)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values=c('lightskyblue', 'orange3')) +
  facet_grid(~ Gene) + ylab('% Samples') +
  theme(text = element_text(size = 15, face='bold')) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))

print(ggarrange(p1, p2, nrow=3, ncol=2, widths=c(4,3)))

dev.off();

