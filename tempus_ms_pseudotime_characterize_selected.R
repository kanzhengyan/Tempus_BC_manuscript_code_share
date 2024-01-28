#!/usr/bin/env Rscript
#
# Zhengyan 'George' Kan
# Jun 28, 2023
#
# characterization analysis of pseudotime (Monocle)
#
source('Code/TEMPUS_BC.R');
library('ggpubr');
library(RColorBrewer);
library(gplots);
library(trend);
library(survival);
library(survminer);

datadir = 'Input';
outdir = 'Output.Manuscript';

# input/output files
plotfile <- paste(outdir, 'tempus_ms_pseudotime_characterize_selected.pdf', sep='/');

# parse PT data
infile <- paste(datadir, '8.sample.clustering.output.batchCorrected.wPseudotime.csv', sep='/');
indata <- read.csv(infile);
ptdata <- indata[,c('X', 'pseudotime')];
colnames(ptdata) <- c('sample_id', 'pseudotime');

# load Tempus rdata
rdata = tempus_bc();
sdata = rdata$sdata;

# add molecular features to sample annotation
mdata <- add_molecular_features(sdata);
mdata <- merge(mdata, ptdata, by='sample_id', sort=F);

#var_characterize_analysis(outdir, 'pseudotime', mdata, dm, 'pseudotime');

# Graphics output file
pdf(file=plotfile, onefile=T, height=10, width=14, pointsize=12);

#-----------------------------------------
# (1) Survival analysis - baseline
#
foo <- mdata[mdata$prepost == 'Pre',];

# median-split
median <- quantile(foo$pseudotime, probs=c(0.5), na.rm=T);
foo$pt <- ifelse(foo$pseudotime > median, 1, 0)
model <- coxph(Surv(pfs, pfs_event) ~ pt, data = foo);
label1 <- paste("HR=", format(coef(summary(model))["pt", "exp(coef)"], digits = 4), 
                ", p=", format(coef(summary(model))["pt", "Pr(>|z|)"], digits = 4, sci = TRUE), sep = "");
labels_legend <- c(paste('Low (<=', round(median, 3), ')', sep=''),
                   paste('High (>', round(median, 3), ')', sep=''))

# KM plot
fit <- survfit(Surv(pfs, pfs_event) ~ pt, data = foo);
plot <- ggsurvplot(fit, data=foo, legend.title="Pseudotime", risk.table=T,
                   fontsize=6, tables.theme=theme_classic(base_size=20),
                   legend.labs=labels_legend, legend='right',
                   pval=FALSE, palette=c('navyblue', 'red'), size=2);
plot$plot  <- plot$plot + 
  ggplot2::annotate("text", x = 1500, y = 0.9, label = label1, size = 7) # x and y coordinates of the text

plot$plot <- ggpar(plot$plot, 
                   font.main = c(20, "bold"),
                   font.x = c(20, "bold"),
                   font.y = c(20, "bold"),
                   font.caption = c(20, "bold"), 
                   font.legend = c(20, "bold"), 
                   font.tickslab = c(20, "bold"))
plot;

#------------------------------------------------------
# (2) - compare vs. key ccategorical features 
# - Pre/Post
# - PAM50 subtype
# - icluster
# - icluster+PAM50
#

# color labels
colors_prepost = c('Pre'='lightgrey', 'Post'='grey20');
colors_pam50 = c("purple", "green", "yellow", "orange", "pink");
colors_icluster = brewer.pal(4, 'Set1')

var <- 'pseudotime';
varlabel <- 'Pseudotime';
comps1 = list(c('Pre', 'Post'));
p1 = ggboxplot(mdata, "prepost", var, xlab='', ylab='', title=varlabel, outlier.shape=NA,
               fill='prepost', palette=colors_prepost,
               order=c('Pre', 'During', 'Post')) +
  stat_compare_means(comparison=comps1, size=4, label='p.format') +
  theme(text = element_text(size = 15, face='bold'), legend.position = "none") +
  geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);

# (2) Examine longitudial change vs. Pre/Post
#
tmp <- table(mdata$patient_id);
pids_paired <- names(tmp[tmp > 1])
mdata_paired = mdata[mdata$subcohorts == '1-Paired' & mdata$patient_id %in% pids_paired,]
p2 <- ggpaired(mdata_paired, x='prepost', y='pseudotime', title='', fill='prepost', 
                       id='patient_id', xlab='', ylab='', point.size=2, 
                       line.color='gray', line.size=.4, palette=colors_prepost, order=c('Pre', 'Post')) +
  stat_compare_means(paired=TRUE, size=4) +
  theme(text = element_text(size = 15, face='bold'), legend.position = "none");

comps1 = list(c('Pre', 'During'), c('During', 'Post'), c('Pre', 'Post'));
idx <- mdata$prepost_pfe %in% c('Pre', 'During', 'Post');
p3 = ggboxplot(mdata[idx,], "prepost_pfe", var, xlab='', ylab='', title=varlabel, outlier.shape=NA,
               fill='prepost_pfe', palette='Set2',
               order=c('Pre', 'During', 'Post')) +
  stat_compare_means(comparison=comps1, size=4, label='p.format') +
  theme(text = element_text(size = 15, face='bold'), legend.position = "none") +
  geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);


#-----------------------------
# PT vs. PAM50
#
comps2 = list(c('LumA', 'LumB'), c('LumB', 'Her2'), c('Basal', 'Her2'), c('LumA', 'Basal'));
idx <- mdata$subtype_pam50 %in% c('LumA', 'LumB', 'Her2', 'Basal');

# trend p-value using mann-kendall test
tmp <- tapply(mdata[idx,]$pseudotime, mdata[idx,]$subtype_pam50, mean);
res <- mk.test(tmp[c('LumA', 'LumB', 'Her2', 'Basal')]);
title <- paste('PAM50 (trend: p=', round(res$p.value, 5), ')', sep='');

p4 = ggboxplot(mdata[idx,], "subtype_pam50", var, xlab='', ylab='Pseudotime', title=title, outlier.shape=NA,
               fill='subtype_pam50', palette=colors_pam50[c(3,4,2,1)], order=c('LumA', 'LumB', 'Her2', 'Basal')) +
  stat_compare_means(comp=comps2, size=4, label='p.format') +
  theme(text = element_text(size = 14, face='bold'), legend.position = "none") +
  geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);

#-----------------------------
# PT vs. icluster
#
comps3 = list(c('IC4', 'IC2'), c('IC2', 'IC3'), c('IC1', 'IC3'), c('IC1', 'IC4'));
idx <- mdata$icluster %in% c('IC1', 'IC2', 'IC3', 'IC4');

# trend p-value using mann-kendall test
tmp <- tapply(mdata[idx,]$pseudotime, mdata[idx,]$icluster, mean);
res <- mk.test(tmp[c('IC4', 'IC3', 'IC2', 'IC1')]);
title <- paste('IC1-4 (trend: p=', round(res$p.value, 5), ')', sep='');

p5 = ggboxplot(mdata[idx,], "icluster", var, xlab='', ylab='Pseudotime', title=title, outlier.shape=NA,
               fill='icluster', palette=colors_icluster[c(4,2,3,1)], order=c('IC4', 'IC2', 'IC3', 'IC1')) +
  stat_compare_means(comp=comps3, size=4, label='p.format') +
  theme(text = element_text(size = 15, face='bold'), legend.position = "none") +
  geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);

#-------------------------
# icluster + PAM50
#

# color subgroup by PAM50
tmp <- matrix(unlist(strsplit(levels(as.factor(mdata$ic_pam50)), ":")), ncol=2, byrow=TRUE)[,2]
collabs <- colors_pam50[as.numeric(as.factor(tmp))];
names(collabs) <- levels(mdata$ic_pam50)

comps4 = list(c('IC1:LumB', 'IC2:LumB'), c('IC2:LumB', 'IC3:LumB'), c('IC1:LumB', 'IC3:LumB'));
p6 = ggboxplot(mdata[!is.na(mdata$ic_pam50),], "ic_pam50", var, xlab='', ylab='Pseudotime', title='', outlier.shape=NA,
               fill='ic_pam50', palette=collabs,
               order=sort(levels(as.factor(mdata$ic_pam50)))) +
  stat_compare_means(comp=comps4, size=4, label='p.format') +
  theme(text = element_text(size = 15, face='bold'), legend.position = "none") +
  geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3) + rotate_x_text(45)

print(ggarrange(p1, p2, p3, p4,
                NULL, NULL, p5, p6, nrow=2, ncol=4, widths=c(1,1.5,2,2)));

#-------------------------------------
# (3) Compare vs. genomic features
#
annot_colors = list(
  mut = c('mut'='firebrick2', 'wt'='deepskyblue2'),
  amp = c('amp'='orange3', 'wt'='deepskyblue2'),
  del = c('wt'='deepskyblue2', 'del'='darkblue'),
  brca_status = c('0'='deepskyblue2', '1'='firebrick2')
)

plist <- alist();
vars = c('brca_somatic_germline',
         'MUT.TP53', 'MUT.ESR1', 'MUT.RB1', 'AMP.MYC', 'MUT.PIK3CA', 'MUT.PTEN',
         'AMP.FGFR1', 'AMP.CCND1', 'DEL.PTEN');
colnames = c('brca_status', 'mut', 'mut', 'mut', 'amp', 'mut', 'mut',
             'amp', 'amp', 'del');
for (i in 1:length(vars)){
  var = vars[i];
  colname = colnames[i]
  if (colname %in% c('mut', 'amp', 'del')) {
    mdata[,'x'] = ifelse(mdata[,var] == 1, colname, 'wt');
    comps = list(c('wt', colname));
  } else {
    mdata[, 'x'] = mdata[,var];
    comps = list(c('0', '1'));
  }

  p<- ggboxplot(mdata, 'x', 'pseudotime', xlab='', ylab='Pseudotime', title=var, outlier.shape=NA, remove=NA,
                fill='x', palette=annot_colors[[colname]]) + 
    stat_compare_means(size=4, label='p.format', comp=comps) + 
    theme(text = element_text(size = 15, face='bold'), legend.position='none') + 
    geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);
  plist[[i]] = p;
}
print(ggarrange(plotlist=plist, nrow=2, ncol=5));

#------------------------------------------
# PT vs. genomic features (multi-panel)
#
vars <- c('MUT.TP53', 'MUT.RB1', 'MUT.ESR1', 'MUT.PIK3CA');
varlabels <- c('TP53', 'RB1', 'ESR1', 'PIK3CA')

# append data into one
for (j in 1:length(vars)){
  var = vars[j];
  df <- data.frame(mdata, value=mdata[,var], varlabel=varlabels[j]);
  if (grepl('MUT', var)) {
    df[,'x'] = ifelse(df[,var] == 1, 'mut', 'wt');
    comps = list(c('wt', 'mut'));
  } else {
    df[,'x'] = ifelse(df[,var] == 1, 'amp', 'wt');
    comps = list(c('wt', 'amp'));
  }
  if (j == 1) { foo <- df; } else { foo <- rbind(foo, df)}
}
foo$varlabel = factor(foo$varlabel, levels=varlabels)

p1 <- ggboxplot(foo, 'x', 'pseudotime', xlab='', ylab='Pseudotime', title='', outlier.shape=NA, remove=NA,
                       facet.by='varlabel', ncol=length(vars),
                       fill='x', palette=annot_colors[['mut']], order=c('wt', 'mut')) + 
  stat_compare_means(size=4, label='p.format', comp=comps) + 
  theme(text = element_text(size = 12, face='bold')) + 
  scale_y_continuous(expand = expansion(mult = c(.1, .15))) +
  geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);

vars <- c('AMP.CCND1', 'AMP.MYC', 'AMP.FGFR1', 'AMP.KAT6A');
varlabels <- c('CCND1', 'MYC', 'FGFR1', 'KAT6A')

# append data into one
for (j in 1:length(vars)){
  var = vars[j];
  df <- data.frame(mdata, value=mdata[,var], varlabel=varlabels[j]);
  if (grepl('MUT', var)) {
    df[,'x'] = ifelse(df[,var] == 1, 'mut', 'wt');
    comps = list(c('wt', 'mut'));
  } else {
    df[,'x'] = ifelse(df[,var] == 1, 'amp', 'wt');
    comps = list(c('wt', 'amp'));
  }
  if (j == 1) { foo <- df; } else { foo <- rbind(foo, df)}
}
foo$varlabel = factor(foo$varlabel, levels=varlabels)

p2 <- ggboxplot(foo, 'x', 'pseudotime', xlab='', ylab='Pseudotime', title='', outlier.shape=NA, remove=NA,
                facet.by='varlabel', ncol=length(vars),
                fill='x', palette=annot_colors[['amp']], order=c('wt', 'amp')) + 
  stat_compare_means(size=4, label='p.format', comp=comps) + 
  theme(text = element_text(size = 12, face='bold')) + 
  scale_y_continuous(expand = expansion(mult = c(.1, .15))) +
  geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);

print(ggarrange(p1, p2, nrow=2, ncol=2));

#------------------------------------------------------
# (4) icluster vs. continuous variables (boxplots)
#
vars = c('proliferative_index', 'tumor_purity', 
         'pam50_cor_luma', 'pam50_cor_lumb', 'pam50_cor_her2', 'pam50_cor_basal',
         'EXPR.ESR1', 'EXPR.PGR', 'EXPR.CCNE1', 'EXPR.CCNE2', 'EXPR.CDK6', 'EXPR.CCND1',
         'h.HALLMARK_ESTROGEN_RESPONSE_EARLY', 'h.HALLMARK_E2F_TARGETS', 'h.HALLMARK_MYC_TARGETS_V1',
         'h.HALLMARK_INFLAMMATORY_RESPONSE', 'h.HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',
         'Paloma3_F1_MYC_E2F_activation', 'Paloma3_F10_estrogen_response', 'Paloma3_F8_EMT');
titles = c('Proliferative Index', 'Tumor Purity (ESTIMATE)', 
           'LumA', 'LumB', 'Her2E', 'Basal',
           'ESR1 (Log2TPM)', 'PGR (Log2TPM)', 'CCNE1 (Log2TPM)', 'CCNE2 (Log2TPM)', 'CDK6 (Log2TPM)', 'CCND1 (Log2TPM)',
           'Estrogen Response (Early)', 'E2F Targets', 'MYC Targets', 'Inflammatory Response', 
           'EPITHELIAL_MESENCHYMAL_TRANSITION', 'MYC/E2F Activation (F1)', 'Estrogen Response (F10)', 'EMT (F8)');
plist = alist();
for (i in 1:length(vars)){
  var = vars[i];
  idx <- mdata$icluster %in% c('IC1', 'IC2', 'IC3', 'IC4');
  ymax <- max(mdata[idx, var]);
  p<- ggscatter(mdata[idx,], 'pseudotime', var, xlab=varlabel, ylab='', title=titles[i],
                add='reg.line', add.params = list(color = "blue", fill = "lightgray"),
                color='pseudotime', alpha=0.6, size=4) +
    gradient_color(c('gray', 'black')) +
    stat_cor(method='spearman', size=5, label.y=ymax*1.15) +
    theme(text = element_text(size = 15, face='bold')) 
  plist[[i]] = p;
}
print(ggarrange(plotlist=plist, nrow=2, ncol=3, legend=F))

dev.off();

