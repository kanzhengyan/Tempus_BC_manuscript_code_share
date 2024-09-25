#!/usr/bin/env Rscript
#
# Zhengyan 'George' Kan
# Jun 28, 2023
# Sept. 20, 2024
#
# characterization analysis of pseudotime (Monocle)
#
setwd('../Fig4.FigS10.pseudotime_characterize');
datadir = 'Input';
outdir = 'Output.Manuscript';

library('ggpubr');
library(RColorBrewer);
library(gplots);
library(trend);
library(reshape2);
library(survival);
library(survminer);

# input/output files
plotfile <- paste(outdir, 'Fig-S10.pseudotime_characterize.pdf', sep='/');

molfeatfile <- paste(datadir, 'Tempus_molecular_features_extended.txt', sep='/');
mdata <- read.table(molfeatfile, header=T, sep='\t');
mdata$pfs <- mdata$pfs1;
mdata$pfs_event <- mdata$pfs1_censored;

# Graphics output file
pdf(file=plotfile, onefile=T, height=10, width=14, pointsize=12);

#-----------------------------------------
# (1) Survival analysis - baseline
#
foo <- mdata[mdata$pre_post == 'Pre',];

# median-split
median <- quantile(foo$pseudotime, probs=c(0.5), na.rm=T);
foo$pt <- ifelse(foo$pseudotime > median, 1, 0)
model <- coxph(Surv(pfs, pfs_event) ~ pt, data = foo);
label1 <- paste("HR=", format(coef(summary(model))["pt", "exp(coef)"], digits = 4), 
                ", p=", format(coef(summary(model))["pt", "Pr(>|z|)"], digits = 4, sci = TRUE), sep = "");
labels_legend <- c(paste('Low (<=', round(median, 3), ')', sep=),
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
colors_prepost_pfe = c('Pre'='lightgrey', 'During'='grey60', 'Post'='grey20');
colors_pam50 = c("purple", "green", "yellow", "orange", "pink");
colors_icluster = c(brewer.pal(4, 'Set1'), 'grey');

# Pre/Post
var <- 'pseudotime';
varlabel <- 'Pseudotime';
comps1 = list(c('Pre', 'Post'));
p1 = ggboxplot(mdata, "pre_post", var, xlab='', ylab='Pseudotime', title=varlabel, outlier.shape=NA,
               fill='pre_post', palette=colors_prepost,
               order=c('Pre', 'Post')) +
  stat_compare_means(comparison=comps1, size=5, label='p.format') +
  theme(text = element_text(size = 15, face='bold'), legend.position = "none") +
  geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);

# longitudial change vs. Pre/Post
tmp <- table(mdata$patient_id);
pids_paired <- names(tmp[tmp > 1])
mdata_paired = mdata[mdata$subcohorts == '1-Paired' & mdata$patient_id %in% pids_paired,]
p2 <- ggpaired(mdata_paired, x='pre_post', y='pseudotime', title='', fill='pre_post', 
                       id='patient_id', xlab='', ylab='', point.size=2, 
                       line.color='gray', line.size=.4, palette=colors_prepost, order=c('Pre', 'Post')) +
  stat_compare_means(paired=TRUE, size=5, label.y=9) +
  theme(text = element_text(size = 15, face='bold'), legend.position = "none");

# Pre/During/Post
ymax <- 12;
ymin <- 0;
comps1 = list(c('Pre', 'During'), c('During', 'Post'), c('Pre', 'Post'));
idx <- mdata$pre_post_pfe %in% c('Pre', 'During', 'Post');
p3 = ggboxplot(mdata[idx,], "pre_post_pfe", var, xlab='', ylab='Pseudotime', title=varlabel, outlier.shape=NA,
               fill='pre_post_pfe', palette=colors_prepost_pfe, ylim=c(ymin, ymax),
               order=c('Pre', 'During', 'Post')) +
  stat_compare_means(comparison=comps1, size=5, label='p.format') +
  theme(text = element_text(size = 15, face='bold'), legend.position = "none") +
  geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);


# longitudial change vs. Pre/During/Post
tmp <- dcast(patient_id~pre_post_pfe, data=mdata, fun.aggregate=length);
pids_paired <- tmp$patient_id[apply(tmp[,c('Pre', 'During', 'Post')], 1, sum) > 1]
mdata_paired = mdata[!is.na(mdata$pre_post_pfe),]

tmp2 <- tapply(mdata$pseudotime, mdata$pre_post_pfe, mean, na.rm=T);
res <- mk.test(tmp2[c('Pre', 'During', 'Post')]);
title <- paste('PT (trend: p=', round(res$p.value, 5), ')', sep='');

p4 <- ggpaired(mdata[!is.na(mdata$pre_post_pfe),], x='pre_post_pfe', y='pseudotime', title=title, fill='pre_post_pfe', 
         id='patient_id', xlab='', ylab='Pseudotime', point.size=2, 
         line.color='gray', line.size=.4, palette=colors_prepost_pfe, ylim=c(ymin, ymax),
         order=c('Pre', 'During', 'Post')) +
  stat_compare_means(paired=TRUE, size=5, label='p.format', label.y=10) +
  theme(text = element_text(size = 15, face='bold'), legend.position = "none");

print(ggarrange(p1, p2, p3, p4, nrow=2, ncol=4, widths=c(1.5, 1.5, 2, 2)));

#-----------------------------
# PT vs. PAM50
#
comps2 = list(c('LumA', 'LumB'), c('LumB', 'Her2'), c('Basal', 'Her2'), c('LumA', 'Basal'));
idx <- mdata$subtype_pam50 %in% c('LumA', 'LumB', 'Her2', 'Basal');

# trend p-value using mann-kendall test
tmp <- tapply(mdata[idx,]$pseudotime, mdata[idx,]$subtype_pam50, mean, na.rm=T);
res <- mk.test(tmp[c('LumA', 'LumB', 'Her2', 'Basal')]);
title <- paste('PAM50 (trend: p=', round(res$p.value, 5), ')', sep='');

p1 = ggboxplot(mdata[idx,], "subtype_pam50", var, xlab='', ylab='Pseudotime', title=title, outlier.shape=NA,
               fill='subtype_pam50', palette=colors_pam50[c(3,4,2,1)], order=c('LumA', 'LumB', 'Her2', 'Basal')) +
  stat_compare_means(comp=comps2, size=4, label='p.format') +
  theme(text = element_text(size = 14, face='bold'), legend.position = "none") +
  geom_smooth(method='lm', se=T, aes(group=1), col='grey') +
  geom_jitter(pch=21, fill='grey', width=0.2, alpha=.6, size=3);

#-----------------------------
# PT vs. icluster
#
comps3 = list(c('IC1', 'IC2'), c('IC1', 'IC3'), c('IC4', 'IC3'), c('IC1', 'IC4'));
idx <- mdata$icluster %in% c('IC1', 'IC2', 'IC3', 'IC4', 'IC5');
orders <- c('IC4', 'IC5', 'IC3', 'IC2', 'IC1');

# trend p-value using mann-kendall test
tmp <- tapply(mdata[idx,]$pseudotime, mdata[idx,]$icluster, mean, na.rm=T);
res <- mk.test(tmp[orders]);
title <- paste('IC1-5 (trend: p=', round(res$p.value, 5), ')', sep='');

p2 = ggboxplot(mdata[idx,], "icluster", var, xlab='', ylab='Pseudotime', title=title, outlier.shape=NA,
               fill='icluster', palette=colors_icluster[c(4,5,2,3,1)], order=orders) +
  stat_compare_means(comp=comps3, size=4, label='p.format') +
  theme(text = element_text(size = 15, face='bold'), legend.position = "none") +
  geom_smooth(method='lm', se=T, aes(group=1), col='grey') +
  geom_jitter(pch=21, fill='grey', width=0.2, alpha=.6, size=3);

comps4 = list(c('IC1', 'IC3'), c('IC1', 'IC2'));
p3 = ggboxplot(mdata[mdata$subtype_pam50 == 'LumB' & mdata$icluster %in% c('IC1', 'IC2', 'IC3'),], "icluster", var, 
               xlab='', ylab='Pseudotime', title='', outlier.shape=NA, facet.by='subtype_pam50',
               fill='icluster', palette=colors_icluster[c(2,3,1)],
               order=c('IC2', 'IC3', 'IC1')) +
  stat_compare_means(comp=comps4, size=4, label='p.format') +
  theme(text = element_text(size = 15, face='bold'), legend.position = "none") +
  geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3) + rotate_x_text(45)

print(ggarrange(p1, p2, p3, nrow=2, ncol=4, widths=c(2,2,1.5)));

#-------------------------------------
# (3) Compare vs. genomic features
#
annot_colors = list(
  mut = c('mut'='firebrick2', 'wt'='lightskyblue'),
  amp = c('amp'='orange3', 'wt'='lightskyblue'),
  del = c('wt'='lightskyblue', 'del'='darkblue'),
  brca_status = c('0'='lightskyblue', '1'='firebrick2')
)

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

vars <- c('AMP.CCND1', 'AMP.MYC', 'AMP.FGFR1');
varlabels <- c('CCND1', 'MYC', 'FGFR1')

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

print(ggarrange(p1, p2, nrow=2, ncol=2, legend=F));

#------------------------------------------------------
# (4) icluster vs. continuous variables (boxplots)
#
vars = c('EXPR.ESR1', 'EXPR.PGR', 'EXPR.CCNE1', 
         'h.HALLMARK_ESTROGEN_RESPONSE_EARLY',
         'Paloma3_F1_MYC_E2F_activation', 'Paloma3_F10_estrogen_response');
titles = c('ESR1 (Log2TPM)', 'PGR (Log2TPM)', 'CCNE1 (Log2TPM)', 'Estrogen Response (Early)', 
           'MYC/E2F Activation (F1)', 'Estrogen Response (F10)');
plist = alist();
for (i in 1:length(vars)){
  var = vars[i];
  idx <- mdata$icluster %in% c('IC1', 'IC2', 'IC3', 'IC4');
  ymax <- max(mdata[idx, var], na.rm=T);
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

