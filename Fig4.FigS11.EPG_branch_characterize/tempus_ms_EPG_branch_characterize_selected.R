#!/usr/bin/env Rscript
#
# Zhengyan 'George' Kan
# June 29, 2023
#
# for EPG 20-node branches
# association analyses - KM plots, boxplot, barplot
#
setwd('../Fig4.FigS11.EPG_branch_characterize')
datadir = 'Input';
outdir = 'Output.Manuscript';

library('ggpubr');
library(ggplot2);
library(reshape2);
library(RColorBrewer);

# input/output files
plotfile <- paste(outdir, 'Fig-4.Fig-S11.EPG_branch_characterize.pdf', sep='/');

#-------------------------------------------------------
# (1) Load updated features, add missing features
#
molfeatfile <- paste(datadir, 'Tempus_molecular_features_extended.txt', sep='/');
mdata <- read.table(molfeatfile, header=T, sep='\t');

# revise molecular feature data
mdata$pfs <- mdata$pfs1;
mdata$pfs_event <- mdata$pfs1_censored;
mdata$branch <- mdata$branch_epg20;
colnames(mdata) <- gsub('^h.HALLMARK', 'HALLMARK', colnames(mdata));

# specifiy color scheme
colors_branch = brewer.pal(4, 'Set1')[c(4, 3,2)]

#-----------------------------------
# (2) KM plot w/ log-rank test
#

# Graphics output file
pdf(file=plotfile, onefile=T, height=10, width=14, pointsize=12);

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

#------------------------------------------------------------
# (3) branch vs. categorical variables (mosaic plots)
#
annot_colors = list(
  time = c('Post'='dimgrey', 'Pre'='white'),
  mut = c('mut'='firebrick2', 'wt'='lightskyblue'),
  amp = c('amp'='orange3', 'wt'='lightskyblue'),
  del = c('wt'='lightskyblue', 'del'='darkblue'),
  subtype = c('Basal'='purple', 'Her2'='green', 'LumA'='yellow', 'LumB'='orange', 'Normal'='pink'),
  brca_status = c('0'='lightskyblue', '1'='firebrick2'),
  icluster = c(brewer.pal(4, 'Set1'), 'grey'),
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
vars = c('time', 'organ', 'subtype_pam50', 'icluster');
colnames = c('time', 'organ', 'subtype', 'icluster');
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

#------------------------------------------------------
# (4) branch vs. continuous variables (boxplots)
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

