#!/usr/bin/env Rscript
#
# Zhengyan 'George' Kan
# May 6, 2022
# June 21, 2023
# Sept. 20, 2024
#
setwd('../Fig2.FigS7-S8.icluster_characterize')
datadir = 'Input';
outdir = 'Output.Manuscript';

library('ggpubr');
library(survival);
library(survminer);
library('reshape2');
library(RColorBrewer);

# input/output files
plotfile <- paste(outdir, 'Fig-2.Fig-S7.icluster_characterize.pdf', sep='/');

#-------------------------------------------------------
# (1) Load updated features, add missing features
#
molfeatfile <- paste(datadir, 'Tempus_molecular_features_extended.txt', sep='/');
mdata <- read.table(molfeatfile, header=T, sep='\t');

# revise molecular feature data
mdata$pfs <- mdata$pfs1;
mdata$pfs_event <- mdata$pfs1_censored;
mdata$prepost <- mdata$pre_post;
colnames(mdata) <- gsub('^h.HALLMARK', 'HALLMARK', colnames(mdata));

# specifiy color scheme
colors_icluster = c(brewer.pal(4, 'Set1'), 'gray')

#-----------------------------------
# (2) KM plot w/ log-rank test
#

# Graphics output file
pdf(file=plotfile, onefile=T, height=10, width=16, pointsize=12);

# % IC, Pre vs. Post
print('Pre:')
print(table(mdata[mdata$pre_post == 'Pre',]$icluster));
print(prop.table(table(mdata[mdata$pre_post == 'Pre',]$icluster)));
print ('Post:')
print(table(mdata[mdata$pre_post == 'Post',]$icluster));
print(prop.table(table(mdata[mdata$pre_post == 'Post',]$icluster)));
print ('All:')
print(table(mdata$icluster))
print(prop.table(table(mdata$icluster)));

# Survival analysis - only baseline features
# HR 95% CI
#
foo <- mdata[mdata$pre_post == 'Pre',];
foo$icluster <- factor(foo$icluster, levels = c('IC2', 'IC1', 'IC3', 'IC4', 'IC5'));
model <- coxph(Surv(pfs, pfs_event) ~ icluster, data = foo);
HR_CI_label1 <- paste('(', round(exp(confint(model))['iclusterIC1', '2.5 %' ], 2), '-', 
                      round(exp(confint(model))['iclusterIC1', '97.5 %' ], 2), ')', sep='')
HR_CI_label2 <- paste('(', round(exp(confint(model))['iclusterIC3', '2.5 %' ], 2), '-', 
                      round(exp(confint(model))['iclusterIC3', '97.5 %' ], 2), ')', sep='')
HR_CI_label3 <- paste('(', round(exp(confint(model))['iclusterIC4', '2.5 %' ], 2), '-', 
                      round(exp(confint(model))['iclusterIC4', '97.5 %' ], 2), ')', sep='')
HR_CI_label4 <- paste('(', round(exp(confint(model))['iclusterIC5', '2.5 %' ], 2), '-', 
                      round(exp(confint(model))['iclusterIC5', '97.5 %' ], 2), ')', sep='')

# different way of calculating HR & p-value
label1 <- paste("IC1: HR=", format(coef(summary(model))["iclusterIC1", "exp(coef)"], digits = 4), 
                ", p=", format(coef(summary(model))["iclusterIC1", "Pr(>|z|)"], digits = 4, sci = TRUE), sep = "");
label2 <- paste("IC3: HR=", format(coef(summary(model))["iclusterIC3", "exp(coef)"], digits = 4), 
                ", p=", format(coef(summary(model))["iclusterIC3", "Pr(>|z|)"], digits = 4, sci = TRUE), sep = "");
label3 <- paste("IC4: HR=", format(coef(summary(model))["iclusterIC4", "exp(coef)"], digits = 4), 
                ", p=", format(coef(summary(model))["iclusterIC4", "Pr(>|z|)"], digits = 4, sci = TRUE), sep = "");
label4 <- paste("IC5: HR=", format(coef(summary(model))["iclusterIC5", "exp(coef)"], digits = 4), 
                ", p=", format(coef(summary(model))["iclusterIC5", "Pr(>|z|)"], digits = 4, sci = TRUE), sep = "");

# KM plot
fit <- survfit(Surv(pfs, pfs_event) ~ icluster, data = foo);
plot <- ggsurvplot(fit, data=foo, legend.title="iCluster",  risk.table=T,
                   ylab='PFS probability',
                   fontsize=6, tables.theme=theme_classic2(base_size=20),
                   legend.labs=levels(as.factor(foo$icluster)),
                   legend="right", pval=FALSE, palette=colors_icluster[c(2,1,3,4,5)], size=2);
plot$plot  <- plot$plot +  
  ggplot2::annotate("text", x = 1500, y = 0.9, label = label1, size = 7) + # x and y coordinates of the text
  ggplot2::annotate("text", x = 1500, y = 0.8, label = label2, size = 7) + # x and y coordinates of the text
  ggplot2::annotate("text", x = 1500, y = 0.7, label = label3, size = 7) + # + x and y coordinates of the text
  ggplot2::annotate("text", x = 1500, y = 0.6, label = label4, size = 7) # + x and y coordinates of the text

plot$plot <- ggpar(plot$plot, 
      font.main = c(20, "bold"),
      font.x = c(20, "bold"),
      font.y = c(20, "bold"),
      font.caption = c(20, "bold"), 
      font.legend = c(20, "bold"), 
      font.tickslab = c(20, "bold"))

plot;

#------------------------------------------------------------
# (3) icluster vs. categorical variables (mosaic plots)
#
annot_colors = list(
  time = c('Post'='grey20', 'Pre'='lightgrey'),
  mut = c('mut'='firebrick2', 'wt'='lightskyblue'),
  amp = c('amp'='orange3', 'wt'='lightskyblue'),
  del = c('wt'='lightskyblue', 'del'='darkblue'),
  subtype = c('Basal'='purple', 'Her2'='green', 'LumA'='yellow', 'LumB'='orange', 'Normal'='pink'),
  brca_status = c('0'='lightskyblue', '1'='firebrick2'),
  epg_branch = brewer.pal(4, 'Set1')[c(4, 3,2)],
  mutsig_status = rainbow(6),
  histology = rainbow(3),
  organ = rainbow(4)
)

# icluster vs. categorical variables
var = 'icluster';
mdata$time = mdata$pre_post;
res <- chisq.test(table(mdata$time, mdata[,var]));
title = paste(var, ': p=', round(res$p.value, 5), sep='');

#--------------------------------------
# mosaicplot: icluster vs. Pre/Post
#
par(mfrow=c(2,3), cex.main=1.5, cex.sub=1.5, cex.lab=1.5, font=2);
mosaicplot(table(mdata$time, mdata[,var])[c(2,1),], main=title, col=colors_icluster, las=1, cex.axis=2, border='grey');

# icluster vs. categorical variables
vars = c('time', 'organ', 'subtype_pam50');
colnames = c('time', 'organ', 'subtype');
for (i in 1:length(vars)){
  var = vars[i];
  colname = colnames[i];
  res <- chisq.test(table(mdata[,var], mdata[,'icluster']));
  title = paste(var, ': p=', round(res$p.value, 5), sep='');
  par(cex.main=1.5, cex.sub=1.5, cex.lab=1.5, font=2);
  if (colname %in% c('mut', 'amp')) { 
    y = ifelse(mdata[,var] == 1, colname, 'wt');
  } else {
    y = mdata[,var];
  }
  mosaicplot(table(mdata[,'icluster'], y), main=title, col=annot_colors[[colname]], ylab='', las=1, cex.axis=1.7, border='grey');
}

#------------------------------------------------------
# (4) icluster vs. continuous variables (boxplots)
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
comps = list(c('IC1', 'IC2'), c('IC1', 'IC3'), c('IC1', 'IC4'), c('IC1', 'IC5'));
orders = c('IC1', 'IC2', 'IC3', 'IC4', 'IC5');
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
  
  plist[[i]]<- ggboxplot(foo, 'icluster', 'value', xlab='', ylab=ylabel, title='', outlier.shape=NA, remove=NA,
                         facet.by='varlabel', scales='free_y', ncol=length(vars),
                         fill='icluster', palette=colors_icluster, order=orders) + 
    stat_compare_means(size=5, label='p.signif', comp=comps) + 
    theme(text = element_text(size = 12, face='bold')) + 
    scale_y_continuous(expand = expansion(mult = c(.1, .15))) +
    geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);
}
print(ggarrange(plotlist=plist[1:2], nrow=2, ncol=2, legend=F))
print(ggarrange(plotlist=plist[3], nrow=2, ncol=2, legend=F, widths=c(2,3)))

#-------------------------------------
# (5) icluster vs. genomic features
# grouped stacked barplots
#

# Mutation
vars = c('MUT.TP53', 'MUT.ESR1', 'MUT.RB1', 'MUT.GATA3');
for (i in 1:length(vars)){
  tmp <- melt(prop.table(table(mdata[,'icluster'], mdata[,vars[i]]), 1))
  tmp$Var3 <- vars[i]
  if (i == 1){
    foo <- tmp;
  } else {
    foo <- rbind(foo, tmp)
  }
}
foo$Var2 <- ifelse(foo$Var2 == 0, 'N', 'Y')
foo$Var3 <- gsub('MUT.', '', foo$Var3);
colnames(foo) <- c('icluster', 'mutation', 'proportion', 'gene')
foo$gene <- factor(foo$gene, levels = c('TP53', 'RB1', 'ESR1', 'GATA3'));

p1 <- ggplot(foo, aes(x = icluster, y = proportion, fill = mutation)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values=c('lightskyblue', 'firebrick2')) +
  facet_grid(~ gene) + ylab('% Samples') +
  theme(text = element_text(size = 15, face='bold')) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))

# Amplification
vars = c('AMP.CCND1', 'AMP.MYC', 'AMP.FGFR1');
for (i in 1:length(vars)){
  tmp <- melt(prop.table(table(mdata[,'icluster'], mdata[,vars[i]]), 1))
  tmp$Var3 <- vars[i]
  if (i == 1){
    foo <- tmp;
  } else {
    foo <- rbind(foo, tmp)
  }
}
foo$Var2 <- ifelse(foo$Var2 == 0, 'N', 'Y')
foo$Var3 <- gsub('AMP.', '', foo$Var3);
colnames(foo) <- c('icluster', 'amp', 'proportion', 'gene')
foo$gene <- factor(foo$gene, levels = c('MYC', 'CCND1', 'FGFR1'));

p2 <- ggplot(foo, aes(x = icluster, y = proportion, fill = amp)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values=c('lightskyblue', 'orange3')) +
  facet_grid(~ gene) + ylab('% Samples') +
  theme(text = element_text(size = 15, face='bold')) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))

print(ggarrange(p1, p2, nrow=3, ncol=2, widths=c(4,3)));

#--------------------------------------
# (6) barplot: icluster vs. Pre/Post
#

# color by IC1-5, use pattern to indicate Pre/Post
par(mfrow=c(2,2), mar=c(4,6,4,4));
tbl <- 100*prop.table(table(mdata$time, mdata$icluster), 1);
barplot(tbl[c(2,1),], beside=T, col=rep(colors_icluster, each=2), ylab='% Samples',
        density=c(30, 300), angle=45, cex.axis=1.5, cex.names=1.5, cex.lab=1.5, font=2,
        ylim=c(0,40))
legend("topright", legend=c('Pre', 'Post'), fill=rep(colors_icluster[1], 2), 
       angle=45, density=c(30, 300), bty='n', cex=1.5)

dev.off();

