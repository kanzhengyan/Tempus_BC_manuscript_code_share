
#=============================================
# (1) Load Tempus BC R data object
# perform routine data revision
#
tempus_bc <- function(){
  
  # load data object
  # rdatafile <- paste(datadir, 'TEMPUS.2023-02-09.RData', sep='/');
  rdatafile <- paste(datadir, 'TEMPUS.2023-05-27.RData', sep='/');
  rdata = get(load(rdatafile));
  sdata = rdata$Annotation$Sample;

  #------------------------------------
  # (1) Revise sample annotation
  #
  # (a) subcohort
  tmp = sdata$subcohorts;
  sdata$subcohorts <- ifelse(tmp %in% c('Pan_Cancer Patients with Biopsy Pre and Post-CDK4_6 (n=22)', 
                                        'Pan_Cancer Patients with Biopsy Pre and Post-CDK4_6 (n=6)'), '1-Paired',
                            ifelse(tmp %in% c('Pan_Cancer Patients with Post-CDK4_6 Biopsy Less Than 1yr from Progression (n=203)',
                                              'Pan_Cancer Patients with Post-CDK4_6 Biopsy Less Than 1yr from Progression (n=3)'), '2-Post', '3-Pre'));

  # (b) PFS status
  sdata$pfs_status <- ifelse(sdata$first_progression_post_cdk_censored_flag == 0, 'Censored',
                           ifelse(sdata$first_progression_post_cdk <= 180, 'PFS<=6m', 'PFS>6m'));

  # (c) Tissue site
  tmp = sdata$DNA_tissue_site_canonical_name;
  tmp[grepl('Lymph node', tmp)] = 'Lymph node';
  tmp[tmp %in% c('Liver', 'Breast', 'Lymph node') == F] = 'Others';
  sdata$organ = tmp;
  
  # (d) Treatment time point
  sdata$time <- ifelse(sdata$PrePostCDKi == 'Pre', 'BL', 'PD');
  sdata$time[sdata$PFEPrePost == 'During'] = 'OT';

  #-----------------------------------
  # (2) create patient data
  #
  foo <- unique(sdata[, c('patient_id', 'subcohorts', 'pfs_status', 'cdk_duration',	'cdk_duration_censored_flag',	'survival_post_cdk',	'cdk_survival_censored_flag',
                          'first_progression_post_cdk', 'first_progression_post_cdk_censored_flag', 'second_progression_post_cdk', 'second_progression_post_cdk_censored_flag'),]);
  tmp <- subset(rdata$Annotation$Patient, select=c('patient_id', 'subcohort_description', 'gender', 'age_at_pdx', 'race_concept_canonical_name', 'medication_concept_canonical_name'));
  pdata <- merge(tmp, foo, by='patient_id', sort=F);
  colnames(pdata)[4] = 'age';
  colnames(pdata)[5] = 'race';
  colnames(pdata)[6] = 'medication';
  
  #-----------------------------------------------
  # (3) Subset & rename sample attributes
  #
  tmp <- subset(sdata, select=c('ID', 'patient_id', 'subcohorts', 'pfs_status', 'PrePostCDKi', 'time',
                                  "RNASeq_analysis_id", "RNASeq_isolate_molecular_assay",           
                                  "DNA_analysis_id", "DNA_isolate_molecular_assay", "DNA_tissue_site_canonical_name", "organ", "DNA_contents_type", "DNA_procedure_type",                       
                                  "DNA_final_tumor_percentage", "DNA_match_type", 
                                  "MSI", "Tempus_Provided_TMB", "Tempus_Provided_TMB_percentile", "TMB", "Indel_load", "Tempus_est_tumor_purity_from_DNA",          
                                  "FACETS_tumor_purity", "ESTIMATE_tumor_purity", "HRD", "m_proliferative_index", "m_cytscore", "pam50",
                                  "biopsy_collection_time_from_index", "PFEPrePost", "time_biopsy_to_cdk_start", "time_biopsy_to_cdk_end",
                                  "first_progression_post_cdk", "first_progression_post_cdk_censored_flag",
                                  "received_palbociclib", "received_ribociclib", "received_abemaciclib",
                                  "fulvestrant_during_cdk", "letrozole_during_cdk", "anastrozole_during_cdk", "cdk_metastatic_lot",
                                  "exposure_before_biopsy", "exposure_adjuvant_endocrine", "exposure_adjuvant_chemo", "exposure_adjuvant_therapy",
                                  "exposure_met_prior_to_cdk_endocrine",  "exposure_met_prior_to_cdk_chemo",  "exposure_met_prior_to_cdk_therapy",
                                  "exposure_met_prior_to_biopsy_endocrine", "exposure_met_prior_to_biopsy_chemo", "exposure_met_prior_to_biopsy_therapy"));
  
  colnames(tmp) <- c('sample_id', 'patient_id', 'subcohorts', 'pfs_status', 'prepost', 'time',
                     "rnaseq_id", "rnaseq_molecular_assay",  "dnaseq_id", "dnaseq_molecular_assay", 
                     "tissue_site",  "organ", "contents_type", "procedure_type", "tumor_purity", "dna_match_type", 
                     "msi", "tmb_tempus", "tmb_percentile", "tmb", "indel_load", "tumor_purity_tempus", 
                     "tumor_purity_facets", "tumor_purity_estimate", "hrd_score", "proliferative_index", "cytscore", "subtype_pam50",
                     "biopsy_collection_time", "prepost_pfe", "time_biopsy_to_cdki_start", "time_biopsy_to_cdki_end", "pfs", "pfs_event",
                     "received_palbo", "received_ribo", "received_abema",
                     "fulvestrant_cdki", "letrozole_cdki", "anastrozole_cdki", "cdki_met_line",
                     "exposure_before_biopsy", "exposure_adjuvant_endocrine", "exposure_adjuvant_chemo", "exposure_adjuvant",
                     "exposure_met_before_cdki_endocrine",  "exposure_met_before_cdki_chemo",  "exposure_met_before_cdki",
                     "exposure_met_before_biopsy_endocrine", "exposure_met_before_biopsy_chemo", "exposure_met_before_biopsy");
  sdata <- tmp;
  
  #-----------------------------------------------  
  # (4) BRCA1/2 pathogenic mutations
  #
  germlinedata <- rdata$Analysis$germline_matrix_8gene;
  midx <- match(sdata$sample_id, rownames(germlinedata))
  germlinedata <- germlinedata[midx,];
  colnames(germlinedata) <- tolower(colnames(germlinedata));
  tmp <- subset(germlinedata, select=c('brca_somatic_germline', 'brca2_germline', 'chek2_germline', 'atm_germline'));
  sdata <- data.frame(sdata, tmp, ddr_germline=as.numeric(apply(tmp, 1, sum) > 0));
  
  #--------------------------------------
  # (5) gene expression data filtering
  #
  expdata_combat = rdata$Expression$log2TPM_corrected;
  means = apply(expdata_combat, 1, mean);
  sds = apply(expdata_combat, 1, sd);
  gidx = (means >= quantile(means, prob=.25, na.rm=T) & sds > quantile(sds, prob=.25, na.rm=T)) | names(means) %in% c('CDK4', 'CDK6');
  expdata = expdata_combat[gidx,];    
  
  # reorder expression data by sample id
  midx = match(sdata$sample_id, colnames(expdata));
  expdata = expdata[,midx];

  # prepare output
  rdata$sdata <- sdata;
  rdata$pdata <- pdata;
  rdata$expdata <- expdata;

  return(rdata)
}

#===================================================
# (2) Extract & annotate & extract GSVA scores
#
tempus_bc_gsva<-function(rdata) {
  
  gsvadata = rdata$Analysis$GSVA;
  
  # Load geneset annotation
  midx <- match(rownames(gsvadata), rdata$Analysis$GSVA_annotation$pathway_name)
  gsannot = rdata$Analysis$GSVA_annotation[midx,];
  colnames(gsannot) = c('geneset', 'collection', 'description');

  # extract GSVA scores
  sdata = rdata$sdata;
  midx = match(sdata$sample_id, colnames(gsvadata));
  gidx = gsannot$collection %in% c('bindea', 'hallmark gene sets', 'C2 CP: REACTOME', 'C2 CP: KEGG', 'C2 CP: WikiPathways', 
                                           'C2 CP: BIOCARTA', 'C2 CP: Pathway Interaction Database', 'C3 TFT: GTRD subset of TFT', 'Pfizer_custom',
                                           'IO_group_curation', 'STING_working_group', 'C2 CGP: chemical and genetic perturbations');
  gsva = gsvadata[gidx,midx];
  gsannot = gsannot[gidx,];
  
  # parse collection names from geneset names
  #cl = unlist(lapply(rownames(gsva), function(x) {rr <- strsplit(x, "\\.")[[1]][1]}));
  #collections = data.frame(geneset=rownames(gsva), collection=cl);
  
  # prepare output
  gsvaobj <- list();
  gsvaobj[['scores']] <- gsva;
  gsvaobj[['annotation']] <- gsannot;
  return(gsvaobj);
}

#-------------------------
# extract CSE GSVA scores
#
tempus_bc_cse_gsva<-function(rdata) {
  
  gsvadata = rdata$Analysis$GSVA_BayesPrism_cancer;
  
  # Load geneset annotation
  midx <- match(rownames(gsvadata), rdata$Analysis$GSVA_annotation$pathway_name)
  gsannot = rdata$Analysis$GSVA_annotation[midx,];
  colnames(gsannot) = c('geneset', 'collection', 'description');
  
  # extract GSVA scores
  sdata = rdata$sdata;
  midx = match(sdata$sample_id, colnames(gsvadata));
  gidx = gsannot$collection %in% c('bindea', 'hallmark gene sets', 'C2 CP: REACTOME', 'C2 CP: KEGG', 'C2 CP: WikiPathways', 
                                   'C2 CP: BIOCARTA', 'C2 CP: Pathway Interaction Database', 'C3 TFT: GTRD subset of TFT', 'Pfizer_custom',
                                   'IO_group_curation', 'STING_working_group', 'C2 CGP: chemical and genetic perturbations');
  gsva = gsvadata[gidx,midx];
  gsannot = gsannot[gidx,];
  
  # prepare output
  gsvaobj <- list();
  gsvaobj[['scores']] <- gsva;
  gsvaobj[['annotation']] <- gsannot;
  return(gsvaobj);
}

#===================================================
# (3) gene expression signature (GSVA) 
# characterization analysis - single
#
gsva_characterize_analysis <- function(outdir, gsname, sdata, gsvadata){
  
  #--------------------------
  # output file names
  #
  plotfile <<- paste(outdir, paste('tempus_gsva_characterize', gsname, 'pdf', sep='.'), sep='/');
  
  #--------------------------------
  # gene expression analysis
  #
  pdf(file=plotfile, onefile=T, height=8, width=10, pointsize=12);
  
  # use adjusted expression data
  gidx = (is.element(rownames(gsvadata), gsname));
  expdata = as.numeric(gsvadata[gidx,]);
  
  # color labels
  colors_organ = rainbow(4);
  colors_prepost = palette('Set1')[1:2]
  colors_pam50 = c("purple", "green", "yellow", "orange", "pink");
  sidx_pfs = grepl('PFS', sdata$pfs_status);
  sidx_lum = sdata$subtype_pam50 %in% c('LumA', 'LumB');
  
  #-------------------------------------------------------------
  # (1) Examine gene expression changes over time 'boxplot'
  # facet by PFS status, PAM50
  #
  comps <- list(c("Pre", "Post"));
  mdata = data.frame(gene=expdata, sdata);
  p1 = ggboxplot(mdata, "prepost", "gene", xlab='', ylab='GSVA', title=gsname, outlier.shape=NA,
                 fill='prepost', palette=colors_prepost,  order=c('Pre', 'Post')) +
    stat_compare_means(comparison=comps, size=4, label='p.format') +
    theme(text = element_text(size = 10, face='bold')) +
    geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);
  
  mdata2 = mdata[sidx_pfs,];
  p2 = ggboxplot(mdata2, "prepost", "gene", xlab='', ylab='GSVA', title=gsname, outlier.shape=NA, facet='pfs_status',
                 fill='prepost', palette=colors_prepost,  order=c('Pre', 'Post')) +
    stat_compare_means(comparison=comps, size=4, label='p.format') +
    theme(text = element_text(size = 10, face='bold'), legend.position = "none") +
    geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);
  
  sidx_subtype = sdata$subtype_pam50 %in% c('LumA', 'LumB', 'Her2');
  mdata3 = mdata[sidx_subtype,]
  p3 = ggboxplot(mdata3, "prepost", "gene", xlab='', ylab='GSVA', title=gsname, outlier.shape=NA, facet='subtype_pam50',
                 fill='prepost', palette=colors_prepost,  order=c('Pre', 'Post')) +
    stat_compare_means(comparison=comps, size=4, label='p.format') +
    theme(text = element_text(size = 10, face='bold'), legend.position = "none") +
    geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);
  
  # Compare gene expression vs. PAM50 subtype, facet by time
  comps2 = list(c('LumA', 'LumB'));
  sidx_lum = sdata$subtype_pam50 %in% c('LumA', 'LumB');
  mdata4 = data.frame(gene=expdata[sidx_lum], sdata[sidx_lum,]);
  p4 = ggboxplot(mdata4, "subtype_pam50", "gene", xlab='', ylab='GSVA', title=gsname, outlier.shape=NA, facet='prepost',
                 fill='subtype_pam50', palette=c("yellow", "orange"),  order=c('LumA', 'LumB')) +
    stat_compare_means(comparison=comps2, size=4, label='p.format') +
    theme(text = element_text(size = 10, face='bold'), legend.position = "none") +
    geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);
  
  print(ggarrange(p1, p2, p3, p4, nrow=2, ncol=2));
  
  #---------------------------------------------------------------------
  # (2) Longitudinal gene expression changes over time
  # for only paired samples (1-Paired)
  # paired 'boxplot', p-value from paired Wilcoxon test
  #

  # data for only sample pairs, determine baseline subtype
  mdata_paired = mdata[mdata$subcohorts == '1-Paired',]
  tmp <- mdata_paired[mdata_paired$prepost == 'Pre', c('patient_id', 'subtype_pam50')];
  colnames(tmp) = c('patient_id', 'subtype_pre');
  mdata_paired <- merge(mdata_paired, tmp, by='patient_id', sort=F, all.x=T);

  p1 <- ggpaired(mdata_paired, x='prepost', y='gene', title=gsname, fill='prepost', id='patient_id', xlab='',
                 line.color='gray', line.size=.4, palette='Set1', order=c('Pre', 'Post')) +
    stat_compare_means(paired=TRUE, size=4) +
    theme(text = element_text(size = 10, face='bold'), legend.position = "none");
  
  mdata_paired2 = mdata_paired[mdata_paired$subtype_pre %in% c('LumA', 'LumB'),];
  p2 <- ggpaired(mdata_paired2, x='prepost', y='gene', title=gsname, facet='subtype_pre', fill='prepost', id='patient_id', xlab='',
                 line.color='gray', line.size=.4, palette='Set1', order=c('Pre', 'Post')) +
    stat_compare_means(paired=TRUE, size=4) +
    theme(text = element_text(size = 10, face='bold'), legend.position = "none");
  
  print(ggarrange(p1, p2, nrow=2, ncol=2, widths=c(1,2)));
  
  #---------------------------------------------------------
  # 5 - compare expression vs. key covariates 
  # (1) tissue origins
  # (2) tissue type - FFPE/FF
  # (3) PAM50 subtype
  #
  comps1 = list(c('Breast', 'Liver'), c('Breast', 'Lymph node'));
  p1 = ggboxplot(mdata, "organ", "gene", xlab='', ylab='GSVA', title=gsname, outlier.shape=NA,
                 fill='organ', palette=colors_organ) +
    stat_compare_means(comparison=comps1, size=4, label='p.format') +
    theme(text = element_text(size = 10, face='bold'), legend.position = "none") +
    geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);

  comps2 = list(c('LumA', 'LumB'), c('LumA', 'Her2'));
  p2 = ggboxplot(mdata, "subtype_pam50", "gene", xlab='', ylab='GSVA', title=gsname, outlier.shape=NA,
                 fill='subtype_pam50', palette=colors_pam50, order=c('Basal', 'Her2', 'LumA', 'LumB', 'Normal')) +
    stat_compare_means(comp=comps2, size=4, label='p.format') +
    theme(text = element_text(size = 10, face='bold'), legend.position = "none") +
    geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);
  
  print(ggarrange(p1, p2, nrow=2, ncol=2));
  
  dev.off();
}

#===================================================
# (3a) single variable (continuous) 
# characterization analysis
# type: gsva, exp, feature
#
var_characterize_analysis <- function(outdir, var, sdata, vardata, vartype){
  
  library('RColorBrewer')
  
  #--------------------------
  # output file names
  #
  plotfile <<- paste(outdir, paste(paste('tempus', vartype, 'characterize', sep='_'), var, 'pdf', sep='.'), sep='/');
  pdf(file=plotfile, onefile=T, height=8, width=10, pointsize=12);
  
  # extract variable data
  gidx = (is.element(rownames(vardata), var));
  expdata = as.numeric(vardata[gidx,]);

  # variable label
  varlabel <- ifelse(vartype == 'gsva', 'GSVA',
                     ifelse(vartype == 'exp', 'Log2TPM', 
                            ifelse(vartype == 'cnv', 'Log2CN', 'Value')));

  # color labels
  colors_organ = rainbow(4);
  colors_prepost = palette.colors(palette='Set1')[c(2,1)];
  colors_time = c(colors_prepost[1], "grey", colors_prepost[2]);
  colors_pam50 = c("purple", "green", "yellow", "orange", "pink");
  sidx_pfs = grepl('PFS', sdata$pfs_status);
  sidx_lum = sdata$subtype_pam50 %in% c('LumA', 'LumB');
  
  #-------------------------------------------------------------
  # (1) Examine variable changes Pre vs. Post 'boxplot'
  # facet by PFS status, PAM50
  #
  comps <- list(c("Pre", "Post"));
  mdata = data.frame(var=expdata, sdata);
  p1 = ggboxplot(mdata, "prepost", "var", xlab='', ylab=varlabel, title=var, outlier.shape=NA,
                 fill='prepost', palette=colors_prepost,  order=c('Pre', 'Post')) +
    stat_compare_means(comparison=comps, size=3, label='p.format') +
    theme(text = element_text(size = 10, face='bold')) +
    geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);
  
  mdata2 = mdata[sidx_pfs,];
  p2 = ggboxplot(mdata2, "prepost", "var", xlab='', ylab=varlabel, title=var, outlier.shape=NA, facet='pfs_status',
                 fill='prepost', palette=colors_prepost,  order=c('Pre', 'Post')) +
    stat_compare_means(comparison=comps, size=3, label='p.format') +
    theme(text = element_text(size = 10, face='bold'), legend.position = "none") +
    geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);
  
  sidx_subtype = sdata$subtype_pam50 %in% c('LumA', 'LumB', 'Her2');
  mdata3 = mdata[sidx_subtype,]
  p3 = ggboxplot(mdata3, "prepost", "var", xlab='', ylab=varlabel, title=var, outlier.shape=NA, facet='subtype_pam50',
                 fill='prepost', palette=colors_prepost,  order=c('Pre', 'Post')) +
    stat_compare_means(comparison=comps, size=3, label='p.format') +
    theme(text = element_text(size = 10, face='bold'), legend.position = "none") +
    geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);
  
  # Compare variable vs. PAM50 subtype, facet by time
  comps2 = list(c('LumA', 'LumB'));
  sidx_lum = sdata$subtype_pam50 %in% c('LumA', 'LumB');
  mdata4 = data.frame(var=expdata[sidx_lum], sdata[sidx_lum,]);
  p4 = ggboxplot(mdata4, "subtype_pam50", "var", xlab='', ylab=varlabel, title=var, outlier.shape=NA, facet='prepost',
                 fill='subtype_pam50', palette=c("yellow", "orange"),  order=c('LumA', 'LumB')) +
    stat_compare_means(comparison=comps2, size=3, label='p.format') +
    theme(text = element_text(size = 10, face='bold'), legend.position = "none") +
    geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);
 
  print(ggarrange(p1, nrow=2, ncol=4)); 
  print(ggarrange(p1, p2, p3, p4, nrow=2, ncol=2));
  
  #-------------------------------------------------------------
  # (2) Examine variable changes over Time 'boxplot'
  # facet by PFS status, PAM50
  #
  comps <- list(c("BL", "PD"), c("BL", "OT"), c("OT", "PD"));
  p1 = ggboxplot(mdata, "time", "var", xlab='', ylab=varlabel, title=var, outlier.shape=NA,
                 fill='time', palette=colors_time) +
    stat_compare_means(comparison=comps, size=3, label='p.format') +
    theme(text = element_text(size = 10, face='bold')) +
    geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);

  p2 = ggboxplot(mdata2, "time", "var", xlab='', ylab=varlabel, title=var, outlier.shape=NA, facet='pfs_status',
                 fill='time', palette=colors_time) +
    stat_compare_means(comparison=comps, size=3, label='p.format') +
    theme(text = element_text(size = 10, face='bold'), legend.position = "none") +
    geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);

  p3 = ggboxplot(mdata3, "time", "var", xlab='', ylab=varlabel, title=var, outlier.shape=NA, facet='subtype_pam50',
                 fill='time', palette=colors_time) +
    stat_compare_means(comparison=comps, size=3, label='p.format') +
    theme(text = element_text(size = 10, face='bold'), legend.position = "none") +
    geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);
  
  # Compare variable vs. PAM50 subtype, facet by time
  p4 = ggboxplot(mdata4, "subtype_pam50", "var", xlab='', ylab=varlabel, title=var, outlier.shape=NA, facet='time',
                 fill='subtype_pam50', palette=c("yellow", "orange"),  order=c('LumA', 'LumB')) +
    stat_compare_means(comparison=comps2, size=3, label='p.format') +
    theme(text = element_text(size = 10, face='bold'), legend.position = "none") +
    geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);
    
  print(ggarrange(p1, p2, p3, p4, nrow=2, ncol=2));
  
  #---------------------------------------------------------------------
  # (3) Longitudinal variable changes Pre vs. Post
  # for only paired samples (1-Paired)
  # paired 'boxplot', p-value from paired Wilcoxon test
  #
  
  # data for only sample pairs, determine baseline subtype
  mdata_paired = mdata[mdata$subcohorts == '1-Paired',]
  tmp <- mdata_paired[mdata_paired$prepost == 'Pre', c('patient_id', 'subtype_pam50')];
  colnames(tmp) = c('patient_id', 'subtype_pre');
  tmp$subtype_pre = paste('Pre:', tmp$subtype_pre);
  mdata_paired <- merge(mdata_paired, tmp, by='patient_id', sort=F, all.x=T);
  
  p1 <- ggpaired(mdata_paired, x='prepost', y='var', title=var, fill='prepost', id='patient_id', xlab='',
                 line.color='gray', line.size=.4, palette='Set1', order=c('Pre', 'Post')) +
    stat_compare_means(paired=TRUE, size=3) +
    theme(text = element_text(size = 10, face='bold'), legend.position = "none");
  
  mdata_paired2 = mdata_paired[mdata_paired$subtype_pre %in% c('Pre: LumA', 'Pre: LumB'),];
  p2 <- ggpaired(mdata_paired2, x='prepost', y='var', title=var, 
                 facet='subtype_pre', fill='prepost', id='patient_id', xlab='',
                 line.color='gray', line.size=.4, palette='Set1', order=c('Pre', 'Post')) +
    stat_compare_means(paired=TRUE, size=3) +
    theme(text = element_text(size = 10, face='bold'), legend.position = "none");
  
  #---------------------------------------------------------------------
  # (4) Longitudinal variable changes over time
  # for only paired samples (1-Paired)
  # paired 'boxplot', p-value from paired Wilcoxon test
  #
  p3 <- ggpaired(mdata_paired, x='time', y='var', title=var, fill='time', id='patient_id', xlab='',
                 line.color='gray', line.size=.4, palette=colors_time, order=c('BL', 'OT', 'PD')) +
    stat_compare_means(paired=TRUE, size=3) +
    theme(text = element_text(size = 10, face='bold'), legend.position = "none");
  
  p4 <- ggpaired(mdata_paired2, x='time', y='var', title=var, facet='subtype_pre', fill='time', id='patient_id', xlab='',
                 line.color='gray', line.size=.4, palette=colors_time, order=c('BL', 'OT', 'PD')) +
    stat_compare_means(paired=TRUE, size=3) +
    theme(text = element_text(size = 10, face='bold'), legend.position = "none");
  
  print(ggarrange(p1, p2, p3, p4, nrow=2, ncol=2, widths=c(1,2)));

  #---------------------------------------------------------
  # 5 - compare variable vs. key covariates 
  # (1) tissue origins
  # (2) tissue type - FFPE/FF
  # (3) PAM50 subtype
  # (4) icluster
  # (5) icluster+PAM50
  #
  comps1 = list(c('Breast', 'Liver'), c('Breast', 'Lymph node'));
  p1 = ggboxplot(mdata, "organ", "var", xlab='', ylab=varlabel, title=var, outlier.shape=NA,
                 fill='organ', palette=colors_organ) +
    stat_compare_means(comparison=comps1, size=3, label='p.format') +
    theme(text = element_text(size = 10, face='bold'), legend.position = "none") +
    geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);
  
  comps2 = list(c('LumA', 'LumB'), c('LumA', 'Her2'), c('LumA', 'Basal'), c('Basal', 'Her2'));
  p2 = ggboxplot(mdata, "subtype_pam50", "var", xlab='', ylab=varlabel, title=var, outlier.shape=NA,
                 fill='subtype_pam50', palette=colors_pam50, order=c('Basal', 'Her2', 'LumA', 'LumB', 'Normal')) +
    stat_compare_means(comp=comps2, size=3, label='p.format') +
    theme(text = element_text(size = 10, face='bold'), legend.position = "none") +
    geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);

  comps3 = list(c('IC1', 'IC2'), c('IC1', 'IC3'), c('IC1', 'IC4'), c('IC2', 'IC3'));
  p3 = ggboxplot(mdata[!is.na(mdata$icluster),], "icluster", "var", xlab='', ylab=varlabel, title=var, outlier.shape=NA,
                 fill='icluster', palette='Set1', order=c('IC1', 'IC2', 'IC3', 'IC4', 'IC5')) +
    stat_compare_means(comp=comps3, size=3, label='p.format') +
    theme(text = element_text(size = 10, face='bold'), legend.position = "none") +
    geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);
  
  #-------------------------
  # icluster + PAM50
  #

  # color subgroup by PAM50
  tmp <- matrix(unlist(strsplit(levels(as.factor(mdata$ic_pam50)), ":")), ncol=2, byrow=TRUE)[,2]
  collabs <- colors_pam50[as.numeric(as.factor(tmp))];
  names(collabs) <- levels(mdata$ic_pam50)
    
  comps4 = list(c('IC1:LumB', 'IC1:Her2'), c('IC1:Basal', 'IC1:LumB'), c('IC1:Basal', 'IC1:Her2'), 
               c('IC2:LumB', 'IC1:LumB'), c('IC2:LumB', 'IC3:LumB'), c('IC1:LumB', 'IC3:LumB'));
  p4 = ggboxplot(mdata[!is.na(mdata$ic_pam50),], "ic_pam50", "var", xlab='', ylab=varlabel, title=var, outlier.shape=NA,
                 fill='ic_pam50', palette=collabs,
                 order=sort(levels(as.factor(mdata$ic_pam50)))) +
    stat_compare_means(comp=comps4, size=3, label='p.format') +
    theme(text = element_text(size = 10, face='bold'), legend.position = "none") +
    geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3) + rotate_x_text(45)
  
  print(ggarrange(p1, p2, p3, p4, nrow=2, ncol=2));
  
  #---------------------------------------------------------
  # (6) compare variable vs. key genomic events (MUT, AMP)
  # facet by Pre/Post
  #
  plist <- list();
  genomics <- c('MUT.ESR1', 'MUT.RB1', 'MUT.TP53', 'MUT.GATA3', 
                'MUT.PTEN', 'MUT.PIK3CA', 'MUT.AKT1', 
                'AMP.KAT6A', 'AMP.CCND1', 'AMP.FGFR1', 'AMP.MYC',
                'DEL.CDKN2A');

  # color labels
  annot_colors = list(
    mut = c('0'='deepskyblue2', '1'='firebrick2'),
    amp = c('0'='deepskyblue2', '1'='orange3'),
    del = c('0'='deepskyblue2', '1'='darkblue')
  )
  mdata$prepost <- factor(mdata$prepost, levels=c('Pre', 'Post'))
  comps = list(c('0', '1'));
  for (i in 1:length(genomics)){
    genomic <- genomics[i];
    if (grepl('MUT', genomic)) {
      colors <- annot_colors$mut;
    } else if (grepl('AMP', genomic)) {
      colors <- annot_colors$amp;
    } else {
      colors <- annot_colors$del;
    }
    plist[[i]] <- ggboxplot(mdata, genomic, "var", xlab=genomic, ylab=varlabel, title=var, outlier.shape=NA, remove=NA,
                            fill=genomic, palette=colors, facet.by='prepost', order=c('0', '1')) + 
      stat_compare_means(size=3, label='p.format', comp=comps) + 
      theme(text = element_text(size = 10, face='bold')) + 
      geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);
  }
  print(ggarrange(plotlist=plist, nrow=2, ncol=3));
  dev.off();
}

#===================================================
# function to turn numbers into color gradient
#
mycolorscale = function(x, color_low, color_high){
  mypal <- colorRampPalette( c(color_low, color_high) )(length(x));
  order = findInterval(x, sort(x));
  return(mypal[order]);
}

#===================================================
# (4) select molecular features (continuous)
#
get_molecular_feature_continuous <- function(rdata) {
  
  # select key features from sample metadata
  sdata = rdata$sdata;
  vars = c('sample_id', 'tmb_tempus', 'tmb_percentile', 'indel_load', 'tumor_purity', 'tumor_purity_tempus', 'hrd_score', 'proliferative_index', 'cytscore',
           'time_biopsy_to_cdki_start', 'time_biopsy_to_cdki_end');
  
  # a) add "infiltration" features
  midx <- match(sdata$sample_id, rownames(rdata$Analysis$infiltration))
  tmp1 <- rdata$Analysis$infiltration[midx,-c(1:2)];
  
  # b) add "ESTIMATE" features
  midx <- match(sdata$sample_id, rownames(rdata$Analysis$ESTIMATE))
  tmp2 <- rdata$Analysis$ESTIMATE[midx,];
  colnames(tmp2) = c('estimate_stromalscore', 'estimate_immunescore', 'estimate_score', 'estimate_purity');
  
  # c) add "FACETS" genomic scar features
  gscars <- rdata$Analysis$FACETS_stat[,-1];
  colnames(gscars) = c('facets_purity', 'facets_ploidy', 'facets_hrdloh', 'facets_tai', 'facets_lst', 'facets_hrdsum')
  midx <- match(sdata$sample_id, rownames(gscars));
  tmp3 <- gscars[midx,]
    
  # d) add PAM50 subtype scores
  pam50 <- subset(rdata$Analysis$PAM50, select=c('pam50.cor.Basal', 'pam50.cor.Her2', 'pam50.cor.LumA', 'pam50.cor.LumB', 'pam50.cor.Normal'))
  colnames(pam50) <- gsub('cor.', '', colnames(pam50))
  midx <- match(sdata$sample_id, rownames(pam50));
  tmp4 <- pam50[midx,];
  
  df <- data.frame(sdata[,vars], tmp1, tmp2, tmp3, tmp4);
  
  # e) add "Mutational_signatures"
  mutsigs <- rdata$Analysis$Mutational_signatures.normalized
  colnames(mutsigs) <- gsub('exp_sig', 'mutsig', colnames(mutsigs));
  mutcounts <- apply(mutsigs, 2, function(x)return(100*sum(x>0.1)/length(x)))
  signames <- names(which(mutcounts > 10));
  mutsigs2 <- mutsigs[,signames];
  fdata <- merge(df, data.frame(sample_id=rownames(mutsigs2), mutsigs2), by='sample_id', sort=F, all.x=T);
  
  # f) add "de novo NMF factors"
  nmf <- t(rdata$Analysis$NMF);
  colnames(nmf) = paste('NMF', colnames(nmf), sep='_')
  fdata <- merge(fdata, data.frame(sample_id=rownames(nmf), nmf), by='sample_id', all.x=T, sort=F);

  # g) add "PALOMA NMF factors"
  nmf2 <- t(rdata$Analysis$projection_to_paloma3_NMF);
  colnames(nmf2) = gsub('Paloma3', 'PALOMA3', colnames(nmf2))
  fdata <- merge(fdata, data.frame(sample_id=rownames(nmf2), nmf2), by='sample_id', all.x=T, sort=F);
  
  rownames(fdata) = fdata$sample_id;
  
  return(fdata);
}

#===================================================
# (5) single variable (discrete) 
# characterization analysis
#
discrete_feature_analysis <- function(outdir, var, sdata, vardata){
  
  #--------------------------
  # output file names
  #
  plotfile <<- paste(outdir, paste('tempus_discrete_feature_analysis', var, 'pdf', sep='.'), sep='/');
  pdf(file=plotfile, onefile=T, height=8, width=10, pointsize=12);
  
  # extract variable data
  gidx = (is.element(rownames(vardata), var));
  vdata = vardata[gidx,];
  
  # color labels
  colors_organ = rainbow(4);
  colors_prepost = palette('Set1')[1:2]
  colors_pfs = palette('Set2')[1:2];
  colors_pam50 = c("green", "yellow", "orange");
  sidx_pre = sdata$prepost == 'Pre';
  sidx_pfs = grepl('PFS', sdata$pfs_status);
  sidx_lum = sdata$subtype_pam50 %in% c('LumA', 'LumB');
  sidx_subtype = sdata$subtype_pam50 %in% c('LumA', 'LumB', 'Her2');
  
  #-------------------------------------------------------------
  # (1) Examine variable vs. key covariates
  # Pre/Post, PFS, PAM50
  #
  par(mfrow=c(2,2));
  mdata = data.frame(var=vdata, sdata);
  tbl <- table(mdata[,c('var', 'prepost')]);
  res <- chisq.test(tbl);
  title <- paste('log10p=', round(-log10(res$p.value), 3), sep='');
  mosaicplot(tbl, col=colors_prepost, main=title, cex=1, xlab=var, las=1);

  # Pre: PFS status
  mdata2 = mdata[sidx_pfs & sidx_pre,]
  tbl <- table(mdata2[,c('var', 'pfs_status')]);
  res <- chisq.test(tbl);
  title <- paste('Pre: log10p=', round(-log10(res$p.value), 3), sep='');
  mosaicplot(tbl, col=colors_pfs, main=title, cex=1, xlab=var, las=1);
  
  # PAM50
  mdata3 = mdata[sidx_subtype,]  
  tbl <- table(mdata3[,c('var', 'subtype_pam50')]);
  res <- chisq.test(tbl);
  title <- paste('log10p=', round(-log10(res$p.value), 3), sep='');
  mosaicplot(tbl, col=colors_pam50, main=title, cex=1, xlab=var, las=1);
  
  # Organ/tissue
  tbl <- table(mdata[,c('var', 'organ')]);
  res <- chisq.test(tbl);
  title <- paste('log10p=', round(-log10(res$p.value), 3), sep='');
  mosaicplot(tbl, col=colors_organ, main=title, cex=1, xlab=var, las=1);
  
  dev.off();
}

#===================================================
# (3a) single variable (continuous) 
# characterization analysis
# type: gsva, exp, feature
#
gene_cnv_analysis <- function(outdir, var, sdata, cnvdata, expdata){
  
  #--------------------------
  # output file names
  #
  plotfile <<- paste(outdir, paste('tempus_gene_cnv_analysis', gene, 'pdf', sep='.'), sep='/');
  pdf(file=plotfile, onefile=T, height=8, width=10, pointsize=12);
  
  # extract variable data
  gidx = (is.element(rownames(cnvdata), var));
  cn = as.numeric(cnvdata[gidx,]);
  exp = as.numeric(expdata[gidx,]);
  
  # variable label
  varlabel <- 'Log2CN-ratio';
  
  # color labels
  colors_organ = rainbow(4);
  colors_prepost = palette('Set1')[1:2]
  sdata$prepost <- factor(sdata$prepost, levels=c('Pre', 'Post'))
  collabs_prepost = colors_prepost[as.numeric(as.factor(sdata$prepost))];
  colors_pam50 = c("purple", "green", "yellow", "orange", "pink");
  collabs_pam50 = colors_pam50[as.numeric(as.factor(sdata$subtype_pam50))];
  sidx_pfs = grepl('PFS', sdata$pfs_status);
  sidx_lum = sdata$subtype_pam50 %in% c('LumA', 'LumB');
  
  #-------------------------------------------------------------
  # (1) Examine variable changes over time 'boxplot'
  # facet by PFS status, PAM50
  #
  comps <- list(c("Pre", "Post"));
  mdata = data.frame(var=cn, exp=exp, sdata);
  p1 = ggboxplot(mdata, "prepost", "var", xlab='', ylab=varlabel, title=var, outlier.shape=NA,
                 fill='prepost', palette=colors_prepost,  order=c('Pre', 'Post')) +
    stat_compare_means(comparison=comps, size=4, label='p.format') +
    theme(text = element_text(size = 10, face='bold')) +
    geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);
  
  mdata2 = mdata[sidx_pfs,];
  p2 = ggboxplot(mdata2, "prepost", "var", xlab='', ylab=varlabel, title=var, outlier.shape=NA, facet='pfs_status',
                 fill='prepost', palette=colors_prepost,  order=c('Pre', 'Post')) +
    stat_compare_means(comparison=comps, size=4, label='p.format') +
    theme(text = element_text(size = 10, face='bold'), legend.position = "none") +
    geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);
  
  sidx_subtype = sdata$subtype_pam50 %in% c('LumA', 'LumB', 'Her2');
  mdata3 = mdata[sidx_subtype,]
  p3 = ggboxplot(mdata3, "prepost", "var", xlab='', ylab=varlabel, title=var, outlier.shape=NA, facet='subtype_pam50',
                 fill='prepost', palette=colors_prepost,  order=c('Pre', 'Post')) +
    stat_compare_means(comparison=comps, size=4, label='p.format') +
    theme(text = element_text(size = 10, face='bold'), legend.position = "none") +
    geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);
  
  # Compare variable vs. PAM50 subtype, facet by time
  comps2 = list(c('LumA', 'LumB'));
  sidx_lum = sdata$subtype_pam50 %in% c('LumA', 'LumB');
  mdata4 = data.frame(var=cn[sidx_lum], sdata[sidx_lum,]);
  p4 = ggboxplot(mdata4, "subtype_pam50", "var", xlab='', ylab=varlabel, title=var, outlier.shape=NA, facet='prepost',
                 fill='subtype_pam50', palette=c("yellow", "orange"),  order=c('LumA', 'LumB')) +
    stat_compare_means(comparison=comps2, size=4, label='p.format') +
    theme(text = element_text(size = 10, face='bold'), legend.position = "none") +
    geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);
  
  print(ggarrange(p1, p2, p3, p4, nrow=2, ncol=2));
  
  #---------------------------------------------------------------------
  # (2) Longitudinal variable changes over time
  # for only paired samples (1-Paired)
  # paired 'boxplot', p-value from paired Wilcoxon test
  #
  
  # data for only sample pairs, determine baseline subtype
  mdata_paired = mdata[mdata$subcohorts == '1-Paired',]
  tmp <- mdata_paired[mdata_paired$prepost == 'Pre', c('patient_id', 'subtype_pam50')];
  colnames(tmp) = c('patient_id', 'subtype_pre');
  mdata_paired <- merge(mdata_paired, tmp, by='patient_id', sort=F, all.x=T);
  
  p1 <- ggpaired(mdata_paired, x='prepost', y='var', title=var, fill='prepost', id='patient_id', xlab='',
                 line.color='gray', line.size=.4, palette='Set1', order=c('Pre', 'Post')) +
    stat_compare_means(paired=TRUE, size=4) +
    theme(text = element_text(size = 10, face='bold'), legend.position = "none");
  
  mdata_paired2 = mdata_paired[mdata_paired$subtype_pre %in% c('LumA', 'LumB'),];
  p2 <- ggpaired(mdata_paired2, x='prepost', y='var', title=var, facet='subtype_pre', fill='prepost', id='patient_id', xlab='',
                 line.color='gray', line.size=.4, palette='Set1', order=c('Pre', 'Post')) +
    stat_compare_means(paired=TRUE, size=4) +
    theme(text = element_text(size = 10, face='bold'), legend.position = "none");
  
  print(ggarrange(p1, p2, nrow=2, ncol=2, widths=c(1,2)));
  
  #---------------------------------------------------------
  # 5 - compare variable vs. key covariates 
  # (1) tissue origins
  # (2) tissue type - FFPE/FF
  # (3) PAM50 subtype
  #
  comps1 = list(c('Breast', 'Liver'), c('Breast', 'Lymph node'));
  p1 = ggboxplot(mdata, "organ", "var", xlab='', ylab=varlabel, title=var, outlier.shape=NA,
                 fill='organ', palette=colors_organ) +
    stat_compare_means(comparison=comps1, size=4, label='p.format') +
    theme(text = element_text(size = 10, face='bold'), legend.position = "none") +
    geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);
  
  comps2 = list(c('LumA', 'LumB'), c('LumA', 'Her2'));
  p2 = ggboxplot(mdata, "subtype_pam50", "var", xlab='', ylab=varlabel, title=var, outlier.shape=NA,
                 fill='subtype_pam50', palette=colors_pam50, order=c('Basal', 'Her2', 'LumA', 'LumB', 'Normal')) +
    stat_compare_means(comp=comps2, size=4, label='p.format') +
    theme(text = element_text(size = 10, face='bold'), legend.position = "none") +
    geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);
  
  print(ggarrange(p1, p2, nrow=2, ncol=2));
  
  #------------------------------------------------------------
  # (6) correlation of gene CN vs. expression
  #

  # scatter plot of gene CN vs. expression
  p1 <- ggscatter(mdata, x = 'var', y = "exp", size = 1, add = "reg.line", conf.int = TRUE, 
                    xlab=varlabel, ylab='log2TPM', main=gene) +
      stat_cor(method = "spearman", col='black', cex=4) + geom_point(pch=21, fill=collabs_prepost, size=3) +
      theme(text = element_text(size = 12, face='bold'));
    
  # gene CN vs. expression, facet by "Pre/Post"  
  p2 <- ggscatter(mdata, x = 'var', y = "exp", size = 1, add = "reg.line", conf.int = TRUE, 
                    xlab=varlabel, ylab='log2TPM', main=gene, facet.by='prepost') +
      stat_cor(method = "spearman", col='black', cex=4) + geom_point(pch=21, fill=collabs_prepost, size=3) +
    theme(text = element_text(size = 12, face='bold'));
    
  print(ggarrange(p1, p2, ncol=2, nrow=2, widths=c(1,2)));
  
  dev.off();
}


#===================================================
# (6) Add a custom gene signature (ASC Epithelial)
#
add_gsva_ASC <- function(rdata, gsvadata){
  
  library('GSVA');
  library('qusage');
  
  # input/output files
  gmtfile = paste(datadir, 'ASC_epithelial.gmt', sep='/');
  geneset = read.gmt(gmtfile);
  gsname = paste('custom', names(geneset)[1], sep='.');
  
  # gene expression data
  expdata = rdata$Expression$log2TPM_corrected;
  
  #------------------------------------------------------------
  # calculate GSVA scores for custom signatures
  #
  gsva_asc = gsva(expr = as.matrix(expdata),
                  gset.idx.list = geneset,
                  kcdf=c("Gaussian"),
                  mx.diff = TRUE,
                  parallel.sz=4,
                  abs.ranking=FALSE);
  
  # merge with existing GSVA scores
  #gsvaobj <- tempus_bc_gsva(rdata);
  gsvadata <- gsvaobj$scores;
  midx = match(colnames(gsvadata), colnames(gsva_asc));
  gsvadata = rbind(gsvadata, gsva_asc[,midx]);
  rownames(gsvadata)[dim(gsvadata)[1]] = gsname;
  
  df <- data.frame(geneset=gsname, collection='Pfizer_custom', description='Smith et al. 2018. Cell Reports')
  annot <- rbind(gsvaobj$annotation, df);
  
  gsvaobj <- list(scores=gsvadata, annotation=annot);
  return(gsvaobj)
}

#=====================================================
# (7) BayesPrism cell-type specific
# gene expression data (CSE)
# VST: Variance Stabilization Transformation
#
BayesPrism_cse <- function(){
  
  library(DESeq2);
  
  #projdir <- '/Users/kanzhengyan/Documents/Clinical Genomics/Tempus_BC';
  #rdatafile <- paste(projdir, 'Data Analysis/Whijae/2022-08-24.BayesPRISM_Output/tempus.ted.Tempus_brca_deconv_by_wu.v2.rdata', sep='/');
  projdir <- '..'
  rdatafile <- paste(projdir, 'Data Analysis/Whijae/2023-01-10.BayesPRISM_Output_Updated/tempus.batch.corrected.ted.Tempus_brca_deconv_by_wu.v2.rdata', sep='/');
  rdata_cse <- get(load(rdatafile));
  cell.type.specific.matrix <- rdata_cse$res$first.gibbs.res$Znkg.merged;
  csedata_cancer <- varianceStabilizingTransformation(round(t(cell.type.specific.matrix[,"Cancer",])) + 1);
  
  return(csedata_cancer);
}

#=====================================================
# (8) Add NMF-based scores for CDK dependencies
#
add_feature_nmf_response_score <- function(fdata){
  fname <- paste(datadir, 'CDK4.CDK6.response.scores.Tempus.tsv', sep='/');
  scoredata <- read.table(fname, sep='\t', header=T);
  tmp <- gsub('\\.', '_', tolower(colnames(scoredata)));
  tmp <- gsub('response_scores_bulk', 'bulk_score', tmp);
  tmp <- gsub('response_scores_tumor_specific', 'cse_score', tmp);
  tmp <- gsub('tumor_specific', 'cse', tmp);
  tmp <- gsub('_tempus', '', tmp);
  tmp <- gsub('subtypes', 'subtype', tmp);
  tmp[5:length(tmp)] <- paste('NMFRS', tmp[5:length(tmp)], sep='_')
  colnames(scoredata) <- tmp;
  scoredata$sample_id <- gsub('^pfizer ', '', scoredata[,1]);
  
  hdrs = c('sample_id', "NMFRS_cdk4_bulk_score", "NMFRS_cdk6_bulk_score", "NMFRS_cdk4_cse_score", "NMFRS_cdk6_cse_score");
  mdata <- merge(fdata, subset(scoredata, select=hdrs), by='sample_id', all.x=T, sort=F);
  return(mdata);
}

#=====================================================
# (9) Add molecular features (e.g. IC1-5) to "sdata"
#
add_molecular_features<-function(sdata){
  molfeatfile <- paste(datadir, 'Tempus_molecular_features_extended_v3.txt', sep='/');
  mfdata <- read.table(molfeatfile, header=T, sep='\t');
  idx <- grepl('^pam50_', colnames(mfdata)) |
    grepl('^nmf_', colnames(mfdata)) |
    grepl('^Paloma3_', colnames(mfdata)) |
    grepl('^MUT', colnames(mfdata)) |
    grepl('^AMP', colnames(mfdata)) |
    grepl('^DEL', colnames(mfdata)) |
    grepl('^EXP', colnames(mfdata)) |
    grepl('^facets', colnames(mfdata)) |
    grepl('^h.HALLMARK_', colnames(mfdata)) |
    colnames(mfdata) %in% c('sample_id', 'icluster');
  
  #-------------------------------------------
  # add icluster+pam50 subgroups
  #
  groups <- paste(mfdata$icluster, mfdata$subtype_pam50, sep=':')
  selected_groups <- names(table(groups))[table(groups) >= 7];
  idx2 <- groups %in% selected_groups & mfdata$icluster != 'IC5' & !is.na(mfdata$icluster)
  mfdata$ic_pam50 <- groups;
  mfdata$ic_pam50[!idx2] <- NA;
  
  mdata <- merge(sdata, mfdata[,idx], by='sample_id', sort=F)
  return(mdata)
}


#=============================================================
# Differential analysis of numerical variables ("data)
# comparing group variable e.g. Pre/Post ("class")
#
DE_class <- function(data, foo){
  
  res <- apply(data,1, function(x){
    fm <- lm(x~class,data=foo);
    return(data.frame(pval=anova(fm)["class","Pr(>F)"], statistic=coef(summary(fm))[2,"Estimate"]));
  });
  
  # FDR q-value
  result <- data.frame(variable=names(res), matrix(unlist(res), ncol=2, byrow=T));
  colnames(result)[2:3] = c('pval', 'coef');
  result$qval <- p.adjust(result$pval, method="BH");
  
  # calculate group means
  means <- t(apply(data, 1, function(x) { fm <- tapply(x, foo$class, mean, na.rm=T) }));
  colnames(means) = c('mean1', 'mean2');
  result = data.frame(result, means);
  
  return(result)
}

