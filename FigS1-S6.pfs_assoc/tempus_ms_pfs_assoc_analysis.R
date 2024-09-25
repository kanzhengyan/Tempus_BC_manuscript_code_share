## Author: Vini Bonato
## Purpose: BC Tempus manuscript
## Created on: 06/27/2023
## Last modified on: 09/24/2024

## Libraries
library(openxlsx)
library(survival)
library(coxphf)
library(qvalue)
library(forestplot)
library(plyr)
library(Cairo)
library(survminer)

sessionInfo()
## R version 4.1.0 (2021-05-18)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 19045)

## Matrix products: default

## locale:
## [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
## [4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    

## attached base packages:
## [1] grid      stats     graphics  grDevices utils     datasets  methods   base     

## other attached packages:
## [1] survminer_0.4.9  ggpubr_0.6.0     ggplot2_3.5.1    Cairo_1.6-0      plyr_1.8.8       forestplot_3.1.3 abind_1.4-5     
## [8] checkmate_2.1.0  qvalue_2.26.0    coxphf_1.13.1    survival_3.5-5   openxlsx_4.2.5.2

## loaded via a namespace (and not attached):
## [1] zoo_1.8-12        tidyselect_1.2.1  xfun_0.45         reshape2_1.4.4    purrr_1.0.1       splines_4.1.0    
## [7] lattice_0.21-8    carData_3.0-5     colorspace_2.1-0  vctrs_0.6.1       generics_0.1.3    utf8_1.2.3       
## [13] survMisc_0.5.6    rlang_1.1.4       pillar_1.9.0      glue_1.6.2        withr_3.0.0       lifecycle_1.0.4  
## [19] stringr_1.5.1     munsell_0.5.1     ggsignif_0.6.4    gtable_0.3.5      zip_2.3.0         knitr_1.47       
## [25] fansi_1.0.4       broom_1.0.6       Rcpp_1.0.13       xtable_1.8-4      backports_1.4.1   scales_1.3.0     
## [31] km.ci_0.5-6       gridExtra_2.3     stringi_1.7.12    rstatix_0.7.2     dplyr_1.1.2       KMsurv_0.1-5     
## [37] cli_3.6.1         tools_4.1.0       magrittr_2.0.3    tibble_3.2.1      tidyr_1.3.0       car_3.1-2        
## [43] pkgconfig_2.0.3   Matrix_1.7-0      data.table_1.14.8 R6_2.5.1          compiler_4.1.0   


## Aux functions

leveling <- function(x) {
  ux <- unique(x)
  ux[order(tabulate(match(x, ux)), decreasing = TRUE)]
}

GetHRs <- function(dataset, covariates, sample.time, feature.type){

	covariates$summ.stats <- NA
	for(i in 1:nrow(covariates)){
		if(covariates$type[i] == "categorical"){
			new.levels                        <- leveling(dataset[, covariates$variable[i]])
			dataset[, covariates$variable[i]] <- addNA(factor(dataset[, covariates$variable[i]], levels = new.levels))
			covariates$summ.stats[i]          <- paste(paste(levels(dataset[, covariates$variable[i]]), " (N=", as.vector(summary(dataset[, covariates$variable[i]])), ")", sep = ""), collapse = "; ", sep = "")
			dataset[, covariates$variable[i]] <- factor(dataset[, covariates$variable[i]], levels = levels(dataset[, covariates$variable[i]])[!is.na(dataset[, covariates$variable[i]])])
		}
		if(covariates$type[i] == "mutation"){
			dataset[, covariates$variable[i]] <- addNA(factor(dataset[, covariates$variable[i]], levels = c("0", "1"), labels = c("WT", "Mut")))
			covariates$summ.stats[i] <- paste(paste(levels(dataset[, covariates$variable[i]]), " (N=", as.vector(summary(dataset[, covariates$variable[i]])), ")", sep = ""), collapse = "; ", sep = "")
			dataset[, covariates$variable[i]] <- factor(dataset[, covariates$variable[i]], levels = levels(dataset[, covariates$variable[i]])[!is.na(dataset[, covariates$variable[i]])])
		}
		if(covariates$type[i] == "numerical"){
			dataset[, covariates$variable[i]] <- as.numeric(as.vector(dataset[, covariates$variable[i]]))
			covariates$summ.stats[i]          <- paste("Mean = ", format(mean(dataset[, covariates$variable[i]], na.rm = TRUE), digits = 3), 
												  	   "; SD = ", format(sd(dataset[, covariates$variable[i]], na.rm = TRUE), digits = 3), 
													   "; Median = ", format(median(dataset[, covariates$variable[i]], na.rm = TRUE), digits = 3), 
													   sep = "")
			dataset[, covariates$variable[i]] <- scale(dataset[, covariates$variable[i]])
		}
	}				

	## CPH-Firth models, unadjusted:
	univ.formulas <- sapply(covariates$variable, function(x) as.formula(paste('Surv(pfs, pfs_event)~ age_at_cdk + cdk_metastatic_lot + ', x)))
	univ.models   <- lapply(univ.formulas, function(x){
										variables <- strsplit(as.character(as.vector(x))[[3]], " + ", fixed = TRUE)[[1]]
										test.data <- dataset[, c("pfs", "pfs_event", variables)]
										test.data <- test.data[complete.cases(test.data), ]
										try(coxphf(x, data = test.data, maxstep = 0.01, maxit = 1000, firth = FALSE), silent = TRUE)
									 })
	univ.results <- lapply(univ.models,
						   function(x){
								if(class(x)[1] != "try-error"){
									variable.name   <- strsplit(strsplit(as.character(x$formula), " ~ ")[[3]], " + ", fixed = TRUE)[[1]][3]
									overall.p.value <- signif(pchisq(as.numeric(-2 * x$loglik[1] + 2 * x$loglik[2]), df = x$df, lower.tail = F), digits = 2)
									res.out         <- NULL
									start.j         <- min(grep(variable.name, names(x$coefficients)))
									for(j in start.j:length(x$coefficients)){
										  contrast     <- paste(gsub(variable.name, "", names(x$coefficients)[j]), " vs. ", levels(dataset[, variable.name])[1], sep = "")
										  beta         <- as.numeric(as.vector(x$coefficients[j]))
										  HR           <- exp(as.numeric(as.vector(x$coefficients[j])))
										  LB           <- log(as.numeric(as.vector(x$ci.lower[j])))
										  UB           <- log(as.numeric(as.vector(x$ci.upper[j])))
										  rawP         <- as.numeric(as.vector(x$prob))[j]
										  signedlog10P <- sign(beta) * (-log10(rawP))
										  res          <- data.frame(contrast = contrast, est = beta, est.LB95 = LB, est.UB95 = UB, HR = HR, pval = rawP, signedLog10p = signedlog10P)
										  res.out[[j]] <- res}
									res.out <- rbind.fill(res.out)
									res.out$variable <- variable.name
									return(res.out)
								}
							})
                           
	results <- rbind.fill(univ.results)
	results <- merge(covariates, results, by = "variable", all = TRUE, sort = FALSE)
	results$contrast[results$contrast == " vs. "] <- "---"
	results$contrast[is.na(results$contrast)] <- "---"
	results$variable.labels <- paste(results$variable, " (", results$contrast, ")", sep = "")
	results$variable.labels <- gsub(" (---)", "", results$variable.labels, fixed = TRUE)
	results$summ.stats      <- gsub("; NA (N=0)", "", results$summ.stats, fixed = TRUE)
	results$Sample.Time     <- sample.time
	results$feature.type    <- feature.type
	results
}

## Paths
workDir      <- "./FigS1-S6.pfs_assoc/"
inputFolder  <- paste(workDir, "Input/", sep = "")
outputFolder <- paste(workDir, "Output.Manuscript/", sep = "")


## Patient Info and Covariates:
patient.info                    <- read.xlsx(paste(inputFolder, "Supplementary_Table_01.sample_info_molecular_features.xlsx", sep = ""), sheet = "molecular features")
patient.info                    <- patient.info[, c("sample_id_anonymized", "patient_id_anonymized", "pre_post", "pre_during_post", "pfs", "pfs_event")]
covariates                      <- read.table(paste(inputFolder, "tempus_ms_patient_covariates.txt", sep = ""), sep = "\t", header = TRUE)
covariates$pre_post             <- ifelse(covariates$PFEPrePost == "Pre" & !is.na(covariates$PFEPrePost), "Pre", "Post")
covariates$sample_id_anonymized <- paste(covariates$patient_id_anonymized, covariates$pre_post, sep = "_")
covariates                      <- covariates[, c("sample_id_anonymized", "age_at_cdk", "cdk_metastatic_lot")]
patient.info                    <- merge(patient.info, covariates, by = "sample_id_anonymized", all.x = TRUE, sort = FALSE)

## Molecular Features:
anal.mol.features.data    <- read.xlsx(paste(inputFolder, "Supplementary_Table_01.sample_info_molecular_features.xlsx", sep = ""), sheet = "molecular features")
anal.mol.features.data    <- merge(anal.mol.features.data, patient.info[, c("sample_id_anonymized", "age_at_cdk", "cdk_metastatic_lot")], by = "sample_id_anonymized", sort = FALSE, all = TRUE)
anal.mol.features.cov.cat <- data.frame(variable = c("subtype_pam50", "mutsig_sigma_status", "brca_pathogenic_mutation", "bc_susceptibility_gene_pathogenic_mutation",
													 "brca2_germline", "chek2_germline", "atm_germline", "ddr_germline"),
										type     = rep("categorical", 8))
anal.mol.features.cov.num <- data.frame(variable = c("pam50_cor_basal", "pam50_cor_her2", "pam50_cor_luma", "pam50_cor_lumb", "pam50_cor_normal", "tmb", "tmb_norm", "cyt_score", 
													 "facets_tumor_purity", "hrd_index", "proliferative_index", "estimate_stromalscore", "estimate_immunescore", "estimate_score", 
													 "estimate_tumor_purity", "facets_ploidy", "facets_loh", "facets_telomeric_ai", "facets_lst", "mutsig_1_sigma", "mutsig_2_sigma", 
													 "mutsig_3_sigma", "mutsig_13_sigma", "nmf_factor1_bone_specific", "nmf_factor2_basal", "nmf_factor3_ER_expr", "nmf_factor4_ER_HER2_expr", 
													 "nmf_factor5_lumA", "nmf_factor6_immune", "nmf_factor7_liver_specific", "nmf_factor8_ER_response", "nmf_factor9_EMT", 
													 "nmf_factor10_stomach_colon_specific", "nmf_factor11_proliferation", "nmf_factor12_HER2_expr", "Paloma3_F1_MYC_E2F_activation", 
													 "Paloma3_F2", "Paloma3_F3", "Paloma3_F4", "Paloma3_F5", "Paloma3_F6", "Paloma3_F7_IFNG_response", "Paloma3_F8_EMT", "Paloma3_F9", 
													 "Paloma3_F10_estrogen_response", "Paloma3_F11", "Paloma3_F12", "Paloma3_F13", "Paloma3_F14_liver_specific", "est_immune_cells", 
													 "est_b_cells", "est_cd4_cells", "est_cd8_cells", "est_mac_cells", "est_nk_cells"),
										type     = rep("numerical", 55))
anal.mol.features.cov      <- rbind(anal.mol.features.cov.cat, anal.mol.features.cov.num)
anal.mol.features.pre.data <- subset(anal.mol.features.data, pre_post == "Pre")
anal.mol.features.pos.data <- subset(anal.mol.features.data, pre_post == "Post")

## Hallmark Signatures:
hallmark.features.data    <- read.xlsx(paste(inputFolder, "Supplementary_Table_01.sample_info_molecular_features.xlsx", sep = ""), sheet = "hallmark gsva")
hallmark.features.data    <- merge(hallmark.features.data, patient.info, by = "sample_id_anonymized", sort = FALSE, all = TRUE, suffix = c("", ".x"))
hallmark.features.data$patient_id_anonymized.x <- NULL
hallmark.features.cov.num <- data.frame(variable = c("TNFA_SIGNALING_VIA_NFKB", "HYPOXIA", "CHOLESTEROL_HOMEOSTASIS", "MITOTIC_SPINDLE", "WNT_BETA_CATENIN_SIGNALING", "TGF_BETA_SIGNALING", 
													 "IL6_JAK_STAT3_SIGNALING", "DNA_REPAIR", "G2M_CHECKPOINT", "APOPTOSIS", "NOTCH_SIGNALING", "ADIPOGENESIS", "ESTROGEN_RESPONSE_EARLY", 
													 "ESTROGEN_RESPONSE_LATE", "ANDROGEN_RESPONSE", "MYOGENESIS", "PROTEIN_SECRETION", "INTERFERON_ALPHA_RESPONSE", "INTERFERON_GAMMA_RESPONSE", 
													 "APICAL_JUNCTION", "APICAL_SURFACE", "HEDGEHOG_SIGNALING", "COMPLEMENT", "UNFOLDED_PROTEIN_RESPONSE", "PI3K_AKT_MTOR_SIGNALING", 
													 "MTORC1_SIGNALING", "E2F_TARGETS", "MYC_TARGETS_V1", "MYC_TARGETS_V2", "EPITHELIAL_MESENCHYMAL_TRANSITION", "INFLAMMATORY_RESPONSE", 
													 "XENOBIOTIC_METABOLISM", "FATTY_ACID_METABOLISM", "OXIDATIVE_PHOSPHORYLATION", "GLYCOLYSIS", "REACTIVE_OXYGEN_SPECIES_PATHWAY", 
													 "P53_PATHWAY", "UV_RESPONSE_UP", "UV_RESPONSE_DN", "ANGIOGENESIS", "HEME_METABOLISM", "COAGULATION", "IL2_STAT5_SIGNALING", "BILE_ACID_METABOLISM", 
													 "PEROXISOME", "ALLOGRAFT_REJECTION", "SPERMATOGENESIS", "KRAS_SIGNALING_UP", "KRAS_SIGNALING_DN", "PANCREAS_BETA_CELLS"),
										type     = rep("numerical", 50))
hallmark.features.cov      <- hallmark.features.cov.num
hallmark.features.pre.data <- subset(hallmark.features.data, pre_post == "Pre")
hallmark.features.pos.data <- subset(hallmark.features.data, pre_post == "Post")

## Genomic Alterations:
gen.alt.features.data    <- read.xlsx(paste(inputFolder, "Supplementary_Table_01.sample_info_molecular_features.xlsx", sep = ""), sheet = "genomic features")
gen.alt.features.data    <- merge(gen.alt.features.data, patient.info, by = "sample_id_anonymized", sort = FALSE, all = TRUE)
gen.alt.features.data$patient_id_anonymized.x <- NULL
gen.alt.features.cov.mut <- data.frame(variable = c("MUT.ABCB1", "MUT.ABCC3", "MUT.ABL1", "MUT.AKT1", "MUT.ALK", "MUT.APC", "MUT.APOB", "MUT.ARID1A", "MUT.ARID1B", "MUT.ARID2", "MUT.ASXL1", "MUT.ATM", 
													"MUT.ATR", "MUT.ATRX", "MUT.BCL11B", "MUT.BCLAF1", "MUT.BCOR", "MUT.BCORL1", "MUT.BCR", "MUT.BRAF", "MUT.BRCA1", "MUT.BRCA2", "MUT.BRD4", "MUT.CARD11", 
													"MUT.CBFB", "MUT.CD40", "MUT.CDH1", "MUT.CDKN2A", "MUT.CEBPA", "MUT.CFTR", "MUT.CHD4", "MUT.CIC", "MUT.CIITA", "MUT.CREBBP", "MUT.CSF3R", "MUT.CUX1", 
													"MUT.DICER1", "MUT.DIS3L2", "MUT.DNM2", "MUT.DNMT3A", "MUT.DOT1L", "MUT.DPYD", "MUT.DYNC2H1", "MUT.EGFR", "MUT.EP300", "MUT.EPHB1", "MUT.EPHB2", "MUT.ERBB2", 
													"MUT.ERBB3", "MUT.ERBB4", "MUT.ERCC4", "MUT.ERCC6", "MUT.ESR1", "MUT.FANCA", "MUT.FANCM", "MUT.FAT1", "MUT.FGFR2", "MUT.FGFR3", "MUT.FLT3", "MUT.FLT4",
													"MUT.FOXA1", "MUT.FOXQ1", "MUT.GALNT12", "MUT.GATA3", "MUT.GATA6", "MUT.GRIN2A", "MUT.HNF1A", "MUT.IFNAR2", "MUT.IKZF1", "MUT.IRS2", "MUT.JAK1", "MUT.JAK3", 
													"MUT.KAT6A", "MUT.KDM5A", "MUT.KIF1B", "MUT.KIT", "MUT.KMT2A", "MUT.KMT2B", "MUT.KMT2C", "MUT.KMT2D", "MUT.KRAS", "MUT.LRP1B", "MUT.MAP2K4", "MUT.MAP3K1", 
													"MUT.MED12", "MUT.MEN1", "MUT.MKI67", "MUT.MLH1", "MUT.MSH3", "MUT.MTOR", "MUT.MYH11", "MUT.NCOR1", "MUT.NCOR2", "MUT.NF1", "MUT.NOTCH1", "MUT.NOTCH2", 
													"MUT.NOTCH3", "MUT.NRG1", "MUT.PALLD", "MUT.PBRM1", "MUT.PDGFRB", "MUT.PIK3CA", "MUT.PIK3R1", "MUT.PIK3R2", "MUT.PLCG2", "MUT.PMS2", "MUT.POLE", "MUT.PREX2", 
													"MUT.PTCH1", "MUT.PTEN", "MUT.PTPN13", "MUT.PTPRD", "MUT.RAD21", "MUT.RANBP2", "MUT.RB1", "MUT.RECQL4", "MUT.RET", "MUT.ROS1", "MUT.RUNX1", "MUT.SETBP1",
													"MUT.SETD2", "MUT.SF3B1", "MUT.SLIT2", "MUT.SLX4", "MUT.SMARCA4", "MUT.SMO", "MUT.SPEN", "MUT.STAG2", "MUT.TAF1", "MUT.TBX3", "MUT.TERT", "MUT.TET2", 
													"MUT.TP53", "MUT.TSC1", "MUT.TSC2", "MUT.ZFHX3", "MUT.ZNF217", "MUT.ZNF750", "MUT.ZNRF3", "AMP.AURKA", "AMP.AXIN2", "AMP.BRIP1", "AMP.CCND1", "AMP.CD79B", 
													"AMP.CDKN1B", "AMP.CKS1B", "AMP.ELF3", "AMP.FGF3", "AMP.FGF4", "AMP.FGFR1", "AMP.FOXA1", "AMP.FRS2", "AMP.GNA13", "AMP.GNAS", "AMP.H3.3A", "AMP.IDO1", 
													"AMP.IKBKE", "AMP.KAT6A", "AMP.MCL1", "AMP.MDM2", "AMP.MDM4", "AMP.MTMR11", "AMP.MYC", "AMP.NDUFC2", "AMP.NTHL1", "AMP.PAK1", "AMP.PIK3C2B", "AMP.RAD21", 
													"AMP.RAD51C", "AMP.RECQL4", "AMP.RNF139", "AMP.RNF43", "AMP.RPS6KB1", "AMP.RSF1", "AMP.SPOP", "AMP.SRSF2", "AMP.UBE2T", "AMP.YEATS4", "AMP.ZNF217", 
													"AMP.ZNF750", "DEL.ALK", "DEL.APLNR", "DEL.CDKN2A", "DEL.CDKN2B", "DEL.CYP1B1", "DEL.CYP2D6", "DEL.EGFR", "DEL.EPHB1", "DEL.FGF14", "DEL.MTAP", "DEL.NKX2.1",
													"DEL.NRG1", "DEL.PDE4D", "DEL.PRSS1"),
										type     = rep("mutation", 194))
gen.alt.features.cov      <- gen.alt.features.cov.mut
gen.alt.features.pre.data <- subset(gen.alt.features.data, pre_post == "Pre")
gen.alt.features.pos.data <- subset(gen.alt.features.data, pre_post == "Post")

anal.mol.features.pre.res <- GetHRs(dataset = anal.mol.features.pre.data, covariates = anal.mol.features.cov, sample.time = "Pre",  feature.type = "Analytical Molecular Features")
anal.mol.features.pos.res <- GetHRs(dataset = anal.mol.features.pos.data, covariates = anal.mol.features.cov, sample.time = "Post", feature.type = "Analytical Molecular Features")
hallmark.features.pre.res <- GetHRs(dataset = hallmark.features.pre.data, covariates = hallmark.features.cov, sample.time = "Pre",  feature.type = "Hallmark Signature Scores")
hallmark.features.pos.res <- GetHRs(dataset = hallmark.features.pos.data, covariates = hallmark.features.cov, sample.time = "Post", feature.type = "Hallmark Signature Scores")
gen.alt.features.pre.res  <- GetHRs(dataset = gen.alt.features.pre.data,  covariates = gen.alt.features.cov,  sample.time = "Pre",  feature.type = "Genomic Alteration Features")
gen.alt.features.pos.res  <- GetHRs(dataset = gen.alt.features.pos.data,  covariates = gen.alt.features.cov,  sample.time = "Post", feature.type = "Genomic Alteration Features")

stacked.pre      <- rbind.fill(anal.mol.features.pre.res, hallmark.features.pre.res, gen.alt.features.pre.res)
stacked.pos      <- rbind.fill(anal.mol.features.pos.res, hallmark.features.pos.res, gen.alt.features.pos.res)
stacked.pre$qval <- qvalue(p = stacked.pre$pval)$qvalues
stacked.pos$qval <- qvalue(p = stacked.pos$pval)$qvalues

out.list <- list("molecular features (Pre-CDK)"  = subset(stacked.pre, feature.type == "Analytical Molecular Features"),
				 "molecular features (Post-CDK)" = subset(stacked.pos, feature.type == "Analytical Molecular Features"),
				 "hallmark gsva (Pre-CDK)"       = subset(stacked.pre, feature.type == "Hallmark Signature Scores"),
				 "hallmark gsva (Post-CDK)"      = subset(stacked.pos, feature.type == "Hallmark Signature Scores"),
				 "genomic features (Pre-CDK)"    = subset(stacked.pre, feature.type == "Genomic Alteration Features"),
				 "genomic features (Post-CDK)"   = subset(stacked.pos, feature.type == "Genomic Alteration Features"))

write.xlsx(out.list, file = paste(outputFolder, "modelResults.xlsx", sep = ""), rowNames = FALSE, colNames = TRUE, keepNA = TRUE, na.string = "")

## Plots:

anal.mol.features.pre.plot <- anal.mol.features.pre.res
anal.mol.features.pre.plot <- subset(anal.mol.features.pre.plot, pval < 0.05 & variable.labels != "subtype_pam50 (Normal vs. LumB)")
anal.mol.features.pre.plot <- anal.mol.features.pre.plot[order(anal.mol.features.pre.plot$HR), ]
anal.mol.features.pre.plot$variable.labels <- c("Paloma3 signature F10 (Estrogen response)", "NMF factor5 (Luminal A)", "Paloma3 signature F5", 
												"PAM50 score (Luminal A)", "Paloma3 signature F2", "PAM50 score (Normal)",
												"estimate_stromalscore", "Paloma3 signature F11 ()", "Paloma3 signature F8 (EMT)",
												"PAM50 score (Luminal B)", "MutSig (13 sigma)", "NMF factor12 (HER2 expr)", 
												"FACETS LOH", "PAM50 score (Basal)", "Paloma3 signature F9 ()",
												"Proliferative Index", "NMF factor11 (proliferation)", "PAM50 score (Her2)", 
												"Paloma3 signature F1 (MYC E2F activation)", "PAM50 subtype (Luminal B vs. Luminal A)", "DDR germline mutation (Y vs. N)",
												"BRCA pathogenic mutation status (Y vs. N)", "PAM50 subtype (Normal vs. Luminal A)")

CairoPDF(paste(outputFolder, "Fig-S1a.forestplot.analmol.pre.pdf", sep = ""), height = 8, width = 16)
    forestplot(labeltext = list(anal.mol.features.pre.plot$variable.labels), 
               mean      = as.numeric(as.vector(anal.mol.features.pre.plot$est)),
               lower     = as.numeric(as.vector(anal.mol.features.pre.plot$est.LB95)),
               upper     = as.numeric(as.vector(anal.mol.features.pre.plot$est.UB95)),
               zero      = 0,
			   lwd.zero  = 1.2,
			   mar       = unit(c(5, 5, 5, 6), "mm"),
			   clip      = c(-5, 1.75), 
               xticks    = seq(-1, 1.5, 0.25),
			   grid      = structure(seq(-1, 1.5, 0.25), gp = gpar(col = "gray75", lty = 3, lwd = 0.5)),
			   col       = fpColors(box = "gray65", line = "gray65", summary = "gray65", zero = "black"),
               new_page  = FALSE,
               txt_gp    = fpTxtGp(label = gpar(cex = 1.7), xlab = gpar(cex = 2, fontface = "bold"), tick = gpar(cex = 1.7, fontface = "bold")),
               xlab      = "CPH model estimate: log(HR)")
dev.off()

hallmark.features.pre.plot <- hallmark.features.pre.res
hallmark.features.pre.plot <- subset(hallmark.features.pre.plot, pval < 0.05)
hallmark.features.pre.plot <- hallmark.features.pre.plot[order(hallmark.features.pre.plot$HR), ]
hallmark.features.pre.plot$variable.labels <- c("WNT Beta Catenin signaling", "KRAS signaling down", "Estrogen response early", "UV response down",
												"Apical surface", "Mitotic spindle", "Xenobiotic metabolism", "DNA repair", "MYC targets V2", 
												"Cholesterol homeostasis", "Oxidative phosphorylation", "Protein secretion", "E2F targets", 
												"Peroxisome", "PI3K AKT MTOR signaling", "G2M checkpoint", "UV response up", "MYC targets V1", 
												"Unfolded protein response", "Glycolysis", "MTORC1 signaling")

CairoPDF(paste(outputFolder, "Fig-S1b.forestplot.hallmarkgs.pre.pdf", sep = ""), height = 8, width = 16)
    forestplot(labeltext = list(hallmark.features.pre.plot$variable.labels), 
               mean      = as.numeric(as.vector(hallmark.features.pre.plot$est)),
               lower     = as.numeric(as.vector(hallmark.features.pre.plot$est.LB95)),
               upper     = as.numeric(as.vector(hallmark.features.pre.plot$est.UB95)),
               zero      = 0,
			   lwd.zero  = 1.2,
			   mar       = unit(c(5, 5, 5, 8), "mm"),
			   clip      = c(-5, 1.75), 
               xticks    = seq(-0.5, 0.75, 0.25),
			   grid      = structure(seq(-0.5, 0.75, 0.25), gp = gpar(col = "gray75", lty = 3, lwd = 0.5)),
			   col       = fpColors(box = "gray65", line = "gray65", summary = "gray65", zero = "black"),
               new_page  = FALSE,
               txt_gp    = fpTxtGp(label = gpar(cex = 2), xlab = gpar(cex = 2, fontface = "bold"), tick = gpar(cex = 1.7, fontface = "bold")),
               xlab      = "CPH model estimate: log(HR)")
dev.off()

################################################################################################################################################################################################
## KM-plots

additional.features.data     <- read.table(paste(inputFolder, "Tempus_molecular_features_extended.annonymized.txt", sep = ""), sep = "\t", header = TRUE)
additional.features.data     <- additional.features.data[, c("sample_id_anonymized", "c2_cgp.SMID_BREAST_CANCER_ERBB2_UP", "EXPR.CCNE1", "EXPR.ESR1", "EXPR.PGR")]
additional.features.data     <- merge(patient.info, additional.features.data, by = "sample_id_anonymized", sort = FALSE)
names(additional.features.data)[names(additional.features.data) == "c2_cgp.SMID_BREAST_CANCER_ERBB2_UP"] <- "SMID_BREAST_CANCER_ERBB2_UP"
additional.features.pre.data <- subset(additional.features.data, pre_post == "Pre")
additional.features.cov.num  <- data.frame(variable = c("SMID_BREAST_CANCER_ERBB2_UP", "EXPR.CCNE1", "EXPR.ESR1", "EXPR.PGR"), type = rep("numerical", 4))
additional.features.pre.res  <- GetHRs(dataset = additional.features.pre.data, covariates = additional.features.cov.num, sample.time = "Pre", feature.type = "Analytical Molecular Features")
stacked.pre                  <- rbind.fill(stacked.pre, additional.features.pre.res)
stacked.pre$qval             <- qvalue(p = stacked.pre$pval)$qvalues

km.features <- c("MUT.TP53", "MUT.RB1", "MUT.ESR1", "proliferative_index", "pam50_cor_basal", "pam50_cor_her2", "pam50_cor_luma", 
				 "EXPR.CCNE1", "Paloma3_F1_MYC_E2F_activation", "E2F_TARGETS", "MYC_TARGETS_V1", "EXPR.ESR1", "EXPR.PGR", 
				 "Paloma3_F10_estrogen_response")
					   
features.set1 <- anal.mol.features.pre.data[, c("sample_id_anonymized", "pfs", "pfs_event", names(anal.mol.features.pre.data)[names(anal.mol.features.pre.data) %in% km.features])]					   
features.set2 <- hallmark.features.pre.data[, c("sample_id_anonymized", names(hallmark.features.pre.data)[names(hallmark.features.pre.data) %in% km.features])]					   
features.set3 <- gen.alt.features.pre.data[, c("sample_id_anonymized", names(gen.alt.features.pre.data)[names(gen.alt.features.pre.data) %in% km.features])]					   
features.set4 <- additional.features.pre.data[, c("sample_id_anonymized", names(additional.features.pre.data)[names(additional.features.pre.data) %in% km.features])]					   

km.data <- merge(features.set1, features.set2, by = "sample_id_anonymized", sort = FALSE)
km.data <- merge(km.data, features.set3, by = "sample_id_anonymized", sort = FALSE)
km.data <- merge(km.data, features.set4, by = "sample_id_anonymized", sort = FALSE)
all(km.features %in% names(km.data))
#[1] TRUE

#----------------------------------------------
# KM plots for selected features
#
pdf(file=paste(outputFolder, 'Fig-S1-S4-S5-S6.kmplots.selectedfeatures.pre.pdf', sep='/'), onefile=T, height=9, width=8, pointsize=12);

for (i in 1:length(km.features)) {
  feature <- km.features[i];
  tmp     <- km.data[,feature];
  r       <- stacked.pre[stacked.pre$variable == feature,]
  # annotate p-value, HR from CPH analysis
  label1 <- paste('p=', format(r$pval, digits=4, sci = TRUE), sep='');
  label2 <- paste('HR=', format(r$HR, digits=4, sci = FALSE), sep='');
  y.label1 <- 1
  y.label2 <- 0.9
   
  if(unique(r$type) == "mutation"){
	  if (grepl('^MUT', feature)){
		vals <- ifelse(tmp == 1, 'MUT', 'WT')
		} else if(grepl('^AMP', feature)) {
		  vals <- ifelse(tmp == 1, 'AMP', 'WT')
		} else if(grepl('^DEL', feature)){
		  vals <- ifelse(tmp == 1, 'DEL', 'WT')
		}
  }

  if(unique(r$type) == "numerical"){
	  vals <- ifelse(tmp > median(tmp, na.rm = T), paste("High (>", round(median(tmp, na.rm = T), 4), ")", sep = ""), paste("Low (<=", round(median(tmp, na.rm = TRUE), 4), ")", sep = "")) 
  }

  if(unique(r$type) == "categorical" & unique(r$variable) == "brca_pathogenic_mutation"){
	  vals <- ifelse(tmp == "Yes", 'MUT', 'WT')
  }
  
  if(unique(r$type) == "categorical" & unique(r$variable) == "subtype_pam50"){
	  vals     <- tmp
	  vals[tmp %in% c("Basal", "Normal")] <- NA
	  label1   <- paste('LumA vs LumB:\n p=', format(subset(r, contrast == "LumA vs. LumB")$pval, digits=3, sci = TRUE), "; HR=", format(subset(r, contrast == "LumA vs. LumB")$HR, digits=3, sci = FALSE), sep = "")
	  label2   <- paste('Her2 vs LumB:\n p=', format(subset(r, contrast == "Her2 vs. LumB")$pval, digits=3, sci = TRUE), "; HR=", format(subset(r, contrast == "Her2 vs. LumB")$HR, digits=3, sci = FALSE), sep = "")
	  y.label1 <- 0.9
	  y.label2 <- 0.6
  }

  vardata = data.frame(sample_id_anonymized = km.data$sample_id_anonymized, value = vals);
  foo <- merge(km.data, vardata, by='sample_id_anonymized', sort=F);
  fit <- survfit(as.formula(paste('Surv(pfs, pfs_event)~', 'value')), data=foo);
  p <- ggsurvplot(fit, conf.int = FALSE, ggtheme = theme_classic(base_size = 19), palette=c('red', 'navyblue', 'green'), 
                  title=paste(feature, paste('(n=', length(vardata$value[!is.na(vardata$value)]), ')', sep='')),
                  legend.labs=names(table(foo$value)), legend.title = '',
                  ylab = "PFS probability", xlab = "Time (days)", 
                  pval=F, risk.table=T, fontsize=7, font.tickslab=22, size=2, censor.size=10);
  
  p$plot  <- p$plot + ggplot2::annotate("text", x = 1500, y = y.label1, label = label1, size = 8) # x and y coordinates of the text
  p$plot  <- p$plot + ggplot2::annotate("text", x = 1500, y = y.label2, label = label2, size = 8) # x and y coordinates of the text
  
  print(p)
}

dev.off()			   