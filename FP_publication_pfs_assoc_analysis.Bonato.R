## Author: Vini Bonato
## Purpose: BC Tempus manuscript
## Created on: 06/27/2023
## Last modified on: 06/27/2023

## Libraries

library(openxlsx)
library(survival)
library(coxphf)
library(qvalue)
library(forestplot)
library(plyr)
library(Cairo)

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
	univ.formulas <- sapply(covariates$variable, function(x) as.formula(paste('Surv(pfs, pfs_event)~', x)))
	univ.models   <- lapply(univ.formulas, function(x){
										try(coxphf(x, data = subset(dataset, !is.na(dataset[, as.character(as.vector(x))[3]])), maxstep = 0.01, maxit = 1000, firth = FALSE), silent = TRUE)
									 })
	univ.results <- lapply(univ.models,
						   function(x){
								if(class(x)[1] != "try-error"){
									variable.name   <- strsplit(as.character(x$formula), " ~ ")[[3]]
									overall.p.value <- signif(pchisq(as.numeric(-2 * x$loglik[1] + 2 * x$loglik[2]), df = x$df, lower.tail = F), digits = 2)
									x.coef          <- as.data.frame(coef(x))
									x.ci            <- as.data.frame(x$conf.int)
									res.out         <- NULL
									for(j in 1:length(x$coefficients)){
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
									#res.out$overall.p.value  <- overall.p.value
									#res.out$mlog10.overall.p <- -log10(overall.p.value)
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

workDir <- "C:/Users/Bonatv/OneDrive - Pfizer/PFS_Assoc_Features/Stats/"
dataDir <- "C:/Users/Bonatv/OneDrive - Pfizer/PFS_Assoc_Features/"

## Start
patient.info              <- read.xlsx("C:/Users/Bonatv/OneDrive - Pfizer/PFS_Assoc_Features/Supplementary_table_02.molecular_features_hallmark_gsva.xlsx", sheet = "patient clinical")

anal.mol.features.data    <- read.xlsx("C:/Users/Bonatv/OneDrive - Pfizer/PFS_Assoc_Features/Supplementary_table_02.molecular_features_hallmark_gsva.xlsx", sheet = "molecular features")
anal.mol.features.data    <- merge(anal.mol.features.data, patient.info, by = "sample_id", sort = FALSE, all = TRUE)
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
anal.mol.features.pre.data <- subset(anal.mol.features.data, prepost == "Pre")
anal.mol.features.pos.data <- subset(anal.mol.features.data, prepost == "Post")

hallmark.features.data    <- read.xlsx("C:/Users/Bonatv/OneDrive - Pfizer/PFS_Assoc_Features/Supplementary_table_02.molecular_features_hallmark_gsva.xlsx", sheet = "hallmark gsva")
hallmark.features.data    <- merge(hallmark.features.data, patient.info, by = "sample_id", sort = FALSE, all = TRUE)
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
hallmark.features.pre.data <- subset(hallmark.features.data, prepost == "Pre")
hallmark.features.pos.data <- subset(hallmark.features.data, prepost == "Post")

gen.alt.features.data    <- read.table("C:/Users/Bonatv/OneDrive - Pfizer/PFS_Assoc_Features/tempus_ms_genomic_features.txt")
gen.alt.features.data$sample_id <- row.names(gen.alt.features.data)
gen.alt.features.data    <- merge(gen.alt.features.data, patient.info, by = "sample_id", sort = FALSE, all = TRUE)
gen.alt.features.cov.mut <- data.frame(variable = c("MUT.ABCB1", "MUT.ABCC3", "MUT.ABL1", "MUT.AKT1", "MUT.ALK", "MUT.APC", "MUT.APOB", "MUT.ARID1A", "MUT.ARID1B", "MUT.ARID2", "MUT.ASXL1", "MUT.ATM", 
													"MUT.ATR", "MUT.ATRX", "MUT.BCL11B", "MUT.BCLAF1", "MUT.BCOR", "MUT.BCORL1", "MUT.BCR", "MUT.BRAF", "MUT.BRCA1", "MUT.BRCA2", "MUT.BRD4", "MUT.CARD11", 
													"MUT.CBFB", "MUT.CD40", "MUT.CDH1", "MUT.CDKN2A", "MUT.CEBPA", "MUT.CFTR", "MUT.CHD4", "MUT.CIC", "MUT.CIITA", "MUT.CREBBP", "MUT.CSF3R", "MUT.CUX1", 
													"MUT.DICER1", "MUT.DIS3L2", "MUT.DNM2", "MUT.DNMT3A", "MUT.DOT1L", "MUT.DPYD", "MUT.DYNC2H1", "MUT.EGFR", "MUT.EP300", "MUT.EPHB1", "MUT.EPHB2", 
													"MUT.ERBB2", "MUT.ERBB3", "MUT.ERBB4", "MUT.ERCC4", "MUT.ERCC6", "MUT.ESR1", "MUT.FANCA", "MUT.FANCM", "MUT.FAT1", "MUT.FGFR2", "MUT.FGFR3", "MUT.FLT3", 
													"MUT.FLT4", "MUT.FOXA1", "MUT.FOXQ1", "MUT.GALNT12", "MUT.GATA3", "MUT.GATA6", "MUT.GRIN2A", "MUT.HNF1A", "MUT.IFNAR2", "MUT.IKZF1", "MUT.IRS2", "MUT.JAK1", 
													"MUT.JAK3", "MUT.KAT6A", "MUT.KDM5A", "MUT.KDM6A", "MUT.KIF1B", "MUT.KIT", "MUT.KMT2A", "MUT.KMT2B", "MUT.KMT2C", "MUT.KMT2D", "MUT.KRAS", "MUT.LRP1B", 
													"MUT.MAP2K4", "MUT.MAP3K1", "MUT.MED12", "MUT.MEN1", "MUT.MKI67", "MUT.MLH1", "MUT.MSH3", "MUT.MTOR", "MUT.MYH11", "MUT.NCOR1", "MUT.NCOR2", "MUT.NF1", 
													"MUT.NOTCH1", "MUT.NOTCH2", "MUT.NOTCH3", "MUT.NRG1", "MUT.PALLD", "MUT.PBRM1", "MUT.PDGFRB", "MUT.PIK3CA", "MUT.PIK3R1", "MUT.PIK3R2", "MUT.PLCG2", 
													"MUT.PMS2", "MUT.POLE", "MUT.PREX2", "MUT.PTCH1", "MUT.PTEN", "MUT.PTPN13", "MUT.PTPRD", "MUT.RAD21", "MUT.RANBP2", "MUT.RB1", "MUT.RECQL4", "MUT.RET", 
													"MUT.ROS1", "MUT.RUNX1", "MUT.SETBP1", "MUT.SETD2", "MUT.SF3B1", "MUT.SLIT2", "MUT.SLX4", "MUT.SMARCA4", "MUT.SMO", "MUT.SPEN", "MUT.STAG2", "MUT.TAF1", 
													"MUT.TBX3", "MUT.TERT", "MUT.TET2", "MUT.TP53", "MUT.TSC1", "MUT.TSC2", "MUT.ZFHX3", "MUT.ZNF217", "MUT.ZNF750", "MUT.ZNRF3", "AMP.AURKA", "AMP.AXIN2", 
													"AMP.BRIP1", "AMP.CCND1", "AMP.CD79B", "AMP.CDKN1B", "AMP.CKS1B", "AMP.ELF3", "AMP.FGF3", "AMP.FGF4", "AMP.FGFR1", "AMP.FOXA1", "AMP.FRS2", "AMP.GNA13", 
													"AMP.GNAS", "AMP.H3.3A", "AMP.IDO1", "AMP.IKBKE", "AMP.KAT6A", "AMP.MCL1", "AMP.MDM2", "AMP.MDM4", "AMP.MTMR11", "AMP.MYC", "AMP.NDUFC2", "AMP.NTHL1", 
													"AMP.PAK1", "AMP.PIK3C2B", "AMP.RAD21", "AMP.RAD51C", "AMP.RECQL4", "AMP.RNF139", "AMP.RNF43", "AMP.RPS6KB1", "AMP.RSF1", "AMP.SPOP", "AMP.SRSF2", 
													"AMP.UBE2T", "AMP.YEATS4", "AMP.ZNF217", "AMP.ZNF750", "DEL.ALK", "DEL.APLNR", "DEL.CDKN2A", "DEL.CDKN2B", "DEL.CYP1B1", "DEL.CYP2D6", "DEL.EGFR", 
													"DEL.EPHB1", "DEL.FGF14", "DEL.MTAP", "DEL.NKX2.1", "DEL.NRG1", "DEL.PDE4D", "DEL.PRSS1"),
										type     = rep("mutation", 195))
gen.alt.features.cov      <- gen.alt.features.cov.mut
gen.alt.features.pre.data <- subset(gen.alt.features.data, prepost == "Pre")
gen.alt.features.pos.data <- subset(gen.alt.features.data, prepost == "Post")

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

write.xlsx(out.list, file = paste(workDir, "modelResults.xlsx", sep = ""), rowNames = FALSE, colNames = TRUE, keepNA = TRUE, na.string = "")

## Plots:

anal.mol.features.pre.plot <- anal.mol.features.pre.res
anal.mol.features.pre.plot <- subset(anal.mol.features.pre.plot, pval < 0.05 & variable.labels != "subtype_pam50 (Normal vs. LumB)")
anal.mol.features.pre.plot <- anal.mol.features.pre.plot[order(anal.mol.features.pre.plot$HR), ]
anal.mol.features.pre.plot$variable.labels <- c("PAM50 subtype (Luminal A vs. Luminal B)", "Paloma3 signature F10 (Estrogen response)", "NMF factor5 (Luminal A)",
												"PAM50 score (Luminal A)", "Paloma3 signature F2", "Paloma3 signature F5", "PAM50 score (Normal)", "NMF factor12 (HER2 expr)", 
												"MutSig (13 sigma)", "Proliferative index", "Paloma3 signature F9", "Paloma3 signature F8 (EMT)", "NMF factor11 (proliferation)",
												"Paloma3 signature F1 (MYC E2F activation)", "PAM50 score (Basal)", "PAM50 score (Her2)", "DDR germline mutation (Y vs. N)",
												"BRCA pathogenic mutation status (Y vs. N)")

CairoPDF(paste(workDir, "FP_publication_AnalMol_Pre.pdf", sep = ""), height = 8, width = 16)
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
               txt_gp    = fpTxtGp(label = gpar(cex = 2), xlab = gpar(cex = 2, fontface = "bold"), tick = gpar(cex = 1.7, fontface = "bold")),
               xlab      = "CPH model estimate: log(HR)")
dev.off()

anal.mol.features.pos.plot <- anal.mol.features.pos.res
anal.mol.features.pos.plot <- subset(anal.mol.features.pos.plot, pval < 0.05 & variable.labels != "subtype_pam50 (Normal vs. LumB)")
anal.mol.features.pos.plot <- anal.mol.features.pos.plot[order(anal.mol.features.pos.plot$HR), ]
anal.mol.features.pos.plot$variable.labels <- c("ATM germline mutation (Y vs. N)", "PAM50 score (Basal)", "PAM50 score (Luminal B)", "HRD index", 
												"FACETS LOH", "FACETS LST", "NMF factor5 (Luminal A)", "Paloma3 signature F4", "Paloma3 signature F6",
												"Paloma3 signature F8 (EMT)", "Paloma3 signature F10 (Estrogen response)", "Paloma3 signature F13",
												"Paloma3 signature F14 (Liver specific)")

CairoPDF(paste(workDir, "FP_publication_AnalMol_Post.pdf", sep = ""), height = 8, width = 16)
    forestplot(labeltext = list(anal.mol.features.pos.plot$variable.labels), 
               mean      = as.numeric(as.vector(anal.mol.features.pos.plot$est)),
               lower     = as.numeric(as.vector(anal.mol.features.pos.plot$est.LB95)),
               upper     = as.numeric(as.vector(anal.mol.features.pos.plot$est.UB95)),
               zero      = 0,
			   lwd.zero  = 1.2,
			   mar       = unit(c(5, 5, 5, 6), "mm"),
			   clip      = c(-5, 5), 
               xticks    = seq(-1, 5, 1),
			   grid      = structure(seq(-1, 5, 1), gp = gpar(col = "gray75", lty = 3, lwd = 0.5)),
			   col       = fpColors(box = "gray65", line = "gray65", summary = "gray65", zero = "black"),
               new_page  = FALSE,
               txt_gp    = fpTxtGp(label = gpar(cex = 2), xlab = gpar(cex = 2, fontface = "bold"), tick = gpar(cex = 1.7, fontface = "bold")),
               xlab      = "CPH model estimate: log(HR)")
dev.off()

hallmark.features.pre.plot <- hallmark.features.pre.res
hallmark.features.pre.plot <- subset(hallmark.features.pre.plot, pval < 0.05)
hallmark.features.pre.plot <- hallmark.features.pre.plot[order(hallmark.features.pre.plot$HR), ]
hallmark.features.pre.plot$variable.labels <- c("Estrogen response early", "Coagulation", "E2F targets", "Oxidative phosphorylation", 
												"MYC targets V2", "Protein secretion", "G2M checkpoint", "Peroxisome", "MYC targets V1", 
												"UV response up", "Cholesterol homeostasis", "PI3K AKT MTOR signaling", "Unfolded protein response", 
												"Glycolysis", "MTORC1 signaling")

CairoPDF(paste(workDir, "FP_publication_Hallmark_Pre.pdf", sep = ""), height = 8, width = 16)
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

hallmark.features.pos.plot <- hallmark.features.pos.res
hallmark.features.pos.plot <- subset(hallmark.features.pos.plot, pval < 0.05)
hallmark.features.pos.plot <- hallmark.features.pos.plot[order(hallmark.features.pos.plot$HR), ]
hallmark.features.pos.plot$variable.labels <- c("Estrogen response early", "Xenobiotic metabolism", "Bile acid metabolism", "Coagulation",
												"Estrogen response late", "TNFA signaling via NFKB")

CairoPDF(paste(workDir, "FP_publication_Hallmark_Post.pdf", sep = ""), height = 8, width = 16)
    forestplot(labeltext = list(hallmark.features.pos.plot$variable.labels), 
               mean      = as.numeric(as.vector(hallmark.features.pos.plot$est)),
               lower     = as.numeric(as.vector(hallmark.features.pos.plot$est.LB95)),
               upper     = as.numeric(as.vector(hallmark.features.pos.plot$est.UB95)),
               zero      = 0,
			   lwd.zero  = 1.2,
			   mar       = unit(c(5, 5, 5, 8), "mm"),
			   clip      = c(-5, 1.75), 
               xticks    = seq(-0.5, 0.5, 0.25),
			   grid      = structure(seq(-0.5, 0.5, 0.25), gp = gpar(col = "gray75", lty = 3, lwd = 0.5)),
			   col       = fpColors(box = "gray65", line = "gray65", summary = "gray65", zero = "black"),
               new_page  = FALSE,
               txt_gp    = fpTxtGp(label = gpar(cex = 2), xlab = gpar(cex = 2, fontface = "bold"), tick = gpar(cex = 1.7, fontface = "bold")),
               xlab      = "CPH model estimate: log(HR)")
dev.off()

gen.alt.features.pre.plot <- gen.alt.features.pre.res
gen.alt.features.pre.plot <- subset(gen.alt.features.pre.plot, pval < 0.05)
gen.alt.features.pre.plot <- gen.alt.features.pre.plot[order(gen.alt.features.pre.plot$HR), ]
gen.alt.features.pre.plot$variable.labels <- c("AKT1 (Mut vs. WT)", "BRCA2 (Mut vs. WT)", "KDM6A (Mut vs. WT)",
											   "NF1 (Mut vs. WT)", "PIK3CA (Mut vs. WT)", "TP53 (Mut vs. WT)",
											   "FGFR1 (Ampl. vs. WT)", "NTHL1 (Ampl. vs. WT)", "RAD51C (Ampl. vs. WT)",
											   "RNF43 (Ampl. vs. WT)", "UBE2T (Ampl. vs. WT)", "YEATS4 (Ampl. vs. WT)")

CairoPDF(paste(workDir, "FP_publication_GenAlterations_Pre.pdf", sep = ""), height = 8, width = 16)
    forestplot(labeltext = list(gen.alt.features.pre.plot$variable.labels), 
               mean      = as.numeric(as.vector(gen.alt.features.pre.plot$est)),
               lower     = as.numeric(as.vector(gen.alt.features.pre.plot$est.LB95)),
               upper     = as.numeric(as.vector(gen.alt.features.pre.plot$est.UB95)),
               zero      = 0,
			   lwd.zero  = 1.2,
			   mar       = unit(c(5, 5, 5, 8), "mm"),
			   clip      = c(-5, 2), 
               xticks    = seq(-2, 2, 0.5),
			   grid      = structure(seq(-2, 2, 0.5), gp = gpar(col = "gray75", lty = 3, lwd = 0.5)),
			   col       = fpColors(box = "gray65", line = "gray65", summary = "gray65", zero = "black"),
               new_page  = FALSE,
               txt_gp    = fpTxtGp(label = gpar(cex = 2), xlab = gpar(cex = 2, fontface = "bold"), tick = gpar(cex = 1.7, fontface = "bold")),
               xlab      = "CPH model estimate: log(HR)")
dev.off()


gen.alt.features.pos.plot <- gen.alt.features.pos.res
gen.alt.features.pos.plot <- subset(gen.alt.features.pos.plot, pval < 0.05)
gen.alt.features.pos.plot <- gen.alt.features.pos.plot[order(gen.alt.features.pos.plot$HR), ]
gen.alt.features.pos.plot$variable.labels <- c("EGFR (Mut vs. WT)", "ESR1 (Mut vs. WT)", "KAT6A (Mut vs. WT)", "KIF1B (Mut vs. WT)",
											   "KMT2D (Mut vs. WT)", "PLCG2 (Mut vs. WT)", "PTEN (Mut vs. WT)", "RECQL4 (Mut vs. WT)", 
											   "TET2 (Mut vs. WT)", "TP53 (Mut vs. WT)", "TSC1 (Mut vs. WT)", "MYC (Ampl. vs. WT)", 
											   "YEATS4 (Ampl. vs. WT)", "APLNR (Del vs. WT)")

CairoPDF(paste(workDir, "FP_publication_GenAlterations_Post.pdf", sep = ""), height = 8, width = 16)
    forestplot(labeltext = list(gen.alt.features.pos.plot$variable.labels), 
               mean      = as.numeric(as.vector(gen.alt.features.pos.plot$est)),
               lower     = as.numeric(as.vector(gen.alt.features.pos.plot$est.LB95)),
               upper     = as.numeric(as.vector(gen.alt.features.pos.plot$est.UB95)),
               zero      = 0,
			   lwd.zero  = 1.2,
			   mar       = unit(c(5, 5, 5, 8), "mm"),
			   clip      = c(-5, 2), 
               xticks    = seq(-1, 2, 0.5),
			   grid      = structure(seq(-1, 2, 0.5), gp = gpar(col = "gray75", lty = 3, lwd = 0.5)),
			   col       = fpColors(box = "gray65", line = "gray65", summary = "gray65", zero = "black"),
               new_page  = FALSE,
               txt_gp    = fpTxtGp(label = gpar(cex = 2), xlab = gpar(cex = 2, fontface = "bold"), tick = gpar(cex = 1.7, fontface = "bold")),
               xlab      = "CPH model estimate: log(HR)")
dev.off()




