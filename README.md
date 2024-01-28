## Overview

This repository contains the source code associated with the manuscript titled ""Real-world clinical genomics analysis revealed a dichotomy between ER-independent and ER-dependent mechanisms of drug resistance to CDK4/6i plus endocrine therapies in HR+/HER2- metastatic breast cancers" 

This project is licensed under the Apache License, Version 2.0 - see the [LICENSE](LICENSE) file for details.

### Contact
- Email: zhengyan.kan@pfizer.com

## Requirements

To run the code in this repository, you will need the following:

R version 4.1.0 or higher

## Descriptions

Order	Display items	Description	Author	Code name
1	Table 1	Cohort summary	George Kan	tempus_ms_cohort_summary_stats.R 
1	Figure 1b, d	Volcano plot showing molecular features with significant changes in Pre vs. Post	George Kan	tempus_ms_prepost_compare_feature.R 
2	Figure 1c, e	Forest plot showing molecular features with significant association with PFS	Vini Bonato	FP_publication_pfs_assoc_analysis.Bonato.R 
3	Figure 1f, h, j	Box plot showing distribution of key gene expressions in Pre vs. Post	George Kan	tempus_ms_prepost_compare_feature.R 
4	Figure 1g, I, k	KM plot showing PFS association of key gene expression (median-split)	George Kan	tempus_ms_pfs_assoc_feature.R 
5	Figure 1i	Bar plot showing top 7 genes with significant changes in GAF in Post vs. Pre	Ji Wen	Figure1i.stackplot.R 
6	Figure 2a	Complex heatmap showing distinct pattersn of molecular features among the four  integrative clusters IC1-4	Ji Wen	Figure2a.R 
7	Figure 2b	KM plot showing PFS association of four integrative clusters IC1-4	George Kan	tempus_ms_icluster_characterize_selected.R 
8	Figure 2c	Bar plot showing changes in the prevalence of IC1-4 in Pre vs. Post	George Kan	tempus_ms_icluster_characterize_selected.R 
9	Figure 2d	Box plot showing distribution of key gene expressions vs. IC1-4	George Kan	tempus_ms_icluster_characterize_selected.R 
10	Figure 2e-f	Bar plot showing distribution of mutation and copy number amplification status vs. IC1-4	George Kan	tempus_ms_icluster_characterize_selected.R 
11	Figure 3a	Complex heatmap showing top ranked hits from DE comparisons of gene CSE signatures vs. IC1-4 and PAM50 subtypes	George Kan	tempus_ms_DE_gsva_analysis.R 
12	Figure 3b-d	Volcano plot showing top hits from DE comparisons of IC1, Basal and Her2 vs. others 	George Kan	tempus_ms_DE_gsva_analysis.R 
13	Figure 3e-h	Box plots showing the distributions of gene CSE signatures vs. IC1-4, PAM50 subtypes and Pre/Post	George Kan	tempus_ms_expr_gsva_multipanel_analysis.R 
14	Figure 4b-c	Box plots showing the distribution of Pseudotime vs. Pre/Post, IC1-4 and PAM50 subtypes	George Kan	tempus_ms_pseudotime_characterize_selected.R 
15	Figure 4e	Complex heatmap showing correlation of molecular features vs. Pseudotime	George Kan	tempus_ms_pseudotime_covariate_analysis.R 
16	Figure 4g	Bar plot showing the distribution of EPG branches vs. IC1-4	George Kan	tempus_ms_EPG_branch_characterize_selected.R 
17	Figure 4h	Box plot showing distribution of key gene expressions vs. EPG branches	George Kan	tempus_ms_EPG_branch_characterize_selected.R 
18	Figure 4i-j	Bar plot showing distribution of mutation and copy number amplification status vs. EPG branches	George Kan	tempus_ms_EPG_branch_characterize_selected.R 
19	Figure 4k-n	Box plot showing the trend of Pseudotime and key gene exprssion vs. EPG trajectory nodes	George Kan	tempus_ms_epg_pseudotime_compare.R 
20	Figure 5b-h	Box plot showing the distribution of Dependency scores vs. IC1-4, Pre/Post and key mutation status	George Kan	tempus_ms_enrs_score_multiplanel_analysis.R 
21	Supplementary Figure 1a	Bar plot showing changes in the prevalence of PAM50 subtypes vs. Pre/Post	George Kan	tempus_ms_prepost_compare_feature.R 
22	Supplementary Figure 1b-f	Box plot showing distribution of molecular features in Pre vs. Post	George Kan	tempus_ms_prepost_compare_feature.R 
23	Supplementary Figure 2a-i	KM plot showing PFS association of molecular features (median-split)	George Kan	tempus_ms_pfs_assoc_feature.R 
24	Supplementary Figure 3a-f	Box plot showing distribution of gene expression signatures in Pre vs. Post	George Kan	tempus_ms_prepost_compare_gsva.R 
25	Supplementary Figure 3g-i	KM plot showing PFS association of gene expression signatures (median-split)	George Kan	tempus_ms_pfs_assoc_gsva.R 
26	Supplementary Figure 4a-b	Oncoprint map showing changes in GAF in Pre vs. Post	Ji Wen	Figure.Supp4a.landscape.R 
27	Supplementary Figure 4c-d	Distribution of gene expression signatures vs. mutation or copy number amplification status	George Kan	tempus_ms_assoc_genomic_vs_molecular_feature.R 
28	Supplementary Figure 6a	Forest plot showing genomic alteration events with significant PFS association. 	Vini Bonato	FP_publication_pfs_assoc_analysis.Bonato.R 
29	Supplementary Figure 6b-d	KM plot showing PFS association of mutation and copy number amplification	George Kan	tempus_ms_pfs_assoc_genomic_feature.R 
30	Supplementary Figure 7a	Bar plot showing the distribution of IC1-4 vs. PAM50 subtypes	George Kan	tempus_ms_icluster_characterize_selected.R 
31	Supplementary Figure 7b-c	Box plot showing the distribution of molecular features vs. IC1-4	George Kan	tempus_ms_icluster_characterize_selected.R 
32	Supplementary Figure 8a-d	Box plot showing the distribution of gene CSE vs. PAM50 subtypes	George Kan	tempus_ms_expr_gsva_multipanel_analysis.R 
33	Supplementary Figure 9a	KM plot showing PFS association of Pseudotime (median-split) at baseline	George Kan	tempus_ms_pseudotime_characterize_selected.R 
34	Supplementary Figure 9b, h	Boxplot showing the distribution of Pseudotime vs. IC-PAM50 subtypes and mutation statuses	George Kan	tempus_ms_pseudotime_characterize_selected.R 
35	Supplementary Figure 9bc-g	Scatter plot showing the correlation of Pseudotime vs. gene expression	George Kan	tempus_ms_pseudotime_characterize_selected.R 
36	Supplementary Figure 10a-b	Bar plot showing the distribution of EPG branches vs. PAM50 subtypes and Pre/Post	George Kan	tempus_ms_EPG_branch_characterize_selected.R 
37	Supplementary Figure 10c-e	Box plot showing distribution of molecular features and Pseudotime vs. EPG branches	George Kan	tempus_ms_epg_pseudotime_compare.R 
38	Supplementary Figure 10f-g	Box plot showing the trend of molecular features vs. EPG trajectory nodes	George Kan	tempus_ms_epg_pseudotime_compare.R ![image](https://github.com/kanzhengyan/Tempus_BC_manuscript_code_share/assets/41809702/5bac01ea-0694-4d66-b32e-21f4f7f350bc)

