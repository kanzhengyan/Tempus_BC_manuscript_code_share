rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE)
options(width=200)
#setwd("Figures/ForQC.20231018")
library(ComplexHeatmap)
library(reshape2)
library(tidyr)
library(circlize)

molfeatfile <- 'Input/Tempus_molecular_features_extended_v3.txt';
outprefix2=paste0("Figure2a.v",format(Sys.time(), "%Y%m%d"))

####################################################
mfdata <- read.table(molfeatfile, header=T, sep='\t');
mfdata = mfdata[ mfdata$icluster%in%c("IC1","IC2","IC3","IC4"),]
mfdata = mfdata[ order(mfdata$icluster) ,]

Matrix.data1=mfdata[,grep("EXPR",colnames(mfdata),value=T)]
Matrix.data1=mfdata[,c('EXPR.ESR1',"EXPR.CCNE1",'EXPR.CCNE2','EXPR.ERBB2','EXPR.CDK6','EXPR.PGR')]
colnames(Matrix.data1)=gsub("EXPR.","EXPR:",colnames(Matrix.data1))
#Matrix.data1 =  scale(  Matrix.data1) 

Matrix.data2=mfdata[,grep("h.HALLMARK_",colnames(mfdata),value=T)]
Matrix.data2=mfdata[,c('h.HALLMARK_E2F_TARGETS',"h.HALLMARK_MYC_TARGETS_V1",'h.HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',
						'h.HALLMARK_APOPTOSIS','h.HALLMARK_INFLAMMATORY_RESPONSE','h.HALLMARK_WNT_BETA_CATENIN_SIGNALING',
						'h.HALLMARK_ESTROGEN_RESPONSE_EARLY')]
colnames(Matrix.data2)=gsub("h.HALLMARK_","",colnames(Matrix.data2))
colnames(Matrix.data2)[colnames(Matrix.data2)=='EPITHELIAL_MESENCHYMAL_TRANSITION']="EMT"
colnames(Matrix.data2)[colnames(Matrix.data2)=='ESTROGEN_RESPONSE_EARLY']="Estrogen response early"
colnames(Matrix.data2)[colnames(Matrix.data2)=='WNT_BETA_CATENIN_SIGNALING']="Wnt/beta-catenin signaling"
colnames(Matrix.data2)[colnames(Matrix.data2)=='INFLAMMATORY_RESPONSE']="Inflammatory response"
colnames(Matrix.data2)[colnames(Matrix.data2)=='E2F_TARGETS']="E2F targets"
colnames(Matrix.data2)[colnames(Matrix.data2)=='MYC_TARGETS_V1']="MYC targets"
colnames(Matrix.data2)[colnames(Matrix.data2)=='APOPTOSIS']="Apoptosis"

Matrix.data3=mfdata[,grep("pam50_cor",colnames(mfdata),value=T)]
colnames(Matrix.data3)=gsub("pam50_cor","PAM50:cor",colnames(Matrix.data3))
Matrix.data3.z= scale( Matrix.data3 )

Matrix.data4=mfdata[,c("estimate_stromalscore","estimate_immunescore" ,"estimate_tumor_purity" )]
colnames(Matrix.data4)=gsub("estimate_","ESTIMATE:",colnames(Matrix.data4))
#Matrix.data4=mfdata[,grep("est",colnames(mfdata),value=T)]
Matrix.data4.z = scale( Matrix.data4 )

Matrix.data5=mfdata[,grep("Paloma3_",colnames(mfdata),value=T)]
sel.col=c("Paloma3_F1_MYC_E2F_activation","Paloma3_F7_IFNG_response","Paloma3_F8_EMT", "Paloma3_F10_estrogen_response"   )
#Matrix.data5=cbind(Matrix.data5[,colnames(Matrix.data5)%in%sel.col], Paloma3_other_factors=rowSums(Matrix.data5[,!colnames(Matrix.data5)%in%sel.col]))
Matrix.data5=  Matrix.data5[,colnames(Matrix.data5)%in%sel.col]  
colnames(Matrix.data5)=gsub("Paloma3_","",colnames(Matrix.data5))
Matrix.data5.z= scale( Matrix.data5 )
Matrix.data5 [Matrix.data5<0]=0
#Matrix.data5 <-  apply(Matrix.data5,1,function(x){as.numeric(x)/sum(as.numeric(x))})
#Matrix.data5[is.na(Matrix.data5)]=0

grep("nmf_factor",colnames(mfdata),value=T)
sel.col=c("nmf_factor2_basal","nmf_factor3_ER_expr","nmf_factor5_lumA","nmf_factor6_immune","nmf_factor8_ER_response", "nmf_factor9_EMT" ,"nmf_factor11_proliferation")
Matrix.data6=  mfdata[,colnames(mfdata)%in%sel.col] 
Matrix.data6.z= scale( Matrix.data6 )
#Matrix.data6 <-  apply(Matrix.data6,1,function(x){as.numeric(x)/sum(as.numeric(x))})
Matrix.data6[is.na(Matrix.data6)]=0

grep("tmb",colnames(mfdata),value=T)
sel.col=c("proliferative_index","cyt_score" ,"hrd_index","tmb","tmb_norm")
Matrix.data7=  mfdata[,colnames(mfdata)%in%sel.col] 
Matrix.data7.z= scale( Matrix.data7 )

grep("mutsig_\\d",colnames(mfdata),value=T)
sel.col=grep("mutsig_\\d",colnames(mfdata),value=T)
Matrix.data8=  mfdata[,colnames(mfdata)%in%sel.col] 
Matrix.data8.z= scale( Matrix.data8 )
Matrix.data8[is.na(Matrix.data8)]=0

grep("pre",colnames(mfdata),value=T)
sel.col=c("pre_post","histology","procedure_type","pfs1","pfs1_censored","pfs2","tissue_site","tissue_origin","subtype_pam50","MUT.RB1","MUT.ESR1","MUT.TP53","MUT.PIK3CA","AMP.MYC","AMP.CCND1","DEL.PTEN","DEL.RB1","icluster"  )
Matrix.data9=  mfdata[,colnames(mfdata)%in%sel.col] 
for(gene in grep("^MUT.|^AMP.|^DEL.",colnames(Matrix.data9),value=TRUE) )
{	Matrix.data9[,gene]=ifelse(Matrix.data9[,gene]==1,'Yes','No')	}
Matrix.data9$tissue_origin[grep("^Lymph node",Matrix.data9$tissue_origin)]="Lymph nodes"
Matrix.data9$tissue_origin[grep("lung",Matrix.data9$tissue_origin,ignore.case=T )]="Lung"
Matrix.data9$tissue_origin[Matrix.data9$tissue_origin=="Pelvic bones, sacrum, coccyx and associated joints"]="Pelvis"
#Matrix.data9$tissue_origin[!Matrix.data9$tissue_origin%in%c("Liver","Breast","Lymph nodes","Thorax","Skin of trunk","Pelvis","Lung")]="Others"
Matrix.data9$tissue_origin[!Matrix.data9$tissue_origin%in%c("Liver","Breast","Lymph nodes")]="Others"
sort(table(Matrix.data9$tissue_origin))

BaseColors = c("turquoise","blue","brown","yellow","green","red","black","pink","magenta",
                "purple","greenyellow","tan","salmon","cyan", "midnightblue", "lightcyan",
                "grey60", "lightgreen", "lightyellow", "royalblue", "darkred", "darkgreen",
                "darkturquoise", "darkgrey",
                "orange", "darkorange", "white", "skyblue", "saddlebrown", "steelblue", 
                "paleturquoise", "violet", "darkolivegreen", "darkmagenta" );
	
columnha <-HeatmapAnnotation(icluster = Matrix.data9$icluster,
					"Treatment status" =  Matrix.data9$pre_post  ,
					"PAM50 subtype" = Matrix.data9$subtype_pam50,
					#"Histology" = Matrix.data9$histology,
					#"procedure_type" = Matrix.data9$procedure_type,
					#"PFS" = Matrix.data9$pfs1,
					#"pfs1_censored" = Matrix.data9$pfs1_censored,
					#"pfs2" = Matrix.data9$pfs2,
					#"tissue_site" = Matrix.data9$tissue_site,
					"Tissue origin" = Matrix.data9$tissue_origin,				
                    col = list(icluster =c("IC1"="#E41A1C","IC2"="#377EB8","IC3"="#4DAF4A","IC4"="#984EA3"),
					  "Treatment status" = c("Pre" =  "lightgrey","Post"="gray20"),
					  "PAM50 subtype" =c("Basal"="purple","Her2" = "green","LumA"="yellow","LumB"="orange", "Normal"="gray"),
					  #"Histology" = c("Infiltrating duct carcinoma"="red","Lobular carcinoma" = "green","Adenocarcinoma"="blue","Carcinoma, metastatic"="orange" ),
					  #"PFS" = colorRamp2(c(0,  1000), c(  "white", "darkblue")),
					  "Tissue origin"= c("Breast"=BaseColors[1],"Skin of trunk"=BaseColors[2],"Liver"=BaseColors[3],"Lymph nodes"=BaseColors[4],"Others"=BaseColors[5],"Pelvis"=BaseColors[6],"Thorax"=BaseColors[7],"Lung"= BaseColors[8])
                    ),
                    annotation_height = unit(c(4,4,4,4,4,4,4,4,4,4,4,4,4,4,4), "mm"),
                    annotation_legend_param = list(  #Time=list(title = "Time"),#annotation_height = unit(c(4,4,4), "mm"),
												     #PFS=list(title = "PFS",at = seq(0,1000,250),labels = c("0","250", "500","750", ">1000")),
													 "PAM50 subtype" = list(title = "PAM50 subtype")  ),
					annotation_name_side = "right",
					show_legend=c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE)
)
	
ht_list <-Heatmap(as.matrix(t(Matrix.data1)),  column_split=mfdata[,'icluster'], column_title_gp = gpar(fontsize = 25,fontface= "bold"), column_dend_height = unit(2, "cm"),
		top_annotation=columnha, show_heatmap_legend = FALSE,rect_gp = gpar(type = "none"),height = -0.5,
		show_column_names = FALSE, cluster_rows = FALSE, cluster_columns = TRUE,cluster_column_slices = FALSE,show_row_names = FALSE ) %v%
HeatmapAnnotation(  "MUT:RB1" = Matrix.data9$MUT.RB1,
					"MUT:PIK3CA" = Matrix.data9$MUT.PIK3CA,
					"MUT:ESR1" = Matrix.data9$MUT.ESR1,
					"MUT:TP53" = Matrix.data9$MUT.TP53,
					"AMP:CCND1" = Matrix.data9$AMP.CCND1,
					#DEL.PTEN = Matrix.data9$DEL.PTEN,
					#DEL.RB1 = Matrix.data9$DEL.RB1,
					col = list( "MUT:RB1"=c("No"="gray75","Yes"="black"),
					  "MUT:PIK3CA"=c("No"="gray75","Yes"="black"),
					  "MUT:ESR1"=c("No"="gray75","Yes"="black"),
					  "MUT:TP53"=c("No"="gray75","Yes"="black"),
					  "AMP:CCND1"=c("No"="gray75","Yes"="red")
					  #DEL.PTEN=c("0"="gray75","1"="black"),
					  #DEL.RB1=c("0"="gray75","1"="black"),
					  ),
                    annotation_height = unit(c(4,4,4,4,4,4,4,4,4,4,4,4,4,4,4), "mm"),
                    annotation_legend_param = list("MUT:RB1"=list(title = "MUT"),
													 "AMP:CCND1" = list(title = "AMP")  ),
					annotation_name_side = "right",
					show_legend=c(TRUE,FALSE,FALSE,FALSE,TRUE)
)%v%
HeatmapAnnotation(	TMB = anno_barplot(Matrix.data7$tmb, ylim = c(0, max(Matrix.data7$tmb, na.rm = TRUE)), axis = TRUE,axis_param = list(side = "left")),
					"Proliferative index" = anno_barplot(Matrix.data7$proliferative_index, ylim = c(0, max(Matrix.data7$proliferative_index, na.rm = TRUE)), extend = 0, axis = TRUE,axis_param = list(side = "left")),
					"Cyt score" = anno_barplot(Matrix.data7$cyt_score, ylim = c(0, max(Matrix.data7$cyt_score, na.rm = TRUE)), axis = TRUE,axis_param = list(side = "left")),
					#hrd_index = anno_barplot(Matrix.data7$hrd_index, ylim = c(0, max(Matrix.data7$hrd_index, na.rm = TRUE)), axis = TRUE,axis_param = list(side = "left")),
					annotation_height = unit(c(8,8,8), "mm"),annotation_name_side = "right")%v%
HeatmapAnnotation("Paloma3 factors" = anno_barplot(as.matrix(  (Matrix.data5) ), show_legend=TRUE,
				gp = gpar(fill = c("green", "yellow","red", "blue","magenta","gray"),col= c("green", "yellow","red", "blue","magenta","gray")),bar_width = 1),height = unit(2, "cm"))%v%
Heatmap(as.matrix(t(Matrix.data4.z)),row_dend_side = "left",row_names_side ="right", 
		show_column_names = FALSE,cluster_columns = TRUE,cluster_column_slices = FALSE,name="ESTIMATE z-score") %v%
Heatmap(as.matrix(t(Matrix.data3.z )),row_dend_side = "left",row_names_side ="right", 
		show_column_names = FALSE,cluster_columns = TRUE,cluster_column_slices = FALSE,name="PAM50 z-score") %v%
Heatmap(as.matrix(t(Matrix.data2)),row_dend_side = "left",row_names_side ="right", 
		show_column_names = FALSE,cluster_columns = TRUE,cluster_column_slices = FALSE,name="GSVA") %v%
Heatmap(as.matrix(t(Matrix.data1)),row_dend_side = "left",row_names_side ="right",column_split=mfdata[,'icluster'], column_title_gp = gpar(fontsize = 23), column_dend_height = unit(2, "cm"),
		show_column_names = FALSE,  cluster_columns = TRUE,cluster_column_slices = FALSE,name="EXPR")  
#Heatmap(as.matrix(t(Matrix.data5)),row_dend_side = "left",row_names_side ="right", col= colorRamp2(c(0, max(Matrix.data5)), c("white", "red")),
#		show_column_names = FALSE,row_names_gp = gpar(fontsize=8),cluster_columns = TRUE,cluster_column_slices = FALSE,name="Paloma3NMF") %v%
#Heatmap(as.matrix(t(Matrix.data6)),row_dend_side = "left",row_names_side ="right", col= colorRamp2(c(0, max(Matrix.data6)), c("white", "red")),
#		show_column_names = FALSE,row_names_gp = gpar(fontsize=8),cluster_columns = TRUE,cluster_column_slices = FALSE,name="denovoNMF") %v%
#HeatmapAnnotation(DeNovoNMF = anno_barplot(as.matrix(  (Matrix.data6) ), 
#				gp = gpar(fill = BaseColors,col= BaseColors),bar_width = 1),height = unit(2, "cm"))%v%
#Heatmap(as.matrix(t(Matrix.data8)),row_dend_side = "left",row_names_side ="right", col= colorRamp2(c(0, max(Matrix.data8)), c("white", "red")),
#		show_column_names = FALSE,row_names_gp = gpar(fontsize=8),cluster_columns = TRUE,cluster_column_slices = FALSE,name="sigMA") %v%
#HeatmapAnnotation(SigMA = anno_barplot(as.matrix(  (Matrix.data8) ), 
#				gp = gpar(fill = BaseColors,col= BaseColors),bar_width = 1),height = unit(2, "cm")) %v%
				  
lgd_Matrix.data5 = Legend(labels = colnames(Matrix.data5), title = "Paloma3 factors", legend_gp = gpar(fill = c("green", "yellow","red", "blue","magenta","gray")  ) )
#lgd_Matrix.data6 = Legend(labels = colnames(Matrix.data6), title = "DeNovoNMF", legend_gp = gpar(fill = BaseColors  ) )
#lgd_Matrix.data8 = Legend(labels = colnames(Matrix.data8), title = "SigmA", legend_gp = gpar(fill = BaseColors  ) )

#dev.off()
pdf(  paste0(outprefix2,".pdf") ,20,10,useDingbats=F)
 
draw(ht_list, merge_legend = TRUE ,
    heatmap_legend_list = list(lgd_Matrix.data5),
    annotation_legend_side="bottom",heatmap_legend_side = "bottom" 
)
dev.off() 
