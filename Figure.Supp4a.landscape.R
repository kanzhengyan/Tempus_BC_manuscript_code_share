rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE)
options(width=200)
#setwd("Figures/ForQC.20231018")

options(stringsAsFactors = FALSE)
library(ComplexHeatmap)
library(reshape2)
library(tidyr)
library(ggplot2)

load("/LJ/CompBio2/share/data/omics/TEMPUS/Pfizer_Pancan-CDK4-6/Pfizer_Advanced_Molecular_Data_Package/TEMPUS.20230527.RData")
 
names(TEMPUS)
str(TEMPUS$Annotation$Sample)

table(TEMPUS$Annotation$Sample$PrePostCDKi)
#Post  Pre
# 227  200
Mutation=TEMPUS$Mutation$Somatic[order(TEMPUS$Mutation$Somatic$sample_id),]
#write.table(SMC_mutect2,file="SMC_mutect2.all.txt",sep="\t")
head(Mutation, 100)[1:10,1:10]

CNgaincut=7 
CNlosscut=0
#outprefix1=paste0("Figure1B.plot_dat_new.",CNgaincut,CNlosscut)
outprefix2=paste0("Figure.Supp4a.landscape.v",format(Sys.time(), "%Y%m%d"))

#################################################################################################################################################################################################################################################################################################

filter1=read.table("Input/tempus_ms_genomic_genelist.txt")
filter1=filter1[,1]
filter2=read.csv("Input/1.2.PostvsPre.TempusCN.stackedbar.70.csv",check.names=FALSE,row.names=1)
filter2=rownames(filter2)[filter2$'p-value'<0.05 & filter2$'# Altered Post'>0.05*227 ]
plot_genes=intersect(filter2,filter1) 
plot_genes=setdiff(plot_genes, "YEATS4")
#plot_genes= unique(c("ESR1", "RB1", "AKT1", "PIK3CA","CDKN1A", "PTEN", "TP53", "EGFR", "ERBB2", "KRAS", "CDKN2A","FGFR1", "FGFR2", "AURKA", "MAP2K4", "CDK4",   "CCND1",  "CDKN2A", "MYC", "ERBB2", "PTEN") )
#length(plot_genes)
 
length(table(Mutation$sample_id ) )
#426
plot_samples= sort( unique(TEMPUS$Annotation$Sample$ID ) )  
table(gsub(".*_","",plot_samples))
# Post  Pre
# 227  200

###cnv
CNV_master=read.csv("Input/master_table_CNV.v1.csv")
plot_dat_cnv=CNV_master[ CNV_master$gene_symbol%in%plot_genes, ]
table(plot_dat_cnv$paired_mutation_status,plot_dat_cnv$pre_post_paired)
table(plot_dat_cnv$paired_mutation_status,plot_dat_cnv$pre_or_post)
plot_dat_cnv$paired_mutation_status[is.na(plot_dat_cnv$paired_mutation_status)]="unpaired"

plot_dat_cnv=plot_dat_cnv[, c("sample_id","gene_symbol","mut_effect","paired_mutation_status"  )]
plot_dat_cnv$mut_effect=paste0("CN_",plot_dat_cnv$mut_effect)
plot_dat_cnv$tag=paste(plot_dat_cnv$gene_symbol,plot_dat_cnv$sample_id)
which(duplicated(plot_dat_cnv[,1:3]))

###somatic mut
plot_dat_mut=Mutation[(Mutation$gene_canonical_name %in% plot_genes &  Mutation$sample_id  %in% plot_samples)  ,   ]
table(plot_dat_mut$Removed_after_manual )
dim(plot_dat_mut)
plot_dat_mut=plot_dat_mut[plot_dat_mut$Removed_after_manual=='',]
dim(plot_dat_mut)
plot_dat_mut$GENE_NAME=plot_dat_mut$gene_canonical_name
plot_dat_mut$SAMPLE_ID=plot_dat_mut$sample_id
plot_dat_mut$CONSEQUENCE=plot_dat_mut$mutation_effect
plot_dat_mut$CONSEQUENCE=gsub("&.*","",plot_dat_mut[,"CONSEQUENCE"])
plot_dat_mut$CONSEQUENCE=gsub("_variant","",plot_dat_mut[,"CONSEQUENCE"])
plot_dat_mut[plot_dat_mut$CONSEQUENCE%in%c("exon_loss","initiator_codon","upstream_gene"),]
plot_dat_mut$CONSEQUENCE=gsub( "disruptive_inframe_deletion", "inframe_indel", plot_dat_mut$CONSEQUENCE  )
plot_dat_mut$CONSEQUENCE=gsub( "disruptive_inframe_insertion", "inframe_indel", plot_dat_mut$CONSEQUENCE  )
plot_dat_mut$CONSEQUENCE=gsub( "inframe_insertion", "inframe_indel", plot_dat_mut$CONSEQUENCE  )
plot_dat_mut$CONSEQUENCE=gsub( "inframe_deletion", "inframe_indel", plot_dat_mut$CONSEQUENCE  )
plot_dat_mut$CONSEQUENCE=gsub( "splice_acceptor", "splice", plot_dat_mut$CONSEQUENCE  )
plot_dat_mut$CONSEQUENCE=gsub( "splice_donor",    "splice", plot_dat_mut$CONSEQUENCE  )
plot_dat_mut$CONSEQUENCE=gsub( "splice_region",    "splice", plot_dat_mut$CONSEQUENCE  )
plot_dat_mut$CONSEQUENCE=gsub( "exon_loss",    "frameshift", plot_dat_mut$CONSEQUENCE  )
plot_dat_mut$CONSEQUENCE=gsub( "initiator_codon",    "start_lost", plot_dat_mut$CONSEQUENCE  )
plot_dat_mut=plot_dat_mut[!plot_dat_mut$mutation_effect%in%c("upstream_gene_variant"),]
table(plot_dat_mut$CONSEQUENCE) 
plot_dat_mut[plot_dat_mut$CONSEQUENCE%in%c("exon_loss","initiator_codon","upstream_gene"),]

plot_dat_mut$tag1=paste(plot_dat_mut$GENE_NAME,plot_dat_mut$SAMPLE_ID,plot_dat_mut$CONSEQUENCE)
tmpdata=plot_dat_mut
table(tmpdata$paired_mutation_status,tmpdata$pre_post_paired)
#             No Yes
#  Acquired    0  21
#  Lost        0   5
#  Maintained  0  44
table(tmpdata$paired_mutation_status,tmpdata$pre_or_post)
#             Post Pre
#  Acquired     21   0
#  Lost          0   5
#  Maintained   22  22
tmpdata$paired_mutation_status[is.na(tmpdata$paired_mutation_status)]="unpaired"
plot_dat_mut_status= aggregate(x = tmpdata["paired_mutation_status"],   
                       by = tmpdata[c("SAMPLE_ID", "GENE_NAME")],
                       FUN = function(market.values){
                               paste(market.values,collapse=";")
                             })
plot_dat_mut_status$tag=paste(plot_dat_mut_status$GENE_NAME,plot_dat_mut_status$SAMPLE_ID)
plot_dat_mut= aggregate(x = tmpdata["CONSEQUENCE"],   
                       by = tmpdata[c("SAMPLE_ID", "GENE_NAME")],
                       FUN = function(market.values){
                               paste(market.values,collapse=";")
                             })
#plot_dat_mut=rbind(plot_dat_mut,c("BRO7F-014_PD","RB1","splice") )
plot_dat_mut$tag=paste(plot_dat_mut$GENE_NAME,plot_dat_mut$SAMPLE_ID)

###germline
#write.csv(SMC_palbo$Germline$SMC_germline[SMC_palbo$Germline$SMC_germline$Gene%in%c("BRCA1","BRCA2"),],file="germline.all.csv")
names(TEMPUS$Germline)
table (TEMPUS$Germline$functional_impact ) 

Germline_From_WonChul=read.delim("Input/g_molecular_master_file_filtered_pathogenic_germline_revised.txt")
 dim(Germline_From_WonChul)
# [1] "variant_type_code"               "variant_type"                    "result_orig_unit"                "secondary_gene_canonical_system" "secondary_gene_canonical_code"   "secondary_gene_canonical_name"
# [7] "secondary_chromosome"            "position_2"                      "chrom_pos_2"                     "log_odds_ratio"                  "copy_number"                     "topology"
#[13] "rearr_num_reads"                 "effect"                          "fusion_caller"                   "tmb_percentile"

#plot_dat_gBRCA=TEMPUS$Germline[ TEMPUS$Germline$functional_impact =='P', ]
plot_dat_gBRCA=cbind(sample_id=TEMPUS$Annotation$Sample$ID[match(Germline_From_WonChul$analysis_id,TEMPUS$Annotation$Sample$DNA_analysis_id)],Germline_From_WonChul)
 
plot_dat_gBRCA=plot_dat_gBRCA[, c("sample_id","gene_canonical_name"  )]
plot_dat_gBRCA=cbind(plot_dat_gBRCA,germline="germline")

plot_dat_gBRCA$tag=paste(plot_dat_gBRCA$gene_canonical_name,plot_dat_gBRCA$sample_id)

###fusion
SV_master_raw=read.csv("SV_master_raw/master_table_SV.csv")
table(SV_master_raw$paired_mutation_status,SV_master_raw$pre_post_paired)
table(SV_master_raw$paired_mutation_status,SV_master_raw$pre_or_post)
SV_master_raw$paired_mutation_status[is.na(SV_master_raw$paired_mutation_status)]="unpaired"

tmp1=SV_master_raw[,c("sample_id","gene_canonical_name","paired_mutation_status")]
tmp2=SV_master_raw[,c("sample_id","secondary_gene_canonical_name","paired_mutation_status")]
colnames(tmp2)=c("sample_id","gene_canonical_name","paired_mutation_status")
SV_master=unique( rbind( tmp1,tmp2))
plot_dat_fusion=SV_master[SV_master$gene_canonical_name!='',]
colnames(plot_dat_fusion)=c("sample","gene","paired_mutation_status")

dim(plot_dat_fusion)
plot_dat_fusion$gene=gsub(":.*","",plot_dat_fusion$gene)

plot_dat_fusion[plot_dat_fusion$gene=="ESR1",]
plot_dat_fusion=unique(plot_dat_fusion)
dim(plot_dat_fusion)
plot_dat_fusion=plot_dat_fusion[plot_dat_fusion$gene%in%plot_genes,]
dim(plot_dat_fusion)
plot_dat_fusion=cbind(plot_dat_fusion,fusion="fusion")
plot_dat_fusion$tag=paste(plot_dat_fusion$gene,plot_dat_fusion$sample)
which(duplicated(plot_dat_fusion[,1:2]))

###mege all data togethoer
plot_dat_new=merge(plot_dat_cnv,plot_dat_mut, by="tag",all=T)
plot_dat_new=cbind( tag=plot_dat_new$tag, Gene=gsub(" .*","",plot_dat_new$tag), Sample=gsub(".* ","",plot_dat_new$tag) ,  CN=plot_dat_new$mut_effect, Mut=as.character(plot_dat_new$CONSEQUENCE) )
plot_dat_new=merge(plot_dat_new,plot_dat_gBRCA, by="tag",all=T)
plot_dat_new=cbind( tag=plot_dat_new$tag, Gene=gsub(" .*","",plot_dat_new$tag), Sample=gsub(".* ","",plot_dat_new$tag) ,  CN=plot_dat_new$CN, Mut= plot_dat_new$Mut , Germline=plot_dat_new$germline )
plot_dat_new=merge(plot_dat_new,plot_dat_fusion, by="tag",all=T)
plot_dat_new=cbind( tag=plot_dat_new$tag, Gene=gsub(" .*","",plot_dat_new$tag), Sample=gsub(".* ","",plot_dat_new$tag) ,  CN=plot_dat_new$CN, Mut= plot_dat_new$Mut , Germline=plot_dat_new$Germline, Fusion=plot_dat_new$fusion )
plot_dat_new=data.frame(plot_dat_new)
plot_dat_new$MUT=apply(plot_dat_new[,4:7], 1, 
              function(x) {
                x = x[x != ""]
				x = x[!is.na(x)]
                y=paste(unique(x), collapse = ";")
				y=unlist(strsplit(y,";") )
				paste(unique(y), collapse = ";")
              }
)
plot_dat_new$MUTN=apply(plot_dat_new[,4:7], 1, 
              function(x) {
                x = x[x != ""]
				x = x[!is.na(x)]
                y=paste(unique(x), collapse = ";")
				y=unlist(strsplit(y,";") )
				length(unique(y))
              }
)
unique(strsplit(paste(names(table(plot_dat_new$MUT)),collapse=";"),";" ) [[1]])
# [1] "CN_Amp"        "frameshift"    "fusion"        "missense"      "inframe_indel" "splice"        "stop_gained"   "CN_Del"        "germline"      "stop_lost"     "start_lost"

plot_dat_status=merge(plot_dat_cnv,plot_dat_mut_status, by="tag",all=T)
plot_dat_status=cbind( tag=plot_dat_status$tag, Gene=gsub(" .*","",plot_dat_status$tag), Sample=gsub(".* ","",plot_dat_status$tag) ,  CN=plot_dat_status$paired_mutation_status.x, Mut=as.character(plot_dat_status$paired_mutation_status.y) )
plot_dat_status=merge(plot_dat_status,plot_dat_gBRCA, by="tag",all=T)
plot_dat_status=cbind( tag=plot_dat_status$tag, Gene=gsub(" .*","",plot_dat_status$tag), Sample=gsub(".* ","",plot_dat_status$tag) ,  CN=plot_dat_status$CN, Mut= plot_dat_status$Mut , Germline=plot_dat_status$germline )
plot_dat_status=merge(plot_dat_status,plot_dat_fusion, by="tag",all=T)
plot_dat_status=cbind( tag=plot_dat_status$tag, Gene=gsub(" .*","",plot_dat_status$tag), Sample=gsub(".* ","",plot_dat_status$tag) ,  CN=plot_dat_status$CN, Mut= plot_dat_status$Mut , Germline=plot_dat_status$Germline, Fusion=plot_dat_status$paired_mutation_status )
plot_dat_status=data.frame(plot_dat_status)
plot_dat_status$MUT=apply(plot_dat_status[,4:7], 1, 
              function(x) {
                x = x[x != ""]
				x = x[!is.na(x)]
                y=paste(unique(x), collapse = ";")
				y=unlist(strsplit(y,";") )
				paste(unique(y), collapse = ";")
              }
)
plot_dat_status$MUTN=apply(plot_dat_status[,4:7], 1, 
              function(x) {
                x = x[x != ""]
				x = x[!is.na(x)]
                y=paste(unique(x), collapse = ";")
				y=unlist(strsplit(y,";") )
				length(unique(y))
              }
)
unique(strsplit(paste(names(table(plot_dat_status$MUT)),collapse=";"),";" ) [[1]])
# [1] "CN_Amp"        "frameshift"    "fusion"        "missense"      "inframe_indel" "splice"        "stop_gained"   "CN_Del"        "germline"      "stop_lost"     "start_lost"
#992     TP53 dea061fa-1f2f-4974-8953-e29239209ee5_Pre   TP53  dea061fa-1f2f-4974-8953-e29239209ee5_Pre CN_Amp                      frameshift     <NA>   <NA>             CN_Amp;frameshift    2
#992     TP53 dea061fa-1f2f-4974-8953-e29239209ee5_Pre   TP53  dea061fa-1f2f-4974-8953-e29239209ee5_Pre       Lost                 Maintained     <NA>       <NA> Lost;Maintained    2



setdiff( plot_dat_new$Sample, plot_samples)
setdiff( plot_dat_status$Sample, plot_samples)
setdiff( Mutation$sample_id, plot_samples)

plot_dat_new=plot_dat_new[plot_dat_new$Sample%in%plot_samples,]
plot_dat_status=plot_dat_status[plot_dat_status$Sample%in%plot_samples,]

table(plot_dat_new$MUT)
table(plot_dat_status$MUT)
#table(plot_dat_new$Sample)

#write.csv(plot_dat_new,file=paste0(outprefix1,".csv") )

mat=acast(data.frame(plot_dat_new), Gene ~ Sample  , value.var="MUT") 
status_mat=acast(data.frame(plot_dat_status), Gene ~ Sample  , value.var="MUT") 

for (sample in setdiff( plot_samples , colnames(mat) ) )
{
	mat=cbind(mat,sample=NA)
	colnames(mat)[ncol(mat)]=sample
	status_mat=cbind(status_mat,sample=NA)
	colnames(status_mat)[ncol(status_mat)]=sample
}

mat=mat[,plot_samples]
mat=mat[rownames(mat)%in%plot_genes ,]
status_mat=status_mat[,plot_samples]
status_mat=status_mat[rownames(status_mat)%in%plot_genes ,]

get_type_fun = function(x) strsplit(x, ";")[[1]]

#mat=gsub("(stop_gained|stop_lost|start_lost)","stop_gained/stop_lost/start_lost",mat)
#plot_dat_new$MUT=gsub("(stop_gained|stop_lost|start_lost)","stop_gained/stop_lost/start_lost",plot_dat_new$MUT)
mat=gsub("(stop_gained|stop_lost)","stop_gained/stop_lost",mat)
plot_dat_new$MUT=gsub("(stop_gained|stop_lost)","stop_gained/stop_lost",plot_dat_new$MUT)
col = c( "frameshift"= "darkgray",
         "missense"= "darkgreen",
         "stop_gained/stop_lost"= "black", #"start_lost" = "black"  
         "splice"= "orange",
		 "CN_Amp"= "red",
		 "CN_Del" = "blue",
         "germline" = "green",
         "fusion" ="purple",
         "inframe_indel"="brown" ,
 		 "mixed" = "turquoise4",
		 "WT" = "white"
		 )

if (length(setdiff(unique(strsplit(paste(names(table(plot_dat_new$MUT)),collapse=";"),";" ) [[1]]), names(col) ) ) ){stop("error")}		 
setdiff(unique(strsplit(paste(names(table(plot_dat_new$MUT)),collapse=";"),";" ) [[1]]), names(col) )
setdiff( names(col),unique(strsplit(paste(names(table(plot_dat_new$MUT)),collapse=";"),";" ) [[1]]) )

col_status = c( "Acquired"= "hotpink",
         "Lost"= "deepskyblue3",
         "Maintained"= "darkblue", 
         "unpaired"= "gray30",
		 "germline" = "green"
		 )
if (length(setdiff(unique(strsplit(paste(names(table(plot_dat_status$MUT)),collapse=";"),";" ) [[1]]), names(col_status) ) ) ){stop("error")}		 
setdiff(unique(strsplit(paste(names(table(plot_dat_status$MUT)),collapse=";"),";" ) [[1]]), names(col_status) )
setdiff( names(col_status),unique(strsplit(paste(names(table(plot_dat_status$MUT)),collapse=";"),";" ) [[1]]) )
		 
alter_fun = function(x, y, w, h, v) {
	n = sum(v)
     grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "white", col = "lightgray" , lwd =0.1))
	 
	 h=h*0.9
	 w=w*0.9
	 if(n) grid.rect(x, y - h*0.5 + 1:n/n*h, w , 1/n*h, 
            gp = gpar(fill = col[names(which(v))], col = NA), just = "top")
}
  
tmp=table(TEMPUS$Annotation$Sample$patient_id)
sel.patient=names(tmp)[(tmp>=2)]
sel.patient=grep("2cda8962-55da-4383-9b36-3340d44f9461",sel.patient,invert=TRUE,value=TRUE)
  

PostPre_freq2= matrix(0, nc = 2,nr = nrow(mat)) 
rownames(PostPre_freq2)=rownames(mat)
colnames(PostPre_freq2)=c("Pre","Post" )
{ 
    pdsamples=grep("_Post",plot_samples,value=TRUE)
    nonpdsamples=grep("_Pre",plot_samples,value=TRUE )	
	PostPre_freq2[ ,'Pre']=rowSums(!is.na(mat[,nonpdsamples,drop=F ]))
	PostPre_freq2[ ,'Post']=rowSums(!is.na(mat[,pdsamples,drop=F ]))
}

gene_order=plot_genes
mat=mat[gene_order,]
status_mat=status_mat[gene_order,]
PostPre_freq2=PostPre_freq2[gene_order, ]

sample_order=colnames(mat)[do.call(order, data.frame(t(mat)))]

#text=paste(pdspecific_freq, paste0("(",round( pdspecific_freq/length(sel.patient)*100,digits=1), "%)") )

row_ha = rowAnnotation(  "N of alterations" = anno_barplot(PostPre_freq2, gp = gpar(fill = c("Pre" =  "lightgrey","Post"="gray20") )), 
					#," " = anno_text(text) 	,
					width = unit(3, "cm")
					)
lgd_rowha = Legend(labels = colnames(PostPre_freq2), title = "Time point", legend_gp = gpar(fill = c("Pre" =  "lightgrey","Post"="gray20")),nrow = 1  )
library(circlize)
#col_fun = colorRamp2(c(0,6.1,6.2, 30), c("red3", "red3","darkgreen","darkgreen"))
#colseq=col_fun(TEMPUS$Annotation$Sample[colnames(mat), "first_progression_post_cdk" ])
rownames(TEMPUS$Annotation$Sample)=TEMPUS$Annotation$Sample$ID
columnha <- HeatmapAnnotation(	which="column",
								empty = anno_block(  gp = gpar(fill = c("Pre" =  "lightgrey","Post"="gray20")) ) 
								#TimePoint =  TEMPUS$Annotation$Sample[colnames(mat), "PrePostCDKi" ]  ,
								#col = list(TimePoint = c("Pre"="lightseagreen","Post"="lightsalmon")),
								#annotation_legend_param= list(title = "Time Point",nrow=1,at = c("Pre","Post"   ))								
							)
#colnames(mat)=gsub(".*-","",colnames(mat))

ht=oncoPrint(mat,
                top_annotation=columnha,
                #row_split = gene_split ,
                right_annotation = row_ha,
                show_pct=FALSE,
             alter_fun=alter_fun,
             get_type=get_type_fun,
             col = col,
             show_column_names = FALSE,
			 remove_empty_columns = FALSE, 
			 remove_empty_rows = FALSE,
             row_order = gene_order,
             column_order=sample_order,
             #column_title=column_title,
             column_title="",
			 column_split=factor(gsub(".*_","",colnames(mat)),levels= unique(gsub(".*_","",colnames(mat))) ),
			 column_gap = unit(0.5, "mm"),
             heatmap_legend_param = list(title = "Alterations",nrow=1,at = c("missense","frameshift","inframe_indel","stop_gained/stop_lost","splice","germline","CN_Amp","CN_Del","fusion" ))
             #heatmap_legend_param = list(title = "Alteration status",nrow=1,at = c("Acquired","Lost","Maintained","unpaired","germline"    ))
		 #top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(),
			#									PD = PDanno,
			#									col = list(PD = c("Patients progressed" = "blue", "Patients has not Progressed" = "green"  ))
			#	 )
		)
 
dev.off()
pdf(  paste0(outprefix2,".pdf") ,10,5,useDingbats=F)

#draw(ht, heatmap_legend_list = list(
#  Legend(labels = c("MISSENSE"),  
#         legend_gp=gpar(col =c("black") ),
#         type = "points" ,
#         pch = c(16, 16, 16, 16, 26 )  
#  )
#))
draw(ht, merge_legend = TRUE, 
    annotation_legend_list = list(lgd_rowha), 
    annotation_legend_side="bottom",heatmap_legend_side = "bottom" )
dev.off() 




q(save="no")










Baseline_freq= matrix(0, nc = 11,nr = length(plot_genes)) 
rownames(Baseline_freq)=plot_genes
colnames(Baseline_freq)=c("missense","frameshift","inframe_indel","start_lost/stop_gained/stop_lost","splice","germline","CN_Amp","CN_Del","fusion","mixed","WT")
for(gene in plot_genes)
{
	for ( muttype in colnames(Baseline_freq) )
	{
		if (muttype == "inframe_indel")
		{
			Baseline_freq[gene,muttype]=sum(plot_dat_new$MUT%in%c("inframe_indel"  ) & plot_dat_new$Gene==gene &  grepl("_Pre",plot_dat_new$Sample ) )
		}else if (muttype == "start_lost/stop_gained/stop_lost")
		{
			Baseline_freq[gene,muttype]=sum(plot_dat_new$MUT%in%c("start_lost","stop_gained","stop_lost") & plot_dat_new$Gene==gene & grepl("_Pre",plot_dat_new$Sample ) )
		}else
		{
			#sum(plot_dat_new$MUT==muttype & plot_dat_new$Gene==gene & !grepl("_PD",plot_dat_new$Sample ) )
			Baseline_freq[gene,muttype]=sum(plot_dat_new$MUT==muttype & plot_dat_new$Gene==gene &  grepl("_Pre",plot_dat_new$Sample ) )
		}
	}
	Baseline_freq[gene,"mixed"]=sum(plot_dat_new$MUTN>1 & plot_dat_new$Gene==gene &  grepl("_Pre",plot_dat_new$Sample ) )
}
rowSums(Baseline_freq)
table(plot_dat_new[ grepl("_Pre",plot_dat_new$Sample ) ,"Gene"])
#BRO7F−008_6weeks was included as baseline sample, so total +1 here
#Baseline_freq[,"WT"]= length(grep("_Baseline",unique(SMC_mutect2$SAMPLE_ID),value=TRUE )) +1 -rowSums(Baseline_freq)
Baseline_freq[,"WT"]= length(grep("_Pre",plot_samples,value=TRUE  )) -rowSums(Baseline_freq)
rownames(Baseline_freq)=paste("Pre",rownames(Baseline_freq) )

PD_freq= matrix(0, nc = 11,nr = length(plot_genes)) 
rownames(PD_freq)=plot_genes
colnames(PD_freq)=c("missense","frameshift","inframe_indel","start_lost/stop_gained/stop_lost","splice","germline","CN_Amp","CN_Del","fusion","mixed","WT")
for(gene in plot_genes)
{
	for ( muttype in colnames(PD_freq) )
	{
		if (muttype == "inframe_indel")
		{
			PD_freq[gene,muttype]=sum(plot_dat_new$MUT%in%c("inframe_indel" ) & plot_dat_new$Gene==gene & grepl("_Post",plot_dat_new$Sample ) )
		}else if (muttype == "start_lost/stop_gained/stop_lost")
		{
			PD_freq[gene,muttype]=sum(plot_dat_new$MUT%in%c("start_lost","stop_gained","stop_lost") & plot_dat_new$Gene==gene & grepl("_Post",plot_dat_new$Sample ) )
		}
		else
		{
			#sum(plot_dat_new$MUT==muttype & plot_dat_new$Gene==gene & grepl("_PD",plot_dat_new$Sample ) )
			PD_freq[gene,muttype]=sum(plot_dat_new$MUT==muttype & plot_dat_new$Gene==gene & grepl("_Post",plot_dat_new$Sample ) )
		}
	}
	PD_freq[gene,"mixed"]=sum(plot_dat_new$MUTN>1 & plot_dat_new$Gene==gene & grepl("_Post",plot_dat_new$Sample ) )
	PD_freq[gene,"WT"]= 0
}
sort(rowSums(PD_freq))
sort(table(plot_dat_new[ grepl("_Post",plot_dat_new$Sample ) ,"Gene"]))
#PD_freq[,"WT"]= length(grep("_PD",unique(SMC_mutect2$SAMPLE_ID),value=TRUE ))-rowSums(PD_freq)
PD_freq[,"WT"]= length(grep("_Post",plot_samples,value=TRUE ))-rowSums(PD_freq)
rownames(PD_freq)=paste("Post",rownames(PD_freq) )

Allsample_freq=rbind(Baseline_freq,PD_freq)
 

Allsample_freq_table=melt(Allsample_freq)
colnames(Allsample_freq_table)[colnames(Allsample_freq_table)=="value"]="Count"
colnames(Allsample_freq_table)[colnames(Allsample_freq_table)=="Var2"]="Alterations"
Allsample_freq_table$sample=gsub(" .*","",Allsample_freq_table$Var1)
Allsample_freq_table$Gene=gsub(".* ","",Allsample_freq_table$Var1)

col = c( "frameshift"= "darkgray",
         "missense"= "darkgreen",
         "start_lost/stop_gained/stop_lost"= "black",
         "splice"= "orange",
		 "CN_Amp"= "red",
		 "CN_Del" = "blue",
         "germline" = "green",
         "fusion" ="purple",
         "inframe_indel"="brown",
		 "mixed" = "turquoise4",
		 "WT" = "white"
		 )
		 

#selected_gene_table=Allsample_freq_table[Allsample_freq_table$Gene%in%c('TP53', 'PIK3CA', 'GATA3', 'KMT2C', 'PTEN', 'ESR1', 'RB1', 'BRCA1','BRCA2', 'FGFR1', 'MYC'),]
##dev.off()
##pdf(  "4.3.PostvsPre.stackedbar.pdf" ,9.5,3,useDingbats=F)
#ggplot(selected_gene_table, aes(x = sample, y = Count , fill = Alterations , colour  = Alterations)) +ylab("Percentage of altered samples")+
# geom_col(position = "fill",colour = "black",size=0.05 ) + 
# scale_y_continuous(labels = scales::percent) +
# scale_fill_manual(values =  col[levels(selected_gene_table$Alterations)] )+
# facet_grid(  ~ Gene   )+theme_classic()+ theme(axis.title.x=element_blank(),legend.position = "bottom",legend.key.size = unit(0.4, "cm") )+  
#  guides(fill = guide_legend(nrow = 1,title.position = "top"))
##dev.off() 

Allsample_freq2=cbind(Baseline_freq,PD_freq)
rownames(Allsample_freq2)=gsub("Pre ","",rownames(Allsample_freq2)) 
colnames(Allsample_freq2)=c(paste(colnames(Baseline_freq),"in Pre") , paste(colnames(PD_freq),"in Post")  )
Allsample_freq2=data.frame(Allsample_freq2,check.names=FALSE)
Allsample_freq2[,"# Altered Pre"]=NA
Allsample_freq2[,"% Altered Pre"]=NA
Allsample_freq2[,"# Altered Post"]=NA
Allsample_freq2[,"% Altered Post"]=NA
Allsample_freq2[,"p-value"]=NA
Allsample_freq2_new=Allsample_freq2
percent <- function(x, digits = 2, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}
for(gene in rownames(Allsample_freq2))
{
	BLalt=grep(" in Pre",colnames(Allsample_freq2),value=T) 
	n1=Allsample_freq2[gene, BLalt[BLalt=="WT in Pre"] ]
	n2=sum(Allsample_freq2[gene, BLalt[BLalt!="WT in Pre"] ])
	PDalt=grep(" in Post",colnames(Allsample_freq2),value=T) 
	n3=Allsample_freq2[gene, PDalt[PDalt=="WT in Post"] ]
	n4=sum(Allsample_freq2[gene, PDalt[PDalt!="WT in Post"] ])
	muttable <- matrix(c(n1, n2, n3, n4), nrow = 2,
	              dimnames =    list(c("WT", "MT"),   c("Pre", "Post")))
	Allsample_freq2_new[gene,"# Altered Pre"]=n2
	Allsample_freq2_new[gene,"% Altered Pre"]=percent(n2/(n2+n1))
	Allsample_freq2_new[gene,"# Altered Post"]=n4
	Allsample_freq2_new[gene,"% Altered Post"]=percent(n4/(n4+n3))
	Allsample_freq2_new[gene,"p-value"]=fisher.test(muttable,alternative="greater"  )$p.value
	
	#Allsample_freq2_new[gene, BLalt ]=paste0(Allsample_freq2[gene, BLalt ] , " (", signif(Allsample_freq2[gene, BLalt ]/sum(Allsample_freq2[gene, BLalt ]),3), ")")
	#Allsample_freq2_new[gene, PDalt ]=paste0(Allsample_freq2[gene, PDalt ] , " (", signif(Allsample_freq2[gene, PDalt ]/sum(Allsample_freq2[gene, PDalt ]),3), ")")
}

Allsample_freq2_new=Allsample_freq2_new[order(Allsample_freq2_new$"p-value"),]
Allsample_freq2_new$FDR = p.adjust(Allsample_freq2_new$"p-value" , method = 'BH')

#write.csv(Allsample_freq2_new,file=paste0(outprefix2,".csv"))
 
pdf(  paste0(outprefix2,".pdf")   ,3.5,3.5,useDingbats=F)
for (i in 1:ceiling(length(plot_genes)/9) )
{ 
	ind=((i-1)*9+1):(i*9)
	GeneList=rownames(Allsample_freq2_new)[ind]
	selected_gene_table=Allsample_freq_table[Allsample_freq_table$Gene%in% GeneList,]
	selected_gene_table$Gene=factor(selected_gene_table$Gene,levels=GeneList)
	selected_gene_table$sample=factor(selected_gene_table$sample,levels=c("Pre","Post"))
	selected_gene_table$Alterations=factor(selected_gene_table$Alterations,levels=rev(levels(selected_gene_table$Alterations)))

	p1=ggplot(selected_gene_table, aes(x = sample, y = Count , fill = Alterations , colour  = Alterations)) +ylab("Percentage of altered samples")+
	 geom_col(position = "fill",colour = "black",size=0.05 ) + 
	 scale_y_continuous(labels = scales::percent) +
	 scale_fill_manual(values =  col[levels(selected_gene_table$Alterations)] )+
	 facet_grid(  ~ Gene   )+theme_classic()+ 
	 theme(axis.title.x=element_blank(),axis.text.x = element_text(angle =  45, hjust = 1),legend.position = "bottom",legend.key.size = unit(0.4, "cm") )+  
	 theme(legend.position = "none") +
	  guides(fill = guide_legend(nrow = 1,title.position = "top"))
	print(p1)
}
dev.off() 




  