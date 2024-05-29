rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE)
options(width=200)
#setwd("Figures/ForQC.20231018")

options(stringsAsFactors = FALSE)
library(ComplexHeatmap)
library(reshape2)
library(tidyr)
library(ggplot2)
library(stringr)

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
outprefix2=paste0("Figure.Supp4b.PDspecific.landscape.v",format(Sys.time(), "%Y%m%d"))

#################################################################################################################################################################################################################################################################################################

filter1=read.table("Input/tempus_ms_genomic_genelist.txt")
filter1=filter1[,1]
filter2=read.csv("Input/3.master_table_gene.PDspecific.Frequency.csv",check.names=FALSE,row.names=1)
filter2=rownames(filter2)[filter2$'Sum'>=3  ]
plot_genes=intersect(filter2,filter1) 
#plot_genes=setdiff(plot_genes, "YEATS4")
#plot_genes= unique(c("ESR1", "RB1", "AKT1", "PIK3CA","CDKN1A", "PTEN", "TP53", "EGFR", "ERBB2", "KRAS", "CDKN2A","FGFR1", "FGFR2", "AURKA", "MAP2K4", "CDK4",   "CCND1",  "CDKN2A", "MYC", "ERBB2", "PTEN") )
#length(plot_genes)
 
length(table(Mutation$sample_id ) )
#426
#plot_samples= sort( unique(TEMPUS$Annotation$Sample$ID ) )  
sel.patient= names(table(TEMPUS$Annotation$Sample$patient_id ))[table(TEMPUS$Annotation$Sample$patient_id )>=2]
plot_samples= sort( unique(c( TEMPUS$Annotation$Sample$ID[ TEMPUS$Annotation$Sample$patient_id%in%sel.patient ]  ) ) )
plot_samples= grep("2cda8962-55da-4383-9b36-3340d44f9461",plot_samples,invert=TRUE,value=TRUE)#2cda8962-55da-4383-9b36-3340d44f9461 pre match, post tumor only
table(gsub(".*_","",plot_samples)) 
# Post  Pre
# 26  26
 
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
SV_master_raw=read.csv("Input/master_table_SV.csv")
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
				#paste(unique(y), collapse = ";")
				paste(y, collapse = ";")
              }
)
plot_dat_new$MUTN=apply(plot_dat_new[,4:7], 1, 
              function(x) {
                x = x[x != ""]
				x = x[!is.na(x)]
                y=paste(unique(x), collapse = ";")
				y=unlist(strsplit(y,";") )
				#length((y))
				length(unique(y))
              }
)

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

#mat=gsub("(stop_gained|stop_lost|start_lost)","stop_gained/stop_lost/start_lost/start_lost",mat)
#plot_dat_new$MUT=gsub("(stop_gained|stop_lost|start_lost)","stop_gained/stop_lost/start_lost/start_lost",plot_dat_new$MUT)
mat=gsub("(stop_gained|stop_lost|start_lost)","stop_gained/stop_lost/start_lost",mat)
plot_dat_new$MUT=gsub("(stop_gained|stop_lost|start_lost)","stop_gained/stop_lost/start_lost",plot_dat_new$MUT)
col = c( "frameshift"= "darkgray",
         "missense"= "darkgreen",
         "stop_gained/stop_lost/start_lost"= "black", #"start_lost" = "black"  
         "splice"= "orange",
		 "CN_Amp"= "red",
		 "CN_Del" = "blue",
         "germline" = "green",
         "fusion" ="purple",
         "inframe_indel"="brown" ,
 		 "mixed" = "turquoise4",
		 "WT" = "white"
		 )
		 
muttypes=unique(strsplit(paste(names(table(plot_dat_new$MUT)),collapse=";"),";" ) [[1]])
muttypes

if (length(setdiff(muttypes, names(col) ) ) ){stop("error")}		 
setdiff(muttypes, names(col) )
setdiff( names(col),muttypes )

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
            gp = gpar(fill = col_status[names(which(v))], col = NA), just = "top")
}
  
tmp=table(TEMPUS$Annotation$Sample$patient_id)
sel.patient=names(tmp)[(tmp>=2)]
sel.patient=grep("2cda8962-55da-4383-9b36-3340d44f9461",sel.patient,invert=TRUE,value=TRUE)
 
tmpmat=mat 
tmpmat[is.na(tmpmat)]=''
pdspecific_freq=rep(0,nrow(mat))
names(pdspecific_freq)=rownames(mat)
for(p in sel.patient)
{
    samples=grep(p,plot_samples,value=TRUE)
    pdsamples=grep("_Post",samples,value=TRUE)
    nonpdsamples=grep("_Pre",samples,value=TRUE )
    if(length(nonpdsamples))
    {
		pdspecific_genextype= matrix(0, nc = 9,nr = nrow(mat)) 
		colnames(pdspecific_genextype)=c("missense","frameshift","inframe_indel","stop_gained/stop_lost/start_lost","splice","germline","CN_Amp","CN_Del","fusion")
		for ( muttype in colnames(pdspecific_genextype) )
		{
			tmp=str_count(tmpmat[,pdsamples,drop=F ] ,muttype) - str_count(tmpmat[,nonpdsamples,drop=F ] ,muttype) 
			tmp[tmp<0]=0
			pdspecific_genextype[,muttype]=pdspecific_genextype[,muttype]+tmp
		}
        tmp = rowSums(pdspecific_genextype)>0
        pdspecific_freq=pdspecific_freq+tmp
    }
}

pdspecific_freq2= matrix(0, nc = 10,nr = nrow(mat)) 
rownames(pdspecific_freq2)=rownames(mat)
colnames(pdspecific_freq2)=c("missense","frameshift","inframe_indel","stop_gained/stop_lost/start_lost","splice","germline","CN_Amp","CN_Del","fusion","mixed" )
for(p in sel.patient)
{
    samples=grep(p,plot_samples,value=TRUE)
    pdsamples=grep("_Post",samples,value=TRUE)
    nonpdsamples=grep("_Pre",samples,value=TRUE )
    if(length(nonpdsamples))
    {
		pdspecific_genextype= matrix(0, nc = 9,nr = nrow(mat)) 
		colnames(pdspecific_genextype)=c("missense","frameshift","inframe_indel","stop_gained/stop_lost/start_lost","splice","germline","CN_Amp","CN_Del","fusion")
		for ( muttype in colnames(pdspecific_genextype) )
		{
			tmp=str_count(tmpmat[,pdsamples,drop=F ] ,muttype) - str_count(tmpmat[,nonpdsamples,drop=F ] ,muttype) 
			tmp[tmp<0]=0
			#tmp[tmp>0]=1
			pdspecific_genextype[,muttype]=pdspecific_genextype[,muttype]+tmp
		}
		ind = rowSums(pdspecific_genextype)>1 #genes with mixed mutaiton comparingt this PD and this Pre of this patient
		pdspecific_genextype[ind,]=0
		pdspecific_genextype=cbind(pdspecific_genextype,mixed=ind)
		pdspecific_freq2 =pdspecific_freq2 + pdspecific_genextype
    }
	else
	{
		print("error")
	}
}
gene_order=names(pdspecific_freq)[order(pdspecific_freq,names(pdspecific_freq)=='RB1',decreasing=TRUE) ]
mat=mat[gene_order,]
status_mat=status_mat[gene_order,]

sample_order=colnames(status_mat)[do.call(order, data.frame(t(status_mat)))]

pdspecific_freq=pdspecific_freq[gene_order ]
pdspecific_freq2=pdspecific_freq2[gene_order, ]
pdspecific_freq2=pdspecific_freq2[,colSums(pdspecific_freq2)!=0]
text=paste(pdspecific_freq, paste0("(",round( pdspecific_freq/length(sel.patient)*100,digits=1), "%)") )

row_ha = rowAnnotation(  "N(%) of PD-\nspecific alterations" = anno_barplot(pdspecific_freq2, gp = gpar(fill = col[colnames(pdspecific_freq2)] )) ,
						 " " = anno_text(text),
						 width = unit(5, "cm")
						 )
lgd_rowha = Legend(labels = colnames(pdspecific_freq2), title = "PD-specific alteration type", legend_gp = gpar(fill = col[colnames(pdspecific_freq2)]),nrow = 1  )
library(circlize)
#col_fun = colorRamp2(c(0,6.1,6.2, 30), c("red3", "red3","darkgreen","darkgreen"))
#colseq=col_fun(TEMPUS$Annotation$Sample[colnames(mat), "first_progression_post_cdk" ])
rownames(TEMPUS$Annotation$Sample)=TEMPUS$Annotation$Sample$ID
columnha <- HeatmapAnnotation(	which="column",
								empty = anno_block(  gp = gpar(fill = c("Pre" =  "lightgrey","Post"="gray20")) ) 
								#TimePoint =  TEMPUS$Annotation$Sample[colnames(status_mat), "PrePostCDKi" ]  ,
								#col = list(TimePoint = c("Pre"="lightseagreen","Post"="lightsalmon")),
								#annotation_legend_param= list(title = "Time Point",nrow=1,at = c("Pre","Post"   ))								
							)
#colnames(mat)=gsub(".*-","",colnames(mat))

ht=oncoPrint(status_mat,
                top_annotation=columnha,
                #row_split = gene_split ,
                right_annotation = row_ha,
                show_pct=FALSE,
             alter_fun=alter_fun,
             get_type=get_type_fun,
             col = col_status,
             show_column_names = FALSE,
			 remove_empty_columns = FALSE, 
			 remove_empty_rows = FALSE,
             row_order = gene_order,
             column_order=sample_order,
             #column_title=column_title,
             column_title="",
			 column_split=factor(gsub(".*_","",colnames(mat)),levels= c("Pre","Post") ),
			 column_gap = unit(0.5, "mm"),
            # heatmap_legend_param = list(title = "Alterations",nrow=1,at = c("missense","frameshift","inframe_indel","stop_gained","splice","germline","CN_Amp","CN_Del","fusion","start_lost"))
             heatmap_legend_param = list(title = "Alteration status",nrow=1,at = c("Acquired","Lost","Maintained","unpaired","germline"    )
				)
		 #top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(),
			#									PD = PDanno,
			#									col = list(PD = c("Patients progressed" = "blue", "Patients has not Progressed" = "green"  ))
			#	 )
		)
 
#dev.off()
pdf(  paste0(outprefix2,".pdf") ,15,5,useDingbats=F)

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



 