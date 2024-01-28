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
#outprefix1=paste0("Figure4a.plot_dat_new.",CNgaincut,CNlosscut)
outprefix2=paste0("Figure1i.stackedbar.",CNgaincut,CNlosscut,'.v',format(Sys.time(), "%Y%m%d"))
  
#################################################################################################################################################################################################################################################################################################
#################copy number 

length(table(Mutation$sample_id ) )
#421
plot_samples= sort( unique(TEMPUS$Annotation$Sample$ID ) )  
table(gsub(".*_","",plot_samples))
# Post  Pre
# 225  197

filter1=read.table("Input/tempus_ms_genomic_genelist.txt")
filter1=filter1[,1]
filter2=read.csv("Input/1.2.PostvsPre.TempusCN.stackedbar.70.csv",check.names=FALSE,row.names=1)
filter2=rownames(filter2)[filter2$'p-value'<0.05 & filter2$'# Altered Post'>0.05*227 ]
plot_genes=intersect(filter2,filter1)
plot_genes=setdiff(plot_genes, "YEATS4")
#plot_genes= unique(c("ESR1", "RB1", "AKT1", "PIK3CA","CDKN1A", "PTEN", "TP53", "EGFR", "ERBB2", "KRAS", "CDKN2A","FGFR1", "FGFR2", "AURKA", "MAP2K4", "CDK4",   "CCND1",  "CDKN2A", "MYC", "ERBB2", "PTEN") )
#plot_genes= unique(c("ESR1", "RB1", "TP53", "FGFR2", "PTEN") )
length(plot_genes)

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
#plot_dat_mut=plot_dat_mut[!duplicated(plot_dat_mut$tag1),]
plot_dat_mut= aggregate(x = plot_dat_mut["CONSEQUENCE"],   
                       by = plot_dat_mut[c("SAMPLE_ID", "GENE_NAME")],
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

names(TEMPUS$SV  )
 tmp1=TEMPUS$SV$SV_from_DNA_filtered[, c("sample_id","gene_canonical_name"  )]
 tmp2=TEMPUS$SV$SV_from_DNA_filtered[, c("sample_id","secondary_gene_canonical_name"  )]
colnames(tmp1)=c("sample","gene"  )
colnames(tmp2)=c("sample","gene"  )
 tmp3=TEMPUS$SV$SV_from_RNASeq_filtered[, c("sample_id","gene_canonical_name"  )]
 tmp4=TEMPUS$SV$SV_from_RNASeq_filtered[, c("sample_id","secondary_gene_canonical_name"  )]
colnames(tmp3)=c("sample","gene"  )
colnames(tmp4)=c("sample","gene"  )

plot_dat_fusion=rbind(tmp1, tmp2, tmp3, tmp4)

dim(plot_dat_fusion)
plot_dat_fusion=plot_dat_fusion[ plot_dat_fusion$gene!='',]
plot_dat_fusion$gene=gsub(":.*","",plot_dat_fusion$gene)

plot_dat_fusion[plot_dat_fusion$gene=="ESR1",]
plot_dat_fusion=unique(plot_dat_fusion)
dim(plot_dat_fusion)
plot_dat_fusion=plot_dat_fusion[plot_dat_fusion$gene%in%plot_genes,]
dim(plot_dat_fusion)
plot_dat_fusion=cbind(plot_dat_fusion,fusion="fusion")
plot_dat_fusion$tag=paste(plot_dat_fusion$gene,plot_dat_fusion$sample)

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

setdiff( plot_dat_new$Sample, plot_samples)
setdiff( Mutation$sample_id, plot_samples)

plot_dat_new=plot_dat_new[plot_dat_new$Sample%in%plot_samples,]
 
table(plot_dat_new$MUT)
table(plot_dat_new$Sample)

#write.csv(plot_dat_new,file=paste0(outprefix1,".csv") )

grep(";",names(table(plot_dat_new$MUT)) ,invert=T ,value=T)
unique(unlist(strsplit(paste( unique(plot_dat_new$MUT),collapse=';'),';' )) )
#[1] "CN_Amp"        "CN_Del"        "frameshift"    "fusion"        "germline"      "inframe_indel" "missense"      "splice"        "start_lost"    "stop_gained"   "stop_lost"

Baseline_freq= matrix(0, nc = 11,nr = length(plot_genes)) 
rownames(Baseline_freq)=plot_genes
colnames(Baseline_freq)=c("missense","frameshift","inframe_indel","stop_gained","splice","germline","CN_Amp","CN_Del","fusion","mixed","WT")
for(gene in plot_genes)
{
	for ( muttype in colnames(Baseline_freq) )
	{
		if (muttype == "inframe_indel")
		{
			Baseline_freq[gene,muttype]=sum(plot_dat_new$MUT%in%c("inframe_indel"  ) & plot_dat_new$Gene==gene &  grepl("_Pre",plot_dat_new$Sample ) )
		}else if (muttype == "stop_gained")
		{
			#Baseline_freq[gene,muttype]=sum(plot_dat_new$MUT%in%c("start_lost","stop_gained","stop_lost") & plot_dat_new$Gene==gene & grepl("_Pre",plot_dat_new$Sample ) )
			Baseline_freq[gene,muttype]=sum(plot_dat_new$MUT%in%c( "stop_gained" ) & plot_dat_new$Gene==gene & grepl("_Pre",plot_dat_new$Sample ) )
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
colnames(PD_freq)=c("missense","frameshift","inframe_indel","stop_gained","splice","germline","CN_Amp","CN_Del","fusion","mixed","WT")
for(gene in plot_genes)
{
	for ( muttype in colnames(PD_freq) )
	{
		if (muttype == "inframe_indel")
		{
			PD_freq[gene,muttype]=sum(plot_dat_new$MUT%in%c("inframe_indel" ) & plot_dat_new$Gene==gene & grepl("_Post",plot_dat_new$Sample ) )
		}else if (muttype == "stop_gained")
		{
			#PD_freq[gene,muttype]=sum(plot_dat_new$MUT%in%c("start_lost","stop_gained","stop_lost") & plot_dat_new$Gene==gene & grepl("_Post",plot_dat_new$Sample ) )
			PD_freq[gene,muttype]=sum(plot_dat_new$MUT%in%c( "stop_gained" ) & plot_dat_new$Gene==gene & grepl("_Post",plot_dat_new$Sample ) )
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
         "stop_gained"= "black",
         "splice"= "orange",
		 "CN_Amp"= "red",
		 "CN_Del" = "blue",
         "germline" = "green",
         "fusion" ="purple",
         "inframe_indel"="brown",
		 "mixed" = "turquoise4",
		 "WT" = "white"
		 )
 
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
 
pdf(  paste0(outprefix2,".pdf")   ,7, 4.5,useDingbats=F)
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
	 geom_col(position = "fill",colour = "NA",size=0.05,show_guide=FALSE ) + 
	 scale_y_continuous(labels = scales::percent) + coord_cartesian(ylim=c(0,0.5))+
	 scale_fill_manual(values =  col[levels(selected_gene_table$Alterations)] ,limits = setdiff(levels(selected_gene_table$Alterations),'WT') ,na.value="white")+
	 facet_grid(  ~ Gene   )+theme_classic()+ 
	 theme(axis.title.x=element_blank(),axis.text.x = element_text(angle =  45, hjust = 1),legend.position = "right",legend.key.size = unit(0.4, "cm") )+  
	  guides(fill = guide_legend(ncol = 1,byrow = TRUE,title.position = "top"))
	print(p1)
 
}
dev.off() 




q()