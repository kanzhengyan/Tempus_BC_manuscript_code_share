rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE)
#options(width=200)
library(ComplexHeatmap)
library(reshape2)
library(tidyr)
library(ggplot2)

setwd('../Fig1.FigS2.stackplot_oncoprint')

molfeatfile="Input/Tempus_molecular_features_extended.annonymized.txt"
mfdata <- read.table(molfeatfile, header=T, sep='\t');
rownames(mfdata)=mfdata$sample_id_anonymized

plot_dat_new<- read.csv(file= "Input/Tempus.alternations.csv")  
plot_dat_new$tag=paste(plot_dat_new$Gene,plot_dat_new$Sample)

outpdf=paste0("Output.Manuscript/Figure1e.stackedbar.v",format(Sys.time(), "%Y%m%d"),".pdf")

####################### 
  
plot_samples= sort( unique(mfdata$sample_id_anonymized  ) )  
table(gsub(".*_","",plot_samples))
 
plot_genes=c( "ESR1","RB1","LRP1B","FGFR2","MDM2","KMT2A","TP53" )
      
setdiff( plot_dat_new$Sample, plot_samples)

plot_dat_new=plot_dat_new[plot_dat_new$Sample%in%plot_samples&plot_dat_new$Gene%in%plot_genes,]

plot_dat_new$MergedAlterations
plot_dat_new$MergedAlterations=sapply(plot_dat_new$MergedAlterations,  
                                     function(x) {
                                       y=unlist(strsplit(x,";") )
                                       paste(unique(y), collapse = ";")
                                       #some sample can have 2 missense mutations. we are intereste in sample level mutation frequency difference between Pre and Post, so if a sample contains two missense mutations for the same gene, count only once
                                     }
)

table(plot_dat_new$MergedAlterations)
table(plot_dat_new$Sample)

grep(";",names(table(plot_dat_new$MergedAlterations)) ,invert=T ,value=T)
unique(unlist(strsplit(paste( unique(plot_dat_new$MergedAlterations),collapse=';'),';' )) )
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
			Baseline_freq[gene,muttype]=sum(plot_dat_new$MergedAlterations%in%c("inframe_indel"  ) & plot_dat_new$Gene==gene &  grepl("_Pre",plot_dat_new$Sample ) )
		}else if (muttype == "stop_gained")
		{
			Baseline_freq[gene,muttype]=sum(plot_dat_new$MergedAlterations%in%c( "stop_gained" ) & plot_dat_new$Gene==gene & grepl("_Pre",plot_dat_new$Sample ) )
		}else
		{
			Baseline_freq[gene,muttype]=sum(plot_dat_new$MergedAlterations==muttype & plot_dat_new$Gene==gene &  grepl("_Pre",plot_dat_new$Sample ) )
		}
	}
	Baseline_freq[gene,"mixed"]=sum(plot_dat_new$UniqueAlterationN>1 & plot_dat_new$Gene==gene &  grepl("_Pre",plot_dat_new$Sample ) )
}
rowSums(Baseline_freq)
table(plot_dat_new[ grepl("_Pre",plot_dat_new$Sample ) ,"Gene"])
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
			PD_freq[gene,muttype]=sum(plot_dat_new$MergedAlterations%in%c("inframe_indel" ) & plot_dat_new$Gene==gene & grepl("_Post",plot_dat_new$Sample ) )
		}else if (muttype == "stop_gained")
		{
			PD_freq[gene,muttype]=sum(plot_dat_new$MergedAlterations%in%c( "stop_gained" ) & plot_dat_new$Gene==gene & grepl("_Post",plot_dat_new$Sample ) )
		}
		else
		{
			PD_freq[gene,muttype]=sum(plot_dat_new$MergedAlterations==muttype & plot_dat_new$Gene==gene & grepl("_Post",plot_dat_new$Sample ) )
		}
	}
	PD_freq[gene,"mixed"]=sum(plot_dat_new$UniqueAlterationN>1 & plot_dat_new$Gene==gene & grepl("_Post",plot_dat_new$Sample ) )
	PD_freq[gene,"WT"]= 0
}
sort(rowSums(PD_freq))
sort(table(plot_dat_new[ grepl("_Post",plot_dat_new$Sample ) ,"Gene"]))
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
}

Allsample_freq2_new=Allsample_freq2_new[order(Allsample_freq2_new$"p-value"),]
Allsample_freq2_new$FDR = p.adjust(Allsample_freq2_new$"p-value" , method = 'BH')
 
pdf(  outpdf  ,7, 4.5,useDingbats=F)
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




#q(save="no")