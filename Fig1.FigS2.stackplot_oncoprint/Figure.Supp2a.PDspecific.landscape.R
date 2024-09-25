rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE)
#options(width=200)

setwd('../Fig1.FigS2.stackplot_oncoprint')

options(stringsAsFactors = FALSE)
library(ComplexHeatmap)
library(reshape2)
library(tidyr)
library(ggplot2)
library(stringr)

molfeatfile="Input/Tempus_molecular_features_extended.annonymized.txt"
mfdata <- read.table(molfeatfile, header=T, sep='\t');
rownames(mfdata)=mfdata$sample_id_anonymized

plot_dat_new<- read.csv(file= "Input/Tempus.alternations.csv")  
plot_dat_new$tag=paste(plot_dat_new$Gene,plot_dat_new$Sample)

outpdf=paste0("Output.Manuscript/SuppFigure2A.PDspecific.landscape.v",format(Sys.time(), "%Y%m%d"),".pdf")

########################################################################################################################################################################################################################################################################################################################################################################################################################################################################
sel.patient= names(table(mfdata$patient_id ))[table(mfdata$patient_id )>=2]
sel.patient= grep("tempus075",sel.patient,invert=TRUE,value=TRUE) #tempus075 pre match, post tumor only
plot_samples= sort( unique(c( mfdata$sample_id_anonymized[ mfdata$patient_id%in%sel.patient ]  ) ) )
table(gsub(".*_","",plot_samples)) 

plot_genes=c("ESR1","KMT2D","BCL11B","CDKN1B","ELF3","FLT3","MYH11","PTPN13","RB1","TSC2"  ) #genes with >=3 PD specific mutations

plot_dat_new=plot_dat_new[plot_dat_new$Sample%in%plot_samples&plot_dat_new$Gene%in%plot_genes,]

table(plot_dat_new$MergedAlterations)

mat=acast(data.frame(plot_dat_new), Gene ~ Sample  , value.var="MergedAlterations") 
for (sample in setdiff( plot_samples , colnames(mat) ) )
{
  mat=cbind(mat,sample=NA)
  colnames(mat)[ncol(mat)]=sample
}
mat=mat[,plot_samples]
mat=mat[rownames(mat)%in%plot_genes ,]
#status_mat=status_mat[rownames(status_mat)%in%plot_genes ,]

get_type_fun = function(x) strsplit(x, ";")[[1]]
column_title = "Mutations specific in PD samples (potential acquired mutation)"

mat=gsub("(stop_gained|stop_lost|start_lost)","stop_gained/stop_lost/start_lost",mat)
plot_dat_new$MergedAlterations=gsub("(stop_gained|stop_lost|start_lost)","stop_gained/stop_lost/start_lost",plot_dat_new$MergedAlterations)
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

muttypes=unique(strsplit(paste(names(table(plot_dat_new$MergedAlterations)),collapse=";"),";" ) [[1]])
muttypes

if (length(setdiff(muttypes, names(col) ) ) ){stop("error")}		 
setdiff(muttypes, names(col) )
setdiff( names(col),muttypes )

if (length(setdiff(unique(strsplit(paste(names(table(plot_dat_new$MergedAlterations)),collapse=";"),";" ) [[1]]),
                   names(col)
) ) ){stop("error")}

alter_fun = function(x, y, w, h, v) {
  n = sum(v)
  grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#EFEFEF", col = "NA" ))
  
  h=h*0.9
  w=w*0.9
  if(n) grid.rect(x, y - h*0.5 + 1:n/n*h, w , 1/n*h, 
                  gp = gpar(fill = col[names(which(v))], col = 'black'), just = "top")
  
}

rownames(mat)
 
#calculate pdspecific_freq. for each gene, this is the number of patient with PD specific mutations of this gene in PD samples
#pdspecific_freq is a vector
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
      #tmp=  grepl(muttype, mat[,pdsamples,drop=F ] ) & !grepl(muttype, mat[,nonpdsamples,drop=F ] )  
      tmp=str_count(tmpmat[,pdsamples,drop=F ] ,muttype) - str_count(tmpmat[,nonpdsamples,drop=F ] ,muttype) 
      tmp[tmp<0]=0
      tmp[tmp>0]=1 #two pd-specific missense for the same patient count once
      pdspecific_genextype[,muttype]=pdspecific_genextype[,muttype]+tmp
    }
    tmp = rowSums(pdspecific_genextype)>0
    pdspecific_freq=pdspecific_freq+tmp
  }
}

#calculate pdspecific_freq2. it is matrix with more detailed info of number of PD specific alternations for each mutation type
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
      tmp[tmp>0]=1 #two pd-specific missense count once
      pdspecific_genextype[,muttype]=pdspecific_genextype[,muttype]+tmp
    }
    ind = rowSums(pdspecific_genextype)>1 #genes with mixed mutaiton comparingt this PD and this Pre of this patient
    pdspecific_genextype[ind,]=0
    pdspecific_genextype=cbind(pdspecific_genextype,mixed=ind)
    pdspecific_freq2 =pdspecific_freq2 + pdspecific_genextype
  }
  else
  {
    print(paste("error" ,p))
  }
}

#gene_order= c("ESR1","RB1", "ERBB2","PTEN" , "CCND1", "ASXL3", "TP53", "PIK3CA", "BRCA1","BRCA2", "PPM1D", "RIF1",  "USP2",  "KDM6A","FGFR1", "FBN2","ATAD2","C5orf42","PCNXL2","PRDM1")
gene_order=names(pdspecific_freq)[order(pdspecific_freq,names(pdspecific_freq)=='RB1',decreasing=TRUE) ]
plotmat=mat[gene_order,]
pdspecific_freq=pdspecific_freq[gene_order ]
pdspecific_freq2=pdspecific_freq2[gene_order, ]
text=paste(pdspecific_freq, paste0("(",round( pdspecific_freq/length(sel.patient)*100,digits=1), "%)") )

row_ha = rowAnnotation(  "N (%) of PD-\nspecific alterations" = anno_barplot(pdspecific_freq2, width = unit(4, "cm"), gp = gpar(fill = col[colnames(pdspecific_freq2)] )) ," " = anno_text(text))

library(circlize)
col_fun = colorRamp2(c(0,6.1,6.2, 30), c("red3", "red3","darkgreen","darkgreen"))
columnha <- HeatmapAnnotation(	which="column",
                               "Pre/Post" =  mfdata[colnames(plotmat), "pre_post" ]  ,
                               col = list("Pre/Post" = c("Pre"="lightgrey","Post"="gray20")),
                               annotation_legend_param= list(title = "Pre/Post",nrow=1,at = c("Pre","Post"   ))								
)

 
ht=oncoPrint(plotmat,
             top_annotation=columnha,
             #top_annotation=NULL,
             #row_split = gene_split ,
             right_annotation = row_ha,
             show_pct=FALSE,
             alter_fun=alter_fun,
             get_type=get_type_fun,
             col = col,
             show_column_names = TRUE,
             remove_empty_columns = FALSE, 
             remove_empty_rows = FALSE,
             row_order = gene_order,
             column_order=ncol(plotmat):1,
             #column_title=column_title,
             column_title=" ",
             column_split=factor(gsub("_.*","",colnames(plotmat)),levels= unique(gsub("_.*","",colnames(plotmat))) ),
             column_gap = unit(3, "mm"),
             column_labels= colnames(plotmat)  ,
             heatmap_legend_param = list(title = "Alterations",nrow=1,at = c("missense","frameshift","inframe_indel","stop_gained/stop_lost/start_lost","splice","germline","CN_Amp","CN_Del","fusion"))
)
 
pdf(  outpdf ,16,5,useDingbats=F)
draw(ht, merge_legend = TRUE, 
     #annotation_legend_list = list(lgd_PFS), 
     annotation_legend_side="bottom",heatmap_legend_side = "bottom" )
dev.off() 



#q(save='no')
