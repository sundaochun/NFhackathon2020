library(synapser)
synLogin(email="@@@", password="XXX")
set.seed('99999')  #set seed for reproducibility


pnf_cell_line <- synTableQuery("SELECT * FROM syn22351884")$asDataFrame()  
table(paste0(pnf_cell_line$specimenID, "_", pnf_cell_line$modelOf))

PNF.tab<-pnf_cell_line %>% dplyr::select(totalCounts,Symbol,zScore,specimenID,sex,tumorType,nf1Genotype,modelOf) %>% unique() %>% dplyr::arrange(specimenID)
PNF.pdata<-PNF.tab %>% dplyr::select(specimenID,sex,tumorType,nf1Genotype,modelOf) %>% unique()
rownames(PNF.pdata)<-PNF.pdata[,1]
table(pnf_cell_line$specimenID)#19098 gene per sample

diffexdata <- PNF.tab%>%dplyr::select(specimenID,Symbol,totalCounts)%>%
  group_by(specimenID,Symbol,totalCounts)%>%
  unique()%>%spread(specimenID,totalCounts)

PNF.mx<-data.frame(diffexdata)
rownames(PNF.mx)<-PNF.mx[,1]
PNF.mx<-PNF.mx[,-1]

library(Biobase)
colnames(PNF.mx)[4:5]<-c("ipNF05.5 (mixed clone)","ipNF05.5 (single clone)")
PNF.pdata<-PNF.tab %>% dplyr::select(specimenID,sex,tumorType,nf1Genotype,modelOf) %>% unique()
rownames(PNF.pdata)<-PNF.pdata[,1]
colnames(PNF.mx)==rownames(PNF.pdata)
PNF.pdata<-new("AnnotatedDataFrame",data=PNF.pdata)
PNF.eSet<-ExpressionSet(assayData=as.matrix(PNF.mx),phenoData=PNF.pdata,annotation="pnf_RNAseq")	
library(limma)
plotMDS(PNF.eSet, col=ifelse(PNF.eSet$modelOf=="normal","red","blue"))


datExpr = as.data.frame(t(exprs(PNF.eSet[,-c(1,2,9,10)]))) #remove the Nf1+/+ cell lines "ipn02.3","ipn97.4","pn02.3","pn97.4"  #Zsocre doesn't work, totalCounts are used
datExpr = as.data.frame(t(exprs(PNF.eSet[,sampleNames(PNF.eSet) %in% rownames(datTraits)]))) #only 5 tumor lines have drug response data 

mads=apply(datExpr,2,mad)
hist(mads)

datExpr.test<-datExpr[,mads>1000]
write.table(t(datExpr.test),"NTAP.5pnf.filtered.txt", quote = FALSE, sep = "\t")

#"GSE41747_expVal.tsv" is downloaded from synapse with ID:syn6130081  
#"GSE41747_phenotype_data.tsv" is downloaded from synapse with ID: syn6130082
GEO<-read.table("GSE41747_expVal.tsv", header=TRUE,row.names=1,stringsAsFactors = FALSE,sep="\t")
GEOanno<-read.table("GSE41747_phenotype_data.tsv",header = TRUE,stringsAsFactors = FALSE,sep = "\t")

GEO.HU<-GEO[,GEOanno$tissue=="pNF"]
GEO.MU<-GEO[,GEOanno$tissue=="Mouse neurofibroma"]

write.table(GEO.HU,"GEO41747.HU.pNF.txt",quote=FALSE,sep = "\t")

