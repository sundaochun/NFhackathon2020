library(dendextend)
library(tidyverse) 
library(synapser)
library(umap) 
library(dbscan)
library(dplyr)
library(pheatmap)
library(synapser)
library(WGCNA)
library(ggplot2)
library(ggdendro)
library(gplots)
library(biomaRt)
library(org.Hs.eg.db)
library(GOSemSim)


synLogin(email="@gmail.com", password="XXX")
set.seed('99999')  #set seed for reproducibility


#######################
#Start analysis 
setwd("~/NFhackathom2020/input")

NTAP.5pnf<-read.table("NTAP.5pnf.filtered.txt",header = TRUE,stringsAsFactors = FALSE,row.names = 1,sep="\t")
GEO41747.pNF<-read.table("GEO41747.HU.pNF.txt",header = TRUE,stringsAsFactors = FALSE, row.names = 1,sep="\t")

pNF.merged<-merge(NTAP.5pnf, GEO41747.pNF, by=0, all=FALSE)
rownames(pNF.merged)<-pNF.merged$Row.names

NTAPm<-pNF.merged[,2:6]
GEO41747m<-pNF.merged[,7:19]



setwd("~/NFhackathom2020/output")
nSets = 2;
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("NTAP", "GEO41747")

multiExpr = vector(mode = "list", length = nSets)
multiExpr[[1]] = list(data = as.data.frame(t(NTAPm)));
multiExpr[[2]] = list(data = as.data.frame(t(GEO41747m)));

names(multiExpr)<-c("NTAP","GEO")
exprSize = checkSets(multiExpr)


##############
# To generate the preserved networks

gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
gsg$allOK

if (!gsg$allOK)
{
  # Print information about the removed genes:
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes], 
                                              collapse = ", ")))
  for (set in 1:exprSize$nSets)
  {
    if (sum(!gsg$goodSamples[[set]]))
      printFlush(paste("In set", setLabels[set], "removing samples",
                       paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
    # Remove the offending genes and samples
    multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
    # Update exprSize
  }
  exprSize = checkSets(multiExpr)
}



sampleTrees = list()
for (set in 1:nSets)
{
  sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}


pdf(file = "ConsensusSampleClustering.pdf", width = 12, height = 12);
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets)
  plot(sampleTrees[[set]], main = paste("Sample Consensus clustering on all genes in", setLabels[set]),
       xlab="", sub="", cex = 0.7);
dev.off();


multiExpr[[2]]$data = multiExpr[[2]]$data[-8, ]#remove GSM1023529 outlier

collectGarbage();
# Check the size of the leftover data
exprSize = checkSets(multiExpr)
exprSize


gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
gsg$allOK

if (!gsg$allOK)
{
  # Print information about the removed genes:
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes], 
                                              collapse = ", ")))
  for (set in 1:exprSize$nSets)
  {
    if (sum(!gsg$goodSamples[[set]]))
      printFlush(paste("In set", setLabels[set], "removing samples",
                       paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
    # Remove the offending genes and samples
    multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
    # Update exprSize
  }
  exprSize = checkSets(multiExpr)
}

sampleTrees = list()
for (set in 1:nSets)
{
  sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}


powers = c(seq(4,10,by=1), seq(12,20, by=2));
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets);
# Call the network topology analysis function for each set in turn
for (set in 1:nSets)
  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,networkType = "signed hybrid", corFnc = "bicor",verbose = 5)[[2]]);
collectGarbage();
# Plot the results:
colors = c("black", "red")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
             "Max connectivity");
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (set in 1:nSets)
{
  for (col in 1:length(plotCols))
  {
    ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
    ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
  }
}


sizeGrWindow(8, 6)
pdf(file = "Consensus.scaleFreeAnalysis.pdf", wi = 8, he = 6)
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
for (col in 1:length(plotCols)) for (set in 1:nSets)
{
  if (set==1)
  {
    plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
         main = colNames[col]);
    addGrid();
  }
  if (col==1)
  {
    text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         labels=powers,cex=cex1,col=colors[set]);
  } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
         labels=powers,cex=cex1,col=colors[set]);
  if (col==1)
  {
    legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
  } else
    legend("topright", legend = setLabels, col = colors, pch = 20) ;
}
dev.off();

# Define data set dimensions
nGenes = exprSize$nGenes;
nSamples = exprSize$nSamples;




net.con = blockwiseConsensusModules(
  multiExpr, power = 8, networkType="signed hybrid",blockSizePenaltyPowwer = "inf",
  maxBlockSize =20000, 
  deepSplit = 2,
  TOMType = "signed Nowick", 
  minModuleSize = 30,
  pamRespectsDendro = FALSE, 
  mergeCutHeight = 0.1, numericLabels = TRUE,
  minKMEtoStay = 0,
  saveTOMs = TRUE, verbose = 5)



consMEs = net.con$multiMEs;
moduleLabels.con = net.con$colors;
# Convert the numeric labels to color labels
moduleColors.con = labels2colors(moduleLabels.con)
consTree = net.con$dendrograms[[1]]; 

sizeGrWindow(8,6);
pdf(file = "ConsensusDendrogram-auto.pdf", wi = 8, he = 6)
plotDendroAndColors(consTree, moduleColors.con,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram and module colors")

dev.off()


# Recalculate consMEs to give them color names
consMEsC = multiSetMEs(multiExpr, universalColors = moduleColors.con);
# We add the weight trait to the eigengenes and order them by consesus hierarchical clustering:
MET.con = consensusOrderMEs(consMEsC);

sizeGrWindow(8,10);
pdf(file = "EigengeneNetworks.NTAPvsGEO41747.pdf", width= 8, height = 10);
par(cex = 0.9)
plotEigengeneNetworks(MET.con, setLabels, marDendro = c(0,2,2,1), marHeatmap = c(3,3,2,1),
                      zlimPreservation = c(0.5, 1), xLabelsAngle = 90)
dev.off();

MET.con = consensusOrderMEs(consMEsC)                        

sizeGrWindow(8,10);
pdf(file = "EigengeneNetworks.NTAPvsGEO41747.pdf", width= 8, height = 10);
par(cex = 0.9)
plotEigengeneNetworks(MET.con, setLabels, letterSubPlots = TRUE,colorLabels = FALSE, marDendro = c(0,2,2,1), marHeatmap = c(3,3,2,1),
                      zlimPreservation = c(0.5, 1), xLabelsAngle = 90)
dev.off() 



pdf("TwoNK_modules.NTAPvsGEO41747.pdf",height=8,width=12)
plotDendroAndColors(geneTreeA1, modulesA1, "Modules", dendroLabels=FALSE, hang=0.03, addGuide=TRUE,
                    guideHang=0.05, main="Gene dendrogram and module colors (NATP)") 
plotDendroAndColors(geneTreeA2, modulesA1, "Modules", dendroLabels=FALSE, hang=0.03, addGuide=TRUE,
                    guideHang=0.05, main="Gene dendrogram and module colors (GEO41747)") 

dev.off()


#To determine which concensus networks are preserved better between NTAP cell lines and primary pNF tumors from GEO41747
multiColor.test = list(NTAP=moduleColors.con)
mp=modulePreservation(multiExpr,multiColor.test,referenceNetworks=1,verbose=3,networkType="signed hybrid",
                      nPermutations=30,maxGoldModuleSize=100,maxModuleSize=400) 

stats.con = mp$preservation$Z$ref.NTAP$inColumnsAlsoPresentIn.GEO 

PreservedConNKs<-stats.con[order(-stats.con[,2]),c(1:2)]%>%filter(Zsummary.pres>7) # Zsummary from 10~7 is believed to be conservative according the package authors.  
write.table(stats.con[order(-stats.con[,2]),c(1:2)],"NetworkConservationBetweenNTAPandGEO41747.txt",quote=FALSE,sep="\t")



drug_data <- synGet("syn20700260")$path %>% read.csv() 
head(drug_data)

pnf <- c("ipNF05.5 (single clone)", "ipNF06.2A", "ipNF95.11b C/T", "ipnNF95.11C", "ipNF95.6", "ipNF05.5 (mixed clone)", "ipNF95.11b C")

table(drug_data$response_type)                        

drug_data_filt_1 <- drug_data %>% 
  filter(response_type == "AUC_Simpson") %>% 
  filter(model_name %in% pnf) %>% 
  group_by(drug_screen_id) %>% 
  filter(n() == 1) %>% 
  ungroup()


drug_data_filt <- drug_data_filt_1 %>% 
  group_by(DT_explorer_internal_id) %>% 
  filter(n() > 3) %>% 
  ungroup() %>% 
  dplyr::select(DT_explorer_internal_id, response) %>%
  group_by(DT_explorer_internal_id) %>% 
  summarize('median_response' = median(response))%>% 
  ungroup() 



targets <- synGet("syn17091507")$path %>% readRDS() %>% 
  filter(mean_pchembl > 6) %>% 
  dplyr::select(internal_id, hugo_gene, std_name) %>% 
  distinct()


drug_data_filt_2 <- drug_data_filt_1 %>% 
  dplyr::select(drug_name,DT_explorer_internal_id, response) %>%
  group_by(drug_name) %>% 
  summarize('median_response' = median(response))%>% 
  ungroup()   

#An abitrary cutoff was set as 50 for median_response, and it can be changed. 
TraitsOfDrug <- drug_data_filt_1%>%dplyr::select(model_name,drug_name,response)%>%
  group_by(drug_name)%>%
  unique()%>%spread(model_name,response)%>%
  filter(drug_name %in% drug_data_filt_2$drug_name[drug_data_filt_2$median_response<50])  



############# use DT_explorer_internal_id for durgs

AnnoDrug <- drug_data_filt_1%>%dplyr::select(DT_explorer_internal_id,drug_name)%>%
  group_by(DT_explorer_internal_id)%>%
  unique()


NewAnno<- TraitsOfDrug %>% 
  left_join(AnnoDrug, by= c("drug_name"="drug_name"))  


finalAnno<- NewAnno %>% left_join(targets, by= c("DT_explorer_internal_id"="internal_id"))%>%
  dplyr::select(drug_name,DT_explorer_internal_id,hugo_gene,std_name)                   
               

datTraits= data.frame(t(column_to_rownames(TraitsOfDrug,"drug_name")))
rownames(datTraits)<-c("ipNF05.5 (mixed clone)","ipNF05.5 (single clone)","ipNF06.2A", "ipNF95.11bC","ipNF95.11bC/T", "ipNF95.6","ipnNF95.11c") 
datTraits<-datTraits[-c(3,5),]

Traits = vector(mode="list", length = nSets)
Traits[[1]]$data <-datTraits                
Traits[[2]]$data <-datTraits



# Set up variables to contain the module-trait correlations
moduleTraitCor = list();
moduleTraitPvalue = list();

moduleTraitCor[[1]] = cor(consMEs[[1]]$data, Traits[[1]]$data, use = "p")
moduleTraitPvalue[[1]] = corPvalueFisher(moduleTraitCor[[1]], exprSize$nSamples[1])


# Convert numerical lables to colors for labeling of modules in the plot
MEColors = labels2colors(as.numeric(substring(names(consMEs[[1]]$data), 3)))
MEColorNames = paste("ME", MEColors, sep="")
Color2NKs<-cbind(MEColors,names(consMEs[[1]]$data))
colnames(Color2NKs)[2]<-"NKs"

rownames(moduleTraitCor[[1]])


library(RColorBrewer)
hmcol <- rev(colorRampPalette(c("blue", "white", "red"))(n = 10000))
DrugDist<-dist(t(moduleTraitCor[[1]]),method = "euclidean")
DrugHC<-hclust(DrugDist,method = "average")


pdf("moduleTraitCor.NTAP.inCon-test_lessdrug.pdf", width = 200, height = 8)

heatmap.2(moduleTraitCor[[1]] 
          #Rowv=as.dendrogram(generow.hc), 
          ,Rowv=NULL
          ,Colv=as.dendrogram(DrugHC)
          #,Colv=NULL
          # ,sepwidth=c(0.01,0)
          # ,sepcolor="white"
          # ,colsep=1:ncol( )
          # ,rowsep=1:(nrow())
          ,main="moduleTraitCor.NTAP.inCon-test_lessdrug"
          ,cexRow=0.8,cexCol=0.8
          ,col= rev(hmcol)
          ,key=TRUE,keysize=1, symkey=TRUE, density.info="none", trace="none",scale="none")	
dev.off()


# Open a suitably sized window 
sizeGrWindow(10,7)
pdf(file = "ModuleTraitRelationships-NTAPnew2.pdf", wi = 200, he = 7);
# Plot the module-trait relationship table for set number 1
set = 1
textMatrix =  paste(signif(moduleTraitCor[[set]], 2), "\n(",signif(moduleTraitPvalue[[set]], 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor[[set]])
par(mar = c(6, 8.8, 3, 2.2));
labeledHeatmap(Matrix = moduleTraitCor[[set]],
               xLabels = names(Traits[[set]]$data),
               yLabels = MEColorNames,
               ySymbols = MEColorNames,
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.2,
               zlim = c(-1,1),
               main = paste("Module--trait relationships in", setLabels[set]))
dev.off();




################
#These are drug_name that define the right boundary of DrugClusters. They are determined manually 
#based on the heatmaps "moduleTraitCor.NTAP.inCon-test_lessdrug.pdf"
BoundaryID<-c(
  "NCGC00346882.01",
  "NCGC00346646.01",
  "NCGC00263218.02",
  "NCGC00346658.01",
  "NCGC00185071.02",
  "NCGC00250378.01",
  "NCGC00242596.01")
##create the Drug groups according to the heatmap clustering
j=1;m=1;
for (i in 1:length(BoundaryID)){
  j<-DrugClustDF$x[DrugClustDF$label == BoundaryID[i]];
  DrugClustDF$cluster[m:j]<-paste0("DrugCluster",i);
  m<-j+1
} 
#make the drug_name consistent among anotations
finalAnno$drug_name<-gsub("-", ".", finalAnno$drug_name)


DrugTargetInCluster<-merge(DrugClustDF,finalAnno, by.x="label",by.y="drug_name",all=FALSE)

#Obtain the ENTREZID for GO analysis
symbols <- DrugTargetInCluster$hugo_gene%>%unique()
DrugTargetEnzID<-mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'SYMBOL')%>%stack()
colnames(DrugTargetEnzID)<-c("EnzID","Symbol")

DrugCandidates<- DrugTargetInCluster %>% left_join(DrugTargetEnzID, by= c("hugo_gene"="Symbol")) %>% group_by(EnzID)%>%unique()

##prepare for GO semantic similarity analysis                 
hsGO.MF <- godata('org.Hs.eg.db', ont="MF")
hsGO.BP <- godata('org.Hs.eg.db', ont="BP")
hsGO.CC <- godata('org.Hs.eg.db', ont="CC")

DrugCluster1<-DrugCandidates$EnzID[DrugCandidates$cluster=="DrugCluster1"]
DrugCluster2<-DrugCandidates$EnzID[DrugCandidates$cluster=="DrugCluster2"]
DrugCluster3<-DrugCandidates$EnzID[DrugCandidates$cluster=="DrugCluster3"]
DrugCluster4<-DrugCandidates$EnzID[DrugCandidates$cluster=="DrugCluster4"]
DrugCluster5<-DrugCandidates$EnzID[DrugCandidates$cluster=="DrugCluster5"]
DrugCluster6<-DrugCandidates$EnzID[DrugCandidates$cluster=="DrugCluster6"]
DrugCluster7<-DrugCandidates$EnzID[DrugCandidates$cluster=="DrugCluster7"]


#Generate a list contains the EntrezID for each consensus networks between NTAP and GSE41747 from above analysis
NTAP27nks<-data.frame(cbind(moduleColors.con,colnames(multiExpr[[1]]$data)))
colnames(NTAP27nks)[2]<-"Symbol"
NTAP27symbols <-NTAP27nks$Symbol%>%unique()
NTAP27EnzID<-mapIds(org.Hs.eg.db, NTAP27symbols, 'ENTREZID', 'SYMBOL')%>%stack()
colnames(NTAP27EnzID)<-c("EnzID","Symbol")
NTAP27Candidates<- NTAP27nks %>% left_join(NTAP27EnzID, by= c("Symbol"="Symbol")) %>% group_by(EnzID)%>%unique()


#Create variables with EntrezID for each consensus networks using their corresponding color names
NTAP.symbols <- colnames(multiExpr[[1]]$data)
NTAPEnzID<-mapIds(org.Hs.eg.db,  NTAP.symbols, 'ENTREZID', 'SYMBOL')%>%stack()
colnames(NTAPEnzID)<-c("EnzID","Symbol")
NTAPgeneNetwork<-data.frame(cbind(moduleColors.con,colnames(multiExpr[[1]]$data)))
colnames(NTAPgeneNetwork)[2]<-"Symbol"
NTAPgeneNKs.EnzID<- NTAPgeneNetwork %>% 
  left_join(NTAPEnzID, by= c("Symbol"="Symbol")) %>% group_by(EnzID)%>%unique()
ConNKgroups<-unique(NTAPgeneNKs.EnzID$moduleColors.con)


for (i in 1:length(ConNKgroups)){
  nam<- paste0("NK_",ConNKgroups[i])
  assign(nam, NTAP27Candidates$EnzID[NTAP27Candidates$moduleColors.con==ConNKgroups[i]])
}

#Create a list of containing EntrezID of genes in consensus networks and EntrezID of hugo_symbol of putative drug targets
#for semantic similarity analysis on GO terms
list4Comp<-list()
for (i in 1:length(ConNKgroups)){
  list4Comp[[ConNKgroups[i]]]<-eval(parse(text=paste0("NK_",ConNKgroups[i])))
}


list4Comp[["DrugCluster1"]]<-DrugCluster1
list4Comp[["DrugCluster2"]]<-DrugCluster2
list4Comp[["DrugCluster3"]]<-DrugCluster3
list4Comp[["DrugCluster4"]]<-DrugCluster4
list4Comp[["DrugCluster5"]]<-DrugCluster5
list4Comp[["DrugCluster6"]]<-DrugCluster6
list4Comp[["DrugCluster7"]]<-DrugCluster7


###########################################################
#each line of the below codes may need 5~20 min to complete       
##########################################################
GOsimilarity.MF <-mclusterSim(list4Comp, semData=hsGO.MF, measure="Wang", combine="BMA")
GOsimilarity.BP <-mclusterSim(list4Comp, semData=hsGO.BP, measure="Wang", combine="BMA")
GOsimilarity.CC <-mclusterSim(list4Comp, semData=hsGO.CC, measure="Wang", combine="BMA")

pdf(file = "GOsimilarity.MF.pdf", width= 10, height = 10)
pheatmap(GOsimilarity.MF)
dev.off()   

pdf(file = "GOsimilarity.BP.pdf", width= 10, height = 10)              
pheatmap(GOsimilarity.BP)
dev.off()                         

pdf(file = "GOsimilarity.CC.pdf", width= 10, height = 10)         
pheatmap(GOsimilarity.CC)     
dev.off()


NK4select<-c(rownames(PreservedConNKs))

DrugCl4select<-c("DrugCluster1","DrugCluster2","DrugCluster3","DrugCluster4","DrugCluster5","DrugCluster6")


#The drug candidates should have higher similarity score (>0.8) in one of the three GO branches
GOssMF<-data.frame(GOsimilarity.MF) %>%rownames_to_column(var = "Name") %>% gather(COlNN, Value, -Name) 
GOssMF<-GOssMF%>%filter(Name %in%  NK4select)%>%filter(COlNN %in% DrugCl4select)%>%filter(Value>0.8)%>%mutate("Candi"=paste0(COlNN,".",Name))


GOssBP<-data.frame(GOsimilarity.BP) %>%rownames_to_column(var = "Name") %>% gather(COlNN, Value, -Name) 
GOssBP<-GOssBP%>%filter(Name %in%  NK4select)%>%filter(COlNN %in% DrugCl4select)%>%filter(Value>0.8)%>%mutate("Candi"=paste0(COlNN,".",Name))

GOssCC<-data.frame(GOsimilarity.CC) %>% rownames_to_column(var = "Name") %>% gather(COlNN, Value, -Name) 
GOssCC<-GOssCC%>%filter(Name %in%  NK4select)%>%filter(COlNN %in% DrugCl4select)%>%filter(Value>0.8)%>%mutate("Candi"=paste0(COlNN,".",Name))

##For the drug mechanisms, both strong positive and negative similarities are interesting, so an arbitrary cutoff is used here filter(median(Value)<(-0.1)) 
DrugGeneCor<-data.frame(t(moduleTraitCor[[1]])) %>%
  rownames_to_column(var = "Name") %>%      
  gather(COlNN, Value, -Name) %>% group_by(COlNN)%>%filter(median(Value)<(-0.1))%>%
  left_join(data.frame(Color2NKs),by=c("COlNN"="NKs"))%>%filter(MEColors %in% NK4select)

DrugNKsCandidate<-DrugTargetInCluster%>%dplyr::select(label,cluster,hugo_gene,std_name)%>%left_join(DrugGeneCor,by=c("label"="Name"))%>%group_by(cluster)%>%mutate("Candi"=paste0(cluster,".",MEColors))

##For the drug candidates, the negative correlation between gene networks and drug response is preferred because we hope the drug targeted genes are highly expressed 
#and essential in PN cells (a higher eigengene expression in the networks), and targeting those genes by the drug candidates will lead to the better responses (a smaller AUC). 
DrugNKsCandidate.refine<-DrugNKsCandidate%>%filter(Candi %in% c(GOssMF$Candi,GOssBP$Candi,GOssCC$Candi))%>%filter(Value<0)%>%group_by(label)%>%
  mutate("Score"=n())%>%dplyr::select(label,std_name,cluster,Score)%>%distinct(label,cluster,std_name,Score) %>% arrange(desc(Score))

write.table(DrugNKsCandidate.refine,"DrugNKsCandidate.refine.txt", quote=FALSE,sep="\t")

## A different format of drug candidate rank with clear drug target information is below. 
drug_data_filt_2$drug_name<-gsub("-",".",drug_data_filt_2$drug_name)
DrugNKsCandidate.refine2<-DrugNKsCandidate.refine%>%left_join(drug_data_filt_2,by=c("label"="drug_name"))%>%left_join(finalAnno,by=c("label"="drug_name"))
write.table(DrugNKsCandidate.refine2,"DrugNKsCandidate.refine2.txt", quote=FALSE,sep="\t")  




