# NFhackathon2020
![logo!](/images/NF-Terminators.jpg "Logo")
## Team: NF Terminators

#Gene network-based drug discovery in plexiform neurofibromas 


##Abstract

* Plexiform neurofibromas (PN) are embryonic tumors predisposed with loss of function of the NF1 gene in children and young people. Besides the recent FDA approved drug, SUNITINIB, treatment options are minimal. Drug screen using the immortalized PN cell lines provides an efficient way to obtain a pool of drug candidates. However, there are multiple challenges to scrutinize the candidates before clinical-related studies: 1) in vitro screen has a different tumor microenvironment from in vivo tests, which yields a low success rate in the latter, 2) a novel chemical will have an unknown mechanism and targets, 3) a “well-studied” drug commonly has numerous targets and side-effects, 4) it is hard to evaluate the possibility of the drug combination. 

* To tackle these issues, 1) we established the consensus gene regulation networks among cell lines and PN tissues by using the transcriptomes from Hackathon data and publicly available data set, 2) the drug response data were correlated to the preserved gene networks and then clustered to reveal the similar pattern among all drug candidates, 3) the known drug target genes within drug clusters and the genes within preserved networks were annotated using gene ontology (GO) analysis, respectively, and the similarities among the GO terms of a gene network and a drug cluster were computed using GO semantic analysis to confirm a biologically reliable drug-gene relationship, 4) drug candidates can be combined according to their drug clusters and correlated gene networks to reduce the toxicities or enhance the effects.    

* We ranked the drug candidates according to drug responses, potential mechanisms of a drug and its targets, and biological consistency between tumor cells and tissue. Considering more criteria such as FDA approval and weighted drug toxicity scores, the drug candidate rank can be further prioritized for a higher success rate in the preclinical and clinical studies.


## Methods

### Data
We used RNA-Seq data of the tumor cell lines in drug screen provided by NFhackathon2020 from synapse and a public PN tissue data (GSE41747) from synapse. 

### Data processing and visualization 
The transcriptomes from different sources were trimed and merged. WGCNA package was used for the network analysis and GOSemSim package was used for GO sematic analysis in R environment. 

## Results

Our approach involved the following steps:
1. Trim the data to satisfy the downstream analysis
2. Define the preserved consensus networks between PN cell lines and primary PN tissue. 
3. Identify the correlations among drug responses and eigengene expression of regulatory networks. 
4. Define the drug clusters according to the pattern of correlations among drugs and networks.
5. Analyze GO of the genes in networks and the known drug targeted genes in drug clusters, and score the similarity according to the GO semantic analysis
6. Rank the drug candidates considering the responses, correlation with the preserved networks, and known drug targets/mechanisms

### Step 1: Dataset trimming and preparation 
Transcriptomes from both PN cell lines (NTAP) and primary PN tissue (GEO41747) should be fit the "Scale free topology model" and reach the model fit score >0.8.   
![ScaleFree!](/images/ScaleFreeTopologyModelFit.JPG  "ScaleFree")

### Step 2: Define the preserved consensus networks between PN cell lines and PN tissue
The discrepancy between in vivo and in vitro systems is expected. Considering the intertumoral heterogeneity among cells, we should narrow down the drug targets to the conserved gene regulatory network. WGCNA package is used to find the networks and define the eigengenes of each network in PN cell lines (NTAP), PN primary tumor (GEO41747), and the preserved ones shared by these two models.

![NTAP eigenene!](/images/EigengeneNetwork_NTAP.jpg "NTAP eigengene")
![GEO eigenene!](/images/EigengeneNetwork_GEO41747.jpg "GEO eigengene")
![Preserved networks!](/images/PreservedNetworks.jpg "Preserved networks")
The color blocks on the x-axis and y-axis of the heatmaps represent the different shared consensus gene networks.

### Step 3: Correlate the drug responses (AUC_Simpson) and eigengene expression of regulatory networks.
Pearson correlations were calculated between the eigengene of preserved networks and the drug responses among the 5 PN cell lines that have the RNAseq data. A part of result is demonstrated below. Each row is preserved regulatory networks and each column is a drug in the screen. To narrow down the candidates, the arbitrary cutoff of median_response<50 was used, and can be modified accordingly. The results provide the unique "fingerprint" of each drug in the screen.

![Correlations drug&networks!](/images/CorrelationDrug&Networks.JPG "Correlations drug&networks")

### Step 4: Define the drug clusters according to the pattern of correlations among drugs and networks.
A further unsupervised clustering on the column was done to group the drugs with similar pattern, which is assumed to have similar biology in PN cells. For those chemicals has unknown targets, the results will provide insights to demonstrate their potential targets or targeted networks according to the well-studied chemicals with the similar pattern (drug fingerprints) in the group. The six durg clusters were defined according to the patterns. In the heatmap, each row is preserved regulatory networks and each column is a drug in the screen.

![Pattern drug&networks!](/images/ModuleDrugResponsetCorrelation.jpg "Pattern drug&networks")

### Step 5: Analyze GO of the genes in networks and the known drug targeted genes in drug clusters, and score the similarity according to the GO semantic analysis.
With the gene ontology annotations, the enriched GO terms within each of the preserved networks were determined. With the drug annotations from the hackathon2020, drugs with known target genes were also analyzed for enriched GO terms. Within each of the three major branches of GO, Molecular Function (MF), Biology Process (BP), and Cell Component (CC), the semantic distances (similarity) among the GO of regulatory networks and GO of drug clusters were calculated. 

The correlation heatmap under the MF branch is demonstrated below. We can see the drug clusters are similar among themselves, partially because that current drug design is to commonly target functional kinases or growth factor receptors. Interestingly, we also observe the similarities between drug clusters and preserved networks, which provides a new layer of insight for the drug mechanisms. Furthermore, the similarity will enhance or verify the biology of the drugs. With the further characterization of the biology of the preserved networks, drugs strongly correlated with distinct functional networks can be used in combination to enhance the inhibitory effects or reduce the toxicities in preclinical studies. Here, to narrow down the candidates, only the similary value more than 0.8 in any of the three GO branches will be used for next step.

![GOsimilarity.MF!](/images/DrugTargets&PresevedNetwork.GOsimilarity.MF.jpg "GOsimilarity.MF")
![GOsimilarity.MF.score!](/images/DrugNetworkScores.jpg "GOsimilarity.score.MF")

### Step 6: Rank the drug candidates considering the responses, correlation with the preserved networks, and known drug targets/mechanisms
To rank the candidates, with the consideration of the above steps, we added the weight to the drugs with known targets that the more targets the drug has, the higher it may rank. We did so because the drugs with known targets were relatively well-studied and even approved by the FDA in different diseases, making the potential toxicity, preclinical experiment design, and mechanism studies more efficient. The code is also easy to be tweaked to prioritize the new drug candidates, and our analysis in Step 3,4 and 5 will extraordinarily helpful for the identification of potential targeted genes or networks. The top 10 candidates are shown below. Surprisingly, the sunitinib ranks the first, which is the first approved drug for PN by FDA this year.

![Top10Rank!](/images/Rank.jpg "Top10Rank")


## Conclusion/Discussion:
* Our strategy works efficiently to rank the candiates, and the top candidate, SUNITINIB, is the first FDA approved drug on PN.
* Using the preserved networks, we enrich candidates targeting both cells and tumors.
* The algorithm currently weights higher on the “well-annotated” drug candidates and can be tweaked to identify novel candidates.
* The gene networks can help researchers to nail down the mechanisms.
* The drug clusters can serve as a start point to explore the targets of a novel candidate.
* Drugs can be combined according to their clusters.
 
### Additional Questions:
#### 1. More available transcriptome data of PN related cell lines used in screening will enhance the prediction of gene regulatory network significantly. One of the immortalized Schwann line with the loss of one NF1 allele were include in the network analysis to increase the sample number. Actually, from the multiple dimension reduction analysis, the transcriptional differences between PN cell lines and "normal" Schwann cells are not dramtically distinct. There is a cell line name discrepancy between the RNAseq data and drug screen data. 

#### 2. The genetic information such as Exome and SNP data were not used in the current analysis because 1) there are limited sample numbers in PN to statistically significant, 2) the PN has been reported to have low mutation burdens, suggesting the individual mutation or SNPs may not significantly influence the tumor biology and drug responses, 3) the principle of this analysis is to identify the preserved networks among both cell lines and tumors, which decreases the weight of gene mutations with in each sample. However, it doesn't mean the mutations are not important in drug responses, and it can be tackled given more time and available samples. 

#### 3. What the next steps to select the candidate drugs for the best potentials for preclinical study?
* Factor in the weighted toxicity scores
* Use FDA drug label database to design a combination of candidates 
* MPNST drug screen can be integrated for candidates to inhibit both 

#### 4. Reproduction: 
* PN cell line transcriptomes and drug responses data are provided by NFhackathon2020 and downloaded from Synapse. Extra approval may be required to access the data.
* "GSE41747_expVal.tsv" is downloaded from Synapse with ID:syn6130081 and "GSE41747_phenotype_data.tsv" is downloaded from Synapse with ID: syn6130082.
* The R code can be found in /R Codes/DrugCandidates.R 
* Docker image: https://hub.docker.com/repository/docker/sundaochun/daochunproject  Note: 1）the docker image crashes at line 545 when I test it on a iMac pro (2017) with MacOS10.14.6, but runs smoothly on a PC server with windows10. Line 545,546 and 547 require huge compuatation resource and could take up to one hour to finish. The R codes can run smoothly in Rstudio Version 1.2.5033 with R-4.0.2 on both of computers without using the docker container image; 2) With unknow reasons, when generating the GOsimilarity plot in pdf, the code lines 551~553, 555~557, 559~561 may need to run twice to output correct plots.

### This is my first time to join NFhackathon and share codes through Github. Please provide feedback and join the efforts ending the NF.   
