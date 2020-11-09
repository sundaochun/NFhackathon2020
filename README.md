# NFhackathon2020
![logo!](/Capture.JPG "Logo")
## Team NFFighters

#Integrated drug target discovery using consensus gene regulation networks and drug response clusters in plexiform neurofibromas 


##Abstract

* Plexiform neurofibromas (PN) are embryonic tumors predisposed with loss of function of the NF1 gene in children and young people. PN are benign tumors but have an increased chance of transforming into malignant peripheral nerve sheath tumors. Besides the recent FDA approved drug, selumetinib, treatment options are very limited. Drug screen using the immortalized PN cell lines provides an efficient way to obtain a pool of drug candidates. However, there are multiple challenges to scrutinize the candidates before clinical trials: 1) in vivo screen has a different tumor microenvironment from in vivo study, which may void the preclinical application, 2) a new chemical has an unknown mechanism and toxicity, 3) a well-studied drug commonly has numerous targets and side-effects, 4) it is hard to evaluate the possibility of the candidate to use as a combination with others. 

* To tackle these issues, 1) we identified the similarities between cell lines and PN tissues by establishing consensus gene regulation networks using the transcriptomes from the Hackathon data and publicly available data set, 2) the drug response data were correlated to the gene networks, and then clustered to reveal the similar pattern among the known and novel drug candidates, 3) the known molecular targets/genes within a drug cluster and the genes within a network were annotated using gene ontology (GO) analysis, respectively, and the similarity between the GOs of a gene network and GOs of a drug cluster were computed using GO semantic analysis to confirm a biologically reliable drug-gene relationship, 4) drug candidates targeting distinct networks can be considered for combination to reduce the toxicities or enhance the effects.  

* We ranked the drug candidates according to drug responses, potential mechanisms of a drug and its targets, and biological consistency between tumor cells and tissue. Considering more criteria such as FDA approval and weighted drug toxicity scores, the drug candidate rank can be further prioritized for a higher success rate in the preclinical and clinical studies.


## Methods

### Data
We used RNA-Seq data of the tumor cell lines in drug screen provided by NFhackathon2020 from synapse and a public PN tissue data (GEO41747) from synapse. 

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
Transcriptomes from both PN cell lines (NTAP) and primary PN tissue (GEO41747) should be fit the "Scale free topology model" and reach the model score >0.8.   
![ScaleFree!](/images/Scale Free Topology Model Fit.JPG "ScaleFree")


