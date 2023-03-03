# awesome-multi-omics

A [community-maintained](https://github.com/mikelove/awesome-multi-omics/graphs/contributors) list of software packages for multi-omics data analysis.

While many of the packages here are marketed for "omics" data (transcriptomics, proteomics, etc.), other more general terms for this type of data analysis are: 

* multi-modal
* multi-table
* multi-way

The common thread among the methods listed here is that the same samples are measured across different assays. The data can be described as multiple matrices/tables with the same number of samples and varying number of features.

The repo is in the style of Sean Davis'
[awesome-single-cell](https://github.com/seandavi/awesome-single-cell)
repo for single-cell analysis methods.

[Contributions welcome](https://github.com/mikelove/awesome-multi-omics/blob/master/CONTRIBUTING.md)...

For brevity, below lists only the first author of multi-omics methods.

## Software packages and methods

### Multi-omics correlation or factor analysis

- 2007 - **SCCA** - Parkhomenko - sparse CCA - [paper 1](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2367499/), [paper 2](https://doi.org/10.2202/1544-6115.1406) 
- 2008 - **PCCA** - Waaijenborg - penalized CCA / CCA-EN - [paper](https://doi.org/10.2202/1544-6115.1329)
- 2009 - [PMA](https://CRAN.r-project.org/package=PMA) - Witten - Sparse Multi CCA - [paper 1](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2697346/), [paper 2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2861323/)
- 2009 - **sPLS** - Lê Cao - sparse PLS - [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2640358/)
- 2009 - [gesca](https://CRAN.r-project.org/package=gesca) - Hwang - RGSCA regularized generalized structured component analysis - [paper](https://doi.org/10.1007/s11336-009-9119-y) 
- 2010 - **Regularized dual CCA** - Soneson - [paper](https://doi.org/10.1186/1471-2105-11-191)
- 2011 - [RGCCA](https://cran.r-project.org/package=RGCCA) - Tenenhaus - Regularized Generalized CCA and Sparse Generalized CCA - [paper 1](https://www.ncbi.nlm.nih.gov/pubmed/28536930), [paper 2](https://www.ncbi.nlm.nih.gov/pubmed/24550197)
- 2011 - **SNMNMF** - Zhang - Sparse Network-regularized Multiple Non-negative Matrix Factorization - [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3117336/)
- 2011 - [scca](https://github.com/tomwhoooo/scca_3.0) - Lee - Sparse Canonical Covariance Analysis for High-throughput Data - [paper](https://doi.org/10.2202/1544-6115.1638)
- 2012 - [STATIS/DiSTATIS](https://github.com/HerveAbdi/DistatisR) - Abdi - structuring three-way statistical tables - [paper](https://doi.org/10.1002/wics.198)
- 2012 - **joint NMF** - Zhang - extension of NMF to multiple datasets - [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3479191/)
- 2012 - **sMBPLS** - Li - sparse MultiBlock Partial Least Squares - [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3463121/)
- 2012 - **Bayesian group factor analysis** - Virtanen - [paper](http://proceedings.mlr.press/v22/virtanen12.html)
- 2012 - [RIMBANET](http://research.mssm.edu/integrative-network-biology/RIMBANET/RIMBANET_overview.html) - Zhu - Reconstructing Integrative Molecular Bayesian Networks - [paper](https://doi.org/10.1371/journal.pbio.1001301)
- 2013 - [FactoMineR](https://cran.r-project.org/package=FactoMineR) - Abdi - MFA: multiple factor analysis - [paper](https://doi.org/10.1002/wics.1246)
- 2013 - [JIVE](https://genome.unc.edu/jive/) - Lock - joint & individual variance explained - [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3671601/)
- 2013 - [pandaR](https://bioconductor.org/packages/release/bioc/html/pandaR.html) - Schlauch - Passing Attributes between Networks for Data Assimilation - [paper](https://doi.org/10.1093/bioinformatics/btx139)
- 2014 - [omicade4](https://bioconductor.org/packages/omicade4) - Meng - MCIA: multiple co-interia analysis - [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4053266/)
- 2014 - [STATegRa](https://bioconductor.org/packages/STATegRa) - Planell - DISCO, JIVE, & O2PLS - [paper](https://doi.org/10.3389/fgene.2021.620453)
- 2014 - **Joint factor model** - Ray - [paper](https://doi.org/10.1093/bioinformatics/btu064)
- 2014 - [GFAsparse](https://research.cs.aalto.fi/pml/software/GFAsparse/) - Khan - group factor analysis sparse [paper 1](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4147909/), [paper 2](https://doi.org/10.1093/bioinformatics/btw207)
- 2015 - **Sparse CCA** - Gao (3rd paper first author is Chen) - [paper 1](https://doi.org/10.1214/15-AOS1332), [paper 2](https://doi.org/10.1214/16-AOS1519), [paper 3](https://arxiv.org/abs/1311.6186)
- 2015 - [CCAGFA](https://cran.r-project.org/package=CCAGFA) - Klami - Bayesian Canonical Correlation Analysis and Group Factor Analysis - [paper 1](https://doi.org/10.1109/TNNLS.2014.2376974), [paper 2](http://www.jmlr.org/papers/v18/16-509.html)
- 2016 - [CMF](https://cran.r-project.org/package=CMF) - Klami - collective matrix factorization - [paper](https://arxiv.org/abs/1312.5921)
- 2016 - [moGSA](https://bioconductor.org/packages/mogsa) - Meng - multi-omics gene set analysis - [paper](https://doi.org/10.1101/046904)
- 2016 - [iNMF](https://github.com/yangzi4/iNMF) - Yang - integrative NMF - [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/pmid/26377073/)
- 2016 - [BASS](https://github.com/judyboon/BASS) - Zhao - Bayesian group factor analysis - [paper](https://arxiv.org/abs/1411.2698)
- 2016 - `imputeMFA` in [missMDA](https://cran.r-project.org/web/packages/missMDA/index.html) - Voillet - multiple imputation for multiple factor analysis (MI-MFA) - [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5048483/)
- 2016 - [PLSCA](https://github.com/derekbeaton/PLSCA_Framework) - Beaton - Partial Least Square Correspondence Analysis - [paper](https://doi.org/10.1037/met0000053)
- 2017 - [mixOmics](https://bioconductor.org/packages/mixOmics) - Rohart - various methods - [paper1](https://doi.org/10.1371/journal.pcbi.1005752), [paper2](https://doi.org/10.1093/bioinformatics/bty1054)
- 2017 - [mixedCCA](https://github.com/irinagain/mixedCCA) - Yoon - sparse CCA for data of mixed types - [paper](https://arxiv.org/abs/1807.05274)
- 2017 - [SLIDE](https://github.com/irinagain/SLIDE_Rpackage) - Gaynanova - Structural Learning and Integrative Decomposition of Multi-View Data - [paper](https://arxiv.org/abs/1707.06573)
- 2017 - [fCCAC](https://github.com/pmb59/fCCAC/) - Madrigal - functional canonical correlation analysis to evaluate covariance - [paper](https://doi.org/10.1093/bioinformatics/btw724)
- 2017 - [TSKCCA](https://github.com/kosyoshida/TSKCCA) - Yoshida - Sparse kernel canonical correlation analysis - [paper](https://doi.org/10.1186/s12859-017-1543-x)
- 2017 - **SMSMA** - Kawaguchi - Supervised multiblock sparse multivariable analysis - [paper](https://doi.org/10.1093/biostatistics/kxx011)
- 2018 - [AJIVE](https://github.com/idc9/r_jive) - Feng - angle-based JIVE - [paper](https://arxiv.org/abs/1704.02060)
- 2018 - [MOFA](https://github.com/bioFAM/MOFA) - Argelaguet - multi-omics factor analysis - [paper 1](https://doi.org/10.15252/msb.20178124), [paper 2](https://www.biorxiv.org/content/10.1101/837104v1), [application](https://doi.org/10.1101/519207)
- 2018 - [PCA+CCA](https://github.com/pachterlab/PCACCA/) - Brown - [paper](https://doi.org/10.1371/journal.pgen.1007841)
- 2018 - [JACA](https://github.com/Pennisetum/JACA) - Zhang - Joint Association and Classification Analysis - [paper](https://arxiv.org/abs/1811.08511)
- 2018 - **iPCA** - Tang - Integrated Principal Components Analysis - [paper](https://arxiv.org/abs/1810.00832)
- 2018 - [pCIA](https://www.med.upenn.edu/long-lab/software.html) - Min - penalized COI - [paper](https://www.ncbi.nlm.nih.gov/pubmed/30165424)
- 2018 - **sSCCA** - Safo - structured sparse CCA - [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5677597/)
- 2018 - **SWCCA** - Min - Sparse Weighted CCA - [paper](https://arxiv.org/abs/1710.04792)
- 2018 - [OmicsPLS](https://github.com/selbouhaddani/OmicsPLS) - Bouhaddani  - O2PLS implemented in R, with an alternative cross-validation scheme - [paper](https://doi.org/10.1186/s12859-018-2371-3)
- 2018 - [SCCA-BC](https://github.com/pimentel/scca-bc) - Pimentel - Biclustering by sparse canonical correlation analysis - [paper](https://doi.org/10.1007/s40484-017-0127-0)
- 2018 - [mixKernel](https://cran.r-project.org/package=mixKernel) - Mariette - kernel method for unsupervised multi-omics integration - [paper 1](http://dx.doi.org/10.1093/bioinformatics/btx682), [paper 2](http://dx.doi.org/10.1093/nargab/lqac014)
- 2019 - [WON-PARAFAC](https://github.com/NKI-CCB/won-parafac) - Kim - weighted orthogonal nonnegative parallel factor analysis - [paper](https://doi.org/10.1038/s41467-019-13027-2) 
- 2019 - [BIDIFAC](https://github.com/lockEF/bidifac) - Park - bidimensional integrative factorization - [paper 1](https://doi.org/10.1111/biom.13141), [paper 2](https://arxiv.org/abs/2002.02601)
- 2019 - [SmCCNet](https://cran.r-project.org/web/packages/SmCCNet/index.html) - Shi - sparse multiple canonical correlation network analysis - [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6931269/)
- 2020 - [msPLS](https://github.com/acsala/sPLSPM) - Csala - multiset sparse partial least squares path modeling - [paper](https://doi.org/10.1186/s12859-019-3286-3)
- 2020 - **MOTA** - Fan - network-based multi-omic data integration for biomarker discovery - [paper](https://doi.org/10.3390/metabo10040144)
- 2020 - [D-CCA](https://github.com/shu-hai/D-CCA) - Shu - Decomposition-based Canonical Correlation Analysis - [paper](https://doi.org/10.1080/01621459.2018.1543599)
- 2020 - [COMBI](https://bioconductor.org/packages/combi) - Hawinkel - Compositional Omics Model-Based Integration - [paper](https://doi.org/10.1093/nargab/lqaa050)
- 2020 - [DPCCA](https://github.com/gwgundersen/dpcca) - Gundersen - Deep Probabilistic CCA - [paper](http://proceedings.mlr.press/v115/gundersen20a.html)
- 2020 - [MEFISTO](https://biofam.github.io/MOFA2/MEFISTO.html) - Velten - spatial or temporal relationships - [preprint](https://doi.org/10.1101/2020.11.03.366674)
- 2020 - [MultiPower](https://github.com/ConesaLab/MultiPower) - Tarazona - Sample size in multi-omic experiments - [paper](https://doi.org/10.1038/s41467-020-16937-8)
- 2020 - [mixedCCA](https://cran.r-project.org/web/packages/mixedCCA/) - Yoon - Sparse semiparametric CCA for data of mixed types - [paper](https://doi.org/10.1093/biomet/asaa007)

### Ecology multi-table literature

- 1994 - **COI** - Doledec - Co‐inertia analysis - [paper](https://doi.org/10.1111/j.1365-2427.1994.tb01741.x)
- 2007 - [ade4](https://CRAN.r-project.org/package=ade4) - Dray - Implementing the Duality Diagram for Ecologists - [paper](http://dx.doi.org/10.18637/jss.v022.i04)

### Chemometrics multi-table literature

- 1987 -    - Wold - Multi‐way principal components‐and PLS‐analysis - [paper](https://doi.org/10.1002/cem.1180010107)
- 1996 -    - Wold - Hierarchical multiblock PLS - [paper](https://doi.org/10.1002/(SICI)1099-128X(199609)10:5/6%3C463::AID-CEM445%3E3.0.CO;2-L)
- 2003 -    - Trygg - O2‐PLS, a two‐block (X–Y) latent variable regression (LVR) - [paper](https://doi.org/10.1002/cem.775)
- 2011 -    - Hanafi - Connections between multiple COI and consensus PCA - [paper](https://doi.org/10.1016/j.chemolab.2010.05.010)
- 2015 - [THEME](https://github.com/ThomData/R_THEME) - Verron - THEmatic Model Exploration - [paper](https://doi.org/10.1002/cem.2759)

### Behavioral research multi-table literature

- 2013 - DISCO SCA - Schouteden - distinctive and common components with simultaneous-component analysis - [paper 1](https://www.ncbi.nlm.nih.gov/pubmed/23361416), [paper 2](https://www.ncbi.nlm.nih.gov/pubmed/24178130)

### Multi-omics clustering / classification / prediction

*Note: I think that prediction of genomic tracks, e.g. ChIP-seq, from other genomic tracks is a large area of research that may deserve a separate repository. Below are methods for clustering / classification of samples into sub-types or prediction of outcomes.*

- 2009 - [iCluster](https://cran.r-project.org/package=iCluster) - Shen - [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2800366/)
- 2011 - [PSDF](https://sites.google.com/site/patientspecificdatafusion/home/patientSpecificDataFusion.tar.gz) - Yuan - Data fusion by Bayesian nonparametric Dirichlet modeling - [website](https://sites.google.com/site/patientspecificdatafusion/), [publication](https://doi.org/10.1371/journal.pcbi.1002227)
- 2012 - [MDI](https://warwick.ac.uk/fac/cross_fac/zeeman_institute/zeeman_research/software/) - Kirk - [paper1](https://academic.oup.com/bioinformatics/article/28/24/3290/244641), [paper2](https://www.degruyter.com/document/doi/10.1515/sagmb-2015-0055/html)
- 2013 - [iClusterPlus](https://bioconductor.org/packages/iClusterPlus) - Mo - [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3600490/)
- 2013 - [BCC](https://github.com/ttriche/bayesCC) - Lock - Bayesian consensus clustering - [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3789539/)
- 2013 - [iBAG](https://github.com/umich-biostatistics/iBAG) - Wang - Integrative Bayesian Analysis of Genomics - [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3546799/)
- 2014 - [SNF](http://compbio.cs.toronto.edu/SNF/SNF/Software.html) - Wang - [paper](https://www.ncbi.nlm.nih.gov/pubmed/24464287)
- 2015 - moCluster - Meng - Derivative of iClusterPlus - [paper](https://doi.org/10.1021/acs.jproteome.5b00824)
- 2017 - [clusternomics](https://cran.r-project.org/web/packages/clusternomics/index.html) - Gabasova - [paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005781)
- 2017 - PFA - Shi - Pattern Fussion analysis - [publication](https://doi.org/10.1093/bioinformatics/btx176)
- 2019 - [IBOOST](http://dlin.web.unc.edu/software/iboost/) - Wong - [paper](https://doi.org/10.1186/s13059-019-1640-4)
- 2019 - [Spectrum](https://cran.r-project.org/web/packages/Spectrum/index.html) - John - [paper](https://doi.org/10.1093/bioinformatics/btz704)
- 2019 - [NEMO](https://github.com/Shamir-Lab/NEMO) - Rappoport - Similarity-based Clustering - [paper](https://doi.org/10.1093/bioinformatics/btz058)
- 2020 - [INF](https://gitlab.fbk.eu/MPBA/INF) - Chierici and Bussola - [paper](https://doi.org/10.1101/2020.04.01.020685)
- 2021 - [ClustOmics](https://github.com/galadrielbriere/ClustOmics) - Brière - Consensus clustering - [paper](https://doi.org/10.1186/s12859-021-04279-1)
- 2021 - [MOGONET](https://github.com/txWang/MOGONET) - Tongxin Wang - Multi-Omics Graph cOnvolutional NETworks - [paper](https://pubmed.ncbi.nlm.nih.gov/34103512/)

### Multi-omics autoencoders

- 2019 - [maui](https://github.com/BIMSBbioinfo/maui) - Ronen - Stacked VAE + clustering predictive of survival - [paper](https://doi.org/10.26508/lsa.201900517)
- 2019 - [IntegrativeVAEs](https://github.com/CancerAI-CL/IntegrativeVAEs) - Simidjievski - Variational autoencoders + classification - [paper](https://doi.org/10.3389/fgene.2019.01205)
- 2019 - [OmiVAE](https://github.com/zhangxiaoyu11/OmiVAE) - Xiaoyu Zhang -  Integrated Multi-omics Analysis Using Variational Autoencoders - [paper](https://arxiv.org/abs/1908.06278)
- 2021 - [DeepProg](https://github.com/lanagarmire/DeepProg) - Poirion - DL and ML ensemble + survival prediction - [paper](https://doi.org/10.1186/s13073-021-00930-x)
- 2021 - [SHAE](https://github.com/BoevaLab/Supervised-hierarchical-autoencoders-for-cancer-survival) - Wissel - Supervised Hierarchical Autoencoder + survival prediction - [preprint](https://doi.org/10.1101/2021.09.16.460589)

### Multi-omics networks

- 2018 - [MolTi-DREAM](https://github.com/gilles-didier/MolTi-DREAM/) - Didier - identifying communities from multiplex networks, and annotated the obtained clusters [article](https://dx.doi.org/10.12688%2Ff1000research.15486.2)
- 2018 [NetICS](https://github.com/cbg-ethz/netics) - Christos Dimitrakopoulos - Network-based integration of multi-omics data for prioritizing cancer genes - [paper](https://pubmed.ncbi.nlm.nih.gov/29547932/)
- 2019 - [RWR-MH](https://github.com/alberto-valdeolivas/RWR-MH) - Valdeolivas - Random walk with restart on multiplex and heterogeneous biological networks [article](https://doi.org/10.1093/bioinformatics/bty637)
- 2020 - [MOGAMUN](https://bioconductor.org/packages/MOGAMUN/) - Novoa-del-toro - A multi-objective genetic algorithm to find active modules in multiplex biological networks [preprint](https://www.biorxiv.org/content/10.1101/2020.05.25.114215v1)
- 2021 - [RWRF](https://github.com/Sepstar/RWRF/) - Wen -  Random Walk with Restart for multi-dimensional data Fusion [paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04029-3) 

### Single cell multi-omics

- 2018 - [cardelino](https://github.com/PMBio/cardelino) -    - gene expression states to clones (SNVs from scRNA-seq + bulk exome data) - 
- 2018 - [clonealign](https://github.com/kieranrcampbell/clonealign) - Campbell - gene expression states to clones (scRNA-seq + scDNA-seq (CNV)) - [paper](https://doi.org/10.1101/344309)
- 2020 - [CiteFuse](https://sydneybiox.github.io/CiteFuse/) - Kim - CITE-seq data analysis [paper](https://doi.org/10.1093/bioinformatics/btaa282)
- 2021 - [CoSpar](https://cospar.readthedocs.io/) - Wang - infer dynamics by integrating state and lineage information - [paper](https://www.biorxiv.org/content/10.1101/2021.05.06.443026v1)

### Multi-study correlation or factor analysis

- 2016 - [MSFA](https://github.com/rdevito/MSFA) - De Vito - multi-study factor analysis: same features, different samples - [paper](https://arxiv.org/abs/1611.06350)

### Multi-omics simulation

- 2016 - [InterSIM](https://cran.r-project.org/package=InterSIM) - Chalise - methylation, gene expression and protein expression - [paper](https://dx.doi.org/10.1016%2Fj.cmpb.2016.02.011)
- 2019 - [MOSim](https://bioconductor.org/packages/MOSim) -  Martinez-Mira - RNA-seq, ATAC-seq (DNase-seq), ChIP-seq, small RNA-seq and Methyl-seq. - [preprint](http://dx.doi.org/10.1101/421834)
- 2019 - [OmicsSIMLA](https://omicssimla.sourceforge.io/) - Chung - DNA, CNV, WGBS, RNAseq, Protein expression - [paper](https://doi.org/10.1093/gigascience/giz045)

## Multi-omics reviews / evaluations

- 2008 - Holmes - [Multivariate data analysis: The French way](https://projecteuclid.org/euclid.imsc/1207580085)
- 2014 - Kohl - [A practical data processing workflow for multi-OMICS projects](https://doi.org/10.1016/j.bbapap.2013.02.029)
- 2016 - Josse - [Measuring multivariate association and beyond](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5658146/)
- 2016 - Ebrahim - [Multi-omic data integration enables discovery of hidden biological regularities](https://doi.org/10.1038/ncomms13091)
- 2016 - Meng - [Dimension reduction techniques for the integrative analysis of multi-omics data](https://doi.org/10.1093/bib/bbv108)
- 2016 - Li - [A review on machine learning principles for multi-view biological data integration](https://doi.org/10.1093/bib/bbw113)
- 2017 - Huang - [More Is Better: Recent Progress in Multi-Omics Data Integration Methods](https://doi.org/10.3389/fgene.2017.00084)
- 2017 - Hasin - [Multi-omics approaches to disease](https://doi.org/10.1186/s13059-017-1215-1)
- 2017 - Allen - [Statistical data integration: Challenges and opportunities](http://www.stat.rice.edu/~gallen/gallen_data_integration_2017.pdf)
- 2018 - Rappoport - [Multi-omic and multi-view clustering algorithms: review and cancer benchmark](https://doi.org/10.1093/nar/gky889)
- 2018 - Bougeard - [Current multiblock methods: Competition or complementarity? A comparative study in a unified framework](https://doi.org/10.1016/j.chemolab.2018.09.003)
- 2018 - Karczewski - [Integrative omics for health and disease](https://doi.org/10.1038/nrg.2018.4)
- 2018 - Yan - [Network approaches to systems biology analysis of complex disease: integrative methods for multi-omics data](https://dx.doi.org/10.1093/bib/bbx066)
- 2019 - Misra - [Integrated omics: tools, advances and future approaches](https://doi.org/10.1530/JME-18-0055)
- 2019 - Chauvel - [Evaluation of integrative clustering methods for the analysis of multi-omics data](https://doi.org/10.1093/bib/bbz015)
- 2019 - McCabe - [Consistency and overfitting of multi-omics methods on experimental data](https://doi.org/10.1093/bib/bbz070) - [code](https://github.com/mccabes292/movie)
- 2019 - Pierre-Jean - [Clustering and variable selection evaluation of 13 unsupervised methods for multi-omics data integration](https://doi.org/10.1093/bib/bbz138)
- 2019 - Pinu - [Systems Biology and Multi-Omics Integration: Viewpoints from the Metabolomics Research Community](https://doi.org/10.3390/metabo9040076)
- 2019 - Wu - [A Selective Review of Multi-Level Omics Data Integration Using Variable Selection](https://doi.org/10.3390/ht8010004)
- 2019 - Sankaran - [Multitable methods for microbiome data integration](https://doi.org/10.3389/fgene.2019.00627) - [code](https://github.com/krisrs1128/multitable_review)
- 2020 - Lee - [Heterogeneous Multi-Layered Network Model for Omics Data Integration and Analysis](https://doi.org/10.3389/fgene.2019.01381)
- 2020 - Herrmann - [Large-scale benchmark study of survival prediction methods using multi-omics data](https://arxiv.org/abs/2003.03621) - [code](https://github.com/HerrMo/multi-omics_benchmark_study)
- 2020 - Nguyen - [Multiview learning for understanding functional multiomics](https://doi.org/10.1371/journal.pcbi.1007677)
- 2020 - Eicher - [Metabolomics and multi-omics integration: a survey of computational methods and resources](https://doi.org/10.3390/metabo10050202)
- 2020 - Cantini - [Benchmarking joint multi-omics dimensionality reduction approaches for cancer study](https://www.biorxiv.org/content/10.1101/2020.01.14.905760v1)
- 2020 - Subramanian - [Multi-omics Data Integration, Interpretation, and Its Application](https://dx.doi.org/10.1177/1177932219899051)
- 2020 - Krassowski - [State of the Field in Multi-Omics Research: From Computational Needs to Data Mining and Sharing](https://doi.org/10.3389%2Ffgene.2020.610798) - [code](https://github.com/krassowski/multi-omics-state-of-the-field)
- 2021 - Espinosa - [Data-Driven Modeling of Pregnancy-Related Complications](https://doi.org/10.1016/j.molmed.2021.01.007)
- 2022 - Jiang - [Uncovering Cross-Cohort Molecular Features with Multi-Omics Integration Analysis](https://doi.org/10.1101/2022.11.10.515908)

## Multi-omics application papers

- 2007 - Fagan - [A multivariate analysis approach to the integration of proteomic and gene expression data](https://doi.org/10.1002/pmic.200600898)
- 2011 - De la Cruz - [The duality diagram in data analysis: Examples of modern applications](https://doi.org/10.1214/10-AOAS408) - [R notebook](http://lbbe-shiny.univ-lyon1.fr/Reproducible_Research/06-AAS.Thioulouse.2011/)
- 2014 - Tomescu - [Integrative omics analysis. A study based on Plasmodium falciparum mRNA and protein data](https://doi.org/10.1186/1752-0509-8-S2-S4)
- 2014 - Costello (NCI/DREAM) - [A community effort to assess and improve drug sensitivity prediction algorithms](https://doi.org/10.1038/nbt.2877)
- 2015 - Wang - [Inferring gene–gene interactions and functional modules using sparse canonical correlation analysis](https://doi.org/10.1214/14-AOAS792)
- 2016 - Wan - [TCGA2STAT: simple TCGA data access for integrated statistical analysis in R](https://doi.org/10.1093/bioinformatics/btv677) - [R notebook](http://www.liuzlab.org/TCGA2STAT/)
- 2017 - Butler - [Integrating single-cell transcriptomic data across different conditions, technologies, and species.](https://www.biorxiv.org/content/10.1101/164889v1)
- 2018 - Skelly - [Reference trait analysis reveals correlations between gene expression and quantitative traits in disjoint samples](https://www.biorxiv.org/content/10.1101/489542v1) - [R notebook](https://daskelly.github.io/reference_traits/reference_trait_analysis_walkthrough.html)
- 2018 - Stuart - [Comprehensive integration of single cell data](https://www.biorxiv.org/content/10.1101/460147v1)
- 2018 - Ash - [Joint analysis of gene expression levels and histological images identifies genes associated with tissue morphology](https://doi.org/10.1101/458711)
- 2019 - Xu - [Identifying subpathway signatures for individualized anticancer drug response by integrating multi-omics data](https://doi.org/10.1186/s12967-019-2010-4)
- 2019 - Ghaemi - [Multiomics modeling of the immunome, transcriptome, microbiome, proteome and metabolome adaptations during human pregnancy](https://doi.org/10.1093/bioinformatics/bty537) Multi-omics in pregnancy using stacked generalization

## Multi-omics data management

- 2017 - [MultiAssayExperiment](https://bioconductor.org/packages/MultiAssayExperiment/) - Ramos - Software for the integration of multi-omics experiments in Bioconductor - [paper](https://doi.org/10.1158/0008-5472.CAN-17-0344).
- 2021 - [muon](https://github.com/pmbio/muon) - Bredikhin - [Multimodal omics analysis framework](https://doi.org/10.1101/2021.06.01.445670)

## Batch effect correction

- 2020 - [MultiBaC](https://github.com/ConesaLab/MultiBaC) - Ugidos - [MultiBaC: A strategy to remove batch effects between different omic data types](https://doi.org/10.1177/0962280220907365) and [R package publication](https://doi.org/10.1093/bioinformatics/btac132)
- [A multivariate method to correct for batch effects in microbiome data]()

## Meetings and workshops

- 2020 - [Mathematical Frameworks for Integrative Analysis of Emerging Biological Data Types](https://www.birs.ca/events/2020/5-day-workshops/20w5197) - [Hackathon details](https://github.com/BIRSBiointegration/Hackathon) - June 14-19, 2020 in Banff, Canada
- 
