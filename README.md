# awesome-multi-omics

List of software packages for multi-omics data analysis. A
shameless copy of Sean Davis'
[awesome-single-cell](https://github.com/seandavi/awesome-single-cell)
repo.

[Contributions welcome](https://github.com/mikelove/awesome-multi-omics/blob/master/CONTRIBUTING.md)...

For brevity, below lists only the first author of multi-omics methods.

## Software packages and methods

### Multi-omics correlation analysis

- 2007 - **SCCA** - Parkhomenko - sparse CCA [paper 1](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2367499/), [paper 2](https://doi.org/10.2202/1544-6115.1406) 
- 2008 - **PCCA** - Waaijenborg - penalized CCA / "CCA-EN" [paper](https://doi.org/10.2202/1544-6115.1329)
- 2009 - [Sparse mCCA](https://CRAN.r-project.org/package=PMA) - Witten - [paper 1](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2697346/), [paper 2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2861323/)
- 2009 - **sPLS** - Lê Cao - sparse PLS [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2640358/)
- 2010 - **Regularized dual CCA** - Soneson - [paper](https://doi.org/10.1186/1471-2105-11-191)
- 2013 - [MFA](https://cran.r-project.org/package=FactoMineR) - Abdi - multiple factor analysis - [paper](https://doi.org/10.1002/wics.1246)
- 2013 - [JIVE](https://genome.unc.edu/jive/) - Lock - joint & individual variance explained - [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3671601/)
- 2014 - [MCIA](https://bioconductor.org/packages/omicade4) - Meng - multiple co-interia analysis - [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4053266/)
- 2014 - [RGCCA](https://cran.r-project.org/package=RGCCA) - Tenenhaus - Regularized generalized CCA - [paper](https://www.ncbi.nlm.nih.gov/pubmed/24550197)
- 2014 - [STATegRa](https://bioconductor.org/packages/STATegRa) - Gomez-Cabrero - DISCO, JIVE, & O2PLS (several papers)
- 2016 - **Bayesian group factor analysis** - [paper](https://arxiv.org/abs/1411.2698)
- 2016 - [MSFA](https://github.com/rdevito/MSFA) - De Vito - multi-study factor analysis: same features, different samples - [paper](https://arxiv.org/abs/1611.06350)
- 2016 - [moGSA](https://bioconductor.org/packages/mogsa) - Meng - multi-omics gene set analysis - [paper](https://www.biorxiv.org/content/10.1101/046904v2)
- 2017 - [mixOmics](https://bioconductor.org/packages/mixOmics) - Rohart - [paper](https://doi.org/10.1371/journal.pcbi.1005752)
- 2018 - [AJIVE](https://github.com/idc9/r_jive) - Feng - angle-based JIVE - [paper](https://arxiv.org/abs/1704.02060)
- 2018 - [MOFA](https://github.com/bioFAM/MOFA) - Argelaguet - multi-omics factor analysis - [paper](http://msb.embopress.org/content/14/6/e8124), [application](https://www.biorxiv.org/content/10.1101/519207v1)
- 2018 - [PCA+CCA](https://github.com/pachterlab/PCACCA/) - Brown - [paper](https://www.biorxiv.org/content/early/2018/07/09/364448)
- ...

### Multi-table methods from ecology

- 1994 - **COI** - Doledec - Co‐inertia analysis - [paper](https://doi.org/10.1111/j.1365-2427.1994.tb01741.x)
- 2007 - [ade4](https://CRAN.r-project.org/package=ade4) - Dray - Implementing the Duality Diagram for Ecologists - [paper](http://dx.doi.org/10.18637/jss.v022.i04)

### Multi-omics classification

- 2009 - [iCluster](https://cran.r-project.org/package=iCluster) - Shen - [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2800366/)
- 2013 - [iCluster+](https://bioconductor.org/packages/iClusterPlus) - Mo - [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3600490/)
- 2014 - [SNF](http://compbio.cs.toronto.edu/SNF/SNF/Software.html) - Wang - [paper](https://www.ncbi.nlm.nih.gov/pubmed/24464287)
- ...

### Integrating multi-omics measured in separate samples

- 2018 - [cardelino](https://github.com/PMBio/cardelino) - assign gene expression states to clones by integrating SNVs from scRNA-seq with bulk exome data
- 2018 - [clonealign](https://github.com/kieranrcampbell/clonealign) - assign gene expression states to clones by integrating scRNA-seq with scDNA-seq (CNV) data - [paper](https://www.biorxiv.org/content/early/2018/06/11/344309)

## Multi-omics reviews / evaluations

- 2014 - [A practical data processing workflow for multi-OMICS projects](https://doi.org/10.1016/j.bbapap.2013.02.029)
- 2016 - [Multi-omic data integration enables discovery of hidden biological regularities](https://www.nature.com/articles/ncomms13091)
- 2016 - [Dimension reduction techniques for the integrative analysis of multi-omics data](https://doi.org/10.1093/bib/bbv108) - from *omicade4* and *moGSA* authors
- 2017 - [More Is Better: Recent Progress in Multi-Omics Data Integration Methods](https://doi.org/10.3389/fgene.2017.00084)
- 2017 - [Multi-omics approaches to disease](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1215-1)
- 2018 - [MOVIE: Multi-Omics VIsualization of Estimated contributions](https://www.biorxiv.org/content/early/2018/07/29/379115) - [code](https://github.com/mccabes292/movie)
- 2018 - [Multi-omic and multi-view clustering algorithms: review and cancer benchmark](https://doi.org/10.1093/nar/gky889)
