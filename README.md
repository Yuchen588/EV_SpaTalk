# EV_SpaTalk
[![R >4.0](https://img.shields.io/badge/R-%3E%3D4.0-brightgreen)](https://www.r-project.org/) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6809147.svg)](https://doi.org/10.5281/zenodo.6809147) [![jupyter](https://img.shields.io/badge/Jupyter--notebook-SpaTalk--tutorial-yellow?logo=jupyter)](https://github.com/multitalk/awesome-cell-cell-communication/blob/main/method/SpaTalk.ipynb) [![DOI](https://img.shields.io/badge/DOI-10.1038%2Fs41467--022--32111--8-yellowgreen)](https://www.nature.com/articles/s41467-022-32111-8)

### An EV-mediated cell-cell interaction inference framework from spatially resolved transcriptomic data

<img src='https://github.com/ZJUFanLab/SpaTalk/blob/main/img/SpaTalk.png'>

[Extracellular vesicles (EVs)](https://pubmed.ncbi.nlm.nih.gov/30914410/), produced by all eukaryotic cells, are essential for intercellular communication. Beyond transporting internal cargo, EVs transfer specific ligands from donor to receptors on recipient cells, playing a crucial role in orchestrating intercellular dynamics and modulating the biological landscape of human ecosystem. However, inferring EV-driven tissue-cellular interactome presents a challenge. To address this, we introduce EV_SpaTalk, an innovative framework designed to infer the significance of [ligand-receptor pairs](https://pubmed.ncbi.nlm.nih.gov/34331449/) in EV-mediated cell-cell crosstalk networks within tumor microenvironments. This framework relies on EV budding power, a key cellular state in [ESCRT+ cells](https://pubmed.ncbi.nlm.nih.gov/31705132/), and the spatial interoperable extrapolation of EV-specific ligand-receptor pairs. EV_SpaTalk models the ligand-receptor signaling network between spatially adjacent sender-receiver niches, utilizing a non-negative linear model and spatial distance mapping to decompose cell components from [spatially resolved transcriptomics (ST)](https://pubmed.ncbi.nlm.nih.gov/31932730/) and paired single-cell RNA-seq data. This framework efficiently deciphers EV-driven ligand-receptor interactions across various cell types, offering an invaluable tool for researchers to investigate intricate EV-mediated spatial cell dynamics in tissue, thereby significantly enhancing our understanding of EVs' roles in cellular communication and tumor biology.

## EV_SpaTalkdb
The EV_SpaTalkdb is a curated database containing 2,997 EV-specific ligand-receptor interaction (LRI) pairs, derived from [cellchatdb] (http://www.cellchat.org/cellchatdb/) and classified into three categories: Cell-Cell Contact, ECM-Receptor, and Secreted Signaling. Validation of these pairs is backed by extensive EV MASS data from [EVpedia] (https://evpedia.info/evpedia2_xe/) and [EV-related studies] (https://pubmed.ncbi.nlm.nih.gov/35918900/; https://pubmed.ncbi.nlm.nih.gov/34817906/), encompassing global EV and EV surface membrane proteomics covers both small (exosomes) and large (e.g. macrovesicles) EVs. Researchers can access the database via the 'EV_spatalk_object@database' command in R, post-loading the EV_spatalk_object, an S4 object that integrates results from EV_SpaTalk workflow. We encourage the scientific community to enhance the EV_SpaTalkdb by contributing verified EV-related LR pairs, thus promoting a cooperative environment for EV research advancements.



# Install

- install dependent packages `devtools` and [`NNLM`](https://github.com/linxihui/NNLM)

```
> install.packages(pkgs = 'devtools')
#Please make sure the following packages were installed (Seurat, tidyverse, future, nnls, parallel, Matrix, parallel, Matrix, VennDiagram, ggplot2, rstatix, dplyr, tidyr, purrr, broom, RobustRankAggreg, gridExtra, circlize)
```

- then install EV_SpaTalk

```
> devtools::install_github("Yuchen588/EV_SpaTalk@master", force = T)

# or download the repository as ZIP
> devtools::install_local("/path/to/EVSpaTalk.tar.gz")
```

# Usage and steps:
EV_spatalk method consists of two components, wherein the first is to use the scRNA-seq profile dissect the cell-type composition of ST data and the second is to infer the spatially resolved EV-mediated cell-cell communications over the decomposed single-cell ST data. Classification and description of EV_spatalk functions are shown in the [tutorial](https://evpedia.info/evpedia2_xe/).

- ### Cell-type deconvolution to reconstruct single-cell ST atlas with known cell types from scRNA-seq data
```
# st_data: A standard Seurat object of ST data (finished raw clustering and dimension reduction)
# sc_data: A standard Seurat object of scRNA-seq data (preferably paired data for ST)
# select.celltype: A character containing the cell types for scRNA-seq data
# sc.anno.id: A character containing the column name corresponding to cell annotations

> EV_spatalk_object@st.seurat.obj <- st_data
> EV_spatalk_object@sc.seurat.obj <- sc_data
> EV.spatalk.results <- st_deco_anno(st = st_data, sc = sc_data, nrand = 2, EV_spatalk_object=EV_spatalk_object, nbin = 5, pval_thresh=0.05, mc.cores=30, sc.anno.id="pop", set.seeds = 1, select.celltype = select.celltype)
```

- ### Screening for spatially approximate sender-receiver niche and ligand-receptor interaction densities
```
# s.cell.type: the vector of celltype ids 
# comm_list= c("Cell-Cell Contact", "ECM-Receptor", "Secreted Signaling")

> EV.spatalk.results <- find_niche_LR(EV_spatalk_object=EV.spatalk.results, prox="inverse", mc.cores=30, s.cell.type = c("Malignant", "T_cells"), comm_list=comm_list, datatype='mean count', method = "pseudocount")
```
- ### Compressing LR pairs that occur with high frequency in each spot/cell
```
> EV.spatalk.results <- find.inter.LR(EV_spatalk_object=EV.spatalk.results, mc.cores=30)
```
- ### Inference of EV-mediate cell-cell communication and ligand-receptor network in space
```
> EV.spatalk.results <- find_EV_spatalk_LR(EV_spatalk_object=EV.spatalk.results, mc.cores=30, seeds = 1)
```
- ### Scaling and integrating the EV-mediated LR interactome into Seurat object
```
> EV.spatalk.results <- add_interaction_score(EV_spatalk_object = EV.spatalk.results, mc.cores=30)
```
- ### Visualization-related features (see our [tutorial](https://evpedia.info/evpedia2_xe/) page for details)
```
> Spatial intensity plot for each EV-mediated LR pair
> LR_spatial_indensity_plot(EV_spatalk_object=EV.spatalk.results, s.LR.pair="CD86_CTLA4")
```
> The pie plot to illustrate the significate of EV-related LR
> LR_pie.plot(EV_spatalk_object=EV.spatalk.results)
```
> Venn diagram summarizing LRs correlated with EV release and spatial distance
> LR_venn.plot(EV_spatalk_object=EV.spatalk.results)
```
> Circos plot representing statistical results of candidate LRs and their enrichment correlation in the EV interactome
> select.LR.id <- c(EV.spatalk.results@inter.LR.results$common_positive_LR, EV.spatalk.results@inter.LR.results$common_negative_LR)
> LR_circos.plot(EV_spatalk_object = EV.spatalk.results, select.LR.id=select.LR.id, specific.LR.id="CD96_PVR")
> circos.clear()
```

# Note
[![CellTalkDB v1.0](https://img.shields.io/badge/CellTalkDB-v1.0-blueviolet)](http://tcm.zju.edu.cn/celltalkdb/) [![KEGG pathway](https://img.shields.io/badge/KEGG-pathway-ff69b4)](https://www.kegg.jp/kegg/pathway.html) [![Reactome pathway](https://img.shields.io/badge/Reactome-pathway-brightgreen)](https://reactome.org/) [![AnimalTFDB v3.0](https://img.shields.io/badge/AnimalTFDB-v3.0-yellowgreen)](http://bioinfo.life.hust.edu.cn/AnimalTFDB/#!/)

- __EV_SpaTalk can be applied to either single-cell (vignette) or spot-based (vignette) ST data__
- __EV_SpaTalk allows to use custom LRIs(wiki), pathways, and TFs database (wiki)__
- __EV_SpaTalk allows to use the parallel processing for st_deco_anno(), find_niche_LR(), find.inter.LR(), find_EV_spatalk_LR(), add_interaction_score()__
- __EV_SpaTalk allows to use other deconvolution methods followed by the inference of cell-cell communications__
  - RCTD, Seurat, SPOTlight, deconvSeq, stereoscope, cell2location, or other methods
- __EV_SpaTalk allows to directly infer cell-cell communications skiping deconvolution (e.g., single-cell based ST data)__
- __EV_SpaTalk can spatially visualize cell-type compositions/distributions (wiki) and cell-cell communications (wiki)__
- LRI and pathways can be download at/[`database/`](https://github.com/ZJUFanLab/SpaTalk/tree/main/data)
- Demo data can be download at /[`demo_data/`](https://github.com/ZJUFanLab/SpaTalk/tree/main/inst/extdata)

__Please refer to the [tutorial vignette](https://raw.githack.com/multitalk/awesome-cell-cell-communication/main/method/tutorial.html) with demo data processing steps. Detailed functions see the [document](https://raw.githack.com/ZJUFanLab/SpaTalk/main/vignettes/SpaTalk.pdf)__

# About
EV_SpaTalk was developed by Yuchen Li. Should you have any questions, please contact Yuchen at ycli16@stanford.edu.

Please cite us as "XXX"
