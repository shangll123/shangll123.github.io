---
layout: paper
title: "Spatially aware dimension reduction method in spatial transcriptomics"
image: /assets/images/papers/SpatialPCA_workflow.png
authors: Lulu Shang, Xiang Zhou
year: 2022
shortref: Shang et al. (2022) 
journal: 
type: statistical
doi: 
---

# Abstract

**Overview**
SpatialPCA is a spatially aware dimension reduction method that aims to infer a low dimensional representation of the gene expression data in spatial transcriptomics. 

**Materials and Methods**
SpatialPCA builds upon the probabilistic version of PCA, incorporates localization information as additional input, and uses a kernel matrix to explicitly model the spatial correlation structure across tissue locations.

**Results**
In the real data applications, SpatialPCA identifies key molecular and immunological signatures in a newly detected tumor surrounding microenvironment, including a tertiary lymphoid structure that shapes the gradual transcriptomic transition during tumorigenesis and metastasis. In addition, SpatialPCA detects the past neuronal developmental history that underlies the current transcriptomic landscape across tissue locations in the cortex.

**Conclusion**
SpatialPCA is computationally efficient, easily scalable to spatial transcriptomics with tens of thousands of spatial locations and thousands of genes. We have illustrated the benefits of SpatialPCA for spatial transcriptomics visualization, spatial domain detection, trajectory inference on the tissue, as well as high-resolution spatial map construction.

<div class="middle">
    <img src="/assets/images/papers/SpatialPCA_workflow.png" alt="photo" width="600"/>
</div>

**Biorxiv**
[https://www.biorxiv.org/content/10.1101/2022.01.19.476966v1.full.pdf](https://www.biorxiv.org/content/10.1101/2022.01.19.476966v1.full.pdf)


**Tutorial of SpatialPCA**
[SpatialPCA Tutorial](https://lulushang.org/SpatialPCA_Tutorial/index.html)

**Analysis Codes**
<br />
[Simulation](https://lulushang.org/docs/Projects/SpatialPCA/Simulation)
<br />
[DLPFC dataset](https://lulushang.org/docs/Projects/SpatialPCA/DLPFC)
<br />
[Slide-Seq mouse cerebellum dataset](https://lulushang.org/docs/Projects/SpatialPCA/Slideseq)
<br />
[Slide-Seq V2 mouse hippocampus dataset](https://lulushang.org/docs/Projects/SpatialPCA/SlideseqV2)
<br />
[Human breast cancer ST dataset](https://lulushang.org/docs/Projects/SpatialPCA/HER2ST)
<br />


