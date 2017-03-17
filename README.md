# The 3D genome organization of *Drosophila melanogaster* through data integration
Genome structures are dynamic and non-randomly organized in the nucleus of higher eukaryotes. To maximize the accuracy and coverage of 3D genome structural models, it is important to integrate all available sources of experimental information about a genome’s organization. It remains a major challenge to integrate such data from various complementary experimental methods. Here, we present an approach for data integration to determine a population of complete 3D genome structures that are statistically consistent with data from both genome-wide chromosome conformation capture (Hi-C) and lamina-DamID experiments. Our structures resolve the genome at the resolution of topological domains, and reproduce simultaneously both sets of experimental data. Importantly, this framework allows for structural heterogeneity between cells, and hence accounts for the expected plasticity of genome structures. As a case study we choose Drosophila melanogaster embryonic cells, for which both data types are available. Our 3D geome structures have strong predictive power for structural features not directly visible in the initial data sets, and reproduce experimental hallmarks of the *D. melanogaster* genome organization from independent and our own imaging experiments. Also they reveal a number of new insights about the genome organization and its functional relevance, including the preferred locations of heterochromatic satellites of different chromosomes, and observations about homologous pairing that cannot be directly observed in the original Hi-C or lamina-DamID data. Our approach allows systematic integration of Hi-C and lamina-DamID data for complete 3D genome structure calculation, while also explicitly considering genome structural variability.

Here is the overview of the population-based genome structure modeling approach and its application to the Drosophila genome
<p align="center">
  <img src="https://github.com/QingjiaoLi/3DDM/blob/master/Fig1.png" width="600" />
</p>

A preprint can be found on [bioRxiv](http://biorxiv.org/content/early/2017/01/15/099911)
Li, Q., et al. (2017). "The 3D genome organization of Drosophila melanogaster through data integration." bioRxiv.
---

This is a collection of modeling codes in python and pbs scripts for submitting jobs. We highly recommend users to run on high performance computing environment (HPC) because of computational resource and running time. 

### Installation
Requirements:

- Python 2.7
- Python packages ``numpy``, ``scipy``
- IMP version 1 ([Integrative Modeling Package](https://integrativemodeling.org/1.0/download/))

How to run:
- Download python codes and shell scripts
- Download the input files from the folder input
- Follow the steps in *instruction.sh*

---
### 
The population of structures (pdb files) generated in the manuscript can be dowonloaded from our webserver ()
