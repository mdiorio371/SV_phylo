[![Codecov test coverage](https://codecov.io/gh/your-username/SV_phylo/branch/main/graph/badge.svg)](https://app.codecov.io/gh/your-username/SV_phylo)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
![Status](https://img.shields.io/badge/status-work--in--progress-orange)


**SV_phylo** is an R package for constructing phylogenies with **structural variant (SV) data** in bacterial genomes.  
It builds on the [`SVMC`](https://github.com/mdiorio371/SVMC) package, which performs SV mapping and characterization, by providing methods to:

- Cluster by SV recurrence 
- Predict systematically impactful SVs
- Compare SV-based phylogenies with SNP/reference trees  

---

## Installation

Both **SV_phylo** and [**SVMC**](https://github.com/mdiorio371/SVMC) are needed: 

```r
# install devtools if not already installed
install.packages("devtools")

# install SVMC first
devtools::install_github("mdiorio371/SVMC")

# then install SV_phylo
devtools::install_github("your-username/SV_phylo")


library(SVMC)
library(SV_phylo)

# 1. Generate structural variant calls with SVMC
sv_table <- SVMC("Salmonella_enterica", n = 200, method = "One-vs-all")

# 2. Cluster genomes by large structural variants
SV_clusters <- sv_cluster(sv_table, minLength = 5e4)

# 3. Subsample clusters and generate parSNP phylogeny
cluster_subsample <- subsample_sv_clusters(SV_clusters, n = 3)
snp_clusters <- snp_phylo(cluster_subsample)

#4. Plot phylo tree with annotated SVs
SV_SNP_tree <- annotate_tree(SV_clusters, snp_clusters)

plot(SV_SNP_tree, main = "SV phylogeny")


#5. 

SV_clusters <- (dist_mat, method = "nj")


# 5. Plot the tree
plot(tree, main = "SV-based phylogeny")

```
