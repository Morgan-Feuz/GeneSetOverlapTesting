# Gene Set Overlap Significance Testing of RNA-seq Data in R

## Background
Two commonly used methods in bioinformatics to determine whether differential gene set overlap is due to random chance or not include the [Hypergeometric disribution](https://en.wikipedia.org/wiki/Hypergeometric_distribution) or [Fisher's exact tests](https://en.wikipedia.org/wiki/Fisher%27s_exact_test). The R software base provides functions such as [phyper](http://stat.ethz.ch/R-manual/R-devel/library/stats/html/Hypergeometric.html) and [fisher.test](http://stat.ethz.ch/R-manual/R-patched/library/stats/html/fisher.test.html), while another good option is the [GeneOverlap](https://bioconductor.org/packages/release/bioc/html/GeneOverlap.html) package available through Bioconductor. 

The classical example used for explaining the hypergeometric distribution is the idea of the probability of sampling white or black balls from an urn without replacement (i.e., [The Urn Model](https://learningds.org/ch/03/theory_urn.html)). 

----------------------------------------------------------------
## Usage 

#### Method 1: Using Base R 
```r
# Define Hypergeometric Distribution parameters
help("phyper")

# Over-representation
p.value <- phyper(overlap_genes-1, group2_DEGs, background-group2_DEGs, group1_DEGs, lower.tail = TRUE, log.p = FALSE)

# or
p.value <- sum(dhyper(overlap_genes:group1_DEGs, group2_DEGs, background-group2_DEGs, group1_DEGs))

# or
fisher.test(matrix(c(overlap_genes, group2_DEGs-overlap_genes, group1_DEGs-overlap_genes, background-group2_DEGs-group1_DEGs+overlap_genes), 2, 2), alternative='greater')$p.value
```
<ins>Note:</ins> Here, essentially we are testing the probablility of sampling the number of `overlap_genes` (or more) in a sample size of `group1_DEGs` from an urn containing `group2_DEGs` and `background-group2_DEGs`. 
+ The `background` or `universe` genes represent the pool of genes used for differential testing analysis.
+ See the additional references section at the bottom of the page for more information, as well as for testing for under-representation using `phyper`. 

<br>

#### Method 2: Using the GeneOverlap Package
To install this package, start R (version "4.4") and enter:
```r
# Installation
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager") BiocManager::install("GeneOverlap")

browseVignettes("GeneOverlap")

# Usage
# Create a GeneOverlap object. Groups = lists and genome.size = int
go.obj <- newGeneOverlap(group1_DEGs$Symbol,
                         group2_DEGs$Symbol,
                         genome.size = universe)

go.obj <- testGeneOverlap(go.obj)
print(go.obj)

# View the results 
head(getIntersection(go.obj))
head(getOddsRatio(go.obj))
head(getContbl(go.obj))
head(getGenomeSize(go.obj))

```

<ins>Note:</ins> Here, the `GeneOverlap` class tests whether the input variables are independent, which can be represented as a contingency table, and then uses Fisher's exact test to find the statistical significance. The results include: 
- The p-value
- Odds ratio, which represents the strength of the association. If an odds ratio is equal to or less than 1, there is no association between the two lists. If the odds ratio is much larger than 1, then the association is strong. 
- Jaccard index, which measures the similarity between the two lists. The Jaccard index varies between 0 and 1, with 0 indicating there is no similarity between the two and 1 meaning that the two lists are identical.

-------------------------------------------------------------------------

#### Additional References
- https://mengnote.blogspot.com/2012/12/calculate-correct-hypergeometric-p.html
- https://stats.stackexchange.com/questions/16247/calculating-the-probability-of-gene-list-overlap-between-an-rna-seq-and-a-chip-c/16259#16259
- https://www.biostars.org/p/485827/#9483835
- https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Hypergeometric.html
- https://bioconductor.org/packages/release/bioc/html/GeneOverlap.html
- https://seqqc.wordpress.com/2019/07/25/how-to-use-phyper-in-r/
- https://www.geeksforgeeks.org/hypergeometric-distribution-in-r-programming/
- https://dputhier.github.io/ASG/practicals/go_statistics_td/go_statistics_td_2015.html
