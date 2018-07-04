# Tutorial: *ABBA* *BABA* statistics
### Simon Martin (shm45@cam.ac.uk)

___
## Requirements
Python 2.7
Numpy 1.10+
R 3.0+

___
## Introduction

ABBA BABA statistics (also called D statistics) provide a simple and powerful test for a deviation from a strict bifurcating evolutionary history. They are therefore frequently used to test for introgression using genome-scale SNP data.

In this practical we will perform an ABBA BABA analysis using a **combination of available software and some code written from scratch in R**. We will analyse genomic data from several populations of *Heliconius* butterflies.

#### Workflow
Starting with genotype data from multiple individuals, we first **infer allele frequencies** at each SNP. We then **compute the *D* statistic** and then use a **block jackknife** method to test for a significant deviation from the null expectation of *D*=0. Finally we **estimate *f* the 'admixture proportion'**.

#### The hypotheses
We will study multiple races from three species: *Heliconius melpomene*, *Heliconius timareta* and *Heliconius cydno*. These species have partially overlapping ranges and they are thought to hybridise where they occur in sympatry. Our sample set includes two pairs of sympatric races of *H. melpomene* and *H. cydno* from Panama and the western slopes of the Andes in Colombia. There are also two pairs of sympatric races of *H. melpomene* and H. timareta from the eastern slopes of the Andes in Colombia and Peru. We hypothesize that hybridisation between species in sympatry will lead to sharing of genetic variation between H. cydno and the sympatric races of H. melpomene from the west, and between H. timareta and the corresponding sympatric races of *H. melpomene* from the east. This should include sharing of wing pattesourcerning alleles between H. melpomene and H. timareta, which are co-mimetics. There is also another race of *H. melpomene* from French Guiana that is allopatric from both *H. timareta* and *H. cydno*, which should have not experienced recent genetic exchange with either species and therefore serves as a control. Finally, there are two samples from an outgroup species *Heliconius numata silvana*.

![Species Map](images/hel92_pops.png)


## Preparation

Open a terminal window and navigate to a folder where you will run the excersise and store all the input and output data files.

This tutorial makes use of a collection of python scripts that must be downloaded from [GitHub](https://github.com/simonhmartin)

```bash
git clone https://github.com/simonhmartin/genomics_general
```

## A genome wide test for introgression

In its simplest formulation, the *ABBA* *BABA* test relies on counts of sites in the genome that match the *ABBA* and *BABA* genotype patterns. That is, given three ingroup populations and an outgroup with the relationship (((P1,P2),P3),O), and given a single genome sequence representing each population (ie, H1, H2 and H3), ***ABBA*** sites are those at which H2 and H3 **share a derived allele ('B')**, while **H1 has the ancestral state ('A')**, as defined by the outgroup sample. Likewise, ***BABA*** represents sites at which **H1 and H3 share the derived allele**.

Ignoring recurrant mutation, the two SNP patterns can only be produced if some parts of the genome have genealogies that do not follow the 'species tree', but instead group H2 with H3 or H1 with H3. If the populations split fairly recently, such 'discordant' genealogies are expected to occur in some parts of the genome due to variation in lineage sorting. In the absence of any deviation from a strict bifurcating topology, **we expect roughly equal proportions of the genome to show the two discordant genealogies** (((H2,H3),H1),O) and (((H1,H3),H2),O). By counting *ABBA* and *BABA* SNPs across the genome (or a large proportion of it), we are therefore **approximating the proportion of the genome represented by the two discordant genealogies**, which means **we expect a 1:1 ratio of *ABBA* and *BABA* SNPs**. A deviation could come about as a result of gene flow between populations P3 and P2 for example, although it could also indicate other phenomena that break our assumptions, such as ancestral population structure, or variable substitution rates.

To quantify the deviation from the expected ratio, we calculate *D*, which is the difference in the sum of *ABBA* and *BABA* patterns across the genome, divided by their sum:

D = \[sum(ABBA) - sub(BABA)\] / \[sum(ABBA) + sub(BABA)\]

**Therefore, D ranges from -1 to 1, and should equal 0 under the null hypothesis. D > 1 indicates and excess of *ABBA*, and D < 1 indicates an excess of *BABA*.**

If we have multiple samples from each population, then counting *ABBA* and *BABA* sites is less straghtforward. One option is to consider only sites at which all samples from the same population share the same allele, but that will discard a large amount of useful data. A preferable option is to use the allele frequencies at each site to quantify the extent to which the genealogy is skewed toward the *ABBA* or *BABA* pattern. This is effectively equivalent to counting *ABBA* and *BABA* SNPs using all possible sets of four haploid genomes at each site. *ABBA* and *BABA* are therefore no longer binary states, but rather numbers between 0 and 1 that represent the frequency of allele combinations matching each genealogy. They are computed based on the frequency of the derived allele (*p*) and ancestral allele (1-*p*) in each population as follows:

*ABBA* = (1-*p1*) x *p2* x *p3* x 1-*pO*
*BABA* = *p1* x (1-*p2*) x *p3* x 1-*pO*

#### Genome wide allele frequencies

To compute these values from population genomic data, we need to first determine the frequency of the derived allele in each populaton at each polymorphic site in the genome. We will compute these from the *Heliconius* genotype data provided using a python script. The input file has already been filtered to contain only bi-allelic sites. The frequencies script requires that we define populations. These are defined in the file `hel92.pop.txt`.

```bash
python genomics_general/freq.py -g data/hel92.DP8MP4BIMAC2HET75dist1K.geno.gz -p mpg -p ros -p vul -p mal -p ama -p chi -p zel -p flo -p txn -p slv --popsFile data/hel92.pop.txt --target derived -o data/hel92.DP8MP4BIMAC2HET75dist1K.derFreq.tsv.gz
```
By setting `--target derived` we obtain the frquency of the derived allele in each population at each site. This is based on using the final population specified (*H. numata silvana*, or '*slv*') as the outgroup. Sites at which this population is not fixed for the ancestral state are discarded.

#### Genome wide ABBA BABA analysis

**(NOTE: here we're working in R, or R Studio if you prefer)**

To learn how the ABBA BABA test works, we will writing the code from scratch to do the test. Therefore, **start a new R script**.

First we define functions for computing the ABBA and BABA proportions at each site. These will take as input the frequency of the derived allele in populations P1, P2 and P3 (i.e. *p1*, *p1* and *p3*). (The frequency of the ancestral allele in the outgroup will be 1 at all sites, by definition, so this can be ignored).

```R
abba <- function(p1, p2, p3) (1 - p1) * p2 * p3

baba <- function(p1, p2, p3) p1 * (1 - p2) * p3
```

These functions will each return a single value if `p1`, `p2` and `p3` are single vlaues, but they will return a vector of values for each site if `p1`, `p2` and `p3` are each vectors.

We then define a function for *D*, which takes the vectors of `ABBA` and `BABA` values for all sites as inputs.


```R
D.stat <- function(ABBA, BABA) (sum(ABBA) - sum(BABA)) / 
                               (sum(ABBA) + sum(BABA))
```


Read in our allele frequency data.

```R
freq_table = read.table("data/hel92.DP8MP4BIMAC2HET75dist1K.derFreq.tsv.gz", header=T, as.is=T)
```

This has created an object called `freq_table` that contains the frequencies for the derived allele at each SNP.

We can check the number of sites in this table, and also look at the first few rows to get a feel for the data.

```
nrow(freq_table)

head(freq_table)
```

Note that the first two columns give the name of the scaffold (i.e. the chromosome) and the position on the chromosome of each site. The remaining columns are the allele frequencies for the different subspecies, as indicated in the figure above.

Now, to compute D, we need to define populations P1, P2 and P3. We will start with an obvious and **previously published test case**:
We will ask whether there is evidence of introgression between ***H. melpomene rosina* (*ros*)** and ***H. cydno chioneus* (*chi*)**. These will be **P2** and **P3** respectively. **P1** will be our **allopatirc** population ***H. melpomene melpomene* from French Guiana (*mpg*)**.

We set these populations and then compute *D* by extracting the the derived allele frequencies for all SNPs for the three populations.

```R
P1 <- "mel_mel"
P2 <- "mel_ros"
P3 <- "cyd_chi"

ABBA <- abba(freq_table[,P1], freq_table[,P2], freq_table[,P3])
BABA <- baba(freq_table[,P1], freq_table[,P2], freq_table[,P3])

D = D.stat(ABBA,BABA)

D
```

We get a **strongly positive D statistic** (remember D varies from -1 to 1), indicating an excess of ABBA over BABA. This indicates that ***H. cydno chioneus* from Panama** (`cyd_chi`) **shares more genetic variation with the sympatric *H. melpomene rosina* from Panama** (`mel_ros`) than with the allopatirc *H. melpomene melpomene* from French Guiana (`mel_mel`). This is consistent with hybridisation and gene flow between the two species where they occur in sympatry.

However, we currently don't know whether this result is statistically robust. In particular, we don't know whether the excess of ABBA is evenly distributed across the genome. If it results from odd ancestry at just one part of teh genome, we would have less confidence that there has been significant intogression.

To test for a consistent genome-wide signal we use a block-jackknife procedure. 

#### Block Jackknife

The Jackknife procedure allows us to compute the variance of *D* despite non-independence among sites. A more conventional bootstrapping approach, where we would randomly resample sites and recalculate *D*, is not appropriate, because nearby sites in the genome have similar ancestry, making them non-independnent observations.

The block jackknife procedure estimates the standard deviation for so-called 'pseudovalues' of the mean genome-wide *D*, where each pseudovalue is computed by excluding a defined block of the genome, taking the difference between the mean genom-wide *D* and *D* computed when the block is omitted. The block size needs to exceed the distance at which autocorrelation occurs. In our case, we will use a block size of 1 Mb. In fact we know that linkage disequilibrium decays to background levels at a distance of well below 1 Mb, so this is a conservative block size.

To define the start and end positions of each block, we need to know the lengths of each chromosome. These are provided in the file `"data/Hmel2_chrom_lengths.txt"`. We load this file and then make a vector of chromosome lengths with the chromosome names ("chr1", chr2" etc. as the names of this vector).

```R
chrom_table <- read.table("data/Hmel2_chrom_lengths.txt")
chrom_lengths <- chrom_table[,2]
names(chrom_lengths) <- chrom_table[,1]

chrom_lengths
```

The code to run the jackknife procedure is fairly simple, but we are not going to write it here, because it requires more advanced R knowledge ('for loops' or the 'apply' methods). Instead, we some R functions for this porpose are provided in a separate script, which we can import now.

```R
source("genomics_general/jackknife.R")
```

The first step in the process is to define the blocks that will be omitted from the genome in each iteration of the jackknife. the function `get_genome_blocks` in the jackknife script will do this. It requires that we specify the block size (here we use 1 Mb) and the lengths of chromosomes.

```R
blocks = get_genome_blocks(block_size=1e6, chrom_lengths=chrom_lengths)
n_blocks = nrow(blocks)
```
This gives a table with a row for each block, giving its start and end position and the chromosome it is on.

```R
head(blocks)
nrow(blocks)
```

In each jackknife iteration, we will recompute the D statistic using all sites in the genome that are not in a given block. We therefore make a list of sites that we will retain in each iteration: those sites that do not fall within each block.

```R
indices <- get_genome_jackknife_indices(chromosome=freq_table$scaffold,
                                        position=freq_table$position,
                                        block_info=blocks)
```

Now we can run the block jackknifing procedure. We provide the the D statistic function we created earlier, which it will use each iteration. We also provide the input ABBA BABA dataframe, and the block indices.

```R
D_sd <- get_jackknife_sd(jackknife_indices=indices, FUN=D.stat, ABBA,BABA)
```


This provides an unbiased estimate of the standard deviation of *D*. From this we can estimate the standard error, Z score and p-value for the test of whether *D* deviates significantly from zero.

```R
D_err <- D_sd/sqrt(n_blocks)
D_Z <- D / D_err
D_p <- 2*pnorm(-abs(D_Z))
```

In this case the deviation is hugely significant.

#### In your own time
What do we find if we change the identity of P1, P2 and P3? For example, we can test for recent admixture by asking whether *H. melpomene rosina* shares more variation with *H. cydno chioneus* than does H. melpomene vulcanus (a close relative of *H. m. rosina* from nearby).


#### Estimating the admixture proportion

The *D* statistic provides a powerful test for introgression, but it does not quantify the proportion of the genome that has been shared. A related method has been developed to estimate *f*, admixture proportion.

The idea behind this approach is that we compare the observed excess of ABBA over BABA sites, to that which would be expected under complete admixture. To approximate the expectation under complete admixture we re-count ABBA and BABA but substituting a second population of the P3 species in the place of P2. If you lack a second population, you can simply split your P3 samples into two. In this case, we have two populations to represent each species, so if we're using H. cydno chioneus (chi) as P3a, we can use H. cydno zelinde (zel) as P3b).

```R
P3a <- "cyd_chi"
P3b <- "cyd_zel"
```

We then compute ABBA and BABA using P1, P2 and P3a, as above, as well as using P1, P3b and P3a.

```R
ABBA_1_2_3a <- abba(freq_table[,P1], freq_table[,P2], freq_table[,P3a])
BABA_1_2_3a <- baba(freq_table[,P1], freq_table[,P2], freq_table[,P3a])

ABBA_1_3b_3a <- abba(freq_table[,P1], freq_table[,P3b], freq_table[,P3a])
BABA_1_3b_3a <- baba(freq_table[,P1], freq_table[,P3b], freq_table[,P3a])
```

Finally, we can estimate *f* as the observed difference between ABBA and BABA to that we see when we substitute P3b for P2.


```R
f <- (sum(ABBA_1_2_3a) - sum(BABA_1_2_3a)) /
     (sum(ABBA_1_3b_3a) - sum(BABA_1_3b_3a))
```

This reveals that nearly 30% of the genome has been shared between *H. melpomene* and *H. cydno* in sympatry. The admixture proportion can be interpreted as the average proportion of foreign ancestry in any given haploid genome. It can also be interpreted as the expected frequency of foreign alleles at any given site in the genome in a population.

