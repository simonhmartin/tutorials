# Tutorial: *ABBA* *BABA* analysis in sliding windows
___
## Requirements
Python 2.7
Numpy 1.10+
R 3.0+

___
## Introduction

ABBA BABA statistics (also called D statistics) provide a simple and powerful test for a deviation from a strict bifurcating evolutionary history. They are therefore frequently used to test for introgression using genome-scale SNP data.

Although originally developed to be employed for genome-wide tests for introgression, they can also be applied in smaller windows, which can allow exploration of the **genomic landscape of introgression**.

In this practical we will perform window-based ABBA BABA analysis using a **available software** and then **write code in R for plotting the results**. We will analyse genomic data from several populations of *Heliconius* butterflies.

Note that the details and theory of ABBA BABA statistics are more fully explored in my other [tutorial on whole genome ABBA BABA analyses](https://github.com/simonhmartin/tutorials/tree/master/ABBA_BABA_whole_genome)

#### Workflow Overview

Starting with genotype data from **whole-genome sequencing** of multiple individuals, we will run a script that computes a measure of the admixture proportion in individual windows across each chromosome. We will then make plots to test hypotheses about adaptive introgression.

#### Data

We will study multiple races from three species: *Heliconius melpomene*, *Heliconius timareta* and *Heliconius cydno*. These species have partially overlapping ranges and they are thought to hybridise where they occur in sympatry. Our sample set includes two pairs of sympatric races of *H. melpomene* and *H. cydno* from Panama and the western slopes of the Andes in Colombia. There are also two pairs of sympatric races of *H. melpomene* and H. timareta from the eastern slopes of the Andes in Colombia and Peru. Finally, there are two samples from an outgroup species *Heliconius numata*, which are necessary for performing the ABBA BABA analyses.

All samples were sequenced using high-depth **whole-genome sequencing**, and genotypes have been called for each individual for each site in the genome using a standard pipeline. The data has been filtered to retain only **bi-allelic** single nucleotide polymorphisms (SNPs). This dataset includes SNP data from chromosome 18, which is known to carry a wing patterning locus of particular interest.

#### Hypotheses

We hypothesize that hybridisation between species in sympatry will lead to sharing of genetic variation between *H. cydno* and the **sympatric** races of *H. melpomene* from the west, and between *H. timareta* and the corresponding sympatric races of *H. melpomene* from the east of the Andes.

However, not all parts fo the genome are expected to be equally affected. In particular, we suspect that the wing patterninging gene *optix* on chromosome 18 has been under strong selection. Differential regulation of *optix* can give rise to different distributions of red pigmentation on the wing, as seen in different subspecies of *H. melpomene*, or the absence of red, as seen in *H. cydno*.

*Heliconius* wing patterns act as warnings to predators that they are toxic. Some species participate in Müllerian mimicry, whereby toxic species have evolved to resemble one-another, which helps to reinfoce predator learning. Mimicry can either be achived through independent convergence on the same wing patterns, or through exchange of wing patterning alleles through **adaptive introgression**. We therefore predict that co-mimetic populations of different species might show an excess signal of introgression in the vicinity of *optix*.

We have an entirely different expectation for populations with different wing patterns. If predators in a given area recognise the most commobn local pattern as toxic, it will be costly to have a foreign wing pattern. Likewise any hybrid individual that has an intermediate wing pattern will also be at risk of higher predation. We therefore predict that between populations with different wing patterns there should be a reduction in the extent of introgression in the vicinity of *optix*.


![Species Map](images/map_and_tree.jpg)


#### Quantifying admixture across the genome

A detailed explanation of the *ABBA BABA* test is given in in my other [tutorial on whole genome ABBA BABA analyses](https://github.com/simonhmartin/tutorials/tree/master/ABBA_BABA_whole_genome).

Briefly, the test uses three populations and an outgroup with the relationship (((P1,P2),P3),O), and investigates whether there is an excess of shared variation between P2 and P3 (compared to that shared between P1 and P3).

This excess can be expressed in terms of the *D statistic*, which ranges from -1 to 1, and should equal 0 under the null hypothesis of no introgression. D > 1 indicates possible introgression between P3 and P2 (or other factors that would result in a deviation from a strict bifurcating species history).

This test was designed to be used at the whole-genome scale. The *D* statistic is not well suited for comparing admixture levels across the genome, because its absolute value depends on factors such as the effective population size, which can vary across the genome.

The *f* estimator descibed in the [other tutorial](https://github.com/simonhmartin/tutorials/tree/master/ABBA_BABA_whole_genome) is better, because it by definition reflects the admixture proportion, but it is highly sensitive to stochastic error at the small scale. A statistic called *f<sub>d</sub>* was therefore developed for this purpose that is more robust to the error introduced by using small numbers of SNPs ([Martin et al. 2015](https://doi.org/10.1093/molbev/msu269)). While the conventional f estimator assumes that P3 is the donor population and P2 the recipient, *f<sub>d</sub>* infers the donor on a site-by-site basis.

#### Selecting populations

The interpretation of these statistics is strongly dependent on the populations selected. Firstly, the test is most sensitive to introgression from **P3 into P2**, rather than the other way around.

Secondly, *f<sub>d</sub>* should be interpreted as a **quantification of excess shared variation between P3 and P2** that is **not also shared with P1**. If there is ongoing gene flow between P1 and P2, then any introgression from P3 to P2 will be underestimated.

Finally, we are only able to quantify introgression that occured **more recently than the split between P1 and P2**.

Therefore, if we want to quantify the **maximum amount of detectable introgression** that has occurred across the genome, we should choose a P1 that is **allopatric and not too closely related to P2**.

However, we can also this feature of the test to our advantage. If we select a P1 that shares ongoing gene flow with P2, then the test will instead be revealing parts of the genome at which **P2 and P3 share variation that is not shared by P1**. This can be useful for identifying wing patterning alleles, as these are often the only genomic regions at which subspecies remain distinct in the face of gene flow.

## Practical

### Preparation

* Open a terminal window and navigate to a folder where you will run the excersise and store all the input and output data files.

* This tutorial makes use of a collection of python scripts that must be downloaded from [GitHub](https://github.com/simonhmartin)

```bash
git clone https://github.com/simonhmartin/genomics_general
```

### Sliding window analysis

* Run the the analysis python script for two separate cases. In both, P1 is the allpatric *H. melpomene melpomene* (`mel_mel`). P2 and P3 are the two populations we expect to be sharing genes. In the first case we are quantifying introgression between *H. melpomene rosina* (`mel_ros`) and *H. cydno chioneus* (`cyd_chi`) both from Panama. In the second we are quantifying introgression between *H. melpomene amaryllis* (`mel_ama`) and *H. timareta thelxinoe* (`tim_txn`) both from Peru.

```bash
python genomics_general/ABBABABAwindows.py \
-g data/hel92.DP8HET75MP9BIminVar2.chr18.geno.gz -f phased \
-o data/hel92.DP8HET75MP9BIminVar2.chr18.ABBABABA_mel_ros_chi_num.w25m250.csv.gz \
-P1 mel_mel -P2 mel_ros -P3 cyd_chi -O num \
--popsFile data/hel92.pop.txt -w 25000 -m 250 --T 2

python genomics_general/ABBABABAwindows.py \
-g data/hel92.DP8HET75MP9BIminVar2.chr18.geno.gz -f phased \
-o data/hel92.DP8HET75MP9BIminVar2.chr18.ABBABABA_mel_ama_txn_num.w25m250.csv.gz \
-P1 mel_mel -P2 mel_ama -P3 tim_txn -O num \
--popsFile data/hel92.pop.txt -w 25000 -m 250 --T 2
```

We provide the scripty it with an input file containing genotype data (`-g`), an output file (`-o`), ingroup populations and outgroup (`-P1`, `-P2`, `-P3` and `-O`), and a file specifying which population each sample is in (`--popsFile`).

We also give parameters for the windows. These ae "coordinate" windows, which means each window is the same length relative to the reference genome, but the number of SNPs per window can vary. The window size (`-w`) will be 25,000 bp. Windows will be required to contain a minimum (`-m`) of 250 SNPs to be considered valid.

Finally, we tell the script to use two threads (`-T`). If you have a multi-core machine, you can increase this value and the script will run faster.

#### Plotting window statistics

* Open R

* Set the working directory to the tutorial directory. You can do this with the `setwd()` command, or in RStudio using the menus.

We need to load each file of window statistics into R. We will make a list containing both datasets. 

* First input the names of teh input files

```R
file_names <- c("data/hel92.DP8HET75MP9BIminVar2.chr18.ABBABABA_mel_ama_txn_num.w25m250.csv.gz",
                "data/hel92.DP8HET75MP9BIminVar2.chr18.ABBABABA_mel_ros_chi_num.w25m250.csv.gz")

stats_tables = lapply(file_names, read.csv)
```

*f<sub>d<sub>* is meaningless when D is negative, as it is designed to quantify the excess of ABBA over BABA only whgen an excess exists.

* We therefore convert all *f<sub>d</sub>* values to 0 at sites where *D* is negative. 

```R
for (x in 1:length(stats_tables)){
stats_tables[[x]]$fd = ifelse(stats_tables[[x]]$D < 0, 0, stats_tables[[x]]$fd)
    }
```

* We can then plot of *f<sub>d</sub>* across the chromosome for the two cases we have analysed.

```R
par(mfrow=c(length(stats_tables), 1), mar = c(4,4,1,1))

for (x in 1:length(stats_tables)){
    plot(stats_tables[[x]]$mid, stats_tables[[x]]$fd,
    type = "l", xlim=c(0,17e6),ylim=c(0,1),ylab="Admixture Proportion",xlab="Position")
    rect(1000000,0,1250000,1, col = rgb(0,0,0,0.2), border=NA)
    }
```

This reveals that there is considerable heterogeneity in the extent of admixture across the chromosome. If we consider the region around optix, we see evidence for reduced introgression between *H. melpomene rosina* and *H. cydno chioneus*, as we predicted. By contrast, we see evidence for elevated introgression between *H. melpomene amaryllis* and *H. timareta thelxinoe*, which suggests that their shared wing patterns might result from adaptive introgression. Given this evidence, it would be recommended to make a phylogeny for the region around optix to test whether the H. timareta allele appears to be 'nested' within the H. melpomene clade. In this case, previous papers have confirmed that that is the case ([Pardo-Diaz et al. 2012](https://doi.org/10.1371/journal.pgen.1002752), [Wallbank et al. 2016](https://doi.org/10.1371/journal.pbio.1002353)).

#### In your own time
What happens when we change the identity of P1, P2 an P3? What happend if we change the window size?

