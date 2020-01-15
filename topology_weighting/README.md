# Tutorial: Topology Weighting
___
## Requirements
* `Python 3+
    * Numpy 1.10+
    * [ete3](http://etetoolkit.org/)
    * [msprime](https://msprime.readthedocs.io/en/stable/) (for part 3 only)
* R 3.0+
    * ape
    * data.table
* [Phyml](http://www.atgc-montpellier.fr/phyml/) (for part 2 only)

___
## Introduction

Topology weighting is a means to quantify relationships among taxa that are not necessarily monophyletic. It provides a summary of a complex genealogy by considering simpler "taxon topologies" and quantifying the proportion of sub-trees that match each taxon topology. The method we use to compute the weightings is called *Twisst*: Topology weighting by iterative sampling of sub-trees.

In this practical we will use a simulation to explore how topology weightings provide a summary of the genealogical history. We will then try to infer topology weights across our simulated chromosome using neighbour joining trees inferred for narrow windows.


#### Workflow Overview

In **Part 1** of the practical we will analyse a set of genealogies that represent the history of a part of chromosome that evolved under a fairly complex history including population subdivision, gene flow and selection. We will compute topology weightings across this genomic region `Twisst` and then explore the results in `R`.

In **Part 2**, we take a step backwards into the real world, in which we ***don't know*** the true genealogical history, but instead we have a set of sequences from which we hope to ***infer*** this history. We will use an unsophisticated approach to do this: making phylogenies for windows across the genome, using a standard phylogenetics tool. By comparing our inferred histories to the truth in `R`, we will gain insights into the **tradeoff between power and resolution** in genealogical inference. We will then apply this understanding by analysing real sequence data from *Heliconius* butteflies.

In **Part 3** we will explore a new way of thinking about genealogies: the **tree sequence**. We will simulate a tree sequence for a chromosome using `msprime` and compute topology weightings for it using `twisst`. We will then also compute the weightings from a tree sequence inferred from human sequence data, and visualise the results. 


## Practical Part 1. Analysis of simulated genelogies

#### Download code and data

The scripts and example data for this part of the practical are in the `twisst` package on github.

```bash
wget https://github.com/simonhmartin/twisst/archive/v0.2.tar.gz
tar -xzf v0.2.tar.gz
rm v0.2.tar.gz
```

* The example data we will use consists of a text file of genealogies coded as 'newick' trees. We can look at the first tree in the file:

```bash
zcat twisst-0.2/examples/msms_4of10_l50k_r500_sweep.trees.gz | head -n 1
```

The numbers before each `:` are the sample names. The numbers after the `:` are the branch lengths. We will only be considering the tree shape and not branch lengths in this tutorial.

* We can also check the total number of distinct genealogies for this region of the chromosome:

```bash
zcat twisst-0.2/examples/msms_4of10_l50k_r500_sweep.trees.gz | wc -l
```

* For plotting, we also need to know where these genealogies occur on the chromosome. This data is provided in a second file with three columns: chromosome, start and end for each genealogy.

```bash
zcat twisst-0.2/examples/msms_4of10_l50k_r500_sweep.data.tsv.gz | head
``` 

As you can see, some genealogies occupy very narrow regions of the chromosome, as small as 1 bp.

#### Compute topology weightings

* We run [`twisst`](https://github.com/simonhmartin/twisst) to compute the weightings for each topology.

The only information *Twisst* requires is which group each sample belongs to. This may be determined by species, phenotype or geography. In our case there are four groups of 10 haploid samples each. Samples 1:10 will be group 'A', 11:20 group 'B' etc.

```bash
python twisst-0.2/twisst.py \
-t twisst-0.2/examples/msms_4of10_l50k_r500_sweep.trees.gz \
-w msms_4of10_l50k_r500_sweep.weights.tsv.gz \
-g A 1,2,3,4,5,6,7,8,9,10 \
-g B 11,12,13,14,15,16,17,18,19,20 \
-g C 21,22,23,24,25,26,27,28,29,30 \
-g D 31,32,33,34,35,36,37,38,39,40
```

This command tells `twisst` to consider all possible combinations samples in which there is one sample per group. For example the first combination examined will be samples `1`, `11`, `21` and `31`, representing groups `A`, `B`, `C` and `D`, respectively. Ignoring all other branches in the tree, Twisst extracts the subtree for the four samples of interest and records its topology, which could have one of three possible shapes: `(((A,B),C),D)`, `(((A,C),B),D)` or `(((B,C),A),D)` (Note that here the trees are represented as rooted, with D as the outgroup, but in reality `twisst` does not consider the rooting and it does not change the result).

Check what the results look like

```bash
#first 10 lines
zcat msms_4of10_l50k_r500_sweep.weights.tsv.gz | head -n 10
#total number of lines
zcat msms_4of10_l50k_r500_sweep.weights.tsv.gz | wc -l
```

The three columns in the weights file represent the three topologies, which are defined in the file too. The numbers are not proportions, because twisst reports the total number of combinations representing each topology. 

You might see that some adjacent lines have identical weightings. This has to do with the fact that some recombination events change the relationships among the samples, but not in a way that influences the weightings. This is something to think about.

#### Analyse the results

* Open `R` or `RStudio` and, if necessary, set the working directory to where you have saved the files. You can use the you can use the `setwd()` command or, in `RStudio`, using the menus.

* Start a new R script to record the commands

* First we will import a set of functions distributed with `twisst` that will help with plotting.

```R
source("twisst-0.2/plot_twisst.R")
```

* Please note the cool name of the above script.

* We define the files containing the weights for each genealogy, and the start and end positions for each block along the chromosome.

```R
#weights file with a column for each topology
weights_file <- 'sim/msms_5of8_l500k_r5k_sweep.weights.tsv.gz'

#coordinates file for each window
window_data_file <- 'sim/msms_5of8_l500k_r5k_sweep.data.tsv.gz'
```
* We already know the structure of these two files, but instead of reading them in and working with them directly, we will use the convenient `import.twisst` function.

```R
twisst_data <- import.twisst(weights_file, window_data_file)
```

* plot the raw weightings using the provided `plot.twisst` function.

```R
plot.twisst(twisst_data)
```

You might need to enlarge this plot to see it clearly. The trees at the top of the plot show the 3 different topologies we have weighted.

The lower plot shows the weightings. You will see coloumns of colour of varying width. Each column corresponds to a single block with a unique genealogy. Some blocks are all one colour and reach a value of 1. That indicates that all subtrees in that block have the same topology, indicating a consistent and completely sorted genealogy. Other columns have two or more colours overlayed, indicating that the genealogy is more complex evolutionary history, with individuals jumping between groups. A completely random genealogy, in which there is no clustring by group, would have equal weightings for all three topologies.

It is often desirable to smooth the weightings so that we can see more clearly how they vary across the chromosome.

* Create smoothed weightings using the `smooth.twisst` function and re-plot.

```R
twisst_data_smooth <- smooth.twisst(twisst_data, span_bp=5000)

plot.twisst(twisst_data_smooth)
```

This averaged the weightings over a 5,000 bp window. You can explore what happens when you change the `span_bp` parameter.

Now we see more clearly that the dominant topologies are `topo1` and `topo3`. In this case, the simulations involved populations spliting acording to `topo1`, but adaptive introgression was simulated from `C` into `B`, which is why `topo3` is more prevalent than topo2, and also why topo3 has a large spike in the middle of the region. This is the location of the selected locus. 

We can look at the overall destribution of weightings too, or just check the mean values.

```R
plot.twisst.summary.boxplot(twisst_data)

twisst_data$weights_mean
```

So if the simulation history followed `topo1`, and introgression created `topo3`, why is `topo2` not zero?

Lineage sorting is often incomplete complete if the taxa split recently, and even when it is complete, we can find genealogies that are discordant with the 'species tree' due to stochasticity in lineage sorting in the past. If you look carefully at the first plot we made, you might find one narrow window in which `topo2` has a weighting of 1. This indicates a *completely sorted*, but *discordant* genealogy. If the difference between incomplete lineage sorting and discordance is not immediately clear to you, you are not alone. In Part 3 we will do our own simulations to look at the condistions under which incomplete sorting and discordance increase or decrease.

___

## Practical Part 2. Infering weightings from sequence data

Above we have used the 'true' genealogies as they were simulated. In most cases, all we have is sequence data, and its evolutionary history has to be inferred. In fact there are two things we do not know:
1. We do not know the genealogical relationship among all individuals
2. We do not know the 'breakpoints' at which recombination has changed the relationship as we move along the chromosome

In this part, we will start from sequence data (the sequences were simulated under the history covered in Part 1, but we pretend that we do not know that at this stage). We will use a fairly straightforward approach in which we infer genealogies in windows along the genome. We will tthen run `twisst` on these to see whether we can recover someting close to the underlying truth.

Note that one of the lessons in this part is that inferring trees in windows is **crude and potentially deeply flawed**. In Part 3 we will touch on alternative approaches.

#### Download code and data

The scripts for this part are in the genomics_general package on github:

```bash
wget https://github.com/simonhmartin/genomics_general/archive/v0.3.tar.gz
tar -xzf v0.3.tar.gz
rm v0.3.tar.gz
```

* To ensure that the libraries are recognisable by python, add the `genomics_general' directory to the Python path

```bash
export PYTHONPATH=$PYTHONPATH:genomics_general-0.3
```

* The sequence file we will use is also part of the `twisst` package, downloaded in Part 1. The file is in simple `.geno` format, which has columns for chromosome, position and genotyope for each individual:

```bash
zcat twisst-0.2/examples/msms_4of10_l50k_r500_sweep.seqgen.SNP.geno.gz | head
```

#### Infering trees for windows

We have a file of SNPs distributed across a 50 kb genomic region. We will infer trees in windows of a defined number of SNPs, such that each window has a similar amount of information, but might differ in its absolute span across the chromosome, depending on the SNP density.

There is an underlying **tradeoff** in this approach. We want to select a window size that is ***large enough*** to provide the necessary ***power*** for tree inference, but ***small enough*** to achive enough ***resolution*** to captue how genealogical histris change across the chromosome.

* We will run the script that reads the SNP file in windows and then infers a tree for each window using [Phyml](http://www.atgc-montpellier.fr/phyml/). Phyml is capable of maximum likelihood inference, but here will will not use optimisation, so the trees output will be Neighbour-Joining trees, inferred using the [BIONJ](http://www.atgc-montpellier.fr/bionj/) algorithm.

* we run the script four times, using a range of different window sizes

```
for x in 20 50 100 500
do
echo "Inferring trees with window size $x"

python genomics_general-0.3/phylo/phyml_sliding_windows.py \
-g twisst-0.2/examples/msms_4of10_l50k_r500_sweep.seqgen.SNP.geno.gz \
--prefix msms_4of10_l50k_r500_sweep.seqgen.SNP.w$x.phyml_bionj \
--windType sites -w $x  --model HKY85 --optimise n

done
```

Each time it ran, the script generated two output files: `.trees.gz` files contain the trees for each window. `.data.tsv` files contain the coordinates for each window.

* We can now compute the weightings across the chromosome using the trees files as input.

```bash
for x in 20 50 100 500
do
echo "Running Twisst for window size $x"

python twisst-0.2/twisst.py \
-t msms_4of10_l50k_r500_sweep.seqgen.SNP.w$x.phyml_bionj.trees.gz \
-w msms_4of10_l50k_r500_sweep.seqgen.SNP.w$x.phyml_bionj.weights.tsv \
-g A 1,2,3,4,5,6,7,8,9,10 \
-g B 11,12,13,14,15,16,17,18,19,20 \
-g C 21,22,23,24,25,26,27,28,29,30 \
-g D 31,32,33,34,35,36,37,38,39,40 \

done
```

#### Plotting inferred weights

As we did in Part 1, we can now plot the weights for these inferred trees across the chromosome. This can be done in the same R script as before.

* **Open R again** (if you have restarted R, you may need to reload the `plot_twisst.R` script).

```R
source("twisst-0.2/plot_twisst.R")
```

* As before we read in the weights and window data files. This time we will load the original ***true*** weights from the simulated genealogies, as well as the four files of ***inferred*** weights that we have just computed.


```R
weights_files <- c('msms_4of10_l50k_r500_sweep.weights.tsv.gz',
                   'msms_4of10_l50k_r500_sweep.seqgen.SNP.w20.phyml_bionj.weights.tsv',
                   'msms_4of10_l50k_r500_sweep.seqgen.SNP.w50.phyml_bionj.weights.tsv',
                   'msms_4of10_l50k_r500_sweep.seqgen.SNP.w100.phyml_bionj.weights.tsv',
                   'msms_4of10_l50k_r500_sweep.seqgen.SNP.w500.phyml_bionj.weights.tsv')

window_data_files <- c('twisst-0.2/examples/msms_4of10_l50k_r500_sweep.data.tsv.gz',
                      'msms_4of10_l50k_r500_sweep.seqgen.SNP.w20.phyml_bionj.data.tsv',
                      'msms_4of10_l50k_r500_sweep.seqgen.SNP.w50.phyml_bionj.data.tsv',
                      'msms_4of10_l50k_r500_sweep.seqgen.SNP.w100.phyml_bionj.data.tsv',
                      'msms_4of10_l50k_r500_sweep.seqgen.SNP.w500.phyml_bionj.data.tsv')
```

* Load in wieghhtings and window data. Note that when given multiple input files, the `import.twisst` function will interpret them as separate chromosomes.

```R
twisst_data <- import.twisst(weights_files, window_data_files)

```

* Now we plot again to compare the true and inferred weightings. (you might need to expand your plot window to display the multiple plots correctly).

```R
plot.twisst(twisst_data, show_topos=FALSE)

```

How well did the inference capture the truth? Which window size is best? Where is there too little power, and where is there too little resolution?


## Part 3 (coming soon!)
