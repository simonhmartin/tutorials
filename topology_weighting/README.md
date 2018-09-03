# Tutorial: Topology Weighting
___
## Requirements
* Python 2.7
* Numpy 1.10+
* R 3.0+
* [ete](http://etetoolkit.org/) 3 (Python module)
* [msms](http://www.mabs.at/ewing/msms/index.shtml)
* [seq-gen](http://tree.bio.ed.ac.uk/software/seqgen/) (for Part 2 only)
* [Phyml](http://www.atgc-montpellier.fr/phyml/) (for Part 2 only)
___
## Introduction

Topology weighting is a means to quantify relationships between taxa that are not necessarily monophyletic. It provides a summary of a complex genealogy by considering simpler "taxon topologies" and quantifying the proportion of sub-trees that match each taxon topology. The method we use to compute the weightings is called *Twisst*: Topology weighting by iterative sampling of sub-trees.

In this practical we will use a simulation to explore how topology weightings provide a summary of the genealogical history. We will then try to infer topology weights across our simulated chromosome using neighbour joining trees inferred for narrow windows.


#### Workflow Overview

In the first part of the practical we will simulate a fairly complex genealogical history including recombination gene flow and selection.
We will then compute topology weightings for all genelogies across this genomic region.
We will plot the weightings in R to see what these tell us about the relationships among the simulated populations.
In the second part of the practical, we will create a simulated sequence alignment and use it as if it was genomic data we had generated.
We will infer genealogies for sequence windows across the region
Finally, we will compute weightings for the inferred trees and compare these with the true weightings.

## Preparation

* Download the collections of python scripts required for this tutorial [GitHub](https://github.com/simonhmartin)

```bash
wget https://github.com/simonhmartin/genomics_general/archive/master.zip
unzip master.zip
rm master.zip

wget https://github.com/simonhmartin/twisst/archive/master.zip
unzip master.zip
rm master.zip
```

* To ensure that the libraries are recognisable by python, add the @genomics_general' directory to the Python path

```bash
export PYTHONPATH=$PYTHONPATH:genomics_general-master
```

* Now make a subdirectory to store the simulation data

```bash
mkdir sim
```

## Practical Part 1. Analysis of simulated Genelogies

#### Coalescent simulation
Our simulation includes five populations each represented by 8 haploid samples. The populations split in the order (((A,B),(C,D)),E). In other words, A and B are sister taxa, C and D are sister taxa, and E is an outgroup to the clade of A, B, C and D. We will add gene flow between C and B. Finally, we will add a beneficial allele to poulation C, which should sweep through population B (i.e. adaptive introgression). The simulated region will be 500 Kb in length, with a population recombination rate of 0.01. The selected locus will be at the centre of the region.

* We use [msms](http://www.mabs.at/ewing/msms/index.shtml) to perform a coalescent simulation with selection.

```bash
msms 40 1 -T -I 5 8 8 8 8 8 -ej 0.5 2 1 -ej 0.5 4 3 -ej 1 3 1 -ej 1.5 5 1 -m 2 3 0.1 -r 5000 500000 \
 -SAA 1000 -SAa 1000 -Sp 0.5 -SI 0.1 5 0 0 1 0 0 -N 100000  |
 grep ";" > sim/msms_5of8_l500k_r5k_sweep.txt
```
The above command will simulate the scenario described above. Note that to modify the evolutionary history, you can change the population joins (`-ej` flags, times are given in coalescent units of 4N generations), migration (`-m` flag), selection strength (`-SAA` and `-SAa`) and selection initiation parameters (`-SI`).

The final part of the command extracts only the tree objects from the output, as these are all we need.

We can view the simulated trees using `less -S sim/msms_5of8_l500k_r5k_sweep.txt`  (type `q` at any time to exit). Each line represents the genealogy for a non-recombining 'block' of sequence, with the the genealogies separated by recombination events. Most blocks are only a few base pairs long (lengths given at the start of each line). You might have noiticed that adjacent genealogies are usually very similar, and sometimes identical. This is because they are usually separated by only a single recombination event, making them highly non-independent. We can check the number of blocks by counting the lines in the file: `wc -l sim/msms_5of8_l500k_r5k_sweep.txt`.

* Next we extract the start and end positions of each genealogy, which will be used for plotting later. We also make a file of trees with the block sizes removed.

```bash
python genomics_general-master/phylo/parse_ms_trees.py --msOutput sim/msms_5of8_l500k_r5k_sweep.txt \
--outTrees sim/msms_5of8_l500k_r5k_sweep.trees.gz --outData sim/msms_5of8_l500k_r5k_sweep.data.tsv.gz
```

The outputs of this command are both zipped with `gzip` to save space. You can view the files using `zcat sim/msms_5of8_l500k_r5k_sweep.data.tsv.gz | less -S` (type `q` to exit).

#### *Twisst*

* We run [*Twisst*](https://github.com/simonhmartin/twisst) to compute the weightings for each topology. The core script is `twisst.py`, but we use the script `run_twisst_parallel.py` to run the script with multi-threading. Adjust the `--threads` argument according to the number of available cores.

The only information *Twisst* requires is which group each sample belongs to. This may be determined by species, phenotype or geography. In our case we know the groupings because we simulated the data. Samples 1:8 will be group 'A', 8:16 group 'B' etc.

```bash
python twisst-master/run_twisst_parallel.py --threads 2 -t sim/msms_5of8_l500k_r5k_sweep.trees.gz --method complete \
-g A 1,2,3,4,5,6,7,8 -g B 9,10,11,12,13,14,15,16 -g C 17,18,19,20,21,22,23,24 \
-g D 25,26,27,28,29,30,31,32 -g E 33,34,35,36,37,38,39,40 |
gzip > sim/msms_5of8_l500k_r5k_sweep.weights.tsv.gz
```
There will be around 70,000 trees to analyse, so this will take some time to run. If you have a multi-core computer, you can increase the number of threads used to increase the speed.

The command above tells twisst to use the `complete` method. This will compute the exact weightings by considering all possible sub-trees. The main output `sim/msms_5of8_l500k_r5k_sweep.weights.tsv.gz` has a line for each input tree with 15 columns, one for each of the possible topologies describing the relationship between the five groups. Each value gives the number of subtrees that match the given topology. Because complete weighting was used, the values in each line should sum to the total number of possible subtrees (8<sup>5</sup> = 32768).

#### Plotting

* Open R and, if necessary, set the working directory to the tutorial directory. In RStudio you can do this using the menus. In the terminal you can use th `setwd()` command.

* Start a new R script to record all the plotting commands

* First we will import the APE library and an additional set of functions from Twisst that will help with plotting.

```R
source("twisst-master/plot_twisst.R")
```

* We define the files containing the weights for each genealogy, and the start and end positions for each block along the chromosome.

```R
#weights file with a column for each topology
weights_file <- 'sim/msms_5of8_l500k_r5k_sweep.weights.tsv.gz'

#coordinates file for each window
window_data_file <- 'sim/msms_5of8_l500k_r5k_sweep.data.tsv.gz'
```

* Read in the weighting data from twosst along with the window positional data

```R
twisst_data <- import.twisst(weights_file, window_data_file)
```

* plot the raw weightings using the provided `plot.twisst` function.

```R
plot.twisst(twisst_data)
```

You might need to enlarge this plot to see it clearly. The trees at the top of the plot show the 15 different topologies we have weighted.

The lower plot shows the weightings. You will see coloumns of colour of varying width. Each column corresponds to a single block with a unique genealogy. Some blocks are all one colour. That indicates that all subtrees in that block have the same topology, indicating a consistent and completely sorted genealogy. Other columns have two or more colours stacked upon one another, indicating that the genealogy has a more complex evolutionary history, with more than one topology represented among the subtrees.

It is often desirable to smooth the weightings so that we can see more clearly how they vary across the chromosome.

* Create smoothed weightings and re-plot.

```R
twisst_data_smooth <- smooth.twisst(twisst_data, span=0.02)

plot.twisst(twisst_data_smooth)
```

Now we see more clearly that the dominant colour is that corresponding to Topology 3, which matches the population branching pattern in the simulation. The fact that other topologies are also pepresented indicates that lineage sorting does not always follow the population branching pattern. In particular, we expect to see some representation of Topology 12, in which B groups with C, and is nested within the C, D clade. This topology is expected to be more common than under null expectations due to the gene flow from population C to B. In particular, we see a major increase in Topology 12 around the selected locu, consistent with adaptive introgression.

Note that Topology 13 also groups B and C, but has a generally low weighting. Gene flow that occurs strictly from C to B will only tend to increase the abundanc of Topology 12, but not 13. Thus, topology weightings carry information about the direction of gene flow.


#### In your own time
How do the weightings change if we make changes to the simulation, such as decreasing the split times, or increasing the rate of migration? (You might need to look at the [msms manual](http://www.mabs.at/ewing/msms/Manual.pdf) to figure out what to change.)

___

## Practical Part 2. Infering weightings from sequence data

Above we have used the 'true' genealogies as they were simulated. In most situations, all we have is sequence data, and its evolutionary history has to be inferred. To demonstrate how this is done for topology weighting, we will use the same simulated genealogies as above to simulate sequence data, and then see whether we can recover the same result using topology weighting based on inferred genealogies.

In other words, we wil **assume that we know nothing about the *true* history of these samples, and try to infer that from simulated sequence data**.

#### Simulating sequences

* Open a new terminal window and navigate to the tutorial folder

* To get sequence data, we will use [[seq-gen](http://tree.bio.ed.ac.uk/software/seqgen/). Seq-gen can take the complete set of genealogies and simulate a single 'recombinant' sequence alignment for our 40 samples. It requires that we provide the number of 'partitions', which corresponds to the number of distinct genealogies for the region. It also requires a model of molecular evolution (we will use Hasegawa, Kishino and Yano 1985) and a scaling factor to convert the branch lengths from coalescent units (4N generations) to mutational distances. We will use 1%.


```bash
partitions=$(wc -l sim/msms_5of8_l500k_r5k_sweep.txt)

seq-gen -mHKY -l 500000 -s 0.01 -p $partitions < sim/msms_5of8_l500k_r5k_sweep.txt |
gzip > sim/msms_5of8_l500k_r5k_sweep.seqgen.phy.gz
```

The output is a sequence alignment in Phylip sequential format. This carries no information about where the recombination breakpoints are, except for the information caried in the sequences themselves.

We will use a simple approach of dividing the sequence into narrow windows and inferring the genealogy for each window. We will use some Python scripts to do this.

* First we need to convert the sequences into a more usable format, in which each row is a site and individual genotypes are given in columns. (NOTE: Conveniently, the simulation provides haploid sequences. Typical genomic data from diploid samples shoudl ideally be phased so that individual haplotypes can be distinguished for tree inference.)

```bash
python genomics_general-master/seqToGeno.py --seqFile sim/msms_5of8_l500k_r5k_sweep.seqgen.phy.gz \
--format phylip --genoFile sim/msms_5of8_l500k_r5k_sweep.seqgen.geno.gz
```

* This file contains both invariant and variant sites (SNPs). We only need the latter to make trees, so we will filter the file using another python script.

```bash
python genomics_general-master/filterGenotypes.py -i sim/msms_5of8_l500k_r5k_sweep.seqgen.geno.gz --minAlleles 2 \
-o sim/msms_5of8_l500k_r5k_sweep.seqgen.SNP.geno.gz  --threads 2
```

#### Infering trees for windows and computing weights

We now have a file of SNPs distributed across out 500 Kb genomic region. We want to infer trees in windows, and we want to select a window size that is **narrow enough to capture the variation in genealogical histroy** across the chromosome, but **contains sufficient SNPs to provide the necessary power for tree inference**. We will use 50 SNP windows in this example.

* The next script reads the SNP file in windows and then infers a tree for each window using [Phyml](http://www.atgc-montpellier.fr/phyml/). Phyml is capable of maximum likelihood inference, but here will will not use optimisation, so the trees output will be Neighbour-Joining trees, inferred using the [BIONJ](http://www.atgc-montpellier.fr/bionj/) algorithm.

```
python genomics_general-master/phylo/phyml_sliding_windows.py -g sim/msms_5of8_l500k_r5k_sweep.seqgen.SNP.geno.gz \
--windType sites -w 50 --prefix sim/msms_5of8_l500k_r5k_sweep.seqgen.SNP.w50.phyml_bionj --model HKY85 --optimise n --threads 1
```

This generates two output files: `sim/msms_5of8_l500k_r5k_sweep.seqgen.SNP.w50.phyml_bionj.trees.gz` contains the trees for each window. `sim/msms_5of8_l500k_r5k_sweep.seqgen.SNP.w50.phyml_bionj.data.tsv` contains the coordinates for each window.

* Finally, we can compute the weightings for each of these inferred window trees.

```bash
python twisst-master/run_twisst_parallel.py --method complete --threads 2 -t sim/msms_5of8_l500k_r5k_sweep.seqgen.SNP.w50.phyml_bionj.trees.gz \
-g A 1,2,3,4,5,6,7,8 -g B 9,10,11,12,13,14,15,16 -g C 17,18,19,20,21,22,23,24 \
-g D 25,26,27,28,29,30,31,32 -g E 33,34,35,36,37,38,39,40 |
gzip > sim/msms_5of8_l500k_r5k_sweep.seqgen.SNP.w50.phyml_bionj.weights.tsv.gz
```
#### Plotting inferred weights

As we did above, we can now plot the weights for these inferred trees across the chromosome. This can be done in the same R script as before.

* **Open R again** (if you have restarted R, you may need to reload the `plot_twisst.R` script).

```R
source("twisst-master/plot_twisst.R")
```

* As before we read in the weights and window data files, and normalise the weights. Here we will load both the ***true*** weights and the ***inferred*** weights that we have just computed.

```R
weights_files = c("sim/msms_5of8_l500k_r5k_sweep.weights.tsv.gz",
                  "sim/msms_5of8_l500k_r5k_sweep.seqgen.SNP.w50.phyml_bionj.weights.tsv.gz")

window_data_files = c("sim/msms_5of8_l500k_r5k_sweep.data.tsv.gz",
                      "sim/msms_5of8_l500k_r5k_sweep.seqgen.SNP.w50.phyml_bionj.data.tsv")
```

* Load in wieghhtings and window data and run smoothing again

```R
twisst_data <- import.twisst(weights_files, window_data_files)

twisst_data_smooth <- smooth.twisst(twisst_data, span=0.02)
```

* Now we plot again to compare the true and inferred weightings.

```R
plot.twisst(twisst_data_smooth)
```
How well did the inference capture the truth?

#### In your own time:
What effect does changing the window size have on inference? For example, what if we inferred treas for windows of 500 SNPs, or just 20 SNPs? Can you explain this behaviour in terms of the tradeoff between power and resolution? 

