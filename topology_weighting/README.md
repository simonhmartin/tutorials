# Tutorial: Topology Weighting
___
## Requirements
* Python 3+
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

In **Part 2**, we take a step backwards into the real world, in which we ***don't know*** the true genealogical history, but instead we have a set of sequences from which we hope to ***infer*** this history. We will use an unsophisticated approach to do this: making phylogenies for windows across the genome, using a standard phylogenetics tool. By comparing our inferred histories to the truth in `R`, we will gain insights into the **tradeoff between power and resolution** in genealogical inference.

In **Part 3** we will explore a new way of thinking about genealogies: the **tree sequence**. We will simulate a tree sequence for a chromosome using `msprime` and compute topology weightings for it using `twisst`. We will then explore what happens when we change parameters of the simulation.

___
## Practical Part 1. Analysis of simulated genealogies

#### Download code and data

The scripts and example data for this part of the practical are in the `twisst` package on github.

```bash
wget https://github.com/simonhmartin/twisst/archive/v0.2.1.tar.gz
tar -xzf v0.2.1.tar.gz
rm v0.2.1.tar.gz
```

* The example data we will use consists of a text file of genealogies coded as 'newick' trees. We can look at the first tree in the file:

```bash
zcat twisst-0.2.1/examples/msms_4of10_l50k_r500_sweep.trees.gz | head -n 1
```

It's pretty ugly, but don't be afraid. The numbers before each `:` are the sample names. The numbers after the `:` are the branch lengths. We will only be considering the tree shape and not branch lengths in this tutorial.

* We can also check the total number of distinct genealogies for this region of the chromosome:

```bash
zcat twisst-0.2.1/examples/msms_4of10_l50k_r500_sweep.trees.gz | wc -l
```

* For plotting, we also need to know where these genealogies occur on the chromosome. This data is provided in a second file with three columns: chromosome, start and end for each genealogy.

```bash
zcat twisst-0.2.1/examples/msms_4of10_l50k_r500_sweep.data.tsv.gz | head
``` 

As you can see, some genealogies occupy very narrow regions of the chromosome, as small as 1 bp.

#### Compute topology weightings

* We run [`twisst`](https://github.com/simonhmartin/twisst) to compute the weightings for each topology.

The only information *Twisst* requires is which group each sample belongs to. This may be determined by species, phenotype or geography. In our case there are four groups of 10 haploid samples each. Samples 1:10 will be group 'A', 11:20 group 'B' etc.

```bash
python twisst-0.2.1/twisst.py \
-t twisst-0.2.1/examples/msms_4of10_l50k_r500_sweep.trees.gz \
-w msms_4of10_l50k_r500_sweep.weights.tsv.gz \
-g A 1,2,3,4,5,6,7,8,9,10 \
-g B 11,12,13,14,15,16,17,18,19,20 \
-g C 21,22,23,24,25,26,27,28,29,30 \
-g D 31,32,33,34,35,36,37,38,39,40 \
--outgroup D
```

This command tells `twisst` to consider all possible combinations samples in which there is one sample per group. For example the first combination examined will be samples `1`, `11`, `21` and `31`, representing groups `A`, `B`, `C` and `D`, respectively. Ignoring all other branches in the tree, `twisst` records the topology of the 'subtree' containing just the four samples of interest, which could have one of three possible shapes: `(((A,B),C),D)`, `(((A,C),B),D)` or `(((B,C),A),D)` (Note that here the trees are represented as rooted, with D as the outgroup. We can tell `twisst` which is the outgroup so it displays the trees as correctly rooted, but this does not affect the results).

Check what the results look like

```bash
#first 10 lines
zcat msms_4of10_l50k_r500_sweep.weights.tsv.gz | head -n 30
#total number of lines
zcat msms_4of10_l50k_r500_sweep.weights.tsv.gz | wc -l
```

The three columns in the weights file represent the three topologies, which are defined in the file too. The numbers are not proportions, because twisst reports the total number of combinations representing each topology. 

You will see that some adjacent lines have identical weightings. This has to do with the fact that some recombination events change the relationships among the samples, but not in a way that influences the weightings. It's worth understanding why this is the case.

#### Analyse the results

* Open `R` or `RStudio` and, if necessary, set the working directory to where you have saved the files. You can use the you can use the `setwd()` command or, in `RStudio`, using the menus.

* Start a new R script to record the commands

* First we will import a set of functions distributed with `twisst` that will help with plotting.

```R
source("twisst-0.2.1/plot_twisst.R")
```

* Please note the cool name of the above script.

* We define the files containing the weights for each genealogy, and the start and end positions for each block along the chromosome.

```R
#weights file with a column for each topology
weights_file <- 'msms_4of10_l50k_r500_sweep.weights.tsv.gz'

#coordinates file for each window
window_data_file <- 'twisst-0.2.1/examples/msms_4of10_l50k_r500_sweep.data.tsv.gz'
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

This averaged the weightings over a 5 kb window. You can explore what happens when you change the `span_bp` parameter.

Now we see more clearly that the dominant topologies are `topo1` and `topo3`. In this case, the simulations involved populations spliting acording to `topo1`, but adaptive introgression was simulated from `C` into `B`, which is why `topo3` is more prevalent than topo2, and also why topo3 has a large spike in the middle of the region. This is the location of the selected locus. 

* We can look at the overall destribution of weightings too, or just check the mean values.

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

* The scripts for this part are in the genomics_general package on github, which we need to download:

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
zcat twisst-0.2.1/examples/msms_4of10_l50k_r500_sweep.seqgen.SNP.geno.gz | head -n 3 | cut -f 1-8
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
-g twisst-0.2.1/examples/msms_4of10_l50k_r500_sweep.seqgen.SNP.geno.gz \
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

python twisst-0.2.1/twisst.py \
-t msms_4of10_l50k_r500_sweep.seqgen.SNP.w$x.phyml_bionj.trees.gz \
-w msms_4of10_l50k_r500_sweep.seqgen.SNP.w$x.phyml_bionj.weights.tsv \
-g A 1,2,3,4,5,6,7,8,9,10 \
-g B 11,12,13,14,15,16,17,18,19,20 \
-g C 21,22,23,24,25,26,27,28,29,30 \
-g D 31,32,33,34,35,36,37,38,39,40 \
--outgroup D

done
```

#### Plotting inferred weights

As we did in Part 1, we can now plot the weights for these inferred trees across the chromosome. This can be done in the same R script as before.

* **Open R again** (if you have restarted R, you may need to reload the `plot_twisst.R` script).

```R
source("twisst-0.2.1/plot_twisst.R")
```

* As before we read in the weights and window data files. This time we will load the original ***true*** weights from the simulated genealogies, as well as the four files of ***inferred*** weights that we have just computed.


```R
weights_files <- c('msms_4of10_l50k_r500_sweep.weights.tsv.gz',
                   'msms_4of10_l50k_r500_sweep.seqgen.SNP.w20.phyml_bionj.weights.tsv',
                   'msms_4of10_l50k_r500_sweep.seqgen.SNP.w50.phyml_bionj.weights.tsv',
                   'msms_4of10_l50k_r500_sweep.seqgen.SNP.w100.phyml_bionj.weights.tsv',
                   'msms_4of10_l50k_r500_sweep.seqgen.SNP.w500.phyml_bionj.weights.tsv')

window_data_files <- c('twisst-0.2.1/examples/msms_4of10_l50k_r500_sweep.data.tsv.gz',
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


## Practical Part 2B: Infering weightings from sequence data

If you're tired of simulated data, here's a chance to work with real data. Note this part of the practical takes a bit longer, so you might want to skip it if you are low on time.

We will now run the same steps as above on real data. The input data are phased sequences from [this paper on *Heliconius* butterflies](https://doi.org/10.1371/journal.pbio.2006288). Fig 1 of the paper shows the makeup of the data set, with 9 populations of 10 samples each and two outgroup individuals.

We are interested in the role of gene flow in shaping the relationships between three species. *H. cydno* and *H. timarata* are sister species that occur on opposite sides of the Andes Mountains. *H. melpomene* occurs on both sides of the Andes, so there is sympatry between *cydno* and *melpomene* to the west and between *timareta* and *melpomene* on the eastern slopes. Hybrids are very rare in the wild, but genomic evidence indicates that gene flow occurs in both areas of sympatry.

The sampling includes two cydno populations (also called races because they have different colour patterns), two timareta populations and five melpomene populations. There are also two outgroup individuals.

#### Infer trees for windows

There are many combinations of populations we could use. Regardless of which combination we use, we need to start by computing trees across the genome for the complete set of samples.

* In this case we have a phased vcf file for a 1 mb region of chromosome 18, so our command to infer the trees has two steps. First we convert the vcf into the correct format, and then pipe it to the `phyml_sliding_windows` script to infer the trees.

```bash
python genomics_general-0.3/VCF_processing/parseVCF.py -i twisst-0.2.1/examples/heliconius92.chr18.500001-1500000.phased.vcf.gz |
python genomics_general-0.3/phylo/phyml_sliding_windows.py --threads 2 \
--prefix heliconius92.chr18.500001-1500000.phyml_bionj \
--windType sites -w 50 --model GTR --optimise n
```
We've set the window size to 50 SNPs, because that seemed a reasonable compromise based on the simulations above.

This could take 10-15 minutes (~1300 trees to make), so it's a good time to get a cup of coffee.

#### Compute topology weights

For `twisst`, we usually select 4 or 5 taxa. Any more and the number of possible topologies becomes large, so you would need to have a clear hypothesis about which particular topologies you want to focus on.

Here we will focus on two sympatric pairs. *H. cydno chioneus* ('chi') and H. melpomene rosina ('ros') from Panama, and *H. timareta thelxinoe* ('txn') and *H. melpomene amaryllis* ('ama') from Peru. Hybridisation occurs in both location, but you will notice that while the Panama pair have divergent colour patterns, the Peru pair are identical. This is hypothesised to have resulted from adaptive introgression.

The 1 Mb region on chromosome 18 that we are targeting contains the gene optix, which is the controler of the red forwing band shared by the pair in Peru.

* We will run `twisst` for this pair of taxa, but instead of specifying all individuals that belong to each group, we just provide a groups file. You can view the format of the populations file.

```bash
head -n 25 twisst-0.2.1/examples/heliconius92.pop.txt
```

* Then run `twisst`

```bash
python twisst-0.2.1/twisst.py -t heliconius92.chr18.500001-1500000.phyml_bionj.trees.gz \
-w heliconius92.chr18.500001-1500000.phyml_bionj.weights.tsv \
-g chi -g txn -g ros -g ama --groupsFile twisst-0.2.1/examples/heliconius92.pop.txt
```

#### plotting the result

* **Open R again** (if you have restarted R, you may need to reload the `plot_twisst.R` script).

```R
source("twisst-0.2.1/plot_twisst.R")
```

```R
weights_file = "heliconius92.chr18.500001-1500000.phyml_bionj.weights.tsv"
data_file = "heliconius92.chr18.500001-1500000.phyml_bionj.data.tsv"

twisst_data <- import.twisst(weights_file, data_file)
```

* Smooth data and plot

```R
twisst_data_smooth <- smooth.twisst(twisst_data, span_bp = 50000)

plot.twisst(twisst_data_smooth, tree_type="unrooted")
```
We see a clear shift from the 'species' topology, to the 'geography' topology between position 1 and 1.2 Mb. *optix* is found at around 1 Mb, so the signature of introgression is mainly confined to the intergenic region near *optix*.

We used unrooted trees because there is no outgroup in this set of four taxa. This means **we cannot determine the direction of introgression**. To clarify which taxa have shared genes, and in which direction we need to include an outgroup. Fortunately we have an outgroup in the form of two samples from the more distant species *H. numata* ('num').

#### Repeating the analysis with an outgroup

* **Return to the regular terminal** and re-run `twisst` specifying the outgroup.

```bash
python twisst-0.2.1/twisst.py -t heliconius92.chr18.500001-1500000.phyml_bionj.trees.gz \
-w heliconius92.chr18.500001-1500000.phyml_bionj.5pops.weights.tsv \
-g chi -g txn -g ros -g ama -g num --groupsFile twisst-0.2.1/examples/heliconius92.pop.txt --outgroup num
```

Now there are 15 different possible topologies!

* **Back in R** We can plot the new results.

```R
weights_file = "heliconius92.chr18.500001-1500000.phyml_bionj.5pops.weights.tsv"
data_file = "heliconius92.chr18.500001-1500000.phyml_bionj.data.tsv"

twisst_data <- import.twisst(weights_file, data_file)
```

```R
twisst_data_smooth <- smooth.twisst(twisst_data, span_bp = 50000)

plot.twisst(twisst_data_smooth)
```

This shows that two topologies dominate. One is the 'species' topology, in which *cydno* (chi) groups with *timareta* (txn) and the two *melpomene* populations (ros and ama) group together. The other is a specific topology in which txn is found grouped with ama (**nested within the *melpomene* pair**). This tells us that the predominant direction of introgression was from melpomene into timareta. This fits with our current understanding of this group in which timareta expanded down the eastern slopes of the andes, hybridisating with melpomene and acquiring the favoured warning patterns in each region.

___

## Practical Part 3: Topology weighting using Tree Sequence format

So far, we have analysed tree files that have a distinct tree for each chromosome 'block' with a unique genealogy, or for each window in the case of infered trees. This format is somewhat wasteful, because adjacent trees are often extremely similar, usually only different by a single recombination event, which moves one branch from one point in the tree to another.

The Tree Sequence format is efficient because it records not only the connections between nodes on the tree, but the the length of the chromosome for which each connection exists. This format is used by [`msprime`](https://msprime.readthedocs.io/en/stable/), where more information can be found.

#### Simulating a tree sequence

We will use `msprime` to simulate a tree sequence. `msprime` is a coalescent simulator, which means it works by computing the probability that any two individuals share a common ancestor at a given time in the past. In a single population, this is determined by the population size. With multiple populations, this is also affected by the rates of migration between populations, and how long ago the descend from a single ancestral population.

We will use the Python interactive environment for working with `msprime` and the tree sequence, and also to analyse it using a function from `twisst`.

* To ensure that we can import the `twisst` module from within puthon, add it to our python path. (This assumes you already downloaded the `twisst` package in section 1.

```bash
export PYTHONPATH=$PYTHONPATH:twisst-0.2.1
```

* Now **open a Python interactive session** (type 'python'), and also open a script in a text editor, because we are going to modify and rerun some of these lines multiple times

* Import the required modules.

```python
import msprime
import twisst
```

* We will start with a simple simulation of 10 samples from a single population. We have to specify the length and rcombination rate, so msprime will give us a sequence of more than one genealogy, separated by recombination. Here we also specify the random seed just to ensure that in this case we all get the same simulation.

```python
ts = msprime.simulate(sample_size=10,
                      Ne=1000,
                      length=10000,
                      recombination_rate=5e-8,
                      random_seed = 1)
```

* We can check how many distinct genealogies are in the tree sequence.

```python
ts.num_trees
```

* And view them using a nice visualisation method provided with tree sequence objects.

```python
for tree in ts.trees():
    print("interval = ", tree.interval)
    print(tree.draw(format="unicode"))
```
Can you tell what is different between the trees? And what has stayed unchanged?

* If you would like to explore further, the related tool [`tskit`](https://tskit.readthedocs.io/en/latest/index.html) has many inbuilt functions to analyse tree sequences.

* We will now set up a larger simulation with multiple populations. This requires a few different components, which we will define separately. First we define the number of samples and populatin size ofr each population.

```python
pop_n = 10
pop_Ne = 10000

population_configurations = [msprime.PopulationConfiguration(sample_size=pop_n, initial_size=pop_Ne),
                             msprime.PopulationConfiguration(sample_size=pop_n, initial_size=pop_Ne),
                             msprime.PopulationConfiguration(sample_size=pop_n, initial_size=pop_Ne),
                             msprime.PopulationConfiguration(sample_size=pop_n, initial_size=pop_Ne)]
```

* Next we set the migration rates between populations. These can be deinfined in the form of a matrix. Here we set a moderate level of migration in both directions between the second and third populations, and no migration between the others.

```python
migration_matrix = [[0,    0,    0,    0],
                    [0,    0,    1e-4, 0],
                    [0,    1e-4, 0,    0],
                    [0,    0,    0,    0]]
```

* Finally, we set split times, which, in the coalescent world view are joins going backwards in time. In `msprime` these are called mass migrations. So the split between the first two populations is modeled as a mass migration of all individuals from the second into the first population. Further back in time, the third and fourth populations also mass migrate into the first population. Once the first two populations have merged (backwards in time) we also turn off migration between the second and third (because the second technically no longer exists).


```python
t_12 = 1000
t_13 = 5000
t_14 = 10000

demographic_events = [msprime.MassMigration(time=t_12, source=1, destination=0, proportion=1.0), # first merge
                      msprime.MigrationRateChange(time=t_12, rate=0, matrix_index=(2, 1)), # mig stop after merge
                      msprime.MigrationRateChange(time=t_12, rate=0, matrix_index=(1, 2)),
                      msprime.MassMigration(time=t_13, source=2, destination=0, proportion=1.0), #next merge
                      msprime.MassMigration(time=t_14, source=3, destination=0, proportion=1.0)] #final merge
```

* Now we are ready to simulate the tree sequence. We set the length to 50 kb and the recombination rate to 5e-8.

```python
ts = msprime.simulate(population_configurations = population_configurations,
                      migration_matrix = migration_matrix,
                      demographic_events = demographic_events,
                      length = 50000,
                      recombination_rate = 5e-8
                      )
```

* Again we can check the number of trees

```python
ts.num_trees
```
* and, if we dare, we can look at the first tree in the tree sequence:

```python
print(ts.first().draw(format="unicode"))
```

* To compute the topology weights, we define the taxa for twisst. This is easy, because each `msprime` names samples with numbers, so samples 1-10 are the first population, 11-20 the second, etc.

```python
taxa = [list(range(x,x+10)) for x in (0,10,20,30)]
taxonNames = ['A','B','C','D']
```

* And then run a function from `twisst` the computes weightings from a tree sequence

```python
weightsData = twisst.weightTrees(ts.trees(), taxa, treeFormat="ts",
                                 taxonNames = taxonNames, outgroup = 'D', verbose=False)
```

* And now write these weights to a file

twisst.writeWeights("msprime.4pop.weights.tsv", weightsData)


* For plotting we also want to record the start and end positions of each tree on the chromosome. We will extract these ourselves (the intervals are all stored within the Tree Sequence object) and write them to a file.

```python
with open("msprime.4pop.data.tsv", "wt") as dataFile:
    dataFile.write("chrom\tstart\tend\n")
    dataFile.write("\n".join(["\t".join(["chr1", str(tree.interval[0]), str(tree.interval[1])]) for tree in ts.trees()]) + "\n")
```

* Now we **switch back to R** to plot the weights as before. (We could have done this directly in Python but we didn't feel like it.)

```R
weights_file <- 'msprime.4pop.weights.tsv'

window_data_file <- 'msprime.4pop.data.tsv'

twisst_data <- import.twisst(weights_file, window_data_file)
```

* We can make a pretty plot acros the chromosome again

```R
plot.twisst(twisst_data)
```
The profile from this `msprime` simulation looks quite different from that we saw at the beginning, simulated with Hudson's `ms`. There seem to be less dramatic switches in topology weights as we move along the chromosome. Don't ask me why.

* Or, perhaps more interestingly, we can look at the weightings of each topology

```R
plot.twisst.summary.boxplot(twisst_data)
```

Now, what do you expect will happen if we:
Increase or decrease population size?
Make the population split times more or less recent?
Increase or decrease migration rates?

* Once you've made your predictions - test them out!
