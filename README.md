# Classiphy
Classify gene trees regards to the processes generating species-tree gene-tree discordance


## Classiphy analysis overview

The basic steps in a Classiphy analysis involve: 

* Simulation of phylogeny under regimes corresponding to different processes that might contribute to discord;
* Calculation of summary statistics on simulated data to obtain a training dataset that can be used to train a classification function;
* Calculation of summary statistics on empirical data
* Construction of a discriminant analysis function based on principal components extracted from the training data set, and Application of the discriminant analysis function to the original data to classify it with respect to the different processes that might underlie gene tree discord.

The current focus of Classiphy is to distinguish two processes-- Horizontal Gene Transfer (HGT) and Incomplete Lineage Sorting (ILS)

## 1.Simulating gene trees for training data

Classiphy provides a wrapper function `simphy.simu` which uses program **Simphy** for simulating species trees and gene trees. 

```{r eval=FALSE}
simphy.simu(simphy.executable,control.file)
```

User need to provide:

* `simphy.executable`: the path to the executable binary file of Simphy
* `control.file`: Simphy's control file, for details check out its [online user manual](https://github.com/adamallo/SimPhy/wiki/Manual)


Below is an example:
```
-rs 2 //2 species trees
-rl f:50 // 50 locus trees (i.e., gene trees)
-sb f:0.000001 
-sd f:0.0000005 //species birth (sb) and death rate (sd)
-st f:5000000 // total specie tree depth
-sl f:20 // number of species
-si f:1 // 1 individual per species
-sp f:500000 // effective population size
-su f:0.00000005 //substitution rate
-sg f:1 // generation time
-lb f:0 // no gene duplication
-ld f:0 // no gene loss
-lg f:0 // no gene conversion
-lk 1 // Distance-dependent HGT
-gt u:0.00000001,0.00000005 // uniform distributed HGT rate, specifying the min and max
-lt f:gt
-V 1 // screen output option 1
-o test //output folder name
-om 0
-on 0
-ol 0 //turn off some output options so only output trees
```

With this example control file, Simphy will generate a main folder `test` to store the simulation results. In the main folder, there are two subfolders `1` and `2`, corresponding to the two species tree replicates we specified in the first line (`-rs 2`) . Heretoafter we refer these two folders as `repfolder`. In each repfolder, there are

* a species tree file `s_tree.trees`, a one-line file with one tree in Newick format
* a locus tree file `l_trees.trees`,  containing 50 locus trees in Newick format (`-rl f:50`), whose topology might be differ from the species tree due to HGT
* a gene tree file `g_trees.trees`, containing 50 gene trees, whose topology might be differ from its locus tree due to ILS

In principle, user can use any program to do the simulation as long as the simulated trees are organized into three files like this. We only considering the situation where one individual is sequenced for each species (like most of phylogenomic studies), so the tip labels in the locus and gene trees need to be consistent with those of the species tree, and currently, branch length is not used in calculating summary statistics, so the trees can be topology-only.

To simulating training data for a empirical dataset (i.e., a set of gene trees), two parameters should set according to the empirical dataset `-rl`: the number of locus, and `-sl`: the number of species. Simphy can simulate species trees with tree depth (`-st`, can be a range) and a speciation/extinction model(`-sb` and `-sb`), or take an estimated species tree as the starting tree. User can decide (or test) how much information from the empirical dataset to brought in this simulation step.


## 2.Calculating summary statistics on simulated data

The `training.data` function prepares the training data matrix that can be use for training the classification function. In this matrix, four (sets) of summary statistics are included by default: 

* RF distance and MDC based Z score 
* Triplet score 
* branch MDC score
* subtree frequency scores


```{r eval=FALSE}
trainingdata<-training.data(repfolder,phylonet.path)
```

User need to provide:

* `repfolder`: the path to one repfolder. In the example above, it could be `test/1` or `test/2`
* `phylonet.path`: the path to the executable file of PhyloNet (which is used to calculate MDC between gene tree and species trees)

This function will return the training data as a data frame, one row per gene tree. Hence, for multiple repfolders, this function can be run sequentially or parallelly across the repfolders, and `rbind` the result into one big training data matrix. 

User can add in additional summary statistics to the data frame of training data using column bind `cbind`. Make sure the name of additional columns start with "predictor." so that the DAPC function in step 4 will recognize these additional columns as the predictor variables that it should use.

Currently, CLASSIPHY specify the true model for training data using RF distance between locus tree and species tree-- larger-than-zero distance is labelled as HGT gene trees, and zero distance as ILS gene trees (i.e., we do not identify HGT events that result in no topological changes on the tree). This information is stored in the `model.id` column of the training data. If needed, user can change the content of `model.id` (e.g., interested in processes other than HGT), DAPC model will be trained to identify the custom categories in user-specified `model.id` .

If HGT events are rare, majority of the gene trees will be ILS trees, which will make the machine learning algorithm preferentially assign gene trees to ILS. We recommend balancing the training dataset using function `balance.training.data`. 

```{r eval=FALSE}
trainingdata<-balance.training.data(trainingdata)
```
However, balancing is not mandatory. As the DAPC analysis does not take much computation time, researchers could try with both balanced and original training data.

## 3.Calculating summary statistics on empirical data

```{r eval=FALSE}
testingdata<-sumstat(sfile,gfile,phylonet.path)
```

Similar to the `training.data` function, here user need to specify the path to the (estimated) species tree file `sfile`, and the gene tree file `gfile`. (If the species tree and gene tree file are named as s_tree.trees, and g_trees.trees, user can also just specify the folder path as `repfolder=`)

This returns a data matrix with summary statistics calculated from the empirical data, which will be the testing data set for next step. Currently, this function does not work with gene trees with more-than-one sequences per species, and gene trees with missing data.

If user added additional summary statistics to the training data in step 2, the same columns should also be added to the testing data.

## 4.Running DAPC classification

The `dapc` function in R package `adegenet` is used for classification. This package uses a wrapper function `classifyData`  written by Jeet Sukumaran originally for the program [archipelago](https://github.com/jeetsukumaran/archipelago)  

```{r eval=FALSE}
result<-classifyData(target.summary.stats = testingdata,training.summary.stats = trainingdata,n.pca = 'optimize',n.da=NULL,n.pca.optimization.penalty.weight = 0)
```

The returning object `result` is list of two:

* `classification.results` A data frame showing the classification result on the testing data, one row per gene tree: the assigned model (i.e., HGT or ILS), and the posterior weight in each model.
* `trained.model` The trained model, which includes the information that can be use to assess how the model performed based on the posterior prediction of the training data set. 
