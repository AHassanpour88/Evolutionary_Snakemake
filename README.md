<img src="https://github.com/AHassanpour88/Evolutionary_Snakemake/tree/main/images/logo.PNG" align="right" alt="Logo" width="300" height="200">


## Authors
- Azadeh Hassanpour, azadeh.hassanpour@uni-goettingen.de
- Johannes Geibel, johannes.geibel@fli.de
- Torsten Pook, Torsten.pook@wur.nl

## Copyright
Copyright © 2020 – 2024

## License
This program falls under a NonCommercial-NoDerivatives-NoDistribution Public License. By using these scripts, I confirm that I represent an academic institute and will use this script only for research purposes. I explicitly acknowledge the terms in the license agreement [here](https://github.com/AHassanpour88/Evolutionary_Snakemake/blob/main/License.md). I understand that any commercial use requires a commercial license from the owner of the script. For more information about a commercial license, please contact Torsten Pook.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

## Manuscript Availability
A manuscript for the optimization framework is available in [BioarXiv](https://academic.oup.com/g3journal/article/10/6/1915/6026363) / in review. The scripts in this directory are for the example from the toy dairy cattle example and adaptation is needed to use these scripts for other purposes.

# Optimization of breeding program designs through an evolutionary algorithm (MoBPSopti Project)

Modern breeding programs have increased significantly in size and complexity, with various highly interdependent parameters and many contrasting breeding goals. As a result, resource allocation in these programs has become more complex, and deriving an optimal breeding strategy has become increasingly challenging. Optimization in breeding programs is crucial, as it helps us make the most efficient use of limited resources such as time, money, and breeding materials, which increases the likelihood of achieving our objectives and ensure adaptability to changing circumstances. The resource allocation problem is a key challenge in optimization theory with various practical uses and many algorithms have been proposed, all sharing a common structure with two parts: simulation and optimization. First, they use simulation to assess performance under different parameter settings. Then, they use optimization to find the best way to allocate these parameters, achieving different objectives while reducing costs or increasing benefits. In this repository, we use an evolutionary optimization framework to optimize resource allocation in a breeding program with the use of stochastic simulation combining with Snakemake workkflow management system, providing you to fine-tune and enhance your breeding program while efficiently allocating your breeding resources to attain various objectives.

This repository contains all the components you need to for optimization including libraries, system tools, and code, where we store, manage, and deploy our optimization pipeline:
```
|-- Snakefile
|-- config
    |-- config.yaml
    |-- iterations.csv
|-- scripts
    |-- sampleScript.r
    |-- simuScript.r
    |-- evoScript.r
    |-- visualizeParameters.r
|-- Functions
    |-- cost_function.R
    |-- generate_new.R
    |-- termination_criteria.R
|-- environment
    |-- installMoBPSOpti.sh
    |-- environment.yaml
    |-- MoBPS_1.10.40.tar.gz
    |-- RandomFieldsUtils_1.0.6.tar.gz
    |-- miraculix_1.0.5.tar.gz
```
In the sections that follow, we provide detailed explanations for each component, ensuring that you have a clear understanding of their functionalities and how they contribute to our optimization process. If needed, we'll also provide broader insights and introductions to different sections to make your experience with this repository as informative and effective as possible.

## Who can use our framework?
Our framework is not just for one type of breeding program; it works for both animals and plants. While there are subtle differences in the theoretical concepts, objectives, tools, and methods between animal and plant breeding, the overall process of designing and optimizing these programs remains intricate. To use our framework, you need to have a clear plan for your breeding program that can be created using stochastic simulation. A well-designed simulation script serves as the backbone of your optimization, allowing you to systematically explore various aspects of breeding programs and draw meaningful conclusions. When you create your simulation plan, remember one important thing: **Flexibility**. 
When designing a simulation script for a breeding program, it's crucial to build flexibility into the plan to handle diverse scenarios. This flexibility allows you to explore various aspects of breeding programs and adapt to different situations within a specified range. Being able to switch between different strategies allows for the investigation of their impact on genetic progress and breeding program efficiency. Below are key components to consider when specifying the details of your breeding program in the simulation script:

  *  Selection and mating: Implement different selection criteria and methods. This could involve choosing between different indices, adjusting selection intensities, or exploring marker-assisted selection/ genomic selection. Design the script to accommodate varying mating designs (random mating, controlled mating, family-based mating), population sizes, and the presence of different subpopulations. 

  *  Genetic and environmental factors: Allow for the adjustment of heritability values, genetic correlations, and other genetic parameters. Consider incorporating environmental variables that influence trait expression.

  *  Breeding Objectives: Define breeding objectives that reflect different breeding goals. For example, the script should allow the user to prioritize traits differently or adjust the economic weights assigned to each trait.

  *  Scenario-specific Settings: Integrate scenario-specific settings into the script. For instance, the script could handle scenarios like population bottlenecks, introduction of new genetic material, or changes in selection pressure. Design the script to record relevant data at each generation and allow for easy adaptation of analysis methods. 

## How to design your breeding program?
To get started, the foundation of our framework relies on the use of stochastic simulation. That's where Modular Breeding Program Simulator (MoBPS) comes in to streamline various aspects of breeding programs. MoBPS is a simulation tool that helps to closely mimic real-world breeding programs. The software takes care of backend-“stuff“ that you have to account for, but is not main part of analysis (like meiosis, trait simulation, and ensuring that your data is efficiently stored and ready for analysis). Moreover, it comes pre-equipped with functions for common breeding actions such as breeding value estimation, phenotyping, and selection, making your breeding program more efficient. With MoBPS, we can look at different traits in plants and animals, deal with changes in the environment, and balance different goals more precisely. Note that, design and simulation of a breeding program are not transferable between breeding programs due to their program-specific nature, but successful breeding actions and methodologies employed in one program can serve as valuable references for implementation in different programs.

### Introduction to the R package MoBPS

Getting started with MoBPS as a stochastic simulator might feel a bit overwhelming initially because there are many options to choose from. Unlike deterministic simulations, here you have to simulate every individual breeding action individually, which means more programming work and a deeper understanding of the breeding system. To start, we provide a short tutorial (created by Johannes Geibel) which lead you through the basic functionality of the R package MoBPS; <https://www.g3journal.org/content/10/6/1915>) which can be found at [GitHub](https://github.com/tpook92/MoBPS "https://github.com/tpook92/MoBPS"). Further help may give the [MoBPS YouTube channel](https://www.youtube.com/channel/UC4LDcBka39NidOF1y_65FFw "https://www.youtube.com/channel/UC4LDcBka39NidOF1y_65FFw"). To explore additional examples and gain a deeper understanding of practical applications, interested readers are encouraged to review materials from previous workshops on MoBPS available at the following GitHub repository: [GitHub](https://github.com/AHassanpour88/MoBPS/blob/master/MoBPS_Workshop_GFT%2526KWS). 

<details><summary><b>Show instructions</b></summary>

**1. Prerequisites:** 

Before you can use MoBPS, you need to install the package on your computer. It is also helpful to install the packages `RandomFieldsUtils` and `miraculix`, as they can speed up calculations. In order to compile the packages from source, **you may need to install [rtools](https://cran.r-project.org/bin/windows/Rtools/) on your Windows computer** first.
Following code installs some packages which are needed packages for the installation from GitHub.

```{r install_packages, echo=TRUE}
pck <- c('devtools','stringr')
install.packages(pck[!pck %in% rownames(installed.packages())],
                 repos="https://cloud.r-project.org")
invisible(lapply(pck,library,character.only = TRUE, quietly = TRUE))
```
We can now install working `RandomFieldsUtils` and `miraculix` versions.
> Note that following code installs the currently supported versions from the MoBPS GitHub page, which may deviate from the latest releases.
> 
```{r install_dependencies, eval=FALSE, exercise=FALSE, include=TRUE}
pck <- character(2)
temp <- readLines('https://github.com/tpook92/MoBPS')
pck[1] <- str_extract(grep('RandomFieldsUtils_[0-9.]+tar.gz',temp, value = TRUE),
                      'RandomFieldsUtils_[0-9.]+tar.gz')[[1]]
pck[2] <- str_extract(grep('miraculix_[0-9.]+tar.gz', temp, value = TRUE),
                      'miraculix_[0-9.]+tar.gz')[[1]]
pck <- paste0('https://github.com/tpook92/MoBPS/raw/master/', pck)
# the installation itselves is commented out to not be accidentially run!
# install.packages(pck, repos = NULL, type = 'source')
```
> **WARNING!** MoBPS uses `miraculix` by default to speed up your code and reduce the memory need, if installed. However, the underlying coding is system dependent and cross platform shipment of the population list object will change the data! To be on the safe side, you can set the option `miraculix = FALSE` in `creating.diploid()`. To check, whether an object was created under usage of `miraculix`, use `object_name$info$miraculix` .

Finally, we install the current `MoBPS` version from GitHub.

```{r install_MoBPS, echo=TRUE, exercise=FALSE}
#devtools::install_github('tpook92/MoBPS/pkg')
library(MoBPS)
```

If you want to use preconfigured SNP maps, also install `MoBPSmaps`

```{r install_MoBPSmaps, eval=FALSE, exercise=FALSE, include=TRUE}
temp <- readLines('https://github.com/tpook92/MoBPS')
pck <- str_extract(grep('MoBPSmaps_[0-9.]+tar.gz',temp, value = TRUE),
                   'MoBPSmaps_[0-9.]+tar.gz')[[1]]
pc <- paste0('https://github.com/tpook92/MoBPS/raw/master/', pck)
# the installation itselves is commented out to not be accidentially run!
# install.packages(pck, repos = NULL, type = 'source')
```

The R package MoBPS is mainly based on two functions: `creating.diploid()` for the creation of random animals and `breeding.diploid()` to mimic breeding actions (e.g. mating and selection of groups). The following video will introduce you to the basic concept, but the content is also captured by this tutorial: <https://www.youtube.com/watch?v=xAl51woiU5s>

**2. Creation of a base population:**

> A short remark on the notation: The example population list will be denoted as `popList`. When just showing exemplary usage without the use for the example program, the results will be stored as `temp` and overwritten at the next example. The exercise population shall be denoted as `ownPop` to stay consistent.

**2.1. Base population without phenotype:**

We will start by simulating a population of 50 individuals and 500 SNPs. The basic function to do so is `creating.diploid()`, which creates the initial population list of founder animals. There also exists a tailored `summary()` function, printing the most necessary information.

```{r show_base_creation, echo=TRUE}
popList <- creating.diploid(nindi = 50, # number of individuals
                            nsnp = 500) # number of SNP
summary(popList)
```

As you noticed, MoBPS automatically divided the population into 50% male and 50% female individuals and assigned them to separate cohorts (a group of individuals which share properties at a special time point). So the population consists of two cohorts.

If you want to change this (e.g. only males or more females than males), add the parameter `sex.quota` which specifies the share of female individuals. So following code creates the same population as before, just all individuals being male. You can additionally name your cohort by the parameter `name.cohort`.

```{r show_base_creation_2, echo=TRUE}
popList <- MoBPS::creating.diploid(nindi =  50,
                                   sex.quota = 0, # share of females
                                   name.cohort = 'founder_male', # name of the cohort
                                   nsnp = 500)
summary(popList)
```

If you want to add additional founder cohorts, use `creating.diploid()` on the existing population (specified via `population`). It is not necessary to also specify the genome info again - will be used from first spacification.

```{r}
popList <- MoBPS::creating.diploid(population = popList, # existing population to which the new individuals will be added
                                   nindi =  50,
                                   sex.quota = 1,
                                   name.cohort = 'founder_female')
summary(popList)
```

If you want to extract the names of all cohorts in the population, use `get.cohorts()`. If you set `extended = TRUE`, you will get a matrix with more information instead of a simple vector.

```{r}
get.cohorts(popList)
```

```{r}
get.cohorts(popList, extended = TRUE)
```

By now, all markers are located on one single chromosome with a default length of 5 Morgan (M) and a physical size of 5 kilo bases (kb). The number of chromosomes can be either changed by the parameter `chr.nr` which leads to the according number of 5M chromosomes with markers evenly distributed on chromosomes, or by directly supplying an vector with the intended lengths in M via `chromosome.length` (the number of markers per chromosome has to be specificed in this case).

```{r}
temp <- creating.diploid(nindi =  50,
                         sex.quota = 0,
                         nsnp = 500,
                         chr.nr = 2) # number of 5M chromosomes 
summary(temp)

temp <- creating.diploid(nindi =  50,
                         sex.quota = 0, 
                         nsnp = c(300,100,100),
                         chromosome.length = c(2,1,0.5)) # chromosome length [M]
summary(temp)
```

It is also possible to create a base population from real data by supplying e.g. a VCF file.

```{r eval=FALSE, include=TRUE}
creating.diploid(vcf = "/path/to/file.vcf")
```
**2.2. Traits:**

**2.2.1. Adding traits to the base population:** 

By now, our base population consists only of individuals with genomes. To simulate breeding programs, we need to additionally simulate phenotypic traits and the according genetic basis. We will start doing so by randomly drawing additive effects from a normal distribution with additive variance of $\sigma_{A;MKG}^2=25^2kg^2$ for 100 SNP which contribute to the trait milk yield [kg] (MKG) via the function `creating.trait()`. Note that it is also possible to directly simulate the trait in `creating.diploid()`.

```{r}
temp <- creating.trait(population = popList,
                          # name of trait
                          trait.name = 'MKG',
                          # number of additive effects
                          n.additive = 100,
                          # additive genetic target variance
                          var.target =  25^2)
summary(temp)
```

```{r}
get.cohorts(temp)
```

```{r echo=FALSE, include=FALSE, eval=FALSE}

hist(get.bv(temp, gen = 1), freq = FALSE,
     main="", xlab= "true breeding value")


lines(density(get.bv(temp, gen = 1)))
lines(seq(0,200,0.1), dnorm(seq(0,200,0.1), mean = 100, sd = 25), col = 'red')
legend('topright',
       legend = c('simulated', 'expected'),
       col = c('black', 'red'),
       lty = 1)
```

It is also possible to use multiple correlated traits. Only additional specification of residual and genetic correlation is needed. We specify a second trait fat percentage (FP) with additive variance ($r_{A;FP}$) is $0.5\%^2$. Our additive genetic correlation ($r_{A;MKG;FP}$) to be 0.5 and the environments are assumed to be uncorrelated for now ($r_E=0$).

```{r}
popList <- creating.trait(population = popList,
                          trait.name = c('MGK','FP'), 
                          n.additive = c(100,200),
                          # target variance
                          var.target = c(MKG = 25^2, FP = 0.5^2),
                          # indices of correlated traits
                          shuffle.traits = c(1,2),
                          # genetic correlation (0.25) matrix
                          trait.cor = matrix(c(1,.25,
                                                 .25,1),2,2),
                          # residual correlation (0) matrix
                          new.residual.correlation = matrix(c(1,0,
                                                            0,1),2,2))
summary(popList)
```

Be aware that using a genetic and residual correlation matrix which results in not defined phenotypic correlations, will throw an error message. Note that such a correlation cannot occur in practice even though pair-wise correlation estimations might lead to such estimates:

```{r error=TRUE}
temp <- creating.trait(population = popList,
                       shuffle.traits = c(1,2),
                       shuffle.cor = matrix(c(1,.5,
                                              .5,1),2,2),
                       new.phenotype.correlation = matrix(c(1,-1,
                                                            -1,1),2,2))
```
**3. Breeding Actions:**

**3.1. Generation of a population structure (random mating):**  

The individuals in the founder population are unrelated. This is normally not the case in practice. Simulations therefore start with a number of burn-in generations which simply apply random mating. The function for all breeding actions is `breeding.diploid()`. It generates new individuals by mating old ones and thereby simulating recombination and mutation events. If not specified differently, the last available males and females are chosen for random mating.

```{r}
temp <- breeding.diploid(population = popList,
                         # number of newly created individuals
                         breeding.size = 100,
                         verbose = TRUE,
                         display.progress = FALSE) 
summary(temp)
```

When applying the function `get.cohorts()` on the updated population, we see that a male and a female cohort, each of size 50, were created in generation 2.

```{r}
get.cohorts(temp, extended = TRUE)

```

As we do not only want to run it for one burn in generation, but for five, we use a loop to save code lines and additionally use a more useful cohort name.

```{r}
n_burnin <- 5
for(i_burnin in 1:n_burnin){
  popList <- breeding.diploid(popList,
                           breeding.size = 100,
                           name.cohort = paste0('burnin_',i_burnin),
                           verbose = FALSE)
}

summary(popList)
```

The population now exists of 600 individuals in 6 generations (1 founder and 5 burnin) and 12 cohorts. MoBPS automatically added "\_M" and "\_F" to the cohort names.

```{r}
tail(get.cohorts(popList, extended = TRUE))
```
**3.1.1. Automate burnin phase:**  
Note that you could also use the function `founder.simulation()` to automate the burnin phase.

The distribution of the allefrequencies is sampled from a beta distribution:

```{r}
par(mfrow = c(1,2))
hist(rbeta(1e4,1,1),
     xlab = 'allele frequency',
     main = 'shape1 = 1, shape2 = 1')
hist(rbeta(1e4,0.5,1),
     xlab = 'allele frequency',
     main = 'shape1 = 0.5, shape2 = 1')
par(mfrow = c(1,1))
```

```{r}
# this is only to show how it works, but not enougth generations e.t.c
# for a realistic set up!!!
temp <- MoBPS::founder.simulation(nindi = 100,
                                  sex.quota = 0.5,
                                  nsnp = 1e3,
                                  beta.shape1 = 0.5,
                                  beta.shape2 = 1,
                                  n.gen = 50,
                                  verbose = FALSE
                                  ) 
```

`founder.simulation()` prints information on the effective population size as well as fixed markers, allele frequencies and LD decay. This allows to fine tune the generation process based on ones needs. The result then is simply a haplotype matrix:

```{r}
str(temp,max.level = 1)
temp[1:5,1:5]
```

Which can be used as founder population in `creating.diploid()`:

```{r}
temp <- creating.diploid(dataset = temp)
summary(temp)
```

If you want to keep more information from the generation of your base population, consider setting `big.output = TRUE`

```{r}
temp <- MoBPS::founder.simulation(nindi = 100,
                                  sex.quota = 0.5,
                                  nsnp = 100,
                                  n.gen = 2,
                                  big.output = TRUE,
                                  verbose = FALSE,
                                  plot = FALSE
                                  ) 
str(temp,max.level = 1)
# haplotype matrix: temp[[1]]
# map: temp[[2]]
# population list: temp[[3]]
# A matrix: temp[[4]]
```
**3.2. Selection:**  

The basic principle of breeding is that we do not use all individuals, but only a selected subset for matings. Let's start with random selecting 5 males from the last generation and assigning an own cohort name to them.

```{r eval=FALSE}
temp <- breeding.diploid(popList,
                         selection.size = 5, # selected number
                         selection.criteria = 'random', # criterium
                         selection.m.cohorts = 'burnin_5_M', # from cohort
                         
                         name.cohort = 'selected_male',
                         # generation number has still to be 6
                         add.gen = n_burnin + 1,
                         # copying individuals instead of creating new ones
                         copy.individual.m = TRUE)
summary(temp)
tail(get.cohorts(temp,extended = TRUE))
```

If we do not only want to apply random selection, but selecting the best males based on estimated breeding values (EBV), we need to additionally simulate phenotypes for our individuals. This needs the additional specification of the environmental variances ($\sigma^2_E$) of the traits. Remember Following connection between heritability ($h^2$), $\sigma^2_A$ and $\sigma^2_E$:

$$\sigma^2_E=\frac{\sigma^2_A}{h^2} - \sigma^2_A = \frac{\sigma^2_A (1-h^2)}{h^2}$$

If we assume heritabilities of 0.8 and 0.3 respectively for our traits, $\sigma^2_{E;MKG}=156.25$ and $\sigma^2_{E;FP}=0.58$. As only cows can give milk, we simulate the phenotypes for all not yet phenotyped females by the following code. Note that you have to set `sigma.e` only once, if you do not assume changes over the generations! `sigma.e` can also be automatically fitted by the use of the `heritability` parameter.

```{r eval=FALSE, include=FALSE}
(25^2*(1-0.8))/0.8
(0.5^2*(1-0.3))/0.3
```

The individuals to phenotype can be specified using the parameter group `phenotyping` :
| Parameter              | Usage                                                                                                                                      |
| ---------------------- | ------------------------------------------------------------------------------------------------------------------------------------------ |
| `phenotyping`          | Quick access to phenotyping for all: "all", non-phenotyped: "non_obs", non-phenotyped male: "non_obs_m", non-phenotyped female: "non_obs_f" |
| `phenotyping.gen`      | Vector of generation for which to generate additional phenotypes                                                                           |
| `phenotyping.cohorts` | Vector of cohorts for which to generate additional phenotypes                                                                              |
| `phenotyping.database` | Matrix of groups for which to generate additional phenotypes                                                                               |


As milk traits can be recorded only in females, we add phenotypes for all females, which were not phenotyped yet:

```{r phenotyping}
popList <- breeding.diploid(popList,
                            phenotyping = "non_obs_f",
                            # heritability = c(0.8,0.3) ; var.a = c(25^5,0.5^2)
                            sigma.e = c(156.25, 0.58))

# alternatively, you could also directly specify the heritability
temp <- breeding.diploid(popList, heritability = c(0.8,0.3), sigma.e.gen = 1)
temp$info$last.sigma.e.value
```

Following code selects 10 out of 50 males from the last cohort. Breeding values for the males are estimated from female phenotypes via genomic breeding value estimation (default).

```{r eval=TRUE}
temp <- breeding.diploid(popList,
                         selection.size = 10,
                         # selection based on EBVs
                         selection.criteria = 'bve',
                         selection.m.cohorts = 'burnin_5_M',
                         selection.m.gen = n_burnin + 1,
                         # estimate BV
                         bve = TRUE,
                         # use all by now available generations
                         bve.gen = c(1:(1+n_burnin)),
                         
                         name.cohort = 'selected_male',
                         add.gen = n_burnin + 1,
                         copy.individual.m = TRUE)
summary(temp)
```

**3.3. Mating and selection**

Of course, we want to use these selected males as sires for a new generation. To do so, we can perform selection and mating at once and again use a loop to repeat this process for a number of generations. As we did not specify which animals to use for reproduction, the last generation will automatically be used - if something else is intended use selection.m.gen/database/cohorts and selection.f.gen/database/cohorts for choose potential paternal/maternal sires:

```{r eval=TRUE}
n_burnin <- 5
n_selection <- 5 # number of selecting generations
for(i_selection in 1:n_selection){
  popList <- breeding.diploid(popList,
                              # create 100 calves each year
                              breeding.size = 100,
                              # select 5 bulls (10%) and all cows
                              selection.size = c(5,50),
                              # select based on EBVs
                              selection.criteria = 'bve',
                              # estimate breeding values
                              bve = TRUE,
                              # use all previous generations 
                              bve.gen = c(1:(n_burnin + i_selection)),
                              # add phenotypes for new females, if not yet present
                              phenotyping = "non_obs_f",
                         
                              name.cohort = paste0('offspring_',i_selection),
                              verbose = FALSE)
}
summary(popList)
```

> Note: Depending on the size of your population, the number of SNPs and the chosen method, the time needed for breeding value estimation is the big problem of most simulations!


**4. Evaluating results**

`MoBPS` contains some utility functions which directly produce simple plots or extract information from a population list object. You already know about `summary` and `get.cohorts`, which you can use to keep track on the simulation. In post simulation evaluation, you are most likely interested in the development of the true breeding values over time. This can be achieved via `bv.development()`.

```{r eval=TRUE}
bv.development(population = popList,
               gen = 1:11, # generations
               development = 1, # plot only true breeding values
               display.cohort.name = TRUE)
```

> Note that by default 100 is added to the breeding values.

As the average breeding value of MKG reaches more than 150, we realize a breeding gain of \> 50 kg milk which is $> 2\sigma_{A;MKG}$. For FP, we reach less than 1% which is $< 2\sigma_{A;FP}$.

We can also use boxplots via `bv.development.box()` either by generation or for selected cohorts. We additionally suppress the return value of this function by `invisible()`.

```{r}
invisible(bv.development.box(popList,
                             gen = 1:11))

invisible(bv.development.box(popList,
                             cohorts = get.cohorts(popList)))

```

The accuracy of the breeding value estimation for a cohort ($r_{BV;EBV}$) and the according variances will be returned by `analyze.bv()`. The first list elements shows the correlations and the second one the additive variances.

```{r warning=FALSE}
analyze.bv(popList,
           cohorts = 'burnin_5_M') # first male generation selection was performed for
```

The higher accuracy for MKG due to a higher heritability explains the stronger gain for MKG. Further note that we did not phenotype on the male side, we can therefore not get any correlation to the phenotypes.

Selection always comes along with inbreeding. `kinship.development()` calculates average between and within individual kinship coefficients by generation:

```{r}
 temp <- kinship.development(popList,
                    gen=1:11,
                    display.cohort.name = TRUE)
```

Note that the output of `kinship.development()` is a matrix which contains the average between-individual kinship coefficients in the first column and the within-individual kinship coefficients in the second column.

```{r}
head(temp)
```

As $2C_{xy}=R_{xy}$ and $2C_{xx}-1=F_X$ , we can use the results to calculate average relationship and inbreeding:

```{r}
# average between individual relationship coefficient
(Rxy <- round(2*temp[,1],3))
```

```{r}
# average inbreeding coefficient
(Fx <- round(2*temp[,2] - 1, 3))
```
> Note that the burn in phase (random mating) led to only minor inbreeding, while the selection phase increased it much more. You should also see, that there are random fluctuations due to the small population size.

</details>

## Our framework as Snakemake Workflow Management System
Snakemake is being used as a workflow management system for creating and executing data analysis pipelines to automate a set of scripts that need to be run in a specific order regularly. It helps to ensure that all of the steps in the process are executed correctly and efficiently every time from start to finish. Snakemake is a powerful tool that allows us to define a workflow in terms of rules, where each rule defines the input and output files and the commands to execute and transform the input files into the output files. One of the key features of Snakemake is that it handles dependencies between the rules automatically. It tracks the status of input and output files and only executes a rule when its input files are newer than its output files. This ensures that the workflow is executed correctly and efficiently, as Snakemake only re-runs necessary rules. This can save a lot of time and resources when dealing with large and complex workflows. Moreover, Snakemake allows us to customize the workflow to integrate it well with other systems and tools. For example, we could specify custom shell commands to run, use different types of cluster resources, or invoke other software programs. Snakemake can execute tasks in parallel on multi-core machines or distribute tasks across a cluster or cloud resources. Snakemake also provides detailed logging and reporting functionality, so we could easily track the progress of our workflow and troubleshoot any issues that may arise. We have developed a Snakemake script available in the [Snakemake](https://github.com/AHassanpour88/Evolutionary_Snakemake/edit/main/Snakefile) specifically for optimizing breeding programs. Snakemake reads various input files from our framework, which we will explain in more detail to help you configure it correctly. Before diving into our Snakemake rules and the vital components of our pipeline, you'll find a brief guide below on how to install Miniconda and set up our optimization framework.

>  **WARNING!** It is important to note that the Snakemake script should not be altered or modified by the user.

### Snakemake environment for optimization

<details><summary><b>Show instructions</b></summary>
 
**1. Get miniconda**

For easy installations of Snakemake, we will further need to install miniconda. Miniconda is a leightweight version of the Conda package management system. The following code will get you the installer and installs you the latest miniconda version with python 3.9. Answer all questions with yes (y). When you need to scroll down for the terms of usage in more, you can do so by hitting return, but be carefull, as this will also answer the acceptance with no, if you are too fast, and you need to repeat the process then.
If you are working on MacOS or do want to have an older python version than 3.9 in your base environment, find your appropriate installer at https://docs.conda.io/en/latest/miniconda.html

> Note that the `.sh` script may be named differently in your case!!

```{sh}
cd
wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.9.2-Linux-x86_64.sh
ls
bash Miniconda3-latest-Linux-x86_64.sh # may to be replaced by another filename
```
If you installed miniconda correctly, you should have a folder miniconda3 in your home directory. Further, there exists a hidden file (hidden files start with `.` and are only shown when you use the option `-a` with `ls`) named  `~/.bashrc`. The bashrc contains code which is executed on system startup. The installation writes some conda specific content to the end of it:

You can close `less` by typing `q`!

You should now trigger the execution of the `.bashrc` by

```{sh}
source ~/.bashrc
```
Then, conda environments can be started with `conda activate` and an activated conda environment is indicated at your prompt.

As long as a conda environment is activated, programs will be used from this environment before searching for global program installations.

Further download the Anaconda cheat sheet to learn about conda usage:
https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf

**2. Installing conda packages**

Conda generally eases the process of installing packages, but it may become slow when setting up larger environments. Therefore first install [mamba](https://github.com/mamba-org/mamba) (C++ re-implementation of the conda installation tools and really faster - couldn’t believe by myself first) via conda:

```{sh}
conda install -n base -c conda-forge mamba
# -n --> name of environment to install to
# -c --> name of channel to install from
```

Installation of conda packages should be done from now on by `mamba` (works without activating the base environment)
```{sh}
mamba install -n <env_name> -c <channel_name> <pck1> <pck2>
```
You can further create additional environments, if package incompatibilities exist:
```{sh}
mamba create -n <env_name> -c <channel_name> <pck1> <pck2>
```
Alternatively, creation of environments based on a environment file in `.yml` format is helpful
```{sh}
name:
  mobpsopti
channels:
  - conda-forge
- bioconda
- defaults
dependencies:
  - python
- pandas
- bioconda:snakemake>=7.14
- r-base>=4.1
- r-yaml
- r-tidyverse
- r-data.table
- r-remotes
- r-stringr
- r-patchwork
- r-purrr
```
The installation then works as following
```{sh}
mamba env create -f=/path/to/environment.yml
# or:
mamba env create -f=/path/to/environment.yml -n <newName> # overwrites default name
```
The new environment can be activated by:
```{sh}
conda activate mobpsopti
```
**3. mobpsopti environment for optimization**

Now that you are familiar how to create special environment, to make our optimization work easier for you, we've created a bash script that helps you install `Miniconda` and a special environment called `mobpsopti`. You can find all necessary file for this part in [environment folder](https://github.com/AHassanpour88/Evolutionary_Snakemake/tree/main/environment). This `mobpsopti` environment is like a toolbox with different tools we need for optimization. It includes things like Python, Pandas, and Snakemake for managing tasks, as well as R and some R packages for data analysis. All of these tools are important for our optimization work and help us get things done more efficiently. Before you start your optimization please make sure that you installed everything correctly. Further, if you plan to distribute the snakemake jobs across a larger [cluster](https://snakemake.readthedocs.io/en/stable/executing/cluster.html), you need to setup a profile for your specific cluster scheduler (e.g. slurm). You can find some default profiles in this [project](https://github.com/snakemake-profiles/doc). 

</details>

### Snakemake rules

<img src="https://github.com/AHassanpour88/Evolutionary_Snakemake/blob/main/images/Snakemake.png" align="right" alt="Snakemake" width="300">

Our Snakemake process is built around four rules, which are translated into four distinct R scripts. For each of these rules, there is a corresponding R script in [script folder](https://github.com/AHassanpour88/Evolutionary_Snakemake/tree/main/scripts). Each rule describes a step in an analysis defining how to obtain output files from input files. Dependencies between rules are determined automatically. Snakemake starts with the first rule and it checks if all the inputs required for that rule are available. If the inputs are ready, it runs the corresponding R script for that rule. After completion, it moves on to the next rule in the sequence. This process repeats until all four rules have been executed. All 4 scripts share a common feature. They uniformly utilize two essential lines of code:

```{r}
args <- commandArgs(TRUE)
```
This line uses the commandArgs() function to capture the command-line arguments passed to the R script. The TRUE argument indicates that we want to include all arguments, including the script name. The result is stored in the variable args, which is now a vector containing the script name, parameter settings and their values for simulation/optimization, random seed, number of replication and any additional arguments we need in each step.

```{r}
output <- args[1]
```
This line selects the first element of the args vector (all used variables) and assigns it to the variable output that can be used as input for the next step.

Let's break down each rule and see how Snakemake uses them one by one for optmization. After formulating your breeding program and the parameters you want to optimize into an optimization problem, our framework starts. Here you see a visualization of the Snakemake workflow for our evolutionary optimization model. Shown are the rule names defined and their input-output relationships. 


**Rule 1: sampleVAl (sampleScript):**
```{sh}
checkpoint sampleVAl:
    """
    input:
        sampleScript = config["sampleScript"]
    output:
        "results/parameters/params_1.csv" 
    """
```
In this rule, the initial script [sampleScript.r](https://github.com/AHassanpour88/Evolutionary_Snakemake/edit/main/scripts/sampleScript.r) will be executed. If the rule successfully completes, it will result in the creation of a **results** folder. Inside this folder, a **parameters** directory will be generated, and the initialization parameters will be stored in a file named `params_1.csv`. Initialization step allows users to define parameter bounds and constraints, such as a budget. The purpose of this script is to generate parameter settings for a simulation script. User should modify this script based on the number of parameters that he/she wants to optimize. In our pipeline we allow two specific types of parameter settings for optimization, which can be varied based on different variables or actions. 

* Continuous variables: represented by generating random samples from a uniform distribution within the specified range of the continuous variable. Continuous parameters may include the number of candidates selected, genotyped, mated, or any other continuous variables that could impact the program's success. By varying these parameters and analyzing their effects on the objective of the breeding program, decision-makers have insights on how to optimize the program's outcomes.
* Class variable: represented by generating random samples with n potential values from a categorical distribution (also called a generalized Bernoulli distribution), where each setting has an equal probability. The class type of action involves making specific decisions throughout the breeding program.  It is crucial to identify when and where key decisions should be made to minimize the impact of unavoidable extensions in terms of time and cost during the implementation of various strategies at different stages.

>  **NOTE!** The script should be customized by the user to accommodate the desired breeding program for optimization.

**Rule 2: runSimulations:**
```{sh}
rule runSimulations:
    input:
        parameters = "results/parameters/params_{iter}.csv",
        script = ancient(config["simuScript"])
    output:
        "results/iteration_{iter}/results_{sim}_{rep}.RData"
```
The **runSimulations** rule requires input data in two main forms. First, it necessitates input parameters retrieved from a file located at `results/parameters/params_{iter}.csv,` where the placeholder {iter} corresponds to the specific iteration being processed.  The input data for this step is automatically generated in the preceding step. The second component of the rule involves [simuScript.r](https://github.com/AHassanpour88/Evolutionary_Snakemake/edit/main/scripts/simuScript.r) execution. The simulation model is used to find the values of the decision variables (i.e. model inputs) that optimize (maximize or minimize) the breeding objectives (i.e. corresponding values of model outputs) and ensure
that any constraints are satisfied. As the simulations are executed, the rule generates an output file that places in the **results** directory as `results/iteration_{iter}/results_{sim}_{rep}.RData`. In this path, {iter} represents the iteration number, while {sim} and {rep} act as placeholders for values that are specific to the simulation and the number of replication that a simulaion shuld be performed being performed. The simulation script serves as the backbone of our pipeline, and it requires careful design tailored to your unique needs and specific breeding programs. It's essential to create this script to be both flexible and adaptable, capable of accommodating various parameter settings that you'll use for optimization. The script should be versatile enough to easily handle these parameter adjustments, which allows you to fine-tune and optimize different aspects of your simulations and efficiently explore a wide range of scenarios.

**Burn-in Phase purpose (input file for simuScript outside of framework):**

To make a complex breeding program work faster, we suggest setting up a special phase called the burn-in separately. During this phase to speed up the process and save time, we simulate initial breeding cycles in the burn-in phase to create a basic population structure. We do these calculations in advance and store them as these computations will be the same accross all simulations and instead of redoing the simulation each time, we load in this precomputed population list after the burn-in phase. You only need to do this burn-in phase once and this file becomes the starting point for your simulation, making sure all your optimization efforts start from the same place for consistency.

> **NOTE!** If your breeding program is not computationaly too expensive to simulate, you can merge the burn-in phase with the simuScript step.

**simuScript (input file for evoScript):**

Once the parameter settings are generated using the sample script, they are passed to the simulation script. In the simulation script, users can design their breeding program and define the **objective function** to be optimized. simuScript is nearly identical to the BASF_burnin_phase.R script, differing only in the starting/ending generation and the inclusion of a target function. 

> **NOTE!** If you've already completed a burn-in phase, please ensure that your simulation script correctly reads and imports the outputs from the burn_in phase script at the begining of your simuScript.

The simulation script generates a population of candidate solutions based on the provided parameters. Our toy example breeding program does not require a burn-in phase. 

>  **NOTE!** The script should be customized by the user to accommodate the desired breeding program for optimization.

**Rule 3: runEvolutionary (evoScript: input file for simuScript till termination criteria is achieved):**
```{sh}
checkpoint runEvolutionary:
    input:
        files = "results/iteration_{iter}/results_{sim}_{rep}.RData",
        script = ancient(config["evoScript"])
    output:
        nextpars = "results/parameters/params_{iter}.csv"
```
The [evoScript.r](https://github.com/AHassanpour88/Evolutionary_Snakemake/edit/main/scripts/evoScript.r) relies on the results derived from the simulation script as its inputs, along with their respective parameter configurations. Upon successful execution, the rule generates a fresh set of parameter settings stored at `results/parameters/params_{iter}.csv`, where the variable {iter} signifies the current iteration number. This step involves evaluating the objective function for each candidate solution and applying various operators iteratively, and the algorithm generates new parameter settings for the next iteration as inputs for the simulation script  by selecting and modifying the best parameters settings. This process continues until the algorithm converges towards an optimized solution for the breeding program. The evolutionary script also provides access to the data file generated in each iteration, enabling users to monitor and evaluate the algorithm's progress. The evolutionary script generates a folder called **EvoStatus** starting from iteration 2. Within this folder, an **EvoStatus.RData** file is updated in each iteration. This file contains the current status of the evolutionary algorithm, including relevant data and variables for monitoring and analysis purposes. By saving the progress in separate files for each iteration, users can track the evolution of the algorithm and access historical data for further analysis and evaluation. This data facilitates informed decisions on parameter settings and termination criteria, allowing adjustments during the algorithm's development and implementation. 

>  **WARNING!** Modifying this script is not advisable and should be avoided by the user.

**Rule 4: gatherResults (visualizeParameters):**
```{sh}
rule gatherResults:
    input:
        files = aggregate_input_final,
        script = config["visualizeParameters"]
    output:
        "rplot.pdf" 
```
The [visualizeParameters.r](https://github.com/AHassanpour88/Evolutionary_Snakemake/edit/main/scripts/visualizeParameters.r) is the final component of the pipeline, responsible for demonstrating the convergence of parameters over iterations. This script takes the final output generated from all iterations and produces visualizations that aid in understanding the level of convergence achieved by the parameters. When the termination criteria is achieved, this script create one plot for each parameter and validate if they converge to an optimal value. The plot will be saved as `rplot.pdf` in the main folder at the end. 

>  **WARNING!** Modifying this script is not advisable and should be avoided by the user.

### Costumized parameter settings for different rules in Snakemake

**Config [folder](https://github.com/AHassanpour88/Evolutionary_Snakemake/tree/main/config)**

**config.yaml**
   
Start designing and optimizing your new breeding program by utilizing the `config.yaml` file. This crucial file enables you to configure the breeding pipeline according to the specific requirements of your program. Specify the number of parameters, differentiate between class and continuous variables, and provide appropriate variable names to optimize your breeding program effectively. 

Plugins for user within config file:

* `@nfactors` Number of parameter that you want to optimize in your simulation process.
* `@name_parameter` Names of the parameters you want to optimize are listed here. These names should be exactly as the names of the parameter you want to optimize in the simulation process. Each parameter is assigned a specific index based on their position in this variable. Please note that the indexing starts from 1 and continues sequentially. Users must accurately refer to the correct indices.
* `@binary_parameter` Existence of class variables by TRUE/FALSE
* `@number_binary_parameter` Number of class variables
* `@Sim_init` Number of initial simulations to perform in first iteration to have a good coverage of the search space
* `@time` Define the time constraint for single simulation
* `@memory` Define the memory constraint for single simulation
> To calculate the approximate time and memory for a single simulation, perform only one simulation (with upper bound of the desired parameter settings) and use resource tracking with **sacct** for this step to provide insights into job execution and resource consumption.
* `@min_range_mu` The mutation operator need a minimum value range for each parameter. This range is essential in cases where the standard deviation of a parameter being mutated is small or when the parameter has a limited set of possible values
* `@mut_offspring` Mutation probability in recombined offspring.
* `@mut_parent` Mutation probability in a single parent. 
* `@Hybrid_Breeding`or `@Line_Breeding` by TRUE/FALSE 
* `@linked_parameter` Designed when optimizing closely linked parameters by TRUE/FALSE. When one parameter is linked to other parameters, a mutation in that parameter can affect the others. This option ensures that any mutation in a linked parameter causes changes in all other linked parameters. Our approach is here to reduce the number of mutations that occur in highly linked parameters, which can otherwise lead to an unstable or inefficient optimization process. 

> **NOTE!** `linked_parameter_female` and `linked_parameter_male` are two plugins and a way to organize groups of parameters, where the first value in the list is considered the main parameter. This main parameter is linked with other parameters listed afterward. For instance, if you have `linked_parameter_female: [1,2,3]`, it means that parameters 2 and 3 are linked to parameter 1. So, any changes or adjustments made to parameter 1 will also affect parameters 2 and 3. It's a way of indicating relationships between parameters and understanding that they work together as a group. The same approach is applied in `@Line_Breeding`

> **NOTE!** The linked_parameter option can be challenging to apply universally across different breeding programs. It's advisable to set `linked_parameter_adjustment` to TRUE in config file and customize the evolutionary script in `linked_parameter_adjustment` part to align with your program's specific requirements.

* `@base_cost_ini` Initial number of each parameter in a base scenario that the cost needs to be calculated for that
* `@cost_par` Only index of parameters that are included in cost functions

> **NOTE!** Class variables are not directly included in the cost and they usually depend on the existance of the other cohorts 

* `@cheapest_unit_female` Identifying the index of cheapest unit in the breeding program on female side
* `@cheapest_unit_male` Identifying the index of cheapest unit in the breeding program on male side
 > **NOTE!** The idea of identifying the cheapest unit in the breeding program is at the end of the program all parameters stay within the base budget and by doing this we are making sure that the cost is not spent in the most 	expensive part but rather to the cheapest unit by identifying the index of cheapest parameter on female and male side	
* `@cost_cheapest_female` Identifying the cost of cheapest unit in the breeding program on female side
* `@cost_cheapest_male` identifying the cost of cheapest unit in the breeding program on male side
* 
**iterations.csv**

The optimization problem's search space varies in size and complexity depending on the number of parameters involved. When dealing with problems that have numerous parameters within a vast search space, it is crucial to carefully select the number of parents and offspring. This selection process ensures that potentially valuable solutions are not lost too early, and that the optimization process does not converge too quickly. To address this, the user can utilize a table to define different parameter settings based on the complexity of the problem. By breaking down the iterations interval, the user can initially generate more parents and offspring, gradually decreasing the numbers as better solutions are found. 

| iteration_start      | iteration_end | n_sel_high_target | n_sel_kernel | n_sel_last2 | n_off_recombination | n_off_mutation | n_off_random | nrep |
|----------------------|---------------|-------------------|--------------|-------------|----------------------|----------------|--------------|------|
| 2                    | 3             | 70                | 30           | 0           | 200                  | 0              | 0            | 1    |
| 4                    | 30            | 30                | 15           | 5           | 170                  | 80             | 0            | 1    |
| 31                   | 50            | 20                | 7            | 3           | 180                  | 90             | 0            | 1    |

These parameter configurations have proven effective in our toy example breeding program.
* `@iteration_start` What parameter settings in which iteration should be used in the interval between iteration_start and iteration_end
* `@iteration_end` What parameter settings in which iteration should be used in the interval between iteration_start and iteration_end
* `@n_sel_high_target` Number of selected parents for the current iteration based on high performance
* `@n_sel_kernel` Number of parents determined using the kernel regression method
* `@n_sel_last2` Number of parents from two cycles before (which should be set to 0 for the second iterations)
* `@n_off_recombination` Number of offspring generated through recombination
* `@n_off_mutation` Number of offspring produced through mutation
* `@nrep` Number of replications, which specifies how many times the same simulation should be performed with different random seeds. As the optimization process progresses in later stages, it is advisable to increase the number of replications. This helps to reduce the variance of the stochastic simulation output, providing more reliable and accurate results.
* `@n_off_random` Number of offspring generated completely from scratch. 

>  **n_off_random** is associated with the `generate_new` function found in the [function folder](https://github.com/AHassanpour88/Evolutionary_Snakemake/edit/main/Functions), where the user has the ability to design a new range for decision variables and to dynamically adjust the range of decision variables based on the ongoing analysis of the optimization process. By examining the results obtained in each iteration, the user can fine-tune and redefine the acceptable ranges for the decision variables. This enables the optimization algorithm to explore different regions of the search space and potentially discover more optimal solutions. Due to the unexplored nature of the search space, it is not advisable to activate this functionality in the early stages.

**Function [folder](https://github.com/AHassanpour88/Evolutionary_Snakemake/tree/main/Functions)**

**cost_function**

The cost_function is an essential function that requires user modification. This function is responsible for generating a cost function which calculates the cost associated with each parameter. In case there is a constraint on the cost scale, the function adjusts the parameters accordingly. It ensures that the parameters remain within the specified cost range. For instance, if the cost exceeds the upper limit, the function automatically scales down all parameters to bring them within the acceptable range. Creating a cost function involves defining a mathematical expression that quantifies the expenses or constraints associated with a particular task or system.  Let's consider a cost function for a breeding program:

```{r}
# Define a cost function for a breeding program
calc_cost <- function(genetic_factors) {
  
  # Genetic factors:
  # genetic_factors[1]: Cost per unit of high-quality seeds
  # genetic_factors[2]: Maintenance cost for optimal growing conditions
  # genetic_factors[3]: Labor cost for genetic analysis
  
  # Calculate the breeding program cost based on the given genetic factors
  breeding_cost <- genetic_factors[1] * 5000 + 
                   genetic_factors[2] * 10000 + 
                   genetic_factors[3] * 3000
  
  return(unname(breeding_cost))
}
# Example usage:
# Suppose the cost per unit of high-quality seeds is 0.50, maintenance cost is 1.00, and labor cost for genetic analysis is 0.30
genetic_factors_example <- c(0.50, 1.00, 0.30)
calc_cost(genetic_factors_example)
```
The breeding program cost is then calculated using a formula that combines these genetic factors. 

In the context of a breeding program, it's common that not all parameters being optimized directly contribute to the cost. Some parameters, even though they play a role in the breeding process, might not incur additional expenses. Therefore, it's important to ensure that only those parameters associated with a cost are considered in the cost function. 

**generate_new**

The generate_new function is responsible for creating new offspring from scratch and requires user modification. This functionality can be activated in the iterations.csv [file](https://github.com/AHassanpour88/Evolutionary_Snakemake/edit/main/config/iterations.csv) by setting the value of 'n_off_random' to the desired number of offspring. Our algorithm incorporates a random search criterion that can be employed to prioritize exploration rather than solely focusing on exploiting the search space. By reducing the influence of other selection criteria and increasing the weight assigned to random search, the algorithm becomes more adept at thoroughly exploring a vast search space. This becomes particularly crucial in later iterations where the likelihood of becoming trapped in local optima is higher. By activateing this criterion, the algorithm gains the ability to uncover unexpected solutions that may have been overlooked during the initialization phase. To enhance flexibility, users have the option to define new boundaries for decision variables based on previous results or extend the search beyond the initial parameter bounds in [sampleScript](https://github.com/AHassanpour88/Evolutionary_Snakemake/edit/main/scripts/sampleScript.r). This enables a more targeted exploration of the solution space and the discovery of new solutions that might have been disregarded in earlier iterations.

**termination_criteria:**
to limit computational time before reaching a maximum number of iterations, we defined dual termination criteria based on the average and standard deviation thresholds for the objective function through a kernel regression applied to the most promising candidate parameters in each iteration. The termination criterion is satisfied when both the mean and standard deviation of the best solutions from the last ten iterations are below predefined threshold values. Determining suitable threshold values for these measurements involves considering the problem-specific characteristics, desired optimization precision, and computational resources available for the optimization process. If the the pipeline does not achieve convergence based on these criteria in evolutionary algorithm, the pipeline will stop when the last number of iterations in the table for the last row in the iteration_end is reached. In general, it is recommended that the decision-makers regularly evaluate and track the progress of the algorithm across different iterations based on **EvoStatus.RData** and adjust the termination criteria during the algorithm's development and implementation. If a user determines that the final termination criterion is insufficient and wishes to continue the process, they should consider the following actions:

1- Remove the final `rplot.pdf` file from the main folder to avoid any conflicting data.

2- Adjust and increase the number of iterations in the [iterations.csv](https://github.com/AHassanpour88/Evolutionary_Snakemake/edit/main/config/iterations.csv).

### Acknowledgment

This research was supported by BASF Belgium Coordination Center CommV.
