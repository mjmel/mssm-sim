# Simulation Code

## Installation

To compile, navigate to the current directory and run `bash compile.sh`.

## Example Usage 

To simulate a population with a single beneficial effect and a single deleterious effect, run `./simulate N NUb NUd NUn sb sd sfs_sample_size num_samples`, where 
* `N` is the size of the population
* `Ub` is the beneficial mutation rate (per-genome, per-generation)
* `Ud` is the deleteriou mutation rate (per-genome, per-generation)
* `Un` is the neutral mutation rate (per-genome, per-generation)
* `sb` is the beneficial fitness effect
* `sd` is the deleterious fitness effect (positive by convention)
* `sfs_sample_size` is the sample size used for measuring the site frequency spectrum
* `num_samples` is the number of roughly independent samples measured

To simulate a population with a beneficial exponential-like DFE and a deleterious exponential-like DFE, run `./simulate_2explike N NUb NUd NUn sigma_b sigma_d beta_b beta_d sfs_sample_size num_samples`. Parameters are the same as above, except
* `sigma_b` is the scale parameter of the beneficial DFE
* `sigma_d` is the scale parameter of the deleterious DFE
* `beta_b` is the steepness parameter of the beneficial DFE
* `beta_d` is the steepness parameter of the deleterious DFE

## Output

Both types of simulation will output the following for each of the `num_samples` samples:
* `Generation`: Generation at which the sample is measured
* `Avg log fitness`: Average log fitness of the population at the sampled generation
* `Pi neutral`: Pairwise heterozygosity of neutral mutations
* `Site frequency spectrum neutral`: Neutral site frequency spectrum, computed using a sample from the population of size `sfs_sample_size`. 
* `Pi nonneutral`: Pairwise heterozygosity of non-neutral mutations
* `Site frequency spectrum neutral`: Non-neutral site frequency spectrum, computed using a sample from the population of size `sfs_sample_size`. 
* `Fixed fitness effects`: Fitness effects of mutations which have fixed within the population since the last sampled generation.
* `Fixed arisal fits`: Relative arisal fitnesses of lineages which have fixed within the population since the last sampled generation.
