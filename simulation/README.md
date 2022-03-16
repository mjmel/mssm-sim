# Simulation Code

## Installation

To compile, navigate to the current directory and run `bash compile.sh`.

## Example Usage 

To simulate a population with a single beneficial effect and a single deleterious effect, run `./simulate N NUb NUd NUn sb sd sfs_sample_size num_samples`, where 
* N is the size of the population
* Ub is the beneficial mutation rate (per-genome, per-generation)
* Ud is the deleteriou mutation rate (per-genome, per-generation)
* Un is the neutral mutation rate (per-genome, per-generation)
* sb is the beneficial fitness effect
* sd>0 is the deleterious fitness effect
* sfs_sample_size is the sample size used for measuring the site frequency spectrum
* num_samples is the number of roughly independent samples measured

To simulate a population with a beneficial exponential-like DFE and a deleterious exponential-like DFE, run `./simulate_2explike N NUb NUd NUn sigma_b sigma_d beta_b beta_d sfs_sample_size num_samples`. Parameters are the same as above, except
* sigma_b is the scale parameter of the beneficial DFE
* sigma_d is the scale parameter of the deleterious DFE
* beta_b is the steepness parameter of the beneficial DFE
* beta_d is the steepness parameter of the deleterious DFE
