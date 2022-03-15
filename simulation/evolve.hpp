#ifndef EVOLVE_HPP
#define EVOLVE_HPP
#include <cmath>

#include "stats.hpp"
#include "individual.hpp"
#include "mutation.hpp"
#include "dfe.hpp"
#include "diversity.hpp"

template<class Type>
void print_vector(const std::vector<Type> & v){
	for (auto & i : v)
		std::cout << i << ' ';
	std::cout << std::endl;
}

class Results {
public:
	int time;
	double avg_log_fitness;
	double pi_neutral;
	std::vector<double> sfs_neutral;
	double pi_nonneutral;
	std::vector<double> sfs_nonneutral;
	std::vector<double> fixed_fitness_effects;
	std::vector<double> fixed_arisal_fits;
   	std::vector<double> relative_fitness;
};

void print_results(Results & results){
    std::cout.precision(10);
    std::cout << "Generation:" << std::endl;
    std::cout << results.time << std::endl;
    std::cout << "Avg log fitness:" << std::endl;
    std::cout<< results.avg_log_fitness << std::endl;
    std::cout << "Pi neutral:" << std::endl;
    std::cout << results.pi_neutral << std::endl;
    std::cout << "Site frequency spectrum neutral:" << std::endl;
    print_vector(results.sfs_neutral); 
    std::cout << "Pi nonneutral:" << std::endl;
    std::cout << results.pi_nonneutral << std::endl;
    std::cout << "Site frequency spectrum nonneutral:" << std:: endl;
    print_vector(results.sfs_nonneutral); 
    std::cout << "Fixed fitness effects:" << std::endl;
    print_vector(results.fixed_fitness_effects);
    std::cout << "Fixed arisal fits:" << std::endl;
    print_vector(results.fixed_arisal_fits);
    // std::cout << "Relative fitnesses:" << std::endl;
    // print_vector(results.relative_fitness);
}

template<class NonneutralDFE, class NeutralDFE>
inline Results evolve(Random & random, LabelGenerator & label_generator, int N, 
						NonneutralDFE & nonneutral_dfe, NeutralDFE & neutral_dfe,
					    int sfs_sample_size, int num_samples){
	Results results;

	results.sfs_neutral.reserve(sfs_sample_size-1);
	results.sfs_nonneutral.reserve(sfs_sample_size-1);
   	results.relative_fitness.reserve(sfs_sample_size);

	Population population;
	Population new_population;
	Population diversity_sample;

	population.reserve(N);
	new_population.reserve(N);
	diversity_sample.reserve(sfs_sample_size);

    double total_fitness = 0.0;
	double pop_mean_fitness = 0.0;
	for (int i = 0; i < N; i++){
	    population.push_back(Individual{total_fitness+=1.0});
	}

	int t(0);
	double log_fitness_baseline = 0.0;
	double log_fitness_baseline_step = 0.0;
	int collected_samples(0);
	bool at_equilibrium = false;
	int output_interval(0);

	FrequencySpectrum freq_spec_tmp_neutral;
	FrequencySpectrum freq_spec_tmp_nonneutral;

	std::cerr << "Collecting samples!" << std::endl;
	while (collected_samples < num_samples)
	{	
		double new_total_fitness = 0.0;
		double temp_mean_fitness = total_fitness/float(N);
		for (int i = 0; i < N; i++){
			new_population.push_back(*std::lower_bound(population.begin(), population.end(), sample_uniform(random)*total_fitness));
			nonneutral_dfe.mutate_individual(random, label_generator, new_population.back(), temp_mean_fitness);
			neutral_dfe.mutate_individual(random, label_generator, new_population.back(), temp_mean_fitness); 
			new_total_fitness += new_population.back().fitness;
			new_population.back().weight = new_total_fitness;
		}

		population.swap(new_population);
		new_population.clear();
		total_fitness = new_total_fitness;
		pop_mean_fitness = total_fitness/float(N); 

		bool found_fixed = true;
		while (found_fixed){
			auto query_for_fixation = population.front().get_earliest_mutation();
			if (query_for_fixation.label != null_label){
				for (auto & individual: population){
					if (individual.get_earliest_mutation().label != query_for_fixation.label){
						found_fixed = false;
						break;
					}
				}	
				if (found_fixed){
					if (!at_equilibrium){    
						output_interval = t;
					}  

					at_equilibrium = true;
					results.fixed_fitness_effects.push_back(std::log(query_for_fixation.fitness_effect));
					results.fixed_arisal_fits.push_back(query_for_fixation.arisal_x);
					for (auto & individual: population){
						individual.remove_earliest_mutation();
					}
				}		
			}
			else 
				found_fixed = false; 
		}

		if (at_equilibrium) {
			if (t % output_interval == 0){
				collected_samples++;
				
				log_fitness_baseline_step = std::log(pop_mean_fitness);
				log_fitness_baseline += log_fitness_baseline_step;
				for (int i = 0; i < N; i++){
					population[i].fitness = population[i].fitness/pop_mean_fitness; 
					population[i].weight = population[i].weight/pop_mean_fitness;	
				}
				total_fitness = total_fitness/pop_mean_fitness;	
			
				results.time = t; 
				results.avg_log_fitness = log_fitness_baseline;
				for (int i = 0; i < sfs_sample_size; i++){ // NOTE: This is a random sample, because individuals were sampled into the population in random order. 
					results.relative_fitness.push_back(std::log(population[i].fitness));
					diversity_sample.push_back(population[i]);
				}
	
				freq_spec_tmp_neutral = calculate_frequency_spectrum(diversity_sample, true);
				freq_spec_tmp_nonneutral = calculate_frequency_spectrum(diversity_sample, false);

				results.sfs_neutral.insert(std::end(results.sfs_neutral), std::begin(freq_spec_tmp_neutral), std::end(freq_spec_tmp_neutral)); 
    			results.pi_neutral = calculate_pi(freq_spec_tmp_neutral);
				results.sfs_nonneutral.insert(std::end(results.sfs_nonneutral), std::begin(freq_spec_tmp_nonneutral), std::end(freq_spec_tmp_nonneutral)); 
    			results.pi_nonneutral = calculate_pi(freq_spec_tmp_nonneutral);

				print_results(results);

				results.sfs_neutral.clear();
				results.sfs_nonneutral.clear();
				results.fixed_fitness_effects.clear();
				results.fixed_arisal_fits.clear();
				results.relative_fitness.clear();
				diversity_sample.clear();

			}
		}	
		t++;   
	}
	std::cerr << "Collected all samples!" << std::endl; 

	return results;
}

#endif
