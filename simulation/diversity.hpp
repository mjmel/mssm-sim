#ifndef DIVERSITY_HPP
#define DIVERSITY_HPP

#include <valarray>
#include <set>
#include <set>
#include <map>
#include "mutation.hpp"
#include "individual.hpp"

typedef std::valarray<double> FrequencySpectrum;

FrequencySpectrum calculate_frequency_spectrum(Population & diversity_sample, bool get_neutral_sfs){  
	
	int n = diversity_sample.size();
	FrequencySpectrum frequency_spectrum(0.0, n-1);

	MutationList mutations;
	for(auto & individual: diversity_sample){
		mutations.insert(mutations.end(),individual.mutations->begin(), individual.mutations->end());
	}

 	std::set<Mutation> unique_mutations(mutations.begin(), mutations.end());
 	std::multiset<Mutation> all_mutations(mutations.begin(), mutations.end());

	double neutral_value = 1.0;
	for (auto mutation : unique_mutations){
		if ( (get_neutral_sfs && mutation.fitness_effect==neutral_value) || (!get_neutral_sfs && mutation.fitness_effect!=neutral_value) ){
			int i = all_mutations.count(mutation);
			if (i < n)
				frequency_spectrum[i-1]++;
		}
	}

	return frequency_spectrum;
}

double calculate_pi(const FrequencySpectrum & frequency_spectrum){
	double pi(0.0);
	double i = 1.0;
	double n = float(frequency_spectrum.size() + 1);
	for (auto f: frequency_spectrum){
		pi += i * (n - i) * f;
		i++;
	}
	pi /= (n * (n - 1.0) / 2.0);
	return pi;
}

#endif
