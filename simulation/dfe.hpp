#ifndef DFE_HPP
#define DFE_HPP

#include <vector>
#include <cmath>
using namespace std;

#include "stats.hpp"
#include "mutation.hpp"

class ExpLikeDFE{
public:
	double sigma;
	double beta;
	double U;

	ExpLikeDFE(double sigma, double beta, double U):sigma(sigma), beta(beta), U(U){
		fitness_distribution = create_gamma(1/beta, 1.0);
	};

	double get_fitness_effect(Random & random){
		double u_rand = fitness_distribution(random);
		return std::exp( sigma*pow(u_rand, 1/beta) ); 
	}

	Mutation get_mutation(Random & random, LabelGenerator & label_generator){
		return Mutation{label_generator.get_next_label(), this->get_fitness_effect(random)};
	}

	void mutate_individual(Random & random, LabelGenerator & label_generator, Individual & individual, double & pop_mean_fitness){
		int no_mutations = sample_poisson(random, U);
		if (no_mutations > 0){
			for (int i = 0; i < no_mutations; ++i){
				individual.add_mutation(get_mutation(random, label_generator), pop_mean_fitness);
			}
		}
	}

private:
	std::gamma_distribution<> fitness_distribution;
};


class DeltaDFE{
public:
	double w;
	double U;

	DeltaDFE(double s, double U): w{std::exp(s)}, U(U) {};

	double get_fitness_effect(Random & random){ return w; }

	Mutation get_mutation(Random & random, LabelGenerator & label_generator){
		return Mutation{label_generator.get_next_label(), this->get_fitness_effect(random)};
	}

	void mutate_individual(Random & random, LabelGenerator & label_generator, Individual & individual, double & pop_mean_fitness){
		int no_mutations = sample_poisson(random, U);
		if (no_mutations > 0){
			for (int i = 0; i < no_mutations; ++i){
				individual.add_mutation(get_mutation(random, label_generator), pop_mean_fitness);
			}
		}
	}

};

template <class DFE1, class DFE2>
class CompositeDFE{
public:
	DFE1 first_dfe;
	DFE2 second_dfe;

	CompositeDFE(DFE1 first_dfe, DFE2 second_dfe): first_dfe(first_dfe), second_dfe(second_dfe) {};

	void mutate_individual(Random & random, LabelGenerator & label_generator, Individual & individual, double & pop_mean_fitness){
		first_dfe.mutate_individual(random, label_generator, individual, pop_mean_fitness);
		second_dfe.mutate_individual(random, label_generator, individual, pop_mean_fitness);
	}
};

#endif
