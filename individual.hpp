#ifndef INDIVIDUAL_HPP
#define INDIVIDUAL_HPP

#include <memory>
#include "mutation.hpp"

class Individual{
public:
	double weight; 
	double fitness;	
	std::shared_ptr<MutationList> mutations; 

	Individual(double weight): fitness(1.0), weight(weight){mutations = std::make_shared<MutationList>();};
	
	bool operator<(const double p) const { return weight < p; } 

	void add_mutation(Mutation && mutation, double & pop_mean_fitness){ // Note: Mutation && mutation is an rvalue reference
		if (!mutations.unique()){ // Checks if other individual(s) in the population share the same mutations shared_ptr. If so, swap the individual's mutations shared_ptr with a newly constructed shared_ptr.
			std::make_shared<MutationList>(*mutations).swap(mutations);
		}
		fitness *= mutation.fitness_effect;
		mutation.arisal_x = std::log(fitness/pop_mean_fitness);
		mutations->push_back(mutation);
	}

	Mutation & get_earliest_mutation(){
		if(mutations->size() > 0){
			return mutations->front();
		}
		else{
			return null_mutation;
		}
	}

	void remove_earliest_mutation(){
		mutations->erase(mutations->begin());
	}
};

typedef std::vector<Individual> Population;

#endif
