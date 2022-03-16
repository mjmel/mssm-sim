/*  Simulates population with specified N, U \rho(s) and records fitness distribution
Code is derived from Ivana Cvijovic (icvijovic@fas.harvard.edu) and Ben Good's code (https://github.com/icvijovic/background-selection).

SUGGESTED COMPILATION:
bash compile.sh
USAGE:
two-effect: ./simulate N NUb NUd NUn sb sd sfs_sample_size num_samples
*/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
using namespace std;

#include "individual.hpp"
#include "mutation.hpp"
#include "dfe.hpp"
#include "evolve.hpp"


void print_params(int argc, char * argv[]){
    std::cout << "Params: " << std::endl;
    for(int i=1; i<argc;i++){
        printf("%s ", argv[i]);
    }
    printf("\n");

}

int simulate_two_effect_DFE(int argc, char* argv[]){
	if (argc != 9)
	{
        std::cout << "USAGE: " << argv[0] << " N NUb NUd NUn sb sd sfs_sample_size num_samples" << std::endl;
        return 1;
    }
    else{
   		std::cout << "Two-effect DFE" << std::endl;
        std::cout << "Param template: " << std::endl << "N NUb NUd NUn sb sd sfs_sample_size num_samples" << std::endl;
    
    	print_params(argc, argv);

		Random random = create_random();
		LabelGenerator label_generator{};

		double N = strtod(argv[1], NULL);
		double Ub = strtod(argv[2], NULL)/N;
        double Ud = strtod(argv[3], NULL)/N;
    	double Un = strtod(argv[4], NULL)/N;
		double sb = strtod(argv[5], NULL);
		double sd = -strtod(argv[6], NULL);
		int sfs_sample_size = stoi(argv[7], NULL);
		int num_samples = stoi(argv[8], NULL);

		DeltaDFE mutation_type_b{sb, Ub};
		DeltaDFE mutation_type_d{sd, Ud};
		CompositeDFE<DeltaDFE, DeltaDFE> non_neutral_dfe{mutation_type_b, mutation_type_d};
		DeltaDFE neutral_dfe{0, Un};

		Results results;
		results = evolve(random, label_generator, N, non_neutral_dfe, neutral_dfe, sfs_sample_size, num_samples);
        
		return 0;
	}
}

int main(int argc, char* argv[])
{
     	return simulate_two_effect_DFE(argc, argv);
}
