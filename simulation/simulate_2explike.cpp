/*  Simulates population with specified N, U \rho(s) and records fitness distribution
Code by Matthew Melissa (matthewmelissa@gmail.com) and Ivana Cvijovic (icvijovic@fas.harvard.edu), inspired by Ben Good's code.

SUGGESTED COMPILATION:
bash compile.sh
USAGE:
two-explike: ./simulate_2explike N NUb NUd NUn sigma_b sigma_d beta_b beta_d sfs_sample_size num_samples
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

int simulate_two_explike_DFE(int argc, char* argv[]){
    if (argc != 11)
    {
        std::cout << "USAGE: " << argv[0] << " N NUb NUd NUn sigma_b sigma_d beta_b beta_d sfs_sample_size num_samples" << std::endl;
        return 1;
    }
    else{
        
		std::cout << "Two-explike DFE" << std::endl;
        std::cout << "Param template: " << std::endl << "N NUb NUd NUn sigma_b sigma_d beta_b beta_d sfs_sample_size num_samples" << std::endl;
        
        print_params(argc, argv);
        
        Random random = create_random();
        LabelGenerator label_generator{};
        
        double N = strtod(argv[1], NULL);
        double Ub = strtod(argv[2], NULL)/N;
        double Ud = strtod(argv[3], NULL)/N;
        double Un = strtod(argv[4], NULL)/N;
        double sigma_b = strtod(argv[5], NULL);
        double sigma_d = -strtod(argv[6], NULL);
        double beta_b = strtod(argv[7], NULL);
        double beta_d = strtod(argv[8], NULL);
        int sfs_sample_size = stoi(argv[9], NULL);
        int num_samples = stoi(argv[10], NULL);
        
        ExpLikeDFE mutation_type_b{sigma_b, beta_b, Ub};
        ExpLikeDFE mutation_type_d{sigma_d, beta_d, Ud};
        CompositeDFE<ExpLikeDFE,ExpLikeDFE> non_neutral_dfe{mutation_type_b, mutation_type_d};
        DeltaDFE neutral_dfe{0, Un};
        
        Results results;
        results = evolve(random, label_generator, N, non_neutral_dfe, neutral_dfe, sfs_sample_size, num_samples);
        
        return 0;
    }
}

int main(int argc, char* argv[])
{
        return simulate_two_explike_DFE(argc, argv);
}