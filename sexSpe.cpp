#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <random>
#include <time.h>


void initiate_pop(std::vector < std::vector < std::vector < float >>> & pop, const unsigned int & N_diplo, const unsigned int & N_SA, const unsigned int & N_SC);
void evolve_pop(std::vector < std::vector < std::vector < float >>> & pop, const unsigned int & N_diplo, const unsigned int & N_SA, const unsigned int & N_SC, const std::vector <float> & param_fitness);
void print_pop(std::vector < std::vector < std::vector < float >>> & pop, const unsigned int & N_diplo, const unsigned int & N_SA, const unsigned int & N_SC);
void get_sexes(const std::vector < std::vector < std::vector < float >>> & pop, std::vector <size_t> & sexes, std::vector <size_t> & males, std::vector <size_t> & females);
void get_fitness(const std::vector < std::vector < std::vector <float>>> & pop, const std::vector <size_t> & males, const std::vector <size_t> & females, std::vector <float> & male_fitness, std::vector <float> & female_fitness, const unsigned int & N_diplo, const unsigned int & N_SA, const unsigned int & N_SC, const std::vector <float> & param_fitness);
void get_parents(std::vector <size_t> & mothers, std::vector <size_t> & females, const std::vector <float> & female_fitness, const unsigned int & N_diplo);
void weighted_random_sampling(const std::vector <float> & weights, const unsigned int & nSamples, std::vector <size_t> & sample);

int main(int argc, char* argv[]){
	const unsigned int N_diplo(std::stoi(argv[1]));
	const unsigned int N_SA(std::stoi(argv[2]));
	const unsigned int N_SC(std::stoi(argv[3]));
	const size_t nGenerations(std::stoi(argv[4]));
	const unsigned int N_reprod_males(std::stoi(argv[5]));

	std::vector <float> param_fitness;
	// optimal level of expression of a gene in a sex:
	const float e_m_opt(1.0); param_fitness.push_back(e_m_opt); // [0]
	const float e_f_opt(1.0); param_fitness.push_back(e_f_opt); // [1]
	// optimum genic values:
	// loci SA:
	const float g_m_opt_SA(1.0); param_fitness.push_back(g_m_opt_SA); // [2]
	const float g_f_opt_SA(0.0); param_fitness.push_back(g_f_opt_SA); // [3]
	// loci SC:
	const float g_m_opt_SC(0.5); param_fitness.push_back(g_m_opt_SC); // [4]
	const float g_f_opt_SC(0.5); param_fitness.push_back(g_f_opt_SC); // [5]
	// sex-specific strength of selection:
	// on expression:
	const float s_m_e(1.0);	param_fitness.push_back(s_m_e); // [6]
	const float s_f_e(1.0);	param_fitness.push_back(s_f_e); // [7]
	// on genic value:
	const float s_m_g(1.0); param_fitness.push_back(s_m_g); // [8]
	const float s_f_g(1.0); param_fitness.push_back(s_f_g); // [9]
	
	size_t generation(0);

	std::vector < std::vector < std::vector < float >>> pop; // [ind][haplotype][alleles at different loci]

	initiate_pop(pop, N_diplo, N_SA, N_SC);
	print_pop(pop, N_diplo, N_SA, N_SC);

	for(generation=0; generation<nGenerations; ++generation){
		evolve_pop(pop, N_diplo, N_SA, N_SC, param_fitness);
	}	
	

	return(0);
}

void initiate_pop(std::vector < std::vector < std::vector < float >>> & pop, const unsigned int & N_diplo, const unsigned int & N_SA, const unsigned int & N_SC){
	/*
	locus 1 = sex determination (0/0 = female, 0/1 = male)
	locus 2 = recombination modifier, 0 <= r <= 1, probability of recombination between haplotypes = r1.r2
	locus 3 = modifier of expression, affecting all genes on the haplotype, is sex-dependent. 
	*/

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<float> distribution_modifierRecomb(0, 1);
	std::uniform_real_distribution<float> distribution_modifierExpression(0, 1);
	std::uniform_real_distribution<float> distribution_locExpression(0, 1);
	std::uniform_real_distribution<float> distribution_locGenicValue(0, 1);
	
	size_t i(0);
	size_t j(0);
	size_t k(0);
	
	std::vector <float> haplotype;
	std::vector < std::vector <float>> individual;

	for( i=0; i<(3 + (1 + N_SA + N_SC)*2); ++i ){
		haplotype.push_back(0);
	}

	for( i=0; i<N_diplo; ++i ){
		pop.push_back(individual);
		for( j=0; j<2; ++j ){
			pop[i].push_back(haplotype);
		}
		for( j=0; j<2; ++j ){
			if( i%2 == 0 ){
				pop[i][0][0] = 0;
				pop[i][1][0] = 0;
				
			}else{
				pop[i][0][0] = 0;
				pop[i][1][0] = 1;
				
			}

			pop[i][j][1] = distribution_modifierRecomb(rd); // locus 2: modifier recomb
			pop[i][j][2] = distribution_modifierExpression(rd); // locus 3: modifier global expression
			pop[i][j][3] = distribution_locExpression(rd); // locus 4: local expression neutral locus
			pop[i][j][4] = distribution_locGenicValue(rd); // locus 5: genic value neutral locus
			
			for( k=0; k<N_SA; ++k){
				pop[i][j][5 + (k*2)] = distribution_locExpression(rd);
				pop[i][j][5 + (k*2 + 1)] = distribution_locGenicValue(rd);
			}
			
			for( k=0; k<N_SC; ++k){
				pop[i][j][5 + N_SA*2 + (k*2)] = distribution_locExpression(rd);
				pop[i][j][5 + N_SA*2 + (k*2 + 1)] = distribution_locGenicValue(rd);
			}
		}
	}
}


void print_pop(std::vector < std::vector < std::vector < float >>> & pop, const unsigned int & N_diplo, const unsigned int & N_SA, const unsigned int & N_SC){
	size_t i(0);
	size_t j(0);
	size_t k(0);
	
	std::cout << "N individuals in pop: " << pop.size() << std::endl;
	for(i=0; i<N_diplo; ++i){
		std::cout << "N alleles: " << pop[i].size() << std::endl;
	
		for(j=0; j<2; ++j){
			std::cout << "loci allele " << j << ": ";
			for(k=0; k<pop[i][j].size(); ++k){
				std::cout << pop[i][j][k] << " ";
			}
			std::cout << std::endl;
		}
	}
}


void evolve_pop(std::vector < std::vector < std::vector < float >>> & pop, const unsigned int & N_diplo, const unsigned int & N_SA, const unsigned int & N_SC, const std::vector <float> & param_fitness){
	std::vector <size_t> males;
	std::vector <size_t> females;
	std::vector <size_t> fathers;
	std::vector <size_t> mothers;
	std::vector <size_t> sexes; // sexes of individuals within the population
	std::vector <float> male_fitness;
	std::vector <float> female_fitness;

	get_sexes(pop, sexes, males, females);
	get_fitness(pop, males, females, male_fitness, female_fitness, N_diplo, N_SA, N_SC, param_fitness);
	get_parents(fathers, males, male_fitness, N_diplo);
	get_parents(mothers, females, female_fitness, N_diplo);
}


void get_sexes(const std::vector < std::vector < std::vector < float >>> & pop, std::vector <size_t> & sexes, std::vector <size_t> & males, std::vector <size_t> & females){
	/* fills the vectors from 'void evolve_pop()':
		sexes: vector of size N_pop. Entry == 0 for female. Entry == 1 for male.
		males: vector of size n_males. Entry == position 0-based in vector 'pop' of male.
		females: vector of size n_females. Entry == position 0-based in vector 'pop' of female.
	*/
	
	size_t i(0);
	
	for(i=0; i<pop.size(); ++i){
		if( pop[i][0][0] + pop[i][1][0] == 0){
			females.push_back(i);
			sexes.push_back(0);
		}else{
			males.push_back(i);
			sexes.push_back(1);
		}
	}
}


void get_fitness(const std::vector < std::vector < std::vector <float>>> & pop, const std::vector <size_t> & males, const std::vector <size_t> & females, std::vector <float> & male_fitness, std::vector <float> & female_fitness, const unsigned int & N_diplo, const unsigned int & N_SA, const unsigned int & N_SC, const std::vector <float> & param_fitness){
	size_t i(0);
	float fitness(0.0);
	/* fills the vectors from 'void evolve_pop()':
		male_fitness: vector of size n_males. Entry == male fitness 
		female_fitness: vector of size n_females. Entry == female fitness 
	*/
	
	for(i=0; i<males.size(); ++i){
		fitness = 0.0;
		
	}

} 


void get_parents(std::vector <size_t> & mothers, std::vector <size_t> & females, const std::vector <float> & female_fitness, const unsigned int & N_diplo){
	
}


void weighted_random_sampling(const std::vector <float> & weights, const unsigned int & nSamples, std::vector <size_t> & sample){
        /*
        returns a vector of positions (size_t) corresponding to the sampled positions.
        */
        size_t i(0);
        size_t j(0);
        float sum(0.0);
        float random(0.0);
        float cumul(0.0);

        for(i=0; i<weights.size(); ++i){
                sum += weights[i];
        }

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<float> distribution_real(0, sum);

        for(i=0; i<nSamples; ++i){
                random = distribution_real(rd);
                cumul = 0.0;
                for(j=0; j<weights.size(); ++j){
                        cumul += weights[j];
                        if( random <= cumul){ sample.push_back(j); break;}
                }
        }
}


