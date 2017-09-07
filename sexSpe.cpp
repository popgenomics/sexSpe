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
#include <assert.h>

static std::random_device rd;
static std::mt19937 rng(rd());

void initiate_pop(std::vector < std::vector < std::vector < float >>> & pop, const unsigned int & N_diplo, const unsigned int & N_SA, const unsigned int & N_SC);
void initiate_newPop(std::vector < std::vector < std::vector < float >>> & new_pop, const unsigned int & N_diplo);
void evolve_pop(std::vector < std::vector < std::vector < float >>> & pop, const unsigned int & N_diplo, const unsigned int & N_reprod_males, const unsigned int & N_SA, const unsigned int & N_SC, const std::vector <float> & param_fitness,   const float & mutation_rate, const float & recombination_rate, const size_t & nLocus);
void print_pop(std::vector < std::vector < std::vector < float >>> & pop, const unsigned int & N_diplo, const unsigned int & N_SA, const unsigned int & N_SC);
void get_sexes(const std::vector < std::vector < std::vector < float >>> & pop, std::vector <size_t> & sexes, std::vector <size_t> & males, std::vector <size_t> & females);
void get_fitness(const std::vector < std::vector < std::vector <float>>> & pop, const std::vector <size_t> & males, const std::vector <size_t> & females, std::vector <float> & male_fitness, std::vector <float> & female_fitness, const unsigned int & N_diplo, const unsigned int & N_SA, const unsigned int & N_SC, const std::vector <float> & param_fitness);
void get_fathers(std::vector <size_t> & fathers, const std::vector <size_t> & males, const std::vector <float> & male_fitness, const unsigned int & N_diplo, const unsigned int & N_reprod_males);
void get_mothers(std::vector <size_t> & mothers, const std::vector <size_t> & females, const std::vector <float> & female_fitness, const unsigned int & N_diplo);
void sampling(const std::vector <float> & weights, const unsigned int & nSamples, std::vector <size_t> & sample, const bool & replace);
size_t getSampledPosition(const std::vector <float> & urn);

void make_babies(const std::vector <size_t> & fathers, const std::vector <size_t> & mothers, std::vector < std::vector < std::vector < float >>> & pop, std::vector < std::vector < std::vector < float >>> & new_pop, const unsigned int & N_diplo, const float & mutation_rate, const float & recombination_rate, const size_t & nLocus);
void is_recombination(const float & recombination_rate, int & test_recombination, size_t & pos_recombination, const size_t & nLocus);


int main(int argc, char* argv[]){
	const unsigned int N_diplo(std::stoi(argv[1]));
	const unsigned int N_SA(std::stoi(argv[2]));
	const unsigned int N_SC(std::stoi(argv[3]));
	const size_t nGenerations(std::stoi(argv[4]));
	const unsigned int N_reprod_males(std::stoi(argv[5]));

	const size_t nLocus(3 + (N_SA+N_SC)*2); // sex_det / recomb / expression / loc_xp_ntrl / genicV_ntrl / etc ...
	const float mutation_rate( (nLocus-1) * 0.00001); // mutation_rate = nLocus x proba_of_mutation_of_a_locus
	const float recombination_rate( (nLocus-1) * 0.00001); // mutation_rate = nLocus x proba_of_mutation_of_a_locus

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

	for(generation=0; generation<nGenerations; ++generation){
	/*	std::cout << "Generation: " << generation << std::endl;
	std::cout << "Before: " << std::endl;
	print_pop(pop, N_diplo, N_SA, N_SC);*/
		evolve_pop(pop, N_diplo, N_reprod_males, N_SA, N_SC, param_fitness, mutation_rate, recombination_rate, nLocus);
/*	std::cout <<  std::endl << "After: " << std::endl;
	print_pop(pop, N_diplo, N_SA, N_SC);
	std::cout << std::endl;*/
	}	

	return(0);
}

//void sample_with(const std::vector <float> & weights, const unsigned int & nSamples, std::vector <size_t> & sample){
//        /*
//        returns a vector of positions (size_t) corresponding to the sampled positions.
//        */
//        size_t i(0);
//        size_t j(0);
//        float sum(0.0);
//        float random(0.0);
//        float cumul(0.0);
//
//        for(i=0; i<weights.size(); ++i){
//                sum += weights[i];
//        }
//
//        std::uniform_real_distribution<float> distribution_real(0, sum);
//
//        for(i=0; i<nSamples; ++i){
//                random = distribution_real(rng);
//                cumul = 0.0;
//                for(j=0; j<weights.size(); ++j){
//                        cumul += weights[j];
//                        if( random <= cumul){ sample.push_back(j); break;}
//                }
//        }
//}

//void sample_without(const std::vector <float> & weights, const unsigned int & nSamples, std::vector <size_t> & sample){
//        /*
//        fills the vector of positions named 'sample' (size_t) with positions corresponding of the sampled weights.
//        */
//	assert(nSamples <= weights.size());
//	std::vector <float> cumul_weights;
//	std::vector <size_t> positions;
//	size_t i(0);
//	size_t j(0);
//	size_t position_tmp(0);
//	float cumul(0.0);
//	float random_value(0.0);
//
//	for(i=0; i<weights.size(); ++i){
//		cumul += weights[i];
//		cumul_weights.push_back(cumul);
//		positions.push_back(i);
//		++j;
//	}
//
//	i=0;
//	do{
//		--j;
//		random_value = getReal(0.0, cumul_weights[j]);
//		position_tmp = getPosition(cumul_weights, random_value);
//
//		assert(position_tmp<positions.size());
//
//		sample.push_back(positions[position_tmp]);
//
//		clean_vectors(weights, cumul_weights, positions, position_tmp);
//		++i;
//	}while(i<nSamples);
//}


/*size_t getPosition(const std::vector <float> & cumul_weights, const float & random_value){
	size_t pos(0);
	float cumul(0.0);

	size_t i(0);
	
	while(random_value > cumul_weights[pos]){
		++pos;
	}

	return(pos);
}
*/


//float getReal(const float & from, const float & to){
	// returns a real in the [from, to] interval
//       	std::uniform_real_distribution<float> distribution_real(from, to);
//	return( distribution_real(rng) );
//}


/*void clean_vectors(const std::vector <float> & weights, std::vector <float> & cumul_weights, std::vector <size_t> & positions, const size_t & position_tmp){
	size_t i(0);
	float cumul(0.0);
	positions.erase(positions.begin()+position_tmp);
	cumul_weights.clear();

	for(i=0; i<positions.size(); ++i){
		cumul += weights[ positions[i] ];
		cumul_weights.push_back(cumul);
	}	
}
*/

void initiate_pop(std::vector < std::vector < std::vector < float >>> & pop, const unsigned int & N_diplo, const unsigned int & N_SA, const unsigned int & N_SC){
	/*
	locus 1 = sex determination (0/0 = female, 0/1 = male)
	locus 2 = recombination modifier, 0 <= r <= 1, probability of recombination between haplotypes = r1.r2
	locus 3 = modifier of expression, affecting all genes on the haplotype, is sex-dependent. 
	*/

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

			pop[i][j][1] = distribution_modifierRecomb(rng); // locus 2: modifier recomb
			pop[i][j][2] = distribution_modifierExpression(rng); // locus 3: modifier global expression
			pop[i][j][3] = distribution_locExpression(rng); // locus 4: local expression neutral locus
			pop[i][j][4] = distribution_locGenicValue(rng); // locus 5: genic value neutral locus
			
			for( k=0; k<N_SA; ++k){
				pop[i][j][5 + (k*2)] = distribution_locExpression(rng);
				pop[i][j][5 + (k*2 + 1)] = distribution_locGenicValue(rng);
			}
			
			for( k=0; k<N_SC; ++k){
				pop[i][j][5 + N_SA*2 + (k*2)] = distribution_locExpression(rng);
				pop[i][j][5 + N_SA*2 + (k*2 + 1)] = distribution_locGenicValue(rng);
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


void evolve_pop(std::vector < std::vector < std::vector < float >>> & pop, const unsigned int & N_diplo, const unsigned int & N_reprod_males, const unsigned int & N_SA, const unsigned int & N_SC, const std::vector <float> & param_fitness, const float & mutation_rate, const float & recombination_rate, const size_t & nLocus){
	std::vector <size_t> males;
	std::vector <size_t> females;
	std::vector <size_t> fathers;
	std::vector <size_t> mothers;
	std::vector <size_t> sexes; // sexes of individuals within the population
	std::vector <float> male_fitness;
	std::vector <float> female_fitness;
	std::vector < std::vector < std::vector < float >>>  new_pop;

	initiate_newPop(new_pop, N_diplo);

	get_sexes(pop, sexes, males, females);
	assert( males.size() != 0);
	assert( females.size() != 0);

	get_fitness(pop, males, females, male_fitness, female_fitness, N_diplo, N_SA, N_SC, param_fitness);
	
	get_mothers(mothers, females, female_fitness, N_diplo);
	get_fathers(fathers, males, male_fitness, N_diplo, N_reprod_males);
	
	make_babies(fathers, mothers, pop, new_pop, N_diplo, mutation_rate, recombination_rate, nLocus);

/*	std::cout << "New pop: " << new_pop[0][0][2] << std::endl;	
	std::cout << pop[0][0][2] << " -Pop- ";	
	//pop = new_pop;
	std::cout << pop[0][0][2] << std::endl;

	size_t i(0);
	std::cout << "males:\t\t";
	for(i=0; i<males.size(); ++i){
		std::cout << males[i] << " "; 
	}
	std::cout << std::endl;
	std::cout << "male fitness:\t";
	for(i=0; i<males.size(); ++i){
		std::cout << male_fitness[i] << " "; 
	}
	std::cout << std::endl;
	std::cout << "fathers:\t";
	for(i=0; i<N_diplo; ++i){
		std::cout << fathers[i] << " "; 
	}
	std::cout << std::endl;
	std::cout << "females:\t";
	for(i=0; i<females.size(); ++i){
		std::cout << females[i] << " "; 
	}
	std::cout << std::endl;
	std::cout << "female fitness:\t";
	for(i=0; i<females.size(); ++i){
		std::cout << female_fitness[i] << " "; 
	}
	std::cout << std::endl;
	std::cout << "mothers:\t";
	for(i=0; i<N_diplo; ++i){
		std::cout << mothers[i] << " "; 
	}
	std::cout << std::endl << std::endl;
*/	
	males.clear();
	females.clear();
	fathers.clear();
	mothers.clear();
	sexes.clear(); // sexes of individuals within the population
	male_fitness.clear();
	female_fitness.clear();
	new_pop.clear();
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
	/*
	optimal level of expression of a gene in a sex:
		e_m_opt [0]
		e_f_opt [1]
	optimum genic values:
		loci SA:
			g_m_opt_SA [2]
			g_f_opt_SA [3]
		loci SC:
			g_m_opt_SC [4]
			g_f_opt_SC [5]
	sex-specific strength of selection:
		on expression:
			s_m_e [6]
			s_f_e [7]
		on genic value:
			s_m_g [8]
			s_f_g [9]
	*/
	
	size_t i(0);
	size_t j(0);
	const float e_m_opt(param_fitness[0]);
	const float e_f_opt(param_fitness[1]);
	const float g_m_opt_SA(param_fitness[2]);
	const float g_f_opt_SA(param_fitness[3]);
	const float g_m_opt_SC(param_fitness[4]);
	const float g_f_opt_SC(param_fitness[5]);
	const float s_m_e(param_fitness[6]);
	const float s_f_e(param_fitness[7]);
	const float s_m_g(param_fitness[8]);
	const float s_f_g(param_fitness[9]);
	
	float e_1_glob(0.0);
	float e_2_glob(0.0);
	float e_1_sex(0.0);
	float e_2_sex(0.0);
	float e_tot_sex(0.0);
	float g_1(0.0);
	float g_2(0.0);

	float fitness(1.0);
	/* fills the vectors from 'void evolve_pop()':
		male_fitness: vector of size n_males. Entry == male fitness 
		female_fitness: vector of size n_females. Entry == female fitness 
	*/
	
	for(i=0; i<males.size(); ++i){ // loop over males
		fitness = 1.0;
		e_1_glob = pop[males[i]][0][2];
		e_2_glob = pop[males[i]][1][2];
		// loop over loci SA:
		for(j=5; j<5+N_SA*2; j+= 2){
			e_1_sex = pop[males[i]][0][j] * e_1_glob;
			e_2_sex = pop[males[i]][1][j] * e_2_glob;
			e_tot_sex = e_1_sex + e_2_sex;
			g_1 = pop[males[i]][0][j+1];
			g_2 = pop[males[i]][1][j+1];
			fitness *= exp(-1 * s_m_g * e_1_sex * pow(g_1 - g_m_opt_SA, 2));
			fitness *= exp(-1 * s_m_g * e_2_sex * pow(g_2 - g_m_opt_SA, 2));
			fitness *= exp(-1 * s_m_e * pow(e_tot_sex - e_m_opt, 2));
		}
		
		// loop over loci SC:
		for(j=5+N_SA*2; j<5+N_SA*2+N_SC*2; j+= 2){
			e_1_sex = pop[males[i]][0][j] * e_1_glob;
			e_2_sex = pop[males[i]][1][j] * e_2_glob;
			e_tot_sex = e_1_sex + e_2_sex;
			g_1 = pop[males[i]][0][j+1];
			g_2 = pop[males[i]][1][j+1];
			fitness *= exp(-1 * s_m_g * e_1_sex * pow(g_1 - g_m_opt_SC, 2));
			fitness *= exp(-1 * s_m_g * e_2_sex * pow(g_2 - g_m_opt_SC, 2));
			fitness *= exp(-1 * s_m_e * pow(e_tot_sex - e_m_opt, 2));
		}
		male_fitness.push_back(fitness);
	} // end of loop over males

	for(i=0; i<females.size(); ++i){ // loop over females
		fitness = 1.0;
		e_1_glob = pop[females[i]][0][2];
		e_2_glob = pop[females[i]][1][2];
		// loop over loci SA:
		for(j=5; j<5+N_SA*2; j+= 2){
			e_1_sex = pop[females[i]][0][j] * e_1_glob;
			e_2_sex = pop[females[i]][1][j] * e_2_glob;
			e_tot_sex = e_1_sex + e_2_sex;
			g_1 = pop[females[i]][0][j+1];
			g_2 = pop[females[i]][1][j+1];
			fitness *= exp(-1 * s_f_g * e_1_sex * pow(g_1 - g_f_opt_SA, 2));
			fitness *= exp(-1 * s_f_g * e_2_sex * pow(g_2 - g_f_opt_SA, 2));
			fitness *= exp(-1 * s_f_e * pow(e_tot_sex - e_f_opt, 2));
		}
		
		// loop over loci SC:
		for(j=5+N_SA*2; j<5+N_SA*2+N_SC*2; j+= 2){
			e_1_sex = pop[females[i]][0][j] * e_1_glob;
			e_2_sex = pop[females[i]][1][j] * e_2_glob;
			e_tot_sex = e_1_sex + e_2_sex;
			g_1 = pop[females[i]][0][j+1];
			g_2 = pop[females[i]][1][j+1];
			fitness *= exp(-1 * s_f_g * e_1_sex * pow(g_1 - g_f_opt_SC, 2));
			fitness *= exp(-1 * s_f_g * e_2_sex * pow(g_2 - g_f_opt_SC, 2));
			fitness *= exp(-1 * s_f_e * pow(e_tot_sex - e_f_opt, 2));
		}
		female_fitness.push_back(fitness);
	} // end of loop over females
} 


void get_fathers(std::vector <size_t> & fathers, const std::vector <size_t> & males, const std::vector <float> & male_fitness, const unsigned int & N_diplo, const unsigned int & N_reprod_males){
	size_t i(0);
	
	std::vector <float> male_fitness_tmp;
	std::vector <size_t> fathers_tmp;
	std::vector <size_t> fathers_tmp2;
	size_t papa;

	// sampling without replacement of reproductive males.
	sampling(male_fitness, N_reprod_males, fathers_tmp, false); // N_reprod_males are sampled from 'male_fitness' to 'fathers_tmp'

	// production of a vector containing the fitnesses of the sampled reproductive males	
	for(i=0; i<fathers_tmp.size(); ++i){
		male_fitness_tmp.push_back(male_fitness[fathers_tmp[i]]);
	}

	// sampling with replacement of fathers among the reproductive males
	sampling(male_fitness_tmp, N_diplo, fathers_tmp2, true); // N_diplo males are sampled into the vector 'fathers_tmp2'

	for(i=0; i<N_diplo; ++i){
		papa = fathers_tmp[fathers_tmp2[i]];
		fathers.push_back(males[papa]); // returns positions of fathers in the whole males+females population
	}
	
/*	for(i=0; i<fathers_tmp2.size(); ++i){
		fathers_tmp3.push_back(fathers_tmp[fathers_tmp2[i]]);
	}
	
	for(i=0; i<fathers_tmp3.size(); ++i){
		fathers.push_back(males[fathers_tmp3[i]]);
	}*/

}


void get_mothers(std::vector <size_t> & mothers, const std::vector <size_t> & females, const std::vector <float> & female_fitness, const unsigned int & N_diplo){
	size_t i(0);
	
	std::vector <size_t> mothers_tmp;
	
	sampling(female_fitness, N_diplo, mothers_tmp, true);
	
	for(i=0; i<N_diplo; ++i){
		mothers.push_back(females[mothers_tmp[i]]);
	}
}


void sampling(const std::vector <float> & weights, const unsigned int & nSamples, std::vector <size_t> & sample, const bool & replace){
        /*
        fills the vector of positions named 'sample' (size_t) with positions corresponding to the sampled weights.
	arg 1: vector of weights (= the urn)
	arg 2: number of balls to sample from the urn
	arg 3: vector of results to 'return'
	arg 4: boolean. If 'false', then sampling without replacement. If 'true', then sampling with replacement
        */
	if( !replace ){
		assert(nSamples <= weights.size());
	}	
	size_t i(0);
	size_t j(0);
	
	std::vector <float> weights_tmp;
	std::vector <size_t> positions_tmp;
	for(i=0; i<weights.size(); ++i){
		weights_tmp.push_back(weights[i]);
		positions_tmp.push_back(i);
	}

	i=0;	
	while( i<nSamples ){
		j = getSampledPosition(weights_tmp);
		sample.push_back(positions_tmp[j]);
		if( !replace ){
			weights_tmp.erase(weights_tmp.begin()+j);
			positions_tmp.erase(positions_tmp.begin()+j);
		}
		++i;
	}
}

size_t getSampledPosition(const std::vector <float> & urn){
	size_t i(0);
	size_t res(urn.size());

	std::discrete_distribution<> distribution_discrete(urn.begin(), urn.end());

	while( res>= urn.size() ){ // because discrete_distribution is able to sample positions greater than the size of the vector
		res = distribution_discrete(rng);
	}
	
	return(res);
}



void initiate_newPop(std::vector < std::vector < std::vector < float >>> & new_pop, const unsigned int & N_diplo){
	size_t i(0);
	size_t j(0);
	std::vector < std::vector <float>> individual;
	std::vector <float> haplotype;
	
	for(i=0; i<N_diplo; ++i){
		new_pop.push_back(individual);
		for(j=0; j<2; ++j){
			new_pop[i].push_back(haplotype);
		}
	}
}


void make_babies(const std::vector <size_t> & fathers, const std::vector <size_t> & mothers, std::vector < std::vector < std::vector < float >>> & pop, std::vector < std::vector < std::vector < float >>> & new_pop, const unsigned int & N_diplo, const float & mutation_rate, const float & recombination_rate, const size_t & nLocus){
	size_t i(0);
	int gamete_id(0);

	int test_mutation(0);
	int test_recombination(0);
	size_t pos_recombination(0);

	std::vector <float> haplotype;

	std::uniform_int_distribution<> distribution_uniform(0, 1); // chose one of the two possible gametes
	std::binomial_distribution<int> distribution(1, mutation_rate); // test if there is a mutation event along the gamete

	for(i=0; i<N_diplo; ++i){
		// within father
		is_recombination(recombination_rate, test_recombination, pos_recombination, nLocus);
		gamete_id = distribution_uniform(rng);	
		new_pop[i][0] = pop[fathers[i]][gamete_id];
		
		// within mother
		gamete_id = distribution_uniform(rng);
		new_pop[i][1] = pop[mothers[i]][gamete_id];
	}

	for(i=0; i<N_diplo; ++i){
		pop[i] = new_pop[i];
	}
}


void is_recombination(const float & recombination_rate, int & test_recombination, size_t & pos_recombination, const size_t & nLocus){
	std::binomial_distribution<int> distribution(1, recombination_rate); // test if there is a recomb event along the gamete
	test_recombination = distribution(rng);
	
	if( test_recombination == 1 ){
		std::uniform_int_distribution<> distribution_position_locus(0, (nLocus-3)/2 - 1); // chose one of the two possible gametes
		pos_recombination =  distribution_position_locus(rng);
		std::cout << "recombination before locus " << pos_recombination << std::endl;
	}
}

