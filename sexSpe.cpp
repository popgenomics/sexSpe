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
void evolve_pop(std::vector < std::vector < std::vector < float >>> & pop, const unsigned int & N_diplo, const unsigned int & N_reprod_males, const unsigned int & N_SA, const unsigned int & N_SC, const std::vector <float> & param_fitness, const float & mutation_rate, const float & inversion_rate, const size_t & nLocus, const std::string & name, const size_t & generation, const size_t & nGenerations);
void write_pop(std::vector < std::vector < std::vector < float >>> & pop, const unsigned int & N_diplo, const unsigned int & N_SA, const unsigned int & N_SC, const std::string & name, const size_t & generation);
void write_pop_sum(std::vector < std::vector < std::vector < float >>> & pop, const unsigned int & N_diplo, const unsigned int & N_SA, const unsigned int & N_SC, const std::string & name, const size_t & generation);
void get_sexes(const std::vector < std::vector < std::vector < float >>> & pop, std::vector <size_t> & sexes, std::vector <size_t> & males, std::vector <size_t> & females);
void get_fitness(const std::vector < std::vector < std::vector <float>>> & pop, const std::vector <size_t> & males, const std::vector <size_t> & females, std::vector <float> & male_fitness, std::vector <float> & female_fitness, const unsigned int & N_diplo, const unsigned int & N_SA, const unsigned int & N_SC, const std::vector <float> & param_fitness);
void get_fathers(std::vector <size_t> & fathers, const std::vector <size_t> & males, const std::vector <float> & male_fitness, const unsigned int & N_diplo, const unsigned int & N_reprod_males);
void get_mothers(std::vector <size_t> & mothers, const std::vector <size_t> & females, const std::vector <float> & female_fitness, const unsigned int & N_diplo);
void sampling(const std::vector <float> & weights, const unsigned int & nSamples, std::vector <size_t> & sample, const bool & replace);
size_t getSampledPosition(const std::vector <float> & urn);

void make_babies(const std::vector <size_t> & fathers, const std::vector <size_t> & mothers, std::vector < std::vector < std::vector < float >>> & pop, std::vector < std::vector < std::vector < float >>> & new_pop, const unsigned int & N_diplo, const float & mutation_rate, const float & inversion_rate, const size_t & nLocus, const unsigned int & N_SA, const unsigned int & N_SC);
void is_recombination(const float & recombination_rate, int & test_recombination, size_t & pos_recombination, const unsigned int & N_SA, const unsigned int & N_SC);
void recombination(const std::vector < std::vector <float >> & parent, const size_t & pos_recombination, const int & gamete_id, std::vector <float> & gamete);
void mutation(std::vector <float> & gamete, const size_t & nLocus);
void write_results(const std::vector < std::vector < std::vector <float>>> & pop, const unsigned & N_diplo, const size_t & generation, const std::string & name, const unsigned int & N_SA, const unsigned int & N_SC, const bool & header);

void write_summary(const std::vector < std::vector < std::vector <float>>> & pop, const std::vector <float> & male_fitness, const std::vector <float> & female_fitness, const std::vector <size_t> & sexes, const unsigned int & N_diplo, const unsigned int & N_SA, const unsigned int & N_SC, const size_t & generation);
float mean(const std::vector <float> & x);

int main(int argc, char* argv[]){
	const unsigned int N_diplo(std::stoi(argv[1]));
	const unsigned int N_SA(std::stoi(argv[2]));
	const unsigned int N_SC(std::stoi(argv[3]));
	const size_t nGenerations(std::stoi(argv[4]));
	const unsigned int N_reprod_males(std::stoi(argv[5]));
	const std::string name(argv[6]);

	const size_t nLocus(4 + (1+N_SA+N_SC)*3); // sex_det / recomb / expression (male and female) / loc_xp_ntrl / genicV_ntrl / etc ...
	const float mutation_rate( (nLocus-1) * 0.001); // mutation_rate = nLocus x proba_of_mutation_of_a_locus
	const float inversion_rate(0.0001);

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

	for(generation=0; generation<=nGenerations; ++generation){
		evolve_pop(pop, N_diplo, N_reprod_males, N_SA, N_SC, param_fitness, mutation_rate, inversion_rate, nLocus, name, generation, nGenerations);
	}	

	return(0);
}


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

	for( i=0; i<(4 + (1 + N_SA + N_SC)*3); ++i ){
		haplotype.push_back(0.0);
	}

	for( i=0; i<N_diplo; ++i ){
		pop.push_back(individual);
		for( j=0; j<2; ++j ){
			pop[i].push_back(haplotype);
		}
		for( j=0; j<2; ++j ){
			if( i%2 == 0 ){
				pop[i][0][0] = 0; // female
				pop[i][1][0] = 0;
				
			}else{
				pop[i][0][0] = 1; // male
				pop[i][1][0] = 0;
				
			}

			pop[i][j][1] = distribution_modifierRecomb(rng); // locus 2: modifier recomb
			pop[i][j][2] = distribution_modifierExpression(rng); // locus 3: modifier global expression (g_e_m)
			pop[i][j][3] = distribution_modifierExpression(rng); // locus 3: modifier global expression (g_e_f)
			pop[i][j][4] = distribution_locExpression(rng); // locus 4: local expression neutral locus (ntrl_loc_e_m)
			pop[i][j][5] = distribution_locExpression(rng); // locus 4: local expression neutral locus (ntrl_loc_e_f)
			pop[i][j][6] = distribution_locGenicValue(rng); // locus 5: genic value neutral locus
			
			for( k=0; k<N_SA; ++k){
				pop[i][j][7 + (k*3)] = distribution_locExpression(rng);
				pop[i][j][7 + (k*3 + 1)] = distribution_locExpression(rng);
				pop[i][j][7 + (k*3 + 2)] = distribution_locGenicValue(rng);
			}
			
			for( k=0; k<N_SC; ++k){
				pop[i][j][7 + N_SA*3 + (k*3)] = distribution_locExpression(rng);
				pop[i][j][7 + N_SA*3 + (k*3 + 1)] = distribution_locExpression(rng);
				pop[i][j][7 + N_SA*3 + (k*3 + 2)] = distribution_locGenicValue(rng);
			}
		}
	}
}


void write_pop(std::vector < std::vector < std::vector < float >>> & pop, const unsigned int & N_diplo, const unsigned int & N_SA, const unsigned int & N_SC, const std::string & name, const size_t & generation){
	std::string outfileName = std::string("population_") + name + "_generation_" + std::to_string(generation) + ".txt";
	std::ofstream outputFlux(outfileName.c_str(), std::ios::out);

	if( outputFlux ){
		size_t i(0);
		size_t j(0);
		size_t k(0);

		// header
		for(i=0; i<N_diplo; ++i){
			for(j=0; j<2; ++j){
				outputFlux << "ind_" << i << "_" << j << "\t";
			}
		}
		outputFlux << std::endl;
		
		// haplotypes 
		for(i=0; i<(4 + (1+N_SA+N_SC)*3); ++i){ // loop over loci
			for(j=0; j<N_diplo; ++j){ // loop over individuals
				for(k=0; k<2; ++k){ // loop over alleles
					outputFlux << pop[j][k][i] << "\t";
				}
			}
			outputFlux << std::endl;
		}
	outputFlux.close();
	}else{
		std::cerr <<  "ERROR: cannot open the file " << outfileName << std::endl;
		exit(EXIT_FAILURE);
	}
}

void write_pop_sum(std::vector < std::vector < std::vector < float >>> & pop, const unsigned int & N_diplo, const unsigned int & N_SA, const unsigned int & N_SC, const std::string & name, const size_t & generation){
	std::string outfileName = std::string("population_") + name + ".txt";
	if(generation == 0){
		std::ofstream outputFlux(outfileName.c_str(), std::ios::out);
		if( outputFlux ){
			outputFlux << "generation\t";
			outputFlux << "recomb_X\trecomb_Y\t";
			outputFlux << "globXP_m_X\tglobXP_m_Y\tglobXP_f_X\tglobXP_f_Y\t";
			outputFlux << "ntrl_XP_m_X\tntrl_XP_m_Y\tntrl_XP_f_X\tntrl_XP_f_Y\tntrl_val_X\tntrl_val_Y\t";
			outputFlux << "SA_XP_m_X\tSA_XP_m_Y\tSA_XP_f_X\tSA_XP_f_Y\tSA_val_X\tSA_val_Y\t";
			outputFlux << "SC_XP_m_X\tSC_XP_m_Y\tSC_XP_f_X\tSC_XP_f_Y\tSC_val_X\tSC_val_Y" << std::endl;
			outputFlux.close();
		}else{
			std::cerr <<  "ERROR: cannot open the file " << outfileName << std::endl;
			exit(EXIT_FAILURE);
		}
	}

	std::ofstream outputFlux(outfileName.c_str(), std::ios::app);

	if( outputFlux ){
		size_t i(0);
		size_t j(0);
		size_t k(0);

		std::vector <float> recomb_X; std::vector <float> recomb_Y;
		std::vector <float> globXP_m_X; std::vector <float> globXP_m_Y; std::vector <float> globXP_f_X; std::vector <float> globXP_f_Y;
		std::vector <float> ntrl_XP_m_X; std::vector <float> ntrl_XP_m_Y; std::vector <float> ntrl_XP_f_X; std::vector <float> ntrl_XP_f_Y; std::vector <float> ntrl_val_X; std::vector <float> ntrl_val_Y;
		std::vector <float> SC_XP_m_X; std::vector <float> SC_XP_m_Y; std::vector <float> SC_XP_f_X; std::vector <float> SC_XP_f_Y; std::vector <float> SC_val_X; std::vector <float> SC_val_Y;
		std::vector <float> SA_XP_m_X; std::vector <float> SA_XP_m_Y; std::vector <float> SA_XP_f_X; std::vector <float> SA_XP_f_Y; std::vector <float> SA_val_X; std::vector <float> SA_val_Y;
		
		for(i=0; i<N_diplo; ++i){
			for(j=0; j<2; ++j){
				if( pop[i][j][0]==1 ){ // if Y chromosome
					recomb_Y.push_back(pop[i][j][1]);
					globXP_m_Y.push_back(pop[i][j][2]);
					globXP_f_Y.push_back(pop[i][j][3]);
					ntrl_XP_m_Y.push_back(pop[i][j][4]);
					ntrl_XP_f_Y.push_back(pop[i][j][5]);
					ntrl_val_Y.push_back(pop[i][j][6]);
					
					for(k=0; k<N_SA; k++){
						SA_XP_m_Y.push_back(pop[i][j][7 + 3*k]);
						SA_XP_f_Y.push_back(pop[i][j][7 + 3*k + 1]);
						SA_val_Y.push_back(pop[i][j][7  + 3*k + 2]);
					}
					
					for(k=0; k<N_SC; k++){
						SC_XP_m_Y.push_back(pop[i][j][7 + 3*(N_SA + k)]);
						SC_XP_f_Y.push_back(pop[i][j][7 + 3*(N_SA + k) + 1]);
						SC_val_Y.push_back(pop[i][j][7  + 3*(N_SA + k) + 2]);
					}
				}else{ // if X chromosome
					recomb_X.push_back(pop[i][j][1]);
					globXP_m_X.push_back(pop[i][j][2]);
					globXP_f_X.push_back(pop[i][j][3]);
					ntrl_XP_m_X.push_back(pop[i][j][4]);
					ntrl_XP_f_X.push_back(pop[i][j][5]);
					ntrl_val_X.push_back(pop[i][j][6]);
					
					for(k=0; k<N_SA; k++){
						SA_XP_m_X.push_back(pop[i][j][7 + 3*k]);
						SA_XP_f_X.push_back(pop[i][j][7 + 3*k + 1]);
						SA_val_X.push_back(pop[i][j][7  + 3*k + 2]);
					}
					
					for(k=0; k<N_SC; k++){
						SC_XP_m_X.push_back(pop[i][j][7 + 3*(N_SA + k)]);
						SC_XP_f_X.push_back(pop[i][j][7 + 3*(N_SA + k) + 1]);
						SC_val_X.push_back(pop[i][j][7  + 3*(N_SA + k) + 2]);
					}
				}
			} // end of loop ove haplotypes
		} // end of loop over individuals

	
		outputFlux << generation << "\t";
		outputFlux << mean(recomb_X) << "\t" << mean(recomb_Y) << "\t";
		outputFlux << mean(globXP_m_X) << "\t" << mean(globXP_m_Y) << "\t" << mean(globXP_f_X) << "\t" << mean(globXP_f_Y) << "\t";
		outputFlux << mean(ntrl_XP_m_X) << "\t" << mean(ntrl_XP_m_Y) << "\t" << mean(ntrl_XP_f_X) << "\t" << mean(ntrl_XP_f_Y) << "\t" << mean(ntrl_val_X) << "\t" << mean(ntrl_val_Y) << "\t";
		outputFlux << mean(SA_XP_m_X) << "\t" << mean(SA_XP_m_Y) << "\t" << mean(SA_XP_f_X) << "\t" << mean(SA_XP_f_Y) << "\t" << mean(SA_val_X) << "\t" << mean(SA_val_Y) << "\t";
		outputFlux << mean(SC_XP_m_X) << "\t" << mean(SC_XP_m_Y) << "\t" << mean(SC_XP_f_X) << "\t" << mean(SC_XP_f_Y) << "\t" << mean(SC_val_X) << "\t" << mean(SC_val_Y) << std::endl;

		outputFlux.close();

	}else{
		std::cerr <<  "ERROR: cannot open the file " << outfileName << std::endl;
		exit(EXIT_FAILURE);
	}
}


void evolve_pop(std::vector < std::vector < std::vector < float >>> & pop, const unsigned int & N_diplo, const unsigned int & N_reprod_males, const unsigned int & N_SA, const unsigned int & N_SC, const std::vector <float> & param_fitness, const float & mutation_rate, const float & inversion_rate, const size_t & nLocus, const std::string & name, const size_t & generation, const size_t & nGenerations){
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
	assert( males.size() != 0 ); // There is no male in the population"
	assert( females.size() != 0 ); // There is no female in the population"

	get_fitness(pop, males, females, male_fitness, female_fitness, N_diplo, N_SA, N_SC, param_fitness);
	
	get_mothers(mothers, females, female_fitness, N_diplo);
	get_fathers(fathers, males, male_fitness, N_diplo, N_reprod_males);
	
	make_babies(fathers, mothers, pop, new_pop, N_diplo, mutation_rate, inversion_rate, nLocus, N_SA, N_SC);
	
	if( generation%100 == 0 || generation==nGenerations ){
		write_pop_sum(pop, N_diplo, N_SA, N_SC, name, generation);
//		write_summary(pop, male_fitness, female_fitness, sexes, N_diplo, N_SA, N_SC, generation);
//		write_results(pop, N_diplo, generation, name, N_SA, N_SC, false);
	}


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
	// 0: sex det
	// 1: global_recomb
	// 2: global_expression in male
	// 3: global_expression in female
	// 4: local_expression in male
	// 5: local_expression in female
	// 6: genic value
	// ...
	
	for(i=0; i<males.size(); ++i){ // loop over males
		fitness = 1.0;
		e_1_glob = pop[males[i]][0][2];
		e_2_glob = pop[males[i]][1][2];
		// loop over loci SA:
		for(j=7; j<7+N_SA*3; j+= 3){ // j=7; j=10 (for N_SA=2)
			e_1_sex = pop[males[i]][0][j] * e_1_glob;
			e_2_sex = pop[males[i]][1][j] * e_2_glob;
			e_tot_sex = e_1_sex + e_2_sex;
			g_1 = pop[males[i]][0][j+2];
			g_2 = pop[males[i]][1][j+2];
			fitness *= exp(-1 * s_m_g * e_1_sex * pow(g_1 - g_m_opt_SA, 2));
			fitness *= exp(-1 * s_m_g * e_2_sex * pow(g_2 - g_m_opt_SA, 2));
			fitness *= exp(-1 * s_m_e * pow(e_tot_sex - e_m_opt, 2));
		}
		
		// loop over loci SC:
		for(j=7+N_SA*3; j<7+(N_SA+N_SC)*3; j+= 3){ //j=13; j=16; j=19 (for N_SA=2 and N_SC=3)
			e_1_sex = pop[males[i]][0][j] * e_1_glob;
			e_2_sex = pop[males[i]][1][j] * e_2_glob;
			e_tot_sex = e_1_sex + e_2_sex;
			g_1 = pop[males[i]][0][j+2];
			g_2 = pop[males[i]][1][j+2];
			fitness *= exp(-1 * s_m_g * e_1_sex * pow(g_1 - g_m_opt_SC, 2));
			fitness *= exp(-1 * s_m_g * e_2_sex * pow(g_2 - g_m_opt_SC, 2));
			fitness *= exp(-1 * s_m_e * pow(e_tot_sex - e_m_opt, 2));
		}
		male_fitness.push_back(fitness);
	} // end of loop over males

	for(i=0; i<females.size(); ++i){ // loop over females
		fitness = 1.0;
		e_1_glob = pop[females[i]][0][3];
		e_2_glob = pop[females[i]][1][3];
		// loop over loci SA:
		for(j=5; j<7+N_SA*3; j+= 3){
			e_1_sex = pop[females[i]][0][j+1] * e_1_glob;
			e_2_sex = pop[females[i]][1][j+1] * e_2_glob;
			e_tot_sex = e_1_sex + e_2_sex;
			g_1 = pop[females[i]][0][j+2];
			g_2 = pop[females[i]][1][j+2];
			fitness *= exp(-1 * s_f_g * e_1_sex * pow(g_1 - g_f_opt_SA, 2));
			fitness *= exp(-1 * s_f_g * e_2_sex * pow(g_2 - g_f_opt_SA, 2));
			fitness *= exp(-1 * s_f_e * pow(e_tot_sex - e_f_opt, 2));
		}
		
		// loop over loci SC:
		for(j=7+N_SA*3; j<7+(N_SA+N_SC)*3; j+= 3){
			e_1_sex = pop[females[i]][0][j+1] * e_1_glob;
			e_2_sex = pop[females[i]][1][j+1] * e_2_glob;
			e_tot_sex = e_1_sex + e_2_sex;
			g_1 = pop[females[i]][0][j+2];
			g_2 = pop[females[i]][1][j+2];
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
	if( !replace ){ // without replacement
		assert(nSamples <= weights.size());// The number of elements to sample is bigger than the urn;
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


void make_babies(const std::vector <size_t> & fathers, const std::vector <size_t> & mothers, std::vector < std::vector < std::vector < float >>> & pop, std::vector < std::vector < std::vector < float >>> & new_pop, const unsigned int & N_diplo, const float & mutation_rate, const float & inversion_rate, const size_t & nLocus, const unsigned int & N_SA, const unsigned int & N_SC){
	size_t i(0);
	int gamete_id(0);

	int test_mutation(0);
	int test_recombination(0);
	int test_inversion(0);

	size_t pos_recombination(0);

	std::vector <float> gamete;

	std::uniform_int_distribution<> distribution_uniform(0, 1); // chose one of the two possible gametes
	std::binomial_distribution<int> distribution(1, mutation_rate); // test if there is a mutation event along the gamete
	std::binomial_distribution<int> distribution_inversion(1, inversion_rate);

	for(i=0; i<N_diplo; ++i){ // loop over parents: gametogenese (recombination + mutation) and then fertilization
		////////////
		// father //
		////////////
		// sampling among 2 haplotypes for a diploid
		gamete_id = distribution_uniform(rng);

		// recombination
		is_recombination(pop[fathers[i]][0][1]*pop[fathers[i]][1][1], test_recombination, pos_recombination, N_SA, N_SC); // first arg = r1.r2 (product of recombination rates coded by the 2 haplotypes)

		if( test_recombination ){ // if recombination occurs
			recombination(pop[fathers[i]], pos_recombination, gamete_id, gamete); // recombination
		}else{
			gamete = pop[fathers[i]][gamete_id]; // no recombination
		}

		// inversion of the Y chromosome
		test_inversion = distribution_inversion(rng);
		if( test_inversion ){
			if(gamete[0]==1){ gamete[1]=0; }
		}

		// mutation
		test_mutation = distribution(rng);
		if( test_mutation ){ // if mutation occurs
			mutation(gamete, nLocus);
		}
		new_pop[i][0] = gamete;
		
		gamete.clear();


		////////////
		// mother //
		////////////
		// sampling among 2 haplotypes for a diploid
		gamete_id = distribution_uniform(rng);

		// recombination
		is_recombination(pop[mothers[i]][0][1]*pop[mothers[i]][1][1], test_recombination, pos_recombination, N_SA, N_SC); // first arg = r1.r2 (product of recombination rates coded by the 2 haplotypes)

		if( test_recombination ){
			recombination(pop[mothers[i]], pos_recombination, gamete_id, gamete); // recombination
		}else{
			gamete = pop[mothers[i]][gamete_id];
		}

		// mutation
		test_mutation = distribution(rng);
		if( test_mutation ){ // if mutation occurs
			mutation(gamete, nLocus);
		}
		new_pop[i][1] = gamete;

		gamete.clear();


	}

	for(i=0; i<N_diplo; ++i){
		pop[i] = new_pop[i];
	}
	new_pop.clear();
}


void is_recombination(const float & recombination_rate, int & test_recombination, size_t & pos_recombination, const unsigned int & N_SA, const unsigned int & N_SC){
	/* returns:
		1) test_recombination (0 or 1): 0 if no recombination; 1 if recombination occured
		2) pos_recombination (an even integer): 0 if it occured before the neutral "locExp - GenicValue", 2 before the first locus, 4 ...
							has to be converted for taking into account the i) sex det, ii) recombination and
							iii) global expression loci. In the current version, such conversion only requires to 
							add +3 to the returned pos_recombination value.
	*/	
	std::binomial_distribution<int> distribution(1, recombination_rate); // test if there is a recomb event along the gamete
	test_recombination = distribution(rng);
	
	if( test_recombination == 1 ){
		std::uniform_int_distribution<> distribution_position_locus(0, N_SA+N_SC);
		pos_recombination =  4 + 3*distribution_position_locus(rng);
	}
}


void recombination(const std::vector < std::vector <float >> & parent, const size_t & pos_recombination, const int & gamete_id, std::vector <float> & gamete){
	assert( parent.size() == 2 ); // The parent is not diploid (?!?)
	
	size_t i(0);

	for( i=0; i<parent[0].size(); ++i ){
		if( i<pos_recombination ){
			gamete.push_back( parent[gamete_id][i] );
		}else{
			gamete.push_back( parent[std::abs(1-gamete_id)][i] );
		}
	}
}


void mutation(std::vector <float> & gamete, const size_t & nLocus){
	std::uniform_int_distribution<> distribution_position_locus(1, nLocus-1);
	const int mutated_locus(distribution_position_locus(rng));	

	// version 1: mutation effect is randomly sampled in [0-1]	
//	std::uniform_real_distribution<float> distribution_mutationEffect(0, 1);
//	if(gamete[mutated_locus] != 0 ){	
//		gamete[mutated_locus] = distribution_mutationEffect(rng) ; // no restauration of a lost of function or inversion
//	}

	// version 2: mutation effect is conditionated by the previous value. new_value = old_value x [0.9 - 1.1]
	std::uniform_real_distribution<float> distribution_mutationEffect(0.9, 1.1);
	gamete[mutated_locus] = gamete[mutated_locus]*distribution_mutationEffect(rng);
	if( gamete[mutated_locus]>1 ){
		gamete[mutated_locus] = 1;
	}

}


void write_results(const std::vector < std::vector < std::vector <float>>> & pop, const unsigned & N_diplo, const size_t & generation, const std::string & name, const unsigned int & N_SA, const unsigned int & N_SC, const bool & header){
	std::string outfileName = std::string("output_") + name + ".txt";

	if( header ){
		std::ofstream outputFlux(outfileName.c_str(), std::ios::out);
		if(outputFlux){
			size_t i(0);
			size_t j(0);
			outputFlux << "Generation";
			for(i=0; i<N_diplo; ++i){
				outputFlux << "\tind_" << i << "\tsex_det_ind" << i;
				outputFlux << "\trecomb_ind" << i;
				outputFlux << "\tglob_exp_ind" << i;
				outputFlux << "\tntrl_exp_1/ntrl_exp_2_ind" << i;
				outputFlux << "\tntrl_gen_1/ntrl_gen_2_ind" << i;
				
				for(j=0; j<(N_SA); ++j){
					outputFlux << "\tSA_" << j << "_exp_1/SA_" << j << "_exp_2_ind" << i;
					outputFlux << "\tSA_" << j << "_gen_1/SA_" << j << "_gen_2_ind" << i;
				}
			
				for(j=0; j<(N_SC); ++j){
					outputFlux << "\tSC_" << j << "_exp_1/SC_" << j << "_exp_2_ind" << i;
					outputFlux << "\tSC_" << j << "_gen_1/SC_" << j << "_gen_2_ind" << i;
				}
			}
			outputFlux << std::endl;
		}else{
			std::cerr <<  "ERROR: cannot open the file " << outfileName << std::endl;
			exit(EXIT_FAILURE);
		}
		
	}else{
		std::ofstream outputFlux(outfileName.c_str(), std::ios::app);
		if(outputFlux){
			size_t i(0);
			size_t j(0);
			outputFlux << "generation_" << generation;
			for(i=0; i<N_diplo; ++i){
				outputFlux << "\tind_" << i << ":";
				for(j=0; j<pop[i][0].size(); ++j){
					outputFlux << "\t" << pop[i][0][j] << "/" << pop[i][1][j];
				}
			}
			outputFlux << std::endl;
		}else{
			std::cerr <<  "ERROR: cannot open the file " << outfileName << std::endl;
			exit(EXIT_FAILURE);
		}
	}	
}


void write_summary(const std::vector < std::vector < std::vector <float>>> & pop, const std::vector <float> & male_fitness, const std::vector <float> & female_fitness, const std::vector <size_t> & sexes, const unsigned int & N_diplo, const unsigned int & N_SA, const unsigned int & N_SC, const size_t & generation){
	// n_females / n_males_noInversion / n_males_inversion / recomb_females / recomb_males_noInversion / globExp_females / globExp_males_noInv / globExp_males_inv / fitness_females / fitness_males_noInv / fitness_males_inv
	size_t i(0);
	size_t j(0);

	size_t nFem(0);
	size_t nMal(0);
	size_t nMal_noInv(0);
	size_t nMal_inv(0);

	
	std::vector <float> fitness_fem;
	std::vector <float> fitness_mal_noInv;
	std::vector <float> fitness_mal_inv;
	
	std::vector <float> recomb_X;
	std::vector <float> recomb_Y_noInv;
	
	std::vector <float> globExp_X;
	std::vector <float> globExp_Y_noInv;
	std::vector <float> globExp_Y_inv;

	std::vector <float> ntrl_local_xp_X;
	std::vector <float> ntrl_genic_X;
	std::vector <float> SA_local_xp_X;
	std::vector <float> SA_genic_X;
	std::vector <float> SC_local_xp_X;
	std::vector <float> SC_genic_X;
	
	std::vector <float> ntrl_local_xp_Y_noInv;
	std::vector <float> ntrl_genic_Y_noInv;
	std::vector <float> SA_local_xp_Y_noInv;
	std::vector <float> SA_genic_Y_noInv;
	std::vector <float> SC_local_xp_Y_noInv;
	std::vector <float> SC_genic_Y_noInv;

	std::vector <float> ntrl_local_xp_Y_inv;
	std::vector <float> ntrl_genic_Y_inv;
	std::vector <float> SA_local_xp_Y_inv;
	std::vector <float> SA_genic_Y_inv;
	std::vector <float> SC_local_xp_Y_inv;
	std::vector <float> SC_genic_Y_inv;
	
	for(i=0; i<N_diplo; ++i){
		if( sexes[i]==0 ){ // if females
			fitness_fem.push_back(female_fitness[nFem]);	
			recomb_X.push_back(pop[i][0][1]);		recomb_X.push_back(pop[i][1][1]);
			globExp_X.push_back(pop[i][0][2]);		globExp_X.push_back(pop[i][1][2]);
			ntrl_local_xp_X.push_back(pop[i][0][3]);	ntrl_local_xp_X.push_back(pop[i][1][3]);
			ntrl_genic_X.push_back(pop[i][0][4]);		ntrl_genic_X.push_back(pop[i][1][4]);
			
			for(j=0; j<N_SA; ++j){
				SA_local_xp_X.push_back( pop[i][0][5 + 2*j]);	SA_local_xp_X.push_back( pop[i][1][5 + 2*j]);
				SA_genic_X.push_back( pop[i][0][6 + 2*j]);	SA_genic_X.push_back( pop[i][1][6 + 2*j]);
			}
			
			for(j=0; j<N_SC; ++j){
				SC_local_xp_X.push_back( pop[i][0][5 + 2*N_SA + 2*j]);	SC_local_xp_X.push_back( pop[i][1][5 + 2*N_SA + 2*j]);
				SC_genic_X.push_back( pop[i][0][6 + 2*N_SA + 2*j]);	SC_genic_X.push_back( pop[i][1][6 + 2*N_SA + 2*j]);
			}
			
			nFem++;
		}else{
			if( pop[i][0][1]*pop[i][1][1] == 0 ){
				// males with inversion
				nMal_inv++;
				fitness_mal_inv.push_back(male_fitness[i]);
				recomb_X.push_back(pop[i][1][1]);
			
				globExp_Y_inv.push_back(pop[i][0][2]);
				globExp_X.push_back(pop[i][1][2]);

				ntrl_local_xp_Y_inv.push_back(pop[i][0][3]);
				ntrl_genic_Y_inv.push_back(pop[i][0][4]);

				for(j=0; j<N_SA; ++j){
					SA_local_xp_Y_inv.push_back( pop[i][0][5 + 2*j]);
					SA_local_xp_X.push_back( pop[i][1][5 + 2*j]);
					
					SA_genic_Y_inv.push_back( pop[i][0][6 + 2*j]);
					SA_genic_X.push_back( pop[i][1][6 + 2*j]);
				}
				
				for(j=0; j<N_SC; ++j){
					SC_local_xp_Y_inv.push_back( pop[i][0][5 + 2*N_SA + 2*j]);
					SC_local_xp_X.push_back( pop[i][1][5 + 2*N_SA + 2*j]);
					
					SC_genic_Y_inv.push_back( pop[i][0][6 + 2*N_SA + 2*j]);
					SC_genic_X.push_back( pop[i][1][6 + 2*N_SA + 2*j]);
				}
			

			}else{
				// males without inversion
				nMal_noInv++;
				fitness_mal_noInv.push_back(male_fitness[i]);
				recomb_Y_noInv.push_back(pop[i][0][1]);
				recomb_X.push_back(pop[i][1][1]);
			
				globExp_Y_noInv.push_back(pop[i][0][2]);
				globExp_X.push_back(pop[i][1][2]);

				ntrl_local_xp_Y_noInv.push_back(pop[i][0][3]);
				ntrl_genic_Y_noInv.push_back(pop[i][0][4]);
				
				for(j=0; j<N_SA; ++j){
					SA_local_xp_Y_noInv.push_back( pop[i][0][5 + 2*j]);
					SA_local_xp_X.push_back( pop[i][1][5 + 2*j]);
					
					SA_genic_Y_noInv.push_back( pop[i][0][6 + 2*j]);
					SA_genic_X.push_back( pop[i][1][6 + 2*j]);
				}
				
				for(j=0; j<N_SC; ++j){
					SC_local_xp_Y_noInv.push_back( pop[i][0][5 + 2*N_SA + 2*j]);
					SC_local_xp_X.push_back( pop[i][1][5 + 2*N_SA + 2*j]);
					
					SC_genic_Y_noInv.push_back( pop[i][0][6 + 2*N_SA + 2*j]);
					SC_genic_X.push_back( pop[i][1][6 + 2*N_SA + 2*j]);
				}
			}
		nMal++;
		}
	}

	float recomb_X_avg(mean(recomb_X));
	float recomb_Y_noInv_avg(mean(recomb_Y_noInv));

	float fitness_fem_avg(mean(fitness_fem));
	float fitness_mal_noInv_avg(mean(fitness_mal_noInv));
	float fitness_mal_inv_avg(mean(fitness_mal_inv));
	
	float globExp_X_avg(mean(globExp_X));
	float globExp_Y_noInv_avg(mean(globExp_Y_noInv));
	float globExp_Y_inv_avg(mean(globExp_Y_inv));

	float ntrl_local_xp_X_avg(mean(ntrl_local_xp_X));
	float ntrl_genic_X_avg(mean(ntrl_genic_X));
	float SA_local_xp_X_avg(mean(SA_local_xp_X));
	float SA_genic_X_avg(mean(SA_genic_X));
	float SC_local_xp_X_avg(mean(SC_local_xp_X));
	float SC_genic_X_avg(mean(SC_genic_X));
	
	float ntrl_local_xp_Y_inv_avg(mean(ntrl_local_xp_Y_inv));
	float ntrl_genic_Y_inv_avg(mean(ntrl_genic_Y_inv));
	float SA_local_xp_Y_inv_avg(mean(SA_local_xp_Y_inv));
	float SA_genic_Y_inv_avg(mean(SA_genic_Y_inv));
	float SC_local_xp_Y_inv_avg(mean(SC_local_xp_Y_inv));
	float SC_genic_Y_inv_avg(mean(SC_genic_Y_inv));

	float ntrl_local_xp_Y_noInv_avg(mean(ntrl_local_xp_Y_noInv));
	float ntrl_genic_Y_noInv_avg(mean(ntrl_genic_Y_noInv));
	float SA_local_xp_Y_noInv_avg(mean(SA_local_xp_Y_noInv));
	float SA_genic_Y_noInv_avg(mean(SA_genic_Y_noInv));
	float SC_local_xp_Y_noInv_avg(mean(SC_local_xp_Y_noInv));
	float SC_genic_Y_noInv_avg(mean(SC_genic_Y_noInv));

	if( generation==0 ){
		std::cout << "generation ";
		std::cout << "nFem nMal_noInv nMal_inv ";
		std::cout << "recomb_X_avg recomb_Y_noInv_avg ";
		std::cout << "globExp_X_avg globExp_Y_noInv_avg globExp_Y_inv_avg ";
		std::cout << "fitness_fem_avg fitness_mal_noInv_avg fitness_mal_inv_avg ";
		std::cout << "ntrl_local_xp_X_avg ntrl_local_xp_Y_inv_avg ntrl_local_xp_Y_noInv_avg ";
		std::cout << "ntrl_genic_X_avg ntrl_genic_Y_inv_avg ntrl_genic_Y_noInv_avg ";
		std::cout << "SA_local_xp_X_avg SA_local_xp_Y_inv_avg SA_local_xp_Y_noInv_avg ";
		std::cout << "SA_genic_X_avg SA_genic_Y_inv_avg SA_genic_Y_noInv_avg ";
		std::cout << "SC_local_xp_X_avg SC_local_xp_Y_inv_avg SC_local_xp_Y_noInv_avg ";
		std::cout << "SC_genic_X_avg SC_genic_Y_inv_avg SC_genic_Y_noInv_avg" << std::endl;
	}

	std::cout << generation << " ";
	std::cout << nFem << " " << nMal_noInv << " " << nMal_inv << " ";
	std::cout << recomb_X_avg << " " << recomb_Y_noInv_avg << " ";
	std::cout << globExp_X_avg << " " << globExp_Y_noInv_avg << " " << globExp_Y_inv_avg << " ";
	std::cout << fitness_fem_avg << " " << fitness_mal_noInv_avg << " " << fitness_mal_inv_avg << " ";
	std::cout << ntrl_local_xp_X_avg << " " << ntrl_local_xp_Y_inv_avg << " " << ntrl_local_xp_Y_noInv_avg << " ";
	std::cout << ntrl_genic_X_avg << " " << ntrl_genic_Y_inv_avg << " " << ntrl_genic_Y_noInv_avg << " ";
	std::cout << SA_local_xp_X_avg << " " << SA_local_xp_Y_inv_avg << " " << SA_local_xp_Y_noInv_avg << " ";
	std::cout << SA_genic_X_avg << " " << SA_genic_Y_inv_avg << " " << SA_genic_Y_noInv_avg << " ";
	std::cout << SC_local_xp_X_avg << " " << SC_local_xp_Y_inv_avg << " " << SC_local_xp_Y_noInv_avg << " ";
	std::cout << SC_genic_X_avg << " " << SC_genic_Y_inv_avg << " " << SC_genic_Y_noInv_avg << std::endl;
}


float mean(const std::vector <float> & x){
	float res(0.0);
	size_t i(0);
	const size_t N_obs(x.size());

	if(N_obs < 10){
		return(0.0);
	}else{	
		for(i=0; i<N_obs; ++i){
			res += x[i];
		}
		return(res/N_obs);
	}
}


