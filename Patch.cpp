#include "Patch.h"
#include "Individuals.h"

using namespace std;

//random generators

std::random_device rd2;
std::mt19937 rdgenPatch(rd2());//constructor: requires x-coordinate, y-coordinate, and number of patch

Patch::Patch(int xx, int yy, int IDPATCH) {

	//Parameters of the grid
	IDPatch = IDPATCH;
	x = xx;
	y = yy;
	env_qual = -9999;
	K = -9999;
	product = 0.0;

	//Counters for the population
	N_adults_init = 0;
	N_immat_init = 0;
	N_popsize = 0;
	N_adults = 0;
	N_succ =0;
	N_fail =0;
	N_juvs=0;
	N_immat_surv = 0;
	N_adult_surv = 0;
	N_juv_surv = 0; 
	N_dead_adults = 0;
	N_dead_juv = 0;
	N_dead_immat = 0;
	N_NewAdults = 0;
	N_immat_final = 0;
	N_adults_predisp = 0;
	N_disp = 0;
	N_stay = 0;
	N_immigr = 0;
	N_Disp_dead = 0;
	N_adults_postdisp = 0;
	N_popsize_surv = 0;
	N_popsize_final = 0;

	//Parameters for the reproduction
	dens_depdce = 0.0;
	ExpOff = 0.0;		//Expected number of offspring after density-dependance

}

//Add individual objects on a patch
void Patch::add_ind(Individual indi) {
	pop_init_vec.push_back(indi);
}

//Get LBS from number of successful and failed breeders
float Patch::get_LocBreedProd(int a, int b) {
	if (b == 0) return 0.0;
	else {
		float r;
		r = (float)a / (float)b;
		return r;
	}
}


//initialize patch quality at t0 and calculate density-dependence factor
void Patch::init_patch_qual(parameters para) {
	std::normal_distribution<float> init_enviro_norm(0.0, para.SD0_enviro);

	IDPatch = (para.NCOL*y) + x;

	env_qual = init_enviro_norm(rdgenPatch);
	if (env_qual < (-1.0)) { env_qual = (-0.9999999999999999); }
	if (env_qual > 1.0) { env_qual = 1.0; }

	K = para.K_baseline * (env_qual + 1.0);

}

//initialize population of individuals on a patch
void Patch::init_pop_onPatch(parameters para, int id) {
		
	for (int n = 0; n < para.N_init; n++) {
		Individual inito = Individual(x, y, IDPatch);
		
		inito.ID = id++;
		inito.init_age(para, n);
		N_popsize_init++;
	}
}


//recalculate environmental quality for a patch based on previous year
void Patch::enviro_quality(parameters para) {

	if (para.enviro == 1) {
		std::normal_distribution<float> enviro_norm(env_qual, para.SD1_enviro);

		env_qual = enviro_norm(rdgenPatch);
	}

	else { uniform_real_distribution<float> enviro_stochast(-1, 1);
			env_qual = enviro_stochast(rdgenPatch);
	}

	if (env_qual < (-1.0)) { env_qual = (-0.9999999999999999); }
	if (env_qual > 1.0) { env_qual = 1.0; }

	K = para.K_baseline * (env_qual + 1.000000);
	
}



//Produce new individuals on the patch
int Patch::reproduction(int j, int id, parameters para) {
	
	//Calculate mean nb of offspring
		dens_depdce = (pow(para.lambda, (1.0 / para.b)) - 1.0) /K;
		ExpOff = para.lambda * pow((1 + dens_depdce * N_adults_init), -para.b);

	
	//Create juveniles based on expected nb of offspring for the patch
	int NOff;
		
			if (ExpOff > 0)
			{
				poisson_distribution<int> pois_off(ExpOff);
				NOff = pois_off(rdgenPatch);
			}
			else { NOff = 0; }


			if (NOff > 0) {
				N_succ++;
				pop_init_vec[j].BreedPerf = true;

				for (int w = 0; w < NOff; w++) {

					id++;
					Individual indi = pop_init_vec[j];		//new ind is created in patch
					indi.init_juv(para, id);
					indi.mutation_prospect(para, id);
					indi.mutation_emigr(para, id);
					N_juvs++;
					
					indi.ind_death(para);

						if (indi.alive == false) {
							N_dead_juv++;
						}

						if (indi.alive == true) { N_juv_surv++;
						pop_init_vec.push_back(indi);
						}
				}
			}	//end of if individual is successful

			if (NOff < 1) {
				 N_fail++;
				pop_init_vec[j].BreedPerf = false;
				}						//otherwise, ind is a failed breeder

	product = get_LocBreedProd(N_succ, N_succ+N_fail);

	return id;
}


//Stage-dependent mortality on the patch
void Patch::mortality(parameters para, int j) {

		pop_init_vec[j].ind_death(para);		//calls death function of individual

	if (pop_init_vec[j].alive == false) {				//if individual is dead, then it is deleted from patch

    if (pop_init_vec[j].stage == "pre-breeder") {
			N_dead_immat++;		
		}

		else  if (pop_init_vec[j].stage == "adult") {
			N_dead_adults++;

		}
	}


}

//Recruitment: pre-breeders become adults
void Patch::recruitment(int j) {

	if (pop_vec[j].age == pop_vec[j].age_Recruit) {
			N_NewAdults++;
		}
}


//Get individual dispersing proba
void Patch::dispersal(parameters para,int j) {

		pop_vec[j].disperse(para, product);

		if (pop_vec[j].disperser == true) {
			N_disp++;
		}

		else  N_stay++;
	}




//reinitialization of the population at the end of the year
void Patch::reinit(void) {

	pop_vec.clear();

	N_adults_init = N_adults_postdisp;
	N_immat_init = N_immat_final;
	N_popsize_init = N_popsize_final;

}


//reinitialization of counters at the end of the year
void Patch::clear_pop (void){

	N_succ = 0;
	N_fail = 0;
	N_juvs = 0;
	N_adults = 0;
	N_popsize = 0;
	N_immat_surv = 0;
	N_adult_surv = 0;
	N_juv_surv = 0;
	N_dead_adults = 0;
	N_dead_juv = 0;
	N_dead_immat = 0;
	N_NewAdults = 0;
	N_immat_final = 0;
	N_adults_predisp = 0;
	N_disp = 0;
	N_stay = 0;
	N_immigr = 0;
	N_Disp_dead = 0;
	N_adults_postdisp = 0;
	N_popsize_surv = 0;
	N_popsize_final = 0;
	succrate = 0;
	ExpOff = -999;

}

//Parameters to write in the output file
void Patch::OutPop(int rep, int years,  std::ofstream *outpop) {

	*outpop << rep << " "
		<< years << " "
		<< IDPatch << " "
		<< env_qual << " "
		<< K << " "
		<< product << " "
		<< ExpOff << " "
		<< N_adults_init << " "
		<< N_immat_init << " "
		<< N_succ << " " 
		<< N_fail << " "
		<< N_juvs << " "
		<< N_adults_predisp << " "
		<< N_disp << " " 
		<< N_stay << " "
		<< N_immigr << " "
		<< N_adults_postdisp << " "
		<< N_immat_final << " "
		<< std::endl;
}
