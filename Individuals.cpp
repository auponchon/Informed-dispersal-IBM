//This is individual.cpp

#include "Individuals.h"


using namespace std;

std::random_device rd3;
std::mt19937 rdgenIND(rd3());



Individual::Individual(int xx, int yy, int idpatchh) {
	ID = 0;
	age = 0;						
	age_Recruit = 0;				
	x = xx;							
	y = yy;
	new_x = xx;
	new_y = yy;
	patchID = idpatchh;
	new_patchID = idpatchh;
	stage = "NA";					
	nb_prosp_patch = -9999;
	BreedPerf = false;
	disperser = false;
	real_emigr = 0;
	emigr_fail_interc = -9999;
	emigr_suc_interc = -9999;
	emigr_suc_slope = -9999;
	emigr_fail_slope = -9999;
	alive = true;

}



//Initialize age at recruitment, age and status
void Individual::init_age(parameters para, int id) {

	ID = id++;

	std::poisson_distribution<int> pois_recruit(para.MeanRecruitAge);
	age_Recruit = pois_recruit(rdgenIND);

	std::uniform_int_distribution<> unif_age(para.min_age, para.max_age);
	age = unif_age(rdgenIND);

	if (age < age_Recruit) {
		stage = "pre-breeder";
	}
	else { stage = "adult"; }

}

//Initialize emigration parameters
void Individual::init_emigr(parameters para, int id) {


	if (para.settl == 3) {		//3=="evol", when emigration is fixed

		if (para.emigr == 1) {		
		emigr_fail_interc = para.fixed_disp;
	}

		if (para.emigr != 1) {
			emigr_suc_interc = para.suc_disp_rate;
			emigr_fail_interc = para.fail_disp_rate;
		}
	}


	if (para.settl == 4) {			//4=="disp_evol", when emigration is evolving along with prospecting

		std::uniform_real_distribution<float> unif_emigr_intercept(0, 1);
		std::uniform_real_distribution<float> unif_emigr_fail_intercept(0, 1);
		std::uniform_real_distribution<float> unif_emigr_suc_intercept(0, 1);
		std::uniform_real_distribution<float> unif_emigr_slope(-1, 0);

		if (para.emigr == 1) {
			emigr_fail_interc = unif_emigr_intercept(rdgenIND);
	}

		if (para.emigr > 1) {
			emigr_fail_interc = unif_emigr_fail_intercept(rdgenIND);
			
			emigr_suc_interc =  unif_emigr_suc_intercept(rdgenIND);
		}

		if (para.emigr == 3) {
			emigr_suc_slope = unif_emigr_slope(rdgenIND);
			emigr_fail_slope = unif_emigr_slope(rdgenIND);
		}
	}
}


//Create a new individual when an adult produce one offspring
void Individual::init_juv(parameters para, int id) {

	std::poisson_distribution<int> pois_recruit(para.MeanRecruitAge);
	
	ID = id;
	age = 0;
	age_Recruit = pois_recruit(rdgenIND);
	real_emigr = 0;
	stage = "juvenile";
	BreedPerf = false;
}


//Initialize nb of prospected patches
void  Individual::init_prospect(parameters para, int id) {
	std::uniform_int_distribution<> unif_patch(0, para.max_prosp_patch);

	if (para.settl == 1) nb_prosp_patch = 0;		//no prospecting = random dispersal
	if (para.settl == 2) nb_prosp_patch = para.max_prosp_patch;		//fully informed dispersal (all avalaible patches prospected)
	if (para.settl == 3 || para.settl == 4) nb_prosp_patch = unif_patch(rdgenIND);		
}


//Mutation of prospecting in the case when settlement is "evol"
void Individual::mutation_prospect(parameters para, int id) {

	if (para.settl == 3 || para.settl == 4) {

	std::uniform_real_distribution<float> rate_mut(0, 1);
	std::uniform_int_distribution<> change_col(0, 1);

		double mutation_rd;
		mutation_rd = rate_mut(rdgenIND);

		if (mutation_rd < para.mutation_rate) {
			int colchange = change_col(rdgenIND);			//decides if it's gonna be + or - 1 COL
			if (colchange == 0) { nb_prosp_patch = nb_prosp_patch - 1; }
			if (colchange == 1) { nb_prosp_patch = nb_prosp_patch + 1; }
			if (nb_prosp_patch > para.max_prosp_patch) { nb_prosp_patch = para.max_prosp_patch; }
			if (nb_prosp_patch < 0) { nb_prosp_patch = 0; }

		}
	}
}


//Mutation of eigration in the case when settlement is "evol"
void Individual::mutation_emigr(parameters para, int id) {
	
	if (para.settl == 4) {

		std::uniform_real_distribution<float> rate_mut(0, 1);
		std::normal_distribution<float> float_normal(0, 0.1);

		float mutation_rd1, mutation_rd2, mutation_rd3, mutation_rd4;
		mutation_rd1 = rate_mut(rdgenIND);
		mutation_rd2 = rate_mut(rdgenIND);
		mutation_rd3 = rate_mut(rdgenIND);
		mutation_rd4 = rate_mut(rdgenIND);

		if (mutation_rd1 < para.mutation_rate) {
			float rd_disp1 = float_normal(rdgenIND);
			emigr_fail_interc += rd_disp1;
		}

		if (para.emigr != 1) {
			if (mutation_rd2 < para.mutation_rate) {
				float rd_disp2 = float_normal(rdgenIND);
				emigr_suc_interc += rd_disp2;

			}
		}

		if (para.emigr == 3) {
			if (mutation_rd3 < para.mutation_rate) {
				float rd_disp3 = float_normal(rdgenIND);
				emigr_suc_slope += rd_disp3;
				emigr_suc_slope = rate_mut(rdgenIND);
			}

			if (mutation_rd4 < para.mutation_rate) {
				float rd_disp4 = float_normal(rdgenIND);
				emigr_fail_slope += rd_disp4;
				emigr_fail_slope = rate_mut(rdgenIND);
			}
		}

	}

}


//Mortality depending on individual stage
bool Individual::ind_death(parameters para) {
	uniform_real_distribution<double> unif_rates(0, 1);

	double rdmSurv = unif_rates(rdgenIND);
	double survstage=0;

	if (stage == "juvenile") { survstage = para.survJuv; }
	else if (stage == "pre-breeder") { survstage = para.survImmat; }
	else if (stage == "adult") { survstage = para.survAdult; }

	if (rdmSurv > survstage) {   // if adults survive
		alive = false;
	}

	return alive;
}

//Aging of individual
void Individual::age_ind(void) {
	age++;

	if (age == 1) {
		stage = "pre-breeder";
	}
	if (age == age_Recruit) {
		stage = "adult";
	}
	if (age > 100) {
		alive = false;
	}

}


//Dispersal probability depending on emigration scenario
void Individual::disperse(parameters para, float patch_LBS) {

	//Define value for dispersal
	float rd_dispers;

	if (para.settl != 4) {

		if (para.emigr == 1) {
			real_emigr = para.fixed_disp;
		}

		if (para.emigr == 2) {
			if (BreedPerf == true) {
				real_emigr = para.suc_disp_rate;
			}
			else { real_emigr = para.fail_disp_rate; }
		}

		if (para.emigr == 3) {
			if (BreedPerf == false) {
				real_emigr = para.fail_disp_rate + (para.suc_disp_rate - para.fail_disp_rate)*patch_LBS;
			}
			else { real_emigr = para.suc_disp_rate; }
		}
	}

	if (para.settl == 4) {

		if (para.emigr == 1) {
			real_emigr = emigr_fail_interc;
		}

		if (para.emigr == 2) {
			if (BreedPerf == false) {
				real_emigr = emigr_fail_interc;
			}
			else {
				real_emigr = emigr_suc_interc;
			}
		}

		if (para.emigr == 3) {
			if (BreedPerf == false) {
				real_emigr = emigr_fail_interc + emigr_fail_slope*patch_LBS;
			}
			if (BreedPerf == true) {
				real_emigr = emigr_suc_interc + emigr_suc_slope*patch_LBS;
		}
		}

	}


		//Disperser or not?
	uniform_real_distribution<float> proba(0, 1);

		rd_dispers = proba(rdgenIND);

		if (rd_dispers < real_emigr ) {
			disperser = true;
			new_x = -999;
			new_y = -999;
			new_patchID = -999;

		}
	

		if (rd_dispers >= real_emigr) {
			disperser = false;
			real_emigr = 0.0;
			new_x = x;
			new_y = y;
			new_patchID = patchID;

	}
}


//Create a vector of patches which are prospected
void Individual::prospecting(parameters para) {
	if (nb_prosp_patch > 1) {
		for (int d = 0; d < para.N_patches; d++) {
			if (d != patchID) {
				patchID_selected.push_back(d);
			}
		}

		shuffle(patchID_selected.begin(), patchID_selected.end(), rdgenIND);
	}
}


//Additional mortality cost due to prospecting
void Individual::prospect_cost(parameters para, int j) {
	if (para.prosp_cost > 0) {

	uniform_real_distribution<double> unif_rates(0, 1);

	double	rdm_prosp_cost = unif_rates(rdgenIND);	
		if (rdm_prosp_cost < para.prosp_cost * nb_prosp_patch) {
			alive = false;
			disperser = false;
			new_patchID = -999;
			new_x = -999;
			new_y = -999;
		}
	}

}	

//reinitialization of individual featuress at the end of the year
void Individual::reinit_ind(int a, int b, int c) {
	patchID = a;
	x = b;
	y = c;
	new_patchID = a;
	new_x = b;
	new_y = c;
	real_emigr = -999;
	disperser = false;
	BreedPerf = false;
}



//Parameters to write in the output file
void Individual::OutInd(int rep, int years, std::ofstream *outind) {
	if (alive == true && stage=="adult") {

		*outind << rep << " "
			<< years << " "
			<< patchID << " "
			<< new_patchID << " "
			<< ID << " "
			<< age << " "
			<< age_Recruit << " "
			<< BreedPerf << " "
			<< disperser << " "
			<< real_emigr << " "
			<< emigr_fail_interc << " " 
			<< emigr_suc_interc << " " 
			<< emigr_fail_slope << " " 
			<< emigr_suc_slope << " " 
			<< nb_prosp_patch << std::endl;
	}
}