//This is parameters.cpp

#include "Parameters.h"


parameters::parameters() {
	ScenNb = 1;				//simulation scenario number
	YEARS = 20001;				//number of years of simulations
	REP = 10;				//number of replicates
	emigr = 2;		// Emigration strategy: can be 1 =random, 2 ="pers" or "3" =pers_social
	settl = 4;		// Settlement strategy: "3"="evol" (only prospecting evolves) or 4="emigr_evol" (both emigration and prospecting evolve)
	selection = 2;	// Patch selection strategy: can be 3="deterministic", 2="accurate", 1="inaccurate"
	cost = true;			//If cost=false, no prospecting cost when prospecting; if cost=true, prospecting cost
	
	if (selection == 2) {
		softmax = 50;
	}    //coefficient for the accurate selection process (only works when selection=2)
	else { softmax = 0; }

	//Parameters for the grid
	NROW = 5;					//Number of cells on x
	NCOL = 5;					// Number of cells in y
	N_patches = NROW * NCOL;				//Total number of patches
	max_prosp_patch = (NROW * NCOL) - 1;		//Maximum number of prospected patches

		//Parameters for the stochastic but auto-correlated environment
	SD_enviro = 1;		// SD of environmental quality
	ac = 0.8;					//auto-correlation coefficient
								
	//Carrying capacity
	N_init=100;			//Number of initial number of individuals in a patch
	K_baseline = 100.0;		//Number of individuals in a medium environment
	K_max = 1000;			//maximum number of individuals in a patch
	b = 1;

	//Mutation in evolutionary model

	if (settl == 3 || settl == 4) {
		mutation_rate = 0.01;	}
	else { mutation_rate = 0; }

	//Cost of prospecting per patch
	if (cost) {
		prosp_cost = 0.01;
	}	
	else { prosp_cost = 0.0; }

	//Settings for ouput files
	seq_out_pop = 1000;			//Output file for populations written every *** years
	seq_out_ind = 20000;			//Output file for individuals written every *** years

	
	//////////////////////////
	// DEMOGRAPHIC PARAMETERS
	//Life history traits
	min_age = 1;				//minimum age for individuals at initialization
	max_age = 15;				//maximum age for individuals at initialization
	MeanRecruitAge = 5;			//Mean age at recruitment
	fixed_disp = 0.5;			// Emigr proba when emigr="random"
	fail_disp_rate = 0.5;		//Emigr proba for failed breeders when emigr="pers"
	suc_disp_rate = 0.05;		//Emigr proba for successful breeders when emigr="pers"
	fail_disp_slope = -0.45;	//Slope of emigr proba for failed breeders when emigr="pers_social"
	suc_disp_slope = 0;			//Slope of emigr proba for successful breeders when emigr="pers_social"

	//Constants for density-dependant breeding success 
	lambda = 2; // Mean number of offspring

	//Survival rates
	survJuv = 0.6;	//Survival of juveniles
	survImmat = 0.7;	//Survival of pre-breeders
	survAdult = 0.85;	//Survival of adults
	
}


//Parameters to write in the parameters output file
void parameters::outpara(std::string name) { //This defines what parameters are going to be output from the model
	std::ofstream out_param;
	out_param.open(name.c_str());

	out_param << "sim " << REP << std::endl;
	out_param << "gen " << YEARS << std::endl;
	out_param << "NROW " << NROW << std::endl;
	out_param << "NCOL " << NCOL << std::endl;
	out_param << "Initial K " << K_baseline << std::endl;
	out_param << "Enviro_init_SD " << SD_enviro << std::endl;
	out_param << "Autocor " << ac << std::endl;
	out_param << "Emigr " << emigr << std::endl;
	out_param << "Settl " << settl << std::endl;
	out_param << "Selection " << selection << std::endl;
	out_param << "Mutation rate " << mutation_rate << std::endl;
	out_param << "Prospecting cost " << prosp_cost << std::endl;

	out_param.close();
}



parameters::~parameters()
{
}
