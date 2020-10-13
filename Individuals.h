//This is Individuals.h

#pragma once

#include "Parameters.h"

#include <stdio.h>
#include <fstream>
#include <math.h>
#include <random>
#include <string>
#include <algorithm>


using namespace std;
//--------------------------------------------
class Individual {
public:
	Individual(int, int, int);//x, y, patchID
	~Individual() {};

	int ID;						//Individual ID
	int age;						//age of individual
	int age_Recruit;				//age at which they become adult and breed for the first time
	int patchID;					//ID of the current patch occupied by the individual
	int x;							//coordinate x of current patch
	int y;						//Coordinate y of current patch
	int new_patchID;			//when individuals disperse, ID of the new chosen patch
	int new_x;		//when individuals disperse, x coordinate of the new chosen patch
	int new_y;		////when individuals disperse, y coordinate of the new chosen patch

	string stage;		    // can be juvenile, pre-breeder or adult
	bool BreedPerf;			//whether the individual have successfully bred (can be true only for adults)
	bool disperser;			//whether the individual disperse
	float real_emigr;		//Actual emigration probability 
	float emigr_suc_interc;		//Emigration probability intercept for successful breeders  when emigr="pers_social" or"pers"
	float emigr_suc_slope;		//Emigration probability slope for successful breeders when emigr=="pers_social"
	float emigr_fail_interc;		//Emigration probability intercept for failed breeders
	float emigr_fail_slope;		//Emigration probability slope for failed breeders when emigr=="pers_social"
	int nb_prosp_patch;		//number of prospected patches
	bool alive;			//whether individual is alive

	//Functions
	void init_age(parameters, int);		//Initialize age at recruitment, age and status when individuals are created
	void init_prospect(parameters, int);		//Initialize nb of prospected patches (can be 0, 24 or randdomly assigned, depending on the settlement scenario)
	void init_emigr(parameters, int);			//Initialize the emigration probabilities depending on the emigration scenario when individuals are created
	void init_juv(parameters, int);		//Create a new individual when an adult produce one offspring
	void mutation_prospect(parameters, int);	//Mutation of prospecting in the case when settlement is "evol" or "evol_disp"
	void mutation_emigr(parameters, int);		//Mutation of emigration parameters when settl = "evol_disp"
	bool ind_death(parameters);				//Mortality depending on individual stage
	void age_ind(void);						//aging of individual if they survive after reproduction
	void disperse(parameters, float);		//Emigration probability depending on emigration scenario
	void prospecting(parameters);		//Create a vector of patches which are prospected
	void prospect_cost(parameters, int);		//Additional mortality cost due to prospecting
	void reinit_ind(int, int, int);			//Update individual position after dispersal/philopatry
	void OutInd(int rep, int years, std::ofstream *outind);		//Parameters to write in the output file

	std::vector<int>patchID_selected;		//Vector of patch objects

};

