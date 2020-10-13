//This is Patch.h

#pragma once

#include <stdio.h>
#include <vector>
#include <random>

#include "Individuals.h"
#include "Parameters.h"

using namespace std;

class Patch {
public:
	Patch(int, int, int);   //x,y,IDPatch

///Parameters for the grid
	int IDPatch;		//ID of the patch
	int x;			//x coordinate
	int y;		//y coordinate
	float env_qual;		//Environmental quality (bteween -1 and 1)
	float K;			//Carrying capacity depending on environmental quality
	float product;		//Breeding success: number of successful breeders / number of failed breeders


	//Parameters for the population installed on each patch
	int N_adults_init;		//Number of adults at the beginning of the annual cycle, before reproduction
	int N_immat_init;		//Number of pre-breeders at the beginning of the annual cycle, before reproduction
	int N_popsize_init;		//Total number of individuals (pre-breeders and adults) at the beginning of the annual cycle, before reproduction
	int N_adults;			//Number of adults for reproduction
	int N_succ;				//Number of individuals that successfullly produce at least one offpsring
	int N_fail;				//Number of individuals that did not produce any offspring
	int N_juvs;				//Number of juveniles (new-born) produced during reproduction
	int N_popsize;		//Number of total individuals after reproduction
	int N_dead_adults;		//Number of adults which die after reproduction
	int N_dead_immat;	//Number of pre-breeders which die after reproduction
	int N_dead_juv;		//Number of juveniles which die after reproduction
	int N_immat_surv;		//Number of pre-breeders which survive after reproduction
	int N_adult_surv;		//Number of pre-breeders which survive after reproduction
	int N_popsize_surv;		//Total number of individuals surviving after reproduction
	int N_juv_surv;		//Number of juveniles which survive after reproduction
	int N_NewAdults;	//Number of pre-breeders which reach age at recruitment and become adults
	int N_immat_final;		//Number of pre-breeders at the end of the annual cycle
	int N_adults_predisp;	//Number of adults in a patch before dispersal (= N_adult_surv)
	int N_disp;		//Number of adults dispersing to a new patch 
	int N_stay;		//Number of adult remaining in their current patch
	int N_immigr;	//Number of immigrants
	int N_adults_postdisp;		//Total number of adults after dispersal (N_stay + N_immigr - N_disp)
	int N_Disp_dead;		//Number of dead dispersers when prospecting involves costs
	int N_popsize_final;		//Number of total individuals (adults + pre-breeders) at the end of the annual cycle

	

	//Parameters for reproduction
	float dens_depdce;		//Density-dependence equation
	float succrate;		//Breeding success: Number of successful breeders / Number of failed breeders
	float ExpOff;					//Expected number of offspring after density-dependance

	std::vector<Individual> pop_init_vec, pop_vec;		//Vectors containing individual objects


//Functions to initialise and make the population live and die on the patches
void add_ind( Individual);		//Add individual objects on a patch
void init_patch_qual(parameters);		//initialize patch quality at t0
void init_pop_onPatch(parameters, int);		//initialize individual features on a patch
void enviro_quality(parameters);		//recalculate environmental quality for a patch based on previous year
int reproduction(int, int, parameters);			//Produce new individuals on the patch
float get_LocBreedProd(int, int);			//Get LBS from number of successful and failed breeders
void mortality(parameters, int);		//Stage-dependent mortality on the patch
void recruitment(int);		//Recruitment: pre-breeders become adults
void dispersal(parameters, int);		//Count of individuals dispersing or staying in the patch
void reinit(void);		//reinitialization of the population at the end of the year
void clear_pop(void);		//reinitialization of counters

void OutPop(int rep, int years,  std::ofstream *outpop);	//population parameters to write in the output file



};