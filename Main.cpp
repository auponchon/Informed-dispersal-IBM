//This is main.cpp

#include <iostream>
#include "MainHeader.h"


//If code has to be run on a HPC
#if CLUSTER
int main(int argc, char* argv[])
{
	// Get the current directory.
	char* buffer = getcwd(NULL, 0);
	dir = buffer;
	free(buffer);
	dir = dir + "/"; //Current directory path
	dirout = dir + "Outputs/"; //Outpus folder path


	para.ScenNb			 = std::atoi(argv[1]);
	para.REP			 = std::atoi(argv[2]);
	para.YEARS	 	  	 = std::atoi(argv[3]);
	para.emigr			 = std::atoi(argv[4]);
	para.settl			 = std::atoi(argv[5]);
	para.selection		 = std::atoi(argv[6]);
	para.enviro			 = std::atoi(argv[7]);
	para.mutation_rate	       = std::atof(argv[8]);
	para.prosp_cost		 = std::atof(argv[9]);
	para.seq_out_pop	       = std::atoi(argv[10]);
	para.seq_out_ind	       = std::atoi(argv[11]);


	RunModel();

	cout << "Simulations completed, YOU ROCK IT, YOU'RE THE BOSS" << endl;

	return 0;
}
#else


//If code is run from a normal computer
int _tmain(int argc, _TCHAR* argv[])
{

	char* buffer = _getcwd(NULL, 0);
	dir = buffer;
	free(buffer);
	dir = dir + "\\"; //Current directory path
	dirout = dir;
	clock_t tStart = clock();

	RunModel();

	printf("Time taken: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
	cout << "Simulations completed, YOU ROCK IT, YOU'RE THE BOSS" << endl;
	system("pause");

return 0;

}	//end of main fuction
#endif

//Functions to convert integers or float to string
//---------------------------------------------------------------------------
string Int2Str(int x)
{
	ostringstream o;
	if (!(o << x)) return "ERROR";
	return o.str();
}
//---------------------------------------------------------------------------
string Float2Str(double x)
{
	ostringstream o;
	if (!(o << x)) return "ERROR";
	return o.str();
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Main model function
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void RunModel(void) {

//	for (SimN = 0; SimN < para.ScenNb; SimN++) {

		
		//generate headers in the output files
		outpop_header();
		outind_header();


		////////REPLICATE SIMULATION
		for (rep = 0; rep < para.REP; rep++) {
			id = 0;

			Initialise();

			///YEARS WITHIN REPLICATES
			for (years = 0; years < para.YEARS; years++) {

				cout << "Replicate " << rep << " year " << years << endl;

				Enviro_Quality();
				Reproduction();
				Mortality();
				Age_and_Recruit();
				Dispersal();
				Get_OutInd(para);
				Settlement();
				Get_OutPop(para);
				Reinitialisation();


			}	//end of years

		}	//end of replicates

		if (outpop.is_open())   outpop.close();   //Checks if the population output file is open and closes it
		if (outind.is_open())   outind.close();
		if (outparam.is_open()) outparam.close();
//	}		//end of simulation Nb


}		//end of running model





/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Initialisation
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void Initialise(void) {


	string name_param = dirout + "Sc" + Int2Str(para.ScenNb) + "_Em_" + Int2Str(para.emigr) + "_Settl_" + Int2Str(para.settl) + "_Select_" + Int2Str(para.selection)  
		+ "_Env_" + Float2Str(para.SD1_enviro) + "_Cost_" + Float2Str(para.prosp_cost) + "_Param.txt";
	para.outpara(name_param);
	
	grid.clear();
	id = 0;

//Initialisation of patches and populations 

	for (y = 0; y < para.NROW; y++) {
		for (x = 0; x < para.NCOL; x++) {
			

			//Create a vector of patches
			idpatch = (para.NCOL*y) + x;		//get patchID
			Patch backup = Patch(x, y, idpatch);		//create patch object 
			backup.init_patch_qual(para);			//initialize patch i

			backup.init_pop_onPatch(para, id);		//initialize population on the patch
			
			//Create a vector of individuals on each patch
			for (int n=0; n < backup.N_popsize_init; n++) {
				id++;
				Individual ind = Individual(x,y, idpatch);
				ind.init_age(para, id);
				ind.init_prospect(para, id);
				ind.init_emigr(para, id);
				backup.add_ind(ind);

				if (ind.stage == "pre-breeder") { backup.N_immat_init++; }
				else { backup.N_adults_init++; }

				}

			grid.push_back(backup);

			}	//end of x
		}//end of y 


}		//end of Initialisation function



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//YEARLY VARIATION OF ENVIRONMENTAL QUALITY OVER THE WHOLE GRID
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Enviro_Quality() {

	for (int p=0; p<para.N_patches;p++) {
			grid[p].enviro_quality(para);
		}
}//end of enviro_qual function



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//REPRODUCTION
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Reproduction() {
	for (int p = 0; p < para.N_patches; p++) {
	
		if( grid[p].N_popsize_init >0){
				for (int n = 0; n < grid[p].N_popsize_init; n++) {

					if (grid[p].pop_init_vec[n].stage == "adult") {
						id = grid[p].reproduction(n, id, para);
					} //end of if ind is adult
			}	//end of n in individuals

			grid[p].N_adults = grid[p].N_succ + grid[p].N_fail;

			grid[p].N_popsize = grid[p].N_adults + grid[p].N_immat_init + grid[p].N_juv_surv;


			}		//end of if pop is empty

		}		//end of grid
}		//end of reproduction



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//MORTALITY
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Mortality() {

	for (int p = 0; p < para.N_patches; p++) {

			if (grid[p].N_popsize > 0) {

	
				for (int n = 0; n < grid[p].N_popsize; n++) {

					if (grid[p].pop_init_vec[n].stage != "juvenile") {
					grid[p].mortality(para, n);
				}

					if (grid[p].pop_init_vec[n].alive == true) {
						grid[p].N_popsize_surv++;
						grid[p].pop_vec.push_back(grid[p].pop_init_vec[n]); }


				} //end of looping over individuals

				grid[p].N_adult_surv = grid[p].N_adults - grid[p].N_dead_adults;
				grid[p].N_immat_surv = grid[p].N_immat_init - grid[p].N_dead_immat;

			}	//end of if pop not empty

			grid[p].pop_init_vec.clear();

	}	//end of p
}	//end of mortality


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//AGING AND RECRUITMENT
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Age_and_Recruit(void) {
	for (int p = 0; p < para.N_patches; p++) {

		if (grid[p].N_popsize_surv > 0) {
			for (int n = 0; n < grid[p].N_popsize_surv; n++) {

				grid[p].pop_vec[n].age_ind();
				grid[p].recruitment(n);
			}

			grid[p].N_adults_predisp = grid[p].N_NewAdults + grid[p].N_adult_surv;
			grid[p].N_immat_final = grid[p].N_immat_surv + grid[p].N_juv_surv - grid[p].N_NewAdults;

		}
	}
}

	

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//DISPERSAL
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Dispersal(void) {
	int  zz;
	int new_patchID_temp = -999;
	std::uniform_real_distribution<float> rate_rd(0, 1);

	for (int p = 0; p < para.N_patches; p++) {

		if (grid[p].N_popsize_surv > 0) {	

			for (int n = 0; n < grid[p].N_popsize_surv; n++) {

				if (grid[p].pop_vec[n].stage == "adult") {

					grid[p].dispersal(para, n);

				}

					if (grid[p].pop_vec[n].disperser == true) {

						//1) IF INDIVIDUALS PROSPECT
						if (grid[p].pop_vec[n].nb_prosp_patch > 1) {
							grid[p].pop_vec[n].prospecting(para);


							float productivity = 0.0;

							//DETERMINISTIC PROCESS IN PATCH SELECTION
							if (para.selection == 3) {		//if patch selection is fixed

								for (int k = 0; k < grid[p].pop_vec[n].nb_prosp_patch; k++) {					// k = Number of prospected colonies
									zz = grid[p].pop_vec[n].patchID_selected[k];

									if (grid[zz].product >= productivity) {
										productivity = grid[zz].product;
										new_patchID_temp = grid[zz].IDPatch;

									}	//end of finding max of LBS in prospected colonies					
								}		//end of k
							}		//end of fixed choice


							//LOTTERY OR SOFTMAX PROCESS FOR PATCH SELECTION
							if (para.selection == 1 || para.selection == 2) {		//if patch selection is inaccurate (1) or accurate (2)
								float cum_sumtot = 0.0;


								//Calculate the total cumulated sum of productivity
								for (int k = 0; k < grid[p].pop_vec[n].nb_prosp_patch; k++) {					// k = Number of prospected colonies
									zz = grid[p].pop_vec[n].patchID_selected[k];

									if (para.selection == 2)
									{
										cum_sumtot += exp(para.softmax * grid[zz].product);
									}
									else { cum_sumtot += grid[zz].product; }

									patch_prosp_vec.push_back(grid[zz]);

								}		//end of k


								//Sort the patch by decreasing productivity and select one

								float rdm_cumsum = rate_rd(rdgen);
								float cumsum = 0.0;

								int l = 0;

								while (cumsum < rdm_cumsum) {

									if (para.selection == 2) { cumsum += (exp(para.softmax * patch_prosp_vec[l].product)) / cum_sumtot; }
									else { cumsum += patch_prosp_vec[l].product / cum_sumtot; }

									new_patchID_temp = patch_prosp_vec[l].IDPatch;
									productivity = patch_prosp_vec[l].product;
									l++;
								}		//end of patch selection				

							}//end of if reading==fixed or lottery



							grid[p].pop_vec[n].patchID_selected.clear();
							patch_prosp_vec.clear();

							//If prospecting entails some mortality cost
							grid[p].pop_vec[n].prospect_cost(para, n);

							// 2) Getting the new coordinates of the patch selected
							grid[p].pop_vec[n].new_patchID = new_patchID_temp;
							grid[p].pop_vec[n].new_x = grid[new_patchID_temp].x;
							grid[p].pop_vec[n].new_y = grid[new_patchID_temp].y;
						}	//end of if individuals prospect 


						//2) if individuals do not prospect, they disperse randomly to a new breeding patch other than their current one		
						if (grid[p].pop_vec[n].nb_prosp_patch < 2) {
							std::uniform_int_distribution<> unif_patch(0, para.max_prosp_patch);

							//Get the ID of a random patch other than the current one
							do {
								new_patchID_temp = unif_patch(rdgen);
							} while (new_patchID_temp == grid[p].IDPatch);

							grid[p].pop_vec[n].new_patchID = new_patchID_temp;
							grid[p].pop_vec[n].new_x = grid[new_patchID_temp].x;
							grid[p].pop_vec[n].new_y = grid[new_patchID_temp].y;
						}		//end of random dispersal	
					}		//end of if individuals disperse

					//3) Individuals do not disperse
					if (grid[p].pop_vec[n].disperser == false) {
						grid[p].pop_vec[n].real_emigr = 0;
						grid[p].pop_vec[n].new_patchID = grid[p].IDPatch;
						grid[p].pop_vec[n].new_x = grid[p].x;
						grid[p].pop_vec[n].new_y = grid[p].y;
					}

			}	//end of n in individuals

		}	//end of if pop empty
	}	//end of patches
}	//end of dispersal fdunction



void Settlement(void) {
		
	for (int p = 0; p < para.N_patches; p++) {

			for (int n = 0; n < grid[p].N_popsize_surv; n++) {
				if (grid[p].pop_vec[n].alive == true) {

					if (grid[p].pop_vec[n].stage == "adult") {

						//SETTLEMENT IN NEW SITE
						if (grid[p].pop_vec[n].disperser == true) {
							grid[grid[p].pop_vec[n].new_patchID].pop_init_vec.push_back(grid[p].pop_vec[n]);
							grid[grid[p].pop_vec[n].new_patchID].N_immigr++;
							grid[grid[p].pop_vec[n].new_patchID].N_adults_postdisp++;
							
						}		//end of if individuals disperse

						else  {
							grid[p].pop_init_vec.push_back(grid[p].pop_vec[n]);
							grid[p].N_adults_postdisp++;
						}
					}		//end of looping on adults
						
					//Keep pre-breeders in the same patch
					else  {
						grid[p].pop_init_vec.push_back(grid[p].pop_vec[n]);
					}

					}	//end of if individuals alive

					//Removal of individuals which died prospecting
				if (grid[p].pop_vec[n].alive == false && grid[p].pop_vec[n].stage == "adult")  {grid[p].N_Disp_dead++;	}

			}		//end of individuals
	}		//end of p
}	//end of function



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//REINITIALIZATION
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Reinitialisation(void) {
	for (int p = 0; p < para.N_patches; p++) {

		grid[p].N_popsize_final = grid[p].N_immat_final + grid[p].N_adults_postdisp;
	
		for (int n = 0; n < grid[p].N_popsize_final; n++) {
			grid[p].pop_init_vec[n].reinit_ind(grid[p].IDPatch, grid[p].x, grid[p].y);

			}

			grid[p].reinit();
			grid[p].clear_pop();

		}
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Outputs
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//Write headers in population output file
void outpop_header(void) {
	std::string name_pop;
	name_pop = dirout + "Sc" + Int2Str(para.ScenNb) + "_Em_" + Int2Str(para.emigr) + "_Settl_" + Int2Str(para.settl) + "_Select_" + Int2Str(para.selection)  + "_Env_" + Float2Str(para.SD1_enviro) + "_Cost_" + Float2Str(para.prosp_cost) + "_pop.txt";
	outpop.open(name_pop.c_str());
	outpop << "sim gen PatchID env.qual K succrate ExpOff Nadult Nimmat Nsuc Nfail Njuv NpreDisp NDisp Nstay Nimmigr Nadfinal Nimmatfin" << endl;
}

//Write headers in individual output file
void outind_header(void) {
	std::string name_ind;
	name_ind = dirout + "Sc" + Int2Str(para.ScenNb) + "_Em_" + Int2Str(para.emigr) + "_Settl_" + Int2Str(para.settl) + "_Select_" + Int2Str(para.selection)  + "_Env_" + Float2Str(para.SD1_enviro) + "_Cost_" + Float2Str(para.prosp_cost) + "_ind.txt";
	outind.open(name_ind.c_str());
	outind << "sim gen PatchID NewPatchID ID Age AgeRecruit BreedPerf Disp DispProba FailInterc SucInterc FailSlope SucSlope NbProspPatch" << endl;
}


//Write population counters in population output file
void Get_OutPop(parameters para) {

	for (int q = 0; q < para.YEARS; q += para.seq_out_pop) {
		if (years == q) {
			for (int p = 0; p < para.N_patches; p++) {
				grid[p].OutPop(rep, years, &outpop);
			}

		}
	}
}


//Write individual features in individual output file
void Get_OutInd(parameters para) {

			for (int q = 0; q < para.YEARS; q += para.seq_out_ind) {
				if (years == q) {
					for (int p = 0; p < para.N_patches; p++) {
						if (grid[p].N_popsize_surv > 0) {

			for (int n = 0; n < grid[p].N_popsize_surv; n++) {
				
						grid[p].pop_vec[n].OutInd(rep, years, &outind);
					}
				}
			}
		}
	}
			}
