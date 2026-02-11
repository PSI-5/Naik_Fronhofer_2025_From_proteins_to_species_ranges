#include<stdio.h>
#include<iostream>
#include <random>
#include <iostream>
#include <cstdlib>                                                              //standard C library
#include <ctime>                                                                //access system time library
#include <fstream>                                                              //file streaming library
#include <string>                                                               //string library included
#include <sstream>                                                              //string streaming for reading numbers from
#include <vector>
#include <cmath>        //standard math library
#include <algorithm>
#include <gsl/gsl_rng.h>                                                //gsl random number generator
#include <gsl/gsl_randist.h>                                    //gsl random distributions
#include <gsl/gsl_statistics.h>                                 //gsl some statistical methods
#include <gsl/gsl_statistics_double.h>                  //gsl some double statistical methods
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_math.h>
#include <chrono>
#include <list>
using namespace std;
using namespace std::chrono;
const gsl_rng *gBaseRand;

//____________________//
class individual{
	public:
              //birth patch
	      int b_patch_no;   
              //Residence patch
	      int patch_no;
	      //dispersal probability
	      float dp;
              //Baseline energy barrier  
	      float DHm;
	      //Enzyme entropy contribution
              float DSm;
	      //Enzyme heat capacity
              float DCp;
              //rate constant for death
	      float Ed;

              //Mutation effect
              float mut_eff;
              float DDCp;
              float DDHm;

	      //Constructor
	      individual() {     
		 b_patch_no=0;     
		 dp=0;
		 patch_no=0;
		 DHm=0;
		 DSm=0;
		 DCp=0;
		 Ed=0;
		 mut_eff=0;

    }
};

//___Global Variables_____//
int Ntot; //total population
int N0; //initial population patch 
float t; //time
float b; //growth
float d; //death rate
float Beta; //competition coefficient
int tot_runs; //total replicates
float c_birth_max; //max birth rate (cb)
float c_death_max; //max death rate (cd)
float c_disp_max; //max dispersal rate
float disp; //dispersal rate
float disp_mort; //dispersal mortality
float disp_ini_sd; //dispersal initial standard deviation
float disp_mut_sd; //dispersal mutation standard deviation
float mut_rate; //mutation rate

bool disp_on; //dispersal window
bool burn_in;
float pe; //patch_extinction probability
int b_count=0; //birth count
const int RS = 2;
int world_size;
int max_pf=-1;
int min_pf=-1;
vector <individual> all_ind;
vector <int> patch_density; //patch vector for density
vector <float> birth_max; //patch vector for maximum birth rate
vector <float> disp_max; //patch vector for maximum dispersal rate
vector <float> death_max; //patch vector for maximum dispersal rate
vector <float> sum_beta; //patch vector for sum of patch beta
vector <double> patch_temp; //patch vector for temperature
vector <int> DHm_m_count; //mutations in dhm
vector <int> DSm_m_count; //mutations in dsm
vector <int> DCp_m_count; //mutations in dcp
vector <int> d_count; //mutations in dp

//Temperature dependance parameters
float Ed; //kJ
float DHm;
float DHm_ini_sd; 
float DHm_mut_sd;
float DSm;
float DCp;
float DCp_ini_sd;
float DCp_mut_sd;
float f0; //birth constant
float D0; //dispersal constant
float beta0; //death constant
float k=8.6e-5; //ev K-1
float h=4.1357e-15; 
float Tref=30+273; //K
float Tm=300;
vector<float> Trange(10000);
//____Functions________//
void specify_rng(unsigned long randSeed)
{
        gBaseRand = gsl_rng_alloc(gsl_rng_rand);

        srand(randSeed);
        unsigned long r = rand();
        gsl_rng_set(gBaseRand, r);
}

double ran()
{
        return gsl_rng_uniform(gBaseRand);
}

//Gaussian Randoms

double gauss(double sd)
{
        return gsl_ran_gaussian(gBaseRand,sd);
}

//Poisson Random

int poisson(double sd)
{
        return gsl_ran_poisson(gBaseRand,sd);
}
//Exponential Random
double expo(double sd)
{
        return gsl_ran_exponential(gBaseRand, 1/sd);
}

//Simplify Lognormal Random

double lognorm(double zeta, double sigma)
{
	double var;															//variance of resulting randoms
	double mu,s_d;														//mean and sd of lognormal distr.
																		//to be calculated by mean and
																		//sigma of resulting randoms
	var = sigma*sigma;
	s_d = sqrt(log((var/(2*zeta))+1));
	mu = log(zeta)-0.5*(s_d*s_d);
	return gsl_ran_lognormal(gBaseRand,mu,s_d);
}

//Maximum element of vector
double max(vector<double> data)
{
        return *max_element(data.begin(), data.end());
}

//Median and IQR
double median(vector<float> v)
{
  sort(v.begin(),v.end());
  int n = v.size()-1;
  if (v.size()<=2){
	  return v.at(v.size()-1);
  }
  else{  
          
	  return v.at((n+1)/2)+double((n+1)%2)*(v.at(((n+1)/2)+1)-v.at((n+1)/2))/double(2);

  }
}



//Read from parameters.in file
void readParameters(){
        ifstream parinfile("parameters.in");
        string buffer;
        istringstream is;
        getline(parinfile,buffer); getline(parinfile,buffer);
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
        is >> world_size;
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
        is >> N0;
        getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
        is >> b;
        getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
        is >> Beta;
        getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
        is >> d;
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
        is >> disp;
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
        is >> disp_mort;
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
        is >> disp_ini_sd;
	 getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
        is >> disp_mut_sd;
        getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
        is >> pe;
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
        is >> Ed;
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
        is >> DHm;
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
        is >> DHm_ini_sd;
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
        is >> DHm_mut_sd;
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
        is >> DSm;
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
        is >> DCp;
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
        is >> DCp_ini_sd;
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
        is >> DCp_mut_sd;
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
        is >> mut_rate;
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
        is >> tot_runs;
        parinfile.close();
}

//Temperature dependance
float temp_drate(float D0,float T,float Tm,float Ed){
     // float d_rate=D0*exp(-(Eb/(k*T)));
      float d_rate=D0*exp(-Ed/(k*T));
      return d_rate;
}


float DH(float T,float Tm,float DHm,float DCp){
    float dH=(DHm + DCp*(T-Tm));
    return dH;
}

float DS(float T,float Tm,float DSm,float DCp){
    float dS=(DSm+ DCp*log(T/Tm));
    return dS;
    }
float s_T(float T,float Tm,float DHm,float DSm,float DCp){
	float temp_scale;
        
	if (DS(T,Tm,DSm,DCp)<0 and DCp<0){
	  temp_scale=(k*T/h)*exp(-DH(T,Tm,DHm,DCp)/(k*T))*exp(DS(T,Tm,DSm,DCp)/(k));
	}
	else{
	   temp_scale=0; 
	}
	return temp_scale;
}
float temp_brate(int patch_no,float f0,float T,float Tm,float DHm,float DSm,float DCp){
	float temp_scale;
        if (DS(T,Tm,DSm,DCp)<0 and DCp<0){
         temp_scale=(k*T/h)*exp(-DH(T,Tm,DHm,DCp)/(k*T))*exp(DS(T,Tm,DSm,DCp)/(k));
        }
	else{
	 temp_scale=0;
	}
	float f_rate=f0*temp_scale/(1+sum_beta[patch_no]);
        return f_rate;
}

void initialise(){
	//Clear time
        t=0;    
        
	//Clear patch properties
	all_ind.clear();  //all individual vector
	patch_density.clear(); //patch-wise population vector
	birth_max.clear();   //patch-wise max birth rate vector
	disp_max.clear();   //patch-wise max dispersal rate vector
	death_max.clear();  //patch-wise max death rate vector
        sum_beta.clear();
	DHm_m_count.clear();
        DSm_m_count.clear();
        DCp_m_count.clear();
	d_count.clear();

	DSm=-(0.9-DHm)/330;
	//Calculate f0 and beta0 and D0
	f0=b*exp(DH(Tref,Tm,DHm,DCp)/(k*Tref))*exp(-DS(Tref,Tm,DSm,DCp)/(k))*(h/(k*Tref));
        D0=d*exp(Ed/(k*Tref));
        beta0=Beta*exp(DH(Tref,Tm,DHm,DCp)/(k*Tref))*exp(-DS(Tref,Tm,DSm,DCp)/(k))*(h/(k*Tref));
	
	//Clear max birth and death rates
	c_birth_max=0;
        c_death_max=0;
	c_disp_max=0;

	//Initial patch density vector
	vector <int> initial_patch_density;

	for (int y=0;y<world_size;++y){
	  if (y<int(world_size/2)+5 and y>=int(world_size/2)-5) {	
 	    initial_patch_density.push_back(N0);
	  }
	  else{
	    initial_patch_density.push_back(0);
	  }
	}
	
	//Initial patch-wise temperature
	
	for (int y=0;y<world_size;++y){
		float T=273+0.5*y;
		patch_temp.push_back(T);

	}
	
	//loop over patches
	for (int x=0;x<world_size;++x){

                int N=initial_patch_density.at(x);  
	       	//Initialise patch-wise population
		patch_density.push_back(N);
                DHm_m_count.push_back(0);
                DSm_m_count.push_back(0);
		DCp_m_count.push_back(0);
                d_count.push_back(0);
                sum_beta.push_back(0);
		
		//Calculate world population
		Ntot=Ntot+N;


		//Max birth rate value tracking 
		birth_max.push_back(0);

		//Max dispersal rate value tracking
		disp_max.push_back(0);

		//Max death rate value tracking
                death_max.push_back(0);


		//Add inidividuals to patch
                for (int i=0;i<patch_density[x];++i){
                    //Initialise single individual
		    individual newind; 
		    newind.patch_no=x;
		    newind.b_patch_no=x;
		    newind.Ed=Ed;
		    newind.dp=disp+gauss(disp_ini_sd);
		    newind.DHm=(DHm+gauss(DHm_ini_sd));
		    newind.DSm=(newind.DHm-0.9)/330;
		    newind.DCp=-abs(DCp+gauss(DCp_ini_sd));
		    sum_beta.at(x)=sum_beta.at(x)+beta0*s_T(patch_temp.at(x),Tm,newind.DHm,newind.DSm,newind.DCp);
                    
		    //Add inidividual to world		 
		    all_ind.push_back(newind);

		    //Track patch-wise max birth rate
		    if (birth_max.at(x)<temp_brate(x,f0,patch_temp.at(x),Tm,newind.DHm,newind.DSm,newind.DCp)){
			 birth_max.at(x)=temp_brate(x,f0,patch_temp.at(x),Tm,newind.DHm,newind.DSm,newind.DCp);
			 
		    }

		    //Track patch-wise max dispersal rate
		    if (disp_max.at(x)<(newind.dp)){
                         disp_max.at(x)=(newind.dp);
                    }

		    //Track patch-wise max death rate
		    if (death_max.at(x)<temp_drate(D0,patch_temp.at(x),Tm,newind.Ed)){
                         death_max.at(x)=temp_drate(D0,patch_temp.at(x),Tm,newind.Ed);
                    }

		}
	}
	//World population
	Ntot=all_ind.size();
			
}


void check_max_rate(){

	// declare variables
	vector <double> cbirth_all;
        cbirth_all.clear();
	vector <double> cdisp_all;
        cdisp_all.clear();
	vector <double> cdeath_all;
        cdeath_all.clear();
	
	// birth rate 
        for (int x = 0; x < world_size; ++x) {
		if (patch_density.at(x)!=0){
	          float T=patch_temp.at(x);		
		  float N=patch_density.at(x);
                  double help_b =birth_max.at(x);
                  cbirth_all.push_back(help_b);
		  
		}
		else{
			cbirth_all.push_back(0);
		}
        }
	c_birth_max=max(cbirth_all);

	// death rate
        for (int x = 0; x < world_size; ++x) {
                if (patch_density.at(x)!=0){
                  float T=patch_temp.at(x);
                  float N=patch_density.at(x);
                  double help_d =death_max.at(x);
                  cdeath_all.push_back(help_d);

                }
                else{
                        cdeath_all.push_back(0);
                }
        }
        c_death_max=max(cdeath_all);

	if (disp_on==1){
  	  // dispersal rate 
          for (int x = 0; x < world_size; ++x) {
		if (patch_density.at(x)!=0){
                  double help_d =disp_max.at(x);
                  cdisp_all.push_back(help_d);
		}
		else{
			cdisp_all.push_back(0);
		}
          }
	
          c_disp_max=max(cdisp_all);
	}
	else{
	       	c_disp_max=0;
	}
    
}


void birth(int x){
	b_count=b_count+1;
	//Inherited properties
        individual baby;
        
	baby.patch_no=all_ind[x].patch_no;
	baby.b_patch_no=all_ind[x].patch_no;
	baby.Ed=all_ind[x].Ed;
        if (max_pf==-1 and min_pf==-1){
                max_pf=baby.patch_no;
                min_pf=baby.patch_no;
        }
        else{
                if (baby.patch_no>max_pf and patch_density[baby.patch_no]>=9)
                {
                        max_pf=baby.patch_no;
                }
                if (baby.patch_no<min_pf and patch_density[baby.patch_no]>=9)
                {
                        min_pf=baby.patch_no;
                }
        }
	if (ran()<mut_rate){
  	   baby.dp=abs(all_ind[x].dp + gauss(disp_mut_sd));
	}
	else{
           baby.dp=all_ind[x].dp;
	}

	if (ran()<mut_rate){
          baby.DCp=(all_ind[x].DCp+gauss(DCp_mut_sd));
	  DCp_m_count[baby.patch_no]=DCp_m_count[baby.patch_no]+1;
	}
	else{
	  baby.DCp=all_ind[x].DCp;
	}

	if (ran()<mut_rate){
	   baby.DHm=(all_ind[x].DHm+gauss(DHm_mut_sd));
	   baby.DSm=((baby.DHm-0.9)/330);
	   DHm_m_count[baby.patch_no]=DHm_m_count[baby.patch_no]+1;
	   DSm_m_count[baby.patch_no]=DSm_m_count[baby.patch_no]+1;
	}
	else{
	   baby.DHm=all_ind[x].DHm;
           baby.DSm=all_ind[x].DSm;
	}

	if (temp_brate(baby.patch_no,f0,patch_temp.at(baby.patch_no),Tm,baby.DHm,baby.DSm,baby.DCp)>0){
        
//	baby.mut_eff=temp_brate(baby.patch_no,f0,patch_temp.at(baby.patch_no),Tm,baby.DHm,baby.DSm,baby.DCp)-temp_brate(all_ind[x].patch_no,f0,patch_temp.at(all_ind[x].patch_no),Tm,all_ind[x].DHm,all_ind[x].DSm,all_ind[x].DCp);
        
        
        	baby.mut_eff=f0*s_T(patch_temp.at(baby.patch_no),Tm,baby.DHm,baby.DSm,baby.DCp)- f0*s_T(patch_temp.at(all_ind[x].patch_no),Tm,all_ind[x].DHm,all_ind[x].DSm,all_ind[x].DCp);	
        	baby.DDCp=baby.DCp-all_ind[x].DCp;
                baby.DDHm=baby.DHm-all_ind[x].DHm;
		//Add individ0ual
        	all_ind.push_back(baby);
                float beta_sT=beta0*s_T(patch_temp.at(baby.patch_no),Tm,baby.DHm,baby.DSm,baby.DCp);
        	sum_beta.at(baby.patch_no)=sum_beta.at(baby.patch_no)+beta_sT;
        	
        	
        	
        	//Check max birth rate
        	if (temp_brate(baby.patch_no,f0,patch_temp.at(baby.patch_no),Tm,baby.DHm,baby.DSm,baby.DCp)>birth_max[baby.patch_no]){
                    birth_max[baby.patch_no]=temp_brate(baby.patch_no,f0,patch_temp.at(baby.patch_no),Tm,baby.DHm,baby.DSm,baby.DCp);
        	    
                }
                //Check max death rate
        	if (temp_drate(D0,patch_temp.at(baby.patch_no),Tm,baby.Ed)>death_max[baby.patch_no]){
                    death_max[baby.patch_no]=temp_drate(D0,patch_temp.at(baby.patch_no),Tm,baby.Ed);
                }
        
        
        	//Update patch and total population
        	patch_density[baby.patch_no]=patch_density[baby.patch_no]+1;     
        	Ntot=Ntot+1;
	}	
        }

void death(int x){
	//Save properties
	int patch_no=all_ind[x].patch_no;
	int N=patch_density[patch_no];
	float T=patch_temp[patch_no];
	float DHm_d=all_ind[x].DHm;
	float DSm_d=all_ind[x].DSm;
	float DCp_d=all_ind[x].DCp;
	float Ed_d=all_ind[x].Ed;
	float dispersal=all_ind[x].dp;

        //Kill
	all_ind.at(x) = all_ind.back();
	all_ind.pop_back();

	//Update population
	patch_density[patch_no]=patch_density[patch_no]-1;
	Ntot=Ntot-1;
        float beta_sT=beta0*s_T(T,Tm,DHm_d,DSm_d,DCp_d);
        sum_beta.at(patch_no)=sum_beta.at(patch_no)-beta_sT;
	//Update range limits
        if (patch_density[patch_no]<10){
           if (max_pf==patch_no){
              while (patch_density[max_pf]<10){
                    max_pf=max_pf-1;
              }
           }
           if (min_pf==patch_no){
              while (patch_density[min_pf]<10){
                     min_pf=min_pf+1;
              }
           }
        }
	//Check patch birth max
	if (birth_max.at(patch_no)==temp_brate(patch_no,f0,T,Tm,DHm_d,DSm_d,DCp_d)){
		birth_max.at(patch_no)=0;
		for (int i=0;i<Ntot;++i){
			if (all_ind.at(i).patch_no==patch_no){
				if (birth_max.at(patch_no)<temp_brate(patch_no,f0,T,Tm,all_ind[i].DHm,all_ind[i].DSm,all_ind[i].DCp)){
					birth_max.at(patch_no)=temp_brate(patch_no,f0,T,Tm,all_ind[i].DHm,all_ind[i].DSm,all_ind[i].DCp);
				}
			}
		}
	}
        
        if (disp_max.at(patch_no)==dispersal){
               // disp_max.at(patch_no)=0;
                for (int y=0;y<Ntot;++y){
                        if (all_ind.at(y).patch_no==patch_no){
				disp_max.at(patch_no)=0;
                                if (disp_max.at(patch_no)<all_ind.at(y).dp){
                                        disp_max.at(patch_no)=all_ind[y].dp;
                                }
                        }
                }
        }

	if (death_max.at(patch_no)==temp_drate(D0,T,Tm,Ed_d)){
                death_max.at(patch_no)=0;
                for (int i=0;i<Ntot;++i){
                        if (all_ind.at(i).patch_no==patch_no){
                                if (death_max.at(patch_no)<temp_drate(D0,T,Tm,all_ind[i].Ed)){
                                        death_max.at(patch_no)=temp_drate(D0,T,Tm,all_ind[i].Ed);
                                }
                        }
                }
        }
}


void patch_extinction(){
     int x=floor(ran()*world_size);
     //if (int(t)%10==0){
     if (ran()<1){
       if (patch_density[x]!=0){	
          int i=0;		   
          while (i<all_ind.size()){
		   if (ran()<0.9){
			   //Kill
                           all_ind.at(i) = all_ind.back();
                           all_ind.pop_back();
			   Ntot=Ntot-1;
			   patch_density[x]=patch_density[x]-1;

		   }
		   else{
			   i=i+1;
		   }
	  }
	}
   }

	
	
}
void dispersal(int x){

	int home_patch=all_ind.at(x).patch_no;

	//Successful dispersal
	if (ran()>disp_mort){
          		
          int target_patch;
          //Select target patch	
          if (burn_in==true){
            if (home_patch==int(world_size/2)-5){
                target_patch=home_patch+1;
            }
            else if(home_patch==int(world_size/2)+5){
                target_patch=home_patch-1;
            }
            else{
                if (ran()>0.5){
                        target_patch=home_patch+1;
                }
                else{
                        target_patch=home_patch-1;
                }
            }
          }
          else{
            if (home_patch==0){
                target_patch=home_patch+1;
            }
            else if(home_patch==(world_size-1)){
                target_patch=home_patch-1;
            }
            else{
                if (ran()>0.5){
                        target_patch=home_patch+1;
                }
                else{
                        target_patch=home_patch-1;
                }
            }
	  }	

	  //Update individual patch property 
          all_ind.at(x).patch_no=target_patch;
	  //Update patch population
	  patch_density[home_patch]=patch_density[home_patch]-1;
	  patch_density[target_patch]=patch_density[target_patch]+1;  
          float beta_sT_home=beta0*s_T(patch_temp.at(home_patch),Tm,all_ind.at(x).DHm,all_ind[x].DSm,all_ind.at(x).DCp);
          sum_beta.at(home_patch)=sum_beta.at(home_patch)-beta_sT_home;
          float beta_sT_target=beta0*s_T(patch_temp.at(target_patch),Tm,all_ind.at(x).DHm,all_ind[x].DSm,all_ind.at(x).DCp);
          sum_beta.at(target_patch)=sum_beta.at(target_patch)+beta_sT_target;

          //Track patch-wise birth rate at target
	  if (birth_max.at(target_patch)<temp_brate(target_patch,f0,patch_temp[target_patch],Tm,all_ind.at(x).DHm,all_ind[x].DSm,all_ind.at(x).DCp)){
		 birth_max.at(target_patch)=temp_brate(target_patch,f0,patch_temp[target_patch],Tm,all_ind.at(x).DHm,all_ind.at(x).DSm,all_ind.at(x).DCp);
	  }

	  //Track patch-wise birth rate at home
          if (birth_max.at(home_patch)==temp_brate(home_patch,f0,patch_temp[home_patch],Tm,all_ind.at(x).DHm,all_ind.at(x).DSm,all_ind.at(x).DCp)){
                birth_max.at(home_patch)=0;
                for (int i=0;i<Ntot;++i){
                        if (all_ind.at(i).patch_no==home_patch){
                                if (birth_max.at(home_patch)<temp_brate(home_patch,f0,patch_temp[home_patch],Tm,all_ind.at(i).DHm,all_ind.at(i).DSm,all_ind.at(i).DCp)){
                                       birth_max.at(home_patch)=temp_brate(home_patch,f0,patch_temp[home_patch],Tm,all_ind.at(i).DHm,all_ind.at(i).DSm,all_ind.at(i).DCp);
				}
                        }
                }
          }

          //Track patch-wise death rate at target
          if (death_max.at(target_patch)<temp_drate(D0,patch_temp[target_patch],Tm,all_ind.at(x).Ed)){
			death_max.at(target_patch)=temp_drate(D0,patch_temp[target_patch],Tm,all_ind.at(x).Ed);
          }

          //Track patch-wise birth rate at home
          if (death_max.at(home_patch)==temp_drate(D0,patch_temp[home_patch],Tm,all_ind.at(x).Ed)){
               death_max.at(home_patch)=0;
                for (int i=0;i<Ntot;++i){
                        if (all_ind.at(i).patch_no==home_patch){
                                if (death_max.at(home_patch)<temp_drate(D0,patch_temp[home_patch],Tm,all_ind.at(i).Ed)){
                                        death_max.at(home_patch)=temp_drate(D0,patch_temp[home_patch],Tm,all_ind.at(i).Ed);
                                }
                        }
                }
          }

	  //Track patch-wise dispersal rate at target
          if (disp_max.at(target_patch)<all_ind.at(x).dp){
                 disp_max.at(target_patch)=all_ind.at(x).dp;
          }

          //Track patch-wise dispersal rate at home
          if (disp_max.at(home_patch)==all_ind.at(x).dp){
                disp_max.at(home_patch)=0;
                for (int i=0;i<Ntot;++i){
                        if (all_ind.at(i).patch_no==home_patch){
                                if (disp_max.at(home_patch)<all_ind.at(i).dp){
                                        disp_max.at(home_patch)=all_ind.at(i).dp;
                                }
                        }
                }
          }
	}

	//Unsuccessful dispersal
	else{
            death(x);
	}
}


vector <float> average_mut_eff(int x){
  float avg=0;
   vector <float> vec;
   vector <float> stats;
   if (patch_density.at(x)==0){
           avg=0;
           vec.push_back(0);
           stats.push_back(0);
           stats.push_back(0);
           stats.push_back(0);
           return stats;

   }
   else{
    for (int i=0;i<Ntot;++i){
      if (all_ind.at(i).patch_no==x){
         int pop_dens=patch_density.at(x);
	 if (all_ind.at(i).mut_eff!=0){
           if (avg==0){
                 avg=all_ind.at(i).mut_eff;
           }
           else{
             avg=(all_ind.at(i).mut_eff+avg)/2;
           }
           vec.push_back(all_ind.at(i).mut_eff);
	 }
      }
    }
    if (vec.size()==1){
      double med=vec.at(0);
      float q25=vec.at(0);
      float q75=vec.at(0);
      stats.push_back(med);
      stats.push_back(q25);
      stats.push_back(q75);
      return stats;
    }
    else if (vec.size()==0){
      stats.push_back(0);
      stats.push_back(0);
      stats.push_back(0);
      return stats;
    }
    else{
      double med=median(vec);
      vector<float> left(vec.begin(), vec.begin() + vec.size() / 2);
      vector<float> right(vec.begin() + vec.size() / 2, vec.end());
      float q25=median(left);
      float q75=median(right);
      stats.push_back(med);
      stats.push_back(q25);
      stats.push_back(q75);
      return stats;
    }
    }

}

vector <float> average_b_patch_no(int x){
  float avg=0;
   vector <float> vec;
   vector <float> stats;
   if (patch_density.at(x)==0){
           avg=0;
           vec.push_back(0);
           stats.push_back(0);
           stats.push_back(0);
           stats.push_back(0);
           return stats;

   }
   else{
    for (int i=0;i<Ntot;++i){
      if (all_ind.at(i).patch_no==x){
         int pop_dens=patch_density.at(x);
         if (avg==0){
                 avg=all_ind.at(i).b_patch_no;
         }
         else{
           avg=(all_ind.at(i).b_patch_no+avg)/2;
         }
         vec.push_back(all_ind.at(i).b_patch_no);

      }
    }
    if (vec.size()==1){
      double med=vec.at(0);
      float q25=vec.at(0);
      float q75=vec.at(0);
      stats.push_back(med);
      stats.push_back(q25);
      stats.push_back(q75);
      return stats;
    }
    else{
      double med=median(vec);
      vector<float> left(vec.begin(), vec.begin() + vec.size() / 2);
      vector<float> right(vec.begin() + vec.size() / 2, vec.end());
      float q25=median(left);
      float q75=median(right);
      stats.push_back(med);
      stats.push_back(q25);
      stats.push_back(q75);
      return stats;
    }
    }

}

float average_disp(int x){
   float avg=0;

   if (Ntot==0){
	   avg=0;

   }   
   else{
    for (int i=0;i<Ntot;++i){
      if (all_ind.at(i).patch_no==x){	   
	 if (avg==0){
		 avg=all_ind.at(i).dp;
	 }	 
	 else{
           avg=(all_ind.at(i).dp+avg)/2;
	 }
      }
    }
   }

return avg;
}
vector <float> average_DCp(int x){
   float avg=0;
   vector <float> vec;
   vector <float> stats;
   if (patch_density.at(x)==0){
           avg=0;
           vec.push_back(0);
           stats.push_back(0);
           stats.push_back(0);
           stats.push_back(0);
           return stats;

   }
   else{
    for (int i=0;i<Ntot;++i){
      if (all_ind.at(i).patch_no==x){
         int pop_dens=patch_density.at(x);
         if (avg==0){
                 avg=all_ind.at(i).DCp;
         }
         else{
           avg=(all_ind.at(i).DCp+avg)/2;
         }
         vec.push_back(all_ind.at(i).DCp);

      }
    }
     if (vec.size()==1){
      double med=vec.at(0);
      float q25=vec.at(0);
      float q75=vec.at(0);
      stats.push_back(med);
      stats.push_back(q25);
      stats.push_back(q75);
      return stats;
    }
    else{
      double med=median(vec);
      sort(vec.begin(),vec.end());
      vector<float> left(vec.begin(), vec.begin() + vec.size() / 2);
      vector<float> right(vec.begin() + vec.size() / 2, vec.end());
      float q25=median(left);
      float q75=median(right);
      stats.push_back(med);
      stats.push_back(q25);
      stats.push_back(q75);
      return stats;
    }
   }

}



vector <float> average_DHm(int x){
   float avg=0;
   vector <float> vec;
   vector <float> stats;
   if (patch_density.at(x)==0){
           avg=0;
	   vec.push_back(0);
	   stats.push_back(0);
           stats.push_back(0);
           stats.push_back(0);
           return stats;

   }
   else{
    for (int i=0;i<Ntot;++i){
      if (all_ind.at(i).patch_no==x){
	 int pop_dens=patch_density.at(x);     
         if (avg==0){
                 avg=all_ind.at(i).DHm;
         }
         else{
           avg=(all_ind.at(i).DHm+avg)/2;
         }
	 vec.push_back(all_ind.at(i).DHm);

      }
    }
     if (vec.size()==1){
      double med=vec.at(0);
      float q25=vec.at(0);
      float q75=vec.at(0);
      stats.push_back(med);
      stats.push_back(q25);
      stats.push_back(q75);
      return stats;
    }
    else{
      sort(vec.begin(),vec.end());	    
      double med=median(vec);
      vector<float> left(vec.begin(), vec.begin() + vec.size() / 2);
      vector<float> right(vec.begin() + vec.size() / 2, vec.end());
      float q25=median(left);
      float q75=median(right);
      stats.push_back(med);
      stats.push_back(q25);
      stats.push_back(q75);
      return stats;
    }
   }
}
    
   
vector <float> average_DSm(int x){
   float avg=0;
   vector <float> vec;
   vector <float> stats;
   if (patch_density.at(x)==0){
           avg=0;
           vec.push_back(0);
           stats.push_back(0);
           stats.push_back(0);
           stats.push_back(0);
           return stats;

   }
   else{
    for (int i=0;i<Ntot;++i){
      if (all_ind.at(i).patch_no==x){
         int pop_dens=patch_density.at(x);
         if (avg==0){
                 avg=all_ind.at(i).DSm;
         }
         else{
           avg=(all_ind.at(i).DSm+avg)/2;
         }
         vec.push_back(all_ind.at(i).DSm);

      }
    }
     if (vec.size()==1){
      double med=vec.at(0);
      float q25=vec.at(0);
      float q75=vec.at(0);
      stats.push_back(med);
      stats.push_back(q25);
      stats.push_back(q75);
      return stats;
    }
    else{
      sort(vec.begin(),vec.end());	    
      double med=median(vec);
      vector<float> left(vec.begin(), vec.begin() + vec.size() / 2);
      vector<float> right(vec.begin() + vec.size() / 2, vec.end());
      float q25=median(left);
      float q75=median(right);
      stats.push_back(med);
      stats.push_back(q25);
      stats.push_back(q75);
      return stats;
    }
   }
}

vector <float> average_dp(int x){
   float avg=0;
   vector <float> vec;
   vector <float> stats;
   if (patch_density.at(x)==0){
           avg=0;
           vec.push_back(0);
           stats.push_back(0);
           stats.push_back(0);
           stats.push_back(0);
           return stats;

   }
   else{
    for (int i=0;i<Ntot;++i){
      if (all_ind.at(i).patch_no==x){
         int pop_dens=patch_density.at(x);
         if (avg==0){
                 avg=all_ind.at(i).dp;
         }
         else{
           avg=(all_ind.at(i).dp+avg)/2;
         }
         vec.push_back(all_ind.at(i).dp);

      }
    }
    if (vec.size()==1){
      double med=vec.at(0);
      float q25=vec.at(0);
      float q75=vec.at(0);
      stats.push_back(med);
      stats.push_back(q25);
      stats.push_back(q75);
      return stats;
    }
    else{
      sort(vec.begin(),vec.end());
      double med=median(vec);
      vector<float> left(vec.begin(), vec.begin() + vec.size() / 2);
      vector<float> right(vec.begin() + vec.size() / 2, vec.end());
      float q25=median(left);
      float q75=median(right);
      stats.push_back(med);
      stats.push_back(q25);
      stats.push_back(q75);
      return stats;
    }
    }
}

void bd_event(){
  //Choose event random number	
  double event_ran=ran();
  //Potential Birth event
  
  if (event_ran < (c_birth_max)/(c_birth_max+c_death_max+c_disp_max)) {
	//Select individual
	int x=floor(ran()*Ntot);
        int patch_no=all_ind[x].patch_no;	
	int N=patch_density[patch_no];
	float T=patch_temp[patch_no];
	float c_birth=temp_brate(patch_no,f0,T,Tm,all_ind.at(x).DHm,all_ind.at(x).DSm,all_ind.at(x).DCp);
	//Birth event
	if (ran()<(c_birth/c_birth_max)){
		 birth(x);
              }

        }

 //Potential Death or Dispersal eventexp(-Ea/(K*T))
  else {
	  
	if (event_ran>(c_birth_max+c_death_max)/(c_birth_max+c_death_max+c_disp_max) and disp_on==1){
	      int x=floor(ran()*Ntot);
	      int patch_no=all_ind[x].patch_no;
              d_count[patch_no]=d_count[patch_no]+1;
	      float c_dispersal=all_ind[x].dp;
	      if (ran()<c_dispersal/c_disp_max){
	         dispersal(x);
	      }
	}
	//Potential Death
        else{
	
	      //Select individual
	      int x=floor(ran()*Ntot);
	      int patch_no=all_ind[x].patch_no;
	      float T=patch_temp[patch_no];
	      float c_death=temp_drate(D0,T,Tm,all_ind.at(x).Ed);
	      if (ran()<(c_death)/c_death_max){
                 death(x); //Update rate constants
	      }
	}
  }
}


int main(int argc,char* argv[]){

  readParameters();  //read parameters for simulation runs
  specify_rng(int(stof(argv[1])+1)*time(NULL));
  
  int run_name = stof(argv[1]);
  
  for (int run_no=0;run_no<tot_runs;++run_no) {
          auto start = high_resolution_clock::now();
      //    int run_name=run_no;
          
          //Initialise world
	  initialise();    
	  
          //Calculate c_birth_max
	  check_max_rate();

          //Output file [result_(run number).out]
          stringstream outputindiv_path_stream;
          outputindiv_path_stream <<"result_"<<run_name<<"_"<<run_no<<".csv";
          string outputindiv_path = outputindiv_path_stream.str();
          ofstream outputindiv(outputindiv_path.c_str());
          // headers
	  outputindiv <<"run"<< " " <<"patch"<<" "<<"time"<< " " <<"cons"<<" "<<"dispersal_events"<<" "<<"DHm_mut_count"<<" "<<"DSm_mut_count"<<" "<<"DCp_mut_count"<<" "<<"average_dp"<<" "<<"q25_dp"<<" "<<"q75_dp"<<" "<<"average_b_patch"<<" "<<"q25_b_patch"<<" "<<"q75_b_patch"<<" "<<"average_DHm"<<" "<<"q25_DHm"<<" "<<"q75_DHm"<<" "<<"average_DSm"<<" "<<"q25_DSm"<<" "<<"q75_DSm"<<" "<<"average_DCp"<<" "<<"q25_DCp"<<" "<<"q75_DCp" <<" "<<"average_mut_eff"<<" "<<"q25_mut_eff"<<" "<<"q75_mut_eff"<<" "<<"sum_beta"<<endl;
          //average dispersal rate
	  vector <float> avg_DHm;
	  vector <float> avg_DSm;
	  vector <float> avg_DCp;
	  vector <float> avg_b_patch_no;
	  vector <float> avg_dp;
	  vector <float> avg_mut_eff;
	  list <vector<individual>> list_ind;
	  float avg_disp; 
	  int iter=0;
	  int bn_iter=1;
	  int ind_iter=0;
	  int j=0;
	  
	  //dispersal event counter
          while (int(t)<5000){
	    if (Ntot!=0 and patch_density[world_size-1]<1 and patch_density[0]<1){ 	  
              j=j+1;		    
              //Update time
 	      double lambda = (c_birth_max+c_death_max+c_disp_max)*Ntot;
	      double time_interval = expo(lambda);
	      t=t+time_interval;
	      
	      //Event
	      int time_step;
	      if (int(t)<2000){
                      burn_in=true;
		      time_step=1;
         //           disp_on=true;
              }
              else{
                      burn_in=false;
		      time_step=2;
        //            disp_on=false;
              }

	      disp_on=1;
	      // Get starting timepoint
	      bd_event();
	      // Get ending timepoint
	      //patch_extinction();
	      check_max_rate();
	      //Write output
              if (int(t)>=iter){
                for (int x=0;x<world_size;++x){	      
		  avg_DHm=average_DHm(x);
		  avg_DSm=average_DSm(x);
		  avg_DCp=average_DCp(x);
		  avg_mut_eff=average_mut_eff(x);
		  avg_b_patch_no=average_b_patch_no(x);
		  avg_dp=average_dp(x);
		  outputindiv <<j<< " " <<x<<" "<<t<< " " <<patch_density[x]<<" "<<d_count[x]<<" "<<DHm_m_count[x]<<" "<<DSm_m_count[x]<<" "<<DCp_m_count[x]<<" "<<avg_dp[0]<<" "<<avg_dp[1]<<" "<<avg_dp[2]<<" "<<avg_b_patch_no[0]<<" "<<avg_b_patch_no[1]<<" "<<avg_b_patch_no[2]<<" "<<avg_DHm[0]<<" "<<avg_DHm[1]<<" "<<avg_DHm[2]<<" "<<avg_DSm[0]<<" "<<avg_DSm[1]<<" "<<avg_DSm[2]<<" "<<avg_DCp[0]<<" "<<avg_DCp[1]<<" "<<avg_DCp[2] <<" "<<avg_mut_eff[0]<<" "<<avg_mut_eff[1]<<" "<<avg_mut_eff[2]<<" "<<sum_beta[x]<<endl; //Update output file
	          b_count=0;
                  DHm_m_count[x]=0;
                  DSm_m_count[x]=0;
                  DCp_m_count[x]=0;
		  d_count[x]=0;
		
		}
		if (iter%1==0 and int(t)>=2000){
                //list_ind.push_back(all_ind);
                     stringstream file_path_stream;
                     file_path_stream <<"individual_trait_"<<iter<<"_"<<run_name<<"_"<<run_no<<".csv";
                     string file_path = file_path_stream.str();
                     ofstream file(file_path.c_str());
                     //range limits
                     file <<"max_pf"<<" "<<max_pf<< " " <<"min_pf"<< " " <<min_pf<<endl;
                     // headers
		     file <<"individual"<< " " <<"DHm"<<" "<<"DSm"<< " " <<"DCp"<<" "<<"disp_rate"<<" "<<"patch_no"<<" "<<"b_patch_no"<<" "<<"mut_eff"<<" "<<"DDHm"<<" "<<"DDCp"<<endl;
		     
                     for (int x=0;x<all_ind.size();++x){
                          if (all_ind[x].patch_no==max_pf or all_ind[x].patch_no==max_pf-1 or all_ind[x].patch_no==max_pf-2 or all_ind[x].patch_no==min_pf or  all_ind[x].patch_no==min_pf+1 or all_ind[x].patch_no==min_pf+2 or all_ind[x].patch_no==int(world_size/2)) {
			     file <<x<< " " <<all_ind[x].DHm<<" "<<all_ind[x].DSm<<" "<<all_ind[x].DCp<<" "<<all_ind[x].dp<<" "<<all_ind[x].patch_no<<" "<<all_ind[x].b_patch_no<<" "<<all_ind[x].mut_eff<<" "<<all_ind[x].DDHm<<" "<<all_ind[x].DDCp<<endl; //Update output file
			  }
		     }
		}
                iter=iter+time_step;

	      }
	    }
	    else {
		  for (int x=0;x<world_size;++x){
                  avg_DHm=average_DHm(x);
                  avg_DSm=average_DSm(x);
                  avg_DCp=average_DCp(x);
		  avg_mut_eff=average_mut_eff(x);
		  avg_b_patch_no=average_b_patch_no(x);
                  avg_dp=average_dp(x);
		  outputindiv <<j<< " " <<x<<" "<<t<< " " <<patch_density[x]<<" "<<d_count[x]<<" "<<DHm_m_count[x]<<" "<<DSm_m_count[x]<<" "<<DCp_m_count[x]<<" "<<avg_dp[0]<<" "<<avg_dp[1]<<" "<<avg_dp[2]<<" "<<avg_b_patch_no[0]<<" "<<avg_b_patch_no[1]<<" "<<avg_b_patch_no[2]<<" "<<avg_DHm[0]<<" "<<avg_DHm[1]<<" "<<avg_DHm[2]<<" "<<avg_DSm[0]<<" "<<avg_DSm[1]<<" "<<avg_DSm[2]<<" "<<avg_DCp[0]<<" "<<avg_DCp[1]<<" "<<avg_DCp[2]<<" "<<avg_mut_eff[0]<<" "<<avg_mut_eff[1]<<" "<<avg_mut_eff[2]<<" "<<sum_beta[x]<<endl; //Update output file

		  }
		  break;
	    }
	  }
          auto stop = high_resolution_clock::now();
	  auto duration = duration_cast<microseconds>(stop - start);
          cout << "Time taken by function: "<< duration.count() << " microseconds" << endl;
	  
}
return 0;
}

		    
    
