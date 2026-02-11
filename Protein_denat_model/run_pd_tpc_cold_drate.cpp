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
              //Birth patch
	      int b_patch_no;
	      //Residence patch
	      int patch_no;
	      //dispersal probability
	      float dp;
              //Baseline energy barrier  
	      float Hdd;
	      //Free-energy
              float dGr;
	      //Entropy 
              float dSr;

	      //Mutation effect
	      float mut_eff;
	      float DdGr;
	      float DdSr;

	      //Neutral loci
	      //Free energy
	      float neut_dGr;
	      //Entropy
	      float neut_dSr;


	      //Constructor
	      individual() {     
		 dp=0;
		 patch_no=0;
		 b_patch_no=0;
		 Hdd=0;
		 dGr=0;
		 dSr=0;
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
float disp_mut_sd; //dispersal mutation standard deviation
float disp_ini_sd;   //dispersal initial standard deviation
bool disp_on; //dispersal window
bool burn_in;
float pe; //patch_extinction probability
float Ed; //death scaling rate
float mut_rate; //mutation rate
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
vector <int> dSr_m_count; //mutations in T0
vector <int> dGr_m_count; //mutations in sigma
vector <int> d_count; //mutations in dp

//Temperature dependance parameters
float Hdd; //kJ
float dSr;
float dSr_ini_sd;
float dSr_mut_sd;
float dGr;
float dGr_ini_sd;
float dGr_mut_sd;
int n =1;
float f0; //birth constant
float D0; 
float beta0; //intraspecific competition constant
float k=8.6e-5; //ev K-1
float h=4.1357e-15; 
float Tref=30+273; //K
float Tm=30+273;
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
  if (v.size()<=2){
          return v.at(v.size()-1);
  }
  else{
          int n = v.size()-1;
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
        is >> Hdd;
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
        is >> dSr;
        getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
        is >> dSr_ini_sd;
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
        is >> dSr_mut_sd;
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
        is >> dGr;
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
        is >> dGr_ini_sd;
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
        is >> dGr_mut_sd;
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
        is >> mut_rate;
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
        is >> tot_runs;
        parinfile.close();
}

//Temperature dependance
float s_T(float T,float Tref, float Hdd,float dGr,float dSr){
	float temp_scale;
        if ((dGr - (T - Tref) * dSr)<0){
          temp_scale=exp(-Hdd/(k*T))/(1+exp((dGr-(T-Tref)*(dSr))/(k*T)));
        }
        else{
          temp_scale=0;
        }
	return temp_scale;
}

float temp_brate(int patch_no,float f0,float T,float Tref, float Hdd,float dGr,float dSr,int n){
      float temp_scale;
      
      if ((dGr - (T - Tref) * dSr)<0){
         temp_scale=exp(-Hdd/(k*T))/(1+exp((dGr-(T-Tref)*(dSr))/(k*T)));
      }
      else{
	 temp_scale=0;
      }

      float f_rate=f0*temp_scale/(1+sum_beta.at(patch_no));
      return f_rate;
}

float temp_drate(float D0,float T,float Tm,float Hdd){
      float d_rate=D0*exp(-Ed/(k*T))+exp(-0.3*(T-273));
      return d_rate;
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
        dGr_m_count.clear();
        dSr_m_count.clear();
        d_count.clear(); 
	//Calculate f0 and beta0 and D0
	f0=b*exp(Hdd/(k*Tref))*(1+exp(dGr/(k*Tref)));
        D0=(d-exp(-0.3*(Tref-273)))*exp(Ed/(k*Tref));
        beta0=Beta*exp(Hdd/(k*Tref))*(1+exp(dGr/(k*Tref)));
	//D0=d*exp(Ed/(k*Tref));
	//Clear max birth and death rates
	c_birth_max=0;
        c_death_max=0;
	c_disp_max=0;

	//Inital patch density vector
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
                sum_beta.push_back(0);
		dGr_m_count.push_back(0);
                dSr_m_count.push_back(0);
                d_count.push_back(0);

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
		    
		    //Initialise traits
		    newind.Hdd=Hdd;
		    newind.dSr=(dSr+(gauss(dSr_ini_sd)));
		    newind.dGr=dGr+gauss(dGr_ini_sd);
		    
		    //Initialise patch/dispersal traits
		    newind.patch_no=x;
		    newind.b_patch_no=x;
		    newind.dp=disp+gauss(disp_ini_sd);

		    //Sum compeition coefficients
                    sum_beta.at(x)=sum_beta.at(x)+beta0*s_T(patch_temp.at(x),Tref,Hdd,newind.dGr,newind.dSr);  
		    
		    //Neutral loci
		    newind.neut_dSr=(dSr+(gauss(dSr_ini_sd)));
                    newind.neut_dGr=dGr+gauss(dGr_ini_sd);
		    
		    //Add inidividual to world		 
		    all_ind.push_back(newind);
		    //Track patch-wise max birth rate
		    if (birth_max.at(x)<temp_brate(x,f0,patch_temp.at(x),Tref,newind.Hdd,newind.dGr,newind.dSr,n)){
			 birth_max.at(x)=temp_brate(x,f0,patch_temp.at(x),Tref,newind.Hdd,newind.dGr,newind.dSr,n);
		    }

		    //Track patch-wise max dispersal rate
		    if (disp_max.at(x)<(newind.dp)){
                         disp_max.at(x)=(newind.dp);
                    }

		    //Track patch-wise max death rate
		    if (death_max.at(x)<temp_drate(D0,patch_temp.at(x),Tm,newind.Hdd)){
                         death_max.at(x)=temp_drate(D0,patch_temp.at(x),Tm,newind.Hdd);
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
	baby.Hdd=all_ind[x].Hdd;
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
  	   baby.dGr=all_ind[x].dGr+gauss(dGr_mut_sd);
	   dGr_m_count[baby.patch_no]=dGr_m_count[baby.patch_no]+1;
        }
	else{
	   baby.dGr=all_ind[x].dGr;
	}
        if (ran()<mut_rate){
                 baby.dSr=(all_ind[x].dSr+gauss(dSr_mut_sd));
		 dSr_m_count[baby.patch_no]=dSr_m_count[baby.patch_no]+1;
        }
        else{
                 baby.dSr=all_ind[x].dSr;
        }
        
	//Update neutral loci
	
	if (ran()<mut_rate){
           baby.neut_dGr=all_ind[x].neut_dGr+gauss(dGr_mut_sd);
        }
        else{
           baby.neut_dGr=all_ind[x].neut_dGr;
        }
        if (ran()<mut_rate){
                 baby.neut_dSr=(all_ind[x].neut_dSr+gauss(dSr_mut_sd));
        }
        else{
                 baby.neut_dSr=all_ind[x].neut_dSr;
        }
        
	if (temp_brate(baby.patch_no,f0,patch_temp.at(baby.patch_no),Tref,baby.Hdd,baby.dGr,baby.dSr,n)>0){
	   //Update mutation effect + corresponding trait changes
 	   baby.mut_eff=f0*s_T(patch_temp.at(baby.patch_no),Tref,baby.Hdd,baby.dGr,baby.dSr) - f0*s_T(patch_temp.at(all_ind[x].patch_no),Tref,all_ind[x].Hdd,all_ind[x].dGr,all_ind[x].dSr);
	   baby.DdGr=baby.dGr-all_ind[x].dGr;
	   baby.DdSr=baby.dSr-all_ind[x].dSr;
	//baby.mut_eff=temp_brate(baby.patch_no,f0,patch_temp.at(baby.patch_no),Tref,baby.Hdd,baby.dGr,baby.dSr,n) - temp_brate(all_ind[x].patch_no,f0,patch_temp.at(all_ind[x].patch_no),Tref,all_ind[x].Hdd,all_ind[x].dGr,all_ind[x].dSr,n);
	
	   float beta_sT=beta0*s_T(patch_temp.at(baby.patch_no),Tref,baby.Hdd,baby.dGr,baby.dSr);
	   sum_beta.at(baby.patch_no)=sum_beta.at(baby.patch_no)+beta_sT;
	
   	   //Add individual
	   all_ind.push_back(baby);
        
           //Check max birth rate
 	   if (temp_brate(baby.patch_no,f0,patch_temp.at(baby.patch_no),Tref,baby.Hdd,baby.dGr,baby.dSr,n)>birth_max[baby.patch_no]){
                birth_max[baby.patch_no]=temp_brate(baby.patch_no,f0,patch_temp.at(baby.patch_no),Tref,baby.Hdd,baby.dGr,baby.dSr,n);
	    
           }
           //Check max death rate
	   if (temp_drate(D0,patch_temp.at(baby.patch_no),Tm,baby.Hdd)>death_max[baby.patch_no]){
                death_max[baby.patch_no]=temp_drate(D0,patch_temp.at(baby.patch_no),Tm,baby.Hdd);
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
	float Hdd_d=all_ind[x].Hdd;
	float dSr_d=all_ind[x].dSr;
	float dGr_d=all_ind[x].dGr;
	float dispersal=all_ind[x].dp;

        //Kill
	all_ind.at(x) = all_ind.back();
	all_ind.pop_back();

	//Update population
	patch_density[patch_no]=patch_density[patch_no]-1;
	Ntot=Ntot-1;
        float beta_sT=beta0*s_T(patch_temp.at(patch_no),Tref,Hdd_d,dGr_d,dSr_d);
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
	if (birth_max.at(patch_no)==temp_brate(patch_no,f0,T,Tref,Hdd_d,dGr_d,dSr_d,n)){
		birth_max.at(patch_no)=0;
		for (int i=0;i<Ntot;++i){
			if (all_ind.at(i).patch_no==patch_no){
				if (birth_max.at(patch_no)<temp_brate(patch_no,f0,T,Tref,all_ind[i].Hdd,all_ind[i].dGr,all_ind[i].dSr,n)){
					birth_max.at(patch_no)=temp_brate(patch_no,f0,T,Tref,all_ind[i].Hdd,all_ind[i].dGr,all_ind[i].dSr,n);
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

	if (death_max.at(patch_no)==temp_drate(D0,T,Tm,Hdd_d)){
                death_max.at(patch_no)=0;
                for (int i=0;i<Ntot;++i){
                        if (all_ind.at(i).patch_no==patch_no){
                                if (death_max.at(patch_no)<temp_drate(D0,T,Tm,all_ind[i].Hdd)){
                                        death_max.at(patch_no)=temp_drate(D0,T,Tm,all_ind[i].Hdd);
                                }
                        }
                }
        }
}


void patch_extinction(){
     int x=floor(ran()*world_size);
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
          //Update intraspecific competition
	  float beta_sT_home=beta0*s_T(patch_temp.at(home_patch),Tref,all_ind.at(x).Hdd,all_ind.at(x).dGr,all_ind.at(x).dSr);
          sum_beta.at(home_patch)=sum_beta.at(home_patch)-beta_sT_home;
          float beta_sT_target=beta0*s_T(patch_temp.at(target_patch),Tref,all_ind.at(x).Hdd,all_ind.at(x).dGr,all_ind.at(x).dSr);
	  sum_beta.at(target_patch)=sum_beta.at(target_patch)+beta_sT_target;
          
	  //Track patch-wise birth rate at target
	  if (birth_max.at(target_patch)<temp_brate(target_patch,f0,patch_temp[target_patch],Tref,all_ind.at(x).Hdd,all_ind[x].dGr,all_ind.at(x).dSr,n)){
		 birth_max.at(target_patch)=temp_brate(target_patch,f0,patch_temp[target_patch],Tref,all_ind.at(x).Hdd,all_ind[x].dGr,all_ind.at(x).dSr,n);
	  }

	  //Track patch-wise birth rate at home
          if (birth_max.at(home_patch)==temp_brate(home_patch,f0,patch_temp[home_patch],Tref,all_ind.at(x).Hdd,all_ind.at(x).dGr,all_ind.at(x).dSr,n)){
               birth_max.at(home_patch)=0;
                for (int i=0;i<Ntot;++i){
                        if (all_ind.at(i).patch_no==home_patch){
                                if (birth_max.at(home_patch)<temp_brate(home_patch,f0,patch_temp[home_patch],Tref,all_ind.at(x).Hdd,all_ind.at(x).dGr,all_ind.at(x).dSr,n)){
                                        birth_max.at(home_patch)=temp_brate(target_patch,f0,patch_temp[target_patch],Tref,all_ind.at(x).Hdd,all_ind[x].dGr,all_ind.at(x).dSr,n);
                                }
                        }
                }
          }

          //Track patch-wise death rate at target
          if (death_max.at(target_patch)<temp_drate(D0,patch_temp[target_patch],Tm,all_ind.at(x).Hdd)){
			death_max.at(target_patch)=temp_drate(D0,patch_temp[target_patch],Tm,all_ind.at(x).Hdd);
          }

          //Track patch-wise birth rate at home
          if (death_max.at(home_patch)==temp_drate(D0,patch_temp[home_patch],Tm,all_ind.at(x).Hdd)){
               death_max.at(home_patch)=0;
                for (int i=0;i<Ntot;++i){
                        if (all_ind.at(i).patch_no==home_patch){
                                if (death_max.at(home_patch)<temp_drate(D0,patch_temp[home_patch],Tm,all_ind.at(i).Hdd)){
                                        death_max.at(home_patch)=temp_drate(D0,patch_temp[home_patch],Tm,all_ind.at(i).Hdd);
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

vector <float> average_neut_dGr(int x){
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
                 avg=all_ind.at(i).neut_dGr;
         }
         else{
           avg=(all_ind.at(i).neut_dGr+avg)/2;
         }
         vec.push_back(all_ind.at(i).neut_dGr);

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

vector <float> average_neut_dSr(int x){
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
                 avg=all_ind.at(i).neut_dSr;
         }
         else{
           avg=(all_ind.at(i).neut_dSr+avg)/2;
         }
         vec.push_back(all_ind.at(i).neut_dSr);

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



vector <float> average_dGr(int x){
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
                 avg=all_ind.at(i).dGr;
         }
         else{
           avg=(all_ind.at(i).dGr+avg)/2;
         }
         vec.push_back(all_ind.at(i).dGr);

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

vector <float> average_dSr(int x){
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
                 avg=all_ind.at(i).dSr;
         }
         else{
           avg=(all_ind.at(i).dSr+avg)/2;
         }
         vec.push_back(all_ind.at(i).dSr);

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
	float c_birth=temp_brate(patch_no,f0,T,Tref,all_ind.at(x).Hdd,all_ind.at(x).dGr,all_ind.at(x).dSr,n);
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
	      float c_death=temp_drate(D0,T,Tm,all_ind.at(x).Hdd);
	      if (ran()<(c_death)/c_death_max){
                 death(x); //Update rate constants
	      }
	}
  }
}


int main(int argc,char* argv[]){

  readParameters();  //read parameters for simulation runs
  specify_rng(int(stof(argv[1])+1)*time(NULL));
  
  int run_name= stof(argv[1]);
  
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
	  outputindiv <<"run"<< " " <<"patch"<<" "<<"time"<< " " <<"cons"<<" "<<"dispersal_events"<<" "<<"dGr_mut_count"<<" "<<"dSr_mut_count"<<" "<<"average_dp"<<" "<<"q25_dp"<<" "<<"q75_dp"<<" "<<"average_dGr"<<" "<<"q25_dGr"<<" "<<"q75_dGr"<<" "<<"average_dSr"<<" "<<"q25_dSr"<<" "<<"q75_dSr"<<" "<<"average_neut_dGr"<<" "<<"q25_neut_dGr"<<" "<<"q75_neut_dGr"<<" "<<"average_neut_dSr"<<" "<<"q25_neut_dSr"<<" "<<"q75_neut_dSr"<<" "<<"average_mut_eff"<<" "<<"q25_mut_eff"<<" "<<"q75_mut_eff"<<" "<<"sum_beta"<<endl;

          //average dispersal rate
	  vector <float> avg_dGr;
	  vector <float> avg_dSr;
          vector <float> avg_neut_dGr;
          vector <float> avg_neut_dSr;
	  vector <float> avg_dp;
	  vector <float> avg_mut_eff;
	  list <vector<individual>> list_ind;
	  float avg_disp; 
	  int iter=0;
	  int bn_iter=1;
	  int j=0;
	  //dispersal event counter
          while (int(t)<5000){
	    if (Ntot!=0 and patch_density[world_size-1]<1 and patch_density[0]<1){ 	  
              j=abs(j+1);		    
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
		  avg_dGr=average_dGr(x);
		  avg_dSr=average_dSr(x);
		  avg_neut_dGr=average_neut_dGr(x);
		  avg_neut_dSr=average_neut_dSr(x);
                  avg_dp=average_dp(x);
		  avg_mut_eff=average_mut_eff(x);
		  outputindiv <<j<< " " <<x<<" "<<t<< " " <<patch_density[x]<<" "<<d_count[x]<<" "<<dGr_m_count[x]<<" "<<dSr_m_count[x]<<" "<<avg_dp[0]<<" "<<avg_dp[1]<<" "<<avg_dp[2]<<" "<<avg_dGr[0]<<" "<<avg_dGr[1]<<" "<<avg_dGr[2]<<" "<<avg_dSr[0]<<" "<<avg_dSr[1]<<" "<<avg_dSr[2]<<" "<<avg_neut_dGr[0]<<" "<<avg_neut_dGr[1]<<" "<<avg_neut_dGr[2]<<" "<<avg_neut_dSr[0]<<" "<<avg_neut_dSr[1]<<" "<<avg_neut_dSr[2]<<" "<<avg_mut_eff[0]<<" "<<avg_mut_eff[1]<<" "<<avg_mut_eff[2]<<" "<<sum_beta[x]<<endl; //Update output file
	          b_count=0;
		  dGr_m_count[x]=0;
                  dSr_m_count[x]=0;
                  d_count[x]=0;
	        } 
		if (iter%1==0 and int(t)>=2000){
	        //list_ind.push_back(all_ind);
                     stringstream file_path_stream;
                     file_path_stream <<"individual_trait_"<<iter<<"_"<<run_name<<"_"<<run_no<<".csv";;
                     string file_path = file_path_stream.str();
                     ofstream file(file_path.c_str());
                     //range limits
                     file <<"max_pf"<<" "<<max_pf<< " " <<"min_pf"<< " " <<min_pf<<endl;
		     // headers
                     file <<"individual"<< " " <<"dGr"<<" "<<"dSr"<< " " <<"disp_rate"<<" "<<"patch_no"<<" "<<"b_patch_no"<<" "<<"mut_eff"<< " " <<"DdGr"<<" "<<"DdSr"<<endl;
                     for (int x=0;x<all_ind.size();++x){
                       if (all_ind[x].patch_no==max_pf or all_ind[x].patch_no==max_pf-1 or all_ind[x].patch_no==max_pf-2 or all_ind[x].patch_no==min_pf or  all_ind[x].patch_no==min_pf+1 or all_ind[x].patch_no==min_pf+2 or all_ind[x].patch_no==int(world_size/2)) {

                         file <<x<< " " <<all_ind[x].dGr<< " " <<all_ind[x].dSr<<" "<<all_ind[x].dp<<" "<<all_ind[x].patch_no<<" "<<all_ind[x].b_patch_no<<" "<<all_ind[x].mut_eff<< " " <<all_ind[x].DdGr<< " " <<all_ind[x].DdSr<<endl; //Update output file

                       }
                     }
	        }
	        iter=iter+time_step;	
	      }
	    }
	    else {
	            for (int x=0;x<world_size;++x){
                      avg_dGr=average_dGr(x);
                      avg_dSr=average_dSr(x);
		      avg_neut_dGr=average_neut_dGr(x);
                      avg_neut_dSr=average_neut_dSr(x);
		      avg_dp=average_dp(x);
		      avg_mut_eff=average_mut_eff(x);
		      outputindiv <<j<< " " <<x<<" "<<t<< " " <<patch_density[x]<<" "<<d_count[x]<<" "<<dGr_m_count[x]<<" "<<dSr_m_count[x]<<" "<<avg_dp[0]<<" "<<avg_dp[1]<<" "<<avg_dp[2]<<" "<<avg_dGr[0]<<" "<<avg_dGr[1]<<" "<<avg_dGr[2]<<" "<<avg_dSr[0]<<" "<<avg_dSr[1]<<" "<<avg_dSr[2]<<" "<<avg_neut_dGr[0]<<" "<<avg_neut_dGr[1]<<" "<<avg_neut_dGr[2]<<" "<<avg_neut_dSr[0]<<" "<<avg_neut_dSr[1]<<" "<<avg_neut_dSr[2]<<" "<<avg_mut_eff[0]<<" "<<avg_mut_eff[1]<<" "<<avg_mut_eff[2]<<" "<<sum_beta[x]<<endl; //Update output file
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

		    
    
