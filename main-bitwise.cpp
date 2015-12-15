//#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <sstream>
//#include "stdafx.h"
#include <malloc.h>
#include <string>
#include <stdio.h>
#include <time.h>
#include <limits.h>
#pragma intrinsic(_BitScanForward64)
//#define N 100
//#define nsteps 10


using namespace std;

unsigned int rand64bits();
unsigned int rand32bits();
unsigned rand256();
int bitCount(unsigned int * arr, unsigned int sizeOfArr);
void exit(string message, int exitcode);
double *vector_double(long m);
int *vector_int(long m);
double **matrix_double(long m, long n);
int **matrix_int(long m, long n);
string convertInt2Str(int number);
void *simulate (double nskip);
double configEnergy(unsigned int *old, int bit, double h, double J);
long nsteps=1000000;
int N=30;
unsigned int unit=1;


int main() {
    long nskip=nsteps/10;
    double u;
    
    
    int i,j;	//Generic integers
    clock_t start, finish; // Time variables
    ofstream myfile; //Output file handle
    
	//Initializing random number generator
    srand(time(NULL));
    
    //Throwing away the first 10000 numbers, as it is customary
    for (i = 0; i < 10000; i++) {
        u = rand();
    }
    
    //Start the clock!
    start = time(NULL);
    printf("Time length: %d \n",nsteps);
    
    simulate(nskip);
        
    //Exiting
	finish = time(NULL);
    printf("Total elapsed time: %i minutes %i seconds\n",(int)((long)(finish - start)/60),(int)((long)(finish - start) % 60));
    printf("Press any key to continue!");
    getchar();
    return 0;
    
}
    
    
//---------------------------------------------------------------------    
// Simulating

void *simulate(double nskip){
  double JperkB=10; //Interaction strength
  double *vector_JperkB;
  long i,j,k,l, iT, ih, u;
  long upspins, isteps;
  double T, h, energydiff, prob, v, sumallconfig, magntemp;
  double kB=1.0;
  double J=JperkB*kB;
    
  double *Tvect; // Temperature values
  int nT=3;
  double hmin=-2.0*J;
  double hmax=2.0*J;
  double dh=0.02*J;
  double *hvect; //Magnetic field values between hmin and hmax values with dh steps
  int nh=(int)((hmax-hmin)/dh+1);
  
  //int *initspinconfig;
  unsigned int *spinconfig; //Array where the spin configurations will be stored for a fixed T, h for every step
  //double *sites; //Spin configuration around the chosen site (the chosen site and its 2 neighbors)
  
  double **finalmagn; //Magnetization averaged on time steps for all T, h
  
  ofstream myfile; //Output file handle
  string filename; //Filename string
  
  
  //Memory allocation
  Tvect = vector_double(nT);
  Tvect[0]=JperkB/2.0;
  Tvect[1]=JperkB;
  Tvect[2]=JperkB*2.0;
  hvect = vector_double(nh);
  for ( i = 0; i < nh; i++ ){
      hvect[i]=hmin+dh*i;
  }
  //initspinconfig = vector_int(N);
  spinconfig = new unsigned int[1];
  //sites = vector_double(3);
  finalmagn = matrix_double(nT,nh);
    
  
  
  //Cycle for different values of temperature
  for ( iT = 0; iT < nT; iT++){
      T=Tvect[iT];
            
      //Cycle for different magnetic field values
      for ( ih = 0; ih < nh; ih++){
          h=hvect[ih];
          magntemp=0;
                    
          //Fixing all spins down
          memset(spinconfig, 0, sizeof(int));
          
          upspins=0; //Number of spins down
          //Fixing the initial spin config to be half up and half down
          while (!(upspins == N/2)){
            spinconfig[0] = rand32bits();
            spinconfig[0] = spinconfig[0] & ((unit<<30) - unit);
            upspins = bitCount(spinconfig, 1);
            
                        
          }
          
          //Cycle for Monte Carlo steps
          for ( isteps=0; (long)isteps<nsteps-1; isteps++){
                       
            //Choosing a site at random
            u = rand() % N; //Random number for choosing site
            
            if ( (double)rand()/RAND_MAX < std::min(1.0, exp((-1)*2.0*configEnergy(spinconfig,u,h,J)/kB/T)) ){
                 spinconfig[0] ^= (unit << (u&31)); //Flip spin
            }
            
            if ( isteps > nskip ){
                 magntemp += bitCount(spinconfig,1);
            }
                        
                        
          }
          
          finalmagn[iT][ih]=(2.0*(double)magntemp-N*(nsteps-nskip))/(nsteps-nskip);
        
        
          
      }
      
      
      
  }
  
  
    filename = "magnetization.txt";
    
    //printf("%s",filename);
    myfile.open(filename.c_str());
    //myfile <<alpha<<"\t"<<beta<<"\n";
    for ( i = 0; i < nh; i++ )
    {
        myfile <<i+1<<"\t"<<finalmagn[0][i]<<"\t"<<finalmagn[1][i]<<"\t"<<finalmagn[2][i]<<"\n";
	}
    myfile.flush();
    myfile.close();
    
    
    free(hvect);
    free(Tvect);
    free(vector_JperkB);
    
    free(spinconfig);


    for (k=0; k<nT; k++){
        free(finalmagn[k]);
    }
    free(finalmagn);


  
    
    
}


double configEnergy(unsigned int *old, int bit, double h, double J){
       int neighbor1, neighbor2;
       int bitvalue, neighvalue;
       
       neighbor1 = bit-1;
       neighbor2 = bit+1;
       if (neighbor1 < 0){ neighbor1 += N; }
       if (neighbor2 > N-1){ neighbor2 -= N; }
       
       bitvalue = ((old[0] & (unit << (bit&31)))>0);
       neighvalue = ((old[0] & (unit << (neighbor1&31)))>0) + ((old[0] & (unit << (neighbor2&31)))>0);
       
       return (double)(2*bitvalue-1)*(J*(double)(2*neighvalue-2) + h);
       
}

    
    
    
    
    
    
    
    
    
    
    
    
    
    
//--------------------------------------------------------------------------------
//General routines, useful for everything

// Borrowed from
// http://stackoverflow.com/questions/8120062/generate-random-64-bit-integer
unsigned rand256()
{
    static unsigned const limit = RAND_MAX - RAND_MAX % 256;
    unsigned result = rand();
    while ( result >= limit ) {
        result = rand();
    }
    return result % 256;
}

//See above.
unsigned int rand64bits()
{
    unsigned long long results = 0ULL;
    for ( int count = 8; count > 0; -- count ) {
        results = 256U * results + rand256();
    }
    return results;
}

unsigned int rand32bits()
{
    unsigned long long results = 0ULL;
    for ( int count = 7; count > 0; -- count ) {
        results = 256U * results + rand256();
    }
    return results;
}

//Counts the 1 bits in an array.
int bitCount(unsigned int * arr, unsigned int sizeOfArr)
{
	unsigned int count = 0;
	for (unsigned int i=0; i<sizeOfArr; i+=1)
	{
		count += __builtin_popcount(arr[i]);
	}
	return count;
}

//Writing exit messages
void exit(string message, int exitcode = 0) {
    cout << "Application exited with error: " << endl << message << endl;
    exit(exitcode);
}

//Converting integer to string
string convertInt2Str(int number) {
    stringstream ss;
    ss << number;
    return ss.str();
}

//Allocating vector of integers, size: m
int *vector_int(long m) {
    int *res = (int *) calloc(m, sizeof (int));
    return (res);
}

//Allocating matrix of doubles, size: m x n
double **matrix_double(long m, long n) {
    double **res = (double **) calloc(m, sizeof (double *));
    long i;
    for (i = 0; i < m; i++)
        res[i] = (double *) calloc(n, sizeof (double));
    return (res);
}

//Allocating matrix of integers, size: m x n
int **matrix_int(long m, long n) {
    int **res = (int **) calloc(m, sizeof (int *));
    long i;
    for (i = 0; i < m; i++)
        res[i] = (int *) calloc((long)n, sizeof (int));
    return (res);
}

//Allocating vector of doubles, size: m
double *vector_double(long m) {

    double *res = (double *) calloc(m, sizeof (double));
    return (res);
}




    
    
