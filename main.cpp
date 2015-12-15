//#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <sstream>
//#define N 100
//#define nsteps 10


using namespace std;

void exit(string message, int exitcode);
double *vector_double(long m);
int *vector_int(long m);
double **matrix_double(long m, long n);
int **matrix_int(long m, long n);
string convertInt2Str(int number);
void *simulate (double nskip);
long nsteps=10000;
int N=30;


int main() {
    int nprocs=1;
    long nskip=nsteps/10;
    double u;
    
    
    int i,j;	//Generic integers
    clock_t start, finish; // Time variables
    ofstream myfile; //Output file handle
    
    //if ( argc != 2 ) // 1 argument is the default (Tlength)
//    {
//        exit("Not appropriate amount of arguments", 1);
//    }
//    else 
//    {
//        Tlength = strtod(argv[1],NULL); // Setting the time (length) of simulation
//	}
	
	//Initializing random number generator
    srand(time(NULL));
    
    //Throwing away the first 10000 numbers, as it is customary
    for (i = 0; i < 10000; i++) {
        u = rand();
    }
    
    //Start the clock!
    start = time(NULL);
    printf("Number of threads: %i \n",nprocs);
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
  long downspins, isteps;
  double T, h, energydiff, prob, v, sumallconfig, magntemp;
  double energy, sumenergy;
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
  double **spinconfig; //Array where the spin configurations will be stored for a fixed T, h for every step
  double *sites; //Spin configuration around the chosen site (the chosen site and its 2 neighbors)
  
  double **finalmagn; //Magnetization averaged on time steps for all T, h
  double **finalenergy;
  
  ofstream myfile; //Output file handle
  string filename; //Filename string
  ofstream myfile2; //Output file handle
  string filename2; //Filename string
  
  
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
  spinconfig = matrix_double(N,nsteps);
  sites = vector_double(3);
  finalmagn = matrix_double(nT,nh);
  finalenergy = matrix_double(nT,nh);
    
  
  
  //Cycle for different values of temperature
  for ( iT = 0; iT < nT; iT++){
      T=Tvect[iT];
            
      //Cycle for different magnetic field values
      for ( ih = 0; ih < nh; ih++){
          h=hvect[ih];
                    
          //Fixing all spins up
          for (i=0; i<N; i++){
              spinconfig[i][0]=1.0;
          }
          
          downspins=0; //Number of spins down
          //Fixing the initial spin config to be half up and half down
          while (downspins<N/2){
            //If the chosen spin is up, change it to down, until half of
            //them will be down
            u = (int)(rand() % N); //Random number for choosing site
            if (spinconfig[u][0]==1.0){
                spinconfig[u][0]=-1.0;
                downspins=downspins+1;
            }
                        
          }
          
          //Cycle for Monte Carlo steps
          for ( isteps=0; (long)isteps<nsteps-1; isteps++){
            //New spin configuration is defined as the previous
            for(i=0; i<N; i++){
                spinconfig[i][isteps+1]=spinconfig[i][isteps];
            }
            
            //Choosing a site at random, and putting its and its neighbors
            //configuration into variable sites
            u = (int)(rand() % N); //Random number for choosing site
            if (u==0 || u==N-1){ //Dealing with boundaries (periodic boundary condition)
                if (u==0){
                    sites[0]=spinconfig[N-1][isteps];
                    sites[1]=spinconfig[0][isteps];
                    sites[2]=spinconfig[1][isteps];
                }
                else{
                    sites[2]=spinconfig[0][isteps];
                    sites[0]=spinconfig[N-2][isteps];
                    sites[1]=spinconfig[N-1][isteps];
                }
            }
            else{
                sites[0]=spinconfig[u-1][isteps];
                sites[1]=spinconfig[u][isteps];
                sites[2]=spinconfig[u+1][isteps];
            }
            
            //Energy change if we were to flip the chosen spin
            energydiff=2.0*J*(sites[0]*sites[1]+sites[1]*sites[2])+2.0*h*sites[1];
            
            //Metropolis algorithm
            if ( 1 > exp(-energydiff/kB/T) ){
                 prob=exp(-energydiff/kB/T);
            }
            else{
                 prob=1;
            }
            
            //If random number is smaller than prob, flip spin
            v = (double)rand()/((double)RAND_MAX + 1);
            if (v < prob){
                spinconfig[u][isteps+1]=-sites[1]; //Flip spin (accept)
            } //Else do nothing (reject)
            
            
                        
          }
          
          
        magntemp=0;
        energy=0;
        //Calculating magnetization averaged on time
        for (j=(long)nskip; j<nsteps; j++){
            sumallconfig=0;
            sumenergy=0;
            for (i=0; i<N-1; i++){
                sumallconfig=sumallconfig+spinconfig[i][j];
                sumenergy=sumenergy-J*(spinconfig[i][j]*spinconfig[i+1][j])-h*spinconfig[i][j];
            }
            sumenergy=sumenergy-J*(spinconfig[N-1][j]*spinconfig[0][j])-h*spinconfig[N-1][j];
            energy=energy+sumenergy;
            magntemp=magntemp+sumallconfig;
        }
        finalmagn[iT][ih]=(double)magntemp/(nsteps-nskip);
        finalenergy[iT][ih]=(double)energy/(nsteps-nskip);
        
        
          
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
    
    filename2 = "energy.txt";
    
    //printf("%s",filename);
    myfile2.open(filename2.c_str());
    //myfile <<alpha<<"\t"<<beta<<"\n";
    for ( i = 0; i < nT; i++ )
    {
        myfile2 <<i+1<<"\t"<<finalenergy[i][0]<<"\t"<<finalenergy[i][1]<<"\t"<<finalenergy[i][2]<<"\n";
	}
    myfile2.flush();
    myfile2.close();
    
    
    free(hvect);
    free(Tvect);
    free(vector_JperkB);
    free(sites);
    
    
    for (l=0; l<N; l++){
        free(spinconfig[l]);
    }
    free(spinconfig);


    for (k=0; k<nT; k++){
        free(finalmagn[k]);
    }
    free(finalmagn);
    
    for (k=0; k<nT; k++){
        free(finalenergy[k]);
    }
    free(finalenergy);


  
    
    
}




    
    
    
    
    
    
    
    
    
    
    
    
    
    
//--------------------------------------------------------------------------------
//General routines, useful for everything

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




    
    
