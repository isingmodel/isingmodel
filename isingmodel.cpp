#include <stdio.h>
#include "twist.h"
#include <math.h>
#include <stdlib.h>

void connect(int L, int* spin, int** nn, int N);
void init(int N, int* spin);
void monte(int* spin, int** nn, int step, int N, int z, double* expo);
void Elecmag(int* spin, int** nn, int N, int z);




int main(int argc, char ** argv) {
    
    //declaration
    printf("1\n");
    init_rnd(gus());
    double E = 0;
    double M = 0;
    int L = 5;//atoi(argv[1]);              //atoi : a to integer
    double T = 1;//atof(argv[2]);           // atof = ascii to float
    int z = 4;                          // the number of netighbor, z=4 for 2d ising model.
    int ensemble = 1;                   //ensemble size
    int N = L*L;                        // total number of nodes
    int step = 1000;
    int* spin = (int*)calloc(N, sizeof(int));
    int** nn = (int**)calloc(N, sizeof(int*));
    double* expo = (double*)calloc(2, sizeof(double));
    expo[0] = exp(-4./T);
    expo[1] = exp(-8./T);
    printf("2\n");
    //initializing//
    for(int i=0; i<N; i++)
    {
        nn[i] = (int*)calloc(z, sizeof(int));
    }
    printf("3\n");
    connect(L, spin, nn, N);
    printf("4\n");
    init(N, spin);
    ///////////////
    printf("5\n");
    monte(spin, nn, step, N, z, expo);

    printf("6\n");
    
    //Elecmag(spin, nn, N, z);
    
    
    
    
    
    
    //printf("%f\n", E);
   // printf("%f\n", M);
    return 0;
}




//////connecting every single nodes with their neighbors.//////
void connect(int L, int* spin, int** nn, int N)
{
    for(int i=0; i<N; i++)
    {
        nn[i][0] = i+1;
        nn[i][1] = i-1;
        nn[i][2] = i+L;
        nn[i][3] = i-L;
        if(i%L == 0)
        {
            nn[i][1] = i+L-1;
        }
        if(i%L == L-1)
        {
            nn[i][0] = i-L+1;
        }
        if(i<L)
        {
            nn[i][3] = N-L+i;
        }
        if(i>N-L-1)
        {
            nn[i][2] = i+L-N;
        }
    }
}




///////////////initializing/////////////////////////////
void init(int N, int* spin)
{
    for(int i = 0; i<N; i++)
    {
        if(drnd()<0.5)
        {
            spin[i] = 1;
        }
        else
        {
            spin[i] = -1;
        }
    }
}





//////////////montecarlo step////////////////////////
void monte(int* spin, int** nn, int step, int N, int z, double* expo)
{ 
    int x;
    double E;
    int EE;
    for (int i = 0; i < N; i++) // repeating 'step' times
    {
        E = 0.; //initializing Energy
        x = (int)(drnd() * N); // choose a random node
        for(int j = 0; j<z; j++)
        {
            E += spin[x] * spin[nn[x][j]];
        }
        
        EE = int(E/4.+0.5)-1;
        if(E <= 0)
        {
            spin[x] = -1 * spin[x];
        }
        else if(drnd() < expo[EE])
        {
            spin[x] = -1 * spin[x];
        }
    }
}




//////////////calculating E and M/////////////////////
void Elecmag(int* spin, int** nn, int N, int z)
{
    float M = 0;   //initializing
    float E = 0;   //initializing
    for(int i = 0; i<N; i++)
    {
        M += spin[i] * 1.0;
        for(int j = 0; j<z; j++)
        {
            E += - spin[i]*spin[nn[i][j]];
        }
    }
    E = E/2;
}