#include <iostream> 
#include <stdio.h>
#include <cstdlib>
#include <fstream>
#include <cmath>
using namespace std;

/*
Program to study the Ising Method through a Monte Carlo approach.
*/

// We choose the values for the general variables

float T = 2.5; //Temperature
int N = 999; //Number of particles
ofstream output; 

int randomizer_positive_or_negative(void)
{     
    if(rand() % 2 == 0)
    {
        return 1;
    }
    else
    {
        return -1;
    }
}

void montecarlo(int data[999][999])
{
    int n = rand()%N, m = rand()%N;

    double p, energy_difference;

    energy_difference = 2*data[n][m]*(data[(n+1)%N][m] + data[(n-1+N)%N][m] + data[n][(m+1)%N] + data[n][(m-1+N)%N]); 

    if(pow(exp(1), -energy_difference/T) > 1)
    {
        p = 1;
    }
    else
    {
        p = pow(exp(1), -energy_difference/T);
    }

    double epsilon = rand()/double(RAND_MAX);

    if(epsilon < p)
    {
        data[n][m] = -data[n][m];
    }
}

int main()
{
    int data[999][999];

    if(1 == 1)
    {
        for(int i=0; i<N; i++)
        {
            for(int j=0; j<N; j++)
            {
                if(T < 2){data[i][j] = randomizer_positive_or_negative();}
                if(T >= 2){data[i][j] = 1;}
            }}}

    output.open("ising_data.dat");

    int pasos = 0;
    int mostrar = 0;

    for(int k=1; k<(pow(10,2))*N*N; k++)
    {
        mostrar++;
        montecarlo(data);
        if(mostrar == N*N)
        {
            for(int i=0; i<N; i++)
            {
                for(int j=0; j<N; j++)
                {
                    output << data[i][j];
                    if(j<N-1){output << ", ";}
                }
                output << endl;
            } 
            output << endl;

            pasos ++;
            mostrar=0;
        }
    }

    output.close();

    cout << pasos;
    return 0;
}