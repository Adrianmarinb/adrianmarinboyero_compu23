#include <iostream> 
#include <stdio.h>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <complex>
using namespace std; 

//General constants
const double N = 100, ciclos = 25, λ = 0.3, n_total = 120;
# define M_PI  3.14159265358979323846  /* pi */
const complex<double> i_m = complex<double>(0.0, 1.0);

double calcular_k_tilda() //Calculamos k_tilda
{
    return 2*M_PI*(ciclos/N);
}

double calcular_s_tilda(double k_tilda) //Calculamos s_tilda
{
    return 0.25 * pow(k_tilda, 2);
}

void generar_V_j(double k_tilda, double V_j[int(N)]) //Definimos V_j
{
    for(int j=0; j<=N; j++)
    {
        if (((j >= 2.0*N/10.0) && (j <= 3.0*N/10.0)) || ((j >= 6.0*N/10.0) && (j <= 7.0*N/10.0)))
        {
            V_j[j] = λ * pow(k_tilda, 2);
        }
        else
        {
           V_j[j] = 0.0;
        }

        V_j[0]=V_j[1]=0.0;
    }
    return;
}

void onda_0(complex<double> onda[int(N)][int(n_total)], double k_tilda)
{
    onda[0][0] = onda[int(N-1)][0] = 0;
    double suma = 0;
    for(int j=1; j<N; j++)
    {
        onda[j][0] = exp(i_m * k_tilda * complex<double>(j,0)) * exp(-8*(4*j - N)*(4*j - N)/(N*N)) * (1/((N/16)*sqrt(2*M_PI)));
        suma += pow(onda[j][0].real(),2) + pow(onda[j][0].imag(),2);
    }

    for(int i=1; i<N; i++)
    {
        onda[i][0] = onda[i][0] / sqrt(suma);
    }

    return;
}


void calcular_alpha(complex<double> alpha[int(N)], double s_tilda, double V_j[int(N)])
{
    alpha[int(N-1)] = 0;
    for(int j=N-1; j>0; j--)
    {
        alpha[j-1] = -pow((complex<double>(-2,0) + complex<double>(0,2)/complex<double>(s_tilda,0) - complex<double>(V_j[j],0) + alpha[j]), -1);
    }
    return;
}

void calcular_beta(double s_tilda, double V_j[int(N)], complex<double> alpha[int(N)], complex<double> onda[int(N)][int(n_total)], int n, complex<double> beta[int(N)][int(n_total)])
{
    complex<double> b[int(N)][int(n_total)];

    beta[int(N-1)][n] = complex<double>(0,0);

    for(int j=N-1; j>0; j--)
    {
        b[j][n] = complex<double>(0,4) * onda[j][n]*pow(complex<double>(s_tilda,0), -1);
        beta[j-1][n] = pow((complex<double>(-2,0) + complex<double>(0,2)/complex<double>(s_tilda,0) - complex<double>(V_j[j],0) + alpha[j]), -1) * (b[j][n] - beta[j][n]);
    }

    return;
}

void generar_chi(complex<double> chi[int(N)][int(n_total)], complex<double> alpha[int(N)], complex<double> beta[int(N)][int(n_total)], int n)
{

    chi[int(N)][n] =  chi[0][n] = 0;

    for(int j=0; j<int(N); j++)
    {
        chi[j+1][n] = alpha[j]*chi[j][n] + beta[j][n];
    }
}

void onda_n(complex<double> onda[int(N)][int(n_total)], complex<double> chi[int(N)][int(n_total)], int n)
{
    for(int j=0; j<N; j++)
    {
        onda[j][n+1]= chi[j][n] - onda[j][n];
    }

    return;
}

void imprimir_resultados(complex<double> onda[int(N)][int(n_total)], ofstream& output, double V_j[int(N)])
{
    for(int n=0; n<n_total; n++)
    {
        for(int j=0; j<=N; j++)
        {
            {
                output << j << ", " << pow(onda[j][n].real(),2) + pow(onda[j][n].imag(),2) << ", " << V_j[j] << endl;
            }
        }
        output << endl;
    }
    return;
}

void calcular_norma(complex<double> onda[int(N)][int(n_total)], int n)
{
    double sum = 0;
    for(int i=0; i<N; i++)
    {
        sum += pow(onda[i][n].real(),2) + pow(onda[i][n].imag(),2);
    }

    cout << sum << " ";
}


int main()
{
    double k_tilda = calcular_k_tilda();
    double s_tilda = calcular_s_tilda(k_tilda);
    double V_j[int(N+1)];
    generar_V_j(k_tilda, V_j);
    complex<double> onda[int(N)][int(n_total)];
    complex<double> alpha[int(N)];
    complex<double> beta[int(N)][int(n_total)];
    complex<double> chi[int(N)][int(n_total)];

    onda_0(onda, k_tilda);
    calcular_alpha(alpha, s_tilda, V_j);
    for(int n=0; n<int(n_total); n++)
    {
        calcular_beta(s_tilda, V_j, alpha, onda, n, beta);
        generar_chi(chi, alpha, beta, n); 
        onda_n(onda, chi, n);
        calcular_norma(onda, n);
    }
    ofstream output;
    output.open("schrodinger_data.dat");
    imprimir_resultados(onda, output, V_j);
    output.close();
    return 0;
}