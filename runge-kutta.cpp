#include <iostream> 
#include <stdio.h>
#include <fstream>
#include <cmath>
using namespace std;

// Valores numéricos para el reescalado, la Tierra y la Luna
const double G = 6.67 * pow(10, -11), M_T = 5.9736 * pow(10, 24);
const double M_L = 0.07349 * pow(10, 24), d_TL = 3.844 * pow(10, 8);
const double ω = 2.6617 * pow(10, -6);
const double R_T = 6.378160 * pow(10, 6), R_L = 1.7375 * pow(10, 6);
const double Δ = G*M_T / pow(d_TL, 3), μ = M_L / M_T; 
const double h = 60; 

void reescalar_vectores(double& r, double& p_r, double& p_φ, double m)
{
    r = r/d_TL;
    p_r = p_r/(m*d_TL);
    p_φ = p_φ/(m*pow(d_TL, 2)); 
    return;
}

double calcular_r(double p_r)
{
    return p_r;
}

double calcular_φ(double r, double p_φ)
{
    return p_φ/pow(r, 2);
}

double calcular_p_r(double r, double φ, double p_φ, double t)
{
    double r_tilda;
    r_tilda = sqrt(1 + pow(r, 2) - 2*r*cos(φ - ω*t));
    return pow(p_φ, 2)/pow(r, 3) - Δ*( 1/pow(r, 2) + (μ/pow(r_tilda, 3))*(r - cos(φ - ω*t)));
}

double calcular_p_φ(double r, double φ, double t)
{
    double r_tilda;
    r_tilda = sqrt(1 + pow(r, 2) - 2*r*cos(φ - ω*t));
    return -(Δ*μ*r)/pow(r_tilda, 3)*sin(φ - ω*t);
}

void repeticiones_runge_kutta(double& r, double& φ, double& p_r, double& p_φ, ofstream& output, string output_file, int pasos)
{
    output.open(output_file);
    double t = 0;
    double k1[5], k2[5], k3[5], k4[5]; // 1 es radio, 2 es φ, 3 es p_r, 4 es p_φ

    output << 0.0 << ", " << 0.0 << endl; // Posiciones iniciales (x,y) Tierra
    output << cos(ω*t) << ", " << sin(ω*t) << endl; // Posiciones iniciales (x,y) Luna
    output << r * cos(φ) << ", " << r * sin(φ) << endl; // Posiciones iniciales (x,y) Cohete
    output << endl;

    for(int i = 1; i<pasos; i++) // Casa del fantasma del código 
    {
        k1[1] = h*calcular_r(p_r);
        k1[2] = h*calcular_φ(r, p_φ);
        k1[3] = h*calcular_p_r(r, φ, p_φ, t);
        k1[4] = h*calcular_p_φ(r, φ, t);  // Aquí vive el fantasma :)

        k2[1] = h*calcular_r(p_r + k1[3]/2);
        k2[2] = h*calcular_φ(r + k1[1]/2, p_φ + k1[4]/2);
        k2[3] = h*calcular_p_r(r + k1[1]/2, φ + k1[2]/2, p_φ + k1[4]/2, t + h/2);
        k2[4] = h*calcular_p_φ(r + k1[1]/2, φ + k1[2]/2, t + h/2);

        k3[1] = h*calcular_r(p_r + k2[3]/2);
        k3[2] = h*calcular_φ(r + k2[1]/2, p_φ + k2[4]/2);
        k3[3] = h*calcular_p_r(r + k2[1]/2, φ + k2[2]/2, p_φ + k2[4]/2, t + h/2);
        k3[4] = h*calcular_p_φ(r + k2[1]/2, φ + k2[2]/2, t + h/2);

        k4[1] = h*calcular_r(p_r + k3[3]);
        k4[2] = h*calcular_φ(r + k3[1], p_φ + k2[4]);
        k4[3] = h*calcular_p_r(r + k3[1], φ + k3[2], p_φ + k3[4], t + h/2);
        k4[4] = h*calcular_p_φ(r + k3[1], φ + k3[2], t);

        r += (1.0/6.0)*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
        φ += (1.0/6.0)*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2]);
        p_r += (1.0/6.0)*(k1[3] + 2*k2[3] + 2*k3[3] + k4[3]);
        p_φ += (1.0/6.0)*(k1[4] + 2*k2[4] + 2*k3[4] + k4[4]);

        if(i%10 == 0)
        {
            output << 0.0 << ", " << 0.0 << endl; // Posicion (x,y) Tierra
            output << cos(ω*t) << ", " << sin(ω*t) << endl; // Posicion (x,y) Luna
            output << r * cos(φ) << ", " << r * sin(φ) << endl; // Posicion (x,y) Cohete
            output << endl;
        }

        t += h;
    }
    output.close();

    return;
}

int main()
{
    ofstream output; string output_file = "mision_lunar.txt";
    double altura, angulo_grados, r, φ, p_r, p_φ, masa_cohete, v_inicial;
    int pasos;
    
    altura = 0; // Altura inical sobre la tierra para la nave (m)
    angulo_grados = 51; // Angulo del cohete (º) respecto al eje x en lanzamiento
    masa_cohete = 7.5 * pow(10, 4); // Masa del cohete (kg)
    v_inicial = 11.088 * pow(10, 3); // Modulo de la velocidad inicial del cohete (m/s)

    pasos = 1.5 * pow(10, 4); //Iteraciones que queremos que de el programa
    
    φ = M_PI*angulo_grados/180; r = R_T + altura; p_r = masa_cohete*v_inicial; p_φ = 0;

    reescalar_vectores(r, p_r, p_φ, masa_cohete);

    repeticiones_runge_kutta(r, φ, p_r, p_φ, output, output_file, pasos);

    return 0;
}