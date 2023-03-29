#include <iostream> 
#include <cmath> 
#include <fstream>
#include <string>
#include <stdio.h>
using namespace std;

/*

Program to simulate the evolution of our Solar System for given magnitudes and 8 planets. It also generates an aproximated calculation of each planet's orbital period, as well as a record of the evolution of the energy of the total system over time. 

- Adrián Marín Boyero, 3ºA Physics

*/  


//Global constants
double c = 1.496*pow(10, 11);
double G = 6.67*pow(10, -11);
double M_s = 1.989*pow(10, 30);
double h = 0.0;
double t = 0.0;

//Function for the module of a distance in cartesian coords
double module(double x_1, double x_2, double y_1, double y_2)
{
    return sqrt(pow((x_1-x_2), 2) + pow((y_1-y_2), 2));
}

//Function to extract intial data into a matriz
void extract_data_into_matrix(double matriz[9][6], string& perro_salchicha)
{
   ifstream fichero;
   int i, j;
   fichero.open("initial_data.txt");
   if(fichero.is_open())
   {
    while(!fichero.eof())
    {   
        getline(fichero, perro_salchicha);
        for(i=0; i<9; i++)
        {
            for(j=0; j<6; j++)
            {fichero >> matriz[i][j];}
        }
    }
    fichero.close();
   }
   else{cout << "ERROR";}
}

//Function to reescale the intial data for an easier use of them
void rescaling_data(double initial_data[9][6], double data[9][6])
{
    for(int i=0; i<9; i++)
    {
       data[i][0] = initial_data[i][0]; //Planet's names are the same - Column 0
       data[i][1] = initial_data[i][1]/M_s; //We rescale the masses - Column 1
       data[i][2] = initial_data[i][2]/c; //We reescale x_0 - Column 2
       data[i][3] = initial_data[i][3]/c; //We reescale y_0 - Column 3
       data[i][4] = ((1/c)*pow(((G*M_s)/(pow(c, 3))), -0.5))*initial_data[i][4]; //We reescale orb_vel_x - Column 4
       data[i][5] = ((1/c)*pow(((G*M_s)/(pow(c, 3))), -0.5))*initial_data[i][5]; //We reescale orb_vel_y - Column 5
    }
    return;
}

//Function that returns the acceleration of an element for a specific cartesian coord
double newton_second_law(double data[9][6], int i, int j, string coord)
{
    double acceleration = 0;
    int column;
    if(coord == "x"){column = 2;}
    if(coord == "y"){column = 3;}
    acceleration += - data[j][1]*(data[i][column] - data[j][column])/(pow(module(data[i][2], data[j][2], data[i][3], data[j][3]), 3)); 
    return acceleration;
}

//Velvet algorithm
void velvet(double data[9][6], double acceleration[9][2], double h)
{
    double omega[9][2];
    
    for(int m=0; m<9; m++) //We empty acceleration[9][2] and omega[9][2]
    {
        for(int n=0; n<2; n++)
        {
            acceleration[m][n] = 0;
            omega[m][n] = 0;
        }}

    for(int i=0; i<9; i++) //We calculate acceleration[9][2]
    {
        for(int j=0; j<9; j++)
        {if(i!=j)
            {
                acceleration[i][0] += newton_second_law(data, i, j, "x");
                acceleration[i][1] += newton_second_law(data, i, j, "y");
            }}}

    for(int i_o=0; i_o<9; i_o++)  //We evaluate omega[9][2] using omega = coord + (h/2)*accel
    {
        for(int j_o=0; j_o<2; j_o++)
        {
            omega[i_o][j_o] = data[i_o][(j_o) + 4] + (h/2)*acceleration[i_o][j_o];
        }}
    
    for(int i=0; i<9; i++) //we calculate new coords by using coord_new = coord + h*omega
    {
        for(int j=0; j<2; j++)
        {
            data[i][j+2] = data[i][j+2] + h*omega[i][j];
            acceleration[i][j] = 0;
        }}

    for(int i=0; i<9; i++)  //We recalculate acceleration by making use of the new coords
    {
        for(int j=0; j<9; j++)
        {
            if(i!=j)
            {
                acceleration[i][0] += newton_second_law(data, i, j, "x");
                acceleration[i][1] += newton_second_law(data, i, j, "y");
            }}}
    
    for(int p=0; p<9; p++) //We recalculate new velocities by using v_new = omega + acceleration_new
    {
        for(int r=0; r<2; r++)
        {
            data[p][r+4] = omega[p][r] + (h/2)*acceleration[p][r];
        }}

}

//Function to calculate each planet's energy
void calculate_energy(double data[9][6], double energy[9], int step)
{
    //We calculate Mechanical energy (potential + kinetic) for each planet
    for(int pot=0; pot<9; pot++)
    {  
        for(int pot_2=0; pot_2<9; pot_2++)
        {
            if(pot!=pot_2)
            {
                energy[pot+1] += - data[pot_2][1]*data[pot][1]*(1/module(data[pot][2], data[pot_2][2], data[pot][3], data[pot_2][3]));
            }}
        energy[pot+1] += 0.5*data[pot][1]*(pow(data[pot][4], 2) + pow(data[pot][5], 2));
    }
}

//Record the orbital period of specific planets, for a given aproximation
void fill_periods(double data[9][6], double periods[8], int step, double h, double y_limit, int initial, int final)
{
    for(int i = initial; i<= final; i++)
    {
        if(periods[i] == 0.0)
        {
            //If the planet is "close" in y_coord to the x axis, and its x_coord is positive we can aproximate periods
            if (data[i+1][3] > 0 && data[i+1][3] < y_limit && data[i+1][2] > 0) 
            {
                periods[i] = 58.1*(h*step);
            }}}
    return;
}

int main()
{
    double initial_data[9][6], data[9][6]; double acceleration[9][2];
    string perro_salchicha, planets[9] = {"Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"};
    ofstream output, energy_file, periods_file; 
    double periods[8] = {}, energy[9] = {};
    int iterations = 15000;

    //We first get values from our file with initial data
    extract_data_into_matrix(initial_data, perro_salchicha);
    rescaling_data(initial_data, data);

    output.open("planets_data_adrian.dat");
    energy_file.open("energy_values.txt");
    
    for(int step=0; step < iterations; step++) //First iterations to have higher precision for closer planets using a higher h
    {
        double energy[9] = {};
        h = 0.01;
        velvet(data, acceleration, h);

        //We fill the periods of closer planers
        double y_limit = 0.01;
        fill_periods(data, periods, step, h, y_limit, 0, 3);

        //We output the energy of each planet
        calculate_energy(data, energy, step);
        if(step % 30 == 0)
        {
            energy_file << 58.1*(h*step) << "   ";
            for(int e=0; e<10; e++)
            {
                energy_file << energy[e] << "   ";
            }
            energy_file << endl;
        }

        //We output the new coords for each planet into the output file, with the correct format
        if(step % 30 == 0)
        {   
            for(int s=0; s<9; s++)
            {
                output << data[s][2] << ", " << data[s][3] << endl; 
            }
            output << endl;
        }}
    
    t = 0;
    rescaling_data(initial_data, data);
    for(int step2=0; step2 < 3*iterations; step2++) //Second iterations with lower h so that further planets complete their orbits
    {   
        h = 0.3;
        double y_limit = 0.04;
        velvet(data, acceleration, h);
        fill_periods(data, periods, step2, h, y_limit, 4, 7); //We fill the periods of closer planers
    }
    output.close();
    energy_file.close();

    //We print the period of the planets in a document
    periods_file.open("periods.txt");
    periods_file << "----------------------------------------------" << endl;
    periods_file << "Diferent orbital periods in our Solar System: " << endl;
    periods_file << "----------------------------------------------" << endl << endl;
    for(int per=0; per<8; per++)
    {
        periods_file << planets[per+1] << ": "<< periods[per] << " days or " << periods[per]/365.24 << " Earth years." << endl;
    }
    periods_file.close();

    return 0;
}