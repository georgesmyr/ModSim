#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <random>
#include <ctime> // for time()
#include <string>
#include "mt19937.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;

// Lattice parameters
const int Nx = 3;
const int Ny = 3;
const int Nz = 3;
const double d = 1.0; // Atomic radius
const double l = 3 * d; // Lattice spacing for FCC lattice

// Monte Carlo Simulation
const int MC_STEPS = 1000;
const int OUTPUT_STEPS = 100;
const double DELTA = 1 * d; // Maximum displacement

// Lennard Jones Dynamics Parameters
const double RCUT = 2.5 * d;
const double BETA = 1;
double DENSITY = 0.5;
const double ECUT {4 * (pow(d / RCUT, 12) - pow(d / RCUT, 6))};
const int NTEST = 20;
const double LAMBDA = 1;

double energy = 0;
double virial = 0;




class ParticleInfo{
    /*
     Description:
     -----------
     This class contains the contribution of the particle in the enegy and virial term of the system.
     */
public:
    double Energy {0};
    double Virial {0};
    
    void CalculateEnergyVirial(unsigned long pid, vector<double> *box, vector<vector<double>> *positions){
        /*
         Description:
         -----------
         This function calculates the contribution of the particle pid to the system's energy and virial term.
         
         Input:
         -----
         pid: The index of the particle of which we want to calculate the contributins to the energy and mu.
         *box: Pointer pointing to the box vector.
         *positions: Pointer pointing to the current positions.
         
         */
        for (auto i {0}; i < (*box)[0]; ++i){
            if (i == pid){
                continue;
            }
            double dist = distance((*positions)[pid], (*positions)[i], (*box));
            if (dist <= RCUT){
                double temp = pow(dist, -6);
                Energy += 4 * temp * (temp - 1) - ECUT;
                Virial += 24 * temp * (2 * temp - 1);
            }
        }
    }
    
};


class Measurement{
    /*
     This class calculates and holds the average pressure and MueExcess when it's called.
     */
public:
    double AveragePressure {0};
    double MuExcess {0};
    
    void CalculateAveragePressure(vector<double> *box, vector<vector<double>> *positions){
        
        double totalVirial {0};
        double dist {0};
        for (int i {0}; i < (*box)[0]; ++i){
            for (int j {0}; j < i; ++j){
                dist = distance((*positions)[i], (*positions)[j], *box);
                totalVirial += 24 * (2 * pow(dist,-12) - pow(dist, -6));
            }
        }
        double volume = (*box)[1] * (*box)[2] * (*box)[3];
        AveragePressure = DENSITY / BETA + totalVirial / (3 * volume) ;
        
    }
    
    void CalculateMuExcess(vector<double> *box, vector<vector<double>> *positions){
        /*
         Description:
         ------------
         This function calculates the MuExcess. It inserts a particle in a random position in the box, calculates
         the change in the systems energy and exponentiates it. For statistical significance, it does this NTEST times
         and averages.
         
         Input:
         ------
         *box: Pointer pointing to the vector describing the box.
         *positions: Pointer pointing to the current positions of the particles.
         */
        
        double totalMuExcess {0};
        for (int i {0}; i < NTEST; ++i){
            // Insert test particle at random position
            vector<double> testParticlePosition {dsfmt_genrand() * (*box)[1], dsfmt_genrand() * (*box)[2], dsfmt_genrand() * (*box)[3]};
            (*positions).push_back(testParticlePosition);
            ParticleInfo testInfo;
            testInfo.CalculateEnergyVirial((*box)[0], box, positions);
            totalMuExcess += exp(-BETA * testInfo.Energy);
            (*positions).pop_back();
        }
        MuExcess = - log(totalMuExcess / NTEST) / BETA;
    }
};
