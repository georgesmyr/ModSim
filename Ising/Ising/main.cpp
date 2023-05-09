#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <random>
#include <ctime> // for time()
#include <string>
#include <unordered_map>
#include "mt19937.h"

#include <chrono>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;

const string path = "/Users/georgesmyridis/Desktop/Trading/ModSim/Ising/measurements/";

class IsingModel{

private:
    
    // Lattice parameters
    int L = 30;
    int spins[30][30] = {0};
    int J = 1;
    double N = L * L;
    double TEMP = 2;
    double BETA = 1 / TEMP;
    
    // Statistical quantities
    int magnetisation = 0;
    int energy = 0;
    
    // Monte Carlo simulation parameters
    const unsigned long STEPS = 200000;
    
    // Output file
    const string NAME = "measurement0";
    

public:
    
    
    double cached_exp(double x) {
        // Define the four possible values of x.
        const double values[] = {J * 4.0, J * 8.0, J * 12.0, J * 16.0};

        // Initialize a static unordered map to cache the exponential values.
        static unordered_map<double, double> cache = {
            {values[0], exp(- BETA * values[0])},
            {values[1], exp(- BETA * values[1])},
            {values[2], exp(- BETA * values[2])},
            {values[3], exp(- BETA * values[3])}
        };

        // Return the cached value if it exists, otherwise calculate and cache it.
        if (cache.find(x) != cache.end()) {
            return cache[x];
        } else {
            double result = exp(x);
            cache[x] = result;
            return result;
        }
    }
    
    
    void initConfig(){
        /*
         Sets the state of the system corresponding to zero temperature.
         */
        for (int i = 0; i < L; ++i){
            for (int j = 0; j < L; ++j){
                spins[i][j] = 1;
            }
        }
        
        // Calculation of energy
        for (int i = 0; i < L; ++i) {
            int ip1 = (i + 1) % L;
            int im1 = (i - 1 + L) % L;
            for (int j = 0; j < L; ++j) {
                int jp1 = (j + 1) % L;
                int jm1 = (j - 1 + L) % L;
                energy += spins[i][j] * (spins[ip1][j] + spins[im1][j] + spins[i][jp1] + spins[i][jm1]
                                      + spins[im1][jm1] + spins[im1][jp1] + spins[ip1][jm1] + spins[ip1][jp1]);
            }
        }
        energy *= -J;
        // Calculation of magnetisation
        magnetisation = 0;
        for (int i = 0; i < L; ++i){
            for (int j = 0; j < L; ++j){
                magnetisation += spins[i][j];
            }
        }
    }
    
    
    int energyChange(int i, int j){
        /*
         Returns the energy change if we flip the spin (i, j).
         */
        int sum = 0;
        for (int di = -1; di <= 1; ++di){
            for (int dj = -1; dj <= 1; ++dj){
                if (di == 0 && dj == 0){
                    continue;
                }
                sum += spins[(i + di + L) % L][(j + dj + L) % L];
            }
        }
        return 2 * spins[i][j] * sum;
    }
    
    
    void run(){
        
        // Initialise first configuration
        initConfig();
        
        // Open output csv file
        ofstream outfile(path + NAME + ".dat");
        if (!outfile.is_open()) {
            cerr << "Error: Unable to open output file" << endl;
        }
        outfile << "sweep,energy,magnetisation" << endl;
        
    
        // Repeat for #STEPS steps
        int i = 0, j = 0, energy_change = 0;
        double r = 0;
        for (auto step = 0; step < STEPS; ++step){
            
            // Choose random spin
            i = int(dsfmt_genrand() * L);
            j = int(dsfmt_genrand() * L);
                            
            // Calculate the change in energe
            energy_change = energyChange(i, j);
            if (energy_change <= 0){
                // If energy change is non-positive
                // change spin, update energy and magnetisation
                spins[i][j] *= -1;
                energy += energy_change;
                magnetisation += 2 * spins[i][j];
            }else{
                r = dsfmt_genrand();
                if (r < cached_exp(energy_change)){
                    spins[i][j] *= -1;
                    energy += energy_change;
                    magnetisation += 2 * spins[i][j];
                }
            }
            outfile << step / N << "," << energy / N << "," << magnetisation / N << endl;
        }
        outfile.close();
    }
    
    
    double autocorrilationTime(){
        
        return 0;
    }
    
    
};

int main(int argc, const char * argv[]) {
    
//    dsfmt_seed(time(NULL)); // Initialise random seed
    dsfmt_seed(10);
    
    IsingModel is;
    is.run();


    return 0;
}
