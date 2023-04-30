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

const string path {"/Users/georgesmyridis/Documents/Physics/Books-Notes/Graduate/Physics/Modeling_Simulations/Scripts/ModSim/LennarJonesNVT"};
const string init_filename = "xyz.dat";

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

void generate_fcc(){
    /*
     Description:
     ------------
     This function creates a configuration of hard spheres on a face-centered cubic crystal.
     
     Parameteres:
     ------------
     Nx, Ny, Nz: The number of hard spheres in each direction.
     d: The diameter of the hard spheres.
     spcacing: The lattice spacing in units of the diameter.
     
     Output:
     -------
     pc_xyz.dat file in the directory (path). The first line has the total number of particles in the box. The
     second line has the coordinates of the box corner in x-direction, and the next for the y- and z- directions
     respectively. The rest of the lines have the x-, y-, and z- coordinates of the hard sheres, and their diameter.
    */
    
    int N {4 * (Nx - 1) * (Ny - 1) * (Nz - 1)};
    double Lx {(Nx - 1) * l};
    double Ly {(Ny - 1) * l};
    double Lz {(Nz - 1) * l};
    
    ofstream outfile( path + "/positions/xyz.dat");
    if (!outfile.is_open()) {
        cerr << "Error: Unable to open output file" << endl;
        return;
    }
    
    outfile << N << endl;
    outfile << Lx << endl;
    outfile << Ly << endl;
    outfile << Lz << endl;
    
    vector<vector<double>> lattice {};
        
    for (int i {0}; i < Nx - 1; i++){
        for (int j {0}; j < Ny - 1; j++){
            for (int k {0}; k < Nz - 1; k++){

                vector<double> lattice_point1 = {i * l + d / 2, j * l + d / 2, k * l + d / 2, d};
                vector<double> lattice_point2 = {i * l + d / 2, (0.5 + j) * l + d / 2, (0.5 + k) * l + d / 2, d};
                vector<double> lattice_point3 = {(0.5 + i) * l + d / 2, (0.5 + j) * l + d / 2, k * l + d / 2, d};
                vector<double> lattice_point4 = {(0.5 + i) * l + d / 2, j * l + d / 2, (0.5 + k) * l + d / 2, d};
                
                lattice.push_back(lattice_point1);
                lattice.push_back(lattice_point2);
                lattice.push_back(lattice_point3);
                lattice.push_back(lattice_point4);
                
                outfile << i * l + d / 2 << " " << j * l + d / 2 << " " << k * l + d / 2<< " " << d << endl;
                outfile << i * l + d / 2 << " " << (0.5 + j) * l + d / 2 << " " << (0.5 + k) * l + d / 2 << " " << d << endl;
                outfile << (0.5 + i) * l + d / 2 << " " << (0.5 + j) * l + d / 2 << " " << k * l + d / 2 << " " << d << endl;
                outfile << (0.5 + i) * l + d / 2 << " " << j * l + d / 2 << " " << (0.5 + k) * l + d / 2 << " " << d << endl;
                
            }
        }
    }
}


void read_data(string file_name, vector<double> *box, vector<vector<double>> *r){
    /*
     Description:
     ------------
     This function is reads the data from a file in the path that is specifically specified with the constant string called path.
     I included the path explicitly because I could not find where XCode saved the file.

     Parameters:
     ----------
     file_name: The name of the file in the directory to read from.
     *box: Pointer, pointing to the location of the vector describing the box.
     *r: Pointer, pointing to the location of the vector containing the position of every particle.

     Outputs:
     --------
     The output is a vector with elements vectors of doubles. The first element (at position 0) contains the information about the box,
     i.e. the number of particles contained, and its dimensions in the three directions. The rest elements are vectors of which the
     first three elements are the x-,y- and z- coordinates of the hard spheres and the fourth one is the diameter of the sphere.
    */

    // Open file to read from.
    ifstream infile;
    infile.open( path + "/positions/" + file_name);
    if (not infile.is_open()){
        cerr << "Could not read the file." << endl;
        return;
    }

    int N {0};
    double Lx {0}, Ly {0}, Lz {0};
    double x {0}, y {0}, z {0}, d {0};
    
    infile >> N;
    infile >> Lx >> Ly >> Lz;
    *box = {double(N), Lx, Ly, Lz};
    
    for (int i {0}; i < N; i++){
        infile >> x >> y >> z >> d;
        vector <double> position {x,y,z,d};
//        cout << "(" << x << "," << y << "," << z << "," << d << ")" << endl;
        (*r).push_back(position);
    }
    infile.close();
    
    assert ((*box)[0] == (*r).size());
}


double distance(vector<double> vector1, vector<double> vector2, vector<double> box) {
    /*
     Description:
     ------------
     This function returns the distance between the centers of two particles, given the boundary conditions of the box.
     
     Input:
     ------
     vector1: The position of one particle.
     vector2: The position of the second particle.
     box: The  vector describing the box.
     
     Output:
     -------
     The distance.
     */
    
    double minD {0}, dist2 {0};
    for (auto i {0}; i < 3; ++i){
        minD = vector1[i] - vector2[i];
        minD -= (int)(2.0 * minD /  pow(box[i+1],2));
        dist2 += minD * minD;
    }
    
    return sqrt(dist2);
}


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


int move_particle(vector<double> box, vector<vector<double>> *positions){
    /*
     Description:
     ------------
     This function chooses a random particle, and randomly displaces it, with displacement chosen uniformly between -delta and delta.
     It then checks if this displacement creates any overlaps between the hard spheres. If there are no overlaps, it updates the positions.
     Otherwise, it rejects it.
     
     Input:
     -----
     box: Vector describing the box.
     positions: Current positions of the the particles.
     */

    double n_particles {box[0]};

    // Choose random particle.
    unsigned long sphere_id = round(dsfmt_genrand() * (n_particles - 1));
    //cout << "The random sphere pickes has id: " << sphere_id << endl;
    
    ParticleInfo info;
    info.CalculateEnergyVirial(sphere_id, &box, positions);

    // Perform random displacement.
    vector<vector<double>> candidate_positions {*positions};
    for (int i {0}; i < 3; ++i){
        candidate_positions[sphere_id][i] += (2 * dsfmt_genrand() - 1) * DELTA;
        
        // Account for periodic boundaries.
        candidate_positions[sphere_id][i] -= int(candidate_positions[sphere_id][i] / box[i+1]) * box[i+1];
        if (candidate_positions[sphere_id][i] <  0){
            candidate_positions[sphere_id][i] += box[i+1];
        }
        assert(candidate_positions[sphere_id][i] >= 0 and candidate_positions[sphere_id][i] <= box[i+1]);
    }
    
    ParticleInfo newInfo;
    newInfo.CalculateEnergyVirial(sphere_id, &box, &candidate_positions);
    
    double dE = newInfo.Energy - info.Energy;
    if (dE < 0 or dsfmt_genrand() < exp(- BETA * dE)){
        energy += dE;
        virial += newInfo.Virial - info.Virial;
        *positions = candidate_positions;
        return 1;
    }
    
    return 0;
}


void set_density(double density, vector<double> *box, vector<vector<double>> *positions){
    
    double volume = (*box)[1] * (*box)[2] * (*box)[3];
    double targetVol = (*box)[0] / density;
    double scaleFactor = pow(targetVol / volume, 1/3);
    
    for (int i {0}; i < (*box)[0]; ++i ){
        for (int j {0}; j < 3; ++j){
            (*positions)[i][j] *= scaleFactor;
        }
    }
    for (int i {1}; i < 4; ++i){
        (*box)[i] *= scaleFactor;
    }
    
}


int main(int argc, const char * argv[]) {
    
    // Set random seed.
    dsfmt_seed(time(NULL));
    
    // Initialise box vector and positions.
    vector<double> box {};
    vector<vector<double>> r {};

    // Generate FCC.
    generate_fcc();
    read_data(init_filename, &box, &r);
    set_density(DENSITY, &box, &r);
    
    //Check that RCUT is smaller than half the box's length.
    for (int d {1}; d < 4; ++d){
        assert(RCUT <= 0.5 * box[d]);
    }
    
    // Calculate the energy and virial term of the initial configuration.
    for (int n {0}; n < box[0]; ++n){
        ParticleInfo info;
        info.CalculateEnergyVirial(n, &box, &r);
        energy += info.Energy;
        virial += info.Virial;
    }
    energy *= 0.5;
    virial *= 0.5;
    
    double volume = box[1] * box[2] * box[3];
    printf("Srarting Volume: %f\n", volume);
    printf("Starting Energy: %f\n", energy);
    printf("Strating Virial: %f\n", virial);
    
    double MuIG = log(pow(LAMBDA, 3.0) * box[0] / volume) / BETA;
    

    // MONTE CARLO SIMULATION
    double *density {nullptr};
    density = &DENSITY;
    for (int i {1}; i < 61; ++i){
        (*density) =  i * 0.02;
        set_density(DENSITY, &box, &r);
        
        ofstream outfile( path + "/measurments/measurments" + to_string(i) + ".csv");
        outfile << "Step,AveragePressure,MuExcess" << endl;
        if (!outfile.is_open()) {
            cerr << "Error: Unable to open output file" << endl;
            return 0;
        }
        
        cout << "DENSITY: " << DENSITY << endl;
        int accepted {0};
        for (unsigned long step {0}; step < MC_STEPS; ++step){
            for (int n {0}; n < box[0]; ++n){
                accepted += move_particle(box, &r);
            }
            
            Measurement ms;
            ms.CalculateAveragePressure(&box, &r);
            ms.CalculateMuExcess(&box, &r);
            
            if (step > 900){
                outfile << step << "," << ms.AveragePressure << "," << ms.MuExcess + MuIG << endl;
            }
            
            if (step % 10 * OUTPUT_STEPS == 0){
                //cout << "Step: " << step << ". Move Acceptance: " << (double)accepted / (box[0] * OUTPUT_STEPS) << endl;
                accepted = 0;
            }
        }
    }
        
    // --------------- START TIMER  ----------------- //
//    auto start = chrono::high_resolution_clock::now();
//
//
//
//
//
//    auto end = chrono::high_resolution_clock::now();
//    auto duration = chrono::duration_cast<chrono::microseconds>(end - start);
//    cout << "Execution time: " << duration.count() << endl;
    
    return 0;
}
