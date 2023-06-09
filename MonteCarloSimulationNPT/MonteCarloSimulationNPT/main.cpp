//  Hard Sphere Monte Carlo Simulation.
//
//  Author: Georgios Smyridis
//

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

const string path {"/Users/georgesmyridis/Documents/Physics/Books-Notes/Graduate/Physics/Modeling_Simulations/Scripts/ModSim/MonteCarloSimulationNPT"};
const string init_filename = "xyz.dat";

// Lattice parameters
const int Nx = 3;
const int Ny = 3;
const int Nz = 3;
const double d = 1.0; // Atomic radius
const double l = 1.7 * d; // Lattice spacing for FCC lattice

// Monte Carlo Simulation
const int mc_steps = 200000;
const int output_steps = 500;

// Maximum Changes
const double delta = 1 * d; // Maximum displacement
const double delta_vol = 1; // Maximum volume change



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


void write_data(int step, vector<double> box, vector<vector<double>> positions){
    /*
     Description:
     -----------
     This function writes the box parameters and the positions of the sphere i a .dat file adding the step of the simulation  in the name.
     
     Input:
     ------
     step: Step of the Monte Carlo simulation.
     box: Vector desrcibing the box.
     positions: Current positions of the spheres.
     */

    ofstream outfile( path + "/positions/xyz" + to_string(step) + ".dat");
    if (!outfile.is_open()) {
        cerr << "Error: Unable to open output file" << endl;
        return;
    }
    
    outfile << box[0] << endl;
    outfile << box[1] << endl;
    outfile << box[2] << endl;
    outfile << box[3] << endl;
    
    for (int i {0}; i < positions.size(); i++){
        outfile << positions[i][0] << " " << positions[i][1] << " " << positions[i][2] << endl;
    }
    outfile.close();
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
    
    double dx = vector2[0] - vector1[0];
    double dy = vector2[1] - vector1[1];
    double dz = vector2[2] - vector1[2];
    double Lx = box[1], Ly = box[2], Lz = box[3];
    
    dx -= round(dx/Lx)*Lx;
    dy -= round(dy/Ly)*Ly;
    dz -= round(dz/Lz)*Lz;
    return sqrt(dx*dx + dy*dy + dz*dz);
}


int check_overlap(vector<double> box, vector<vector<double>> positions){
    /*
     Description:
     ------------
     This function checks if there is any overlap between the particles, with their current positions.
     
     Input:
     ------
     box: Vector describing the box.
     positions: Vector containing the positions of all particles.
     
     Output:
     -------
     1: If there is overlap.
     0: It there is no overlap.
     */
    for (int i {0}; i < positions.size(); i++){
        for (int j {0}; j < i; j++){
            if (distance((positions)[i],(positions[j]), box) < d){
                return 1;
            }
        }
    }
    return 0;
}


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

    double n_particles {box[0]}, Lx {box[1]}, Ly {box[2]}, Lz {box[3]};

    // Choose random particle.
    unsigned long sphere_id = round(dsfmt_genrand() * (n_particles - 1));
    //cout << "The random sphere pickes has id: " << sphere_id << endl;

    // Generate random displacement for each direction.
    double dx {(2 * dsfmt_genrand() - 1) * delta};
    double dy {(2 * dsfmt_genrand() - 1) * delta};
    double dz {(2 * dsfmt_genrand() - 1) * delta};
    
    
    // Perform displacement.
    vector<vector<double>> candidate_positions {*positions};
    candidate_positions[sphere_id][0] += dx;
    candidate_positions[sphere_id][1] += dy;
    candidate_positions[sphere_id][2] += dz;

    // Account for periodic boundaries.
    candidate_positions[sphere_id][0] -= int(candidate_positions[sphere_id][0]/Lx)*Lx;
    candidate_positions[sphere_id][1] -= int(candidate_positions[sphere_id][1]/Ly)*Ly;
    candidate_positions[sphere_id][2] -= int(candidate_positions[sphere_id][2]/Lz)*Lz;

    if (candidate_positions[sphere_id][0] < 0)
        candidate_positions[sphere_id][0] += Lx;
    if (candidate_positions[sphere_id][1] < 0)
        candidate_positions[sphere_id][1] += Ly;
    if (candidate_positions[sphere_id][2] < 0)
        candidate_positions[sphere_id][2] += Lz;

    // Check if there is overlap.
    int overlap {check_overlap(box, candidate_positions)};

    if (overlap == 1){
        // If there is overlap, reject the displacement, return 0.
        return 0;
    }
    else{
        // If there isn't overlap, accept the displacement, update positions and return 1.
        *positions = candidate_positions;
        return 1;
    }
}


int change_volume(vector<double> *box, vector<vector<double>> *positions, double betaP){
    /*
     Description:
     ------------
     This function initially chooses a random volume change. Then, it rescales the positions of the particles and checks for overlaps.
     If there are no overlaps, the box vector and the positions of the particles are updated (rescaled).
     
     Input:
     ------
     *box: Pointer that points to the vector that describes the box.
     *positions: Pointer that points to the vector of particle positions.
     
     Output:
     ------
     1: If there are no overlaps.
     0: If there are overlaps.
     */
    long unsigned n_particles {(*positions).size()};
    double dv {(2 * dsfmt_genrand() - 1) * delta_vol};
    double volume {(*box)[1] * (*box)[2] * (*box)[3]};
    double new_volume {volume + dv};
    double scale_factor {pow(new_volume / volume, 1.0/3.0)};
    
    // Rescale the positions
    vector<vector<double>> candidate_positions {*positions};
    for (int i {0}; i < n_particles; i++){
        for (int j {0}; j < 3; j++){
            candidate_positions[i][j] *= scale_factor;
        }
    }
    
    // Check if there is particle overlap
    if (check_overlap(*box, candidate_positions) == 0){
        
        // If there is no overlap, check the accepting criteria
        double acc = exp( - betaP * (new_volume - volume) + n_particles * log( new_volume / volume) );
        
        // If acceptance criterion larger than a random number from 0 to 1,
        // accept volume change, and update box and positions.
        if (acc > dsfmt_genrand()){
            // Accept volume change, update positions and box vector
            *positions = candidate_positions;
            for (int i {1}; i < (*box).size(); i++){
                (*box)[i] *= scale_factor;
            }
            return 1;
        }
        // If acceptance criterion smaller, do nothing and return 0.
        else{
            return 0;
        }
    }
    // If there is overlap, then do nothing and return 0.
    else{
        return 0;
    }
}

int main(int argc, const char * argv[]) {
    
    // Open file to write volumes
    ofstream outfile( path + "/results/results.csv");
    if (!outfile.is_open()) {
        cerr << "Error: Unable to open output file" << endl;
        return 0;
    }
    outfile << "betaP,packingfrac,deltaV,rate" << endl;
    
    // Initialisation of parameters
    dsfmt_seed( time (NULL)); // Initialise random seed.
    
    for (int i {0}; i < 250; i++){
        
        double betaP {0.1 + i * 0.4};
        
        // Initialise box vector and positions.
        vector<double> box {};
        vector<vector<double>> r {};

        // Generate FCC.
        //generate_fcc();
        read_data(init_filename, &box, &r);
        
        //Perform Monte Carlo Simulation.
        long unsigned n_particles = r.size();
        int accepted_moves {0};
        int accepted_vols {0};
        double rate {0}; //Rate of accepting volume change
        double deltaV {1};
        double packing_fraction {0};
        double packing_fraction_sum {0};

        for (int step {0}; step < mc_steps; step ++){
            
            for (int n {0}; n < n_particles; n++){
                accepted_moves += move_particle(box, &r);
            }
            accepted_vols += change_volume(&box, &r, betaP);
            
            if (step % (output_steps) == 0 and step != 0){
                
                rate = double(accepted_vols) / output_steps;
                
                double rate = double(accepted_vols)/output_steps;
                
                cout << "================ Step: " << step << " ================" << endl;
                cout << "Move acceptance rate: " << accepted_moves / (n_particles * output_steps) << endl;
                cout << "Volume change acceptance rate: " << rate << endl;
                cout << "DeltaV: " << deltaV << endl;
                //write_data(step, box, r);
                
//                if(rate < 0.4){
//                    deltaV *= 0.9;
//                } // To ensure that volume acceptance rate remains around half
//                if(rate > 0.6){
//                    deltaV *= 1.1;
//                } // If acceptance is too high or low, we change the deltaV
//                accepted_vols = 0;
                
                if(step >= 150000 && step % (output_steps) == 0){ // This is to ensure that the system has reached equillibrium before we start measuring packing fractions
                                                           // For all cases of P, the system reaches equillibrium by 150,000 steps, so this doesn't need to change.
                    packing_fraction = n_particles * M_PI * pow(d, 3) /(6 * box[0] * box[1] * box[2]);
                    packing_fraction_sum += packing_fraction;
                }
                
                accepted_moves = 0;
                
            }

        }
        outfile << betaP << "," << packing_fraction_sum / 30 << "," << deltaV << "," << rate << endl;

    }
    return 0;

}


