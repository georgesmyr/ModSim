
//  PairCorrelationNVT
//
//  Created by Georgios Smyridis
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

const string path {"/Users/georgesmyridis/Documents/Physics/Books-Notes/Graduate/Physics/Modeling_Simulations/Scripts/ModSim/PairCorrelationNVT"};
const string init_filename = "xyz.dat";

// Lattice parameters
const int Nx = 3;
const int Ny = 3;
const int Nz = 3;
const double d = 1.0; // Atomic radius
const double l = 2 * d; // Lattice spacing for FCC lattice

// Monte Carlo Simulation
const int mc_steps = 1000;
const int output_steps = 100;

// Maximum Changes
const double delta = 0.2 * d; // Maximum displacement


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
    
    cout << "Generating FCC..." << endl;
    
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
    cout << "FCC generated." << endl;
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
    
    cout << "Reading Data..." << endl;

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
    
    cout << "Data read." << endl;
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


void save_RAD(vector<double> histogram, double bin_size, int num){
    /*
     Description:
     -----------
     This function saves the radial distribution function histogram in a csv file.
     
     Input:
     ------
     step: Step of the Monte Carlo simulation.
     box: Vector desrcibing the box.
     positions: Current positions of the spheres.
     */
    
    cout << "Saving RAD..." << endl;

    ofstream outfile( path + "/results/rad" + to_string(num) + ".csv");
    if (!outfile.is_open()) {
        cerr << "Error: Unable to open output file" << endl;
        return;
    }
    
    outfile << "density,bin" << endl;
    for (int i {0}; i < histogram.size(); i++){
        outfile << histogram.at(i) << "," << (i+1) * bin_size << endl;
    }
    cout << "RAD saved." << endl;
}


void calculate_distances(vector<double> box, vector<vector<double>> positions, int num){

    cout << "Calculating distances..." << endl;

    ofstream outfile( path + "/distances/dist" + to_string(num) + ".csv");
    if (!outfile.is_open()) {
        cerr << "Error: Unable to open output file" << endl;
        return;
    }
    outfile << "distance" << endl;

    long unsigned n_particles = box[0];
    for (int i {0}; i < n_particles; i++){
        for (int j {0}; j < n_particles; j++){
            if (i != j){
                outfile << distance(positions[i], positions[j], box) << endl;
            }
        }
    }
    
}

void calculate_RAD(vector<double> box, vector<vector<double>> positions, int num){
    /*
     Description:
     -----------
     This function calculates the radial distribution of a particle configuration inside a box with positions given as an input.
     It places the values in bins of a histogram and then saves it as a csv file.
     
     Input:
     ------
     box: Vector describing the box
     positions: The current positions of all the particles.
     */
    
    cout << "Calculating RAD..." << endl;
    
    // Parameter Initialisation
    double max_distance {box[1] * sqrt(3)}; // The maximum distance between particles is the diagonal of the box
    int N {int(box[0])}; // Total number of particles in the box
    double density {N / (box[1] * box[2] * box[3])};
    int num_bins {50};
    double bin_size {max_distance / num_bins};
    
    // Histogram and normalising factors initialisation
    vector<double> histogram(num_bins, 0);
    vector<double> hist_factors(num_bins, 1.0);
    
    // We skip the bins that represent distances exclusively inside the hard ball, i.e. smaller than the radius
    int skip {0};
    skip = floor(d / (2 * bin_size)) - 1;
    // Set the normalising parameters for the rest entries
    for (int i {0}; i > skip; i++){
        hist_factors[i] = (M_PI * density / 6) * (pow(1 + i * bin_size, 3) - pow(1 + (i-1) * bin_size,3));
    }

    // Calculate all the distances and place them in the correct bin
    for (int i {0}; i < box[0]; i++){
        for (int j {0}; j < i; j++){

            double dist {0};
            dist = distance(positions[i], positions[j], box);
            
            int bin {0};
            bin = floor(dist / bin_size) - 1;
            histogram.at(bin) += 1;

        }
    }
    
    // Normalise
    for (int i {0}; i < num_bins; i++){
        histogram.at(i) /=  N * hist_factors.at(i);
        cout << "Bin: " << i << " Value: " << histogram.at(i) << endl;
    }
    
    cout << "RAD calculated." << endl;
    
    // Save the histogram in a csv file
    save_RAD(histogram, bin_size, num);
}


int main(int argc, const char * argv[]) {
    
    //Initialies random seed
    dsfmt_seed(time(NULL));
    
    for (int sample {0}; sample < 10; sample++){

        // Initialise box and positions
        vector<double> box {};
        vector<vector<double>> r {};
        
        //Generate FCC and read positions
        generate_fcc();
        read_data(init_filename, &box, &r);
        
        // Run Monte Carlo Simulation
        int n_particles {0};
        n_particles = int(box[0]);
        int accepted_moves {0};
        for (int step {0}; step < mc_steps; step++){
            for (int n {0}; n < n_particles; n++){
                accepted_moves += move_particle(box, &r);
                }
            
            if (step % output_steps  == 0){
                cout << "================ STEP: " << step << " ================" << endl;
                cout << "Move acceptance rate: " << double(accepted_moves) / (output_steps * n_particles) << endl;
                accepted_moves = 0;
            }
        }
        
        // Calculate and save the radial distribution function
        calculate_RAD(box, r, sample);
        calculate_distances(box, r, sample);
    }
    return 0;
}
