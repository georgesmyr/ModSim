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

const string path {"/Users/georgesmyridis/Desktop/Trading/ModSim/StructDet/SS2D/positions2D/"};

// FCC PARAMETERS
const int NDIM = 2;
int Nx {10}, Ny {10}; //change also POSITIONS size!!!
int N {2 * Nx * Ny};
const double RADIUS = 1.0;
const double DIAMETER = 2 * RADIUS;
double LATTICE_SPACING = 4.0 * RADIUS;

//SYSTEM VARIABLES
double POSITIONS[200][4] = {{0}};
double BOX[4] = {0};

//POTENTIAL PARAMETERS
const double THRESH1 {2 * RADIUS};
const double LAMBDA {1.5};
const double THRESH2 {LAMBDA * THRESH1};
const double POTENTIAL {2};

//MONTE CARLO SIMULATION PARAMETERS
int STEPS {30000};
double BETA {1};
const double DISPLACEMENT_FRACTION {0.05};
double DELTA_MAX {DISPLACEMENT_FRACTION * LATTICE_SPACING};
const unsigned int OUTPUT_STEPS {500};


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
     xy.dat file in the directory (path). The first line has the total number of particles in the box. The
     second line has the coordinates of the box corner in x-direction, and the next for the y- and z- directions
     respectively. The rest of the lines have the x-, y-, and z- coordinates of the hard sheres, and their diameter.
    */
    
    cout << "Generating FCC..." << endl;
    int N {2 * Nx * Ny};
    cout << N << endl;
    double Lx {Nx * LATTICE_SPACING};
    double Ly {Ny * LATTICE_SPACING};
    
    cout << Lx << endl;
    cout << Ly << endl;
    
    ofstream outfile( path + "xy.dat");
    if (!outfile.is_open()) {
        cerr << "Error: Unable to open output file" << endl;
        return;
    }
    outfile << N << endl;
    outfile << Lx << endl;
    outfile << Ly << endl;
    outfile << 0 << endl;
        
    for (int i {0}; i < Nx; i++){
        for (int j {0}; j < Ny; j++){
            outfile << i * LATTICE_SPACING << " " << j * LATTICE_SPACING << " " << 0 << " " << DIAMETER << endl;
            outfile << (i + 0.5) * LATTICE_SPACING << " " << (j + 0.5) * LATTICE_SPACING << " " << 0 << " " << DIAMETER << endl;
        }
    }
    cout << "FCC generated." << endl;
    outfile.close();
}

void read_data(){
    /*
     DESCRIPTION: Reads all the initial positions from a dat file.
     */
    
    ifstream infile;
    infile.open(path + "xy.dat");
    if (not infile.is_open()){
        cerr << "Could not read the file." << endl;
        return;
    }
    infile >> BOX[3] >> BOX[0] >> BOX[1] >> BOX[2];
    for(int i = 0; i < N; ++i){
        for(int j = 0; j < 4; ++j){
            infile >> POSITIONS[i][j];
        }
    }
}

double distance(unsigned int index1, unsigned int index2){
    /*
     DESCRIPTION: Returns the distance between the centers of two particles, given the boundary conditions of the box.
     INPUT: index1, index2: The indices of the two particles.
     OUTPUT: The distance.
     */

    double minD {0}, dist2 {0}, factor {0};
    for (auto i {0}; i < NDIM; ++i){
        minD = POSITIONS[index1][i] - POSITIONS[index2][i];
        factor = (int)(2.0 * minD /  BOX[i]);
        minD -= factor * BOX[i];
        dist2 += minD * minD;
    }
    
    return sqrt(dist2);
}

double potential(unsigned int idx1, unsigned int idx2){
    /*
     DESCRIPTION: Returns the potential energy between the particles with idx1 and idx2.
     */
    double dist {distance(idx1, idx2)};
    if(dist < THRESH1){
        return INFINITY;
    }
    else if(dist < THRESH2){
        return POTENTIAL;
    }
    else{
        return 0;
    }
}

double energy_contr(unsigned int idx){
    /*
     DESCRIPTION: Calculates the energy contribution for the particle idx.
     */
    double contribution {0};
    for(auto i = 0; i < N; ++i){
        if(i != idx){
            contribution += potential(idx, i);
        }
        if(contribution == INFINITY){
            return contribution;
        }
    }
    return contribution;
}


int move_particle(){
    
    
    //CHOOSE RANDOM PARTICLE
    int idx = (int)(dsfmt_genrand() * N);
    
    //CALCULATE CURRENT ENERGY CONTRIBUTION FRO MPARTICLE idx
    double contr_before {energy_contr(idx)};
    
    //CALCULATE RANDOM DISPLACEMENTS AND CHANGE POSITION
    double displacements[NDIM] = {0,0};
    for(auto i {0}; i < NDIM; ++i){
        displacements[i] = (2 * dsfmt_genrand() - 1) * DELTA_MAX;
        POSITIONS[idx][i] += displacements[i];
    
        //Account for periodic boundary conditions
        if(POSITIONS[idx][i] < 0){
            POSITIONS[idx][i] += BOX[i];
        }else if(POSITIONS[idx][i] > BOX[i]){
            POSITIONS[idx][i] -= BOX[i];
        }
        assert(POSITIONS[idx][i] >= 0 and POSITIONS[idx][i] <= BOX[i]);
    }
    
    //CALCULATE CHANGE IN ENERGY
    double contr_after {energy_contr(idx)};
    double dE {contr_after - contr_before};
        
    //ACCEPT OR REJECT MOVE
    //Accept and return 1
    if(dE != INFINITY and (dE <= 0 or dsfmt_genrand() < exp(- BETA * dE))){
        return 1;
    }
    //Reject, reverse the displacements and return 0.
    for(int i = 0; i < NDIM; ++i){
        POSITIONS[idx][i] -= displacements[i];
        
        //Account for periodic boundary conditions
        if(POSITIONS[idx][i] < 0){
            POSITIONS[idx][i] += BOX[i];
        }else if(POSITIONS[idx][i] > BOX[i]){
            POSITIONS[idx][i] -= BOX[i];
        }
        assert(POSITIONS[idx][i] >= 0 and POSITIONS[idx][i] <= BOX[i]);
    }
    
    return 0;
}

void check_overlaps(){
    /*
     DESCRIPTION: The function checks if there are overlaps between the particles.
     */
    for(int i = 0; i < N - 1; ++i){
        for(int j = i + 1; j < N; ++j){
            assert(distance(i,j) > 2 * RADIUS);
        }
    }
}

void write_data(unsigned int num){
    /*
     DESCRIPTION: This function writes the positions of the particles in a dat file named xy(num).dat.
     */
    ofstream outfile(path + "xy" + to_string(num) + ".dat");
    if (!outfile.is_open()){
        cerr << "Error: Unable to open output file." << endl;
        return;
    }
        outfile << BOX[3] << endl;
    for(int i = 0; i < 3; ++i){
        outfile << BOX[i] << endl;
    }
    for(auto i = 0; i < N; ++i){
        for(int j = 0; j < 4; ++j){
            outfile << POSITIONS[i][j] << " ";
        }
        outfile << endl;
    }
    outfile.close();
}

void print_box(){
    for(int i = 0; i < 4; ++i){
        cout << BOX[i] << " ";
    }
    cout << endl;
}

void print_positions(){
    for(int i = 0; i < N; ++i){
        cout << i << ") ";
        for(int j = 0; j < 4; ++j){
            cout << POSITIONS[i][j] << " ";
        }
        cout << endl;
    }
}

void set_packing_fraction(double packing_fraction){

    double area {BOX[0] * BOX[1]};
    double particle_area {M_PI * pow(RADIUS,2)};

    double target_area = (BOX[3] * particle_area)/ packing_fraction;
    double scale_factor = pow(target_area/area, 1.0/2.0);
    cout << "Scale factor: " << scale_factor << endl;
    for (int i {0}; i < BOX[3]; i++){
        for (int j{0}; j < NDIM; j++){
            POSITIONS[i][j] *= scale_factor;
        }
    }
    for (int i{0}; i < NDIM; i++){
        BOX[i] *= scale_factor;
    }
    LATTICE_SPACING *= scale_factor;
    DELTA_MAX = DISPLACEMENT_FRACTION * LATTICE_SPACING;
}


int main(int argc, const char * argv[]) {
    
    dsfmt_seed(time(NULL));
    //GENERATE FCC
    generate_fcc();
    read_data();
    set_packing_fraction(0.5);
    check_overlaps();

    float accepted {0};
    for(auto step = 0; step < STEPS; ++step){
        accepted += move_particle();
        check_overlaps();
        if(step % OUTPUT_STEPS == 0){
            printf("Acceptance ratio: %f\n", accepted / OUTPUT_STEPS);
            write_data(step / OUTPUT_STEPS);
            accepted = 0;
        }
    }

    
    
    
    return 0;
}




