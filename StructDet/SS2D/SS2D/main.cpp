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
int Nx {3}, Ny {3};
const int N {2 * Nx * Ny};
const double RADIUS = 1.0, LATTICE_SPACING = 4.0 * RADIUS;

//SYSTEM VARIABLES
double POSITIONS[18][4] = {{0}};
double BOX[4] = {0};

//POTENTIAL PARAMETERS
const double THRESH1 {2 * RADIUS};
const double LAMBDA {1.5};
const double THRESH2 {LAMBDA * THRESH1};
const double POTENTIAL {1};

//MONTE CARLO SIMULATION PARAMETERS
int STEPS {10000};
double DELTA_MAX {0.1 * LATTICE_SPACING};


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
    double Lx {Nx * LATTICE_SPACING};
    double Ly {Ny * LATTICE_SPACING};
    
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
            outfile << i * LATTICE_SPACING << " " << j * LATTICE_SPACING << " " << 0 << " " << RADIUS << endl;
            outfile << (i + 0.5) * LATTICE_SPACING << " " << (j + 0.5) * LATTICE_SPACING << " " << 0 << " " << RADIUS << endl;
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

void print_positions(){
    /*
     DESCRIPTION: Prints out all the positions.
     */
    cout << "---- POSITIONS ----" << endl;
    for(int i = 0; i < N; ++i){
        cout << i << ") ";
        for(int j = 0; j < NDIM; ++j){
            cout << POSITIONS[i][j] << " ";
        }
        cout << endl;
    }
    cout << "----------------" << endl;
}

void print_box(){
    /*
     DESCRIPTION: Prints out the information of the box, i.e., dimensions and number of particles.
     */
    cout << "BOX: ";
    for(int i = 0; i < 4; ++i){
        cout << BOX[i] << " ";
    }
}

double distance(unsigned int index1, unsigned int index2){
    /*
     DESCRIPTION: Returns the distance between the centers of two particles, given the boundary conditions of the box.
     INPUT: index1, index2: The indices of the two particles.
     OUTPUT: The distance.
     */

    double minD {0}, dist2 {0}, periodicFactor {0};
    for (auto i {0}; i < NDIM; ++i){
        minD = POSITIONS[index1][i] - POSITIONS[index2][i];
        periodicFactor = round(minD / BOX[i]);
        minD -= periodicFactor * BOX[i];
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
        printf("(%d,%d) Distace = %f\n", idx1, idx2, dist);
        return INFINITY;
    }
    else if(dist < THRESH2){
        return POTENTIAL;
    }
    else{
        return 0;
    }
}

double calc_energy_contr(unsigned int idx){
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

double calc_energy(){
    /*
     DESCRIPTION: Calculates the enrgy of the system.
     */
    double energy {0};
    for(int i = 0; i < BOX[3] - 1; ++i){
        for(int j = i + 1; j < BOX[3]; ++j){
            energy += potential(i,j);
        }
        if(energy == INFINITY){
            return energy;
        }
    }
    return energy;
}


int move_particle(){
    
    
    //CHOOSE RANDOM PARTICLE
    int idx = (int)(dsfmt_genrand() * N);
    cout << "Index: " << idx << endl;
    for(int i = 0; i < NDIM; ++i){
        cout << POSITIONS[idx][i] << " ";
    }
    
    //CALCULATE CANDIDATE POSITION FOR THE CHOSEN PARTICLE
    double new_pos[4] = {0,0,0,1};
    for(auto i {0}; i < NDIM; ++i){
        new_pos[i] = POSITIONS[idx][i] + (2 * dsfmt_genrand() - 1) * DELTA_MAX;
    
        // ACCOUNT FOR PERIODIC BOUNDARY CONDITIONS
        if(new_pos[i] < 0){
            new_pos[i] += BOX[i];
        }else if(new_pos[i] > BOX[i]){
            new_pos[i] -= BOX[i];
        }
        assert(new_pos[i] >= 0 and new_pos[i] <= BOX[i]);
    }
    
    //CALCULATE CHANGE IN ENERGY
    double contr1 {calc_energy_contr(idx)};
    double eb {calc_energy()};
    //Change the position
    for(int i = 0; i < NDIM; ++i){
        POSITIONS[idx][i] = new_pos[i];
    }
    
    double ea {calc_energy()};
    double contr2 {calc_energy_contr(idx)};
    cout << "Energy before: " << contr1 << endl;
    cout << "Energy after: " << contr2 << endl;
    cout << "DEE: " << contr2-contr1 << endl;
    cout << "DE: " << ea - eb << endl;
    
    for(int i = 0; i < NDIM; ++i){
        cout << POSITIONS[idx][i] << " ";
    }cout << endl;
            
    
    

    
    return 1;
}


int main(int argc, const char * argv[]) {
    
    dsfmt_seed(time(NULL));
    //GENERATE FCC
//    generate_fcc();
    
//    vector<vector<double>> positions = {};
//    vector<double> box = {};
//    read_data(&positions, &box);
    
    read_data();
    print_positions();
    move_particle();
    print_positions();
    
    
    
    
    
    return 0;
}




