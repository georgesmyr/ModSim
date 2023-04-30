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

const string path {"/Users/georgesmyridis/Documents/Physics/Books-Notes/Graduate/Physics/Modeling_Simulations/Scripts/ModSim/MolecularDynamics"};
const string init_filename = "xyz.dat";

const double SIGMA = 1;
const double EPSILON = 1;
const int NDIM = 3;

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
        (*r).push_back(position);
    }
    infile.close();
    
    assert ((*box)[0] == (*r).size());
}

void print_vector(vector<vector<double>> pos){
    for (auto i = 0; i < pos.size(); ++i){
        cout << "Vector " << i << ": ";
        for (auto j = 0; j < pos[i].size(); ++j){
            cout << pos[i][j] <<  " ";
        }
        cout << endl;
    }
}


double distance(vector<double> *vector1, vector<double> *vector2, vector<double> *box){
    /*
     Description:
     ------------
     This function returns the distance between the centers of two particles, given the boundary conditions of the box.
     
     Input:
     ------
     *vector1: The position of one particle.
     *vector2: The position of the second particle.
     *box: Specs of the box.
     
     Output:
     -------
     The distance.
     */
    
    double minD {0}, dist2 {0};
    for (auto i {0}; i < 3; ++i){
        minD = (*vector1)[i] - (*vector2)[i];
        minD -= (int)(2.0 * minD /  pow((*box)[i+1],2));
        dist2 += minD * minD;
    }
    
    return sqrt(dist2);
}


class MolecularDynamics{
    
private:
    
    double TEMPERATURE = 1; // Instantaneous temperature
    double DELTA_T = 1; // Time step
    double MASS = 1;
    double time = 0; // Current time
    double TMAX = 1000; // Maximum time
    
public:
    
    vector<double> box = {}; // Box parameters
    double Nparticles = 0;
    
    vector<vector<double>> positions = {}; // Current positions
    
    vector<vector<double>> velocities = {}; // Current velocities
    vector<double> total_momentum = {0, 0, 0}; // Total momentum of the system
    
    vector<vector<double>> forces = {}; // Total force acting on each particle at current position
    vector<vector<double>> forcesP = {}; // Total force acting on each particle at future position

    
    double kinetic_energy = 0; // Aggregate kinetic energy of all particles
    double potential_energy = 0; // Aggregate potential energy of all particles
    double total_energy_pp = 0; // Total energy per particle
    
    
    //INITIALISE CONFIGURATION//
    void InitConfig(){
        /*
         Description:
         -----------
         Initialises the configuration of the system. It reads the positions from a '_.dat' file with name init_filename.
         For the velocities (momentums), it randomly chooses them from [-1, 1].
         */
        
        cout << "Initialising configuration.." << endl;
        
        
        // Initialise Positions: Read them from dat file.
        cout << "Importing initial positions.." << endl;
        read_data(init_filename, &box, &positions);
        cout << "Positions initialised." << endl;
        Nparticles = box[0];
        
        // Initialise velocities: Random value from -1 to +1 for each direction
        cout << "Initialising velocities.." << endl;
        for (auto i = 0; i < Nparticles; ++i){
            vector<double> velocity = {2 * (dsfmt_genrand() - 0.5), 2 * (dsfmt_genrand() - 0.5), 2 * (dsfmt_genrand() - 0.5)};
            velocities.push_back(velocity);
            
            for (auto j = 0; j < 3; ++j){
                // Calculate total momentum
                total_momentum[j] += MASS * velocity[j];
                // Calculate total energy
                kinetic_energy += MASS * pow(velocity[j], 2);
            }
        }
        cout << "Velocities initialised." << endl;
        
        // Normalisation
        kinetic_energy /= Nparticles;
        for (unsigned long i = 0; i < 3; ++i){
            total_momentum[i]  /= Nparticles;
        }
        double scale_factor = sqrt(3 * TEMPERATURE / kinetic_energy);
        
        // Shift velocities so that total momentum is zero
        for (auto i = 0; i < Nparticles; ++i){
            for (auto j = 0; j < 3; ++j){
                velocities[i][j] = scale_factor * (velocities[i][j] - total_momentum[j]);
            }
        }
        
        // Initialise list of total force acting on all particles to vanishing.
        cout << "Initialising forces.." << endl;
        for (auto i = 0; i < Nparticles; ++i){
            forces.push_back({0, 0, 0});
            forcesP.push_back({0, 0, 0});
        }
        
        calculate_forces_potential();
        for (auto i = 0; i < Nparticles; ++i){
            for (auto j = 0; j < NDIM; ++j){
                forces[i][j] = forcesP[i][j];
            }
        }
        
        cout << "Forces initialised." << endl;
                
        assert(velocities.size() == Nparticles);
        assert(forces.size() == Nparticles);
        assert(forces.size() == Nparticles);
        
        cout << "Configuration initialised!" << endl;
    }
    
    
     //CALCULATE FORCES & TOTAL ENERGY PER PARTICLE//
    void calculate_forces_potential(){
        /*
         Description:
         -----------
         Calculates the total force acting on all particles.
         */
        
        // Set future forces to zero
        for (auto i = 0; i < Nparticles; ++i){
            for (auto j = 0; j < NDIM; ++j){
                forcesP[i][j] = 0;
            }
        }
        
        for(auto index1 = 0; index1 < Nparticles - 1; ++index1){
            for(auto index2 = index1 + 1; index2 < Nparticles; ++index2){
                
                // Caclulate displacement vector
                vector<double> dist_vec = {0, 0, 0};
                for (auto i = 0; i < NDIM; ++i){
                    dist_vec[i] = positions[index1][i] - positions[index2][i];
                }
                // Calculate unitary vector in displacement direction
                double dist = distance(&positions[index1], &positions[index2], &box);
                for (auto i = 0; i < NDIM; ++i){
                    dist_vec[i] /= dist;
                }
                
                // Calculate force vector
                double force_magnitude = (24 * EPSILON / SIGMA) * (2 * pow(dist, 13) - pow(dist,7));
                // Add it to the total forces
                for(auto i = 0; i < NDIM; ++i){
                    forcesP[index1][i] += force_magnitude * dist_vec[i];;
                    forcesP[index2][i] -= force_magnitude * dist_vec[i];;
                }
                
                // Add potential energy
                potential_energy += 4 * EPSILON * (pow(SIGMA / dist, 12) - pow(SIGMA / dist, 6));
            }
        }
        
        assert(forcesP.size() == Nparticles);
    }
    
    
    //INTEGRATE EQUATIONS OF MOTION//
    void integrate(bool NVT = false){
        /*
         Description:
         -----------
         Integrates Newton's equations of motion with the 'velocity Verlet' algorithm.
         */
        
        // Set kinetic energy and total momentum to 0
        kinetic_energy = 0;
        for (auto i = 0; i < Nparticles; ++i){
            total_momentum[i] = 0;
        }
        if (NVT == false){
            for (auto i = 0; i < Nparticles; ++i){
                for (auto j = 0; j < NDIM; ++j){
                    // Update positions
                    positions[i][j] += velocities[i][j] * DELTA_T + pow(DELTA_T, 2) * forces[i][j] / (2 * MASS);
                    // Calculate future forces and potential energy
                    calculate_forces_potential();
                    // Update velocities
                    velocities[i][j] += (forces[i][j] + forcesP[i][j]) / (2 * MASS);
                    // Add in kinetic energy
                    kinetic_energy += pow(velocities[i][j], 2) / 2;
                    // Add total momentum
                    total_momentum[j] += MASS * velocities[i][j];
                    
                    // Update forces
                    forces[i][j] = forcesP[i][j];
                }
            }
        } else {
            for (auto i = 0; i < Nparticles; ++i){
                for (auto j = 0; j < Nparticles; ++j){
                    velocities[i][j] += DELTA_T * forces[i][j] / 2;
                }
            }
            
        }
    }
    
    
    void save_data(){
        /*
         Description:
         ------------
         Saves the positions in a .dat file.
         Saves the total energy per particle in a csv file
         */
        
        
    }
    
    void run_simulation(){
        /*
         Description:
         ------------
         Runs the simulation.
         */
        
        InitConfig();
        while (time < TMAX){
            integrate(NVT = false);
            time += DELTA_T;
        }
        
        
        
        
        
        
        
        
    }
};

int main(int argc, const char * argv[]) {
    dsfmt_seed(time(NULL));

    MolecularDynamics md;
    md.InitConfig();
    print_vector(md.forcesP);
    

    return 0;
}
