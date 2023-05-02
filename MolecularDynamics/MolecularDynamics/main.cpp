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

const string path {"/Users/georgesmyridis/Desktop/Trading/ModSim/MolecularDynamics"};
const string init_filename = "xyz.dat";

void print_vector(vector<vector<double>> pos){
    for (auto i = 0; i < pos.size(); ++i){
        cout << "Vector " << i << ": ";
        for (auto j = 0; j < pos[i].size(); ++j){
            cout << pos[i][j] <<  " ";
        }
        cout << endl;
    }
}

// Lattice parameters
const int Nx = 2;
const int Ny = 2;
const int Nz = 2;
const double d = 1.0; // Atomic radius
const double l = 2 * d; // Lattice spacing for FCC lattice


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
        
    for (int i {0}; i < Nx - 1; i++){
        for (int j {0}; j < Ny - 1; j++){
            for (int k {0}; k < Nz - 1; k++){
                outfile << i * l + d / 2 << " " << j * l + d / 2 << " " << k * l + d / 2<< " " << d << endl;
                outfile << i * l + d / 2 << " " << (0.5 + j) * l + d / 2 << " " << (0.5 + k) * l + d / 2 << " " << d << endl;
                outfile << (0.5 + i) * l + d / 2 << " " << (0.5 + j) * l + d / 2 << " " << k * l + d / 2 << " " << d << endl;
                outfile << (0.5 + i) * l + d / 2 << " " << j * l + d / 2 << " " << (0.5 + k) * l + d / 2 << " " << d << endl;
                
            }
        }
    }
    cout << "FCC generated." << endl;
}

class MolecularDynamics{
    
private:
    
    const int NDIM = 3; // Number of spatial dimensions
    double TEMPERATURE = 1; // Instantaneous temperature
    const double DENSITY = 0.1; // Density of the system

    const double TMAX = 1000; // Maximum time
    const double DELTA_T = 0.001; // Time step
    double time = 0; // Current time
    
    double RCUT = 0;
    double ECUT = 0;
    
    double MAX_DISTANCE = 0; // Max distance in the box, equal to the diagonal
    
    const string IMPORT_PATH = path + "/positions/" + init_filename;
    const string STATISTICS_PATH = path + "/statistics/";
    
    
    vector<double> box = {}; // Box parameters
    double NPARTICLES = 0;
    
public:

    vector<vector<double>> positions = {}; // Current positions
    
    vector<vector<double>> velocities = {}; // Current velocities
    vector<double> total_momentum = {0, 0, 0}; // Total momentum of the system
    vector<double> total_momentum_pp = {0, 0, 0}; // Total momentum per particle
    
    vector<vector<double>> forces = {}; // Total force acting on each particle at current position
    vector<vector<double>> forcesP = {}; // Total force acting on each particle at future position

    
    double kinetic_energy = 0; // Aggregate kinetic energy of all particles
    double kinetic_energy_pp = 0; // Kinetic energy per particle
    double potential_energy = 0; // Aggregate potential energy of all particles
    double potential_energy_pp = 0; // Potential energy per particle
    double total_energy = 0; // Aggregate total energy
    double total_energy_pp = 0; // Total energy per particle
    
    void print_energies(){
        
        cout << "-----------------------------" << endl;
        cout << "Kinetic energy: " << kinetic_energy << endl;
        cout << "Potential energy: " << potential_energy << endl;
        cout << "Total energy: " << total_energy << endl;
    }
    
    void import_positions(){
        /*
         Description: Imports positions from dat file in FILE_PATH.
        */

        ifstream infile;
        infile.open(IMPORT_PATH);
        if (not infile.is_open()){
            cerr << "Could not read the file." << endl;
            return;
        }
        double Lx {0}, Ly {0}, Lz {0}, x {0}, y {0}, z {0}, d {0};
        infile >> NPARTICLES;
        infile >> Lx >> Ly >> Lz;
        box = {Lx, Ly, Lz};
        for (int i {0}; i < NPARTICLES; i++){
            infile >> x >> y >> z >> d;
            vector <double> position {x,y,z,d};
            positions.push_back(position);
        }
        infile.close();
        assert (NPARTICLES == positions.size());
    }
    
    //INITIALISE CONFIGURATION//
    void set_RCUT(){
        /*
         DESCRIPTION: Sets RCUT and ECUT
         */
        RCUT = MAX_DISTANCE / 2;
        ECUT = 4 * (pow(RCUT, -12) - pow(RCUT, -6));
        cout << " - RCUT: " << RCUT << endl;
        cout << " - ECUT: " << ECUT << endl;
    }
    
    void set_density(){
        /*
         DESCRIPTION: Sets the density of the system by scaling the box and the positions.
         */
        double volume = 1;
        cout << " - Box dimensions: ";
        for(auto i = 0; i < NDIM; ++i){
            cout << box[i] << "  ";
            volume *= box[i];
        }
        cout << endl;
        cout << " - Volume: " << volume << endl;
        double target_volume = NPARTICLES / DENSITY;
        cout << " - Target volume: " << target_volume << endl;
        double scale_factor = pow(target_volume / volume, 1.0 / NDIM);
        cout << " - Scale factor: " << scale_factor << endl;
        
        for (auto i = 0; i < NPARTICLES; ++i){
            for(auto j = 0; j < NDIM; ++j){
                positions[i][j] *= scale_factor;
            }
        }
        MAX_DISTANCE = 0;
        cout << " - New box dimensions: ";
        for (auto i = 0; i < NDIM; ++i){
            box[i] *= scale_factor;
            MAX_DISTANCE += pow(box[i], 2);
            cout << box[i] << "  ";
        }
        cout << endl;
        MAX_DISTANCE = sqrt(MAX_DISTANCE);
    }
    
    void InitConfig(){
        /*
         DESCRIPTION: Initialises the configuration of the system. It reads the positions from a '_.dat' file with name init_filename.
         For the velocities (momentums), it randomly chooses them from [-1, 1].
         */
        
        cout << "INITIALISING CONFIGURATION.." << endl;
        cout << "-----------------------------" << endl;
        
        // Initialise Positions: Read them from dat file.
        cout << "Initialising positions.." << endl;
        cout << " - Importing.." << endl;
        import_positions();
        cout << "Positions initialised." << endl;
        cout << "-----------------------------" << endl;
        // Set Density
        cout << "Setting density.." << endl;
        set_density();
        cout << "Density set." << endl;
        cout << "-----------------------------" << endl;
        
        cout << "Setting RCUT and ECUT.." << endl;
        set_RCUT();
        cout << "Set." << endl;
        cout << "-----------------------------" << endl;

        // Initialise velocities: Random value from -1 to +1 for each direction
        cout << "Initialising velocities.." << endl;
        for (auto i = 0; i < NPARTICLES; ++i){
            vector<double> velocity = {2 * (dsfmt_genrand() - 0.5), 2 * (dsfmt_genrand() - 0.5), 2 * (dsfmt_genrand() - 0.5)};
            velocities.push_back(velocity);

            for (auto j = 0; j < 3; ++j){
                // Calculate total momentum
                total_momentum_pp[j] += velocity[j] / NPARTICLES;
            }
        }
        // Shift velocities so that total momentum is zero
        for (auto i = 0; i < NPARTICLES; ++i){
            for (auto j = 0; j < NDIM; ++j){
                velocities[i][j] -= total_momentum_pp[j];
            }
        }
        total_momentum_pp = {0, 0, 0};
        
        // Scale velocities so that kinetic energy matches the temperature
        for (auto i = 0; i < NPARTICLES; ++i){
            for (auto j = 0; j < NDIM; ++j){
                kinetic_energy += pow(velocities[i][j], 2);
            }
        }
        kinetic_energy /= 2;
        double scale_factor = sqrt(1.5 * TEMPERATURE / kinetic_energy);
        for (auto i = 0; i < NPARTICLES; ++i){
            for (auto j = 0; j < NDIM; ++j){
                velocities[i][j] *= scale_factor;
            }
        }
        kinetic_energy = 1.5 * TEMPERATURE;
             
        cout << "Velocities initialised." << endl;
        cout << "-----------------------------" << endl;

        // Initialise list of total force acting on all particles to vanishing.
        cout << "Initialising forces.." << endl;
        for (auto i = 0; i < NPARTICLES; ++i){
            forces.push_back({0, 0, 0});
            forcesP.push_back({0, 0, 0});
        }

        update_forces();
        for (auto i = 0; i < NPARTICLES; ++i){
            for (auto j = 0; j < NDIM; ++j){
                forces[i][j] = forcesP[i][j];
            }
        }
        cout << "Forces initialised." << endl;
        cout << "-----------------------------" << endl;
        assert(velocities.size() == NPARTICLES);
        assert(forces.size() == NPARTICLES);
        assert(forcesP.size() == NPARTICLES);
        cout << "CONFIGURATION INITIALISED!" << endl;
        
        total_energy = potential_energy + kinetic_energy;
        print_energies();
    }
    
    double distance(unsigned int index1, unsigned int index2){
        /*
         DESCRIPTION: Returns the distance between the centers of two particles, given the boundary conditions of the box.
         
         INPUT: index1, index2: The indices of the two particles.
         
         OUTPUT: The distance.
         */
        
        double minD {0}, dist2 {0};
        for (auto i {0}; i < 3; ++i){
            minD = positions[index1][i] - positions[index2][i];
            minD -= (int)(2.0 * minD /  pow(box[i],2));
            dist2 += minD * minD;
        }
        return sqrt(dist2);
    }
    
   void update_forces(){
       /*
        DESCRIPTION: Calculates the total force acting on all particles.
        */
       
       // Set potential energy and future forces to zero
       potential_energy = 0;
       for (auto i = 0; i < NPARTICLES; ++i){
           for (auto j = 0; j < NDIM; ++j){
               forcesP[i][j] = 0;
           }
       }
       
       double dist = 0;
       double force_magnitude = 0;
       vector<double> dist_vec = {0, 0, 0};
       for(auto index1 = 0; index1 < NPARTICLES - 1; ++index1){
           for(auto index2 = index1 + 1; index2 < NPARTICLES; ++index2){
               
               // Caclulate displacement vector
               dist_vec = {0, 0, 0};
               dist = 0;
               for (auto i = 0; i < NDIM; ++i){
                   dist_vec[i] = positions[index1][i] - positions[index2][i];
                   // Nearest image convention
                   if (dist_vec[i] > box[i] / 2){
                       dist_vec[i] -= box[i];
                   }else if (dist_vec[i] < - box[i] / 2){
                       dist_vec[i] += box[i];
                   }
                   dist += pow(dist_vec[i], 2);
                   assert(dist_vec[i] >= - box[i] and dist_vec[i] <= box[i]);
               }
               dist = sqrt(dist);
               assert(dist >= 0 and dist <= MAX_DISTANCE);
               
            //Calculate unitary vector in displacement direction
               for (auto i = 0; i < NDIM; ++i){
                   dist_vec[i] /= dist;
               }
               // Calculate force vector
               force_magnitude = 24 * (2 * pow(dist, -13) - pow(dist, -7));
               // Add it to the total forces
               for(auto i = 0; i < NDIM; ++i){
                   forcesP[index1][i] += force_magnitude * dist_vec[i];;
                   forcesP[index2][i] -= force_magnitude * dist_vec[i];;
               }

               // Add potential energy
               potential_energy += 4 * (pow(dist, -12) - pow(dist, -6)) - ECUT;
           }
       }
       assert(forcesP.size() == NPARTICLES);
   }
    
    
    //INTEGRATE EQUATIONS OF MOTION//
    void update_positions(){
        /*
         DESCRIPTION: Updates positions with the 'velocity Verlet' algorithm.
         */
        
        // Set kinetic energy and total momentum to 0
        kinetic_energy = 0;
        for (auto i = 0; i < NPARTICLES; ++i){
            total_momentum[i] = 0;
        }

        double x_j = 0; // Temporary value for the position component along j direction
        unsigned long multiple {0};
        for (auto i = 0; i < NPARTICLES; ++i){
            for (auto j = 0; j < NDIM; ++j){
                // Update positions
                positions[i][j] += velocities[i][j] * DELTA_T + pow(DELTA_T, 2) * forces[i][j] / 2;
                x_j = positions[i][j];
                if (positions[i][j] > box[j]){
                    multiple = (unsigned long)(x_j / box[j]);
                    positions[i][j] = x_j - multiple * box[j];
                }else if(positions[i][j] < 0){
                    multiple = (unsigned long)(abs(x_j) / box[j]);
                    positions[i][j] = x_j + (multiple + 1) * box[j];

                }
                assert(positions[i][j] >= 0 and positions[i][j] <= box[j]);
            }
        }
    }
    
    void update_velocities(){
        /*
         DESCRIPTION: Update velocities, kinetic energy and total linear momentum
         for updated positions and forces.
         */
        
        kinetic_energy = 0;
        total_momentum = {0, 0, 0};
        for (auto i = 0; i < NPARTICLES; ++i){
            for (auto j = 0; j < NDIM; ++j){
                // Updata velocities
                velocities[i][j] += DELTA_T * (forces[i][j] + forcesP[i][j]) / 2;
                // Aggregate kinetic energy
                kinetic_energy += pow(velocities[i][j], 2);
                // Aggregate total linear momentum
                total_momentum[j] += velocities[i][j];
            }
        }
        kinetic_energy /= 2;
    }
    
    void run_simulation(){
        
        // Open energies.csv file
        ofstream outfile_en(STATISTICS_PATH + "energies.csv");
        if (!outfile_en.is_open()) {
            cerr << "Error: Unable to open output file" << endl;
            return;
        }
        // Open diffusion.csv file
        ofstream outfile_dif(STATISTICS_PATH + "diffusion.csv");
        if (!outfile_dif.is_open()) {
            cerr << "Error: Unable to open output file" << endl;
            return;
        }
        
        outfile_en << "time,kinetic_energy_pp,potential_energy_pp,total_energy_pp,temperature" << endl;
        outfile_dif << "time," << endl;
        
        // Initialise configuration
        time = 0;
        InitConfig();
        while(time < TMAX){
            // Update positions
            update_positions();
            // Update forces and potential for updated positions
            update_forces();
            // Update velocities and kinetic energy
            update_velocities();
            // Update total energy
            total_energy = kinetic_energy + potential_energy;
            // Make future forces current
            for (auto i = 0; i < NPARTICLES; ++i){
                for (auto j = 0; j < NDIM; ++j){
                    forces[i][j] = forcesP[i][j];
                }
            }
            
            // Save energies
            outfile_en << time << "," << kinetic_energy / NPARTICLES << "," << potential_energy / NPARTICLES << "," << total_energy / NPARTICLES;
            outfile_en << "," << TEMPERATURE << endl;
            
            // Calculate and save self diffusion coefficient
            
            
            
            // Update time
            time += DELTA_T;
        }
        
        outfile_en.close();
        outfile_dif.close();

    }
    
};



int main(int argc, const char * argv[]) {
    dsfmt_seed(time(NULL));
//    generate_fcc();
    MolecularDynamics md;
    md.run_simulation();
    md.print_energies();

    return 0;
}
