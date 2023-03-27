

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <algorithm>
#include <unordered_set>

#include "mt19937.h"

using namespace std;

const int N = 9;
const int K = sqrt(9);

class SudokuSolver{
    
private:
    
    double BETA {1};
    const double COOLING_RATE {0.000001};
    int board[N][N] = {
        {0, 0, 9, 0, 7, 0, 0, 0, 2},
        {0, 0, 0, 0, 0, 8, 4, 0, 0},
        {0, 0, 0, 2, 0, 0, 5, 7, 0},
        {0, 0, 0, 7, 0, 0, 9, 0, 0},
        {5, 0, 0, 0, 0, 0, 0, 0, 3},
        {0, 0, 6, 0, 0, 1, 0, 0, 0},
        {0, 8, 2, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 3, 0, 0, 0, 0, 0},
        {6, 0, 0, 0, 8, 0, 0, 0, 0}
    };
        
    vector<vector<int>> init_ConstPos(){
        /*
         Description:
         ------------
         Gets the positions of non-zero elements of the starting board.
         These elements cannot be altered during the solution of the game.
         */
        
        vector<vector<int>> constPos = {};
        for (auto row = 0; row < N; ++row){
            for (auto col = 0; col < N; ++col){
                if (board[row][col] != 0){
                    constPos.push_back({row,col});
                }
            }
        }
        return constPos;
    }
    
    const vector<vector<int>> constPos = init_ConstPos();
    int costFunction = 0;
    vector<int> randomPos = {};
    int randomEntry = 0;
    bool solved = false;
    
public:
    
    void print_board() {
        for (int i = 0; i < N; i++) {
            if (i % K == 0 && i != 0) {
                cout << "---------------------" << endl;
            }
            for (int j = 0; j < N; j++) {
                if (j % K == 0 && j != 0) {
                    cout << "| ";
                }
                cout << board[i][j] << " ";
            }
            cout << endl;
        }
    }
    
    void GenRandomPosition(){
        
        // Pick a random position, not one of the initial ones, and put a random integer in it.
        do {
            randomPos = {int (dsfmt_genrand() * N), int (dsfmt_genrand() * N)};
        } while(find(constPos.begin(), constPos.end(), randomPos) != constPos.end());
        
        assert (randomPos[0] >= 0 and randomPos[0] <= 8);
        assert (randomPos[1] >= 0 and randomPos[1] <= 8);
    }
    
    void GenRandomEntry(){
        
        randomEntry = int (dsfmt_genrand() * N + 1);
        assert (randomEntry >= 1 and randomEntry <= 9);
    }
    
    int getCostFunction(){
        
        int duplicates  = 0;
        int value = 0;
        // Find duplicates in all rows.
        for (int row = 0; row < N; ++row){
            unordered_set<int> values;
            for (int col = 0; col < N; ++col){
                value = board[row][col];
                if (value != 0 and !values.insert(value).second){
                    duplicates += 1;
                }
            }
        }
        
        // Find duplicates in all columns
        for (int col = 0; col < N; ++col){
            unordered_set<int> values;
            for (int row = 0; row < N; ++row){
                value = board[row][col];
                if (value != 0 and !values.insert(value).second){
                    duplicates += 1;
                }
            }
        }
        return duplicates;
    }
    
    void AcceptOrReject(){
        
        int oldEntry = 0;
        int randomX = randomPos[0];
        int randomY = randomPos[1];
        oldEntry = board[randomX][randomY];
        board[randomX][randomY] = randomEntry;
        
        int newCostFunction = getCostFunction();
        // Accept or reject the new solution
        int dCF = newCostFunction - costFunction;
        //cout << "Cost function difference: " << dCF << endl;
        if (dCF < 0 or dsfmt_genrand() < exp(- BETA * dCF)){
            //cout << "Accepted!" << endl;
            costFunction += dCF;
        }else{
            //cout << "Rejected!" << endl;
            board[randomX][randomY] = oldEntry;
        }
    }
    
    void cooldown(){
        BETA *= (1 - COOLING_RATE);
    }
    
    bool is_solved(){
        
        // Check if the cost function is zero
        if (costFunction != 0){
            return false;
        }
        
        // Check if there's any zero entry
        for (int row = 0; row < N; ++row){
            for (int col = 0; col < N; ++col){
                if (board[row][col] == 0){
                    return false;
                }
            }
        }
        int sum = 0;
        // Check if the sum in all rows is equal to 45
        for (int row = 0; row < N; ++row){
            sum = 0;
            for (int col = 0; col < N; ++col){
                sum += board[row][col];
            }
            if (sum != 45){
                return false;
            }
        }
        // Check if the sum in all columns is equal to 45
        for (int col = 0; col < N; ++col){
            sum = 0;
            for (int row = 0; row < N; ++row){
                sum += board[row][col];
            }
            if (sum != 45){
                return false;
            }
        }
        return true;
    }
    
    void solve(){
//        int i = 0;
        do{
            GenRandomPosition();
            GenRandomEntry();
            AcceptOrReject();
            cooldown();
//            i += 1;
        }while (not is_solved());
        
        cout << "Cost function: " << costFunction <<  endl;
        print_board();
    }
    
    
};

//vector<list<int>> init(){
//    vector<int[2]> temp = {};
//    temp.push_back([1,2]);
//    return temp;
//}

int main(int argc, const char * argv[]) {
    
    dsfmt_seed( time (NULL)); // Initialise random seed.

    SudokuSolver sudoku;
    sudoku.print_board();
    
    sudoku.solve();

    
    
    
    return 0;
}
