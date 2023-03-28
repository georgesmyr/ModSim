

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <random>

#include "mt19937.h"

using namespace std;

const int BOARD_SIZE = 9;
const int N = sqrt(BOARD_SIZE);

class SudokuSolver{
    
private:
    
    double BETA {1};
    const double COOLING_RATE {0.99};
    int board[BOARD_SIZE][BOARD_SIZE] = {
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
    
    unordered_set<int> GetAllPossibleValues(){
        /*
         Description:
         ------------
         For given size of sudoku instance, it returns an unordered set
         with all the possible values that one can place in the cells.
         */
        unordered_set<int> values = {};
        for (int i = 1; i < BOARD_SIZE + 1; ++i){
            assert(i >= 1 && i <= BOARD_SIZE);
            values.emplace(i);
        }
        return values;
    }
    unordered_set<int> AllPossibleValues = GetAllPossibleValues();
    
    
    vector<vector<int>> GetFixedPositions(){
        /*
         Description:
         ------------
         Gets the positions of non-zero elements of the starting board.
         These elements cannot be altered during the solution of the game.
         */
        
        vector<vector<int>> FixedPos = {};
        for (auto row = 0; row < BOARD_SIZE; ++row){
            for (auto col = 0; col < BOARD_SIZE; ++col){
                if (board[row][col] != 0){
                    FixedPos.push_back({row,col});
                }
            }
        }
        return FixedPos;
    }
    const vector<vector<int>> FixedPositions = GetFixedPositions();

    
public:
    
    void print_board() {
        /*
         Description:
         ------------
         Prints the grid with the current values of each cell.
         */
        for (int i = 0; i < BOARD_SIZE; i++) {
            if (i % N == 0 && i != 0) {
                cout << "---------------------" << endl;
            }
            for (int j = 0; j < BOARD_SIZE; j++) {
                if (j % N == 0 && j != 0) {
                    cout << "| ";
                }
                cout << board[i][j] << " ";
            }
            cout << endl;
        }
    }
    
    void PreFill(){}
    
    vector<int> GetSquare(vector<int> cell){
        /*
         Description:
         ------------
         It returns the box that the given cell belongs to.
         */
        int i = cell[0] / N;
        int j = cell[1] / N;
        assert(i >=0 && i <= N);
        assert(j >=0 && j <= N);
        return {i, j};
    }
    
    int GetRandomElement(unordered_set<int> set){
        /*
         Description:
         ------------
         Given an unordered set as input, it selects and returns one of its
         elements randomly.
         */
        if (set.size() == 1){
            for(auto it = set.begin(); it != set.end();){
                return *it;
            }
        }
        random_device rd;
        mt19937 gen(rd());
        uniform_int_distribution<> dis(0, set.size() - 1);
        auto it = set.begin();
        advance(it, dis(gen));
        int randomElement = *it;
        return randomElement;
    }

    void InitialiseSolution(){
        /*
         Description:
         ------------
         It generates a random solution. It places random numbers in the non-fixed cells of each
         box, so that all the numbers are contained exactly once in each box.
         */

        int rowLow = 0, rowHigh = 0;
        int colLow = 0, colHigh = 0;

        // Iteration over all squares
        for (int i = 0; i < N; ++i){
            rowLow = i * N;
            rowHigh = (i + 1) * N - 1;
            for (int j = 0; j < N; ++j){
                colLow = j * N;
                colHigh = (j + 1) * N - 1;
                // Iteration over all cells inside the square
                // Find all the allowed values we can assign in the box.
                unordered_set<int> AllowedValues = AllPossibleValues;
                for (int row = rowLow; row < rowHigh + 1; ++row){
                    for (int col = colLow; col < colHigh + 1; ++col){
                        AllowedValues.erase(board[row][col]);
                    }
                }

                // Iteration over all cells inside the square.
                // Put random values in the non-fixed cells
                // so that all of them are unique.
                for (int row = rowLow; row < rowHigh + 1; ++row){
                    for (int col = colLow; col < colHigh + 1; ++col){
                        //cout << row << " - " << col << endl;
                        // If the cell is fixed, go to the next.
                        vector<int> pos = {row, col};
                        if (find(FixedPositions.begin(), FixedPositions.end(), pos) != FixedPositions.end()){
                            //cout << "Cannot change entry of pre-filled cell (" << row << "," << col << "). It remains: " << board[row][col] << endl;
                            continue;
                        }

                        // Otherwise, randomly assign a number to it so that it's unique in the box.
                        if (AllowedValues.size() > 0){
                            int randomInteger = GetRandomElement(AllowedValues);
                            AllowedValues.erase(randomInteger);
                            board[row][col] = randomInteger;
                            //cout << "Setting cell (" << row << "," << col << ") = " << randomInteger << endl;
                        }
                    }
                }
            }
        }
        print_board();
    }
    
    vector<int> GenRanPosInGrid(){
        /*
         Description:
         -----------
         Generates random position from all elements
         */
        
        // Pick a random position, not one of the initial ones, and put a random integer in it.
        vector<int> randomPos = {};
        do {
            randomPos = {(int) round(dsfmt_genrand() * BOARD_SIZE), (int) round(dsfmt_genrand() * BOARD_SIZE)};
        } while(find(FixedPositions.begin(), FixedPositions.end(), randomPos) != FixedPositions.end());
        
        assert (randomPos[0] >= 0 and randomPos[0] <= 8);
        assert (randomPos[1] >= 0 and randomPos[1] <= 8);
        
        return randomPos;
    }
    
    vector<int> GenRanPosInSquare(vector<int> square){
        
        int rowLow = square[0] * N, rowHigh = (square[0] + 1) * N - 1;
        int colLow = square[1] * N, colHigh = (square[1] + 1) * N - 1;
        
        int row = rowLow + round(dsfmt_genrand() * (rowHigh - rowLow));
        int col = colLow + round(dsfmt_genrand() * (colHigh - colLow));
        vector<int> pos = {row, col};
        
        return pos;
    }

    int GenRandomEntry(){
        int randomEntry = 0;
        /*
         Description:
         ------------
         Generates random entry, i.e. a number from 1 to BOARD_SIZE.
         */
        randomEntry = 1 + round(dsfmt_genrand() * (BOARD_SIZE - 1));
        assert (randomEntry >= 1 and randomEntry <= BOARD_SIZE);
        return randomEntry;
    }

//    int getCostFunction(){
//
//        int duplicates  = 0;
//        int value = 0;
//        // Find duplicates in all rows.
//        for (int row = 0; row < N; ++row){
//            unordered_set<int> values;
//            for (int col = 0; col < N; ++col){
//                value = board[row][col];
//                if (value != 0 and !values.insert(value).second){
//                    duplicates += 1;
//                }
//            }
//        }
//
//        // Find duplicates in all columns
//        for (int col = 0; col < N; ++col){
//            unordered_set<int> values;
//            for (int row = 0; row < N; ++row){
//                value = board[row][col];
//                if (value != 0 and !values.insert(value).second){
//                    duplicates += 1;
//                }
//            }
//        }
//        return duplicates;
//    }
//
//    void AcceptOrReject(){
//
//        int oldEntry = 0;
//        int randomX = randomPos[0];
//        int randomY = randomPos[1];
//        oldEntry = board[randomX][randomY];
//        board[randomX][randomY] = randomEntry;
//
//        int newCostFunction = getCostFunction();
//        // Accept or reject the new solution
//        int dCF = newCostFunction - costFunction;
//        //cout << "Cost function difference: " << dCF << endl;
//        if (dCF < 0 or dsfmt_genrand() < exp(- BETA * dCF)){
//            //cout << "Accepted!" << endl;
//            costFunction += dCF;
//        }else{
//            //cout << "Rejected!" << endl;
//            board[randomX][randomY] = oldEntry;
//        }
//    }
//
//    void cooldown(){
//        BETA *= (1 - COOLING_RATE);
//    }
//
//    bool is_solved(){
//
//        // Check if the cost function is zero
//        if (costFunction != 0){
//            return false;
//        }
//
//        // Check if there's any zero entry
//        for (int row = 0; row < N; ++row){
//            for (int col = 0; col < N; ++col){
//                if (board[row][col] == 0){
//                    return false;
//                }
//            }
//        }
//        int sum = 0;
//        // Check if the sum in all rows is equal to 45
//        for (int row = 0; row < N; ++row){
//            sum = 0;
//            for (int col = 0; col < N; ++col){
//                sum += board[row][col];
//            }
//            if (sum != 45){
//                return false;
//            }
//        }
//        // Check if the sum in all columns is equal to 45
//        for (int col = 0; col < N; ++col){
//            sum = 0;
//            for (int row = 0; row < N; ++row){
//                sum += board[row][col];
//            }
//            if (sum != 45){
//                return false;
//            }
//        }
//        return true;
//    }
//
//    void solve(){
////        int i = 0;
//        do{
//            GenRandomPosition();
//            GenRandomEntry();
//            AcceptOrReject();
//            cooldown();
////            i += 1;
//        }while (not is_solved());
//
//        cout << "Cost function: " << costFunction <<  endl;
//        print_board();
//    }
    
    
};

//vector<list<int>> init(){
//    vector<int[2]> temp = {};
//    temp.push_back([1,2]);
//    return temp;
//}

int main(int argc, const char * argv[]) {
    
    dsfmt_seed( time (NULL)); // Initialise random seed.

    SudokuSolver sudoku;
    //sudoku.print_board();
    //sudoku.InitialiseSolution();
    //sudoku.GenRandomPosition();
    
    unordered_set<int> vals = {1,2,3,4,5,6,7,8,9};
    int entry = 0;
    do{
        entry = sudoku.GenRandomEntry();
        cout << entry << endl;
        vals.erase(entry);
    }while(vals.size() > 0);
    
    return 0;
}
