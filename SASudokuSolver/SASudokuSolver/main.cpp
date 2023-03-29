

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
    
    double Temperature {100};
    const double COOLING_RATE {0.99};
    const int MARKOV_CHAIN_REPS = 10;
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
    
    int CostFunction = 0;
    
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

    void InitSolution(){
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
            randomPos = {(int) dsfmt_genrand() * BOARD_SIZE, (int) dsfmt_genrand() * BOARD_SIZE};
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
    
    
    int GetCostContribution(int row, int col){
        /*
         Description:
         ------------
         Calculates the contribution of row and col to the cost function.
         That is, it counts the number of integers that are not present
         in row and col.
         */
        int CostContribution = 0;
        
        unordered_set<int> NotPresentInts = AllPossibleValues;
        for (int i = 0; i < BOARD_SIZE; ++i){
            NotPresentInts.erase(board[row][i]);
        }
        CostContribution += NotPresentInts.size();
        
        NotPresentInts = AllPossibleValues;
        for (int i = 0; i < BOARD_SIZE; ++i){
            NotPresentInts.erase(board[i][col]);
        }
        CostContribution += NotPresentInts.size();
        
        return CostContribution;
    }
    
    
    int GetCostFunction(){
        /*
         Desctiption:
         ------------
         Calculates the cost function, or objective function to optimise.
         It is simply the total number of integers that are not present in
         each row and columns.
         */

        int CostFunction = 0;
        for (int i = 0; i < BOARD_SIZE; ++i){
            CostFunction += GetCostContribution(i, i);
        }
        return CostFunction;
    }
    
    
    void Swap(vector<int> cell1, vector<int> cell2){
        /*
         Description:
         ------------
         Swaps the values of two cells.
         */
        
        int row1 = cell1[0], col1 = cell1[1];
        int row2 = cell2[0], col2 = cell2[1];
        int temp = board[row1][col1];
        board[row1][col1] = board[row2][col2];
        board[row2][col2] = temp;
        
    }
        
    void AttemptSwap(vector<int> cell1, vector<int> cell2){
        /*
         Description:
         ------------
         Attempts a swap between the values of twe cells. If the cost function
         is lower the swap is accepted. If the swap leads to higher cost fucntion
         it is accepted with probability exp(-Δ/T).
         */
        
        int row1 = cell1[0], col1 = cell1[1];
        int row2 = cell2[0], col2 = cell2[1];
        
        // Calculate starting cost function contribution.
        int CostBefore = GetCostContribution(row1, col1) + GetCostContribution(row2, col2);
        //cout << "Cost before: " << CostBefore << endl;
        
        // Make the swap.
        Swap(cell1, cell2);
        
        // Calculate new cost function contribution and change (delta).
        int CostAfter = GetCostContribution(row1, col1) + GetCostContribution(row2, col2);
        int DeltaCost = CostAfter - CostBefore;
        
        //cout << "Cost after: " << CostAfter << endl;
        //cout << "Delta: " << DeltaCost << endl;
        
        // Accept or reject.
        double randNum = dsfmt_genrand();
        //cout << "Random num: " << randNum << endl;
        //cout << "Threshold: " << exp(-DeltaCost/Temperature) << endl;
        if (DeltaCost > 0 and randNum > exp(-DeltaCost/Temperature)){
            Swap(cell1, cell2);
        }else{
            CostFunction += DeltaCost;
        }
    }
    
    
    void Cooldown(){
        Temperature *= COOLING_RATE;
    }
    
    
    

    void Solve(){

        // 1. Initialise random solution.
        InitSolution();
        CostFunction = GetCostFunction();

        while (CostFunction != 0 && Temperature > 0){
            
            for (int i = 0; i < 14; ++i){
                // 2. Choose random non-fixed cell in grid.
                vector<int> cell1 = {};
                do{
                    cell1 = GenRanPosInGrid();
                }while(find(FixedPositions.begin(), FixedPositions.end(), cell1) != FixedPositions.end());

                // 3. Choose random non-fixed cell in the same square
                vector<int> square = GetSquare(cell1);
                vector<int> cell2 = {};
                do{
                    cell2 = GenRanPosInSquare(square);
                }while(find(FixedPositions.begin(), FixedPositions.end(), cell2) != FixedPositions.end());
                
                AttemptSwap(cell1, cell2);
            }
            cout << "Cost Function: " << CostFunction << endl;
            Cooldown();
        }
    }
};



int main(int argc, const char * argv[]) {
    
    dsfmt_seed( time (NULL)); // Initialise random seed.

    SudokuSolver sudoku;
    //sudoku.print_board();
    sudoku.InitSolution();
    
    sudoku.Solve();

    
    return 0;
}
