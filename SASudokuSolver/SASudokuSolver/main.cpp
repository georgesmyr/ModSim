

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <random>

#include "mt19937.h"

using namespace std;

string path = "/Users/georgesmyridis/Documents/Physics/Books-Notes/Graduate/Physics/Modeling_Simulations/Scripts/ModSim/SASudokuSolver/performance/";
const int BOARD_SIZE = 9;
const int N = sqrt(BOARD_SIZE);

double get_std(vector<int> set){
    /*
     Description:
     -----------
     Calculates the standard diviatioin for a set of integers.
     */
    double mean = 0, std = 0, sum = 0;
    
    for (int i = 0; i < set.size(); ++i){
        sum += set[i];
    }
    mean = sum / set.size();
    sum = 0;
    for (int i = 0; i < set.size(); ++i){
        sum += pow((mean - set[i]), 2);
    }
    std = sum / set.size();
    return std;
}

class SudokuSolver{
    
private:
    
    double Temperature {20};
    double InitTemperature = Temperature;
    double COOLING_RATE {0.99};
    int MARKOV_CHAIN_REPS = 10;
    vector<vector<int>> board = {
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
    vector<vector<int>> FixedPositions = GetFixedPositions();

    
public:
    
    void SetBoard(vector<vector<int>> boardd){
        board = boardd;
    }
    
    void Cooldown(){
        Temperature *= COOLING_RATE;
    }
        
    void print_board() {
        /*
         Description:
         ------------
         Prints the grid with the current values of each cell.
         */
        cout << "---------------------" << endl;
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
        cout << "---------------------" << endl;
    }
    
    void PreFill(string level){
        
        if (level == "easy"){
            board = {
                {0, 7, 0, 0, 0, 0, 6, 0, 0},
                {9, 0, 0, 0, 0, 3, 0, 0, 1},
                {0, 6, 0, 0, 0, 8, 7, 3, 0},
                {0, 0, 0, 0, 9, 0, 0, 0, 0},
                {0, 8, 0, 0, 0, 0, 0, 4, 0},
                {0, 0, 0, 0, 5, 0, 0, 0, 0},
                {0, 1, 5, 4, 0, 0, 0, 9, 0},
                {4, 0, 0, 0, 0, 0, 0, 0, 7},
                {0, 0, 3, 0, 0, 0, 0, 0, 0}
            };
            
        }else if (level == "medium"){
            board = {
                {2, 0, 0, 7, 0, 8, 0, 0, 5},
                {0, 8, 4, 0, 0, 0, 1, 7, 0},
                {0, 0, 0, 0, 0, 0, 0, 0, 0},
                {0, 6, 0, 0, 0, 0, 0, 1, 0},
                {0, 0, 0, 0, 0, 0, 0, 0, 0},
                {0, 2, 0, 0, 0, 0, 0, 4, 0},
                {0, 0, 0, 0, 0, 0, 0, 0, 0},
                {0, 5, 7, 0, 0, 0, 3, 8, 0},
                {3, 0, 0, 1, 0, 9, 0, 0, 6}
            };
            
        }else if (level == "hard"){
            board = {
                {0, 0, 0, 0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 3, 0, 8, 5},
                {0, 0, 1, 0, 2, 0, 0, 0, 0},
                {0, 0, 0, 5, 0, 7, 0, 0, 0},
                {0, 0, 4, 0, 0, 0, 1, 0, 0},
                {0, 9, 0, 0, 0, 0, 0, 0, 0},
                {5, 0, 0, 0, 0, 0, 0, 7, 3},
                {0, 0, 2, 0, 1, 0, 0, 0, 0},
                {0, 0, 0, 0, 4, 0, 0, 0, 9}
            };
        }else{
            return;
        }
            FixedPositions = GetFixedPositions();
    }
    
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

    vector<vector<int>> InitialSolution = {};
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

                        // If the cell is fixed, go to the next.
                        vector<int> pos = {row, col};
                        if (find(FixedPositions.begin(), FixedPositions.end(), pos) != FixedPositions.end()){
                            continue;
                        }

                        // Otherwise, randomly assign a number to it so that it's unique in the box.
                        if (AllowedValues.size() > 0){
                            int randomInteger = GetRandomElement(AllowedValues);
                            AllowedValues.erase(randomInteger);
                            board[row][col] = randomInteger;
                        }
                    }
                }
            }
        }
        InitialSolution = board;
        //print_board();
    }
    
    vector<int> GenRanPosInGrid(){
        /*
         Description:
         -----------
         Generates random position from all elements
         */
        
        // Pick a random position, not one of the initial ones, and put a random integer in it.
        vector<int> randomPos = {};
        int row = 0, col = 0;
        do {
            row = int(dsfmt_genrand() * BOARD_SIZE);
            col = int(dsfmt_genrand() * BOARD_SIZE);
            randomPos = {row, col};
        } while(find(FixedPositions.begin(), FixedPositions.end(), randomPos) != FixedPositions.end());
        
        assert (randomPos[0] >= 0 and randomPos[0] <= 8);
        assert (randomPos[1] >= 0 and randomPos[1] <= 8);
        
        return randomPos;
    }
    
    vector<int> GenRanPosInSquare(vector<int> square){
        /*
         Description:
         ------------
         Generates random position of cell inside a specific square.
         */
        
        int rowLow = square[0] * N, rowHigh = (square[0] + 1) * N - 1;
        int colLow = square[1] * N, colHigh = (square[1] + 1) * N - 1;
        int row = -1, col = -1;
        vector<int> pos = {row, col};
        do{
            row = rowLow + round(dsfmt_genrand() * (rowHigh - rowLow));
            col = colLow + round(dsfmt_genrand() * (colHigh - colLow));
            pos = {row, col};
        }while(find(FixedPositions.begin(), FixedPositions.end(), pos) != FixedPositions.end());
        
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
    
    int CostFunction = 0;
    vector<int> CostFunctions = {};
    
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
        
        // Make the swap.
        Swap(cell1, cell2);
        
        // Calculate new cost function contribution and change (delta).
        int CostAfter = GetCostContribution(row1, col1) + GetCostContribution(row2, col2);
        int DeltaCost = CostAfter - CostBefore;
        
        // Accept or reject.
        double randNum = dsfmt_genrand();
        if (DeltaCost > 0 and randNum > exp(-DeltaCost/Temperature)){
            Swap(cell1, cell2);
        }else{
            CostFunction += DeltaCost;
        }
    }
    
    
    void SetControlParameters(double InitTemp, double CoolingRate, int MarkovChainReps){
        /*
         Description:
         ------------
         Sets the control parameters.
         */
        Temperature = InitTemp;
        InitTemperature = InitTemp;
        COOLING_RATE = CoolingRate;
        MARKOV_CHAIN_REPS = MarkovChainReps;
    }
    
    
    void SetIdealControlParameters(){
        /*
         Description:
         ------------
         Sets ideal ini parameters, i.e. Temperature so that approximately 80% of the swaps are accepted,
         cooling rate 0.9, and Markov chain reps so that there's a goood chance of every non-fixed pair
         be selected at least once in each chain.
         */
        cout << "Setting ideal control parameters.." << endl;
        
        // Find ideal temperature.
        vector<int> cost_contributions = {};
        const int TRIALS = 100;
        vector<int> cell1 = {}, square = {}, cell2 = {};
        int temp_entry = 0;
        for (int trial = 0; trial < TRIALS; ++trial){
            cell1 = GenRanPosInGrid();
            square = GetSquare(cell1);
            do{
                cell2 = GenRanPosInSquare(square);
            }while (cell1 == cell2);
        
            temp_entry = board[cell1[0]][cell1[1]];
            board[cell1[0]][cell1[1]] = board[cell2[0]][cell2[1]];
            board[cell2[0]][cell2[1]] = temp_entry;
            
            cost_contributions.push_back(GetCostFunction());
        }
        double std = get_std(cost_contributions);
        Temperature = std;
        cout << std << endl;
        // Revert back to the initial random solution
        board = InitialSolution;
        
        // Set ideal Markov Chain Reps
        int NumNonFixed = pow(BOARD_SIZE, 2) - FixedPositions.size(); // Number of non-fixed cells.
        MARKOV_CHAIN_REPS = pow(NumNonFixed, 2);
        cout << "Markov chain length: " << MARKOV_CHAIN_REPS << endl;
        
        // Set ideal cooling rate
        COOLING_RATE = 0.99;
        
    }
    
    
    vector<double> PerformanceState = {};


    void Solve(){
        int steps = 0;
        // 1. Initialise random solution.
        SetBoard(InitialSolution);
        CostFunction = GetCostFunction();

        while (CostFunction != 0 && Temperature > 0 && steps < 100000){
            // Markov chain iteration
            steps += 1;
            // 2. Choose random non-fixed cell in grid.
            vector<int> cell1 = GenRanPosInGrid();
            
            // 3. Choose random non-fixed cell in the same square
            vector<int> square = GetSquare(cell1);
            vector<int> cell2 = GenRanPosInSquare(square);
             
            // 4. Apply neighbourhood operator and accept or reject.
            AttemptSwap(cell1, cell2);
            CostFunctions.push_back(CostFunction);
            //cout << "Cost Function: " << CostFunction << endl;
            if (steps % MARKOV_CHAIN_REPS == 0){
                Cooldown();
            }
        }
        cout << "Temperature: " << InitTemperature << "|Cooldown Rate: " << COOLING_RATE <<
                "|Markov Length: " << MARKOV_CHAIN_REPS << "|Steps: " << steps << "|Cost: " << CostFunction << endl;
        PerformanceState = {InitTemperature, COOLING_RATE, (double) MARKOV_CHAIN_REPS, (double) steps, (double) CostFunction};
    }
    
};



int main(int argc, const char * argv[]) {
    
    dsfmt_seed( time (NULL)); // Initialise random seed.
    
    SudokuSolver s;
    s.PreFill("medium");
    s.print_board();
    s.InitSolution();
    s.print_board();
    //s.SetIdealControlParameters();
    s.SetControlParameters(50, 0.99, 10);
    s.Solve();

    ofstream outfile(path + "energies_medium.csv");
    if (!outfile.is_open()) {
        cerr << "Error: Unable to open output file" << endl;
        return  0;
    }
    outfile << "Energy" << endl;
    for (auto& cost : s.CostFunctions){
        outfile << cost << endl;
    }


//    vector<double> temperatures = {1, 2.5, 5, 10, 25};
//    vector<double> cooling_rates = {0.25, 0.5, 0.75, 0.9, 0.99};
//    vector<double> mc_lengths = {10, 25, 50, 100};
//
//    for (int i = 0; i < 10; ++i){
//        ofstream outfile(path + "performance" + to_string(i) + ".csv");
//        if (!outfile.is_open()) {
//            cerr << "Error: Unable to open output file" << endl;
//            return  0;
//        }
//        outfile << "InitTemp,CoolingRate,MarkovLength,Steps,CostFunction" << endl;
//
//        SudokuSolver sudoku;
//        sudoku.InitSolution();
//        for (auto& temp: temperatures){
//            for (auto& cr : cooling_rates){
//                for (auto& ml: mc_lengths){
//                    sudoku.SetControlParameters(temp, cr, ml);
//                    sudoku.Solve();
//                    vector<double> perf = sudoku.PerformanceState;
//                    outfile << perf[0] << "," << perf[1] << "," << perf[2] << "," << perf[3] << "," << perf[4] << endl;
//                }
//            }
//        }
//    }
    

    return 0;
}
