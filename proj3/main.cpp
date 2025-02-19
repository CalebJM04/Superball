#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cctype>
#include "types.h"
#include "board_utils.h"
#include "move_generator.h"

using namespace std;

int main(int argc, char** argv) {
    // Initialize random seed
    // srand(time(NULL));

    // Check command line arguments
    if (argc < 5) {
        cerr << "Usage: " << argv[0] << " rows cols min_score_size colors\n";
        cerr << "Example: " << argv[0] << " 8 8 3 pbyRg\n";
        cerr << "Colors: p=purple, b=blue, y=yellow, r=red, g=green\n";
        cerr << "Board format:\n";
        cerr << "  . = empty cell\n";
        cerr << "  * = empty goal cell\n";
        cerr << "  lowercase letter = regular piece\n";
        cerr << "  uppercase letter = piece on goal\n";
        return 1;
    }

    // Parse command line arguments
    int rows = atoi(argv[1]);
    int cols = atoi(argv[2]);
    int minScore = atoi(argv[3]);
    string colors = argv[4];

    // Validate arguments
    if (rows <= 0 || cols <= 0 || minScore <= 0) {
        cerr << "Error: Invalid dimensions or minimum score size\n";
        return 1;
    }

    if (colors.empty()) {
        cerr << "Error: Must specify at least one color\n";
        return 1;
    }

    int seed = (argc >= 6) ? atoi(argv[5]) : time(NULL);
    srand(seed);

    // Read board (one line per row)
    vector<string> board;
    string line;
    for (int i = 0; i < rows; i++) {
        if (!getline(cin, line)) {
            cerr << "Error: Failed to read row " << i << " of the board\n";
            return 1;
        }
        // Pad short lines with spaces
        while ((int)line.size() < cols) {
            line.push_back(' ');
        }
        // Trim long lines
        if ((int)line.size() > cols) {
            line = line.substr(0, cols);
        }
        board.push_back(line);
    }

    // Build isGoal matrix (goal cells are '*' or uppercase letters)
    vector<vector<bool>> isGoal(rows, vector<bool>(cols, false));
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            char ch = board[i][j];
            if (ch == '*' || (isalpha(ch) && isupper(ch))) {
                isGoal[i][j] = true;
            }
        }
    }

    // Create initial board state
    BoardState state;
    state.board = board;

    // Choose and output the best move
    // Choose and output the best move
Move chosen = chooseBestMove(state, isGoal, minScore, rows, cols, colors); // 500 was best
    
    // Output the chosen move
    if (chosen.type == SCORE) {
        cout << "SCORE " << chosen.r1 << " " << chosen.c1 << "\n";
    } else if (chosen.type == SWAP) {
        cout << "SWAP " << chosen.r1 << " " << chosen.c1 << " " 
             << chosen.r2 << " " << chosen.c2 << "\n";
    }

    return 0;
}