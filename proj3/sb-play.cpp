#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cctype>
#include <sstream>
#include <algorithm>
#include <limits>
#include <ctime>
#include "disjoint.h" // Uses your provided disjoint set implementation
using namespace std;

// Return the canonical (lowercase) color.
char pieceColor(char ch) {
    return tolower(ch);
}

// Structure to hold information about a connected component.
struct GroupInfo {
    int root;
    int size;
    bool hasGoal;
    int goalR, goalC; // a representative goal cell (if any)
};

enum MoveType { SWAP, SCORE };

struct Move {
    MoveType type;
    // For SWAP: use (r1,c1) and (r2,c2)
    // For SCORE: only r1 and c1 are used.
    int r1, c1, r2, c2;
};

// Board state structure.
struct BoardState {
    vector<string> board;
};


// This is essentially the same analysis as before but we now also only report the groups at their root.
vector<GroupInfo> analyzeBoard(const vector<string>& board, const vector<vector<bool>>& isGoal,
                               int minScore, int rows, int cols) {
    // Create disjoint set (one element per cell).
    DisjointSetByRankWPC ds(rows * cols);

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            char ch = board[i][j];
            if (ch == '.' || ch == '*') continue; //Corrected *
            int index = i * cols + j;

            // Only look right and down to avoid duplicate unions.
            if (j + 1 < cols) {
                char nbr = board[i][j + 1];
                if (nbr != '.' && nbr != '*' && pieceColor(nbr) == pieceColor(ch)) { //Corrected *
                    int r1 = ds.Find(index);
                    int r2 = ds.Find(i * cols + (j + 1));
                    if (r1 != r2)
                        ds.Union(r1, r2);
                }
            }
            if (i + 1 < rows) {
                char nbr = board[i + 1][j];
                if (nbr != '.' && nbr != '*' && pieceColor(nbr) == pieceColor(ch)) { //Corrected *
                    int r1 = ds.Find(index);
                    int r2 = ds.Find((i + 1) * cols + j);
                    if (r1 != r2)
                        ds.Union(r1, r2);
                }
            }
        }
    }

    vector<int> compSize(rows * cols, 0);
    vector<bool> compHasGoal(rows * cols, false);
    vector<int> compRep(rows * cols, -1);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            char ch = board[i][j];
            if (ch == '.' || ch == '*') continue; //Corrected *
            int index = i * cols + j;
            int root = ds.Find(index);
            compSize[root]++;
            if (isGoal[i][j]) {
                compHasGoal[root] = true;
                compRep[root] = index;
            }
        }
    }

    vector<GroupInfo> groups;
    for (int i = 0; i < rows * cols; i++) {
        if (compSize[i] >= minScore && compHasGoal[i] && ds.Find(i) == i) { // only consider roots
            GroupInfo info;
            info.root = i;
            info.size = compSize[i];
            info.hasGoal = true;
            if (compRep[i] != -1) {
                info.goalR = compRep[i] / cols;
                info.goalC = compRep[i] % cols;
            } else {
                info.goalR = 0;
                info.goalC = 0;
            }
            groups.push_back(info);
        }
    }
    return groups;
}

// A simple heuristic: sum the sizes of scoreable groups.
double evaluateState(const BoardState & state, const vector<vector<bool>>& isGoal,
                     int minScore, int rows, int cols) {
    vector<GroupInfo> groups = analyzeBoard(state.board, isGoal, minScore, rows, cols);
    double val = 0;
    
    // Score existing groups based on their color values
    for (auto &g : groups) {
        int r = g.goalR, c = g.goalC;
        char color = pieceColor(state.board[r][c]);
        int pointValue = 0;
        switch(color) {
            case 'p': pointValue = 2; break;
            case 'b': pointValue = 3; break;
            case 'y': pointValue = 4; break;
            case 'r': pointValue = 5; break;
            case 'g': pointValue = 6; break;
        }
        // Value = size * point value * (1 + bonus for very large groups)
        val += g.size * pointValue * (1.0 + max(0, g.size - minScore) * 0.1);
    }
    
    // Look for potential future groups
    vector<vector<bool>> visited(rows, vector<bool>(cols, false));
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (visited[i][j] || state.board[i][j] == '.' || state.board[i][j] == '*') 
                continue;
                
            // Count connected pieces
            char color = pieceColor(state.board[i][j]);
            int size = 0;
            bool nearGoal = false;
            
            // DFS to count group size and check proximity to goals
            vector<pair<int,int>> stack({{i,j}});
            visited[i][j] = true;
            
            while (!stack.empty()) {
                auto [r, c] = stack.back();
                stack.pop_back();
                size++;
                
                // Check if this piece is adjacent to a goal
                for (int dr = -1; dr <= 1; dr++) {
                    for (int dc = -1; dc <= 1; dc++) {
                        int nr = r + dr, nc = c + dc;
                        if (nr >= 0 && nr < rows && nc >= 0 && nc < cols) {
                            if (isGoal[nr][nc]) nearGoal = true;
                        }
                    }
                }
                
                // Add neighbors of same color
                int dr[] = {1, -1, 0, 0};
                int dc[] = {0, 0, 1, -1};
                for (int k = 0; k < 4; k++) {
                    int nr = r + dr[k], nc = c + dc[k];
                    if (nr >= 0 && nr < rows && nc >= 0 && nc < cols &&
                        !visited[nr][nc] && 
                        state.board[nr][nc] != '.' && state.board[nr][nc] != '*' &&
                        pieceColor(state.board[nr][nc]) == color) {
                        visited[nr][nc] = true;
                        stack.push_back({nr, nc});
                    }
                }
            }
            
            // Add value for potential future groups
            if (size >= (minScore - 1) && nearGoal) {
                int pointValue = 0;
                switch(color) {
                    case 'p': pointValue = 2; break;
                    case 'b': pointValue = 3; break;
                    case 'y': pointValue = 4; break;
                    case 'r': pointValue = 5; break;
                    case 'g': pointValue = 6; break;
                }
                val += (size * pointValue * 0.3); // Discount future potential
            }
        }
    }
    
    return val;
}

// After a move, new random pieces are added to empty cells.
// For SWAP moves, 5 new pieces are added; for SCORE moves, 3 are added.
BoardState addRandomCells(const BoardState & state, int numCells, const string & colors, int rows, int cols) {
    BoardState newState = state;
    vector<pair<int, int>> emptyCells;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            char ch = newState.board[i][j];
            if (ch == '.' || ch == '*')
                emptyCells.push_back({ i, j });
        }
    }
    int count = min((int)emptyCells.size(), numCells);
    random_shuffle(emptyCells.begin(), emptyCells.end());
    for (int k = 0; k < count; k++) {
        int r = emptyCells[k].first;
        int c = emptyCells[k].second;
        int idx = rand() % colors.size();
        newState.board[r][c] = colors[idx];
    }
    return newState;
}

// Simulate a SWAP move.
BoardState simulateSwap(const BoardState & state, int r1, int c1, int r2, int c2) {
    BoardState newState = state;
    swap(newState.board[r1][c1], newState.board[r2][c2]);
    return newState;
}

// Simulate a SCORE move. Starting from cell (startR, startC) (which must be in a scoreable group),
// remove all connected pieces (of the same color) by setting them to '.' (or '*' if on a goal).  Corrected.
BoardState simulateScoreMove(const BoardState & state, int startR, int startC, int rows, int cols,
                              const vector<vector<bool>>& isGoal) {
    BoardState newState = state;
    vector<vector<bool>> visited(rows, vector<bool>(cols, false));
    char color = pieceColor(state.board[startR][startC]);

    vector<pair<int, int>> stack;
    stack.push_back({ startR, startC });
    visited[startR][startC] = true;

    while (!stack.empty()) {
        auto [r, c] = stack.back();
        stack.pop_back();

        // Remove the piece.
        newState.board[r][c] = (isGoal[r][c] ? '*' : '.'); //Corrected this

        int dr[4] = { 1, -1, 0, 0 };
        int dc[4] = { 0, 0, 1, -1 };
        for (int k = 0; k < 4; k++) {
            int nr = r + dr[k], nc = c + dc[k];
            if (nr < 0 || nr >= rows || nc < 0 || nc >= cols) continue;
            if (visited[nr][nc]) continue;
            char ch = state.board[nr][nc];
            if (ch == '.' || ch == '*') continue;
            if (pieceColor(ch) == color) {
                visited[nr][nc] = true;
                stack.push_back({ nr, nc });
            }
        }
    }
    return newState;
}

// This expectimax decision function (one-ply lookahead with sampling) examines all legal moves,
// simulates the random new cell additions, and chooses the move with highest average evaluation.
Move chooseMoveExpectimax(const BoardState& state, const vector<vector<bool>>& isGoal,
                         int minScore, int rows, int cols, const string& colors, int samples = 10) {
    // Use more samples by default for better accuracy
    constexpr int MIN_GROUP_SIZE_FOR_IMMEDIATE = 5;  // Immediately take large groups
    constexpr double PRUNING_THRESHOLD = -1000.0;    // Prune clearly bad moves
    
    struct MoveCandidate {
        Move move;
        double value;
        BoardState nextState;
    };
    vector<MoveCandidate> candidates;
    
    // First pass: Generate and evaluate all possible moves
    // 1. Score moves
    vector<GroupInfo> groups = analyzeBoard(state.board, isGoal, minScore, rows, cols);
    for (const auto& g : groups) {
        // Immediately return very large scoring groups
        if (g.size >= MIN_GROUP_SIZE_FOR_IMMEDIATE) {
            return Move{SCORE, g.goalR, g.goalC, 0, 0};
        }
        
        BoardState nextState = simulateScoreMove(state, g.goalR, g.goalC, rows, cols, isGoal);
        // Quick evaluation without sampling
        double quickVal = evaluateState(nextState, isGoal, minScore, rows, cols);
        
        if (quickVal > PRUNING_THRESHOLD) {
            candidates.push_back({
                Move{SCORE, g.goalR, g.goalC, 0, 0},
                quickVal,
                nextState
            });
        }
    }
    
    // 2. Swap moves
    vector<pair<int, int>> validPieces;
    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            if (state.board[r][c] != '.' && state.board[r][c] != '*') {
                validPieces.emplace_back(r, c);
            }
        }
    }
    
    // Generate swap moves in parallel
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < validPieces.size(); i++) {
        for (size_t j = i + 1; j < validPieces.size(); j++) {
            int r1 = validPieces[i].first, c1 = validPieces[i].second;
            int r2 = validPieces[j].first, c2 = validPieces[j].second;
            
            char cell1 = state.board[r1][c1];
            char cell2 = state.board[r2][c2];
            
            if (pieceColor(cell1) == pieceColor(cell2)) continue;
            
            BoardState swapState = simulateSwap(state, r1, c1, r2, c2);
            double quickVal = evaluateState(swapState, isGoal, minScore, rows, cols);
            
            if (quickVal > PRUNING_THRESHOLD) {
                #pragma omp critical
                candidates.push_back({
                    Move{SWAP, r1, c1, r2, c2},
                    quickVal,
                    swapState
                });
            }
        }
    }
    
    // If no candidates found, return default move
    if (candidates.empty()) {
        return Move{SWAP, 0, 0, 0, 0};
    }
    
    // Sort candidates by quick evaluation
    sort(candidates.begin(), candidates.end(),
         [](const MoveCandidate& a, const MoveCandidate& b) {
             return a.value > b.value;
         });
         
    // Keep only top N candidates for detailed evaluation
    constexpr int MAX_DETAILED_CANDIDATES = 10;
    if (candidates.size() > MAX_DETAILED_CANDIDATES) {
        candidates.resize(MAX_DETAILED_CANDIDATES);
    }
    
    // Second pass: Detailed evaluation with sampling
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < candidates.size(); i++) {
        double totalVal = 0;
        int actualSamples = candidates[i].move.type == SCORE ? samples : samples * 2;
        
        for (int s = 0; s < actualSamples; s++) {
            int addCells = candidates[i].move.type == SCORE ? 3 : 5;
            BoardState outcome = addRandomCells(candidates[i].nextState, addCells, colors, rows, cols);
            
            // Use weighted evaluation
            double eval = evaluateState(outcome, isGoal, minScore, rows, cols);
            totalVal += eval * (s < samples ? 1.0 : 0.5);  // Weight early samples more
        }
        
        candidates[i].value = totalVal / (samples * 1.5);  // Adjust for weighted samples
    }
    
    // Find best move after detailed evaluation
    auto bestCandidate = max_element(candidates.begin(), candidates.end(),
                                   [](const MoveCandidate& a, const MoveCandidate& b) {
                                       return a.value < b.value;
                                   });
                                   
    return bestCandidate->move;
}


int main(int argc, char** argv) {
    srand(time(NULL));
    if (argc < 5) {
        cerr << "Usage: " << argv[0] << " rows cols min_score_size colors\n";
        return 1;
    }
    int rows = atoi(argv[1]);
    int cols = atoi(argv[2]);
    int minScore = atoi(argv[3]);
    string colors = argv[4];

    // Read board (one line per row).
    vector<string> board;
    string line;
    for (int i = 0; i < rows; i++) {
        getline(cin, line);
        while ((int)line.size() < cols)
            line.push_back(' ');
        board.push_back(line);
    }

    // Build isGoal matrix (goal cells are '*' or uppercase letters).
    vector<vector<bool>> isGoal(rows, vector<bool>(cols, false));
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            char ch = board[i][j];
            if (ch == '*' || (isalpha(ch) && isupper(ch)))
                isGoal[i][j] = true;
        }
    }

    BoardState state;
    state.board = board;
    Move chosen = chooseMoveExpectimax(state, isGoal, minScore, rows, cols, colors, 5);
    if (chosen.type == SCORE)
        cout << "SCORE " << chosen.r1 << " " << chosen.c1 << "\n";
    else if (chosen.type == SWAP)
        cout << "SWAP " << chosen.r1 << " " << chosen.c1 << " " << chosen.r2 << " " << chosen.c2 << "\n";
    return 0;
}