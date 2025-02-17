#include "move_generator.h"
#include "board_utils.h"
#include <algorithm>
#include <cstdlib>
#include <utility>
#include <vector>
#include <iostream>
using namespace std;

//------------------------------------------------------------------------------
// Aggregate score evaluation function.
// Computes sub-scores from connected groups on the board using DFS,
// then combines them with weighted constants.
double evaluateState(const BoardState& state, 
                     const vector<vector<bool>>& isGoal,
                     int minScore, int rows, int cols) {
    // Sub-score variables
    int numSets = 0;               // Count of scoring sets (group on goal with size >= minScore)
    int numScorableTiles = 0;      // Total number of tiles in scoring sets
    double scoringSetScore = 0;    // Sum of squares for scoring sets
    double nearScoringSetScore = 0;// Sum of squares for near-scoring sets (>= 3 on a goal)
    double nonScoringSetScore = 0; // Sum of squares for all other groups
    double adjacencyScore = 0;     // (Not computed in this example)
    double positioningScore = 0;   // (Not computed in this example)

    // Track visited cells so that each connected group is processed once.
    vector<vector<bool>> visited(rows, vector<bool>(cols, false));

    // Loop over every cell.
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (visited[i][j])
                continue;
            char ch = state.board[i][j];
            if (ch == '.' || ch == '*')
                continue;
            
            // DFS to gather the connected group (same color)
            vector<pair<int,int>> group;
            vector<pair<int,int>> stack;
            stack.push_back({i, j});
            visited[i][j] = true;
            while (!stack.empty()) {
                auto [r, c] = stack.back();
                stack.pop_back();
                group.push_back({r, c});
                const int dr[4] = {1, -1, 0, 0};
                const int dc[4] = {0, 0, 1, -1};
                for (int k = 0; k < 4; k++) {
                    int nr = r + dr[k], nc = c + dc[k];
                    if (nr < 0 || nr >= rows || nc < 0 || nc >= cols)
                        continue;
                    if (visited[nr][nc])
                        continue;
                    char neighbor = state.board[nr][nc];
                    if (neighbor == '.' || neighbor == '*')
                        continue;
                    if (pieceColor(neighbor) == pieceColor(ch)) {
                        visited[nr][nc] = true;
                        stack.push_back({nr, nc});
                    }
                }
            }
            
            int groupSize = group.size();
            bool hasGoal = false;
            for (auto& cell : group) {
                int r = cell.first, c = cell.second;
                if (isGoal[r][c]) {
                    hasGoal = true;
                    break;
                }
            }
            
            // Classify the group into one of three buckets.
            if (hasGoal) {
                if (groupSize >= minScore) {
                    numSets++;
                    numScorableTiles += groupSize;
                    scoringSetScore += groupSize * groupSize;
                } else if (groupSize >= 3) {
                    nearScoringSetScore += groupSize * groupSize;
                } else {
                    nonScoringSetScore += groupSize * groupSize;
                }
            } else {
                nonScoringSetScore += groupSize * groupSize;
            }
        }
    }
    
    // Combine sub-scores with empirical weightings.
    double finalScore = (numSets + numScorableTiles / 3.0) * 8000.0
                        + (scoringSetScore * 1.3 + nearScoringSetScore * 1.3 + nonScoringSetScore * 1.1) * 100.0
                        + (adjacencyScore * 20.0 + positioningScore) * 10.0;
                        
    return finalScore;
}

//------------------------------------------------------------------------------
// Add random cells to fill the board (or parts of it) after a move.
BoardState addRandomCells(const BoardState& state, int numCells, 
                           const string& colors, int rows, int cols) {
    BoardState newState = state;
    vector<pair<int, int>> emptyCells;
    
    // Locate empty cells.
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            char ch = newState.board[i][j];
            if (ch == '.' || ch == '*')
                emptyCells.push_back({i, j});
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

//------------------------------------------------------------------------------
// Simulate swapping two cells.
BoardState simulateSwap(const BoardState& state, int r1, int c1, int r2, int c2) {
    BoardState newState = state;
    swap(newState.board[r1][c1], newState.board[r2][c2]);
    return newState;
}

//------------------------------------------------------------------------------
// Simulate a scoring move: remove a connected group from the board.
BoardState simulateScoreMove(const BoardState& state, int startR, int startC, 
                             int rows, int cols,
                             const vector<vector<bool>>& isGoal) {
    BoardState newState = state;
    vector<vector<bool>> visited(rows, vector<bool>(cols, false));
    char color = pieceColor(state.board[startR][startC]);

    // DFS to remove the connected group.
    vector<pair<int, int>> stack;
    stack.push_back({startR, startC});
    visited[startR][startC] = true;

    while (!stack.empty()) {
        auto [r, c] = stack.back();
        stack.pop_back();

        newState.board[r][c] = (isGoal[r][c] ? '*' : '.');

        const int dr[] = {1, -1, 0, 0};
        const int dc[] = {0, 0, 1, -1};
        for (int k = 0; k < 4; k++) {
            int nr = r + dr[k], nc = c + dc[k];
            if (nr < 0 || nr >= rows || nc < 0 || nc >= cols) continue;
            if (visited[nr][nc]) continue;
            char ch = state.board[nr][nc];
            if (ch == '.' || ch == '*') continue;
            if (pieceColor(ch) == color) {
                visited[nr][nc] = true;
                stack.push_back({nr, nc});
            }
        }
    }
    return newState;
}

//------------------------------------------------------------------------------
// Move selection function using the final aggregate score.
// It generates candidate moves (both score and swap moves), evaluates each
// outcome using evaluateState(), and returns the move with the best score.
Move chooseMoveExpectimax(const BoardState& state, 
                          const vector<vector<bool>>& isGoal,
                          int minScore, int rows, int cols, 
                          const string& colors, int samples) {
    constexpr int MIN_GROUP_SIZE_FOR_IMMEDIATE = 5;  // Immediately take large groups
    constexpr double PRUNING_THRESHOLD = -1000.0;    // Prune clearly bad moves
    
    struct MoveCandidate {
        Move move;
        double value;
        BoardState nextState;
    };
    vector<MoveCandidate> candidates;
    
    // --- 1. Consider score moves ---
    vector<GroupInfo> groups = analyzeBoard(state.board, isGoal, minScore, rows, cols);
    for (const auto& g : groups) {
        // If a very large scoring group is found, return immediately.
        if (g.size >= MIN_GROUP_SIZE_FOR_IMMEDIATE) {
            return Move{SCORE, g.goalR, g.goalC, 0, 0};
        }
        
        BoardState nextState = simulateScoreMove(state, g.goalR, g.goalC, rows, cols, isGoal);
        double quickVal = evaluateState(nextState, isGoal, minScore, rows, cols);
        
        if (quickVal > PRUNING_THRESHOLD) {
            candidates.push_back({
                Move{SCORE, g.goalR, g.goalC, 0, 0},
                quickVal,
                nextState
            });
        }
    }
    
    // --- 2. Consider swap moves ---
    vector<pair<int, int>> validPieces;
    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            if (state.board[r][c] != '.' && state.board[r][c] != '*')
                validPieces.emplace_back(r, c);
        }
    }
    
    for (size_t i = 0; i < validPieces.size(); i++) {
        for (size_t j = i + 1; j < validPieces.size(); j++) {
            int r1 = validPieces[i].first, c1 = validPieces[i].second;
            int r2 = validPieces[j].first, c2 = validPieces[j].second;
            
            char cell1 = state.board[r1][c1];
            char cell2 = state.board[r2][c2];
            if (pieceColor(cell1) == pieceColor(cell2))
                continue;
            
            BoardState swapState = simulateSwap(state, r1, c1, r2, c2);
            double quickVal = evaluateState(swapState, isGoal, minScore, rows, cols);
            
            if (quickVal > PRUNING_THRESHOLD) {
                candidates.push_back({
                    Move{SWAP, r1, c1, r2, c2},
                    quickVal,
                    swapState
                });
            }
        }
    }
    
    // Fallback: If no candidates, try any valid swap.
    if (candidates.empty()) {
        for (size_t i = 0; i < validPieces.size(); i++) {
            for (size_t j = i + 1; j < validPieces.size(); j++) {
                int r1 = validPieces[i].first, c1 = validPieces[i].second;
                int r2 = validPieces[j].first, c2 = validPieces[j].second;
                return Move{SWAP, r1, c1, r2, c2};
            }
        }
        return Move{SWAP, 0, 0, 0, 0}; // Last resort
    }
    
    // Sort candidates by their quick evaluation scores.
    sort(candidates.begin(), candidates.end(),
         [](const MoveCandidate& a, const MoveCandidate& b) {
             return a.value > b.value;
         });
         
    // Retain only a limited number of top candidates for further evaluation.
    constexpr int MAX_DETAILED_CANDIDATES = 10;
    if (candidates.size() > MAX_DETAILED_CANDIDATES)
        candidates.resize(MAX_DETAILED_CANDIDATES);
    
    // --- 3. Detailed evaluation using random sampling ---
    for (auto& candidate : candidates) {
        double totalVal = 0;
        int actualSamples = (candidate.move.type == SCORE) ? samples : samples * 2;
        
        for (int s = 0; s < actualSamples; s++) {
            int addCells = (candidate.move.type == SCORE) ? 3 : 5;
            BoardState outcome = addRandomCells(candidate.nextState, addCells, colors, rows, cols);
            double eval = evaluateState(outcome, isGoal, minScore, rows, cols);
            totalVal += eval * (s < samples ? 1.0 : 0.5);  // Weight early samples more
        }
        candidate.value = totalVal / (samples * 1.5);  // Normalize the evaluation
    }
    
    // --- 4. Choose the best move based on the detailed evaluation ---
    auto bestCandidate = max_element(candidates.begin(), candidates.end(),
                                   [](const MoveCandidate& a, const MoveCandidate& b) {
                                       return a.value < b.value;
                                   });
                                   
    return bestCandidate->move;
}
