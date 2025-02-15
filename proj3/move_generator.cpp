#include "move_generator.h"
#include "board_utils.h"
#include <algorithm>
#include <cstdlib>
#include <utility>
#include <vector>

using namespace std;

BoardState addRandomCells(const BoardState& state, int numCells, 
                         const string& colors, int rows, int cols) {
    BoardState newState = state;
    vector<pair<int, int>> emptyCells;
    
    // Find all empty cells
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            char ch = newState.board[i][j];
            if (ch == '.' || ch == '*')
                emptyCells.push_back({i, j});
        }
    }
    
    // Add random pieces to empty cells
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

BoardState simulateSwap(const BoardState& state, int r1, int c1, int r2, int c2) {
    BoardState newState = state;
    swap(newState.board[r1][c1], newState.board[r2][c2]);
    return newState;
}

BoardState simulateScoreMove(const BoardState& state, int startR, int startC, 
                            int rows, int cols,
                            const vector<vector<bool>>& isGoal) {
    BoardState newState = state;
    vector<vector<bool>> visited(rows, vector<bool>(cols, false));
    char color = pieceColor(state.board[startR][startC]);

    // Use DFS to find and remove connected pieces
    vector<pair<int, int>> stack;
    stack.push_back({startR, startC});
    visited[startR][startC] = true;

    while (!stack.empty()) {
        auto [r, c] = stack.back();
        stack.pop_back();

        // Remove the piece
        newState.board[r][c] = (isGoal[r][c] ? '*' : '.');

        // Check neighbors
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
    
    // Generate swap moves
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
        // Try to find any valid swap as a fallback
        for (size_t i = 0; i < validPieces.size(); i++) {
            for (size_t j = i + 1; j < validPieces.size(); j++) {
                int r1 = validPieces[i].first, c1 = validPieces[i].second;
                int r2 = validPieces[j].first, c2 = validPieces[j].second;
                return Move{SWAP, r1, c1, r2, c2};
            }
        }
        return Move{SWAP, 0, 0, 0, 0}; // Last resort
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
    for (auto& candidate : candidates) {
        double totalVal = 0;
        int actualSamples = candidate.move.type == SCORE ? samples : samples * 2;
        
        for (int s = 0; s < actualSamples; s++) {
            int addCells = candidate.move.type == SCORE ? 3 : 5;
            BoardState outcome = addRandomCells(candidate.nextState, addCells, colors, rows, cols);
            
            // Use weighted evaluation
            double eval = evaluateState(outcome, isGoal, minScore, rows, cols);
            totalVal += eval * (s < samples ? 1.0 : 0.5);  // Weight early samples more
        }
        
        candidate.value = totalVal / (samples * 1.5);  // Adjust for weighted samples
    }
    
    // Find best move after detailed evaluation
    auto bestCandidate = max_element(candidates.begin(), candidates.end(),
                                   [](const MoveCandidate& a, const MoveCandidate& b) {
                                       return a.value < b.value;
                                   });
                                   
    return bestCandidate->move;
}