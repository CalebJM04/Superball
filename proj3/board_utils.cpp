#include "board_utils.h"
#include "disjoint.h"
#include <cctype>
#include <algorithm>
#include <vector>
#include <utility>

using namespace std;

char pieceColor(char ch) {
    return tolower(ch);
}

vector<GroupInfo> analyzeBoard(const vector<string>& board, 
                             const vector<vector<bool>>& isGoal,
                             int minScore, int rows, int cols) {
    // Create disjoint set (one element per cell)
    DisjointSetByRankWPC ds(rows * cols);

    // Connect adjacent same-colored pieces
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            char ch = board[i][j];
            if (ch == '.' || ch == '*') continue;
            int index = i * cols + j;

            // Only look right and down to avoid duplicate unions
            if (j + 1 < cols) {
                char nbr = board[i][j + 1];
                if (nbr != '.' && nbr != '*' && pieceColor(nbr) == pieceColor(ch)) {
                    int r1 = ds.Find(index);
                    int r2 = ds.Find(i * cols + (j + 1));
                    if (r1 != r2)
                        ds.Union(r1, r2);
                }
            }
            if (i + 1 < rows) {
                char nbr = board[i + 1][j];
                if (nbr != '.' && nbr != '*' && pieceColor(nbr) == pieceColor(ch)) {
                    int r1 = ds.Find(index);
                    int r2 = ds.Find((i + 1) * cols + j);
                    if (r1 != r2)
                        ds.Union(r1, r2);
                }
            }
        }
    }

    // Analyze connected components
    vector<int> compSize(rows * cols, 0);
    vector<bool> compHasGoal(rows * cols, false);
    vector<int> compRep(rows * cols, -1);
    
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            char ch = board[i][j];
            if (ch == '.' || ch == '*') continue;
            int index = i * cols + j;
            int root = ds.Find(index);
            compSize[root]++;
            if (isGoal[i][j]) {
                compHasGoal[root] = true;
                compRep[root] = index;
            }
        }
    }

    // Build group info for qualifying components
    vector<GroupInfo> groups;
    for (int i = 0; i < rows * cols; i++) {
        if (compSize[i] >= minScore && compHasGoal[i] && ds.Find(i) == i) {
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

double evaluateState(const BoardState& state, 
                    const vector<vector<bool>>& isGoal,
                    int minScore, int rows, int cols) {
    vector<GroupInfo> groups = analyzeBoard(state.board, isGoal, minScore, rows, cols);
    double val = 0.0;
    
    // Score existing groups based on their color values
    for (const auto& g : groups) {
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
                const int dr[] = {1, -1, 0, 0};
                const int dc[] = {0, 0, 1, -1};
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