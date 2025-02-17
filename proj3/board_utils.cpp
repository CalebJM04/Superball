#include "board_utils.h"
#include "disjoint.h"
#include <cctype>
#include <algorithm>
#include <vector>
#include <utility>
#include <cmath>

using namespace std;

char pieceColor(char ch) {
    return tolower(ch);
}

static int colorToPoints(char color) {
    switch (tolower(color)) {
        case 'p': return 2;
        case 'b': return 3;
        case 'y': return 4;
        case 'r': return 5;
        case 'g': return 6;
        default:  return 1;
    }
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
                     int minScore, int rows, int cols)
{
    // --- 1) Find actual scoreable groups via union-find (existing code) ---
    vector<GroupInfo> groups = analyzeBoard(state.board, isGoal, minScore, rows, cols);

    double val = 0.0;
    // Award points for groups that already meet minScore and are on a goal
    for (const auto& g : groups) {
        int r = g.goalR, c = g.goalC;
        char color = tolower(state.board[r][c]);
        int pointValue = colorToPoints(color);

        // Existing logic: bigger groups get increasing bonus
        // e.g. size * pointValue * (1 + (size - minScore)*0.1)
        val += g.size * pointValue * (1.0 + max(0, g.size - minScore) * 0.1);
    }

    // --- 2) Compute "future potential" from standard DFS approach (existing) ---
    vector<vector<bool>> visited(rows, vector<bool>(cols, false));
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            char ch = state.board[i][j];
            if (ch == '.' || ch == '*' || visited[i][j]) continue;

            // DFS to find how many connected pieces of the same color,
            // and whether any piece is near a goal.
            char color = tolower(ch);
            int pointValue = colorToPoints(color);

            int size = 0;
            bool nearGoal = false;

            vector<pair<int,int>> stack{{i,j}};
            visited[i][j] = true;

            while (!stack.empty()) {
                auto [r, c] = stack.back();
                stack.pop_back();
                size++;

                // Check if near a goal
                for (int dr = -1; dr <= 1; dr++) {
                    for (int dc = -1; dc <= 1; dc++) {
                        int nr = r + dr, nc = c + dc;
                        if (nr >= 0 && nr < rows && nc >= 0 && nc < cols) {
                            if (isGoal[nr][nc]) {
                                nearGoal = true;
                            }
                        }
                    }
                }

                // Add neighbors of the same color
                static const int DR[4] = {1, -1, 0, 0};
                static const int DC[4] = {0, 0, 1, -1};
                for (int k = 0; k < 4; k++) {
                    int nr = r + DR[k], nc = c + DC[k];
                    if (nr >= 0 && nr < rows && nc >= 0 && nc < cols) {
                        if (!visited[nr][nc] && tolower(state.board[nr][nc]) == color &&
                            state.board[nr][nc] != '.' && state.board[nr][nc] != '*') {
                            visited[nr][nc] = true;
                            stack.push_back({nr, nc});
                        }
                    }
                }
            }

            // If this group is just under minScore but near a goal,
            // give a partial bonus to encourage building it up.
            if (size >= (minScore - 1) && nearGoal) {
                val += (size * pointValue * 0.3);
            }
        }
    }

    // --- 3) New: "Combo potential" check for near-mergeable groups ---
    // For each color, see if we can find two separate groups close enough
    // that a single swap would merge them into a bigger group.
    // We'll do a quick pass checking adjacency of boundaries of different groups
    // of the same color. This is just an example heuristic.
    {
        // 3a) Build a union-find or BFS again to identify each distinct group
        //     (for the same color).
        // We can reuse a modified version of 'analyzeBoard' or just do a quick pass:
        struct CellInfo {
            int compId;
            char color;
        };
        vector<vector<CellInfo>> comp(rows, vector<CellInfo>(cols, {-1, '.'}));
        int compCounter = 0;

        // Mark connected components of same color
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                if (state.board[i][j] == '.' || state.board[i][j] == '*') {
                    comp[i][j].compId = -1;
                } else {
                    comp[i][j].color = tolower(state.board[i][j]);
                }
            }
        }

        // BFS/DFS to label each group with compId
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                if (comp[i][j].compId == -1 && comp[i][j].color != '.' && comp[i][j].color != '*') {
                    // new group
                    char ccolor = comp[i][j].color;
                    vector<pair<int,int>> stack{{i,j}};
                    comp[i][j].compId = compCounter;

                    while (!stack.empty()) {
                        auto [r, c] = stack.back();
                        stack.pop_back();
                        static const int DR[4] = {1, -1, 0, 0};
                        static const int DC[4] = {0, 0, 1, -1};
                        for (int k = 0; k < 4; k++) {
                            int nr = r + DR[k], nc = c + DC[k];
                            if (nr >= 0 && nr < rows && nc >= 0 && nc < cols) {
                                if (comp[nr][nc].color == ccolor && comp[nr][nc].compId == -1) {
                                    comp[nr][nc].compId = compCounter;
                                    stack.push_back({nr, nc});
                                }
                            }
                        }
                    }
                    compCounter++;
                }
            }
        }

        // 3b) For each pair of distinct components of the same color, check if
        // they are near each other (within 1 or 2 cells). If so, add a small
        // bonus that encourages merging them.
        vector<int> compSizes(compCounter, 0);
        vector<char> compColors(compCounter, '.');

        // Count sizes and track color of each component
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                int cid = comp[i][j].compId;
                if (cid >= 0) {
                    compSizes[cid]++;
                    compColors[cid] = comp[i][j].color;
                }
            }
        }

        // Check adjacency between components of the same color
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                int cid1 = comp[i][j].compId;
                if (cid1 < 0) continue;
                // For each neighbor in a small radius
                for (int dr = -2; dr <= 2; dr++) {
                    for (int dc = -2; dc <= 2; dc++) {
                        if (dr == 0 && dc == 0) continue;
                        int nr = i + dr, nc = j + dc;
                        if (nr < 0 || nr >= rows || nc < 0 || nc >= cols) continue;
                        int cid2 = comp[nr][nc].compId;
                        if (cid2 < 0 || cid2 == cid1) continue;
                        // If same color, and close enough, add a bonus
                        if (compColors[cid1] == compColors[cid2]) {
                            int totalSize = compSizes[cid1] + compSizes[cid2];
                            // If merging them would exceed minScore (or become significantly bigger):
                            if (totalSize >= minScore) {
                                int colorPoints = colorToPoints(compColors[cid1]);
                                // Example formula: give a partial bonus for the possibility of merging
                                // The closer they are, the bigger the bonus. You can refine the distance measure
                                double dist = sqrt(dr*dr + dc*dc);
                                double proximityFactor = (dist <= 1.0 ? 1.0 : 0.5);

                                // Weighted by how big the combined group could be
                                double comboBonus = totalSize * colorPoints * 0.05 * proximityFactor;
                                val += comboBonus;
                            }
                        }
                    }
                }
            }
        }
    }

    return val;
}