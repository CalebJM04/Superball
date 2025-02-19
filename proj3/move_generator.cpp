#include "board_utils.h"
#include "move_generator.h"
#include "disjoint.h"
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <utility>
#include <vector>
#include <unordered_map>
using namespace std;

// Helper functions for board evaluation

bool isCorner(int row, int col, int rows, int cols) {
  return (row == 2 &&
          (col == 0 || col == 1 || col == cols - 2 || col == cols - 1)) ||
         (row == 5 &&
          (col == 0 || col == 1 || col == cols - 2 || col == cols - 1));
}

bool isEdge(int row, int col, int rows, int cols) {
  return (row >= 2 && row <= 5 &&
          (col == 0 || col == 1 || col == cols - 2 || col == cols - 1));
}

int getEmptyCellCount(const BoardState &state, int rows, int cols) {
  int count = 0;
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      if (state.board[i][j] == '.' || state.board[i][j] == '*') {
        count++;
      }
    }
  }
  return count;
}

bool hasAdjacentDifferentColorGroup(const BoardState &state, int row, int col,
                                    int rows, int cols) {
  const int dr[] = {1, -1, 0, 0};
  const int dc[] = {0, 0, 1, -1};
  char currentColor = pieceColor(state.board[row][col]);

  for (int k = 0; k < 4; k++) {
    int nr = row + dr[k], nc = col + dc[k];
    if (nr < 0 || nr >= rows || nc < 0 || nc >= cols)
      continue;
    char neighborColor = pieceColor(state.board[nr][nc]);
    if (neighborColor != '.' && neighborColor != '*' &&
        neighborColor != currentColor) {
      return true;
    }
  }
  return false;
}

double computeConnectivityScore(const BoardState &state,
                                const vector<pair<int, int>> &group, int rows,
                                int cols) {
  double score = 0;
  for (const auto &cell : group) {
    if (hasAdjacentDifferentColorGroup(state, cell.first, cell.second, rows,
                                       cols)) {
      score += 2.9; // 2.9 worked best, 1763 avg
    }
  }
  return score;
}

// Utility function to convert 2D coordinates to a 1D index
int cellToIndex(int row, int col, int cols) {
  return row * cols + col;
}

// Utility function to convert 1D index to 2D coordinates
pair<int, int> indexToCell(int index, int cols) {
  return {index / cols, index % cols};
}

// Using disjoint sets to identify connected groups
vector<vector<pair<int, int>>> findConnectedGroups(const BoardState &state, int rows, int cols) {
  // Total number of cells in the board
  int totalCells = rows * cols;
  
  // Create disjoint set data structure
  DisjointSetByRankWPC dset(totalCells);
  
  // Map to store group representatives and their corresponding cells
  unordered_map<int, vector<pair<int, int>>> groupMap;
  
  // First pass: Union adjacent cells of the same color
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      char currentPiece = state.board[i][j];
      if (currentPiece == '.' || currentPiece == '*') continue;
      
      int currentIdx = cellToIndex(i, j, cols);
      char currentColor = pieceColor(currentPiece);
      
      // Check all four neighbors
      const int dr[] = {0, 1, 0, -1};
      const int dc[] = {1, 0, -1, 0};
      
      for (int k = 0; k < 4; k++) {
        int ni = i + dr[k];
        int nj = j + dc[k];
        
        if (ni < 0 || ni >= rows || nj < 0 || nj >= cols) continue;
        
        char neighborPiece = state.board[ni][nj];
        if (neighborPiece == '.' || neighborPiece == '*') continue;
        
        char neighborColor = pieceColor(neighborPiece);
        if (neighborColor == currentColor) {
          int neighborIdx = cellToIndex(ni, nj, cols);
          dset.Union(dset.Find(currentIdx), dset.Find(neighborIdx));
        }
      }
    }
  }
  
  // Second pass: Collect cells belonging to each group
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      char currentPiece = state.board[i][j];
      if (currentPiece == '.' || currentPiece == '*') continue;
      
      int currentIdx = cellToIndex(i, j, cols);
      int groupRep = dset.Find(currentIdx);
      
      groupMap[groupRep].push_back({i, j});
    }
  }
  
  // Convert the map to a vector of groups
  vector<vector<pair<int, int>>> groups;
  for (const auto &entry : groupMap) {
    groups.push_back(entry.second);
  }
  
  return groups;
}

double computeFinalScore(const BoardState &state,
                         const vector<vector<bool>> &isGoal, int minScore,
                         int rows, int cols) {
  int numSets = 0;
  int numScorableTiles = 0;
  double scoringSetScore = 0;
  double nearScoringSetScore = 0;
  double nonScoringSetScore = 0;
  double positioningScore = 0;
  double connectivityScore = 0;

  // Find all connected groups using disjoint sets
  vector<vector<pair<int, int>>> groups = findConnectedGroups(state, rows, cols);

  for (const auto &group : groups) {
    int groupSize = group.size();
    bool hasGoal = false;
    char groupColor = pieceColor(state.board[group[0].first][group[0].second]);

    // Calculate position-based scoring for the group
    double groupPositionScore = 0;
    for (auto &cell : group) {
      int r = cell.first, c = cell.second;
      if (isGoal[r][c])
        hasGoal = true;
      if (isCorner(r, c, rows, cols))
        groupPositionScore += 5.0; // 5.0 - avg 1782
      else if (isEdge(r, c, rows, cols))
        groupPositionScore += 9.5; // 9.5 - 1932 avg
    }

    // Calculate connectivity score for the group
    double groupConnectivityScore =
        computeConnectivityScore(state, group, rows, cols);

    // Classify and score the group
    if (hasGoal) {
      if (groupSize >= minScore) {
        numSets++;
        numScorableTiles += groupSize;
        scoringSetScore += groupSize * groupSize * pow(1.1, max(0, groupSize - minScore)); //                  scoringSetScore += groupSize * groupSize;
        positioningScore += groupPositionScore * 1.8; // Higher weight for scoring groups - 1932 avg
        connectivityScore += groupConnectivityScore * 1.3; // 1.3 - 1966 avg
      } else if (groupSize >= 3) {
        nearScoringSetScore += groupSize * groupSize;
        positioningScore += groupPositionScore * 1.2;
        connectivityScore += groupConnectivityScore * 1.2;
      } else {
        nonScoringSetScore += groupSize * groupSize;
        positioningScore += groupPositionScore;
        connectivityScore += groupConnectivityScore;
      }
    } else {
      nonScoringSetScore += groupSize * groupSize;
      positioningScore += groupPositionScore * 0.8;
      connectivityScore += groupConnectivityScore * 0.8;
    }
  }

  // Combine sub-scores with empirical weightings
  double finalScore = (numSets + numScorableTiles / 3.0) * 8000.0 +
                      (scoringSetScore * 1.3 + nearScoringSetScore * 1.3 +
                       nonScoringSetScore * 1.1) *
                          100.0 +
                      positioningScore * 20.0 + connectivityScore * 10.0;

  return finalScore;
}

// Simulate swapping two cells.
BoardState simulateSwap(const BoardState &state, int r1, int c1, int r2,
                        int c2) {
  BoardState newState = state;
  swap(newState.board[r1][c1], newState.board[r2][c2]);
  return newState;
}

// Simulate a scoring move using disjoint sets
BoardState simulateScoreMove(const BoardState &state, int startR, int startC,
                           int rows, int cols,
                           const vector<vector<bool>> &isGoal) {
  BoardState newState = state;
  char targetColor = pieceColor(state.board[startR][startC]);
  
  // Create disjoint set
  DisjointSetByRankWPC dset(rows * cols);
  
  // First pass: build the disjoint set of connected components
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      if (pieceColor(state.board[i][j]) != targetColor) continue;
      
      int currentIdx = cellToIndex(i, j, cols);
      
      // Check all four neighbors
      const int dr[] = {0, 1, 0, -1};
      const int dc[] = {1, 0, -1, 0};
      
      for (int k = 0; k < 4; k++) {
        int ni = i + dr[k];
        int nj = j + dc[k];
        
        if (ni < 0 || ni >= rows || nj < 0 || nj >= cols) continue;
        
        if (pieceColor(state.board[ni][nj]) == targetColor) {
          int neighborIdx = cellToIndex(ni, nj, cols);
          dset.Union(dset.Find(currentIdx), dset.Find(neighborIdx));
        }
      }
    }
  }
  
  // Find the representative of the group containing the startR, startC cell
  int startIdx = cellToIndex(startR, startC, cols);
  int groupRep = dset.Find(startIdx);
  
  // Remove all cells in the same group
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      if (pieceColor(state.board[i][j]) != targetColor) continue;
      
      int currentIdx = cellToIndex(i, j, cols);
      if (dset.Find(currentIdx) == groupRep) {
        newState.board[i][j] = (isGoal[i][j] ? '*' : '.');
      }
    }
  }
  
  return newState;
}

// Update the analyzeBoard function to use disjoint sets
vector<GroupInfo> analyzeBoard(const vector<vector<char>> &board,
                              const vector<vector<bool>> &isGoal,
                              int minScore, int rows, int cols) {
  vector<GroupInfo> result;
  
  // Create disjoint set
  DisjointSetByRankWPC dset(rows * cols);
  unordered_map<int, vector<pair<int, int>>> groupMap;
  unordered_map<int, char> groupColors;
  
  // First pass: Union adjacent cells of the same color
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      char currentPiece = board[i][j];
      if (currentPiece == '.' || currentPiece == '*') continue;
      
      int currentIdx = cellToIndex(i, j, cols);
      char currentColor = pieceColor(currentPiece);
      
      // Check all four neighbors
      const int dr[] = {0, 1, 0, -1};
      const int dc[] = {1, 0, -1, 0};
      
      for (int k = 0; k < 4; k++) {
        int ni = i + dr[k];
        int nj = j + dc[k];
        
        if (ni < 0 || ni >= rows || nj < 0 || nj >= cols) continue;
        
        char neighborPiece = board[ni][nj];
        if (neighborPiece == '.' || neighborPiece == '*') continue;
        
        char neighborColor = pieceColor(neighborPiece);
        if (neighborColor == currentColor) {
          int neighborIdx = cellToIndex(ni, nj, cols);
          dset.Union(dset.Find(currentIdx), dset.Find(neighborIdx));
        }
      }
    }
  }
  
  // Second pass: Collect cells belonging to each group and record group colors
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      char currentPiece = board[i][j];
      if (currentPiece == '.' || currentPiece == '*') continue;
      
      int currentIdx = cellToIndex(i, j, cols);
      int groupRep = dset.Find(currentIdx);
      
      groupMap[groupRep].push_back({i, j});
      groupColors[groupRep] = pieceColor(currentPiece);
    }
  }
  
  // Third pass: Analyze each group
  for (const auto &entry : groupMap) {
    const vector<pair<int, int>> &group = entry.second;
    int groupSize = group.size();
    int groupRep = entry.first;
    char groupColor = groupColors[groupRep];
    
    // Check if group contains a goal cell
    bool hasGoal = false;
    int goalR = -1, goalC = -1;
    
    for (const auto &cell : group) {
      int r = cell.first, c = cell.second;
      if (isGoal[r][c]) {
        hasGoal = true;
        goalR = r;
        goalC = c;
        break;
      }
    }
    
    // Only include groups that have a goal cell and meet the minimum score requirement
    if (hasGoal && groupSize >= minScore) {
      GroupInfo info;
      // Note: Assuming GroupInfo only contains size and position, not color
      info.size = groupSize;
      info.goalR = goalR;
      info.goalC = goalC;
      result.push_back(info);
    }
  }
  
  return result;
}

// Simplified move selection without future simulation
Move chooseMove(const BoardState &state, const vector<vector<bool>> &isGoal, 
                int minScore, int rows, int cols) {
  // Dynamic threshold based on game phase
  int emptyCells = getEmptyCellCount(state, rows, cols);
  int MIN_GROUP_SIZE_FOR_IMMEDIATE = 5;
  
  struct MoveCandidate {
    Move move;
    double value;
  };
  vector<MoveCandidate> candidates;
  
  // First, consider all scoring moves
  vector<GroupInfo> groups = analyzeBoard(state.board, isGoal, minScore, rows, cols);
  for (const auto &g : groups) {
    if (g.size >= MIN_GROUP_SIZE_FOR_IMMEDIATE) {
      // Take immediate large scores without further evaluation
      return Move{SCORE, g.goalR, g.goalC, 0, 0};
    }
    
    // Evaluate the move by seeing how good the board would be after scoring
    BoardState nextState = simulateScoreMove(state, g.goalR, g.goalC, rows, cols, isGoal);
    double moveValue = computeFinalScore(nextState, isGoal, minScore, rows, cols);
    candidates.push_back({Move{SCORE, g.goalR, g.goalC, 0, 0}, moveValue});
  }
  
  // Consider swap moves (prioritizing corner/edge pieces)
  vector<pair<int, int>> validPieces;
  for (int r = 0; r < rows; r++) {
    for (int c = 0; c < cols; c++) {
      if (state.board[r][c] != '.' && state.board[r][c] != '*')
        validPieces.emplace_back(r, c);
    }
  }
  
  // Prioritize corner and edge pieces
  sort(validPieces.begin(), validPieces.end(),
       [rows, cols](const pair<int, int> &a, const pair<int, int> &b) {
         bool aCorner = isCorner(a.first, a.second, rows, cols);
         bool bCorner = isCorner(b.first, b.second, rows, cols);
         if (aCorner != bCorner) return aCorner > bCorner;
         
         bool aEdge = isEdge(a.first, a.second, rows, cols);
         bool bEdge = isEdge(b.first, b.second, rows, cols);
         return aEdge > bEdge;
       });
  
  // Evaluate promising swaps (limit to reasonable number to avoid excessive computation)
  const int MAX_SWAP_CANDIDATES = 500000; //50k high score
  int swapsChecked = 0;
  
  for (size_t i = 0; i < validPieces.size() && swapsChecked < MAX_SWAP_CANDIDATES; i++) {
    for (size_t j = i + 1; j < validPieces.size() && swapsChecked < MAX_SWAP_CANDIDATES; j++) {
      int r1 = validPieces[i].first, c1 = validPieces[i].second;
      int r2 = validPieces[j].first, c2 = validPieces[j].second;
      
      char cell1 = state.board[r1][c1];
      char cell2 = state.board[r2][c2];
      if (pieceColor(cell1) == pieceColor(cell2))
        continue;
      
      swapsChecked++;
      BoardState swapState = simulateSwap(state, r1, c1, r2, c2);
      double moveValue = computeFinalScore(swapState, isGoal, minScore, rows, cols);
      
      candidates.push_back({Move{SWAP, r1, c1, r2, c2}, moveValue});
    }
  }
  
  // Fallback handling for empty board or no good moves
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
  
  // Find best move
  auto bestCandidate = max_element(candidates.begin(), candidates.end(),
                                  [](const MoveCandidate &a, const MoveCandidate &b) {
                                    return a.value < b.value;
                                  });
  
  return bestCandidate->move;
}

// Main function that interfaces with the game system
Move chooseBestMove(const BoardState &state, const vector<vector<bool>> &isGoal,
                    int minScore, int rows, int cols, const string &colors) {
  return chooseMove(state, isGoal, minScore, rows, cols);
}