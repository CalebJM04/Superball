#include "board_utils.h"
#include "move_generator.h"
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <utility>
#include <vector>
using namespace std;

//------------------------------------------------------------------------------
// Aggregate score evaluation function renamed to computeFinalScore.
// Computes sub-scores from connected groups on the board using DFS,
// then combines them with weighted constants.

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
      score += 0.5;
    }
  }
  return score;
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

  vector<vector<bool>> visited(rows, vector<bool>(cols, false));

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      if (visited[i][j])
        continue;
      char ch = state.board[i][j];
      if (ch == '.' || ch == '*')
        continue;

      vector<pair<int, int>> group;
      vector<pair<int, int>> stack;
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

      // Calculate position-based scoring for the group
      double groupPositionScore = 0;
      for (auto &cell : group) {
        int r = cell.first, c = cell.second;
        if (isGoal[r][c])
          hasGoal = true;
        if (isCorner(r, c, rows, cols))
          groupPositionScore += 2.0;
        else if (isEdge(r, c, rows, cols))
          groupPositionScore += 1.0;
      }

      // Calculate connectivity score for the group
      double groupConnectivityScore =
          computeConnectivityScore(state, group, rows, cols);

      // Classify and score the group
      if (hasGoal) {
        if (groupSize >= minScore) {
          numSets++;
          numScorableTiles += groupSize;
          scoringSetScore += groupSize * groupSize;
          positioningScore +=
              groupPositionScore * 1.6; // Higher weight for scoring groups // 1.5
          connectivityScore += groupConnectivityScore * 1.5;
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
        positioningScore +=
            groupPositionScore * 0.8; // Lower weight for non-goal groups
        connectivityScore += groupConnectivityScore * 0.8;
      }
    }
  }

  // Combine sub-scores with empirical weightings
  double finalScore = (numSets + numScorableTiles / 3.0) * 10000.0 +
                      (scoringSetScore * 1.5 + nearScoringSetScore * 1.4 +
                       nonScoringSetScore * 1.2) *
                          150.0 +
                      positioningScore * 75.0 + connectivityScore * 50.0;

  return finalScore;
}

//------------------------------------------------------------------------------
// Add random cells to fill the board (or parts of it) after a move.
BoardState addRandomCells(const BoardState &state, int numCells,
                          const string &colors, int rows, int cols) {
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
BoardState simulateSwap(const BoardState &state, int r1, int c1, int r2,
                        int c2) {
  BoardState newState = state;
  swap(newState.board[r1][c1], newState.board[r2][c2]);
  return newState;
}

//------------------------------------------------------------------------------
// Simulate a scoring move: remove a connected group from the board.
BoardState simulateScoreMove(const BoardState &state, int startR, int startC,
                             int rows, int cols,
                             const vector<vector<bool>> &isGoal) {
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
      if (nr < 0 || nr >= rows || nc < 0 || nc >= cols)
        continue;
      if (visited[nr][nc])
        continue;
      char ch = state.board[nr][nc];
      if (ch == '.' || ch == '*')
        continue;
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
// outcome using computeFinalScore(), and returns the move with the best score.
Move chooseMoveExpectimax(const BoardState &state,
                          const vector<vector<bool>> &isGoal, int minScore,
                          int rows, int cols, const string &colors,
                          int samples) {
  // Dynamic threshold based on game phase
  int emptyCells = getEmptyCellCount(state, rows, cols);
  int MIN_GROUP_SIZE_FOR_IMMEDIATE = (emptyCells < 20) ? 5 : // :
                                         (emptyCells < 40) ? 5
                                     : (emptyCells < 60)
                                         ? 6
                                         : 7; // Adjusted for late game
  constexpr double PRUNING_THRESHOLD = -1000.0;

  struct MoveCandidate {
    Move move;
    double value;
    BoardState nextState;
  };
  vector<MoveCandidate> candidates;

  // Score moves consideration
  vector<GroupInfo> groups =
      analyzeBoard(state.board, isGoal, minScore, rows, cols);
  for (const auto &g : groups) {
    if (g.size >= MIN_GROUP_SIZE_FOR_IMMEDIATE) {
      // Add extra verification for very early game
      if (emptyCells > 60 && g.size <= 5) {
        // Don't immediately take small groups in early game
        continue;
      }
      return Move{SCORE, g.goalR, g.goalC, 0, 0};
    }

    BoardState nextState =
        simulateScoreMove(state, g.goalR, g.goalC, rows, cols, isGoal);
    double quickVal =
        computeFinalScore(nextState, isGoal, minScore, rows, cols);

    if (quickVal > PRUNING_THRESHOLD) {
      candidates.push_back(
          {Move{SCORE, g.goalR, g.goalC, 0, 0}, quickVal, nextState});
    }
  }

  // Swap moves consideration
  vector<pair<int, int>> validPieces;
  for (int r = 0; r < rows; r++) {
    for (int c = 0; c < cols; c++) {
      if (state.board[r][c] != '.' && state.board[r][c] != '*')
        validPieces.emplace_back(r, c);
    }
  }

  // Prioritize corner and edge pieces for swaps
  sort(validPieces.begin(), validPieces.end(),
       [rows, cols](const pair<int, int> &a, const pair<int, int> &b) {
         bool aCorner = isCorner(a.first, a.second, rows, cols);
         bool bCorner = isCorner(b.first, b.second, rows, cols);
         if (aCorner != bCorner)
           return aCorner > bCorner;

         bool aEdge = isEdge(a.first, a.second, rows, cols);
         bool bEdge = isEdge(b.first, b.second, rows, cols);
         return aEdge > bEdge;
       });

  for (size_t i = 0; i < validPieces.size(); i++) {
    for (size_t j = i + 1; j < validPieces.size(); j++) {
      int r1 = validPieces[i].first, c1 = validPieces[i].second;
      int r2 = validPieces[j].first, c2 = validPieces[j].second;

      char cell1 = state.board[r1][c1];
      char cell2 = state.board[r2][c2];
      if (pieceColor(cell1) == pieceColor(cell2))
        continue;

      BoardState swapState = simulateSwap(state, r1, c1, r2, c2);
      double quickVal =
          computeFinalScore(swapState, isGoal, minScore, rows, cols);

      if (quickVal > PRUNING_THRESHOLD) {
        candidates.push_back({Move{SWAP, r1, c1, r2, c2}, quickVal, swapState});
      }
    }
  }

  // Fallback handling remains the same
  if (candidates.empty()) {
    for (size_t i = 0; i < validPieces.size(); i++) {
      for (size_t j = i + 1; j < validPieces.size(); j++) {
        int r1 = validPieces[i].first, c1 = validPieces[i].second;
        int r2 = validPieces[j].first, c2 = validPieces[j].second;
        return Move{SWAP, r1, c1, r2, c2};
      }
    }
    return Move{SWAP, 0, 0, 0, 0};
  }

  // Sort and evaluate top candidates
  sort(candidates.begin(), candidates.end(),
       [](const MoveCandidate &a, const MoveCandidate &b) {
         return a.value > b.value;
       });

  constexpr int MAX_DETAILED_CANDIDATES = 10;
  if (candidates.size() > MAX_DETAILED_CANDIDATES)
    candidates.resize(MAX_DETAILED_CANDIDATES);

  // Monte Carlo sampling with adaptive sample sizes
  for (auto &candidate : candidates) {
    double totalVal = 0;
    int actualSamples =
        (candidate.move.type == SCORE)
            ? samples
            : (emptyCells < 30 ? samples * 3
                               : samples * 2); // More samples in late game

    for (int s = 0; s < actualSamples; s++) {
      int addCells = (candidate.move.type == SCORE) ? 3 : 5;
      BoardState outcome =
          addRandomCells(candidate.nextState, addCells, colors, rows, cols);
      double eval = computeFinalScore(outcome, isGoal, minScore, rows, cols);
      // Progressive sampling weights
      double weight = (s < samples) ? 1.0 : (s < samples * 2) ? 0.7 : 0.5;
      totalVal += eval * weight;
    }
    candidate.value = totalVal / (actualSamples * 0.8); // Normalized evaluation
  }

  auto bestCandidate =
      max_element(candidates.begin(), candidates.end(),
                  [](const MoveCandidate &a, const MoveCandidate &b) {
                    return a.value < b.value;
                  });

  return bestCandidate->move;
}