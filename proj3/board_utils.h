#pragma once

#include <vector>
#include <string>
#include "types.h"

/**
 * Returns the canonical (lowercase) color of a piece.
 * @param ch The character representing the piece
 * @return The lowercase version of the character
 */
char pieceColor(char ch);

/**
 * Analyzes the board to find connected components that meet scoring criteria.
 * @param board The current game board
 * @param isGoal Matrix indicating which cells are goal cells
 * @param minScore Minimum group size required for scoring
 * @param rows Number of rows in the board
 * @param cols Number of columns in the board
 * @return Vector of GroupInfo structures for qualifying groups
 */
std::vector<GroupInfo> analyzeBoard(const std::vector<std::string>& board, 
                                   const std::vector<std::vector<bool>>& isGoal,
                                   int minScore, int rows, int cols);

/**
 * Evaluates the current board state using a heuristic function.
 * The evaluation considers:
 * - Existing scoreable groups
 * - Potential future groups
 * - Color values (p=2, b=3, y=4, r=5, g=6)
 * - Group sizes and proximity to goals
 * 
 * @param state Current board state
 * @param isGoal Matrix indicating which cells are goal cells
 * @param minScore Minimum group size required for scoring
 * @param rows Number of rows in the board
 * @param cols Number of columns in the board
 * @return A double value representing the board's evaluation score
 */
double evaluateState(const BoardState& state, 
                    const std::vector<std::vector<bool>>& isGoal,
                    int minScore, int rows, int cols);