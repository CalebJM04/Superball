#pragma once

#include "types.h"
#include <vector>
#include <string>

/**
 * Adds random pieces to empty cells on the board.
 * @param state Current board state
 * @param numCells Number of random cells to add
 * @param colors String containing possible colors to choose from
 * @param rows Number of rows in the board
 * @param cols Number of columns in the board
 * @return New board state with random cells added
 */
BoardState addRandomCells(const BoardState& state, 
                         int numCells, 
                         const std::string& colors, 
                         int rows, 
                         int cols);

/**
 * Simulates a swap move between two positions on the board.
 * @param state Current board state
 * @param r1 Row of first piece
 * @param c1 Column of first piece
 * @param r2 Row of second piece
 * @param c2 Column of second piece
 * @return New board state after the swap
 */
BoardState simulateSwap(const BoardState& state, 
                       int r1, 
                       int c1, 
                       int r2, 
                       int c2);

/**
 * Simulates scoring a group of connected pieces.
 * Starting from the given position, removes all connected pieces of the same color.
 * Goal cells become '*' when cleared, regular cells become '.'.
 * 
 * @param state Current board state
 * @param startR Starting row position
 * @param startC Starting column position
 * @param rows Number of rows in the board
 * @param cols Number of columns in the board
 * @param isGoal Matrix indicating which cells are goal cells
 * @return New board state after the scoring move
 */
BoardState simulateScoreMove(const BoardState& state, 
                            int startR, 
                            int startC, 
                            int rows, 
                            int cols,
                            const std::vector<std::vector<bool>>& isGoal);

/**
 * Chooses the best move based on current board evaluation.
 * Considers both swap and score moves and selects the one that 
 * maximizes the board position quality without simulating random future pieces.
 * 
 * Strategy:
 * 1. Immediately scores groups above threshold size
 * 2. Evaluates each possible move based on resulting board quality
 * 3. Prioritizes pieces in strategic positions (corners/edges)
 * 
 * @param state Current board state
 * @param isGoal Matrix indicating which cells are goal cells
 * @param minScore Minimum group size required for scoring
 * @param rows Number of rows in the board
 * @param cols Number of columns in the board
 * @param colors String containing possible colors
 * @return The chosen Move
 */
Move chooseBestMove(const BoardState& state, 
                   const std::vector<std::vector<bool>>& isGoal,
                   int minScore, 
                   int rows, 
                   int cols, 
                   const std::string& colors);

/**
 * Legacy function for backward compatibility.
 * Calls chooseBestMove internally.
 */
Move chooseMoveExpectimax(const BoardState& state, 
                         const std::vector<std::vector<bool>>& isGoal,
                         int minScore, 
                         int rows, 
                         int cols, 
                         const std::string& colors, 
                         int samples = 10);