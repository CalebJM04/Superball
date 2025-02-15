#pragma once

#include <vector>
#include <string>

// Defines the type of move that can be made
enum MoveType { 
    SWAP,   // Swapping two pieces
    SCORE   // Scoring a group of pieces 
};

// Structure to represent a move in the game
struct Move {
    MoveType type;
    // For SWAP: use (r1,c1) and (r2,c2)
    // For SCORE: only r1 and c1 are used
    int r1, c1;     // First position (required for both move types)
    int r2, c2;     // Second position (only used for SWAP)
};

// Structure to represent the game board state
struct BoardState {
    std::vector<std::string> board;
};

// Structure to hold information about a connected component/group
struct GroupInfo {
    int root;           // Root cell of the group in the disjoint set
    int size;           // Number of pieces in the group
    bool hasGoal;       // Whether the group contains a goal cell
    int goalR, goalC;   // Position of a representative goal cell (if any)
};