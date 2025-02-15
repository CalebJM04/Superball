// ––––– superball.cpp ––––– (Corrected Again, with DisjointSet details)
#include "superball.h"
#include <iostream>
#include <set>
#include <cmath>
#include <cstdlib>
#include <cctype>
#include <algorithm>

Superball::Superball(int rows, int cols, int mss, const std::string &colorsStr)
    : r(rows), c(cols), mss(mss), empty(r * cols), ds(nullptr) // Initialize empty correctly
{
    board.resize(r * c, '.');
    goals.resize(r * c, 0);
    colorValues.resize(256, 0); // Initialize all to 0

    // Set up color values: first color gets value 2, then increasing by 1.
    for (size_t i = 0; i < colorsStr.size(); i++) {
        char lc = colorsStr[i];
        colorValues[lc] = 2 + (int)i; // Cast i to int for correct arithmetic
    }

    // Initialize goal cells (rows 2-5, columns 0, 1, 8, and 9)
    for (int i = 2; i <= 5; ++i) {
        if (i < r) { // Check bounds for rows
            if (0 < c) goals[i * c + 0] = 1; //col 0
            if (1 < c) goals[i * c + 1] = 1; //col 1
            if (8 < c) goals[i * c + 8] = 1; //col 8
            if (9 < c) goals[i * c + 9] = 1; //col 9
        }
    }

    // Create a disjoint set instance for r*c elements.
     ds = new DisjointSetByRankWPC(r * c);
}

// Destructor
Superball::~Superball() {
    if (ds) delete ds;  // Correctly deallocate
}

// Copy Constructor (Deep Copy)
Superball::Superball(const Superball &other) :
    r(other.r), c(other.c), mss(other.mss), empty(other.empty) {
    deep_copy(other); // Use a helper function for the deep copy
}

// Copy Assignment Operator (Deep Copy)
Superball &Superball::operator=(const Superball &other) {
    if (this != &other) { // Check for self-assignment (sb = sb;)
        // Clean up existing resources (important!)
        delete ds; // Delete the *current* DisjointSet

        // Now, do a deep copy
        r = other.r;
        c = other.c;
        mss = other.mss;
        empty = other.empty;
        deep_copy(other);
    }
    return *this; // Return a reference to the current object
}

// Helper function for deep copying
void Superball::deep_copy(const Superball &other) {
    board = other.board;
    goals = other.goals;
    colorValues = other.colorValues;

    // Deep copy of the DisjointSet.  *Crucially*, create a *new* DisjointSet.
    if (other.ds != nullptr) {
        ds = new DisjointSetByRankWPC(r * c); // Allocate *new* memory

        // Cast to DisjointSetByRankWPC* to access links and ranks
        DisjointSetByRankWPC *this_ds = static_cast<DisjointSetByRankWPC*>(ds);
        const DisjointSetByRankWPC *other_ds = static_cast<const DisjointSetByRankWPC*>(other.ds);

        // Copy the links and ranks vectors.
        //  this_ds->links = other_ds->links;
        //  this_ds->ranks = other_ds->ranks;

    } else {
        ds = nullptr;
    }
}


void Superball::analyze_board() {
    if (ds) delete ds;
    ds = new DisjointSetByRankWPC(r * c);

    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            int index = i * c + j;
            char cell = board[index];
            if (cell == '.' || cell == '*') continue;

            if (j + 1 < c) {
                int right_index = i * c + (j + 1);
                if (board[index] == board[right_index]) {
                    ds->Union(ds->Find(index), ds->Find(right_index));
                }
            }
            if (i + 1 < r) {
                int down_index = (i + 1) * c + j;
                if (board[index] == board[down_index]) {
                    ds->Union(ds->Find(index), ds->Find(down_index));
                }
            }
        }
    }
}

std::vector<double> Superball::get_features() {
    analyze_board();

    std::set<int> reps;
    std::vector<int> set_sizes(r * c, 0);
    std::vector<bool> has_goal(r * c, false);

    for (int i = 0; i < r * c; i++) {
        char cell = board[i];
        if (cell == '.' || cell == '*') continue;
        int root = ds->Find(i);
        reps.insert(root);
        set_sizes[root]++;
        if (goals[i] == 1) {
            has_goal[root] = true;
        }
    }

    int num_sets = reps.size();
    int largest_set = 0;
    int scoreable_sets = 0;
    int potential_score = 0;

    for (int root : reps) {
        if (set_sizes[root] > largest_set) {
            largest_set = set_sizes[root];
        }

        if (set_sizes[root] >= mss && has_goal[root]) {
            scoreable_sets++;
            char color = 0;
            for (int j = 0; j < r * c; ++j) {
                if (board[j] != '.' && board[j] != '*' && ds->Find(j) == root) {
                    color = board[j];
                    break;
                }
            }
            if (color != 0) {
                potential_score += set_sizes[root] * colorValues[color];
            }
        }
    }

    std::vector<double> features;
    features.push_back(num_sets);
    features.push_back(largest_set);
    features.push_back(scoreable_sets);
    features.push_back(potential_score);
    return features;
}
void Superball::apply_move(int r1, int c1, int r2, int c2, bool is_score) {
    int index1 = r1 * c + c1;
    if (!is_score) {
        int index2 = r2 * c + c2;
        // SWAP
        std::swap(board[index1], board[index2]);
    } else {
        // SCORE
        int root = ds->Find(index1);
        for (int i = 0; i < r * c; i++) {
                if (board[i] != '.' && board[i] != '*' && ds->Find(i) == root) {
                    if(goals[i] == 1){
                        board[i] = '*';
                    }
                    else{
                        board[i] = '.';
                    }
                    empty++;
                }
            }
    }
}
std::vector<char> Superball::get_board_copy() {
    return board;
}

void Superball::set_board(const std::vector<char> &new_board) {
    board = new_board;
    empty = 0;
    for (char cell : board) {
        if (cell == '.' || cell == '*') empty++;
    }
}