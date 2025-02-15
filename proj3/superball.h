// ––––– superball.h –––––
#ifndef SUPERBALL_H
#define SUPERBALL_H

#include <vector>
#include <string>
#include "disjoint.h"

class Superball {
public:
    // Constructor
    Superball(int rows, int cols, int mss, const std::string &colorsStr);

    // Destructor
    ~Superball();

    // Copy Constructor (NEW!)
    Superball(const Superball &other);

    // Copy Assignment Operator (NEW!)
    Superball &operator=(const Superball &other);


    int r, c, mss;
    int empty;
    std::vector<char> board;
    std::vector<int> goals;
    std::vector<int> colorValues;
    DisjointSet *ds;

    void analyze_board();
    std::vector<double> get_features();
    void apply_move(int r1, int c1, int r2, int c2, bool is_score);
    std::vector<char> get_board_copy();
    void set_board(const std::vector<char> &new_board);

private:
    // Helper function for deep copying (NEW!)
    void deep_copy(const Superball &other);
};

#endif