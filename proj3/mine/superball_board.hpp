#pragma once
#include <vector>
#include <map>
#include <iostream>

class SuperballBoard
{
public:
    int r, c;
    std::vector<char> cells;
    SuperballBoard(int rows, int cols);

    bool Swap(int idxA, int idxB);
    int ScoreGreen();
    bool NoMovesLeft();
    void PrintBoard();
};

// Mapping function
inline int EncodeCell(char c)
{
    if (c == '.' || c == '*') return 0;
    else if (c == 'g' || c == 'G') return 1;
    else if (c == 'r' || c == 'R') return 2;
    else if (c == 'b' || c == 'B') return 3;
    else if (c == 'p' || c == 'P') return 4;
    return 5;
}
