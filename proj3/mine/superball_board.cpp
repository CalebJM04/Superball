#include "superball_board.hpp"
#include <cstdlib>

SuperballBoard::SuperballBoard(int rows, int cols) : r(rows), c(cols), cells(rows * cols, '.') {}

bool SuperballBoard::Swap(int idxA, int idxB)
{
    if (cells[idxA] == '.' || cells[idxB] == '.') return false;
    std::swap(cells[idxA], cells[idxB]);
    return true;
}

int SuperballBoard::ScoreGreen()
{
    int scored = 0;
    for (size_t i = 0; i < cells.size(); i++)
    {
        if (cells[i] == 'g' || cells[i] == 'G')
        {
            cells[i] = '.';
            scored++;
        }
    }
    return scored;
}

bool SuperballBoard::NoMovesLeft()
{
    return false; // For now, assume we always have moves.
}

void SuperballBoard::PrintBoard()
{
    for (int i = 0; i < r; i++)
    {
        for (int j = 0; j < c; j++)
        {
            std::cout << cells[i * c + j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}
