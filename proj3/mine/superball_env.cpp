#include "superball_env.hpp"
#include <cassert>

SuperballEnv::SuperballEnv() : board(ROWS, COLS)
{
    ActionSize = MAX_ACTIONS;
    for (int i = 0; i < ROWS * COLS; i++)
    {
        if (i % COLS < COLS - 1) swapList.emplace_back(i, i + 1);
        if (i / COLS < ROWS - 1) swapList.emplace_back(i, i + COLS);
    }
}

size_t SuperballEnv::ActionSize = 0;

arma::colvec SuperballEnv::InitialSample()
{
    return BoardToState();
}

bool SuperballEnv::IsTerminal(const arma::colvec& state)
{
    return board.NoMovesLeft();
}

double SuperballEnv::Step(arma::colvec& state, const size_t action, arma::colvec& nextState)
{
    double reward = PerformAction(action);
    nextState = BoardToState();
    return reward;
}

double SuperballEnv::PerformAction(size_t action)
{
    if (action == 0) return board.ScoreGreen();
    auto [a, b] = swapList[action - 1];
    return board.Swap(a, b) ? 0.1 : -0.1;
}

arma::colvec SuperballEnv::BoardToState()
{
    arma::colvec st(ROWS * COLS);
    for (size_t i = 0; i < ROWS * COLS; i++) st(i) = EncodeCell(board.cells[i]);
    return st;
}
