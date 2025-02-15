#pragma once

#include <prereqs.hpp>
#include <environment.hpp>
#include "superball_board.hpp"

class SuperballEnv
{
public:
    static constexpr size_t ROWS = 8;
    static constexpr size_t COLS = 10;
    static constexpr size_t MAX_ACTIONS = 1 + 30; // Score + 30 swaps

    SuperballEnv();
    static size_t ActionSize;
    
    arma::colvec InitialSample();
    bool IsTerminal(const arma::colvec& state);
    double Step(arma::colvec& state, const size_t action, arma::colvec& nextState);

private:
    SuperballBoard board;
    std::vector<std::pair<int,int>> swapList;
    
    arma::colvec BoardToState();
    double PerformAction(size_t action);
};
