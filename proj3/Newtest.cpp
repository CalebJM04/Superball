#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include "disjoint.h"

using namespace std;

class SuperballGame {
public:
    int rows, cols;
    vector<vector<char>> board;
    set<pair<int, int>> goalCells;
    DisjointSetByRankWPC ds;

    SuperballGame(vector<vector<char>> initBoard) 
        : board(initBoard), rows(initBoard.size()), cols(initBoard[0].size()), ds(rows * cols) {
        initializeGoalCells();
    }

    void initializeGoalCells() {
        for (int r = 2; r <= 5; ++r) {
            goalCells.insert({r, 0});
            goalCells.insert({r, 1});
            goalCells.insert({r, 8});
            goalCells.insert({r, 9});
        }
    }

    int getIndex(int r, int c) { return r * cols + c; }

    void findGroups() {
        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                if (board[r][c] == '.' || board[r][c] == '*') continue;
                
                if (r + 1 < rows && board[r][c] == board[r + 1][c]) 
                    ds.Union(getIndex(r, c), getIndex(r + 1, c));
                if (c + 1 < cols && board[r][c] == board[r][c + 1]) 
                    ds.Union(getIndex(r, c), getIndex(r, c + 1));
            }
        }
    }

    vector<pair<string, vector<int>>> generateMoves() {
        vector<pair<string, vector<int>>> moves;

        // Generate score moves
        for (auto [r, c] : goalCells) {
            int groupRoot = ds.Find(getIndex(r, c));
            int count = countGroupSize(groupRoot);
            if (count >= 5) moves.push_back({"SCORE", {r, c}});
        }

        // Generate swap moves
        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                if (board[r][c] == '.') continue;
                
                vector<pair<int, int>> directions = {{1, 0}, {0, 1}};
                for (auto [dr, dc] : directions) {
                    int nr = r + dr, nc = c + dc;
                    if (nr < rows && nc < cols && board[nr][nc] != '.') 
                        moves.push_back({"SWAP", {r, c, nr, nc}});
                }
            }
        }
        return moves;
    }

    int countGroupSize(int root) {
        int size = 0;
        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                if (ds.Find(getIndex(r, c)) == root) size++;
            }
        }
        return size;
    }

    int evaluateBoard() {
        int scoreOpportunities = 0;
        int nearScoreOpportunities = 0;
        
        for (auto [r, c] : goalCells) {
            int groupRoot = ds.Find(getIndex(r, c));
            int count = countGroupSize(groupRoot);
            if (count >= 5) scoreOpportunities++;
            else if (count >= 3) nearScoreOpportunities++;
        }
        return scoreOpportunities * 10 + nearScoreOpportunities * 5;
    }

    pair<string, vector<int>> bestMove() {
        vector<pair<string, vector<int>>> moves = generateMoves();
        
        for (auto move : moves) {
            if (move.first == "SCORE") return move;
        }

        pair<string, vector<int>> bestMove;
        int bestScore = -1;

        for (auto move : moves) {
            if (move.first == "SWAP") {
                int r1 = move.second[0], c1 = move.second[1];
                int r2 = move.second[2], c2 = move.second[3];
                
                swap(board[r1][c1], board[r2][c2]);
                findGroups();
                int score = evaluateBoard();
                swap(board[r1][c1], board[r2][c2]);

                if (score > bestScore) {
                    bestScore = score;
                    bestMove = move;
                }
            }
        }
        return bestMove;
    }
};

int main() {
    vector<vector<char>> board = {
        {'.', '.', '.', '.', '.', '.', '.', '.', '.', '.'},
        {'.', '.', '.', '.', 'p', '.', '.', '.', '.', '.'},
        {'B', '*', '.', '.', '.', '.', '.', '.', '*', '*'},
        {'*', '*', '.', '.', 'b', '.', '.', '.', '*', '*'},
        {'*', '*', '.', '.', 'g', '.', '.', '.', '*', '*'},
        {'*', '*', '.', '.', '.', '.', '.', '.', '*', '*'},
        {'.', '.', 'g', '.', '.', '.', '.', '.', '.', '.'},
        {'.', '.', '.', '.', '.', '.', '.', '.', '.', '.'},
    };

    SuperballGame game(board);
    game.findGroups();
    pair<string, vector<int>> move = game.bestMove();
    
    cout << move.first;
    for (int i : move.second) cout << " " << i;
    cout << endl;
    
    return 0;
}
