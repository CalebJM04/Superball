// ––––– ga.cpp –––––
#include "ga.h"
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <iostream>

GeneticAlgorithm::GeneticAlgorithm(int population_size, int num_generations, double mutation_rate,
                                   double crossover_rate, int num_features, int num_runs_per_fitness)
    : population_size(population_size), num_generations(num_generations),
      mutation_rate(mutation_rate), crossover_rate(crossover_rate),
      num_features(num_features), num_runs_per_fitness(num_runs_per_fitness)
{
    std::srand(std::time(0)); // Seed the random number generator

    // Initialize population with random weights (range -1 to 1)
    for (int i = 0; i < population_size; i++) {
        Chromosome chrom;
        chrom.weights.resize(num_features, 0.0); // Initialize weights
        for (int j = 0; j < num_features; j++) {
            chrom.weights[j] = ((double)std::rand() / RAND_MAX) * 2 - 1;
        }
        chrom.fitness = 0; // Initialize fitness
        population.push_back(chrom);
    }
}

GeneticAlgorithm::~GeneticAlgorithm() {} // Nothing to do in the destructor

double GeneticAlgorithm::evaluate_fitness(Chromosome &chromosome, Superball &sb) {
    double total_score = 0;
    int moves = 50; // Simulate up to this number of moves.  Adjust as needed.

    for (int run = 0; run < num_runs_per_fitness; run++) {
        Superball sb_copy = sb; // Create a deep copy

        for (int m = 0; m < moves; m++) { // Simulate up to 'moves' moves.
            sb_copy.analyze_board();
            std::vector<double> features = sb_copy.get_features();
            double best_move_score = -1e9; // Initialize to a very low value
            bool best_is_score = false;
            int best_r1 = 0, best_c1 = 0, best_r2 = 0, best_c2 = 0, best_score_r = 0, best_score_c = 0;

            // Try SWAP moves (adjacent swaps)
            for (int r1 = 0; r1 < sb_copy.r; r1++) {
                for (int c1 = 0; c1 < sb_copy.c; c1++) {
                    // Swap with right
                    if (c1 + 1 < sb_copy.c) {
                        Superball temp_sb = sb_copy;
                        temp_sb.ds = new DisjointSetByRankWPC(sb_copy.r * sb_copy.c); // Create new DS
                        temp_sb.apply_move(r1, c1, r1, c1 + 1, false);
                        temp_sb.analyze_board();
                        std::vector<double> temp_features = temp_sb.get_features();
                        double score = 0;
                        for (size_t k = 0; k < chromosome.weights.size(); k++) {
                            score += chromosome.weights[k] * temp_features[k];
                        }

                        if (score > best_move_score) {
                            best_move_score = score;
                            best_is_score = false;
                            best_r1 = r1;
                            best_c1 = c1;
                            best_r2 = r1;
                            best_c2 = c1 + 1;
                        }
                        delete temp_sb.ds; // Clean up temporary DS
                    }
                    // Swap with down
                    if (r1 + 1 < sb_copy.r) {
                        Superball temp_sb = sb_copy;
                         temp_sb.ds = new DisjointSetByRankWPC(sb_copy.r * sb_copy.c); // Create new DS
                        temp_sb.apply_move(r1, c1, r1 + 1, c1, false);
                        temp_sb.analyze_board();
                        std::vector<double> temp_features = temp_sb.get_features();
                        double score = 0;
                        for (size_t k = 0; k < chromosome.weights.size(); k++) {
                            score += chromosome.weights[k] * temp_features[k];
                        }

                        if (score > best_move_score) {
                            best_move_score = score;
                            best_is_score = false;
                            best_r1 = r1;
                            best_c1 = c1;
                            best_r2 = r1 + 1;
                            best_c2 = c1;
                        }
                        delete temp_sb.ds; //Clean up temporary DS
                    }
                }
            }

            // Try SCORE moves
                for (int score_r = 0; score_r < sb_copy.r; score_r++) {
                    for (int score_c = 0; score_c < sb_copy.c; score_c++) {
                        int index = score_r * sb_copy.c + score_c;
                        if (sb_copy.goals[index] == 1 && sb_copy.board[index] != '.' && sb_copy.board[index] != '*') {
                             Superball temp_sb = sb_copy;
                            temp_sb.ds = new DisjointSetByRankWPC(sb_copy.r*sb_copy.c);
                            temp_sb.analyze_board(); // Analyze before scoring
                            int root = temp_sb.ds->Find(index);
                            std::vector<int> set_sizes(temp_sb.r * temp_sb.c, 0);
                            for(int i = 0; i < temp_sb.r * temp_sb.c; i++){ //Calculate set sizes
                                if(temp_sb.board[i] != '.' && temp_sb.board[i] != '*'){
                                    set_sizes[temp_sb.ds->Find(i)]++;
                                }
                            }
                            if (set_sizes[root] >= temp_sb.mss) { // Ensure it's a valid score
                                double score_before = temp_sb.get_features()[3];  //Potential score before
                                temp_sb.apply_move(score_r, score_c, 0, 0, true); // 0,0 are dummies
                                double score_after = temp_sb.get_features()[3]; //Potential score after

                                double move_score = score_before - score_after; //Difference is score gained
                                if (move_score > best_move_score) {
                                    best_move_score = move_score;
                                    best_is_score = true;
                                    best_score_r = score_r;
                                    best_score_c = score_c;
                                }
                            }
                            delete temp_sb.ds;
                        }
                    }
                }


            // Apply the best move (to the copy!)
            if (best_is_score) {
                sb_copy.apply_move(best_score_r, best_score_c, 0, 0, true);
                total_score += best_move_score; // Add to the total score *after* applying the score.
            } else {
                sb_copy.apply_move(best_r1, best_c1, best_r2, best_c2, false);
            }

            if(sb_copy.empty <= 5) break; //End if there are less than 5 moves left
        }
        total_score += sb_copy.get_features()[3]; // Add final potential score if any is left
    }
    return total_score / num_runs_per_fitness;
}



Chromosome GeneticAlgorithm::selection() {
    // Tournament selection: pick two random individuals and choose the fitter.
    int idx1 = std::rand() % population_size;
    int idx2 = std::rand() % population_size;
    return (population[idx1].fitness > population[idx2].fitness) ? population[idx1] : population[idx2];
}

Chromosome GeneticAlgorithm::crossover(const Chromosome &parent1, const Chromosome &parent2) {
    Chromosome child;
    child.weights.resize(num_features);

    // Single-point crossover (you can experiment with other methods)
    int crossover_point = std::rand() % num_features;
    for (int i = 0; i < num_features; i++) {
        if (i < crossover_point) {
            child.weights[i] = parent1.weights[i];
        } else {
            child.weights[i] = parent2.weights[i];
        }
    }
    //Alternative average crossover
    /*for (int i = 0; i < num_features; i++) {
        child.weights[i] = (parent1.weights[i] + parent2.weights[i])/2;
    }*/

    return child;
}

void GeneticAlgorithm::mutation(Chromosome &chromosome) {
    for (int i = 0; i < num_features; i++) {
        if (((double)std::rand() / RAND_MAX) < mutation_rate) {
            // Add a small random value (e.g., Gaussian noise)
            chromosome.weights[i] += ((double)std::rand() / RAND_MAX) * 0.2 - 0.1; // Adjust range as needed
        }
    }
}

Chromosome GeneticAlgorithm::run(Superball &sb) {
    // Evaluate initial population.
    for (auto &chrom : population) {
        chrom.fitness = evaluate_fitness(chrom, sb);
    }

    for (int gen = 0; gen < num_generations; gen++) {
        std::vector<Chromosome> new_population;
        while ((int)new_population.size() < population_size) {
            Chromosome parent1 = selection();
            Chromosome parent2 = selection();
            Chromosome child;
            if (((double)std::rand() / RAND_MAX) < crossover_rate) {
                child = crossover(parent1, parent2);
            } else {
                // If no crossover, clone a parent (e.g., parent1)
                child = parent1;
            }
            mutation(child);
            child.fitness = evaluate_fitness(child, sb);  // Evaluate *after* mutation
            new_population.push_back(child);
        }
        population = new_population;

        // Optionally, output the best fitness of this generation.
        double best_fit = -1e9;
        for (auto &chrom : population) {
            if (chrom.fitness > best_fit) {
                best_fit = chrom.fitness;
            }
        }
        std::cout << "Generation " << gen << " best fitness: " << best_fit << std::endl;
    }

    // Return the best chromosome.
    Chromosome best = population[0];
    for (auto &chrom : population) {
        if (chrom.fitness > best.fitness) {
            best = chrom;
        }
    }
    return best;
}

Chromosome GeneticAlgorithm::get_best_chromosome() const {
    Chromosome best = population[0];
    for (auto &chrom : population) {
        if (chrom.fitness > best.fitness) {
            best = chrom;
        }
    }
    return best;
}