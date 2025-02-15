// ––––– ga.h –––––
#ifndef GA_H
#define GA_H

#include <vector>
#include "superball.h" // Include Superball class definition

struct Chromosome {
    std::vector<double> weights; // Heuristic weights for features
    double fitness;
};

class GeneticAlgorithm {
public:
    GeneticAlgorithm(int population_size, int num_generations, double mutation_rate,
                     double crossover_rate, int num_features, int num_runs_per_fitness);
    ~GeneticAlgorithm();

    // Runs the GA on a given Superball board and returns the best chromosome.
    Chromosome run(Superball &sb);

    // Returns a copy of the best chromosome found so far.
    Chromosome get_best_chromosome() const;

private:
    std::vector<Chromosome> population;
    int population_size;
    int num_generations;
    double mutation_rate;
    double crossover_rate;
    int num_features;
    int num_runs_per_fitness;

    double evaluate_fitness(Chromosome &chromosome, Superball &sb);
    Chromosome selection();
    Chromosome crossover(const Chromosome &parent1, const Chromosome &parent2);
    void mutation(Chromosome &chromosome);
};

#endif