#include <mlpack/methods/reinforcement_learning/q_learning.hpp>
#include "superball_env.hpp"

using namespace mlpack::rl;

int main()
{
    SuperballEnv env;
    SuperballEnv::ActionSize = 1 + env.swapList.size();

    TrainingConfig config;
    config.ExplorationSteps() = 1000;
    config.ExplorationDecay() = 0.999;
    config.ExplorationMin() = 0.05;
    config.StepSize() = 0.01;
    config.Discount() = 0.95;

    const size_t stateDim = SuperballEnv::ROWS * SuperballEnv::COLS; 
    const size_t actionDim = SuperballEnv::ActionSize;

    arma::mat network;
    network.randu(stateDim, actionDim); // Random weight initialization

    QLearning<SuperballEnv> agent(config, env, network);

    for (size_t episode = 0; episode < 1000; ++episode)
    {
        agent.Episode();
    }

    network.save("trained_model.bin");
    return 0;
}
