#include "common.h"

#include <chrono>
#include <random>
#include <numeric>
#include <utility>
#include <unordered_set>

void PrintUsage(const char* Exe, std::FILE* File){
    std::fprintf(File, "Usage: %s <edge-list> -s <subset-size> -k <seed_size> -p <infection-prob> -a <MC runs> -r <retries> -m <mutation-rate> -w <- weigh-by-degree\n", Exe);
}

template <typename RngT>
bool Infect(double p, RngT &Rng){
    if (p == 1.0) {
        return true;
    }
    std::uniform_real_distribution<> Dist(0, 1);
    return Dist(Rng) < p;
}

// It seems there may be different weights in each direction... I guess that's fine?
template <typename RngT>
std::pair<node, bool> Infect(const digraph &G, const edge_probs &E, node x, RngT &Rng){
    std::vector<double> weights(G[x].Out.size());
    std::generate(weights.begin(), weights.end(), [i=0, &E, &G, &x]() mutable {
        node n = G[x].Out[i];
        auto e = edge{x, n};
        auto weight = E.at(e);
        i++;
        return weight;
    });
    weights.push_back(1 - std::accumulate(weights.begin(), weights.end(), decltype(weights)::value_type(0.0)));
    std::discrete_distribution<> Dist(weights.begin(), weights.end());
    auto n = Dist(Rng);
    if(n+1 > G[x].Out.size()){
        return std::pair<node, bool> (n, false);
    }
    return std::pair<node, bool> (G[x].Out[n], true);
}

nodes AllNeighbors(const digraph &G, node x){
    auto neighbors = G[x].In;
    neighbors.insert(neighbors.end(), G[x].Out.begin(), G[x].Out.end());
    return neighbors;
}

// Treat it like an undirected graph
template <typename RngT>
int RandomNeighbor(const digraph &G, node x, RngT &Rng){
    auto neighbors = AllNeighbors(G, x);
    std::uniform_int_distribution<> Dist(0, neighbors.size() - 1);
    return neighbors[Dist(Rng)];
}

std::unordered_set<node> AllOutNeighbors(const digraph &G, const nodes &S){
    std::unordered_set<node> N;
    for (auto &x : S)
    {
        for (auto &n : G[x].Out)
        {
            N.insert(n);
        }
    }
    return N;
}

using TimeLapse = std::vector<std::pair<int,int>>;
struct MCResult {
    bool success;
    int influence_extreme=0;
    TimeLapse time_lapse;
};
using MCResults = std::vector<MCResult>;
// (probability, start-set, simulation average time, mutations to convergence)
struct MoranResult {
    double probability = 0.0;
    nodes seed;
    double sim_time = 0.0;
    int mutations = 0;
    MCResults monte;
};
using MoranResults = std::vector<MoranResult>;

// Print statistics of nodes
void MoranKStats(const digraph &G, const std::vector<double> &P, MoranResult &i, double &mean_degree, double &mean_p){
    for (auto &n : i.seed){
        mean_degree += G[n].Out.size() + G[n].In.size();
        mean_p += P[n];
    }
    mean_degree = mean_degree/i.seed.size();
    mean_p = mean_p/i.seed.size();
}

// Print statistics of nodes
void MoranKStats(const digraph &G, const edge_probs &E, MoranResult &i, double &mean_degree, double &mean_p){
    for (auto &n : i.seed){
        // Only include out-degree for this model
        mean_degree += G[n].Out.size();
        for (auto &v : G[n].Out){
            mean_p += E.at(edge{n, v});
        }
    }
    mean_degree = mean_degree/i.seed.size();
    mean_p = mean_p/i.seed.size();
}

// No fixed nodes.
template <typename RngT>
double MonteCarloEstimate(const digraph &G, const nodes &Infected, const bitvector &States, const std::vector<double> P, int k, int accuracy, RngT &Rng)
{
    int giveup, win;

    // Give up when we're at less than 1%... or we get stuck long time
    giveup = std::floor(G.size() * 0.01);
    // Win when we're within 1% of domination... or we get stuck long time
    win = std::floor(G.size() * 0.99);
    if (k <= giveup) {
        giveup = 0;
    } if (k >= win) {
        throw std::runtime_error("k is too big!");
    }

    int successes = 0;
    bitvector TrialStates;
    std::uniform_int_distribution<> Dist(0, G.size() - 1);

    for (int i=0; i < accuracy; i++)
    {
        TrialStates = States;
        auto progress = k;

        while (progress < win && progress > giveup)
        {
            auto x = Dist(Rng);

            if (Infect(P[x], Rng))
            {
                auto n = RandomNeighbor(G, x, Rng);

                // We've made progress
                if (TrialStates[n] == true && TrialStates[x] != true)
                {
                    progress++;
                    TrialStates.set(x, true);
                }
                // We've lost progress
                if (TrialStates[n] == false && TrialStates[x] == true)
                {
                    progress--;
                    TrialStates.set(x, false);
                }
            }
        }

        if (progress > giveup){
            successes++;
        }

    }
    return successes/(double)accuracy;
}


template <typename RngT>
double MonteCarloEstimate(const digraph &G, const nodes &Infected, const bitvector &States, const std::vector<double> P, int k, int accuracy, MoranResults &Results, RngT &Rng)
{
    MCResults mc_res;
    int giveup = 0, win = G.size();
    int successes = 0;
    bitvector TrialStates;

    for (int i=0; i < accuracy; i++)
    {
        MCResult mc;
        TrialStates = States;
        int progress = k, max = k, min = k, iteration=0;
        TimeLapse lapse = {std::pair<int,int>(iteration, progress)};
        std::uniform_int_distribution<> Dist(0, G.size() - 1);

        while (progress < win && progress > giveup)
        {
            auto x = Dist(Rng);

            if (Infect(P[x], Rng))
            {
                auto n = RandomNeighbor(G, x, Rng);

                // We've made progress
                if (TrialStates[n] == true && TrialStates[x] != true)
                {
                    progress++;
                    TrialStates.set(x, true);
                    if (progress > max){
                        max++;
                    }
                }
                // We've lost progress
                if (TrialStates[n] == false && TrialStates[x] == true)
                {
                    progress--;
                    TrialStates.set(x, false);
                    if (progress < min){
                        min--;
                    }
                }
            }
            if (lapse.back().second != progress) {
                lapse.push_back(std::pair<int,int>(iteration, progress));
            }
            iteration++;
        }

        if (progress > giveup){
            mc.success = true;
            mc.influence_extreme = min;
            successes++;
        }else {
            mc.success = false;
            mc.influence_extreme = max;
        }
        mc.time_lapse = lapse;
        mc_res.push_back(mc);

    }

    MCResults& current = Results.back().monte;
    current.insert(current.end(), mc_res.begin(), mc_res.end());
    return successes/(double)accuracy;
}

// No fixed nodes.
template <typename RngT>
double MonteCarloEstimate(const digraph &G, const nodes &Infected, const bitvector &States, const edge_probs E, int k, int accuracy, RngT &Rng)
{
    int giveup, win;

    // Give up when we're at less than 1%... or we get stuck long time
    giveup = 0;//std::floor(G.size() * 0.05);
    // Win when we're within ...% of domination... or we get stuck long time
    win = std::floor(G.size() * 0.9);
    if (k <= giveup) {
        giveup = 0;
    } if (k >= win) {
        throw std::runtime_error("k is too big!");
    }

    int successes = 0;
    bitvector TrialStates;
    std::uniform_int_distribution<> Dist(0, G.size() - 1);

    for (int i=0; i < accuracy; i++)
    {
        TrialStates = States;
        auto progress = k;

        while (progress < win && progress > giveup)
        {
            auto x = Dist(Rng);

            auto Inf = Infect(G, E, x, Rng);
            if (Inf.second)
            {
                auto n = Inf.first;

                // We've made progress
                if (TrialStates[n] == true && TrialStates[x] != true)
                {
                    progress++;
                    TrialStates.set(x, true);
                }
                // We've lost progress
                if (TrialStates[n] == false && TrialStates[x] == true)
                {
                    progress--;
                    TrialStates.set(x, false);
                }
            }
        }

        if (progress > giveup){
            successes++;
        }

    }
    return successes/(double)accuracy;
}


template <typename RngT>
double MonteCarloEstimate(const digraph &G, const nodes &Infected, const bitvector &States, const edge_probs E, int k, int accuracy, MoranResults &Results, RngT &Rng)
{
    MCResults mc_res;
    int giveup = 0, win = G.size();
    int successes = 0;
    bitvector TrialStates;

    for (int i=0; i < accuracy; i++)
    {
        MCResult mc;
        TrialStates = States;
        int progress = k, max = k, min = k, iteration=0;
        TimeLapse lapse = {std::pair<int,int>(iteration, progress)};
        std::uniform_int_distribution<> Dist(0, G.size() - 1);

        while (progress < win && progress > giveup){
            auto x = Dist(Rng);

            auto Inf = Infect(G, E, x, Rng);
            if (Inf.second){
                auto n = Inf.first;

                // We've made progress
                if (TrialStates[n] == true && TrialStates[x] != true){
                    progress++;
                    TrialStates.set(x, true);
                    if (progress > max){
                        max++;
                    }
                }
                // We've lost progress
                if (TrialStates[n] == false && TrialStates[x] == true){
                    progress--;
                    TrialStates.set(x, false);
                    if (progress < min){
                        min--;
                    }
                }
            }
            if (lapse.back().second != progress) {
                lapse.push_back(std::pair<int,int>(iteration, progress));
            }
            iteration++;
        }

        if (progress > giveup){
            mc.success = true;
            mc.influence_extreme = min;
            successes++;
        }else {
            mc.success = false;
            mc.influence_extreme = max;
        }
        mc.time_lapse = lapse;
        mc_res.push_back(mc);

    }

    MCResults& current = Results.back().monte;
    current.insert(current.end(), mc_res.begin(), mc_res.end());
    return successes/(double)accuracy;
}


nodes TopD(const digraph &G, int k) {
    nodes N(G.size());
    std::iota(N.begin(), N.end(), node(0));
    std::sort(N.begin(), N.end(), [&](const node a, const node b) {
            return G[a].In.size() + G[a].Out.size() > G[b].In.size() + G[b].Out.size();
        });
    N = nodes(N.begin(), N.begin()+k);
    return N;
}

class RandomIndex {
public:
    template <typename RngT>
    node operator()(const digraph &G, const nodes &N, const std::vector<double> P, RngT &Rng){
        std::uniform_int_distribution<> KDist(0, N.size() - 1);
        return KDist(Rng);
    }

    template <typename RngT>
    node operator()(const digraph &G, const nodes &N, const edge_probs &E, RngT &Rng){
        std::uniform_int_distribution<> KDist(0, N.size() - 1);
        return KDist(Rng);
    }
};

class PFitnessIndex {
public:
    template <typename RngT>
    node operator()(const digraph &G, const nodes &N, const std::vector<double> P, RngT &Rng){
        auto it = std::max_element(N.begin(), N.end(), [&](const node a, const node b) {
            return P[a] < P[b];
        });

        return std::distance(N.begin(), it);
    }
};

class DFitnessIndex {
public:
    template <typename RngT>
    node operator()(const digraph &G, const nodes &N, const std::vector<double> P, RngT &Rng){
        auto it = std::min_element(N.begin(), N.end(), [&](const node a, const node b) {
            return G[a].In.size() + G[a].Out.size() < G[b].In.size() + G[b].Out.size();
        });

        return std::distance(N.begin(), it);
    }
};

class D_Init {
    public:
        template <typename RngT>
        nodes operator()(const digraph &G, const std::vector<double> &P, int k, RngT Rng) {
            return TopD(G, k);
        }

};

class P_Init {
    public:
        template <typename RngT>
        nodes operator()(const digraph &G, const std::vector<double> &P, int k, RngT Rng) {
            nodes N(G.size());
            std::iota(N.begin(), N.end(), node(0));
            std::sort(N.begin(), N.end(), [&](const node a, const node b) {
                    return  P[a] < P[b];
                });
            N = nodes(N.begin(), N.begin()+k);
            return N;
        }

};

class R_Init {
    public:
        template <typename RngT>
        nodes operator()(const digraph &G, const std::vector<double> &P, int k, RngT Rng) {
            std::uniform_int_distribution<> Dist(0, G.size() - 1);
            std::unordered_set<node> InfectedSet;
            nodes Infected;
            // Select k nodes
            while(InfectedSet.size() < k)
            {
                auto x = Dist(Rng);
                InfectedSet.insert(x);
            }
            Infected.insert(Infected.begin(), InfectedSet.begin(), InfectedSet.end());

            return Infected;
        }

        template <typename RngT>
        nodes operator()(const digraph &G, const edge_probs &E, int k, RngT Rng) {
            std::uniform_int_distribution<> Dist(0, G.size() - 1);
            std::unordered_set<node> InfectedSet;
            nodes Infected;
            // Select k nodes
            while(InfectedSet.size() < k)
            {
                auto x = Dist(Rng);
                InfectedSet.insert(x);
            }
            Infected.insert(Infected.begin(), InfectedSet.begin(), InfectedSet.end());

            return Infected;
        }

};

template <typename FitT, typename RngT>
void Mutate(const digraph &G, nodes &Infected, bitvector &States, const std::vector<double> P, int accuracy, double mutateRate, MoranResults &Results, RngT &Rng, FitT F, bool weighted=false, bool graph = false)
{
    double average;
    int divisor = 1;

    auto start = std::chrono::steady_clock::now();
    auto currentProb = graph ? MonteCarloEstimate(G, Infected, States, P, Infected.size(), accuracy, Results, Rng) : MonteCarloEstimate(G, Infected, States, P, Infected.size(), accuracy, Rng);
    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration<double>(end - start).count();

    average = average + (duration)/divisor;
    divisor++;

    // Increase every time we fail to find a better set. Reset when we do. assume convergence if > 10 and abort.
    int stuck = 0;
    // number of (succesful) mutations before abort
    int mutations = 0;

    int swaps = std::floor(Infected.size() * mutateRate);

    std::uniform_int_distribution<> AllDist(0, G.size() - 1);

    for (;;)
    {
        int i = 0;
        nodes NewNodes, OldNodes, Mutation = Infected;
        bitvector NewStates = States;
        // Swap without replacement, order is important!
        while (i < swaps-1)
        {

            auto x = AllDist(Rng);
            if (weighted) {
                x = RandomNeighbor(G, x, Rng);
            }

            if (NewStates[x] == true)
            {
                continue;
            }
            NewStates.flip(x);

            // Dynamic as vector changes size
            auto index = F(G, Mutation, P, Rng);

            // Need to keep so we can flip state afterwards (doing it now would mean we could add it right back in)
            OldNodes.push_back(Mutation[index]);

            // Remove old node so it can't accidentally be swapped twice
            std::swap(Mutation[index], Mutation.back());
            Mutation.pop_back();

            // Need to keep so we can add to infected afterwards (doing it now would would mean we could remove it again immediately)
            NewNodes.push_back(x);

            i++;
        }

        Mutation.insert(Mutation.end(), NewNodes.begin(), NewNodes.end());

        for (auto &n : OldNodes) {
            NewStates.flip(n);
        }

        start = std::chrono::steady_clock::now();
        auto mutationProb = graph ? MonteCarloEstimate(G, Mutation, NewStates, P, Mutation.size(), accuracy, Results, Rng) : MonteCarloEstimate(G, Mutation, NewStates, P, Mutation.size(), accuracy, Rng);
        end = std::chrono::steady_clock::now();
        duration = std::chrono::duration<double>(end - start).count();
        average = average + (duration)/divisor;
        divisor++;

        if (mutationProb > currentProb)
        {
            std::swap(Infected, Mutation);
            std::swap(States, NewStates);

            currentProb = mutationProb;

            stuck = 0;
            mutations++;
        } else if (stuck > 9) {
            break;
        }
        else {
            stuck++;
        }
    }

    Results.back().probability = currentProb;
    Results.back().sim_time = average/accuracy;
    Results.back().mutations = mutations;
    Results.back().seed = Infected;
}

template <typename FitT, typename RngT>
void Mutate(const digraph &G, nodes &Infected, bitvector &States, const edge_probs E, int accuracy, double mutateRate, MoranResults &Results, RngT &Rng, FitT F, bool weighted=false, bool graph = false)
{
    double average;
    int divisor = 1;

    auto start = std::chrono::steady_clock::now();
    auto currentProb = graph ? MonteCarloEstimate(G, Infected, States, E, Infected.size(), accuracy, Results, Rng) : MonteCarloEstimate(G, Infected, States, E, Infected.size(), accuracy, Rng);
    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration<double>(end - start).count();

    average = average + (duration)/divisor;
    divisor++;

    // Increase every time we fail to find a better set. Reset when we do. assume convergence if > 10 and abort.
    int stuck = 0;
    // number of (succesful) mutations before abort
    int mutations = 0;

    int swaps = std::floor(Infected.size() * mutateRate);

    std::uniform_int_distribution<> AllDist(0, G.size() - 1);

    for (;;)
    {
        int i = 0;
        nodes NewNodes, OldNodes, Mutation = Infected;
        bitvector NewStates = States;
        // Swap without replacement, order is important!
        while (i < swaps-1)
        {

            auto x = AllDist(Rng);
            if (weighted) {
                x = RandomNeighbor(G, x, Rng);
            }

            if (NewStates[x] == true)
            {
                continue;
            }
            NewStates.flip(x);

            // Dynamic as vector changes size
            auto index = F(G, Mutation, E, Rng);

            // Need to keep so we can flip state afterwards (doing it now would mean we could add it right back in)
            OldNodes.push_back(Mutation[index]);

            // Remove old node so it can't accidentally be swapped twice
            std::swap(Mutation[index], Mutation.back());
            Mutation.pop_back();

            // Need to keep so we can add to infected afterwards (doing it now would would mean we could remove it again immediately)
            NewNodes.push_back(x);

            i++;
        }

        Mutation.insert(Mutation.end(), NewNodes.begin(), NewNodes.end());

        for (auto &n : OldNodes) {
            NewStates.flip(n);
        }

        start = std::chrono::steady_clock::now();
        auto mutationProb = graph ? MonteCarloEstimate(G, Mutation, NewStates, E, Mutation.size(), accuracy, Results, Rng) : MonteCarloEstimate(G, Mutation, NewStates, E, Mutation.size(), accuracy, Rng);
        end = std::chrono::steady_clock::now();
        duration = std::chrono::duration<double>(end - start).count();
        average = average + (duration)/divisor;
        divisor++;

        if (mutationProb > currentProb)
        {
            std::swap(Infected, Mutation);
            std::swap(States, NewStates);

            currentProb = mutationProb;

            stuck = 0;
            mutations++;
        } else if (stuck > 9) {
            break;
        }
        else {
            stuck++;
        }
    }

    Results.back().probability = currentProb;
    Results.back().sim_time = average/accuracy;
    Results.back().mutations = mutations;
    Results.back().seed = Infected;
}


template <typename FitT, typename InitT, typename RngT>
MoranResults MoranGenetic(const digraph &G, const std::vector<double> &P, int Retries, int accuracy, int k, double mutateRate, RngT &Rng, FitT F, InitT I, bool weighted=false, bool graph=false)
{
    MoranResults Results;
    bitvector States(G.size());

    for (int i = 0; i < Retries; i++)
    {
        if (i > 0) {
            States.clear();
        }

        auto Infected = I(G, P, k, Rng);

        // Weigh probabilities by degree
        if (weighted) {
            for (int j=0; j < Infected.size()-1; j++)
            {
                auto x = RandomNeighbor(G, Infected[j], Rng);
                States.set(x);
                Infected[j] = x;
            }
        }

        for (auto &n : Infected) {
            States.set(n);
        }

        Results.push_back(MoranResult());

        Mutate(G, Infected, States, P, accuracy, mutateRate, Results, Rng, F, graph);
    }

    std::sort(Results.begin(), Results.end(), [](MoranResult a, MoranResult b) {
        return a.probability > b.probability;
    });

    return Results;
}

template <typename FitT, typename InitT, typename RngT>
MoranResults MoranGenetic(const digraph &G, const edge_probs &E, int Retries, int accuracy, int k, double mutateRate, RngT &Rng, FitT F, InitT I, bool weighted=false, bool graph=false)
{
    MoranResults Results;
    bitvector States(G.size());

    for (int i = 0; i < Retries; i++)
    {
        if (i > 0) {
            States.clear();
        }

        auto Infected = I(G, E, k, Rng);

        // Weigh probabilities by degree
        if (weighted) {
            for (int j=0; j < Infected.size()-1; j++){
                auto x = RandomNeighbor(G, Infected[j], Rng);
                States.set(x);
                Infected[j] = x;
            }
        }

        for (auto &n : Infected) {
            States.set(n);
        }

        Results.push_back(MoranResult());

        Mutate(G, Infected, States, E, accuracy, mutateRate, Results, Rng, F, graph);
    }

    std::sort(Results.begin(), Results.end(), [](MoranResult a, MoranResult b) {
        return a.probability > b.probability;
    });

    return Results;
}

template <typename RngT>
MoranResult MoranTopD(const digraph &G, int k, int a, std::vector<double> &P, RngT &Rng) {
    MoranResult res;
    auto N = TopD(G, k);

    bitvector States(G.size());
    for (auto &n : N){
        States.set(n);
    }

    res.probability = MonteCarloEstimate(G, N, States, P, k, a, Rng);
    res.seed = N;
    return res;
}

template <typename RngT>
MoranResult MoranTopD(const digraph &G, int k, int a, edge_probs &E, RngT &Rng) {
    MoranResult res;
    auto N = TopD(G, k);

    bitvector States(G.size());
    for (auto &n : N){
        States.set(n);
    }

    res.probability = MonteCarloEstimate(G, N, States, E, k, a, Rng);
    res.seed = N;
    return res;
}

template <typename RngT>
MoranResult MoranBotP(const digraph &G, int k, int a, std::vector<double> &P, RngT &Rng) {
    MoranResult res;
    nodes N(G.size());
    std::iota(N.begin(), N.end(), node(0));
    std::sort(N.begin(), N.end(), [&](const auto a, const auto b){
        return P[a] < P[b];
    });
    N = nodes(N.begin(), N.begin()+k);

    bitvector States(G.size());
    for (auto &n : N){
        States.set(n);
    }

    res.probability = MonteCarloEstimate(G, N, States, P, k, a, Rng);
    res.seed = N;
    return res;
}

// Produce random edge-weights per a random (out-degree size) Dirichlet Distribution
template <typename RngT>
edge_probs InitMoranEdgeProbabilities(const digraph &G, RngT &Rng){
    edge_probs E = {};
    std::vector<int> Weights = {1, 10, 100, 1000};
    for (node u = 0; u < G.size(); ++u){
        auto size = G[u].Out.size();
        double sum = 0.0;
        int weight;
        std::vector<double> Gamma(size);
        std::uniform_int_distribution<> WeightDist(0, Weights.size()-1);
        for (size_t i=0; i<size; i++){
            weight = Weights[WeightDist(Rng)];
            std::gamma_distribution<double> GamDist(weight,1);
            auto n = GamDist(Rng);
            Gamma[i] = n;
            sum += n;
        }
        // Produce 1 extra and ignore it to make room for 'no node' choice (this value is then recovered later, which is admittedly a little silly)
        weight = Weights[WeightDist(Rng)];
        std::gamma_distribution<double> GamDist(weight,1);
        sum += GamDist(Rng);

        for (size_t i=0; i<size; i++){
            E.insert({edge{u, G[u].Out[i]}, Gamma[i]/sum});
        }
    }
    return E;
}

int main(int argc, char** argv){
    if(argc == 1){
        PrintUsage(argv[0],stdout);
        return 0;
    }

    try {
        command_line Cmd(argc,argv);

        digraph G;
        int k, acc, Retries;
        double p = 0.5, MutateRate;
        std::vector<double> P;
        edge_probs E;

        //std::random_device rd;
        //std::mt19937_64 Rng(rd());
        std::mt19937_64 Rng;

        G = ReadGraph(OpenIfstream(Cmd.ReadString()));

        if(Cmd.IsFlag('s')){
            std::vector<node> Nodes(Cmd.ReadUint<size_t>());
            std::iota(Nodes.begin(),Nodes.end(),node(0));

            InduceSubset(Nodes,&G,nullptr,nullptr,nullptr,nullptr);
        }

        if(Cmd.IsFlag('k')){
            // Is this okay?
            k = Cmd.ReadUint<size_t>();
            if (k >= G.size()){
                throw std::runtime_error("k is too big!");
            }
        }else{
            k = G.size()/2;
        }

        if (Cmd.IsFlag('P')) {
            P = ReadProbabilities(OpenIfstream(Cmd.ReadString()));
        } else if(Cmd.IsFlag('p')){
            p = Cmd.ReadDouble();
            P.assign(G.size(), p);
        } else {
            P.assign(G.size(), p);
        }
        if (Cmd.IsFlag('E')) {
            E = InitMoranEdgeProbabilities(G, Rng);
        }

        if(Cmd.IsFlag('a')){
            acc = Cmd.ReadUint<size_t>();
        }else{
            acc = 1000;
        }

        if(Cmd.IsFlag('r')){
            Retries = Cmd.ReadUint<size_t>();
        }else{
            Retries = 10;
        }

        if(Cmd.IsFlag('m')){
            MutateRate = Cmd.ReadDouble();
        }else{
            MutateRate = 0.2;
        }

        bool weighted = false;

        if(Cmd.IsFlag('w')){
            weighted = true;
        }

        bool deg_init = false, p_init = false;

        if(Cmd.IsString("degree")){
            deg_init = true;
            // No point
        } else if (Cmd.IsString("probability")){
            p_init = true;
        }

        bool graph = false;

        if(Cmd.IsFlag('g')){
            graph = true;
        }

        MoranResults results;

        if(Cmd.IsFlag('d')){
            double mean_degree, mean_p;
            MoranResult res;
            if (!E.empty()){
                res = MoranTopD(G, k, acc, E, Rng);

                MoranKStats(G, E, res, mean_degree, mean_p);
            } else {
                res = MoranTopD(G, k, acc, P, Rng);
                double mean_degree, mean_p;
                MoranKStats(G, P, res, mean_degree, mean_p);

            }

            std::cout << k << "," << MutateRate << "," << Retries << "," << "Top-Degree," << "Top-Degree," << "Top-Degree," << res.probability << "," << "Top-Degree," << mean_degree << "," << mean_p << "\n";
        }else if (Cmd.IsFlag('e')) {
            auto res = MoranBotP(G, k, acc, P, Rng);
            double mean_degree, mean_p;
            MoranKStats(G, P, res, mean_degree, mean_p);

            std::cout << k << "," << MutateRate << "," << Retries << "," << "Bot-P," << "Bot-P," << "Bot-P," << res.probability << "," << "Bot-P," << mean_degree << "," << mean_p << "\n";
        }else {
            if(deg_init) {
                results = MoranGenetic(G, P, Retries, acc, k, MutateRate, Rng, PFitnessIndex(), D_Init(), weighted, graph);
            } else if (p_init) {
                results = MoranGenetic(G, P, Retries, acc, k, MutateRate, Rng, DFitnessIndex(), P_Init(), weighted, graph);
            } else if (!E.empty()){
                results = MoranGenetic(G, E, Retries, acc, k, MutateRate, Rng, RandomIndex(), R_Init(), weighted, graph);
            } else {
                results = MoranGenetic(G, P, Retries, acc, k, MutateRate, Rng, RandomIndex(), R_Init(), weighted, graph);
            }


            // Janky
            for (auto i : results)
            {
                std::cout << k << ",";
                //std::cout << p << ",";
                std::cout << MutateRate << ",";
                std::cout << Retries << ",";
                if (deg_init){
                    std::cout << "top-degree,p-fitness,";
                } else if (p_init) {
                    std::cout << "bot-p,k-fitness,";
                } else {
                    std::cout << "none,none,";
                }
                //std::cout << (weighted ? "true" : "false") << ",";

                std::cout << i.probability << ",";
                std::cout << i.sim_time << ",";
                std::cout << i.mutations << ",";

                double mean_degree, mean_p;
                if (E.empty()){
                    MoranKStats(G, P, i, mean_degree, mean_p);
                } else {
                    MoranKStats(G, E, i, mean_degree, mean_p);
                }
                std::cout << mean_degree << ",";
                std::cout << mean_p << "\n";


                if (graph) {
                    // Warning! large files ahead! would suggest using small values for r and a.

                    std::ofstream fail_data, succ_data, summary;
                    fail_data.open("results/MC-data-fail.csv", std::ios::trunc);
                    fail_data << "#This contains the progression of a (failed) monte carlo simulation of the most recent run of the program (hope commenting works)\n";
                    fail_data << "iteration,progress\n";

                    succ_data.open("results/MC-data-success.csv", std::ios::trunc);
                    succ_data << "#This contains the progression of a (successful) monte carlo simulation of the most recent run of the program (hope commenting works)\n";
                    succ_data << "iteration,progress\n";

                    summary.open("results/MC-summary.csv", std::ios::trunc);
                    summary << "#This contains the summary of progression of monte carlo simulations of the most recent run of the program (hope commenting works)\n";
                    summary << "success,influence extrema\n";

                    int maxmax, minmin = INT_MAX;
                    for (auto &mc : i.monte){
                        if (!mc.success){
                            maxmax = std::max(mc.influence_extreme, maxmax);
                        } else {
                            minmin = std::min(mc.influence_extreme, minmin);
                        }
                    }

                    for (auto &mc : i.monte){
                        if (mc.influence_extreme == maxmax){
                            for (int j = 0; j < mc.time_lapse.size(); j++){
                                fail_data << mc.time_lapse[j].first << "," << mc.time_lapse[j].second << "\n";
                            }
                        } else if (mc.influence_extreme == minmin){
                            for (int j = 0; j < mc.time_lapse.size(); j++){
                                succ_data << mc.time_lapse[j].first << "," << mc.time_lapse[j].second << "\n";
                            }
                        }

                        summary << (mc.success ? "true" : "false") << "," << mc.influence_extreme << "\n";

                    }

                    summary.close();
                    fail_data.close();
                    succ_data.close();

                }
            }
        }

        return 0;
    }catch(command_line_exception& e){
        std::cerr << e.what() << '\n';

        PrintUsage(argv[0], stderr);
    }catch(std::exception& e){
        std::cerr << e.what() << '\n';
    }

    return 1;
}
