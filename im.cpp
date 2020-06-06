#include "common.h"

#include <chrono>
#include <random>
#include <numeric>
#include <functional>

using namespace std::literals;

// TODO: Different moran algorithms
void PrintUsage(const char* Exe, std::FILE* File){
    std::fprintf(File,
        "Usage: %s <edge-list> <b-probs> (-b <q-scale>|<q-probs>) <features> (<seed>|-g <group>) (-r <subset-size>|-s <subset-size>) <base-truth> <normalizer> {<algorithms>}\n\n"

        "To specify algorithms you first need to specify a base truth spread estimator. The list of spread estimators is:\n"
        " - mc -i <repetitions=1000>             | Monte Carlo Simulation\n"
        " - simpath -t <pruning-threshold=0.001> | Simpath-Spread\n\n"

        "After that you need to specify a normalizer:\n"
        " - identity       | n_uv(p) = p\n"
        " - sigmoid        | n_uv(p) = tanh(p)\n"
        " - degree         | n_uv(p) = p/deg^in(v)\n"
        " - sigmoid_degree | n_uv(p) = tanh(p)/deg^in(v)\n\n"

        "After this you can specify a list of algorithms to run sequentially:\n"
        " - greedy -k <features> <spread-estimator> | Greedy selection algorithm using given spread estimator (see above).\n"
        " - top -k <features> <heuristic>           | Pick the best features according to the heuristic.\n\n"

        "Heuristics can be any of the spread estimators (see above) or one of the following:\n"
        " - nodes | Select the features associated with most nodes.\n"
        " - edges | Select the features associated with most nodes, weighted by in-edges.\n", Exe);
}

struct identity_normalize {
    double operator()(double x) const {
        return x;
    }

    double operator()(const digraph&, const edge_probs& P, edge e) const {
        return P.at(e);
    }

    double Inverse(double x) const {
        return x;
    }

    double Inverse(const digraph&, const edge_probs& P, edge e) const {
        return P.at(e);
    }
};

struct degree_normalize {
    double operator()(const digraph& G, const edge_probs& P, edge e) const {
        return P.at(e)/static_cast<double>(G[e.v].In.size());
    }

    double Inverse(const digraph& G, const edge_probs& P, edge e) const {
        return P.at(e)*static_cast<double>(G[e.v].In.size());
    }
};

struct sigmoid_normalize {
    double operator()(double x) const {
        return std::tanh(x);
    }

    double operator()(const digraph&, const edge_probs& P, edge e) const {
        return operator()(P.at(e));
    }

    double Inverse(double x) const {
        return std::atanh(x);
    }

    double Inverse(const digraph&, const edge_probs& P, edge e) const {
        return Inverse(P.at(e));
    }
};

struct sigmoid_degree_normalize {
    double operator()(const digraph& G, const edge_probs& P, edge e) const {
        return std::tanh(P.at(e))/static_cast<double>(G[e.v].In.size());
    }

    double Inverse(const digraph& G, const edge_probs& P, edge e) const {
        return std::atanh(P.at(e)*static_cast<double>(G[e.v].In.size()));
    }
};

template <typename RngT>
class linear_threshold_activator {
public:
    RngT Rng;

private:
    std::vector<double> Theta;

public:
    linear_threshold_activator(const digraph& G, RngT Rng):Rng(Rng), Theta(G.size()) {}

    template <typename NormalizeT>
    double Init(const digraph&, const edge_probs&, const nodes& S, NormalizeT){
        std::uniform_real_distribution<double> Dist(0.0, 1.0);
        for(auto& theta:Theta){
            theta = Dist(Rng);
        }

        for(auto& u:S){
            Theta[u] = -1.0;
        }

        return 1.0;
    }

    template <typename NormalizeT>
    bool operator()(const digraph& G, const edge_probs& P, node u, node v, NormalizeT Normalize){
        auto& theta = Theta[v];
        if(theta >= 0.0){
            theta -= Normalize(G, P, {u, v});

            if(theta < 0.0){
                return true;
            }
        }

        return false;
    }
};

template <typename ActivatesT>
class monte_carlo_spread {
public:
    std::size_t I;

    ActivatesT Activates;

private:
    const nodes* S;

    nodes Q;

public:
    monte_carlo_spread()=default;

    template <typename... Ts>
    monte_carlo_spread(const digraph& G, const nodes& S, std::size_t I, Ts&&... ts)
        :I(I), Activates(G, std::forward<Ts>(ts)...), S(&S), Q() {}

    // It is probably a mistake to ever copy this. Hence we only allow moves.
    monte_carlo_spread(monte_carlo_spread&& rhs)=default;
    monte_carlo_spread& operator=(monte_carlo_spread&& rhs)=default;

    template <typename NormalizeT>
    double operator()(const digraph& G, const edge_probs& P, NormalizeT Normalize){
        double Spread = 0.0;

        for(std::size_t i = 0; i < I; ++i){
            Q = *S;

            auto Prob = Activates.Init(G, P, *S, Normalize);

            Spread += Prob*static_cast<double>(S->size());

            while(!Q.empty()){
                auto u = Q.back();
                Q.pop_back();

                for(auto& v:G[u].Out){
                    if(Activates(G, P, u, v, Normalize)){
                        Q.push_back(v);

                        Spread += Prob;
                    }
                }
            }
        }

        return Spread/static_cast<double>(I);
    }
};

class simpath_spread {
public:
    double eta;

private:
    const nodes* S;

    struct stack_frame {
        node x;
        nodes::const_iterator it;
        double pp;
    };

    // These are simply stored to avoid having to always reallocate them per invocation.
    // They are always reset to their default state whenever `operator()` has finished.
    std::vector<stack_frame> Q;
    bitvector V;

public:
    simpath_spread()=default;

    simpath_spread(const digraph& G, const nodes& S, double eta):eta(eta), S(&S), Q(), V(G.size()){
        for(auto& u:S){
            V.set(u);
        }
    }

    // It is probably a mistake to ever copy this. Hence we only allow moves.
    simpath_spread(simpath_spread&& rhs)=default;
    simpath_spread& operator=(simpath_spread&& rhs)=default;

    template <typename NormalizeT>
    double operator()(const digraph& G, const edge_probs& P, NormalizeT Normalize){
        double Spread = 0.0;

        for(auto& u:*S){
            Spread += 1.0;

            Q.push_back({u, G[u].Out.begin(), 1.0});

            for(;;){
                auto& q = Q.back();
                auto x = q.x;

                if(q.it == G[x].Out.end()){
                    Q.pop_back();

                    if(Q.empty()){
                        break;
                    }

                    V.clear(x);

                    continue;
                }

                auto y = *q.it++;
                if(V[y]){
                    // $y \in Q or y \in D[x] or y \not\in W$. Either way skip it.
                    continue;
                }

                auto pp = q.pp*Normalize(G, P, {x, y});
                if(pp < eta){
                    continue;
                }

                V.set(y);

                Spread += pp;

                Q.push_back({y, G[y].Out.begin(), pp});
            }
        }

        return Spread;
    }
};

class nodes_spread {
public:
    double operator()(const digraph&, const edge_probs&, const nodes& Nodes){
        return static_cast<double>(Nodes.size());
    }
};

class in_edges_spread {
public:
    double operator()(const digraph& G, const edge_probs&, const nodes& Nodes){
        double Spread = 0.0;
        for(auto& v:Nodes){
            Spread += static_cast<double>(G[v].In.size());
        }

        return Spread;
    }
};

class weighted_in_edges_spread {
public:
    double operator()(const digraph& G, const edge_probs& Q, const nodes& Nodes){
        double Spread = 0.0;
        for(auto& v:Nodes){
            for(auto& u:G[v].In){
                // TODO: Normalize?
                Spread += Q.at({u, v});
            }
        }

        return Spread;
    }
};

class out_edges_spread {
public:
    double operator()(const digraph& G, const edge_probs&, const nodes& Nodes){
        double Spread = 0.0;
        for(auto& v:Nodes){
            Spread += static_cast<double>(G[v].Out.size());
        }

        return Spread;
    }
};

class weighted_out_edges_spread {
public:
    double operator()(const digraph& G, const edge_probs& Q, const nodes& Nodes){
        double Spread = 0.0;
        for(auto& v:Nodes){
            for(auto& u:G[v].Out){
                // TODO: Normalize?
                Spread += Q.at({u, v});
            }
        }

        return Spread;
    }
};

void IncreaseProbabilities(const digraph& G, edge_probs* Pf, const edge_probs& Q, const nodes& Nodes){
    for(auto& v:Nodes){
        for(auto& u:G[v].In){
            edge e = {u, v};
            Pf->at(e) += Q.at(e);
        }
    }
}

void IncreaseProbabilities(const digraph& G, edge_probs* Pf, const edge_probs& Q, const nodes& Nodes, edge_probs* P){
    for(auto& v:Nodes){
        for(auto& u:G[v].In){
            edge e = {u, v};
            P->at(e) = Pf->at(e) += Q.at(e);
        }
    }
}

void ResetProbabilities(const digraph& G, edge_probs* Pf, const edge_probs& P, const nodes& Nodes){
    for(auto& v:Nodes){
        for(auto& u:G[v].In){
            edge e = {u, v};
            Pf->at(e) = P.at(e);
        }
    }
}

template <typename EstimateSpreadT, typename NormalizeT>
std::size_t GreedyStep(const digraph& G, edge_probs* P, const edge_probs& Q, const groups& Groups, const features& Phi,
                       EstimateSpreadT EstimateSpread, NormalizeT Normalize){
    assert(!Phi.empty());

    auto Pf = *P;

    auto MaxSpread = -std::numeric_limits<double>::infinity();
    auto MaxIndex = std::numeric_limits<std::size_t>::max();

    for(std::size_t i = 0; i < Phi.size(); ++i){
        auto& Group = Groups.at(Phi[i]);

        IncreaseProbabilities(G, &Pf, Q, Group);

        auto Spread = EstimateSpread(G, Pf, Normalize);

        ResetProbabilities(G, &Pf, *P, Group);

        if(Spread > MaxSpread){
            MaxSpread = Spread;
            MaxIndex = i;
        }
    }

    return MaxIndex;
}

template <typename TimeT>
struct timed_feature {
    timed_feature()=default;
    timed_feature(feature f, TimeT t):f(f), t(t){}

    feature f;
    TimeT t;

    operator feature() const {
        return f;
    }
};

template <typename TimeT>
using timed_features = std::vector<timed_feature<TimeT>>;

template <typename TimeT, typename EstimateSpreadT, typename NormalizeT>
timed_features<TimeT> Greedy(const digraph& G, edge_probs* P, const edge_probs& Q, const groups& Groups, features Phi,
                      std::size_t K, EstimateSpreadT EstimateSpread, NormalizeT Normalize){
    if(K > Phi.size()){
        throw std::runtime_error("|Phi| < K.");
    }

    timed_features<TimeT> F;
    F.reserve(K);

    auto Pf = *P;

    while(F.size() < K){
        auto i = GreedyStep(G, P, Q, Groups, Phi, EstimateSpread, Normalize);

        auto& f = Phi[i];

        IncreaseProbabilities(G, &Pf, Q, Groups.at(f), P);

        F.push_back({f, TimeT::clock::now()});

        std::swap(f, Phi.back());
        Phi.pop_back();
    }

    return F;
}

template <typename EstimateSpreadT, typename NormalizeT, typename TimeT>
timed_features<TimeT> Greedy(const digraph& G, edge_probs* P, const edge_probs& Q, const groups& Groups, features Phi,
                             TimeT MinTime, EstimateSpreadT EstimateSpread, NormalizeT Normalize){
    if(Phi.empty()){
        throw std::runtime_error("|Phi| == 0.");
    }

    timed_features<TimeT> F = {};

    auto Pf = *P;

    TimeT t;

    do{
        auto i = GreedyStep(G, P, Q, Groups, Phi, EstimateSpread, Normalize);

        auto& f = Phi[i];

        IncreaseProbabilities(G, &Pf, Q, Groups.at(f), P);

        t = TimeT::clock::now();

        F.push_back({f, t});

        std::swap(f, Phi.back());
        Phi.pop_back();
    }while(!Phi.empty() && t < MinTime);

    return F;
}

struct feature_spread {
    feature_spread()=default;
    feature_spread(feature f, double Spread = 0.0):f(f), Spread(Spread){}

    feature f;
    double Spread;

    friend bool operator<(feature_spread lhs, feature_spread rhs){
        // We want to sort in increasing order.
        return lhs.Spread > rhs.Spread;
    }

    operator feature() const {
        return f;
    }
};

using feature_spreads = std::vector<feature_spread>;

template <typename EstimateSpreadT, typename NormalizeT>
feature_spreads EstimateEach(const digraph& G, const edge_probs& P, const edge_probs& Q, const groups& Groups,
                             const features& Phi, EstimateSpreadT EstimateSpread, NormalizeT Normalize){
    feature_spreads F(Phi.begin(), Phi.end());

    auto Pf = P;
    for(auto& f:F){
        auto& Group = Groups.at(f.f);

        IncreaseProbabilities(G, &Pf, Q, Group);

        f.Spread = EstimateSpread(G, Pf, Normalize);

        ResetProbabilities(G, &Pf, P, Group);
    }

    return F;
}

template <typename EstimateSpreadT>
feature_spreads EstimateEach(const digraph& G, const edge_probs& Q, const groups& Groups,
                             const features& Phi, EstimateSpreadT EstimateSpread){
    feature_spreads F(Phi.begin(), Phi.end());

    for(auto& f:F){
        auto& Group = Groups.at(f.f);

        f.Spread = EstimateSpread(G, Q, Group);
    }

    return F;
}

template <typename TimeT, typename... Ts, typename EstimateSpreadT, typename NormalizeT>
timed_features<TimeT> Top(const digraph& G, const edge_probs& P, const edge_probs& Q, const groups& Groups,
                          const features& Phi, std::size_t K, EstimateSpreadT EstimateSpread, NormalizeT Normalize){
    auto F = EstimateEach(G, P, Q, Groups, Phi, EstimateSpread, Normalize);
    // TODO: Use `std::partial_sort`?
    std::sort(F.begin(), F.end());

    auto t = TimeT::clock::now();

    timed_features<TimeT> r;
    r.reserve(K);

    for(std::size_t i = 0; i < K; ++i){
        r.push_back({F[i], t});
    }

    return r;
}

template <typename TimeT, typename... Ts, typename EstimateSpreadT>
timed_features<TimeT> Top(const digraph& G, const edge_probs& Q, const groups& Groups, const features& Phi,
                          std::size_t K, EstimateSpreadT EstimateSpread){
    auto F = EstimateEach(G, Q, Groups, Phi, EstimateSpread);
    std::sort(F.begin(), F.end());

    auto t = TimeT::clock::now();

    timed_features<TimeT> r;
    r.reserve(K);

    for(std::size_t i = 0; i < K; ++i){
        r.push_back({F[i], t});
    }

    return r;
}

template <typename NormalizeT>
bool ValidateProbabilities(const digraph& G, const edge_probs& P, NormalizeT Normalize){
    auto Valid = true;
    for(node v = 0; v < G.size(); ++v){
        double p = 0.0;

        for(auto& u:G[v].In){
            p += Normalize(G, P, {u, v});
        }

        if(p > 1.0){
            std::cerr << "p[" << v << "] = " << p << " > 1.\n";
            Valid = false;
        }
    }

    return Valid;
}

struct result {
    feature Feature;
    double Spread;
    std::chrono::duration<double> Time;
};

struct results_column {
    std::string Name;
    std::vector<result> Results;
};

using results = std::vector<results_column>;

template <typename ComputeSpreadT, typename NormalizeT, typename TimeT>
void AddResults(results* Results, const digraph& G, edge_probs* P, const edge_probs& Q, const groups& Groups,
                const timed_features<TimeT>& F, ComputeSpreadT ComputeSpread, NormalizeT Normalize, TimeT t,
                std::string Name){
    Results->emplace_back();
    auto& Column = Results->back();

    Column.Name = std::move(Name);

    Column.Results.reserve(F.size()+1);

    Column.Results.push_back({0, ComputeSpread(G, *P, Normalize), 0.0s});

    for(auto& f:F){
        IncreaseProbabilities(G, P, Q, Groups.at(f));

        Column.Results.push_back({f.f, ComputeSpread(G, *P, Normalize), f.t-t});
    }
}

void PrintResults(const results& Results){
    bool First = true;
    for(auto& Column:Results){
        if(First){
            First = false;
        }else{
            std::cout << ',';
        }

        std::cout << Column.Name << "_features," << Column.Name << "_spread," << Column.Name << "_time";
    }

    std::cout << '\n';

    auto Rows = std::max_element(Results.begin(), Results.end(), [](auto& lhs, auto& rhs){
        return lhs.Results.size() < rhs.Results.size();
    })->Results.size();

    for(std::size_t i = 0; i < Rows; ++i){
        First = true;
        for(auto& Column:Results){
            if(First){
                First = false;
            }else{
                std::cout << ',';
            }

            if(i >= Column.Results.size()){
                std::cout << ",,";
                continue;
            }

            auto& Result = Column.Results[i];

            std::cout << Result.Feature << "," << Result.Spread << "," << Result.Time.count();
        }

        std::cout << '\n';
    }
}

template <typename ComputeSpreadT, typename NormalizeT, typename TimeT>
void FinalizeResults(results* Results, const digraph& G, edge_probs P, const edge_probs& Q, const groups& Groups,
                     const timed_features<TimeT>& F, ComputeSpreadT ComputeSpread, NormalizeT Normalize, TimeT t,
                     std::string Name){
    AddResults(Results, G, &P, Q, Groups, F, ComputeSpread, Normalize, t, std::move(Name));

    ValidateProbabilities(G, P, Normalize);
}

template <typename CallbackT>
bool ReadSpreadHeuristic(command_line* Cmd, CallbackT Callback, std::string* Name){
    if(Cmd->IsString("nodes")){
        if(Name){
            *Name += "_nodes";
        }

        Callback(nodes_spread());
    }else if(Cmd->IsString("edges")){
        if(Name){
            *Name += "_edges";
        }

        Callback(in_edges_spread());
    }else{
        return false;
    }

    return true;
}

template <typename CallbackT, typename RngT>
bool ReadSpreadEstimator(command_line* Cmd, const digraph& G, const nodes& S, RngT& Rng, CallbackT Callback,
                         std::string* Name){
    if(Cmd->IsString("simpath")){
        if(Name){
            *Name += "_simpath";
        }

        double eta = 0.001;

        if(Cmd->IsFlag('t')){
            eta = Cmd->ReadDouble();

            if(Name){
                *Name += "_t";
                *Name += Cmd->PreviousString();
            }
        }

        simpath_spread EstimateSpread(G, S, eta);

        Callback(std::ref(EstimateSpread));
    }else if(Cmd->IsString("mc")){
        if(Name){
            *Name += "_mc";
        }

        std::size_t I = 1000;

        if(Cmd->IsFlag('i')){
            I = Cmd->ReadUint<std::size_t>();

            if(Name){
                *Name += "_i";
                *Name += Cmd->PreviousString();
            }
        }

        monte_carlo_spread<linear_threshold_activator<std::mt19937_64&>> EstimateSpread(G, S, I, Rng);

        Callback(std::ref(EstimateSpread));
    }else{
        return false;
    }

    return true;
}

template <typename CallbackT>
bool ReadNormalizer(command_line* Cmd, CallbackT Callback){
    if(Cmd->IsString("identity")){
        Callback(identity_normalize());
    }else if(Cmd->IsString("degree")){
        Callback(degree_normalize());
    }else if(Cmd->IsString("sigmoid")){
        Callback(sigmoid_normalize());
    }else if(Cmd->IsString("sigmoid_degree")){
        Callback(sigmoid_degree_normalize());
    }else{
        return false;
    }

    return true;
}

template <typename ClockT, typename RngT, typename ComputeSpreadT, typename NormalizeT>
void RunTests(command_line* Cmd, results* Results, const digraph& G, const edge_probs& B, const edge_probs& Q,
              const groups& Groups, const features& Phi, const nodes& S, RngT& Rng, ComputeSpreadT ComputeSpread,
              NormalizeT Normalize){
    std::string Name;
    while(!Cmd->IsDone()){
        if(Cmd->IsString("greedy")){
            Name = "greedy";

            if(Cmd->IsFlag('k')){
                auto K = Cmd->ReadUint<std::size_t>();

                if(!ReadSpreadEstimator(Cmd, G, S, Rng, [&](auto EstimateSpread){
                    auto P = B;

                    auto t = ClockT::now();

                    auto F = Greedy<decltype(t)>(G, &P, Q, Groups, Phi, K, EstimateSpread, Normalize);

                    FinalizeResults(Results, G, B, Q, Groups, F, ComputeSpread, Normalize, t, Name);
                }, &Name)){
                    Cmd->InvalidArgument();
                }
            }else if(Cmd->IsFlag('t')){
                char* Suffix;
                auto T = std::chrono::duration<double>(Cmd->ReadDouble(&Suffix));

                if(std::strcmp(Suffix,"d") == 0){
                    T *= 60.0f*60.0f*24.0f;
                }else if(std::strcmp(Suffix,"h") == 0){
                    T *= 60.0f*60.0f;
                }else if(std::strcmp(Suffix,"m") == 0){
                    T *= 60.0f;
                }else if(std::strcmp(Suffix,"s") != 0){
                    throw command_line_exception(Cmd->GetIndex()-1, Suffix);
                }

                if(!ReadSpreadEstimator(Cmd, G, S, Rng, [&](auto EstimateSpread){
                    auto P = B;

                    auto t = ClockT::now();
                    auto EndTime = std::chrono::time_point_cast<typename ClockT::duration>(t+T);

                    auto F = Greedy(G, &P, Q, Groups, Phi, EndTime, EstimateSpread, Normalize);

                    FinalizeResults(Results, G, B, Q, Groups, F, ComputeSpread, Normalize, t, Name);
                }, &Name)){
                    Cmd->InvalidArgument();
                }
            }else{
                Cmd->InvalidArgument();
            }
        }else if(Cmd->IsString("top")){
            Name = "top";

            if(Cmd->IsFlag('k')){
                auto K = Cmd->ReadUint<std::size_t>();

                if(!ReadSpreadEstimator(Cmd, G, S, Rng, [&](auto EstimateSpread){
                    auto t = ClockT::now();

                    auto F = Top<typename ClockT::time_point>(G, B, Q, Groups, Phi, K, EstimateSpread, Normalize);

                    FinalizeResults(Results, G, B, Q, Groups, F, ComputeSpread, Normalize, t, Name);
                }, &Name) && !ReadSpreadHeuristic(Cmd, [&](auto EstimateSpread){
                    auto t = ClockT::now();

                    auto F = Top<typename ClockT::time_point>(G, Q, Groups, Phi, K, EstimateSpread);

                    FinalizeResults(Results, G, B, Q, Groups, F, ComputeSpread, Normalize, t, Name);
                }, &Name)){
                    Cmd->InvalidArgument();
                }
            }else{
                Cmd->InvalidArgument();
            }
        }else{
            Cmd->InvalidArgument();
        }
    }
}

template <typename ClockT, typename RngT, typename ComputeSpreadT>
void RunTests(command_line* Cmd, results* Results, const digraph& G, const edge_probs& B, const edge_probs& Q,
              const groups& Groups, const features& Phi, const nodes& S, RngT& Rng, ComputeSpreadT ComputeSpread){
    if(!ReadNormalizer(Cmd, [&](auto Normalize){
        RunTests<ClockT>(Cmd, Results, G, B, Q, Groups, Phi, S, Rng, ComputeSpread, Normalize);
    })){
        Cmd->InvalidArgument();
    }
}

template <typename ClockT, typename RngT>
void RunTests(command_line* Cmd, results* Results, const digraph& G, const edge_probs& B, const edge_probs& Q,
              const groups& Groups, const features& Phi, const nodes& S, RngT& Rng){
    if(!ReadSpreadEstimator(Cmd, G, S, Rng, [&](auto ComputeSpread){
        RunTests<ClockT>(Cmd, Results, G, B, Q, Groups, Phi, S, Rng, ComputeSpread);
    }, nullptr)){
        Cmd->InvalidArgument();
    }
}

int main(int argc, char** argv){
    if(argc == 1){
        PrintUsage(argv[0],stdout);
        return 0;
    }

    try {
        command_line Cmd(argc,argv);

        digraph G;
        edge_probs B, Q;
        node_features Nf;
        groups Groups;
        features Phi;
        nodes S;

        std::mt19937_64 Rng;

        G = ReadGraph(OpenIfstream(Cmd.ReadString()));

        B = ReadProbabilities(OpenIfstream(Cmd.ReadString()), G);

        if(Cmd.IsFlag('b')){
            Q = B;

            auto Scale = Cmd.ReadDouble();
            for(auto& q:Q){
                q.second *= Scale;
            }
        }else{
            Q = ReadProbabilities(OpenIfstream(Cmd.ReadString()), G);
        }

        if(Cmd.IsFlag('s')){
            auto Scale = Cmd.ReadDouble();
            for(auto& b:B){
                b.second *= Scale;
            }
        }

        Nf = ReadFeatures(OpenIfstream(Cmd.ReadString()), G);

        Groups = ExtractGroups(Nf);

        Phi = ExtractPhi(Groups);

        if(Cmd.IsFlag('g')){
            auto SeedGroup = Cmd.ReadUint<feature>();

            S = Groups.at(SeedGroup);
        }else{
            S = ReadSeeds(OpenIfstream(Cmd.ReadString()), G);
        }

        if(Cmd.IsFlag('r')){
            auto Limit = Cmd.ReadUint<std::size_t>();

            if(Limit < G.size()){
                std::vector<node> Nodes(G.size());

                {
                    std::iota(Nodes.begin(), Nodes.end(), node(0));

                    std::shuffle(Nodes.begin(), Nodes.end(), Rng);

                    Nodes.resize(Limit);
                }

                InduceSubset(Nodes, &G, &B, &Q, &Groups, &S);
            }
        }else if(Cmd.IsFlag('s')){
            auto Limit = Cmd.ReadUint<std::size_t>();

            auto Nodes = S;
            bitvector NodeSet(G.size());
            for(auto& Node:Nodes){
                NodeSet.set(Node);
            }

            for(std::size_t i = 0;i < Nodes.size() && Nodes.size() < Limit;++i){
                for(auto& v:G[Nodes[i]].Out){
                    if(!NodeSet[v]){
                        NodeSet.set(v);
                        Nodes.push_back(v);

                        if(Nodes.size() >= Limit){
                            break;
                        }
                    }
                }
            }

            InduceSubset(Nodes, &G, &B, &Q, &Groups, &S);
        }

        using my_clock = std::chrono::steady_clock;

        results Results = {};

        RunTests<std::chrono::steady_clock>(&Cmd, &Results, G, B, Q, Groups, Phi, S, Rng);

        PrintResults(Results);

        return 0;
    }catch(command_line_exception& e){
        std::cerr << e.what() << '\n';

        PrintUsage(argv[0], stderr);
    }catch(std::exception& e){
        std::cerr << e.what() << '\n';
    }

    return 1;
}
