#ifndef COMMON_H_INCLUDED
    #define COMMON_H_INCLUDED

#include <vector>
#include <string>
#include <unordered_map>

#include <algorithm>

#include <fstream>
#include <sstream>
#include <iostream>

#include <cassert>
#include <cstddef>
#include <cstring>

// Similar to `std::bitset`, but dynamically sized.
class bitvector {
public:
    using word_type = std::size_t;
    static constexpr std::size_t word_size = CHAR_BIT*sizeof(word_type);

    bitvector():Words(){}
    explicit bitvector(std::size_t Size):Words((Size+word_size-1)/word_size), Size(Size) {}

    bool operator[](std::size_t i) const {
        return (Words[i/word_size]&Mask(i)) != 0;
    }

    void set(std::size_t i, bool b = true){
        auto& w = Words[i/word_size];
        auto  m = Mask(i);

        if(b){
            w |=  m;
        }else{
            w &= ~m;
        }
    }
    void clear(std::size_t i){
        set(i, false);
    }
    void flip(std::size_t i){
        Words[i/word_size] ^= Mask(i);
    }

    void clear(){
        for(auto& Word:Words){
            Word = 0;
        }
    }

    std::size_t size() const noexcept {
        return Size;
    }

private:
    std::vector<word_type> Words;

    std::size_t Size;

    static word_type Mask(std::size_t i){
        return word_type(1) << i%word_size;
    }
};

using node = std::size_t;
using feature = std::size_t;

struct edge {
    node u, v;

    friend bool operator==(edge lhs, edge rhs) noexcept {
        return lhs.u == rhs.u && lhs.v == rhs.v;
    }
    friend bool operator!=(edge lhs, edge rhs) noexcept {
        return lhs.u != rhs.u || lhs.v != rhs.v;
    }
};

template <std::size_t>
struct fnv_constants;

template<>
struct fnv_constants<32> {
    static constexpr std::uint32_t init  = 0x811C9DC5;
    static constexpr std::uint32_t prime = 0x1000193;
};

template<>
struct fnv_constants<64> {
    static constexpr std::uint64_t init  = 0xCBF29CE484222325;
    static constexpr std::uint64_t prime = 0x100000001B3;
};

template<>
struct std::hash<edge> {
    std::size_t operator()(edge e) const noexcept {
        using constants = fnv_constants<CHAR_BIT*sizeof(std::size_t)>;

        // FNV-1a hash
        std::size_t r = constants::init;
        for(auto it = reinterpret_cast<unsigned char*>(&e), End = reinterpret_cast<unsigned char*>(&e+1); it != End; ++it){
            r ^= *it;
            r *= constants::prime;
        }

        return r;
    }
};

using nodes = std::vector<node>;

struct edge_list {
    nodes In;
    nodes Out;
};

using digraph = std::vector<edge_list>;
using edge_probs = std::unordered_map<edge, double>;

using features = std::vector<feature>;
using node_features = std::vector<features>;
// Features can have many values, but as far as I can tell there is no reason for them to not just have linearly
// increasing values like nodes. If they have linearly increasing values we could just use a vector.
// TODO: Fix this?
using groups = std::unordered_map<feature, nodes>;

using influence = std::vector<double>;

inline void Require(std::istream& In){
    if(!In){
        throw std::runtime_error("IO error.");
    }
}

template <typename T>
inline void Grow(std::vector<T>* v, std::size_t MaxId){
    auto Size = MaxId+1;
    if(Size == 0){
        throw std::runtime_error("Overflow when growing vector.");
    }

    if(Size > v->size()){
        v->resize(Size);
    }
}

template <typename T>
inline bool InitializeSet(std::vector<T>* v){
    std::sort(v->begin(), v->end());

    // TODO: Duplicate entries are apparently a thing. Should we allow that?
    auto NewEnd = std::unique(v->begin(), v->end());
    if(NewEnd == v->end()){
        return true;
    }else{
        v->erase(NewEnd, v->end());

        return false;
    }
}

inline std::ifstream OpenIfstream(const char* File){
    std::ifstream r(File);

    if(!r){
        throw std::runtime_error("Unable to open file.");
    }

    return r;
}

inline digraph ReadGraph(std::istream& In){
    digraph G = {};

    node u, v;
    while(In >> u){
        Require(In >> v);

        Grow(&G, std::max(u, v));

        G[u].Out.push_back(v);
        G[v].In.push_back(u);
    }

    for(auto& E:G){
        InitializeSet(&E.In);
        InitializeSet(&E.Out);

        if(E.In.empty() && E.Out.empty()){
            throw std::runtime_error("Node has no edges.");
        }
    }

    return G;
}

inline digraph ReadGraph(std::istream&& In){
    return ReadGraph(In);
}

inline edge_probs InitProbabilities(const digraph& G){
    edge_probs B = {};

    for(node u = 0; u < G.size(); ++u){
        for(auto& v:G[u].Out){
            B.emplace(edge{u, v}, 0.0);
        }
    }

    return B;
}

inline std::vector<double> ReadProbabilities(std::istream&& In){
    double p;
    std::vector<double> P;
    while(In >> p) {
        P.push_back(p);
    }

    return P;
}

inline edge_probs ReadProbabilities(std::istream& In, const digraph& G){
    auto B = InitProbabilities(G);

    edge e;
    while(In >> e.u){
        Require(In >> e.v);
        Require(In >> B.at(e));
    }

    return B;
}

inline edge_probs ReadProbabilities(std::istream&& In, const digraph& G){
    return ReadProbabilities(In, G);
}

inline node_features ReadFeatures(std::istream& In, const digraph& G){
    node_features Nf(G.size());

    std::string LineStr;
    while(std::getline(In, LineStr)){
        std::stringstream Line(LineStr);
        node u;
        if(!(Line >> u)){
            continue;
        }

        if(u >= G.size()){
            throw std::runtime_error("Invalid node when parsing features.");
        }

        feature f;

        auto& F = Nf[u];
        while(Line >> f){
            F.push_back(f);
        }
    }

    for(auto& F:Nf){
        if(!InitializeSet(&F)){
            throw std::runtime_error("Duplicate features.");
        }
    }

    return Nf;
}

inline node_features ReadFeatures(std::istream&& In, const digraph& G){
    return ReadFeatures(In, G);
}

inline groups ExtractGroups(const node_features& Nf){
    groups Groups = {};

    for(node u = 0; u < Nf.size(); ++u){
        auto& F = Nf[u];
        for(auto& f:F){
            Groups[f].push_back(u);
        }
    }

    return Groups;
}

inline features ExtractPhi(groups& Groups){
    features Phi = {};

    for(auto& F:Groups){
        Phi.push_back(F.first);
    }

    return Phi;
}

inline nodes ReadSeeds(std::istream& In, const digraph& G){
    nodes S = {};

    node u;
    while(In >> u){
        if(u >= G.size()){
            throw std::runtime_error("Invalid node index when reading seed.");
        }

        S.push_back(u);
    }

    if(!InitializeSet(&S)){
        throw std::runtime_error("Duplicate features.");
    }

    return S;
}

inline nodes ReadSeeds(std::istream&& In, const digraph& G){
    return ReadSeeds(In, G);
}

inline void InduceSubset(const nodes& Nodes, digraph* G, edge_probs* B, edge_probs* Q, groups* Groups, nodes* S){
    std::unordered_map<node, node> NodeMap;

    {
        digraph Gp = {};
        Gp.reserve(Nodes.size());

        for(node u = 0; u < Nodes.size(); ++u){
            NodeMap[Nodes[u]] = u;
            Gp.push_back(std::move((*G)[Nodes[u]]));
        }

        *G = std::move(Gp);
    }

    auto FixNodes = [&NodeMap](nodes* Vec){
        auto it = Vec->begin();
        for(auto jt = Vec->begin(); jt != Vec->end(); ++jt){
            auto u = NodeMap.find(*jt);
            if(u == NodeMap.end()){
                continue;
            }

            *it++ = u->second;
        }

        Vec->erase(it, Vec->end());

        std::sort(Vec->begin(),Vec->end());
    };

    for(auto& E:*G){
        FixNodes(&E.In);
        FixNodes(&E.Out);
    }

    if(Groups){
        for(auto& Group:*Groups){
            FixNodes(&Group.second);
        }
    }

    if(S){
        FixNodes(S);
    }

    auto FixEdges = [&NodeMap](const edge_probs& P){
        edge_probs Pp = {};
        for(auto& p:P){
            auto u = NodeMap.find(p.first.u);
            auto v = NodeMap.find(p.first.v);
            if(u == NodeMap.end() || v == NodeMap.end()){
                continue;
            }

            Pp[{u->second, v->second}] = p.second;
        }

        return Pp;
    };

    if(B){
        *B = FixEdges(*B);
    }

    if(Q){
        *Q = FixEdges(*Q);
    }
}

class command_line_exception:public std::runtime_error {
public:
    command_line_exception(int i, const char* Arg)
        :std::runtime_error("Invalid command line at argument "+std::to_string(i)+": "+Arg) {}

    command_line_exception(int i)
        :std::runtime_error("Expected more than "+std::to_string(i-1)+" command line arguments.") {}
};

class command_line {
public:
    command_line(int argc, char** argv):argc(argc), argv(argv), i(1) {
        assert(i <= argc);
    }

    bool IsDone() const noexcept {
        return i == argc;
    }

    int GetIndex() const noexcept {
        return i;
    }

    const char* PreviousString() const noexcept {
        return argv[i-1];
    }

    const char* ReadString(){
        Validate();
        return argv[i++];
    }

    bool IsString(const char* Str){
        if(IsDone()){
            return false;
        }

        Validate();

        if(std::strcmp(argv[i], Str) == 0){
            ++i;
            return true;
        }else{
            return false;
        }
    }

    bool IsFlag(char c){
        if(IsDone()){
            return false;
        }

        auto Arg = argv[i];
        if(Arg[0] == '-' && Arg[1] == c && Arg[2] == '\0'){
            ++i;
            return true;
        }else{
            return false;
        }
    }

    double ReadDouble(char** EndPtr){
        Validate();

        errno = 0;
        auto r = std::strtod(argv[i], EndPtr);

        if(errno != 0){
            throw std::system_error(errno, std::generic_category());
        }

        ++i;

        return r;
    }

    double ReadDouble(){
        Validate();

        char* EndPtr;
        auto r = ReadDouble(&EndPtr);

        if(*EndPtr != '\0'){
            --i;
            InvalidArgument();
        }

        return r;
    }

    void InvalidArgument(){
        throw command_line_exception(i, argv[i]);
    }

    template <typename T>
    T ReadUint(){
        static_assert(std::is_unsigned<T>::value,"`T` must be an unsigned integer type.");

        Validate();

        char* EndPtr;
        errno = 0;
        auto r = std::strtoull(argv[i], &EndPtr, 10);

        if(errno != 0){
            throw std::system_error(errno, std::generic_category());
        }

        if(*EndPtr != '\0'){
            throw command_line_exception(i, argv[i]);
        }

        if(static_cast<T>(r) != r){
            throw std::runtime_error("Seed group is too large.");
        }

        ++i;

        return static_cast<T>(r);
    }

private:
    int argc;
    char** argv;
    int i;

    void Validate(){
        if(IsDone()){
            throw command_line_exception(i);
        }
    }
};

#endif // COMMON_H_INCLUDED
