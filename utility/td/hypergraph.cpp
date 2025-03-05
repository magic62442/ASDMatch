#ifndef HYPERGRAPH_HZY_CPP
#define HYPERGRAPH_HZY_CPP
#include "hypergraph.h"

using namespace std;
HyperG::HyperG(size_t N_, vector <set<string> > & v, map <string, size_t> & f) : N(N_) {
    /*
            construct Hypergraph via a group of string set with number mapping given
            N_: number of element
            v: string set vector
            f: number mapping
    */
    M = v.size();
    e.clear();

    for(size_t i = 0; i < v.size(); ++i) {
        VertexSet temp;
        for(auto it = v[i].begin(); it != v[i].end(); ++it) 
            temp.Set(f[*it]);
        e.push_back(temp);    
    }
}

HyperG::HyperG(VertexSet S, vector<VertexSet> & v) {
    N = S.size();
    M = v.size();
    e.clear();

    size_t numbered[MAXBITNUMBER];
    vector <size_t> V;
    S.getelement(V);

    for(size_t i = 0; i < V.size(); ++i)
        numbered[V[i]] = i;

    for(size_t i = 0; i < v.size(); ++i) {
        V.clear();
        v[i].getelement(V);
        VertexSet tmp;
        for(size_t j = 0; j < V.size(); ++j)
            tmp.Set(numbered[V[j]]);
        e.push_back(tmp);
    }
}

vector <VertexSet> & HyperG::PrimalG() {
    if(__PrimalG.size() == N) {
        // the primal graph has been initialized
        return __PrimalG;
    }

    __PrimalG.clear();
    for(size_t i = 0; i < N; ++i) {
        VertexSet temp;
        for(size_t j = 0; j < e.size(); ++j)
            if(e[j].test(i)) {
                temp.merge(e[j]);
            }
        __PrimalG.push_back(temp);
    }
    return __PrimalG;
}

HyperG HyperG::induced(VertexSet S) {
    vector <VertexSet> e_;
    VertexSet S_;
    for(size_t i = 0; i < e.size(); ++i) {
        VertexSet tmp = e[i] & S;
            if(tmp.size() >= 2) {
                e_.push_back(tmp);
                S_.merge(tmp);
            }
    }
    S = S_;

    return HyperG(S, e_);    
}

#endif