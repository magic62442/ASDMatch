#ifndef TWBOUND_HZY_CPP
#define TWBOUND_HZY_CPP

#include "twbound.h"
#include "hypergraph.h"
#include "SCsolver.h"
#include "td_utils.hpp"
#include <queue>

size_t MinFHW(VertexSet S, vector <VertexSet> E, SCsolver & Sv) {
    vector <size_t> V;
    S.getelement(V);

    db res = inf; 
    size_t ans;

    for(size_t i = 0; i < V.size(); ++i) {
        VertexSet ei = Aggregation(V[i], E);
        ei.intersect(S);
        db tmp = Sv.solve(ei);
        if(tmp <= res) {
            res = tmp;
            ans = V[i];
        }
    }

    return ans;
}

db GreedyMinFHW(VertexSet S, vector <VertexSet> E, SCsolver & Sv, Order & o) {
    vector <size_t> V;
    S.getelement(V);

    if(V.size() < 2)
        return 1.;

    VertexSet Ne[MAXELENUM];
    db width[MAXELENUM];
    priority_queue < pair<db, size_t> > Q;

    for(size_t i = 0; i < V.size(); ++i) {
        size_t u = V[i];
        Ne[u] = Aggregation(u, E);
        Ne[u].intersect(S);
        width[u] = Sv.solve(Ne[u]);
        Q.push(make_pair(-width[u], u));
    }

    db res = 1.;
    // for(size_t i = 0; i < V.size() - 1; ++i) {
    //     size_t v = MinFHW(S, E, Sv);
    //     VertexSet ei = Aggregation(v, E);
    //     ei.intersect(S);
    //     res = max(res, Sv.solve(ei));
    //     ei.reset(v);
    //     S.reset(v);
    //     E.push_back(ei);
    // }

    for(size_t i = 0; i < V.size() - 1; ++i) {
        pair<db, size_t> tmp = Q.top(); Q.pop();
        while(!S.test(tmp.second) || -tmp.first != width[tmp.second]) {
            tmp = Q.top();
            Q.pop();
        }
        size_t u = tmp.second;
        o.push_back(u);
        res = max(res, -tmp.first);
        Ne[u].reset(u);
        S.reset(u);

        /*
            eliminate u
        */

        vector <size_t> Ve;
        Ne[u].getelement(Ve);
        for(size_t j = 0; j < Ve.size(); ++j) {
            size_t v = Ve[j];
            Ne[v] = Ne[v] | Ne[u];
            Ne[v].reset(u);
            db w = Sv.solve(Ne[v]);
            if(w != width[v])
                Q.push(make_pair(-w, v));
            width[v] = w;
        }
    }
    S.getelement(V);
    o.push_back(V[0]);

    return res;
}

db FHW_ub(VertexSet S, vector <VertexSet> E, SCsolver & Sv, Order & o) {
#ifdef UP_MINW
    return GreedyMinFHW(S, E, Sv, o);
#else
    return Sv.solve(S);
#endif
}   

db MMD_plus(VertexSet S, vector <VertexSet> E, SCsolver & Sv) {
    vector <size_t> V;
    S.getelement(V);

    if(V.size() < 2)
        return 1.;

    VertexSet Ne[MAXELENUM];
    double width[MAXELENUM];
    priority_queue < pair <db, size_t> > Q;
    for(size_t i = 0; i < V.size(); ++i) {
        size_t v = V[i];
        Ne[v] = Aggregation(v, E);
        Ne[v].intersect(S);
        width[v] = Sv.solve(Ne[v]);
        Q.push( make_pair(-width[v], v) );
    }

    db res = 0;
    for(size_t i = 0; i < V.size() - 1; ++i) {
        pair <db, size_t> tmp = Q.top();
        Q.pop();
        size_t u = tmp.second;
        while(!S.test(u) || -width[u] != tmp.first) { // find min-width vertex
            tmp = Q.top();
            u = tmp.second;
            Q.pop();
        }
        res = max(res, -tmp.first);
        S.reset(u);

        vector <size_t> Ve;
        Ne[u].reset(u);
        Ne[u].getelement(Ve);

        size_t v;
        // pair<db, size_t> min_w, max_w ;
        pair<size_t, size_t> least_c; // select a neighbour of (minimum width, maximum width, minmum common neighbours)
        // min_w = make_pair(inf, -1);
        // max_w = make_pair(0, -1);
        least_c = make_pair(MAXELENUM, -1);
        for(size_t j = 0; j < Ve.size(); ++j) {
            v = Ve[j];
            // if(Sv.solve(Ne[v]) < min_w.first)
            //     min_w = make_pair(Sv.solve(Ne[v]), v);
            // if(Sv.solve(Ne[v]) > max_w.first)
            //     max_w = make_pair(Sv.solve(Ne[v]), v);
            if((Ne[u] | Ne[v]).size() < least_c.first)
                least_c = make_pair((Ne[u] | Ne[v]).size(), v);
        }

        v = least_c.second; // choose lease_c strategy

        /*
            contract u -> v
        */

        Ne[v] = Ne[u] | Ne[v];

        for(size_t j = 0; j < E.size(); ++j) // update edge set
            if( (int) E[j].test(v) + E[j].test(u) == 1)
                E[j].Set(u), E[j].Set(v);
        
        Ne[v].getelement(Ve);
        Ne[v].reset(u);

        for(size_t j = 0; j < Ve.size(); ++j) { // recompute width
            size_t w = Ve[j];
            Ne[w].reset(u), Ne[w].Set(v);
            double tmp = FSC_LP(E, Ne[w]); // Edge set changed, could not use previous result
            if(tmp != width[w]) {
                width[w] = tmp;
                Q.push(make_pair(-tmp, w));
            }
        }
    }
    return res;
}

db MIN2FSC(VertexSet S, vector <VertexSet> E, SCsolver & Sv) {
    vector <size_t> V;
    S.getelement(V);

    if(V.size() < 2)
        return 1.;

    db minfsc1 = inf, minfsc2 = inf;

    // VertexSet Ne[V.size()];
    for(size_t i = 0; i < V.size(); ++i) {
        VertexSet ei = Aggregation(V[i], E);
        ei.intersect(S);
        db fsc = Sv.solve(ei);
        if(fsc < minfsc1) {
            minfsc2 = minfsc1;
            minfsc1 = fsc;
        }
        else if(fsc < minfsc2)
            minfsc2 = fsc;
    }
    return minfsc2;
}

db GammaR(VertexSet S, vector <VertexSet> E, SCsolver & Sv) {
    vector <size_t> V;
    S.getelement(V);
    size_t n = V.size();

    if(V.size() < 2)
        return 1.;

    pair<db, pair<size_t, VertexSet> > width_Ne[n];

    for(size_t i = 0; i < n; ++i) {
        VertexSet ei = Aggregation(V[i], E);
        ei.intersect(S);
        db fsc = Sv.solve(ei);
        width_Ne[i] = make_pair(fsc, make_pair(V[i], ei));
    }
    sort(width_Ne, width_Ne + n);

    size_t c[n], p[MAXELENUM];
    memset(c, 0, sizeof c);
    for(size_t i = 0; i < n; ++i)
        p[width_Ne[i].second.first] = i;

    size_t i = 0;
    while(i == c[i] && i < n - 1) {
        vector <size_t> tmp;
        width_Ne[i].second.second.getelement(tmp);
        for(auto it = tmp.begin(); it != tmp.end(); it++)
            c[p[*it]]++;
        i++;
    }
    return width_Ne[i].first;
}


db FHW_lb(VertexSet S, vector <VertexSet> E, SCsolver & Sv) {
#ifdef LOW_DELTA2
    return MIN2FSC(S, E, Sv);
#else
#ifdef LOW_GAMMAR
    return GammaR(S, E, Sv);
#else
#ifdef LOW_MMDP
    return MMD_plus(S, E, Sv);
#else
    return 1.;
#endif
#endif
#endif
}

#endif