#ifndef SCSOLVER_HZY_H
#define SCSOLVER_HZY_H

#include "hypergraph.h"

struct SCsolver {
    std::vector < VertexSet > E;
    std::map <VertexSet, double> res;

    SCsolver() {}
    SCsolver(std::vector < VertexSet > & _E) {
        E = _E;
    }

    db solve(VertexSet S);  
};

double calc_width(HyperG & H, Order o, std::vector <VertexSet> E, double * w);
double fractional_set_cover(size_t n, std::vector< VertexSet > & E, VertexSet e);
double FSC_LP(std::vector< VertexSet > & E, VertexSet e);

#endif