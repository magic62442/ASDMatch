#ifndef TWBOUND_HZY_H
#define TWBOUND_HZY_H

#include "hypergraph.h"
#include "SCsolver.h"

db FHW_ub(VertexSet V, std::vector <VertexSet> E, SCsolver & Sv, Order & o);
db FHW_lb(VertexSet V, std::vector <VertexSet> E, SCsolver & Sv);

#endif