#ifndef PREPROCESSING_HZY_H
#define PREPROCESSING_HZY_H

#include "hypergraph.h"
using namespace std;

bool InducedHyperG(HyperG& H, vector<size_t>& remain_vertex, map<size_t, size_t>& Vres_map);
bool DelDegree1Vertex(HyperG& H, Order& prefix_o, map<size_t, size_t>& Vres_map);
bool DelCoveredEdge(HyperG& H);
bool DelISOVertex(HyperG& H, Order& prefix_o, map<size_t, size_t>& Vres_map);
void Preprocessing(HyperG& H, Order& prefix_o, map<size_t, size_t>& Vres_map);

#endif // PREPROCESSING_HZY_H
