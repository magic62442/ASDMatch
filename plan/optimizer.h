//
// Created by Qiyan LI on 2024/5/22.
//

#ifndef IN_MEMORY_JOIN_OPTIMIZER_H
#define IN_MEMORY_JOIN_OPTIMIZER_H

#include "subset_structure.h"

void
simplePlan(const Graph &query, HyperTree &t, CandidateSpace &cs, PrefixNode *&pt,
           const std::vector<SubsetStructure> &dpStructures, bool *visited, VertexID *partMatch,
           VertexID **candidates, ui *candCount, double &minCost, PrefixNode *bestPT, bool share);
void
simplePlan(const Graph &query, HyperTree &t, CandidateSpace &cs, PrefixNode *&pt,
           double &minCost);
void setTDExtention(HyperTree &t, const Graph &query);
void bagPermutationDP(const Graph &query, HyperTree &t, CandidateSpace &cs, PrefixNode *&pt,
                      vector<SubsetStructure> &dpStructures, double &bestCost, int reorder = 1, bool connected = true,
                      bool globalOrderType = false);
void optCostPlan(const Graph &query, HyperTree &t, CandidateSpace &cs, PrefixNode *&pt,
                 vector<SubsetStructure> &dpStructures, double &bestCost, int reorder = 1, bool connected = true,
                 bool globalOrderType = false);
void reorderBags(HyperTree &t, PrefixNode *pt, const std::vector<SubsetStructure> &dpStructures, int heuristic,
                 const Graph &query, CandidateSpace &cs);
#endif //IN_MEMORY_JOIN_OPTIMIZER_H
