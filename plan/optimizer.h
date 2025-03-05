//
// Created by anonymous authors on 2024/5/22.
//

#ifndef IN_MEMORY_JOIN_OPTIMIZER_H
#define IN_MEMORY_JOIN_OPTIMIZER_H

#include "subset_structure.h"

void simpleTwoStepPlan(const Graph &query, HyperTree &t, const CandidateSpace &cs);
void dpTwoStepPlan(const Graph &query, HyperTree &t, CandidateSpace &cs, bool *visited, VertexID *partMatch,
                   VertexID **candidates, ui *candCount, const std::vector<SubsetStructure> &dpStructures);
void dpTwoBagPlan(const Graph &query, HyperTree &t, CandidateSpace &cs, std::vector<double> &cost1,
                  std::vector<double> &cost2, bool &refine, bool *visited, VertexID *partMatch, VertexID **candidates,
                  ui *candCount, const SubsetStructure &s0, const SubsetStructure &s1);
void
twoBagMaxSharePlan(const Graph &query, HyperTree &t, CandidateSpace &cs, std::vector<double> &cost, bool *visited,
                   VertexID *partMatch, VertexID **candidates, ui *candCount, const SubsetStructure &s0,
                   const SubsetStructure &s1);
void twoBagHeuristicPlan(const Graph &query, HyperTree &t, CandidateSpace &cs, std::vector<double> &cost, bool *visited, VertexID *partMatch,
                         VertexID **candidates, ui *candCount);
void
optPrefixTreeCost(PrefixNode *root, uint64_t prevID, HyperTree &t, const std::vector<VertexID> &nIDs,
                  const Graph &query, const std::vector<SubsetStructure> &dpStructures, double &totalCost,
                  std::vector<double> &costs, const vector<vector<VertexID>> &matchedAttrs,
                  std::vector<std::vector<VertexID>> &nodeOrders, CandidateSpace &cs,
                  const std::vector<double> &factors);
void
optSharePlan(const Graph &query, HyperTree &t, CandidateSpace &cs, PrefixNode *&pt, bool *visited, VertexID *partMatch,
             VertexID **candidates, ui *candCount, std::vector<double> &costs, std::vector<double> &noShareCosts,
             std::vector<HyperTree> &trees, std::vector<PrefixNode *> &prefixTrees, std::ofstream &outStream,
             const std::vector<SubsetStructure> &dpStructures);
void
optInitPlan(const Graph &query, HyperTree &t, CandidateSpace &cs, PrefixNode *&pt, bool *visited, VertexID *partMatch,
            VertexID **candidates, ui *candCount, double &minCost, size_t budget, bool adaptive, bool materialize,
            bool optimizeGlobal, std::ofstream &outStream, std::vector<SubsetStructure> &dpStructures);
void
simplePlan(const Graph &query, HyperTree &t, CandidateSpace &cs, PrefixNode *&pt, bool *visited, VertexID *partMatch,
           VertexID **candidates, ui *candCount, double &minCost, PrefixNode *bestPT, bool share);
void heuristicInitPlan(const Graph &query, HyperTree &t, CandidateSpace &cs, PrefixNode *&pt, bool *visited,
                       VertexID *partMatch, VertexID **candidates, ui *candCount, double &minCost, size_t budget,
                       bool adaptive, bool share, bool materialize, const std::vector<SubsetStructure> &dpStructures);
void adjustSubTree(HyperTree &t, PrefixNode *&root, VertexID nID);
void
adjustSubTree(const Graph &query, HyperTree &t, PrefixNode *&root, VertexID nID, CandidateSpace &cs,
              const std::vector<VertexID> &prefix, const std::vector<int> &mappingSizes);
void subTreePlan(const Graph &query, HyperTree &t, CandidateSpace &cs, const vector<VertexID> &nIDs,
                 const std::vector<SubsetStructure> &dpStructures, uint64_t prevID,
                 std::map<uint64_t, PrefixNode *> &bestPTs, std::map<uint64_t, double> &ptCosts,
                 const vector<vector<VertexID>> &matchedAttrs, std::map<uint64_t, std::vector<std::vector<VertexID>>> &nodeOrders,
                 std::vector<double> &factors);
void
newPlanHeuristic(const Graph &query, const HyperTree &t, CandidateSpace &cs, std::priority_queue<CandidatePlan> &plans,
                 const PrefixNode *current, const std::vector<VertexID> &attrsInPath,
                 const std::vector<VertexID> &siblingAttr, int depth, std::vector<ui> &childPoses,
                 const CandidatePlan &cp, double parentCost,
                 std::unordered_set<PrefixNode *, PrefixNodePtrHash, PrefixNodePtrEqual> &visitedPT,
                 const std::vector<ui> &nodePriority, const std::vector<std::vector<VertexID>> &attrIntersections,
                 const std::vector<std::vector<VertexID>> &sharedAttrs, const std::vector<uint64_t> &cc);

void
subTreePlan(const Graph &query, HyperTree &t, PrefixNode *&bestPT, CandidateSpace &cs, const vector<VertexID> &nIDs,
            const vector<VertexID> &prevAttrs, const vector<vector<VertexID>> &matchedAttrs,
            vector<vector<VertexID>> &localOrders, std::vector<double> &factors, bool *visited, VertexID *partMatch,
            VertexID **candidates, ui *candCount, double &minCost, bool share,
            const std::vector<SubsetStructure> &dpStructures);
void setTDExtention(HyperTree &t, const Graph &query);
void
extendPrefixTree(int depth, const Graph &query, HyperTree &t, const std::vector<VertexID> &partitionOrder, VertexID nID,
                 const std::vector<ui> &lengths, PrefixNode *&newPrefix, HyperTree *&newTree, const PrefixNode *oldPT,
                 const std::vector<ui> &childPos, CandidateSpace &cs);
double computeCostDP(const PrefixNode *pt, const Graph &query, CandidateSpace &cs, bool *visited, VertexID *partMatch,
                     VertexID **candidates, ui *candCount, std::vector<std::vector<VertexID>> &localOrders,
                     const std::vector<SubsetStructure> &dpStructures, const std::vector<double> &factors,
                     const std::vector<std::vector<VertexID>> &matchedAttrs, const std::vector<VertexID> &prevAttrs);
void twoBagDP(const Graph &query, HyperTree &t, PrefixNode *&bestPT, CandidateSpace &cs, const vector<VertexID> &nIDs,
              const vector<VertexID> &prevAttrs, vector<vector<VertexID>> &localOrders, double &minCost,
              const std::vector<SubsetStructure> &dpStructures, const vector<VertexID> &sharedAttrs);
void addOneBagDP(const Graph &query, HyperTree &t, PrefixNode *&pt, CandidateSpace &cs, VertexID nID,
                 vector<vector<VertexID>> &localOrders, double &minCost,
                 const std::vector<SubsetStructure> &dpStructures,
                 const std::vector<std::vector<VertexID>> &attrIntersections);
void bagDPHeuristic(const Graph &query, HyperTree &t, CandidateSpace &cs, PrefixNode *&pt,
                    const std::vector<SubsetStructure> &dpStructures, double &minCost);
void
bagSetDP(const Graph &query, HyperTree &t, CandidateSpace &cs, PrefixNode *&pt, vector<SubsetStructure> &dpStructures,
         double &bestCost, bool connected = true);

void bagPermutationDP(const Graph &query, HyperTree &t, CandidateSpace &cs, PrefixNode *&pt,
                      vector<SubsetStructure> &dpStructures, double &bestCost, int reorder = 1, bool connected = true,
                      bool globalOrderShare = true);
void maxSharePermutationDP(const Graph &query, HyperTree &t, CandidateSpace &cs, PrefixNode *&pt,
                           vector<SubsetStructure> &dpStructures, double &bestCost, bool connected = true);
void reorderBags(HyperTree &t, PrefixNode *pt, const std::vector<SubsetStructure> &dpStructures, int heuristic,
                 const Graph &query, CandidateSpace &cs);
#endif //IN_MEMORY_JOIN_OPTIMIZER_H
