//
// Created by anonymous authors on 2024/2/27.
//

#ifndef IN_MEMORY_JOIN_ESTIMATOR_H
#define IN_MEMORY_JOIN_ESTIMATOR_H

#include "relation.h"
#include "candidate_space.h"
#include "decomposition.h"
#include "utility"
#include <list>

extern std::unordered_map<uint64_t, double> subsetToCard;

struct VertexPriority {
//    bool independent;   // dynamically updated. whether neighbors of the vertex are all in the front
    ui num;        // number of occurrences in all tree nodes
    ui backwardNbr;
    ui degree;
    ui candSize;

//    VertexPriority() : independent(true), num(0), backwardNbr(0), degree(0), candSize(0) {}
    VertexPriority() : num(0), backwardNbr(0), degree(0), candSize(std::numeric_limits<ui>::max()) {}

    bool operator<(const VertexPriority &rhs) const {
//        if (independent != rhs.independent)
//            return rhs.independent;  // false < true
        if (num != rhs.num)
            return num > rhs.num;    // more occurrences should come first
        if (backwardNbr != rhs.backwardNbr)
            return backwardNbr > rhs.backwardNbr;  // more backward neighbors come first
        if (degree != rhs.degree)
            return degree > rhs.degree;  // larger degree comes first
        return candSize < rhs.candSize;  // smaller candidate size comes first
    }
};

struct CandidatePlan {
    PrefixNode *pt;
    std::vector<std::vector<VertexID>> orders;
    std::vector<std::vector<ui>> numBackNbr;
    std::vector<VertexID> nodePriority;
    std::vector<double> nodeCosts;
    int shareNum;
    double parentCost;
    CandidatePlan(PrefixNode *pt, const std::vector<std::vector<VertexID>> &orders, const vector<std::vector<ui>> &numBackNbr,
                  const vector<VertexID> &nodePriority, int shareNum, double parentCost) :
            pt(pt), orders(orders), numBackNbr(numBackNbr), nodePriority(nodePriority), shareNum(shareNum), parentCost(parentCost) {}

    virtual ~CandidatePlan() = default;

    bool operator<(const CandidatePlan &rhs) const {
        if (parentCost < rhs.parentCost) return true;
        else if (parentCost == rhs.parentCost) {
            for (VertexID nID = 0; nID < nodePriority.size(); ++nID) {
                if (this->numBackNbr[nID] > rhs.numBackNbr[nID])
                    return true;
                else if (this->numBackNbr[nID] == rhs.numBackNbr[nID]) continue;
                else return false;
            }
            return shareNum > rhs.shareNum;
        }
        return false;
    }
};

void initPoses(const std::vector<VertexID> &order, const Graph &query, std::vector<std::vector<VertexID>> &vertexParents,
               std::vector<VertexID> &cartesianParent, const std::vector<std::vector<size_t>> &dist);

void
optimizedCartesianProduct(CandidateSpace &cs, VertexID v, const std::vector<VertexID> &path, VertexID *candidate,
                          ui &candCount);
void
optimizedCartesianProduct(CandidateSpace &cs, VertexID v, const std::vector<VertexID> &path, DynamicArray<TrieNode *> *edgeColumn);

void
cardEstimateLabeled(const std::vector<VertexID> &order, const std::vector<std::vector<VertexID>> &vertexParents,
                    const std::vector<VertexID> &cartesianParent, CandidateSpace &cs, bool *visited,
                    VertexID *partMatch, VertexID **candidates, ui *candCount, std::vector<ui> &poses,
                    std::vector<double> &estimation);

void
cardEstimateLabeled(VertexID u, bool lastLevel, const std::vector<VertexID> &prevOrder, const Graph &query,
                    CandidateSpace &cs, std::vector<std::vector<VertexID>> &prevMatches,
                    std::vector<std::vector<VertexID>> &nextMatches, std::vector<double> &weights,
                    std::vector<double> &nextWeights, double &estimation);

double estimateCartesian(VertexID u1, VertexID u2, CandidateSpace &cs);

std::vector<VertexID> RIOrder(const Graph &query);
std::vector<VertexID> GQLOrder(const Graph &query, const CandidateSpace &cs);
std::vector<VertexID> readLastLineAsVector(const std::string& filename);
std::vector<VertexID> simpleOrder(const Graph &query, const CandidateSpace &cs, const std::vector<VertexID> &vertices,
                                  const std::vector<int> &repetitions);
std::vector<VertexID>
simpleOrder(const Graph &query, const CandidateSpace &cs, const std::vector<VertexID> &prefix,
            std::vector<ui> &numBackWard, const std::vector<VertexID> &prevOrder,
            int noChangePos);

std::vector<ui>
computeNumBackWard(const Graph &query, const std::vector<VertexID> &prefix, const std::vector<VertexID> &localOrder);
std::vector<double> computeCost(const std::vector<VertexID> &order, const Graph &query, CandidateSpace &cs,
                                bool *visited = nullptr, VertexID *partMatch = nullptr, VertexID **candidates = nullptr,
                                ui *candCount = nullptr);
double computeCost(const std::vector<VertexID> &order, const std::vector<std::vector<VertexID>> &vertexParents,
                   const std::vector<VertexID> &cartesianParent, const std::vector<double> &cards, CandidateSpace &cs);
double computeCost(const std::vector<VertexID> &prefix, const std::vector<VertexID> &order, const Graph &query,
                   CandidateSpace &cs, bool *visited, VertexID *partMatch, VertexID **candidates, ui *candCount);
double computeCost(const PrefixNode *pt, const std::vector<std::vector<VertexID>> &localOrders, const Graph &query,
                   CandidateSpace &cs, bool *visited, VertexID *partMatch, VertexID **candidates, ui *candCount,
                   const vector<vector<VertexID>> &matchedAttrs, const std::vector<VertexID> &prevAttrs,
                   const std::vector<double> &factors);
ui maxNumBackWard(const std::vector<VertexID> &order, const Graph &query);
bool saveSubsetToCard(std::ofstream& ofs);
bool loadSubsetToCard(std::ifstream& ifs);

#endif //IN_MEMORY_JOIN_ESTIMATOR_H
