//
// Created by Qiyan LI on 2024/9/19.
//

#ifndef IN_MEMORY_JOIN_SUBSET_STRUCTURE_H
#define IN_MEMORY_JOIN_SUBSET_STRUCTURE_H

#include "estimator.h"
#include <unordered_map>
#include <bitset>

// using 165 bits to encode a permutation
typedef std::bitset<165> Permutation;

// a data structure for the subset lattice, used for dp-based order selection
struct SubsetStructure {
    std::vector<VertexID> elements;
    ui n;
    std::vector<uint64_t> cc;
    std::vector<std::vector<uint64_t>> subsets;
    std::vector<std::vector<std::vector<int>>> subsetOf;
    std::vector<std::vector<std::vector<int>>> supersetOf;
    std::vector<std::vector<std::vector<VertexID>>> orders;
    std::vector<std::vector<double>> subsetToCost;
    std::vector<std::vector<std::vector<VertexID>>> remainingOrders;
    std::vector<std::vector<double>> remainingCost;
    std::map<uint64_t, std::map<uint64_t, double>> extendCost;
    std::map<uint64_t, std::map<uint64_t, std::vector<VertexID>>> extendOrder;

    SubsetStructure(const std::vector<VertexID> &elements, const Graph &query, bool connected);
    SubsetStructure() = default;
    int getPosition(uint64_t id, ui k) const {
        auto it = std::lower_bound(subsets[k].begin(), subsets[k].end(), id);
        if (it != subsets[k].end() && *it == id) {
            return std::distance(subsets[k].begin(), it);
        }
        return -1; // Not found
    }

    // Generate relationships between subsets and supersets
    void generateSubsetSupersetRelationships();
    void optimalPlanDP(const Graph &query, CandidateSpace &cs, bool *visited, VertexID *partMatch,
                       VertexID **candidates, ui *candCount);
    void reverseDP(const Graph &query, CandidateSpace &cs);
    void buildExtendCost(const Graph &query, CandidateSpace &cs);
    void prefixPlan(uint64_t prefixID, ui k, std::vector<VertexID> &optOrder, double &optCost) const;
    void extendPrefixPlan(const Graph &query, CandidateSpace &cs, uint64_t id1, ui k1, uint64_t id2, ui k2,
                          std::vector<VertexID> &optOrder, double &optCost) const;
    std::vector<VertexID> getOptOrder() const { return orders[n][0]; }
    std::vector<VertexID> getOptOrder(uint64_t subsetID, ui k) const {
        if (k == 0) return std::vector<VertexID>();
        int pos = getPosition(subsetID, k);
        return orders[k][pos];
    }
    std::vector<VertexID> getOptOrder(uint64_t prefixID, ui k1, uint64_t subsetID, ui k2, const Graph &query,
                                      CandidateSpace &cs) const;
    double getOptCost() const { return subsetToCost[n][0]; }
    double getOptCost(uint64_t subsetID, ui k) const {
        if (k == 0) return 0.0;
        int pos = getPosition(subsetID, k);
        return subsetToCost[k][pos];
    }
    double getRemainingCost(uint64_t subsetID, ui k) const {
        if (k == 0) return getOptCost();
        int pos = getPosition(subsetID, k);
        return remainingCost[k][pos];
    }
    std::vector<uint64_t> getSupersets(uint64_t subsetID, ui k) const {
        int pos = getPosition(subsetID, k);
        std::vector<uint64_t> supersets;
        for (int supersetPos: supersetOf[k][pos])
            supersets.push_back(subsets[k + 1][supersetPos]);
        return supersets;
    }
    bool saveToSteam(std::ofstream &ofs) const;
    bool readFromStream(std::ifstream &ifs);
};

struct bagSetStructure {
    std::map<uint64_t, std::map<uint64_t,double>> indepCosts;
    std::map<uint64_t, std::map<uint64_t,std::vector<std::vector<VertexID>>>> localOrders;
    std::map<uint64_t, uint64_t> intersections;
    void init(const std::vector<std::vector<VertexID>> &sharedAttrs,
              const std::vector<SubsetStructure> &dpStructures, const Graph &query, bool connected = true);
};

struct bagPermuteStructure{
    std::map<uint64_t, std::unordered_map<Permutation ,double>> indepCosts;
    std::map<uint64_t, std::unordered_map<Permutation ,double>> subtreeCosts;
    std::map<uint64_t, std::unordered_map<Permutation ,std::vector<std::vector<VertexID>>>> localOrders;
    std::map<uint64_t, std::unordered_map<Permutation, std::vector<int>>> orderIDs;
    std::vector<std::vector<std::vector<VertexID>>> basicOrders;
    std::map<uint64_t, uint64_t> intersections;
    void init(const std::vector<std::vector<VertexID>> &sharedAttrs, const std::vector<SubsetStructure> &dpStructures,
              const Graph &query, bool maxShare = false, bool connected = true);
    void initSubtreeCost(const std::vector<std::vector<VertexID>> &sharedAttrs,
                         const std::vector<SubsetStructure> &dpStructures, const Graph &query,
                         CandidateSpace &cs, bool connected = true);
    void initBasicOrder(const std::vector<std::vector<VertexID>> &sharedAttrs, const std::vector<SubsetStructure> &dpStructures,
                        const Graph &query, bool maxShare = false, bool connected = true);
};

VertexID findExtraElement(uint64_t idSubsetK, uint64_t idSubsetKMinus1);
ui subsetSize(uint64_t subsetID, VertexID n);
void generateSubsetsK(const std::vector<VertexID> &set, int size, int index, std::vector<VertexID> &current,
                      vector<uint64_t> &subsets, const std::vector<uint64_t> &cc, const Graph &query, bool connected);
void generateExtendedSubsets(const std::vector<VertexID>& set, uint64_t subsetID, std::vector<uint64_t>& extendedSubsets);
void
generateSubsets(const std::vector<VertexID> &set, vector<uint64_t> &subsets);
std::vector<VertexID> getSubsetFromID(uint64_t id, ui maxVID);
void smallCard(const Graph &query, CandidateSpace &cs);
bool isSubset(uint64_t subset1ID, uint64_t subset2ID);
Permutation encodePermutation(const std::vector<uint32_t>& elements);
std::vector<uint32_t> decodePermutation(const Permutation & encoded);
std::vector<Permutation>
extensions(const Permutation &pi, const std::vector<VertexID> &elements, const std::vector<uint64_t> &cc,
           const Graph &query, bool connected);
SubsetStructure buildSubStructure(const SubsetStructure& original, uint64_t subsetID, ui k);
SubsetStructure buildSubStructure(const SubsetStructure& original, uint64_t id1, ui k1, uint64_t id2, ui k2);


#endif //IN_MEMORY_JOIN_SUBSET_STRUCTURE_H
