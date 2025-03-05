//
// Created by anonymous authors on 2024/6/29.
//

#include "join.h"
#include "optimizer.h"

#ifndef IN_MEMORY_JOIN_ADAPTIVE_H
#define IN_MEMORY_JOIN_ADAPTIVE_H

extern double gGlobalTime;
extern int gNumFixedCase;
extern int gMaxNewFixed;

// Custom hash function for std::vector<ui>
struct VectorHash {
    std::size_t operator()(const std::vector<ui>& v) const {
        std::size_t hash = 0;
        std::hash<int> hasher;
        for (ui i : v) {
            hash ^= hasher(i) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
};

// Custom equality comparator for std::vector<ui>
struct VectorEqual {
    bool operator()(const std::vector<ui>& v1, const std::vector<ui>& v2) const {
        return v1 == v2;
    }
};

void
initialPlan(const Graph &query, HyperTree &t, CandidateSpace &cs, PrefixNode *&pt, bool *visited, VertexID *partMatch,
            VertexID **candidates, ui *candCount, std::vector<double> &costs, std::vector<VertexID> &largeNodes,
            size_t budget);
void
safePlan(const Graph &query, HyperTree &t, CandidateSpace &cs, PrefixNode *&pt, bool *visited, VertexID *partMatch,
         VertexID **candidates, ui *candCount, std::vector<double> &costs, size_t budget, bool share = false);

void
largeNodePlan(const Graph &query, HyperTree &t, const vector<VertexID> &largeNodes, CandidateSpace &cs, PrefixNode *&pt,
              bool *visited, VertexID *partMatch, VertexID ***candidates, ui **candCount, size_t budget,
              std::vector<double> &costs, std::vector<VertexID> &defaultPartition, std::vector<ui> &lengths);
bool adaptiveNodeJoin(const HyperTree &t, VertexID nID, CandidateSpace &cs, bool *visited,
                      const std::vector<VertexID> &notInBag, VertexID *partMatch, int mappingSize,
                      VertexID **candidates, ui *candCount, std::vector<ui> &poses,
                      std::vector<std::vector<VertexID>> &tuples, size_t &budget);
void adaptiveGlobalJoin(std::vector<std::vector<VertexID>> &result, size_t &budget, size_t &count, const Graph &query,
                        HyperTree &t, CandidateSpace &cs, VertexID *partMatch, int mappingSize, int pathLength,
                        std::vector<TrieNode *> &nodes, bool *visited, const std::vector<std::vector<VertexID>> &nIDs,
                        const std::vector<std::vector<VertexID>> &vertexParents,
                        const std::vector<VertexID> &cartesianParent, ui ***iters, ui *iterSize, bool traversal,
                        std::vector<std::vector<DynamicArray<TrieNode *> *>> &edgeColumns, TrieNode *empty, bool skip);
void materializeJoin(HyperTree &t, PrefixNode *pt, CandidateSpace &cs, bool *visited,
                     vector<vector<std::vector<VertexID>>> &tuples, std::vector<VertexID> &largeNodes, size_t &budget,
                     ui height, VertexID *partMatch, VertexID **pCandidates, ui *pCandCount, VertexID **neighbors,
                     ui *neighborCount, std::vector<std::vector<ui>> &nodePoses, std::vector<ui> &pPoses,
                     std::vector<ui> &childPoses, ui ***nodeCandidates, ui **nodeCandCount);
void adaptiveMaterializeJoin(const Graph &query, HyperTree &t, PrefixNode *pt, CandidateSpace &cs,
                             std::vector<TrieNode *> &trieNodes, bool *visited, std::vector<std::vector<VertexID>> &result,
                             size_t &count, bool traverse, size_t budget, bool skip);
std::vector<ui>
resetLengths(HyperTree &t, const std::vector<ui> &lengths, VertexID nID, const PrefixNode *oldPT, int depth,
             const HyperTree &defaultT, const vector<vector<int>> &attrIDMap);
void
resetMaterialize(const Graph &query, HyperTree &t, const HyperTree &defaultT, VertexID nID, int depth,
                 const std::vector<ui> &lengths, PrefixNode *&newPrefix, HyperTree *&newTree,
                 const PrefixNode *oldPT, CandidateSpace &cs, std::vector<bool> &fixed,
                 std::vector<int> &mappingSizes, std::vector<std::pair<int, VertexID>> &depend,
                 vector<vector<std::vector<VertexID>>> &tuples,
                 std::map<PrefixNode *, std::vector<VertexID>> &bagsBelow,
                 const std::vector<std::vector<int>> &attrIDMap, const std::vector<PrefixNode *> &dynamicPartition,
                 VertexID *partMatch);
void
resetPoses(const HyperTree &t, const HyperTree &defaultT, VertexID nID, VertexID **nodeCandidates, ui *nodeCandCount,
           std::vector<ui> &nodePoses, VertexID **pCandidates, ui *pCandCount, std::vector<ui> &pPoses, ui newLength,
           ui pathLength, const std::vector<VertexID> &resetTuple, bool *visited, CandidateSpace &cs,
           VertexID *partMatch, const std::vector<std::vector<int>> &attrIDMap);
void updateSkipMatch(std::vector<PrefixNode *> &path, int depth, ui newLength, std::vector<bool> &skipMatch,
                     const std::vector<std::vector<int>> &attrIDMap,
                     map<PrefixNode *, vector<VertexID>> &bagsBelow,
                     std::vector<PrefixNode *> &dynamicPartition);
void makeBackup(ui oldLength, ui newLength, int depth, const HyperTree &newT,
                std::vector<std::vector<std::vector<VertexID>>> &tuples,
                std::vector<std::vector<std::vector<VertexID>>> &backup, const vector<PrefixNode *> &path, VertexID nID,
                VertexID *partMatch, const std::vector<VertexID> &bagToOrder, size_t &budget);
bool buildFromExist(const HyperTree &t, VertexID nID, VertexID *partMatch, std::vector<std::vector<VertexID>> &tuples,
                    std::vector<std::vector<VertexID>> &backup, size_t &budget, std::vector<TrieNode *> &trieNodes,
                    int pathLength);
void attributeOpen(std::vector<bool> &skipMatch, CandidateSpace &cs, int attrID, int depth, PrefixNode *current,
                   std::vector<TrieNode *> &lastNodes,
                   std::vector<std::vector<DynamicArray<TrieNode *> *>> &edgeColumns, ui ***pIters, ui *pIterSizes,
                   VertexID **pCandidates, ui *pCandCount, std::vector<ui> &pPoses, DynamicArray<TrieNode *> **children,
                   VertexID *partMatch, std::vector<bool> &currentFixed, VertexID **neighbors, ui *neighborCount,
                   bool parentNewVersion);
void jump(const std::vector<PrefixNode *> &dynamicPartition, PrefixNode *&current, PrefixNode *&pn,
          std::vector<ui> &childPoses, std::vector<int> &ids, vector<PrefixNode *> &path,
          std::vector<bool> &nextCands, int &depth, const std::vector<std::vector<int>> &attrIDMap, PrefixNode *tmp);
void rebuildTrie(const HyperTree &oldT, const HyperTree &newT, std::vector<TrieNode *> &lastNodes,
                 std::vector<PrefixNode *> &path, ui extendSize, std::vector<bool> &skipBuild,
                 std::vector<std::vector<std::vector<VertexID>>> &tuples, size_t &budget,
                 std::vector<size_t> &tupleSizes, const std::vector<bool> &oldFixed);
bool moveToVertex(PrefixNode *current, std::vector<TrieNode *> &lastNodes,
                  std::vector<std::vector<TrieNode *>> &traversedNodes, TrieNode *empty, VertexID v);
bool
moveTrie(const HyperTree &t, std::vector<std::vector<TrieNode *>> &traversedNodes, std::vector<TrieNode *> &lastNodes,
         PrefixNode *attribute, VertexID *partMatch, TrieNode *empty, const std::vector<std::vector<int>> &attrIDMap,
         map<PrefixNode *, vector<VertexID>> &bagsBelow, std::vector<bool> &skipMatch, int pathLength);
void adaptiveShareJoin(const Graph &query, HyperTree &t, PrefixNode *pt, CandidateSpace &cs,
                       std::vector<TrieNode *> &trieNodes, bool *visited, std::vector<std::vector<VertexID>> &result,
                       size_t &count, bool traverse, size_t budget, bool skip);
#endif //IN_MEMORY_JOIN_ADAPTIVE_H
