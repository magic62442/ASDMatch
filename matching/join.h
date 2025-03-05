//
// Created by anonymous authors on 2024/3/8.
//

#ifndef IN_MEMORY_JOIN_JOIN_H
#define IN_MEMORY_JOIN_JOIN_H

#include "estimator.h"

extern size_t gNumResult;
extern size_t gNumCartesian;
extern size_t gNumTraverse;

ui binarySearch(const DynamicArray<TrieNode *> &child, const ui begin, const ui end, const ui target);
// for joining unary relation (set intersection)
ui leapfrogSeek(const DynamicArray<TrieNode *> &child, ui begin, ui end, ui target);
void leapFrogJoin(DynamicArray<TrieNode *> **children, ui num, ui **iters, ui &iterSize);

void
nodeJoinScope(const HyperTree &t, VertexID nID, CandidateSpace &cs, TrieNode *root, bool *visited, VertexID *partMatch,
              VertexID **candidates, ui *candCount, std::vector<ui> &poses, bool skip);
void nodeJoinScope(const HyperTree &t, VertexID nID, CandidateSpace &cs, bool *visited, VertexID *partMatch,
                   VertexID **candidates, ui *candCount, std::vector<ui> &poses);
void nodeJoin(const HyperTree &t, VertexID nID, CandidateSpace &cs, TrieNode *root, bool *visited, VertexID *partMatch,
              int mappingSize, VertexID **candidates, ui *candCount, std::vector<ui> &poses,
              std::vector<std::vector<VertexID>> &tuples);
void treeJoinScope(std::vector<std::vector<VertexID>> &result, CandidateSpace &cs, const HyperTree &t,
                   std::vector<TrieNode *> &nodes, bool *visited, bool skip);
void emptyHeadedJoin(const HyperTree &t, CandidateSpace &cs, std::vector<TrieNode *> &nodes, bool *visited,
                     std::vector<std::vector<VertexID>> &result, std::vector<double> &times);
void sharedJoin(const HyperTree &t, const PrefixNode *pt, const Graph &query, CandidateSpace &cs,
                std::vector<TrieNode *> &trieNodes, bool *visited, std::vector<std::vector<VertexID>> &result,
                size_t &count, bool traverse);
void
leapFrogTrieJoin(std::vector<std::vector<VertexID>> &result, const HyperTree &t, CandidateSpace &cs,
                 const std::vector<VertexID> &order, VertexID *partMatch, int mappingSize,
                 const std::vector<TrieNode *> &nodes, const std::vector<ui> &compressionSize, bool *visited,
                 int extendLevel, const std::vector<std::vector<VertexID>> &nIDs, ui ***iters, ui *iterSize);
void globalJoin(std::vector<std::vector<VertexID>> &result, size_t &count, const Graph &query, const HyperTree &t,
                CandidateSpace &cs, VertexID *partMatch, int mappingSize, int pathLength,
                std::vector<TrieNode *> &nodes, bool *visited, const std::vector<std::vector<VertexID>> &nIDs,
                const std::vector<std::vector<VertexID>> &vertexParents, const std::vector<VertexID> &cartesianParent,
                ui ***iters, ui *iterSize, bool traversal,
                std::vector<std::vector<DynamicArray<TrieNode *> *>> &edgeColumns);
void traverse(std::vector<std::vector<VertexID>> &result, const HyperTree &t, const std::vector<VertexID> &order,
              VertexID *partMatch, int mappingSize, const std::vector<TrieNode *> &nodes, bool *visited);
void
traverse(size_t &count, const Graph &query, const HyperTree &t, const std::vector<VertexID> &order, VertexID *partMatch,
         int mappingSize, const std::vector<TrieNode *> &nodes, bool *visited, int extendLevel);
void
produceResult(std::vector<std::vector<VertexID>> &result, const HyperTree &t, VertexID *partMatch,
              const std::vector<TrieNode *> &nodes, bool *visited);
void storeMatches(const std::vector<std::vector<VertexID>> &result, std::ofstream &outStream);
void buildTrie(std::vector<std::vector<VertexID>> &tuples, TrieNode *root, const std::vector<VertexID> &order);
void executeTwoBag(std::vector<std::vector<VertexID>> &result, CandidateSpace &cs, const HyperTree &t,
                   std::vector<TrieNode *> &nodes, bool *visited);
void getExample(const Graph &data, std::map<uint64_t, ui> &distribution);

#endif //IN_MEMORY_JOIN_JOIN_H
