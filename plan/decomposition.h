//
// Created by anonymous authors on 2024/2/27.
//

#ifndef IN_MEMORY_JOIN_DECOMPOSITION_H
#define IN_MEMORY_JOIN_DECOMPOSITION_H

#include "config.h"
#include "graph.h"
#include "candidate_space.h"
#include <stack>

struct HyperNode {
    VertexID *attributes;   // sorted by the matching order
    ui numAttributes;
    std::vector<std::vector<VertexID>> attributesBefore;
    std::vector<VertexID> cartesianParent;
    std::vector<std::vector<VertexID>> nIDs; // used for the global join
    VertexID *prefix;
    ui prefixSize;
    double width;

    HyperNode(): attributes(nullptr), prefix(nullptr), numAttributes(0), prefixSize(0), width(0.0) {}

    ~HyperNode() {
        delete[] attributes;
        delete[] prefix;
    }

    void initPoses(const Graph &query, const std::vector<std::vector<size_t>> &dist);
    void initPoses(const vector<std::vector<VertexID>> &sharedAttrs, const Graph &query, const std::vector<std::vector<size_t>> &dist, bool global);
    void copyTo(HyperNode &other) const;
};

struct HyperTree {
    HyperNode *nodes;
    std::vector<VertexID>* v2n;           // for each query vertex, store the nodes that contains it.
    ui numAttributes;
    ui numNodes;
    std::vector<std::vector<VertexID>> edges;
    std::vector<ui> compressionSizes;
    std::vector<std::vector<VertexID>> nIDs;
    int extendLevel;
    std::vector<VertexID> globalOrder;
    std::vector<std::vector<VertexID>> trieOrder;
    /*******only used for scope-style join********/
    std::vector<std::vector<VertexID>> nodesAtStep;
    // dimension 1 is nID, dimension 2 is mapping size, dimension 3 are the children to call
    std::vector<VertexID> defaultPartition;
    std::vector<std::vector<VertexID>> attributesBefore;
    bool newGlobalNode;    // whether there is a virtual global node
    std::vector<VertexID> cartesianParent;
    /*******only used for global join and traverse********/
    std::vector<VertexID> depthToNID;
    std::vector<ui> levels;
    std::vector<std::vector<VertexID>> groups;
    std::map<int, std::vector<VertexID>> adaptiveDepthToNID;
    std::map<int, std::vector<ui>> adaptiveLevels;
    std::map<int, std::vector<std::vector<VertexID>>> adaptiveGroups;

    /*********************************************/

    HyperTree(): nodes(nullptr), v2n(nullptr), numAttributes(0), numNodes(0), extendLevel(0), newGlobalNode(false) {};

    ~HyperTree() {
        delete[] nodes;
        delete[] v2n;
    }

    void initPoses(const Graph &query, const CandidateSpace &cs, bool handleNode);
    void buildTrieOrder();
    void buildTraverseStruct(const Graph &query);
    void buildTraverseAdaptive(const Graph &query, int currentExtendLevel);
    void selectRoot(VertexID &root, std::vector<std::vector<VertexID>> &cohesion,
                    std::vector<std::vector<VertexID>> &children, std::vector<VertexID> &parents);
    void buildFromTD(FHD &fhd);
    void print(const Graph &query) const;
    void addGlobalNode(const std::vector<VertexID> &globalAttrs);
    void writeToStream(ostream &outStream);
};

// the attribute tree structure
struct PrefixNode {
    VertexID u;
    std::vector<VertexID> attributesBefore;
    VertexID cartesianParent;
    // hyper nodes to call in this node
    std::vector<VertexID> nIDsToCall;
    // hyper nodes to join, for the prefix nodes in the path to the global join
    std::vector<VertexID> nIDsToJoin;
    // hyper nodes to build tries
    std::vector<VertexID> nIDsToBuild;
    std::vector<PrefixNode *> children;
    bool pathToGlobal;

    PrefixNode(VertexID id) : u(id), cartesianParent(99), pathToGlobal(true) {}
    int findChildIndex(VertexID id);

    PrefixNode(const PrefixNode& other) = delete;
    PrefixNode& operator=(const PrefixNode &other) = delete;
    ~PrefixNode() {
        for (auto & c : children)
            delete c;
    }

    ui getHeight() const;
    std::vector<PrefixNode *> locate(VertexID nID) const;
    void locate(PrefixNode *attribute, std::vector<PrefixNode *> &path, std::vector<ui> &childPoses);
    int locate(PrefixNode *attribute, VertexID firstID);
    void refine(const std::vector<std::vector<VertexID>> &sharedAttrs);
    void initPoses(const std::vector<std::vector<VertexID>> &bagAttrs, const Graph &query, const std::vector<std::vector<size_t>> &dist);
    void initNIDsToBuild(ui numNodes);
    void print() const;
    void checkCallAll(ui numNodes) const;
    std::vector<VertexID> getBagsBelow() const;
    ui numAttr(const HyperTree &t) const;

    bool operator==(const PrefixNode &rhs) const;

    bool operator!=(const PrefixNode &rhs) const;
    PrefixNode *clone() const;
    void mergeToRight(std::vector<std::vector<VertexID>> &localOrders, std::vector<VertexID> &rightCall);
    void getTraverseOrder(vector<PrefixNode *> &attributes, vector<VertexID> &nIDs, const HyperTree &t);
    void addBagsBelow(std::map<PrefixNode *, std::vector<VertexID>> &bagsBelow);
};

struct PrefixNodePtrHash {
    std::size_t operator()(const PrefixNode* ptr) const {
        return std::hash<int>()(ptr->u);
    }
};

struct PrefixNodePtrEqual {
    bool operator()(const PrefixNode* lhs, const PrefixNode* rhs) const {
        return *lhs == *rhs;
    }
};

PrefixNode *buildPrefixTree(std::vector<std::vector<VertexID>> &orders, const Graph &query,
                            const std::vector<std::vector<size_t>> &dist,
                            std::unordered_set<PrefixNode *, PrefixNodePtrHash, PrefixNodePtrEqual> &visitedPT,
                            bool &exist);
PrefixNode *fullAttributeTree(const HyperTree &t, std::vector<PrefixNode *> &attributeOrder,
                              map<PrefixNode *, vector<VertexID>> &bagsBelow, bool share);
PrefixNode *fullAttributeTree(PrefixNode *attrTree, std::vector<PrefixNode *> &attributeOrder,
                              map<PrefixNode *, vector<VertexID>> &bagsBelow, const HyperTree &t);
std::vector<std::vector<int>>
buildAttrIDMap(const vector<PrefixNode *> &attributes, const map<PrefixNode *, vector<VertexID>> &bagsBelow,
               const HyperTree &t);
std::vector<std::pair<int, VertexID>> buildDependent(const vector<PrefixNode *> &attributes);
PrefixNode* cloneTree(const PrefixNode* root);
std::vector<PrefixNode *> allPrefixTree(const vector<vector<VertexID>> &nodeAttrs, const vector<VertexID> &globalAttrs,
                                        const Graph &query, const std::vector<std::vector<size_t>> &dist);

HyperTree twoBagSubQuery(const HyperTree &t, const Graph &query);

void buildFromPrefixTree(PrefixNode *prefixTree, const std::vector<std::vector<VertexID>> &nodeOrder, HyperTree &t,
                         const std::vector<std::vector<VertexID>> &sharedAttrs, const Graph &query, CandidateSpace &cs);
void buildFromPrefixTree(const std::vector<PrefixNode *> &prefixTrees,
                         const std::vector<std::vector<std::vector<VertexID>>> &bestOrders,
                         std::vector<HyperTree> &trees, const HyperTree &reference,
                         const std::vector<std::vector<VertexID>> &sharedAttrs, const Graph &query, CandidateSpace &cs);
void removeNID(VertexID nID, PrefixNode *root);
void removeNID(VertexID nID, std::vector<PrefixNode *> &nodes, int depth, PrefixNode *&root);
std::vector<int> getMappingSizes(const HyperTree &t, const PrefixNode *pt);
uint64_t getSubsetID(const std::vector<VertexID>& subset);
bool subsetConnectivity(const Graph &query, const std::vector<uint64_t> &cc, const vector<VertexID> &subset);
bool orderConnectivity(const Graph &query, const vector<VertexID> &order, const std::vector<uint64_t> &cc);
bool orderConnectivity(const Graph &query, const vector<VertexID> &order);
void
genAllPrefixTree(const HyperTree &t, const Graph &query, CandidateSpace &cs, std::vector<PrefixNode *> &prefixTrees,
                 bool global);
void extendPrefixTree(const Graph &query, const HyperTree &t, queue<PrefixNode *> &plans, const PrefixNode *pt,
                      const PrefixNode *current, const std::vector<VertexID> &attrsInPath,
                      const std::vector<VertexID> &siblingAttr, int depth, std::vector<ui> &childPoses,
                      std::unordered_set<PrefixNode *, PrefixNodePtrHash, PrefixNodePtrEqual> &visitedPT,
                      const std::vector<std::vector<VertexID>> &attrIntersections,
                      const std::vector<std::vector<VertexID>> &sharedAttrs, const std::vector<uint64_t> &cc);
std::vector<VertexID> globalOrder(const Graph &query, const HyperTree &t, const CandidateSpace &cs, const std::vector<VertexID> &prefix);
void reorderBags(const Graph &query, HyperTree &t);
void changeNIDsToCall(const Graph &query, HyperTree &t, PrefixNode *pt);
void addGlobal(const Graph &query, HyperTree &t, PrefixNode *pt, VertexID nID, VertexID maxCostID, CandidateSpace &cs,
               vector<vector<VertexID>> &nodeOrders, std::vector<ui> &prefixSizes, bool globalOrderShare);
std::vector<int>
getPrefixAttrID(PrefixNode *pt, const std::vector<VertexID> &prefix, const std::vector<std::vector<int>> &attrIDMap,
                const std::vector<PrefixNode *> &dynamicPartition);
int convertToScopePlan(const Graph &query, HyperTree &t, PrefixNode *&pt, CandidateSpace &cs);

#endif //IN_MEMORY_JOIN_DECOMPOSITION_H
