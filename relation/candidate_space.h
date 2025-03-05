//
// Created by anonymous authors on 2024/2/29.
//

#ifndef IN_MEMORY_JOIN_CANDIDATE_SPACE_H
#define IN_MEMORY_JOIN_CANDIDATE_SPACE_H

#include "graph.h"
#include "compute_set_intersection.h"
#include <queue>
#include <set>

struct Dag {
    VertexID root;
    std::vector<std::vector<VertexID>> edges;
    std::vector<std::vector<VertexID>> children;
    std::vector<std::vector<VertexID>> parents;
    std::vector<VertexID> bfsSequence;

    Dag() = default;
    explicit Dag(const Graph &query);
    virtual ~Dag() = default;
    Dag &operator=(const Dag &rhs);
    void buildDAG(const Graph &query, const Graph &data);
    void selectRoot(const Graph &query, const Graph &data);
};

class CandidateSpace {
private:
    bool **_bitsetVertex;
    bool **_bitsetEdge;
    // required for candidate filtering
    std::map<LabelID, Triangle *> *_queryTriangle;
    std::map<LabelID, FourCycle *> *_queryFourCycle;
    std::map<LabelID, ui> *_queryNumTriangle;
    std::map<LabelID, ui> *_queryNumFourCycle;
    std::map<LabelID, Triangle *> *_dataTriangle;
    std::map<LabelID, FourCycle *> *_dataFourCycle;
    std::map<LabelID, ui> *_dataNumTriangle;
    std::map<LabelID, ui> *_dataNumFourCycle;

    const Graph &_query;
    const Graph &_data;
    Dag _dag;

public:
    std::vector<std::vector<VertexID>> candidateSet;
    // index candidate edge by u1, v(a candidate of u1), u2
    std::vector<std::vector<std::vector<std::vector<VertexID>>>> candidateEdge;
    // all source shortest path for the query graph
    std::vector<std::vector<VertexID>> next; // Next vertex on the shortest path, used for cartesian product
    std::vector<std::vector<size_t>> dist;

    CandidateSpace(const Graph &query, const Graph &data, std::map<LabelID, Triangle *> *queryTriangles=nullptr,
                   std::map<LabelID, Triangle *> *dataTriangles=nullptr,
                   std::map<LabelID, FourCycle *> *queryFourCycles=nullptr,
                   std::map<LabelID, FourCycle *> *dataFourCycles=nullptr,
                   std::map<LabelID, ui> *queryNumTriangles=nullptr,
                   std::map<LabelID, ui> *queryNumFourCycles=nullptr,
                   std::map<LabelID, ui> *dataNumTriangles=nullptr, std::map<LabelID, ui> *dataNumFourCycles=nullptr);

    ~CandidateSpace() = default;

    bool buildCandVEQ();
    bool init();
    bool filter();
    bool filter(VertexID u, VertexID v);
    void construct();
    void buildCandidateEdge(bool checkEdge);
    bool TriangleSafety(EdgeID queryEdge, EdgeID dataEdge);
    bool FourCycleSafety(EdgeID queryEdge, EdgeID dataEdge);
    bool edgeBipartiteSafety(VertexID u, VertexID v);
    bool buildCandCFL();
    VertexID selectCFLRoot();
    void buildCFLBFSTree(VertexID root, std::vector<std::vector<VertexID>> &levels, std::vector<VertexID> &order,
                         std::vector<std::vector<VertexID>> &before, std::vector<std::vector<VertexID>> &lowerLevelAfter,
                         std::vector<std::vector<VertexID>> &sameLevelAfter,
                         std::vector<std::vector<VertexID>> &children);
    void writeToStream(std::ofstream &outStream);
    void setQueryGraphWeights(Graph &query);
    size_t getDist(VertexID u1, VertexID u2) const;
    std::vector<VertexID> reconstructPath(VertexID i, VertexID j) const;
    ui getMaxSize() const;
    bool checkExists(VertexID v1, VertexID v2) const {
        return _data.getEdgeID(v1, v2) != -1;
    };
};

struct FilterVertex {
    double penalty;
    int filteredTime;
    VertexID u;

    FilterVertex(double penalty, int filteredTimes, VertexID u) : penalty(penalty), filteredTime(filteredTimes),
                                                                  u(u) {}

    bool operator<(const FilterVertex &rhs) const {
        return penalty > rhs.penalty;
    }
    bool operator>(const FilterVertex &rhs) const {
        return penalty < rhs.penalty;
    }
};

#endif //IN_MEMORY_JOIN_CANDIDATE_SPACE_H
