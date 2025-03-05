//
// Created by anonymous authors on 2024/2/27.
//

#ifndef IN_MEMORY_JOIN_GRAPH_H
#define IN_MEMORY_JOIN_GRAPH_H


#include "config.h"
#include "utils.h"
#include <map>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <queue>
#ifndef __APPLE__
#include <immintrin.h>
#include <x86intrin.h>
#endif

// for candidate filtering
struct Triangle {
    EdgeID secondEdge, thirdEdge;
};

struct FourCycle {
    // 1 -[first edge]- 2 -[second edge]- 3 -[third edge]- 4 -[fourth edge]- 1
    // there may aslo exists 1 - 3 and 2 -4 edges
    EdgeID secondEdge, thirdEdge, fourthEdge, oneThreeEdge, twoFourEdge;
};

class Graph {
private:
    ui _numVertices;           // number of vertices in the graph
    ui _numEdges;              // number of directed edges in the graph.
    ui _numLabels;
    LabelID *_labels;
    // undirected graph has _numEdges / 2 undirected edges.
    EdgeID *_offsets;          // vertex v's edges are [offset[v], offset[v + 1])
    EdgeID *_reverseID;
    VertexID *_nbrs;          // incoming or outgoing neighbors
    std::map<LabelID, ui> *_nlc;
    std::map<LabelID, EdgeID> *_labelOffsets; // vertex v's edges with label l starts from [offset[labelOffset[l]]
    std::map<LabelID, ui> _labelCount;
    ui *_offset2;   // offset of 'verticesByLabel'
    VertexID *_verticesByLabel; // vertices grouped by label
    // only used for query graphs
    ui *_nodeWeight;
    ui *_edgeWeight;

public:
    Graph();

    virtual ~Graph();

    const ui &getNumVertices() const {
        return _numVertices;
    }

    const ui &getNumEdges() const {
        return _numEdges;
    }

    ui getDegree(VertexID v) const {
        return _offsets[v + 1] - _offsets[v];
    }

    ui getMaxDegree() const {
        ui maxDegree = 0;
        for (VertexID i = 0; i < _numVertices; ++i) {
            ui degree = _offsets[i + 1] - _offsets[i];
            if (degree > maxDegree)
                maxDegree = degree;
        }

        return maxDegree;
    }

    const EdgeID *getOffsets() const {
        return _offsets;
    }

    const VertexID *getNbrs() const {
        return _nbrs;
    }

    const ui * getVerticesByLabel(const LabelID l, ui& count) const {
        count = _offset2[l + 1] - _offset2[l];
        return _verticesByLabel + _offset2[l];
    }

    VertexID *const getNeighbors(const VertexID v, ui& count) const {
        count = _offsets[v + 1] - _offsets[v];
        return _nbrs + _offsets[v];
    }

    LabelID getVertexLabel(const VertexID v) const {
        return _labels[v];
    }

    LabelID *getLabels() const {
        return _labels;
    }

    EdgeID getReverseID(const EdgeID e) const {
        return _reverseID[e];
    }

    EdgeID getEdgesByLabel(const VertexID v, const LabelID label, ui& count) const {
        count = _nlc[v][label];
        return _labelOffsets[v][label];
    }

    const std::map<LabelID, ui> &getNLC(const VertexID u) const {
        return _nlc[u];
    }

    ui getVertexWeight(const VertexID u) const {
        return _nodeWeight[u];
    }

    ui getEdgeWeight(const EdgeID e) const {
        return _edgeWeight[e];
    }

    void setVertexWeight(VertexID u, ui weight) {
        _nodeWeight[u] = weight;
    }

    void setEdgeWeight(EdgeID e, ui weight) {
        _edgeWeight[e] = weight;
    }

    void setLabel(VertexID u, LabelID l) {
        _labels[u] = l;
    }

    void loadGraphFromTextFile(const std::string& file);
    void buildReverseID();
    void buildNLC();
    EdgeID getEdgeID(VertexID v, VertexID w) const;
    EdgeID getUndirectedEID(VertexID v, VertexID w) const;
    void indexTriangles(std::map<LabelID, Triangle *> *&triangle, std::map<LabelID, ui> *&numTriangle);
    void indexFourCycles(std::map<LabelID, FourCycle *> *&fourCycle, std::map<LabelID, ui> *&numFourCycle);
    void readFromStream(std::ifstream &inFile);
    void writeToStream(std::ofstream &outFile) const;
    void initWeights();
    void allSourcesShortestPaths(std::vector<std::vector<size_t>> &dist, std::vector<std::vector<VertexID>> &next) const;
    void computeConnectedComponents(const std::vector<VertexID> &vertices, std::vector<std::vector<VertexID>> &components) const;
    bool isConnected(const std::vector<VertexID> &vertices) const;
    std::vector<VertexID> bfsForCC(VertexID start, set<VertexID> &visited, const std::vector<VertexID> &vertices) const;
    void buildHyperGraph(HyperG &h) const;
    void buildFHD(FHD &fhd) const;
    void sortNeighborsUnlabeled() {
        for (ui i = 0; i < _numVertices; ++i) {
            std::sort(_nbrs + _offsets[i], _nbrs + _offsets[i + 1]);
        }
    }
    void generateRandomGraph(ui numVertices, ui numEdges);
    Graph generateRandomSubgraph(ui numVertices, ui numEdges) const;
    void writeToTextFile(const std::string &filename) const;
};

void readGraphBinary(const std::string &file, Graph &g, std::map<LabelID, Triangle *> *&triangle, std::map<LabelID, FourCycle *> *&fourCycle, std::map<LabelID, ui> *&numTriangle, std::map<LabelID, ui> *&numFourCycle);
void writeGraphBinary(const std::string &file, const Graph &g, std::map<LabelID, Triangle *> *&triangle, std::map<LabelID, FourCycle *> *&fourCycle, std::map<LabelID, ui> *&numTriangle, std::map<LabelID, ui> *&numFourCycle);
void release(std::map<LabelID, Triangle *> *queryTriangle, std::map<LabelID, FourCycle *> *queryFourCycle,
             std::map<LabelID, ui> *queryNumTriangle, std::map<LabelID, ui> *queryNumFourCycle,
             std::map<LabelID, Triangle *> *dataTriangle, std::map<LabelID, FourCycle *> *dataFourCycle,
             std::map<LabelID, ui> *dataNumTriangle, std::map<LabelID, ui> *dataNumFourCycle, ui qm, ui gm);
bool canAddVertex(const Graph& query, const std::vector<VertexID>& currentSet, VertexID newVertex);
void maxIndependentSet(const Graph& graph, const std::vector<VertexID>& subset, std::vector<VertexID>& maxSet);

#endif //IN_MEMORY_JOIN_GRAPH_H
