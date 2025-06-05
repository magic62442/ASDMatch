//
// Created by Qiyan LI on 2024/2/29.
//

#include "candidate_space.h"

BipartiteMaximumMatching BPSolver;

Dag::Dag(const Graph &query) {
    ui querySize = query.getNumVertices();
    children.resize(querySize);
    parents.resize(querySize);
    bfsSequence.resize(querySize);
    root = -1;
}

void Dag::buildDAG(const Graph &query, const Graph &data) {
    bool *visited = new bool[query.getNumVertices()];
    bool *popped = new bool[query.getNumVertices()];
    memset(visited, false, sizeof(bool) * query.getNumVertices());
    memset(popped, false, sizeof(bool) * query.getNumVertices());
    if (root == -1) selectRoot(query, data);
    bfsSequence[0] = root;
    visited[root] = true;
    ui begin = 0, end = 1;
    while (begin < end) {
        ui oldEnd = end;
        while (begin < oldEnd) {
            VertexID parent = bfsSequence[begin];
            ++begin;
            popped[parent] = true;
            ui numNeighbor;
            const VertexID *neighbors = query.getNeighbors(parent, numNeighbor);
            for (int i = 0; i < numNeighbor; ++i) {
                VertexID child = neighbors[i];
                if (!popped[child]) {
                    children[parent].push_back(child);
                    parents[child].push_back(parent);
                    if (!visited[child]) {
                        visited[child] = true;
                        bfsSequence[end] = child;
                        end += 1;
                    }
                }
            }
        }
    }
    delete[] visited;
    delete[] popped;
}

void Dag::selectRoot(const Graph &query, const Graph &data) {
    root = 0;
    double minRatio = static_cast<double>(data.getNumVertices());
    for (VertexID u = 0; u < query.getNumVertices(); ++u) {
        const std::map<LabelID, ui> &uNLC = query.getNLC(u);
        ui candSize = 0;
        ui totalCount;
        LabelID l = query.getVertexLabel(u);
        const VertexID *dataV = data.getVerticesByLabel(l, totalCount);
        for (int i = 0; i < totalCount; ++i) {
            VertexID v = dataV[i];
            const std::map<LabelID, ui> &vNLC = data.getNLC(v);
            if (nlcValid(uNLC, vNLC))
                ++candSize;
        }
        double ratio = static_cast<double>(candSize) / static_cast<double>(query.getDegree(u));
        if (ratio < minRatio) {
            minRatio = ratio;
            root = u;
        }
    }
}

Dag &Dag::operator=(const Dag &rhs) {
    // Self-assignment check
    if (this == &rhs) return *this;
    root = rhs.root;
    edges = rhs.edges;
    children = rhs.children;
    parents = rhs.parents;
    bfsSequence = rhs.bfsSequence;

    return *this;
}

CandidateSpace::CandidateSpace(const Graph &query, const Graph &data, std::map<LabelID, Triangle *> *queryTriangles,
                               std::map<LabelID, Triangle *> *dataTriangles,
                               std::map<LabelID, FourCycle *> *queryFourCycles, std::map<LabelID, FourCycle *> *dataFourCycles, std::map<LabelID, ui> *queryNumTriangles,
                               std::map<LabelID, ui> *queryNumFourCycles, std::map<LabelID, ui> *dataNumTriangles, std::map<LabelID, ui> *dataNumFourCycles)
        : _query(query), _data(data), _queryTriangle(queryTriangles), _dataTriangle(dataTriangles),
        _queryFourCycle(queryFourCycles), _dataFourCycle(dataFourCycles), _queryNumTriangle(queryNumTriangles),
        _queryNumFourCycle(queryNumFourCycles), _dataNumTriangle(dataNumTriangles), _dataNumFourCycle(dataNumFourCycles) {
    _bitsetVertex = new bool*[query.getNumVertices()];
    for (int i = 0; i < query.getNumVertices(); ++i) {
        _bitsetVertex[i] = new bool[data.getNumVertices()];
        memset(_bitsetVertex[i], false, data.getNumVertices());
    }
    _bitsetEdge = new bool*[query.getNumEdges()];
    for (int i = 0; i < query.getNumEdges(); ++i) {
        _bitsetEdge[i] = new bool[data.getNumEdges()];
        memset(_bitsetEdge[i], false, data.getNumEdges());
    }
    _dag = Dag(query);
    candidateSet.resize(query.getNumVertices());
    candidateEdge.resize(query.getNumVertices());
}

bool CandidateSpace::buildCandVEQ() {
    if (!init()) return false;
    if (!filter()) return false;
    construct();
    return true;
}

bool CandidateSpace::init() {
    _dag.buildDAG(_query, _data);
//    BPSolver.global_initialize(_query.getMaxDegree(), _data.getMaxDegree());
    // the root candidate : nlc
    VertexID root = _dag.bfsSequence[0];
    LabelID rootLabel = _query.getVertexLabel(root);
    ui count;
    const VertexID *dataV = _data.getVerticesByLabel(rootLabel, count);
    const std::map<LabelID, ui> &rootNLC = _query.getNLC(root);
    for (int i = 0; i < count; ++i) {
        VertexID v = dataV[i];
        const std::map<LabelID, ui> &vNLC = _data.getNLC(v);
        if (nlcValid(rootNLC, vNLC)) {
            candidateSet[root].push_back(v);
            _bitsetVertex[root][v] = true;
        }
    }
    ui *numVisitedParent = new ui[_data.getNumVertices()];
    memset(numVisitedParent, 0, sizeof(ui) * _data.getNumVertices());
    const VertexID *dataNbr = _data.getNbrs();
    const EdgeID *queryOffset = _query.getOffsets();
    const VertexID *queryNbr = _query.getNbrs();
    for (int i = 1; i < _query.getNumVertices(); ++i) {
        VertexID u = _dag.bfsSequence[i];
        LabelID l = _query.getVertexLabel(u);
        const std::map<LabelID, ui> &uNLC = _query.getNLC(u);
        ui numParent = 0;
        for (int j = 0; j < _dag.parents[u].size(); ++j) {
            VertexID parentU = _dag.parents[u][j];
            EdgeID queryEdge = _query.getEdgeID(parentU, u);
            for (VertexID parentV : candidateSet[parentU]) {
                ui numNbr;
                EdgeID startEdge = _data.getEdgesByLabel(parentV, l, numNbr);
                for (int k = 0; k < numNbr; ++k) {
                    EdgeID dataEdge = startEdge + k;
                    VertexID v = dataNbr[dataEdge];
                    if (numVisitedParent[v] < numParent) continue;
                    const std::map<LabelID, ui> &vNLC = _data.getNLC(v);
                    if (!nlcValid(uNLC, vNLC)) continue;
                    if (numVisitedParent[v] == numParent) {
                        numVisitedParent[v] += 1;
                        if (numVisitedParent[v] == 1) {
                            candidateSet[u].push_back(v);
                            _bitsetVertex[u][v] = true;
                        }
                    }
                }
            }
            ++numParent;
        }
        for (int j = 0; j < candidateSet[u].size(); ++j) {
            VertexID v = candidateSet[u][j];
            _bitsetVertex[u][v] = false;
            if (numVisitedParent[v] == numParent) _bitsetVertex[u][v] = true;
            else {
                candidateSet[u][j] = candidateSet[u].back();
                candidateSet[u].pop_back();
                --j;
            }
            numVisitedParent[v] = 0;
        }
    }
    delete[] numVisitedParent;

    for (VertexID u = 0; u < _query.getNumVertices(); ++u) {
        if (candidateSet[u].empty())
            return false;
    }
    for (VertexID u1 = 0; u1 < _query.getNumVertices(); ++u1) {
        for (EdgeID queryEdge = queryOffset[u1]; queryEdge < queryOffset[u1 + 1]; ++queryEdge) {
            VertexID u2 = queryNbr[queryEdge];
            LabelID l2 = _query.getVertexLabel(u2);
            for (VertexID v1 : candidateSet[u1]) {
                ui numNbr;
                EdgeID startEdge = _data.getEdgesByLabel(v1, l2, numNbr);
                for (int k = 0; k < numNbr; ++k) {
                    EdgeID dataEdge = startEdge + k;
                    VertexID v2 = dataNbr[dataEdge];
                    if (_bitsetVertex[u2][v2]) {
                        _bitsetEdge[queryEdge][dataEdge] = true;
                    }
                }
            }
        }
    }

    return true;
}

bool CandidateSpace::TriangleSafety(EdgeID queryEdge, EdgeID dataEdge) {
    std::map<LabelID, Triangle *> qTriangle = _queryTriangle[queryEdge];
    std::map<LabelID, Triangle *> dTriangle = _dataTriangle[dataEdge];
    for (const auto &iter : qTriangle) {
        LabelID l = iter.first;
        ui qCount = _queryNumTriangle[queryEdge][l];
        ui dCount = _dataNumTriangle[dataEdge][l];
        if (dCount < qCount) return false;
        if (dCount > 2 * qCount) return true;
        for (int i = 0; i < qCount; ++i) {
            Triangle qT = qTriangle[l][i];
            bool found = false;
            for (int j = 0; j < dCount; ++j) {
                Triangle dT = dTriangle[l][j];
                bool valid = true;
                valid &= _bitsetEdge[qT.secondEdge][dT.secondEdge];
                if (!valid) continue;
                valid &= _bitsetEdge[qT.thirdEdge][dT.thirdEdge];
                if (!valid) continue;
                found = true;
                break;
            }
            if (!found) return false;
        }
    }

    return true;
}

bool CandidateSpace::FourCycleSafety(EdgeID queryEdge, EdgeID dataEdge) {
    if (_dataFourCycle == nullptr) return true;
    std::map<LabelID, FourCycle *> qFourCycle = _queryFourCycle[queryEdge];
    std::map<LabelID, FourCycle *> dFourCycle = _dataFourCycle[dataEdge];
    for (const auto &iter : qFourCycle) {
        LabelID l = iter.first;
        ui qCount = _queryNumFourCycle[queryEdge][l];
        ui dCount = _dataNumFourCycle[dataEdge][l];
        if (dCount < qCount) return false;
        if (dCount > 2 * qCount) return true;
        for (int i = 0; i < qCount; ++i) {
            FourCycle qF = qFourCycle[l][i];
            bool found = false;
            for (int j = 0; j < dCount; ++j) {
                FourCycle dF = dFourCycle[l][j];
                bool valid = true;
                valid &= _bitsetEdge[qF.secondEdge][dF.secondEdge];
                if (!valid) continue;
                valid &= _bitsetEdge[qF.thirdEdge][dF.thirdEdge];
                if (!valid) continue;
                valid &= _bitsetEdge[qF.fourthEdge][dF.fourthEdge];
                if (!valid) continue;
                if (qF.oneThreeEdge != -1)
                    valid &= (dF.oneThreeEdge != -1 && _bitsetEdge[qF.oneThreeEdge][dF.oneThreeEdge]);
                if (qF.twoFourEdge != -1)
                    valid &= (dF.twoFourEdge != -1 && _bitsetEdge[qF.twoFourEdge][dF.twoFourEdge]);
                if (valid) {
                    found = true;
                    break;
                }
            }
            if (!found) return false;
        }
    }


    return true;
}

bool CandidateSpace::edgeBipartiteSafety(VertexID u, VertexID v) {
    ui uDegree = _query.getDegree(u), vDegree = _data.getDegree(v);
    if (_query.getDegree(u) == 1) return true;
    const EdgeID *queryOffset = _query.getOffsets();
    const VertexID *queryNbrs = _query.getNbrs();
    const EdgeID *dataOffset = _data.getOffsets();
    const VertexID *dataNbrs = _data.getNbrs();
    BPSolver.reset();
    std::vector<std::pair<int, int>> poses;
    for (int i = 0; i < uDegree; ++i) {
        EdgeID queryEdgeID = queryOffset[u] + i;
        for (int j = 0; j < vDegree; ++j) {
            EdgeID dataEdgeID = dataOffset[v] + j;
            if (_bitsetEdge[queryEdgeID][dataEdgeID]) {
                BPSolver.add_edge(i, j);
                poses.emplace_back(i, j);
            }
        }
    }
    bool flag = BPSolver.FindUnmatchableEdges(uDegree);
    if (!flag) return false;
    for (auto &p: poses) {
        int i = p.first, j = p.second;
        if (!BPSolver.matchable[i][j]) {
            EdgeID queryEdgeID = queryOffset[u] + i;
            EdgeID dataEdgeID = dataOffset[v] + j;
            _bitsetEdge[queryEdgeID][dataEdgeID] = false;
            _bitsetEdge[_query.getReverseID(queryEdgeID)][_data.getReverseID(dataEdgeID)] = false;
        }
    }

    return true;
}

bool CandidateSpace::filter() {
    const EdgeID *queryOffset = _query.getOffsets();
    const VertexID *queryNbr = _query.getNbrs();
    const VertexID *dataNbr = _data.getNbrs();
    std::vector<int> filteredTimes(_query.getNumVertices(), 0);
    std::vector<double> penalties(_query.getNumVertices(), 0.5);
    std::priority_queue<FilterVertex> vq;
    for (VertexID u = 0; u < _query.getNumVertices(); u++)
        vq.emplace(0.5, 0, u);
    ui selectedDegreeSum = 0;
    ui maxSum = 5 * _query.getNumVertices();
    std::vector<int> vertexStep(_query.getNumVertices(), 0);
    int currentStep = 0;
    while (!vq.empty() && selectedDegreeSum <= maxSum) {
        FilterVertex fv = vq.top();
        vq.pop();
        double penalty = fv.penalty;
        int filteredTime = fv.filteredTime;
        VertexID u = fv.u;
        if (filteredTime < vertexStep[u]) continue;
        ++currentStep;
        selectedDegreeSum += _query.getDegree(u);
        int size1 = candidateSet[u].size();
        for (int i = 0; i < candidateSet[u].size(); ++i) {
            VertexID v = candidateSet[u][i];
            if (!filter(u, v)) {
                for (EdgeID queryEdge = queryOffset[u]; queryEdge < queryOffset[u + 1]; ++queryEdge) {
                    VertexID u2 = queryNbr[queryEdge];
                    LabelID l2 = _query.getVertexLabel(u2);
                    ui numNbr;
                    EdgeID startEdge = _data.getEdgesByLabel(v, l2, numNbr);
                    for (int j = 0; j < numNbr; ++j) {
                        EdgeID dataEdge = startEdge + j;
                        if (_bitsetEdge[queryEdge][dataEdge]) {
                            _bitsetEdge[queryEdge][dataEdge] = false;
                            _bitsetEdge[_query.getReverseID(queryEdge)][_data.getReverseID(dataEdge)] = false;
                        }
                    }
                }
                candidateSet[u][i] = candidateSet[u].back();
                candidateSet[u].pop_back();
                --i;
                _bitsetVertex[u][v] = false;
            }
        }
        int size2 = candidateSet[u].size();
        penalties[u] = 1.0;
        if (size2 == size1) continue;
        for (EdgeID queryEdge = queryOffset[u]; queryEdge < queryOffset[u + 1]; ++queryEdge) {
            VertexID u2 = queryNbr[queryEdge];
            penalties[u2] = penalties[u2] * static_cast<double>(size2) / static_cast<double>(size1);
            if (penalties[u2] > MAX_PENALTY) continue;
            vertexStep[u2] = currentStep;
            vq.emplace(penalties[u2], currentStep, u2);
        }
    }

    for (VertexID u = 0; u < _query.getNumVertices(); ++u) {
        if (candidateSet[u].empty())
            return false;
    }

    return true;
}

bool CandidateSpace::filter(VertexID u, VertexID v) {
    const EdgeID *queryOffset = _query.getOffsets();
    const VertexID *queryNbr = _query.getNbrs();
    for (EdgeID queryEdge = queryOffset[u]; queryEdge < queryOffset[u + 1]; ++queryEdge) {
        VertexID u2 = queryNbr[queryEdge];
        LabelID l2 = _query.getVertexLabel(u2);
        ui numNbr;
        EdgeID startEdge = _data.getEdgesByLabel(v, l2, numNbr);
        bool valid = false;
        for (int j = 0; j < numNbr; ++j) {
            EdgeID dataEdge = startEdge + j;
            if (!_bitsetEdge[queryEdge][dataEdge]) continue;
//            if (!TriangleSafety(queryEdge, dataEdge) || !FourCycleSafety(queryEdge, dataEdge)) {
//                _bitsetEdge[queryEdge][dataEdge] = false;
//                _bitsetEdge[_query.getReverseID(queryEdge)][_data.getReverseID(dataEdge)] = false;
//                continue;
//            }
            valid = true;
        }
        if (!valid) return false;
    }

//    return edgeBipartiteSafety(u, v);
    return true;
}

void CandidateSpace::construct() {
    buildCandidateEdge(true);
    std::vector<std::set<VertexID>> deleted(_query.getNumVertices());
    // clear candidate vertices that do not have edge
    for (VertexID u = 0; u < _query.getNumVertices(); ++u) {
        ui count;
        const VertexID *neighbors = _query.getNeighbors(u, count);
        for (int i = 0; i < candidateSet[u].size(); ++i) {
            VertexID v = candidateSet[u][i];
            for (int j = 0; j < count; ++j) {
                VertexID u2 = neighbors[j];
                if (candidateEdge[u][v][u2].empty()) {
                    deleted[u].insert(v);
                    candidateSet[u][i] = candidateSet[u].back();
                    candidateSet[u].pop_back();
                    --i;
                    break;
                }
            }
        }
        std::sort(candidateSet[u].begin(), candidateSet[u].end());
    }
    // clear candidate edges for deleted vertices
    for (VertexID u = 0; u < _query.getNumVertices(); ++u) {
        ui count;
        const VertexID *neighbors = _query.getNeighbors(u, count);
        for (int i = 0; i < candidateSet[u].size(); ++i) {
            VertexID v = candidateSet[u][i];
            if (std::find(deleted[u].begin(), deleted[u].end(), v) != deleted[u].end()) {
                candidateEdge[u][v].clear();
                continue;
            }
            for (int j = 0; j < count; ++j) {
                VertexID u2 = neighbors[j];
                bool updated = false;
                for (int k = 0; k < candidateEdge[u][v][u2].size(); ++k) {
                    VertexID v2 = candidateEdge[u][v][u2][k];
                    if (std::find(deleted[u2].begin(), deleted[u2].end(), v2) != deleted[u2].end()) {
                        candidateEdge[u][v][u2][k] = candidateEdge[u][v][u2].back();
                        candidateEdge[u][v][u2].pop_back();
                        --k;
                        updated = true;
                    }
                }
                if (updated) std::sort(candidateEdge[u][v][u2].begin(), candidateEdge[u][v][u2].end());
            }
        }
    }

    // delete bitset
    for (int i = 0; i < _query.getNumEdges(); ++i) {
        delete[] _bitsetEdge[i];
    }
    for (int i = 0; i < _query.getNumVertices(); ++i) {
        delete[] _bitsetVertex[i];
    }
    delete[] _bitsetVertex;
    delete[] _bitsetEdge;
}

void CandidateSpace::buildCandidateEdge(bool checkEdge) {
    const EdgeID *queryOffset = _query.getOffsets();
    const VertexID *queryNbr = _query.getNbrs();
    const VertexID *dataNbr = _data.getNbrs();
    for (VertexID u = 0; u < _query.getNumVertices(); ++u) {
        candidateEdge[u].clear();
        candidateEdge[u].resize(_data.getNumVertices());
        for (EdgeID queryEdge = queryOffset[u]; queryEdge < queryOffset[u + 1]; ++queryEdge) {
            VertexID u2 = queryNbr[queryEdge];
            LabelID l2 = _query.getVertexLabel(u2);
            for (int i = 0; i < candidateSet[u].size(); ++i) {
                VertexID v = candidateSet[u][i];
                ui numNbr;
                EdgeID startEdge = _data.getEdgesByLabel(v, l2, numNbr);
                candidateEdge[u][v].resize(_query.getNumVertices());
                candidateEdge[u][v][u2].resize(numNbr);
                size_t index = 0;
                for (int j = 0; j < numNbr; ++j) {
                    EdgeID dataEdge = startEdge + j;
                    VertexID v2 = dataNbr[dataEdge];
                    if (!_bitsetVertex[u2][v2]) continue;
                    if (checkEdge && !_bitsetEdge[queryEdge][dataEdge]) continue;
                    candidateEdge[u][v][u2][index] = v2;
                    ++index;
                }
                candidateEdge[u][v][u2].resize(index);
            }
        }
    }
}

void CandidateSpace::writeToStream(std::ofstream &outStream) {
    for (VertexID u1 = 0; u1 < _query.getNumVertices(); u1++) {
        ui numNbr;
        const VertexID *neighbors = _query.getNeighbors(u1, numNbr);
        for (ui i = 0; i < numNbr; ++i) {
            VertexID u2 = neighbors[i];
            outStream << "# " << u1 << " " << u2 << std::endl;
            for (int j = 0; j < candidateSet[u1].size(); ++j) {
                VertexID v = candidateSet[u1][j];
                for (int k = 0; k < candidateEdge[u1][v][u2].size(); ++k) {
                    VertexID v2 = candidateEdge[u1][v][u2][k];
                    outStream << v << " " << v2 << std::endl;
                }
            }
        }
    }
}

void CandidateSpace::setQueryGraphWeights(Graph &query) {
    query.initWeights();
    const EdgeID *queryOffset = query.getOffsets();
    const VertexID *queryNbr = query.getNbrs();
    for (VertexID u = 0; u < query.getNumVertices(); ++u) {
        query.setVertexWeight(u, candidateSet[u].size());
        for (EdgeID queryEdge = queryOffset[u]; queryEdge < queryOffset[u + 1]; ++queryEdge) {
            VertexID u2 = queryNbr[queryEdge];
            ui weight = 0;
            for (int i = 0; i < candidateSet[u].size(); ++i) {
                VertexID v = candidateSet[u][i];
                weight += candidateEdge[u][v][u2].size();
            }
            query.setEdgeWeight(queryEdge, weight);
        }
    }
    query.allSourcesShortestPaths(dist, next);
}

size_t CandidateSpace::getDist(VertexID u1, VertexID u2) const {
    return dist[u1][u2];
}

std::vector<VertexID> CandidateSpace::reconstructPath(VertexID i, VertexID j) const {
    if (next[i][j] == std::numeric_limits<VertexID>::max()) return {}; // Path does not exist
    std::vector<VertexID> path = {i};
    while (i != j) {
        i = next[i][j];
        path.push_back(i);
    }
    return path;
}

ui CandidateSpace::getMaxSize() const {
    ui maxSize = 0;
    for (VertexID u = 0; u < candidateSet.size(); ++u) {
        ui size = candidateSet[u].size();
        if (size > maxSize) maxSize = size;
    }

    return maxSize;
}

bool CandidateSpace::buildCandCFL() {
    VertexID root = selectCFLRoot();
    std::vector<std::vector<VertexID>> levels, before, lowerLevelAfter, sameLevelAfter, children;
    std::vector<VertexID> order;
    buildCFLBFSTree(root, levels, order, before, lowerLevelAfter, sameLevelAfter, children);
    for (VertexID u = 0; u < _query.getNumVertices(); ++u) {
        ui count;
        _data.getVerticesByLabel(_query.getVertexLabel(u), count);
        candidateSet[u].resize(count);
    }
//     candidates for root
    LabelID rootLabel = _query.getVertexLabel(root);
    ui count;
    const VertexID *dataV = _data.getVerticesByLabel(rootLabel, count);
    const std::map<LabelID, ui> &rootNLC = _query.getNLC(root);
    ui rootDegree = _query.getDegree(root);
    std::vector<ui> candSize(_query.getNumVertices(), 0);
    for (int i = 0; i < count; ++i) {
        VertexID v = dataV[i];
        if (_data.getDegree(v) < rootDegree) continue;
        const std::map<LabelID, ui> &vNLC = _data.getNLC(v);
        if (!nlcValid(rootNLC, vNLC)) continue;
        candidateSet[root][candSize[root]++] = v;
        _bitsetVertex[root][v] = true;
    }
    const VertexID *dataNbr = _data.getNbrs();
    std::vector<ui> numVisitedParent(_data.getNumVertices(), 0);
    std::vector<ui> updatedVertices(_data.getNumVertices());
    for (int i = 0; i < levels.size(); ++i) {
        // candidates for u: 1. union neighbors of candidates 2. intersecting unions
        for (int j = 0; j < levels[i].size(); ++j) {
            VertexID u = levels[i][j];
            ui degree = _query.getDegree(u);
            LabelID l = _query.getVertexLabel(u);
            const std::map<LabelID, ui> &uNLC = _query.getNLC(u);
            ui numParent = 0, numUpdated = 0;
            for (int j2 = 0; j2 < before[u].size(); ++j2) {
                VertexID parentU = before[u][j2];
                for (VertexID parentV : candidateSet[parentU]) {
                    ui numNbr;
                    EdgeID startEdge = _data.getEdgesByLabel(parentV, l, numNbr);
                    for (int k = 0; k < numNbr; ++k) {
                        EdgeID dataEdge = startEdge + k;
                        VertexID v = dataNbr[dataEdge];
                        if (_data.getDegree(v) < degree) continue;
                        const std::map<LabelID, ui> &vNLC = _data.getNLC(v);
                        if (!nlcValid(uNLC, vNLC)) continue;
                        if (numVisitedParent[v] == numParent) {
                            numVisitedParent[v] += 1;
                            if (numVisitedParent[v] == before[u].size()) {
                                candidateSet[u][candSize[u]++] = v;
                                _bitsetVertex[u][v] = true;
                            }
                            if (numParent == 0) {
                                updatedVertices[numUpdated++] = v;
                            }
                        }
                    }
                }
                ++numParent;
            }
            for (int k = 0; k < numUpdated; ++k) {
                numVisitedParent[updatedVertices[k]] = 0;
            }
        }
        // backward pruning: intersect neighbors of v with candidates of u' in 'sameLevelAfter'
        for (int j = levels[i].size() - 1; j >= 0; --j) {
            VertexID u = levels[i][j];
            ui degree = _query.getDegree(u);
            LabelID l = _query.getVertexLabel(u);
//            const std::map<LabelID, ui> &uNLC = _query.getNLC(u);
            ui numParent = 0, numUpdated = 0;
            for (int j2 = 0; j2 < sameLevelAfter[u].size(); ++j2) {
                VertexID parentU = sameLevelAfter[u][j2];
                for (VertexID parentV : candidateSet[parentU]) {
                    ui numNbr;
                    EdgeID startEdge = _data.getEdgesByLabel(parentV, l, numNbr);
                    for (int k = 0; k < numNbr; ++k) {
                        EdgeID dataEdge = startEdge + k;
                        VertexID v = dataNbr[dataEdge];
                        if (_data.getDegree(v) < degree) continue;
//                        const std::map<LabelID, ui> &vNLC = _data.getNLC(v);
//                        if (!nlcValid(uNLC, vNLC)) continue;
                        if (numVisitedParent[v] == numParent) {
                            numVisitedParent[v] += 1;
                            if (numParent == 0) {
                                updatedVertices[numUpdated++] = v;
                            }
                        }
                    }
                }
                ++numParent;
            }
            for (VertexID v: candidateSet[u]) {
                if (!_bitsetVertex[u][v]) continue;
                if (numVisitedParent[v] != sameLevelAfter[u].size())
                    _bitsetVertex[u][v] = false;
            }
            for (int k = 0; k < numUpdated; ++k) {
                numVisitedParent[updatedVertices[k]] = 0;
            }
        }
    }
    // bottom-up pruning: intersect neighbors of v with candidates of u' in 'after'
    for (int i = levels.size() - 2; i >= 0; --i) {
        for (int j = levels[i].size() - 1; j >= 0; --j) {
            VertexID u = levels[i][j];
            LabelID l = _query.getVertexLabel(u);
            ui degree = _query.getDegree(u);
//            const std::map<LabelID, ui> &uNLC = _query.getNLC(u);
            ui numParent = 0, numUpdated = 0;
            for (int j2 = 0; j2 < lowerLevelAfter[u].size(); ++j2) {
                VertexID parentU = lowerLevelAfter[u][j2];
                for (VertexID parentV : candidateSet[parentU]) {
                    ui numNbr;
                    EdgeID startEdge = _data.getEdgesByLabel(parentV, l, numNbr);
                    for (int k = 0; k < numNbr; ++k) {
                        EdgeID dataEdge = startEdge + k;
                        VertexID v = dataNbr[dataEdge];
                        if (_data.getDegree(v) < degree) continue;
//                        const std::map<LabelID, ui> &vNLC = _data.getNLC(v);
//                        if (!nlcValid(uNLC, vNLC)) continue;
                        if (numVisitedParent[v] == numParent) {
                            numVisitedParent[v] += 1;
                            if (numParent == 0) {
                                updatedVertices[numUpdated++] = v;
                            }
                        }
                    }
                }
                ++numParent;
            }
            for (VertexID v: candidateSet[u]) {
                if (!_bitsetVertex[u][v]) continue;
                if (numVisitedParent[v] != lowerLevelAfter[u].size())
                    _bitsetVertex[u][v] = false;
            }
            for (int k = 0; k < numUpdated; ++k) {
                numVisitedParent[updatedVertices[k]] = 0;
            }
        }
    }
    // compact
    for (VertexID u = 0; u < _query.getNumVertices(); ++u) {
        ui pos = 0;
        for (int i = 0; i < candSize[u]; ++i) {
            VertexID v = candidateSet[u][i];
            if (_bitsetVertex[u][v]) {
                candidateSet[u][pos++] = v;
            }
        }
        candidateSet[u].resize(pos);
        std::sort(candidateSet[u].begin(), candidateSet[u].end());
    }
    for (VertexID u = 0; u < _query.getNumVertices(); ++u) {
        if (candidateSet[u].empty())
            return false;
    }
    // build edges
    buildCandidateEdge(false);
    // delete bitset
    for (int i = 0; i < _query.getNumEdges(); ++i) {
        delete[] _bitsetEdge[i];
    }
    for (int i = 0; i < _query.getNumVertices(); ++i) {
        delete[] _bitsetVertex[i];
    }
    delete[] _bitsetVertex;
    delete[] _bitsetEdge;

    return true;
}

VertexID CandidateSpace::selectCFLRoot() {
    auto compare = [](const std::pair<VertexID, double> &l, const std::pair<VertexID, double> &r) {
        return l.second < r.second;
    };
    std::priority_queue<std::pair<VertexID, double>, std::vector<std::pair<VertexID, double>>, decltype(compare)> rankQ(compare);
    for (VertexID u = 0; u < _query.getNumVertices(); ++u) {
        ui degree = _query.getDegree(u);
        LabelID label = _query.getVertexLabel(u);
        ui count;
        _data.getVerticesByLabel(label, count);
        double rank = count / (double) degree;
        rankQ.emplace(u, rank);
    }
    while (rankQ.size() > 3) rankQ.pop();
    VertexID root = 0;
    double minScore = _data.getNumVertices();
    while (!rankQ.empty()) {
        VertexID u = rankQ.top().first;
        const std::map<LabelID, ui> &uNLC = _query.getNLC(u);
        ui degree = _query.getDegree(u);
        LabelID label = _query.getVertexLabel(u);
        ui count;
        const VertexID *cand = _data.getVerticesByLabel(label, count);
        ui num = 0;
        for (int i = 0; i < count; ++i) {
            VertexID v = cand[i];
            if (_data.getDegree(v) < degree) continue;
            const std::map<LabelID, ui> &vNLC = _data.getNLC(v);
            if (!nlcValid(uNLC, vNLC)) continue;
            ++num;
        }
        double score = num / (double) degree;
        if (score < minScore) {
            root = u;
            minScore = score;
        }
        rankQ.pop();
    }

    return root;
}

void
CandidateSpace::buildCFLBFSTree(VertexID root, std::vector<std::vector<VertexID>> &levels, std::vector<VertexID> &order,
                                std::vector<std::vector<VertexID>> &before, std::vector<std::vector<VertexID>> &lowerLevelAfter,
                                std::vector<std::vector<VertexID>> &sameLevelAfter,
                                std::vector<std::vector<VertexID>> &children) {
    std::queue<VertexID> q;
    q.push(root);
    std::vector<bool> visited(_query.getNumVertices(), false);
    visited[root] = true;
    before = std::vector<std::vector<VertexID>>(_query.getNumVertices());
    lowerLevelAfter = std::vector<std::vector<VertexID>>(_query.getNumVertices());
    sameLevelAfter = std::vector<std::vector<VertexID>>(_query.getNumVertices());
    children = std::vector<std::vector<VertexID>>(_query.getNumVertices());
    while (!q.empty()) {
        std::vector<VertexID> level;
        while (!q.empty()) {
            VertexID u = q.front();
            level.push_back(u);
            order.push_back(u);
            q.pop();
        }
        levels.push_back(level);
        for (int i = 0; i < level.size(); ++i) {
            VertexID u1 = level[i];
            for (int j = i + 1; j < level.size(); ++j) {
                VertexID u2 = level[j];
                if (_query.getEdgeID(u1, u2) != -1) sameLevelAfter[u1].push_back(u2);
            }
            ui count;
            const VertexID *neighbors = _query.getNeighbors(u1, count);
            for (int j = 0; j < count; ++j) {
                if (!visited[neighbors[j]]) {
                    visited[neighbors[j]] = true;
                    q.push(neighbors[j]);
                    children[u1].push_back(neighbors[j]);
                }
            }
        }
    }
    for (int i = 0; i < order.size(); ++i) {
        VertexID u1 = order[i];
        for (int j = 0; j < i; ++j) {
            VertexID u2 = order[j];
            if (_query.getEdgeID(u1, u2) != -1) before[u1].push_back(u2);
        }
        for (int j = i + 1; j < order.size(); ++j) {
            VertexID u2 = order[j];
            if (std::find(sameLevelAfter[u1].begin(), sameLevelAfter[u1].end(), u2) != sameLevelAfter[u1].end()) continue;
            if (_query.getEdgeID(u1, u2) != -1) lowerLevelAfter[u1].push_back(u2);
        }
    }
}
