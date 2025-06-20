//
// Created by Qiyan LI on 2024/2/27.
//

#include "graph.h"

Graph::Graph() {
    _numVertices = 0;
    _numEdges = 0;
    _numLabels = 0;
    _labels = nullptr;
    _offsets = nullptr;
    _reverseID = nullptr;
    _nbrs = nullptr;
    _nlc = nullptr;
    _labelOffsets = nullptr;
    _offset2 = nullptr;
    _verticesByLabel = nullptr;
    _nodeWeight = nullptr;
    _edgeWeight = nullptr;
}

Graph::~Graph() {
    delete[] _labels;
    delete[] _offsets;
    delete[] _reverseID;
    delete[] _nbrs;
    delete[] _nlc;
    delete[] _labelOffsets;
    delete[] _offset2;
    delete[] _verticesByLabel;
    delete[] _nodeWeight;
    delete[] _edgeWeight;
}

void Graph::buildReverseID() {
    _reverseID = new EdgeID[_numEdges];
    for (VertexID v1 = 0; v1 < _numVertices; ++v1) {
        for (EdgeID e12 = _offsets[v1]; e12 < _offsets[v1 + 1]; ++e12) {
            VertexID v2 = _nbrs[e12];
            if (v2 < v1) continue;
            EdgeID e21 = getUndirectedEID(v2, v1);
            _reverseID[e12] = e21;
            _reverseID[e21] = e12;
        }
    }
}

void Graph::buildNLC() {
    _nlc = new std::map<LabelID, ui>[_numVertices];
    _labelOffsets = new std::map<LabelID, EdgeID>[_numVertices];
    for (ui i = 0; i < _numVertices; ++i) {
        ui count;
        const VertexID * neighbors = getNeighbors(i, count);
        for (ui j = 0; j < count; ++j) {
            VertexID u = neighbors[j];
            LabelID label = getVertexLabel(u);
            if (_nlc[i].find(label) == _nlc[i].end()) _nlc[i][label] = 1;
            else _nlc[i][label] += 1;
            if (_labelOffsets[i].find(label) == _labelOffsets[i].end()) {
                EdgeID e = _offsets[i] + j;
                _labelOffsets[i][label] = e;
            }
        }
    }
}

void Graph::loadGraphFromTextFile(const std::string &file) {
    std::ifstream infile(file);

    if (!infile.is_open()) {
        std::cout << "Can not open file " << file << " ." << std::endl;
        exit(-1);
    }

    char type;
    infile >> type >> _numVertices >> _numEdges;
    _numEdges *= 2;
    _offsets = new ui[_numVertices +  1];
    _offsets[0] = 0;

    _nbrs = new VertexID[_numEdges];
    _labels = new LabelID[_numVertices];
    _numLabels = 0;

    LabelID max_label_id = 0;
    std::vector<ui> _nbrsoffset(_numVertices, 0);

    while (infile >> type) {
        if (type == 'v') { // Read vertex.
            VertexID id;
            LabelID  label;
            ui degree;
            infile >> id >> label >> degree;

            _labels[id] = label;
            _offsets[id + 1] = _offsets[id] + degree;

            if (_labelCount.find(label) == _labelCount.end()) {
                _labelCount[label] = 0;
                if (label > max_label_id)
                    max_label_id = label;
            }

            _labelCount[label] += 1;
        }
        else if (type == 'e') { // Read edge.
            VertexID begin;
            VertexID end;
            infile >> begin >> end;

            ui offset = _offsets[begin] + _nbrsoffset[begin];
            _nbrs[offset] = end;

            offset = _offsets[end] + _nbrsoffset[end];
            _nbrs[offset] = begin;

            _nbrsoffset[begin] += 1;
            _nbrsoffset[end] += 1;
        }
    }
    infile.close();
    _numLabels = max_label_id + 1;
    for (ui i = 0; i < _numVertices; ++i) {
        std::sort(_nbrs + _offsets[i], _nbrs + _offsets[i + 1],
                  [this](const VertexID &a, const VertexID &b) {
                      if (_labels[a] == _labels[b])
                          return a < b;
                      return _labels[a] < _labels[b];
                  });
    }
    _verticesByLabel = new VertexID[_numVertices];
    _offset2 = new ui[_numLabels + 1];
    _offset2[0] = 0;
    _verticesByLabel[0] = 0;
    ui total = 0;
    for (LabelID l = 0; l < _numLabels; ++l) {
        _offset2[l + 1] = total;
        total += _labelCount[l];
    }
    for (VertexID u = 0; u < _numVertices; ++u) {
        LabelID label = _labels[u];
        _verticesByLabel[_offset2[label + 1]++] = u;
    }

    buildNLC();
    buildReverseID();
}

EdgeID Graph::getEdgeID(VertexID v, VertexID w) const {
    if (v > w) std::swap(v, w);
    LabelID l = _labels[w];
    auto iter = _nlc[v].find(l);
    if (iter == _nlc[v].end()) return -1;
    EdgeID labelOffset = _labelOffsets[v][l];
    int low = (int)labelOffset;
    int high = low + int(iter->second) - 1;
    int mid;

    while (high - low >= 16) {
        mid = low + ((high - low) >> 1);
#ifndef __APPLE__
        _mm_prefetch((char *) &_nbrs[(mid + 1 + high) / 2], _MM_HINT_T0);
        _mm_prefetch((char *) &_nbrs[(mid - 1 + low) / 2], _MM_HINT_T0);
#endif
        if (_nbrs[mid] == w) return mid;
        if (_nbrs[mid] > w) high = mid - 1;
        else low = mid + 1;
    }

    for (int i = low; i <= high; ++i) {
        if (_nbrs[i] == w) return i;
    }

    return -1;
}

EdgeID Graph::getUndirectedEID(VertexID v, VertexID w) const {
    LabelID l = _labels[w];
    auto iter = _nlc[v].find(l);
    if (iter == _nlc[v].end()) return -1;
    EdgeID labelOffset = _labelOffsets[v][l];
    int low = (int)labelOffset;
    int high = low + int(iter->second) - 1;
    int mid;

    while (high - low >= 16) {
        mid = low + ((high - low) >> 1);
#ifndef __APPLE__
        _mm_prefetch((char *) &_nbrs[(mid + 1 + high) / 2], _MM_HINT_T0);
        _mm_prefetch((char *) &_nbrs[(mid - 1 + low) / 2], _MM_HINT_T0);
#endif
        if (_nbrs[mid] == w) return mid;
        if (_nbrs[mid] > w) high = mid - 1;
        else low = mid + 1;
    }

    for (int i = low; i <= high; ++i) {
        if (_nbrs[i] == w) return i;
    }

    return -1;
}

void Graph::indexTriangles(std::map<LabelID, Triangle *> *&triangle, std::map<LabelID, ui> *&numTriangle) {
    ui totalNum = 0;
    triangle = new std::map<LabelID, Triangle *>[_numEdges];
    numTriangle = new std::map<LabelID, ui>[_numEdges];
    for (VertexID v1 = 0; v1 < _numVertices; ++v1) {
        for (EdgeID e12 = _offsets[v1]; e12 < _offsets[v1 + 1]; ++e12) {
            VertexID v2 = _nbrs[e12];
            std::map<LabelID, Triangle *> tri;
            for (EdgeID e23 = _offsets[v2]; e23 < _offsets[v2 + 1]; ++e23) {
                VertexID v3 = _nbrs[e23];
                LabelID l3 = _labels[v3];
                if (v3 == v1) continue;
                EdgeID e13 = getUndirectedEID(v1, v3);
                if (e13 != -1) {
                    if (numTriangle[e12].find(l3) == numTriangle[e12].end()) {
                        tri[l3] = new Triangle[_labelCount[l3]];
                        numTriangle[e12][l3] = 0;
                    }
                    tri[l3][numTriangle[e12][l3]].secondEdge = e23;
                    tri[l3][numTriangle[e12][l3]].thirdEdge = e13;
                    ++numTriangle[e12][l3];
                    ++totalNum;
                }
            }
            for (auto it : tri) {
                LabelID l3 = it.first;
                triangle[e12][l3] = new Triangle[numTriangle[e12][l3]];
                memcpy(triangle[e12][l3], tri[l3], sizeof(Triangle) * numTriangle[e12][l3]);
                delete[] it.second;
            }
        }
    }
    std::cout << "total num triangle " << totalNum << std::endl;
}

void Graph::indexFourCycles(std::map<LabelID, FourCycle *> *&fourCycle, std::map<LabelID, ui> *&numFourCycle) {
    ui totalNum = 0;
    fourCycle = new std::map<LabelID, FourCycle *>[_numEdges];
    numFourCycle = new std::map<LabelID, ui>[_numEdges];
    for (VertexID v1 = 0; v1 < _numVertices; ++v1) {
        for (EdgeID e12 = _offsets[v1]; e12 < _offsets[v1 + 1]; ++e12) {
            VertexID v2 = _nbrs[e12];
            ui num = 0;
            std::map<LabelID, FourCycle *> fc;
            for (EdgeID e23 = _offsets[v2]; e23 < _offsets[v2 + 1]; ++e23) {
                VertexID v3 = _nbrs[e23];
                LabelID l3 = _labels[v3];
                if (v3 == v1) continue;
                for (EdgeID e34 = _offsets[v3]; e34 < _offsets[v3 + 1]; ++e34) {
                    VertexID v4 = _nbrs[e34];
                    LabelID l4 = _labels[v4];
                    if (v4 == v1 || v4 == v2) continue;
                    EdgeID e14 = getUndirectedEID(v1, v4);
                    if (e14 != -1) {
                        LabelID l = l3 * MAX_NUM_LABEL + l4;
                        if (numFourCycle[e12].find(l) == numFourCycle[e12].end()) {
                            fc[l] = new FourCycle[_numEdges];
                            numFourCycle[e12][l] = 0;
                        }
                        fc[l][numFourCycle[e12][l]].secondEdge = e23;
                        fc[l][numFourCycle[e12][l]].thirdEdge = e34;
                        fc[l][numFourCycle[e12][l]].fourthEdge = e14;
                        fc[l][numFourCycle[e12][l]].oneThreeEdge = getUndirectedEID(v1, v3);
                        fc[l][numFourCycle[e12][l]].twoFourEdge = getUndirectedEID(v2, v4);
                        ++numFourCycle[e12][l];
                        ++totalNum;
                    }
                }
            }
            for (auto it : fc) {
                LabelID l = it.first;
                fourCycle[e12][l] = new FourCycle[numFourCycle[e12][l]];
                memcpy(fourCycle[e12][l], fc[l], sizeof(FourCycle) * numFourCycle[e12][l]);
                delete[] it.second;
            }
            if (totalNum > MAX_FOURCYCLE_NUM) {
                std::cout << "too many four cycles. should not index them" << std::endl;
                delete[] fourCycle;
                delete[] numFourCycle;
                fourCycle = nullptr;
                return;
            }
        }
    }

    std::cout << "total num four cycle " << totalNum << std::endl;
}

void Graph::readFromStream(std::ifstream &inFile) {
    if (!inFile.is_open()) {
        std::cout << "Can not open the binary graph file." << std::endl;
        exit(-1);
    }
    inFile.read(reinterpret_cast<char*>(&_numVertices), sizeof(_numVertices));
    inFile.read(reinterpret_cast<char*>(&_numEdges), sizeof(_numEdges));
    inFile.read(reinterpret_cast<char*>(&_numLabels), sizeof(_numLabels));
    readArrayFromStream(inFile, _labels, _numVertices);
    readArrayFromStream(inFile, _offsets, _numVertices + 1);
    readArrayFromStream(inFile, _reverseID, _numEdges);
    readArrayFromStream(inFile, _nbrs, _numEdges);
    _nlc = new std::map<LabelID, ui>[_numVertices];
    for (VertexID u = 0; u < _numVertices; ++u)
        readMapFromStream(inFile, _nlc[u]);
    _labelOffsets = new std::map<LabelID, EdgeID>[_numVertices];
    for (VertexID u = 0; u < _numVertices; ++u)
        readMapFromStream(inFile, _labelOffsets[u]);
    readMapFromStream(inFile, _labelCount);
    readArrayFromStream(inFile, _offset2, _numLabels + 1);
    _offset2[0] = 0;
    readArrayFromStream(inFile, _verticesByLabel, _numVertices);
}

void Graph::writeToStream(std::ofstream &outFile) const {
    outFile.write(reinterpret_cast<const char*>(&_numVertices), sizeof(_numVertices));
    outFile.write(reinterpret_cast<const char*>(&_numEdges), sizeof(_numEdges));
    outFile.write(reinterpret_cast<const char*>(&_numLabels), sizeof(_numLabels));
    writeArrayToStream(outFile, _labels, _numVertices);
    writeArrayToStream(outFile, _offsets, _numVertices + 1);
    writeArrayToStream(outFile, _reverseID, _numEdges);
    writeArrayToStream(outFile, _nbrs, _numEdges);
    for (VertexID u = 0; u < _numVertices; ++u)
        writeMapToStream(outFile, _nlc[u]);
    for (VertexID u = 0; u < _numVertices; ++u)
        writeMapToStream(outFile, _labelOffsets[u]);
    writeMapToStream(outFile, _labelCount);
    writeArrayToStream(outFile, _offset2, _numLabels + 1);
    writeArrayToStream(outFile, _verticesByLabel, _numVertices);
}

void Graph::initWeights() {
    _nodeWeight = new ui[_numVertices];
    memset(_nodeWeight, 0, sizeof(ui) * _numVertices);
    _edgeWeight = new ui[_numEdges];
    memset(_edgeWeight, 0, sizeof(ui) * _numEdges);
}

void Graph::allSourcesShortestPaths(std::vector<std::vector<size_t>> &dist, std::vector<std::vector<VertexID>> &next) const {
    ui n = _numVertices;
    dist = std::vector<std::vector<size_t>>(n, std::vector<size_t>(n, std::numeric_limits<size_t>::max()));
    next = std::vector<std::vector<VertexID>>(n, std::vector<VertexID>(n, std::numeric_limits<VertexID>::max()));

    for (ui i = 0; i < n; ++i) {
        dist[i][i] = std::numeric_limits<size_t>::max();
        next[i][i] = i; // Self-loops are initialized to the vertex itself
    }

    // Initialize distances and paths for direct edges
    for (VertexID v = 0; v < n; ++v) {
        for (EdgeID e = _offsets[v]; e < _offsets[v + 1]; ++e) {
            VertexID u = _nbrs[e];
            dist[v][u] = getEdgeWeight(e);
            next[v][u] = u;
        }
    }

    // Floyd-Warshall algorithm to update distances and paths
    for (ui k = 0; k < n; ++k) {
        for (ui i = 0; i < n; ++i) {
            for (ui j = 0; j < n; ++j) {
                if (dist[i][k] == std::numeric_limits<size_t>::max() || dist[k][j] == std::numeric_limits<size_t>::max()) continue;
                if (dist[i][k] * dist[k][j] < dist[i][j]) {
                    dist[i][j] = dist[i][k] * dist[k][j];
                    next[i][j] = next[i][k];
                }
            }
        }
    }
}

void Graph::computeConnectedComponents(const std::vector<VertexID> &vertices, std::vector<std::vector<VertexID>> &components) const {
    std::set<VertexID> visited;
    for (VertexID v : vertices) {
        if (visited.find(v) == visited.end()) {
            std::vector<VertexID> component = bfsForCC(v, visited, vertices);
            components.push_back(component);
        }
    }
}

bool Graph::isConnected(const std::vector<VertexID> &vertices) const {
    std::set<VertexID> visited;
    return bfsForCC(vertices[0], visited, vertices).size() == vertices.size();
}

std::vector<VertexID> Graph::bfsForCC(VertexID start, set<VertexID> &visited, const std::vector<VertexID>& vertices) const {
    std::vector<VertexID> component;
    std::queue<VertexID> queue;
    queue.push(start);

    while (!queue.empty()) {
        VertexID current = queue.front();
        queue.pop();

        if (visited.find(current) == visited.end()) {
            visited.insert(current);
            component.push_back(current);

            // Iterate through all neighbors of current vertex
            for (EdgeID i = _offsets[current]; i < _offsets[current + 1]; ++i) {
                VertexID nbr = _nbrs[i];
                // Check if this neighbor is part of the subgraph we are interested in
                if (visited.find(nbr) == visited.end() && std::find(vertices.begin(), vertices.end(), nbr) != vertices.end()) {
                    queue.push(nbr);
                }
            }
        }
    }

    return component;
}

void Graph::buildHyperGraph(HyperG &h) const {
    h.N = _numVertices;
    h.M = _numEdges / 2;
    h.e.clear();
    for (VertexID u = 0; u < _numVertices; ++u) {
        for (EdgeID e = _offsets[u]; e < _offsets[u + 1]; ++e) {
            VertexID u2 = _nbrs[e];
            if (u2 < u) continue;
            VertexSet temp;
            temp.Set(u);
            temp.Set(u2);
            h.e.push_back(temp);
        }
    }
}

void Graph::buildFHD(FHD &fhd) const {
    HyperG H;
    buildHyperGraph(H);
    HyperG H_old = H;
    Order prefix_o;
    std::map<size_t, size_t> Vres_map;
    for(size_t i = 0; i < H.N; ++i)
        Vres_map[i] = i;
    Preprocessing(H, prefix_o, Vres_map);
    Order elim_o;
    double ans = DPFHD(H, elim_o);
    for(size_t i = 0; i < elim_o.size(); ++i)
        elim_o[i] = Vres_map[elim_o[i]];
    elim_o.insert(elim_o.begin(), prefix_o.begin(), prefix_o.end());
    H = H_old;
    fhd = FHD(H, elim_o);
    fhd.Refine();
}

void Graph::writeToTextFile(const std::string &filename) const {
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        throw std::runtime_error("Unable to open file");
    }
    outFile << "t " << _numVertices << " " << _numEdges / 2 << "\n";
    for (ui i = 0; i < _numVertices; ++i) {
        ui degree = _offsets[i + 1] - _offsets[i];
        outFile << "v " << i << " " << _labels[i] << " " << degree << "\n";
    }
    for (ui i = 0; i < _numVertices; ++i) {
        for (ui j = _offsets[i]; j < _offsets[i + 1]; ++j) {
            ui neighbor = _nbrs[j];
            if (i < neighbor) {
                outFile << "e " << i << " " << neighbor << "\n";
            }
        }
    }

    outFile.close();
}
