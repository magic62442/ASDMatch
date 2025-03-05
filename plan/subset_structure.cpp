//
// Created by anonymous authors on 2024/9/19.
//

#include "subset_structure.h"

// Function to generate all subsets of a given size
void generateSubsetsK(const std::vector<VertexID> &set, int size, int index, std::vector<VertexID> &current,
                      vector<uint64_t> &subsets, const std::vector<uint64_t> &cc, const Graph &query, bool connected) {
    if (size == 0) {
        if (connected) {
            if (!subsetConnectivity(query, cc, current))
                return;
        }
        subsets.push_back(getSubsetID(current));
        return;
    }
    for (int i = index; i <= set.size() - size; ++i) {
        current.push_back(set[i]);
        generateSubsetsK(set, size - 1, i + 1, current, subsets, cc, query, connected);
        current.pop_back();
    }
}

void
generateSubsets(const std::vector<VertexID> &set, vector<uint64_t> &subsets) {
    size_t n = set.size();
    uint64_t totalSubsets = 1ULL << n;
    for (uint64_t subsetID = 1; subsetID < totalSubsets; ++subsetID) {
        std::vector<VertexID> subset;
        for (size_t i = 0; i < n; ++i) {
            // Check if the i-th bit is set in subsetID
            if (subsetID & (1ULL << i)) {
                subset.push_back(set[i]);
            }
        }
        subsets.push_back(getSubsetID(subset));
    }
}

void generateExtendedSubsets(const std::vector<VertexID>& set, uint64_t subsetID, std::vector<uint64_t>& extendedSubsets) {
    size_t n = set.size();
    uint64_t totalSetID = getSubsetID(set);
    uint64_t remainingSet = totalSetID & ~subsetID;
    uint64_t totalRemainingSubsets = 1ULL << __builtin_popcountll(remainingSet);
    for (uint64_t additionalSubsetID = 0; additionalSubsetID < totalRemainingSubsets; ++additionalSubsetID) {
        uint64_t extendedSubset = subsetID;
        uint64_t mask = remainingSet;
        int bitIndex = 0;
        while (mask) {
            if (additionalSubsetID & (1ULL << bitIndex)) {
                extendedSubset |= mask & (-mask); // Add the least significant set bit from remainingSet
            }
            mask &= (mask - 1);  // Clear the least significant set bit in remainingSet
            ++bitIndex;
        }

        extendedSubsets.push_back(extendedSubset);
    }
}

VertexID findExtraElement(uint64_t idSubsetK, uint64_t idSubsetKMinus1) {
    // Compute XOR to find the differing bit
    uint64_t differingBits = idSubsetK ^ idSubsetKMinus1;

    // Find the position of the differing bit
    int position = 0;
    while (differingBits > 1) {
        differingBits >>= 1;
        position++;
    }

    // Return the element corresponding to the position
    return position;
}

ui subsetSize(uint64_t subsetID, VertexID n) {
    ui count = 0;
    for (int i = 0; i < n; ++i) {
        count += subsetID & 1;
        subsetID >>= 1;
    }
    return count;
}

// Function to convert an ID to a subset
std::vector<VertexID> getSubsetFromID(uint64_t id, ui maxVID) {
    std::vector<VertexID> subset;
    for (VertexID i = 0; i <= maxVID; ++i) {
        if (id & (1ULL << i)) {
            subset.push_back(i);
        }
    }
    return subset;
}

void smallCard(const Graph &query, CandidateSpace &cs) {
    for (VertexID u = 0; u < query.getNumVertices(); ++u) {
        uint64_t id = 1ULL << u;
        subsetToCard[id] = double(cs.candidateSet[u].size());
        for (VertexID u2 = 0; u2 < query.getNumVertices(); ++u2) {
            if (u == u2) continue;
            uint64_t edgeID = id + (1ULL << u2);
            if (subsetToCard.find(edgeID) != subsetToCard.end()) continue;
            if (query.getEdgeID(u, u2) != -1) {
                double card = 0.0;
                for (VertexID v: cs.candidateSet[u]) {
                    card += cs.candidateEdge[u][v][u2].size();
                }
                subsetToCard[edgeID] = card;
            }
        }
    }
}

bool isSubset(uint64_t subset1ID, uint64_t subset2ID) {
    return (subset1ID & subset2ID) == subset1ID;
}

SubsetStructure::SubsetStructure(const std::vector<VertexID> &elements, const Graph &query, bool connected)
        : elements(elements), n(elements.size()), subsets(n + 1) {
    std::sort(this->elements.begin(), this->elements.end());
    std::vector<std::vector<VertexID>> components;
    query.computeConnectedComponents(elements, components);
    for (auto & component : components) {
        cc.push_back(getSubsetID(component));
    }
    for (int k = 1; k <= n; ++k) {
        std::vector<VertexID> current;
        generateSubsetsK(elements, k, 0, current, subsets[k], cc, query, connected);
        std::sort(subsets[k].begin(), subsets[k].end());
    }
    for (int i = 0; i < n + 1; ++i) {
        subsetToCost.emplace_back(subsets[i].size(), 0);
        remainingCost.emplace_back(subsets[i].size(), 0);
        orders.emplace_back(subsets[i].size(), std::vector<VertexID>());
        remainingOrders.emplace_back(subsets[i].size(), std::vector<VertexID>());
    }
    generateSubsetSupersetRelationships();
}

void SubsetStructure::generateSubsetSupersetRelationships() {
    subsetSupersetRelationships(n, elements, subsets, subsetOf, supersetOf);
}

void SubsetStructure::optimalPlanDP(const Graph &query, CandidateSpace &cs, bool *visited, VertexID *partMatch,
                                    VertexID **candidates, ui *candCount) {
    for (int i = 0; i < subsets[1].size(); ++i) {
        orders[1][i].push_back(elements[i]);
        subsetToCost[1][i] = subsetToCard[1 << elements[i]];
    }
    for (int level = 2; level < n + 1; ++level) {
        for (int i = 0; i < subsets[level].size(); ++i) {
            uint64_t id = subsets[level][i];
            // traverse subsets
            int minPos = 0;
            double minCost = std::numeric_limits<double>::max();
            for (int j = 0; j < subsetOf[level][i].size(); ++j) {
                int pos = subsetOf[level][i][j];
                uint64_t childID = subsets[level - 1][pos];
                VertexID u = findExtraElement(id, childID);
                std::vector<VertexID> subsetContent = getSubsetFromID(childID, query.getNumVertices());
                double cost = subsetToCost[level - 1][pos];
                double listSize = 0.0;
                bool cartesian = true;
                for (VertexID u2: subsetContent) {
                    if (query.getEdgeID(u2, u) != -1) {
                        cartesian = false;
                        uint64_t edgeID = (1 << u) + (1 << u2);
                        listSize += subsetToCard[edgeID] / subsetToCard[1 << u2];
                    }
                }
                if (cartesian) {
                    VertexID cartesianParent = subsetContent[0];
                    size_t minDist = std::numeric_limits<size_t>::max();
                    for (VertexID u2: subsetContent) {
                        if (cs.dist[u][u2] < minDist) {
                            minDist = cs.dist[u][u2];
                            cartesianParent = u2;
                        }
                    }
                    uint64_t id2 = (1 << u) + (1 << cartesianParent);
                    if (subsetToCard.find(id2) == subsetToCard.end()) {
                        subsetToCard[id2] = estimateCartesian(cartesianParent, u, cs);
                    }
                    listSize = subsetToCard[id2] / subsetToCard[1 << cartesianParent];
                }
                cost += subsetToCard[childID] * listSize;
                if (cost < minCost) {
                    minPos = pos;
                    minCost = cost;
                }
            }
            uint64_t childID = subsets[level - 1][minPos];
            VertexID u = findExtraElement(id, childID);
            subsetToCost[level][i] = minCost;
            orders[level][i] = orders[level -1][minPos];
            orders[level][i].push_back(u);
            if (subsetToCard.find(id) == subsetToCard.end()) {
                std::vector<std::vector<VertexID>> vertexParents(level);
                std::vector<VertexID> cartesianParent(level);
                for (int j = 0; j < level; ++j) {
                    VertexID u2 = orders[level][i][j];
                    for (int k = 0; k < j; ++k) {
                        VertexID u3 = orders[level][i][k];
                        if (query.getEdgeID(u2, u3) != -1)
                            vertexParents[j].push_back(u3);
                    }
                }
                for (int j = 1; j < level; ++j) {
                    VertexID u2 = orders[level][i][j];
                    if (vertexParents[j].empty()) {
                        size_t minDist = std::numeric_limits<size_t>::max();
                        for (int k = 0; k < j; ++k) {
                            VertexID u3 = orders[level][i][k];
                            if (cs.dist[u3][u2] < minDist) {
                                minDist = cs.dist[u3][u2];
                                cartesianParent[j] = u3;
                            }
                        }
                    }
                }
                std::vector<ui> poses(query.getNumVertices(), 0);
                std::vector<double> estimations(level, 0.0);
                cardEstimateLabeled(orders[level][i], vertexParents, cartesianParent, cs, visited, partMatch, candidates,
                                    candCount, poses, estimations);
                double estimation = estimations.back();
                if (estimation == 0) estimation = 1;
                subsetToCard[id] = estimation;
            }
        }
    }
}

void SubsetStructure::reverseDP(const Graph &query, CandidateSpace &cs) {
    remainingOrders[n][0] = {};
    remainingCost[n][0] = 0;
    for (int level = n - 1; level >= 1; --level) {
        for (int i = 0; i < subsets[level].size(); ++i) {
            uint64_t id = subsets[level][i];
            // traverse subsets
            int minPos = 0;
            double minCost = std::numeric_limits<double>::max();
            for (int j = 0; j < supersetOf[level][i].size(); ++j) {
                int pos = supersetOf[level][i][j];
                uint64_t parentID = subsets[level + 1][pos];
                VertexID u = findExtraElement(parentID, id);
                std::vector<VertexID> subsetContent = getSubsetFromID(id, query.getNumVertices());
                double cost = remainingCost[level + 1][pos];
                double listSize = 0.0;
                bool cartesian = true;
                for (VertexID u2 : subsetContent) {
                    if (query.getEdgeID(u2, u) != -1) {
                        cartesian = false;
                        uint64_t edgeID = (1 << u) + (1 << u2);
                        listSize += subsetToCard[edgeID] / subsetToCard[1 << u2];
                    }
                }
                if (cartesian) {
                    VertexID cartesianParent = subsetContent[0];
                    size_t minDist = std::numeric_limits<size_t>::max();
                    for (VertexID u2: subsetContent) {
                        if (cs.dist[u][u2] < minDist) {
                            minDist = cs.dist[u][u2];
                            cartesianParent = u2;
                        }
                    }
                    uint64_t id2 = (1 << u) + (1 << cartesianParent);
                    if (subsetToCard.find(id2) == subsetToCard.end()) {
                        subsetToCard[id2] = estimateCartesian(cartesianParent, u, cs);
                    }
                    listSize = subsetToCard[id2] / subsetToCard[1 << cartesianParent];
                }
                cost += subsetToCard[id] * listSize;
                if (cost < minCost) {
                    minPos = pos;
                    minCost = cost;
                }
            }
            uint64_t parentID = subsets[level + 1][minPos];
            VertexID u = findExtraElement(parentID, id);
            remainingCost[level][i] = minCost;
            remainingOrders[level][i] = {u};
            for (VertexID u2 : remainingOrders[level + 1][minPos]) remainingOrders[level][i].push_back(u2);
        }
    }
}

void SubsetStructure::prefixPlan(uint64_t prefixID, ui k, std::vector<VertexID> &optOrder, double &optCost) const {
    if (k == 0) {
        optOrder = getOptOrder();
        optCost = getOptCost();
        return;
    }
    int pos = getPosition(prefixID, k);
    optOrder = remainingOrders[k][pos];
    optCost = remainingCost[k][pos];
}

void SubsetStructure::extendPrefixPlan(const Graph &query, CandidateSpace &cs, uint64_t id1, ui k1, uint64_t id2,
                                       ui k2, std::vector<VertexID> &optOrder, double &optCost) const {
    if (k1 == 0) {
        optOrder = getOptOrder(id2, k2);
        optCost = getOptCost(id2, k2);
        return;
    }
    SubsetStructure newS = buildSubStructure(*this, id1, k1, id2, k2);
    for (int level = k1 + 1; level < k2 + 1; ++level) {
        for (int i = 0; i < newS.subsets[level].size(); ++i) {
            uint64_t id = newS.subsets[level][i];
            // traverse subsets
            int minPos = 0;
            double minCost = std::numeric_limits<double>::max();
            for (int j = 0; j < newS.subsetOf[level][i].size(); ++j) {
                int pos = newS.subsetOf[level][i][j];
                uint64_t childID = newS.subsets[level - 1][pos];
                VertexID u = findExtraElement(id, childID);
                std::vector<VertexID> subsetContent = getSubsetFromID(childID, query.getNumVertices());
                double cost = newS.subsetToCost[level - 1][pos];
                double listSize = 0.0;
                bool cartesian = true;
                for (VertexID u2: subsetContent) {
                    if (query.getEdgeID(u2, u) != -1) {
                        cartesian = false;
                        uint64_t edgeID = (1 << u) + (1 << u2);
                        listSize += subsetToCard[edgeID] / subsetToCard[1 << u2];
                    }
                }
                if (cartesian) {
                    VertexID cartesianParent = subsetContent[0];
                    size_t minDist = std::numeric_limits<size_t>::max();
                    for (VertexID u2: subsetContent) {
                        if (cs.dist[u][u2] < minDist) {
                            minDist = cs.dist[u][u2];
                            cartesianParent = u2;
                        }
                    }
                    uint64_t id2 = (1 << u) + (1 << cartesianParent);
                    if (subsetToCard.find(id2) == subsetToCard.end()) {
                        subsetToCard[id2] = estimateCartesian(cartesianParent, u, cs);
                    }
                    listSize = subsetToCard[id2] / subsetToCard[1 << cartesianParent];
                }
                cost += subsetToCard[childID] * listSize;
                if (cost < minCost) {
                    minPos = pos;
                    minCost = cost;
                }
            }
            uint64_t childID = newS.subsets[level - 1][minPos];
            VertexID u = findExtraElement(id, childID);
            newS.subsetToCost[level][i] = minCost;
            newS.orders[level][i] = newS.orders[level -1][minPos];
            newS.orders[level][i].push_back(u);
        }
    }

    optOrder = newS.getOptOrder();
    optCost = newS.getOptCost();
}

std::vector<VertexID>
SubsetStructure::getOptOrder(uint64_t prefixID, ui k1, uint64_t subsetID, ui k2, const Graph &query,
                             CandidateSpace &cs) const {
    SubsetStructure newStructure;
    newStructure.elements = elements;
    newStructure.n = n;
    newStructure.subsets.resize(n + 1);
    newStructure.subsetOf.resize(n + 1);
    newStructure.supersetOf.resize(n + 1);
    newStructure.orders.resize(n + 1);
    newStructure.subsetToCost.resize(n + 1);
    std::set<uint64_t> visited;
    std::queue<std::pair<uint64_t, int>> subsetQ;
    subsetQ.emplace(prefixID, k1);
    visited.insert(prefixID);
    while (!subsetQ.empty()) {
        uint64_t id = subsetQ.front().first;
        int num = subsetQ.front().second;
        subsetQ.pop();
        int originalPos = getPosition(id, num);
        newStructure.subsets[num].push_back(id);
        if (num < k2) {
            for (auto pPos: supersetOf[num][originalPos]) {
                uint64_t pID = subsets[num + 1][pPos];
                if (visited.find(pID) == visited.end() && isSubset(pID, subsetID)) {
                    subsetQ.emplace(pID, num + 1);
                    visited.insert(pID);
                }
            }
        }
    }
    for (ui i = k1; i <= k2; ++i) {
        std::sort(newStructure.subsets[i].begin(), newStructure.subsets[i].end());
    }
    newStructure.generateSubsetSupersetRelationships();
    for (ui i = k1; i <= k2; ++i) {
        newStructure.subsetToCost[i].resize(newStructure.subsets[i].size());
        newStructure.orders[i].resize(newStructure.subsets[i].size());
    }
    newStructure.subsetToCost[k1][0] = 0.0;
    for (int level = k1 + 1; level < k2 + 1; ++level) {
        for (int i = 0; i < newStructure.subsets[level].size(); ++i) {
            uint64_t id = newStructure.subsets[level][i];
            // traverse subsets
            int minPos = 0;
            double minCost = std::numeric_limits<double>::max();
            for (int j = 0; j < newStructure.subsetOf[level][i].size(); ++j) {
                int pos = newStructure.subsetOf[level][i][j];
                uint64_t childID = newStructure.subsets[level - 1][pos];
                VertexID u = findExtraElement(id, childID);
                std::vector<VertexID> subsetContent = getSubsetFromID(childID, query.getNumVertices());
                double cost = newStructure.subsetToCost[level - 1][pos];
                double listSize = 0.0;
                bool cartesian = true;
                for (VertexID u2: subsetContent) {
                    if (query.getEdgeID(u2, u) != -1) {
                        cartesian = false;
                        uint64_t edgeID = (1 << u) + (1 << u2);
                        listSize += subsetToCard[edgeID] / subsetToCard[1 << u];
                    }
                }
                if (cartesian) {
                    VertexID cartesianParent = subsetContent[0];
                    size_t minDist = std::numeric_limits<size_t>::max();
                    for (VertexID u2: subsetContent) {
                        if (cs.dist[u][u2] < minDist) {
                            minDist = cs.dist[u][u2];
                            cartesianParent = u2;
                        }
                    }
                    uint64_t id2 = (1 << u) + (1 << cartesianParent);
                    if (subsetToCard.find(id2) == subsetToCard.end()) {
                        subsetToCard[id2] = estimateCartesian(cartesianParent, u, cs);
                    }
                    listSize = subsetToCard[id2] / subsetToCard[1 << u];
                }
                cost += subsetToCard[childID] * listSize;
                if (cost < minCost) {
                    minPos = pos;
                    minCost = cost;
                }
            }
            uint64_t childID = newStructure.subsets[level - 1][minPos];
            VertexID u = findExtraElement(id, childID);
            newStructure.subsetToCost[level][i] = minCost;
            newStructure.orders[level][i] = newStructure.orders[level -1][minPos];
            newStructure.orders[level][i].push_back(u);
        }
    }

    return newStructure.orders[k2][0];
}

bool SubsetStructure::saveToSteam(std::ofstream &ofs) const {
    if (!ofs) return false;
    uint64_t size = elements.size();
    ofs.write(reinterpret_cast<const char*>(&size), sizeof(size));
    if (size > 0) {
        ofs.write(reinterpret_cast<const char*>(elements.data()), size * sizeof(VertexID));
    }
    ofs.write(reinterpret_cast<const char*>(&n), sizeof(n));
    size = cc.size();
    ofs.write(reinterpret_cast<const char*>(&size), sizeof(size));
    if (size > 0) {
        ofs.write(reinterpret_cast<const char*>(cc.data()), size * sizeof(uint64_t));
    }
    size = subsets.size();
    ofs.write(reinterpret_cast<const char*>(&size), sizeof(size));
    for (const auto &subset : subsets) {
        uint64_t subset_size = subset.size();
        ofs.write(reinterpret_cast<const char*>(&subset_size), sizeof(subset_size));
        if (subset_size > 0) {
            ofs.write(reinterpret_cast<const char*>(subset.data()), subset_size * sizeof(uint64_t));
        }
    }
    size = subsetOf.size();
    ofs.write(reinterpret_cast<const char*>(&size), sizeof(size));
    for (const auto &vec2d : subsetOf) {
        uint64_t vec2d_size = vec2d.size();
        ofs.write(reinterpret_cast<const char*>(&vec2d_size), sizeof(vec2d_size));
        for (const auto &vec : vec2d) {
            uint64_t vec_size = vec.size();
            ofs.write(reinterpret_cast<const char*>(&vec_size), sizeof(vec_size));
            if (vec_size > 0) {
                ofs.write(reinterpret_cast<const char*>(vec.data()), vec_size * sizeof(int));
            }
        }
    }
    size = supersetOf.size();
    ofs.write(reinterpret_cast<const char*>(&size), sizeof(size));
    for (const auto &vec2d : supersetOf) {
        uint64_t vec2d_size = vec2d.size();
        ofs.write(reinterpret_cast<const char*>(&vec2d_size), sizeof(vec2d_size));
        for (const auto &vec : vec2d) {
            uint64_t vec_size = vec.size();
            ofs.write(reinterpret_cast<const char*>(&vec_size), sizeof(vec_size));
            if (vec_size > 0) {
                ofs.write(reinterpret_cast<const char*>(vec.data()), vec_size * sizeof(int));
            }
        }
    }
    size = orders.size();
    ofs.write(reinterpret_cast<const char*>(&size), sizeof(size));
    for (const auto &vec2d : orders) {
        uint64_t vec2d_size = vec2d.size();
        ofs.write(reinterpret_cast<const char*>(&vec2d_size), sizeof(vec2d_size));
        for (const auto &vec : vec2d) {
            uint64_t vec_size = vec.size();
            ofs.write(reinterpret_cast<const char*>(&vec_size), sizeof(vec_size));
            if (vec_size > 0) {
                ofs.write(reinterpret_cast<const char*>(vec.data()), vec_size * sizeof(VertexID));
            }
        }
    }
    size = subsetToCost.size();
    ofs.write(reinterpret_cast<const char*>(&size), sizeof(size));
    for (const auto &vec : subsetToCost) {
        uint64_t vec_size = vec.size();
        ofs.write(reinterpret_cast<const char*>(&vec_size), sizeof(vec_size));
        if (vec_size > 0) {
            // 处理 double 精度，只保留两位小数
            std::vector<int64_t> scaled;
            scaled.reserve(vec_size);
            for (double d : vec) {
                scaled.push_back(static_cast<int64_t>(d * 100.0));
            }
            ofs.write(reinterpret_cast<const char*>(scaled.data()), vec_size * sizeof(int64_t));
        }
    }
    size = remainingOrders.size();
    ofs.write(reinterpret_cast<const char*>(&size), sizeof(size));
    for (const auto &vec2d : remainingOrders) {
        uint64_t vec2d_size = vec2d.size();
        ofs.write(reinterpret_cast<const char*>(&vec2d_size), sizeof(vec2d_size));
        for (const auto &vec : vec2d) {
            uint64_t vec_size = vec.size();
            ofs.write(reinterpret_cast<const char*>(&vec_size), sizeof(vec_size));
            if (vec_size > 0) {
                ofs.write(reinterpret_cast<const char*>(vec.data()), vec_size * sizeof(VertexID));
            }
        }
    }
    size = remainingCost.size();
    ofs.write(reinterpret_cast<const char*>(&size), sizeof(size));
    for (const auto &vec : remainingCost) {
        uint64_t vec_size = vec.size();
        ofs.write(reinterpret_cast<const char*>(&vec_size), sizeof(vec_size));
        if (vec_size > 0) {
            std::vector<int64_t> scaled;
            scaled.reserve(vec_size);
            for (double d : vec) {
                scaled.push_back(static_cast<int64_t>(d * 100.0));
            }
            ofs.write(reinterpret_cast<const char*>(scaled.data()), vec_size * sizeof(int64_t));
        }
    }

    return ofs.good();
}

bool SubsetStructure::readFromStream(std::ifstream &ifs) {
    if (!ifs) return false;
    uint64_t size;
    ifs.read(reinterpret_cast<char*>(&size), sizeof(size));
    elements.resize(size);
    if (size > 0) {
        ifs.read(reinterpret_cast<char*>(elements.data()), size * sizeof(VertexID));
    }
    ifs.read(reinterpret_cast<char*>(&n), sizeof(n));
    ifs.read(reinterpret_cast<char*>(&size), sizeof(size));
    cc.resize(size);
    if (size > 0) {
        ifs.read(reinterpret_cast<char*>(cc.data()), size * sizeof(uint64_t));
    }
    ifs.read(reinterpret_cast<char*>(&size), sizeof(size));
    subsets.resize(size);
    for (auto &subset : subsets) {
        uint64_t subset_size;
        ifs.read(reinterpret_cast<char*>(&subset_size), sizeof(subset_size));
        subset.resize(subset_size);
        if (subset_size > 0) {
            ifs.read(reinterpret_cast<char*>(subset.data()), subset_size * sizeof(uint64_t));
        }
    }
    ifs.read(reinterpret_cast<char*>(&size), sizeof(size));
    subsetOf.resize(size);
    for (auto &vec2d : subsetOf) {
        uint64_t vec2d_size;
        ifs.read(reinterpret_cast<char*>(&vec2d_size), sizeof(vec2d_size));
        vec2d.resize(vec2d_size);
        for (auto &vec : vec2d) {
            uint64_t vec_size;
            ifs.read(reinterpret_cast<char*>(&vec_size), sizeof(vec_size));
            vec.resize(vec_size);
            if (vec_size > 0) {
                ifs.read(reinterpret_cast<char*>(vec.data()), vec_size * sizeof(int));
            }
        }
    }
    ifs.read(reinterpret_cast<char*>(&size), sizeof(size));
    supersetOf.resize(size);
    for (auto &vec2d : supersetOf) {
        uint64_t vec2d_size;
        ifs.read(reinterpret_cast<char*>(&vec2d_size), sizeof(vec2d_size));
        vec2d.resize(vec2d_size);
        for (auto &vec : vec2d) {
            uint64_t vec_size;
            ifs.read(reinterpret_cast<char*>(&vec_size), sizeof(vec_size));
            vec.resize(vec_size);
            if (vec_size > 0) {
                ifs.read(reinterpret_cast<char*>(vec.data()), vec_size * sizeof(int));
            }
        }
    }
    ifs.read(reinterpret_cast<char*>(&size), sizeof(size));
    orders.resize(size);
    for (auto &vec2d : orders) {
        uint64_t vec2d_size;
        ifs.read(reinterpret_cast<char*>(&vec2d_size), sizeof(vec2d_size));
        vec2d.resize(vec2d_size);
        for (auto &vec : vec2d) {
            uint64_t vec_size;
            ifs.read(reinterpret_cast<char*>(&vec_size), sizeof(vec_size));
            vec.resize(vec_size);
            if (vec_size > 0) {
                ifs.read(reinterpret_cast<char*>(vec.data()), vec_size * sizeof(VertexID));
            }
        }
    }
    ifs.read(reinterpret_cast<char*>(&size), sizeof(size));
    subsetToCost.resize(size);
    for (auto &vec : subsetToCost) {
        uint64_t vec_size;
        ifs.read(reinterpret_cast<char*>(&vec_size), sizeof(vec_size));
        std::vector<int64_t> scaled(vec_size);
        if (vec_size > 0) {
            ifs.read(reinterpret_cast<char*>(scaled.data()), vec_size * sizeof(int64_t));
            vec.resize(vec_size);
            for (size_t i = 0; i < vec_size; ++i) {
                vec[i] = static_cast<double>(scaled[i]) / 100.0;
            }
        }
    }
    ifs.read(reinterpret_cast<char*>(&size), sizeof(size));
    remainingOrders.resize(size);
    for (auto &vec2d : remainingOrders) {
        uint64_t vec2d_size;
        ifs.read(reinterpret_cast<char*>(&vec2d_size), sizeof(vec2d_size));
        vec2d.resize(vec2d_size);
        for (auto &vec : vec2d) {
            uint64_t vec_size;
            ifs.read(reinterpret_cast<char*>(&vec_size), sizeof(vec_size));
            vec.resize(vec_size);
            if (vec_size > 0) {
                ifs.read(reinterpret_cast<char*>(vec.data()), vec_size * sizeof(VertexID));
            }
        }
    }
    ifs.read(reinterpret_cast<char*>(&size), sizeof(size));
    remainingCost.resize(size);
    for (auto &vec : remainingCost) {
        uint64_t vec_size;
        ifs.read(reinterpret_cast<char*>(&vec_size), sizeof(vec_size));
        std::vector<int64_t> scaled(vec_size);
        if (vec_size > 0) {
            ifs.read(reinterpret_cast<char*>(scaled.data()), vec_size * sizeof(int64_t));
            vec.resize(vec_size);
            for (size_t i = 0; i < vec_size; ++i) {
                vec[i] = static_cast<double>(scaled[i]) / 100.0;
            }
        }
    }

    return ifs.good();
}

void SubsetStructure::buildExtendCost(const Graph &query, CandidateSpace &cs) {
    extendCost[0][0] = 0.0;
    extendOrder[0][0] = std::vector<VertexID>();
    std::map<uint64_t, std::set<uint64_t>> allSubsets;
    for (int level = 1; level < n + 1; ++level) {
        for (int i = 0; i < subsets[level].size(); ++i) {
            uint64_t id = subsets[level][i];
            allSubsets[id].insert(id);
            allSubsets[id].insert(0);
            extendCost[id][id] = getOptCost(id, level);
            extendOrder[id][id] = std::vector<VertexID>();
            if (level == 1) {
                VertexID u = elements[i];
                extendOrder[0][id] = {u};
                extendCost[0][id] = cs.candidateSet[u].size();
                continue;
            }
            for (int j = 0; j < subsetOf[level][i].size(); ++j) {
                int pos = subsetOf[level][i][j];
                uint64_t childID = subsets[level - 1][pos];
                VertexID u = findExtraElement(id, childID);
                std::vector<VertexID> subsetContent = getSubsetFromID(childID, query.getNumVertices());
                double listSize = 0.0;
                bool cartesian = true;
                for (VertexID u2: subsetContent) {
                    if (query.getEdgeID(u2, u) != -1) {
                        cartesian = false;
                        uint64_t edgeID = (1 << u) + (1 << u2);
                        listSize += subsetToCard[edgeID] / subsetToCard[1 << u2];
                    }
                }
                if (cartesian) {
                    VertexID cartesianParent = subsetContent[0];
                    size_t minDist = std::numeric_limits<size_t>::max();
                    for (VertexID u2: subsetContent) {
                        if (cs.dist[u][u2] < minDist) {
                            minDist = cs.dist[u][u2];
                            cartesianParent = u2;
                        }
                    }
                    uint64_t id2 = (1 << u) + (1 << cartesianParent);
                    if (subsetToCard.find(id2) == subsetToCard.end()) {
                        subsetToCard[id2] = estimateCartesian(cartesianParent, u, cs);
                    }
                    listSize = subsetToCard[id2] / subsetToCard[1 << cartesianParent];
                }
                double newCost = subsetToCard[childID] * listSize;
                for (uint64_t subsetID: allSubsets[childID]) {
                    allSubsets[id].insert(subsetID);
                    double cost = newCost + extendCost[subsetID][childID];
                    if (extendCost[subsetID].find(id) == extendCost[subsetID].end() || extendCost[subsetID][id] > cost) {
                        extendCost[subsetID][id] = cost;
                        extendOrder[subsetID][id] = extendOrder[subsetID][childID];
                        extendOrder[subsetID][id].push_back(u);
                    }
                }
            }
        }
    }
}

// Function to build the sub-structure of direct and indirect supersets
SubsetStructure buildSubStructure(const SubsetStructure& original, uint64_t subsetID, ui k) {
    if (k == 0) return original;
    std::set<uint64_t> visited;
    // Construct the new sub-structure
    SubsetStructure newStructure;
    newStructure.elements = original.elements;
    ui n = original.n;
    newStructure.n = n;
    newStructure.subsets.resize(n + 1);
    newStructure.subsetOf.resize(n + 1);
    newStructure.supersetOf.resize(n + 1);
    newStructure.orders.resize(n + 1);
    newStructure.subsetToCost.resize(n + 1);
    newStructure.generateSubsetSupersetRelationships();
    std::queue<std::pair<uint64_t, int>> subsetQ;
    subsetQ.emplace(subsetID, k);
    visited.insert(subsetID);
    while (!subsetQ.empty()) {
        uint64_t id = subsetQ.front().first;
        int num = subsetQ.front().second;
        subsetQ.pop();
        int originalPos = original.getPosition(id, num);
        newStructure.subsets[num].push_back(id);
        for (auto pPos: original.supersetOf[num][originalPos]) {
            uint64_t pID = original.subsets[num + 1][pPos];
            if (visited.find(pID) == visited.end()) {
                subsetQ.emplace(pID, num + 1);
                visited.insert(pID);
            }
        }
    }
    for (ui i = k; i <= n; ++i) {
        std::sort(newStructure.subsets[i].begin(), newStructure.subsets[i].end());
    }
    newStructure.generateSubsetSupersetRelationships();
    for (ui i = k; i <= n; ++i) {
        newStructure.subsetToCost[i].resize(newStructure.subsets[i].size());
        newStructure.orders[i].resize(newStructure.subsets[i].size());
    }
    newStructure.subsetToCost[k][0] = 0.0;

    return newStructure;
}

SubsetStructure buildSubStructure(const SubsetStructure& original, uint64_t id1, ui k1, uint64_t id2, ui k2) {
    std::vector<VertexID> elements2 = getSubsetFromID(id2, original.elements.back());
    SubsetStructure newStructure;
    newStructure.elements = elements2;
    ui n = elements2.size();
    newStructure.n = n;
    newStructure.subsets.resize(n + 1);
    newStructure.subsetOf.resize(n + 1);
    newStructure.supersetOf.resize(n + 1);
    newStructure.orders.resize(n + 1);
    newStructure.subsetToCost.resize(n + 1);
    newStructure.generateSubsetSupersetRelationships();
    std::queue<std::pair<uint64_t, int>> subsetQ;
    subsetQ.emplace(id1, k1);
    std::set<uint64_t> visited;
    visited.insert(id1);
    while (!subsetQ.empty()) {
        uint64_t id = subsetQ.front().first;
        int num = subsetQ.front().second;
        subsetQ.pop();
        int originalPos = original.getPosition(id, num);
        newStructure.subsets[num].push_back(id);
        for (auto pPos: original.supersetOf[num][originalPos]) {
            uint64_t pID = original.subsets[num + 1][pPos];
            VertexID u = findExtraElement(pID, id);
            if (std::find(elements2.begin(), elements2.end(), u) == elements2.end()) continue;
            if (visited.find(pID) == visited.end()) {
                subsetQ.emplace(pID, num + 1);
                visited.insert(pID);
            }
        }
    }
    for (ui i = k1; i <= k2; ++i) {
        std::sort(newStructure.subsets[i].begin(), newStructure.subsets[i].end());
    }
    newStructure.generateSubsetSupersetRelationships();
    for (ui i = k1; i <= k2; ++i) {
        newStructure.subsetToCost[i].resize(newStructure.subsets[i].size());
        newStructure.orders[i].resize(newStructure.subsets[i].size());
    }
    newStructure.subsetToCost[k1][0] = 0.0;

    return newStructure;
}

void bagSetStructure::init(const std::vector<std::vector<VertexID>> &sharedAttrs,
                           const std::vector<SubsetStructure> &dpStructures, const Graph &query, bool connected) {
    ui n = query.getNumVertices();
    // costs and orders for k = 1
    for (int i = 0; i < sharedAttrs.size(); ++i) {
        uint64_t id = 1 << i;
        const SubsetStructure &s = dpStructures[i];
        std::vector<std::vector<VertexID>> localOrder(dpStructures.size());
        localOrder[i] = dpStructures[i].getOptOrder();
        intersections[id] = getSubsetID(sharedAttrs[i]);
        for (ui k = 0; k <= sharedAttrs[i].size(); ++k) {
            std::vector<VertexID> current;
            std::vector<uint64_t> subsets;
            generateSubsetsK(sharedAttrs[i], k, 0, current, subsets, s.cc, query, connected);
            for (uint64_t sharedID: subsets) {
                std::vector<VertexID> remainingOrder;
                double remainingCost;
                s.prefixPlan(sharedID, subsetSize(sharedID, query.getNumVertices()), remainingOrder, remainingCost);
                localOrder[i] = remainingOrder;
                localOrders[id][sharedID] = localOrder;
                indepCosts[id][sharedID] = remainingCost;
            }
        }
    }
}

Permutation encodePermutation(const std::vector<uint32_t>& elements) {
    std::bitset<165> encoded;
    size_t length = elements.size();

    // Encode the length in the first 5 bits
    for (size_t bit = 0; bit < 5; ++bit) {
        if (length & (1 << bit)) {
            encoded.set(bit);
        }
    }

    // Encode the elements starting from the 6th bit
    for (size_t i = 0; i < length; ++i) {
        uint32_t value = elements[i] & 0x1F;  // Mask to 5 bits
        size_t startBit = 5 + i * 5;
        for (size_t bit = 0; bit < 5; ++bit) {
            if (value & (1 << bit)) {
                encoded.set(startBit + bit);
            }
        }
    }
    return encoded;
}

std::vector<uint32_t> decodePermutation(const Permutation & encoded) {
    // Decode the length from the first 5 bits
    size_t length = 0;
    for (size_t bit = 0; bit < 5; ++bit) {
        if (encoded.test(bit)) {
            length |= (1 << bit);
        }
    }

    // Decode the elements
    std::vector<uint32_t> elements(length);
    for (size_t i = 0; i < length; ++i) {
        uint32_t value = 0;
        size_t startBit = 5 + i * 5;
        for (size_t bit = 0; bit < 5; ++bit) {
            if (encoded.test(startBit + bit)) {
                value |= (1 << bit);
            }
        }
        elements[i] = value;
    }
    return elements;
}

void
bagPermuteStructure::init(const std::vector<std::vector<VertexID>> &sharedAttrs,
                          const std::vector<SubsetStructure> &dpStructures,
                          const Graph &query, bool maxShare, bool connected) {
    ui n = query.getNumVertices();
    // costs and orders for k = 1
    for (int i = 0; i < sharedAttrs.size(); ++i) {
        uint64_t id = 1 << i;
        const SubsetStructure &s = dpStructures[i];
        std::vector<std::vector<VertexID>> localOrder(dpStructures.size());
        localOrder[i] = dpStructures[i].getOptOrder();
        intersections[id] = getSubsetID(sharedAttrs[i]);
        for (ui k = 0; k <= sharedAttrs[i].size(); ++k) {
            std::vector<VertexID> current;
            std::vector<uint64_t> subsets;
            generateSubsetsK(sharedAttrs[i], k, 0, current, subsets, s.cc, query, connected);
            for (uint64_t sharedID: subsets) {
                std::vector<VertexID> remainingOrder;
                double remainingCost;
                s.prefixPlan(sharedID, subsetSize(sharedID, query.getNumVertices()), remainingOrder, remainingCost);
                localOrder[i] = remainingOrder;
                std::vector<VertexID> sharedSet = getSubsetFromID(sharedID, query.getNumVertices());
                do {
                    if (connected && !orderConnectivity(query, sharedSet)) continue;
                    Permutation pi = encodePermutation(sharedSet);
                    localOrders[id][pi] = localOrder;
                    if (!maxShare) indepCosts[id][pi] = remainingCost;
                    else indepCosts[id][pi] = (double)remainingOrder.size();
                } while(std::next_permutation(sharedSet.begin(), sharedSet.end()));
            }
        }
    }
}
