//
// Created by Qiyan LI on 2024/2/27.
//

#include "estimator.h"

std::unordered_map<uint64_t, double> subsetToCard = std::unordered_map<uint64_t, double>();

void initPoses(const std::vector<VertexID> &order, const Graph &query, std::vector<std::vector<VertexID>> &vertexParents,
               std::vector<VertexID> &cartesianParent, const std::vector<std::vector<size_t>> &dist) {
    ui numAttributes = order.size();
    vertexParents = std::vector<std::vector<VertexID>>(numAttributes);
    cartesianParent = std::vector<VertexID>(numAttributes);
    bool cartesian = false;
    for (int i = 0; i < numAttributes; ++i) {
        VertexID u = order[i];
        for (int j = 0; j < i; ++j) {
            VertexID u2 = order[j];
            if (query.getEdgeID(u, u2) != -1)
                vertexParents[i].push_back(u2);
        }
        if (i > 0 && vertexParents[i].empty()) cartesian = true;
    }
    if (cartesian) {
        for (ui i = 1; i < numAttributes; ++i) {
            VertexID u = order[i];
            if (vertexParents[i].empty()) {
                size_t minDist = std::numeric_limits<size_t>::max();
                for (VertexID j = 0; j < i; ++j) {
                    VertexID u2 = order[j];
                    if (dist[u][u2] < minDist) {
                        minDist = dist[u][u2];
                        cartesianParent[i] = u2;
                    }
                }
            }
        }
    }
}

void
optimizedCartesianProduct(CandidateSpace &cs, VertexID v, const std::vector<VertexID> &path, VertexID *candidate,
                          ui &candCount) {
    VertexID u = path[0], uPrime = path.back();
    VertexID *result = new VertexID[cs.candidateSet[uPrime].size()];
    ui numResult = 0;
    candCount = 0;
    VertexID **pathCand = new VertexID *[path.size()];
    pathCand[0] = new VertexID[1];
    pathCand[0][0] = v;
    ui *pathCandCount = new ui[path.size()];
    ui *poses = new ui[path.size()];
    memset(pathCandCount, 0, sizeof(ui) * path.size());
    memset(poses, 0, sizeof(ui) * path.size());
    pathCandCount[0] = 1;
    int depth = 0;
    while (depth >= 0) {
        while (poses[depth] < pathCandCount[depth]) {
            VertexID u2 = path[depth];
            VertexID v2 = pathCand[depth][poses[depth]];
            ++poses[depth];
            VertexID next = path[depth + 1];
            VertexID *nextCand = cs.candidateEdge[u2][v2][next].data();
            ui nextCount = cs.candidateEdge[u2][v2][next].size();
            if (depth == path.size() - 2) {
                VertexID *temp = new VertexID[cs.candidateSet[uPrime].size()];
                VertexID i = 0, j = 0, k = 0;
                while (i < numResult && j < nextCount) {
                    if (result[i] < nextCand[j]) {
                        temp[k++] = result[i++];
                    } else if (result[i] > nextCand[j]) {
                        temp[k++] = nextCand[j++];
                    } else {
                        temp[k++] = result[i++];
                        ++j;
                    }
                }
                while (i < numResult) temp[k++] = result[i++];
                while (j < nextCount) temp[k++] = nextCand[j++];
                for (i = 0; i < k; ++i) result[i] = temp[i];
                numResult = k;
                delete[] temp;
            }
            else {
                ++depth;
                poses[depth] = 0;
                pathCand[depth] = nextCand;
                pathCandCount[depth] = nextCount;
            }
        }
        --depth;
    }

    for (int i = 0; i < numResult; ++i) {
        if (result[i] != v) {
            candidate[candCount] = result[i];
            ++candCount;
        }
    }
    delete[] result;
    delete[] pathCand[0];
    delete[] pathCand;
    delete[] pathCandCount;
    delete[] poses;
}

void
optimizedCartesianProduct(CandidateSpace &cs, VertexID v, const std::vector<VertexID> &path, DynamicArray<TrieNode *> *edgeColumn) {
    VertexID u = path[0], uPrime = path.back();
    VertexID *result = new VertexID[cs.candidateSet[uPrime].size()];
    ui numResult = 0;
    ui candCount = 0;
    VertexID **pathCand = new VertexID *[path.size()];
    pathCand[0] = new VertexID[1];
    pathCand[0][0] = v;
    ui *pathCandCount = new ui[path.size()];
    ui *poses = new ui[path.size()];
    memset(pathCandCount, 0, sizeof(ui) * path.size());
    memset(poses, 0, sizeof(ui) * path.size());
    pathCandCount[0] = 1;
    int depth = 0;
    while (depth >= 0) {
        while (poses[depth] < pathCandCount[depth]) {
            VertexID u2 = path[depth];
            VertexID v2 = pathCand[depth][poses[depth]];
            ++poses[depth];
            VertexID next = path[depth + 1];
            VertexID *nextCand = cs.candidateEdge[u2][v2][next].data();
            ui nextCount = cs.candidateEdge[u2][v2][next].size();
            if (depth == path.size() - 2) {
                VertexID *temp = new VertexID[cs.candidateSet[uPrime].size()];
                VertexID i = 0, j = 0, k = 0;
                while (i < numResult && j < nextCount) {
                    if (result[i] < nextCand[j]) {
                        temp[k++] = result[i++];
                    } else if (result[i] > nextCand[j]) {
                        temp[k++] = nextCand[j++];
                    } else {
                        temp[k++] = result[i++];
                        ++j;
                    }
                }
                while (i < numResult) temp[k++] = result[i++];
                while (j < nextCount) temp[k++] = nextCand[j++];
                for (i = 0; i < k; ++i) result[i] = temp[i];
                numResult = k;
                delete[] temp;
            }
            else {
                ++depth;
                poses[depth] = 0;
                pathCand[depth] = nextCand;
                pathCandCount[depth] = nextCount;
            }
        }
        --depth;
    }

    for (int i = 0; i < numResult; ++i) {
        if (result[i] != v) {
            edgeColumn[0][candCount]->value = result[i];
            ++candCount;
        }
    }
    edgeColumn->setSize(candCount);
    delete[] result;
    delete[] pathCand[0];
    delete[] pathCand;
    delete[] pathCandCount;
    delete[] poses;
}

void
cardEstimateLabeled(const std::vector<VertexID> &order, const std::vector<std::vector<VertexID>> &vertexParents,
                    const std::vector<VertexID> &cartesianParent, CandidateSpace &cs, bool *visited,
                    VertexID *partMatch, VertexID **candidates, ui *candCount, std::vector<ui> &poses,
                    std::vector<double> &estimation) {
    if (order.empty()) {
        estimation.push_back(0.0);
        return;
    }
    if (order.size() == 1) {
        estimation.push_back(cs.candidateSet[order[0]].size());
        return;
    }
    if (order.size() == 2) {
        VertexID u1 = order[0], u2 = order[1];
        estimation.push_back(cs.candidateSet[u1].size());
        bool connected = true;
        double sz = 0.0;
        for (VertexID v: cs.candidateSet[u1]) {
            sz += cs.candidateEdge[u1][v][u2].size();
        }
        if (sz == 0.0) sz = cs.candidateSet[u1].size() * cs.candidateSet[u2].size();
    }
    int depth = 0;
    VertexID u = order[0];
    VertexID **neighbors = new VertexID *[order.size()];
    ui maxSize = cs.getMaxSize();
    for (int i = 0; i < order.size(); ++i) {
        neighbors[i] = new VertexID[maxSize];
    }
    ui *neighborCount = new ui[order.size()];
    std::vector<ui> numSamples(order.size());
    std::vector<ui> numUsed(order.size(), 0);
    std::vector<std::vector<double>> estimate(order.size());
    for (int i = 0; i < order.size(); ++i) estimate[i] = std::vector<double>(order.size() - i, 0.0);
    numSamples[0] = SAMPLE_SIZE;
    std::vector<ui> totalCandCount(order.size(), 0);
    totalCandCount[0] = cs.candidateSet[u].size();
    // randomly sample a subset of candidates
    if (totalCandCount[0] < SAMPLE_SIZE) candCount[0] = totalCandCount[0];
    else candCount[0] = SAMPLE_SIZE;
    memcpy(candidates[0], cs.candidateSet[u].data(), cs.candidateSet[u].size() * sizeof(VertexID));
    if (candCount[0] < cs.candidateSet[u].size())
        sampleKElements(candidates[0], cs.candidateSet[u].size(), candCount[0]);
    ui totalCount = 0;
    numSamples[1] = SAMPLE_SIZE / candCount[0];
    while (depth >= 0) {
        while (poses[depth] < candCount[depth]) {
            VertexID v = candidates[depth][poses[depth]];
            ++poses[depth];
            visited[v] = true;
            partMatch[order[depth]] = v;
            ++depth;
            numUsed[depth] = 0;
            for (int i = 0; i < estimate[depth].size(); ++i) {
                estimate[depth][i] = 0.0;
            }
            poses[depth] = 0;
            u = order[depth];
            const std::vector<VertexID> &parents = vertexParents[depth];
            if (!parents.empty()) {
                for (int i = 0; i < parents.size(); ++i) {
                    VertexID pU = parents[i];
                    neighborCount[i] = 0;
                    for (int j = 0; j < cs.candidateEdge[pU][partMatch[pU]][u].size(); ++j) {
                        VertexID v2 = cs.candidateEdge[pU][partMatch[pU]][u][j];
                        if (visited[v2]) continue;
                        neighbors[i][neighborCount[i]] = v2;
                        ++neighborCount[i];
                    }
                }
                ComputeSetIntersection::LeapfrogJoin(neighbors, neighborCount, parents.size(), candidates[depth], totalCandCount[depth]);
            }
            else {
                const std::vector<VertexID> &path = cs.reconstructPath(cartesianParent[depth], u);
                optimizedCartesianProduct(cs, partMatch[cartesianParent[depth]], path, candidates[depth], totalCandCount[depth]);
            }
            if (depth == order.size() - 1) {
                --depth;
                estimate[depth][0] += 1.0;
                estimate[depth][1] += static_cast<double>(totalCandCount[depth + 1]);
                totalCount += totalCandCount[depth + 1];
                ++numUsed[depth];
                visited[partMatch[order[depth]]] = false;
            }
            else {
                // randomly sample a subset of candidates. set the number of samples
                ui sampleSize;
                sampleSize = totalCandCount[depth] / SAMPLE_PORTION_LABELED;
                sampleSize = std::min(numSamples[depth], sampleSize);
                if (totalCandCount[depth] != 0 && sampleSize == 0) sampleSize = 1;
                if (sampleSize < totalCandCount[depth])
                    sampleKElements(candidates[depth], totalCandCount[depth], sampleSize);
                candCount[depth] = sampleSize;
                if (sampleSize != 0) numSamples[depth + 1] = numSamples[depth] / sampleSize;
                else numSamples[depth + 1] = 0;
            }
        }
        --depth;
        if (depth >= 0) {
            // update estimations at the current depth
            estimate[depth][0] += 1.0;
            for (int i = 1; i < estimate[depth].size(); ++i) {
                estimate[depth][i] += static_cast<double>(estimate[depth + 1][i - 1]) * totalCandCount[depth] / candCount[depth];
            }
            numUsed[depth] += numUsed[depth + 1];
            if (poses[depth] < candCount[depth])
                numSamples[depth + 1] = (numSamples[depth] - numUsed[depth]) / (candCount[depth] - poses[depth]);
            visited[partMatch[order[depth]]] = false;
        }
    }

    estimation = estimate[0];
    for (int i = 0; i < order.size(); ++i) {
        delete[] neighbors[i];
    }
    delete[] neighbors;
    delete[] neighborCount;
}

void
cardEstimateLabeled(VertexID u, bool lastLevel, const std::vector<VertexID> &prevOrder, const Graph &query,
                    CandidateSpace &cs, std::vector<std::vector<VertexID>> &prevMatches,
                    std::vector<std::vector<VertexID>> &nextMatches, std::vector<double> &weights,
                    std::vector<double> &nextWeights, double &estimation) {
    std::vector<VertexID> attributesBefore;
    VertexID cartesianParent = 99;
    for (int i = 0 ; i < prevOrder.size(); ++i) {
        VertexID u2 = prevOrder[i];
        if (query.getEdgeID(u, u2) != -1)
            attributesBefore.push_back(u2);
    }
    if (attributesBefore.empty()) {
        size_t minDist = std::numeric_limits<size_t>::max();
        for (int i = 0; i < prevOrder.size(); ++i) {
            VertexID u2 = prevOrder[i];
            if (cs.dist[u][u2] < minDist) {
                minDist = cs.dist[u][u2];
                cartesianParent = u2;
            }
        }
    }
    if (prevMatches.empty() && !prevOrder.empty()) {
        estimation = 1.0;
        return;
    }
    VertexID *candidates = new VertexID[cs.candidateSet[u].size()];
    ui totalCandCount = 0, candCount = 0;
    if (prevMatches.empty() && prevOrder.empty()) {
        memcpy(candidates, cs.candidateSet[u].data(), cs.candidateSet[u].size() * sizeof(VertexID));
        totalCandCount = candCount = cs.candidateSet[u].size();
        if (candCount >= SAMPLE_SIZE) candCount = SAMPLE_SIZE;
        sampleKElements(candidates, cs.candidateSet[u].size(), candCount);
        for (int i = 0; i < candCount; ++i) {
            std::vector<VertexID> nextMatch(query.getNumVertices());
            nextMatch[u] = candidates[i];
            nextMatches.push_back(nextMatch);
            nextWeights.push_back(double(totalCandCount) / double(candCount));
        }
        delete[] candidates;
        return;
    }
    VertexID **neighbors = new VertexID *[prevOrder.size()];
    for (int i = 0; i < prevOrder.size(); ++i) {
        neighbors[i] = new VertexID[cs.candidateSet[u].size()];
    }
    ui *neighborCount = new ui[prevOrder.size()];
    estimation = 0.0;
    ui numSamples = SAMPLE_SIZE / prevMatches.size();
    for (int i = 0; i < prevMatches.size(); ++i) {
        std::set<VertexID> visited;
        const std::vector<VertexID> prevMatch = prevMatches[i];
        for (VertexID u2: prevOrder) visited.insert(prevMatch[u2]);
        if (!attributesBefore.empty()) {
            for (int j = 0; j < attributesBefore.size(); ++j) {
                VertexID pU = attributesBefore[j];
                neighborCount[j] = 0;
                for (int k = 0; k < cs.candidateEdge[pU][prevMatch[pU]][u].size(); ++k) {
                    VertexID v2 = cs.candidateEdge[pU][prevMatch[pU]][u][k];
                    if (visited.find(v2) != visited.end()) continue;
                    neighbors[j][neighborCount[j]] = v2;
                    ++neighborCount[j];
                }
            }
            ComputeSetIntersection::LeapfrogJoin(neighbors, neighborCount, attributesBefore.size(), candidates, totalCandCount);
        }
        else {
            const std::vector<VertexID> &path = cs.reconstructPath(cartesianParent, u);
            optimizedCartesianProduct(cs, prevMatch[cartesianParent], path, candidates, totalCandCount);
        }
        if (!lastLevel) {
            ui sampleSize = totalCandCount / SAMPLE_PORTION_LABELED;
            sampleSize = std::min(sampleSize, numSamples);
            if (totalCandCount != 0 && sampleSize == 0) sampleSize = 1;
            sampleKElements(candidates, totalCandCount, sampleSize);
            candCount = sampleSize;
            for (int j = 0; j < candCount; ++j) {
                std::vector<VertexID> nextMatch = prevMatch;
                nextMatch[u] = candidates[j];
                nextMatches.push_back(nextMatch);
                nextWeights.push_back(weights[i] * totalCandCount / candCount);
            }
            numSamples = (numSamples - sampleSize) / (prevMatches.size() - i);
        }

        estimation += weights[i] * totalCandCount;
    }

    for (int i = 0; i < prevOrder.size(); ++i) {
        delete[] neighbors[i];
    }

    delete[] neighbors;
    delete[] neighborCount;
    delete[] candidates;
}

double estimateCartesian(VertexID u1, VertexID u2, CandidateSpace &cs) {
    double card = 0.0;
    const std::vector<VertexID> &path = cs.reconstructPath(u1, u2);
    VertexID *candidates = new VertexID [cs.candidateSet[u2].size()];
    ui candCount = 0;
    // sample vertices for cartesian product
    VertexID *vSamples = new VertexID [cs.candidateSet[u1].size()];
    memcpy(vSamples, cs.candidateSet[u1].data(), cs.candidateSet[u1].size() * sizeof(VertexID));
    ui numSamples = cs.candidateSet[u1].size() / SAMPLE_PORTION_LABELED;
    if (numSamples >= SAMPLE_SIZE) numSamples = SAMPLE_SIZE;
    if (numSamples == 0) numSamples = 1;
    sampleKElements(vSamples, cs.candidateSet[u1].size(), numSamples);
    for (int i = 0; i < numSamples; ++i) {
        optimizedCartesianProduct(cs, vSamples[i], path, candidates, candCount);
        card += double(candCount) * cs.candidateSet[u1].size() / numSamples;
    }
    delete[] candidates;
    delete[] vSamples;
    if (card < cs.candidateSet[u1].size()) card = static_cast<double>(cs.candidateSet[u1].size());
    if (card < cs.candidateSet[u2].size()) card = static_cast<double>(cs.candidateSet[u2].size());
    return card;
}

std::vector<VertexID> RIOrder(const Graph &query) {
    std::vector<VertexID> order;
    VertexID maxDegreeU = 0;
    ui maxDegree = 0;
    for (VertexID u = 1; u < query.getNumVertices(); ++u) {
        if (query.getDegree(u) > maxDegree) {
            maxDegree = query.getDegree(u);
            maxDegreeU = u;
        }
    }
    order.push_back(maxDegreeU);
    while (order.size() != query.getNumVertices()) {
        std::vector<int> numBackWard(query.getNumVertices(), -1);
        std::vector<int> breakTieNum1(query.getNumVertices(), 0);
        std::vector<int> breakTieNum2(query.getNumVertices(), 0);
        std::vector<VertexID> remaining;
        for (VertexID u = 0; u < query.getNumVertices(); ++u) {
            if (std::find(order.begin(), order.end(), u) != order.end())
                continue;
            else remaining.push_back(u);
            numBackWard[u] = 0;
            for (VertexID u2 : order) {
                if (query.getEdgeID(u, u2) != -1) ++numBackWard[u];
            }
        }
        for (VertexID u: remaining) {
            for (VertexID u2: order) {
                for (VertexID u3: remaining) {
                    if (query.getEdgeID(u, u3) != -1 && query.getEdgeID(u2, u3) != -1) {
                        ++breakTieNum1[u];
                        break;
                    }
                }
            }
        }
        for (VertexID u: remaining) {
            for (VertexID u2: remaining) {
                if (query.getEdgeID(u, u2) != -1) {
                    bool notConnectedToPrev = true;
                    for (VertexID u3: order) {
                        if (query.getEdgeID(u2, u3) != -1)
                            notConnectedToPrev = false;
                    }
                    if (notConnectedToPrev) ++breakTieNum2[u2];
                }
            }
        }
        int maxNumBack = -1, maxTie1 = 0, maxTie2 = 0;
        for (int i = 0; i < numBackWard.size(); ++i) {
            if (numBackWard[i] > maxNumBack) {
                maxNumBack = numBackWard[i];
            }
        }
        VertexID nextU = remaining[0];
        for (VertexID u: remaining) {
            if (numBackWard[u] == maxNumBack && (breakTieNum1[u] > maxTie1 || breakTieNum1[u] == maxTie1 &&
            breakTieNum2[u] >= maxTie2)) {
                nextU = u;
                maxTie1 = breakTieNum1[u];
                maxTie2 = breakTieNum2[u];
            }
        }
        order.push_back(nextU);
    }

    return order;
}

std::vector<VertexID> GQLOrder(const Graph &query, const CandidateSpace &cs) {
    std::vector<VertexID> order;
    ui minSize = cs.candidateSet[0].size();
    for (VertexID u = 0; u < query.getNumVertices(); ++u) {
        if (cs.candidateSet[u].size() < minSize) minSize = cs.candidateSet[u].size();
    }
    for (VertexID u = 0; u < query.getNumVertices(); ++u) {
        if (cs.candidateSet[u].size() == minSize) {
            order.push_back(u);
            break;
        }
    }
    while (order.size() != query.getNumVertices()) {
        std::vector<VertexID> candidates;
        for (VertexID u = 0; u < query.getNumVertices(); ++u) {
            if (std::find(order.begin(), order.end(), u) != order.end()) continue;
            bool connected = false;
            for (VertexID u2: order) {
                if (query.getEdgeID(u, u2) != -1) {
                    connected = true;
                    break;
                }
            }
            if (connected) candidates.push_back(u);
        }
        minSize = cs.candidateSet[candidates[0]].size();
        VertexID nextU = candidates[0];
        for (VertexID u: candidates) {
            if (cs.candidateSet[u].size() < minSize) {
                minSize = cs.candidateSet[u].size();
                nextU = u;
            }
        }
        order.push_back(nextU);
    }

    return order;
}

std::vector<VertexID> simpleOrder(const Graph &query, const CandidateSpace &cs, const std::vector<VertexID> &vertices,
                                  const std::vector<int> &repetitions) {
    ui n = query.getNumVertices();
    std::map<VertexID, VertexPriority> vset;
    std::vector<bool> visited(query.getNumVertices(), true);
    for (VertexID u: vertices) visited[u] = false;
    std::vector<std::vector<VertexID>> components;
    query.computeConnectedComponents(vertices, components);
    std::vector<size_t> componentWeight(components.size());
    for (int i = 0; i < components.size(); ++i) {
        componentWeight[i] = 1;
        for (int j = 0; j < components[i].size(); ++j) {
            componentWeight[i] *= cs.candidateSet[components[i][j]].size();
        }
    }
    for (int i = 0; i < components.size(); ++i) {
        for (int j = i + 1; j < components.size(); ++j) {
            if (componentWeight[j] > componentWeight[i]) {
                std::swap(components[i], components[j]);
            }
        }
    }
    std::vector<ui> prefixSum(components.size());
    prefixSum[0] = 0;
    for (int i = 1; i < prefixSum.size(); ++i) {
        prefixSum[i] = prefixSum[i - 1] + components[i - 1].size();
    }
    std::vector<VertexID> order;
    std::vector<ui> degrees(query.getNumVertices(), 0);
    for (int i = 0; i < vertices.size(); ++i) {
        VertexID u = vertices[i];
        for (int j = 0; j < vertices.size(); ++j) {
            VertexID u2 = vertices[j];
            if (query.getEdgeID(u, u2) != -1) {
                ++degrees[u];
                ++degrees[u2];
            }
        }
    }
    for (int i = 0; i < components.size(); ++i) {
        const std::vector<VertexID> &component = components[i];
        for (VertexID u: component) {
            VertexPriority vp;
            vp.num = repetitions[u];
            vp.backwardNbr = 0;
            vp.degree = degrees[u];
            vp.candSize = cs.candidateSet[u].size();
            vset[u] = vp;
        }
        VertexID next;
        VertexPriority min;

        while (!vset.empty()) {
            for (auto &it : vset) {
                VertexID u = it.first;
                VertexPriority priority = it.second;
                bool connected = false;
                if (order.size() == prefixSum[i]) connected = true;
                else {
                    for (VertexID u2: order) {
                        if (query.getEdgeID(u, u2) != -1) {
                            connected = true;
                            break;
                        }
                    }
                }
                if (!connected) continue;
                if (priority < min) {
                    min = priority;
                    next = u;
                }
            }
            order.push_back(next);
            vset.erase(next);
            visited[next] = true;
            ui numNbr;
            const VertexID *neighbors = query.getNeighbors(next, numNbr);
            for (int j = 0; j < numNbr; ++j) {
                VertexID u2 = neighbors[j];
                if (!visited[u2]) ++vset[u2].backwardNbr;
            }

            min.num = min.backwardNbr = min.degree = 0;
            min.candSize = std::numeric_limits<ui>::max();
        }
    }

    return order;
}

std::vector<VertexID>
simpleOrder(const Graph &query, const CandidateSpace &cs, const std::vector<VertexID> &prefix,
            std::vector<ui> &numBackWard, const std::vector<VertexID> &prevOrder, int noChangePos) {
    std::map<VertexID, VertexPriority> vset;
    std::vector<VertexID> vertices;
    for (int i = 0; i < noChangePos; ++i) {
        vertices.push_back(prevOrder[i]);
    }
    std::vector<VertexID> order;
    if (vertices.empty()) return order;
    std::vector<bool> visited(query.getNumVertices(), true);
    for (VertexID u: vertices) visited[u] = false;
    std::vector<std::vector<VertexID>> components;
    query.computeConnectedComponents(vertices, components);
    std::vector<size_t> componentWeight(components.size());
    for (int i = 0; i < components.size(); ++i) {
        componentWeight[i] = 1;
        for (int j = 0; j < components[i].size(); ++j) {
            componentWeight[i] *= cs.candidateSet[components[i][j]].size();
        }
    }
    for (int i = 0; i < components.size(); ++i) {
        for (int j = i + 1; j < components.size(); ++j) {
            if (componentWeight[j] > componentWeight[i]) {
                std::swap(components[i], components[j]);
            }
        }
    }
    std::vector<ui> prefixSum(components.size());
    prefixSum[0] = 0;
    for (int i = 1; i < prefixSum.size(); ++i) {
        prefixSum[i] = prefixSum[i - 1] + components[i - 1].size();
    }
    std::vector<ui> degrees(query.getNumVertices(), 0);
    for (int i = 0; i < vertices.size(); ++i) {
        VertexID u = vertices[i];
        for (int j = 0; j < vertices.size(); ++j) {
            VertexID u2 = vertices[j];
            if (query.getEdgeID(u, u2) != -1) {
                ++degrees[u];
                ++degrees[u2];
            }
        }
    }
    for (int i = 0; i < components.size(); ++i) {
        const std::vector<VertexID> &component = components[i];
        for (VertexID u: component) {
            VertexPriority vp;
            vp.num = 0;
            vp.backwardNbr = 0;
            for (VertexID u2: prefix) {
                if (query.getEdgeID(u, u2) != -1)
                    ++vp.backwardNbr;
            }
            vp.degree = degrees[u];
            vp.candSize = cs.candidateSet[u].size();
            vset[u] = vp;
        }
        VertexID next;
        VertexPriority min;

        while (!vset.empty()) {
            for (auto &it : vset) {
                VertexID u = it.first;
                VertexPriority priority = it.second;
                bool connected = false;
                if (order.size() == prefixSum[i]) connected = true;
                else {
                    for (VertexID u2: order) {
                        if (query.getEdgeID(u, u2) != -1) {
                            connected = true;
                            break;
                        }
                    }
                }
                if (!connected) continue;
                if (priority < min) {
                    min = priority;
                    next = u;
                }
            }
            order.push_back(next);
            numBackWard.push_back(vset[next].backwardNbr);
            vset.erase(next);
            visited[next] = true;
            ui numNbr;
            const VertexID *neighbors = query.getNeighbors(next, numNbr);
            for (int j = 0; j < numNbr; ++j) {
                VertexID u2 = neighbors[j];
                if (!visited[u2]) ++vset[u2].backwardNbr;
            }

            min.num = min.backwardNbr = min.degree = 0;
            min.candSize = std::numeric_limits<ui>::max();
        }
    }
    for (int i = noChangePos + 1; i < prevOrder.size(); ++i) {
        order.push_back(prevOrder[i]);
    }
    return order;
}

std::vector<ui>
computeNumBackWard(const Graph &query, const std::vector<VertexID> &prefix, const std::vector<VertexID> &localOrder) {
    std::vector<ui> numBackWard(localOrder.size(), 0);
    for (int i = 0; i < localOrder.size(); ++i) {
        VertexID u = localOrder[i];
        for (int j = 0; j < prefix.size(); ++j) {
            if (query.getEdgeID(prefix[j], u) != -1)
                ++numBackWard[i];
        }
        for (int j = 0; j < i; ++j) {
            if (query.getEdgeID(localOrder[j], u) != -1)
                ++numBackWard[i];
        }
    }
    return numBackWard;
}

std::vector<double> computeCost(const std::vector<VertexID> &order, const Graph &query, CandidateSpace &cs,
                                bool *visited, VertexID *partMatch, VertexID **candidates, ui *candCount) {
    std::vector<std::vector<VertexID>> vertexParents;
    std::vector<VertexID> cartesianParent;
    initPoses(order, query, vertexParents, cartesianParent, cs.dist);
    std::vector<double> cards;
    uint64_t id = 0;
    int maxPos = -1;
    for (int i = 0; i < order.size(); ++i) {
        id += 1 << order[i];
        if (subsetToCard.find(id) == subsetToCard.end()) maxPos = i;
    }
    std::vector<VertexID> subsequence(order.begin(), order.begin() + maxPos + 1);
    std::vector<ui> poses(query.getNumVertices(), 0);
    if (visited != nullptr)
        cardEstimateLabeled(subsequence, vertexParents, cartesianParent, cs, visited, partMatch, candidates, candCount, poses, cards);
    id = 0;
    for (int i = 0; i <= maxPos; ++i) {
        id += 1 << subsequence[i];
        subsetToCard[id] = cards[i];
    }
    for (int i = maxPos + 1; i < order.size(); ++i) {
        id += 1 << order[i];
        cards.push_back(subsetToCard[id]);
    }
    std::vector<double> costs(order.size());
    for (int i = 0; i < order.size(); ++i) {
        VertexID u = order[i];
        if (i == 0) {
            costs[i] = subsetToCard[1 << u];
            continue;
        }
        double listSize = 0;
        if (vertexParents[i].empty()) {
            VertexID u2 = cartesianParent[i];
            uint64_t id2 = (1 << u) + (1 << u2);
            if (subsetToCard.find(id2) == subsetToCard.end()) {
                subsetToCard[id2] = estimateCartesian(u2, u, cs);
            }
            listSize += subsetToCard[id2] / cs.candidateSet[u2].size();
        }
        else {
            for (ui u2: vertexParents[i]) {
                uint64_t edgeID = (1 << u) + (1 << u2);
                listSize += subsetToCard[edgeID] / subsetToCard[1 << u2];
            }
        }
        costs[i] += listSize * cards[i - 1];
    }

    return costs;
}

double computeCost(const std::vector<VertexID> &order, const std::vector<std::vector<VertexID>> &vertexParents,
                   const std::vector<VertexID> &cartesianParent, const std::vector<double> &cards, CandidateSpace &cs) {
    double cost = 0.0;
    for (int i = 0; i < order.size(); ++i) {
        VertexID u = order[i];
        if (i == 0) {
            cost += cs.candidateSet[u].size();
            continue;
        }
        double listSize = 0.0;
        if (vertexParents[i].empty()) {
            VertexID u2 = cartesianParent[i];
            uint64_t id2 = (1 << u) + (1 << u2);
            if (subsetToCard.find(id2) == subsetToCard.end()) {
                subsetToCard[id2] = estimateCartesian(u2, u, cs);
            }
            listSize = subsetToCard[id2] / cs.candidateSet[u2].size();
        }
        else {
            for (ui u2: vertexParents[i]) {
                uint64_t edgeID = (1 << u) + (1 << u2);
                listSize += subsetToCard[edgeID] / subsetToCard[1 << u2];
            }
        }
        cost += listSize * cards[i - 1];
    }

    return cost;
}

// order should extend prefix. The cost does not include the prefix cost
double computeCost(const std::vector<VertexID> &prefix, const std::vector<VertexID> &order, const Graph &query,
                   CandidateSpace &cs, bool *visited, VertexID *partMatch, VertexID **candidates, ui *candCount) {
    double cost = 0.0;
    uint64_t prefixID = 0;
    int maxCardPos = 0;
    for (VertexID u : prefix) prefixID += 1 << u;
    uint64_t id = prefixID;
    for (int i = 0; i < order.size(); ++i) {
        id += 1 << order[i];
        if (subsetToCard.find(id) == subsetToCard.end())
            maxCardPos = i;
    }
    std::vector<std::vector<VertexID>> vertexParents;
    std::vector<VertexID> cartesianParent;
    initPoses(order, query, vertexParents, cartesianParent, cs.dist);
    std::vector<VertexID> subsequence(maxCardPos);
    for (int i = 0; i < maxCardPos; ++i) subsequence[i] = order[i];
    std::vector<double> estimation;
    std::vector<ui> poses(query.getNumVertices(), 0);
    if (!subsequence.empty())
        cardEstimateLabeled(subsequence, vertexParents, cartesianParent, cs, visited, partMatch, candidates, candCount, poses, estimation);
    id = prefixID;
    for (int i = 0; i < subsequence.size(); ++i) {
        id += 1 << order[i];
        if (subsetToCard.find(id) == subsetToCard.end()) subsetToCard[id] = estimation[i];
    }
    id = prefixID;
    std::vector<VertexID> nodeOrder = prefix;
    for (VertexID u: order) nodeOrder.push_back(u);
    initPoses(nodeOrder, query, vertexParents, cartesianParent, cs.dist);
    for (int i = 0; i < order.size(); ++i) {
        double card = 1.0;
        if (id != 0) card = subsetToCard[id];
        VertexID u = order[i];
        if (prefix.empty() && i == 0) {
            cost += cs.candidateSet[u].size();
            id += 1 << order[i];
            continue;
        }
        double listSize = std::numeric_limits<double>::max();
        if (vertexParents[i + prefix.size()].empty()) {
            VertexID u2 = cartesianParent[i + prefix.size()];
            uint64_t id2 = (1 << u) + (1 << u2);
            if (subsetToCard.find(id2) == subsetToCard.end()) {
                subsetToCard[id2] = estimateCartesian(u2, u, cs);
            }
            listSize = subsetToCard[id2] / cs.candidateSet[u2].size();
        }
        else {
            for (VertexID u2: vertexParents[i + prefix.size()]) {
                uint64_t edgeID = (1 << u) + (1 << u2);
                listSize += subsetToCard[edgeID] / subsetToCard[1 << u2];
            }
        }
        cost += listSize * card;
        id += 1 << order[i];
    }

    return cost;
}

double computeCost(const PrefixNode *pt, const std::vector<std::vector<VertexID>> &localOrders, const Graph &query,
                   CandidateSpace &cs, bool *visited, VertexID *partMatch, VertexID **candidates, ui *candCount,
                   const vector<vector<VertexID>> &matchedAttrs, const std::vector<VertexID> &prevAttrs,
                   const std::vector<double> &factors) {
    // pt is the virtual root
    std::vector<VertexID> attrsInPath;
    int depth = 0;
    const PrefixNode *pn = pt;
    ui height = pn -> getHeight();
    std::vector<ui> childPoses(height, 0);
    std::vector<PrefixNode *> nodes(height, nullptr);
    // generate root-to-leaf paths
    std::vector<std::vector<VertexID>> paths;
    std::vector<bool> computed(query.getNumVertices(), false);
    double cost = 0.0;
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            ++childPoses[depth];
            VertexID u = current -> u;
            attrsInPath.push_back(u);
            nodes[depth] = current;
            if (!current -> children.empty()) {
                ++depth;
                childPoses[depth] = 0;
            }
            else {
                std::vector<VertexID> path = prevAttrs;
                for (VertexID u2 : attrsInPath) path.push_back(u2);
                paths.push_back(path);
                if (childPoses[depth] < pn -> children.size())
                    attrsInPath.pop_back();
            }
            if (depth > 0) pn = nodes[depth - 1];
        }
        --depth;
        if (depth >= 0) {
            if (depth == 0) pn = pt;
            else pn = nodes[depth - 1];
            attrsInPath.pop_back();
            if (childPoses[depth] < pn -> children.size()) attrsInPath.pop_back();
        }
    }
    for (const auto &path : paths) {
        std::vector<double> costs = computeCost(path, query, cs, visited, partMatch, candidates, candCount);
        for (int i = 0; i < costs.size(); ++i) {
            VertexID u = path[i];
            if (!computed[u]) {
                computed[u] = true;
                cost += costs[i];
            }
        }
    }
    pn = pt;
    depth = 0;
    childPoses = std::vector<ui>(height, 0);
    attrsInPath.clear();
    for (VertexID nID : pt -> nIDsToCall)
        cost += computeCost(matchedAttrs[nID], localOrders[nID], query, cs, visited, partMatch, candidates, candCount) * factors[nID];
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            ++childPoses[depth];
            VertexID u = current -> u;
            attrsInPath.push_back(u);
            nodes[depth] = current;
            for (VertexID nID: current -> nIDsToCall) {
                std::vector<VertexID> prefix = matchedAttrs[nID];
                for (VertexID u2 : attrsInPath) prefix.push_back(u2);
                cost += computeCost(prefix, localOrders[nID], query, cs, visited, partMatch, candidates, candCount) * factors[nID];
            }
            if (!current -> children.empty()) {
                ++depth;
                childPoses[depth] = 0;
            }
            else if (childPoses[depth] < pn -> children.size())
                attrsInPath.pop_back();
            if (depth > 0) pn = nodes[depth - 1];
        }
        --depth;
        if (depth >= 0) {
            if (depth == 0) pn = pt;
            else pn = nodes[depth - 1];
            attrsInPath.pop_back();
            if (childPoses[depth] < pn -> children.size()) attrsInPath.pop_back();
        }
    }

    return cost;
}

ui maxNumBackWard(const std::vector<VertexID> &order, const Graph &query) {
    ui maxNum = 1;
    for (int i = 1; i < order.size(); ++i) {
        ui num = 0;
        VertexID u = order[i];
        for (int j = 0; j < i; ++j) {
            VertexID u2 = order[j];
            if (query.getEdgeID(u, u2) != -1) ++num;
        }

        if (maxNum < num) maxNum = num;
    }

    return maxNum;
}

bool saveSubsetToCard(std::ofstream& ofs) {
    if (!ofs) return false;
    uint64_t map_size = subsetToCard.size();
    ofs.write(reinterpret_cast<const char*>(&map_size), sizeof(map_size));
    for (const auto& pair : subsetToCard) {
        ofs.write(reinterpret_cast<const char*>(&pair.first), sizeof(pair.first));
        double scaled = static_cast<double>(static_cast<int64_t>(pair.second * 100.0)) / 100.0;
        ofs.write(reinterpret_cast<const char*>(&scaled), sizeof(scaled));
    }

    return ofs.good();
}

bool loadSubsetToCard(std::ifstream& ifs) {
    if (!ifs) return false;
    subsetToCard.clear();
    subsetToCard[0] = 0;
    uint64_t map_size;
    ifs.read(reinterpret_cast<char*>(&map_size), sizeof(map_size));
    for (uint64_t i = 0; i < map_size; ++i) {
        uint64_t key;
        double value;
        ifs.read(reinterpret_cast<char*>(&key), sizeof(key));
        ifs.read(reinterpret_cast<char*>(&value), sizeof(value));
        subsetToCard[key] = value;
    }

    return ifs.good();
}