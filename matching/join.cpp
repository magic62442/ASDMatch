//
// Created by anonymous authors on 2024/3/8.
//

#include "join.h"

size_t gNumResult = 0;
size_t gNumCartesian = 0;
size_t gNumTraverse = 0;

ui binarySearch(const DynamicArray<TrieNode *> &child, const ui begin, const ui end, const ui target) {
    ui offset_begin = begin;
    ui offset_end = end;
    while (offset_end - offset_begin >= 16) {
        auto mid = (offset_begin + offset_end) / 2;
#ifndef __APPLE__
        _mm_prefetch((char *) &child[(mid + 1 + offset_end) / 2], _MM_HINT_T0);
        _mm_prefetch((char *) &child[(mid - 1 + offset_begin) / 2], _MM_HINT_T0);
#endif
        if (child[mid]->value == target) {
            return mid;
        } else if (child[mid]->value < target) {
            offset_begin = mid + 1;
        } else {
            offset_end = mid;
        }
    }

    // linear search fallback
    for (auto offset = offset_begin; offset < offset_end; ++offset) {
        if (child[offset]->value >= target) {
            return offset;
        }
    }

    return offset_end;
}

// for the amortized cost O(1 + log(N/m) in the leapfrog paper, we first run an exponential search
// then run a binary search
ui leapfrogSeek(const DynamicArray<TrieNode *> &child, ui begin, ui end, ui target) {
    int bound = 1, halfBound;
    ui offset = begin + bound;
    while (offset < end && child[offset]->value < target) {
        bound *= 2;
        offset = begin + bound;
    }
    halfBound = bound / 2;
    if (offset + bound < end) end = offset + bound;
    return binarySearch(child, begin + halfBound, end, target);
}

void leapFrogJoin(DynamicArray<TrieNode *> **children, ui num, ui **iters, ui &iterSize) {
    if (num == 1) {
        for (int i = 0; i < (*children[0]).size(); ++i)
            iters[i][0] = i;
        iterSize = (*children[0]).size();
        return;
    }
#ifdef COLLECT_STATISTICS
    ++gNumInterSection;
#endif
    iterSize = 0;
    for (int i = 0; i < num; ++i) {
        if (children[i]->empty())
            return;
    }
    // Create a vector of pairs for sorting
    std::vector<std::pair<DynamicArray<TrieNode*>*, ui>> childPairs(num);
    for (ui i = 0; i < num; ++i) {
        childPairs[i] = std::make_pair(children[i], i); // Pair with original index
    }

    // Sort children by the first value
    std::sort(childPairs.begin(), childPairs.end(), [](const auto &a, const auto &b) {
        return (*a.first)[0]->value < (*b.first)[0]->value;
    });
    ui *iter = new ui[num];
    memset(iter, 0, sizeof(ui) * num);
    int p = 0;
    VertexID xPrime = (*childPairs[num - 1].first)[0]->value;
    while (true) {
        VertexID x = (*childPairs[p].first)[iter[p]]->value;
        if (x == xPrime) {
            for (ui j = 0; j < num; ++j) {
                // Map sorted iter positions back to original indices
                iters[iterSize][childPairs[j].second] = iter[j];
            }
            ++iterSize;
            ++iter[p];
        }
        else {
            iter[p] = leapfrogSeek(*childPairs[p].first, iter[p], (*childPairs[p].first).size(), xPrime);
        }
        if (iter[p] == (*childPairs[p].first).size()) break;
        xPrime = (*childPairs[p].first)[iter[p]]->value;
        p = (++p) % num;
    }
    delete[] iter;
}

/*
 * scope-style node join algorithm
 * for a given partition specified by the partial match, join all remaining attributes and store them in tau
 * order: attribute order, vertexParents[u]: N(u) that are before u in order
 * reverse: the reverse mapping of the attribute order
 * */

void
nodeJoinScope(const HyperTree &t, VertexID nID, CandidateSpace &cs, TrieNode *root, bool *visited, VertexID *partMatch,
              VertexID **candidates, ui *candCount, std::vector<ui> &poses, bool skip) {
    const HyperNode &tau = t.nodes[nID];
    const VertexID *order = tau.attributes;
    int mappingSize = int(tau.prefixSize);
    if (mappingSize == tau.numAttributes) {
        root->addMatch(partMatch, tau.attributes, tau.numAttributes, 0, nullptr, nullptr, 0);
        return;
    }
    ui compressionSize = t.compressionSizes[nID];
    const std::vector<std::vector<VertexID>> &vertexParents = tau.attributesBefore;
    // handle the first level
    VertexID u = order[mappingSize];
    VertexID **neighbors = new VertexID *[tau.numAttributes];
    ui *neighborCount = new ui[tau.numAttributes];
    if (mappingSize == 0) {
        memcpy(candidates[0], cs.candidateSet[u].data(), cs.candidateSet[u].size() * sizeof(VertexID));
        candCount[0] = cs.candidateSet[u].size();
    }
    else {
        if (mappingSize + compressionSize == tau.numAttributes) {
            // compute compressions
            VertexID **cData = new VertexID *[compressionSize];
            ui *length = new VertexID[compressionSize];
            for (int i = 0; i < compressionSize; ++i) {
                u = order[mappingSize + i];
                cData[i] = new VertexID[cs.candidateSet[u].size()];
                const std::vector<VertexID> &parents = vertexParents[mappingSize + i];
                for (int j = 0; j < parents.size(); ++j) {
                    VertexID pU = parents[j];
                    neighbors[j] = cs.candidateEdge[pU][partMatch[pU]][u].data();
                    neighborCount[j] = cs.candidateEdge[pU][partMatch[pU]][u].size();
                }
                ComputeSetIntersection::LeapfrogJoin(neighbors, neighborCount, parents.size(), cData[i], length[i]);
            }
            // add match and compressions to trie
            if (!skip) root->addMatch(partMatch, tau.attributes, tau.numAttributes, 0, cData, length, compressionSize);
            else ++gNumResult;
            delete[] cData;
            delete[] length;
            return;
        }
        const std::vector<VertexID> &parents = vertexParents[mappingSize];
        if (!parents.empty()) {
            for (int i = 0; i < parents.size(); ++i) {
                VertexID pU = parents[i];
                neighbors[i] = cs.candidateEdge[pU][partMatch[pU]][u].data();
                neighborCount[i] = cs.candidateEdge[pU][partMatch[pU]][u].size();
            }
            ComputeSetIntersection::LeapfrogJoin(neighbors, neighborCount, parents.size(), candidates[mappingSize], candCount[mappingSize]);
        }
        else {
            VertexID cartesianParent = tau.cartesianParent[mappingSize];
            const std::vector<VertexID> &path = cs.reconstructPath(cartesianParent, u);
            optimizedCartesianProduct(cs, partMatch[cartesianParent], path, candidates[mappingSize], candCount[mappingSize]);
        }
    }
    int depth = mappingSize;
    while (depth >= mappingSize) {
        while (poses[depth] < candCount[depth]) {
            VertexID v = candidates[depth][poses[depth]];
            ++poses[depth];
            if (visited[v]) continue;
            visited[v] = true;
            partMatch[order[depth]] = v;
            if (depth + compressionSize + 1 == tau.numAttributes) {
                // compute compressions
                VertexID **cData = new VertexID *[compressionSize];
                ui *length = new VertexID[compressionSize];
                for (int i = 0; i < compressionSize; ++i) {
                    u = order[depth + i + 1];
                    cData[i] = new VertexID[cs.candidateSet[u].size()];
                    const std::vector<VertexID> &parents = vertexParents[depth + i + 1];
                    for (int j = 0; j < parents.size(); ++j) {
                        VertexID pU = parents[j];
                        neighbors[j] = cs.candidateEdge[pU][partMatch[pU]][u].data();
                        neighborCount[j] = cs.candidateEdge[pU][partMatch[pU]][u].size();
                    }
                    ComputeSetIntersection::LeapfrogJoin(neighbors, neighborCount, parents.size(), cData[i], length[i]);
                }
                // add match and compressions to trie
                if (!skip) root->addMatch(partMatch, tau.attributes, tau.numAttributes, 0, cData, length, compressionSize);
                else ++gNumResult;
                delete[] cData;
                delete[] length;
                visited[partMatch[order[depth]]] = false;
            }
            else {
                ++depth;
                poses[depth] = 0;
                u = order[depth];
                const std::vector<VertexID> &parents = vertexParents[depth];
                if (!parents.empty()) {
                    for (int i = 0; i < parents.size(); ++i) {
                        VertexID pU = parents[i];
                        neighbors[i] = cs.candidateEdge[pU][partMatch[pU]][u].data();
                        neighborCount[i] = cs.candidateEdge[pU][partMatch[pU]][u].size();
                    }
                    ComputeSetIntersection::LeapfrogJoin(neighbors, neighborCount, parents.size(), candidates[depth], candCount[depth]);
                }
                else {
                    VertexID cartesianParent = tau.cartesianParent[depth];
                    const std::vector<VertexID> &path = cs.reconstructPath(cartesianParent, u);
                    optimizedCartesianProduct(cs, partMatch[cartesianParent], path, candidates[depth], candCount[depth]);
                }
            }
        }
        --depth;
        if (depth >= mappingSize) visited[partMatch[order[depth]]] = false;
    }

    delete[] neighbors;
    delete[] neighborCount;
}

void nodeJoinScope(const HyperTree &t, VertexID nID, CandidateSpace &cs, bool *visited, VertexID *partMatch,
                   VertexID **candidates, ui *candCount, std::vector<ui> &poses) {
    const HyperNode &tau = t.nodes[nID];
    const VertexID *order = tau.attributes;
    int mappingSize = int(tau.prefixSize);
    if (mappingSize == tau.numAttributes) {
        return;
    }
    ui compressionSize = t.compressionSizes[nID];
    const std::vector<std::vector<VertexID>> &vertexParents = tau.attributesBefore;
    // handle the first level
    VertexID u = order[mappingSize];
    VertexID **neighbors = new VertexID *[tau.numAttributes];
    ui *neighborCount = new ui[tau.numAttributes];
    if (mappingSize == 0) {
        memcpy(candidates[0], cs.candidateSet[u].data(), cs.candidateSet[u].size() * sizeof(VertexID));
        candCount[0] = cs.candidateSet[u].size();
    }
    else {
        if (mappingSize + compressionSize == tau.numAttributes) {
            // compute compressions
            VertexID **cData = new VertexID *[compressionSize];
            ui *length = new VertexID[compressionSize];
            for (int i = 0; i < compressionSize; ++i) {
                u = order[mappingSize + i];
                cData[i] = new VertexID[cs.candidateSet[u].size()];
                const std::vector<VertexID> &parents = vertexParents[mappingSize + i];
                for (int j = 0; j < parents.size(); ++j) {
                    VertexID pU = parents[j];
                    neighbors[j] = cs.candidateEdge[pU][partMatch[pU]][u].data();
                    neighborCount[j] = cs.candidateEdge[pU][partMatch[pU]][u].size();
                }
                ComputeSetIntersection::LeapfrogJoin(neighbors, neighborCount, parents.size(), cData[i], length[i]);
            }
            delete[] cData;
            delete[] length;
            return;
        }
        const std::vector<VertexID> &parents = vertexParents[mappingSize];
        if (!parents.empty()) {
            for (int i = 0; i < parents.size(); ++i) {
                VertexID pU = parents[i];
                neighbors[i] = cs.candidateEdge[pU][partMatch[pU]][u].data();
                neighborCount[i] = cs.candidateEdge[pU][partMatch[pU]][u].size();
            }
            ComputeSetIntersection::LeapfrogJoin(neighbors, neighborCount, parents.size(), candidates[mappingSize], candCount[mappingSize]);
        }
        else {
            VertexID cartesianParent = tau.cartesianParent[mappingSize];
            const std::vector<VertexID> &path = cs.reconstructPath(cartesianParent, u);
            optimizedCartesianProduct(cs, partMatch[cartesianParent], path, candidates[mappingSize], candCount[mappingSize]);
        }
    }
    int depth = mappingSize;
    while (depth >= mappingSize) {
        while (poses[depth] < candCount[depth]) {
            VertexID v = candidates[depth][poses[depth]];
            ++poses[depth];
            if (visited[v]) continue;
            visited[v] = true;
            partMatch[order[depth]] = v;
            if (depth + compressionSize + 1 == tau.numAttributes) {
                // compute compressions
                VertexID **cData = new VertexID *[compressionSize];
                ui *length = new VertexID[compressionSize];
                for (int i = 0; i < compressionSize; ++i) {
                    u = order[depth + i + 1];
                    cData[i] = new VertexID[cs.candidateSet[u].size()];
                    const std::vector<VertexID> &parents = vertexParents[depth + i + 1];
                    for (int j = 0; j < parents.size(); ++j) {
                        VertexID pU = parents[j];
                        neighbors[j] = cs.candidateEdge[pU][partMatch[pU]][u].data();
                        neighborCount[j] = cs.candidateEdge[pU][partMatch[pU]][u].size();
                    }
                    ComputeSetIntersection::LeapfrogJoin(neighbors, neighborCount, parents.size(), cData[i], length[i]);
                }
                delete[] cData;
                delete[] length;
                visited[partMatch[order[depth]]] = false;
            }
            else {
                ++depth;
                poses[depth] = 0;
                u = order[depth];
                const std::vector<VertexID> &parents = vertexParents[depth];
                if (!parents.empty()) {
                    for (int i = 0; i < parents.size(); ++i) {
                        VertexID pU = parents[i];
                        neighbors[i] = cs.candidateEdge[pU][partMatch[pU]][u].data();
                        neighborCount[i] = cs.candidateEdge[pU][partMatch[pU]][u].size();
                    }
                    ComputeSetIntersection::LeapfrogJoin(neighbors, neighborCount, parents.size(), candidates[depth], candCount[depth]);
                }
                else {
                    VertexID cartesianParent = tau.cartesianParent[depth];
                    const std::vector<VertexID> &path = cs.reconstructPath(cartesianParent, u);
                    optimizedCartesianProduct(cs, partMatch[cartesianParent], path, candidates[depth], candCount[depth]);
                }
            }
        }
        --depth;
        if (depth >= mappingSize) visited[partMatch[order[depth]]] = false;
    }

    delete[] neighbors;
    delete[] neighborCount;
}

void nodeJoin(const HyperTree &t, VertexID nID, CandidateSpace &cs, TrieNode *root, bool *visited, VertexID *partMatch,
              int mappingSize, VertexID **candidates, ui *candCount, std::vector<ui> &poses,
              std::vector<std::vector<VertexID>> &tuples) {
    const HyperNode &tau = t.nodes[nID];
    const VertexID *order = tau.attributes;
    if (mappingSize == tau.numAttributes) {
        std::vector<VertexID> tuple(tau.numAttributes - tau.prefixSize);
        for (int i = 0; i < tau.numAttributes - tau.prefixSize; ++i) {
            tuple[i] = partMatch[tau.attributes[i + tau.prefixSize]];
        }
        tuples.push_back(tuple);
        return;
    }
    const std::vector<std::vector<VertexID>> &vertexParents = tau.attributesBefore;
    // handle the first level
    VertexID u = order[mappingSize];
    VertexID **neighbors = new VertexID *[tau.numAttributes];
    ui *neighborCount = new ui[tau.numAttributes];
    if (mappingSize == 0) {
        memcpy(candidates[0], cs.candidateSet[u].data(), cs.candidateSet[u].size() * sizeof(VertexID));
        candCount[0] = cs.candidateSet[u].size();
    }
    else {
        const std::vector<VertexID> &parents = vertexParents[mappingSize];
        if (!parents.empty()) {
            for (int i = 0; i < parents.size(); ++i) {
                VertexID pU = parents[i];
                neighbors[i] = cs.candidateEdge[pU][partMatch[pU]][u].data();
                neighborCount[i] = cs.candidateEdge[pU][partMatch[pU]][u].size();
            }
            ComputeSetIntersection::LeapfrogJoin(neighbors, neighborCount, parents.size(), candidates[mappingSize], candCount[mappingSize]);
        }
        else {
            VertexID cartesianParent = tau.cartesianParent[mappingSize];
            const std::vector<VertexID> &path = cs.reconstructPath(cartesianParent, u);
            optimizedCartesianProduct(cs, partMatch[cartesianParent], path, candidates[mappingSize], candCount[mappingSize]);
        }
    }
    int depth = mappingSize;
    while (depth >= mappingSize) {
        while (poses[depth] < candCount[depth]) {
            VertexID v = candidates[depth][poses[depth]];
            ++poses[depth];
            if (visited[v]) continue;
            visited[v] = true;
            partMatch[order[depth]] = v;
            if (depth + 1 == tau.numAttributes) {
                std::vector<VertexID> tuple(tau.numAttributes - tau.prefixSize);
                for (int i = 0; i < tau.numAttributes - tau.prefixSize; ++i) {
                    tuple[i] = partMatch[tau.attributes[i + tau.prefixSize]];
                }
                tuples.push_back(tuple);
                visited[partMatch[tau.attributes[depth]]] = false;
            }
            else {
                ++depth;
                poses[depth] = 0;
                u = order[depth];
                const std::vector<VertexID> &parents = vertexParents[depth];
                if (!parents.empty()) {
                    for (int i = 0; i < parents.size(); ++i) {
                        VertexID pU = parents[i];
                        neighbors[i] = cs.candidateEdge[pU][partMatch[pU]][u].data();
                        neighborCount[i] = cs.candidateEdge[pU][partMatch[pU]][u].size();
                    }
                    ComputeSetIntersection::LeapfrogJoin(neighbors, neighborCount, parents.size(), candidates[depth], candCount[depth]);
                }
                else {
                    VertexID cartesianParent = tau.cartesianParent[depth];
                    const std::vector<VertexID> &path = cs.reconstructPath(cartesianParent, u);
                    optimizedCartesianProduct(cs, partMatch[cartesianParent], path, candidates[depth], candCount[depth]);
                }
            }
        }
        --depth;
        if (depth >= mappingSize) visited[partMatch[order[depth]]] = false;
    }
    delete[] neighbors;
    delete[] neighborCount;
}

// nodes are roots of tries, and the children are the 0th attribute of each node
void treeJoinScope(std::vector<std::vector<VertexID>> &result, CandidateSpace &cs, const HyperTree &t,
                   std::vector<TrieNode *> &nodes, bool *visited, bool skip) {
    VertexID *partMatch = new VertexID[t.numAttributes];
    VertexID ***nodeCandidates = new VertexID **[nodes.size()];
    ui **nodeCandCount = new ui *[nodes.size()];
    ui ***iters = new ui **[t.extendLevel];
    ui *iterSizes = new ui [t.extendLevel];
    std::vector<std::vector<ui>> nodePoses(nodes.size());
    for (VertexID nID = 0; nID < nodes.size(); ++nID) {
        const HyperNode &tau = t.nodes[nID];
        nodeCandidates[nID]= new VertexID *[tau.numAttributes];
        nodeCandCount[nID] = new ui[tau.numAttributes];
        nodePoses[nID] = std::vector<ui>(tau.numAttributes, 0);
        for (int i = 0; i < tau.numAttributes; ++i) {
            VertexID u = tau.attributes[i];
            nodeCandidates[nID][i] = new VertexID [cs.candidateSet[u].size()];
            nodeCandCount[nID][i] = 0;
        }
    }
    for (int i = 0; i < t.extendLevel; ++i) {
        ui depth = t.defaultPartition.size() + i;
        VertexID u = t.globalOrder[depth];
        iters[i] = new ui *[cs.candidateSet[u].size()];
        for (int j = 0; j < cs.candidateSet[u].size(); ++j) {
            iters[i][j] = new ui[t.nIDs[depth].size()];
        }
        iterSizes[i] = 0;
    }
    if (t.defaultPartition.empty()) {
        for (VertexID nID = 0; nID < t.numNodes; ++nID) {
            nodeJoinScope(t, nID, cs, nodes[nID], visited, partMatch, nodeCandidates[nID], nodeCandCount[nID],
                          nodePoses[nID], skip);
        }
        if (!skip) {
            leapFrogTrieJoin(result, t, cs, t.globalOrder, partMatch, 0, nodes,
                             t.compressionSizes, visited, t.extendLevel, t.nIDs, iters, iterSizes);
        }
        delete[] partMatch;
        for (int i = 0; i < nodes.size(); ++i) {
            delete[] nodeCandCount[i];
            for (int j = 0; j < t.nodes[i].numAttributes; ++j) {
                delete[] nodeCandidates[i][j];
            }
            delete[] nodeCandidates[i];
        }
        delete[] nodeCandidates;
        delete[] nodeCandCount;
        for (int i = 0; i < t.extendLevel; ++i) {
            int depth = t.defaultPartition.size() + i;
            VertexID u = t.globalOrder[depth];
            for (int j = 0; j < cs.candidateSet[u].size(); ++j) {
                delete iters[i][j];
            }
            delete iters[i];
        }
        delete[] iters;
        delete[] iterSizes;
        return;
    }
    VertexID **pCandidates = new VertexID *[t.defaultPartition.size()];
    ui *pCandCount = new ui[t.defaultPartition.size()];
    const std::vector<std::vector<VertexID>> &vertexParents = t.attributesBefore;
    for (int i = 0; i < t.defaultPartition.size(); ++i) {
        VertexID u = t.defaultPartition[i];
        pCandidates[i] = new VertexID[cs.candidateSet[u].size()];
        pCandCount[i] = 0;
    }
    std::vector<ui> pPoses(t.defaultPartition.size(), 0);

    VertexID u = t.defaultPartition[0];
    pCandidates[0] = cs.candidateSet[u].data();
    pCandCount[0] = cs.candidateSet[u].size();
    VertexID **neighbors = new VertexID *[t.defaultPartition.size()];
    ui *neighborCount = new ui[t.defaultPartition.size()];
    std::vector<std::vector<TrieNode *>> traversedNodes(nodes.size());
    std::vector<TrieNode *> lastNodes(nodes.size());
    for (VertexID nID = 0; nID < nodes.size(); ++nID) {
        traversedNodes[nID].push_back(nodes[nID]);
    }
    int depth = 0;
    while (depth >= 0) {
        while (pPoses[depth] < pCandCount[depth]) {
            VertexID v = pCandidates[depth][pPoses[depth]];
            ++pPoses[depth];
            if (visited[v]) continue;
            visited[v] = true;
            partMatch[t.defaultPartition[depth]] = v;
            for (VertexID nID: t.nodesAtStep[depth]) {
                for (int i = 0; i < nodePoses[nID].size(); ++i)
                    nodePoses[nID][i] = 0;
                nodeJoinScope(t, nID, cs, nodes[nID], visited, partMatch, nodeCandidates[nID], nodeCandCount[nID],
                              nodePoses[nID], false);
                TrieNode *pointer = traversedNodes[nID].back();
                while (traversedNodes[nID].size() < t.nodes[nID].prefixSize + 1) {
                    pointer = pointer->nodeChild[0];
                    traversedNodes[nID].push_back(pointer);
                }
            }
            if (depth < t.defaultPartition.size() - 1) {
                // match the next partition attribute
                ++depth;
                pPoses[depth] = 0;
                u = t.defaultPartition[depth];
                const std::vector<VertexID> &parents = vertexParents[depth];
                if (!parents.empty()) {
                    for (int i = 0; i < parents.size(); ++i) {
                        VertexID pU = parents[i];
                        neighbors[i] = cs.candidateEdge[pU][partMatch[pU]][u].data();
                        neighborCount[i] = cs.candidateEdge[pU][partMatch[pU]][u].size();
                    }
                    ComputeSetIntersection::LeapfrogJoin(neighbors, neighborCount, parents.size(), pCandidates[depth], pCandCount[depth]);
                }
                else {
                    VertexID cartesianParent = t.cartesianParent[depth];
                    const std::vector<VertexID> &path = cs.reconstructPath(cartesianParent, u);
                    optimizedCartesianProduct(cs, partMatch[cartesianParent], path, pCandidates[depth], pCandCount[depth]);
                }
            }
            else {
                // the tree join part
                for (VertexID nID = 0; nID < nodes.size(); ++nID) {
                    lastNodes[nID] = traversedNodes[nID].back();
                }
                leapFrogTrieJoin(result, t, cs, t.globalOrder, partMatch, t.defaultPartition.size(),
                                 lastNodes, t.compressionSizes, visited, t.extendLevel, t.nIDs, iters, iterSizes);
                visited[partMatch[t.defaultPartition[depth]]] = false;
                // delete the trie rooted at this partial match
                for (VertexID nID: t.nIDs[depth]) {
                    traversedNodes[nID].pop_back();
                    traversedNodes[nID].back()->nodeChild.pop_front();
                }
            }
        }
        --depth;
        if (depth >= 0) {
            visited[partMatch[t.defaultPartition[depth]]] = false;
            // delete the trie rooted at this partial match
            for (VertexID nID: t.nIDs[depth]) {
                traversedNodes[nID].pop_back();
                traversedNodes[nID].back()->nodeChild.pop_front();
            }
        }
    }

    delete[] partMatch;
    delete[] neighbors;
    for (int i = 0; i < nodes.size(); ++i) {
        delete[] nodeCandCount[i];
        for (int j = 0; j < t.nodes[i].numAttributes; ++j) {
            delete[] nodeCandidates[i][j];
        }
        delete[] nodeCandidates[i];
    }
    delete[] nodeCandidates;
    delete[] nodeCandCount;
    for (int i = 0; i < t.extendLevel; ++i) {
        depth = t.defaultPartition.size() + i;
        u = t.globalOrder[depth];
        for (int j = 0; j < cs.candidateSet[u].size(); ++j) {
            delete[] iters[i][j];
        }
        delete[] iters[i];
    }
    delete[] iters;
    delete[] iterSizes;
}


// for each tree node, set the global join attribute order as: 1. parent shared attributes
// 2. children shared attributes 3. remaining attributes
void emptyHeadedJoin(const HyperTree &t, CandidateSpace &cs, std::vector<TrieNode *> &nodes, bool *visited,
                     std::vector<std::vector<VertexID>> &result, std::vector<double> &times) {
    VertexID *partMatch = new VertexID[cs.candidateSet.size()];
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        auto start = std::chrono::steady_clock::now();
        std::vector<std::vector<VertexID>> tuples;
        const HyperNode &tau = t.nodes[nID];
        // local join
        const std::vector<std::vector<VertexID>> &vertexParents = tau.attributesBefore;
        VertexID **candidates = new VertexID *[tau.numAttributes];
        ui *candCount = new ui[tau.numAttributes];
        ui *poses = new ui[tau.numAttributes];
        for (int i = 0; i < tau.numAttributes; ++i) {
            candidates[i] = new VertexID[cs.candidateSet[tau.attributes[i]].size()];
            candCount[i] = 0;
            poses[i] = 0;
        }
        VertexID **neighbors = new VertexID *[tau.numAttributes];
        ui *neighborCount = new ui[tau.numAttributes];
        VertexID u = tau.attributes[0];
        memcpy(candidates[0], cs.candidateSet[u].data(), cs.candidateSet[u].size() * sizeof(VertexID));
        candCount[0] = cs.candidateSet[u].size();
        int depth = 0;
        while (depth >= 0) {
            while (poses[depth] < candCount[depth]) {
                VertexID v = candidates[depth][poses[depth]];
                ++poses[depth];
                if (visited[v]) continue;
                visited[v] = true;
                partMatch[tau.attributes[depth]] = v;
                if (depth == tau.numAttributes - 1) {
                    std::vector<VertexID> tuple(tau.numAttributes);
                    for (int i = 0; i < tau.numAttributes; ++i) {
                        tuple[i] = partMatch[tau.attributes[i]];
                    }
                    tuples.push_back(tuple);
                    visited[partMatch[tau.attributes[depth]]] = false;
                }
                else {
                    ++depth;
                    poses[depth] = 0;
                    u = tau.attributes[depth];
                    const std::vector<VertexID> &parents = vertexParents[depth];
                    if (!parents.empty()) {
                        for (int i = 0; i < parents.size(); ++i) {
                            VertexID pU = parents[i];
                            neighbors[i] = cs.candidateEdge[pU][partMatch[pU]][u].data();
                            neighborCount[i] = cs.candidateEdge[pU][partMatch[pU]][u].size();
                        }
                        ComputeSetIntersection::LeapfrogJoin(neighbors, neighborCount, parents.size(), candidates[depth], candCount[depth]);
                    }
                    else {
                        VertexID cartesianParent = tau.cartesianParent[depth];
                        const std::vector<VertexID> &path = cs.reconstructPath(cartesianParent, u);
                        optimizedCartesianProduct(cs, partMatch[cartesianParent], path, candidates[depth], candCount[depth]);
                    }
                }
            }
            --depth;
            if (depth >= 0) visited[partMatch[tau.attributes[depth]]] = false;
        }
        for (int i = 0; i < tau.numAttributes; ++i) {
            delete[] candidates[i];
        }
        delete[] candidates;
        delete[] candCount;
        delete[] poses;
        delete[] neighbors;
        delete[] neighborCount;
        buildTrie(tuples, nodes[nID], t.trieOrder[nID]);
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsedSeconds = end - start;
        times.push_back(elapsedSeconds.count());
    }
    ui ***iters = new ui **[t.extendLevel];
    ui *iterSizes = new ui [t.extendLevel];
    for (int i = 0; i < t.extendLevel; ++i) {
        VertexID u = t.globalOrder[i];
        iters[i] = new ui *[cs.candidateSet[u].size()];
        for (int j = 0; j < cs.candidateSet[u].size(); ++j) {
            iters[i][j] = new ui[t.nIDs[i].size()];
        }
        iterSizes[i] = 0;
    }
    leapFrogTrieJoin(result, t, cs, t.globalOrder, partMatch, t.defaultPartition.size(),
                     nodes, t.compressionSizes, visited, t.extendLevel, t.nIDs, iters, iterSizes);
    for (int i = 0; i < t.extendLevel; ++i) {
        VertexID u = t.globalOrder[i];
        for (int j = 0; j < cs.candidateSet[u].size(); ++j) {
            delete[] iters[i][j];
        }
        delete[] iters[i];
    }
    delete[] iters;
    delete[] iterSizes;
    delete[] partMatch;
}

void sharedJoin(const HyperTree &t, const PrefixNode *pt, const Graph &query, CandidateSpace &cs,
                std::vector<TrieNode *> &trieNodes, bool *visited, std::vector<std::vector<VertexID>> &result,
                size_t &count, bool traverse) {
    std::vector<std::vector<std::vector<VertexID>>> tuples(t.numNodes);
    VertexID *partMatch = new VertexID[t.numAttributes];
    VertexID ***nodeCandidates = new VertexID **[trieNodes.size()];
    ui **nodeCandCount = new ui *[trieNodes.size()];
    const HyperNode &globalNode = t.nodes[t.numNodes - 1];
    const PrefixNode *pn = pt;
    ui maxSize = cs.getMaxSize();
    ui ***iters = new ui **[t.numAttributes];
    ui *iterSizes = new ui [t.numAttributes];
    for (int i = 0; i < t.numAttributes; ++i) {
        iters[i] = new ui *[maxSize];
        for (int j = 0; j < maxSize; ++j) {
            iters[i][j] = new ui[t.numAttributes];
        }
        iterSizes[i] = 0;
    }
    std::vector<std::vector<ui>> nodePoses(trieNodes.size());
    for (VertexID nID = 0; nID < trieNodes.size(); ++nID) {
        const HyperNode &tau = t.nodes[nID];
        nodeCandidates[nID]= new VertexID *[tau.numAttributes];
        nodeCandCount[nID] = new ui[tau.numAttributes];
        nodePoses[nID] = std::vector<ui>(tau.numAttributes, 0);
        for (int i = 0; i < tau.numAttributes; ++i) {
            VertexID u = tau.attributes[i];
            nodeCandidates[nID][i] = new VertexID [cs.candidateSet[u].size()];
            nodeCandCount[nID][i] = 0;
        }
    }
    ui height = pt->getHeight();
    VertexID **pCandidates = new VertexID *[height];
    ui *pCandCount = new ui[height];
    for (int i = 0; i < height; ++i) {
        pCandidates[i] = new VertexID[maxSize];
        pCandCount[i] = 0;
    }
    std::vector<ui> pPoses(height, 0);
    std::vector<ui> childPoses(height, 0);
    VertexID **neighbors = new VertexID *[height];
    ui *neighborCount = new ui[height];
    DynamicArray<TrieNode *> ** children = new DynamicArray<TrieNode *> *[t.numAttributes];
    // join relational columns together with edge columns
    ui maxEdgeSize = maxNumBackWard(t.globalOrder, query);
    ui maxNum = height + globalNode.numAttributes - globalNode.prefixSize;
    std::vector<std::vector<DynamicArray<TrieNode *> *>> edgeColumns(maxNum);
    for (int i = 0 ; i < maxNum; ++i) {
        edgeColumns[i] = std::vector<DynamicArray<TrieNode *> *>(maxEdgeSize);
        for (int j = 0; j < maxEdgeSize; ++j) {
            edgeColumns[i][j] = new DynamicArray<TrieNode *>();
            for (int k = 0; k < maxSize; ++k) {
                TrieNode * pointer = new TrieNode();
                edgeColumns[i][j]->push_back(pointer);
            }
        }
    }
    std::vector<std::vector<TrieNode *>> traversedNodes(trieNodes.size());
    std::vector<TrieNode *> lastNodes(trieNodes.size());
    for (VertexID nID = 0; nID < trieNodes.size(); ++nID) {
        traversedNodes[nID].push_back(trieNodes[nID]);
        lastNodes[nID] = trieNodes[nID];
    }
    pn = pt;
    std::vector<const PrefixNode *> nodes(height, nullptr);
    std::vector<int> mappingSizes = getMappingSizes(t, pt);
    int depth = 0;
    for (VertexID nID: pn -> nIDsToCall) {
        if (nID != t.numNodes - 1) {
            for (int i = 0; i < nodePoses[nID].size(); ++i)
                nodePoses[nID][i] = 0;
            nodeJoin(t, nID, cs, trieNodes[nID], visited, partMatch, 0, nodeCandidates[nID], nodeCandCount[nID],
                     nodePoses[nID], tuples[nID]);
        }
    }
    // create candidates for the first prefixNode
    if (!pn -> children.empty()) {
        if (pn -> children[0] -> pathToGlobal) {
            for (VertexID nID : pn -> nIDsToBuild) {
                if (tuples[nID].empty()) return;
                buildTrie(tuples[nID], trieNodes[nID], t.trieOrder[nID]);
            }
            if (pn-> children[0] -> nIDsToJoin.empty()) {
                VertexID u = pn -> children[0] -> u;
                edgeColumns[0][0] ->setSize(cs.candidateSet[u].size());
                for (int i = 0; i < cs.candidateSet[u].size(); ++i) {
                    edgeColumns[0][0][0][i]->value = cs.candidateSet[u][i];
                    iters[0][i][0] = i;
                }
                iterSizes[0] = cs.candidateSet[pn -> children[0] -> u].size();
            }
            else {
                for (int i = 0; i < pn -> children[0] -> nIDsToJoin.size(); ++i) {
                    VertexID id = pn -> children[0] -> nIDsToJoin[i];
                    children[i] = &trieNodes[id] -> nodeChild;
                }
                leapFrogJoin(children, pn -> children[0] -> nIDsToJoin.size(), iters[0], iterSizes[0]);
            }
        }
        else {
            pCandidates[0] = cs.candidateSet[pn -> children[0] -> u].data();
            pCandCount[0] = cs.candidateSet[pn -> children[0] -> u].size();
        }
    }
    else {
        for (VertexID nID: pn -> nIDsToCall) {
            if (nID != t.numNodes - 1) {
                if (tuples[nID].empty()) return;
                buildTrie(tuples[nID], trieNodes[nID], t.trieOrder[nID]);
            }
            else {
                globalJoin(result, count, query, t, cs, partMatch, 0, 0, trieNodes, visited,
                           globalNode.nIDs, globalNode.attributesBefore, globalNode.cartesianParent, iters, iterSizes,
                           traverse, edgeColumns);
            }
        }
    }
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            const PrefixNode *current = pn -> children[childPoses[depth]];
            nodes[depth] = current;
            VertexID u = current -> u;
            bool nextLevel = false;
            if (!current -> pathToGlobal) {
                while (pPoses[depth] < pCandCount[depth]) {
                    VertexID v = pCandidates[depth][pPoses[depth]];
                    ++pPoses[depth];
                    if (visited[v]) continue;
                    visited[v] = true;
                    partMatch[u] = v;
                    for (VertexID nID: current -> nIDsToCall) {
                        for (int i = 0; i < nodePoses[nID].size(); ++i)
                            nodePoses[nID][i] = 0;
                        nodeJoin(t, nID, cs, trieNodes[nID], visited, partMatch, mappingSizes[nID], nodeCandidates[nID],
                                 nodeCandCount[nID], nodePoses[nID], tuples[nID]);
                    }
                    if (!current -> children.empty()) {
                        nextLevel = true;
                        break;
                    }
                    else visited[v] = false;
                }
            }
            else {
                while (pPoses[depth] < iterSizes[depth]) {
                    VertexID v;
                    ui childPos = iters[depth][pPoses[depth]][0];
                    if (current -> nIDsToJoin.empty()) {
                        v = edgeColumns[depth][0][0][childPos]->value;
                    }
                    else {
                        VertexID tmp = current -> nIDsToJoin[0];
                        v = traversedNodes[tmp].back()->nodeChild[childPos]->value;
                    }
                    ++pPoses[depth];
                    if (visited[v]) continue;
                    visited[v] = true;
                    partMatch[u] = v;
                    for (VertexID nID: current -> nIDsToCall) {
                        if (nID != t.numNodes - 1) {
                            for (int i = 0; i < nodePoses[nID].size(); ++i)
                                nodePoses[nID][i] = 0;
                            nodeJoin(t, nID, cs, trieNodes[nID], visited, partMatch, mappingSizes[nID], nodeCandidates[nID],
                                     nodeCandCount[nID],
                                     nodePoses[nID], tuples[nID]);
                        }
                    }
                    // all joined relations extend one level
                    for (int i = 0; i < current -> nIDsToJoin.size(); ++i) {
                        VertexID nID = current -> nIDsToJoin[i];
                        childPos = iters[depth][pPoses[depth] - 1][i];
                        traversedNodes[nID].push_back(traversedNodes[nID].back()->nodeChild[childPos]);
                        lastNodes[nID] = traversedNodes[nID].back();
                    }
                    if (!current -> children.empty()) {
                        nextLevel = true;
                        break;
                    }
                    else {
                        for (VertexID nID: current -> nIDsToCall) {
                            if (nID == t.numNodes - 1) {
                                bool flag = true;
                                for (VertexID nID2 : current -> nIDsToBuild) {
                                    if (tuples[nID2].empty()) {
                                        flag = false;
                                        break;
                                    }
                                }
                                if (flag) {
                                    for (VertexID nID2 : current -> nIDsToBuild)
                                        buildTrie(tuples[nID2], trieNodes[nID2], t.trieOrder[nID2]);
                                    globalJoin(result, count, query, t, cs, partMatch, mappingSizes[nID], depth + 1,
                                               lastNodes,
                                               visited, globalNode.nIDs, globalNode.attributesBefore,
                                               globalNode.cartesianParent, iters, iterSizes, traverse,
                                               edgeColumns);
                                }
                                else {
                                    for (VertexID nID2 : current -> nIDsToBuild)
                                        tuples[nID2].clear();
                                }
                            }
                        }
                        visited[v] = false;
                        for (VertexID nID : current -> nIDsToJoin) {
                            traversedNodes[nID].pop_back();
                            lastNodes[nID] = traversedNodes[nID].back();
                        }
                    }
                }
            }
            if (!nextLevel) {
                ++childPoses[depth];
                pPoses[depth] = 0;
                if (childPoses[depth] == pn->children.size()) break;
                current = pn->children[childPoses[depth]];
                u = current->u;
            }
            else {
                ++depth;
                pPoses[depth] = 0;
                childPoses[depth] = 0;
                pn = current;
                current = pn -> children[0];
                u = current->u;
            }
            if (current -> pathToGlobal) {
                bool flag = true;
                for (VertexID nID : pn -> nIDsToBuild) {
                    if (tuples[nID].empty()) {
                        flag = false;
                        break;
                    }
                }
                if (flag) {
                    for (VertexID nID : pn -> nIDsToBuild)
                        buildTrie(tuples[nID], trieNodes[nID], t.trieOrder[nID]);
                }
                else {
                    for (VertexID nID2 : pn -> nIDsToBuild) tuples[nID2].clear();
                    visited[partMatch[pn->u]] = false;
                    for (VertexID nID : pn -> nIDsToJoin) {
                        traversedNodes[nID].pop_back();
                        lastNodes[nID] = traversedNodes[nID].back();
                    }
                    --depth;
                    if (depth == 0) pn = pt;
                    else pn = nodes[depth - 1];
                    continue;
                }
                if (current -> cartesianParent != 99) {
                    VertexID pU = current -> cartesianParent;
                    const std::vector<VertexID> &path = cs.reconstructPath(pU, u);
                    optimizedCartesianProduct(cs, partMatch[pU], path, edgeColumns[depth][0]);
                    iterSizes[depth] = edgeColumns[depth][0]->size();
                    for (int i = 0; i < edgeColumns[depth][0]->size(); ++i) {
                        iters[depth][i][0] = i;
                    }
                }
                else if (current ->nIDsToJoin.empty() && current -> attributesBefore.empty()) {
                    edgeColumns[0][0] ->setSize(cs.candidateSet[u].size());
                    for (int i = 0; i < cs.candidateSet[u].size(); ++i) {
                        edgeColumns[0][0][0][i]->value = cs.candidateSet[u][i];
                        iters[0][i][0] = i;
                    }
                    iterSizes[0] = cs.candidateSet[u].size();
                }
                for (int i = 0; i < current->nIDsToJoin.size(); ++i) {
                    VertexID nID = current->nIDsToJoin[i];
                    children[i] = &lastNodes[nID]->nodeChild;
                }
                for (int i = 0; i < current->attributesBefore.size(); ++i) {
                    VertexID pU = current->attributesBefore[i];
                    VertexID pV = partMatch[pU];
                    edgeColumns[depth][i]->setSize(cs.candidateEdge[pU][pV][u].size());
                    for (int j = 0; j < cs.candidateEdge[pU][pV][u].size(); ++j)
                        edgeColumns[depth][i][0][j]->value = cs.candidateEdge[pU][pV][u][j];
                    children[current->nIDsToJoin.size() + i] = edgeColumns[depth][i];
                }
                if (current->nIDsToJoin.size() + current->attributesBefore.size() > 0)
                    leapFrogJoin(children, current->nIDsToJoin.size() + current->attributesBefore.size(), iters[depth], iterSizes[depth]);
            }
            else {
                if (depth != 0) {
                    const std::vector<VertexID> &parents = current->attributesBefore;
                    if (!parents.empty()) {
                        for (int i = 0; i < parents.size(); ++i) {
                            VertexID pU = parents[i];
                            neighbors[i] = cs.candidateEdge[pU][partMatch[pU]][u].data();
                            neighborCount[i] = cs.candidateEdge[pU][partMatch[pU]][u].size();
                        }
                        ComputeSetIntersection::LeapfrogJoin(neighbors, neighborCount, parents.size(), pCandidates[depth], pCandCount[depth]);
                    }
                    else {
                        VertexID cartesianParent = current->cartesianParent;
                        const std::vector<VertexID> &path = cs.reconstructPath(cartesianParent, u);
                        optimizedCartesianProduct(cs, partMatch[cartesianParent], path, pCandidates[depth], pCandCount[depth]);
                    }
                }
                else {
                    pCandidates[0] = cs.candidateSet[u].data();
                    pCandCount[0] = cs.candidateSet[u].size();
                }
            }
        }
        --depth;
        if (depth >= 0) {
            const PrefixNode *current = pn;
            if (depth == 0) pn = pt;
            else pn = nodes[depth - 1];
            VertexID u = current -> u;
            VertexID v = partMatch[u];
            for (VertexID nID: current -> nIDsToCall) {
                if (nID == t.numNodes - 1) {
                    bool flag = true;
                    for (VertexID nID2 : current -> nIDsToBuild) {
                        if (tuples[nID2].empty()) {
                            flag = false;
                            break;
                        }
                    }
                    if (flag) {
                        for (VertexID nID2 : current -> nIDsToBuild)
                            buildTrie(tuples[nID2], trieNodes[nID2], t.trieOrder[nID2]);
                        globalJoin(result, count, query, t, cs, partMatch, mappingSizes[nID], depth + 1, lastNodes,
                                   visited, globalNode.nIDs, globalNode.attributesBefore,
                                   globalNode.cartesianParent, iters, iterSizes, traverse, edgeColumns);
                    }
                    else {
                        for (VertexID nID2 : current -> nIDsToBuild)
                            tuples[nID2].clear();
                    }
                }
            }
            visited[v] = false;
            if (current -> pathToGlobal) {
                for (VertexID nID : current -> nIDsToJoin) {
                    traversedNodes[nID].pop_back();
                    lastNodes[nID] = traversedNodes[nID].back();
                }
            }
        }
    }
    for (VertexID nID: pn -> nIDsToCall) {
        if (nID == t.numNodes - 1) {
            bool flag = true;
            for (VertexID nID2 : pn -> nIDsToBuild) {
                if (tuples[nID2].empty()) {
                    flag = false;
                    break;
                }
            }
            if (flag) {
                for (VertexID nID2 : pn -> nIDsToBuild)
                    buildTrie(tuples[nID2], trieNodes[nID2], t.trieOrder[nID2]);
                globalJoin(result, count, query, t, cs, partMatch, mappingSizes[nID], depth + 1, lastNodes, visited,
                           globalNode.nIDs, globalNode.attributesBefore,
                           globalNode.cartesianParent, iters, iterSizes, traverse, edgeColumns);
            }
            else {
                for (VertexID nID2 : pn -> nIDsToBuild)
                    tuples[nID2].clear();
            }
        }
    }
    delete[] partMatch;
    delete[] neighbors;
    for (int i = 0; i < trieNodes.size(); ++i) {
        delete[] nodeCandCount[i];
        for (int j = 0; j < t.nodes[i].numAttributes; ++j) {
            delete[] nodeCandidates[i][j];
        }
        delete[] nodeCandidates[i];
    }
    delete[] nodeCandidates;
    delete[] nodeCandCount;
    for (int i = 0; i < t.numAttributes; ++i) {
        for (int j = 0; j < maxSize; ++j) {
            delete[] iters[i][j];
        }
        delete[] iters[i];
    }
    delete[] iters;
    delete[] iterSizes;
    delete[] children;
    for (int i = 0; i < maxNum; ++i) {
        for (int j = 0; j < maxEdgeSize; ++j) {
            DynamicArray<TrieNode *> *dynamicArrayPtr = edgeColumns[i][j];
            for (ui k = 0; k < dynamicArrayPtr->size(); ++k) {
                delete (*dynamicArrayPtr)[k];
            }
            delete dynamicArrayPtr;
        }
    }
}

/*
 * for a given partition specified by the partial match, join all shared attributes
 * nIDsToCall: at step i, the nodes to join are nIDsToCall[i]
 * extendLevel: the total number of levels to extend
 * */

void
leapFrogTrieJoin(std::vector<std::vector<VertexID>> &result, const HyperTree &t, CandidateSpace &cs,
                 const std::vector<VertexID> &order, VertexID *partMatch, int mappingSize,
                 const std::vector<TrieNode *> &nodes, const std::vector<ui> &compressionSize, bool *visited,
                 int extendLevel, const std::vector<std::vector<VertexID>> &nIDs, ui ***iters, ui *iterSize) {
    if (extendLevel == 0) {
        traverse(result, t, order, partMatch, mappingSize, nodes, visited);
        return;
    }
    std::vector<ui> poses(extendLevel, 0);
    std::vector<TrieNode *> lastNodes(nodes.size());
    std::vector<std::vector<TrieNode *>> traversedNodes(nodes.size());
    for (VertexID nID = 0; nID < nodes.size(); ++nID) {
        traversedNodes[nID].push_back(nodes[nID]);
        lastNodes[nID] = nodes[nID];
    }
    // for each level and each nID, the trie node is traversedNodes[]
    DynamicArray<TrieNode *> ** children = new DynamicArray<TrieNode *> *[nodes.size()];
    // initialize the first level
    for (int i = 0; i < nIDs[mappingSize].size(); ++i) {
        VertexID id = nIDs[mappingSize][i];
        children[i] = &nodes[id]->nodeChild;
    }
    leapFrogJoin(children, nIDs[mappingSize].size(), iters[0], iterSize[0]);
    int depth = 0;
    while (depth >= 0) {
        while (poses[depth] < iterSize[depth]) {
            VertexID nID = nIDs[mappingSize + depth][0];
            ui childPos = iters[depth][poses[depth]][0];
            VertexID v = traversedNodes[nID].back()->nodeChild[childPos]->value;
            ++poses[depth];
            if (visited[v]) continue;
            visited[v] = true;
            partMatch[order[mappingSize + depth]] = v;
            // all joined relations extend one level
            for (int i = 0; i < nIDs[mappingSize + depth].size(); ++i) {
                nID = nIDs[mappingSize + depth][i];
                childPos = iters[depth][poses[depth] - 1][i];
                traversedNodes[nID].push_back(traversedNodes[nID].back()->nodeChild[childPos]);
                lastNodes[nID] = traversedNodes[nID].back();
            }
            if (depth == extendLevel - 1) {
                // when a match of shared attributes are found, traverse the remaining attributes and compute the cartesian products
                traverse(result, t, order, partMatch, mappingSize + depth + 1, lastNodes, visited);
                for (int i = 0; i < nIDs[mappingSize + depth].size(); ++i) {
                    nID = nIDs[mappingSize + depth][i];
                    traversedNodes[nID].pop_back();
                    lastNodes[nID] = traversedNodes[nID].back();
                }
                visited[v] = false;
            }
            else {
                ++depth;
                poses[depth] = 0;
                for (int i = 0; i < nIDs[mappingSize + depth].size(); ++i) {
                    nID = nIDs[mappingSize + depth][i];
                    children[i] = &lastNodes[nID]->nodeChild;
                }
                leapFrogJoin(children, nIDs[mappingSize + depth].size(), iters[depth], iterSize[depth]);
            }
        }
        --depth;
        if (depth >= 0) {
            VertexID v = partMatch[t.globalOrder[mappingSize + depth]];
            for (int i = 0; i < nIDs[mappingSize + depth].size(); ++i) {
                VertexID nID = nIDs[mappingSize + depth][i];
                traversedNodes[nID].pop_back();
                lastNodes[nID] = traversedNodes[nID].back();
            }
            visited[v] = false;
        }
    }

    delete[] children;
}

void globalJoin(std::vector<std::vector<VertexID>> &result, size_t &count, const Graph &query, const HyperTree &t,
                CandidateSpace &cs, VertexID *partMatch, int mappingSize, int pathLength,
                std::vector<TrieNode *> &nodes, bool *visited, const std::vector<std::vector<VertexID>> &nIDs,
                const std::vector<std::vector<VertexID>> &vertexParents, const std::vector<VertexID> &cartesianParent,
                ui ***iters, ui *iterSize, bool traversal,
                std::vector<std::vector<DynamicArray<TrieNode *> *>> &edgeColumns) {
    const HyperNode &globalNode = t.nodes[t.numNodes - 1];
    if (globalNode.numAttributes == mappingSize) {
#ifdef COLLECT_RESULT
        if (traversal) traverse(result, t, t.globalOrder, partMatch, mappingSize, nodes, visited);
#else
        traverse(count, query, t, t.globalOrder, partMatch, mappingSize, nodes, visited, t.extendLevel);
#endif
        return;
    }
    std::vector<VertexID> order(globalNode.attributes, globalNode.attributes + globalNode.numAttributes);
    std::vector<ui> poses(globalNode.numAttributes, 0);
    std::vector<TrieNode *> lastNodes(nodes.size());
    std::vector<std::vector<TrieNode *>> traversedNodes(nodes.size());
    for (VertexID nID = 0; nID < nodes.size(); ++nID) {
        traversedNodes[nID].push_back(nodes[nID]);
        lastNodes[nID] = nodes[nID];
    }
    // for each level and each nID, the trie node is traversedNodes[]
    DynamicArray<TrieNode *> ** children = new DynamicArray<TrieNode *> *[t.numAttributes];
    // initialize the first level
    if (mappingSize == 0 && nIDs[0].empty()) {
        VertexID u = order[0];
        edgeColumns[0][0] ->setSize(cs.candidateSet[u].size());
        for (int i = 0; i < cs.candidateSet[u].size(); ++i) {
            edgeColumns[0][0][0][i] -> value = cs.candidateSet[u][i];
            iters[0][i][0] = i;
        }
        iterSize[0] = cs.candidateSet[u].size();
    }
    if (cartesianParent[mappingSize] != 99) {
        VertexID u = order[mappingSize];
        VertexID pU = cartesianParent[mappingSize];
        const std::vector<VertexID> &path = cs.reconstructPath(pU, u);
        optimizedCartesianProduct(cs, partMatch[pU], path, edgeColumns[pathLength][0]);
        iterSize[pathLength] = edgeColumns[pathLength][0]->size();
        for (int i = 0; i < edgeColumns[pathLength][0]->size(); ++i) {
            iters[pathLength][i][0] = i;
        }
    }
    for (int i = 0; i < nIDs[mappingSize].size(); ++i) {
        VertexID id = nIDs[mappingSize][i];
        children[i] = &nodes[id]->nodeChild;
    }
    for (int i = 0; i < vertexParents[mappingSize].size(); ++i) {
        VertexID pU = vertexParents[mappingSize][i];
        VertexID pV = partMatch[pU];
        VertexID u = order[mappingSize];
        edgeColumns[pathLength][i]->setSize(cs.candidateEdge[pU][pV][u].size());
        for (int j = 0; j < cs.candidateEdge[pU][pV][u].size(); ++j)
            edgeColumns[pathLength][i][0][j]->value = cs.candidateEdge[pU][pV][u][j];
        children[nIDs[mappingSize].size() + i] = edgeColumns[pathLength][i];
    }
    if (nIDs[mappingSize].size() + vertexParents[mappingSize].size() != 0)
        leapFrogJoin(children, nIDs[mappingSize].size() + vertexParents[mappingSize].size(), iters[pathLength], iterSize[pathLength]);
    std::vector<std::vector<VertexID>> tuples;
    int depth = 0;
    while (depth >= 0) {
        while (poses[mappingSize + depth] < iterSize[pathLength + depth]) {
            VertexID v;
            ui childPos = iters[pathLength + depth][poses[mappingSize + depth]][0];
            if (nIDs[mappingSize + depth].empty()) v = edgeColumns[pathLength + depth][0][0][childPos]->value;
            else {
                VertexID nID = nIDs[mappingSize + depth][0];
                v = traversedNodes[nID].back()->nodeChild[childPos]->value;
            }
            ++poses[mappingSize + depth];
            if (visited[v]) continue;
            visited[v] = true;
            partMatch[order[mappingSize + depth]] = v;
            // all joined relations extend one level
            for (int i = 0; i < nIDs[mappingSize + depth].size(); ++i) {
                VertexID nID = nIDs[mappingSize + depth][i];
                childPos = iters[pathLength + depth][poses[mappingSize + depth] - 1][i];
                traversedNodes[nID].push_back(traversedNodes[nID].back()->nodeChild[childPos]);
                lastNodes[nID] = traversedNodes[nID].back();
            }
            if (depth + 1 + mappingSize == globalNode.numAttributes) {
#ifdef COLLECT_RESULT
                if (traversal) traverse(result, t, t.globalOrder, partMatch, mappingSize + depth + 1, lastNodes, visited);
#else
                if (t.extendLevel != globalNode.numAttributes) {
                    std::vector<VertexID> tuple(globalNode.numAttributes - t.extendLevel);
                    for (int i = 0; i < globalNode.numAttributes - t.extendLevel; ++i) {
                        tuple[i] = partMatch[globalNode.attributes[i + t.extendLevel]];
                    }
                    tuples.push_back(tuple);
                }
                else {
                    // when a match of shared attributes is found, traverse the remaining attributes
                    traverse(count, query, t, t.globalOrder, partMatch, mappingSize + depth + 1, lastNodes, visited, t.extendLevel);
                }
#endif
                for (int i = 0; i < nIDs[mappingSize + depth].size(); ++i) {
                    VertexID nID = nIDs[mappingSize + depth][i];
                    traversedNodes[nID].pop_back();
                    lastNodes[nID] = traversedNodes[nID].back();
                }
                visited[v] = false;
            }
            else {
                ++depth;
                poses[mappingSize + depth] = 0;
                if (cartesianParent[mappingSize + depth] != 99) {
                    VertexID u = order[mappingSize + depth];
                    VertexID pU = cartesianParent[mappingSize + depth];
                    const std::vector<VertexID> &path = cs.reconstructPath(pU, u);
                    optimizedCartesianProduct(cs, partMatch[pU], path, edgeColumns[pathLength + depth][0]);
                }
                for (int i = 0; i < nIDs[mappingSize + depth].size(); ++i) {
                    VertexID nID = nIDs[mappingSize + depth][i];
                    children[i] = &lastNodes[nID]->nodeChild;
                }
                for (int i = 0; i < vertexParents[mappingSize + depth].size(); ++i) {
                    VertexID pU = vertexParents[mappingSize + depth][i];
                    VertexID pV = partMatch[pU];
                    VertexID u = order[mappingSize + depth];
                    edgeColumns[pathLength + depth][i]->setSize(cs.candidateEdge[pU][pV][u].size());
                    for (int j = 0; j < cs.candidateEdge[pU][pV][u].size(); ++j)
                        edgeColumns[pathLength + depth][i][0][j]->value = cs.candidateEdge[pU][pV][u][j];
                    children[nIDs[mappingSize + depth].size() + i] = edgeColumns[pathLength + depth][i];
                }

                if (nIDs[mappingSize + depth].size() + vertexParents[mappingSize + depth].size() > 0)
                    leapFrogJoin(children, nIDs[mappingSize + depth].size() + vertexParents[mappingSize + depth].size(),
                             iters[pathLength + depth], iterSize[pathLength + depth]);
            }
        }
        --depth;
#ifndef COLLECT_RESULT
        if (depth + 1 + mappingSize == t.extendLevel && !tuples.empty()) {
            buildTrie(tuples, lastNodes[t.numNodes - 1], t.trieOrder[t.numNodes - 1]);
            traverse(count, query, t, t.globalOrder, partMatch, mappingSize + depth + 1, lastNodes, visited, t.extendLevel);
            lastNodes[t.numNodes - 1]->nodeChild.clear();
        }
#endif
        if (depth >= 0) {
            VertexID v = partMatch[order[mappingSize + depth]];
            for (int i = 0; i < nIDs[mappingSize + depth].size(); ++i) {
                VertexID nID = nIDs[mappingSize + depth][i];
                traversedNodes[nID].pop_back();
                lastNodes[nID] = traversedNodes[nID].back();
            }
            visited[v] = false;
        }
    }

    delete[] children;
}

// call this function when a match of shared attributes is found
void traverse(std::vector<std::vector<VertexID>> &result, const HyperTree &t, const std::vector<VertexID> &order,
              VertexID *partMatch, int mappingSize, const std::vector<TrieNode *> &nodes, bool *visited) {
//#ifdef COLLECT_STATISTICS
//    auto start = std::chrono::steady_clock::now();
//#endif
    if (t.defaultPartition.size() > mappingSize) mappingSize = t.defaultPartition.size();
    ui numNodes = nodes.size();
    ui totalLevel = 0;
    std::vector<ui> levels(numNodes, 0);
    for (int i = 0; i < numNodes; ++i) {
        ui level = 0;
        const TrieNode *tau = nodes[i];
        while (!tau->nodeChild.empty()) {
            tau = tau->nodeChild[0];
            ++level;
        }
        totalLevel += level;
        levels[i] = level;
    }
    if (totalLevel == 0) {
        produceResult(result, t, partMatch, nodes, visited);
        return;
    }
    std::vector<TrieNode *> lastNodes(numNodes);
    std::vector<ui> poses(totalLevel, 0);
    std::vector<ui> counts(totalLevel, 0);
    std::vector<std::vector<TrieNode *>> traversedNodes(nodes.size());
    for (VertexID nID = 0; nID < nodes.size(); ++nID) {
        traversedNodes[nID].push_back(nodes[nID]);
        lastNodes[nID] = nodes[nID];
    }
    // initialize the first level
    VertexID firstLevelNID = t.nIDs[mappingSize][0];
    counts[0] = traversedNodes[firstLevelNID][0]->nodeChild.size();
    int depth = 0;
    while (depth >= 0) {
        while (poses[depth] < counts[depth]) {
            VertexID nID = t.nIDs[mappingSize + depth][0];
            VertexID v = traversedNodes[nID].back()->nodeChild[poses[depth]]->value;
            ui childPos = poses[depth];
            ++poses[depth];
            if (visited[v]) continue;
            visited[v] = true;
            partMatch[order[mappingSize + depth]] = v;
            traversedNodes[nID].push_back(traversedNodes[nID].back()->nodeChild[childPos]);
            lastNodes[nID] = traversedNodes[nID].back();
            if (depth == totalLevel - 1) {
                produceResult(result, t, partMatch, lastNodes, visited);
                traversedNodes[nID].pop_back();
                lastNodes[nID] = traversedNodes[nID].back();
                visited[v] = false;
            }
            else {
                ++depth;
                poses[depth] = 0;
                nID = t.nIDs[mappingSize + depth][0];
                counts[depth] = lastNodes[nID]->nodeChild.size();
            }
        }
        --depth;
        if (depth >= 0) {
            VertexID nID = t.nIDs[mappingSize + depth][0];
            traversedNodes[nID].pop_back();
            lastNodes[nID] = traversedNodes[nID].back();
            visited[partMatch[t.globalOrder[mappingSize + depth]]] = false;
        }
    }
}


void
traverse(size_t &count, const Graph &query, const HyperTree &t, const std::vector<VertexID> &order, VertexID *partMatch,
         int mappingSize, const std::vector<TrieNode *> &nodes, bool *visited, int extendLevel) {
    ++gNumTraverse;
    std::vector<VertexID> depthToNID = t.adaptiveDepthToNID.at(extendLevel);
    std::vector<ui> levels = t.adaptiveLevels.at(extendLevel);
    std::vector<std::vector<VertexID>> groups = t.adaptiveGroups.at(extendLevel);
    if (t.defaultPartition.size() > mappingSize) mappingSize = t.defaultPartition.size();
    ui numNodes = nodes.size();
    ui totalLevel = 0;
    for (VertexID nID = 0; nID < nodes.size(); ++nID) {
        totalLevel += levels[nID];
    }
    std::vector<TrieNode *> lastNodes(numNodes);
    std::vector<ui> poses(totalLevel, 0);
    std::vector<ui> numBranches(totalLevel, 0);
    std::vector<std::vector<TrieNode *>> traversedNodes(nodes.size());
    for (VertexID nID = 0; nID < nodes.size(); ++nID) {
        traversedNodes[nID].push_back(nodes[nID]);
        lastNodes[nID] = nodes[nID];
    }
    if (totalLevel == 0) {
        size_t num = 1;
        for (auto &group : groups) {
            if (group.size() == 1) {
                VertexID nID2 = group[0];
                num *= lastNodes[nID2] -> numTuples(visited);
            }
            else {
                std::vector<std::vector<VertexID>> vsets(group.size());
                for (int i = 0; i < group.size(); ++i) {
                    VertexID nID2 = group[i];
                    vsets[i].reserve(lastNodes[nID2]->nodeChild.size());
                    for (int j = 0; j < lastNodes[nID2]->nodeChild.size(); ++j) {
                        VertexID v2 = lastNodes[nID2]->nodeChild[j]->value;
                        if (!visited[v2]) vsets[i].push_back(v2);
                    }
                }

                num *= numUniqueCombine(vsets);
            }
        }

        count += num;
    }
    else {
        VertexID firstLevelNID = depthToNID[0];
        numBranches[0] = traversedNodes[firstLevelNID][0]->nodeChild.size();
        int depth = 0;
        while (depth >= 0) {
            while (poses[depth] < numBranches[depth]) {
                VertexID nID = depthToNID[depth];
                VertexID v = traversedNodes[nID].back()->nodeChild[poses[depth]]->value;
                ui childPos = poses[depth];
                ++poses[depth];
                if (visited[v]) continue;
                visited[v] = true;
//                partMatch[order[mappingSize + depth]] = v;
                traversedNodes[nID].push_back(traversedNodes[nID].back()->nodeChild[childPos]);
                lastNodes[nID] = traversedNodes[nID].back();
                if (depth + 1 < totalLevel) {
                    ++depth;
                    poses[depth] = 0;
                    nID = depthToNID[depth];
                    numBranches[depth] = lastNodes[nID]->nodeChild.size();
                }
                else {
                    size_t num = 1;
                    for (auto &group : groups) {
                        if (group.size() == 1) {
                            VertexID nID2 = group[0];
                            num *= lastNodes[nID2] -> numTuples(visited);
                        }
                        else {
                            std::vector<std::vector<VertexID>> vsets(group.size());
                            for (int i = 0; i < group.size(); ++i) {
                                VertexID nID2 = group[i];
                                vsets[i].reserve(lastNodes[nID2]->nodeChild.size());
                                for (int j = 0; j < lastNodes[nID2]->nodeChild.size(); ++j) {
                                    VertexID v2 = lastNodes[nID2]->nodeChild[j]->value;
                                    if (!visited[v2]) vsets[i].push_back(v2);
                                }
                            }

                            num *= numUniqueCombine(vsets);
                        }
                    }

                    count += num;
                    traversedNodes[nID].pop_back();
                    lastNodes[nID] = traversedNodes[nID].back();
                    visited[v] = false;
                }
            }
            --depth;
            if (depth >= 0) {
                VertexID nID = depthToNID[depth];
                VertexID v = lastNodes[nID] -> value;
                traversedNodes[nID].pop_back();
                lastNodes[nID] = traversedNodes[nID].back();
                visited[v] = false;
            }
        }
    }
}

void
produceResult(std::vector<std::vector<VertexID>> &result, const HyperTree &t, VertexID *partMatch,
              const std::vector<TrieNode *> &nodes, bool *visited) {
    ui numNodes = nodes.size();
    ++gNumCartesian;
    // First, collect all sets from the compression of each node
    std::vector<std::vector<VertexID>> sets;
//    for (ui i = 0; i < numNodes; ++i) {
//        if (nodes[i]->compression != nullptr) {
//            for (ui j = 0; j < nodes[i]->compression->num; ++j) {
//                sets.emplace_back(nodes[i]->compression->data[j], nodes[i]->compression->data[j] + nodes[i]->compression->length[j]);
//                if (nodes[i]->compression->length[j] == 0) return;
//            }
//        }
//    }
    if (sets.empty()) {
#ifdef COLLECT_RESULT
        result.emplace_back(partMatch, partMatch + t.numAttributes);
#endif
#ifdef COLLECT_STATISTICS
        ++gNumResult;
#endif
        return;
    }
    std::vector<ui> indices(sets.size(), 0);
    std::vector<VertexID> attributes;
    for (VertexID nID = 0; nID < numNodes; ++nID) {
        const HyperNode &tau = t.nodes[nID];
        for (ui i = tau.numAttributes - t.compressionSizes[nID]; i < tau.numAttributes; ++i) {
            attributes.push_back(tau.attributes[i]);
        }
    }
    // Cartesian product loop
    while (true) {
        bool iso = true;
        for (size_t i = 0; i < sets.size(); ++i) {
            if (!visited[sets[i][indices[i]]]) {
                visited[sets[i][indices[i]]] = true;
                partMatch[attributes[i]] = sets[i][indices[i]];
            }
            else {
                iso = false;
                for (int j = 0; j < i; ++j) {
                    visited[partMatch[attributes[j]]] = false;
                }
                break;
            }
        }

        if (iso) {
            for (int i = 0; i < sets.size(); ++i) {
                visited[partMatch[attributes[i]]] = false;
            }
#ifdef COLLECT_RESULT
            result.emplace_back(partMatch, partMatch + t.numAttributes);
#endif
#ifdef COLLECT_STATISTICS
            ++gNumResult;
#endif
        }

        // Move to the next element in the Cartesian product
        int updateIndex = sets.size() - 1;
        while (updateIndex >= 0) {
            if (++indices[updateIndex] < sets[updateIndex].size()) {
                break;
            }
            indices[updateIndex] = 0;
            updateIndex--;
        }
        if (updateIndex < 0) {
            break;
        }
    }
}

void storeMatches(const std::vector<std::vector<VertexID>> &result, std::ofstream &outStream) {
    outStream << result.size() << std::endl;
    for (int i = 0; i < result.size(); ++i) {
        for (int j = 0; j < result[i].size(); ++j) {
            outStream << result[i][j] << " ";
        }
        outStream << std::endl;
    }
}

void buildTrie(std::vector<std::vector<VertexID>> &tuples, TrieNode *root, const std::vector<VertexID> &order) {
    for (int i = 0; i < root->nodeChild.size(); ++i) {
        delete root->nodeChild[i];
    }
    root -> nodeChild.clear();
    for (int i = 0; i < tuples.size(); ++i) {
        std::vector<VertexID> tmp = tuples[i];
        for (int j = 0; j < order.size(); ++j) {
            tuples[i][j] = tmp[order[j]];
        }
    }
    std::sort(tuples.begin(), tuples.end());
    for (const auto &tuple : tuples)
        root->addMatch(tuple, 0);
    tuples.clear();
}

void executeTwoBag(std::vector<std::vector<VertexID>> &result, CandidateSpace &cs, const HyperTree &t,
                   std::vector<TrieNode *> &nodes, bool *visited) {
    VertexID *partMatch = new VertexID[t.numAttributes];
    VertexID ***nodeCandidates = new VertexID **[nodes.size()];
    ui **nodeCandCount = new ui *[nodes.size()];
    std::vector<std::vector<ui>> nodePoses(nodes.size());
    for (VertexID nID = 0; nID < nodes.size(); ++nID) {
        const HyperNode &tau = t.nodes[nID];
        nodeCandidates[nID]= new VertexID *[tau.numAttributes];
        nodeCandCount[nID] = new ui[tau.numAttributes];
        nodePoses[nID] = std::vector<ui>(tau.numAttributes, 0);
        for (int i = 0; i < tau.numAttributes; ++i) {
            VertexID u = tau.attributes[i];
            nodeCandidates[nID][i] = new VertexID [cs.candidateSet[u].size()];
            nodeCandCount[nID][i] = 0;
        }
    }
    if (t.defaultPartition.empty()) {
        for (VertexID nID = 0; nID < t.numNodes; ++nID) {
            nodeJoinScope(t, nID, cs, visited, partMatch, nodeCandidates[nID], nodeCandCount[nID], nodePoses[nID]);
        }
        return;
    }
    VertexID **pCandidates = new VertexID *[t.defaultPartition.size()];
    ui *pCandCount = new ui[t.defaultPartition.size()];
    const std::vector<std::vector<VertexID>> &vertexParents = t.attributesBefore;
    for (int i = 0; i < t.defaultPartition.size(); ++i) {
        VertexID u = t.defaultPartition[i];
        pCandidates[i] = new VertexID[cs.candidateSet[u].size()];
        pCandCount[i] = 0;
    }
    std::vector<ui> pPoses(t.defaultPartition.size(), 0);

    VertexID u = t.defaultPartition[0];
    pCandidates[0] = cs.candidateSet[u].data();
    pCandCount[0] = cs.candidateSet[u].size();
    VertexID **neighbors = new VertexID *[t.defaultPartition.size()];
    ui *neighborCount = new ui[t.defaultPartition.size()];
    int depth = 0;
    while (depth >= 0) {
        while (pPoses[depth] < pCandCount[depth]) {
            VertexID v = pCandidates[depth][pPoses[depth]];
            ++pPoses[depth];
            if (visited[v]) continue;
            visited[v] = true;
            partMatch[t.defaultPartition[depth]] = v;
            for (VertexID nID: t.nodesAtStep[depth]) {
                for (int i = 0; i < nodePoses[nID].size(); ++i)
                    nodePoses[nID][i] = 0;
                nodeJoinScope(t, nID, cs, visited, partMatch, nodeCandidates[nID], nodeCandCount[nID], nodePoses[nID]);
            }
            if (depth < t.defaultPartition.size() - 1) {
                // match the next partition attribute
                ++depth;
                pPoses[depth] = 0;
                u = t.defaultPartition[depth];
                const std::vector<VertexID> &parents = vertexParents[depth];
                if (!parents.empty()) {
                    for (int i = 0; i < parents.size(); ++i) {
                        VertexID pU = parents[i];
                        neighbors[i] = cs.candidateEdge[pU][partMatch[pU]][u].data();
                        neighborCount[i] = cs.candidateEdge[pU][partMatch[pU]][u].size();
                    }
                    ComputeSetIntersection::LeapfrogJoin(neighbors, neighborCount, parents.size(), pCandidates[depth], pCandCount[depth]);
                }
                else {
                    VertexID cartesianParent = t.cartesianParent[depth];
                    const std::vector<VertexID> &path = cs.reconstructPath(cartesianParent, u);
                    optimizedCartesianProduct(cs, partMatch[cartesianParent], path, pCandidates[depth], pCandCount[depth]);
                }
            }
            else {
                visited[v] = false;
            }
        }
        --depth;
        if (depth >= 0) {
            visited[partMatch[t.defaultPartition[depth]]] = false;
        }
    }
    delete[] partMatch;
    delete[] neighbors;
    for (int i = 0; i < nodes.size(); ++i) {
        delete[] nodeCandCount[i];
        for (int j = 0; j < t.nodes[i].numAttributes; ++j) {
            delete[] nodeCandidates[i][j];
        }
        delete[] nodeCandidates[i];
    }
    delete[] nodeCandidates;
    delete[] nodeCandCount;
}

void getExample(const Graph &data, std::map<uint64_t, ui> &distribution) {
    ui n = data.getNumVertices();
    VertexID **candidates = new VertexID *[6];
    VertexID **neighbors = new VertexID *[6];
    for (int i = 0; i < 6; ++i) {
        candidates[i] = new VertexID [n];
    }
    ui *candCount = new ui[6];
    ui *neighborCounts = new ui[6];
    bool *visited = new bool[n];
    memset(visited, false, sizeof(bool) * n);
    for (VertexID v1 = 0; v1 < n; ++v1) {
        visited[v1] = true;
        ui numV1Nbr;
        VertexID *v1Nbr = data.getNeighbors(v1, numV1Nbr);
        for (int i2 = 0; i2 < numV1Nbr; ++i2) {
            VertexID v2 = v1Nbr[i2];
            if (v2 < v1) continue;
            visited[v2] = true;
            ui numV2Nbr;
            VertexID *v2Nbr = data.getNeighbors(v2, numV2Nbr);
            neighbors[0] = v1Nbr;
            neighbors[1] = v2Nbr;
            neighborCounts[0] = numV1Nbr;
            neighborCounts[1] = numV2Nbr;
            ComputeSetIntersection::LeapfrogJoin(neighbors, neighborCounts, 2, candidates[2], candCount[2]);
            for (int i3 = 0; i3 < candCount[2]; ++i3) {
                VertexID v3 = candidates[2][i3];
                if (visited[v3]) continue;
                visited[v3] = true;
                for (int i4 = i3 + 1; i4 < candCount[2]; ++i4) {
                    VertexID v4 = candidates[2][i4];
                    if (visited[v4]) continue;
                    visited[v4] = true;
                    ui numV3Nbr, numV4Nbr;
                    VertexID *v3Nbr = data.getNeighbors(v3, numV3Nbr);
                    VertexID *v4Nbr = data.getNeighbors(v4, numV4Nbr);
                    neighbors[0] = v1Nbr;
                    neighbors[1] = v3Nbr;
                    neighbors[2] = v4Nbr;
                    neighborCounts[0] = numV1Nbr;
                    neighborCounts[1] = numV3Nbr;
                    neighborCounts[2] = numV4Nbr;
                    ComputeSetIntersection::LeapfrogJoin(neighbors, neighborCounts, 3, candidates[4], candCount[4]);
                    for (int i5 = 0; i5 < candCount[4]; ++i5) {
                        VertexID v5 = candidates[4][i5];
                        if (visited[v5]) continue;
                        visited[v5] = true;
                        neighbors[0] = v2Nbr;
                        neighbors[1] = v3Nbr;
                        neighbors[2] = v4Nbr;
                        neighborCounts[0] = numV2Nbr;
                        neighborCounts[1] = numV3Nbr;
                        neighborCounts[2] = numV4Nbr;
                        ComputeSetIntersection::LeapfrogJoin(neighbors, neighborCounts, 3, candidates[5], candCount[5]);
                        if (candCount[5] == 4 && candidates[5][2] == candidates[5][3]) {
                            printArray(neighbors, neighborCounts, 3);
                        }
                        for (int i6 = 0; i6 < candCount[5]; ++i6) {
                            VertexID v6 = candidates[5][i6];
                            if (visited[v6]) continue;
                            uint64_t id1 = data.getVertexLabel(v3) + data.getVertexLabel(v1) * 100 +
                                    data.getVertexLabel(v2) * 1e4 + data.getVertexLabel(v4) * 1e6 +
                                    data.getVertexLabel(v5) * 1e8 + data.getVertexLabel(v6) * 1e10;
                            if (distribution.find(id1) == distribution.end()) distribution[id1] = 1;
                            else distribution[id1] += 1;
                            if (id1 == 0) ++gNumResult;
                            uint64_t id2 = data.getVertexLabel(v4) + data.getVertexLabel(v1) * 100 +
                                           data.getVertexLabel(v2) * 1e4 + data.getVertexLabel(v3) * 1e6 +
                                           data.getVertexLabel(v5) * 1e8 + data.getVertexLabel(v6) * 1e10;
                            if (distribution.find(id2) == distribution.end()) distribution[id2] = 1;
                            else distribution[id2] += 1;
                            uint64_t id3 = data.getVertexLabel(v3) + data.getVertexLabel(v2) * 100 +
                                           data.getVertexLabel(v1) * 1e4 + data.getVertexLabel(v4) * 1e6 +
                                           data.getVertexLabel(v6) * 1e8 + data.getVertexLabel(v5) * 1e10;
                            if (distribution.find(id3) == distribution.end()) distribution[id3] = 1;
                            else distribution[id3] += 1;
                            uint64_t id4 = data.getVertexLabel(v4) + data.getVertexLabel(v2) * 100 +
                                           data.getVertexLabel(v1) * 1e4 + data.getVertexLabel(v3) * 1e6 +
                                           data.getVertexLabel(v6) * 1e8 + data.getVertexLabel(v5) * 1e10;
                            if (distribution.find(id4) == distribution.end()) distribution[id4] = 1;
                            else distribution[id4] += 1;
                        }
                        visited[v5] = false;
                    }
                    visited[v4] = false;
                }
                visited[v3] = false;
            }
            visited[v2] = false;
        }
        visited[v1] = false;
    }
    for (int i = 0; i < 6; ++i) delete[] candidates[i];
    delete[] neighbors;
    delete[] candidates;
    delete[] candCount;
    delete[] visited;
}