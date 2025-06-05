//
// Created by Qiyan LI on 2024/5/22.
//

#include "optimizer.h"

void
simplePlan(const Graph &query, HyperTree &t, CandidateSpace &cs, PrefixNode *&pt,
           const std::vector<SubsetStructure> &dpStructures, bool *visited, VertexID *partMatch,
           VertexID **candidates, ui *candCount, double &minCost, PrefixNode *bestPT, bool share) {
    std::vector<std::vector<VertexID>> sharedAttrs(t.numNodes);
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        for (int i = 0; i < t.nodes[nID].numAttributes; ++i) {
            VertexID u = t.nodes[nID].attributes[i];
            if (t.v2n[u].size() > 1) {
                sharedAttrs[nID].push_back(u);
            }
        }
        std::sort(sharedAttrs[nID].begin(), sharedAttrs[nID].end());
    }
    ui numNodes = t.numNodes;
    if (t.newGlobalNode) --numNodes;
    std::vector<std::vector<VertexID>> bestOrders;
    minCost = 0.0;
    for (VertexID nID = 0; nID < numNodes; ++nID) {
        minCost += dpStructures[nID].getOptCost();
        bestOrders.push_back(dpStructures[nID].getOptOrder());
    }
//    if (t.newGlobalNode) bestOrders.push_back(globalOrder(query, t, cs, std::vector<VertexID>()));
    if (t.newGlobalNode && bestPT != nullptr) {
        std::vector<VertexID> lastShared;
        VertexID lastID = pt->getBagsBelow().back();
        for (int i = 0; i < t.nodes[lastID].numAttributes; ++i) {
            VertexID u = t.nodes[lastID].attributes[i];
            if (t.v2n[u].size() >= 2) lastShared.push_back(u);
            else break;
        }
    }
    if (!share) {
        pt = new PrefixNode(99);
        for (VertexID nID = 0; nID < t.numNodes; ++nID) {
            pt->nIDsToCall.push_back(nID);
            if (nID != t.numNodes - 1) pt->nIDsToBuild.push_back(nID);
        }
        if (t.newGlobalNode) {
            bestOrders.push_back(globalOrder(query, t, cs, std::vector<VertexID>()));
        }
    }
    else {
        std::unordered_set<PrefixNode *, PrefixNodePtrHash, PrefixNodePtrEqual> visitedPT;
        bool exist;
        pt = buildPrefixTree(bestOrders, query, cs.dist, visitedPT, exist);
        if (t.newGlobalNode) {
            pt->nIDsToCall.push_back(t.numNodes - 1);
            bestOrders.push_back(globalOrder(query, t, cs, std::vector<VertexID>()));
        }
        pt->initNIDsToBuild(t.numNodes);
        pt->initPoses(sharedAttrs, query, cs.dist);
        std::vector<VertexID> attrsInPath;
        int depth = 0;
        const PrefixNode *pn = pt;
        ui height = pn -> getHeight();
        std::vector<ui> childPoses(height, 0);
        std::vector<PrefixNode *> nodes(height, nullptr);
        std::vector<bool> computed(query.getNumVertices(), false);
        while (depth >= 0) {
            while (childPoses[depth] < pn -> children.size()) {
                PrefixNode *current = pn -> children[childPoses[depth]];
                ++childPoses[depth];
                VertexID u = current -> u;
                uint64_t id = getSubsetID(attrsInPath);
                double cost = 0.0;
                ui num = current->getBagsBelow().size() - 1;
                for (VertexID nID : current -> getBagsBelow()) {
                    if (t.newGlobalNode && nID == t.numNodes - 1)
                        --num;
                }
                if (num != 0) {
                    if (depth == 0) cost = cs.candidateSet[u].size();
                    else {
                        double card = subsetToCard[id];
                        double listSize = 0.0;
                        std::vector<VertexID> vertexParents;
                        for (VertexID u2: attrsInPath) {
                            if (query.getEdgeID(u, u2) != -1) vertexParents.push_back(u2);
                        }
                        if (!vertexParents.empty()) {
                            for (VertexID u2: vertexParents) {
                                uint64_t edgeID = (1 << u) + (1 << u2);
                                listSize += subsetToCard[edgeID] / subsetToCard[1 << u2];
                            }
                        }
                        else {
                            std::vector<VertexID> subsetContent = attrsInPath;
                            std::sort(subsetContent.begin(), subsetContent.end());
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
                        cost = card * listSize;
                    }
                    minCost -= cost * num;
                }
                attrsInPath.push_back(u);
                nodes[depth] = current;
                if (!current -> children.empty()) {
                    ++depth;
                    childPoses[depth] = 0;
                }
                else {
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
    }
    buildFromPrefixTree(pt, bestOrders, t, sharedAttrs, query, cs);
    setTDExtention(t, query);
}

void setTDExtention(HyperTree &t, const Graph &query) {
#ifdef ALL_LEVEL
    t.extendLevel = t.defaultPartition.size();
    if (!t.newGlobalNode) {
        VertexID nID = t.numNodes - 1;
        for (int i = t.nodes[nID].prefixSize; i < t.nodes[nID].numAttributes; ++i) {
            VertexID u = t.nodes[nID].attributes[i];
            if (t.v2n[u].size() > 1)
                t.extendLevel = i + 1;
        }
        for (int i = 0; i < t.extendLevel; ++i) {
            VertexID u = t.nodes[nID].attributes[i];
            if (std::find(t.defaultPartition.begin(), t.defaultPartition.end(), u) == t.defaultPartition.end())
                t.defaultPartition.push_back(u);
        }
    }
    else {
        t.extendLevel = t.nodes[t.numNodes - 1].numAttributes;
        for (int i = 0; i < t.extendLevel; ++i) {
            VertexID u = t.nodes[t.numNodes - 1].attributes[i];
            if (std::find(t.defaultPartition.begin(), t.defaultPartition.end(), u) == t.defaultPartition.end())
                t.defaultPartition.push_back(u);
        }
    }
#endif
    t.buildTraverseStruct(query);
}

void
simplePlan(const Graph &query, HyperTree &t, CandidateSpace &cs, PrefixNode *&pt,
           double &minCost) {
    std::vector<std::vector<VertexID>> sharedAttrs(t.numNodes), nodeAttrs(t.numNodes);
    std::vector<int> repetitions(query.getNumVertices(), 1);
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        for (int i = 0; i < t.nodes[nID].numAttributes; ++i) {
            VertexID u = t.nodes[nID].attributes[i];
            nodeAttrs[nID].push_back(u);
            if (t.v2n[u].size() > 1) {
                sharedAttrs[nID].push_back(u);
            }
        }
        std::sort(sharedAttrs[nID].begin(), sharedAttrs[nID].end());
    }
    ui numNodes = t.numNodes;
    if (t.newGlobalNode) --numNodes;
    std::vector<std::vector<VertexID>> bestOrders;
    minCost = 0.0;
    for (VertexID nID = 0; nID < numNodes; ++nID) {
        std::vector<VertexID> order = simpleOrder(query, cs, nodeAttrs[nID], repetitions);
        bestOrders.push_back(order);
    }
    pt = new PrefixNode(99);
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        pt->nIDsToCall.push_back(nID);
        if (nID != t.numNodes - 1) pt->nIDsToBuild.push_back(nID);
    }
    if (t.newGlobalNode) {
        bestOrders.push_back(globalOrder(query, t, cs, std::vector<VertexID>()));
    }
    buildFromPrefixTree(pt, bestOrders, t, sharedAttrs, query, cs);
    setTDExtention(t, query);
}

void bagPermutationDP(const Graph &query, HyperTree &t, CandidateSpace &cs, PrefixNode *&pt,
                      vector<SubsetStructure> &dpStructures, double &bestCost, int reorder, bool connected,
                      bool globalOrderType) {
    ui n = query.getNumVertices();
    ui numNodes = t.numNodes;
    if (t.newGlobalNode) --numNodes;
    std::vector<std::vector<VertexID>> sharedAttrs(numNodes);
    for (VertexID nID = 0; nID < numNodes; ++nID) {
        for (int i = 0; i < t.nodes[nID].numAttributes; ++i) {
            VertexID u = t.nodes[nID].attributes[i];
            if (t.v2n[u].size() > 1) {
                sharedAttrs[nID].push_back(u);
            }
        }
        std::sort(sharedAttrs[nID].begin(), sharedAttrs[nID].end());
    }
    bagPermuteStructure bs;
    bs.init(sharedAttrs, dpStructures, query, false, connected);
    std::vector<VertexID> bagIDs;
    for (VertexID nID = 0; nID < sharedAttrs.size(); ++nID) bagIDs.push_back(nID);
    std::unordered_map<Permutation, double> piToCost;
    for (ui k = 2; k <= sharedAttrs.size(); ++k) {
        std::vector<uint64_t> bagSets, cc;
        std::vector<VertexID> current;
        generateSubsetsK(bagIDs, k, 0, current, bagSets, cc, query, false);
        for (uint64_t bagSet: bagSets) {
            uint64_t copyID = bagSet;
            VertexID oneBag = 0;
            for (ui i = 0; i < n; ++i) {
                if (copyID & 1) {
                    oneBag = i;
                    break;
                }
                copyID >>= 1;
            }
            uint64_t remainingID = bagSet - (1 << oneBag);
            bs.intersections[bagSet] = bs.intersections[remainingID] & bs.intersections[1 << oneBag];
        }
        for (uint64_t bagSet: bagSets) {
            std::vector<VertexID> maxPrefix = getSubsetFromID(bs.intersections[bagSet], n);
            std::vector<VertexID> bags = getSubsetFromID(bagSet, numNodes);
            std::vector<uint64_t> subsets1;
            generateSubsets(bags, subsets1);
            for (ui k2 = 0; k2 <= maxPrefix.size(); ++k2) {
                std::vector<uint64_t> prefixIDs;
                current.clear();
                generateSubsetsK(maxPrefix, k2, 0, current, prefixIDs, dpStructures[bags[0]].cc, query, connected);
                for (uint64_t prefixID: prefixIDs) {
                    std::vector<VertexID> prefix = getSubsetFromID(prefixID, query.getNumVertices());
                    do {
                        if (connected && !orderConnectivity(query, prefix)) continue;
                        std::vector<double> costs, prefixSum(prefix.size());
                        costs = computeCost(prefix, query, cs);
                        double sharedCost = 0.0;
                        for (int i = 0; i < prefix.size(); ++i) {
                            if (i == 0) prefixSum[0] = costs[0];
                            else prefixSum[i] = prefixSum[i - 1] + costs[i];
                            sharedCost += costs[i];
                        }
                        uint64_t minID1, minID2;
                        double minCost = std::numeric_limits<double>::max();
                        Permutation fullPi = encodePermutation(prefix);
                        for (uint64_t subset1: subsets1) {
                            uint64_t subset2 = bagSet - subset1;
                            if (subset1 == 0 || subset2 == 0 || subset1 > subset2) continue;
                            double cost = bs.indepCosts[subset1][fullPi] + bs.indepCosts[subset2][fullPi];
                            if (cost < minCost) {
                                minCost = cost;
                                minID1 = subset1;
                                minID2 = subset2;
                            }
                        }
                        minCost += sharedCost;
                        Permutation pi;
                        if (k == sharedAttrs.size()) piToCost[pi] = 0.0;
                        std::vector<std::vector<VertexID>> localOrder(numNodes);
                        if (bs.indepCosts[bagSet].find(pi) == bs.indepCosts[bagSet].end() || minCost < bs.indepCosts[bagSet][pi]) {
                            std::vector<VertexID> bags1 = getSubsetFromID(minID1, numNodes);
                            std::vector<VertexID> bags2 = getSubsetFromID(minID2, numNodes);
                            const std::vector<std::vector<VertexID>> &localOrder1 = bs.localOrders[minID1][fullPi];
                            const std::vector<std::vector<VertexID>> &localOrder2 = bs.localOrders[minID2][fullPi];
                            for (VertexID nID: bags1) {
                                localOrder[nID] = prefix;
                                for (VertexID u2 : localOrder1[nID]) localOrder[nID].push_back(u2);
                            }
                            for (VertexID nID: bags2) {
                                localOrder[nID] = prefix;
                                for (VertexID u2 : localOrder2[nID]) localOrder[nID].push_back(u2);
                            }
                            bs.indepCosts[bagSet][pi] = minCost;
                            bs.localOrders[bagSet][pi] = localOrder;
                        }
                        for (int i = 0; i < prefix.size(); ++i) {
                            VertexID u = prefix[i];
                            uint32_t value = u & 0x1F;
                            size_t startBit = i * 5 + 5;
                            for (size_t bit = 0; bit < 5; ++bit) {
                                if ((i + 1) & (1 << bit)) pi.set(bit);
                                else pi.reset(bit);
                            }
                            for (size_t bit = 0; bit < 5; ++bit) {
                                if (value & (1 << bit)) {
                                    pi.set(startBit + bit);
                                }
                            }
                            minCost -= costs[i];
                            if (k == sharedAttrs.size()) piToCost[pi] = prefixSum[i];
                            if (bs.indepCosts[bagSet].find(pi) == bs.indepCosts[bagSet].end() || minCost < bs.indepCosts[bagSet][pi]) {
                                std::vector<VertexID> bags1 = getSubsetFromID(minID1, numNodes);
                                std::vector<VertexID> bags2 = getSubsetFromID(minID2, numNodes);
                                const std::vector<std::vector<VertexID>> &localOrder1 = bs.localOrders[minID1][fullPi];
                                const std::vector<std::vector<VertexID>> &localOrder2 = bs.localOrders[minID2][fullPi];
                                for (VertexID nID: bags1) {
                                    localOrder[nID].clear();
                                    for (int j = i + 1; j < prefix.size(); ++j) localOrder[nID].push_back(prefix[j]);
                                    for (VertexID u2 : localOrder1[nID]) localOrder[nID].push_back(u2);
                                }
                                for (VertexID nID: bags2) {
                                    localOrder[nID].clear();
                                    for (int j = i + 1; j < prefix.size(); ++j) localOrder[nID].push_back(prefix[j]);
                                    for (VertexID u2 : localOrder2[nID]) localOrder[nID].push_back(u2);
                                }
                                bs.indepCosts[bagSet][pi] = minCost;
                                bs.localOrders[bagSet][pi] = localOrder;
                            }
                        }
                    } while (std::next_permutation(prefix.begin(), prefix.end()));
                }
            }
        }
    }
    Permutation bestPi;
    bestCost = std::numeric_limits<double>::max();
    uint64_t allBag = getSubsetID(bagIDs);
    for (auto it = bs.indepCosts[allBag].begin(); it != bs.indepCosts[allBag].end(); ++it) {
        Permutation pi = it->first;
        double prefixCost = piToCost[pi];
        double cost = it->second + prefixCost;
        if (cost < bestCost) {
            bestCost = cost;
            bestPi = pi;
        }
    }
    std::vector<std::vector<VertexID>> bestOrders = bs.localOrders[allBag][bestPi];
    std::vector<VertexID> prefixOrder = decodePermutation(bestPi);
    std::vector<std::vector<VertexID>> old = bestOrders;
    for (VertexID nID = 0; nID < dpStructures.size(); ++nID) {
        bestOrders[nID] = prefixOrder;
        for (VertexID u : old[nID]) bestOrders[nID].push_back(u);
    }
    std::unordered_set<PrefixNode *, PrefixNodePtrHash, PrefixNodePtrEqual> visitedPT;
    bool exist;
    pt = buildPrefixTree(bestOrders, query, cs.dist, visitedPT, exist);
    buildFromPrefixTree(pt, bestOrders, t, sharedAttrs, query, cs);
    changeNIDsToCall(query, t, pt);
    reorderBags(t, pt, dpStructures, reorder, query, cs);
    if (t.newGlobalNode) {
        bestOrders.emplace_back();
        VertexID maxCostID;
//        double maxCost = 0.0;
//        for (VertexID nID2 = 0; nID2 < numNodes; ++nID2) {
//            if (dpStructures[nID2].getOptCost() > maxCost) {
//                maxCost = dpStructures[nID2].getOptCost();
//                maxCostID = nID2;
//            }
//        }
        maxCostID = pt->getBagsBelow().back();
        std::vector<ui> prefixSizes;
        addGlobal(query, t, pt, t.numNodes - 1, maxCostID, cs, bestOrders, prefixSizes, globalOrderType);
        sharedAttrs.push_back(bestOrders.back());
        std::sort(sharedAttrs.back().begin(), sharedAttrs.back().end());
    }
    pt->initNIDsToBuild(t.numNodes);
    pt->initPoses(sharedAttrs, query, cs.dist);
    buildFromPrefixTree(pt, bestOrders, t, sharedAttrs, query, cs);
    setTDExtention(t, query);
}

void optCostPlan(const Graph &query, HyperTree &t, CandidateSpace &cs, PrefixNode *&pt,
                 vector<SubsetStructure> &dpStructures, double &bestCost, int reorder, bool connected,
                 bool globalOrderType) {
    ui n = query.getNumVertices();
    ui numNodes = t.numNodes;
    if (t.newGlobalNode) --numNodes;
    std::vector<std::vector<VertexID>> sharedAttrs(numNodes);
    for (VertexID nID = 0; nID < numNodes; ++nID) {
        for (int i = 0; i < t.nodes[nID].numAttributes; ++i) {
            VertexID u = t.nodes[nID].attributes[i];
            if (t.v2n[u].size() > 1) {
                sharedAttrs[nID].push_back(u);
            }
        }
        std::sort(sharedAttrs[nID].begin(), sharedAttrs[nID].end());
    }
    bagPermuteStructure bs;
    bs.initBasicOrder(sharedAttrs, dpStructures, query, false, connected);
    std::vector<VertexID> bagIDs;
    for (VertexID nID = 0; nID < sharedAttrs.size(); ++nID) bagIDs.push_back(nID);
    std::unordered_map<Permutation, double> piToCost;
    for (ui k = 2; k <= sharedAttrs.size(); ++k) {
        std::vector<uint64_t> bagSets, cc;
        std::vector<VertexID> current;
        generateSubsetsK(bagIDs, k, 0, current, bagSets, cc, query, false);
        for (uint64_t bagSet: bagSets) {
            uint64_t copyID = bagSet;
            VertexID oneBag = 0;
            for (ui i = 0; i < n; ++i) {
                if (copyID & 1) {
                    oneBag = i;
                    break;
                }
                copyID >>= 1;
            }
            uint64_t remainingID = bagSet - (1 << oneBag);
            bs.intersections[bagSet] = bs.intersections[remainingID] & bs.intersections[1 << oneBag];
        }
        for (uint64_t bagSet: bagSets) {
            std::vector<VertexID> maxPrefix = getSubsetFromID(bs.intersections[bagSet], n);
            std::vector<VertexID> bags = getSubsetFromID(bagSet, numNodes);
            std::vector<uint64_t> subsets1;
            generateSubsets(bags, subsets1);
            for (ui k2 = 0; k2 <= maxPrefix.size(); ++k2) {
                std::vector<uint64_t> prefixIDs;
                current.clear();
                generateSubsetsK(maxPrefix, k2, 0, current, prefixIDs, dpStructures[bags[0]].cc, query, connected);
                for (uint64_t prefixID: prefixIDs) {
                    std::vector<VertexID> prefix = getSubsetFromID(prefixID, query.getNumVertices());
                    do {
                        if (connected && !orderConnectivity(query, prefix)) continue;
                        std::vector<double> costs, prefixSum(prefix.size());
                        costs = computeCost(prefix, query, cs);
                        double sharedCost = 0.0;
                        for (int i = 0; i < prefix.size(); ++i) {
                            if (i == 0) prefixSum[0] = costs[0];
                            else prefixSum[i] = prefixSum[i - 1] + costs[i];
                            sharedCost += costs[i];
                        }
                        uint64_t minID1, minID2;
                        double minCost = std::numeric_limits<double>::max();
                        Permutation fullPi = encodePermutation(prefix);
                        for (uint64_t subset1: subsets1) {
                            uint64_t subset2 = bagSet - subset1;
                            if (subset1 == 0 || subset2 == 0 || subset1 > subset2) continue;
                            double cost = bs.indepCosts[subset1][fullPi] + bs.indepCosts[subset2][fullPi];
                            if (cost < minCost) {
                                minCost = cost;
                                minID1 = subset1;
                                minID2 = subset2;
                            }
                        }
                        minCost += sharedCost;
                        Permutation pi;
                        if (k == sharedAttrs.size()) piToCost[pi] = 0.0;
                        std::vector<int> localOrder(numNodes, -1);
                        if (bs.indepCosts[bagSet].find(pi) == bs.indepCosts[bagSet].end() || minCost < bs.indepCosts[bagSet][pi]) {
                            std::vector<VertexID> bags1 = getSubsetFromID(minID1, numNodes);
                            std::vector<VertexID> bags2 = getSubsetFromID(minID2, numNodes);
                            const std::vector<int> &localOrder1 = bs.orderIDs[minID1][fullPi];
                            const std::vector<int> &localOrder2 = bs.orderIDs[minID2][fullPi];
                            for (VertexID nID: bags1) {
                                localOrder[nID] = localOrder1[nID];
                            }
                            for (VertexID nID: bags2) {
                                localOrder[nID] = localOrder2[nID];
                            }
                            bs.indepCosts[bagSet][pi] = minCost;
                            bs.orderIDs[bagSet][pi] = localOrder;
                        }
                        for (int i = 0; i < prefix.size(); ++i) {
                            VertexID u = prefix[i];
                            uint32_t value = u & 0x1F;
                            size_t startBit = i * 5 + 5;
                            for (size_t bit = 0; bit < 5; ++bit) {
                                if ((i + 1) & (1 << bit)) pi.set(bit);
                                else pi.reset(bit);
                            }
                            for (size_t bit = 0; bit < 5; ++bit) {
                                if (value & (1 << bit)) {
                                    pi.set(startBit + bit);
                                }
                            }
                            minCost -= costs[i];
                            if (k == sharedAttrs.size()) piToCost[pi] = prefixSum[i];
                            if (bs.indepCosts[bagSet].find(pi) == bs.indepCosts[bagSet].end() || minCost < bs.indepCosts[bagSet][pi]) {
                                std::vector<VertexID> bags1 = getSubsetFromID(minID1, numNodes);
                                std::vector<VertexID> bags2 = getSubsetFromID(minID2, numNodes);
                                const std::vector<int> &localOrder1 = bs.orderIDs[minID1][fullPi];
                                const std::vector<int> &localOrder2 = bs.orderIDs[minID2][fullPi];
                                for (VertexID nID: bags1) {
                                    localOrder[nID] = localOrder1[nID];
                                }
                                for (VertexID nID: bags2) {
                                    localOrder[nID] = localOrder2[nID];
                                }
                                bs.indepCosts[bagSet][pi] = minCost;
                                bs.orderIDs[bagSet][pi] = localOrder;
                            }
                        }
                    } while (std::next_permutation(prefix.begin(), prefix.end()));
                }
            }
        }
    }
    Permutation bestPi;
    bestCost = std::numeric_limits<double>::max();
    uint64_t allBag = getSubsetID(bagIDs);
    for (auto it = bs.indepCosts[allBag].begin(); it != bs.indepCosts[allBag].end(); ++it) {
        Permutation pi = it->first;
        double prefixCost = piToCost[pi];
        double cost = it->second + prefixCost;
        if (cost < bestCost) {
            bestCost = cost;
            bestPi = pi;
        }
    }
    std::vector<int> bestOrderIDs = bs.orderIDs[allBag][bestPi];
    std::vector<std::vector<VertexID>> bestOrders = std::vector<std::vector<VertexID>>(bestOrderIDs.size());
    for (VertexID nID = 0; nID < dpStructures.size(); ++nID) {
        bestOrders[nID] = bs.basicOrders[nID][bestOrderIDs[nID]];
    }
    std::unordered_set<PrefixNode *, PrefixNodePtrHash, PrefixNodePtrEqual> visitedPT;
    bool exist;
    pt = buildPrefixTree(bestOrders, query, cs.dist, visitedPT, exist);
    buildFromPrefixTree(pt, bestOrders, t, sharedAttrs, query, cs);
    changeNIDsToCall(query, t, pt);
    reorderBags(t, pt, dpStructures, reorder, query, cs);
    if (t.newGlobalNode) {
        bestOrders.emplace_back();
        VertexID maxCostID;
//        double maxCost = 0.0;
//        for (VertexID nID2 = 0; nID2 < numNodes; ++nID2) {
//            if (dpStructures[nID2].getOptCost() > maxCost) {
//                maxCost = dpStructures[nID2].getOptCost();
//                maxCostID = nID2;
//            }
//        }
        maxCostID = pt->getBagsBelow().back();
        std::vector<ui> prefixSizes;
        addGlobal(query, t, pt, t.numNodes - 1, maxCostID, cs, bestOrders, prefixSizes, globalOrderType);
        sharedAttrs.push_back(bestOrders.back());
        std::sort(sharedAttrs.back().begin(), sharedAttrs.back().end());
    }
    pt->initNIDsToBuild(t.numNodes);
    pt->initPoses(sharedAttrs, query, cs.dist);
    buildFromPrefixTree(pt, bestOrders, t, sharedAttrs, query, cs);
    setTDExtention(t, query);
}

void reorderBags(HyperTree &t, PrefixNode *pt, const std::vector<SubsetStructure> &dpStructures, int heuristic,
                 const Graph &query, CandidateSpace &cs) {
    std::map<PrefixNode *, bool> switchable;
    std::map<PrefixNode *, double> score;
    std::vector<PrefixNode *> attributes;
    std::vector<VertexID> nIDs;
    pt->getTraverseOrder(attributes, nIDs, t);
    for (PrefixNode *attr: attributes) {
        switchable[attr] = true;
        score[attr] = 0;
    }
    if (!t.newGlobalNode) {
        PrefixNode *attr = pt;
        switchable[attr] = false;
        while (!attr->children.empty()) {
            attr = attr->children.back();
            switchable[attr] = false;
        }
    }
    bool smallerBetter = true;
    std::map<PrefixNode *, std::vector<VertexID>> attrsInPath;
    std::vector<PrefixNode *> path;
    std::vector<VertexID> attrs;
    int depth = 0;
    PrefixNode *pn = pt;
    ui height = pn -> getHeight();
    std::vector<ui> childPoses = std::vector<ui>(height, 0);
    std::set<VertexID> visitedNID;
    std::map<PrefixNode *, int> numBagsAfter;
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            path.push_back(current);
            attrs.push_back(current->u);
            attrsInPath[current] = attrs;
            for (VertexID nID: current->nIDsToCall) visitedNID.insert(nID);
            const std::vector<VertexID> &below = current->getBagsBelow();
            for (VertexID nID = 0; nID < t.numNodes; ++nID) {
                if (std::find(below.begin(), below.end(), nID) == below.end() && visitedNID.find(nID) ==
                                                                                 visitedNID.end()) {
                    bool includes = false;
                    for (int i = 0; i < t.nodes[nID].numAttributes; ++i) {
                        VertexID u = t.nodes[nID].attributes[i];
                        if (u == current->u) includes = true;
                    }
                    if (includes) {
                        if (numBagsAfter.find(current) == numBagsAfter.end()) numBagsAfter[current] = 1;
                        else ++numBagsAfter[current];
                    }
                }
            }
            if (!current -> children.empty()) {
                ++depth;
                childPoses[depth] = 0;
            }
            else {
                ++childPoses[depth];
                path.pop_back();
                attrs.pop_back();
            }
            if (depth > 0) pn = path[depth - 1];
            else pn = pt;
        }
        --depth;
        if (depth >= 0) {
            ++childPoses[depth];
            path.pop_back();
            attrs.pop_back();
            if (depth == 0) pn = pt;
            else pn = path[depth - 1];
        }
    }
    if (heuristic == 1) {
        smallerBetter = false;
        for (auto it = attributes.rbegin(); it != attributes.rend(); ++it) {
            PrefixNode *attr = *it;
            score[attr] = (double)numBagsAfter[attr];
        }
    }
    if (heuristic == 2) {
        smallerBetter = true;
        for (auto it = attributes.rbegin(); it != attributes.rend(); ++it) {
            PrefixNode *attr = *it;
            for (VertexID nID: attr->nIDsToCall) {
                VertexID id = 0;
                for (int i = 0; i < t.nodes[nID].numAttributes; ++i)
                    id += 1 << t.nodes[nID].attributes[i];
                if (subsetToCard[id] > score[attr]) score[attr] = subsetToCard[id];
            }
            for (PrefixNode *c: attr->children) {
                if (score[c] > score[attr]) score[attr] = score[c];
            }
        }
    }
    if (heuristic == 3) {
        smallerBetter = false;
        for (auto it = attributes.rbegin(); it != attributes.rend(); ++it) {
            PrefixNode *attr = *it;
            const std::vector<VertexID> &attrInPath = attrsInPath[attr];
            for (PrefixNode *c: attr->children) score[attr] += score[c];
            if (attr->u != 99)
                score[attr] += computeCost(attrInPath, query, cs).back();
            for (VertexID nID: attr->nIDsToCall) {
                uint64_t id = getSubsetID(attrInPath);
                score[attr] += dpStructures[nID].getRemainingCost(id, attrInPath.size());
            }
        }
    }
    if (heuristic == 4) {
        smallerBetter = true;
        for (auto it = attributes.rbegin(); it != attributes.rend(); ++it) {
            PrefixNode *attr = *it;
            for (int i = 0; i < attr->children.size(); ++i) {
                PrefixNode *c = attr->children[i];
                score[c] = i;
            }
        }
    }
    // reorder
    for (auto it = attributes.rbegin(); it != attributes.rend(); ++it) {
        PrefixNode *attr = *it;
        std::vector<PrefixNode *> &children = attr->children;
        std::sort(children.begin(), children.end(), [&smallerBetter, &score](PrefixNode *&a, PrefixNode *&b) {
            if (smallerBetter) return score[a] < score[b];
            else return score[a] > score[b];
        });
        for (int i = 0; i < children.size(); ++i) {
            if (!switchable[children[i]]) {
                std::swap(children.back(), children[i]);
                break;
            }
        }
    }
}