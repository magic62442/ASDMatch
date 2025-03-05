//
// Created by anonymous authors on 2024/5/22.
//

#include "optimizer.h"

void simpleTwoStepPlan(const Graph &query, HyperTree &t, const CandidateSpace &cs) {
    t.numAttributes = query.getNumVertices();
    std::vector<VertexID> globalOrder, globalVertices(query.getNumVertices());
    std::vector<int> repetitions(t.numAttributes, 0);
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        std::vector<VertexID> localVertex, localOrder;
        for (int i = 0; i < t.nodes[nID].numAttributes; ++i) {
            localVertex.push_back(t.nodes[nID].attributes[i]);
        }
        localOrder = simpleOrder(query, cs, localVertex, repetitions);
        t.nodes[nID].prefixSize = 0;
        for (int i = 0; i < t.nodes[nID].numAttributes; ++i) {
            t.nodes[nID].attributes[i] = localOrder[i];
        }
    }
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        for (int i = 0; i < t.nodes[nID].numAttributes; ++i) {
            ++repetitions[t.nodes[nID].attributes[i]];
        }
    }
    for (VertexID i = 0; i < query.getNumVertices(); ++i) globalVertices[i] = i;
    t.globalOrder = simpleOrder(query, cs, globalVertices, repetitions);
    t.initPoses(query, cs, true);
}

void dpTwoStepPlan(const Graph &query, HyperTree &t, CandidateSpace &cs, bool *visited, VertexID *partMatch,
                   VertexID **candidates, ui *candCount, const std::vector<SubsetStructure> &dpStructures) {
    t.numAttributes = query.getNumVertices();
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        const SubsetStructure &s = dpStructures[nID];
        std::vector<VertexID> localOrder = s.getOptOrder();
        t.nodes[nID].prefixSize = 0;
        for (int i = 0; i < t.nodes[nID].numAttributes; ++i) {
            t.nodes[nID].attributes[i] = localOrder[i];
        }
    }
    std::vector<VertexID> glocalVertex = globalOrder(query, t, cs, std::vector<VertexID>());
    t.globalOrder = glocalVertex;
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        for (int i = 0; i < t.nodes[nID].numAttributes; ++i) {
            VertexID u = t.nodes[nID].attributes[i];
            if (t.v2n[u].size() == 1) {
                t.globalOrder.push_back(u);
            }
        }
    }
    t.initPoses(query, cs, true);
}

void dpTwoBagPlan(const Graph &query, HyperTree &t, CandidateSpace &cs, std::vector<double> &cost1,
                  std::vector<double> &cost2, bool &refine, bool *visited, VertexID *partMatch, VertexID **candidates,
                  ui *candCount, const SubsetStructure &s0, const SubsetStructure &s1) {
    t.numAttributes = query.getNumVertices();
    std::vector<VertexID> localVertices0(t.nodes[0].attributes, t.nodes[0].attributes + t.nodes[0].numAttributes);
    std::vector<VertexID> localVertices1(t.nodes[1].attributes, t.nodes[1].attributes + t.nodes[1].numAttributes);
    std::vector<VertexID> bestGlobal, bestOrder0, bestOrder1;
    bestOrder0 = s0.getOptOrder();
    bestOrder1 = s1.getOptOrder();
    double minCost = s0.getOptCost() + s1.getOptCost();
    cost1 = {s0.getOptCost(), s1.getOptCost()};
    std::vector<VertexID> sharedAttrs;
    for (VertexID u = 0; u < query.getNumVertices(); ++u) {
        if (t.v2n[u].size() > 1) {
            sharedAttrs.push_back(u);
        }
    }
    std::vector<uint64_t> tmpSubsetIDs;
    std::vector<std::vector<VertexID>> components;
    query.computeConnectedComponents(sharedAttrs, components);
    for (int size = 1; size <= sharedAttrs.size(); ++size) {
        std::vector<VertexID> current;
        generateSubsetsK(sharedAttrs, size, 0, current, tmpSubsetIDs, s0.cc, query, true);
    }
    std::vector<uint64_t> subsetIDs;
    for (uint64_t subsetID: tmpSubsetIDs) {
        if (subsetConnectivity(query, s1.cc, getSubsetFromID(subsetID, query.getNumVertices())))
            subsetIDs.push_back(subsetID);
    }
    cost2 = cost1;
    cost2.push_back(0.0);
    std::map<uint64_t, double> idToCost;
    for (uint64_t subsetID: subsetIDs) {
        ui k = getSubsetFromID(subsetID, query.getNumVertices()).size();
        std::vector<VertexID> order0, order1, prefixOrder;
        double costNode0, costNode1, prefixCost;
        prefixCost = s0.getOptCost(subsetID, k);
        prefixOrder = s0.getOptOrder(subsetID, k);
        s0.prefixPlan(subsetID, k, order0, costNode0);
        s1.prefixPlan(subsetID, k, order1, costNode1);
        double cost = costNode0 + costNode1 + prefixCost;
        idToCost[subsetID] = cost;
        if (cost < minCost) {
            minCost = costNode0 + costNode1 + prefixCost;
            bestGlobal = prefixOrder;
            bestOrder0 = order0;
            bestOrder1 = order1;
            cost2 = {costNode0, costNode1, prefixCost};
        }
    }

    for (int i = 0; i < bestGlobal.size(); ++i) {
        t.nodes[0].attributes[i] = bestGlobal[i];
        t.nodes[1].attributes[i] = bestGlobal[i];
    }
    t.nodes[0].prefixSize = t.nodes[1].prefixSize = bestGlobal.size();
    t.nodes[0].prefix = new VertexID [bestGlobal.size()];
    t.nodes[1].prefix = new VertexID [bestGlobal.size()];
    for (int i = 0; i < bestGlobal.size(); ++i) {
        t.nodes[0].prefix[i] = bestGlobal[i];
        t.nodes[1].prefix[i] = bestGlobal[i];
    }
    for (int i = 0; i < bestOrder0.size(); ++i)
        t.nodes[0].attributes[i + bestGlobal.size()] = bestOrder0[i];
    for (int i = 0; i < bestOrder1.size(); ++i)
        t.nodes[1].attributes[i + bestGlobal.size()] = bestOrder1[i];
    t.defaultPartition = bestGlobal;
    t.globalOrder = bestGlobal;
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        for (int i = 0; i < t.nodes[nID].numAttributes; ++i) {
            VertexID u = t.nodes[nID].attributes[i];
            if (t.v2n[u].size() == 1) {
                t.globalOrder.push_back(u);
            }
        }
    }
    t.initPoses(query, cs, true);
    // check whether this plan can be obtained by refining
//    if (bestGlobal.size() < 2) {
//        refine = true;
//        return;
//    }
//    std::vector<VertexID> copy = bestGlobal;
//    refine = false;
//    do {
//        uint64_t id = 0;
//        double prevCost = cost1[0] + cost1[1];
//        for (int i = 0; i < copy.size(); ++i) {
//            id += 1 << copy[i];
//            if (idToCost[id] > prevCost)
//                break;
//            else if (i == copy.size() - 1) {
//                refine = true;
//                return;
//            }
//            else prevCost = idToCost[id];
//        }
//    } while (std::next_permutation(copy.begin(), copy.end()));
}

void
twoBagMaxSharePlan(const Graph &query, HyperTree &t, CandidateSpace &cs, std::vector<double> &cost, bool *visited,
                   VertexID *partMatch, VertexID **candidates, ui *candCount, const SubsetStructure &s0,
                   const SubsetStructure &s1) {
    t.numAttributes = query.getNumVertices();
    std::vector<VertexID> localVertices0(t.nodes[0].attributes, t.nodes[0].attributes + t.nodes[0].numAttributes);
    std::vector<VertexID> localVertices1(t.nodes[1].attributes, t.nodes[1].attributes + t.nodes[1].numAttributes);
    std::vector<VertexID> bestGlobal, bestOrder0, bestOrder1;
    double costGlobal, cost0, cost1;
    bestOrder0 = s0.getOptOrder();
    bestOrder1 = s1.getOptOrder();
    std::vector<VertexID> sharedAttrs;
    for (VertexID u = 0; u < query.getNumVertices(); ++u) {
        if (t.v2n[u].size() > 1) {
            sharedAttrs.push_back(u);
        }
    }
    std::vector<std::vector<VertexID>> components;
    query.computeConnectedComponents(sharedAttrs, components);
    int largestCCPos = 0;
    ui largestSize = 0;
    for (int i = 0; i < components.size(); ++i) {
        if (components.size() > largestSize) {
            largestSize = components.size();
            largestCCPos = i;
        }
    }
    std::vector<VertexID> largestCC = components[largestCCPos];
    uint64_t sharedID = getSubsetID(largestCC);
    bestGlobal = s0.getOptOrder(sharedID, largestCC.size());
    costGlobal = s0.getOptCost(sharedID, largestCC.size());
    s0.prefixPlan(sharedID, largestCC.size(), bestOrder0, cost0);
    s1.prefixPlan(sharedID, largestCC.size(), bestOrder1, cost1);
    cost = {cost0, cost1, costGlobal};
    for (int i = 0; i < bestGlobal.size(); ++i) {
        t.nodes[0].attributes[i] = bestGlobal[i];
        t.nodes[1].attributes[i] = bestGlobal[i];
    }
    t.nodes[0].prefixSize = t.nodes[1].prefixSize = bestGlobal.size();
    t.nodes[0].prefix = new VertexID [bestGlobal.size()];
    t.nodes[1].prefix = new VertexID [bestGlobal.size()];
    for (int i = 0; i < bestGlobal.size(); ++i) {
        t.nodes[0].prefix[i] = bestGlobal[i];
        t.nodes[1].prefix[i] = bestGlobal[i];
    }
    for (int i = 0; i < bestOrder0.size(); ++i)
        t.nodes[0].attributes[i + bestGlobal.size()] = bestOrder0[i];
    for (int i = 0; i < bestOrder1.size(); ++i)
        t.nodes[1].attributes[i + bestGlobal.size()] = bestOrder1[i];
    t.defaultPartition = bestGlobal;
    t.globalOrder = bestGlobal;
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        for (int i = 0; i < t.nodes[nID].numAttributes; ++i) {
            VertexID u = t.nodes[nID].attributes[i];
            if (t.v2n[u].size() == 1) {
                t.globalOrder.push_back(u);
            }
        }
    }

    t.initPoses(query, cs, true);
}

void twoBagHeuristicPlan(const Graph &query, HyperTree &t, CandidateSpace &cs, std::vector<double> &cost, bool *visited, VertexID *partMatch,
                         VertexID **candidates, ui *candCount) {
    t.numAttributes = query.getNumVertices();
    std::vector<VertexID> localVertices0(t.nodes[0].attributes, t.nodes[0].attributes + t.nodes[0].numAttributes);
    std::vector<VertexID> localVertices1(t.nodes[1].attributes, t.nodes[1].attributes + t.nodes[1].numAttributes);
    std::vector<VertexID> bestGlobal;
    std::vector<ui> num0, num1;
    double bestCost0 = 0.0, bestCost1 = 0.0, bestCostShared = 0.0, minCost = 0.0;
    std::vector<VertexID> bestOrder0 = simpleOrder(query, cs, bestGlobal, num0, localVertices0, t.nodes[0].numAttributes);
    std::vector<VertexID> bestOrder1 = simpleOrder(query, cs, bestGlobal, num1, localVertices1, t.nodes[1].numAttributes);
    std::vector<ui> poses(query.getNumVertices());
    bestCost0 = computeCost(bestGlobal, bestOrder0, query, cs, visited, partMatch, candidates, candCount);
    bestCost1 = computeCost(bestGlobal, bestOrder1, query, cs, visited, partMatch, candidates, candCount);
    minCost = bestCost0 + bestCost1;
    std::vector<ui> nodePriority;
    if (bestCost0 < bestCost1) nodePriority = {1, 0};
    else nodePriority = {0, 1};
    PrefixNode *empty = new PrefixNode(99);
    std::vector<VertexID> sharedAttrs;
    for (VertexID u = 0; u < query.getNumVertices(); ++u) {
        if (t.v2n[u].size() > 1) {
            sharedAttrs.push_back(u);
        }
    }
    std::priority_queue<CandidatePlan> plans;
    for (VertexID u : sharedAttrs) {
        std::vector<VertexID> prefix = {u};
        PrefixNode *newPT = empty->clone();
        PrefixNode *newChild = new PrefixNode(u);
        newPT->children.push_back(newChild);
        int unChangePos0 = 0, unChangePos1 = 0;
        while (unChangePos0 < bestOrder0.size()) {
            if (bestOrder0[unChangePos0] == u) break;
            else ++unChangePos0;
        }
        while (unChangePos1 < bestOrder1.size()) {
            if (bestOrder1[unChangePos1] == u) break;
            else ++unChangePos1;
        }
        std::vector<VertexID> order0, order1;
        std::vector<ui> numBackWard0 = {0}, numBackWard1 = {0};
        if (unChangePos0 == 0) {
            numBackWard0 = num0;
            order0.assign(bestOrder0.begin() + 1, bestOrder0.end());
        }
        else {
            order0 = simpleOrder(query, cs, prefix, numBackWard0, bestOrder0, unChangePos0);
            for (int i = unChangePos0 + 1; i < bestOrder0.size(); ++i)
                numBackWard0.push_back(num0[i]);
        }
        if (unChangePos1 == 0) {
            numBackWard1 = num1;
            order1.assign(bestOrder1.begin() + 1, bestOrder1.end());
        }
        else {
            order1 = simpleOrder(query, cs, prefix, numBackWard1, bestOrder1, unChangePos1);
            for (int i = unChangePos1 + 1; i < bestOrder1.size(); ++i)
                numBackWard1.push_back(num1[i]);
        }
        std::vector<std::vector<VertexID>> orders = {order0, order1};

        std::vector<std::vector<ui>> backwards = {numBackWard0, numBackWard1};
        plans.emplace(newPT, orders, backwards, nodePriority, 2, minCost);
    }
    int numCheckedPlan = 0;
    std::vector<std::vector<VertexID>> components;
    std::vector<uint64_t> cc;
    query.computeConnectedComponents(sharedAttrs, components);
    for (auto & component : components) {
        cc.push_back(getSubsetID(component));
    }
    while (!plans.empty()) {
        CandidatePlan cp = plans.top();
        plans.pop();
        // evaluate the true cost of the candidate plan
        std::vector<VertexID> prefix;
        PrefixNode *pn = cp.pt -> children[0];
        prefix.push_back(pn -> u);
        while (!pn -> children.empty()) {
            pn = pn -> children[0];
            prefix.push_back(pn -> u);
        }
        double costShared = computeCost(std::vector<VertexID>(), prefix, query, cs, visited,
                                        partMatch, candidates, candCount);
        double cost0 = computeCost(prefix, cp.orders[0], query, cs, visited, partMatch, candidates, candCount);
        double cost1 = computeCost(prefix, cp.orders[1], query, cs, visited, partMatch, candidates, candCount);
        double totalCost = costShared + cost0 + cost1;
        if (totalCost < minCost) {
            minCost = totalCost;
            bestCost0 = cost0;
            bestCost1 = cost1;
            bestCostShared = costShared;
            bestOrder0 = cp.orders[0];
            bestOrder1 = cp.orders[1];
            bestGlobal = prefix;
        }
        ++numCheckedPlan;
        if (numCheckedPlan == MAX_NUM_PLAN) break;
        if (totalCost > cp.parentCost) continue;
        // share one more vertex
        std::vector<VertexID> newVertices;
        for (VertexID u: sharedAttrs) {
            bool flag = true;
            for (VertexID u2 : prefix) {
                if (u <= u2) {
                    flag = false;
                    break;
                }
            }
            if (flag) newVertices.push_back(u);
        }
        for (VertexID u : newVertices) {
            std::vector<VertexID> newPrefix = prefix;
            newPrefix.push_back(u);
            if (!subsetConnectivity(query, cc, newPrefix)) continue;
            std::vector<ui> newNumBackNbr0;
            newPrefix = simpleOrder(query, cs, std::vector<VertexID>(), newNumBackNbr0, newPrefix, newPrefix.size());
            std::vector<ui> newNumBackNbr1 = newNumBackNbr0;
            int unChangePos0 = 0;
            while (unChangePos0 < cp.orders[0].size()) {
                if (cp.orders[0][unChangePos0] == u) break;
                else ++unChangePos0;
            }
            int unChangePos1 = 0;
            while (unChangePos1 < cp.orders[1].size()) {
                if (cp.orders[1][unChangePos1] == u) break;
                else ++unChangePos1;
            }
            std::vector<VertexID> newOrder0, newOrder1;
            if (unChangePos0 == 0) {
                newOrder0.assign(cp.orders[0].begin() + 1, cp.orders[0].end());
                newNumBackNbr0 = cp.numBackNbr[0];
            }
            else {
                newOrder0 = simpleOrder(query, cs, newPrefix, newNumBackNbr0, cp.orders[0], unChangePos0);
                for (int i = unChangePos0 + 1; i < cp.numBackNbr[0].size(); ++i)
                    newNumBackNbr0.push_back(cp.numBackNbr[0][i]);
            }
            if (unChangePos1 == 0) {
                newOrder1.assign(cp.orders[1].begin() + 1, cp.orders[1].end());
                newNumBackNbr1 = cp.numBackNbr[1];
            }
            else {
                newOrder1 = simpleOrder(query, cs, newPrefix, newNumBackNbr1, cp.orders[1], unChangePos1);
                for (int i = unChangePos1 + 1; i < cp.numBackNbr[1].size(); ++i)
                    newNumBackNbr1.push_back(cp.numBackNbr[1][i]);
            }
            PrefixNode *newTree = new PrefixNode(99);
            PrefixNode *node = newTree;
            for (int i = 0; i < newPrefix.size(); ++i) {
                PrefixNode *newNode = new PrefixNode(newPrefix[i]);
                node -> children.push_back(newNode);
                node = node -> children[0];
            }
            CandidatePlan newPlan(newTree, {newOrder0, newOrder1}, {newNumBackNbr0, newNumBackNbr1}, nodePriority, cp.shareNum + 2, totalCost);
            plans.push(newPlan);
        }
        delete cp.pt;
    }
    cost = {bestCost0, bestCost1, bestCostShared};
    for (int i = 0; i < bestGlobal.size(); ++i) {
        t.nodes[0].attributes[i] = bestGlobal[i];
        t.nodes[1].attributes[i] = bestGlobal[i];
    }
    t.nodes[0].prefixSize = t.nodes[1].prefixSize = bestGlobal.size();
    t.nodes[0].prefix = new VertexID [bestGlobal.size()];
    t.nodes[1].prefix = new VertexID [bestGlobal.size()];
    for (int i = 0; i < bestGlobal.size(); ++i) {
        t.nodes[0].prefix[i] = bestGlobal[i];
        t.nodes[1].prefix[i] = bestGlobal[i];
    }
    for (int i = 0; i < bestOrder0.size(); ++i)
        t.nodes[0].attributes[i + bestGlobal.size()] = bestOrder0[i];
    for (int i = 0; i < bestOrder1.size(); ++i)
        t.nodes[1].attributes[i + bestGlobal.size()] = bestOrder1[i];
    t.defaultPartition = bestGlobal;
    t.globalOrder = bestGlobal;
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        for (int i = 0; i < t.nodes[nID].numAttributes; ++i) {
            VertexID u = t.nodes[nID].attributes[i];
            if (t.v2n[u].size() == 1) {
                t.globalOrder.push_back(u);
            }
        }
    }
    t.initPoses(query, cs, true);
}

// costs: the last one is the cost for the attribute tree; others are local costs for bags
void
optPrefixTreeCost(PrefixNode *root, uint64_t prevID, HyperTree &t, const std::vector<VertexID> &nIDs,
                  const Graph &query, const std::vector<SubsetStructure> &dpStructures, double &totalCost,
                  std::vector<double> &costs, const vector<vector<VertexID>> &matchedAttrs,
                  std::vector<std::vector<VertexID>> &nodeOrders, CandidateSpace &cs,
                  const std::vector<double> &factors) {
    std::vector<PrefixNode *> nodes(query.getNumVertices(), nullptr);
    std::vector<ui> childPoses(query.getNumVertices(), 0);
    uint64_t id = 0;
    std::vector<VertexID> attrsInPath, globalAttrsInPath;
    PrefixNode *pn = root;
    std::vector<uint64_t> matchedIDs(t.numNodes, 0);
    for (VertexID nID : nIDs) {
        for (int i = 0; i < matchedAttrs[nID].size(); ++i) {
            matchedIDs[nID] += 1 << matchedAttrs[nID][i];
        }
    }
    if (pn -> children.empty()) {
        costs.back() = 0;
        for (VertexID nID : nIDs) {
            dpStructures[nID].prefixPlan(matchedIDs[nID], matchedAttrs[nID].size(), nodeOrders[nID],
                                         costs[nID]);
        }
        totalCost = 0.0;
        for (double c: costs) totalCost += c;
        return;
    }
    for (VertexID nID: root -> nIDsToCall) {
        dpStructures[nID].prefixPlan(matchedIDs[nID], matchedAttrs[nID].size(), nodeOrders[nID], costs[nID]);
    }
    int depth = 0;
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            ++childPoses[depth];
            VertexID u = current -> u;
            if (id == 0) costs.back() += cs.candidateSet[u].size();
            else {
                double card = subsetToCard[id];
                double listSize = 0;
                std::vector<VertexID> vertexParents;
                for (VertexID u2: attrsInPath) {
                    if (query.getEdgeID(u, u2) != -1) vertexParents.push_back(u2);
                }
                if (!vertexParents.empty()) {
                    for (VertexID u2: vertexParents) {
                        uint64_t edgeID = (1 << u) + (1 << u2);
                        listSize += subsetToCard[edgeID] / subsetToCard[1 << u2];
                    }
                } else {
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
                costs.back() += card * listSize;
            }
            nodes[depth] = current;
            id += 1 << u;
            attrsInPath.push_back(u);
            for (VertexID nID: current -> nIDsToCall) {
                nodeOrders[nID].clear();
                for (VertexID u2: attrsInPath) nodeOrders[nID].push_back(u2);
                std::vector<VertexID> remainingOrder;
                dpStructures[nID].prefixPlan(id + matchedIDs[nID], depth + 1 + matchedAttrs[nID].size(),
                                             remainingOrder, costs[nID]);
                for (VertexID u2: remainingOrder) nodeOrders[nID].push_back(u2);
            }
            if (!current -> children.empty()) {
                ++depth;
                childPoses[depth] = 0;
            }
            else if (childPoses[depth] < pn -> children.size()) {
                id -= 1 << u;
                attrsInPath.pop_back();
            }
            if (depth > 0) pn = nodes[depth - 1];
        }
        --depth;
        if (depth >= 0) {
            if (depth == 0) pn = root;
            else pn = nodes[depth - 1];
            id -= 1 << attrsInPath.back();
            attrsInPath.pop_back();
            if (childPoses[depth] < pn -> children.size()) {
                id -= 1 << attrsInPath.back();
                attrsInPath.pop_back();
            }
        }
    }
    totalCost = 0.0;
    for (double c: costs) totalCost += c;
}

void
optSharePlan(const Graph &query, HyperTree &t, CandidateSpace &cs, PrefixNode *&pt, bool *visited, VertexID *partMatch,
             VertexID **candidates, ui *candCount, std::vector<double> &costs, std::vector<double> &noShareCosts,
             std::vector<HyperTree> &trees, std::vector<PrefixNode *> &prefixTrees, std::ofstream &outStream,
             const std::vector<SubsetStructure> &dpStructures) {
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
    auto start = std::chrono::steady_clock::now();
    ui numNodes = t.numNodes;
    if (t.newGlobalNode) --numNodes;
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsedSeconds = end - start;
    std::unordered_set<PrefixNode*, PrefixNodePtrHash, PrefixNodePtrEqual> visitedPT;
    start = std::chrono::steady_clock::now();
    genAllPrefixTree(t, query, cs, prefixTrees, false);
    end = std::chrono::steady_clock::now();
    elapsedSeconds = end - start;
    double enumTreeTime = elapsedSeconds.count();
    outStream << enumTreeTime << " " << std::flush;
    PrefixNode *bestT;
    std::vector<std::vector<VertexID>> bestOrders(t.numNodes);
    std::vector<ui> prefixSizes(t.numNodes, 0);
    double minCost = std::numeric_limits<double>::max();
    std::vector<double> minCosts(t.numNodes + 1, std::numeric_limits<double>::max());
//    std::vector<std::vector<std::vector<VertexID>>> currentOrders;
    for (int i = 0; i < prefixTrees.size(); ++i) {
        PrefixNode *prefixNode = prefixTrees[i];
        std::vector<double> currentCost(t.numNodes + 1, 0.0);
        std::vector<std::vector<VertexID>> currentOrder(t.numNodes);
        std::vector<ui> currentPrefixSize(t.numNodes, 0);
        std::vector<PrefixNode *> nodes(query.getNumVertices(), nullptr);
        int depth = 0;
        std::vector<ui> childPoses(query.getNumVertices(), 0);
        uint64_t id = 0;
        std::vector<VertexID> attrsInPath;
        PrefixNode *pn = prefixNode;
        if (pn -> children.empty()) {
            currentCost.back() = 0;
            for (VertexID nID = 0; nID < dpStructures.size(); ++nID) {
                currentOrder[nID] = dpStructures[nID].getOptOrder();
                currentCost[nID] = dpStructures[nID].getOptCost();
            }
//            currentOrders.push_back(currentOrder);
            double totalCost = 0.0;
            for (double c: currentCost) totalCost += c;
            if (totalCost < minCost) {
                minCost = totalCost;
                minCosts = currentCost;
                prefixSizes = currentPrefixSize;
                bestOrders = currentOrder;
                bestT = prefixTrees[i];
            }
            noShareCosts = currentCost;
            continue;
        }
        for (VertexID nID: prefixNode -> nIDsToCall) {
            currentPrefixSize[nID] = 0;
            currentOrder[nID] = dpStructures[nID].getOptOrder();
            currentCost[nID] = dpStructures[nID].getOptCost();
        }
        while (depth >= 0) {
            while (childPoses[depth] < pn -> children.size()) {
                PrefixNode *current = pn -> children[childPoses[depth]];
                ++childPoses[depth];
                VertexID u = current -> u;
                if (depth == 0) currentCost.back() += cs.candidateSet[u].size();
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
                    currentCost.back() += card * listSize;
                }
                nodes[depth] = current;
                id += 1 << u;
                attrsInPath.push_back(u);
                for (VertexID nID: current -> nIDsToCall) {
                    for (int j = 0; j <= depth; ++j) {
                        if (nodes[j]->pathToGlobal) ++currentPrefixSize[nID];
                        else break;
                    }
                    currentOrder[nID] = attrsInPath;
                    std::vector<VertexID> remainingOrder;
                    dpStructures[nID].prefixPlan(id, depth + 1, remainingOrder, currentCost[nID]);
                    for (VertexID u2: remainingOrder) currentOrder[nID].push_back(u2);
                }
                if (!current -> children.empty()) {
                    ++depth;
                    childPoses[depth] = 0;
                }
                else if (childPoses[depth] < pn -> children.size()) {
                    id -= 1 << u;
                    attrsInPath.pop_back();
                }
                if (depth > 0) pn = nodes[depth - 1];
            }
            --depth;
            if (depth >= 0) {
                if (depth == 0) pn = prefixNode;
                else pn = nodes[depth - 1];
                id -= 1 << attrsInPath.back();
                attrsInPath.pop_back();
                if (childPoses[depth] < pn -> children.size()) {
                    id -= 1 << attrsInPath.back();
                    attrsInPath.pop_back();
                }
            }
        }
        double totalCost = 0.0;
        for (double c: currentCost) totalCost += c;
//        currentOrders.push_back(currentOrder);
        if (totalCost < minCost) {
            minCost = totalCost;
            minCosts = currentCost;
            prefixSizes = currentPrefixSize;
            bestOrders = currentOrder;
            bestT = prefixTrees[i];
        }
    }
    if (t.newGlobalNode) {
        VertexID maxCostID;
        double maxCost = 0.0;
        for (VertexID nID2 = 0; nID2 < numNodes; ++nID2) {
            if (dpStructures[nID2].getOptCost() > maxCost) {
                maxCost = dpStructures[nID2].getOptCost();
                maxCostID = nID2;
            }
        }
        addGlobal(query, t, bestT, t.numNodes - 1, maxCostID, cs, bestOrders, prefixSizes, true);
        bestT->initNIDsToBuild(t.numNodes);
        bestT->initPoses(sharedAttrs, query, cs.dist);
    }

//    buildFromPrefixTree(prefixTrees, currentOrders, trees, t, sharedAttrs, query, cs);
    pt = bestT->clone();
    costs = minCosts;
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        delete[] t.nodes[nID].prefix;
        t.nodes[nID].prefixSize = prefixSizes[nID];
        t.nodes[nID].prefix = new VertexID [prefixSizes[nID]];
        for (int i = 0; i < prefixSizes[nID]; ++i) {
            t.nodes[nID].prefix[i] = bestOrders[nID][i];
        }
        for (int i = 0; i < t.nodes[nID].numAttributes; ++i) {
            t.nodes[nID].attributes[i] = bestOrders[nID][i];
        }
        t.nodes[nID].initPoses(sharedAttrs, query, cs.dist, nID == t.numNodes - 1);
    }
    t.globalOrder = bestOrders.back();
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        for (int i = 0; i < bestOrders[nID].size(); ++i) {
            VertexID u = bestOrders[nID][i];
            if (std::find(t.globalOrder.begin(), t.globalOrder.end(), u) == t.globalOrder.end())
                t.globalOrder.push_back(u);
        }
    }
    t.extendLevel = bestOrders.back().size();
    t.nIDs = std::vector<std::vector<VertexID>>(t.globalOrder.size());
    for (int j = 0; j < t.globalOrder.size(); ++j) {
        VertexID u = t.globalOrder[j];
        for (VertexID nID = 0; nID < bestOrders.size() - 1; ++nID) {
            if (std::find(bestOrders[nID].begin(), bestOrders[nID].end(), u) != bestOrders[nID].end())
                t.nIDs[j].push_back(nID);
        }
    }
    t.buildTrieOrder();
}

void
optInitPlan(const Graph &query, HyperTree &t, CandidateSpace &cs, PrefixNode *&pt, bool *visited, VertexID *partMatch,
            VertexID **candidates, ui *candCount, double &minCost, size_t budget, bool adaptive, bool materialize,
            bool optimizeGlobal, std::ofstream &outStream, std::vector<SubsetStructure> &dpStructures) {
    std::vector<HyperTree> trees;
    std::vector<PrefixNode *> prefixTrees;
    std::vector<double> noShareCosts;
    std::vector<double> costs;
    optSharePlan(query, t, cs, pt, visited, partMatch, candidates, candCount, costs,
                 noShareCosts, trees, prefixTrees, outStream, dpStructures);
    const HyperNode &last = t.nodes[t.numNodes - 1];
    setTDExtention(t, query);
    minCost = 0.0;
    for (auto c: costs) minCost += c;
    for (auto *prefixNode : prefixTrees) delete prefixNode;
}

void
simplePlan(const Graph &query, HyperTree &t, CandidateSpace &cs, PrefixNode *&pt, bool *visited, VertexID *partMatch,
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
    std::vector<SubsetStructure> dpStructures(numNodes);
    std::vector<std::vector<VertexID>> bestOrders;
    minCost = 0.0;
    for (VertexID nID = 0; nID < numNodes; ++nID) {
        std::vector<VertexID> localVertices(t.nodes[nID].attributes, t.nodes[nID].attributes + t.nodes[nID].numAttributes);
        SubsetStructure s(localVertices, query, false);
        s.optimalPlanDP(query, cs, visited, partMatch, candidates, candCount);
        dpStructures[nID] = s;
        minCost += s.getOptCost();
        bestOrders.push_back(s.getOptOrder());
    }
//    if (t.newGlobalNode) bestOrders.push_back(globalOrder(query, t, cs, std::vector<VertexID>()));
    if (t.newGlobalNode && bestPT != nullptr) {
        std::vector<VertexID> lastShared;
        std::vector<VertexID> bagOrder;
        std::vector<PrefixNode *> attributeOrder;
        bestPT->getTraverseOrder(attributeOrder, bagOrder, t);
        VertexID lastID = bagOrder[bagOrder.size() - 2];
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

void heuristicInitPlan(const Graph &query, HyperTree &t, CandidateSpace &cs, PrefixNode *&pt, bool *visited,
                       VertexID *partMatch, VertexID **candidates, ui *candCount, double &minCost, size_t budget,
                       bool adaptive, bool share, bool materialize, const std::vector<SubsetStructure> &dpStructures) {
    std::vector<std::vector<VertexID>> sharedAttrs(t.numNodes);
    std::vector<VertexID> globalAttr;
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        for (int i = 0; i < t.nodes[nID].numAttributes; ++i) {
            VertexID u = t.nodes[nID].attributes[i];
            if (t.v2n[u].size() > 1) {
                if (std::find(globalAttr.begin(), globalAttr.end(), u) == globalAttr.end())
                    globalAttr.push_back(u);
                sharedAttrs[nID].push_back(u);
            }
        }
        std::sort(sharedAttrs[nID].begin(), sharedAttrs[nID].end());
    }
    std::sort(globalAttr.begin(), globalAttr.end());
    ui numNodes = t.numNodes;
    if (t.newGlobalNode) --numNodes;
    std::vector<VertexID> smallNodes, largeNodes;
    std::vector<int> repetitions(query.getNumVertices(), 1);
    size_t currentMem = 0;
    std::vector<VertexID> partitionOrder;
    if (materialize) {
        for (VertexID nID = 0; nID < numNodes; ++nID) {
            if (nID == t.numNodes - 1 && t.newGlobalNode) continue;
            const HyperNode &tau = t.nodes[nID];
            std::vector<double> estimation;
            std::vector<VertexID> localOrder(tau.attributes, tau.attributes + tau.numAttributes);
            localOrder = simpleOrder(query, cs, localOrder, repetitions);
            std::vector<ui> poses(query.getNumVertices(), 0);
            size_t memUsage = 0;
            std::vector<std::vector<VertexID>> attributesBefore;
            std::vector<VertexID> cartesianParent;
            initPoses(localOrder, query, attributesBefore, cartesianParent, cs.dist);
            if (nID != t.numNodes - 1) {
                cardEstimateLabeled(localOrder, attributesBefore, cartesianParent, cs, visited, partMatch,
                                    candidates, candCount, poses, estimation);
                memUsage = estimation.back() * tau.numAttributes * sizeof(VertexID);
            }
            else {
                if (adaptive) {
                    cardEstimateLabeled(localOrder, attributesBefore, cartesianParent, cs, visited, partMatch,
                                        candidates, candCount, poses, estimation);
                    double totalCard = estimation.back();
                    attributesBefore.clear();
                    cartesianParent.clear();
                    initPoses(globalAttr, query, attributesBefore, cartesianParent, cs.dist);
                    cardEstimateLabeled(globalAttr, attributesBefore, cartesianParent, cs, visited, partMatch,
                                        candidates, candCount, poses, estimation);
                    memUsage = totalCard / estimation.back() * tau.numAttributes * sizeof(VertexID);
                }
                else {
                    std::vector<VertexID> materializedAttrs;
                    for (VertexID u : localOrder) {
                        if (t.v2n[u].size() == 1) materializedAttrs.push_back(u);
                    }
                    if (!materializedAttrs.empty())
                        materializedAttrs = simpleOrder(query, cs, materializedAttrs, repetitions);
                    attributesBefore.clear();
                    cartesianParent.clear();
                    initPoses(materializedAttrs, query, attributesBefore, cartesianParent, cs.dist);
                    cardEstimateLabeled(materializedAttrs, attributesBefore, cartesianParent, cs, visited, partMatch,
                                        candidates, candCount, poses, estimation);
                    memUsage = estimation.back() * tau.numAttributes * sizeof(VertexID);
                }

            }
            if (currentMem + memUsage > budget) largeNodes.push_back(nID);
            else {
                currentMem += memUsage;
                if (nID != t.numNodes - 1) smallNodes.push_back(nID);
            }
        }
        if (largeNodes.empty()) smallNodes.push_back(t.numNodes - 1);
        globalAttr.clear();
        std::vector<VertexID> localAttr;
        for (VertexID u = 0; u < query.getNumVertices(); ++u) {
            bool exists = false;
            const HyperNode &globalNode = t.nodes[t.numNodes - 1];
            for (int i = 0; i < globalNode.numAttributes; ++i) {
                if (u == globalNode.attributes[i]) {
                    exists = true;
                    globalAttr.push_back(u);
                    break;
                }
            }
            if (!exists) localAttr.push_back(u);
        }
        // an order for all attributes in large nodes
        for (VertexID nID : largeNodes) {
            for (int i = 0; i < t.nodes[nID].numAttributes; ++i) {
                VertexID u = t.nodes[nID].attributes[i];
                ++repetitions[u];
            }
        }
        globalAttr = simpleOrder(query, cs, globalAttr, repetitions);
        std::vector<ui> numBackWard;
        localAttr = simpleOrder(query, cs, globalAttr, numBackWard, localAttr, localAttr.size());
        partitionOrder = globalAttr;
        for (VertexID u: localAttr) partitionOrder.push_back(u);
    }
    else {
        for (VertexID nID = 0; nID < t.numNodes; ++nID) {
            if (t.newGlobalNode && nID == t.numNodes - 1) continue;
            smallNodes.push_back(nID);
        }
    }

    // the partition Order is projected to each large node
    std::vector<std::vector<VertexID>> orders(t.numNodes);
    std::vector<std::vector<VertexID>> posToNID(partitionOrder.size());
    std::vector<std::vector<double>> estimations(t.numNodes);
    for (VertexID nID : largeNodes) {
        for (int pos = 0; pos < partitionOrder.size(); ++pos) {
            VertexID u = partitionOrder[pos];
            bool exists = false;
            for (int i = 0; i < t.nodes[nID].numAttributes; ++i) {
                if (t.nodes[nID].attributes[i] == u) {
                    exists = true;
                    break;
                }
            }
            if (exists) {
                orders[nID].push_back(u);
                posToNID[pos].push_back(nID);
            }
        }
        std::vector<std::vector<VertexID>> vertexParents;
        std::vector<VertexID> cartesianParent;
        initPoses(orders[nID], query, vertexParents, cartesianParent, cs.dist);
        std::vector<ui> poses(query.getNumVertices(), 0);
        cardEstimateLabeled(orders[nID], vertexParents, cartesianParent, cs, visited, partMatch, candidates,
                            candCount, poses, estimations[nID]);
    }
    size_t totalMem = 0;
    std::vector<ui> lengths = std::vector<ui>(t.numNodes, 0);
    std::vector<std::vector<VertexID>> prefixes(t.numNodes);
    std::vector<size_t> currentMems(t.numNodes, 0);
    std::vector<bool> next(t.numNodes, true);
    std::vector<std::vector<VertexID>> nIDsToCall(partitionOrder.size());
    std::vector<VertexID> defaultPartition;
    if (!largeNodes.empty()) {
        for (VertexID nID : largeNodes) {
            currentMems[nID] = estimations[nID].back() * t.nodes[nID].numAttributes * sizeof(VertexID);
            totalMem += currentMems[nID];
        }
        for (int i = 0; i < partitionOrder.size(); ++i) {
            for (VertexID nID : posToNID[i]) {
                if (!next[nID]) continue;
                totalMem -= currentMems[nID];
                double prefixCard = estimations[nID][lengths[nID]];
                if (nID == t.numNodes && lengths[nID] < t.extendLevel) prefixCard = estimations[nID][t.extendLevel];
                ++lengths[nID];
                currentMems[nID] = estimations[nID].back() / prefixCard * (t.nodes[nID].numAttributes - lengths[nID]) * sizeof(VertexID);
                prefixes[nID].push_back(partitionOrder[i]);
                if (currentMems[nID] < budget / largeNodes.size()) {
                    next[nID] = false;
                    nIDsToCall[i].push_back(nID);
                }
                totalMem += currentMems[nID];
                if (totalMem < budget) break;
            }
            if (totalMem < budget) {
                budget -= totalMem;
                defaultPartition.assign(partitionOrder.begin(), partitionOrder.begin() + i + 1);
                break;
            }
        }
    }
    for (VertexID nID : largeNodes) {
        if (lengths[nID] == 0) {
            for (int i = 0; i < partitionOrder.size(); ++i) {
                if (std::find(posToNID[i].begin(), posToNID[i].end(), nID) != posToNID[i].end()) {
                    lengths[nID] = 1;
                    prefixes[nID].push_back(partitionOrder[i]);
                    nIDsToCall[i].push_back(nID);
                    if (i + 1 > defaultPartition.size()) defaultPartition.assign(partitionOrder.begin(), partitionOrder.begin() + i + 1);
                    break;
                }
            }
        }
    }
    if (!defaultPartition.empty()) {
        if (std::find(largeNodes.begin(), largeNodes.end(), t.numNodes - 1 ) == largeNodes.end())
            nIDsToCall[defaultPartition.size() - 1].push_back(t.numNodes - 1);
    }
    // the local orders for materialized nodes is built elsewhere
    // large nodes
    pt = new PrefixNode(99);
    PrefixNode *pn = pt;
    std::vector<std::vector<VertexID>> matchedAttrs(t.numNodes);
    std::vector<std::vector<VertexID>> localOrders(t.numNodes);
    std::vector<double> factors(t.numNodes, 1.0);
    for (int i = 0; i < defaultPartition.size(); ++i) {
        VertexID u = defaultPartition[i];
        PrefixNode *current = new PrefixNode(u);
        pn->children.push_back(current);
        pn = current;
    }
    pn = pt;
    std::vector<VertexID> prevAttr;
    PrefixNode *subTree = nullptr;
    minCost = 0.0;
    double subtreeCost = 0.0;
    std::vector<std::map<uint64_t, double>> prefixCost(t.numNodes);
    std::vector<std::map<uint64_t, std::vector<VertexID>>> prefixOrder(t.numNodes);
    subTreePlan(query, t, subTree, cs, smallNodes, prevAttr, matchedAttrs, localOrders, factors, visited,
                partMatch, candidates, candCount, subtreeCost, share, dpStructures);
    minCost += subtreeCost;
    if (subTree != nullptr) {
        std::vector<PrefixNode *> old = pn -> children;
        pn -> children.clear();
        for (PrefixNode *c : subTree->children) pn->children.push_back(c->clone());
        for (PrefixNode *c : old) pn->children.push_back(c);
        pn->nIDsToCall = subTree -> nIDsToCall;
        if (!largeNodes.empty()) pn->mergeToRight(localOrders, nIDsToCall[0]);
    }
    else {
        pn -> nIDsToCall = smallNodes;
        if (!largeNodes.empty()) pn->mergeToRight(localOrders, nIDsToCall[0]);
    }
    delete subTree;
    subTree = nullptr;
    uint64_t prevID = 0;
    std::vector<double> costs = computeCost(defaultPartition, query, cs, visited, partMatch, candidates, candCount);
    for (double c: costs) minCost += c;
    for (int i = 0; i < defaultPartition.size(); ++i) {
        if (!pn -> children.empty()) pn = pn -> children.back();
        VertexID u = defaultPartition[i];
        prevID += 1 << u;
        prevAttr.push_back(u);
        for (VertexID nID : nIDsToCall[i]) {
            for (int j = 0; j <= i; ++j) {
                VertexID u2 = defaultPartition[j];
                bool exists = false;
                for (int k = 0; k < t.nodes[nID].numAttributes; ++k) {
                    if (t.nodes[nID].attributes[k] == u2) {
                        exists = true;
                        break;
                    }
                }
                if (exists) matchedAttrs[nID].push_back(u2);
            }
            uint64_t matchID = 0;
            for (VertexID u2 : matchedAttrs[nID]) matchID += 1 << u2;
            factors[nID] = subsetToCard[prevID]  / subsetToCard[matchID];
        }
        subTreePlan(query, t, subTree, cs, nIDsToCall[i], prevAttr, matchedAttrs, localOrders, factors, visited,
                    partMatch, candidates, candCount, subtreeCost, share, dpStructures);
        minCost += subtreeCost;
        if (subTree != nullptr) {
            std::vector<PrefixNode *> old = pn -> children;
            pn -> children.clear();
            for (PrefixNode *c : subTree->children) pn->children.push_back(c->clone());
            for (PrefixNode *c : old) pn->children.push_back(c);
            pn->nIDsToCall = subTree -> nIDsToCall;
            if (i != defaultPartition.size() - 1) pn->mergeToRight(localOrders, nIDsToCall[i + 1]);
        }
        else {
            pn -> nIDsToCall = nIDsToCall[i];
            if (i != defaultPartition.size() - 1) pn->mergeToRight(localOrders, nIDsToCall[i + 1]);
        }
        delete subTree;
    }
    std::vector<std::vector<VertexID>> nodeOrders = matchedAttrs;
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        for (VertexID u: localOrders[nID])
            nodeOrders[nID].push_back(u);
    }
    if (t.newGlobalNode) {
        VertexID maxCostID;
        double maxCost = 0.0;
        for (VertexID nID2 = 0; nID2 < numNodes; ++nID2) {
            if (dpStructures[nID2].getOptCost() > maxCost) {
                maxCost = dpStructures[nID2].getOptCost();
                maxCostID = nID2;
            }
        }
        std::vector<ui> prefixSizes;
        addGlobal(query, t, pt, t.numNodes - 1, maxCostID, cs, nodeOrders, prefixSizes, true);
    }
    pt->initNIDsToBuild(t.numNodes);
    pt ->initPoses(sharedAttrs, query, cs.dist);
    buildFromPrefixTree(pt, nodeOrders, t, sharedAttrs, query, cs);
    t.defaultPartition = defaultPartition;
    // set the extend level
#ifdef COLLECT_RESULT
    t.extendLevel = t.defaultPartition.size();
    const HyperNode &last = t.nodes[t.numNodes - 1];
    for (int i = last.prefixSize; i < last.numAttributes; ++i) {
        ++t.extendLevel;
        t.defaultPartition.push_back(last.attributes[i]);
    }
#else
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
        t.buildTraverseStruct(query);
        if (materialize) {
            uint64_t prefixID = 0, allID = 0;
            for (int i = 0; i < t.extendLevel; ++i) prefixID += 1 << nodeOrders[nID][i];
            for (VertexID u : nodeOrders[nID]) allID += 1 << u;
            double allCard = subsetToCard[allID], prefixCard = subsetToCard[prefixID];
            while (sizeof(VertexID) * (t.nodes[nID].numAttributes - t.extendLevel) * allCard / prefixCard > budget) {
                prefixID += 1 << (nodeOrders[nID][t.extendLevel]);
                prefixCard = subsetToCard[prefixID];
                ++t.extendLevel;
                t.buildTraverseStruct(query);
            }
        }
    }
    else {
        t.extendLevel = t.nodes[t.numNodes - 1].numAttributes;
        for (int i = 0; i < t.extendLevel; ++i) {
            VertexID u = t.nodes[t.numNodes - 1].attributes[i];
            if (std::find(t.defaultPartition.begin(), t.defaultPartition.end(), u) == t.defaultPartition.end())
                t.defaultPartition.push_back(u);
        }
        t.buildTraverseStruct(query);
    }
#endif
}

void adjustSubTree(HyperTree &t, PrefixNode *&root, VertexID nID) {
    root -> nIDsToBuild.push_back(nID);
    root -> nIDsToCall.push_back(nID);
    std::sort(root -> nIDsToCall.begin(), root->nIDsToCall.end());
}

// the heuristic version
void
adjustSubTree(const Graph &query, HyperTree &t, PrefixNode *&root, VertexID nID, CandidateSpace &cs,
              const std::vector<VertexID> &prefix, const std::vector<int> &mappingSizes) {
    root -> nIDsToBuild.push_back(nID);
    VertexID u = root -> u;
    std::vector<VertexID> localVertices;
    for (int i = t.nodes[nID].prefixSize; i < t.nodes[nID].numAttributes; ++i) {
        if (t.nodes[nID].attributes[i] != u)
            localVertices.push_back(t.nodes[nID].attributes[i]);
    }
    std::vector<ui> numBackWard;
    std::vector<VertexID> nodeOrder = simpleOrder(query, cs, prefix, numBackWard, localVertices, localVertices.size());
    // check whether it can merge with existing paths in the subtree
    PrefixNode *pn = root;
    std::vector<VertexID> path = prefix;
    for (int i = 0; i < t.nodes[nID].numAttributes; ++i) {
        u = nodeOrder[i];
        for (PrefixNode *c : pn -> children) {
            if (c -> u == u) {
                path.push_back(u);
                pn = c;
            }
        }
    }
    PrefixNode *last = pn;
    bool share = false;
    for (VertexID nID2 : last -> nIDsToCall) {
        int length = 0;
        for (int i = 0; i < t.nodes[nID].numAttributes - path.size() && i < t.nodes[nID2].numAttributes - mappingSizes[nID2]; ++i) {
            VertexID u1 = nodeOrder[path.size() + i];
            VertexID u2 = t.nodes[nID2].attributes[i + mappingSizes[nID2]];
            if (u1 != u2) break;
            else ++length;
        }
        for (int i = 0; i < length; ++i) {
            VertexID u1 = nodeOrder[path.size() + i];
            PrefixNode *c = new PrefixNode(u1);
            if (nID2 == t.numNodes - 1) {
                c->pathToGlobal = true;
                pn -> children.push_back(c);
            }
            else {
                c->pathToGlobal = false;
                std::vector<PrefixNode *> old = pn -> children;
                pn -> children.clear();
                pn -> children.push_back(c);
                for (PrefixNode *c2 : old) pn -> children.push_back(c2);
            }
            pn = c;
        }
        if (length != 0) {
            share = true;
            pn -> nIDsToCall = {nID, nID2};
            for (auto it = last->nIDsToCall.begin(); it != last->nIDsToCall.end(); ++it) {
                if (*it == nID2) {
                    last->nIDsToCall.erase(it);
                    break;
                }
            }
            break;
        }
    }
    if (!share) {
        last->nIDsToCall.push_back(nID);
        std::sort(last->nIDsToCall.begin(), last->nIDsToCall.end());
    }
}

// the brute-force version
// set prefix sizes before calling it
void subTreePlan(const Graph &query, HyperTree &t, CandidateSpace &cs, const vector<VertexID> &nIDs,
                 const std::vector<SubsetStructure> &dpStructures, uint64_t prevID,
                 std::map<uint64_t, PrefixNode *> &bestPTs, std::map<uint64_t, double> &ptCosts,
                 const vector<vector<VertexID>> &matchedAttrs, std::map<uint64_t, std::vector<std::vector<VertexID>>> &nodeOrders,
                 std::vector<double> &factors) {
    // attributes shared with other nodes under the same root
    std::vector<std::vector<VertexID>> v2n(query.getNumVertices());
    for (VertexID nID : nIDs) {
        for (int i = 0; i < t.nodes[nID].numAttributes; ++i) {
            VertexID u = t.nodes[nID].attributes[i];
            if ((prevID & (1 << u)) != 0) continue;
            v2n[u].push_back(nID);
        }
    }
    std::vector<std::vector<VertexID>> sharedAttrs(t.numNodes);
    for (VertexID nID : nIDs) {
        for (int i = 0; i < t.nodes[nID].numAttributes; ++i) {
            VertexID u = t.nodes[nID].attributes[i];
            if (v2n[u].size() > 1) {
                sharedAttrs[nID].push_back(u);
            }
        }
        std::sort(sharedAttrs[nID].begin(), sharedAttrs[nID].end());
    }
    std::vector<std::vector<std::vector<VertexID>>> permutations(nIDs.size());
    for (VertexID i = 0; i < nIDs.size(); ++i) {
        VertexID nID = nIDs[i];
        do {
            if (orderConnectivity(query, sharedAttrs[nID], dpStructures[nID].cc))
                permutations[i].push_back(sharedAttrs[nID]);
        } while (std::next_permutation(sharedAttrs[nID].begin(), sharedAttrs[nID].end()));
    }
    for (int i = 0; i < nIDs.size(); ++i) {
        VertexID nID = nIDs[i];
        std::vector<VertexID> noShare = {99};
        for (auto u: sharedAttrs[nID]) {
            noShare.push_back(u);
        }
        permutations[i].push_back(noShare);
    }

    std::vector<ui> poses(t.numNodes, 0);
    std::unordered_set<PrefixNode*, PrefixNodePtrHash, PrefixNodePtrEqual> visitedPT;
    std::vector<std::vector<VertexID>> orders(t.numNodes);
    int depth = 0;
    while (depth >= 0) {
        while (poses[depth] < permutations[depth].size()) {
            orders[nIDs[depth]] = permutations[depth][poses[depth]];
            ++poses[depth];
            if (depth == nIDs.size() - 1) {
                bool exist = true;
                PrefixNode *tree = buildPrefixTree(orders, query, cs.dist, visitedPT, exist);
                if (!exist) {
                    visitedPT.insert(tree);
                    // compute the cost for this prefix tree
                    double totalCost;
                    std::vector<double> costs(t.numNodes);
                    std::vector<std::vector<VertexID>> localOrder(t.numNodes), fullOrder(t.numNodes);
                    optPrefixTreeCost(tree, prevID, t, nIDs, query, dpStructures, totalCost, costs, matchedAttrs,
                                      localOrder, cs, factors);
                    for (VertexID nID : nIDs) {
                        fullOrder[nID] = matchedAttrs[nID];
                        for (VertexID u: localOrder[nID]) fullOrder[nID].push_back(u);
                    }
                    // generate subsets of nIDs and update costs
                    std::sort(tree->nIDsToCall.begin(), tree->nIDsToCall.end());
                    std::vector<VertexID> sharedNIDs;
                    std::set_difference(nIDs.begin(), nIDs.end(), tree->nIDsToCall.begin(), tree->nIDsToCall.end(), std::back_inserter(sharedNIDs));
                    uint64_t sharedID = 0;
                    for (VertexID nID : sharedNIDs) sharedID += 1 << nID;
                    // considering shared bags
                    for (int i = 0; i <= tree->nIDsToCall.size(); ++i) {
                        if (tree->nIDsToCall.empty()) {
                            if (ptCosts.find(sharedID) == ptCosts.end() || totalCost < ptCosts[sharedID]) {
                                delete bestPTs[sharedID];
                                bestPTs[sharedID] = tree->clone();
                                ptCosts[sharedID] = totalCost;
                                nodeOrders[sharedID] = fullOrder;
                            }
                            continue;
                        }
                        std::vector<std::vector<bool>> choices = chooseK(tree->nIDsToCall.size(), i);
                        for (auto choice: choices) {
                            std::vector<VertexID> selected;
                            double cost = totalCost;
                            uint64_t bagID = sharedID;
                            std::vector<std::vector<VertexID>> nodeOrder = fullOrder;
                            for (int j = 0; j < tree->nIDsToCall.size(); ++j) {
                                VertexID nID = tree->nIDsToCall[j];
                                if (choice[j]) {
                                    selected.push_back(nID);
                                    bagID += 1 << nID;
                                }
                                else {
                                    nodeOrder[nID].clear();
                                    cost -= costs[nID];;
                                }
                            }
                            if (ptCosts.find(bagID) == ptCosts.end() || cost < ptCosts[bagID]) {
                                delete bestPTs[bagID];
                                bestPTs[bagID] = tree->clone();
                                bestPTs[bagID]->nIDsToCall = selected;
                                ptCosts[bagID] = cost;
                                nodeOrders[bagID] = nodeOrder;
                            }
                        }
                    }
                    // not considering shared bags
                    if (sharedID != 0) {
                        for (int i = 0; i <= tree->nIDsToCall.size(); ++i) {
                            if (tree->nIDsToCall.empty()) continue;
                            std::vector<std::vector<bool>> choices = chooseK(tree->nIDsToCall.size(), i);
                            for (auto choice: choices) {
                                std::vector<VertexID> selected;
                                double cost = 0.0;
                                uint64_t bagID = 0;
                                std::vector<std::vector<VertexID>> nodeOrder = fullOrder;
                                for (VertexID nID : sharedNIDs) nodeOrder[nID].clear();
                                for (int j = 0; j < tree->nIDsToCall.size(); ++j) {
                                    VertexID nID = tree->nIDsToCall[j];
                                    if (choice[j]) {
                                        selected.push_back(nID);
                                        cost += costs[nID];
                                        bagID += 1 << tree->nIDsToCall[j];
                                    }
                                    else nodeOrder[nID].clear();
                                }
                                if (ptCosts.find(bagID) == ptCosts.end() || cost < ptCosts[bagID]) {
                                    delete bestPTs[bagID];
                                    bestPTs[bagID] = new PrefixNode(99);
                                    bestPTs[bagID] -> nIDsToCall = selected;
                                    ptCosts[bagID] = cost;
                                    nodeOrders[bagID] = nodeOrder;
                                }
                            }
                        }
                    }
                }
            }
            else {
                ++depth;
                poses[depth] = 0;
            }
        }
        --depth;
    }

    for (auto it: visitedPT) delete it;
}

// the heuristic version
void
subTreePlan(const Graph &query, HyperTree &t, PrefixNode *&bestPT, CandidateSpace &cs, const vector<VertexID> &nIDs,
            const vector<VertexID> &prevAttrs, const vector<vector<VertexID>> &matchedAttrs,
            vector<vector<VertexID>> &localOrders, std::vector<double> &factors, bool *visited, VertexID *partMatch,
            VertexID **candidates, ui *candCount, double &minCost, bool share,
            const std::vector<SubsetStructure> &dpStructures) {
    std::vector<ui> poses(query.getNumVertices(), 0);
    std::vector<std::vector<VertexID>> v2n(query.getNumVertices());
    std::vector<std::vector<VertexID>> noShareOrders(t.numNodes);
    std::vector<std::vector<ui>> numBackWards(t.numNodes);
    for (VertexID nID : nIDs) {
        for (int i = 0; i < t.nodes[nID].numAttributes; ++i) {
            VertexID u = t.nodes[nID].attributes[i];
            if (std::find(prevAttrs.begin(), prevAttrs.end(), u) != prevAttrs.end()) continue;
            v2n[u].push_back(nID);
            noShareOrders[nID].push_back(u);
        }
    }
    std::vector<std::vector<VertexID>> sharedAttrs(t.numNodes);
    std::vector<VertexID> internalAttrs;
    std::vector<std::vector<VertexID>> attrIntersections(t.numNodes * t.numNodes);
    for (VertexID nID : nIDs) {
        for (int i = 0; i < t.nodes[nID].numAttributes; ++i) {
            VertexID u = t.nodes[nID].attributes[i];
            if (v2n[u].size() > 1) {
                sharedAttrs[nID].push_back(u);
            }
        }
        std::sort(sharedAttrs[nID].begin(), sharedAttrs[nID].end());
    }
    for (VertexID u = 0; u < query.getNumVertices(); ++u) {
        if (v2n[u].size() > 1)
            internalAttrs.push_back(u);
    }
    for (int i = 0; i < nIDs.size(); ++i) {
        VertexID nID1 = nIDs[i];
        for (int j = 0; j < i; ++j) {
            VertexID nID2 = nIDs[j];
            if (nID1 < nID2) std::swap(nID1, nID2);
            std::set_intersection(sharedAttrs[nID1].begin(), sharedAttrs[nID1].end(), sharedAttrs[nID2].begin(),
                                  sharedAttrs[nID2].end(), std::back_inserter(attrIntersections[nID2 * t.numNodes + nID1]));
        }
    }
    minCost = 0.0;
    std::vector<double> nodeCosts(t.numNodes, 0.0);
    for (VertexID nID : nIDs) {
        noShareOrders[nID] = dpStructures[nID].getOptOrder();
        nodeCosts[nID] = factors[nID] * dpStructures[nID].getOptCost();
        minCost += nodeCosts[nID];
        numBackWards[nID] = computeNumBackWard(query, matchedAttrs[nID], noShareOrders[nID]);
    }
    if (!share || nIDs.size() < 2) {
        for (VertexID nID: nIDs)
            localOrders[nID] = noShareOrders[nID];
        return;
    }
    // prioritize bags with larger costs
    std::vector<std::pair<double, ui>> costPairs;
    for (VertexID nID : nIDs) costPairs.emplace_back(nodeCosts[nID], nID);
    std::sort(costPairs.begin(), costPairs.end(), [](const auto &a, const auto &b) {
        return a.first > b.first;
    });
    std::vector<ui> nodePriority;
    for (auto &c : costPairs) nodePriority.push_back(c.second);
    std::vector<std::vector<VertexID>> components;
    std::vector<uint64_t> cc;
    query.computeConnectedComponents(internalAttrs, components);
    for (auto &c : components) cc.push_back(getSubsetID(c));
    std::priority_queue<CandidatePlan> plans;
    PrefixNode *empty = new PrefixNode(99);
    empty -> nIDsToCall = nIDs;
    std::vector<std::vector<VertexID>> bestOrders = noShareOrders;
    plans.emplace(empty, noShareOrders, numBackWards, nodePriority, 0, minCost);
    bestPT = empty -> clone();
    int numPlan = 0;
    std::unordered_set<PrefixNode*, PrefixNodePtrHash, PrefixNodePtrEqual> visitedPT;
    while (!plans.empty()) {
        CandidatePlan cp = plans.top();
        plans.pop();
        // evaluate the true cost of the candidate plan
//        double cost = computeCost(cp.pt, cp.orders, query, cs, visited, partMatch, candidates, candCount,
//                                  matchedAttrs, prevAttrs, factors);
        double cost = computeCostDP(cp.pt, query, cs, visited, partMatch, candidates, candCount, localOrders,
                                    dpStructures, factors, matchedAttrs, prevAttrs);
        if (cost < minCost) {
            minCost = cost;
            delete bestPT;
            bestPT = cp.pt -> clone();
//            bestOrders = cp.orders;
            bestOrders = localOrders;
        }
        ++numPlan;
        if (numPlan == MAX_NUM_PLAN) {
            while (!plans.empty()) {
                cp = plans.top();
                plans.pop();
            }
            break;
        }
        // two ways to extend a plan
        // 1. extend a prefix node by one attribute, select two bags to push down
        // 2. push down one bag in a non-leaf
        PrefixNode *pn = cp.pt;
        ui height = pn -> getHeight();
        std::vector<ui> childPoses(height, 0);
        std::vector<PrefixNode *> nodes(height, nullptr);
        int depth = 0;
        std::vector<VertexID> attrsInPath;
        std::vector<VertexID> siblingAttr;
        newPlanHeuristic(query, t, cs, plans, pn, attrsInPath, siblingAttr, -1, childPoses, cp, cost, visitedPT,
                         nodePriority, attrIntersections, sharedAttrs, cc);
        while (depth >= 0) {
            while (childPoses[depth] < pn -> children.size()) {
                PrefixNode *current = pn -> children[childPoses[depth]];
                VertexID u = current -> u;
                attrsInPath.push_back(u);
                nodes[depth] = current;
                siblingAttr.clear();
                if (depth == 0) {
                    for (PrefixNode *sibling: cp.pt -> children)
                        if (sibling != current) siblingAttr.push_back(sibling -> u);
                }
                else {
                    for (PrefixNode *sibling: nodes[depth - 1] -> children)
                        if (sibling != current) siblingAttr.push_back(sibling -> u);
                }
                ++childPoses[depth];
                newPlanHeuristic(query, t, cs, plans, current, attrsInPath, siblingAttr, depth, childPoses, cp, cost,
                                 visitedPT, nodePriority, attrIntersections, sharedAttrs, cc);
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
                if (depth == 0) pn = cp.pt;
                else pn = nodes[depth - 1];
                attrsInPath.pop_back();
                if (childPoses[depth] < pn -> children.size()) attrsInPath.pop_back();
            }
        }
    }
    int depth = 0;
    std::vector<VertexID> attrsInPath;
    PrefixNode *pn = bestPT;
    ui height = pn -> getHeight();
    std::vector<ui> childPoses(height, 0);
    std::vector<PrefixNode *> nodes(height, nullptr);
    for (VertexID nID : pn -> nIDsToCall) localOrders[nID] = bestOrders[nID];
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            ++childPoses[depth];
            VertexID u = current -> u;
            attrsInPath.push_back(u);
            nodes[depth] = current;
            for (VertexID nID: current -> nIDsToCall) {
                localOrders[nID] = attrsInPath;
                for (VertexID u2 : bestOrders[nID]) localOrders[nID].push_back(u2);
            }
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
            if (depth == 0) pn = bestPT;
            else pn = nodes[depth - 1];
            attrsInPath.pop_back();
            if (childPoses[depth] < pn -> children.size()) attrsInPath.pop_back();
        }
    }

    for (PrefixNode *pt : visitedPT) delete pt;
}

void setTDExtention(HyperTree &t, const Graph &query) {
#ifdef COLLECT_RESULT
    t.extendLevel = t.defaultPartition.size();
//    const HyperNode &last = t.nodes[t.numNodes - 1];
//    for (int i = last.prefixSize; i < last.numAttributes; ++i) {
//        ++t.extendLevel;
//        t.defaultPartition.push_back(last.attributes[i]);
//    }
#else
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
newPlanHeuristic(const Graph &query, const HyperTree &t, CandidateSpace &cs, std::priority_queue<CandidatePlan> &plans,
                 const PrefixNode *current, const std::vector<VertexID> &attrsInPath,
                 const std::vector<VertexID> &siblingAttr, int depth, std::vector<ui> &childPoses,
                 const CandidatePlan &cp, double parentCost,
                 std::unordered_set<PrefixNode *, PrefixNodePtrHash, PrefixNodePtrEqual> &visitedPT,
                 const std::vector<ui> &nodePriority, const std::vector<std::vector<VertexID>> &attrIntersections,
                 const std::vector<std::vector<VertexID>> &sharedAttrs, const std::vector<uint64_t> &cc) {
    VertexID u = current -> u;
    if (current -> nIDsToCall.size() >= 2) {
        std::vector<std::vector<bool>> choices = chooseK(current->nIDsToCall.size(), 2);
        for (auto choice: choices) {
            std::vector<VertexID> selected;
            for (int j = 0; j < current->nIDsToCall.size(); ++j) {
                VertexID nID = current->nIDsToCall[j];
                if (choice[j]) selected.push_back(nID);
            }
            ui pos = selected[0] * t.numNodes + selected[1];
            if (selected[0] > selected[1]) pos = selected[1] * t.numNodes + selected[0];
            std::vector<VertexID> attrIntersect = attrIntersections[pos];
            for (VertexID u2: attrIntersect) {
                if (std::find(siblingAttr.begin(), siblingAttr.end(), u2) != siblingAttr.end()) continue;
                if (std::find(attrsInPath.begin(), attrsInPath.end(), u2) != attrsInPath.end()) continue;
                std::vector<VertexID> prefix = attrsInPath;
                prefix.push_back(u2);
                if (!subsetConnectivity(query, cc, prefix)) continue;
                PrefixNode *newPT = cp.pt -> clone();
                PrefixNode *newPN = newPT;
                for (int i = 0; i <= depth; ++i) {
                    newPN = newPN -> children[childPoses[i] - 1];
                }
                PrefixNode *newChild = new PrefixNode(u2);
                newChild -> nIDsToCall = selected;
                if (selected[0] == t.numNodes - 1 || selected[1] == t.numNodes - 1) newPN->children.push_back(newChild);
                else {
                    std::vector<PrefixNode *> old = newPN -> children;
                    newPN -> children = {newChild};
                    for (PrefixNode * c: old) newPN -> children.push_back(c);
                }
                newPN -> nIDsToCall.clear();
                for (int i = 0; i < current->nIDsToCall.size(); ++i) {
                    if (!choice[i]) newPN -> nIDsToCall.push_back(current -> nIDsToCall[i]);
                }
                if (visitedPT.find(newPT) != visitedPT.end()) {
                    delete newPT;
                    continue;
                }
                else visitedPT.insert(newPT);
                std::vector<std::vector<VertexID>> newOrders = cp.orders, newNumBackWard = cp.numBackNbr;
                std::vector<ui> newNumBackNbr0(cp.numBackNbr[selected[0]].begin(), cp.numBackNbr[selected[0]].begin() + depth + 1);
                newNumBackNbr0.push_back(0);
                for (int i = 0; i <= depth; ++i) {
                    if (query.getEdgeID(u2, attrsInPath[i]) != -1)
                        ++newNumBackNbr0.back();
                }
                std::vector<ui> newNumBackNbr1 = newNumBackNbr0;
                int unChangePos0 = 0;
                while (unChangePos0 < cp.orders[selected[0]].size()) {
                    if (cp.orders[selected[0]][unChangePos0] == u2) break;
                    else ++unChangePos0;
                }
                int unChangePos1 = 0;
                while (unChangePos1 < cp.orders[selected[1]].size()) {
                    if (cp.orders[selected[1]][unChangePos1] == u2) break;
                    else ++unChangePos1;
                }
                std::vector<VertexID> newOrder0, newOrder1;
                if (unChangePos0 == 0) {
                    newOrder0.assign(cp.orders[selected[0]].begin() + 1, cp.orders[selected[0]].end());
                    newNumBackNbr0 = cp.numBackNbr[selected[0]];
                }
                else {
                    newOrder0 = simpleOrder(query, cs, prefix, newNumBackNbr0, cp.orders[selected[0]], unChangePos0);
                    for (int i = unChangePos0 + 1; i < cp.numBackNbr[selected[0]].size(); ++i)
                        newNumBackNbr0.push_back(cp.numBackNbr[selected[0]][i]);
                }
                if (unChangePos1 == 0) {
                    newOrder1.assign(cp.orders[selected[1]].begin() + 1, cp.orders[selected[1]].end());
                    newNumBackNbr1 = cp.numBackNbr[selected[1]];
                }
                else {
                    newOrder1 = simpleOrder(query, cs, prefix, newNumBackNbr1, cp.orders[selected[1]], unChangePos1);
                    for (int i = unChangePos1 + 1; i < cp.numBackNbr[selected[1]].size(); ++i)
                        newNumBackNbr1.push_back(cp.numBackNbr[selected[1]][i]);
                }
                newOrders[selected[0]] = newOrder0;
                newOrders[selected[1]] = newOrder1;
                newNumBackWard[selected[0]] = newNumBackNbr0;
                newNumBackWard[selected[1]] = newNumBackNbr1;
                plans.emplace(newPT, newOrders, newNumBackWard, nodePriority, cp.shareNum + 2, parentCost);
            }
        }
    }
    if (!current -> children.empty()) {
        for (int i = 0; i < current -> children.size(); ++i) {
            PrefixNode *c = current -> children[i];
            VertexID u2 = c -> u;
            for (VertexID nID : current -> nIDsToCall) {
                if (std::find(sharedAttrs[nID].begin(), sharedAttrs[nID].end(), u2) == sharedAttrs[nID].end()) continue;
                PrefixNode *newPT = cp.pt -> clone();
                PrefixNode *newPN = newPT;
                for (int j = 0; j <= depth; ++j) {
                    newPN = newPN -> children[childPoses[j] - 1];
                }
                newPN -> nIDsToCall.clear();
                for (VertexID nID2 : current -> nIDsToCall) {
                    if (nID != nID2) newPN -> nIDsToCall.push_back(nID2);
                }
                newPN -> children[i] -> nIDsToCall.push_back(nID);
                std::sort(newPN -> children[i] -> nIDsToCall.begin(), newPN->children[i]->nIDsToCall.end());
                if (visitedPT.find(newPT) != visitedPT.end()) {
                    delete newPT;
                    continue;
                }
                else visitedPT.insert(newPT);
                std::vector<VertexID> prefix = attrsInPath;
                prefix.push_back(u2);
                std::vector<std::vector<VertexID>> newOrders = cp.orders, newNumBackWard = cp.numBackNbr;
                std::vector<ui> newNumBackNbr(cp.numBackNbr[nID].begin(), cp.numBackNbr[nID].begin() + depth + 1);
                newNumBackNbr.push_back(0);
                std::vector<VertexID> newOrder;
                for (int j = 0; j <= depth; ++j) {
                    if (query.getEdgeID(u2, attrsInPath[j]) != -1)
                        ++newNumBackNbr.back();
                }
                int unChangePos = 0;
                while (unChangePos < cp.orders[nID].size()) {
                    if (cp.orders[nID][unChangePos] == u2) break;
                    else ++unChangePos;
                }
                if (unChangePos == 0) {
                    newOrder.assign(cp.orders[nID].begin() + 1, cp.orders[nID].end());
                    newNumBackNbr = cp.numBackNbr[nID];
                }
                else {
                    newOrder = simpleOrder(query, cs, prefix, newNumBackNbr, cp.orders[nID], unChangePos);
                    for (int j = unChangePos + 1; j < cp.numBackNbr[nID].size(); ++j)
                        newNumBackNbr.push_back(cp.numBackNbr[nID][j]);
                }
                newOrders[nID] = newOrder;
                newNumBackWard[nID] = newNumBackNbr;
                plans.emplace(newPT, newOrders, newNumBackWard, nodePriority, cp.shareNum + 1, parentCost);
            }
        }
    }
}

void
extendPrefixTree(int depth, const Graph &query, HyperTree &t, const std::vector<VertexID> &partitionOrder, VertexID nID,
                 const std::vector<ui> &lengths, PrefixNode *&newPrefix, HyperTree *&newTree, const PrefixNode *oldPT,
                 const std::vector<ui> &childPos, CandidateSpace &cs) {
    newPrefix = oldPT -> clone();
    newTree = new HyperTree();
    PrefixNode *pn = newPrefix;
    std::vector<PrefixNode *> nodes(depth + 1);
    nodes[0] = pn;
    PrefixNode *current = newPrefix;
    for (int i = 1; i < depth + 1; ++i) {
        current = current -> children[childPos[i]];
        if (current -> pathToGlobal) pn = current;
        nodes[i] = current;
    }

    // traverse the rightmost child to find the prefix node to put NID
    ui length = 0;
    VertexID u;
    std::vector<VertexID> prefix;
    std::vector<VertexID> newGlobalPrefix;
    int pos = 0;
    for (; pos < partitionOrder.size(); ++pos) {
        VertexID u2 = partitionOrder[pos];
        if (std::find(t.v2n[u2].begin(), t.v2n[u2].end(), nID) != t.v2n[u2].end()) {
            ++length;
            u = u2;
            prefix.push_back(u2);
            if (length == lengths[nID]) break;
        }
    }
    if (pos < t.defaultPartition.size() - 1) pos = t.defaultPartition.size() - 1;
    current = pn;
    while (current -> pathToGlobal && current -> u != u) {
        if (current -> children.empty()) break;
        if (!current->children.back()->pathToGlobal) break;
        current = current->children.back();
    }
    bool globalChange = true;
    if (current -> u == u) {
        adjustSubTree(t, current, nID);
        globalChange = false;
    }
    else {
        // extend the rightmost branch of the prefix tree
        // current must call the global join. move it to the new node
        for (auto it = current -> nIDsToBuild.begin(); it != current -> nIDsToBuild.end();) {
            if (*it == t.numNodes - 1) {
                current -> nIDsToBuild.erase(it);
                break;
            }
            else ++it;
        }
        int start = -1, end = 0;
        for (int i = 0; i < partitionOrder.size(); ++i) {
            if (partitionOrder[i] == current -> u) start = i;
            if (partitionOrder[i] == u) {
                end = i;
                break;
            }
        }
        removeNID(t.numNodes - 1, current);
        for (int i = start + 1; i <= end; ++i) {
            PrefixNode *newNode = new PrefixNode(partitionOrder[i]);
            current->children.push_back(newNode);
            current = newNode;
            if (t.v2n[newNode->u].size() > 1) newGlobalPrefix.push_back(newNode -> u);
        }
        current->nIDsToCall.push_back(nID);
        current->nIDsToCall.push_back(t.numNodes - 1);
        current->nIDsToBuild.push_back(nID);
    }
    // remove the nID in the original position
    if (depth != -1) removeNID(nID, nodes, depth, newPrefix);
    else {
        for (auto it = pn -> nIDsToBuild.begin(); it != pn -> nIDsToBuild.end();) {
            if (*it == nID) {
                pn -> nIDsToBuild.erase(it);
                break;
            }
            else ++it;
        }
        for (auto it = pn -> nIDsToCall.begin(); it != pn -> nIDsToCall.end();) {
            if (*it == nID) {
                pn -> nIDsToCall.erase(it);
                break;
            }
            else ++it;
        }
    }
    // build a new hypertree to store the orders
    // change the prefix and local order for nID and the global node.
    newTree->newGlobalNode = t.newGlobalNode;
    newTree->numNodes = t.numNodes;
    newTree->nodes = new HyperNode[newTree->numNodes];
    newTree->numAttributes = t.numAttributes;
    for (VertexID nID2 = 0; nID2 < t.numNodes; ++nID2) {
        if (nID2 != nID && nID2 != t.numNodes - 1) {
            t.nodes[nID2].copyTo(newTree->nodes[nID2]);
        }
        else if (nID2 == nID) {
            newTree->nodes[nID2].prefixSize = lengths[nID2];
            newTree->nodes[nID2].prefix = new VertexID [lengths[nID2]];
            newTree->nodes[nID2].numAttributes = t.nodes[nID2].numAttributes;
            newTree->nodes[nID2].attributes = new VertexID [newTree->nodes[nID2].numAttributes];
            for (int i = 0; i < prefix.size(); ++i) {
                newTree->nodes[nID2].attributes[i] = newTree->nodes[nID2].prefix[i] = prefix[i];
            }
            std::vector<VertexID> localVertices;
            for (int i = t.nodes[nID2].prefixSize; i < t.nodes[nID2].numAttributes; ++i) {
                if (t.nodes[nID2].attributes[i] != u)
                    localVertices.push_back(t.nodes[nID2].attributes[i]);
            }
            std::vector<ui> numBackWard;
            localVertices = simpleOrder(query, cs, prefix, numBackWard, localVertices, localVertices.size());
            for (int i = 0; i < localVertices.size(); ++i) {
                newTree->nodes[nID2].attributes[i + lengths[nID2]] = localVertices[i];
            }
        }
        else {
            if (!globalChange) {
                t.nodes[nID2].copyTo(newTree->nodes[nID2]);
            }
            else {
                HyperNode &globalNode = newTree->nodes[nID2];
                globalNode.prefixSize = t.nodes[nID2].prefixSize + newGlobalPrefix.size();
                globalNode.prefix = new VertexID [globalNode.prefixSize];
                globalNode.numAttributes = t.nodes[nID2].numAttributes;
                globalNode.attributes = new VertexID [globalNode.numAttributes];
                for (int i = 0; i < t.nodes[nID2].prefixSize; ++i) {
                    globalNode.attributes[i] = globalNode.prefix[i] = t.nodes[nID2].attributes[i];
                }
                for (int i = 0; i < newGlobalPrefix.size(); ++i) {
                    globalNode.attributes[i + t.nodes[nID2].prefixSize] = globalNode.prefix[i + t.nodes[nID2].prefixSize]
                            = newGlobalPrefix[i];
                }
                std::vector<VertexID> localVertices;
                for (int i = 0; i < t.nodes[nID2].numAttributes; ++i) {
                    bool exists = false;
                    for (int j = 0; j < globalNode.prefixSize; ++j) {
                        if (t.nodes[nID2].attributes[i] == globalNode.prefix[j]) {
                            exists = true;
                            break;
                        }
                    }
                    if (!exists) localVertices.push_back(t.nodes[nID2].attributes[i]);
                }
                for (int i = 0; i < localVertices.size(); ++i) {
                    globalNode.attributes[i + globalNode.prefixSize] = localVertices[i];
                }
            }
        }
        newTree->nodes[nID2].initPoses(query, cs.dist);
    }
    newTree->globalOrder = partitionOrder;
    newTree->defaultPartition.assign(partitionOrder.begin(), partitionOrder.begin() + pos + 1);
    newTree->initPoses(query, cs, false);
    setTDExtention(*newTree, query);
    std::vector<std::vector<VertexID>> bagAttrs(t.numNodes);
    for (VertexID nID2 = 0; nID2 < t.numNodes; ++nID2) {
        bagAttrs[nID2].assign(t.nodes[nID2].attributes, t.nodes[nID2].attributes + t.nodes[nID2].numAttributes);
        std::sort(bagAttrs[nID2].begin(), bagAttrs[nID2].end());
    }
    newPrefix->initPoses(bagAttrs, query, cs.dist);
}

double computeCostDP(const PrefixNode *pt, const Graph &query, CandidateSpace &cs, bool *visited, VertexID *partMatch,
                     VertexID **candidates, ui *candCount, std::vector<std::vector<VertexID>> &localOrders,
                     const std::vector<SubsetStructure> &dpStructures, const std::vector<double> &factors,
                     const std::vector<std::vector<VertexID>> &matchedAttrs, const std::vector<VertexID> &prevAttrs) {
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
    std::vector<uint64_t> ids(matchedAttrs.size());
    for (VertexID nID = 0; nID < matchedAttrs.size(); ++nID) ids[nID] = getSubsetID(matchedAttrs[nID]);
    for (VertexID nID : pt -> nIDsToCall) {
        VertexID id = ids[nID];
        double c;
        dpStructures[nID].prefixPlan(id, matchedAttrs[nID].size(), localOrders[nID], c);
        cost += c * factors[nID];
    }
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            ++childPoses[depth];
            VertexID u = current -> u;
            attrsInPath.push_back(u);
            nodes[depth] = current;
            for (VertexID nID: current -> nIDsToCall) {
                uint64_t id = ids[nID];
                for (VertexID u2 : attrsInPath) id += 1 << u2;
                double c;
                dpStructures[nID].prefixPlan(id, matchedAttrs[nID].size() + attrsInPath.size(),
                                             localOrders[nID], c);
                cost += c * factors[nID];
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

void twoBagDP(const Graph &query, HyperTree &t, PrefixNode *&bestPT, CandidateSpace &cs, const vector<VertexID> &nIDs,
              const vector<VertexID> &prevAttrs, vector<vector<VertexID>> &localOrders, double &minCost,
              const std::vector<SubsetStructure> &dpStructures, const vector<VertexID> &sharedAttrs) {
    std::vector<uint64_t> tmpSubsetIDs;
    const SubsetStructure &s0 = dpStructures[nIDs[0]], &s1 = dpStructures[nIDs[1]];
    std::vector<VertexID> newAttrs;
    for (VertexID u: sharedAttrs) {
        if (std::find(prevAttrs.begin(), prevAttrs.end(), u) == prevAttrs.end()) newAttrs.push_back(u);
    }
    std::sort(newAttrs.begin(), newAttrs.end());
    for (int size = 1; size <= newAttrs.size(); ++size) {
        std::vector<VertexID> current;
        generateSubsetsK(newAttrs, size, 0, current, tmpSubsetIDs, s0.cc, query, true);
    }
    uint64_t prevID = getSubsetID(prevAttrs);
    std::vector<uint64_t> subsetIDs;
    for (uint64_t subsetID: tmpSubsetIDs) {
        std::vector<VertexID> path = prevAttrs;
        const std::vector<VertexID> &newPath = getSubsetFromID(subsetID, query.getNumVertices());
        for (VertexID u: newPath) path.push_back(u);
        std::vector<std::vector<VertexID>> components;
        query.computeConnectedComponents(path, components);
        if (components.size() == 1) subsetIDs.push_back(subsetID + prevID);
    }
    std::vector<VertexID> bestGlobal, bestOrder0, bestOrder1;
    minCost = 0.0;
    double cost0, cost1;
    s0.prefixPlan(prevID, prevAttrs.size(), bestOrder0, cost0);
    s1.prefixPlan(prevID, prevAttrs.size(), bestOrder1, cost1);
    minCost = cost0 + cost1;
    for (uint64_t subsetID: subsetIDs) {
        ui k = getSubsetFromID(subsetID, query.getNumVertices()).size();
        std::vector<VertexID> order0, order1, sharedOrder;
        double costNode0, costNode1, sharedCost;
        s0.extendPrefixPlan(query, cs, prevID, prevAttrs.size(), subsetID, k, sharedOrder, sharedCost);
        s0.prefixPlan(subsetID, k, order0, costNode0);
        s1.prefixPlan(subsetID, k, order1, costNode1);
        double cost = costNode0 + costNode1 + sharedCost;
        if (cost < minCost) {
            minCost = costNode0 + costNode1 + sharedCost;
            bestGlobal = sharedOrder;
            bestOrder0 = order0;
            bestOrder1 = order1;
        }
    }
    bestPT = new PrefixNode(99);
    localOrders[nIDs[0]] = bestOrder0;
    localOrders[nIDs[1]] = bestOrder1;
    PrefixNode *pn = bestPT;
    for (int i = 0; i < bestGlobal.size(); ++i) {
        PrefixNode *child = new PrefixNode(bestGlobal[i]);
        pn->children.push_back(child);
        pn = child;
    }
    pn->nIDsToCall = nIDs;
}

void addOneBagDP(const Graph &query, HyperTree &t, PrefixNode *&pt, CandidateSpace &cs, VertexID nID,
                 vector<vector<VertexID>> &localOrders, double &minCost,
                 const std::vector<SubsetStructure> &dpStructures,
                 const std::vector<std::vector<VertexID>> &attrIntersections) {
    // traverse the attribute tree and try to put nID to each attribute
    const SubsetStructure &s = dpStructures[nID];
    std::vector<VertexID> newShareAttr;
    double newCost = dpStructures[nID].getOptCost();
    localOrders[nID] = dpStructures[nID].getOptOrder();
    PrefixNode *pn = pt;
    ui height = pn -> getHeight();
    std::vector<ui> childPoses(height, 0);
    std::vector<PrefixNode *> nodes(height, nullptr);
    std::vector<VertexID> bagsToMerge;
    std::vector<std::vector<VertexID>> prev;
    std::vector<PrefixNode *> originalAttrs;
    std::vector<PrefixNode *> subTrees;
    for (VertexID nID2 : pt -> nIDsToCall) {
        bagsToMerge.push_back(nID2);
        prev.emplace_back();
        originalAttrs.push_back(pt);
        subTrees.push_back(nullptr);
    }
    PrefixNode *bestAttr = pt;
    std::vector<VertexID> attrsInPath;
    int depth = 0;
    std::vector<std::vector<VertexID>> bestOrders = localOrders;
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            ++childPoses[depth];
            VertexID u = current -> u;
            attrsInPath.push_back(u);
            bool exists = false;
            for (int i = 0; i < t.nodes[nID].numAttributes; ++i) {
                if (t.nodes[nID].attributes[i] == u) {
                    exists = true;
                    break;
                }
            }
            if (!exists) {
                if (childPoses[depth] < pn -> children.size()) attrsInPath.pop_back();
                continue;
            }
            double cost;
            std::vector<VertexID> localOrder;
            s.prefixPlan(getSubsetID(attrsInPath), attrsInPath.size(), localOrder, cost);
            if (cost < newCost) {
                newCost = cost;
                bestOrders[nID] = localOrder;
                bestAttr = current;
            }
            nodes[depth] = current;
            // check whether bag nID can share with one existing bag
            for (VertexID nID2 : current -> nIDsToCall) {
                bagsToMerge.push_back(nID2);
                prev.push_back(attrsInPath);
                originalAttrs.push_back(current);
                subTrees.push_back(nullptr);
            }
            if (!current -> children.empty()) {
                ++depth;
                childPoses[depth] = 0;
            }
            else if (childPoses[depth] < pn -> children.size()) {
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
    int bestPos = -1;
    std::vector<std::vector<VertexID>> orders(t.numNodes);
    for (int i = 0; i < bagsToMerge.size(); ++i) {
        VertexID nID2 = bagsToMerge[i];
        std::vector<VertexID> prevAttr = prev[i];
        uint64_t id = getSubsetID(prevAttr);
        double prevCost = dpStructures[nID2].getRemainingCost(id, prevAttr.size());
        double twoCost;
        std::vector<VertexID> sharedAttrs;
        if (nID < nID2) sharedAttrs = attrIntersections[nID * t.numNodes + nID2];
        else sharedAttrs = attrIntersections[nID2 * t.numNodes + nID];
        twoBagDP(query, t, subTrees[i], cs, {nID2, nID}, prevAttr, orders, twoCost,
                 dpStructures, sharedAttrs);
        twoCost -= prevCost;
        if (twoCost < newCost + 0.1) {
            newCost = twoCost;
            bestPos = i;
            bestOrders = localOrders;
            bestOrders[nID] = orders[nID];
            bestOrders[nID2] = orders[nID2];
        }
    }
    localOrders = bestOrders;
    if (bestPos != -1) {
        pn = originalAttrs[bestPos];
        PrefixNode *subTree = subTrees[bestPos];
        if (subTree->getHeight() == 1) pn->nIDsToCall.push_back(nID);
        else {
            for (PrefixNode *c: subTree->children) pn->children.push_back(c->clone());
            std::vector<VertexID> newCall;
            const std::vector<VertexID> &oldCall = originalAttrs[bestPos]->nIDsToCall;
            for (VertexID nID2: oldCall) {
                if (nID2 != bagsToMerge[bestPos]) newCall.push_back(nID2);
            }
            originalAttrs[bestPos]->nIDsToCall = newCall;
        }

    }
    else {
        bestAttr -> nIDsToCall.push_back(nID);
    }
    for (PrefixNode *attribute : subTrees) delete attribute;
    minCost += newCost;
}

void bagDPHeuristic(const Graph &query, HyperTree &t, CandidateSpace &cs, PrefixNode *&pt,
                    const std::vector<SubsetStructure> &dpStructures, double &minCost) {
    ui numNodes = t.numNodes;
    if (t.newGlobalNode) --numNodes;
    std::vector<std::vector<uint64_t>> subsets(numNodes + 1);
    std::vector<std::vector<std::vector<int>>> subsetOf(numNodes + 1);
    std::vector<std::vector<std::vector<int>>> supersetOf(numNodes + 1);
    std::vector<std::vector<std::vector<std::vector<VertexID>>>> localOrders(numNodes + 1);
    std::vector<std::vector<PrefixNode *>> attributeTrees(numNodes + 1);
    std::vector<std::vector<double>> subsetToCost(numNodes + 1);
    std::vector<VertexID> elements;
    for (VertexID nID = 0; nID < numNodes; ++nID) elements.push_back(nID);
    std::vector<uint64_t> cc;
    for (int k = 1; k <= numNodes; ++k) {
        std::vector<VertexID> current;
        generateSubsetsK(elements, k, 0, current, subsets[k], cc, query, false);
        std::sort(subsets[k].begin(), subsets[k].end());
    }
    subsetSupersetRelationships(numNodes, elements, subsets, subsetOf, supersetOf);
    for (int k = 0; k <= numNodes; ++k) {
        localOrders[k].resize(subsets[k].size());
        for (int i = 0;i < subsets[k].size(); ++i) localOrders[k][i].resize(t.numNodes);
        attributeTrees[k].resize(subsets[k].size());
        subsetToCost[k] = std::vector<double>(subsets[k].size());
    }
    std::vector<std::vector<VertexID>> sharedAttrs(t.numNodes);
    std::vector<std::vector<VertexID>> attrIntersections(t.numNodes * t.numNodes);
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        for (int i = 0; i < t.nodes[nID].numAttributes; ++i) {
            VertexID u = t.nodes[nID].attributes[i];
            if (t.v2n[u].size() > 1) {
                sharedAttrs[nID].push_back(u);
            }
        }
        std::sort(sharedAttrs[nID].begin(), sharedAttrs[nID].end());
    }
    for (VertexID nID1 = 0; nID1 < numNodes; ++nID1) {
        for (VertexID nID2 = nID1 + 1; nID2 < numNodes; ++nID2) {
            std::set_intersection(sharedAttrs[nID1].begin(), sharedAttrs[nID1].end(), sharedAttrs[nID2].begin(),
                                  sharedAttrs[nID2].end(), std::back_inserter(attrIntersections[nID1 * t.numNodes + nID2]));
        }
    }
    std::vector<VertexID> prevAttrs;
    // starting from two bags
    for (int i = 0; i < subsets[2].size(); ++i) {
        uint64_t id = subsets[2][i];
        std::vector<VertexID> nIDs = getSubsetFromID(id, numNodes);
        int pos = nIDs[0] * t.numNodes + nIDs[1];
        twoBagDP(query, t, attributeTrees[2][i], cs, nIDs, prevAttrs, localOrders[2][i],
                 subsetToCost[2][i], dpStructures, attrIntersections[pos]);
    }
    for (int level = 3; level < numNodes + 1; ++level) {
        for (int i = 0; i < subsets[level].size(); ++i) {
            uint64_t id = subsets[level][i];
            int minPos = 0;
            double minC = std::numeric_limits<double>::max();
            for (int j = 0; j < subsetOf[level][i].size(); ++j) {
                int pos = subsetOf[level][i][j];
                uint64_t childID = subsets[level - 1][pos];
                VertexID newBag = findExtraElement(id, childID);
                PrefixNode *old = attributeTrees[level - 1][pos];
                double cost = subsetToCost[level - 1][pos];
                PrefixNode *newPT = old->clone();
                std::vector<std::vector<VertexID>> localOrder = localOrders[level - 1][pos];
                addOneBagDP(query, t, newPT, cs, newBag, localOrder, cost, dpStructures, attrIntersections);
                if (cost < minC) {
                    minC = cost;
                    subsetToCost[level][i] = cost;
                    localOrders[level][i] = localOrder;
                    attributeTrees[level][i] = newPT;
                }
                else delete newPT;
            }
        }
    }
    pt = attributeTrees[numNodes][0];
    std::vector<std::vector<VertexID>> nodeOrders(t.numNodes);
    std::vector<VertexID> attrsInPath;
    int depth = 0;
    const PrefixNode *pn = pt;
    ui height = pn -> getHeight();
    std::vector<ui> childPoses(height, 0);
    std::vector<PrefixNode *> nodes(height, nullptr);
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            ++childPoses[depth];
            VertexID u = current -> u;
            attrsInPath.push_back(u);
            nodes[depth] = current;
            for (VertexID nID: current -> nIDsToCall) {
                nodeOrders[nID] = attrsInPath;
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
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        for (VertexID u: localOrders[numNodes][0][nID])
            nodeOrders[nID].push_back(u);
    }
    if (t.newGlobalNode) {
        VertexID maxCostID;
        double maxCost = 0.0;
        for (VertexID nID2 = 0; nID2 < numNodes; ++nID2) {
            if (dpStructures[nID2].getOptCost() > maxCost) {
                maxCost = dpStructures[nID2].getOptCost();
                maxCostID = nID2;
            }
        }
        std::vector<ui> prefixSizes;
        addGlobal(query, t, pt, t.numNodes - 1, maxCostID, cs, nodeOrders, prefixSizes, true);
    }
    pt->initNIDsToBuild(t.numNodes);
    pt ->initPoses(sharedAttrs, query, cs.dist);
    buildFromPrefixTree(pt, nodeOrders, t, sharedAttrs, query, cs);
    setTDExtention(t, query);
    minCost = subsetToCost[numNodes][0];
    for (int level = 2; level < numNodes; ++level) {
        for (PrefixNode *prefixNode : attributeTrees[level])
            delete prefixNode;
    }
}

void
bagSetDP(const Graph &query, HyperTree &t, CandidateSpace &cs, PrefixNode *&pt, vector<SubsetStructure> &dpStructures,
         double &bestCost, bool connected) {
    for (SubsetStructure &s: dpStructures) s.buildExtendCost(query, cs);
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
    bagSetStructure bs;
    bs.init(sharedAttrs, dpStructures, query, connected);
    // intersections of bags
    std::vector<VertexID> bagIDs;
    for (VertexID nID = 0; nID < sharedAttrs.size(); ++nID) bagIDs.push_back(nID);
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
                    double prefixCost = dpStructures[bags[0]].getOptCost(prefixID, k2);
                    std::vector<uint64_t> extensions;
                    if (k != sharedAttrs.size()) generateExtendedSubsets(maxPrefix, prefixID, extensions);
                    else extensions = {prefixID};
                    std::vector<std::vector<VertexID>> localOrder(numNodes);
                    for (uint64_t sharedID: extensions) {
                        if (connected && !subsetConnectivity(query, dpStructures[bags[0]].cc, getSubsetFromID(sharedID, query.getNumVertices()))) continue;
                        double minCost = std::numeric_limits<double>::max();
                        uint64_t minID1, minID2;
                        std::vector<VertexID> extendOrder = dpStructures[bags[0]].extendOrder[prefixID][sharedID];
                        double sharedCost = dpStructures[bags[0]].extendCost[prefixID][sharedID];
                        if (bs.indepCosts.find(bagSet) == bs.indepCosts.end() ||
                            bs.indepCosts[bagSet].find(prefixID) == bs.indepCosts[bagSet].end())
                            bs.indepCosts[bagSet][prefixID] = minCost;
                        for (uint64_t subset1: subsets1) {
                            uint64_t subset2 = bagSet - subset1;
                            if (subset1 == 0 || subset2 == 0 || subset1 > subset2) continue;
                            double cost = bs.indepCosts[subset1][sharedID] + bs.indepCosts[subset2][sharedID];
                            if (cost < minCost) {
                                minCost = cost;
                                minID1 = subset1;
                                minID2 = subset2;
                            }
                        }
                        minCost = minCost + sharedCost - prefixCost;
                        if (minCost < bs.indepCosts[bagSet][prefixID]) {
                            std::vector<VertexID> bags1 = getSubsetFromID(minID1, numNodes);
                            std::vector<VertexID> bags2 = getSubsetFromID(minID2, numNodes);
                            const std::vector<std::vector<VertexID>> &localOrder1 = bs.localOrders[minID1][sharedID];
                            const std::vector<std::vector<VertexID>> &localOrder2 = bs.localOrders[minID2][sharedID];
                            for (VertexID nID: bags1) {
                                localOrder[nID] = extendOrder;
                                for (VertexID u : localOrder1[nID]) localOrder[nID].push_back(u);
                            }
                            for (VertexID nID: bags2) {
                                localOrder[nID] = extendOrder;
                                for (VertexID u : localOrder2[nID]) localOrder[nID].push_back(u);
                            }
                            bs.indepCosts[bagSet][prefixID] = minCost;
                            bs.localOrders[bagSet][prefixID] = localOrder;
                        }
                    }
                }
            }
        }
    }
    uint64_t bestPrefix = 0;
    bestCost = std::numeric_limits<double>::max();
    uint64_t allBag = getSubsetID(bagIDs);
    for (auto it = bs.indepCosts[allBag].begin(); it != bs.indepCosts[allBag].end(); ++it) {
        uint64_t prefixID = it->first;
        double prefixCost = dpStructures[0].getOptCost(prefixID, subsetSize(prefixID, query.getNumVertices()));
        double cost = it->second + prefixCost;
        if (cost < bestCost) {
            bestCost = cost;
            bestPrefix = it->first;
        }
    }
    std::vector<std::vector<VertexID>> bestOrders = bs.localOrders[allBag][bestPrefix];
    if (bestPrefix != 0) {
        std::vector<VertexID> prefixOrder = dpStructures[0].getOptOrder(bestPrefix, subsetSize(bestPrefix, query.getNumVertices()));
        std::vector<std::vector<VertexID>> old = bestOrders;
        for (VertexID nID = 0; nID < dpStructures.size(); ++nID) {
            bestOrders[nID] = prefixOrder;
            for (VertexID u : old[nID]) bestOrders[nID].push_back(u);
        }
    }
    std::unordered_set<PrefixNode *, PrefixNodePtrHash, PrefixNodePtrEqual> visitedPT;
    bool exist;
    pt = buildPrefixTree(bestOrders, query, cs.dist, visitedPT, exist);
    if (t.newGlobalNode) {
        bestOrders.emplace_back();
        VertexID maxCostID;
        double maxCost = 0.0;
        for (VertexID nID2 = 0; nID2 < numNodes; ++nID2) {
            if (dpStructures[nID2].getOptCost() > maxCost) {
                maxCost = dpStructures[nID2].getOptCost();
                maxCostID = nID2;
            }
        }
        std::vector<ui> prefixSizes;
        addGlobal(query, t, pt, t.numNodes - 1, maxCostID, cs, bestOrders, prefixSizes, true);
        sharedAttrs.push_back(bestOrders.back());
        std::sort(sharedAttrs.back().begin(), sharedAttrs.back().end());
    }
    pt->initNIDsToBuild(t.numNodes);
    pt->initPoses(sharedAttrs, query, cs.dist);
    buildFromPrefixTree(pt, bestOrders, t, sharedAttrs, query, cs);
    setTDExtention(t, query);
}

void bagPermutationDP(const Graph &query, HyperTree &t, CandidateSpace &cs, PrefixNode *&pt,
                      vector<SubsetStructure> &dpStructures, double &bestCost, int reorder, bool connected,
                      bool globalOrderShare) {
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
        addGlobal(query, t, pt, t.numNodes - 1, maxCostID, cs, bestOrders, prefixSizes, globalOrderShare);
        sharedAttrs.push_back(bestOrders.back());
        std::sort(sharedAttrs.back().begin(), sharedAttrs.back().end());
    }
    pt->initNIDsToBuild(t.numNodes);
    pt->initPoses(sharedAttrs, query, cs.dist);
    buildFromPrefixTree(pt, bestOrders, t, sharedAttrs, query, cs);
    setTDExtention(t, query);
}

void maxSharePermutationDP(const Graph &query, HyperTree &t, CandidateSpace &cs, PrefixNode *&pt,
                           vector<SubsetStructure> &dpStructures, double &bestCost, bool connected) {
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
    bs.init(sharedAttrs, dpStructures, query, true, connected);
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
                        costs = std::vector<double>(prefix.size(), 1.0);
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
    if (t.newGlobalNode) {
        bestOrders.emplace_back();
        VertexID maxCostID;
        double maxCost = 0.0;
        for (VertexID nID2 = 0; nID2 < numNodes; ++nID2) {
            if (dpStructures[nID2].getOptCost() > maxCost) {
                maxCost = dpStructures[nID2].getOptCost();
                maxCostID = nID2;
            }
        }
        std::vector<ui> prefixSizes;
        addGlobal(query, t, pt, t.numNodes - 1, maxCostID, cs, bestOrders, prefixSizes, true);
        sharedAttrs.push_back(bestOrders.back());
        std::sort(sharedAttrs.back().begin(), sharedAttrs.back().end());
    }
    pt->initNIDsToBuild(t.numNodes);
    pt->initPoses(sharedAttrs, query, cs.dist);
    // the real cost of pt
    bestCost = 0.0;
    int depth = 0;
    PrefixNode *pn = pt;
    ui height = pn -> getHeight();
    std::vector<ui> childPoses = std::vector<ui>(height, 0);
    std::vector<PrefixNode *> path;
    std::vector<VertexID> attrs;
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            path.push_back(current);
            attrs.push_back(current->u);
            bestCost += computeCost(attrs, query, cs).back();
            for (VertexID nID: current->nIDsToCall) {
                if (t.newGlobalNode && nID == t.numNodes - 1) continue;
                uint64_t id = getSubsetID(attrs);
                bestCost += dpStructures[nID].getRemainingCost(id, attrs.size());
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
    buildFromPrefixTree(pt, bestOrders, t, sharedAttrs, query, cs);
    setTDExtention(t, query);
}

// heuristic: 1: appear in later bags; 2: cardinality smaller; 3: cost larger, 4: do not reorder
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