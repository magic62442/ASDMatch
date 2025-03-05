//
// Created by anonymous authors on 2024/6/29.
//

#include "adaptive.h"
double gGlobalTime = 0.0;
int gNumFixedCase = 0;
int gMaxNewFixed = 0;

void
initialPlan(const Graph &query, HyperTree &t, CandidateSpace &cs, PrefixNode *&pt, bool *visited, VertexID *partMatch,
            VertexID **candidates, ui *candCount, std::vector<double> &costs, std::vector<VertexID> &largeNodes,
            size_t budget) {
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
    bool globalBag = true;
    std::vector<VertexID> globalCand;
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        if (std::includes(sharedAttrs[nID].begin(), sharedAttrs[nID].end(), globalAttr.begin(), globalAttr.end())) {
            globalBag = false;
            globalCand.push_back(nID);
        }
    }
    if (!globalBag) {
        VertexID nID = globalCand[0];
        ui maxNum = t.nodes[nID].numAttributes;
        for (int i = 1; i < globalCand.size(); ++i) {
            ui num = t.nodes[globalCand[i]].numAttributes;
            if (num > maxNum) {
                maxNum = num;
                nID = globalCand[i];
            }
        }
        if (nID != t.numNodes - 1) {
            HyperNode tmp;
            tmp.numAttributes = t.nodes[nID].numAttributes;
            tmp.attributes = new VertexID[tmp.numAttributes];
            memcpy(tmp.attributes, t.nodes[nID].attributes, sizeof(VertexID) * tmp.numAttributes);
            t.nodes[nID].numAttributes = t.nodes[t.numNodes - 1].numAttributes;
            delete[] t.nodes[nID].attributes;
            t.nodes[nID].attributes = new VertexID [t.nodes[nID].numAttributes];
            memcpy(t.nodes[nID].attributes, t.nodes[t.numNodes - 1].attributes, sizeof(VertexID) * t.nodes[t.numNodes - 1].numAttributes);
            t.nodes[t.numNodes - 1].numAttributes = tmp.numAttributes;
            delete[] t.nodes[t.numNodes - 1].attributes;
            t.nodes[t.numNodes - 1].attributes = new VertexID [t.nodes[t.numNodes - 1].numAttributes];
            memcpy(t.nodes[t.numNodes - 1].attributes, tmp.attributes, sizeof(VertexID) * tmp.numAttributes);
            std::swap(sharedAttrs[nID], sharedAttrs[t.numNodes - 1]);
        }
    }
    else {
        sharedAttrs.push_back(globalAttr);
        t.addGlobalNode(globalAttr);
    }
    delete[] t.v2n;
    t.v2n = new std::vector<VertexID> [query.getNumVertices()];
    for (ui i = 0; i < t.numNodes; ++i) {
        for (ui j = 0; j < t.nodes[i].numAttributes; ++j) {
            t.v2n[t.nodes[i].attributes[j]].push_back(i);
        }
    }
    std::vector<int> repetitions(query.getNumVertices(), 0);
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
        bool global = nID == t.numNodes - 1;
        t.nodes[nID].initPoses(sharedAttrs, query, cs.dist, global);
    }
    std::vector<VertexID> vertices(query.getNumVertices());
    for (VertexID u = 0; u < query.getNumVertices(); ++u) {
        repetitions[u] += t.v2n[u].size();
        vertices[u] = u;
    }
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
    globalAttr = simpleOrder(query, cs, globalAttr, repetitions);
    std::vector<ui> numBackWard;
    localAttr = simpleOrder(query, cs, globalAttr, numBackWard, localAttr, localAttr.size());
    t.globalOrder = globalAttr;
    for (VertexID u : localAttr)
        t.globalOrder.push_back(u);
    std::vector<VertexID> materializeNodes;
    size_t currentMem = 0;
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        if (nID == t.numNodes - 1 && t.newGlobalNode) continue;
        const HyperNode &tau = t.nodes[nID];
        std::vector<double> estimation;
        std::vector<VertexID> localOrder(tau.attributes, tau.attributes + tau.numAttributes);
        std::vector<ui> poses(query.getNumVertices(), 0);
        if (nID != t.numNodes - 1) {
            cardEstimateLabeled(localOrder, tau.attributesBefore, tau.cartesianParent, cs, visited, partMatch,
                                candidates, candCount, poses, estimation);
        }
        else {
            std::vector<VertexID> materializedAttrs;
            materializedAttrs.assign(localOrder.begin() + globalAttr.size(), localOrder.end());
            cardEstimateLabeled(materializedAttrs, tau.attributesBefore, tau.cartesianParent, cs, visited, partMatch,
                                candidates, candCount, poses, estimation);
        }
        size_t memUsage = estimation.back() * tau.numAttributes * sizeof(VertexID);
        if (currentMem + memUsage > budget) largeNodes.push_back(nID);
        else {
            currentMem += memUsage;
            if (nID != t.numNodes - 1) materializeNodes.push_back(nID);
        }
    }

    // delay the large node plan after materializing small nodes
    pt = new PrefixNode(99);
    pt->nIDsToCall = materializeNodes;
    pt->nIDsToBuild = materializeNodes;
}

void
safePlan(const Graph &query, HyperTree &t, CandidateSpace &cs, PrefixNode *&pt, bool *visited, VertexID *partMatch,
         VertexID **candidates, ui *candCount, std::vector<double> &costs, size_t budget, bool share) {
    delete[] t.v2n;
    t.v2n = new std::vector<VertexID> [query.getNumVertices()];
    for (ui i = 0; i < t.numNodes; ++i) {
        for (ui j = 0; j < t.nodes[i].numAttributes; ++j) {
            t.v2n[t.nodes[i].attributes[j]].push_back(i);
        }
    }
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
    bool globalBag = true;
    std::vector<VertexID> globalCand;
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        if (std::includes(sharedAttrs[nID].begin(), sharedAttrs[nID].end(), globalAttr.begin(), globalAttr.end())) {
            globalBag = false;
            globalCand.push_back(nID);
        }
    }
    if (!globalBag) {
        VertexID nID = globalCand[0];
        ui maxNum = t.nodes[nID].numAttributes;
        for (int i = 1; i < globalCand.size(); ++i) {
            ui num = t.nodes[globalCand[i]].numAttributes;
            if (num > maxNum) {
                maxNum = num;
                nID = globalCand[i];
            }
        }
        if (nID != t.numNodes - 1) {
            HyperNode tmp;
            tmp.numAttributes = t.nodes[nID].numAttributes;
            tmp.attributes = new VertexID[tmp.numAttributes];
            memcpy(tmp.attributes, t.nodes[nID].attributes, sizeof(VertexID) * tmp.numAttributes);
            t.nodes[nID].numAttributes = t.nodes[t.numNodes - 1].numAttributes;
            delete[] t.nodes[nID].attributes;
            t.nodes[nID].attributes = new VertexID [t.nodes[nID].numAttributes];
            memcpy(t.nodes[nID].attributes, t.nodes[t.numNodes - 1].attributes, sizeof(VertexID) * t.nodes[t.numNodes - 1].numAttributes);
            t.nodes[t.numNodes - 1].numAttributes = tmp.numAttributes;
            delete[] t.nodes[t.numNodes - 1].attributes;
            t.nodes[t.numNodes - 1].attributes = new VertexID [t.nodes[t.numNodes - 1].numAttributes];
            memcpy(t.nodes[t.numNodes - 1].attributes, tmp.attributes, sizeof(VertexID) * tmp.numAttributes);
            std::swap(sharedAttrs[nID], sharedAttrs[t.numNodes - 1]);
        }
    }
    else {
        sharedAttrs.push_back(globalAttr);
        t.addGlobalNode(globalAttr);
    }

    std::vector<int> repetitions(query.getNumVertices(), 0);
    std::vector<VertexID> vertices(query.getNumVertices());
    for (VertexID u = 0; u < query.getNumVertices(); ++u) {
        repetitions[u] += t.v2n[u].size();
        vertices[u] = u;
    }
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
    globalAttr = simpleOrder(query, cs, globalAttr, repetitions);
    std::vector<ui> numBackWard;
    localAttr = simpleOrder(query, cs, globalAttr, numBackWard, localAttr, localAttr.size());
    t.globalOrder = globalAttr;
    for (VertexID u : localAttr)
        t.globalOrder.push_back(u);
    std::vector<std::vector<VertexID>> posToNID(t.globalOrder.size());
    std::vector<ui> lengths(t.numNodes, 0);
    // node orders
    for (VertexID nID = 0; nID < t.numNodes - 1; ++nID) {
        HyperNode &tau = t.nodes[nID];
        // project the global order to the local order
        std::vector<VertexID> nodeOrder;
        for (int pos = 0; pos < t.numAttributes; ++pos) {
            VertexID u = t.globalOrder[pos];
            bool exists = false;
            for (int i = 0; i < tau.numAttributes; ++i) {
                if (tau.attributes[i] == u) {
                    exists = true;
                    break;
                }
            }
            if (exists) {
                nodeOrder.push_back(u);
                posToNID[pos].push_back(nID);
            }
        }
        std::vector<VertexID> reverseOrder = nodeOrder;
        std::reverse(reverseOrder.begin(), reverseOrder.end());
        std::vector<std::vector<VertexID>> vertexParents;
        std::vector<VertexID> cartesianParent;
        initPoses(reverseOrder, query, vertexParents, cartesianParent, cs.dist);
        std::vector<ui> poses(query.getNumVertices(), 0);
        std::vector<double> estimation;
        cardEstimateLabeled(reverseOrder,vertexParents, cartesianParent, cs, visited, partMatch,
                            candidates, candCount, poses, estimation);
        for (ui i = tau.numAttributes - 1; i >= 0; --i) {
            size_t memUsage = estimation[i] * tau.numAttributes * sizeof(VertexID);
            if (t.newGlobalNode && memUsage < budget / (t.numNodes - 1)) break;
            else if (!t.newGlobalNode && memUsage < budget / t.numNodes) break;
            ++lengths[nID];
        }
        std::vector<VertexID> prefix, localOrder;
        prefix.assign(nodeOrder.begin(), nodeOrder.begin() + lengths[nID]);
        localOrder.assign(nodeOrder.begin() + lengths[nID], nodeOrder.end());
        tau.prefixSize = lengths[nID];
        delete[] tau.prefix;
        tau.prefix = new VertexID [tau.prefixSize];
        std::vector<ui> numBackWard;
        localOrder = simpleOrder(query, cs, prefix, numBackWard, localOrder, localOrder.size());
        for (int i = 0; i < tau.prefixSize; ++i) {
            tau.attributes[i] = tau.prefix[i] = prefix[i];
        }
        for (int i = 0; i < tau.numAttributes - tau.prefixSize; ++i) {
            tau.attributes[i + tau.prefixSize] = localOrder[i];
        }
        tau.initPoses(query, cs.dist);
    }
    std::vector<std::vector<VertexID>> nIDsToCall(t.globalOrder.size());
    std::vector<ui> currentLength(t.numNodes, 0);
    // global node order
    for (int pos = 0; pos < t.globalOrder.size(); ++pos) {
        for (VertexID nID : posToNID[pos]) {
            ++currentLength[nID];
            if (currentLength[nID] == lengths[nID]) {
                nIDsToCall[pos].push_back(nID);
            }
        }
    }
    int defaultLen = 0;
    for (int i = 0; i < nIDsToCall.size(); ++i) {
        if (!nIDsToCall[i].empty())
            defaultLen = i + 1;
    }
    if (!t.newGlobalNode) {
        HyperNode &globalNode = t.nodes[t.numNodes - 1];
        std::vector<VertexID> remainingAttrs;
        for (int j = defaultLen; j < t.numAttributes; ++j) {
            VertexID u = t.globalOrder[j];
            for (int k = 0; k < globalNode.numAttributes; ++k) {
                if (globalNode.attributes[k] == u) {
                    remainingAttrs.push_back(u);
                    break;
                }
            }
        }
        if (!remainingAttrs.empty()) {
            std::vector<VertexID> reverseOrder = remainingAttrs;
            std::reverse(reverseOrder.begin(), reverseOrder.end());
            std::vector<std::vector<VertexID>> vertexParents;
            std::vector<VertexID> cartesianParent;
            initPoses(reverseOrder, query, vertexParents, cartesianParent, cs.dist);
            std::vector<ui> poses(query.getNumVertices(), 0);
            std::vector<double> estimation;
            cardEstimateLabeled(reverseOrder,vertexParents, cartesianParent, cs, visited, partMatch,
                                candidates, candCount, poses, estimation);
            for (ui i = globalNode.numAttributes - 1; i >= 0; --i) {
                size_t memUsage = estimation[i] * globalNode.numAttributes * sizeof(VertexID);
                if (memUsage < budget / t.numNodes) break;
                ++defaultLen;
            }
        }
    }
    std::vector<VertexID> defaultPartition(t.globalOrder.begin(), t.globalOrder.begin() + defaultLen);
    HyperNode &globalNode = t.nodes[t.numNodes - 1];
    delete[] globalNode.prefix;
    globalNode.prefix = new VertexID[globalNode.numAttributes];
    globalNode.prefixSize = 0;
    for (int i = 0; i < defaultPartition.size(); ++i) {
        VertexID u = defaultPartition[i];
        for (int j = 0; j < globalNode.numAttributes; ++j) {
            if (globalNode.attributes[j] == u) {
                globalNode.prefix[globalNode.prefixSize] = u;
                ++globalNode.prefixSize;
                break;
            }
        }
    }
    lengths[t.numNodes - 1] = globalNode.prefixSize;
    // project t.globalOrder to the order of global node
    for (int i = 0; i < globalNode.numAttributes; ++i)
        globalNode.attributes[i] = t.globalOrder[i];
    globalNode.initPoses(sharedAttrs, query, cs.dist, true);
    t.defaultPartition = defaultPartition;
    t.initPoses(query, cs, false);
    // build the prefix tree structure
    pt = new PrefixNode(99);
    PrefixNode *pn = pt;
    for (VertexID nID = 0; nID < t.numNodes - 1; ++nID) {
        if (lengths[nID] == 0) {
            pt -> nIDsToCall.push_back(nID);
            if (nID != t.numNodes - 1) pt -> nIDsToBuild.push_back(nID);
        }
    }
    for (int i = 0; i < defaultPartition.size(); ++i) {
        VertexID u = defaultPartition[i];
        PrefixNode *current = new PrefixNode(u);
        pn->children.push_back(current);
        for (VertexID nID : nIDsToCall[i]) {
            current -> nIDsToCall.push_back(nID);
            current -> nIDsToBuild.push_back(nID);
        }
        pn = current;
    }
    pn->nIDsToCall.push_back(t.numNodes - 1);
    pt->initPoses(sharedAttrs, query, cs.dist);
}

void
largeNodePlan(const Graph &query, HyperTree &t, const vector<VertexID> &largeNodes, CandidateSpace &cs, PrefixNode *&pt,
              bool *visited, VertexID *partMatch, VertexID ***candidates, ui **candCount, size_t budget,
              std::vector<double> &costs, std::vector<VertexID> &defaultPartition, std::vector<ui> &lengths) {
    std::vector<VertexID> globalAttr, localAttr;
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
    std::vector<int> repetitions(query.getNumVertices(), 1);
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
    std::vector<VertexID> partitionOrder = globalAttr;
    for (VertexID u: localAttr) partitionOrder.push_back(u);
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
        cardEstimateLabeled(orders[nID], vertexParents, cartesianParent, cs, visited, partMatch, candidates[nID],
                            candCount[nID], poses, estimations[nID]);
    }
    size_t totalMem = 0;
    lengths = std::vector<ui>(t.numNodes, 0);
    std::vector<std::vector<VertexID>> prefixes(t.numNodes);
    std::vector<size_t> currentMem(t.numNodes, 0);
    std::vector<bool> next(t.numNodes, true);
    std::vector<std::vector<VertexID>> nIDsToCall(partitionOrder.size());
    if (!largeNodes.empty()) {
        for (VertexID nID : largeNodes) {
            currentMem[nID] = estimations[nID].back() * t.nodes[nID].numAttributes * sizeof(VertexID);
            totalMem += currentMem[nID];
        }
        for (int i = 0; i < partitionOrder.size(); ++i) {
            for (VertexID nID : posToNID[i]) {
                if (!next[nID]) continue;
                totalMem -= currentMem[nID];
                double prefixCard = estimations[nID][lengths[nID]];
                if (nID == t.numNodes && lengths[nID] < t.extendLevel) prefixCard = estimations[nID][t.extendLevel];
                ++lengths[nID];
                currentMem[nID] = estimations[nID].back() / prefixCard * (t.nodes[nID].numAttributes - lengths[nID]) * sizeof(VertexID);
                prefixes[nID].push_back(partitionOrder[i]);
                if (currentMem[nID] < budget / largeNodes.size()) {
                    next[nID] = false;
                    nIDsToCall[i].push_back(nID);
                }
                totalMem += currentMem[nID];
                if (totalMem < budget) break;
            }
            if (totalMem < budget) {
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
    // the local orders for materialized nodes is built elsewhere
    // large nodes
    PrefixNode *pn = pt;
    for (int i = 0; i < defaultPartition.size(); ++i) {
        VertexID u = defaultPartition[i];
        PrefixNode *current = new PrefixNode(u);
        pn->children.push_back(current);
        for (VertexID nID : nIDsToCall[i]) {
            if (nID == t.numNodes - 1) continue;
            current -> nIDsToCall.push_back(nID);
            current -> nIDsToBuild.push_back(nID);
            std::vector<ui> numBackWard;
            std::vector<VertexID> localOrder;
            for (int j = 0; j < t.nodes[nID].numAttributes; ++j) {
                VertexID u2 = t.nodes[nID].attributes[j];
                if (std::find(prefixes[nID].begin(), prefixes[nID].end(), u2) == prefixes[nID].end())
                    localOrder.push_back(u2);
            }
            localOrder = simpleOrder(query, cs, prefixes[nID], numBackWard, localOrder, localOrder.size());
            t.nodes[nID].prefixSize = lengths[nID];
            delete[] t.nodes[nID].prefix;
            t.nodes[nID].prefix = new VertexID [lengths[nID]];
            for (int j = 0; j < lengths[nID]; ++j) {
                t.nodes[nID].prefix[j] = t.nodes[nID].attributes[j] = prefixes[nID][j];
            }
            for (int j = 0; j < localOrder.size(); ++j) {
                t.nodes[nID].attributes[j + lengths[nID]] = localOrder[j];
            }
            t.nodes[nID].initPoses(query, cs.dist);
        }
        pn = current;
    }
    // global node
    t.globalOrder = partitionOrder;
    t.defaultPartition = defaultPartition;
    pn->nIDsToCall.push_back(t.numNodes - 1);
    std::vector<std::vector<VertexID>> sharedAttrs(t.numNodes);
    for (VertexID nID2 = 0; nID2 < t.numNodes; ++nID2) {
        for (int i = 0; i < t.nodes[nID2].numAttributes; ++i) {
            VertexID u = t.nodes[nID2].attributes[i];
            if (t.v2n[u].size() > 1) sharedAttrs[nID2].push_back(u);
        }
        std::sort(sharedAttrs[nID2].begin(), sharedAttrs[nID2].end());
    }
    HyperNode &globalNode = t.nodes[t.numNodes - 1];
    delete[] globalNode.prefix;
    globalNode.prefix = new VertexID[globalNode.numAttributes];
    globalNode.prefixSize = 0;
    for (int i = 0; i < defaultPartition.size(); ++i) {
        VertexID u = defaultPartition[i];
        for (int j = 0; j < globalNode.numAttributes; ++j) {
            if (globalNode.attributes[j] == u) {
                globalNode.prefix[globalNode.prefixSize] = u;
                ++globalNode.prefixSize;
                break;
            }
        }
    }
    lengths[t.numNodes - 1] = globalNode.prefixSize;
    // project t.globalOrder to the order of global node
    for (int i = 0; i < globalNode.numAttributes; ++i)
        globalNode.attributes[i] = t.globalOrder[i];
    globalNode.initPoses(sharedAttrs, query, cs.dist, true);
    t.initPoses(query, cs, false);
    pt->initPoses(sharedAttrs, query, cs.dist);
}

bool adaptiveNodeJoin(const HyperTree &t, VertexID nID, CandidateSpace &cs, bool *visited,
                      const std::vector<VertexID> &notInBag, VertexID *partMatch, int mappingSize,
                      VertexID **candidates, ui *candCount, std::vector<ui> &poses,
                      std::vector<std::vector<VertexID>> &tuples, size_t &budget) {
    for (VertexID u2: notInBag) visited[partMatch[u2]] = false;
    const HyperNode &tau = t.nodes[nID];
    const VertexID *order = tau.attributes;
    if (mappingSize == tau.numAttributes) {
        std::vector<VertexID> tuple(tau.numAttributes - tau.prefixSize);
        for (int i = 0; i < tau.numAttributes - tau.prefixSize; ++i) {
            tuple[i] = partMatch[tau.attributes[i + tau.prefixSize]];
        }
        tuples.push_back(tuple);
        budget -= sizeof(VertexID) * (tau.numAttributes - tau.prefixSize);
        for (VertexID u2: notInBag) visited[partMatch[u2]] = true;
        return true;
    }
    const std::vector<std::vector<VertexID>> &vertexParents = tau.attributesBefore;
    // handle the first level
    VertexID u = order[mappingSize];
    VertexID **neighbors = new VertexID *[tau.numAttributes];
    ui *neighborCount = new ui[tau.numAttributes];
    if (poses[mappingSize] == 0) {
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
    }
    int depth = mappingSize;
    while (depth >= mappingSize) {
        while (poses[depth] < candCount[depth]) {
            VertexID v = candidates[depth][poses[depth]];
            if (visited[v]) {
                ++poses[depth];
                continue;
            }
            visited[v] = true;
            partMatch[order[depth]] = v;
            if (depth + 1 == tau.numAttributes) {
                if (budget <= (tau.numAttributes - tau.prefixSize) * sizeof(VertexID)) {
                    delete[] neighbors;
                    delete[] neighborCount;
                    for (int i = tau.prefixSize; i < tau.numAttributes; ++i) {
                        visited[partMatch[order[i]]] = false;
                    }
                    for (VertexID u2: notInBag) visited[partMatch[u2]] = true;
                    return false;
                }
                std::vector<VertexID> tuple(tau.numAttributes - tau.prefixSize);
                for (int i = 0; i < tau.numAttributes - tau.prefixSize; ++i) {
                    tuple[i] = partMatch[tau.attributes[i + tau.prefixSize]];
                }
                tuples.push_back(tuple);
                budget -= sizeof(VertexID) * (tau.numAttributes - tau.prefixSize);
                visited[partMatch[tau.attributes[depth]]] = false;
                ++poses[depth];
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
        if (depth >= mappingSize) {
            ++poses[depth];
            visited[partMatch[order[depth]]] = false;
        }
    }
    delete[] neighbors;
    delete[] neighborCount;
    for (VertexID u2: notInBag) visited[partMatch[u2]] = true;
    return true;
}

void adaptiveGlobalJoin(std::vector<std::vector<VertexID>> &result, size_t &budget, size_t &count, const Graph &query,
                        HyperTree &t, CandidateSpace &cs, VertexID *partMatch, int mappingSize, int pathLength,
                        std::vector<TrieNode *> &nodes, bool *visited, const std::vector<std::vector<VertexID>> &nIDs,
                        const std::vector<std::vector<VertexID>> &vertexParents,
                        const std::vector<VertexID> &cartesianParent, ui ***iters, ui *iterSize, bool traversal,
                        std::vector<std::vector<DynamicArray<TrieNode *> *>> &edgeColumns, TrieNode *empty, bool skip) {
    if (skip) return;
#ifdef COLLET_GLOBAL_TIME
    auto start = std::chrono::steady_clock::now();
#endif
    for (VertexID nID = 0; nID < t.numNodes - 1; ++nID) {
        if (nodes[nID] == empty) return;
    }
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
    pathLength += 1;
    // for each level and each nID, the trie node is traversedNodes[]
    DynamicArray<TrieNode *> ** children = new DynamicArray<TrieNode *> *[t.numAttributes];
    // initialize the first level
    if (mappingSize == 0 && nIDs[0].empty()) {
        VertexID u = order[0];
        edgeColumns[pathLength][0] ->setSize(cs.candidateSet[u].size());
        for (int i = 0; i < cs.candidateSet[u].size(); ++i) {
            edgeColumns[pathLength][0][0][i] -> value = cs.candidateSet[u][i];
            iters[pathLength][i][0] = i;
        }
        iterSize[pathLength] = cs.candidateSet[u].size();
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
    std::vector<int> depthToExtendLevels(globalNode.numAttributes, t.extendLevel);
    std::map<int, std::vector<VertexID>> extendLevelToTrieOrder;
    extendLevelToTrieOrder[t.extendLevel] = t.trieOrder.back();
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
            int currentExtendLevel = t.extendLevel;
            if (depth > 0) currentExtendLevel = depthToExtendLevels[depth - 1];
            depthToExtendLevels[depth] = currentExtendLevel;
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
                if (currentExtendLevel != globalNode.numAttributes) {
                    if (budget > (globalNode.numAttributes - currentExtendLevel) * sizeof(VertexID)) {
                        std::vector<VertexID> tuple(globalNode.numAttributes - currentExtendLevel);
                        for (int i = 0; i < globalNode.numAttributes - currentExtendLevel; ++i) {
                            tuple[i] = partMatch[globalNode.attributes[i + currentExtendLevel]];
                        }
                        tuples.push_back(tuple);
                        budget -= (globalNode.numAttributes - currentExtendLevel) * sizeof(VertexID);
                    }
                    else {
                        // change the extend level
                        budget += tuples.size() * tuples[0].size() * sizeof(VertexID);
                        tuples.clear();
                        // backtrack to a depth
                        for (int i = currentExtendLevel; i < globalNode.numAttributes; ++i) {
                            visited[partMatch[order[i]]] = false;
                        }
                        if (currentExtendLevel > mappingSize)
                            depthToExtendLevels[currentExtendLevel - mappingSize - 1] = currentExtendLevel + 1;
                        else ++t.extendLevel;
                        t.buildTraverseAdaptive(query, currentExtendLevel + 1);
                        depth = currentExtendLevel - mappingSize;
                        poses[mappingSize + depth] = 0;
                    }
                }
                else {
                    // when a match of shared attributes is found, traverse the remaining attributes
                    traverse(count, query, t, t.globalOrder, partMatch, mappingSize + depth + 1, lastNodes, visited, currentExtendLevel);
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
        int currentExtendLevel = t.extendLevel;
        if (depth >= 0) currentExtendLevel = depthToExtendLevels[depth];
        if (depth + 1 + mappingSize == currentExtendLevel && !tuples.empty()) {
            budget += tuples.size() * (globalNode.numAttributes - currentExtendLevel) * sizeof(VertexID);
            std::vector<VertexID> &trieOrder = t.trieOrder.back();
            if (extendLevelToTrieOrder.find(currentExtendLevel) == extendLevelToTrieOrder.end()) {
                 std::vector<VertexID> newOrder;
                 for (ui j = 0; j < t.globalOrder.size(); ++j) {
                    VertexID u = t.globalOrder[j];
                     ui k = currentExtendLevel;
                     for (; k < globalNode.numAttributes; ++k) {
                         if (globalNode.attributes[k] == u) {
                             break;
                         }
                     }
                     if (k != globalNode.numAttributes) newOrder.push_back(k - currentExtendLevel);
                 }
                 trieOrder = newOrder;
                 extendLevelToTrieOrder[currentExtendLevel] = newOrder;
            }
            else trieOrder = extendLevelToTrieOrder[currentExtendLevel];
            buildTrie(tuples, lastNodes[t.numNodes - 1], trieOrder);
            // when a match of shared attributes is found, traverse the remaining attributes
            traverse(count, query, t, t.globalOrder, partMatch, mappingSize + depth + 1, lastNodes, visited, currentExtendLevel);
            for (int i = 0; i < lastNodes[t.numNodes - 1]->nodeChild.size(); ++i) {
                delete lastNodes[t.numNodes - 1]->nodeChild[i];
            }
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
#ifdef COLLET_GLOBAL_TIME
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsedSeconds = end - start;
    gGlobalTime += elapsedSeconds.count();
#endif
}

void materializeJoin(HyperTree &t, PrefixNode *pt, CandidateSpace &cs, bool *visited,
                     vector<vector<std::vector<VertexID>>> &tuples, std::vector<VertexID> &largeNodes, size_t &budget,
                     ui height, VertexID *partMatch, VertexID **pCandidates, ui *pCandCount, VertexID **neighbors,
                     ui *neighborCount, std::vector<std::vector<ui>> &nodePoses, std::vector<ui> &pPoses,
                     std::vector<ui> &childPoses, ui ***nodeCandidates, ui **nodeCandCount) {
    PrefixNode *pn = pt;
    std::vector<PrefixNode *> nodes(height, nullptr);
    int depth = 0;
    std::vector<VertexID> notInBag;
    for (int i = 0; i < pt -> nIDsToCall.size(); ++i) {
        VertexID nID = pt -> nIDsToCall[i];
        if (nID != t.numNodes - 1) {
            for (int j = 0; j < nodePoses[nID].size(); ++j)
                nodePoses[nID][j] = 0;
            if (!adaptiveNodeJoin(t, nID, cs, visited, notInBag, partMatch, 0, nodeCandidates[nID],
                                  nodeCandCount[nID],
                                  nodePoses[nID], tuples[nID], budget)) {
                if (!tuples[nID].empty())
                    budget += tuples[nID].size() * tuples[nID][0].size() * sizeof(VertexID);
                tuples[nID].clear();
                largeNodes.push_back(nID);
                // remove nID in the branch: the nIDsToBuild for the global node, and the nIDsToCall for the current node
                removeNID(nID, pt);
                --i;
                for (auto it = pt -> nIDsToBuild.begin(); it != pt -> nIDsToBuild.end();) {
                    if (*it == nID) it = pt -> nIDsToBuild.erase(it);
                    else ++it;
                }
            }
        }
    }
    if (!pn -> children.empty()) {
        pCandidates[0] = cs.candidateSet[pn -> children[0] -> u].data();
        pCandCount[0] = cs.candidateSet[pn -> children[0] -> u].size();
    }
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            nodes[depth] = current;
            VertexID u = current -> u;
            bool nextLevel = false;
            while (pPoses[depth] < pCandCount[depth]) {
                VertexID v = pCandidates[depth][pPoses[depth]];
                ++pPoses[depth];
                if (visited[v]) continue;
                visited[v] = true;
                partMatch[u] = v;
                for (int i = 0; i < current -> nIDsToCall.size(); ++i) {
                    VertexID nID = current -> nIDsToCall[i];
                    for (int j = 0; j < nodePoses[nID].size(); ++j)
                        nodePoses[nID][j] = 0;
                    if (!adaptiveNodeJoin(t, nID, cs, visited, notInBag, partMatch, depth + 1,
                                          nodeCandidates[nID],
                                          nodeCandCount[nID], nodePoses[nID], tuples[nID], budget)) {
                        if (!tuples[nID].empty())
                            budget += tuples[nID].size() * tuples[nID][0].size() * sizeof(VertexID);
                        tuples[nID].clear();
                        largeNodes.push_back(nID);
                        removeNID(nID, nodes, depth, pt);
                    }
                }
                if (!current -> children.empty()) {
                    nextLevel = true;
                    break;
                }
                else visited[v] = false;
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
        --depth;
        if (depth >= 0) {
            const PrefixNode *current = pn;
            if (depth == 0) pn = pt;
            else pn = nodes[depth - 1];
            VertexID u = current -> u;
            VertexID v = partMatch[u];
            visited[v] = false;
        }
    }
}

void adaptiveMaterializeJoin(const Graph &query, HyperTree &t, PrefixNode *pt, CandidateSpace &cs,
                             std::vector<TrieNode *> &trieNodes, bool *visited, std::vector<std::vector<VertexID>> &result,
                             size_t &count, bool traverse, size_t budget, bool skip) {
    t.numAttributes = query.getNumVertices();
    delete[] t.v2n;
    t.v2n = new std::vector<VertexID> [query.getNumVertices()];
    for (ui i = 0; i < t.numNodes; ++i) {
        for (ui j = 0; j < t.nodes[i].numAttributes; ++j) {
            t.v2n[t.nodes[i].attributes[j]].push_back(i);
        }
    }
    std::unordered_map<std::vector<ui>, PrefixNode *, VectorHash, VectorEqual> lengthToPT;
    std::unordered_map<std::vector<ui>, HyperTree *, VectorHash, VectorEqual> lengthToTree;
    std::unordered_map<std::vector<ui>, std::vector<int>, VectorHash, VectorEqual> lengthToMappingSizes;
    std::vector<std::vector<std::vector<VertexID>>> tuples(t.numNodes);
    VertexID *partMatch = new VertexID[t.numAttributes];
    TrieNode *empty = new TrieNode();
    ui maxSize = cs.getMaxSize();
    ui height = t.numAttributes + 1;
    VertexID **pCandidates = new VertexID *[height];
    ui *pCandCount = new ui[height];
    for (int i = 0; i < height; ++i) {
        pCandidates[i] = new VertexID[maxSize];
        pCandCount[i] = 0;
    }
    std::vector<double> costs;
    std::vector<VertexID> largeNodes;
    std::vector<ui> pPoses(height, 0);
    std::vector<ui> childPoses(height, 0);
    VertexID **neighbors = new VertexID *[height];
    ui *neighborCount = new ui[height];
    DynamicArray<TrieNode *> ** children = new DynamicArray<TrieNode *> *[t.numAttributes];
    std::vector<std::vector<TrieNode *>> traversedNodes(trieNodes.size());
    std::vector<TrieNode *> lastNodes(trieNodes.size());
    for (VertexID nID = 0; nID < trieNodes.size(); ++nID) {
        traversedNodes[nID].push_back(trieNodes[nID]);
        lastNodes[nID] = trieNodes[nID];
    }
    VertexID ***nodeCandidates = new VertexID **[trieNodes.size()];
    ui **nodeCandCount = new ui *[trieNodes.size()];
    std::vector<std::vector<ui>> nodePoses(trieNodes.size());
    for (VertexID nID = 0; nID < trieNodes.size(); ++nID) {
        const HyperNode &tau = t.nodes[nID];
        nodeCandidates[nID]= new VertexID *[tau.numAttributes];
        nodeCandCount[nID] = new ui[tau.numAttributes];
        nodePoses[nID] = std::vector<ui>(tau.numAttributes, 0);
        ui maxSize = 0;
        for (int i = 0; i < tau.numAttributes; ++i) {
            if (maxSize < cs.candidateSet[tau.attributes[i]].size())
                maxSize = cs.candidateSet[tau.attributes[i]].size();
        }
        for (int i = 0; i < tau.numAttributes; ++i) {
            nodeCandidates[nID][i] = new VertexID [maxSize];
            nodeCandCount[nID][i] = 0;
        }
    }
    std::vector<VertexID> defaultLengths = std::vector<ui>(t.numNodes, 0);
    for (VertexID nID = 0; nID < t.numNodes; ++nID) defaultLengths[nID] = t.nodes[nID].prefixSize;
    ui ***iters = new ui **[height];
    ui *iterSizes = new ui [height];
    for (int i = 0; i < height; ++i) {
        iters[i] = new ui *[maxSize];
        for (int j = 0; j < maxSize; ++j) {
            iters[i][j] = new ui[t.numAttributes];
        }
        iterSizes[i] = 0;
    }
    iterSizes[0] = 1;
    const HyperNode &globalNode = t.nodes[t.numNodes - 1];
    ui maxEdgeSize = query.getMaxDegree();
    // join relational columns together with edge columns
    std::vector<std::vector<DynamicArray<TrieNode *> *>> edgeColumns(height);
    for (int i = 0 ; i < height; ++i) {
        edgeColumns[i] = std::vector<DynamicArray<TrieNode *> *>(maxEdgeSize);
        for (int j = 0; j < maxEdgeSize; ++j) {
            edgeColumns[i][j] = new DynamicArray<TrieNode *>();
            for (int k = 0; k < maxSize; ++k) {
                TrieNode * pointer = new TrieNode();
                edgeColumns[i][j]->push_back(pointer);
            }
        }
    }
    std::vector<size_t> tupleSizes(t.numNodes, 0);
    std::vector<PrefixNode *> nodes(height, nullptr);
    std::vector<HyperTree *> trees(height);
    std::vector<std::vector<int>> mappingSizes(height);
    std::vector<std::vector<ui>> depthToLengths(height, defaultLengths);
    std::vector<int> defaultSizes = getMappingSizes(t, pt);
    int depth = 0, globalDepth = 0;
    PrefixNode *tmp = new PrefixNode(99);
    tmp->children.push_back(pt);
    PrefixNode *pn = tmp;
    if (lengthToPT.empty()) {
        lengthToPT[defaultLengths] = pt;
        lengthToTree[defaultLengths] = &t;
        defaultSizes = getMappingSizes(t, pt);
    }
    std::vector<VertexID> notInBag;
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            nodes[depth] = current;
            HyperTree *currentT;
            std::vector<int> currentSize;
            if (globalDepth == 0) {
                currentT = &t;
                currentSize = defaultSizes;
                depthToLengths[0] = defaultLengths;
            }
            else {
                currentT = trees[globalDepth - 1];
                currentSize = mappingSizes[globalDepth - 1];
                depthToLengths[globalDepth] = depthToLengths[globalDepth - 1];
            }
            trees[depth] = currentT;
            mappingSizes[depth] = currentSize;
            VertexID u = current -> u;
            bool nextLevel = false, newBranch = false;
            if (!current -> pathToGlobal) {
                while (pPoses[depth] < pCandCount[depth]) {
                    VertexID v = pCandidates[depth][pPoses[depth]];
                    ++pPoses[depth];
                    if (visited[v]) continue;
                    visited[v] = true;
                    partMatch[u] = v;
                    bool deadend = false;
                    for (VertexID nID: current -> nIDsToCall) {
                        for (int i = 0; i < nodePoses[nID].size(); ++i)
                            nodePoses[nID][i] = 0;
                        ui oldSize = tuples[nID].size();
                        if (!adaptiveNodeJoin(*currentT, nID, cs, visited, notInBag, partMatch, currentSize[nID],
                                              nodeCandidates[nID],
                                              nodeCandCount[nID], nodePoses[nID], tuples[nID], budget)) {
                            std::vector<ui> length = depthToLengths[globalDepth];
                            PrefixNode *old = lengthToPT[length];
                            ++length[nID];
                            if (!tuples[nID].empty())
                                budget += tuples[nID].size() * tuples[nID][0].size() * sizeof(VertexID);
                            tuples[nID].clear();
                            // replace the subtree rooted at nodes[globalDepth]
                            PrefixNode *newPrefix;
                            HyperTree *newT;
                            std::vector<int> newMappingSize;
                            if (lengthToPT.find(length) != lengthToPT.end()) {
                                newPrefix = lengthToPT[length];
                                newT = lengthToTree[length];
                                newMappingSize = lengthToMappingSizes[length];
                            }
                            else {
                                extendPrefixTree(depth, query, *currentT, t.globalOrder, nID, length,
                                                 newPrefix, newT, old, childPoses, cs);
                                lengthToPT[length] = newPrefix;
                                lengthToTree[length] = newT;
                                newMappingSize = getMappingSizes(*newT, newPrefix);
                                lengthToMappingSizes[length] = newMappingSize;
                            }
                            // backtracking the prefix tree
                            int newDepth = -1;
                            PrefixNode *pn1 = old, *pn2 = newPrefix;
                            while (pn1 -> u == pn2 -> u) {
                                ++newDepth;
                                if (!pn1->children.empty() && pn2->children.size() > childPoses[newDepth + 1]) {
                                    pn1 = pn1 -> children[childPoses[newDepth + 1]];
                                    pn2 = pn2 -> children[childPoses[newDepth + 1]];
                                }
                                else break;
                                if (newDepth == depth) break;
                            }
                            for (int i = newDepth + 1; i <= depth; ++i) {
                                childPoses[i] = 0;
                                visited[partMatch[nodes[i]->u]] = false;
                            }
                            if (newDepth != depth) newBranch = true;
                            depth = newDepth;
                            // replacing the prefix tree and orders
                            currentT = trees[globalDepth] = newT;
                            depthToLengths[globalDepth] = length;
                            currentSize = mappingSizes[globalDepth] = newMappingSize;
                            pn2 = newPrefix;
                            for (int i = 0; i < globalDepth; ++i) pn2 = pn2 -> children[childPoses[i + 1]];
                            for (int i = globalDepth; i <= depth; ++i) {
                                nodes[i] = pn2;
                                if (i != depth) pn2 = pn2 -> children[childPoses[i + 1]];
                            }
                            if (depth != -1) current = nodes[depth];
                            else current = newPrefix;
                            for (int i = 0; i <= depth; ++i) {
                                if (nodes[i]->pathToGlobal) globalDepth = i;
                            }
                            if (u != current->u) newBranch = true;
                            if (globalDepth == 0) {
                                lengthToTree.erase(defaultLengths);
                                defaultLengths = length;
                                pt = lengthToPT[defaultLengths];
                                tmp->children[0] = pt;
                                t = *lengthToTree[defaultLengths];
                                defaultSizes = lengthToMappingSizes[defaultLengths];
                            }
                            if (depth == 0) pn = tmp;
                            else pn = nodes[depth - 1];
                        }
                        else if (tuples[nID].size() == oldSize) {
                            deadend = true;
                            break;
                        }
                    }
                    if (deadend) {
                        visited[v] = false;
                        continue;
                    }
                    if (newBranch) {
                        pPoses[depth] = 0;
                        break;
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
                    if (depth != 0) {
                        ui childPos = iters[depth][pPoses[depth]][0];
                        if (current -> nIDsToJoin.empty()) {
                            v = edgeColumns[depth][0][0][childPos]->value;
                        }
                        else {
                            VertexID u2 = current -> nIDsToJoin[0];
                            v = traversedNodes[u2].back()->nodeChild[childPos]->value;
                        }
                        ++pPoses[depth];
                        if (visited[v]) continue;
                        visited[v] = true;
                        partMatch[u] = v;
                    }
                    else ++pPoses[depth];
                    bool deadend = false;
                    for (VertexID nID: current -> nIDsToCall) {
                        if (nID != currentT->numNodes - 1) {
                            for (int i = 0; i < nodePoses[nID].size(); ++i)
                                nodePoses[nID][i] = 0;
                            ui oldSize = tuples[nID].size();
                            if (!adaptiveNodeJoin(*currentT, nID, cs, visited, notInBag, partMatch,
                                                  currentSize[nID],
                                                  nodeCandidates[nID],
                                                  nodeCandCount[nID], nodePoses[nID], tuples[nID], budget)) {
                                std::vector<ui> length = depthToLengths[globalDepth];
                                PrefixNode *old = lengthToPT[length];
                                ++length[nID];
                                if (!tuples[nID].empty())
                                    budget += tuples[nID].size() * tuples[nID][0].size() * sizeof(VertexID);
                                tuples[nID].clear();
                                // replace the subtree rooted at nodes[globalDepth]
                                PrefixNode *newPrefix;
                                HyperTree *newT;
                                std::vector<int> newMappingSize;
                                if (lengthToPT.find(length) != lengthToPT.end()) {
                                    newPrefix = lengthToPT[length];
                                    newT = lengthToTree[length];
                                    newMappingSize = lengthToMappingSizes[length];
                                }
                                else {
                                    extendPrefixTree(depth, query, *currentT, t.globalOrder, nID, length,
                                                     newPrefix, newT, old, childPoses, cs);
                                    lengthToPT[length] = newPrefix;
                                    lengthToTree[length] = newT;
                                    newMappingSize = getMappingSizes(*newT, newPrefix);
                                    lengthToMappingSizes[length] = newMappingSize;
                                }
                                currentT = trees[globalDepth] = newT;
                                depthToLengths[globalDepth] = length;
                                currentSize = mappingSizes[globalDepth] = newMappingSize;
                                current = newPrefix;
                                for (int i = 1; i <= globalDepth; ++i)
                                    current = current->children[childPoses[i]];
                                nodes[depth] = current;
                                if (globalDepth == 0) {
                                    lengthToTree.erase(defaultLengths);
                                    defaultLengths = length;
                                    pt = lengthToPT[defaultLengths];
                                    tmp->children[0] = pt;
                                    t = *lengthToTree[defaultLengths];
                                    defaultSizes = lengthToMappingSizes[defaultLengths];
                                }
                            }
                            else if (tuples[nID].size() == oldSize) {
                                for (VertexID nID2 : current -> nIDsToCall) {
                                    if (!tuples[nID2].empty()) budget += tuples[nID2][0].size() * tuples[nID2].size() * sizeof(VertexID);
                                    tuples[nID2].clear();
                                }
                                deadend = true;
                                break;
                            }
                        }
                    }
                    if (depth != 0 && deadend) {
                        visited[v] = false;
                        continue;
                    }
                    // all joined relations extend one level
                    for (int i = 0; i < current -> nIDsToJoin.size(); ++i) {
                        VertexID nID = current -> nIDsToJoin[i];
                        ui childPos = iters[depth][pPoses[depth] - 1][i];
                        traversedNodes[nID].push_back(traversedNodes[nID].back()->nodeChild[childPos]);
                        lastNodes[nID] = traversedNodes[nID].back();
                    }
                    if (!current -> children.empty()) {
                        nextLevel = true;
                        break;
                    }
                    else {
                        for (VertexID nID: current -> nIDsToCall) {
                            if (nID == currentT->numNodes - 1) {
                                bool flag = true;
                                for (VertexID nID2 : current -> nIDsToBuild) {
                                    if (tuples[nID2].empty()) {
                                        flag = false;
                                        break;
                                    }
                                }
                                if (flag) {
                                    for (VertexID nID2 : current -> nIDsToBuild) {
                                        budget += tupleSizes[nID2];
                                        tupleSizes[nID2] = tuples[nID2].size() * tuples[nID2][0].size() * sizeof(VertexID);
                                        buildTrie(tuples[nID2], trieNodes[nID2], currentT->trieOrder[nID2]);
                                    }
                                    adaptiveGlobalJoin(result, budget, count, query, *currentT, cs, partMatch,
                                                       currentSize[nID], depth,
                                                       lastNodes, visited, globalNode.nIDs, globalNode.attributesBefore,
                                                       globalNode.cartesianParent, iters, iterSizes, traverse,
                                                       edgeColumns, empty, skip);
                                }
                                else {
                                    for (VertexID nID2 : current -> nIDsToBuild) {
                                        if (!tuples[nID2].empty()) {
                                            budget += tuples[nID2].size() * tuples[nID2][0].size() * sizeof(VertexID);
                                            tuples[nID2].clear();
                                        }
                                    }
                                }
                            }
                        }
                        if (current->u != 99) visited[v] = false;
                        for (VertexID nID : current -> nIDsToJoin) {
                            traversedNodes[nID].pop_back();
                            lastNodes[nID] = traversedNodes[nID].back();
                        }
                    }
                }
            }
            if (nextLevel) {
                ++depth;
                pPoses[depth] = 0;
                childPoses[depth] = 0;
                depthToLengths[depth] = depthToLengths[depth - 1];
                pn = current;
                current = pn -> children[0];
            }
            else if (!newBranch) {
                ++childPoses[depth];
                pPoses[depth] = 0;
                if (childPoses[depth] == pn->children.size()) break;
                current = pn->children[childPoses[depth]];
            }
            u = current->u;
            nodes[depth] = current;
            if (current -> u == 99) {
                pPoses[0] = 0;
                iterSizes[0] = 1;
            }
            else if (current -> pathToGlobal) {
                if (!newBranch) ++globalDepth;
                bool flag = true;
                for (VertexID nID : pn -> nIDsToBuild) {
                    if (tuples[nID].empty()) {
                        flag = false;
                        break;
                    }
                }
                if (flag) {
                    for (VertexID nID : pn -> nIDsToBuild) {
                        budget += tupleSizes[nID];
                        tupleSizes[nID] = tuples[nID].size() * tuples[nID][0].size() * sizeof(VertexID);
                        buildTrie(tuples[nID], trieNodes[nID], currentT->trieOrder[nID]);
                    }
                }
                else {
                    for (VertexID nID : pn -> nIDsToBuild) {
                        if (!tuples[nID].empty())
                            budget += tuples[nID].size() * tuples[nID][0].size() * sizeof(VertexID);
                        tuples[nID].clear();
                    }
                    if (pn->u == 99) return;
                    visited[partMatch[pn->u]] = false;
                    for (VertexID nID : pn -> nIDsToJoin) {
                        traversedNodes[nID].pop_back();
                        lastNodes[nID] = traversedNodes[nID].back();
                    }
                    --depth;
                    --globalDepth;
                    if (depth == 0) pn = tmp;
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
                    edgeColumns[depth][0] ->setSize(cs.candidateSet[u].size());
                    for (int i = 0; i < cs.candidateSet[u].size(); ++i) {
                        edgeColumns[depth][0][0][i]->value = cs.candidateSet[u][i];
                        iters[depth][i][0] = i;
                    }
                    iterSizes[depth] = cs.candidateSet[u].size();
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
                if (iterSizes[depth] == 0)
                    break;
            }
            else if (current->u != 99){
                if (depth != 1) {
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
                    pCandidates[1] = cs.candidateSet[u].data();
                    pCandCount[1] = cs.candidateSet[u].size();
                }
                if (pCandCount[depth] == 0)
                    break;
            }
        }
        if (nodes[depth]->pathToGlobal) --globalDepth;
        --depth;
        if (depth >= 0) {
            const PrefixNode *current = pn;
            HyperTree *currentT = trees[globalDepth];
            std::vector<int> currentSize = mappingSizes[globalDepth];
            if (depth == 0) pn = tmp;
            else pn = nodes[depth - 1];
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
                        for (VertexID nID2 : current -> nIDsToBuild) {
                            budget += tupleSizes[nID2];
                            tupleSizes[nID2] = tuples[nID2].size() * tuples[nID2][0].size() * sizeof(VertexID);
                            buildTrie(tuples[nID2], trieNodes[nID2], currentT->trieOrder[nID2]);
                        }
                        adaptiveGlobalJoin(result, budget, count, query, *currentT, cs, partMatch, currentSize[nID],
                                           depth, lastNodes, visited, globalNode.nIDs, globalNode.attributesBefore,
                                           globalNode.cartesianParent, iters, iterSizes, traverse, edgeColumns, empty,
                                           skip);
                    }
                    else {
                        for (VertexID nID2 : current -> nIDsToBuild) {
                            if (!tuples[nID2].empty()) {
                                budget += tuples[nID2].size() * tuples[nID2][0].size() * sizeof(VertexID);
                                tuples[nID2].clear();
                            }
                        }
                    }
                }
            }
            if (depth != 0) {
                VertexID u = current -> u;
                VertexID v = partMatch[u];
                visited[v] = false;
            }
            if (current -> pathToGlobal) {
                for (VertexID nID : current -> nIDsToJoin) {
                    traversedNodes[nID].pop_back();
                    lastNodes[nID] = traversedNodes[nID].back();
                }
            }
        }
    }
    tmp->children.clear();
    delete tmp;
    delete[] partMatch;

    for (int i = 0; i < trieNodes.size(); ++i) {
        delete[] nodeCandCount[i];
        for (int j = 0; j < t.nodes[i].numAttributes; ++j) {
            delete[] nodeCandidates[i][j];
        }
        delete[] nodeCandidates[i];
    }
    delete[] nodeCandidates;
    delete[] nodeCandCount;
    delete[] pCandidates;
    delete[] pCandCount;
    delete[] neighbors;
    delete[] neighborCount;
    for (int i = 0; i < t.numAttributes; ++i) {
        for (int j = 0; j < maxSize; ++j) {
            delete[] iters[i][j];
        }
        delete[] iters[i];
    }
    delete[] iters;
    delete[] iterSizes;
    delete[] children;
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < maxEdgeSize; ++j) {
            DynamicArray<TrieNode *> *dynamicArrayPtr = edgeColumns[i][j];
            for (ui k = 0; k < dynamicArrayPtr->size(); ++k) {
                delete (*dynamicArrayPtr)[k];
            }
            dynamicArrayPtr->clear();
            delete dynamicArrayPtr;
        }
    }

    for (auto & it : lengthToPT) {
        if (it.first != defaultLengths)
            delete it.second;
    }

    for (auto & it : lengthToTree) {
        if (it.first != defaultLengths)
            delete it.second;
    }
}

std::vector<ui>
resetLengths(HyperTree &t, const std::vector<ui> &lengths, VertexID nID, const PrefixNode *oldPT, int depth,
             const HyperTree &defaultT, const vector<vector<int>> &attrIDMap) {
    PrefixNode *newPrefix = oldPT -> clone();
    std::vector<PrefixNode *> path = {newPrefix};
    std::vector<PrefixNode *> tmp = newPrefix->locate(nID);
    for (PrefixNode *attribute: tmp) path.push_back(attribute);
    ui newPrefixSize = lengths[nID];
    std::vector<VertexID> newPrefixes(t.nodes[nID].attributes + t.nodes[nID].prefixSize, t.nodes[nID].attributes + newPrefixSize);
    int extendSize = depth;
    if (newPrefixSize <= depth) extendSize = newPrefixSize;
    else {
        bool includes = true;
        for (VertexID u: newPrefixes) {
            bool exists = false;
            for (int i = 0; i < extendSize; ++i) {
                if (defaultT.nodes[nID].attributes[i] == u) {
                    exists = true;
                    break;
                }
            }
            if (!exists) includes = false;
        }
        while (!includes) {
            ++extendSize;
            includes = true;
            for (VertexID u: newPrefixes) {
                bool exists = false;
                for (int i = 0; i < extendSize; ++i) {
                    if (defaultT.nodes[nID].attributes[i] == u) {
                        exists = true;
                        break;
                    }
                }
                if (!exists) includes = false;
            }
        }
        if (extendSize > depth) {
            PrefixNode *lastAttr = path[depth];
            std::vector<VertexID> below = lastAttr->getBagsBelow();
            PrefixNode *newAttr, *attr;
            std::vector<VertexID> oldNIDsToCall = lastAttr->nIDsToCall;
            lastAttr->nIDsToCall.clear();
            std::vector<int> ids(lastAttr->children.size());
            for (int j = 0; j < lastAttr->children.size(); ++j) {
                VertexID firstBag = lastAttr->children[j]->getBagsBelow()[0];
                ids[j] = attrIDMap[firstBag][depth];
            }
            // extend nID and sibling bags by one level
            int pos = oldNIDsToCall.size();
            for (int i = 0; i < oldNIDsToCall.size(); ++i) {
                VertexID nID2 = oldNIDsToCall[i];
                if (nID2 == t.numNodes - 1) {
                    lastAttr->nIDsToCall.push_back(nID2);
                    continue;
                }
                if (nID2 == nID) pos = i;
                if (i >= pos) {
                    VertexID u = defaultT.nodes[nID2].attributes[depth];
                    newAttr = new PrefixNode(u);
                    newAttr->pathToGlobal = false;
                    for (int j = 0; j < t.nodes[nID2].numAttributes; ++j) {
                        if (t.nodes[nID2].attributes[j] == u) {
                            newAttr->attributesBefore = t.nodes[nID2].attributesBefore[j];
                        }
                    }
                    lastAttr->children.push_back(newAttr);
                    ids.push_back(attrIDMap[nID2][depth]);
                    if (i == pos) {
                        attr = newAttr;
                        path.push_back(attr);
                    }
                    else newAttr->nIDsToCall.push_back(nID2);
                }
                else lastAttr->nIDsToCall.push_back(nID2);
            }
            for (int j1 = 0; j1 < lastAttr->children.size() - 1; ++j1) {
                for (int j2 = j1 + 1; j2 < lastAttr->children.size(); ++j2) {
                    if (ids[j2] < ids[j1]) {
                        std::swap(lastAttr->children[j1], lastAttr->children[j2]);
                        std::swap(ids[j1], ids[j2]);
                    }
                }
            }
            newAttr = attr;
            for (int i = depth + 1; i < extendSize; ++i) {
                VertexID u = defaultT.nodes[nID].attributes[i];
                newAttr = new PrefixNode(u);
                newAttr->pathToGlobal = false;
                for (int j = 0; j < t.nodes[nID].numAttributes; ++j) {
                    if (t.nodes[nID].attributes[j] == u) {
                        newAttr->attributesBefore = t.nodes[nID].attributesBefore[j];
                    }
                }
                attr->children.push_back(newAttr);
                attr = newAttr;
                path.push_back(newAttr);
            }
            newAttr->nIDsToCall = {nID};
        }
    }
    std::vector<ui> prefixSizes(t.numNodes);
    for (VertexID nID2 = 0; nID2 < t.numNodes; ++nID2) {
        prefixSizes[nID2] = t.nodes[nID2].prefixSize;
    }
    int newAttrPos = -1, buildPos = 0;
    for (int i = 0; i < path.size(); ++i) {
        PrefixNode *attr = path[i];
        if (attr->u == newPrefixes[0]) {
            newAttrPos = i;
        }
        if (std::find(attr->nIDsToBuild.begin(), attr->nIDsToBuild.end(), nID) != attr->nIDsToBuild.end()) {
            buildPos = i;
        }
    }
    VertexID firstID = nID;
    if (newAttrPos != -1) firstID = path[newAttrPos]->getBagsBelow()[0];
    std::vector<VertexID> bagOrder, bagToOrder(t.numNodes);
    std::vector<PrefixNode *> attributeOrder;
    newPrefix->getTraverseOrder(attributeOrder, bagOrder, t);
    for (ui i = 0; i < bagOrder.size(); ++i) bagToOrder[bagOrder[i]] = i;
    for (int i1 = bagToOrder[firstID]; i1 < bagOrder.size(); ++i1) {
        VertexID nID2 = bagOrder[i1];
        for (int i = 0; i < newPrefixes.size(); ++i) {
            VertexID u2 = newPrefixes[i];
            if (std::find(t.v2n[u2].begin(), t.v2n[u2].end(), nID2) != t.v2n[u2].end()) {
                bool prefixExists = false;
                for (int j = 0; j < t.nodes[nID2].prefixSize; ++j) {
                    if (t.nodes[nID2].attributes[j] == u2) {
                        prefixExists = true;
                        break;
                    }
                }
                if (!prefixExists) ++prefixSizes[nID2];
            }
        }
    }

    delete newPrefix;
    return prefixSizes;
}

void
resetMaterialize(const Graph &query, HyperTree &t, const HyperTree &defaultT, VertexID nID, int depth,
                 const std::vector<ui> &lengths, PrefixNode *&newPrefix, HyperTree *&newTree,
                 const PrefixNode *oldPT, CandidateSpace &cs, std::vector<bool> &fixed,
                 std::vector<int> &mappingSizes, std::vector<std::pair<int, VertexID>> &depend,
                 vector<vector<std::vector<VertexID>>> &tuples,
                 std::map<PrefixNode *, std::vector<VertexID>> &bagsBelow,
                 const std::vector<std::vector<int>> &attrIDMap, const std::vector<PrefixNode *> &dynamicPartition,
                 VertexID *partMatch) {
    newPrefix = oldPT -> clone();
    newTree = new HyperTree();
    // path[0] should be the virtual node
    std::vector<PrefixNode *> path = {newPrefix};
    std::vector<PrefixNode *> tmp = newPrefix->locate(nID);
    for (PrefixNode *attribute: tmp) path.push_back(attribute);
    ui newPrefixSize = lengths[nID];
    VertexID firstID = nID;
    std::vector<VertexID> newPrefixes(t.nodes[nID].attributes + t.nodes[nID].prefixSize, t.nodes[nID].attributes + newPrefixSize);
    int extendSize = depth;
    if (newPrefixSize <= depth) extendSize = newPrefixSize;
    else {
        bool includes = true;
        for (VertexID u: newPrefixes) {
            bool exists = false;
            for (int i = 0; i < extendSize; ++i) {
                if (defaultT.nodes[nID].attributes[i] == u) {
                    exists = true;
                    break;
                }
            }
            if (!exists) includes = false;
        }
        while (!includes) {
            ++extendSize;
            includes = true;
            for (VertexID u: newPrefixes) {
                bool exists = false;
                for (int i = 0; i < extendSize; ++i) {
                    if (defaultT.nodes[nID].attributes[i] == u) {
                        exists = true;
                        break;
                    }
                }
                if (!exists) includes = false;
            }
        }
        if (extendSize > depth) {
            PrefixNode *lastAttr = path[depth];
            PrefixNode *newAttr, *attr;
            std::vector<VertexID> oldNIDsToCall = lastAttr->nIDsToCall;
            lastAttr->nIDsToCall.clear();
            std::vector<int> ids(lastAttr->children.size());
            for (int j = 0; j < lastAttr->children.size(); ++j) {
                VertexID firstBag = lastAttr->children[j]->getBagsBelow()[0];
                ids[j] = attrIDMap[firstBag][depth];
            }
            // extend nID and sibling bags by one level
            int pos = oldNIDsToCall.size();
            for (int i = 0; i < oldNIDsToCall.size(); ++i) {
                VertexID nID2 = oldNIDsToCall[i];
                if (nID2 == t.numNodes - 1) {
                    lastAttr->nIDsToCall.push_back(nID2);
                    continue;
                }
                if (nID2 == nID) pos = i;
                if (i >= pos) {
                    VertexID u = defaultT.nodes[nID2].attributes[depth];
                    newAttr = new PrefixNode(u);
                    newAttr->pathToGlobal = false;
                    for (int j = 0; j < t.nodes[nID2].numAttributes; ++j) {
                        if (t.nodes[nID2].attributes[j] == u) {
                            newAttr->attributesBefore = t.nodes[nID2].attributesBefore[j];
                        }
                    }
                    lastAttr->children.push_back(newAttr);
                    ids.push_back(attrIDMap[nID2][depth]);
                    if (i == pos) {
                        attr = newAttr;
                        path.push_back(attr);
                    }
                    else newAttr->nIDsToCall.push_back(nID2);
                    if (nID2 != t.numNodes - 1 && i > pos && !fixed[u]) {
                        ++mappingSizes[nID2];
                        bool labelExists = false;
                        LabelID l1 = query.getVertexLabel(t.nodes[nID].attributes[depth]);
                        for (int j = depth; j < t.nodes[nID2].numAttributes; ++j) {
                            LabelID l2 = query.getVertexLabel(t.nodes[nID2].attributes[j]);
                            if (l1 == l2) labelExists = true;
                        }
                        if (labelExists) {
                            int id1 = attrIDMap[nID][depth], id2 = attrIDMap[nID2][depth];
                            depend[id2].first = id1;
                            depend[id2].second = -1;
                        }
                    }
                }
                else lastAttr->nIDsToCall.push_back(nID2);
            }
            for (int j1 = 0; j1 < lastAttr->children.size() - 1; ++j1) {
                for (int j2 = j1 + 1; j2 < lastAttr->children.size(); ++j2) {
                    if (ids[j2] < ids[j1]) {
                        std::swap(lastAttr->children[j1], lastAttr->children[j2]);
                        std::swap(ids[j1], ids[j2]);
                    }
                }
            }
            newAttr = attr;
            for (int i = depth + 1; i < extendSize; ++i) {
                VertexID u = defaultT.nodes[nID].attributes[i];
                newAttr = new PrefixNode(u);
                newAttr->pathToGlobal = false;
                for (int j = 0; j < t.nodes[nID].numAttributes; ++j) {
                    if (t.nodes[nID].attributes[j] == u) {
                        newAttr->attributesBefore = t.nodes[nID].attributesBefore[j];
                    }
                }
                attr->children.push_back(newAttr);
                attr = newAttr;
                path.push_back(newAttr);
            }
            newAttr->nIDsToCall = {nID};
        }
    }
    int newAttrPos = 0, buildPos = 0;
    for (int i = 0; i < path.size(); ++i) {
        PrefixNode *attr = path[i];
        if (attr->u == newPrefixes[0]) {
            newAttrPos = i;
        }
        if (std::find(attr->nIDsToBuild.begin(), attr->nIDsToBuild.end(), nID) != attr->nIDsToBuild.end()) {
            buildPos = i;
        }
    }
    firstID = path[newAttrPos]->getBagsBelow()[0];
    int lastPos = extendSize;
    while (std::find(newPrefixes.begin(), newPrefixes.end(), path[lastPos]->u) == newPrefixes.end()) --lastPos;
    std::set<VertexID> built;
    for (int i = lastPos; i >= buildPos; --i) {
        PrefixNode *attr = path[i];
        std::vector<VertexID> below;
        if (i != buildPos) below = attr->getBagsBelow();
        else below = attr->nIDsToBuild;
        attr -> nIDsToBuild.clear();
        for (VertexID bag: below) {
            if (bag == t.numNodes - 1) continue;
            if (built.find(bag) == built.end()) {
                attr -> nIDsToBuild.push_back(bag);
                built.insert(bag);
            }
        }
    }
    newPrefix->addBagsBelow(bagsBelow);
    std::vector<VertexID> bagOrder, bagToOrder(t.numNodes);
    std::vector<PrefixNode *> attributeOrder;
    newPrefix->getTraverseOrder(attributeOrder, bagOrder, t);
    for (ui i = 0; i < bagOrder.size(); ++i) bagToOrder[bagOrder[i]] = i;
    std::vector<VertexID> prevFixed;
    for (VertexID u = 0; u < query.getNumVertices(); ++u) {
        if (fixed[u] && std::find(t.v2n[u].begin(), t.v2n[u].end(), nID) == t.v2n[u].end())
            prevFixed.push_back(u);
    }
    for (int i = buildPos; i < extendSize; ++i) {
        PrefixNode *newAttr = path[i + 1];
        VertexID u = newAttr->u;
        if (i == -1) continue;
        if (i + 1 > newAttrPos) firstID = newAttr->getBagsBelow()[0];
        if (fixed[u]) continue;
        for (int j = 0; j < bagToOrder[firstID]; ++j) {
            VertexID nID2 = bagOrder[j];
            if (std::find(t.v2n[u].begin(), t.v2n[u].end(), nID2) != t.v2n[u].end()) {
                newAttr->nIDsToJoin.push_back(nID2);
                std::vector<VertexID> oldBefore = newAttr->attributesBefore;
                newAttr->attributesBefore.clear();
                for (VertexID u2: oldBefore) {
                    if (std::find(t.nodes[nID2].attributes, t.nodes[nID2].attributes + t.nodes[nID2].numAttributes,
                                  u) == t.nodes[nID2].attributes + t.nodes[nID2].numAttributes) {
                        newAttr->attributesBefore.push_back(u2);
                    }
                }
                if (newAttr->attributesBefore.empty())
                    newAttr->attributesBefore = oldBefore;
            }
        }
        for (VertexID u2: prevFixed) {
            if (query.getEdgeID(u, u2) != -1) {
                bool check = true;
                const std::vector<VertexID> &u2Nodes = t.v2n[u2];
                for (VertexID nID2: u2Nodes) {
                    if (bagToOrder[nID2] < bagToOrder[nID] &&
                        std::find(t.nodes[nID2].attributes, t.nodes[nID2].attributes + t.nodes[nID2].numAttributes,
                                  u) != t.nodes[nID2].attributes + t.nodes[nID2].numAttributes) check = false;
                }
                if (check)
                    newAttr->attributesBefore.push_back(u2);
            }
        }
    }
    for (VertexID u: newPrefixes) fixed[u] = true;
    int offset = t.nodes[nID].prefixSize;
    for (int i = 0; i < depth; ++i) {
        bool exists = false;
        VertexID u = t.nodes[nID].attributes[i];
        for (int j = 0; j <= depth; ++j) {
            if (u == path[j]->u) exists = true;
        }
        if (!exists) --offset;
    }
    int startPos1 = 0;
    for (int i = 0; i < attributeOrder.size(); ++i) {
        if (attributeOrder[i] == path.back()) startPos1 = i;
    }
    for (int i = startPos1 + 1; i < attributeOrder.size(); ++i) {
        VertexID u = attributeOrder[i]->u;
        for (int j = 0; j < newPrefixes.size(); ++j) {
            VertexID u2 = newPrefixes[j];
            int u2ID = attrIDMap[nID][j + offset];
            if (u2 == u) {
                attributeOrder[i]->nIDsToJoin.clear();
            }
            else {
                if (query.getEdgeID(u, u2) != -1) {
                    bool newAttrBefore = true;
                    if (std::find(attributeOrder[i]->attributesBefore.begin(), attributeOrder[i]->attributesBefore.end(),
                                  u2) != attributeOrder[i]->attributesBefore.end()) newAttrBefore = false;
                    for (VertexID nID2: attributeOrder[i]->nIDsToJoin) {
                        if (std::find(t.v2n[u2].begin(), t.v2n[u2].end(), nID2) != t.v2n[u2].end()) {
                            newAttrBefore = false;
                            break;
                        }
                    }
                    if (newAttrBefore) {
                        attributeOrder[i]->attributesBefore.push_back(u2);
//                        int uID = 0;
//                        VertexID uNID = bagsBelow[attributeOrder[i]][0];
//                        for (int k = 0; k < t.nodes[uNID].numAttributes; ++k) {
//                            if (t.nodes[uNID].attributes[k] == u) {
//                                uID = attrIDMap[uNID][k];
//                                break;
//                            }
//                        }
//                        if (u2ID > depend[uID].first) {
//                            depend[uID].first = u2ID;
//                            depend[uID].second = -1;
//                        }
                    }
                }
            }
        }
    }
    std::vector<ui> prefixSizes(t.numNodes);
    std::vector<std::vector<VertexID>> nodeOrders(t.numNodes);
    for (VertexID nID2 = 0; nID2 < t.numNodes; ++nID2) {
        prefixSizes[nID2] = t.nodes[nID2].prefixSize;
        nodeOrders[nID2].assign(t.nodes[nID2].attributes, t.nodes[nID2].attributes + t.nodes[nID2].numAttributes);
    }
    for (int i = 0; i < newPrefixes.size(); ++i) {
        firstID = path[newAttrPos + i]->getBagsBelow()[0];
        VertexID u2 = newPrefixes[i];
        LabelID l2 = query.getVertexLabel(u2);
        for (int i1 = bagToOrder[firstID]; i1 < bagOrder.size(); ++i1) {
            VertexID nID2 = bagOrder[i1];
            int oldMappingSize = mappingSizes[nID2];
            int u2ID = attrIDMap[nID][i + offset];
            int dependAttr = 0;
            bool labelExists = false;
            for (int j = 0; j < t.nodes[nID2].numAttributes; ++j) {
                if (attrIDMap[nID2][j] > u2ID) {
                    dependAttr = attrIDMap[nID2][j];
                    for (int k = j; k < t.nodes[nID2].numAttributes; ++k) {
                        if (query.getVertexLabel(t.nodes[nID2].attributes[k]) == l2)
                            labelExists = true;
                    }
                    break;
                }
            }
            if (nID2 != t.numNodes - 1 && labelExists) {
                depend[dependAttr].first = u2ID;
                depend[dependAttr].second = -1;
            }
            if (std::find(t.v2n[u2].begin(), t.v2n[u2].end(), nID2) != t.v2n[u2].end()) {
                bool prefixExists = false, sharedExists = false;
                for (int j = 0; j < oldMappingSize; ++j) {
                    if (nodeOrders[nID2][j] == u2) {
                        sharedExists = true;
                        break;
                    }
                }
                for (int j = 0; j < t.nodes[nID2].prefixSize; ++j) {
                    if (t.nodes[nID2].attributes[j] == u2) {
                        prefixExists = true;
                        break;
                    }
                }
                std::vector<VertexID> prefixOrder;
                prefixOrder.assign(nodeOrders[nID2].begin(), nodeOrders[nID2].begin() + prefixSizes[nID2]);
                std::vector<VertexID> sortedPrefixOrder;
                if (!prefixExists) ++prefixSizes[nID2];
                if (nID != nID2) {
                    std::vector<int> prefixAttrID = getPrefixAttrID(const_cast<PrefixNode *>(oldPT), prefixOrder,
                                                                    attrIDMap, dynamicPartition);
                    if (!prefixExists) {
                        prefixOrder.push_back(u2);
                        prefixAttrID.push_back(u2ID);
                    }
                    else {
                        for (int j = 0; j < prefixSizes[nID2]; ++j) {
                            if (prefixOrder[j] == u2) {
                                prefixAttrID[j] = u2ID;
                            }
                        }
                    }
                    sortedPrefixOrder = std::vector<VertexID>(prefixOrder.size());
                    std::vector<size_t> indices(prefixOrder.size());
                    for (size_t j = 0; j < indices.size(); ++j) indices[j] = j;
                    std::sort(indices.begin(), indices.end(), [&prefixAttrID](size_t i, size_t j) {
                        return prefixAttrID[i] < prefixAttrID[j];
                    });
                    for (size_t j = 0; j < indices.size(); ++j) {
                        sortedPrefixOrder[j] = prefixOrder[indices[j]];
                    }
                }
                else {
                    sortedPrefixOrder = prefixOrder;
                    sortedPrefixOrder.push_back(u2);
                }
                const std::vector<VertexID> &oldNodeOrder = nodeOrders[nID2];
                std::vector<VertexID> newNodeOrder = sortedPrefixOrder;
                for (VertexID u: oldNodeOrder) {
                    if (std::find(sortedPrefixOrder.begin(), sortedPrefixOrder.end(), u) == sortedPrefixOrder.end())
                        newNodeOrder.push_back(u);
                }
                nodeOrders[nID2] = newNodeOrder;
                if (!sharedExists) ++mappingSizes[nID2];
            }
        }
    }
    // build a new hypertree to store the orders
    newTree->newGlobalNode = t.newGlobalNode;
    newTree->numNodes = t.numNodes;
    newTree->nodes = new HyperNode[newTree->numNodes];
    newTree->numAttributes = t.numAttributes;
    for (VertexID nID2 = 0; nID2 < t.numNodes; ++nID2) {
        t.nodes[nID2].copyTo(newTree->nodes[nID2]);
        newTree->nodes[nID2].prefixSize = prefixSizes[nID2];
        if (prefixSizes[nID2] != 0 && t.nodes[nID2].prefixSize == 0) newTree->nodes[nID2].prefix = new VertexID[t.numAttributes];
        for (int i = 0; i < prefixSizes[nID2]; ++i)
            newTree->nodes[nID2].prefix[i] = nodeOrders[nID2][i];
        for (int i = 0; i < nodeOrders[nID2].size(); ++i)
            newTree->nodes[nID2].attributes[i] = nodeOrders[nID2][i];
    }
    newTree->defaultPartition.clear();
    for (VertexID nID2: bagOrder) {
        for (int i = 0; i < prefixSizes[nID2]; ++i) {
            VertexID u = newTree->nodes[nID2].attributes[i];
            if (std::find(newTree->defaultPartition.begin(), newTree->defaultPartition.end(), u) == newTree->defaultPartition.end())
                newTree->defaultPartition.push_back(u);
        }
    }
    const HyperNode &last = newTree->nodes[newTree->numNodes - 1];
    for (int i = 0; i < last.numAttributes; ++i) {
        VertexID u = last.attributes[i];
        if (std::find(newTree->defaultPartition.begin(), newTree->defaultPartition.end(), u) == newTree->defaultPartition.end())
            newTree->defaultPartition.push_back(u);
    }
    newTree->v2n = new std::vector<VertexID> [newTree->numAttributes];
    for (ui i = 0; i < newTree->numNodes; ++i) {
        for (ui j = 0; j < newTree->nodes[i].numAttributes; ++j) {
            newTree->v2n[newTree->nodes[i].attributes[j]].push_back(i);
        }
    }
    std::vector<std::vector<VertexID>> sharedAttrs(t.numNodes);
    for (VertexID nID2 = 0; nID2 < newTree->numNodes; ++nID2) {
        for (int i = 0; i < newTree->nodes[nID2].numAttributes; ++i) {
            VertexID u = newTree->nodes[nID2].attributes[i];
            if (newTree->v2n[u].size() > 1) {
                sharedAttrs[nID2].push_back(u);
            }
        }
        std::sort(sharedAttrs[nID2].begin(), sharedAttrs[nID2].end());
    }
    for (VertexID nID2 = 0; nID2 < t.numNodes; ++nID2) {
        newTree->nodes[nID2].initPoses(sharedAttrs, query, cs.dist, nID2 == t.numNodes - 1);
    }
    newTree->globalOrder = newTree->defaultPartition;
    for (VertexID nID2 = 0; nID2 < t.numNodes; ++nID2) {
        for (int i = prefixSizes[nID2]; i < newTree->nodes[nID2].numAttributes; ++i) {
            VertexID u = newTree->nodes[nID2].attributes[i];
            if (std::find(newTree->globalOrder.begin(), newTree->globalOrder.end(), u) == newTree->globalOrder.end())
                newTree->globalOrder.push_back(u);
        }
    }
    newTree->initPoses(query, cs, false);
//    if (newTree->extendLevel > t.nodes[t.numNodes - 1].numAttributes)
//        newTree->extendLevel = t.nodes[t.numNodes - 1].numAttributes;
//    setTDExtention(*newTree, query);
    newTree->extendLevel = t.nodes[t.numNodes - 1].numAttributes;
    newTree->buildTraverseStruct(query);
    ++gNumFixedCase;
    int numFixedAttr = 0;
    for (VertexID u = 0; u < query.getNumVertices(); ++u) {
        if (fixed[u]) ++numFixedAttr;
    }
    if (numFixedAttr > gMaxNewFixed) gMaxNewFixed = numFixedAttr;
}

void
resetPoses(const HyperTree &t, const HyperTree &defaultT, VertexID nID, VertexID **nodeCandidates, ui *nodeCandCount,
           std::vector<ui> &nodePoses, VertexID **pCandidates, ui *pCandCount, std::vector<ui> &pPoses, ui newLength,
           ui pathLength, const std::vector<VertexID> &resetTuple, bool *visited, CandidateSpace &cs,
           VertexID *partMatch, const std::vector<std::vector<int>> &attrIDMap) {
    bool recompute = false;
    VertexID **neighbors = new VertexID *[t.numAttributes];
    ui *neighborCount = new ui[t.numAttributes];
    int attrID;
    for (int i = 0; i < newLength - t.nodes[nID].prefixSize; ++i) {
        ui mappingSize = t.nodes[nID].prefixSize + i;
        VertexID u = t.nodes[nID].attributes[mappingSize];
        VertexID v = resetTuple[i];
        partMatch[u] = v;
        visited[v] = true;
        int originalPos = mappingSize;
        for (int j = 0; j < defaultT.nodes[nID].numAttributes; ++j) {
            if (defaultT.nodes[nID].attributes[j] == u) {
                attrID = attrIDMap[nID][j];
                originalPos = j;
                break;
            }
        }
        VertexID *nCandidate = nodeCandidates[mappingSize], *pCandidate = pCandidates[attrID];
        // copy nCandidate to pCandidate if exceeds pathLength
        if (originalPos >= pathLength) {
            memcpy(pCandidates[attrID], nCandidate, sizeof(VertexID) * nodeCandCount[mappingSize]);
            pCandCount[attrID] = nodeCandCount[mappingSize];
        }
        if (!recompute) {
            pPoses[attrID] = std::lower_bound(pCandidate, pCandidate + pCandCount[attrID], v) - pCandidate;
            recompute = true;
        }
        else {
            const std::vector<VertexID> &parents = t.nodes[nID].attributesBefore[mappingSize];
            if (!parents.empty()) {
                for (int j = 0; j < parents.size(); ++j) {
                    VertexID pU = parents[j];
                    neighbors[j] = cs.candidateEdge[pU][partMatch[pU]][u].data();
                    neighborCount[j] = cs.candidateEdge[pU][partMatch[pU]][u].size();
                }
                ComputeSetIntersection::LeapfrogJoin(neighbors, neighborCount, parents.size(), pCandidate, pCandCount[attrID]);
            }
            else {
                VertexID cartesianParent = t.nodes[nID].cartesianParent[mappingSize];
                const std::vector<VertexID> &path = cs.reconstructPath(cartesianParent, u);
                optimizedCartesianProduct(cs, partMatch[cartesianParent], path, pCandidate, pCandCount[attrID]);
            }
            pPoses[attrID] = std::lower_bound(pCandidate, pCandidate + pCandCount[attrID], v) - pCandidate;
        }
    }
    delete[] neighbors;
    delete[] neighborCount;
}

void updateSkipMatch(std::vector<PrefixNode *> &path, int depth, ui newLength, std::vector<bool> &skipMatch,
                     const std::vector<std::vector<int>> &attrIDMap,
                     map<PrefixNode *, vector<VertexID>> &bagsBelow,
                     std::vector<PrefixNode *> &dynamicPartition) {
    PrefixNode *current = path[depth];
    for (PrefixNode *c : current -> children) {
        if (bagsBelow[c].size() > 1) continue;
        VertexID nID = bagsBelow[c][0];
        int id = attrIDMap[nID][depth];
        bool exists = false;
        for (PrefixNode *attr: dynamicPartition) {
            if (attr->u == c->u) exists = true;
        }
        if (exists) skipMatch[id] = true;
    }
    VertexID nID = bagsBelow[path[depth + 1]][0];
    for (int i = depth + 1; i <= newLength; ++i) {
        int id = attrIDMap[nID][i - 1];
        bool exists = false;
        for (PrefixNode *attr: dynamicPartition) {
            if (attr->u == path[i]->u) exists = true;
        }
        if (exists) skipMatch[id] = true;
    }
}

void makeBackup(ui oldLength, ui newLength, int depth, const HyperTree &newT,
                std::vector<std::vector<std::vector<VertexID>>> &tuples,
                std::vector<std::vector<std::vector<VertexID>>> &backup, const vector<PrefixNode *> &path, VertexID nID,
                VertexID *partMatch, const std::vector<VertexID> &bagToOrder, size_t &budget) {
    std::vector<VertexID> newPrefix;
    int offset = oldLength;
    for (int i = 0; i < oldLength; ++i) {
        bool exists = false;
        VertexID u = newT.nodes[nID].attributes[i];
        for (int j = 0; j <= depth; ++j) {
            if (u == path[j]->u) exists = true;
        }
        if (!exists) --offset;
    }
    std::vector<std::vector<VertexID>> newPrefixes;
    std::vector<VertexID> nIDs;

    for (int i = 0; i < newLength - oldLength; ++i) {
        ui mappingSize = offset + i + 1;
        if (mappingSize <= depth) {
            PrefixNode *attribute = path[mappingSize];
            VertexID u = attribute->u;
            VertexID v = partMatch[u];
            newPrefix.push_back(v);
            for (VertexID nID2: attribute->getBagsBelow()) {
                if (newT.nodes[nID2].prefixSize >= mappingSize) {
                    if (nID2 != nID) {
                        for (int j = 0; j < nIDs.size(); ++j) {
                            if (nIDs[j] == nID2) {
                                newPrefixes.erase(newPrefixes.begin() + j);
                                nIDs.erase(nIDs.begin() + j);
                                break;
                            }
                        }
                        nIDs.push_back(nID2);
                        newPrefixes.push_back(newPrefix);
                    }
                }
            }
        }
        else break;
    }
    for (int i = 0; i < nIDs.size(); ++i) {
        VertexID nID2 = nIDs[i];
        newPrefix = newPrefixes[i];
        int pos1 = findFirstExtension(tuples[nID2], newPrefix);
        int pos2 = findFirstGreater(tuples[nID2], newPrefix);
        budget += (pos2 - pos1) * newPrefix.size() * sizeof(VertexID);
        if (!backup[nID2].empty())
            budget += backup[nID2].size() * backup[nID2][0].size() * sizeof(VertexID);
        backup[nID2].assign(tuples[nID2].begin() + pos2, tuples[nID2].end());
        tuples[nID2].erase(tuples[nID2].begin() + pos2, tuples[nID2].end());
        if (!tuples[nID2].empty())
            budget += tuples[nID2][0].size() * pos1 * sizeof(VertexID);
        tuples[nID2].erase(tuples[nID2].begin(), tuples[nID2].begin() + pos1);
        for (auto &tp : tuples[nID2]) {
            tp.assign(tp.begin() + newPrefix.size(), tp.end());
        }
    }
}

bool buildFromExist(const HyperTree &t, VertexID nID, VertexID *partMatch, std::vector<std::vector<VertexID>> &tuples,
                    std::vector<std::vector<VertexID>> &backup, size_t &budget, std::vector<TrieNode *> &trieNodes,
                    int pathLength) {
    if (backup.empty()) return false;
    const HyperNode &bag = t.nodes[nID];
    int oldSize = backup[0].size(), newSize = bag.numAttributes - bag.prefixSize;
    if (oldSize == newSize) return false;
    int oldPrefixSize = bag.numAttributes - oldSize;
    std::vector<VertexID> newPrefix;
    for (int i = oldPrefixSize; i < pathLength; ++i) {
        newPrefix.push_back(partMatch[bag.attributes[i]]);
    }
    bool flag = false;
    int pos1 = findFirstExtension(backup, newPrefix);
    int pos2 = findFirstGreater(backup, newPrefix);
    if (pos1 != backup.size() && pos1 != pos2) {
        budget += oldSize * pos2 * sizeof(VertexID);
        flag = true;
        for (int i = pos1; i < pos2; ++i) {
            const std::vector<VertexID> &tp = backup[i];
            tuples.emplace_back(tp.end() - newSize, tp.end());
        }
        budget -= newSize * (pos2 - pos1) * sizeof(VertexID);
        backup.assign(backup.begin() + pos2, backup.end());
    }
    else {
        budget += oldSize * backup.size() * sizeof(VertexID);
        backup.clear();
    }
    return flag;
}

void attributeOpen(std::vector<bool> &skipMatch, CandidateSpace &cs, int attrID, int depth, PrefixNode *current,
                   std::vector<TrieNode *> &lastNodes,
                   std::vector<std::vector<DynamicArray<TrieNode *> *>> &edgeColumns, ui ***pIters, ui *pIterSizes,
                   VertexID **pCandidates, ui *pCandCount, std::vector<ui> &pPoses, DynamicArray<TrieNode *> **children,
                   VertexID *partMatch, std::vector<bool> &currentFixed, VertexID **neighbors, ui *neighborCount,
                   bool parentNewVersion) {
    VertexID u = current->u;
    if (current->pathToGlobal) {
        if (skipMatch[attrID]) {
            pPoses[attrID] = -1;
            pIterSizes[attrID] = 1;
            bool checkBefore = true;
            for (VertexID u2: current->attributesBefore) {
                if (!cs.checkExists(partMatch[u], partMatch[u2])) checkBefore = false;
            }
            if (checkBefore) {
                for (int i = 0; i < current->nIDsToJoin.size(); ++i) {
                    VertexID nID = current->nIDsToJoin[i];
                    children[i] = &lastNodes[nID]->nodeChild;
                    pIters[attrID][i][0] = binarySearch(*children[i], partMatch[u]);
                }
                if (current->nIDsToJoin.empty()) {
                    pIters[attrID][0][0] = 0;
                    edgeColumns[depth][0][0][0]->value = partMatch[u];
                }
            }
            else pIterSizes[attrID] = 0;
        }
        else {
            pPoses[attrID] = -1;
            if (current -> cartesianParent != 99) {
                VertexID pU = current -> cartesianParent;
                const std::vector<VertexID> &path = cs.reconstructPath(pU, u);
                optimizedCartesianProduct(cs, partMatch[pU], path, edgeColumns[depth][0]);
                pIterSizes[attrID] = edgeColumns[depth][0]->size();
                for (int i = 0; i < edgeColumns[depth][0]->size(); ++i) {
                    pIters[attrID][i][0] = i;
                }
            }
            else if (current ->nIDsToJoin.empty() && current -> attributesBefore.empty()) {
                edgeColumns[depth][0] ->setSize(cs.candidateSet[u].size());
                for (int i = 0; i < cs.candidateSet[u].size(); ++i) {
                    edgeColumns[depth][0][0][i]->value = cs.candidateSet[u][i];
                    pIters[attrID][i][0] = i;
                }
                pIterSizes[attrID] = cs.candidateSet[u].size();
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
                leapFrogJoin(children, current->nIDsToJoin.size() + current->attributesBefore.size(), pIters[attrID], pIterSizes[attrID]);
        }
    }
    else {
        if (!skipMatch[attrID]) {
            if (currentFixed[u] && parentNewVersion) {

            }
            else if (depth != 1 || !current->attributesBefore.empty()) {
                pPoses[attrID] = -1;
                const std::vector<VertexID> &parents = current->attributesBefore;
                if (!parents.empty()) {
                    for (int i = 0; i < parents.size(); ++i) {
                        VertexID pU = parents[i];
                        neighbors[i] = cs.candidateEdge[pU][partMatch[pU]][u].data();
                        neighborCount[i] = cs.candidateEdge[pU][partMatch[pU]][u].size();
                    }
                    ComputeSetIntersection::LeapfrogJoin(neighbors, neighborCount, parents.size(), pCandidates[attrID], pCandCount[attrID]);
                }
                else {
                    VertexID cartesianParent = current->cartesianParent;
                    const std::vector<VertexID> &path = cs.reconstructPath(cartesianParent, u);
                    optimizedCartesianProduct(cs, partMatch[cartesianParent], path, pCandidates[attrID], pCandCount[attrID]);
                }
            }
            else {
                pPoses[attrID] = -1;
                if (pCandCount[attrID] != cs.candidateSet[u].size()) {
                    pCandCount[attrID] = cs.candidateSet[u].size();
                    memcpy(pCandidates[attrID], cs.candidateSet[u].data(), sizeof(VertexID) * pCandCount[attrID]);
                }
            }
        }
        else {
            pCandidates[attrID][0] = partMatch[u];
            pPoses[attrID] = -1;
            pCandCount[attrID] = 1;
        }
    }
}

void jump(const std::vector<PrefixNode *> &dynamicPartition, PrefixNode *&current, PrefixNode *&pn,
          std::vector<ui> &childPoses, std::vector<int> &ids, vector<PrefixNode *> &path,
          std::vector<bool> &nextCands, int &depth, const std::vector<std::vector<int>> &attrIDMap, PrefixNode *tmp) {
    // current becomes the last attribute in dynamicPartition. change current, depth, path, childPose
    PrefixNode *next = dynamicPartition.back();
    std::vector<PrefixNode *> newPath;
    std::vector<ui> newChildPose;
    current->locate(next, newPath, newChildPose);
    VertexID firstBag = next->getBagsBelow()[0];
    for (int i = 0; i < newPath.size(); ++i) {
        ++depth;
        path[depth] = newPath[i];
        childPoses[depth] = newChildPose[i];
        ids[depth] = attrIDMap[firstBag][depth - 1];
    }
    current = next;
    while (path[depth] != next) {
        --depth;
        if (depth == -1) break;
    }
    if (depth == -1) {
        newPath.clear();
        newChildPose.clear();
        path[0]->locate(next, newPath, newChildPose);
        depth = 0;
        for (int i = 0; i < newPath.size(); ++i) {
            ++depth;
            path[depth] = newPath[i];
            childPoses[depth] = newChildPose[i];
            ids[depth] = attrIDMap[firstBag][depth - 1];
        }
    }
    if (depth == 0) pn = tmp;
    else pn = path[depth - 1];
    nextCands[ids[depth]] = true;
}

void rebuildTrie(const HyperTree &oldT, const HyperTree &newT, std::vector<TrieNode *> &lastNodes,
                 std::vector<PrefixNode *> &path, ui extendSize, std::vector<bool> &skipBuild,
                 std::vector<std::vector<std::vector<VertexID>>> &tuples, size_t &budget,
                 std::vector<size_t> &tupleSizes, const std::vector<bool> &oldFixed) {
    std::vector<bool> rebuilt(tuples.size(), false);
    for (ui i = 0; i < extendSize; ++i) {
        PrefixNode *current = path[i + 1];
        if (oldFixed[current->u]) continue;
        for (VertexID nID: current->nIDsToJoin) {
            if (rebuilt[nID]) continue;
            int level = 0, offset;
            std::vector<std::vector<VertexID>> bagTuples;
            if (!skipBuild[nID]) {
                level = oldT.trieOrder[nID].size();
                offset = 0;
            }
            else {
                TrieNode *tau = lastNodes[nID];
                while (!tau->nodeChild.empty()) {
                    tau = tau->nodeChild[0];
                    ++level;
                }
                if (level == 0) continue;
                offset = oldT.trieOrder[nID].size() - level;
            }
            std::vector<VertexID> prevOrder, newOrder;
            for (int j = 0; j < level; ++j) {
                prevOrder.push_back(oldT.trieOrder[nID][j + offset]);
                newOrder.push_back(newT.trieOrder[nID][j + offset]);
            }
            makeContinuous(prevOrder);
            makeContinuous(newOrder);
            if (!skipBuild[nID]) {
                skipBuild[nID] = true;
//                if (tuples[nID].empty()) continue;
                budget += tupleSizes[nID];
                tupleSizes[nID] = tuples[nID].size() * tuples[nID][0].size() * sizeof(VertexID);
                buildTrie(tuples[nID], lastNodes[nID], newOrder);
                rebuilt[nID] = true;
                continue;
            }
            // traverse to reconstruct tuples
            int depth = 0;
            std::vector<int> poses(level, 0);
            std::vector<ui> candCount(level, 0);
            candCount[0] = lastNodes[nID]->nodeChild.size();
            std::vector<VertexID> tuple(level);
            std::vector<TrieNode *> triePath(level);
            while (depth >= 0) {
                TrieNode *parent;
                if (depth != 0) parent = triePath[depth - 1];
                else parent = lastNodes[nID];
                while (poses[depth] < candCount[depth]) {
                    TrieNode *tn = parent->nodeChild[poses[depth]];
                    triePath[depth] = tn;
                    tuple[depth] = tn->value;
                    ++poses[depth];
                    if (depth < level - 1) {
                        ++depth;
                        parent = tn;
                        poses[depth] = 0;
                        candCount[depth] = tn->nodeChild.size();
                    }
                    else bagTuples.push_back(tuple);
                }
                --depth;
            }
            // reorder tuples as local orders
            for (int j1 = 0; j1 < bagTuples.size(); ++j1) {
                std::vector<VertexID> tmp = bagTuples[j1];
                for (int j2 = 0; j2 < prevOrder.size(); ++j2) {
                    bagTuples[j1][prevOrder[j2]] = tmp[j2];
                }
            }
            // rebuild new trie
            buildTrie(bagTuples, lastNodes[nID], newOrder);
            rebuilt[nID] = true;
        }
    }
}

bool moveToVertex(PrefixNode *current, std::vector<TrieNode *> &lastNodes,
                  std::vector<std::vector<TrieNode *>> &traversedNodes, TrieNode *empty, VertexID v) {
    bool flag = false;
    for (VertexID nID: current->nIDsToJoin) {
        int pos = binarySearch(lastNodes[nID]->nodeChild, v);
        if (pos != -1) {
            traversedNodes[nID].push_back(lastNodes[nID]->nodeChild[pos]);
            lastNodes[nID] = traversedNodes[nID].back();
        }
        else {
            flag = true;
            traversedNodes[nID].push_back(empty);
            lastNodes[nID] = empty;
        }
    }

    return flag;
}

bool
moveTrie(const HyperTree &t, std::vector<std::vector<TrieNode *>> &traversedNodes, std::vector<TrieNode *> &lastNodes,
         PrefixNode *attribute, VertexID *partMatch, TrieNode *empty, const std::vector<std::vector<int>> &attrIDMap,
         map<PrefixNode *, vector<VertexID>> &bagsBelow, std::vector<bool> &skipMatch, int pathLength) {
    // dfs traverse subtree rooted at current
    // for each attribute, extend related tries based on partMatch
    PrefixNode *pn = attribute;
    ui height = pn -> getHeight();
    std::vector<ui> childPoses(height, 0);
    std::vector<PrefixNode *> nodes(height, nullptr);
    int depth = 0;
    bool deadEnd = false;
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            ++childPoses[depth];
            if (current->pathToGlobal) continue;
            nodes[depth] = current;
            VertexID v = partMatch[current->u];
            if (moveToVertex(current, lastNodes, traversedNodes, empty, v)) deadEnd = true;
            if (!current -> children.empty()) {
                ++depth;
                childPoses[depth] = 0;
            }
            if (depth > 0) pn = nodes[depth - 1];
            else pn = attribute;
        }
        --depth;
        if (depth >= 0) {
            if (depth == 0) pn = attribute;
            else pn = nodes[depth - 1];
        }
    }

    return deadEnd;
}

void adaptiveShareJoin(const Graph &query, HyperTree &t, PrefixNode *pt, CandidateSpace &cs,
                       std::vector<TrieNode *> &trieNodes, bool *visited, std::vector<std::vector<VertexID>> &result,
                       size_t &count, bool traverse, size_t budget, bool skip) {
    size_t oldBudget = budget;
    std::unordered_map<std::vector<ui>, PrefixNode *, VectorHash, VectorEqual> lengthToPT;
    std::unordered_map<std::vector<ui>, HyperTree *, VectorHash, VectorEqual> lengthToTree;
    std::unordered_map<std::vector<ui>, std::vector<int>, VectorHash, VectorEqual> lengthToMappingSizes;
    std::unordered_map<std::vector<ui>, std::vector<bool>, VectorHash, VectorEqual> lengthToFixed;
    std::unordered_map<std::vector<ui>, std::vector<std::pair<int, VertexID>>, VectorHash, VectorEqual> lengthToDepend;
    std::vector<std::vector<std::vector<VertexID>>> tuples(t.numNodes);
    std::vector<std::vector<std::vector<VertexID>>> backup(t.numNodes);
    VertexID *partMatch = new VertexID[t.numAttributes];
    TrieNode *empty = new TrieNode();
    ui maxSize = cs.getMaxSize();
    ui height = t.numAttributes + 1;
    std::vector<PrefixNode *> attributeOrder;
    std::map<PrefixNode *, std::vector<VertexID>> bagsBelow;
    PrefixNode *attributeTree = fullAttributeTree(pt, attributeOrder, bagsBelow, t);
//    pt->addBagsBelow(bagsBelow);
    std::vector<VertexID> bagOrder, bagToOrder(t.numNodes);
    std::vector<PrefixNode *> vec1;
    pt->getTraverseOrder(vec1, bagOrder, t);
    for (int i = 0; i < bagOrder.size(); ++i) bagToOrder[bagOrder[i]] = i;
    std::vector<std::vector<int>> attrIDMap = buildAttrIDMap(attributeOrder, bagsBelow, t);
    std::vector<std::vector<int>> uToAttrID(query.getNumVertices());
    for (int i = 0; i < attributeOrder.size(); ++i) {
        PrefixNode *attribute = attributeOrder[i];
        if (i > 0) uToAttrID[attribute->u].push_back(i);
    }
    VertexID **pCandidates = new VertexID *[attributeOrder.size()];
    ui *pCandCount = new ui[attributeOrder.size()];
    for (int i = 0; i < attributeOrder.size(); ++i) {
        PrefixNode *attribute = attributeOrder[i];
        if (!attribute -> pathToGlobal) {
            pCandidates[i] = new VertexID[cs.candidateSet[attribute->u].size()];
            pCandCount[i] = 0;
        }
        else pCandidates[i] = nullptr;
    }
    std::vector<ui> pPoses(attributeOrder.size(), 0);
    pPoses[0] = -1;
    std::vector<ui> childPoses(height, 0);
    VertexID **neighbors = new VertexID *[height];
    ui *neighborCount = new ui[height];
    DynamicArray<TrieNode *> ** children = new DynamicArray<TrieNode *> *[t.numAttributes];
    std::vector<std::vector<TrieNode *>> traversedNodes(trieNodes.size());
    std::vector<TrieNode *> lastNodes(trieNodes.size());
    for (VertexID nID = 0; nID < trieNodes.size(); ++nID) {
        traversedNodes[nID].push_back(trieNodes[nID]);
        lastNodes[nID] = trieNodes[nID];
    }
    VertexID ***nodeCandidates = new VertexID **[trieNodes.size()];
    ui **nodeCandCount = new ui *[trieNodes.size()];
    std::vector<std::vector<ui>> nodePoses(trieNodes.size());
//    for (VertexID nID = 0; nID < trieNodes.size(); ++nID) {
//        const HyperNode &tau = t.nodes[nID];
//        nodeCandidates[nID]= new VertexID *[tau.numAttributes];
//        nodeCandCount[nID] = new ui[tau.numAttributes];
//        nodePoses[nID] = std::vector<ui>(tau.numAttributes, 0);
//        for (int i = 0; i < tau.numAttributes; ++i) {
//            int attrID = attrIDMap[nID][i];
//            nodeCandidates[nID][i] = pCandidates[attrID];
//            nodeCandCount[nID][i] = 0;
//        }
//    }
    for (VertexID nID = 0; nID < trieNodes.size(); ++nID) {
        const HyperNode &tau = t.nodes[nID];
        nodeCandidates[nID]= new VertexID *[tau.numAttributes];
        nodeCandCount[nID] = new ui[tau.numAttributes];
        nodePoses[nID] = std::vector<ui>(tau.numAttributes, 0);
        ui mSize = 0;
        for (int i = 0; i < tau.numAttributes; ++i) {
            if (mSize < cs.candidateSet[tau.attributes[i]].size())
                mSize = cs.candidateSet[tau.attributes[i]].size();
        }
        for (int i = 0; i < tau.numAttributes; ++i) {
            nodeCandidates[nID][i] = new VertexID [mSize];
            nodeCandCount[nID][i] = 0;
        }
    }
    std::vector<VertexID> defaultLengths = std::vector<ui>(t.numNodes, 0);
    for (VertexID nID = 0; nID < t.numNodes; ++nID) defaultLengths[nID] = t.nodes[nID].prefixSize;
    ui ***nIters = new ui **[t.nodes[t.numNodes - 1].numAttributes + 1];
    ui *nIterSizes = new ui [t.nodes[t.numNodes - 1].numAttributes + 1];
    for (int i = 0; i < t.nodes[t.numNodes - 1].numAttributes + 1; ++i) {
        nIters[i] = new ui *[maxSize];
        for (int j = 0; j < maxSize; ++j) {
            nIters[i][j] = new ui[t.numAttributes];
        }
        nIterSizes[i] = 0;
    }
    ui ***pIters = new ui **[attributeOrder.size()];
    ui *pIterSizes = new ui[attributeOrder.size()];
    pIterSizes[0] = 1;
    pIters[0] = nullptr;
    for (int i = 1; i < attributeOrder.size(); ++i) {
        PrefixNode *attribute = attributeOrder[i];
        if (attribute -> pathToGlobal) {
            pIters[i] = new ui *[cs.candidateSet[attribute->u].size()];
            for (int j = 0; j < cs.candidateSet[attribute->u].size(); ++j) {
                pIters[i][j] = new ui[t.numAttributes];
            }
            pIterSizes[i] = 0;
        }
        else pIters[i] = nullptr;
    }
    ui lastID = t.numNodes - 1;
    ui maxEdgeSize = query.getMaxDegree();
    // join relational columns together with edge columns
    std::vector<std::vector<DynamicArray<TrieNode *> *>> edgeColumns(t.nodes[lastID].numAttributes + 1);
    for (int i = 0 ; i < t.nodes[lastID].numAttributes + 1; ++i) {
        edgeColumns[i] = std::vector<DynamicArray<TrieNode *> *>(maxEdgeSize);
        for (int j = 0; j < maxEdgeSize; ++j) {
            edgeColumns[i][j] = new DynamicArray<TrieNode *>();
            for (int k = 0; k < maxSize; ++k) {
                TrieNode * pointer = new TrieNode();
                edgeColumns[i][j]->push_back(pointer);
            }
        }
    }
    std::vector<bool> defaultFixed(query.getNumVertices(), false);
    std::vector<bool> skipMatch(attributeOrder.size(), false);
    std::vector<bool> changeMaterialize(attributeOrder.size(), false);
    std::vector<bool> skipBuild(t.numNodes, false);
    std::vector<bool> nextCands(attributeOrder.size(), true);
    std::vector<size_t> tupleSizes(t.numNodes, 0);
    std::vector<PrefixNode *> path(height, nullptr);
    std::vector<int> ids(height);
    std::vector<HyperTree *> trees(query.getNumVertices());
    std::vector<std::vector<int>> mappingSizes(query.getNumVertices());
    std::vector<std::vector<bool>> fixedAttrs(query.getNumVertices(), defaultFixed);
    std::vector<std::vector<ui>> prefixLengths(query.getNumVertices(), defaultLengths);
//    std::vector<std::vector<std::pair<int, VertexID>>> depends(query.getNumVertices());
    std::vector<VertexID> versionToNID(query.getNumVertices(), 0);
    std::vector<std::pair<int, VertexID>> depend = buildDependent(attributeOrder);
    lengthToDepend[defaultLengths] = depend;
    std::vector<std::vector<VertexID>> notInBag(t.numNodes);
    std::vector<int> defaultSizes = getMappingSizes(t, pt);
    int depth = 0;
    PrefixNode *tmp = new PrefixNode(99);
    tmp->children.push_back(pt);
    PrefixNode *pn = tmp;
    defaultSizes = getMappingSizes(t, pt);
    lengthToPT[defaultLengths] = pt;
    lengthToTree[defaultLengths] = &t;
    lengthToMappingSizes[defaultLengths] = defaultSizes;
    lengthToFixed[defaultLengths] = defaultFixed;
    int version = 0;
    std::vector<PrefixNode *> dynamicPartition;
    dynamicPartition.push_back(pt);
    trees[0] = &t;
    mappingSizes[0] = defaultSizes;
    prefixLengths[0] = defaultLengths;
    fixedAttrs[0] = defaultFixed;
//    depends[0] = lengthToDepend[defaultLengths];
    std::vector<bool> attrDeadEnd(attributeOrder.size(), false);
    bool haveNewVersion = false;
    while (depth >= 0) {
        bool jumpFlag = false;
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            VertexID u = current -> u;
            int attrID;
            if (depth == 0) attrID = 0;
            else {
                VertexID firstBag = bagsBelow[current][0];
                attrID = attrIDMap[firstBag][depth - 1];
            }
            if (pPoses[attrID] == -1) nextCands[attrID] = true;
            ids[depth] = attrID;
            path[depth] = current;
            bool nextLevel = false, newBranch = false;
            std::vector<bool> currentFixed = fixedAttrs[version];
            HyperTree *currentT = trees[version];
            std::vector<int> currentSize = mappingSizes[version];
//            std::vector<std::pair<int, VertexID>> currentDepend = depends[version];
            if (!skipMatch[attrID]) haveNewVersion = false;
            if (!nextCands[attrID]) goto finishNextCand;
            for (PrefixNode *c: current->children) {
                int id = attrIDMap[bagsBelow[c][0]][depth];
                depend[id].first = attrID;
                depend[id].second = -1;
            }
            if (!current -> pathToGlobal) {
                if (!skipMatch[attrID] && currentFixed[u] &&
                    dynamicPartition.back()->u != current->u) dynamicPartition.push_back(current);
                ++pPoses[attrID];
                if (pPoses[attrID] != 0 && !skipMatch[attrID]) visited[partMatch[u]] = false;
                while (pPoses[attrID] < pCandCount[attrID]) {
                    for (VertexID nID: bagsBelow[current]) skipBuild[nID] = false;
                    VertexID v = pCandidates[attrID][pPoses[attrID]];
                    if (visited[v] && !skipMatch[attrID]) {
                        ++pPoses[attrID];
                        continue;
                    }
                    visited[v] = true;
                    partMatch[u] = v;
                    bool deadend = moveToVertex(current, lastNodes, traversedNodes, empty, v);
                    if (deadend) {
                        ++pPoses[attrID];
                        if (!skipMatch[attrID]) visited[v] = false;
                        for (VertexID nID: current->nIDsToJoin) {
                            traversedNodes[nID].pop_back();
                            lastNodes[nID] = traversedNodes[nID].back();
                        }
                        continue;
                    }
                    for (PrefixNode *c: current->children) {
                        int cID = attrIDMap[bagsBelow[c][0]][depth];
                        attrDeadEnd[cID] = true;
                    }
                    int oldVersion = version;
                    for (VertexID nID: current -> nIDsToCall) {
                        if (version != oldVersion) continue;
                        for (int i = 0; i < nodePoses[nID].size(); ++i)
                            nodePoses[nID][i] = 0;
                        ui oldSize = tuples[nID].size();
                        if (buildFromExist(*currentT, nID, partMatch, tuples[nID], backup[nID], budget, trieNodes, depth)) {
                            attrDeadEnd[attrID] = false;
                            continue;
                        }
                        else if (!adaptiveNodeJoin(*currentT, nID, cs, visited, notInBag[nID], partMatch,
                                                   currentSize[nID],
                                                   nodeCandidates[nID], nodeCandCount[nID], nodePoses[nID], tuples[nID],
                                                   budget)) {
                            haveNewVersion = true;
                            std::vector<ui> length = prefixLengths[version];
                            PrefixNode *old = lengthToPT[length];
                            std::vector<std::vector<VertexID>> newTuples;
                            // find the new prefix and build a smaller trie
                            ui oldLength = currentT->nodes[nID].prefixSize;
                            ui newLength = oldLength + 1;
                            if (tuples[nID].size() <= 1) {
                                newLength = currentT->nodes[nID].numAttributes;
                                std::vector<VertexID> resetTuple;
                                if (tuples[nID].size() == 1) {
                                    resetTuple = tuples[nID][0];
                                    budget += tuples[nID][0].size() * sizeof(VertexID);
                                    tuples[nID].clear();
                                }
                                else {
                                    for (int i = oldLength; i < newLength; ++i) {
                                        resetTuple.push_back(partMatch[currentT->nodes[nID].attributes[i]]);
                                    }
                                }
                                resetPoses(*currentT, t, nID, nodeCandidates[nID],
                                           nodeCandCount[nID],
                                           nodePoses[nID],
                                           pCandidates, pCandCount, pPoses, newLength, depth, resetTuple,
                                           visited, cs, partMatch, attrIDMap);
                            }
                            else {
                                const std::vector<VertexID> &tuple = tuples[nID][0];
                                for (int i = 0; i < tuple.size(); ++i) {
                                    bool flag = false;
                                    std::vector<VertexID> newPrefix(tuple.begin(), tuple.begin() + i + 1);
                                    int pos = findFirstGreater(tuples[nID], newPrefix);
                                    if (pos != tuples[nID].size()) {
                                        flag = true;
                                        ui newSize = newLength - currentT->nodes[nID].prefixSize;
                                        resetPoses(*currentT, t, nID, nodeCandidates[nID],
                                                   nodeCandCount[nID],
                                                   nodePoses[nID],
                                                   pCandidates, pCandCount, pPoses, newLength, depth, tuple,
                                                   visited, cs, partMatch, attrIDMap);
                                        budget += tuples[nID][0].size() * (tuples[nID].size() - pos) * sizeof(VertexID);
                                        budget += pos * (i + 1) * sizeof(VertexID);
                                        tuples[nID].erase(tuples[nID].begin() + pos, tuples[nID].end());
                                        tuples[nID].shrink_to_fit();
                                        for (auto &tp : tuples[nID]) {
                                            tp.assign(tp.begin() + i + 1, tp.end());
                                        }
                                    }
                                    if (flag) break;
                                    else ++newLength;
                                }
                            }
                            for (int i = oldLength; i < newLength; ++i) {
                                VertexID u2 = currentT->nodes[nID].attributes[i];
                                int pos = bagToOrder[nID];
                                for (int j = pos + 1; j < bagOrder.size(); ++j) {
                                    VertexID nID2 = bagOrder[j];
                                    if (std::find(t.v2n[u2].begin(), t.v2n[u2].end(), nID2) == t.v2n[u2].end()) {
                                        notInBag[nID2].push_back(u2);
                                    }
                                }
                            }
                            length[nID] = newLength;
                            PrefixNode *newPrefix;
                            HyperTree *newT;
                            std::vector<int> newMappingSize = currentSize;
                            std::vector<bool> newFixed = currentFixed;
//                            std::vector<std::pair<int, VertexID>> newDepend = currentDepend;
                            if (lengthToPT.find(length) != lengthToPT.end()) {
                                newPrefix = lengthToPT[length];
                                newT = lengthToTree[length];
                                newMappingSize = lengthToMappingSizes[length];
                                newFixed = lengthToFixed[length];
                                depend = lengthToDepend[length];
                            }
                            else {
                                resetMaterialize(query, *currentT, t, nID, depth, length, newPrefix, newT,
                                                 old, cs, newFixed,
                                                 newMappingSize, depend, tuples, bagsBelow, attrIDMap,
                                                 dynamicPartition,
                                                 partMatch);
                                lengthToPT[length] = newPrefix;
                                lengthToTree[length] = newT;
                                lengthToMappingSizes[length] = newMappingSize;
                                lengthToFixed[length] = newFixed;
                                lengthToDepend[length] = depend;
                            }
                            makeBackup(oldLength, newLength, depth, *newT, tuples, backup, path, nID, partMatch,
                                       bagToOrder, budget);
                            PrefixNode *pn2 = newPrefix;
                            for (int i = 0; i <= depth; ++i) {
                                path[i] = pn2;
                                if (i != depth) pn2 = pn2 -> children[childPoses[i + 1]];
                            }
                            current = path[depth];
                            if (depth > 0) pn = path[depth - 1];
                            int extendSize = newLength;
                            if (newLength > depth) {
                                extendSize = pn2->locate(nID).size() + depth;
                                // one more level
                                for (int i = 0; i < pn2->children.size(); ++i) {
                                    PrefixNode *child = pn2->children[i];
                                    if (std::find(bagsBelow[child].begin(), bagsBelow[child].end(), nID) != bagsBelow[child].end()) {
                                        pn2 = child;
                                        childPoses[depth + 1] = i;
                                        break;
                                    }
                                }
                                // remaining levels
                                for (int i = depth + 1; i <= extendSize; ++i) {
                                    path[i] = pn2;
                                    ids[i] = attrIDMap[nID][i - 1];
                                    if (i != extendSize) {
                                        childPoses[i + 1] = 0;
                                        pn2 = pn2 -> children[0];
                                    }
                                }
                            }
                            else {
                                for (int i = 1; i <= depth; ++i) {
                                    if (newFixed[path[i]->u]) extendSize = i;
                                }
                                if (extendSize != depth) newBranch = true;
                                depth = extendSize;
                            }
                            rebuildTrie(*currentT, *newT, lastNodes, path, extendSize,
                                        skipBuild, tuples, budget, tupleSizes, currentFixed);
                            bool setChangeMaterialize = false;
                            bool valid = true;
                            for (int i = 0; i < extendSize; ++i) {
                                if (!(newFixed[path[i + 1]->u] && !currentFixed[path[i + 1]->u])) continue;
                                if (!setChangeMaterialize) {
                                    setChangeMaterialize = true;
                                    changeMaterialize[ids[i + 1]] = true;
                                }
                                // shrink candidates, starting from next cand
                                for (VertexID u2: path[i + 1] -> attributesBefore) {
                                    if (std::find(newT->nodes[nID].attributes, newT->nodes[nID].attributes + newT->nodes[nID].numAttributes,
                                                  u2) == newT->nodes[nID].attributes + newT->nodes[nID].numAttributes) {
                                        std::vector<VertexID> newCandidates(pCandCount[ids[i + 1]], 0);
                                        int remaining = 0;
                                        for (int j = pPoses[ids[i + 1]] + 1; j < pCandCount[ids[i + 1]]; ++j) {
                                            VertexID v1 = pCandidates[ids[i + 1]][j];
                                            bool connected = true;
                                            for (VertexID u3: path[i + 1] -> attributesBefore) {
                                                if (!cs.checkExists(v1, partMatch[u3])) connected = false;
                                            }
                                            if (connected) newCandidates[remaining++] = v1;
                                        }
                                        pCandCount[ids[i + 1]] = pPoses[ids[i + 1]] + remaining + 1;
                                        for (int j = 0; j < remaining; ++j) {
                                            pCandidates[ids[i + 1]][pPoses[ids[i + 1]] + j + 1] = newCandidates[j];
                                        }
                                        break;
                                    }
                                }
                                if (moveToVertex(path[i + 1], lastNodes, traversedNodes, empty, partMatch[path[i + 1]->u]))
                                    valid = false;
                                for (VertexID u2: path[i + 1] -> attributesBefore) {
                                    VertexID v1 = partMatch[path[i + 1] -> u];
                                    if (!cs.checkExists(v1, partMatch[u2])) {
                                        valid = false;
                                        break;
                                    }
                                }
                                for (VertexID nID2 : path[i + 1] -> nIDsToBuild) {
                                    if (skipBuild[nID2]) continue;
                                    budget += tupleSizes[nID2];
                                    if (!tuples[nID2].empty()) {
                                        tupleSizes[nID2] = tuples[nID2].size() * tuples[nID2][0].size() * sizeof(VertexID);
                                    }
                                    else tupleSizes[nID2] = 0;
                                    if (!valid) {
                                        budget += tupleSizes[nID2];
                                        tupleSizes[nID2] = 0;
                                        tuples[nID2].clear();
                                    }
                                    buildTrie(tuples[nID2], trieNodes[nID2], newT->trieOrder[nID2]);
                                    skipBuild[nID2] = true;
                                }
                                nextCands[ids[i + 1]] = false;
                                dynamicPartition.push_back(path[i + 1]);
                                VertexID u2 = currentT->nodes[nID].attributes[i];
                                for (int id: uToAttrID[u2]) {
                                    if (id > ids[i + 1])
                                        skipMatch[id] = true;
                                }
                            }
                            ++version;
                            versionToNID[version] = nID;
                            currentT = trees[version] = newT;
                            prefixLengths[version] = length;
                            currentSize = mappingSizes[version] = newMappingSize;
                            currentFixed = fixedAttrs[version] = newFixed;
//                            currentDepend = depends[version] = newDepend;
                            if (newBranch) break;
                        }
                        else if (tuples[nID].size() == oldSize) {
                            deadend = true;
                            break;
                        }
                        else attrDeadEnd[attrID] = false;
                    }
                    if (currentFixed[u] && !skipMatch[attrID]) nextCands[attrID] = false;
                    if (deadend) {
                        ++pPoses[attrID];
                        for (VertexID nID : current -> nIDsToJoin) {
                            traversedNodes[nID].pop_back();
                            lastNodes[nID] = traversedNodes[nID].back();
                        }
                        if (!skipMatch[attrID]) visited[v] = false;
                        if (!nextCands[attrID]) {
                            for (VertexID nID: current -> nIDsToCall) {
                                if (!tuples[nID].empty()) {
                                    budget += tuples[nID].size() * tuples[nID][0].size() * sizeof(VertexID);
                                    tuples[nID].clear();
                                }
                            }
                        }
                        continue;
                    }
                    if (newBranch) {
                        break;
                    }
                    if (!current -> children.empty()) {
                        nextLevel = true;
                        break;
                    }
                    if (nextCands[attrID]) {
                        if (!skipMatch[attrID]) visited[v] = false;
                        ++pPoses[attrID];
                    }
                    else {
                        for (VertexID nID2 : current -> nIDsToCall) {
                            if (skipBuild[nID2]) continue;
                            budget += tupleSizes[nID2];
                            if (!tuples[nID2].empty())
                                tupleSizes[nID2] = tuples[nID2].size() * tuples[nID2][0].size() * sizeof(VertexID);
                            else tupleSizes[nID2] = 0;
                            buildTrie(tuples[nID2], trieNodes[nID2], currentT->trieOrder[nID2]);
                            skipBuild[nID2] = true;
                        }
                        break;
                    }
                }
            }
            else {
                if (!skipMatch[attrID] && dynamicPartition.back()->u != current->u)
                    dynamicPartition.push_back(current);
                ++pPoses[attrID];
                if (pPoses[attrID] != 0 && !skipMatch[attrID] && depth != 0) visited[partMatch[u]] = false;
                while (pPoses[attrID] < pIterSizes[attrID]) {
                    for (VertexID nID: bagsBelow[current]) skipBuild[nID] = false;
                    VertexID v;
                    if (depth != 0) {
                        ui childPos = pIters[attrID][pPoses[attrID]][0];
                        if (current -> nIDsToJoin.empty()) {
                            v = edgeColumns[depth][0][0][childPos]->value;
                        }
                        else {
                            VertexID u2 = current -> nIDsToJoin[0];
                            v = traversedNodes[u2].back()->nodeChild[childPos]->value;
                        }
                        if (visited[v] && !skipMatch[attrID]) {
                            ++pPoses[attrID];
                            continue;
                        }
                        visited[v] = true;
                        partMatch[u] = v;
                    }
                    // all joined relations extend one level
                    for (int i = 0; i < current -> nIDsToJoin.size(); ++i) {
                        VertexID nID = current -> nIDsToJoin[i];
                        ui childPos = pIters[attrID][pPoses[attrID]][i];
                        traversedNodes[nID].push_back(traversedNodes[nID].back()->nodeChild[childPos]);
                        lastNodes[nID] = traversedNodes[nID].back();
                    }
                    bool deadend = false;
                    int oldVersion = version;
                    for (PrefixNode *c: current->children) {
                        int cID = attrIDMap[bagsBelow[c][0]][depth];
                        attrDeadEnd[cID] = true;
                    }
                    for (VertexID nID: current -> nIDsToCall) {
                        if (nID == currentT->numNodes - 1) continue;
                        if (version != oldVersion) continue;
                        for (int i = 0; i < nodePoses[nID].size(); ++i)
                            nodePoses[nID][i] = 0;
                        ui oldSize = tuples[nID].size();
                        if (buildFromExist(*currentT, nID, partMatch, tuples[nID], backup[nID], budget, trieNodes,
                                           depth)) {
                            attrDeadEnd[attrID] = false;
                            continue;
                        }
                        else if (!adaptiveNodeJoin(*currentT, nID, cs, visited, notInBag[nID], partMatch,
                                                   currentSize[nID],
                                                   nodeCandidates[nID], nodeCandCount[nID], nodePoses[nID],
                                                   tuples[nID], budget)) {
                            haveNewVersion = true;
                            std::vector<ui> length = prefixLengths[version];
                            PrefixNode *old = lengthToPT[length];
                            std::vector<std::vector<VertexID>> newTuples;
                            // find the new prefix and build a smaller trie
                            ui oldLength = currentT->nodes[nID].prefixSize;
                            ui newLength = oldLength+ 1;
                            if (tuples[nID].size() <= 1) {
                                newLength = currentT->nodes[nID].numAttributes;
                                std::vector<VertexID> resetTuple;
                                if (tuples[nID].size() == 1) {
                                    resetTuple = tuples[nID][0];
                                    budget += tuples[nID][0].size() * sizeof(VertexID);
                                    tuples[nID].clear();
                                }
                                else {
                                    for (int i = oldLength; i < newLength; ++i) {
                                        resetTuple.push_back(partMatch[currentT->nodes[nID].attributes[i]]);
                                    }
                                }
                                resetPoses(*currentT, t, nID, nodeCandidates[nID],
                                           nodeCandCount[nID],
                                           nodePoses[nID],
                                           pCandidates, pCandCount, pPoses, newLength, depth, resetTuple,
                                           visited, cs, partMatch, attrIDMap);
                            }
                            else {
                                const std::vector<VertexID> &tuple = tuples[nID][0];
                                for (int i = 0; i < tuple.size(); ++i) {
                                    bool flag = false;
                                    std::vector<VertexID> newPrefix(tuple.begin(), tuple.begin() + i + 1);
                                    int pos = findFirstGreater(tuples[nID], newPrefix);
                                    if (pos != tuples[nID].size()) {
                                        flag = true;
                                        ui newSize = newLength - currentT->nodes[nID].prefixSize;
                                        resetPoses(*currentT, t, nID, nodeCandidates[nID],
                                                   nodeCandCount[nID],
                                                   nodePoses[nID],
                                                   pCandidates, pCandCount, pPoses, newLength, depth, tuple,
                                                   visited, cs, partMatch, attrIDMap);
                                        budget += tuples[nID][0].size() * (tuples[nID].size() - pos) * sizeof(VertexID);
                                        budget += pos * (i + 1) * sizeof(VertexID);
                                        tuples[nID].erase(tuples[nID].begin() + pos, tuples[nID].end());
                                        tuples[nID].shrink_to_fit();
                                        for (auto &tp : tuples[nID]) {
                                            tp.assign(tp.begin() + i + 1, tp.end());
                                        }
                                    }
                                    if (flag) break;
                                    else ++newLength;
                                }
                            }
                            for (int i = oldLength; i < newLength; ++i) {
                                VertexID u2 = currentT->nodes[nID].attributes[i];
                                int pos = bagToOrder[nID];
                                for (int j = pos + 1; j < bagOrder.size(); ++j) {
                                    VertexID nID2 = bagOrder[j];
                                    if (std::find(t.v2n[u2].begin(), t.v2n[u2].end(), nID2) == t.v2n[u2].end()) {
                                        notInBag[nID2].push_back(u2);
                                    }
                                }
                            }
                            length[nID] = newLength;
                            PrefixNode *newPrefix;
                            HyperTree *newT;
                            std::vector<int> newMappingSize = currentSize;
                            std::vector<bool> newFixed = currentFixed;
//                            std::vector<std::pair<int, VertexID>> newDepend = currentDepend;
                            if (lengthToPT.find(length) != lengthToPT.end()) {
                                newPrefix = lengthToPT[length];
                                newT = lengthToTree[length];
                                newMappingSize = lengthToMappingSizes[length];
                                newFixed = lengthToFixed[length];
                                depend = lengthToDepend[length];
                            }
                            else {
                                resetMaterialize(query, *currentT, t, nID, depth, length, newPrefix, newT,
                                                 old, cs, newFixed,
                                                 newMappingSize, depend, tuples, bagsBelow, attrIDMap,
                                                 dynamicPartition,
                                                 partMatch);
                                lengthToPT[length] = newPrefix;
                                lengthToTree[length] = newT;
                                lengthToMappingSizes[length] = newMappingSize;
                                lengthToFixed[length] = newFixed;
                                lengthToDepend[length] = depend;
                            }
                            makeBackup(oldLength, newLength, depth, *newT, tuples, backup, path, nID, partMatch,
                                       bagToOrder, budget);
                            PrefixNode *pn2 = newPrefix;
                            for (int i = 0; i <= depth; ++i) {
                                path[i] = pn2;
                                if (i != depth) pn2 = pn2 -> children[childPoses[i + 1]];
                            }
                            current = path[depth];
                            if (depth > 0) pn = path[depth - 1];
                            int extendSize = newLength;
                            if (newLength > depth) {
                                extendSize = pn2->locate(nID).size() + depth;
                                // one more level
                                for (int i = 0; i < pn2->children.size(); ++i) {
                                    PrefixNode *child = pn2->children[i];
                                    if (std::find(bagsBelow[child].begin(), bagsBelow[child].end(), nID) != bagsBelow[child].end()) {
                                        pn2 = child;
                                        childPoses[depth + 1] = i;
                                        break;
                                    }
                                }
                                // remaining levels
                                for (int i = depth + 1; i <= extendSize; ++i) {
                                    path[i] = pn2;
                                    ids[i] = attrIDMap[nID][i - 1];
                                    if (i != extendSize) {
                                        childPoses[i + 1] = 0;
                                        if (!pn2->children.empty()) pn2 = pn2 -> children[0];
                                    }
                                }
//                                updateSkipMatch(path, depth, newLength, skipMatch, attrIDMap, bagsBelow, dynamicPartition);
                                rebuildTrie(*currentT, *newT, lastNodes, path, extendSize,
                                            skipBuild, tuples, budget, tupleSizes, currentFixed);
                            }
                            bool setChangeMaterialize = false;
                            bool valid = true;
                            for (int i = 0; i < extendSize; ++i) {
                                if (!(newFixed[path[i + 1]->u] && !currentFixed[path[i + 1]->u])) continue;
                                if (!setChangeMaterialize) {
                                    setChangeMaterialize = true;
                                    changeMaterialize[ids[i + 1]] = true;
                                }
                                // shrink candidates, starting from next cand
                                for (VertexID u2: path[i + 1] -> attributesBefore) {
                                    if (std::find(newT->nodes[nID].attributes, newT->nodes[nID].attributes + newT->nodes[nID].numAttributes,
                                                  u2) == newT->nodes[nID].attributes + newT->nodes[nID].numAttributes) {
                                        std::vector<VertexID> newCandidates(pCandCount[ids[i + 1]], 0);
                                        int remaining = 0;
                                        for (int j = pPoses[ids[i + 1]] + 1; j < pCandCount[ids[i + 1]]; ++j) {
                                            VertexID v1 = pCandidates[ids[i + 1]][j];
                                            bool connected = true;
                                            for (VertexID u3: path[i + 1] -> attributesBefore) {
                                                if (!cs.checkExists(v1, partMatch[u3])) connected = false;
                                            }
                                            if (connected) newCandidates[remaining++] = v1;
                                        }
                                        pCandCount[ids[i + 1]] = pPoses[ids[i + 1]] + remaining + 1;
                                        for (int j = 0; j < remaining; ++j) {
                                            pCandidates[ids[i + 1]][pPoses[ids[i + 1]] + j + 1] = newCandidates[j];
                                        }
                                        break;
                                    }
                                }
                                if (moveToVertex(path[i + 1], lastNodes, traversedNodes, empty, partMatch[path[i + 1]->u]))
                                    valid = false;
                                for (VertexID u2: path[i + 1] -> attributesBefore) {
                                    VertexID v1 = partMatch[path[i + 1] -> u];
                                    if (!cs.checkExists(v1, partMatch[u2])) {
                                        valid = false;
                                        break;
                                    }
                                }
                                for (VertexID nID2 : path[i + 1] -> nIDsToBuild) {
                                    if (skipBuild[nID2]) continue;
                                    budget += tupleSizes[nID2];
                                    if (!tuples[nID2].empty()) {
                                        tupleSizes[nID2] = tuples[nID2].size() * tuples[nID2][0].size() * sizeof(VertexID);
                                    }
                                    else tupleSizes[nID2] = 0;
                                    if (!valid) {
                                        budget += tupleSizes[nID2];
                                        tupleSizes[nID2] = 0;
                                        tuples[nID2].clear();
                                    }
                                    buildTrie(tuples[nID2], trieNodes[nID2], newT->trieOrder[nID2]);
                                    skipBuild[nID2] = true;
                                }
                                nextCands[ids[i + 1]] = false;
                                dynamicPartition.push_back(path[i + 1]);
                                VertexID u2 = currentT->nodes[nID].attributes[i];
                                for (int id: uToAttrID[u2]) {
                                    if (id > ids[i + 1])
                                        skipMatch[id] = true;
                                }
                            }
                            ++version;
                            versionToNID[version] = nID;
                            currentT = trees[version] = newT;
                            prefixLengths[version] = length;
                            currentSize = mappingSizes[version] = newMappingSize;
                            currentFixed = fixedAttrs[version] = newFixed;
//                            currentDepend = depends[version] = newDepend;
                        }
                        else if (tuples[nID].size() == oldSize) {
                            for (VertexID nID2 : current -> nIDsToCall) {
                                if (!tuples[nID2].empty()) budget += tuples[nID2][0].size() * tuples[nID2].size() * sizeof(VertexID);
                                tuples[nID2].clear();
                            }
                            deadend = true;
                            break;
                        }
                        else attrDeadEnd[attrID] = false;
                    }
                    if (depth != 0 && deadend) {
                        if (!skipMatch[attrID]) visited[v] = false;
                        ++pPoses[attrID];
                        for (VertexID nID : current -> nIDsToJoin) {
                            traversedNodes[nID].pop_back();
                            lastNodes[nID] = traversedNodes[nID].back();
                        }
                        for (VertexID nID: current->nIDsToCall) {
                            if (!tuples[nID].empty() && currentT->nodes[nID].prefixSize == depth) {
                                budget += tuples[nID].size() * tuples[nID][0].size() * sizeof(VertexID);
                                tuples[nID].clear();
                            }
                        }
                        continue;
                    }
                    if (!current -> children.empty()) {
                        nextLevel = true;
                        break;
                    }
                    else {
                        VertexID nID = currentT->numNodes - 1;
                        bool flag = true;
                        for (VertexID nID2 : current -> nIDsToBuild) {
                            if (skipBuild[nID2]) continue;
                            if (tuples[nID2].empty()) {
                                flag = false;
                                break;
                            }
                        }
                        if (flag) {
                            for (VertexID nID2 : current -> nIDsToBuild) {
                                if (skipBuild[nID2]) continue;
                                budget += tupleSizes[nID2];
                                tupleSizes[nID2] = tuples[nID2].size() * tuples[nID2][0].size() * sizeof(VertexID);
                                buildTrie(tuples[nID2], trieNodes[nID2], currentT->trieOrder[nID2]);
                                skipBuild[nID2] = true;
                            }
                            adaptiveGlobalJoin(result, budget, count, query, *currentT, cs, partMatch,
                                               currentSize[nID], depth,
                                               lastNodes, visited, currentT->nodes[lastID].nIDs,
                                               currentT->nodes[lastID].attributesBefore,
                                               currentT->nodes[lastID].cartesianParent, nIters, nIterSizes, traverse,
                                               edgeColumns, empty, skip);
                        }
                        else {
                            for (VertexID nID2 : current -> nIDsToBuild) {
                                if (!tuples[nID2].empty()) {
                                    budget += tuples[nID2].size() * tuples[nID2][0].size() * sizeof(VertexID);
                                    tuples[nID2].clear();
                                }
                            }
                        }
                        ++pPoses[attrID];
                        if (current->u != 99 && !skipMatch[attrID]) {
                            visited[v] = false;
                            for (VertexID nID : current -> nIDsToJoin) {
                                traversedNodes[nID].pop_back();
                                lastNodes[nID] = traversedNodes[nID].back();
                            }
                        }
                    }
                }
            }
            finishNextCand:
            if (nextLevel) {
                ++depth;
                childPoses[depth] = 0;
                pn = current;
                current = pn -> children[0];
            }
            else if (newBranch) {
                pn = path[depth - 1];
                continue;
            }
            else if (depth != 0 && currentFixed[u] && !skipMatch[attrID] && pPoses[attrID] == pCandCount[attrID])
                break;
            else {
                if (attrDeadEnd[attrID]) {
                    if (pn == dynamicPartition.back() && pn->pathToGlobal || currentFixed[pn->u]) {
                        for (VertexID nID2: bagsBelow[pn]) {
                            if (!tuples[nID2].empty()) {
                                budget += tuples[nID2].size() * tuples[nID2][0].size() * sizeof(VertexID);
                                tuples[nID2].clear();
                            }
                        }
                        break;
                    }
                    else if (!pn->pathToGlobal && !currentFixed[pn->u])
                        break;
                }
                else if (depth > 0) {
                    attrDeadEnd[ids[depth - 1]] = false;
                    if (depend[ids[depth]].first != 0)
                        depend[ids[depth]].second = partMatch[attributeOrder[depend[ids[depth]].first]->u];
                }
                else depend[ids[depth]].second = 0;
//                ++childPoses[depth];
                while (true) {
                    ++childPoses[depth];
                    if (childPoses[depth] == pn->children.size()) break;
                    if (pn->children[childPoses[depth]]->pathToGlobal) break;
                    int id = attrIDMap[bagsBelow[pn->children[childPoses[depth]]][0]][depth - 1];
                    if (pn->u != attributeOrder[depend[id].first]->u ||
                        (pn->u != 99 && partMatch[pn->u] != depend[id].second) || (pn->u == 99 && depend[id].second != 0)) break;
                }
                if (childPoses[depth] == pn->children.size()) break;
                current = pn->children[childPoses[depth]];
            }
            u = current->u;
            path[depth] = current;
            if (depth != 0) {
                VertexID firstBag = bagsBelow[current][0];
                attrID = attrIDMap[firstBag][depth - 1];
                ids[depth] = attrID;
                if (!currentFixed[current->u]) pPoses[attrID] = -1;
            }
            if (u == 99) {
                attrID = 0;
                pPoses[0] = 0;
                pIterSizes[0] = 1;
            }
            else if (current -> pathToGlobal) {
                bool flag = true;
                for (VertexID nID : pn -> nIDsToBuild) {
                    if (skipBuild[nID]) continue;
                    if (tuples[nID].empty()) {
                        flag = false;
                        break;
                    }
                }
                if (flag) {
                    for (VertexID nID : pn -> nIDsToBuild) {
                        if (skipBuild[nID]) continue;
                        budget += tupleSizes[nID];
                        tupleSizes[nID] = tuples[nID].size() * tuples[nID][0].size() * sizeof(VertexID);
                        buildTrie(tuples[nID], trieNodes[nID], currentT->trieOrder[nID]);
                        skipBuild[nID] = true;
                    }
                }
                else {
                    for (VertexID nID : pn -> nIDsToBuild) {
                        if (!tuples[nID].empty())
                            budget += tuples[nID].size() * tuples[nID][0].size() * sizeof(VertexID);
                        tuples[nID].clear();
                    }
                    // go back to the last partition
                    pPoses[attrID] = pIterSizes[attrID] = 0;
                    if (!skipMatch[attrID]) dynamicPartition.push_back(current);
                    break;
                }
                attributeOpen(skipMatch, cs, attrID, depth, current, lastNodes, edgeColumns, pIters, pIterSizes,
                              pCandidates, pCandCount, pPoses, children, partMatch, currentFixed, neighbors,
                              neighborCount, haveNewVersion);
                if (pIterSizes[attrID] == 0 && !skipMatch[attrID]) {
                    dynamicPartition.push_back(current);
                    pPoses[attrID] = 0;
                    break;
                }
            }
            else {
                attributeOpen(skipMatch, cs, attrID, depth, current, lastNodes, edgeColumns, pIters, pIterSizes,
                              pCandidates, pCandCount, pPoses, children, partMatch, currentFixed, neighbors,
                              neighborCount, haveNewVersion);
                if (pCandCount[attrID] == 0) {
                    for (VertexID nID: current->nIDsToCall) {
                        for (int i = 0; i < trieNodes[nID]->nodeChild.size(); ++i) {
                            delete trieNodes[nID]->nodeChild[i];
                        }
                        trieNodes[nID] -> nodeChild.clear();
                    }
                    pPoses[attrID] = 0;
                    if (currentFixed[u])
                        dynamicPartition.push_back(current);
                    else {
                        if (childPoses[depth] > 0 && currentFixed[pn->children[childPoses[depth] - 1]->u]) {
                            jumpFlag = true;
                        }
                    }
                    break;
                }
            }
        }
        // for backtracking, may not go back to parent.
        PrefixNode *current = path[depth];
        if (depth > 0) pn = path[depth - 1];
        else pn = tmp;
        VertexID u = path[depth]->u;
        HyperTree *currentT = trees[version];
        std::vector<int> currentSize = mappingSizes[version];
        std::vector<bool> currentFixed = fixedAttrs[version];
        int attrID = ids[depth];
        if (current->pathToGlobal || (currentFixed[current->u] && !skipMatch[attrID])) {
            if (!current->pathToGlobal && pPoses[attrID] == pCandCount[attrID]) {
                dynamicPartition.pop_back();
                depend[attrID].second = -1;
                jumpFlag = true;
                for (VertexID nID = 0; nID < t.numNodes; ++nID) {
                    if (std::find(notInBag[nID].begin(), notInBag[nID].end(), u) != notInBag[nID].end()) {
                        std::vector<VertexID> newNotIn;
                        for (VertexID u2: notInBag[nID]) {
                            if (u2 != u)
                                newNotIn.push_back(u2);
                        }
                        notInBag[nID] = newNotIn;
                    }
                }
                if (changeMaterialize[attrID]) {
                    changeMaterialize[attrID] = false;
                    VertexID nID = versionToNID[version];
                    std::vector<PrefixNode *> pathToNID = current->locate(nID);
                    for (int i = 0; i < pathToNID.size(); ++i) {
                        path[depth + i + 1] = pathToNID[i];
                    }
                    int extendSize = depth + pathToNID.size();
                    const HyperTree &oldT = *currentT;
                    --version;
                    currentT = trees[version];
                    const HyperTree &newT = *currentT;
                    currentSize = mappingSizes[version];
                    std::vector<bool> oldFixed = currentFixed;
                    currentFixed = fixedAttrs[version];
                    rebuildTrie(oldT, newT, lastNodes, path, extendSize, skipBuild, tuples, budget, tupleSizes, currentFixed);
                    PrefixNode *pn2 = lengthToPT[prefixLengths[version]];
                    for (int i = 0; i <= depth - 1; ++i) {
                        path[i] = pn2;
                        if (pn2->children.empty()) {
                            depth = i + 1;
                            break;
                        }
                        if (i != depth - 1) pn2 = pn2 -> children[childPoses[i + 1]];
                    }
                    pn = path[depth - 1];
                    std::vector<VertexID> newPrefix;
                    for (VertexID u2 = 0; u2 < query.getNumVertices(); ++u2) {
                        if (oldFixed[u2] && !currentFixed[u2])
                            newPrefix.push_back(u2);
                    }
                    for (int i = 0; i < attributeOrder.size(); ++i) {
                        if (skipMatch[i] && std::find(newPrefix.begin(), newPrefix.end(), attributeOrder[i]->u) != newPrefix.end()) {
                            skipMatch[i] = false;
                        }
                    }
                }
            }
            else if (current->pathToGlobal && pPoses[attrID] == pIterSizes[attrID]) {
                jumpFlag = true;
                if (!skipMatch[attrID]) dynamicPartition.pop_back();
                else {
                    pPoses[attrID] = -1;
                    nextCands[attrID] = true;
                }
            }
        }
        --depth;
        bool deadEnd = attrDeadEnd[attrID];
        if (depth >= 0) {
            current = pn;
            if (depth == 0) pn = tmp;
            else pn = path[depth - 1];
            attrID = ids[depth];
            if (jumpFlag) {
                jump(dynamicPartition, current, pn, childPoses, ids, path, nextCands, depth, attrIDMap, tmp);
            }
            else {
                for (VertexID nID: current -> nIDsToCall) {
                    if (nID == lastID) {
                        bool flag = true;
                        for (VertexID nID2 : current -> nIDsToBuild) {
                            if (skipBuild[nID2]) continue;
                            if (tuples[nID2].empty()) {
                                flag = false;
                            }
                            budget += tupleSizes[nID2];
                            if (!tuples[nID2].empty())
                                tupleSizes[nID2] = tuples[nID2].size() * tuples[nID2][0].size() * sizeof(VertexID);
                            else tupleSizes[nID2] = 0;
                            buildTrie(tuples[nID2], trieNodes[nID2], currentT->trieOrder[nID2]);
                            skipBuild[nID2] = true;
                        }
                        if (flag) {
                            adaptiveGlobalJoin(result, budget, count, query, *currentT, cs, partMatch,
                                               currentSize[nID],
                                               depth, lastNodes, visited, currentT->nodes[lastID].nIDs,
                                               currentT->nodes[lastID].attributesBefore,
                                               currentT->nodes[lastID].cartesianParent, nIters, nIterSizes,
                                               traverse, edgeColumns, empty, skip);
                        }
                        // go back to the last partition
                        jump(dynamicPartition, current, pn, childPoses, ids, path, nextCands, depth, attrIDMap, tmp);
                    }
                }
            }
            if (depth != 0) {
                u = current -> u;
                attrID = ids[depth];
                if (!deadEnd) attrDeadEnd[attrID] = false;
                if (nextCands[attrID] && !skipMatch[attrID]) {
                    VertexID v = partMatch[u];
                    visited[v] = false;
                    for (VertexID nID : current -> nIDsToJoin) {
                        traversedNodes[nID].pop_back();
                        lastNodes[nID] = traversedNodes[nID].back();
                    }
                    if (currentFixed[u] || current->pathToGlobal) {
                        for (VertexID nID: current->nIDsToBuild) {
                            if (!tuples[nID].empty()) budget += tuples[nID].size() * tuples[nID][0].size() * sizeof(VertexID);
                            tuples[nID].clear();
                        }
                    }
                }
                else if (!nextCands[attrID] && !current->pathToGlobal) {
                    for (VertexID nID2 : current -> nIDsToBuild) {
                        if (skipBuild[nID2]) continue;
                        budget += tupleSizes[nID2];
                        if (!tuples[nID2].empty()) {
                            tupleSizes[nID2] = tuples[nID2].size() * tuples[nID2][0].size() * sizeof(VertexID);
                        }
                        else tupleSizes[nID2] = 0;
                        buildTrie(tuples[nID2], trieNodes[nID2], currentT->trieOrder[nID2]);
                        skipBuild[nID2] = true;
                    }
                }
            }
        }
    }
    tmp->children.clear();
    delete tmp;
    delete[] partMatch;

    for (int i = 0; i < trieNodes.size(); ++i) {
        delete[] nodeCandCount[i];
        delete[] nodeCandidates[i];
    }
    delete[] nodeCandidates;
    delete[] nodeCandCount;
    delete[] pCandidates;
    delete[] pCandCount;
    delete[] neighbors;
    delete[] neighborCount;
    for (int i = 0; i < attributeOrder.size(); ++i) {
        PrefixNode *attribute = attributeOrder[i];
        if (attribute -> pathToGlobal && attribute->u != 99) {
            for (int j = 0; j < cs.candidateSet[attribute->u].size(); ++j) delete[] pIters[i][j];
            delete[] pIters[i];
        }
    }
    delete[] pIters;
    delete[] pIterSizes;
    for (int i = 0; i < t.nodes[lastID].numAttributes + 1; ++i) {
        for (int j = 0; j < maxSize; ++j) {
            delete[] nIters[i][j];
        }
        delete[] nIters[i];
    }
    delete[] nIters;
    delete[] nIterSizes;
    delete[] children;
    for (int i = 0; i < t.nodes[lastID].numAttributes + 1; ++i) {
        for (int j = 0; j < maxEdgeSize; ++j) {
            DynamicArray<TrieNode *> *dynamicArrayPtr = edgeColumns[i][j];
            for (ui k = 0; k < dynamicArrayPtr->size(); ++k) {
                delete (*dynamicArrayPtr)[k];
            }
            dynamicArrayPtr->clear();
            delete dynamicArrayPtr;
        }
    }

    for (auto & it : lengthToPT) {
        if (it.first != defaultLengths)
            delete it.second;
    }

    for (auto & it : lengthToTree) {
        if (it.first != defaultLengths)
            delete it.second;
    }
    for (size_t b: tupleSizes) budget += b;
    for (VertexID nID = 0; nID < tuples.size(); ++nID) {
        if (!tuples[nID].empty()) budget += tuples[nID].size() * tuples[nID][0].size() * sizeof(VertexID);
    }
    if (budget != oldBudget)
        exit(1);
}
