//
// Created by anonymous authors on 2024/2/27.
//

#include "decomposition.h"

void HyperNode::initPoses(const Graph &query, const std::vector<std::vector<size_t>> &dist) {
    attributesBefore = std::vector<std::vector<VertexID>>(numAttributes);
    cartesianParent = std::vector<VertexID>(numAttributes);
    bool cartesian = false;
    for (int i = 0; i < numAttributes; ++i) {
        VertexID u = attributes[i];
        for (int j = 0; j < i; ++j) {
            VertexID u2 = attributes[j];
            if (query.getEdgeID(u, u2) != -1)
                attributesBefore[i].push_back(u2);
        }
        if (i > 0 && attributesBefore[i].empty()) cartesian = true;
    }
    if (cartesian) {
        for (ui i = 1; i < numAttributes; ++i) {
            VertexID u = attributes[i];
            if (attributesBefore[i].empty()) {
                size_t minDist = std::numeric_limits<size_t>::max();
                for (VertexID j = 0; j < i; ++j) {
                    VertexID u2 = attributes[j];
                    if (dist[u][u2] < minDist) {
                        minDist = dist[u][u2];
                        cartesianParent[i] = u2;
                    }
                }
            }
        }
    }
}

void HyperNode::initPoses(const vector<std::vector<VertexID>> &sharedAttrs, const Graph &query,
                          const vector<std::vector<size_t>> &dist, bool global) {
    attributesBefore = std::vector<std::vector<VertexID>>(numAttributes);
    cartesianParent = std::vector<VertexID>(numAttributes, 99);
    bool cartesian = false;
    for (int i = 0; i < numAttributes; ++i) {
        VertexID u = attributes[i];
        for (int j = 0; j < i; ++j) {
            VertexID u2 = attributes[j];
            if (query.getEdgeID(u, u2) != -1)
                attributesBefore[i].push_back(u2);
        }
        if (i > 0 && attributesBefore[i].empty()) cartesian = true;
    }
    if (global) {
        nIDs = std::vector<std::vector<VertexID>>(numAttributes);
        // remove attributes that are covered by materialized bags
        for (int i = 0; i < numAttributes; ++i) {
            std::set<VertexID> coveredAttrs;
            VertexID u = attributes[i];
            for (VertexID nID = 0; nID < sharedAttrs.size() - 1; ++nID) {
                if (std::find(sharedAttrs[nID].begin(), sharedAttrs[nID].end(), u) != sharedAttrs[nID].end()) {
                    nIDs[i].push_back(nID);
                    for (VertexID u2: sharedAttrs[nID])
                        coveredAttrs.insert(u2);
                }
            }
            std::vector<ui> attrsToJoin;
            for (VertexID u2: attributesBefore[i]) {
                if (coveredAttrs.find(u2) == coveredAttrs.end())
                    attrsToJoin.push_back(u2);
            }
            attributesBefore[i] = attrsToJoin;
        }
    }
    if (cartesian) {
        for (ui i = 1; i < numAttributes; ++i) {
            VertexID u = attributes[i];
            if (attributesBefore[i].empty()) {
                size_t minDist = std::numeric_limits<size_t>::max();
                for (VertexID j = 0; j < i; ++j) {
                    VertexID u2 = attributes[j];
                    if (dist[u][u2] < minDist) {
                        minDist = dist[u][u2];
                        cartesianParent[i] = u2;
                    }
                }
            }
        }
    }
}

void HyperNode::copyTo(HyperNode &other) const {
    delete[] other.attributes;
    delete[] other.prefix;
    other.numAttributes = numAttributes;
    other.prefixSize = prefixSize;
    other.width = width;
    if (numAttributes != 0) {
        other.attributes = new VertexID[numAttributes];
        std::copy(attributes, attributes + numAttributes, other.attributes);
    } else {
        other.attributes = nullptr;
    }
    if (prefixSize != 0) {
        other.prefix = new VertexID[numAttributes];
        std::copy(prefix, prefix + prefixSize, other.prefix);
    } else {
        other.prefix = nullptr;
    }
    other.attributesBefore = attributesBefore;
    other.cartesianParent = cartesianParent;
    other.nIDs = nIDs;
}

void HyperTree::initPoses(const Graph &query, const CandidateSpace &cs, bool handleNode) {
    numAttributes = query.getNumVertices();
    delete[] v2n;
    v2n = new std::vector<VertexID> [numAttributes];
    for (ui i = 0; i < numNodes; ++i) {
        for (ui j = 0; j < nodes[i].numAttributes; ++j) {
            v2n[nodes[i].attributes[j]].push_back(i);
        }
    }
    if (handleNode) {
        for (int i = 0; i < numNodes; ++i) {
            nodes[i].initPoses(query, cs.dist);
        }
    }

    // set the compression sizes
    compressionSizes = std::vector<ui>(numNodes, 0);
    extendLevel = 0;
    for (int i = 0; i < globalOrder.size(); ++i) {
        if (v2n[globalOrder[i]].size() > 1)
            extendLevel = i + 1;
    }
    buildTrieOrder();
    if (!newGlobalNode && extendLevel < nodes[numNodes - 1].prefixSize) extendLevel = nodes[numNodes - 1].prefixSize;
    nodesAtStep = std::vector<std::vector<VertexID>>(defaultPartition.size());
    for (VertexID nID = 0; nID < numNodes; ++nID) {
        for (int i = 0; i < defaultPartition.size(); ++i) {
            if (std::includes(defaultPartition.begin(), defaultPartition.begin() + i + 1,
                              nodes[nID].prefix, nodes[nID].prefix + nodes[nID].prefixSize)) {
                nodesAtStep[i].push_back(nID);
                break;
            }
        }
    }
    nIDs = std::vector<std::vector<VertexID>>(globalOrder.size());
    for (int i = 0; i < globalOrder.size(); ++i) {
        VertexID u = globalOrder[i];
        nIDs[i] = v2n[u];
    }
    attributesBefore = std::vector<std::vector<VertexID>>(globalOrder.size());
    cartesianParent = std::vector<VertexID>(globalOrder.size());
    for (int i = 0; i < globalOrder.size(); ++i) {
        VertexID u = globalOrder[i];
        for (int j = 0; j < i; ++j) {
            VertexID u2 = globalOrder[j];
            if (query.getEdgeID(u, u2) != -1)
                attributesBefore[i].push_back(u2);
        }
    }
    for (ui i = 1; i < globalOrder.size(); ++i) {
        VertexID u = globalOrder[i];
        if (attributesBefore[i].empty()) {
            size_t minDist = std::numeric_limits<size_t>::max();
            for (VertexID j = 0; j < i; ++j) {
                VertexID u2 = globalOrder[j];
                if (cs.dist[u][u2] < minDist) {
                    minDist = cs.dist[u][u2];
                    cartesianParent[i] = u2;
                }
            }
        }
    }
    adaptiveDepthToNID.clear();
    adaptiveGroups.clear();
    adaptiveLevels.clear();
//    buildTraverseStruct(query);
}

void HyperTree::selectRoot(VertexID &root, std::vector<std::vector<VertexID>> &cohesion,
                           std::vector<std::vector<VertexID>> &children, std::vector<VertexID> &parents) {
    std::queue<std::pair<VertexID, VertexID>> q;
    root = 0;
    double minWidth = std::numeric_limits<double>::max();
    for (VertexID nID = 0; nID < numNodes; ++nID) {
        if (nodes[nID].width < minWidth) {
            root = nID;
            minWidth = nodes[nID].width;
        }
    }
    q.emplace(MAX_PATTERN_SIZE, root);
    children.reserve(numNodes);
    parents.reserve(numNodes);
    while (!q.empty()) {
        VertexID parent = q.front().first, current = q.front().second;
        q.pop();
        if (parent != MAX_PATTERN_SIZE) children[parent].push_back(current);
        parents[current] = parent;
        for (int i = 0; i < edges[current].size(); ++i) {
            VertexID next = edges[current][i];
            if (parent != next) {
                q.emplace(current, next);
            }
        }
    }
    for (VertexID nID = 0; nID < numNodes; ++nID) {
        VertexID parent = parents[nID];
        if (parent == MAX_PATTERN_SIZE) continue;
        for (int i = 0; i < nodes[nID].numAttributes; ++i) {
            for (int j = 0; j < nodes[parent].numAttributes; ++j) {
                VertexID u1 = nodes[nID].attributes[i], u2 = nodes[nID].attributes[j];
                if (u1 == u2)
                    cohesion[nID].push_back(u1);
            }
        }
    }
}

void HyperTree::buildFromTD(FHD &fhd) {
    numNodes = fhd.X.size();
    nodes = new HyperNode[numNodes + 1];
    for (VertexID nID = 0; nID < fhd.X.size(); ++nID) {
        std::vector<size_t> tmpV;
        fhd.X[nID].getelement(tmpV);
        nodes[nID].numAttributes = tmpV.size();
        nodes[nID].attributes = new VertexID[nodes[nID].numAttributes];
        for (int i = 0; i < tmpV.size(); ++i) {
            nodes[nID].attributes[i] = tmpV[i];
        }
    }
    edges = std::vector<std::vector<VertexID>>(numNodes);
    for (auto ei: fhd.eg)
        edges[ei.first].push_back(ei.second);
}

void HyperTree::print(const Graph &query) const {
    for (int i = 0; i < numNodes; ++i) {
        std::cout << "node " << i << " : ";
        for (int j = 0; j < nodes[i].numAttributes; ++j) {
            std::cout << nodes[i].attributes[j] << " ";
        }
        std::cout << std::endl;
        for (int j = 0; j < nodes[i].numAttributes; ++j) {
            for (int k = j + 1; k < nodes[i].numAttributes; ++k) {
                if (query.getEdgeID(nodes[i].attributes[j], nodes[i].attributes[k]) != -1) {
                    std::cout << nodes[i].attributes[j] << " " << nodes[i].attributes[k] << std::endl;
                }
            }
        }
        std::cout << std::endl;
    }
}

void HyperTree::addGlobalNode(const vector<VertexID> &globalAttrs) {
    nodes[numNodes] = HyperNode();
    nodes[numNodes].attributes = new VertexID[globalAttrs.size()];
    nodes[numNodes].numAttributes = globalAttrs.size();
    for (int i = 0; i < globalAttrs.size(); ++i) {
        VertexID u = globalAttrs[i];
        nodes[numNodes].attributes[i] = u;
        if (v2n[u].back() != numNodes) v2n[u].push_back(numNodes);
    }

    ++numNodes;
    newGlobalNode = true;
}

void HyperTree::buildTrieOrder() {
    trieOrder = std::vector<std::vector<VertexID>>(numNodes);
    for (ui i = 0; i < numNodes - 1; ++i) {
        for (ui j = 0; j < globalOrder.size(); ++j) {
            VertexID u = globalOrder[j];
            ui k = nodes[i].prefixSize;
            for (; k < nodes[i].numAttributes; ++k) {
                if (nodes[i].attributes[k] == u) {
                    break;
                }
            }
            if (k != nodes[i].numAttributes) trieOrder[i].push_back(k - nodes[i].prefixSize);
        }
    }
    for (ui j = 0; j < globalOrder.size(); ++j) {
        VertexID u = globalOrder[j];
        ui k = extendLevel;
        for (; k < nodes[numNodes - 1].numAttributes; ++k) {
            if (nodes[numNodes - 1].attributes[k] == u) {
                break;
            }
        }
        if (k != nodes[numNodes - 1].numAttributes) trieOrder[numNodes - 1].push_back(k - extendLevel);
    }
}

void HyperTree::buildTraverseStruct(const Graph &query) {
    buildTraverseAdaptive(query, extendLevel);
    depthToNID = adaptiveDepthToNID[extendLevel];
    levels = adaptiveLevels[extendLevel];
    groups = adaptiveGroups[extendLevel];
    if (!newGlobalNode && extendLevel > nodes[numNodes - 1].prefixSize) {
        trieOrder[numNodes - 1].clear();
        for (int i = 0; i < globalOrder.size(); ++i) {
            VertexID u = globalOrder[i];
            ui k = extendLevel;
            for (; k < nodes[numNodes - 1].numAttributes; ++k) {
                if (nodes[numNodes - 1].attributes[k] == u)
                    break;
            }
            if (k != nodes[numNodes - 1].numAttributes) trieOrder[numNodes - 1].push_back(k - extendLevel);
        }
    }
}

void HyperTree::buildTraverseAdaptive(const Graph &query, int currentExtendLevel) {
    if (adaptiveDepthToNID.find(currentExtendLevel) != adaptiveDepthToNID.end()) return;
    std::vector<VertexID> toExtend;
    int pos = currentExtendLevel;
    if (defaultPartition.size() > nodes[numNodes - 1].numAttributes) pos = defaultPartition.size();
    std::vector<ui> currentLevels = std::vector<ui>(numNodes, 0);
    std::vector<std::vector<VertexID>> localAttr(numNodes);
    for (int i = pos; i < numAttributes; ++i) {
        VertexID u = globalOrder[i];
        VertexID nID = v2n[u][0];
        ++currentLevels[nID];
        localAttr[nID].push_back(u);
    }
    for (VertexID nID = 0; nID < numNodes; ++nID) {
        if (currentLevels[nID] != 0) toExtend.push_back(nID);
    }
    std::sort(toExtend.begin(), toExtend.end(), [currentLevels](VertexID a, VertexID b) {
        return currentLevels[a] > currentLevels[b];
    });
    // divide bags into groups based on labels
    // for each bag, if the last levels have distinct labels, directly multiply the num
    // else, extend to the last level and do intersection for each label
    std::map<LabelID, std::vector<VertexID>> labelToNodes;
    std::set<VertexID> independents;
    // all bags start from one level
    for (VertexID nID : toExtend) {
        VertexID u = localAttr[nID].back();
        LabelID l = query.getVertexLabel(u);
        labelToNodes[l].push_back(nID);
    }
    std::vector<std::vector<VertexID>> currentGroups;
    for (auto &item : labelToNodes) {
        if (item.second.size() > 1) currentGroups.push_back(item.second);
        else independents.insert(item.second[0]);
    }
    for (VertexID nID : toExtend) {
        --currentLevels[nID];
        if (independents.find(nID) == independents.end()) continue;
        ui numAttrs = nodes[nID].numAttributes - nodes[nID].prefixSize;
        VertexID u;
        LabelID l;
        if (currentLevels[nID] > 0) {
            u = localAttr[nID][currentLevels[nID] - 1];
            l = query.getVertexLabel(u);
        }
        std::vector<VertexID> group = {nID};
        currentGroups.push_back(group);
        // expend the level of nID
        while (currentLevels[nID] > 0 && (labelToNodes.find(l) == labelToNodes.end() || labelToNodes[l][0] == nID)) {
            if (labelToNodes.find(l) == labelToNodes.end()) labelToNodes[l].push_back(nID);
            --currentLevels[nID];
            if (currentLevels[nID] > 0) {
                u = localAttr[nID][currentLevels[nID] - 1];
                l = query.getVertexLabel(u);
            }
        }
    }
    std::vector<VertexID> currentdepthToNID;
    for (VertexID nID = 0; nID < numNodes; ++nID) {
        for (int i = 0; i < currentLevels[nID]; ++i) {
            currentdepthToNID.push_back(nID);
        }
    }
    adaptiveDepthToNID[currentExtendLevel] = currentdepthToNID;
    adaptiveLevels[currentExtendLevel] = currentLevels;
    adaptiveGroups[currentExtendLevel] = currentGroups;
}

void HyperTree::writeToStream(ostream &outStream) {
    outStream << "num bags: " << numNodes << std::endl;
    for (VertexID nID = 0; nID < numNodes; ++nID) {
        for (int i = 0; i < nodes[nID].numAttributes; ++i) {
            outStream << nodes[nID].attributes[i];
            if (i != nodes[nID].numAttributes - 1) outStream << ",";
        }
        outStream << std::endl;
    }
}

PrefixNode *buildPrefixTree(std::vector<std::vector<VertexID>> &orders, const Graph &query,
                            const std::vector<std::vector<size_t>> &dist,
                            std::unordered_set<PrefixNode *, PrefixNodePtrHash, PrefixNodePtrEqual> &visitedPT,
                            bool &exist) {
    PrefixNode *root = new PrefixNode(99);
    for (VertexID nID = 0; nID < orders.size(); ++nID) {
        if (orders[nID].empty()) continue;
        const std::vector<VertexID> &order = orders[nID];
        PrefixNode *currentNode = root;
        if (order[0] == 99) {
            root -> nIDsToCall.push_back(nID);
            if (nID == orders.size() - 1) currentNode -> pathToGlobal = true;
            continue;
        }
        for (int i = 0; i < order.size(); ++i) {
            VertexID u = order[i];
            int childIndex = currentNode->findChildIndex(u);
            currentNode->nIDsToCall.push_back(nID);
            if (nID != orders.size() - 1) currentNode -> pathToGlobal = false;
            else currentNode -> pathToGlobal = true;
            if (childIndex == -1) {
                // If the child does not exist, create a new child node
                PrefixNode *newNode = new PrefixNode(u);
                currentNode->children.push_back(newNode);
                currentNode = currentNode->children.back();
            } else {
                // If the child exists, move to the child node
                currentNode = currentNode->children[childIndex];
            }
            if (i == order.size() - 1) {
                currentNode->nIDsToCall.push_back(nID);
                if (nID != orders.size() - 1) currentNode -> pathToGlobal = false;
                else currentNode -> pathToGlobal = true;
            }
        }
    }
    std::queue<PrefixNode *> Q;
    Q.push(root);
    while (!Q.empty()) {
        PrefixNode *pn = Q.front();
        Q.pop();
        for (int i = 0; i < pn -> children.size(); ++i) {
            if (i != pn -> children.size() - 1 && pn -> children[i]->pathToGlobal) {
                std::swap(pn -> children[i], pn -> children.back());
                break;
            }
        }
        for (const auto &c : pn->children)
            Q.push(c);
    }
    // refine the prefix tree
    root -> refine(orders);
    if (visitedPT.find(root) == visitedPT.end()) exist = false;
    else {
        exist = true;
        delete root;
        return nullptr;
    }
    root -> initPoses(orders, query, dist);
    return root;
}

PrefixNode *fullAttributeTree(const HyperTree &t, std::vector<PrefixNode *> &attributeOrder,
                              map<PrefixNode *, vector<VertexID>> &bagsBelow, bool share) {
    PrefixNode *root = new PrefixNode(99);
    root->pathToGlobal = true;
    for (VertexID nID = 0; nID < t.numNodes - 1; ++nID) bagsBelow[root].push_back(nID);
    attributeOrder.push_back(root);
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        bool global = nID == t.numNodes - 1;
        const HyperNode &bag = t.nodes[nID];
        PrefixNode *currentNode = root;
        for (int i = 0; i < bag.numAttributes; ++i) {
            VertexID u = bag.attributes[i];
            int childIndex = -1;
            if (share) childIndex = currentNode->findChildIndex(u);
            if (childIndex == -1) {
                if (global) break;
                // If the child does not exist, create a new child node
                PrefixNode *newNode = new PrefixNode(u);
                newNode->pathToGlobal = false;
                currentNode->children.push_back(newNode);
                currentNode = currentNode->children.back();
                attributeOrder.push_back(newNode);
                bagsBelow[newNode].push_back(nID);
            } else {
                // If the child exists, move to the child node
                currentNode = currentNode->children[childIndex];
                bagsBelow[currentNode].push_back(nID);
                if (global) currentNode->pathToGlobal = true;
            }
        }
    }

    return root;
}

PrefixNode *fullAttributeTree(PrefixNode *attrTree, std::vector<PrefixNode *> &attributeOrder,
                              map<PrefixNode *, vector<VertexID>> &bagsBelow, const HyperTree &t) {
    PrefixNode *root = attrTree->clone();
    PrefixNode *pn1 = attrTree;
    PrefixNode *pn2 = root;
    ui height = pn1 -> getHeight();
    std::vector<ui> childPoses(height, 0);
    std::vector<PrefixNode *> nodes1, nodes2;
    int depth = 0;
    attributeOrder.push_back(pn2);
    for (VertexID nID: pn1->nIDsToCall) {
        if (depth + 1 < t.nodes[nID].numAttributes && nID != t.numNodes - 1) {
            PrefixNode *attr = pn2;
            for (int i = 0; i < t.nodes[nID].numAttributes; ++i) {
                PrefixNode *newAttr = new PrefixNode(t.nodes[nID].attributes[i]);
                newAttr->pathToGlobal = false;
                attributeOrder.push_back(newAttr);
                attr->children.push_back(newAttr);
                attr = newAttr;
            }
            attr->nIDsToCall.push_back(nID);
        }
    }
    while (depth >= 0) {
        while (childPoses[depth] < pn1 -> children.size()) {
            PrefixNode *current1 = pn1 -> children[childPoses[depth]];
            PrefixNode *current2 = pn2 -> children[childPoses[depth]];
            ++childPoses[depth];
            nodes1.push_back(current1);
            nodes2.push_back(current2);
            bool exists = false;
            for (PrefixNode *attr : attributeOrder) {
                if (attr == current2) {
                    exists = true;
                    break;
                }
            }
            if (!exists) attributeOrder.push_back(current2);
            for (VertexID nID: current1->nIDsToCall) {
                if (depth + 1 < t.nodes[nID].numAttributes && nID != t.numNodes - 1) {
                    PrefixNode *attr = current2;
                    for (int i = depth + 1; i < t.nodes[nID].numAttributes; ++i) {
                        PrefixNode *newAttr = new PrefixNode(t.nodes[nID].attributes[i]);
                        newAttr->pathToGlobal = false;
                        attributeOrder.push_back(newAttr);
                        attr->children.push_back(newAttr);
                        attr = newAttr;
                    }
                    attr->nIDsToCall.push_back(nID);
                }
            }
            if (!current1 -> children.empty()) {
                ++depth;
                childPoses[depth] = 0;
            }
            else {
                nodes1.pop_back();
                nodes2.pop_back();
            }
            if (depth > 0) {
                pn1 = nodes1[depth - 1];
                pn2 = nodes2[depth - 1];
            }
            else {
                pn1 = attrTree;
                pn2 = root;
            }
        }
        --depth;
        if (depth >= 0) {
            nodes1.pop_back();
            nodes2.pop_back();
            if (depth == 0) {
                pn1 = attrTree;
                pn2 = root;
            }
            else {
                pn1 = nodes1[depth - 1];
                pn2 = nodes2[depth - 1];
            }
        }
    }
    std::vector<PrefixNode *> path = root->locate(t.numNodes - 1);
    PrefixNode *attr = root;
    if (!path.empty()) attr = path.back();
    for (int i = path.size(); i < t.nodes[t.numNodes - 1].numAttributes; ++i) {
        PrefixNode *newAttr = new PrefixNode(t.nodes[t.numNodes - 1].attributes[i]);
        newAttr->pathToGlobal = true;
        attributeOrder.push_back(newAttr);
        attr->children.push_back(newAttr);
        attr = newAttr;
    }
    attr->nIDsToCall.push_back(t.numNodes - 1);
    attrTree->addBagsBelow(bagsBelow);
    root->addBagsBelow(bagsBelow);

    return root;
}

std::vector<std::vector<int>>
buildAttrIDMap(const vector<PrefixNode *> &attributes, const map<PrefixNode *, vector<VertexID>> &bagsBelow,
               const HyperTree &t) {
    std::vector<std::vector<int>> attrIDMap(t.numNodes);
    for (VertexID nID = 0; nID < t.numNodes; ++nID) attrIDMap[nID].resize(t.nodes[nID].numAttributes);
    int depth = 0;
    const PrefixNode *pn = attributes[0];
    ui height = pn -> getHeight();
    std::vector<ui> childPoses(height, 0);
    std::vector<PrefixNode *> nodes;
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            ++childPoses[depth];
            nodes.push_back(current);
            int id = 0;
            for (int i = 0; i < attributes.size(); ++i) {
                if (attributes[i] == current) {
                    id = i;
                    break;
                }
            }
            for (VertexID nID : bagsBelow.at(current))
                attrIDMap[nID][depth] = id;
            if (!current -> children.empty()) {
                ++depth;
                childPoses[depth] = 0;
            }
            else nodes.pop_back();
            if (depth > 0) pn = nodes[depth - 1];
            else pn = attributes[0];
        }
        --depth;
        if (depth >= 0) {
            nodes.pop_back();
            if (depth == 0) pn = attributes[0];
            else pn = nodes[depth - 1];
        }
    }

    return attrIDMap;
}

std::vector<std::pair<int, VertexID>> buildDependent(const vector<PrefixNode *> &attributes) {
    std::vector<std::pair<int, VertexID>> depend(attributes.size());
    depend[0].first = -1;
    depend[0].second = -1;
    int depth = 0;
    const PrefixNode *pn = attributes[0];
    ui height = pn -> getHeight();
    std::vector<ui> childPoses(height, 0);
    std::vector<PrefixNode *> nodes;
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            ++childPoses[depth];
            nodes.push_back(current);
            int id = 0, pid = 0;
            for (int i = 0; i < attributes.size(); ++i) {
                if (attributes[i] == current) {
                    id = i;
                }
                if (attributes[i] == pn) {
                    pid = i;
                }
            }
            depend[id].first = pid;
            depend[id].second = -1;
            if (!current -> children.empty()) {
                ++depth;
                childPoses[depth] = 0;
            }
            else nodes.pop_back();
            if (depth > 0) pn = nodes[depth - 1];
            else pn = attributes[0];
        }
        --depth;
        if (depth >= 0) {
            nodes.pop_back();
            if (depth == 0) pn = attributes[0];
            else pn = nodes[depth - 1];
        }
    }

    return depend;
}

// Function to clone a PrefixNode tree
PrefixNode* cloneTree(const PrefixNode* root) {
    if (!root) return nullptr;
    PrefixNode* newRoot = new PrefixNode(root->u);
    for (const auto& child : root->children) {
        newRoot->children.push_back(cloneTree(child));
    }
    return newRoot;
}

std::vector<PrefixNode *> allPrefixTree(const vector<vector<VertexID>> &nodeAttrs, const vector<VertexID> &globalAttrs,
                                        const Graph &query, const std::vector<std::vector<size_t>> &dist) {
    std::vector<PrefixNode *> allTrees;
    if (globalAttrs.empty()) return allTrees;

    std::queue<std::pair<PrefixNode *, std::vector<VertexID>>> Q;

    // Initialize queue with each element as a root
    for (VertexID rootElement : globalAttrs) {
        std::vector<VertexID> remainingElements = globalAttrs;
        remainingElements.erase(std::remove(remainingElements.begin(), remainingElements.end(), rootElement), remainingElements.end());
        PrefixNode *root = new PrefixNode(rootElement);
        Q.push({root, remainingElements});
    }

    while (!Q.empty()) {
        auto [currentTree, remainingElements] = Q.front();
        Q.pop();

        if (remainingElements.empty()) {
            allTrees.push_back(currentTree);
        } else {
            for (VertexID nextElement : remainingElements) {
                PrefixNode *newTree = cloneTree(currentTree);
                PrefixNode *newNode = new PrefixNode(nextElement);

                std::queue<PrefixNode *> nodes;
                nodes.push(newTree);
                while (!nodes.empty()) {
                    PrefixNode *currentNode = nodes.front();
                    nodes.pop();
                    currentNode->children.push_back(newNode);
                    Q.push({cloneTree(newTree), remainingElements});
                    currentNode->children.pop_back();
                    for (PrefixNode *child : currentNode->children) {
                        nodes.push(child);
                    }
                }
                delete newNode;
            }
            delete currentTree;
        }
    }

    return allTrees;
}

int PrefixNode::findChildIndex(VertexID id) {
    for (int i = 0; i < children.size(); ++i) {
        if (children[i]->u == id) {
            return i;
        }
    }
    return -1;
}

ui PrefixNode::getHeight() const {
    std::queue<const PrefixNode *> Q;
    Q.push(this);
    ui height = 0;
    while (!Q.empty()) {
        ++height;
        std::queue<const PrefixNode *> newQ;
        while (!Q.empty()) {
            const PrefixNode *pn = Q.front();
            Q.pop();
            for (const auto &c : pn->children)
                newQ.push(c);
        }
        Q = newQ;
    }

    return height;
}

std::vector<PrefixNode *> PrefixNode::locate(VertexID nID) const {
    int depth = 0;
    const PrefixNode *pn = this;
    ui height = pn -> getHeight();
    std::vector<ui> childPoses(height, 0);
    std::vector<PrefixNode *> nodes;
    for (VertexID nID2: pn -> nIDsToCall) {
        if (nID2 == nID) return nodes;
    }
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            ++childPoses[depth];
            nodes.push_back(current);
            for (VertexID nID2: current -> nIDsToCall) {
                if (nID2 == nID) return nodes;
            }
            if (!current -> children.empty()) {
                ++depth;
                childPoses[depth] = 0;
            }
            else nodes.pop_back();
            if (depth > 0) pn = nodes[depth - 1];
            else pn = this;
        }
        --depth;
        nodes.pop_back();
        if (depth >= 0) {
            if (depth == 0) pn = this;
            else pn = nodes[depth - 1];
        }
    }

    return nodes;
}

void PrefixNode::locate(PrefixNode *attribute, std::vector<PrefixNode *> &path, std::vector<ui> &childPoses) {
    int depth = 0;
    PrefixNode *pn = this;
    ui height = pn -> getHeight();
    childPoses = std::vector<ui>(height, 0);
    if (pn == attribute) return;
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            path.push_back(current);
            if (current == attribute) {
                return;
            }
            if (!current -> children.empty()) {
                ++depth;
                childPoses[depth] = 0;
            }
            else {
                ++childPoses[depth];
                path.pop_back();
            }
            if (depth > 0) pn = path[depth - 1];
            else pn = this;
        }
        --depth;
        if (depth >= 0) {
            ++childPoses[depth];
            path.pop_back();
            if (depth == 0) pn = this;
            else pn = path[depth - 1];
        }
    }
}

int PrefixNode::locate(PrefixNode *attribute, VertexID firstID){
    int depth = 0;
    PrefixNode *pn = this;
    ui height = pn -> getHeight();
    std::vector<ui> childPoses = std::vector<ui>(height, 0);
    std::vector<PrefixNode *> path;
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            path.push_back(current);
            if (current->u == attribute->u) {
                std::vector<VertexID> below = current->getBagsBelow();
                if (std::find(below.begin(), below.end(), firstID) != below.end())
                    return depth;
            }
            if (!current -> children.empty()) {
                ++depth;
                childPoses[depth] = 0;
            }
            else {
                ++childPoses[depth];
                path.pop_back();
            }
            if (depth > 0) pn = path[depth - 1];
            else pn = this;
        }
        --depth;
        if (depth >= 0) {
            ++childPoses[depth];
            path.pop_back();
            if (depth == 0) pn = this;
            else pn = path[depth - 1];
        }
    }

    return -1;
}

ui PrefixNode::numAttr(const HyperTree &t) const {
    ui num = 0;
    std::queue<const PrefixNode *> Q;
    std::queue<int> heights;
    Q.push(this);
    heights.push(0);
    while (!Q.empty()) {
        const PrefixNode *pn = Q.front();
        int height = heights.front();
        Q.pop();
        heights.pop();
        ++num;
        for (VertexID nID : pn -> nIDsToCall)
            num += t.nodes[nID].numAttributes - height;
        for (const auto &c : pn->children) {
            Q.push(c);
            heights.push(height + 1);
        }
    }

    return num - 1;
}

void PrefixNode::refine(const std::vector<std::vector<VertexID>> &sharedAttrs) {
    std::queue<PrefixNode *> Q;
    Q.push(this);
    while (!Q.empty()) {
        PrefixNode *pn = Q.front();
        Q.pop();
        std::vector<VertexID> newNID;
        std::set<VertexID> childCallNID;
        std::set<VertexID> buildNID;
        std::vector<int> validPoses;
        std::vector<PrefixNode *> newChildren;
        if (pn -> pathToGlobal) {
            for (VertexID nID = 0; nID < sharedAttrs.size() - 1; ++nID) {
                if (std::find(pn->nIDsToCall.begin(), pn->nIDsToCall.end(), nID) == pn->nIDsToCall.end()
                    && std::find(sharedAttrs[nID].begin(), sharedAttrs[nID].end(), pn->u) != sharedAttrs[nID].end()) {
                    pn->nIDsToJoin.push_back(nID);
                }
            }
        }
        for (int i = 0 ; i < pn->children.size(); ++i) {
            PrefixNode *c = pn->children[i];
            if (c->nIDsToCall.size() > 1) {
                Q.push(c);
                childCallNID.insert(c->nIDsToCall.begin(), c->nIDsToCall.end());
                validPoses.push_back(i);
                if (pn -> pathToGlobal && !c -> pathToGlobal) buildNID.insert(c->nIDsToCall.begin(), c->nIDsToCall.end());
            }
            else {
                delete c;
            }
        }
        for (auto pos: validPoses) {
            newChildren.push_back(pn -> children[pos]);
        }
        for (VertexID nID: pn -> nIDsToCall) {
            if (childCallNID.find(nID) == childCallNID.end())
                newNID.push_back(nID);
        }
        pn -> nIDsToCall = newNID;
        pn -> children = newChildren;
        if (pn -> pathToGlobal) {
            for (VertexID nID: pn -> nIDsToCall) {
                if (nID != sharedAttrs.size() - 1)
                    pn -> nIDsToBuild.push_back(nID);
            }
            for (VertexID nID : buildNID)
                pn -> nIDsToBuild.push_back(nID);
        }
    }
}

void PrefixNode::initPoses(const vector<std::vector<VertexID>> &bagAttrs, const Graph &query,
                           const vector<std::vector<size_t>> &dist) {
    std::vector<VertexID> attributes(query.getNumVertices());
    std::vector<PrefixNode *> nodes(query.getNumVertices(), nullptr);
    int depth = 0;
    std::vector<ui> childPoses(query.getNumVertices(), 0);
    PrefixNode *pn = this;
    // rebuild nIDsToJoin
    std::vector<VertexID> materializedNodes = pn -> nIDsToBuild;
    while (!pn->children.empty()) {
        PrefixNode *current = pn -> children.back();
        if (!(current -> pathToGlobal)) break;
        current -> nIDsToJoin.clear();
        for (VertexID nID : materializedNodes) {
            if (std::find(bagAttrs[nID].begin(), bagAttrs[nID].end(), current->u) != bagAttrs[nID].end()) {
                current -> nIDsToJoin.push_back(nID);
            }
        }
        for (VertexID nID : current -> nIDsToBuild) {
            materializedNodes.push_back(nID);
        }
        pn = current;
    }
    depth = 0;
    pn = this;
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            ++childPoses[depth];
            VertexID u1 = current -> u;
            attributes[depth] = u1;
            current -> attributesBefore.clear();
            for (int i = 0; i < depth; ++i) {
                VertexID u2 = attributes[i];
                if (query.getEdgeID(u1, u2) != -1)
                    current -> attributesBefore.push_back(u2);
            }
            if (current -> pathToGlobal) {
                std::set<VertexID> coveredAttrs;
                for (VertexID nID : current -> nIDsToJoin) {
                    for (VertexID u2 : bagAttrs[nID])
                        coveredAttrs.insert(u2);
                }
                std::vector<ui> attrsToJoin;
                for (VertexID u2: current -> attributesBefore) {
                    if (coveredAttrs.find(u2) == coveredAttrs.end())
                        attrsToJoin.push_back(u2);
                }
                current -> attributesBefore = attrsToJoin;
            }
            bool cartesian = (depth > 0 && current -> attributesBefore.empty() && current -> nIDsToJoin.empty());
            if (cartesian) {
                size_t minDist = std::numeric_limits<size_t>::max();
                for (VertexID i = 0; i < depth; ++i) {
                    VertexID u2 = attributes[i];
                    if (dist[u1][u2] < minDist) {
                        minDist = dist[u1][u2];
                        current -> cartesianParent = u2;
                    }
                }
            }
            nodes[depth] = current;
            if (!current -> children.empty()) {
                ++depth;
                childPoses[depth] = 0;
                pn = nodes[depth - 1];
            }
        }
        --depth;
        if (depth == 0) pn = this;
        else if (depth > 0) pn = nodes[depth - 1];
    }
}

void PrefixNode::print() const {
    std::queue<const PrefixNode *> Q;
    Q.push(this);
    while (!Q.empty()) {
        const PrefixNode *pn = Q.front();
        Q.pop();
        if (!pn -> children.empty()) std::cout << pn -> u << " -- ";
        for (int i = 0; i < pn -> children.size(); ++i) {
            std::cout << " " << pn -> children[i] -> u << " ";
        }
        if (!pn -> children.empty()) std::cout << std::endl;
        for (const auto &c : pn->children)
            Q.push(c);
    }
}

bool PrefixNode::operator==(const PrefixNode &rhs) const {
    if (children.size() != rhs.children.size()) return false;
    for (int i = 0; i < children.size(); ++i) {
        if (*children[i] != *rhs.children[i]) return false;
    }
    return u == rhs.u && nIDsToCall == rhs.nIDsToCall;
}

bool PrefixNode::operator!=(const PrefixNode &rhs) const {
    return !(rhs == *this);
}

void PrefixNode::checkCallAll(ui numNodes) const {
    std::set<VertexID> nIDs;
    std::queue<const PrefixNode *> Q;
    Q.push(this);
    while (!Q.empty()) {
        const PrefixNode *pn = Q.front();
        Q.pop();
        nIDs.insert(pn -> nIDsToCall.begin(), pn->nIDsToCall.end());
        for (const auto &c : pn->children)
            Q.push(c);
    }
//    assert(nIDs.size() == numNodes);
}

std::vector<VertexID> PrefixNode::getBagsBelow() const {
    std::vector<VertexID> nIDs;
    std::stack<const PrefixNode *> S;
    S.push(this);
    while (!S.empty()) {
        const PrefixNode *pn = S.top();
        S.pop();
        for (VertexID nID : pn -> nIDsToCall)
            nIDs.push_back(nID);
        for (auto it = pn->children.rbegin(); it != pn->children.rend(); ++it)
            S.push(*it);
    }

    return nIDs;
}

PrefixNode *PrefixNode::clone() const {
    PrefixNode* newNode = new PrefixNode(u);
    newNode->attributesBefore = attributesBefore;
    newNode->cartesianParent = cartesianParent;
    newNode->nIDsToCall = nIDsToCall;
    newNode->nIDsToJoin = nIDsToJoin;
    newNode->nIDsToBuild = nIDsToBuild;
    newNode->pathToGlobal = pathToGlobal;

    for (const auto& child : children) {
        newNode->children.push_back(child->clone());
    }

    return newNode;
}

void PrefixNode::initNIDsToBuild(ui numNodes) {
    PrefixNode *pn = this;
    // find the path to global
    std::queue<PrefixNode *> q1;
    q1.push(pn);
    ui height = getHeight();
    std::vector<PrefixNode *> path, nodes(height, nullptr);
    int depth = 0;
    std::vector<ui> childPoses(height, 0);
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            ++childPoses[depth];
            nodes[depth] = current;
            current -> pathToGlobal = false;
            current -> nIDsToBuild.clear();
            current -> nIDsToJoin.clear();
            for (VertexID nID: current -> nIDsToCall) {
                if (nID == numNodes - 1) {
                    path.assign(nodes.begin(), nodes.begin() + depth + 1);
                    break;
                }
            }
            if (!current -> children.empty()) {
                ++depth;
                childPoses[depth] = 0;
            }
            if (depth > 0) pn = nodes[depth - 1];
        }
        --depth;
        if (depth >= 0) {
            if (depth == 0) pn = this;
            else pn = nodes[depth - 1];
        }
    }
    this -> pathToGlobal = true;
    for (PrefixNode *node : path) node -> pathToGlobal = true;
    pn = this;
    std::queue<PrefixNode *> Q;
    Q.push(this);
    while (!Q.empty()) {
        pn = Q.front();
        Q.pop();
        for (int i = 0; i < pn -> children.size(); ++i) {
            if (i != pn -> children.size() - 1 && pn -> children[i]->pathToGlobal) {
                std::swap(pn -> children[i], pn -> children.back());
                break;
            }
        }
        for (const auto &c : pn->children)
            Q.push(c);
    }
    pn = this;
    while (true) {
        pn -> nIDsToBuild.clear();
        for (VertexID nID : pn -> nIDsToCall)
            if (nID != numNodes - 1)
                pn -> nIDsToBuild.push_back(nID);
        std::queue<PrefixNode *>q;
        for (PrefixNode *c : pn -> children) {
            if (!c->pathToGlobal)
                q.push(c);
        }
        while (!q.empty()) {
            PrefixNode *current = q.front();
            q.pop();
            for (VertexID nID : current -> nIDsToCall) pn -> nIDsToBuild.push_back(nID);
            for (PrefixNode *c : current -> children) q.push(c);
        }
        std::sort(pn -> nIDsToBuild.begin(), pn->nIDsToBuild.end());
        if (pn -> children.empty() || !pn->children.back()->pathToGlobal) break;
        else pn = pn -> children.back();
    }
}

// check whether the left subtree can be merged to the right path
void PrefixNode::mergeToRight(std::vector<std::vector<VertexID>> &localOrders, std::vector<VertexID> &rightCall) {
    if (this -> children.empty()) return;
    PrefixNode *pn = this;
    PrefixNode *right = this -> children.back();
    // merge bags to the right
    std::vector<VertexID> nIDsToRemove;
    std::sort(nIDsToCall.begin(), nIDsToCall.end());
    for (int i = 0; i < nIDsToCall.size(); ++i) {
        VertexID nID = nIDsToCall[i];
        if (localOrders[nID][0] == right->u) {
            nIDsToRemove.push_back(nID);
            rightCall.push_back(nID);
            std::sort(rightCall.begin(), rightCall.end());
        }
    }
    std::vector<VertexID> old = nIDsToCall;
    nIDsToCall.clear();
    std::set_difference(old.begin(), old.end(), nIDsToRemove.begin(), nIDsToRemove.end(), std::back_inserter(nIDsToCall));
    // merge subtree to the right
    int pos = -1;
    for (int i = 0; i < children.size() - 1; ++i) {
        if (children[i] -> u == right->u) {
            pos = i;
            break;
        }
    }
    if (pos == -1) return;
    PrefixNode *left = this -> children[pos];
    std::queue<PrefixNode *> q;
    q.push(left);
    while (!q.empty()) {
        PrefixNode *subtree = q.front();
        q.pop();
        for (VertexID nID: subtree->nIDsToCall) rightCall.push_back(nID);
        for (PrefixNode *c: subtree->children) q.push(c);
    }
    std::sort(rightCall.begin(), rightCall.end());
    std::vector<PrefixNode *> newChild;
    for (int i = 0; i < children.size(); ++i) {
        if (i != pos) newChild.push_back(children[i]);
    }
    children = newChild;
}

void PrefixNode::getTraverseOrder(vector<PrefixNode *> &attributes, vector<VertexID> &nIDs, const HyperTree &t) {
    int depth = 0;
    const PrefixNode *pn = this;
    ui height = pn -> getHeight();
    std::vector<ui> childPoses(height, 0);
    std::vector<PrefixNode *> nodes;
    for (VertexID nID2: pn -> nIDsToCall) {
        if (nID2 !=t.numNodes - 1)
            nIDs.push_back(nID2);
    }
    attributes.push_back(this);
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            attributes.push_back(current);
            ++childPoses[depth];
            nodes.push_back(current);
            for (VertexID nID2: current -> nIDsToCall) {
                if (nID2 !=t.numNodes - 1) nIDs.push_back(nID2);
            }
            if (!current -> children.empty()) {
                ++depth;
                childPoses[depth] = 0;
            }
            else nodes.pop_back();
            if (depth > 0) pn = nodes[depth - 1];
            else pn = this;
        }
        --depth;
        if (depth >= 0) {
            nodes.pop_back();
            if (depth == 0) pn = this;
            else pn = nodes[depth - 1];
        }
    }
    nIDs.push_back(t.numNodes - 1);
}

void PrefixNode::addBagsBelow(map<PrefixNode *, std::vector<VertexID>> &bagsBelow) {
    int depth = 0;
    PrefixNode *pn = this;
    ui height = pn -> getHeight();
    std::vector<ui> childPoses(height, 0);
    std::vector<PrefixNode *> nodes;
    bagsBelow[pn] = pn->getBagsBelow();
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            bagsBelow[current] = current->getBagsBelow();
            ++childPoses[depth];
            nodes.push_back(current);
            if (!current -> children.empty()) {
                ++depth;
                childPoses[depth] = 0;
            }
            else nodes.pop_back();
            if (depth > 0) pn = nodes[depth - 1];
            else pn = this;
        }
        --depth;
        if (depth >= 0) {
            nodes.pop_back();
            if (depth == 0) pn = this;
            else pn = nodes[depth - 1];
        }
    }
}

HyperTree twoBagSubQuery(const HyperTree &t, const Graph &query) {
    HyperTree twoBagT;
    twoBagT.numNodes = 2;
    twoBagT.nodes = new HyperNode[2];
    twoBagT.edges = std::vector<std::vector<VertexID>>(2);
    twoBagT.edges[0].push_back(1);
    twoBagT.edges[1].push_back(0);
    // select the pair of bags with the maximum number of shared attributes in t
    ui maxShared = 0;
    ui pos1, pos2;
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        std::set<VertexID> vertices1(t.nodes[nID].attributes, t.nodes[nID].attributes + t.nodes[nID].numAttributes);
        for (VertexID nID2: t.edges[nID]) {
            if (nID2 < nID) continue;
            ui shared = 0;
            for (int i = 0; i < t.nodes[nID2].numAttributes; ++i) {
                if (vertices1.find(t.nodes[nID2].attributes[i]) != vertices1.end()) ++shared;
            }
            if (shared > maxShared) {
                maxShared = shared;
                pos1 = nID;
                pos2 = nID2;
            }
        }
    }

    twoBagT.nodes[0] = t.nodes[pos1];
    twoBagT.nodes[1] = t.nodes[pos2];
    twoBagT.numAttributes = query.getNumVertices();
    twoBagT.nodes[0].attributes = new VertexID [t.nodes[pos1].numAttributes];
    memcpy(twoBagT.nodes[0].attributes, t.nodes[pos1].attributes, t.nodes[pos1].numAttributes * sizeof(VertexID));
    twoBagT.nodes[1].attributes = new VertexID [t.nodes[pos2].numAttributes];
    memcpy(twoBagT.nodes[1].attributes, t.nodes[pos2].attributes, t.nodes[pos2].numAttributes * sizeof(VertexID));
    twoBagT.v2n = new std::vector<VertexID> [query.getNumVertices()];
    for (int i = 0; i < twoBagT.nodes[0].numAttributes; ++i) {
        VertexID u = twoBagT.nodes[0].attributes[i];
        twoBagT.v2n[u].push_back(0);
    }
    for (int i = 0; i < twoBagT.nodes[1].numAttributes; ++i) {
        VertexID u = twoBagT.nodes[1].attributes[i];
        twoBagT.v2n[u].push_back(1);
    }

    return twoBagT;
}

void buildFromPrefixTree(PrefixNode *prefixTree, const std::vector<std::vector<VertexID>> &nodeOrder, HyperTree &t,
                         const std::vector<std::vector<VertexID>> &sharedAttrs, const Graph &query, CandidateSpace &cs) {
    PrefixNode *pn = prefixTree;
    ui height = pn -> getHeight();
    std::vector<ui> childPoses(height, 0);
    int depth = 0;
    std::vector<PrefixNode *> nodes(query.getNumVertices(), nullptr);
    t.globalOrder.clear();
    for (VertexID nID = 0; nID < t.numNodes; ++nID) t.nodes[nID].prefixSize = 0;
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            ++childPoses[depth];
            VertexID u = current -> u;
            nodes[depth] = current;
            if (current -> pathToGlobal && u != 99) t.globalOrder.push_back(u);
            for (VertexID nID: current -> nIDsToCall) {
                t.nodes[nID].numAttributes = nodeOrder[nID].size();
                t.nodes[nID].attributes = new VertexID [nodeOrder[nID].size()];
                for (int j = 0; j <= depth; ++j) {
                    if (nodes[j]->pathToGlobal) {
                        VertexID u2 = nodes[j] -> u;
                        if (std::find(t.v2n[u2].begin(), t.v2n[u2].end(), nID) != t.v2n[u2].end())
                            ++t.nodes[nID].prefixSize;
                    }
                    else break;
                }
                t.nodes[nID].prefix = new VertexID [t.nodes[nID].prefixSize];
                for (int j = 0; j < t.nodes[nID].prefixSize; ++j) {
                    t.nodes[nID].prefix[j] = nodeOrder[nID][j];
                }
                for (int j = 0; j < t.nodes[nID].numAttributes; ++j) {
                    t.nodes[nID].attributes[j] = nodeOrder[nID][j];
                }
                t.nodes[nID].initPoses(sharedAttrs, query, cs.dist, nID == t.numNodes - 1);
            }
            if (!current -> children.empty()) {
                ++depth;
                childPoses[depth] = 0;
            }
            if (depth > 0) pn = nodes[depth - 1];
        }
        --depth;
        if (depth >= 0) {
            if (depth == 0) pn = prefixTree;
            else pn = nodes[depth - 1];
        }
    }
    for (VertexID nID: prefixTree -> nIDsToCall) {
        t.nodes[nID].numAttributes = nodeOrder[nID].size();
        t.nodes[nID].attributes = new VertexID [nodeOrder[nID].size()];
        t.nodes[nID].prefixSize = 0;
        for (int j = 0; j < t.nodes[nID].numAttributes; ++j) {
            t.nodes[nID].attributes[j] = nodeOrder[nID][j];
        }
        t.nodes[nID].initPoses(sharedAttrs, query, cs.dist, nID == t.numNodes - 1);
    }
    for (int i = 0; i < nodeOrder.back().size(); ++i) {
        VertexID u = nodeOrder.back()[i];
        if (std::find(t.globalOrder.begin(), t.globalOrder.end(), u) == t.globalOrder.end())
            t.globalOrder.push_back(u);
    }
    for (VertexID nID = 0; nID < t.numNodes - 1; ++nID) {
        for (int j = 0; j < nodeOrder[nID].size(); ++j) {
            VertexID u = nodeOrder[nID][j];
            if (std::find(t.globalOrder.begin(), t.globalOrder.end(), u) == t.globalOrder.end())
                t.globalOrder.push_back(u);
        }
    }
    t.extendLevel = t.nodes[t.numNodes - 1].prefixSize;
    t.nIDs = std::vector<std::vector<VertexID>>(t.globalOrder.size());
    for (int j = 0; j < t.globalOrder.size(); ++j) {
        VertexID u = t.globalOrder[j];
        for (VertexID nID = 0; nID < nodeOrder.size(); ++nID) {
            if (std::find(nodeOrder[nID].begin(), nodeOrder[nID].end(), u) != nodeOrder[nID].end())
                t.nIDs[j].push_back(nID);
        }
    }
    t.buildTrieOrder();
}

void buildFromPrefixTree(const std::vector<PrefixNode *> &prefixTrees,
                         const std::vector<std::vector<std::vector<VertexID>>> &bestOrders,
                         std::vector<HyperTree> &trees, const HyperTree &reference,
                         const std::vector<std::vector<VertexID>> &sharedAttrs, const Graph &query,
                         CandidateSpace &cs) {
    trees.resize(prefixTrees.size());
    for (int i = 0; i < prefixTrees.size(); ++i) {
        HyperTree &t = trees[i];
        t.numAttributes = reference.numAttributes;
        t.numNodes = reference.numNodes;
        t.nodes = new HyperNode[t.numNodes];
        t.v2n = new std::vector<VertexID> [t.numAttributes];
        for (VertexID u = 0; u < t.numAttributes; ++u)
            t.v2n[u] = reference.v2n[u];
        buildFromPrefixTree(prefixTrees[i], bestOrders[i], t, sharedAttrs, query, cs);
    }
}

void removeNID(VertexID nID, PrefixNode *root) {
    for (auto it = root -> nIDsToBuild.begin(); it != root -> nIDsToBuild.end();) {
        if (*it == nID) {
            root -> nIDsToBuild.erase(it);
            break;
        }
        else ++it;
    }
    for (auto it = root -> nIDsToCall.begin(); it != root -> nIDsToCall.end();) {
        if (*it == nID) {
            root -> nIDsToCall.erase(it);
            break;
        }
        else ++it;
    }
}

// remove nID in the branch: the nIDsToBuild for the global node, and the nIDsToCall for the current node
// for nodes in the path, if there is no bags to call, it is also removed
void removeNID(VertexID nID, std::vector<PrefixNode *> &nodes, int depth, PrefixNode *&root) {
    int pos = 0;
    for (; pos < depth + 1; ++pos) {
        if (!nodes[pos]->pathToGlobal)
            break;
    }
    PrefixNode *pn = root;
    if (pos != 0 ) pn = nodes[pos - 1];
    PrefixNode *current = nodes[depth];
    for (auto it = pn -> nIDsToBuild.begin(); it != pn -> nIDsToBuild.end();) {
        if (*it == nID) {
            pn -> nIDsToBuild.erase(it);
            break;
        }
        else ++it;
    }
    for (auto it = current -> nIDsToCall.begin(); it != current -> nIDsToCall.end();) {
        if (*it == nID) {
            current -> nIDsToCall.erase(it);
            break;
        }
        else ++it;
    }
    if (current -> children.empty() && current -> nIDsToCall.empty()) {
        pos = depth - 1;
        PrefixNode *child = current;
        current = root;
        if (pos != -1) current = nodes[pos];
        for (auto it = current -> children.begin(); it != current -> children.end(); ++it) {
            if (*it == child) {
                current -> children.erase(it);
                break;
            }
        }
        while (current -> children.empty() && current -> nIDsToCall.empty()) {
            child = current;
            --pos;
            current = root;
            if (pos != -1) current = nodes[pos];
            for (auto it = current -> children.begin(); it != current -> children.end(); ++it) {
                if (*it == child) {
                    current -> children.erase(it);
                    break;
                }
            }
        }
    }
}

std::vector<int> getMappingSizes(const HyperTree &t, const PrefixNode *pt) {
    std::vector<int> mappingSizes(t.numNodes);
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
                for (VertexID u2: attrsInPath) {
                    if (std::find(t.v2n[u2].begin(), t.v2n[u2].end(), nID) != t.v2n[u2].end())
                        ++mappingSizes[nID];
                }
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

    return mappingSizes;
}

// Function to convert a subset to an ID
uint64_t getSubsetID(const std::vector<VertexID>& subset) {
    uint64_t id = 0;
    for (VertexID vertex : subset) {
        id |= (1ULL << vertex);
    }
    return id;
}

bool subsetConnectivity(const Graph &query, const std::vector<uint64_t> &cc, const vector<VertexID> &subset) {
    if (subset.empty()) return true;
    std::vector<std::vector<VertexID>> subsetCC;
    query.computeConnectedComponents(subset, subsetCC);
    if (subsetCC.size() == 1) return true;
    else return false;
//    ui num = cc.size();
//    for (int i = 0; i < subsetCC.size(); ++i) {
//        for (int j = 0; j < cc.size(); ++j) {
//            if (getSubsetID(subsetCC[i]) == cc[j]) --num;
//        }
//    }
//    return num <= cc.size() - 1;
}

bool orderConnectivity(const Graph &query, const vector<VertexID> &order, const std::vector<uint64_t> &cc) {
    if (order.empty()) return false;
    uint64_t ccID = 1 << order[0];
    for (int i = 1; i < order.size(); ++i) {
        bool connected = false;
        VertexID u1 = order[i];
        for (int j = 0; j < i; ++j) {
            VertexID u2 = order[j];
            if (query.getEdgeID(u1, u2) != -1) {
                connected = true;
                ccID += 1 << u1;
                break;
            }
        }
        if (!connected) {
            bool valid = false;
            for (auto queryCC: cc) {
                if (queryCC == ccID) {
                    valid = true;
                    break;
                }
            }
            if (!valid) return false;
            ccID = 1 << u1;
        }
    }
    return true;
}

bool orderConnectivity(const Graph &query, const vector<VertexID> &order) {
    for (int i = 1; i < order.size(); ++i) {
        VertexID u1 = order[i];
        bool connected = false;
        for (int j = 0; j < i; ++j) {
            VertexID u2 = order[j];
            if (query.getEdgeID(u1, u2) != -1) {
                connected = true;
                break;
            }
        }
        if (!connected) return false;
    }

    return true;
}

void
genAllPrefixTree(const HyperTree &t, const Graph &query, CandidateSpace &cs, std::vector<PrefixNode *> &prefixTrees,
                 bool global) {
    std::vector<ui> poses(query.getNumVertices(), 0);
    std::vector<std::vector<VertexID>> v2n(query.getNumVertices());
    ui numNodes = t.numNodes;
    if (t.newGlobalNode && !global) --numNodes;
    std::vector<std::vector<VertexID>> noShareOrders(numNodes);
    std::vector<std::vector<ui>> numBackWards(numNodes);
    for (VertexID nID = 0; nID < numNodes; ++nID) {
        for (int i = 0; i < t.nodes[nID].numAttributes; ++i) {
            VertexID u = t.nodes[nID].attributes[i];
            v2n[u].push_back(nID);
            noShareOrders[nID].push_back(u);
        }
    }
    std::vector<std::vector<VertexID>> sharedAttrs(t.numNodes);
    std::vector<VertexID> internalAttrs;
    std::vector<std::vector<VertexID>> attrIntersections(t.numNodes * t.numNodes);
    for (VertexID nID = 0; nID < numNodes; ++nID) {
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
    std::vector<std::vector<VertexID>> components;
    std::vector<uint64_t> cc;
    query.computeConnectedComponents(internalAttrs, components);
    for (auto &c : components) cc.push_back(getSubsetID(c));
    for (VertexID nID1 = 0; nID1 < numNodes; ++nID1) {
        for (VertexID nID2 = nID1 + 1; nID2 < numNodes; ++nID2) {
            std::set_intersection(sharedAttrs[nID1].begin(), sharedAttrs[nID1].end(), sharedAttrs[nID2].begin(),
                                  sharedAttrs[nID2].end(), std::back_inserter(attrIntersections[nID1 * t.numNodes + nID2]));
        }
    }
    std::unordered_set<PrefixNode*, PrefixNodePtrHash, PrefixNodePtrEqual> visitedPT;
    PrefixNode *empty = new PrefixNode(99);
    for (VertexID nID = 0; nID < numNodes; ++nID) empty->nIDsToCall.push_back(nID);
    std::queue<PrefixNode *> plans;
    plans.push(empty);
    visitedPT.insert(empty);
    while (!plans.empty()) {
        PrefixNode *pt = plans.front();
        plans.pop();
        prefixTrees.push_back(pt);
        PrefixNode *pn = pt;
        ui height = pn -> getHeight();
        std::vector<ui> childPoses(height, 0);
        std::vector<PrefixNode *> nodes(height, nullptr);
        int depth = 0;
        std::vector<VertexID> attrsInPath;
        std::vector<VertexID> siblingAttr;
        extendPrefixTree(query, t, plans, pt, pn, attrsInPath, siblingAttr, -1, childPoses,
                         visitedPT, attrIntersections, sharedAttrs, cc);
        while (depth >= 0) {
            while (childPoses[depth] < pn -> children.size()) {
                PrefixNode *current = pn -> children[childPoses[depth]];
                VertexID u = current -> u;
                attrsInPath.push_back(u);
                nodes[depth] = current;
                siblingAttr.clear();
                if (depth == 0) {
                    for (PrefixNode *sibling: pt -> children)
                        if (sibling != current) siblingAttr.push_back(sibling -> u);
                }
                else {
                    for (PrefixNode *sibling: nodes[depth - 1] -> children)
                        if (sibling != current) siblingAttr.push_back(sibling -> u);
                }
                ++childPoses[depth];
                extendPrefixTree(query, t, plans, pt, current, attrsInPath, siblingAttr, depth, childPoses,
                                 visitedPT, attrIntersections, sharedAttrs, cc);
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
    }
    for (PrefixNode *pt : prefixTrees) {
        pt->initNIDsToBuild(t.numNodes);
        pt ->initPoses(sharedAttrs, query, cs.dist);
    }
}

void extendPrefixTree(const Graph &query, const HyperTree &t, queue<PrefixNode *> &plans, const PrefixNode *pt,
                      const PrefixNode *current, const std::vector<VertexID> &attrsInPath,
                      const std::vector<VertexID> &siblingAttr, int depth, std::vector<ui> &childPoses,
                      std::unordered_set<PrefixNode *, PrefixNodePtrHash, PrefixNodePtrEqual> &visitedPT,
                      const std::vector<std::vector<VertexID>> &attrIntersections,
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
                bool connected = prefix.empty();
                for (VertexID up: prefix) {
                    if (query.getEdgeID(up, u2) != -1) {
                        connected = true;
                        break;
                    }
                }
                if (!connected) continue;
                prefix.push_back(u2);
                PrefixNode *newPT = pt -> clone();
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
                else {
                    visitedPT.insert(newPT);
                    plans.push(newPT);
                }
            }
        }
    }
    if (!current -> children.empty()) {
        for (int i = 0; i < current -> children.size(); ++i) {
            PrefixNode *c = current -> children[i];
            VertexID u2 = c -> u;
            for (VertexID nID : current -> nIDsToCall) {
                if (std::find(sharedAttrs[nID].begin(), sharedAttrs[nID].end(), u2) == sharedAttrs[nID].end()) continue;
                PrefixNode *newPT = pt -> clone();
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
                else {
                    visitedPT.insert(newPT);
                    plans.push(newPT);
                }
            }
        }
    }
}

std::vector<VertexID> globalOrder(const Graph &query, const HyperTree &t, const CandidateSpace &cs, const std::vector<VertexID> &prefix) {
    std::vector<VertexID> localOrder;
    std::vector<ui> repetitions(query.getNumVertices(), 0);
    std::vector<ui> sizes(query.getNumVertices(), 0);
    const HyperNode &global = t.nodes[t.numNodes - 1];
    for (int i = 0; i < global.numAttributes; ++i) {
        VertexID u = global.attributes[i];
        if (std::find(prefix.begin(), prefix.end(), u) == prefix.end())
            localOrder.push_back(u);
    }
    for (VertexID u : localOrder) {
        repetitions[u] = t.v2n[u].size();
        sizes[u] = cs.candidateSet[u].size();
    }

    std::sort(localOrder.begin(), localOrder.end(), [&repetitions, &sizes](const VertexID &a, const VertexID &b) {
        if (repetitions[a] > repetitions[b]) return true;
        if (repetitions[a] == repetitions[b] && sizes[a] < sizes[b]) return true;
        return false;
    });
    std::vector<VertexID> totalOrder = prefix;
    for (VertexID u : localOrder) totalOrder.push_back(u);
    return totalOrder;
}

void reorderBags(const Graph &query, HyperTree &t) {
    t.numAttributes = query.getNumVertices();
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
    delete[] t.v2n;
    t.v2n = new std::vector<VertexID> [query.getNumVertices()];
    for (ui i = 0; i < t.numNodes; ++i) {
        for (ui j = 0; j < t.nodes[i].numAttributes; ++j) {
            t.v2n[t.nodes[i].attributes[j]].push_back(i);
        }
    }
}

void changeNIDsToCall(const Graph &query, HyperTree &t, PrefixNode *pt) {
    int depth = 0;
    PrefixNode *pn = pt;
    ui height = pn -> getHeight();
    std::vector<ui> childPoses(height, 0);
    std::vector<PrefixNode *> nodes;
    std::map<PrefixNode *, std::vector<PrefixNode *>> newChild;
    // root
    for (VertexID nID: pn ->nIDsToCall) {
        if (!t.newGlobalNode && nID == t.numNodes - 1) continue;
        VertexID u = t.nodes[nID].attributes[0];
        PrefixNode *child = new PrefixNode(u);
        child->pathToGlobal = false;
        child->attributesBefore = t.nodes[nID].attributesBefore[0];
        child->nIDsToCall = {nID};
        newChild[pn].push_back(child);
    }
    while (depth >= 0) {
        while (childPoses[depth] < pn -> children.size()) {
            PrefixNode *current = pn -> children[childPoses[depth]];
            ++childPoses[depth];
            nodes.push_back(current);
            for (VertexID nID: current -> nIDsToCall) {
                if (!t.newGlobalNode && nID == t.numNodes - 1) continue;
                VertexID u = t.nodes[nID].attributes[depth + 1];
                PrefixNode *child = new PrefixNode(u);
                child->pathToGlobal = false;
                child->attributesBefore = t.nodes[nID].attributesBefore[depth + 1];
                child->nIDsToCall = {nID};
                newChild[current].push_back(child);
            }
            if (!current -> children.empty()) {
                ++depth;
                childPoses[depth] = 0;
            }
            else nodes.pop_back();
            if (depth > 0) pn = nodes[depth - 1];
            else pn = pt;
        }
        --depth;
        if (depth >= 0) {
            nodes.pop_back();
            if (depth == 0) pn = pt;
            else pn = nodes[depth - 1];
        }
    }
    for (auto &item : newChild) {
        PrefixNode *attr = item.first;
        std::vector<PrefixNode *> bags = item.second;
        std::vector<PrefixNode *> newChild = bags;
        for (PrefixNode *c: attr->children) newChild.push_back(c);
        attr->children = newChild;
        if (!t.newGlobalNode && std::find(attr->nIDsToCall.begin(), attr->nIDsToCall.end(),
                                         t.numNodes - 1) != attr->nIDsToCall.end()) {
            attr->nIDsToCall = {t.numNodes - 1};
        }
        else attr->nIDsToCall.clear();
    }
}

void addGlobal(const Graph &query, HyperTree &t, PrefixNode *pt, VertexID nID, VertexID maxCostID, CandidateSpace &cs,
               vector<vector<VertexID>> &nodeOrders, std::vector<ui> &prefixSizes, bool globalOrderShare) {
    std::vector<PrefixNode *> attributes2 = pt->locate(maxCostID);
    std::vector<PrefixNode *> attributes;
    for (PrefixNode *attribute: attributes2) {
        if (t.v2n[attribute->u].size() >= 2) attributes.push_back(attribute);
        else break;
    }
//    while (attributes.back()->nIDsToCall.size() + attributes.back()->children.size() < 2) {
//        attributes.pop_back();
//    }
    std::vector<VertexID> prefix;
    for (int i = 0; i < attributes.size(); ++i) prefix.push_back(attributes[i]->u);
    std::vector<VertexID> newShareAttr;
    int pos = attributes.size();
    for (; pos < nodeOrders[maxCostID].size(); ++pos) {
        VertexID u = nodeOrders[maxCostID][pos];
        if (t.v2n[u].size() > 1) {
            newShareAttr.push_back(u);
            prefix.push_back(u);
        }
        else break;
    }
    nodeOrders[nID] = globalOrder(query, t, cs, prefix);
    PrefixNode *attr = pt;
    if (!attributes.empty()) attr = attributes.back();
    for (VertexID u : newShareAttr) {
        PrefixNode *newAttr = new PrefixNode(u);
        attr -> children.push_back(newAttr);
        attr = newAttr;
    }
    if (!prefixSizes.empty()) {
        for (int i = 0; i < attributes.size(); ++i) {
            for (VertexID nID2 : attributes[i]->getBagsBelow())
                prefixSizes[nID2] = i + 1;
        }
        prefixSizes[nID] = prefixSizes[maxCostID] = prefix.size();
    }
    if (!globalOrderShare) {
        newShareAttr.clear();
        prefix.clear();
        nodeOrders[nID] = globalOrder(query, t, cs, prefix);
    }
    if (newShareAttr.empty()) attr->nIDsToCall.push_back(nID);
    else {
        PrefixNode *pn = pt;
        if (!attributes.empty()) pn = attributes.back();
        for (auto it = pn->nIDsToCall.begin(); it != pn->nIDsToCall.end(); ++it) {
            if (*it == maxCostID) {
                pn->nIDsToCall.erase(it);
                break;
            }
        }
        attr->nIDsToCall = {maxCostID, nID};
    }
//    std::vector<PrefixNode *> attributes;
//    std::vector<VertexID> prefix;
//    nodeOrders[nID] = globalOrder(query, t, cs, prefix);
//    int pos = 0;
//    PrefixNode *attr = pt;
//    while (true) {
//        bool flag = false;
//        for (PrefixNode *c : attr->children)
//            if (c->u == nodeOrders[nID][pos]) {
//                flag = true;
//                attr = c;
//                attributes.push_back(attr);
//                ++pos;
//                break;
//            }
//        if (!flag) break;
//    }
//    if (!attributes.empty()) attr = attributes.back();
//    attr->nIDsToCall.push_back(nID);
//    if (!prefixSizes.empty()) {
//        for (int i = 0; i < attributes.size(); ++i) {
//            for (VertexID nID2 : attributes[i]->getBagsBelow())
//                prefixSizes[nID2] = i + 1;
//        }
//    }
}

std::vector<int>
getPrefixAttrID(PrefixNode *pt, const std::vector<VertexID> &prefix, const std::vector<std::vector<int>> &attrIDMap,
                const std::vector<PrefixNode *> &dynamicPartition) {
    std::vector<int> result;
    for (VertexID u : prefix) {
        bool exists = false;
        for (int i = 1; i < dynamicPartition.size(); ++i) {
            if (dynamicPartition[i]->u == u) {
                PrefixNode *current = dynamicPartition[i];
                VertexID firstBag = current->getBagsBelow()[0];
                int depth = pt->locate(current, firstBag);
                exists = true;
                result.push_back(attrIDMap[firstBag][depth]);
                break;
            }
        }
        if (!exists) break;
    }
    ui trivial = 99;
    while (result.size() < prefix.size()) {
        result.push_back(trivial);
        ++trivial;
    }

    return result;
}

int convertToScopePlan(const Graph &query, HyperTree &t, PrefixNode *&pt, CandidateSpace &cs) {
    ui numNodes = t.numNodes;
    if (t.newGlobalNode) --numNodes;
    std::vector<std::vector<VertexID>> prefixOrders(t.numNodes), nodeOrders(t.numNodes);
    for (VertexID nID2 = 0; nID2 < t.numNodes; ++nID2) {
        prefixOrders[nID2].assign(t.nodes[nID2].attributes, t.nodes[nID2].attributes + t.nodes[nID2].prefixSize);
        nodeOrders[nID2].assign(t.nodes[nID2].attributes, t.nodes[nID2].attributes + t.nodes[nID2].numAttributes);
    }
    std::vector<VertexID> defaultPartition;
    std::vector<PrefixNode *> attributeOrder;
    std::vector<VertexID> bagOrder;
    pt->getTraverseOrder(attributeOrder, bagOrder, t);
    for (int i = 0; i < numNodes; ++i) {
        VertexID nID = bagOrder[i];
        ui oldPrefixSize = prefixOrders[nID].size();
        const HyperNode &tau = t.nodes[nID];
        ui newPrefixSize = tau.numAttributes - 2;
        VertexID u_1 = nodeOrders[nID][tau.numAttributes - 2], u_2 = nodeOrders[nID][tau.numAttributes - 1];
        if (query.getEdgeID(u_1, u_2) == -1) ++newPrefixSize;
        if (newPrefixSize < oldPrefixSize) newPrefixSize = oldPrefixSize;
        for (int j = 0; j < oldPrefixSize; ++j) {
            VertexID u = nodeOrders[nID][j];
            if (std::find(defaultPartition.begin(), defaultPartition.end(), u) == defaultPartition.end())
                defaultPartition.push_back(u);
        }
        for (int j = 0; j < newPrefixSize - oldPrefixSize; ++j) {
            VertexID u = nodeOrders[nID][j + oldPrefixSize];
            if (std::find(defaultPartition.begin(), defaultPartition.end(), u) == defaultPartition.end())
                defaultPartition.push_back(u);
            prefixOrders[nID].push_back(u);
        }
        for (int j = i + 1; j < bagOrder.size(); ++j) {
            VertexID nID2 = bagOrder[j];
            for (int k = oldPrefixSize; k < newPrefixSize; ++k) {
                VertexID u = nodeOrders[nID][k];
                if (std::find(nodeOrders[nID2].begin(), nodeOrders[nID2].end(), u) != nodeOrders[nID2].end() &&
                        std::find(prefixOrders[nID2].begin(), prefixOrders[nID2].end(), u) == prefixOrders[nID2].end()) {
                    prefixOrders[nID2].push_back(u);
                    std::vector<VertexID> newNodeOrder = prefixOrders[nID2];
                    for (VertexID u2: nodeOrders[nID2]) {
                        if (std::find(prefixOrders[nID2].begin(), prefixOrders[nID2].end(), u2) == prefixOrders[nID2].end())
                            newNodeOrder.push_back(u2);
                    }
                    nodeOrders[nID2] = newNodeOrder;
                }
            }
        }
    }
    PrefixNode *root = new PrefixNode(99);
    PrefixNode *pn = root;
    std::vector<VertexID> materializedBags;
    for (VertexID nID = 0; nID < t.numNodes; ++nID) {
        if (prefixOrders[nID].size() == 0 && nID != t.numNodes - 1) {
            root->nIDsToCall.push_back(nID);
            root->nIDsToBuild.push_back(nID);
            materializedBags.push_back(nID);
        }
    }
    for (int i = 0; i < defaultPartition.size(); ++i) {
        PrefixNode *attr = new PrefixNode(defaultPartition[i]);
        attr -> pathToGlobal = true;
        pn->children.push_back(attr);
        // find all bags that the prefix is included and adjust the prefix order
        std::vector<VertexID> below;
        for (VertexID nID = 0; nID < t.numNodes; ++nID) {
            if (std::find(materializedBags.begin(), materializedBags.end(), nID) != materializedBags.end()) continue;
            std::vector<std::pair<VertexID, int>> poses;
            for (int j = 0; j < prefixOrders[nID].size(); ++j) {
                int pos = -1;
                VertexID u = prefixOrders[nID][j];
                for (int k = 0; k <= i; ++k) {
                    if (defaultPartition[k] == u)
                        pos = k;
                }
                poses.emplace_back(u, pos);
            }
            bool flag = true;
            for (auto &item : poses) {
                if (item.second == -1) flag = false;
            }
            if (!flag) continue;
            std::sort(poses.begin(), poses.end(), [](const std::pair<VertexID, int>& a, const std::pair<VertexID, int>& b) {
                return a.second < b.second;
            });
            for (int j = 0; j < prefixOrders[nID].size(); ++j) {
                nodeOrders[nID][j] = poses[j].first;
            }
            if (nID != t.numNodes - 1 || i == defaultPartition.size() - 1) {
                below.push_back(nID);
                materializedBags.push_back(nID);
            }
        }
        for (VertexID nID: below) {
            if (nID != t.numNodes - 1) attr->nIDsToBuild.push_back(nID);
        }
        // create one child if they can share one attribute
        std::map<VertexID, std::vector<VertexID>> firstAttrToBags;
        for (VertexID nID : below) {
            if (nID == t.numNodes - 1) {
                attr->nIDsToCall.push_back(nID);
                continue;
            }
            if (prefixOrders[nID].size() == nodeOrders[nID].size() - 2) {
                VertexID u = nodeOrders[nID][prefixOrders[nID].size()];
                firstAttrToBags[u].push_back(nID);
            }
            else attr->nIDsToCall.push_back(nID);
        }
        for (auto &item: firstAttrToBags) {
            if (item.second.size() > 1) {
                VertexID u = item.first;
                PrefixNode *child = new PrefixNode(u);
                child->pathToGlobal = false;
                child->nIDsToCall = item.second;
                std::vector<PrefixNode *> oldChild = attr->children;
                attr->children = {child};
                for (PrefixNode *c: oldChild) attr->children.push_back(c);
            }
            else attr->nIDsToCall.push_back(item.second[0]);
        }
        for (int j = 0; j < attr->nIDsToCall.size(); ++j) {
            if (j == attr->nIDsToCall.size() - 1) continue;
            if (attr->nIDsToCall[j] == t.numNodes - 1) std::swap(attr->nIDsToCall[j], attr->nIDsToCall.back());
        }
        pn = attr;
    }
    root->initPoses(nodeOrders, query, cs.dist);
    delete pt;
    pt = root;
    HyperTree *newTree = new HyperTree();
    newTree->newGlobalNode = t.newGlobalNode;
    newTree->numNodes = t.numNodes;
    newTree->nodes = new HyperNode[newTree->numNodes];
    newTree->numAttributes = t.numAttributes;
    for (VertexID nID2 = 0; nID2 < t.numNodes; ++nID2) {
        t.nodes[nID2].copyTo(newTree->nodes[nID2]);
        newTree->nodes[nID2].prefixSize = prefixOrders[nID2].size();
        if (prefixOrders[nID2].size() != 0 && t.nodes[nID2].prefixSize == 0) newTree->nodes[nID2].prefix = new VertexID[t.numAttributes];
        for (int i = 0; i < prefixOrders[nID2].size(); ++i)
            newTree->nodes[nID2].prefix[i] = nodeOrders[nID2][i];
        for (int i = 0; i < nodeOrders[nID2].size(); ++i)
            newTree->nodes[nID2].attributes[i] = nodeOrders[nID2][i];
    }
    int returnValue = defaultPartition.size();
    newTree->defaultPartition = defaultPartition;
    for (VertexID u: nodeOrders.back()) {
        if (std::find(newTree->defaultPartition.begin(), newTree->defaultPartition.end(), u) == newTree->defaultPartition.end())
            newTree->defaultPartition.push_back(u);
    }
    const HyperNode &last = newTree->nodes[newTree->numNodes - 1];
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
    for (int i = prefixOrders[t.numNodes - 1].size(); i < nodeOrders[t.numNodes - 1].size(); ++i) {
        VertexID u = nodeOrders[t.numNodes - 1][i];
        if (std::find(newTree->globalOrder.begin(), newTree->globalOrder.end(), u) == newTree->globalOrder.end())
            newTree->globalOrder.push_back(u);
    }
    for (VertexID nID2 = 0; nID2 < numNodes; ++nID2) {
        for (int i = prefixOrders[nID2].size(); i < newTree->nodes[nID2].numAttributes; ++i) {
            VertexID u = newTree->nodes[nID2].attributes[i];
            if (std::find(newTree->globalOrder.begin(), newTree->globalOrder.end(), u) == newTree->globalOrder.end())
                newTree->globalOrder.push_back(u);
        }
    }
    newTree->initPoses(query, cs, false);
    newTree->extendLevel = t.nodes[t.numNodes - 1].numAttributes;
    newTree->buildTraverseStruct(query);
    t = *newTree;

    return returnValue;
}