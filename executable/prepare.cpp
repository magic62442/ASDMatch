//
// Created by Qiyan LI on 2024/3/4.
//

#include <chrono>
#include "candidate_space.h"
#include "join.h"
#include "command.h"
#include "optimizer.h"
#include "adaptive.h"

void storeCard(int argc, char **argv) {
    Command cmd(argc, argv);
    std::string queryGraphPath = cmd.getQueryGraphPath();
    std::string dataGraphPath = cmd.getDataGraphPath();
    std::string resultPath = cmd.getResultPath();
    std::string dpStructPath = cmd.getDPStructPath();
    std::ifstream inFile(dpStructPath);
//    if (inFile.good()) return;
    int memKB = cmd.getMemoryBudget();
    if (memKB == 0) memKB = 16e5;
    int baselineType = cmd.getBaselineType();
    size_t budget = (size_t)1024 * memKB;
    Graph q, g;
    q.loadGraphFromTextFile(queryGraphPath);
    g.loadGraphFromTextFile(dataGraphPath);
    CandidateSpace cs(q, g);
    auto start = std::chrono::steady_clock::now();
    if (!cs.buildCandCFL()) {
        return;
    }
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsedSeconds = end - start;
    cs.setQueryGraphWeights(q);
    FHD fhd;
    q.buildFHD(fhd);
    HyperTree t;
    t.buildFromTD(fhd);
    smallCard(q, cs);
    std::vector<TrieNode *> nodes(t.numNodes);
    for (int i = 0; i < t.numNodes; ++i) nodes[i] = new TrieNode();
    std::vector<std::vector<VertexID>> result1, result2, difference;
    bool *visited = new bool[g.getNumVertices()];
    memset(visited, false, sizeof(bool) * g.getNumVertices());
    VertexID *partMatch = new VertexID[q.getNumVertices()];
    VertexID **candidates = new VertexID * [q.getNumVertices()];
    ui maxSize = cs.getMaxSize();
    for (int i = 0; i < q.getNumVertices(); ++i) {
        candidates[i] = new VertexID[maxSize];
    }
    ui *candCount = new ui[q.getNumVertices()];
    memset(candCount, 0, sizeof(ui) * q.getNumVertices());
    reorderBags(q, t);
    ui numNodes = t.numNodes;
    if (t.newGlobalNode) --numNodes;
    start = std::chrono::steady_clock::now();
    std::vector<SubsetStructure> dpStructures(numNodes);
    for (VertexID nID = 0; nID < numNodes; ++nID) {
        std::vector<VertexID> localVertices(t.nodes[nID].attributes, t.nodes[nID].attributes + t.nodes[nID].numAttributes);
        SubsetStructure s(localVertices, q, false);
        s.optimalPlanDP(q, cs, visited, partMatch, candidates, candCount);
        s.reverseDP(q, cs);
        dpStructures[nID] = s;
    }
    std::ofstream ofs(dpStructPath);
    for (const auto &s: dpStructures) s.saveToSteam(ofs);
    saveSubsetToCard(ofs);
    end = std::chrono::steady_clock::now();
    elapsedSeconds = end - start;
    int planTime = static_cast<int>(elapsedSeconds.count() * 10000);
    ofs.write(reinterpret_cast<const char*>(&planTime), sizeof(planTime));
    ofs << std::flush;
}

int main(int argc, char** argv) {
    storeCard(argc, argv);
    return 0;
}