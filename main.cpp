#include <chrono>
#include "candidate_space.h"
#include "join.h"
#include "command.h"
#include "optimizer.h"
#include "adaptive.h"
#include "estimator.h"
#include <cassert>
#include <iomanip>

void clearCache(int &simulation) {
    std::vector<int> data(10000000);
    std::generate(data.begin(), data.end(), [](){ return rand() % 100; });
    simulation = 0;
    for (auto &d : data) simulation += d;
}

void ASDMatch(int argc, char **argv) {
    bool skip = false;
    Command cmd(argc, argv);
    std::string queryGraphPath = cmd.getQueryGraphPath();
    std::string dataGraphPath = cmd.getDataGraphPath();
    std::string resultPath = cmd.getResultPath();
    std::string dpStructPath = cmd.getDPStructPath();
    int memKB = cmd.getMemoryBudget();
    if (memKB == 0) memKB = 16e5;
    int baselineType = cmd.getBaselineType();
    size_t budget = (size_t)1024 * memKB;
    Graph q, g;
    q.loadGraphFromTextFile(queryGraphPath);
    g.loadGraphFromTextFile(dataGraphPath);
    CandidateSpace cs(q, g);
    std::ofstream outFile;
    std::ostream &outStream = resultPath.empty() ? std::cout : outFile ;
    if (!resultPath.empty()) outFile.open(resultPath);
    auto start = std::chrono::steady_clock::now();
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsedSeconds = end - start;
    start = std::chrono::steady_clock::now();
    if (!cs.buildCandCFL()) {
        return;
    }
    end = std::chrono::steady_clock::now();
    elapsedSeconds = end - start;
    outStream << "Filtering time (s): " << elapsedSeconds.count() << std::endl;
    cs.setQueryGraphWeights(q);
    start = std::chrono::steady_clock::now();
    FHD fhd;
    q.buildFHD(fhd);
    HyperTree t;
    t.buildFromTD(fhd);
    end = std::chrono::steady_clock::now();
    elapsedSeconds = end - start;
    double fhdTime = elapsedSeconds.count();
    outStream << "FHD time (s): " << fhdTime << std::endl;
    std::vector<TrieNode *> nodes(t.numNodes);
    for (int i = 0; i < t.numNodes; ++i) nodes[i] = new TrieNode();
    std::vector<std::vector<VertexID>> result;
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
    PrefixNode *pt;
    std::vector<double> costs;
    double minCost;
    size_t count = 0;
    double cardTime = 0.0;
    reorderBags(q, t);
    std::vector<SubsetStructure> dpStructures;
    ui numNodes = t.numNodes;
    if (t.newGlobalNode) --numNodes;
    if (dpStructPath.empty()) {
        smallCard(q, cs);
        dpStructures.resize(numNodes);
        start = std::chrono::steady_clock::now();
        for (VertexID nID = 0; nID < numNodes; ++nID) {
            std::vector<VertexID> localVertices(t.nodes[nID].attributes, t.nodes[nID].attributes + t.nodes[nID].numAttributes);
            SubsetStructure s(localVertices, q, false);
            s.optimalPlanDP(q, cs, visited, partMatch, candidates, candCount);
            s.reverseDP(q, cs);
            dpStructures[nID] = s;
        }
        end = std::chrono::steady_clock::now();
        elapsedSeconds = end - start;
        cardTime = elapsedSeconds.count();
    }
    else {
        std::ifstream ifs(dpStructPath, std::ios::binary);
        dpStructures.resize(numNodes);
        for (VertexID nID = 0; nID < numNodes; ++nID) {
            dpStructures[nID].readFromStream(ifs);
        }
        if (!loadSubsetToCard(ifs)) return;
        int scaled;
        ifs.read(reinterpret_cast<char*>(&scaled), sizeof(int));
        cardTime = static_cast<double>(scaled) / 10000.0;
        outStream << "dp time: " << cardTime << std::endl;
    }
    cardTime += fhdTime;
    int simulation;
    clearCache(simulation);
    bool traverse = false;
    int reorder = 1;
    start = std::chrono::steady_clock::now();
    if (numNodes != 1) bagPermutationDP(q, t, cs, pt, dpStructures, minCost, reorder, true, true);
    else simplePlan(q, t, cs, pt, visited, partMatch, candidates, candCount, minCost, nullptr, false);
    end = std::chrono::steady_clock::now();
    elapsedSeconds = end - start;
    outStream << "Planning Time: " << elapsedSeconds.count() + cardTime << std::endl;
    outStream << "cost: " << minCost << std::endl;
    start = std::chrono::steady_clock::now();
    adaptiveShareJoin(q, t, pt, cs, nodes, visited, result, count, traverse, budget, skip);
    end = std::chrono::steady_clock::now();
    elapsedSeconds = end - start;
    for (int i = 0; i < t.numNodes; ++i) {
        if (t.newGlobalNode && i == t.numNodes - 1) continue;
        delete nodes[i];
        nodes[i] = new TrieNode();
    }
    outStream << "Execution Time: " << elapsedSeconds.count() << std::endl;
    outStream << "Number of matches: " << count << std::endl;
    outStream << "Number of fixed attributes: " << gMaxNewFixed << std::endl;
    outStream << "Number of intersections: " << gNumInterSection << std::endl;
    outStream << "Number of fix cases: " << gNumFixedCase << std::endl;
#ifdef COLLET_GLOBAL_TIME
    outStream << "global join time: " << gGlobalTime << std::endl;
#endif
    t.writeToStream(outStream);
}

int main(int argc, char **argv) {
    ASDMatch(argc, argv);
    return 0;
}
