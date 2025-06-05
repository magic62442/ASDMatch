//
// Created by Qiyan LI on 2024/3/2.
//

#ifndef IN_MEMORY_JOIN_UTILS_H
#define IN_MEMORY_JOIN_UTILS_H

#include "config.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <chrono>
#include <cstring>
#include <random>
#include <unistd.h>
#include "td/hypergraph.h"
#include "td/preprocessing.h"
#include "td/td.h"
#include "td/td_utils.hpp"

static std::random_device rd;
static std::mt19937 gen(rd());

bool nlcValid(const std::map<LabelID, ui> &queryNLC, const std::map<LabelID, ui> &dataNLC);
void sampleKElements(VertexID *array, ui length, ui k);
void leftShit(std::vector<bool> &arr, int k);
std::vector<std::vector<bool>> chooseK(ui n, int k);
int getPosition(uint64_t id, ui k, const std::vector<uint64_t> &allSets);
void subsetSupersetRelationships(ui n, const std::vector<VertexID> &elements,
                                 const std::vector<std::vector<uint64_t>> &subsets, std::vector<std::vector<std::vector<int>>> &subsetOf,
                                 std::vector<std::vector<std::vector<int>>> &supersetOf);
size_t numUniqueCombine(std::vector<std::vector<VertexID>> &vsets);
size_t bfUniqueCombine(std::vector<std::vector<VertexID>> &vsets);
int findFirstExtension(const std::vector<std::vector<uint32_t>>& vec, const std::vector<uint32_t>& prefix);
int findFirstGreater(const std::vector<std::vector<uint32_t>>& vec, const std::vector<uint32_t>& prefix);
void makeContinuous(std::vector<VertexID>& nums);
bool checkBudget(size_t oldBudget, size_t budget, const std::vector<size_t> &tupleSizes, const std::vector<std::vector<std::vector<VertexID>>> &tuples,
                 const std::vector<std::vector<std::vector<VertexID>>> &backup);

double getMemoryUsageGB();

template <typename T>
void printArray(T ** arrays, ui *counts, ui num) {
    for (int i = 0; i < num; ++i) {
        for (int j = 0; j < counts[i]; ++j) {
            std::cout << arrays[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

template <typename T>
void writeArrayToStream(std::ofstream& outFile, const T* array, size_t size) {
    if (array != nullptr && size > 0) {
        outFile.write(reinterpret_cast<const char*>(array), sizeof(T) * size);
    }
}

template <typename T>
void readArrayFromStream(std::ifstream& inFile, T*& array, size_t size) {
    if (size > 0) {
        array = new T[size];
        inFile.read(reinterpret_cast<char*>(array), sizeof(T) * size);
    } else {
        array = nullptr;
    }
}

template <typename T>
void write2DArrayToStream(std::ofstream& outFile, T** array2D, size_t dim1Size, size_t dim2Size) {
    if (array2D != nullptr && dim1Size > 0) {
        for (size_t i = 0; i < dim1Size; ++i) {
            writeArrayToStream(outFile, array2D[i], dim2Size);
        }
    }
}

template <typename T>
void read2DArrayFromStream(std::ifstream& inFile, T**& array2D, size_t dim1Size, size_t dim2Size) {
    array2D = new T*[dim1Size];
    for (size_t i = 0; i < dim1Size; ++i) {
        readArrayFromStream(inFile, array2D[i], dim2Size);
    }
}

template <typename T1, typename T2>
void writeMapToStream(std::ofstream& outFile, const std::map<T1, T2>& map) {
    size_t size = map.size();
    outFile.write(reinterpret_cast<const char*>(&size), sizeof(size));
    for (const auto& item : map) {
        outFile.write(reinterpret_cast<const char*>(&(item.first)), sizeof(T1));
        outFile.write(reinterpret_cast<const char*>(&(item.second)), sizeof(T2));
    }
}

template <typename T1, typename T2>
void readMapFromStream(std::ifstream& inFile, std::map<T1, T2>& map) {
    size_t size;
    inFile.read(reinterpret_cast<char*>(&size), sizeof(size));
    map.clear();
    T1 key;
    T2 value;
    for (size_t i = 0; i < size; ++i) {
        inFile.read(reinterpret_cast<char*>(&key), sizeof(T1));
        inFile.read(reinterpret_cast<char*>(&value), sizeof(T2));
        map[key] = value;
    }
}

struct BipartiteMaximumMatching {
    int *left, *right;
    int left_len, right_len;
    bool *used;
    int **adj, *adj_size;
    int **adj_index;
    bool **matchable;


    bool *bfs_visited;
    int *right_order, *inverse_right_order;
    int **lower_graph, *lower_graph_size;
    int **upper_graph, *upper_graph_size;

    int *Q, *S;
    int qright = 0, qleft = 0;
    int stkright = 0;
    int dfsn[MAX_PATTERN_SIZE], scch[MAX_PATTERN_SIZE], scc_idx[MAX_PATTERN_SIZE], ord, found_scc;
    int scc_dfs(int v) {
        S[stkright++] = v;
        int ret=dfsn[v]=++ord;
        for(int i = 0; i < upper_graph_size[v]; i++){
            int n = upper_graph[v][i];
            if(!dfsn[n])	ret=std::min(ret, scc_dfs(n));
            else if(!scch[n])	ret=std::min(ret, dfsn[n]);
        }
        if(ret==dfsn[v]){
            int u;
            do {
                u=S[--stkright];
                scch[u]=1;
                scc_idx[u] = found_scc;
            } while(u!=v);
            found_scc++;
        }
        return ret;
    }
    bool FindUnmatchableEdges(int required) {
        ord = 0;
        found_scc = 0;
        memset(dfsn, 0, sizeof (dfsn));
        memset(scch, 0, sizeof (scch));
        memset(scc_idx, -1, sizeof (scc_idx));
        int num_matched_ans = solve();
        if (num_matched_ans != required) return false;
        for (int i = 0; i < num_matched_ans; i++) {
            for (int j = 0; j < adj_size[i]; j++) {
                matchable[i][adj[i][j]] = false;
            }
        }
        for (int i = 0; i < num_matched_ans; i++) {
            matchable[i][left[i]] = true;
        }
        int num_rights = 0;
        for (int i = 0; i < num_matched_ans; i++) {
            inverse_right_order[num_rights] = left[i];
            right_order[left[i]] = num_rights++;
        }
        for (int i = 0; i < num_matched_ans; i++) {
            for (int j = 0; j < adj_size[i]; j++) {
                int r = right_order[adj[i][j]];
                if (r == -1) {
                    inverse_right_order[num_rights] = adj[i][j];
                    right_order[adj[i][j]] = num_rights++;
                    r = right_order[adj[i][j]];
                }
                if (r != i) {
                    lower_graph[r][lower_graph_size[r]++] = i;
                }
                if (r != i and r < num_matched_ans) {
                    upper_graph[i][upper_graph_size[i]++] = r;
                }
            }
        }
        for(int i=0;i<num_matched_ans;i++) {
            if (!dfsn[i])
                scc_dfs(i);
        }
        for (int i = 0; i < num_matched_ans; i++) {
            for (int j = 0; j < upper_graph_size[i]; j++) {
                int r = upper_graph[i][j];
                if (scc_idx[i] == scc_idx[r]) {
                    matchable[i][inverse_right_order[r]] = true;
                }
            }
        }

        std::memset(bfs_visited, 0, sizeof(bool) * right_len);
        for (int i = num_matched_ans; i < num_rights; i++) {
            Q[qright++] = i;
            bfs_visited[i] = true;
        }
        while (qright != qleft) {
            int u = Q[qleft++];
            for (int j = 0; j < lower_graph_size[u]; j++) {
                int v = lower_graph[u][j];
                matchable[v][inverse_right_order[u]] = true;
                if (!bfs_visited[v]) {
                    Q[qright++] = v;
                    bfs_visited[v] = true;
                }
            }
        }
        return true;
    }

    void global_initialize(int max_left, int max_right) {
        Q = new int[max_right];
        S = new int[max_right];
        qleft = qright = stkright = 0;
        left = new int[max_left];
        right = new int[max_right];
        left_len = max_left;
        right_len = max_right;
        used = new bool[max_left];
        adj = new int *[max_left];
        adj_index = new int *[max_left];
        matchable = new bool *[max_left];

        lower_graph = new int *[max_right];
        lower_graph_size = new int[max_right];

        upper_graph = new int *[max_left];
        upper_graph_size = new int[max_left];
        right_order = new int[max_right];
        inverse_right_order = new int[max_right];
        for (int i = 0; i < max_left; i++) {
            adj[i] = new int[max_right];
            adj_index[i] = new int[max_right];
            matchable[i] = new bool[max_right];
            upper_graph[i] = new int[max_left];
        }
        for (int i = 0; i < max_right; i++) {
            lower_graph[i] = new int[max_left];
        }
        adj_size = new int[max_left];
        bfs_visited = new bool[max_right];
    }

    ~BipartiteMaximumMatching() {
        delete[] left;
        delete[] right;
        delete[] used;
        for (int i = 0; i < left_len; i++) {
            delete[] adj[i];
            delete[] matchable[i];
            delete[] upper_graph[i];
        }
        delete[] matchable;
        delete[] adj;
        delete[] adj_size;
        delete[] lower_graph;
        delete[] upper_graph;
        delete[] lower_graph_size;
        delete[] upper_graph_size;
        delete[] right_order;
        delete[] inverse_right_order;
        delete[] bfs_visited;
    }

    void reset(bool reset_edges = true) {
        std::memset(left, -1, sizeof(int) * left_len);
        std::memset(right, -1, sizeof(int) * right_len);
        std::memset(used, false, sizeof(bool) * left_len);
        if (reset_edges) {
            std::memset(adj_size, 0, sizeof(int) * left_len);
        }
        std::memset(lower_graph_size, 0, sizeof(int) * right_len);
        std::memset(upper_graph_size, 0, sizeof(int) * left_len);
        std::memset(right_order, -1, sizeof(int) * right_len);
        std::memset(inverse_right_order, -1, sizeof(int) * right_len);
        qleft = qright = stkright = 0;
    }

    void add_edge(int u, int v) {
        adj_index[u][v] = adj_size[u];
        adj[u][adj_size[u]++] = v;
    }

    int remove_edge(int u, int v) {
        if (adj_size[u] > 1) {
            adj_index[u][adj[u][adj_size[u] - 1]] = adj_index[u][v];
            std::swap(adj[u][adj_size[u] - 1], adj[u][adj_index[u][v]]);
        }
        return --adj_size[u];
    }

    void revert(int *tmp_left) {
        for (int i = 0; i < left_len; i++) {
            if (left[i] == -1) continue;
            right[left[i]] = -1;
        }
        std::memcpy(left, tmp_left, sizeof(int) * left_len);
        for (int i = 0; i < left_len; i++) {
            if (left[i] == -1) continue;
            right[left[i]] = i;
        }
    }

    int solve(int ignore = -1) {
        std::memset(left, -1, sizeof(int) * left_len);
        int ans = 0;
        for (VertexID u = 0; u < left_len; u++) {
            if (u == ignore) continue;
            if (left[u] == -1) {
                std::memset(used, false, sizeof(bool) * left_len);
                if (dfs(u)) {
                    ans++;
                }
            }
        }
        return ans;
    }

    bool dfs(VertexID r) {
        if (used[r]) return false;
        used[r] = true;
        for (int i = 0; i < adj_size[r]; i++) {
            VertexID c = adj[r][i];
            VertexID k = right[c];
            if (k == -1 or dfs(k)) {
                left[r] = c;
                right[c] = r;
                return true;
            }
        }
        return false;
    }

    bool single_dfs(VertexID r) {
        if (used[r]) return false;
        used[r] = true;
        for (int i = 0; i < adj_size[r]; i++) {
            VertexID c = adj[r][i];
            VertexID k = right[c];
            if (k == -1 or single_dfs(k)) {
                return true;
            }
        }
        return false;
    }
};

template<typename T>
size_t calculate4DVectorSize(const std::vector<std::vector<std::vector<std::vector<T>>>>& vec) {
    size_t totalCapacity = 0;
    totalCapacity += sizeof(std::vector<T>) * vec.capacity();
    for (const auto& secondLevelVec : vec) {
        totalCapacity += sizeof(std::vector<T>) * secondLevelVec.capacity();
        for (const auto& thirdLevelVec : secondLevelVec) {
            totalCapacity += sizeof(std::vector<T>) * thirdLevelVec.capacity();
            for (const auto& fourthLevelVec : thirdLevelVec) {
                totalCapacity += sizeof(T) * fourthLevelVec.capacity(); // Total capacity of the vector
            }
        }
    }

    return totalCapacity;
}

class UnionFind {
public:
    // Constructor initializes the disjoint sets
    UnionFind(int n) : parent(n), rank(n, 0) {
        for (int i = 0; i < n; ++i) {
            parent[i] = i;  // Each element is its own parent (self root)
        }
    }

    // Find the root of the element x with path compression
    int find(int x) {
        if (parent[x] != x) {
            parent[x] = find(parent[x]);  // Path compression
        }
        return parent[x];
    }

    // Union by rank
    void unite(int x, int y) {
        int rootX = find(x);
        int rootY = find(y);
        if (rootX != rootY) {
            if (rank[rootX] > rank[rootY]) {
                parent[rootY] = rootX;
            } else if (rank[rootX] < rank[rootY]) {
                parent[rootX] = rootY;
            } else {
                parent[rootY] = rootX;
                rank[rootX]++;
            }
        }
    }

    // Check if two elements are in the same set
    bool connected(int x, int y) {
        return find(x) == find(y);
    }

    // Get the disjoint sets as a 2D vector
    std::vector<std::vector<int>> getDisjointSets() {
        std::map<int, std::vector<int>> setsMap;
        for (int i = 0; i < parent.size(); ++i) {
            int root = find(i);
            setsMap[root].push_back(i);
        }
        std::vector<std::vector<int>> disjointSets;
        for (auto& set : setsMap) {
            disjointSets.push_back(set.second);
        }
        return disjointSets;
    }

private:
    std::vector<int> parent;  // Parent of each element
    std::vector<int> rank;    // Rank (approximate depth) of each tree
};

#endif //IN_MEMORY_JOIN_UTILS_H
