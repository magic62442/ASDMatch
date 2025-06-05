//
// Created by Qiyan LI on 2024/3/2.
//

#include "utils.h"

double getMemoryUsageGB() {
    std::ifstream statm("/proc/self/statm");
    long pageSize = sysconf(_SC_PAGESIZE); // Get system page size in bytes

    std::string line;
    double totalProgramSizeGB = 0.0;
    double residentSetSizeGB = 0.0;

    if (std::getline(statm, line)) {
        std::istringstream iss(line);
        long totalProgramPages, residentSetPages;
        if (iss >> totalProgramPages >> residentSetPages) {
            totalProgramSizeGB = totalProgramPages * pageSize / 1024.0 / 1024.0 / 1024.0;
            residentSetSizeGB = residentSetPages * pageSize / 1024.0 / 1024.0 / 1024.0;
        }
    }
    statm.close();

//    std::cout << "Total Program Size: " << totalProgramSizeGB << " GB" << std::endl;
//    std::cout << "Resident Set Size: " << residentSetSizeGB << " GB" << std::endl;
    return residentSetSizeGB;
}

bool nlcValid(const std::map<LabelID, ui> &queryNLC, const std::map<LabelID, ui> &dataNLC) {
    if (queryNLC.size() > dataNLC.size()) return false;
    for (auto it : queryNLC) {
        VertexID uL = it.first;
        ui uCount = it.second;
        if (dataNLC.find(uL) == dataNLC.end() || dataNLC.at(uL) < uCount)
            return false;
    }

    return true;
}

void sampleKElements(VertexID *array, ui length, ui k) {
    std::uniform_int_distribution<ui> dis(0, length - 1);
    std::set<ui> sampledIndices;
    while (sampledIndices.size() < k) {
        ui r = dis(gen);
        sampledIndices.insert(r);
    }
    ui i = 0;
    for (ui index: sampledIndices) {
        array[i] = array[index];
        ++i;
    }
}

void leftShit(std::vector<bool> &arr, int k) {
    int count = 0;
    for (int j = k; j >= 0; --j)
        if ( arr[j] == 1)  ++count;
    if (count == 0)  return;
    for (int j = k; j >= 0; --j)
        arr[j] = false;
    for ( int j=0; j<count;++j)
        arr[j] = true;
}

std::vector<std::vector<bool>> chooseK(ui n, int k) {
    if (k > n) {
        std::vector<std::vector<bool>> result;
        result.emplace_back(n, true);
        return result;
    }
    std::vector<bool> arr(n, false);
    for (int i = 0; i < k; ++i)
        arr[i] = true;
    std::vector<std::vector<bool>> result;
    result.push_back(arr);
    bool flag = true;
    while (flag) {
        flag = false;
        for (int i = 0; i < n - 1; ++i) {
            if (arr[i] && !arr[i + 1]){
                flag = true;
                arr[i] = false;
                arr[i + 1] = true;
                leftShit(arr,i);
                result.push_back(arr);
                break;
            }
        }
    }

    return result;
}

int getPosition(uint64_t id, ui k, const std::vector<uint64_t> &allSets) {
    auto it = std::lower_bound(allSets.begin(), allSets.end(), id);
    if (it != allSets.end() && *it == id) {
        return std::distance(allSets.begin(), it);
    }
    return -1; // Not found
}

void subsetSupersetRelationships(ui n, const std::vector<VertexID> &elements,
                                 const std::vector<std::vector<uint64_t>> &subsets, std::vector<std::vector<std::vector<int>>> &subsetOf,
                                 std::vector<std::vector<std::vector<int>>> &supersetOf) {
    subsetOf.resize(n + 1);
    supersetOf.resize(n + 1);
    for (int k = 1; k <= n; ++k) {
        subsetOf[k].resize(subsets[k].size());
        if (k > 1) {
            for (size_t i = 0; i < subsets[k].size(); ++i) {
                uint64_t subsetID = subsets[k][i];
                for (const auto& element : elements) {
                    uint64_t elementBit = 1ULL << (element % 64);
                    if (subsetID & elementBit) {
                        uint64_t smallerSubsetID = subsetID & ~elementBit;
                        int pos = getPosition(smallerSubsetID, k - 1, subsets[k - 1]);
                        if (pos != -1) {
                            subsetOf[k][i].push_back(pos);
                        }
                    }
                }
            }
        }
        supersetOf[k].resize(subsets[k].size());
        if (k < n) {
            for (size_t i = 0; i < subsets[k].size(); ++i) {
                uint64_t subsetID = subsets[k][i];
                for (const auto& element : elements) {
                    uint64_t elementBit = 1ULL << (element % 64);
                    if (!(subsetID & elementBit)) {
                        uint64_t superSetID = subsetID | elementBit;
                        int pos = getPosition(superSetID, k + 1, subsets[k + 1]);
                        if (pos != -1) {
                            supersetOf[k][i].push_back(pos);
                        }
                    }
                }
            }
        }
    }
}

size_t numUniqueCombine(std::vector<std::vector<VertexID>> &vsets) {
    size_t num = 1;
    for (auto &vset : vsets) num *= vset.size();
    std::map<uint64_t, std::vector<VertexID>> id2Intersection;
    std::vector<std::pair<int, int>> pairs;
    // pairwise intersection
    for (int i = 0; i < vsets.size(); ++i) {
        uint64_t id = 1 << i;
        id2Intersection[id] = vsets[i];
        for (int j = i + 1; j < vsets.size(); ++j) {
            id  = (1 << i) + (1 << j);
            std::vector<VertexID> intersect;
            std::set_intersection(vsets[i].begin(), vsets[i].end(), vsets[j].begin(),
                                  vsets[j].end(), std::back_inserter(intersect));
            id2Intersection[id] = intersect;
            pairs.emplace_back(i, j);
        }
    }
    int n = vsets.size();
    for (int k = 1; k <= pairs.size(); ++k) {
        // choose k sets in n(n-1)/2 pairwises
        std::vector<std::vector<bool>> choices = chooseK(pairs.size(), k);
        for (auto &choice : choices) {
            UnionFind uf(n);
            for (int i = 0; i < pairs.size(); ++i) {
                if (choice[i])
                    uf.unite(pairs[i].first, pairs[i].second);
            }
            std::vector<std::vector<int>> disjointSets = uf.getDisjointSets();
            size_t term = 1;
            for (auto &disjointSet : disjointSets) {
                uint64_t id = 0;
                for (int pos : disjointSet) id |= (1 << pos);
                if (id2Intersection.find(id) != id2Intersection.end())
                    term *= id2Intersection[id].size();
                else {
                    for (int pos : disjointSet) {
                        uint64_t subsetID = id - (1 << pos);
                        if (id2Intersection.find(subsetID) != id2Intersection.end()) {
                            std::vector<VertexID> oldIntersect = id2Intersection[subsetID];
                            std::vector<VertexID> intersect;
                            std::set_intersection(oldIntersect.begin(), oldIntersect.end(), vsets[pos].begin(),
                                                  vsets[pos].end(), std::back_inserter(intersect));
                            id2Intersection[id] = intersect;
                            term *= intersect.size();
                            break;
                        }
                    }
                }
            }
            if (k % 2 == 0) num += term;
            else num -= term;
        }
    }

    return num;
}

size_t bfUniqueCombine(std::vector<std::vector<VertexID>> &vsets) {
    size_t num = 0;
    ui n = vsets.size();
    std::vector<int> poses(n, 0);
    int depth = 0;
    std::vector<VertexID> tuple(n);
    while (depth >= 0) {
        while (poses[depth] < vsets[depth].size()) {
            VertexID v = vsets[depth][poses[depth]];
            tuple[depth] = v;
            ++poses[depth];
            bool exists = false;
            for (int i = 0; i < depth; ++i) {
                if (tuple[i] == v) {
                    exists = true;
                    break;
                }
            }
            if (exists) continue;
            if (depth == n - 1)
                ++num;
            else {
                ++depth;
                poses[depth] = 0;
            }
        }
        --depth;
    }

    return num;
}

int findFirstExtension(const std::vector<std::vector<uint32_t>>& vec, const std::vector<uint32_t>& prefix) {
    int left = 0;
    int right = vec.size();

    while (left < right) {
        int mid = left + (right - left) / 2;
        bool vecLarger = true;
        for (int i = 0; i < prefix.size(); ++i) {
            if (vec[mid][i] < prefix[i]) {
                vecLarger = false;
                break;
            }
            else if(vec[mid][i] > prefix[i]) break;
        }
        if (vecLarger) {
            right = mid;
        } else {
            left = mid + 1;
        }
    }

    return left;
}

int findFirstGreater(const std::vector<std::vector<uint32_t>>& vec, const std::vector<uint32_t>& prefix) {
    int left = 0;
    int right = vec.size();

    while (left < right) {
        int mid = left + (right - left) / 2;
        bool vecSmaller = true;
        for (int i = 0; i < prefix.size(); ++i) {
            if (vec[mid][i] > prefix[i]) {
                vecSmaller = false;
                break;
            }
            else if (vec[mid][i] < prefix[i]) break;
        }
        if (vecSmaller) {
            left = mid + 1;
        } else {
            right = mid;
        }
    }

    return left;
}

void makeContinuous(std::vector<VertexID>& nums) {
    // Copy and sort the original vector
    std::vector<VertexID> sorted_nums = nums;
    std::sort(sorted_nums.begin(), sorted_nums.end());

    // Map each unique value to its new continuous value
    std::map<int, int> value_map;
    int new_value = 0;
    for (int num : sorted_nums) {
        if (value_map.find(num) == value_map.end()) {
            value_map[num] = new_value++;
        }
    }

    // Replace each value in the original vector with its new value
    for (auto& num : nums) {
        num = value_map[num];
    }
}

bool checkBudget(size_t oldBudget, size_t budget, const std::vector<size_t> &tupleSizes, const std::vector<std::vector<std::vector<VertexID>>> &tuples,
                 const std::vector<std::vector<std::vector<VertexID>>> &backup) {
    for (size_t sz: tupleSizes) budget += sz;
    for (int i = 0; i < tuples.size(); ++i) {
        if (!tuples[i].empty()) budget += tuples[i].size() * tuples[i][0].size() * sizeof(VertexID);
        if (!backup[i].empty()) budget += backup[i].size() * backup[i][0].size() * sizeof(VertexID);
    }

    return oldBudget == budget;
}