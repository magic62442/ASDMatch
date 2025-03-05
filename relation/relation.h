//
// Created by anonymous authors on 2024/2/27.
//

#ifndef IN_MEMORY_JOIN_RELATION_H
#define IN_MEMORY_JOIN_RELATION_H

#include "graph.h"
#include "dynamic_array.h"

struct Compression {
    VertexID **data;
    VertexID *length;
    ui num;

    Compression(): data(nullptr), length(nullptr), num(0) {};
    Compression(VertexID **data, VertexID *length, ui num);
    Compression(const Compression &other) = delete;
    Compression& operator=(const Compression &other) = delete;
    ~Compression() {
        for (int i = 0; i < num; ++i)
            delete[] data[i];
        delete[] data;
        delete[] length;
    }

    size_t memoryCost() const;

    friend bool operator==(const Compression& lhs, const Compression& rhs) {
        if (lhs.num != rhs.num) return false;

        for (ui i = 0; i < lhs.num; ++i) {
            if (lhs.length[i] != rhs.length[i]) return false;

            for (VertexID j = 0; j < lhs.length[i]; ++j) {
                if (lhs.data[i][j] != rhs.data[i][j]) return false;
            }
        }

        return true;
    }
};

struct TrieNode {
    VertexID value;
    DynamicArray<TrieNode *> nodeChild;
//    Compression *compression;

//    TrieNode(): value(0), nodeChild(), compression(nullptr) {}
    TrieNode(): value(0), nodeChild() {}
    TrieNode(const TrieNode& other) = delete;
    TrieNode& operator=(const TrieNode& other) = delete;
    ~TrieNode() {
        for (int i = 0; i < nodeChild.size(); ++i) {
            delete nodeChild[i];
        }
        nodeChild.clear();
    }

    size_t memoryCost() const;
    size_t numTuples(bool *visited) const;
    void addMatch(VertexID *match, const VertexID *order, ui matchSize, ui startPos, VertexID **data, VertexID *length,
                  ui num);
    void addMatch(const std::vector<VertexID> &match, ui startPos);
    friend bool operator==(const TrieNode& lhs, const TrieNode& rhs) {
        if (lhs.value != rhs.value) return false;

//        // If both nodes have a compression, compare the compression objects
//        if (lhs.compression != nullptr && rhs.compression != nullptr) {
//            if (!(*lhs.compression == *rhs.compression)) return false;
//        } else if (lhs.compression != nullptr || rhs.compression != nullptr) {
//            // One node has a compression and the other doesn't
//            return false;
//        }

        // If no compression, compare children
        if (lhs.nodeChild.size() != rhs.nodeChild.size()) return false;
        for (int i = 0; i < lhs.nodeChild.size(); ++i) {
            if (!(*lhs.nodeChild[i] == *rhs.nodeChild[i])) return false;
        }

        return true;
    }
};

int binarySearch(const DynamicArray<TrieNode *> &array, VertexID v);

#endif //IN_MEMORY_JOIN_RELATION_H
