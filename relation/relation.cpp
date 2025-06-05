//
// Created by Qiyan LI on 2024/2/27.
//

#include "relation.h"

size_t Compression::memoryCost() const {
    size_t mem = sizeof(*this);
    for (int i = 0; i < num; ++i) {
        mem += length[i] * sizeof(VertexID);
        mem += sizeof(VertexID *);
    }

    return mem;
}

Compression::Compression(VertexID **data, ui *length, ui num) {
    this->data = new VertexID *[num];
    this->length = new VertexID [num];
    this->num = num;
    for (int i = 0; i < num; ++i) {
        this->data[i] = new VertexID [length[i]];
        memcpy(this->data[i], data[i], sizeof(VertexID) * length[i]);
        this->length[i] = length[i];
    }
}


size_t TrieNode::memoryCost() const {
    if (nodeChild.empty()) return 0;
    size_t mem = sizeof(*this);
    mem += nodeChild.memoryCost();
//    if (compression) mem += compression->memoryCost();
    for (int i = 0; i < nodeChild.size(); ++i) {
        mem += nodeChild[i]->memoryCost();
    }
    return mem;
}

void TrieNode::addMatch(VertexID *match, const VertexID *order, ui matchSize, ui startPos, VertexID **data, VertexID *length,
                        ui num) {
    TrieNode* current = this;
    while (true) {
        if (startPos + num == matchSize) {
//            if (num != 0) {
//                current->compression = new Compression(data, length, num);
//            }
            break;
        } else {
            VertexID v = match[order[startPos]];
            if (current->nodeChild.empty() || current->nodeChild.back()->value < v) {
                TrieNode *newChild = new TrieNode();
                newChild->value = v;
                current->nodeChild.push_back(newChild);
                current = newChild; // Continue with the new child
            } else {
                current = current->nodeChild.back(); // Move to the last child
            }
            startPos += 1; // Prepare for the next iteration
        }
    }
}

void TrieNode::addMatch(const std::vector<VertexID> &match, ui startPos) {
    TrieNode* current = this;
    while (startPos < match.size()) {
        VertexID v = match[startPos];
        if (current->nodeChild.empty() || current->nodeChild.back()->value < v) {
            TrieNode *newChild = new TrieNode();
            newChild->value = v;
            current->nodeChild.push_back(newChild);
            current = newChild; // Move to the new child
        } else {
            current = current->nodeChild.back(); // Move to the last child
        }
        startPos++; // Move to the next position
    }
}

size_t TrieNode::numTuples(bool *visited) const {
    size_t num = 0;
    if (nodeChild[0]->nodeChild.empty()) {
        for (int i = 0; i < nodeChild.size(); ++i) {
            if (!visited[nodeChild[i]->value])
                ++num;
        }
    }
    else {
        for (int i = 0; i < nodeChild.size(); ++i) {
            if (!visited[nodeChild[i]->value])
                num += nodeChild[i]->numTuples(visited);
        }
    }
    return num;
}

int binarySearch(const DynamicArray<TrieNode *> &array, VertexID v) {
    if (array.size() == 0) return -1;
    int left = 0;
    int right = array.size() - 1;
    while (left <= right) {
        int mid = left + (right - left) / 2;
        if (array[mid]->value == v) {
            return mid;
        } else if (array[mid]->value < v) {
            left = mid + 1;
        } else {
            right = mid - 1;
        }
    }
    return -1; // Element not found
}
