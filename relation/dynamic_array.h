//
// Created by anonymous authors on 2024/3/13.
//

#ifndef IN_MEMORY_JOIN_DYNAMIC_ARRAY_H
#define IN_MEMORY_JOIN_DYNAMIC_ARRAY_H

#include <cstdlib>
#include <algorithm>
#include "config.h"

template <typename T>
class DynamicArray {
private:
    T* _data;
    ui _capacity;
    ui _start; // Logical _start index
    ui _size;  // Logical _size of the array

    void grow() {
        ui newCapacity = 2 * (_capacity - _start);
        if (newCapacity == 0) newCapacity = 1;
        T* newData = new T[newCapacity];
        std::copy(_data + _start, _data + _start + _size, newData);
        delete[] _data; // Deallocate old array
        _data = newData;
        _capacity = newCapacity;
        _start = 0; // Reset _start as we've compacted the array
    }

public:
    DynamicArray(): _data(nullptr), _capacity(0), _start(0), _size(0) {}

    ~DynamicArray() {
        clear();
    }

    void push_back(const T& value) {
        if (_size + _start == _capacity) {
            grow();
        }
        _data[_start + _size] = value;
        _size++;
    }

    void pop_front() {
        if (_size == 0) return;
        // Call destructor for the element being "popped" if it's not a trivially destructible type
        _data[_start].~T();
        _start++;
        _size--;
    }

    void clear() {
        delete[] _data;
        _start = _size = _capacity = 0;
        _data = nullptr;
    }

    T& operator[](ui index) const {
        return _data[_start + index];
    }

    T& back() const {
        return _data[_start + _size - 1];
    }

    ui size() const { return _size; }

    void setSize(ui size) { _size = size; }

    bool empty() const { return _size == 0; }

    const T* data() const {return _data; }

    size_t memoryCost() const {
        size_t objectMemory = sizeof(*this); // Memory used by the DynamicArray object itself
        size_t elementsMemory = _capacity * sizeof(T); // Memory allocated for the elements
        return objectMemory + elementsMemory;
    }

    // Deleted copy constructor and copy assignment to prevent copying
    DynamicArray(const DynamicArray& other) = delete;
    DynamicArray& operator=(const DynamicArray& other) = delete;
    // Friend declaration of the == operator
    friend bool operator==(const DynamicArray<T>& lhs, const DynamicArray<T>& rhs) {
        // First, check if the sizes are the same
        if (lhs.size() != rhs.size()) return false;

        // Then compare each element
        for (uint32_t i = 0; i < lhs.size(); ++i) {
            if (lhs[i] != rhs[i]) return false;
        }

        return true;
    }
};


#endif //IN_MEMORY_JOIN_DYNAMIC_ARRAY_H
