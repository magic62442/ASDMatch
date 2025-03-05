#ifndef HYPERGRAPH_HZY_H
#define HYPERGRAPH_HZY_H

#include <cstdlib>
#include <string.h>
#include <map>
#include <vector>
// #include <bitmapset.h>
#include <set>
#include <bitset>
#include <algorithm>
#include <unordered_map>
#include "global_setting.h"

#define db double
#define LOG2(X) ((unsigned) (8*sizeof (unsigned long long) - __builtin_clzll((X)) - 1))

// using namespace std;


const uint64_t m1  = 0x5555555555555555; //binary: 0101...
const uint64_t m2  = 0x3333333333333333; //binary: 00110011..
const uint64_t m4  = 0x0f0f0f0f0f0f0f0f; //binary:  4 zeros,  4 ones ...
const uint64_t m8  = 0x00ff00ff00ff00ff; //binary:  8 zeros,  8 ones ...
const uint64_t m16 = 0x0000ffff0000ffff; //binary: 16 zeros, 16 ones ...
const uint64_t m32 = 0x00000000ffffffff; //binary: 32 zeros, 32 ones
const uint64_t h01 = 0x0101010101010101; //the sum of 256 to the power of 0,1,2,3...

struct VertexSet ;
struct Width;
struct HyperG;
struct Hvalue;

typedef std::pair<VertexSet, VertexSet> Hstate;
typedef std::vector<size_t> Order;

struct Width {
    db width;
    std::vector <size_t> o;
    std::vector < std::pair<VertexSet, db> > Bag;

    Width(): width(1e7) {}
    Width(db w): width(w) {}
    Width(db w, const std::vector <size_t> & o_) {
        width = w;
        o = o_;
    }

    void push(size_t x) {
        o.push_back(x);
    }
    void pop() {
        o.pop_back();
    }
    void set_width(db _w) {
        width = _w;
    }

    bool operator < (const Width & w) const {
        return width < w.width;
    }

    Width operator - (db & x) const {
        return Width(width - x);
    }

    Width operator + (Width & x) const {
        Width res(std::max(width, x.width));
        res.o = o;
        res.o.insert(res.o.end(), x.o.begin(), x.o.end());
#ifdef OUTPUT_BAG
        res.Bag = Bag;
        res.Bag.insert(res.Bag.end(), x.Bag.begin(), x.Bag.end());
#endif
        return res;
    }

    size_t size() { return o.size(); }
};

struct VertexSet {
    uint64_t S[UINT64NUM];

    VertexSet() {
        memset(S, 0, sizeof S);
    }
    VertexSet(uint64_t * _s) {
        for(size_t i = 0; i < UINT64NUM; ++i)
            S[i] = _s[i];
    }
    VertexSet(std::vector <size_t> & v) {
        memset(S, 0, sizeof S);
        for(auto it = v.begin(); it != v.end(); ++it)
            this->Set(*it);
    }
    VertexSet(std::set <size_t> & v) {
        memset(S, 0, sizeof S);
        for(auto it = v.begin(); it != v.end(); ++it)
            this->Set(*it);
    }

    void intersect (const VertexSet & V) {
        for(size_t i = 0; i < UINT64NUM; ++i)
            S[i] &= V.S[i];
    }

    void merge (const VertexSet & V) {
        for(size_t i = 0; i < UINT64NUM; ++i)
            S[i] |= V.S[i];
    }

    void Xor (const VertexSet & V) {
        for(size_t i = 0; i < UINT64NUM; ++i)
            S[i] ^= V.S[i];
    }

    void reset(const VertexSet & V) {
        (*this).merge(V);
        (*this).Xor(V);
    }

    void clear() {
        for(size_t i = 0; i < UINT64NUM; ++i)
            S[i] = 0;
    }

    bool test(size_t pos) {
        size_t div = pos / 64, mod = pos % 64;
        if(div >= UINT64NUM) return false;
        return (S[div] & (1ULL << mod)) ? 1 : 0;
    }

    void Set(size_t pos) {
        size_t div = pos / 64, mod = pos % 64;
        if(div >= UINT64NUM) return ;
        S[div] |= (1ULL << mod);
    }

    void reset(size_t pos) {
        size_t div = pos / 64, mod = pos % 64;
        if(div >= UINT64NUM) return ;
        S[div] = (S[div] | (1ULL << mod)) ^ (1ULL << mod);
    }

    void getelement(std::vector <size_t> & V) {
        V.clear();
        for(size_t i = 0; i < UINT64NUM; ++i) {
            uint64_t x = S[i], y;
            for(; x; x = y) {
                y = (x & (x - 1));
                V.push_back(i * 64 + LOG2(x ^ y));
            }
        }
    }

    void rev() {
        for(size_t i = 0; i < UINT64NUM; ++i)
            S[i] = ~S[i];
    }

    bool subset(VertexSet V) {
        bool flag = 1;
        for(size_t i = 0; i < UINT64NUM; ++i)
            flag &= ((V.S[i] & S[i]) == V.S[i]);
        return flag;
    }

    bool none() {
        for(size_t i = 0; i < UINT64NUM; ++i)
            if(S[i] != 0)
                return false;
        return true;
    }

    size_t size() {
        size_t sum = 0;
        for(size_t i = 0; i < UINT64NUM; ++i) {
            uint64_t x = S[i];
            x -= (x >> 1) & m1;             //put count of each 2 bits into those 2 bits
            x = (x & m2) + ((x >> 2) & m2); //put count of each 4 bits into those 4 bits 
            x = (x + (x >> 4)) & m4;        //put count of each 8 bits into those 8 bits 
            sum += (x * h01) >> 56;  //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... 
        }
        return sum;
    }

    VertexSet operator & (const VertexSet & V) {
        VertexSet temp = S;
        temp.intersect(V);
        return temp;
    }
    VertexSet operator | (const VertexSet & V) {
        VertexSet temp = S;
        temp.merge(V);
        return temp;
    }
    VertexSet operator ^ (const VertexSet & V) {
        VertexSet temp = S;
        temp.Xor(V);
        return temp;
    }
    bool operator == (const VertexSet & V) {
        for(size_t i = 0; i < UINT64NUM; ++i)
            if(S[i] != V.S[i])
                return false;
        return true;
    }
    bool operator != (const VertexSet & V) {
        for(size_t i = 0; i < UINT64NUM; ++i)
            if(S[i] != V.S[i])
                return true;
        return false;
    }
    bool operator == (const uint64_t * x) {
        for(size_t i = 0; i < UINT64NUM; ++i)
            if(S[i] != x[i])
                return false;
        return true;
    }
    bool operator < (const VertexSet & V) const {
        for(size_t i = 0; i < UINT64NUM; ++i)
            if(S[i] < V.S[i])
                return true;
            else if(S[i] > V.S[i])
                return false;
        return false;
    }
};


struct Hvalue {
    Width optimal_ans;
    uint64_t comb;

    Hvalue() {}
    Hvalue(db width, uint64_t _comb) : optimal_ans(width), comb(_comb) {}
    Hvalue(Width width, uint64_t _comb) : optimal_ans(width), comb(_comb) {}
};
struct HyperG {
    /*
        use bitset to maintain hyperedge, assume the element has been numbered from 0~N-1
        N: number of element (constant)
        M: number of hyperedge
        e: hyperedge set
        __PrimalG: primal graph 
    */
    size_t N, M;
    std::vector <VertexSet> e;
    std::vector <VertexSet> __PrimalG;

    HyperG() : N(0) {
        M = 0;
    }
    HyperG(size_t N_): N(N_) {}
    HyperG(size_t N_, std::vector <std::set<std::string> > & v, std::map <std::string, size_t> & f);
    HyperG(VertexSet S, std::vector<VertexSet> & v);

    ~HyperG() {
        e.clear();
        __PrimalG.clear();
    }

    std::vector <VertexSet> & PrimalG();
    HyperG induced(VertexSet S);
};



#endif