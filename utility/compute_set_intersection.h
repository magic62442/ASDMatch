//
// Created by Qiyan LI on 2024/3/8.
//

#ifndef IN_MEMORY_JOIN_COMPUTE_SET_INTERSECTION_H
#define IN_MEMORY_JOIN_COMPUTE_SET_INTERSECTION_H


#include "config.h"
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#ifndef __APPLE__
#include <immintrin.h>
#include <x86intrin.h>
#endif

extern size_t gNumInterSection;

/*
 * Because the set intersection is designed for computing common neighbors, the target is uieger.
 */

class ComputeSetIntersection {
public:

    static void ComputeCandidates(const VertexID* larray, ui l_count, const VertexID* rarray,
                                  ui r_count, VertexID* cn, ui &cn_count);
    static void ComputeCandidates(const VertexID* larray, ui l_count, const VertexID* rarray,
                                  ui r_count, ui &cn_count);
    // actually an algorithm for sorted set intersection. LeapFrogTrieJoin is the join algorithm.
    // See Section 3.1 in the LeapFrog paper. only used for joining edge relations.
    static ui LeapfrogSeek(const VertexID *src, ui begin, ui end, ui target);
    static void LeapfrogJoin(VertexID **arrays, ui *counts, ui num, VertexID *cn, ui &cn_count);
    static const ui GallopingSearch(const VertexID *src, ui begin, ui end, ui target);
    static ui BinarySearch(const VertexID *src, ui begin, ui end, ui target);
    static ui BinarySearch(const std::vector<VertexID> &array, VertexID target);

#if SI == 0
    static void ComputeCNGallopingAVX2(const VertexID* larray, ui l_count,
                                       const VertexID* rarray, ui r_count, VertexID* cn,
                                       ui &cn_count);
    static void ComputeCNGallopingAVX2(const VertexID* larray, ui l_count,
                                       const VertexID* rarray, ui r_count, ui &cn_count);

    static void ComputeCNMergeBasedAVX2(const VertexID* larray, ui l_count, const VertexID* rarray,
                                        ui r_count, VertexID* cn, ui &cn_count);
    static void ComputeCNMergeBasedAVX2(const VertexID* larray, ui l_count, const VertexID* rarray,
                                        ui r_count, ui &cn_count);
    static const ui BinarySearchForGallopingSearchAVX2(const VertexID*  array, ui offset_beg, ui offset_end, ui val);
    static const ui GallopingSearchAVX2(const VertexID*  array, ui offset_beg, ui offset_end, ui val);
#elif SI == 1

    static void ComputeCNGallopingAVX512(const VertexID* larray, const ui l_count,
                                         const VertexID* rarray, const ui r_count, VertexID* cn,
                                         ui &cn_count);
    static void ComputeCNGallopingAVX512(const VertexID* larray, const ui l_count,
                                         const VertexID* rarray, const ui r_count, ui &cn_count);

    static void ComputeCNMergeBasedAVX512(const VertexID* larray, const ui l_count, const VertexID* rarray,
                                          const ui r_count, VertexID* cn, ui &cn_count);
    static void ComputeCNMergeBasedAVX512(const VertexID* larray, const ui l_count, const VertexID* rarray,
                                          const ui r_count, ui &cn_count);

#elif SI == 2

    static void ComputeCNNaiveStdMerge(const VertexID* larray, ui l_count, const VertexID* rarray,
                                       ui r_count, VertexID* cn, ui &cn_count);
    static void ComputeCNNaiveStdMerge(const VertexID* larray, ui l_count, const VertexID* rarray,
                                       ui r_count, ui &cn_count);
    static void ComputeCNGalloping(const VertexID * larray, ui l_count, const VertexID * rarray,
                                   ui r_count, VertexID * cn, ui& cn_count);
    static void ComputeCNGalloping(const VertexID * larray, ui l_count, const VertexID * rarray,
                                   ui r_count, ui& cn_count);
#endif
};


#endif //IN_MEMORY_JOIN_COMPUTE_SET_INTERSECTION_H
