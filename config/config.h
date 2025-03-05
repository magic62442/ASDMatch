//
// Created by anonymous authors on 2024/2/27.
//

#ifndef IN_MEMORY_JOIN_CONFIG_H
#define IN_MEMORY_JOIN_CONFIG_H

#include <cstdlib>
#include <utility>
#include <cstdint>

/**
 * Set intersection method.
 * 0: Hybrid method; 1: Merge based set intersections.
 */
#define HYBRID 0

/**
 * Accelerate set intersection with SIMD instructions.
 * 0: AVX2; 1: AVX512; 2: Basic;
 * If building in mac AMD chip, can only use 2.
 * For machines using Intel CPU, use 0.
 */

#ifdef __APPLE__
#define SI 2
#else
#define SI 0
#endif

#define MAX_PATTERN_SIZE 32
#define MAX_NUM_LABEL 1000
#define MAX_FOURCYCLE_NUM 1e10
#define SAMPLE_SIZE 10000
#define SAMPLE_PORTION_LABELED 6
#define MAX_PENALTY 0.95
#define COLLECT_STATISTICS
//#define COLLET_GLOBAL_TIME
#define MAX_NUM_PLAN 200
//#define COLLECT_RESULT

typedef uint32_t ui;
typedef uint32_t VertexID;
typedef uint32_t EdgeID;
typedef uint32_t LabelID;
typedef uint64_t Count;

#endif //IN_MEMORY_JOIN_CONFIG_H
