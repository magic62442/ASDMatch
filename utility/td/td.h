#ifndef TD_HZY_H
#define TD_HZY_H
#define db double

#include "hypergraph.h"

struct FHD {
    std::vector< std::pair<size_t, size_t> > eg;
    std::vector< VertexSet > X;

    void Refine();

    FHD() {}
    FHD(HyperG & H, std::vector <size_t> o);
};

void EnumOrderFHD(HyperG & H, std::vector <size_t> & V, Width & current_ans, Width & optimal_ans, std::vector <VertexSet> & E, VertexSet suffix);

/*
    Divide and conquer set method for fhd:
        H: given hypergraph
        Hstate: prefix and suffix, a state of divide and conquer
        V: divided vertex set, each iteration will be divided up
        optimal_ans: maintain the optimal answer from now on
        E: edge set, maintain connectivity and aggregation operation
        ckf: check process or not, if true then break when find a feasible solution

*/

void DivOrderFHD(HyperG & H, Hstate S, std::vector <size_t> & V, Width & optimal_ans, std::vector <VertexSet> & E, bool ckf, clock_t start_time);

db DPFHD(HyperG & H, FILE * f, Order & up_o);
db DPFHD(HyperG & H, Order & up_o);


#endif