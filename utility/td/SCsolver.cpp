#ifndef SCSOLVER_HZY_CPP
#define SCSOLVER_HZY_CPP

#include "SCsolver.h"
#include "td_utils.hpp"
#include "glpk.h"

/*
    FSC solve by LP with glpk library:
        variation range: [l, r]
        set collection: E[][l~r]
        cover set: e[l~r]
    
    LP:
        target: 
            min ||x_i||_1
        subject to:
            E[][l~r] x >= e
            x >-= 0
        
*/

double fractional_set_cover(size_t n, vector< VertexSet > & E, VertexSet e) {
    // return 1.;
    glp_term_out(GLP_OFF);

    glp_prob *lp;
    size_t m = E.size(); // number of xi, col
    // size_t n = r + 1 - l; // number of row
    int ia[n * m + 1], ja[n * m + 1];
    double ar[n * m + 1], z/*, x[m + 1]*/;

    lp = glp_create_prob();
    glp_set_obj_dir(lp, GLP_MIN);
    glp_add_rows(lp, n);

    for(size_t i = 0; i < n; ++i)
        if(e.test(i))
            glp_set_row_bnds(lp, i + 1, GLP_LO, 1., inf);
        else
            glp_set_row_bnds(lp, i + 1, GLP_LO, 0., inf);
    glp_add_cols(lp, m);
    for(size_t i = 0; i < m; ++i) {
        glp_set_col_bnds(lp, i + 1, GLP_LO, 0., 0.);
        glp_set_obj_coef(lp, i + 1, 1.);
    }

    int cnt = 0;
    for(size_t i = 0; i < n; ++i)
        for(size_t j = 0; j < m; ++j)
            ia[++cnt] = i + 1, ja[cnt] = j + 1, ar[cnt] = (E[j].test(i) ? 1. : 0.);
    glp_load_matrix(lp, cnt, ia, ja, ar);
    glp_simplex(lp, NULL);

    z = glp_get_obj_val(lp);
    glp_delete_prob(lp);

    return z;
}

double calc_width(HyperG & H, vector <size_t> o, vector <VertexSet> E, double * w) {
    VertexSet VS(o);
    db res(0);
    int cnt = 0;

    for(auto it = o.begin(); it != o.end(); ++it) {
        VertexSet ei = Aggregation(*it, E);
        ei.intersect(VS);
        db temp = fractional_set_cover(H.N, H.e, ei);
        w[cnt++] = temp;
        res = max(res, temp);
        ei.reset(*it);
        VS.reset(*it);
        E.push_back(ei);
    }
    return res;
}

double FSC_LP(vector< VertexSet > & E, VertexSet e) {
    if(e.none())
        return 0.;

    vector < VertexSet > newE, temp;
    VertexSet tmp;

    for(auto it = E.begin(); it != E.end(); it++) {
        tmp = *it;
        tmp.intersect(e);
        temp.push_back(tmp);
    }
    sort(temp.begin(), temp.end());
    for(auto it = temp.begin(); it != temp.end(); ++it) {
        auto itt = it;
        bool flag = 1;
        itt++;
        while(itt != temp.end()) {
            if((*itt).subset(*it)) {
                flag = 0;
                break;
            }
            itt++;
        }
        if(flag)
            newE.push_back(*it);
    }

    glp_term_out(GLP_OFF);

    glp_prob *lp;
    size_t m = newE.size(); 
    vector < size_t > V;
    e.getelement(V);
    size_t n = V.size();
    int ia[n * m + 1], ja[n * m + 1];
    double ar[n * m + 1], z;

    lp = glp_create_prob();
    glp_set_obj_dir(lp, GLP_MIN);
    glp_add_rows(lp, n);

    for(size_t i = 0; i < n; ++i)
        glp_set_row_bnds(lp, i + 1, GLP_LO, 1., inf);
    glp_add_cols(lp, m);
    for(size_t i = 0; i < m; ++i) {
        glp_set_col_bnds(lp, i + 1, GLP_LO, 0., 0.);
        glp_set_obj_coef(lp, i + 1, 1.);
    }

    int cnt = 0;
    for(size_t i = 0; i < n; ++i)
        for(size_t j = 0; j < m; ++j)
            ia[++cnt] = i + 1, ja[cnt] = j + 1, ar[cnt] = (newE[j].test(V[i]) ? 1. : 0.);
    glp_load_matrix(lp, cnt, ia, ja, ar);
    glp_simplex(lp, NULL);

    z = glp_get_obj_val(lp);
    glp_delete_prob(lp);

    return z;
}

double SCsolver::solve(VertexSet S) {
    if(res.count(S))
        return res[S];
    double ans = FSC_LP(E, S);
    res[S] = ans;
    return ans;
}
#endif