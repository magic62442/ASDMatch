#ifndef TD_HZY_CPP
#define TD_HZY_CPP

#include "hypergraph.h"
#include "td.h"
#include "twbound.h"
#include "td_utils.hpp"
#include "SCsolver.h"
#include "glpk.h"

using namespace std;

FHD::FHD(HyperG & H, vector <size_t> o) {
    vector <VertexSet> E = H.e;
    vector <VertexSet> e;
    for(size_t i = 0; i < o.size(); ++i) {
        VertexSet ei;
        for(size_t j = i + 1; j < o.size(); ++j) {
            for(auto e_: E) 
                if(e_.test(o[i]) && e_.test(o[j])) {
                    ei.Set(o[j]);
                    break;
                }
        }
        e.push_back(ei);
        E.push_back(ei);
        ei.Set(o[i]);
        X.push_back(ei);
    }
    for(size_t i = 0; i < o.size(); ++i)
        for(size_t j = i + 1; j < o.size(); ++j)
            if(X[j].subset(e[i])) {
                eg.push_back(make_pair(i, j));
                break;
            }
}

void FHD::Refine() {
    bool flag = true;
    while(flag) {
        flag = false;
        for(auto ei: eg) {
            size_t p = ei.first, q = ei.second;
            if(X[p].none() || X[q].none() || p == q)
                continue;
            if(X[p].subset(X[q])) { // merge q to p
                flag = true;
                for(auto & ej: eg) {
                    if(ej.first == q)
                        ej.first = p;
                    if(ej.second == q)
                        ej.second = p;
                }
                X[q].clear();
                break;
            }
            else if(X[q].subset(X[p])) {
                swap(p, q);
                flag = true;
                for(auto & ej: eg) {
                    if(ej.first == q)
                        ej.first = p;
                    if(ej.second == q)
                        ej.second = p;
                }
                X[q].clear();
                break;
            }
        }
    }

    vector< pair<size_t, size_t> > _eg;
    vector< VertexSet > _X;
    map<size_t, size_t> remap;
    size_t cnt = 0;

    for(size_t i = 0; i < X.size(); ++i)
        if(!X[i].none()) {
            _X.push_back(X[i]);
            remap[i] = cnt++;
        }
    for(size_t i = 0; i < eg.size(); ++i)
        if(eg[i].first != eg[i].second) {
            pair<size_t, size_t> tmp = make_pair(remap[eg[i].first], remap[eg[i].second]);
            if(tmp.first == tmp.second)
                continue;
            if(tmp.first > tmp.second)
                swap(tmp.first, tmp.second);
            _eg.push_back(tmp);
        }
    
    sort(_eg.begin(), _eg.end()); // de-duplicate edges
    auto last = unique(_eg.begin(), _eg.end());
    _eg.erase(last, _eg.end());

    X = _X;
    eg = _eg;
}

void FineTuneOrder(HyperG & H, Width & ans, vector <VertexSet> & E) {
    vector <size_t> o = ans.o;
    db width = ans.width, w[o.size()];

    while(1) {
        bool updated = 0; // updated flag

#ifdef SWAPNEAR
        for(size_t i = 0; i < o.size() - 1; ++i) {
            swap(o[i], o[i + 1]);
            db tmp = calc_width(H, o, E, w);
            if(tmp < width - eps) {
                width = tmp;
                updated = 1;
            }
            else 
                swap(o[i], o[i + 1]);
        }
#endif

        for(size_t i = 0; i < o.size() || updated ; ++i)
            for(size_t j = i + 1; j < o.size() || updated; ++j) {
                swap(o[i], o[j]);
                db tmp = calc_width(H, o, E, w);
                if(tmp < width - eps) {
                    width = tmp;
                    updated = 1;
                }
                else 
                    swap(o[i], o[j]);
            }

        if(!updated)
            break;
    }
}

// void DivOrderFHD(HyperG & H, Hstate S, vector <size_t> & V, Width & optimal_ans, vector <VertexSet> & E, bool ckf, clock_t start_time) {

// #ifdef TIMEOUT
//     if(clock() - start_time > time_out)
//         return ;
// #endif

//     size_t n = V.size(), k = n / 2;
//     uint64_t comb = (1ULL << k) - 1, x, y;
//     db w_lim = optimal_ans.width;

// #ifdef MEMORIZATION
//     Hvalue temp;
//     if(H.FindState(S)) {
//         temp = H.GetValue(S);
//         comb = temp.comb;
//         x = comb & -comb, y = comb + x;
//         comb = (((comb & ~y) / x ) >> 1) | y;

//         if(temp.optimal_ans < optimal_ans.width - eps) {
//             optimal_ans = temp.optimal_ans;
//             if(ckf) // find a feasible solution
//                 return;
//         }

//         if(n <= tau || ! ( comb < (1ULL << n) ) ) // the state finished search 
//             return ;
//     } 
// #endif

//     if(n <= tau) {
//         Width cur_res(0);
//         EnumOrderFHD(H, V, cur_res, optimal_ans, E, S.second);
//     }
//     else {

// #ifdef DIVAPPROX

//     if(n >= lambda) {
// #ifdef BYCOVER
//         do {
//             vector < pair <db, size_t> > w;
//             VertexSet VS(V);
//             VS.merge(S.second);

//             for(auto it = V.begin(); it != V.end(); it++) {
//                 VertexSet ei = Aggregation(*it, E);
            
//                 ei.intersect(VS);
//                 db current_width = fractional_set_cover(H.N, H.e, ei);
//                 w.push_back( make_pair(current_width, *it) );
//             }

//             sort(w.begin(), w.end());

//             vector <size_t> V_left, V_right;
        
//             for(size_t i = 0; i < k; ++i)
//                 V_left.push_back(w[i].second);
//             for(size_t i = k; i < w.size(); ++i)
//                 V_right.push_back(w[i].second);

//             Width W_left(optimal_ans.width);
//             Hstate _S = S;
//             _S.second = _S.second | (VertexSet) V_right;
//             DivOrderFHD(H, _S, V_left, W_left, E, ckf);

//             if(W_left.size() == 0) // cannot get better solution from left, skip
//                 break;

//             for(auto it = V_left.begin(); it != V_left.end(); it++) { 
//                 // any order of left set will get same connetivity of E
//                 // maybe a Disjoint-set can accelerate here?
//                 VertexSet ei = Aggregation(*it, E);
//                 E.push_back(ei);
//             }

//             Width W_right(W_left.width); // when find a solution not bad than left, pause
//             _S = S;
//             _S.first = _S.first | (VertexSet) V_left;
//             DivOrderFHD(H, _S, V_right, W_right, E, true);

//             if(W_right.size() != 0)
//                 optimal_ans = W_left + W_right;
            
//             for(size_t i = 0; i < k; ++i)
//                 E.pop_back();
//         } while(0);        
// #endif
//     return ;
//     }
// #endif

//         for(; comb < (1ULL << n); comb = (((comb & ~y) / x ) >> 1) | y) {
// #ifdef TIMEOUT
//             if(clock() - start_time > time_out)
//                 return ;
// #endif


//             x = comb & -comb, y = comb + x;
//             VertexSet l(comb); // get left set
//             VertexSet r(((1ULL << n) - 1) ^ comb); // get right set

//             vector <size_t> V_left, V_right;
//             l.getelement(V_left);
//             r.getelement(V_right);
            
//             for(auto it = V_left.begin(); it != V_left.end(); it++)
//                 *it = V[*it];
//             for(auto it = V_right.begin(); it != V_right.end(); it++)
//                 *it = V[*it];
            
//             // Width W_left(optimal_ans.width);
//             Width W_left(w_lim);
//             Hstate _S = S;
//             _S.second = _S.second | (VertexSet) V_right;
//             DivOrderFHD(H, _S, V_left, W_left, E, ckf, start_time);

//             if(W_left.size() == 0) // cannot get better solution from left, skip
//                 continue;

//             for(auto it = V_left.begin(); it != V_left.end(); it++) { 
//                 // any order of left set will get same connetivity of E
//                 // maybe a Disjoint-set can accelerate here?
//                 VertexSet ei = Aggregation(*it, E);
//                 E.push_back(ei);
//             }

//             // Width W_right(W_left.width + 2 * eps); // when find a solution not bad than left, pause
//             // if(ckf)
//             //     W_right.width = optimal_ans.width;
//             // Width W_right(optimal_ans.width);
//             Width W_right(w_lim);

//             _S = S;
//             _S.first = _S.first | (VertexSet) V_left;
//             DivOrderFHD(H, _S, V_right, W_right, E, ckf, start_time); // true due to only need to check on right
//             for(size_t i = 0; i < k; ++i)
//                 E.pop_back();
            
//             if(W_right.size() == 0)
//                 continue;

//             if(W_right.size() != 0 && max(W_left.width, W_right.width) < optimal_ans.width - eps) {
//                 optimal_ans = W_left + W_right;
//                 if(ckf)
//                     break;
//             }
            
            
            
//         }
//     }

//     if(optimal_ans.size() == 0)
//         optimal_ans.width = inf; // There is no solution under now constraint

// #ifdef MEMORIZATION
//     H.AddState(S, Hvalue(optimal_ans, comb));
// #endif
// }

// void EnumOrderFHD(HyperG & H, vector <size_t> & V, Width & current_ans, Width & optimal_ans, vector <VertexSet> & E, VertexSet suffix) {
//     if(V.size() == 0) {
//         optimal_ans = current_ans;
//         return ;
//     }

//     VertexSet VS(V);
//     VS.merge(suffix);

//     for(auto it = V.begin(); it != V.end(); ++it) {
//         size_t v = *it;
//         VertexSet ei = Aggregation(v, E);

//         ei.intersect(VS);

//         db current_width = fractional_set_cover(H.N, H.e, ei);

//         if(current_width < (optimal_ans.width - eps)) {
//             ei.reset(v);
//             E.push_back(ei);
//             swap(*it, V[V.size() - 1]);
//             V.pop_back();

//             db pre_width = current_ans.width;
//             if(pre_width < current_width)
//                 current_ans.set_width(current_width);
//             current_ans.push(v);
// #ifdef OUTPUT_BAG
//             current_ans.Bag.push_back(make_pair(ei, current_width));
// #endif

//             EnumOrderFHD(H, V, current_ans, optimal_ans, E, suffix);

//             E.pop_back();
//             V.push_back(v);
//             swap(*it, V[V.size() - 1]);

//             current_ans.set_width(pre_width);
//             current_ans.pop();
// #ifdef OUTPUT_BAG
//             current_ans.Bag.pop_back();
// #endif
//         }
//     }
// }

void CleanEdge(VertexSet S, vector <VertexSet> & E) {
    vector <VertexSet> res;
    bool flag = 0;

    for(size_t i = 0; i < E.size(); ++i)
        E[i].intersect(S);

    for(size_t i = 0; i < E.size(); ++i) {
        bool fg = 1;
        for(size_t j = 0; j < E.size(); ++j) 
            if(i != j && ((E[i] & E[j]) == E[i])) { // covered
                if((E[i] == E[j]) && i > j)
                    continue;
                fg = 0;
                flag = 1;
                break;
            }
        if(fg)
            res.push_back(E[i]);
    }
    E = res;
}

db FindSimplicalVertex(VertexSet S, vector <VertexSet> & E, SCsolver & Sv, vector <size_t> & res) {
    res.clear();
    vector <size_t> V;
    S.getelement(V);

    VertexSet Ne[MAXELENUM];

    for(size_t i = 0; i < V.size(); ++i) {
        size_t u = V[i];
        Ne[u] = Aggregation(u, E);
        Ne[u].intersect(S);
        Ne[u].reset(u);
    }

    for(size_t i = 0; i < V.size(); ++i) { // find simplical vertex
        size_t u = V[i];
        vector <size_t> Vi;
        Ne[u].getelement(Vi);
        bool flag = 1;
        for(size_t j = 0; j < Vi.size(); ++j)
            if(!Ne[Vi[j]].subset(Ne[u])) {
                flag = 0;
                break;
            }
        if(flag)
            res.push_back(u);
    }

    db low = 0;
    for(size_t i = 0; i < res.size(); ++i)
        low = max(low, Sv.solve(Ne[res[i]]));
    return low;
}

db DPFHD(HyperG & H, FILE * f, Order & up_o) {
    map <VertexSet, db> FHW[H.N + 1];
    map <VertexSet, Order> o[H.N + 1];
    map <VertexSet, pair<db, db> > Bnd[H.N + 1];
    map <VertexSet, Order> UB_o[H.N + 1];
    VertexSet Uniset;
    vector <size_t> SimplicalVertex;
    // FHW[0][Uniset] = 0;
    for(size_t i = 0; i < H.N; ++i)
        Uniset.Set(i);

    SCsolver Sv(H.e);

    db up = Sv.solve(Uniset), bt;
    for(size_t i = 0; i < H.N; ++i)
        up_o.push_back(i);
    // db bt = Sv.solve(C);
    // Uniset.Xor(C); // arrange clique C to the end


#ifdef PRINTSTATENUMBER
    size_t total_state = 0;
#endif

#ifdef SLMPV
    db low = FindSimplicalVertex(Uniset, H.e, Sv, SimplicalVertex);
    VertexSet InitS(SimplicalVertex);
    size_t cnt = SimplicalVertex.size();
#else
    db low = 0;
    VertexSet InitS;
    size_t cnt = 0;
#endif

    VertexSet C = MaximalClique(H, InitS);
    bt = max(1., max(Sv.solve(C), low));

    FHW[cnt][InitS] = max(low, 1.);
    o[cnt][InitS] = Order();

#ifdef TIMEOUT
    clock_t begin = clock();
#endif

    for(size_t i = cnt; i < H.N; ++i) {
#ifdef PRINTSTATENUMBER
        total_state += FHW[i].size();
#endif
        if(i) {
            FHW[i - 1].clear(); // shrink memory
            o[i - 1].clear();
            Bnd[i - 1].clear();
            UB_o[i - 1].clear();
        }
        for(auto it = FHW[i].begin(); it != FHW[i].end(); it++) {
            pair <db, db> bnd;
            Order tmp_o;
            if(Bnd[i].count(it->first))
                bnd = Bnd[i][it->first];
            else
                bnd = make_pair(0, H.N + 1);

            if(it->second >= up - eps || it->second >= bnd.second + eps)
                continue;

#ifdef TIMEOUT
            if(clock() - begin >= time_out) {
                fprintf(f, "Timeout terminate\n");
#ifdef PRINTSTATENUMBER
                fprintf(f, "total_state_number = %lu\n", total_state);
#endif
                return up;
            }
#endif

            VertexSet S = it->first;
            vector <size_t> V;
            vector <VertexSet> E = H.e;

            for(size_t j = 0; j < H.N; ++j) 
                if(S.test(j)) {
                    VertexSet ei = Aggregation(j, E);
                    ei.reset(j);
                    E.push_back(ei);
                }
                else if(!C.test(j))
                    V.push_back(j);

            CleanEdge(S ^ Uniset, E);
#ifdef SLMPV
            low = FindSimplicalVertex(S ^ Uniset ^ C, E, Sv, SimplicalVertex);
            if(low > 0) {
                VertexSet tmpS = S;
                tmpS.merge(SimplicalVertex);
                low = max(low, it->second);
                cnt = SimplicalVertex.size();  

                if(low >= up - eps)
                    continue;
                if(FHW[i + cnt].count(tmpS) && FHW[i + cnt][tmpS] <= low + eps)
                    continue;

                FHW[i + cnt][tmpS] = low;
                continue;
            }
#endif
            
            for(auto vit = V.begin(); vit != V.end(); ++vit) {
                VertexSet tmpS = S;
                tmpS.Set(*vit);
                VertexSet CS = tmpS ^ Uniset;
                tmp_o.clear();

                VertexSet ei = Aggregation(*vit, E);
                ei.intersect(CS);
                ei.Set(*vit);
                db w = Sv.solve(ei);
                if(w >= up - eps)
                    continue;

                if(FHW[i + 1].count(tmpS) && FHW[i + 1][tmpS] <= w + eps)
                    continue;
                
                ei.reset(*vit);
                E.push_back(ei);

                bnd = make_pair(0, H.N + 1);
                if(Bnd[i + 1].count(tmpS))
                    bnd = Bnd[i + 1][tmpS];
                if(UB_o[i + 1].count(tmpS))
                    tmp_o = UB_o[i + 1][tmpS];

                do {
                    if(bnd.first < 1)
                        bnd.first = FHW_lb(CS, E, Sv);

                    if(bnd.first >= up - eps) {
                        Bnd[i + 1][tmpS] = bnd;
                        break;
                    }

                    if(bnd.second > H.N)
                        bnd.second = FHW_ub(CS, E, Sv, tmp_o);

                    Bnd[i + 1][tmpS] = bnd;
                    UB_o[i + 1][tmpS] = tmp_o;

                    if(max(w, bnd.second) <= up - eps) {
                        up = max(w, bnd.second);
                        up_o = o[i][S];
                        up_o.push_back(*vit);
                        up_o.insert(up_o.end(), tmp_o.begin(), tmp_o.end());
                        fprintf(f, "update ub = %lf update_time: %lfs\n", up, (db) (clock() - begin) / CLOCKS_PER_SEC);
                        // fprintf(f, "%lu\n", i);
                        fprintf(f, "update elim order: ");
                        for(auto v: up_o) 
                            fprintf(f, "%lu ", v);
                        fprintf(f, "\n");
#ifdef TIMEOUT
                        begin = clock();
#endif
                        if(up + eps <= bt)
                            return up;
                    }

                    // if(bnd.second <= w + eps) 
                    //     break;

                    FHW[i + 1][tmpS] = w;
                    o[i + 1][tmpS] = o[i][S];
                    o[i + 1][tmpS].push_back(*vit);
                } while(0);
                E.pop_back();
            }
        }
    }
#ifdef PRINTSTATENUMBER
    fprintf(f, "total_state_number = %lu\n", total_state);
#endif

#ifdef TIMEOUT
    fprintf(f, "running_time: %lfs\n", (db) (clock() - begin) / CLOCKS_PER_SEC);
#endif
    return up;
}

db DPFHD(HyperG & H, Order & up_o) {
    map <VertexSet, db> FHW[H.N + 1];
    map <VertexSet, Order> o[H.N + 1];
    map <VertexSet, pair<db, db> > Bnd[H.N + 1];
    map <VertexSet, Order> UB_o[H.N + 1];
    VertexSet Uniset;
    vector <size_t> SimplicalVertex;
    // FHW[0][Uniset] = 0;
    for(size_t i = 0; i < H.N; ++i)
        Uniset.Set(i);

    SCsolver Sv(H.e);

    db up = Sv.solve(Uniset), bt;
    for(size_t i = 0; i < H.N; ++i)
        up_o.push_back(i);
    // db bt = Sv.solve(C);
    // Uniset.Xor(C); // arrange clique C to the end


#ifdef PRINTSTATENUMBER
    size_t total_state = 0;
#endif

#ifdef SLMPV
    db low = FindSimplicalVertex(Uniset, H.e, Sv, SimplicalVertex);
    VertexSet InitS(SimplicalVertex);
    size_t cnt = SimplicalVertex.size();
#else
    db low = 0;
    VertexSet InitS;
    size_t cnt = 0;
#endif

    VertexSet C = MaximalClique(H, InitS);
    bt = max(1., max(Sv.solve(C), low));

    FHW[cnt][InitS] = max(low, 1.);
    o[cnt][InitS] = Order();

#ifdef TIMEOUT
    clock_t begin = clock();
#endif

    for(size_t i = cnt; i < H.N; ++i) {
#ifdef PRINTSTATENUMBER
        total_state += FHW[i].size();
#endif
        if(i) {
            FHW[i - 1].clear(); // shrink memory
            o[i - 1].clear();
            Bnd[i - 1].clear();
            UB_o[i - 1].clear();
        }
        for(auto it = FHW[i].begin(); it != FHW[i].end(); it++) {
            pair <db, db> bnd;
            Order tmp_o;
            if(Bnd[i].count(it->first))
                bnd = Bnd[i][it->first];
            else
                bnd = make_pair(0, H.N + 1);

            if(it->second >= up - eps || it->second >= bnd.second + eps)
                continue;

#ifdef TIMEOUT
            if(clock() - begin >= time_out) {
                return up;
            }
#endif

            VertexSet S = it->first;
            vector <size_t> V;
            vector <VertexSet> E = H.e;

            for(size_t j = 0; j < H.N; ++j)
                if(S.test(j)) {
                    VertexSet ei = Aggregation(j, E);
                    ei.reset(j);
                    E.push_back(ei);
                }
                else if(!C.test(j))
                    V.push_back(j);

            CleanEdge(S ^ Uniset, E);
#ifdef SLMPV
            low = FindSimplicalVertex(S ^ Uniset ^ C, E, Sv, SimplicalVertex);
            if(low > 0) {
                VertexSet tmpS = S;
                tmpS.merge(SimplicalVertex);
                low = max(low, it->second);
                cnt = SimplicalVertex.size();

                if(low >= up - eps)
                    continue;
                if(FHW[i + cnt].count(tmpS) && FHW[i + cnt][tmpS] <= low + eps)
                    continue;

                FHW[i + cnt][tmpS] = low;
                continue;
            }
#endif

            for(auto vit = V.begin(); vit != V.end(); ++vit) {
                VertexSet tmpS = S;
                tmpS.Set(*vit);
                VertexSet CS = tmpS ^ Uniset;
                tmp_o.clear();

                VertexSet ei = Aggregation(*vit, E);
                ei.intersect(CS);
                ei.Set(*vit);
                db w = Sv.solve(ei);
                if(w >= up - eps)
                    continue;

                if(FHW[i + 1].count(tmpS) && FHW[i + 1][tmpS] <= w + eps)
                    continue;

                ei.reset(*vit);
                E.push_back(ei);

                bnd = make_pair(0, H.N + 1);
                if(Bnd[i + 1].count(tmpS))
                    bnd = Bnd[i + 1][tmpS];
                if(UB_o[i + 1].count(tmpS))
                    tmp_o = UB_o[i + 1][tmpS];

                do {
                    if(bnd.first < 1)
                        bnd.first = FHW_lb(CS, E, Sv);

                    if(bnd.first >= up - eps) {
                        Bnd[i + 1][tmpS] = bnd;
                        break;
                    }

                    if(bnd.second > H.N)
                        bnd.second = FHW_ub(CS, E, Sv, tmp_o);

                    Bnd[i + 1][tmpS] = bnd;
                    UB_o[i + 1][tmpS] = tmp_o;

                    if(max(w, bnd.second) <= up - eps) {
                        up = max(w, bnd.second);
                        up_o = o[i][S];
                        up_o.push_back(*vit);
                        up_o.insert(up_o.end(), tmp_o.begin(), tmp_o.end());
#ifdef TIMEOUT
                        begin = clock();
#endif
                        if(up + eps <= bt)
                            return up;
                    }

                    // if(bnd.second <= w + eps)
                    //     break;

                    FHW[i + 1][tmpS] = w;
                    o[i + 1][tmpS] = o[i][S];
                    o[i + 1][tmpS].push_back(*vit);
                } while(0);
                E.pop_back();
            }
        }
    }
    return up;
}

#endif