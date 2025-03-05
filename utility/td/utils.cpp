#ifndef UTILS_HZY_CPP
#define UTILS_HZY_CPP

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <set>
#include <stack>
#include <algorithm>
#include <pthread.h>
#include <signal.h>
#include <glpk.h>
#include "hypergraph.h"

#include <sys/types.h>
#include <sys/wait.h>

#if __linux__ || __APPLE__
#include <unistd.h>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h> 
#elif __WIN32__
#include <direct.h>
#include <io.h>
#endif

using namespace std;

struct thread_queue {
    pthread_t tids[worker_num];

    thread_queue() {memset(tids, 0, sizeof tids);}

    pthread_t * alloc() {
        for(size_t i = 0; i < worker_num; ++i)
            if(tids[i] == 0)
                return tids + i;
            else {
                int kill_rc = pthread_kill(tids[i], 0);
                if(kill_rc == ESRCH || kill_rc == EINVAL)
                    return tids + i;
            }
        return NULL;
    }
    bool finish() {
        for(size_t i = 0; i < worker_num; ++i)
            if(tids[i] != 0) {
                int kill_rc = pthread_kill(tids[i], 0);
                if(kill_rc != ESRCH && kill_rc != EINVAL)
                    return 0;
            }
        return 1;
    }
};

inline bool is_empty_file(std::ifstream & pFile) {
    return pFile.peek() == std::ifstream::traits_type::eof();
}

inline bool is_timeout(std::ifstream &pFile) {
    std::string temp;
    while(!pFile.eof()) {
        pFile >> temp;
        if(temp == "Timeout")
            return true;
    }
    return false;
}

inline vector <string> data_loader(const char * data_dir, const char * load_type) {
    vector <string> file_name;

#if __linux__ || __APPLE__
    DIR *dir;
    struct dirent * ptr;
    // char base[1000];
    bool append_mode = false, retest_mode = false, fulltest_mode = false;
    if(strcmp(load_type, "append") == 0)
        append_mode = true;
    if(strcmp(load_type, "retest") == 0)
        retest_mode = true;
    if(strcmp(load_type, "full") == 0)
        fulltest_mode = true;
    

    if( (dir = opendir(data_dir) ) == NULL ) {
        perror("Open dir error...");
        exit(1);
    }

    while( (ptr = readdir(dir)) != NULL ) {
        if(strcmp(ptr->d_name, ".") == 0 || strcmp(ptr->d_name, "..") == 0)
            continue;
        else if(ptr->d_type == 8) {// file
            if(fulltest_mode)
                file_name.push_back(ptr->d_name);
            else if(append_mode) {
                std::ifstream pfile(ptr->d_name);
                if(is_empty_file(pfile))
                    file_name.push_back(ptr->d_name);
            }
            else if(retest_mode) {
                std::ifstream pfile(ptr->d_name);
                if(is_empty_file(pfile) || is_timeout(pfile))
                    file_name.push_back(ptr->d_name);
            }
        }
        else if(ptr->d_type == 10) // link file
            continue;
        else if(ptr->d_type == 4) {// dir
            // file_name.push_back(ptr->d_name);
            continue; // do not handle recursive
        }
    }
    closedir(dir);

#endif
    sort(file_name.begin(), file_name.end());
    return file_name;
}

inline vector < set <size_t> > BiconnectComp(vector < vector <bool> > & AdjM);
inline vector < vector <bool> > GetAdjM(size_t n, vector < set < size_t > > & E);
inline vector < vector <bool> > GetAdjM(size_t n, vector <VertexSet > & E);
inline vector < set < size_t > > Numbered(vector < set < string > > & v, map < string, size_t > & f);
inline vector < set < size_t > > InducedHyperG(vector < set < size_t > > & E, set < size_t > V) {
    map < size_t, size_t > M;
    size_t i = 0;

    for(auto it = V.begin(); it != V.end(); ++it)
        M[*it] = i++;

    vector < set < size_t > > IndE;

    for(auto it = E.begin(); it != E.end(); it++) {
        set < size_t > temp, tmp;
        set_intersection((*it).begin(), (*it).end(), V.begin(), V.end(), inserter(temp, temp.begin()));
        if(temp.size()) {
            for(auto itt = temp.begin(); itt != temp.end(); ++itt)
                tmp.insert(M[*itt]);
            IndE.push_back(tmp);
        }
    }

    return IndE ;
}

inline HyperG BuildHyperGraph(const char * file_name, map<size_t, string> & vname) {
    /*
        the data should be organized as edge_name(element 1, element 2, ...), line by line, the detail can refer Hyperbench
    */
    ifstream fin;
    // ofstream fout;

    fin.open(file_name, ios::in);

    string str, s;
    map<string, size_t> f;
    size_t cnt = 0;
    vector <set <string> > v;

    while(getline(fin, str))
        s = s + str;

    // while(getline(fin, str)) {
    while(1) {
        size_t p_L = s.find_first_of("(");
        size_t p_R = s.find_first_of(")");
        if(p_L == str.npos || p_R == str.npos || p_R  <= p_L + 1)
            // continue; // skip this line
            break;
            
        str = s.substr(p_L + 1, p_R - p_L - 1);
        s = s.substr(p_R + 1);
        string curr = "";
        set <string> temp;
        for(size_t i = 0; i < str.size(); ++i) {
            if(str[i] == ' ' && curr == "")
                continue;
            if(str[i] == ',') {
                if(curr != "") {
                    if(!f.count(curr))
                        f[curr] = cnt++;
                    temp.insert(curr);
                    curr = "";
                }
            }
            else
                curr += str[i];
        }
        if(curr != "") {
            temp.insert(curr);
            if(!f.count(curr))
                f[curr] = cnt++;
        }
        v.push_back(temp);
    }

    // vector < set < size_t > > E = Numbered(v, f);
    // vector < vector <bool> > AdjM = GetAdjM(cnt, E);
    // vector < set <size_t> > BiComp = BiconnectComp(AdjM);

    for(auto it = f.begin(); it != f.end(); ++it) {
        // cerr << it->second << " " << it->first << endl;
        vname[it->second] = it->first;
    }
    
    return HyperG(cnt, v, f);
}

inline vector < set < size_t > > Numbered(vector < set < string > > & v, map < string, size_t > & f) {
    vector < set < size_t > > E;

    for(auto vit = v.begin(); vit != v.end(); ++vit) {
        set < size_t > e;
        for(auto it = (*vit).begin(); it != (*vit).end(); ++it)
            e.insert(f[*it]);
        E.push_back(e);
    }
    
    return E;
}

inline vector < vector <bool> > GetAdjM(size_t n, vector < set < size_t > > & E) {
    vector < vector <bool> > AdjM;
    vector <bool> temp(n, 0);
    for(size_t i = 0; i < n; ++i)
        AdjM.push_back(temp);

    for(auto Eit = E.begin(); Eit != E.end(); ++Eit) 
        for(auto it = (*Eit).begin(); it != (*Eit).end(); ++it) 
            for(auto itt = (*Eit).begin(); itt != (*Eit).end(); ++itt)
                AdjM[*it][*itt] = 1;

    return AdjM;
}

inline vector < vector <bool> > GetAdjM(size_t n, vector <VertexSet> & E) {
    vector < vector <bool> > AdjM;
    vector <bool> temp(n, 0);
    for(size_t i = 0; i < n; ++i)
        AdjM.push_back(temp);

    for(size_t i = 0; i < E.size(); ++i) {
        vector <size_t> V;
        E[i].getelement(V);

        for(size_t j = 0; j < V.size(); ++j)
            for(size_t k = j + 1; k < V.size(); ++k)
                AdjM[V[j]][V[k]] = AdjM[V[k]][V[j]] = 1;
    }   
    return AdjM;
}

inline void dfs(size_t x, size_t fa, size_t & id, 
        vector < vector <bool> > & AdjM,
        vector <size_t> & dfn, vector <size_t> & low,
        stack < pair <size_t, size_t> > & S,
        vector < set <size_t> > & BiComp
        )
{
    dfn[x] = low[x] = ++id;
    size_t n = AdjM.size();

    for(size_t i = 0; i < n; ++i) {
        if(i == x || i == fa || !AdjM[x][i])
            continue;

        if(!dfn[i]) {
            S.push( make_pair(x, i) );
            dfs(i, x, id, AdjM, dfn, low, S, BiComp);
            if(low[i] < low[x]) 
                low[x] = low[i];
            if(low[i] >= dfn[x]) {
                pair<size_t, size_t> e;
                set <size_t> Comp;
                do {
                    e = S.top(); S.pop();
                    Comp.insert(e.first);
                    Comp.insert(e.second);
                } while(e.first != x || e.second != i);
                BiComp.push_back(Comp);
            }
        }
        else if(dfn[i] < low[x]) 
            low[x] = dfn[i];
    }
}

inline vector < set <size_t> > BiconnectComp(vector < vector <bool> > & AdjM) {
    size_t n = AdjM.size(), id = 0;

    vector <size_t> dfn(n), low(n);
    stack < pair <size_t, size_t> > S;

    vector < set <size_t> > BiComp;

    dfs(0, (size_t) -1, id, AdjM, dfn, low, S, BiComp);

    return BiComp;
}

inline void MCsearch(vector < VertexSet> & G, VertexSet & curr, VertexSet cand, VertexSet & res) {
    if(curr.size() + cand.size() <= res.size())
        return ;
    if(cand.size() == 0) {
        res = curr;
        return ;
    }

    vector < size_t > V;
    cand.getelement(V);

    size_t v = V[0];

    VertexSet new_cand = cand & G[v];
    curr.Set(v);
    new_cand.reset(v);
    MCsearch(G, curr, new_cand, res);

    curr.reset(v);
    cand.reset(v);
    MCsearch(G, curr, cand, res);
}

inline VertexSet MaximalClique(HyperG & H, VertexSet & elimV) {
    vector <VertexSet> G = H.PrimalG();

    VertexSet res, curr, cand;
    for(size_t i = 0; i < H.N; ++i)
        if(!elimV.test(i)) {
            cand.Set(i);
            G[i].reset(elimV);
        }
    // printf("enter\n");
    MCsearch(G, curr, cand, res);

    return res; 
}

inline vector <HyperG> DivBiCCP(HyperG & H) {
    vector <HyperG> V;
    vector< vector <bool> > AdjM;
    vector < set <size_t> > BiComp;

    if(H.N <= 1 || H.M <= 1) {
        V.push_back(H);
        return V;
    }

    AdjM = GetAdjM(H.N, H.e);
    BiComp = BiconnectComp(AdjM);

    for(size_t i = 0; i < BiComp.size(); ++i)
        V.push_back(H.induced(BiComp[i]));
    return V;
}

inline VertexSet Aggregation(size_t v, vector <VertexSet> & E) {
    VertexSet ei;
    for(auto itt = E.begin(); itt != E.end(); ++itt)
        if((*itt).test(v))
            ei.merge(*itt);
    // ei.reset(v);
    return ei;
}




#endif