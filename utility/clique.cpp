//
// adapted from https://github.com/Gawssin/kCliqueListing/blob/master/DDegCol/DDegCol.c
//

#include "clique.h"

int cmp(const void* a, const void* b)
{
    iddegree *x = (iddegree*)a, *y = (iddegree*)b;
    return y->degree - x->degree;
}

int cmpadj(const void* a, const void* b)
{
    int *x = (int*)a, *y = (int*)b;
    return color[Index[*y]] - color[Index[*x]];
}

void freespecialsparse(specialsparse *g, unsigned char k) {
    unsigned char i;
    free(g->ns);
    for (i = 2; i < k + 1; i++) {
        free(g->d[i]);
        free(g->sub[i]);
    }
    free(g->d);
    free(g->sub);
    free(g->edges);
    free(g->lab);
    free(g->cd);
    free(g->adj);
    free(g);
}

void freesub(specialsparse *g, unsigned char k) {
    // keep g->edges
    unsigned char i;
    free(g->ns);
    for (i = 2; i < k + 1; i++) {
        free(g->d[i]);
        free(g->sub[i]);
    }
    free(g->d);
    free(g->sub);
    free(g->lab);
    free(g->cd);
    free(g->adj);
}

//Compute the maximum of three unsigned integers.
inline unsigned int max3(unsigned int a, unsigned int b, unsigned int c) {
    a = (a > b) ? a : b;
    return (a > c) ? a : c;
}

specialsparse* readedgelist(char* edgelist) {
    unsigned e1 = NLINKS;
    specialsparse *g = (specialsparse *)malloc(sizeof(specialsparse));
    FILE *file;

    g->n = 0;
    g->e = 0;
    file = fopen(edgelist, "r");
    g->edges = (edge *)malloc(e1 * sizeof(edge));
    while (fscanf(file, "%u %u", &(g->edges[g->e].s), &(g->edges[g->e].t)) == 2) {//Add one edge
        g->n = max3(g->n, g->edges[g->e].s, g->edges[g->e].t);
        g->e++;
        if (g->e == e1) {
            e1 += NLINKS;
            g->edges = (edge *)realloc(g->edges, e1 * sizeof(edge));
        }
    }
    fclose(file);
    g->n++;

    g->edges = (edge *)realloc(g->edges, g->e * sizeof(edge));

    return g;
}

void relabel(specialsparse *g) {
    unsigned i, source, target, tmp;

    for (i = 0; i < g->e; i++) {
        source = g->rank[g->edges[i].s];
        target = g->rank[g->edges[i].t];
        if (source < target) {
            tmp = source;
            source = target;
            target = tmp;
        }
        g->edges[i].s = source;
        g->edges[i].t = target;
    }

}

bheap *construct(unsigned n_max) {
    unsigned i;
    bheap *heap = (bheap *)malloc(sizeof(bheap));

    heap->n_max = n_max;
    heap->n = 0;
    heap->pt = (unsigned *)malloc(n_max * sizeof(unsigned));
    for (i = 0; i < n_max; i++) heap->pt[i] = -1;
    heap->kv = (keyvalue *)malloc(n_max * sizeof(keyvalue));
    return heap;
}

void swap(bheap *heap, unsigned i, unsigned j) {
    keyvalue kv_tmp = heap->kv[i];
    unsigned pt_tmp = heap->pt[kv_tmp.key];
    heap->pt[heap->kv[i].key] = heap->pt[heap->kv[j].key];
    heap->kv[i] = heap->kv[j];
    heap->pt[heap->kv[j].key] = pt_tmp;
    heap->kv[j] = kv_tmp;
}

void bubble_up(bheap *heap, unsigned i) {
    unsigned j = (i - 1) / 2;
    while (i > 0) {
        if (heap->kv[j].value > heap->kv[i].value) {
            swap(heap, i, j);
            i = j;
            j = (i - 1) / 2;
        }
        else break;
    }
}

void bubble_down(bheap *heap) {
    unsigned i = 0, j1 = 1, j2 = 2, j;
    while (j1 < heap->n) {
        j = ((j2 < heap->n) && (heap->kv[j2].value < heap->kv[j1].value)) ? j2 : j1;
        if (heap->kv[j].value < heap->kv[i].value) {
            swap(heap, i, j);
            i = j;
            j1 = 2 * i + 1;
            j2 = j1 + 1;
            continue;
        }
        break;
    }
}

void insert(bheap *heap, keyvalue kv) {
    heap->pt[kv.key] = (heap->n)++;
    heap->kv[heap->n - 1] = kv;
    bubble_up(heap, heap->n - 1);
}

void update(bheap *heap, unsigned key) {
    unsigned i = heap->pt[key];
    if (i != -1) {
        ((heap->kv[i]).value)--;
        bubble_up(heap, i);
    }
}

keyvalue popmin(bheap *heap) {
    keyvalue min = heap->kv[0];
    heap->pt[min.key] = -1;
    heap->kv[0] = heap->kv[--(heap->n)];
    heap->pt[heap->kv[0].key] = 0;
    bubble_down(heap);
    return min;
}

//Building the heap structure with (key,value)=(node,degree) for each node
bheap* mkheap(unsigned n, unsigned *v) {
    unsigned i;
    keyvalue kv;
    bheap* heap = construct(n);
    for (i = 0; i < n; i++) {
        kv.key = i;
        kv.value = v[i];
        insert(heap, kv);
    }
    return heap;
}

void freeheap(bheap *heap) {
    free(heap->pt);
    free(heap->kv);
    free(heap);
}

//computing degeneracy ordering and core value
void ord_core(specialsparse* g) {
    unsigned i, j, r = 0, n = g->n;
    keyvalue kv;
    bheap *heap;

    unsigned *d0 = (unsigned *)calloc(g->n, sizeof(unsigned));
    unsigned *cd0 = (unsigned *)malloc((g->n + 1) * sizeof(unsigned));
    unsigned *adj0 = (unsigned *)malloc(2 * g->e * sizeof(unsigned));
    for (i = 0; i < g->e; i++) {
        d0[g->edges[i].s]++;
        d0[g->edges[i].t]++;
    }
    cd0[0] = 0;
    for (i = 1; i < g->n + 1; i++) {
        cd0[i] = cd0[i - 1] + d0[i - 1];
        d0[i - 1] = 0;
    }
    for (i = 0; i < g->e; i++) {
        adj0[cd0[g->edges[i].s] + d0[g->edges[i].s]++] = g->edges[i].t;
        adj0[cd0[g->edges[i].t] + d0[g->edges[i].t]++] = g->edges[i].s;
    }

    heap = mkheap(n, d0);

    g->rank = (unsigned *)malloc(g->n * sizeof(unsigned));
    for (i = 0; i < g->n; i++) {
        kv = popmin(heap);
        g->rank[kv.key] = n - (++r);
        for (j = cd0[kv.key]; j < cd0[kv.key + 1]; j++) {
            update(heap, adj0[j]);
        }
    }
    freeheap(heap);
    free(d0);
    free(cd0);
    free(adj0);
}

//////////////////////////
//Building the special graph structure
void mkspecial(specialsparse *g, unsigned char k) {
    unsigned i, ns, max;
    unsigned *d, *sub;
    unsigned char *lab;

    d = (unsigned *)calloc(g->n, sizeof(unsigned));

    for (i = 0; i < g->e; i++) {
        d[g->edges[i].s]++;
    }

    g->cd = (unsigned *)malloc((g->n + 1) * sizeof(unsigned));
    ns = 0;
    g->cd[0] = 0;
    max = 0;
    sub = (unsigned *)malloc(g->n * sizeof(unsigned));
    lab = (unsigned char *)malloc(g->n * sizeof(unsigned char));
    for (i = 1; i < g->n + 1; i++) {
        g->cd[i] = g->cd[i - 1] + d[i - 1];
        max = (max > d[i - 1]) ? max : d[i - 1];
        sub[ns++] = i - 1;
        d[i - 1] = 0;
        lab[i - 1] = k;
    }

    subg = (specialsparse *)malloc(sizeof(specialsparse));
    subg->edges = (edge *)malloc((max*(max-1)/2) * sizeof(edge));
    subg->cd = (unsigned *)malloc((max + 1) * sizeof(unsigned));
    subg->adj = (unsigned *)malloc((max*(max - 1) / 2) * sizeof(unsigned));
    subg->ns = (unsigned *)malloc((k + 1) * sizeof(unsigned));

    subg->d = (unsigned **)malloc((k + 1) * sizeof(unsigned*));
    subg->sub = (unsigned **)malloc((k + 1) * sizeof(unsigned*));
    subg->lab = (unsigned char *)malloc(max * sizeof(unsigned char));

    for (int i = 2; i <= k; i++) {
        subg->d[i] = (unsigned *)malloc(max * sizeof(unsigned));
        subg->sub[i] = (unsigned *)malloc(max * sizeof(unsigned));
    }

    C = (int *)malloc(max * sizeof(int));
    ig = (iddegree *)malloc(max * sizeof(iddegree));
    adj0 = (unsigned  *)malloc(2 * (max*(max-1)/2) * sizeof(unsigned));
    ind = (int *)malloc(g->n * sizeof(int));
    memset(ind, -1, g->n * sizeof(int));
    loc = (int *)malloc(max * sizeof(int));
    dsub = (unsigned *)calloc(max, sizeof(unsigned));
    Index = (unsigned *)malloc(max * sizeof(unsigned));
    color = (int *)malloc(max * sizeof(int));
    cd0 = (unsigned *)malloc((max + 1) * sizeof(unsigned));

    g->adj = (unsigned *)malloc(g->e * sizeof(unsigned));

    for (i = 0; i < g->e; i++) {
        g->adj[g->cd[g->edges[i].s] + d[g->edges[i].s]++] = g->edges[i].t;
    }

    g->ns = (unsigned *)malloc((k + 1) * sizeof(unsigned));
    g->ns[k] = ns;

    g->d = (unsigned **)malloc((k + 1) * sizeof(unsigned*));
    g->sub = (unsigned **)malloc((k + 1) * sizeof(unsigned*));
    for (i = 2; i < k; i++) {
        g->d[i] = (unsigned *)malloc(g->n * sizeof(unsigned));
        g->sub[i] = (unsigned *)malloc(max * sizeof(unsigned));
    }
    g->d[k] = d;
    g->sub[k] = sub;

    g->lab = lab;
}


void mkspecial_sub(specialsparse *g, unsigned char k) {
    unsigned i, ns, max;
    unsigned *d, *sub;
    unsigned char *lab;
    for (int i = 0; i < g->n; i++)
        g->d[k][i] = 0;
    for (i = 0; i < g->e; i++) {
        g->d[k][g->edges[i].s]++;

    }

    ns = 0;
    g->cd[0] = 0;
    max = 0;

    for (i = 1; i < g->n + 1; i++) {
        g->cd[i] = g->cd[i - 1] + g->d[k][i - 1];
        max = (max > g->d[k][i - 1]) ? max : g->d[k][i - 1];
        g->sub[k][ns++] = i - 1;
        g->d[k][i - 1] = 0;
        g->lab[i - 1] = k;
    }

    for (i = 0; i < g->e; i++) {
        g->adj[g->cd[g->edges[i].s] + g->d[k][g->edges[i].s]++] = g->edges[i].t;
    }

    g->ns[k] = ns;
}

void kclique(int l, int K, specialsparse *g, unsigned *vertices, const LabelID *labels,
             std::vector<std::map<LabelID, int>> &labelFreq) {
    unsigned i, j, k, end, u, v, w;

    if (l == 2) {
        for (i = 0; i < g->ns[2]; i++) {//list all edges
            u = g->sub[2][i];
            // (*n)+=g->d[2][u];
            Count cnt = g->d[2][u];
            if (cnt == 0) continue;
            end = g->cd[u] + cnt;
            std::map<LabelID, int> freq;
            // the clique is: loc[u], vertices[j] for j \in [3, k], loc[g->adj[j]]
            VertexID v = loc[u];
            LabelID label = labels[v];
            if (freq.find(label) == freq.end()) freq[label] = 1;
            else ++freq[label];
            for (j = 3; j <= K; ++j) {
                v = vertices[j];
                label = labels[v];
                if (freq.find(label) == freq.end()) freq[label] = 1;
                else ++freq[label];
            }
            for (j = g->cd[u]; j < end; j++) {
                std::map<LabelID, int> fullFreq = freq;
                v = loc[g->adj[j]];
                label = labels[v];
                if (fullFreq.find(label) == fullFreq.end()) fullFreq[label] = 1;
                else ++fullFreq[label];
                bool exists = false;
                for (auto &item : labelFreq) {
                    if (item == fullFreq) {
                        exists = true;
                        break;
                    }
                }
                if (!exists) labelFreq.push_back(fullFreq);
            }
        }
        return;
    }

    if (l == K)
    {

        for (i = 0; i < g->ns[l]; i++) {
            u = g->sub[l][i];
            vertices[l] = u;
            g->ns[l - 1] = 0;
            end = g->cd[u] + g->d[l][u];
            for (j = g->cd[u]; j < end; j++) {//relabeling nodes and forming U'.
                v = g->adj[j];
                if (g->lab[v] == l) {		//equal to if(1)
                    g->lab[v] = l - 1;
                    g->sub[l - 1][g->ns[l - 1]++] = v;

                    g->d[l - 1][v] = 0;//new degrees
                }
            }

            if (g->ns[l - 1] < 2)
            {
                for (j = 0; j < g->ns[l - 1]; j++) {//restoring labels
                    v = g->sub[l - 1][j];
                    g->lab[v] = l;
                }
                continue;
            }


            int cnt = -1, edge_num = 0;
            for (j = 0; j < g->ns[l - 1]; j++)
            {//reodering adjacency list and computing new degrees

                v = g->sub[l - 1][j];
                if (ind[v] == -1)
                {
                    ind[v] = ++cnt;
                    loc[cnt] = v;
                    dsub[cnt] = 0;
                }
                end = g->cd[v] + g->d[l][v];
                for (k = g->cd[v]; k < end; k++)
                {
                    w = g->adj[k];
                    if (g->lab[w] == l - 1)
                    {
                        if (ind[w] == -1)
                        {
                            ind[w] = ++cnt;
                            loc[cnt] = w;
                            dsub[cnt] = 0;
                        }
                        edge_num++;
                        dsub[ind[v]]++;
                        dsub[ind[w]]++;
                    }
                }


            }

            cd0[0] = 0;
            for (int i = 1; i < g->ns[l - 1] + 1; i++) {
                cd0[i] = cd0[i - 1] + dsub[i - 1];
                ig[i - 1].id = i - 1;
                ig[i - 1].degree = dsub[i - 1];
                dsub[i - 1] = 0;
            }

            for (j = 0; j < g->ns[l - 1]; j++)
            {
                color[j] = -1;
                v = g->sub[l - 1][j];
                end = g->cd[v] + g->d[l][v];
                for (k = g->cd[v]; k < end; k++)
                {
                    w = g->adj[k];
                    if (g->lab[w] == l - 1)
                    {
                        adj0[cd0[ind[v]] + dsub[ind[v]]++] = ind[w];
                        adj0[cd0[ind[w]] + dsub[ind[w]]++] = ind[v];
                    }
                }
            }

            qsort(ig, g->ns[l - 1], sizeof(ig[0]), cmp);

            for (int i = 0; i < g->ns[l - 1]; i++)
            {
                Index[ig[i].id] = i;
                C[i] = 0;
            }

            color[0] = 0;;
            int colorNum = 0;

            for (int i = 1; i < g->ns[l - 1]; i++)
            {
                int tmpdegree = ig[i].degree, tmpid = ig[i].id;

                for (int j = 0; j < tmpdegree; j++)
                {
                    int now = Index[adj0[cd0[tmpid] + j]];
                    if (color[now] != -1)
                        C[color[now]] = 1;
                }
                for (int j = 0; j < ig[0].degree + 1; j++)
                    if (C[j] == 0)
                    {
                        color[i] = j;
                        colorNum = j > colorNum ? j : colorNum;
                        break;
                    }

                for (int j = 0; j < tmpdegree; j++)
                {
                    int now = Index[adj0[cd0[tmpid] + j]];
                    if (color[now] != -1)
                        C[color[now]] = 0;
                }

            }

            int e_num = 0;
            for (j = 0; j < g->ns[l - 1]; j++)
            {

                v = g->sub[l - 1][j];
                end = g->cd[v] + g->d[l][v];
                for (k = g->cd[v]; k < end; k++)
                {
                    w = g->adj[k];
                    if (g->lab[w] == l - 1)
                    {
                        if (color[Index[ind[v]]] < color[Index[ind[w]]])
                        {

                            subg->edges[e_num].s = ind[w];
                            subg->edges[e_num++].t = ind[v];
                        }
                        else if (color[Index[ind[v]]] == color[Index[ind[w]]])
                        {
                            if (ig[Index[ind[v]]].id < ig[Index[ind[w]]].id)
                            {
                                subg->edges[e_num].s = ind[v];
                                subg->edges[e_num++].t = ind[w];
                            }
                            else
                            {
                                subg->edges[e_num].s = ind[w];
                                subg->edges[e_num++].t = ind[v];
                            }
                        }
                        else if (color[Index[ind[v]]] > color[Index[ind[w]]])
                        {
                            subg->edges[e_num].s = ind[v];
                            subg->edges[e_num++].t = ind[w];
                        }
                    }
                }
            }

            subg->n = g->ns[l - 1];
            subg->e = edge_num;
            mkspecial_sub(subg, l - 1);

            kclique(l - 1, K, subg, vertices, labels, labelFreq);
            if (labelFreq.size() > 100) return;
            for (j = 0; j < g->ns[l - 1]; j++) {//restoring labels
                ind[loc[j]] = -1;
                v = g->sub[l - 1][j];
                g->lab[v] = l;
            }

        }
    }

    else
    {

        if (l > g->ns[l])
            return;
        for (int i = 0; i < g->ns[l]; i++) {
            u = g->sub[l][i];
            vertices[l] = loc[u];
            if (color[Index[u]] < l - 1)
                continue;
            g->ns[l - 1] = 0;
            end = g->cd[u] + g->d[l][u];
            for (j = g->cd[u]; j < end; j++) {//relabeling nodes and forming U'.
                v = g->adj[j];
                if (g->lab[v] == l) {
                    g->lab[v] = l - 1;
                    g->sub[l - 1][g->ns[l - 1]++] = v;
                    g->d[l - 1][v] = 0;//new degrees
                }
            }
            for (j = 0; j < g->ns[l - 1]; j++) {//reodering adjacency list and computing new degrees
                v = g->sub[l - 1][j];
                end = g->cd[v] + g->d[l][v];
                for (k = g->cd[v]; k < end; k++) {
                    w = g->adj[k];
                    if (g->lab[w] == l - 1) {
                        g->d[l - 1][v]++;
                    }
                    else {
                        g->adj[k--] = g->adj[--end];
                        g->adj[end] = w;
                    }

                }
            }
            kclique(l - 1, K, g, vertices, labels, labelFreq);
            if (labelFreq.size() > 100) return;
            for (j = 0; j < g->ns[l - 1]; j++) {//restoring labels
                v = g->sub[l - 1][j];
                g->lab[v] = l;
            }
        }
    }
}