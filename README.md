# ASDMatch

## Compile

1. Compile and link to the glpk library.  The glpk library is used to compute fractional hypertree decomposition.

```shell
cd utility/td/glpk-4.55
./configure -prefix={ABSOLUTE_PATH_TO_GLPK_LIB} # relative path: ../glpk_lib
make -j
make install
cd .. # at the 'td' directory

```

2. Build the project.

```shell
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

## Input format

We use the same format as in [RapidMatch](https://github.com/RapidsAtHKUST/RapidMatch).

Both the input query graph and data graph are vertex-labeled. Each graph starts with 't N M' where N is the number of vertices and M is the number of edges. A vertex and an edge are formatted as 'v VertexID LabelId Degree' and 'e VertexId VertexId' respectively. Note that we require that the vertex id is started from 0 and the range is [0,N - 1] where V is the vertex set. The following is an input sample. 

Example:

```
t 5 6
v 0 0 2
v 1 1 3
v 2 2 3
v 3 1 2
v 4 2 2
e 0 1
e 0 2
e 1 2
e 1 3
e 2 4
e 3 4
```

## Execution and output

### ASDMatch.out:

| Option | Description                                          |
| ------ | ---------------------------------------------------- |
| -q     | the query graph path                                 |
| -d     | the data graph path                                  |
| -s     | the cardinalities used in query optimizer (optional) |
| -m     | the memory budget by KB. The default value is 16GB   |
| -r     | the result path (optional)                           |

Example:

```
./ASDMatch.out -q ../dataset/yeast/query_graph/query_8_1.graph -d ../dataset/yeast/data_graph/yeast.graph -s ../dataset/yeast/query_graph/query_8_1.card -r ../result.txt
```

We can also use 'prepare.out' to precompute cardinalities.
```shell
./prepare.out -q {query_graph_path} -d {data_graph_path} -s {card_result_path}
```
