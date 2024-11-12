#include <igraph.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h> /* unlink */

int readGraphMLFile(FILE* file, igraph_t* graph) ;

int getNumberOfConnectedVertices(const igraph_t* graph, igraph_integer_t vertex_id) ;

void getConnectedVertices(const igraph_t* graph, igraph_integer_t vertex_id, igraph_vector_int_t* result) ;

void getWeightMatrix(const igraph_t* graph, double** wmatrix, size_t hasW) ;

void getRepulsionMatrix(const igraph_t* graph, double** wmatrix, size_t hasW) ;
