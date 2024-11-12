#include "readgraphmllib.h"

int readGraphMLFile(FILE* file, igraph_t* graph) {
    if (igraph_read_graph_graphml(graph, file, 0) != IGRAPH_SUCCESS) {
        fprintf(stderr, "Error reading GraphML file.\n");
        return 0;  // Return an error code or use another error handling method.
    }

    return 1;  // Successful graph read.
}

int getNumberOfConnectedVertices(const igraph_t* graph, igraph_integer_t vertex_id) {
    igraph_vector_int_t result;
    igraph_vector_int_init(&result, 0);
    igraph_neighbors(graph, &result, vertex_id, IGRAPH_ALL);
    return igraph_vector_int_size(&result);
}

void getConnectedVertices(const igraph_t* graph, igraph_integer_t vertex_id, igraph_vector_int_t* result) {
    igraph_neighbors(graph, result, vertex_id, IGRAPH_ALL);
}

void getWeightMatrix(const igraph_t* graph, double** wmatrix, size_t hasW) {
  // Get the edge attribute "weight"
  igraph_integer_t num_vertices = igraph_vcount(graph);
  if(hasW) {
    igraph_strvector_t weights;
    igraph_strvector_init(&weights, 0);

    igraph_cattribute_EASV(graph, "EdgeWeight", igraph_ess_all(IGRAPH_EDGEORDER_ID), &weights);

    for (int i = 0; i < igraph_strvector_size(&weights); i++) {
        //printf("Edge %d weight: %s\n", i, VECTOR(weights)[i]);
        igraph_integer_t from;
        igraph_integer_t to;
        igraph_edge(graph, (igraph_integer_t)i, &from, &to);
        //printf("Edge %d (%d -> %d) weight: %d\n", i, from, to, atoi(VECTOR(weights)[i]));
        wmatrix[from][to] = atof(VECTOR(weights)[i]);
        wmatrix[to][from] = atof(VECTOR(weights)[i]);
    }
    //print_matrix_d(wmatrix, wrows, wcols);

    // Free the memory
    igraph_strvector_destroy(&weights);
  }
  else {
    for (int i = 0; i < num_vertices; i++) {
      //printf("Edge %d weight: %s\n", i, VECTOR(weights)[i]);
      igraph_vector_int_t orbital_id_allowed;
      igraph_vector_int_init(&orbital_id_allowed, 0);
      getConnectedVertices(graph, (igraph_integer_t)i, &orbital_id_allowed);
      for (size_t j = 0; j < igraph_vector_int_size(&orbital_id_allowed); ++j) {
        size_t orbital_id = VECTOR(orbital_id_allowed)[j];
        igraph_integer_t from;
        igraph_integer_t to;
        from = i;
        to = orbital_id;
        wmatrix[from][to] = 1.0;
        wmatrix[to][from] = 1.0;
      }
    }
    //print_matrix_d(wmatrix, wrows, wcols);
  }
}

void getRepulsionMatrix(const igraph_t* graph, double** wmatrix, size_t hasW) {
  // Get the edge attribute "weight"
  igraph_integer_t num_vertices = igraph_vcount(graph);
  if(hasW) {
    igraph_strvector_t weights;
    igraph_strvector_init(&weights, 0);

    igraph_cattribute_EASV(graph, "EdgeRepulsion", igraph_ess_all(IGRAPH_EDGEORDER_ID), &weights);

    for (int i = 0; i < igraph_strvector_size(&weights); i++) {
        //printf("Edge %d weight: %s\n", i, VECTOR(weights)[i]);
        igraph_integer_t from;
        igraph_integer_t to;
        igraph_edge(graph, (igraph_integer_t)i, &from, &to);
        //printf("Edge %d (%d -> %d) weight: %d\n", i, from, to, atoi(VECTOR(weights)[i]));
        wmatrix[from][to] = atof(VECTOR(weights)[i]);
        wmatrix[to][from] = atof(VECTOR(weights)[i]);
    }
    //print_matrix_d(wmatrix, wrows, wcols);

    // Free the memory
    igraph_strvector_destroy(&weights);
  }
  else {
    for (int i = 0; i < num_vertices; i++) {
      //printf("Edge %d weight: %s\n", i, VECTOR(weights)[i]);
      igraph_vector_int_t orbital_id_allowed;
      igraph_vector_int_init(&orbital_id_allowed, 0);
      getConnectedVertices(graph, (igraph_integer_t)i, &orbital_id_allowed);
      for (size_t j = 0; j < igraph_vector_int_size(&orbital_id_allowed); ++j) {
        size_t orbital_id = VECTOR(orbital_id_allowed)[j];
        igraph_integer_t from;
        igraph_integer_t to;
        from = i;
        to = orbital_id;
        wmatrix[from][to] = 0.0;
        wmatrix[to][from] = 0.0;
      }
    }
    //print_matrix_d(wmatrix, wrows, wcols);
  }
}
