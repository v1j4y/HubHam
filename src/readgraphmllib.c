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

int print_attributes(const igraph_t *g) {

    igraph_vector_t gtypes, vtypes, etypes;
    igraph_strvector_t gnames, vnames, enames;
    long int i;

    igraph_vector_t vec;
    igraph_strvector_t svec;
    long int j;

    igraph_vector_init(&gtypes, 0);
    igraph_vector_init(&vtypes, 0);
    igraph_vector_init(&etypes, 0);
    igraph_strvector_init(&gnames, 0);
    igraph_strvector_init(&vnames, 0);
    igraph_strvector_init(&enames, 0);

    igraph_cattribute_list(g, &gnames, &gtypes, &vnames, &vtypes,
                           &enames, &etypes);

    /* Graph attributes */
    for (i = 0; i < igraph_strvector_size(&gnames); i++) {
        printf("%s=", STR(gnames, i));
        if (VECTOR(gtypes)[i] == IGRAPH_ATTRIBUTE_NUMERIC) {
            igraph_real_printf(GAN(g, STR(gnames, i)));
            putchar(' ');
        } else {
            printf("\"%s\" ", GAS(g, STR(gnames, i)));
        }
    }
    printf("\n");

    for (i = 0; i < igraph_vcount(g); i++) {
        long int j;
        printf("Vertex %li: ", i);
        for (j = 0; j < igraph_strvector_size(&vnames); j++) {
            printf("%s=", STR(vnames, j));
            if (VECTOR(vtypes)[j] == IGRAPH_ATTRIBUTE_NUMERIC) {
                igraph_real_printf(VAN(g, STR(vnames, j), i));
                putchar(' ');
            } else {
                printf("\"%s\" ", VAS(g, STR(vnames, j), i));
            }
        }
        printf("\n");
    }

    for (i = 0; i < igraph_ecount(g); i++) {
        long int j;
        printf("Edge %li (%i-%i): ", i, (int)IGRAPH_FROM(g, i), (int)IGRAPH_TO(g, i));
        for (j = 0; j < igraph_strvector_size(&enames); j++) {
            printf("%s=", STR(enames, j));
            if (VECTOR(etypes)[j] == IGRAPH_ATTRIBUTE_NUMERIC) {
                igraph_real_printf(EAN(g, STR(enames, j), i));
                putchar(' ');
            } else {
                printf("\"%s\" ", EAS(g, STR(enames, j), i));
            }
        }
        printf("\n");
    }

    /* Check vector-based query functions */
    igraph_vector_init(&vec, 0);
    igraph_strvector_init(&svec, 0);

    for (j = 0; j < igraph_strvector_size(&vnames); j++) {
        if (VECTOR(vtypes)[j] == IGRAPH_ATTRIBUTE_NUMERIC) {
            igraph_cattribute_VANV(g, STR(vnames, j), igraph_vss_all(), &vec);
            for (i = 0; i < igraph_vcount(g); i++) {
                igraph_real_t num = VAN(g, STR(vnames, j), i);
                if (num != VECTOR(vec)[i] &&
                    (!isnan(num) || !isnan(VECTOR(vec)[i]))) {
                    exit(51);
                }
            }
        } else {
            igraph_cattribute_VASV(g, STR(vnames, j), igraph_vss_all(), &svec);
            for (i = 0; i < igraph_vcount(g); i++) {
                const char *str = VAS(g, STR(vnames, j), i);
                if (strcmp(str, STR(svec, i))) {
                    exit(52);
                }
            }
        }
    }

    for (j = 0; j < igraph_strvector_size(&enames); j++) {
        if (VECTOR(etypes)[j] == IGRAPH_ATTRIBUTE_NUMERIC) {
            igraph_cattribute_EANV(g, STR(enames, j),
                                   igraph_ess_all(IGRAPH_EDGEORDER_ID), &vec);
            for (i = 0; i < igraph_ecount(g); i++) {
                igraph_real_t num = EAN(g, STR(enames, j), i);
                if (num != VECTOR(vec)[i] &&
                    (!isnan(num) || !isnan(VECTOR(vec)[i]))) {
                    exit(53);
                }
            }
        } else {
            igraph_cattribute_EASV(g, STR(enames, j),
                                   igraph_ess_all(IGRAPH_EDGEORDER_ID), &svec);
            for (i = 0; i < igraph_ecount(g); i++) {
                const char *str = EAS(g, STR(enames, j), i);
                if (strcmp(str, STR(svec, i))) {
                    exit(54);
                }
            }
        }
    }

    igraph_strvector_destroy(&svec);
    igraph_vector_destroy(&vec);

    igraph_strvector_destroy(&enames);
    igraph_strvector_destroy(&vnames);
    igraph_strvector_destroy(&gnames);
    igraph_vector_destroy(&etypes);
    igraph_vector_destroy(&vtypes);
    igraph_vector_destroy(&gtypes);

    return 0;
}
