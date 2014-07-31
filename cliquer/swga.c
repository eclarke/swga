
#include <stdio.h>
#include <stdlib.h>
#include "cliquer.h"

int main(int argc, char *argv[]) {
  graph_t *g;
  set_t s;
  long n_sets;
  int max_weight;
  int min_weight;

  if (argc!=4) {
    fprintf(stderr,"%s <dimacs_file> <min_weight> <max_weight>\n",argv[0]);
    return 1;
  }
  g=graph_read_dimacs_file(argv[1]);
  if (g==NULL)
    return 1;

  ASSERT(graph_test(g,stderr));

  min_weight = atoi(argv[2]);
  max_weight = atoi(argv[3]);
  fprintf(stdout, "Min weight %i, max weight %i\n", min_weight, max_weight);
  clique_default_options->output=stdout;
  clique_default_options->
  clique_unweighted_find_all(g, min_weight, max_weight, FALSE, clique_default_options);


  return 0;
}
