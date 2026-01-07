#include <stdio.h>
#include <stdlib.h>
#include <trexio.h>

int main(int argc, char** argv) {

  /* Check that the user provided the TREXIO file as argument */
  if (argc != 2) {
    fprintf(stderr, "Usage: %s <path_to_trexio_file>\n", argv[0]);
    return 1;
  }

  const char* filename = argv[1];

  trexio_exit_code rc;

  /* Open the TREXIO file in read-only mode */
  trexio_t* trexio_file = trexio_open(filename, 'r', TREXIO_AUTO, &rc);
  if (rc != TREXIO_SUCCESS) {
    printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
    exit(1);
  }

  /* Read the nuclear repulsion energy */
  double energy;
  rc = trexio_read_nucleus_repulsion(trexio_file, &energy);
  if (rc != TREXIO_SUCCESS) {
    printf("TREXIO Error reading nuclear repulsion energy:\n%s\n",
           trexio_string_of_error(rc));
    exit(1);
  }

  /* Close the TREXIO file */
  rc = trexio_close(trexio_file);
  if (rc != TREXIO_SUCCESS) {
    printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
    exit(1);
  }
  trexio_file = NULL;

  /* Print the result */
  printf("Nuclear repulsion energy (E_NN) = %.5f\n", energy);

  return 0;
}


