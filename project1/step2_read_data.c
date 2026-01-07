#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <trexio.h>

int main(int argc, char** argv) {

  if (argc != 2) {
    fprintf(stderr, "Usage: %s <path_to_trexio_file>\n", argv[0]);
    return 1;
  }
  const char* filename = argv[1];

  trexio_exit_code rc;
  trexio_t* trexio_file = trexio_open(filename, 'r', TREXIO_AUTO, &rc);
  if (rc != TREXIO_SUCCESS) {
    printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
    exit(1);
  }

  /* 1) Read number of up-spin electrons */
  int32_t n_up = 0;
  rc = trexio_read_electron_up_num(trexio_file, &n_up);
  if (rc != TREXIO_SUCCESS) {
    printf("TREXIO Error reading electron_up_num:\n%s\n",
           trexio_string_of_error(rc));
    exit(1);
  }

  /* 2) Read number of molecular orbitals */
  int32_t mo_num = 0;
  rc = trexio_read_mo_num(trexio_file, &mo_num);
  if (rc != TREXIO_SUCCESS) {
    printf("TREXIO Error reading mo_num:\n%s\n", trexio_string_of_error(rc));
    exit(1);
  }

  /* 3) Allocate and read one-electron integrals h[p,q] */
  size_t n_h = (size_t) mo_num * (size_t) mo_num;
  double* h = (double*) malloc(n_h * sizeof(double));
  if (h == NULL) {
    fprintf(stderr, "Malloc failed for h\n");
    exit(1);
  }

  rc = trexio_read_mo_1e_int_core_hamiltonian(trexio_file, h);
  if (rc != TREXIO_SUCCESS) {
    printf("TREXIO Error reading core Hamiltonian (h_pq):\n%s\n",
           trexio_string_of_error(rc));
    exit(1);
  }

  /* Close TREXIO file */
  rc = trexio_close(trexio_file);
  if (rc != TREXIO_SUCCESS) {
    printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
    exit(1);
  }
  trexio_file = NULL;

  /* Print results */
  printf("E_NN = %.5f\n", e_nn);
  printf("n_up (=> Nocc) = %d\n", (int) n_up);
  printf("mo_num = %d\n", (int) mo_num);

  /* Show a couple of matrix elements (h is stored row-major: h[p*mo_num + q]) */
  printf("h[0,0] = %.10f\n", h[0 * mo_num + 0]);

  if (mo_num > 1) {
    printf("h[0,1] = %.10f\n", h[0 * mo_num + 1]);
    printf("h[1,0] = %.10f\n", h[1 * mo_num + 0]);
    printf("h[1,1] = %.10f\n", h[1 * mo_num + 1]);
  }

  free(h);
  h = NULL;

  return 0;
}

