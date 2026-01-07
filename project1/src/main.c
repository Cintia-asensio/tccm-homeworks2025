/*  hw_act_mp2 — Step 2 (cumulative)
 *
 *  This program is based on the Step 1 and Step 2 examples from the repo,
 *  but everything is put together in a single file so it can be extended later.
 *
 *  What is done here:
 *   - Open a TREXIO file
 *   - Read basic information (E_nn, mo_num, nocc)
 *   - Read one-electron integrals in MO basis (core Hamiltonian h[p,q])
 *
 *  The matrix h is stored as a dense mo_num x mo_num matrix.
 *  In C it is stored in row-major order: h[p*mo_num + q].
 */

#include <stdio.h>    /* printf, fprintf */
#include <stdlib.h>   /* malloc, free, exit */
#include <trexio.h>   /* TREXIO API */

/* Simple helper function to handle TREXIO errors */
static void die_trexio(const char* msg, trexio_exit_code rc)
{
  fprintf(stderr, "%s\n%s\n", msg, trexio_string_of_error(rc));
  exit(1);
}

int main(int argc, char** argv)
{
  /* Check that the user provided the TREXIO file as argument */
  if (argc != 2) {
    fprintf(stderr, "Usage: %s file.h5\n", argv[0]);
    return 1;
  }

  /* TREXIO filename given by the user */
  const char* filename = argv[1];

  trexio_exit_code rc;

  /* Open the TREXIO file in read-only mode */
  trexio_t* f = trexio_open(filename, 'r', TREXIO_AUTO, &rc);
  if (rc != TREXIO_SUCCESS) {
    die_trexio("TREXIO open failed:", rc);
  }

  /* Read the nuclear repulsion energy (E_nn) */
  double e_nn = 0.0;
  rc = trexio_read_nucleus_repulsion(f, &e_nn);
  if (rc != TREXIO_SUCCESS) {
    die_trexio("Cannot read E_nn (nucleus_repulsion):", rc);
  }

  /* Read the number of molecular orbitals */
  int32_t mo_num = 0;
  rc = trexio_read_mo_num(f, &mo_num);
  if (rc != TREXIO_SUCCESS) {
    die_trexio("Cannot read mo_num:", rc);
  }

  /* Read number of occupied orbitals (closed-shell case) */
  int32_t nocc = 0;
  rc = trexio_read_electron_up_num(f, &nocc);
  if (rc != TREXIO_SUCCESS) {
    die_trexio("Cannot read nocc (electron_up_num):", rc);
  }

  /* Step 2: read one-electron integrals h[p,q] (core Hamiltonian) */

  /* Allocate memory for h (mo_num x mo_num matrix) */
  double* h = (double*) malloc((size_t)mo_num * (size_t)mo_num * sizeof(double));
  if (h == NULL) {
    fprintf(stderr, "Error: malloc failed for h (mo_num=%d)\n", (int)mo_num);
    exit(1);
  }

  /* Read the core Hamiltonian matrix from TREXIO */
  rc = trexio_read_mo_1e_int_core_hamiltonian(f, h);
  if (rc != TREXIO_SUCCESS) {
    die_trexio("Cannot read mo_1e_int_core_hamiltonian:", rc);
  }

  /* Print a few values just to check that it looks reasonable */
  printf("h(0,0) = %.10f\n", h[0*mo_num + 0]);
  printf("h(0,1) = %.10f\n", h[0*mo_num + 1]);
  printf("h(1,0) = %.10f\n", h[1*mo_num + 0]);
  printf("h(1,1) = %.10f\n", h[1*mo_num + 1]);

  /* Free allocated memory */
  free(h);
  h = NULL;

  /* Close the TREXIO file */
  rc = trexio_close(f);
  if (rc != TREXIO_SUCCESS) {
    die_trexio("TREXIO close failed:", rc);
  }
  f = NULL;

  /* Print some basic information */
  printf("File:   %s\n", filename);
  printf("E_nn:   %.10f\n", e_nn);
  printf("mo_num: %d\n", (int) mo_num);
  printf("nocc:   %d\n", (int) nocc);

  return 0;
}

