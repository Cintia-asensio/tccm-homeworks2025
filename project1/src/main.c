/*  hw_act_mp2 — Step 3
 *
 *  This is Step 1 + Step 2 + Step 3:
 *   - Read the 2-electron integrals (ERI) in "sparse" format
 *
 *  In TREXIO (MO basis) the data is stored as:
 *   - mo_2e_int_eri_size : number of stored (non-zero) integrals
 *   - mo_2e_int_eri      : read a chunk of integrals (here we just read all at once)
 *
 *  The format is:
 *   - eri_idx : int32_t array of length 4*eri_size
 *       for each k:
 *         eri_idx[4*k+0] = p
 *         eri_idx[4*k+1] = q
 *         eri_idx[4*k+2] = r
 *         eri_idx[4*k+3] = s
 *   - eri_val : double array of length eri_size
 *       eri_val[k] = <pq|rs>
 */

#include <stdio.h>     /* printf, fprintf */
#include <stdlib.h>    /* malloc, free, exit */
#include <inttypes.h>  /* PRId64 for printing int64_t */
#include <trexio.h>    /* TREXIO API */

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

  /* Step 3: read ERI in sparse format */

  /* Read how many ERI entries are stored */
  int64_t eri_size = 0;
  rc = trexio_read_mo_2e_int_eri_size(f, &eri_size);
  if (rc != TREXIO_SUCCESS) {
    die_trexio("Cannot read mo_2e_int_eri_size:", rc);
  }

  printf("eri_size: %" PRId64 "\n", eri_size);

  /* If there is nothing stored, we skip reading the arrays */
  if (eri_size > 0) {

    /* Allocate buffers for indices and values */
    int32_t* eri_idx = (int32_t*) malloc((size_t)(4 * eri_size) * sizeof(int32_t));
    double*  eri_val = (double*)  malloc((size_t)eri_size * sizeof(double));

    if (eri_idx == NULL || eri_val == NULL) {
      fprintf(stderr, "Error: malloc failed for ERI buffers (eri_size=%" PRId64 ")\n", eri_size);
      exit(1);
    }

    /* Read everything in one call (TREXIO also supports chunk reading) */
    int64_t offset = 0;
    int64_t buffer_size = eri_size;

    rc = trexio_read_mo_2e_int_eri(f, offset, &buffer_size, eri_idx, eri_val);
    if (rc != TREXIO_SUCCESS) {
      die_trexio("Cannot read mo_2e_int_eri:", rc);
    }

    /* Quick check: print the first integral from the list */
    printf("ERI[0]: (%d,%d,%d,%d) = %.10f\n",
           (int)eri_idx[0], (int)eri_idx[1], (int)eri_idx[2], (int)eri_idx[3], eri_val[0]);

    /* Free buffers (later steps will store ERI in a more useful way) */
    free(eri_idx);
    free(eri_val);
    eri_idx = NULL;
    eri_val = NULL;
  }

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

  /* Free allocated memory */
  free(h);
  h = NULL;

  return 0;
}

