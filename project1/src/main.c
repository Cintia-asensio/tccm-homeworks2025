/*  hw_act_mp2 — Step 4
 *
 *  This is Step 3 accumulative + Step 4:
 *   - Convert the ERI list (sparse) into a 4D array eri[p,q,r,s]
 *
 *  We do this because later (for HF / MP2) it is easier to access integrals with direct indexing.
 *
 *  In this step we:
 *   - Read the ERI in sparse format (eri_idx + eri_val)
 *   - Allocate a dense array of size mo_num^4 (initialized to 0)
 *   - Fill it using the 8-fold symmetry of (pq|rs)
 */

#include <stdio.h>     /* printf, fprintf */
#include <stdlib.h>    /* malloc, calloc, free, exit */
#include <inttypes.h>  /* PRId64 for printing int64_t */
#include <trexio.h>    /* TREXIO API */

/* Simple helper function to handle TREXIO errors */
static void die_trexio(const char* msg, trexio_exit_code rc)
{
  fprintf(stderr, "%s\n%s\n", msg, trexio_string_of_error(rc));
  exit(1);
}

/* Step 4 helper: convert (p,q,r,s) to a single index for a 1D array */
static size_t idx4(int p, int q, int r, int s, int mo_num)
{
  return (((size_t)p * (size_t)mo_num + (size_t)q) * (size_t)mo_num + (size_t)r) * (size_t)mo_num + (size_t)s;
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

  /* Step 4: build a dense ERI array (initialized to 0.0) */
  double* eri = (double*) calloc((size_t)mo_num * (size_t)mo_num * (size_t)mo_num * (size_t)mo_num, sizeof(double));
  if (eri == NULL) {
    fprintf(stderr, "Error: calloc failed for dense ERI (mo_num=%d)\n", (int)mo_num);
    exit(1);
  }

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

    /* Fill the dense tensor using symmetry.
       Here we just copy each (p,q,r,s) into the 8 equivalent positions. */
    for (int64_t k = 0; k < buffer_size; ++k) {
      int p = eri_idx[4*k + 0];
      int q = eri_idx[4*k + 1];
      int r = eri_idx[4*k + 2];
      int s = eri_idx[4*k + 3];
      double v = eri_val[k];

      /* 8-fold symmetry */
      eri[idx4(p,q,r,s,mo_num)] = v;
      eri[idx4(r,q,p,s,mo_num)] = v;
      eri[idx4(p,s,r,q,mo_num)] = v;
      eri[idx4(r,s,p,q,mo_num)] = v;

      eri[idx4(q,p,s,r,mo_num)] = v;
      eri[idx4(s,p,q,r,mo_num)] = v;
      eri[idx4(q,r,s,p,mo_num)] = v;
      eri[idx4(s,r,q,p,mo_num)] = v;
    }

    /* Quick check using two equivalent permutations */
    printf("ERI dense check: (0,0,1,0)=%.10f  (1,0,0,0)=%.10f\n",
           eri[idx4(0,0,1,0,mo_num)], eri[idx4(1,0,0,0,mo_num)]);

    /* Free buffers (we keep only the dense eri array for the next steps) */
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

  free(eri);
  eri = NULL;

  return 0;
}

