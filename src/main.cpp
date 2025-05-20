#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <string>

#include <mpivars_cpp.h>
#include <simulation_library.h>

// Declare OP2
#include "op_seq.h"

static const char help[] = "HyPar - A finite-difference algorithm for solving hyperbolic-parabolic PDEs";

/*!
 * \brief Main driver
 *
 * The main driver function that calls the initialization, solving, and cleaning up functions.
 */
int main(int argc, char **argv)
{
  int ierr = 0, d, n;
  struct timeval main_start, solve_start;
  struct timeval main_end, solve_end;

  MPI_Comm world;
  int rank, nproc;
  MPI_Init(&argc, &argv);
  MPI_Comm_dup(MPI_COMM_WORLD, &world);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  op_init(argc, argv, 1); // initialize OP2

  if (!rank)
    printf("HyPar - Parallel (MPI) version with %d processes\n", nproc);

#ifdef with_petsc
  PetscInitialize(&argc, &argv, (char *)0, help);
  if (!rank)
    printf("Compiled with PETSc time integration.\n");
#endif

#ifdef with_python
  initializePython(rank);
  initializePythonPlotting(rank);
#endif

  op_printf("Hello from OP2!\n");

  gettimeofday(&main_start, NULL);

  int sim_type = -1;
  Simulation *sim = NULL;

  if (!rank)
  {

    std::string ensemble_sim_fname(_ENSEMBLE_SIM_INP_FNAME_);
    std::string sparsegrids_sim_fname(_SPARSEGRIDS_SIM_INP_FNAME_);

    FILE *f_ensemble_sim = fopen(ensemble_sim_fname.c_str(), "r");
    FILE *f_sparsegrids_sim = fopen(sparsegrids_sim_fname.c_str(), "r");

    if (f_ensemble_sim && f_sparsegrids_sim)
    {

      fprintf(stderr, "Error: Cannot have both %s and %s input files.\n",
              _ENSEMBLE_SIM_INP_FNAME_, _SPARSEGRIDS_SIM_INP_FNAME_);
      fprintf(stderr, "Remove one or both of them depending on the kind of simulation you want to run.\n");
      fclose(f_ensemble_sim);
      fclose(f_sparsegrids_sim);
    }
    else if (f_ensemble_sim)
    {

      sim_type = _SIM_TYPE_ENSEMBLE_;
      fclose(f_ensemble_sim);
    }
    else if (f_sparsegrids_sim)
    {

      sim_type = _SIM_TYPE_SPARSE_GRIDS_;
      fclose(f_sparsegrids_sim);
    }
    else
    {

      sim_type = _SIM_TYPE_SINGLE_;
    }
  }

  if (sim_type == _SIM_TYPE_SINGLE_)
  {
    sim = new SingleSimulation;
    if (!rank)
      printf("-- Single Simulation --\n");
    sim = new EnsembleSimulation;
  }
  else if (sim_type == _SIM_TYPE_ENSEMBLE_)
  {
    if (!rank)
      printf("-- Ensemble Simulation --\n");
    sim = new EnsembleSimulation;
  }
  else if (sim_type == _SIM_TYPE_SPARSE_GRIDS_)
  {
    if (!rank)
      printf("-- Sparse Grids Simulation --\n");
    sim = new SparseGridsSimulation;
  }
  else
  {
    fprintf(stderr, "ERROR: invalid sim_type (%d) on rank %d.\n",
            sim_type, rank);
  }

  if (sim == NULL)
  {
    fprintf(stderr, "ERROR: unable to create sim on rank %d.\n",
            rank);
    return 1;
  }

  /* Allocate simulation objects */
  ierr = sim->define(rank, nproc);
  if (!sim->isDefined())
  {
    printf("Error: Simulation::define() failed on rank %d\n",
           rank);
    return 1;
  }
  if (ierr)
  {
    printf("Error: Simulation::define() returned with status %d on process %d.\n",
           ierr, rank);
    return (ierr);
  }

  /* Read Inputs */
  ierr = sim->ReadInputs();
  if (ierr)
  {
    printf("Error: Simulation::ReadInputs() returned with status %d on process %d.\n", ierr, rank);
    return (ierr);
  }

  /* Initialize and allocate arrays */
  ierr = sim->Initialize();
  if (ierr)
  {
    printf("Error: Simulation::Initialize() returned with status %d on process %d.\n", ierr, rank);
    return (ierr);
  }

  /* read and set grid & initial solution */
  ierr = sim->InitialSolution();
  if (ierr)
  {
    printf("Error: Simulation::InitialSolution() returned with status %d on process %d.\n", ierr, rank);
    return (ierr);
  }

  /* Initialize domain boundaries */
  ierr = sim->InitializeBoundaries();
  if (ierr)
  {
    printf("Error: Simulation::InitializeBoundaries() returned with status %d on process %d.\n", ierr, rank);
    return (ierr);
  }

  /* Initialize immersed boundaries */
  ierr = sim->InitializeImmersedBoundaries();
  if (ierr)
  {
    printf("Error: Simulation::InitializeImmersedBoundaries() returned with status %d on process %d.\n", ierr, rank);
    return (ierr);
  }

  /* Initialize solvers */
  ierr = sim->InitializeSolvers();
  if (ierr)
  {
    printf("Error: Simulation::InitializeSolvers() returned with status %d on process %d.\n", ierr, rank);
    return (ierr);
  }

  /* Initialize physics */
  ierr = sim->InitializePhysics();
  if (ierr)
  {
    printf("Error: Simulation::InitializePhysics() returned with status %d on process %d.\n", ierr, rank);
    return (ierr);
  }

  /* Initialize physics data */
  ierr = sim->InitializePhysicsData();
  if (ierr)
  {
    printf("Error: Simulation::InitializePhysicsData() returned with status %d on process %d.\n", ierr, rank);
    return (ierr);
  }

  /* Wrap up initializations */
  ierr = sim->InitializationWrapup();
  if (ierr)
  {
    printf("Error: Simulation::InitializationWrapup() returned with status %d on process %d.\n", ierr, rank);
    return (ierr);
  }

  /* Initializations complete */

  /* Run the solver */
  gettimeofday(&solve_start, NULL);

  /* Use native time-integration */
  ierr = sim->Solve();
  if (ierr)
  {
    printf("Error: Simulation::Solve() returned with status %d on process %d.\n", ierr, rank);
    return (ierr);
  }

  gettimeofday(&solve_end, NULL);
  gettimeofday(&main_end, NULL);

  /* calculate solver and total runtimes */
  long long walltime;
  walltime = ((main_end.tv_sec * 1000000 + main_end.tv_usec) - (main_start.tv_sec * 1000000 + main_start.tv_usec));
  double main_runtime = (double)walltime / 1000000.0;
  ierr = MPIMax_double(&main_runtime, &main_runtime, 1, &world);
  if (ierr)
    return (ierr);
  walltime = ((solve_end.tv_sec * 1000000 + solve_end.tv_usec) - (solve_start.tv_sec * 1000000 + solve_start.tv_usec));
  double solver_runtime = (double)walltime / 1000000.0;
  ierr = MPIMax_double(&solver_runtime, &solver_runtime, 1, &world);
  if (ierr)
    return (ierr);

  /* Write errors and other data */
  sim->WriteErrors(solver_runtime, main_runtime);

  /* Cleaning up */
  delete sim;
  if (!rank)
    printf("Finished.\n");

  MPI_Comm_free(&world);
  MPI_Finalize();

  op_exit(); // Finalise OP2

  return (0);
}