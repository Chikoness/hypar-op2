/*! @file FirstDerivativeSecondOrderNoGhosts.c
    @author Ping-Hsuan Tsai 
    @brief Second order finite-difference approximation to first derivative
*/

#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <firstderivative.h>

#ifdef with_omp
#include <omp.h>
#endif

/*! Computes the second-order finite-difference approximation to the first derivative 
    (\b Note: not divided by the grid spacing):
    \f{equation}{
      \left(\partial f\right)_i = \left\{ \begin{array}{ll} \frac{1}{2}\left(-3f_i+4f_{i+1}-f_{i+2}\right) & i = -g \\ \frac{1}{2}\left( f_{i+1} - f_{i-1} \right) & -g+1 \leq i \leq N+g-2 \\ \frac{1}{2}\left( f_{i-2} -4f_{i-1}+3f_i \right) & i = N+g-1 \end{array}\right.
    \f}
    where \f$i\f$ is the grid index along the spatial dimension of the derivative, \f$g\f$ is the number of ghost points, and \f$N\f$ is the number of grid points (not including the ghost points) in the spatial dimension of the derivative.
    \n\n
    Notes:
    + The first derivative is computed at the grid points or the cell centers.
    + \b Df and \b f are 1D arrays containing the function and its computed derivatives on a multi-
      dimensional grid. The derivative along the specified dimension \b dir is computed by looping
      through all grid lines along \b dir.
    + \b Df and \b f must have the same number of ghost points (including the case with 0 ghost points).
    + Biased stencils are used at the boundary and ghost point values are not used. So this function can be used
      when the ghost values are not filled.
*/
int FirstDerivativeSecondOrderCentralNoGhosts(
                                                double  *Df,    /*!< Array to hold the computed first derivative (with ghost points) */
                                                double  *f,     /*!< Array containing the grid point function values whose first 
                                                                     derivative is to be computed (with ghost points) */
                                                int     dir,    /*!< The spatial dimension along which the derivative is computed */
                                                int     bias,   /*!< Forward or backward differencing for non-central 
                                                                     finite-difference schemes (-1: backward, 1: forward)*/
                                                int     ndims,  /*!< Number of spatial/coordinate dimensions */
                                                int     *dim,   /*!< Local dimensions */
                                                int     ghosts, /*!< Number of ghost points */
                                                int     nvars   /*!< Number of vector components at each grid points */ )
{
  int  i, j, v;

  if ((!Df) || (!f)) {
    fprintf(stderr, "Error in FirstDerivativeSecondOrder(): input arrays not allocated.\n");
    return(1);
  }

  /* create index and bounds for the outer loop, i.e., to loop over all 1D lines along
     dimension "dir"                                                                    */
  int indexC[ndims], index_outer[ndims], bounds_outer[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  int N_outer;  _ArrayProduct1D_(bounds_outer,ndims,N_outer);

#pragma omp parallel for schedule(auto) default(shared) private(i,j,v,index_outer,indexC)
  for (j=0; j<N_outer; j++) {
    _ArrayIndexnD_(ndims,j,bounds_outer,index_outer,0);
    _ArrayCopy1D_(index_outer,indexC,ndims);
    /* left boundary */
    {
      i = 0;
      int qC, qR, qR2;
      indexC[dir] = i  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qC );
      indexC[dir] = i+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qR);
      indexC[dir] = i+2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qR2);
      for (v=0; v<nvars; v++)  Df[qC*nvars+v] = (-0.5*f[qR2*nvars+v]+2*f[qR*nvars+v]-1.5*f[qC*nvars+v]);
    }
    /* interior */
    for (i = 1; i < dim[dir]-1; i++) {
      int qC, qL, qR;
      indexC[dir] = i-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qL);
      indexC[dir] = i  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qC );
      indexC[dir] = i+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qR);
      for (v=0; v<nvars; v++)  Df[qC*nvars+v] = 0.5 * (f[qR*nvars+v]-f[qL*nvars+v]);
    }
    /* right boundary */
    {
      i = dim[dir] - 1;
      int qC, qL, qL2;
      indexC[dir] = i-2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qL2);
      indexC[dir] = i-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qL);
      indexC[dir] = i  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qC );
      for (v=0; v<nvars; v++)  Df[qC*nvars+v] = (0.5*f[qL2*nvars+v]-2*f[qL*nvars+v]+1.5*f[qC*nvars+v]);
    }
  }
  
  return 0;
}
