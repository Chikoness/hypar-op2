'''
Python script to create plots from the solution of a
HyPar simulation.

- if op_overwrite is set to "no", a plot is generated
for each variable (solution vector component) and each
simulation time for which the solution is available.
- if op_overwrite is set to "yes", a single plot is
created for for each variable (solution vector component).
- solution must be in binary format

Make sure the environment variable "HYPAR_DIR" is set
and points to the correct location (/path/to/hypar)
'''

import os
import time
import sys

import numpy as np
import numpy.linalg as LA
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

hypar_dir = os.environ.get('HYPAR_DIR')
hypar_dir_python = hypar_dir + '/Examples/Python'
sys.path.append(hypar_dir_python)

import modHyParUtils as hyparutils

font = {'size':22}
matplotlib.rc('font', **font)
colormap='Spectral'

fig_nv=1
fig_nh=3
figsize=(1.4*fig_nh*10,fig_nv*10)
plt_dir_name='plots'

'''
Set up the simulation parameters
'''
sim_path = '.'
sim_inp_data = hyparutils.readHyParInpFile(sim_path+'/simulation.inp')
solver_inp_data = hyparutils.readHyParInpFile(sim_path+'/solver.inp')

if not sim_inp_data:
    nsims = 1
else:
    nsims = int(sim_inp_data['nsims'][0])

print('Number of simulations: ', nsims)

if solver_inp_data['op_file_format'] != 'binary':
    raise Exception("op_file_format must be 'binary' in solver.inp.")

ndims = int(solver_inp_data['ndims'][0])
nvars = int(solver_inp_data['nvars'][0])
size = np.int32(solver_inp_data['size']).reshape(nsims,ndims)

if not os.path.exists(plt_dir_name):
      os.makedirs(plt_dir_name)

if solver_inp_data['op_overwrite'] == 'no':

  niter = int(solver_inp_data['n_iter'][0])
  dt = float(solver_inp_data['dt'][0])
  t_final = dt*niter

  op_write_iter = int(solver_inp_data['file_op_iter'][0])
  dt_snapshots = op_write_iter*dt
  n_snapshots = int(niter/op_write_iter) + 1

  print('Simulation parameters:')
  print('  ndims = ', ndims)
  print('  nvars = ', nvars)
  print('  grid size = ', size)
  print('  simulation dt = ', dt)
  print('  niter = ', niter)
  print('  final time = ', t_final)
  print('  snapshot dt = ', dt_snapshots)
  print('  number of snapshots = ', n_snapshots)

  for s in range(nsims):

    if nsims > 1:
      size_sim = size[s,:]
    else:
      size_sim = size

    '''
    Load simulation data (solution snapshots)
    '''
    grid, sol_fom_snapshots = hyparutils.getSimulationSnapshots(  sim_path,
                                                                  nsims,
                                                                  s,
                                                                  n_snapshots,
                                                                  ndims,
                                                                  nvars,
                                                                  size_sim )
    sol_fom_snapshots = np.float32(sol_fom_snapshots)

    grid, sol_rom_snapshots = hyparutils.getSimulationSnapshots(  sim_path,
                                                                  nsims,
                                                                  s,
                                                                  n_snapshots,
                                                                  ndims,
                                                                  nvars,
                                                                  size_sim,
                                                                  op_root='op_rom')
    sol_rom_snapshots = np.float32(sol_rom_snapshots)

    x = grid[:size_sim[0]]
    y = grid[size_sim[0]:]
    print('2D Domain:');
    print(' x: ', np.min(x), np.max(x))
    print(' x.shape: ', x.shape)
    print(' y: ', np.min(y), np.max(y))
    print(' y.shape: ', y.shape)
    y2d, x2d = np.meshgrid(y, x)

    for var in range(nvars):
      for i in range(n_snapshots):
        fig, axes = plt.subplots(fig_nv,fig_nh,figsize=figsize)
        sol_fom_2d = np.transpose(sol_fom_snapshots.reshape(n_snapshots,size_sim[1],size_sim[0],nvars))[var,:,:,i]
        sol_rom_2d = np.transpose(sol_rom_snapshots.reshape(n_snapshots,size_sim[1],size_sim[0],nvars))[var,:,:,i]
        sol_diff_2d = sol_fom_2d - sol_rom_2d
        fom_sol_norm = LA.norm(sol_fom_2d)
        diff_norm = LA.norm(sol_diff_2d)
        if fom_sol_norm > 1e-14:
          diff_norm = diff_norm / fom_sol_norm
        # Plot solutions
        plot_fom = axes[0].pcolor(x2d, y2d, sol_fom_2d, cmap=colormap)
        plot_rom = axes[1].pcolor(x2d, y2d, sol_rom_2d, cmap=colormap)
        plot_diff = axes[2].pcolor(x2d, y2d, sol_diff_2d, cmap=colormap)
        axes[0].set_title('FOM: var {:}, t={:.3}'.format(var,i*dt_snapshots))
        axes[1].set_title('ROM: var {:}, t={:.3}'.format(var,i*dt_snapshots))
        axes[2].set_title('Diff: var {:}, t={:.3}'.format(var,i*dt_snapshots))
        fig.colorbar(plot_fom, ax=axes[0])
        fig.colorbar(plot_rom, ax=axes[1])
        fig.colorbar(plot_diff, ax=axes[2])
        if nsims > 1:
          plt_fname = plt_dir_name+'/fig_'+f'{s:02d}'+'_'+f'{var:02d}'+'_'+f'{i:05d}'+'.png'
        else:
          plt_fname = plt_dir_name+'/fig_'+f'{var:02d}'+'_'+f'{i:05d}'+'.png'
        print('Saving %s' % plt_fname)
        plt.savefig(plt_fname)
        plt.close()

else:

  niter = int(solver_inp_data['n_iter'][0])
  dt = float(solver_inp_data['dt'][0])
  t_final = dt*niter

  n_snapshots = 1

  print('Simulation parameters:')
  print('  ndims = ', ndims)
  print('  nvars = ', nvars)
  print('  grid size = ', size)
  print('  final time = ', t_final)
  print('  number of snapshots = ', n_snapshots)

  for s in range(nsims):

    if nsims > 1:
      size_sim = size[s,:]
    else:
      size_sim = size

    '''
    Load simulation data (solution snapshots)
    '''
    grid,sol_fom_snapshots = hyparutils.getSimulationSnapshots( sim_path,
                                                                nsims,
                                                                s,
                                                                n_snapshots,
                                                                ndims,
                                                                nvars,
                                                                size_sim )
    grid,sol_rom_snapshots = hyparutils.getSimulationSnapshots( sim_path,
                                                                nsims,
                                                                s,
                                                                n_snapshots,
                                                                ndims,
                                                                nvars,
                                                                size_sim,
                                                                op_root='op_rom')
    sol_rom_snapshots = np.float32(sol_rom_snapshots)
    sol_fom_snapshots = np.float32(sol_fom_snapshots)

    x = grid[:size_sim[0]]
    y = grid[size_sim[0]:]
    print('2D Domain:');
    print(' x: ', np.min(x), np.max(x))
    print(' x.shape: ', x.shape)
    print(' y: ', np.min(y), np.max(y))
    print(' y.shape: ', y.shape)
    y2d, x2d = np.meshgrid(y, x)

    for var in range(nvars):
      fig, axes = plt.subplots(fig_nv,fig_nh,figsize=figsize)
      sol_fom_2d = np.transpose(sol_fom_snapshots.reshape(size_sim[1],size_sim[0],nvars))[var,:,:]
      sol_rom_2d = np.transpose(sol_rom_snapshots.reshape(size_sim[1],size_sim[0],nvars))[var,:,:]
      sol_diff_2d = sol_fom_2d - sol_rom_2d
      fom_sol_norm = LA.norm(sol_fom_2d)
      diff_norm = LA.norm(sol_diff_2d)
      if fom_sol_norm > 1e-14:
        diff_norm = diff_norm / fom_sol_norm
      # Plot solutions
      plot_fom = axes[0].pcolor(x2d, y2d, sol_fom_2d, cmap=colormap)
      plot_rom = axes[1].pcolor(x2d, y2d, sol_rom_2d, cmap=colormap)
      plot_diff = axes[2].pcolor(x2d, y2d, sol_diff_2d, cmap=colormap)
      axes[0].set_title('FOM: var {:}, t={:.3}'.format(var,t_final))
      axes[1].set_title('ROM: var {:}, t={:.3}'.format(var,t_final))
      axes[2].set_title('Diff: var {:}, t={:.3}'.format(var,t_final))
      fig.colorbar(plot_fom, ax=axes[0])
      fig.colorbar(plot_rom, ax=axes[1])
      fig.colorbar(plot_diff, ax=axes[2])
      if nsims > 1:
        plt_fname = plt_dir_name+'/fig_'+f'{s:02d}'+'_'+f'{var:02d}'+'.png'
      else:
          plt_fname = plt_dir_name+'/fig_'+f'{var:02d}'+'.png'
      print('Saving %s' % plt_fname)
      plt.savefig(plt_fname)
      plt.close()
