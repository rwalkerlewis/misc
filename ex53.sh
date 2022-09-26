#!/bin/bash

echo "Problem: $1"
test=$1

if [ "$test" = "2dql" ]; then
  # 2D Quad Linear
  ./ex53  -sol_type quadratic_linear -dm_refine 2 \
      -displacement_petscspace_degree 2 -tracestrain_petscspace_degree 1 -pressure_petscspace_degree 1 \
      -dmts_check .0001 -ts_max_steps 5 -ts_monitor_extreme -dm_view hdf5:./output/sol_2dql.h5 -monitor_solution \
        hdf5:./output/sol_2dql.h5::append

elif [ "$test" = "3dql" ]; then
  # 3D Quad Linear
  ./ex53  -dm_plex_dim 3 -sol_type quadratic_linear -dm_refine 1 \
      -displacement_petscspace_degree 2 -tracestrain_petscspace_degree 1 -pressure_petscspace_degree 1 \
      -dmts_check .0001 -ts_max_steps 5 -ts_monitor_extreme -dm_view hdf5:./output/sol_3dql.h5 -monitor_solution \
        hdf5:./output/sol_3dql.h5::append

if [ "$test" = "2dtl" ]; then
  # 2D Trig Linear
  ./ex53  -sol_type trig_linear -dm_refine 2 \
      -displacement_petscspace_degree 2 -tracestrain_petscspace_degree 1 -pressure_petscspace_degree 1 \
      -dmts_check .0001 -ts_max_steps 5 -ts_monitor_extreme -dm_view hdf5:./output/sol_2dtl.h5 -monitor_solution \
        hdf5:./output/sol_2dtl.h5::append

elif [ "$test" = "3dtl" ]; then
  # 3D Trig Linear
  ./ex53  -dm_plex_dim 3 -sol_type trig_linear -dm_refine 1 \
      -displacement_petscspace_degree 2 -tracestrain_petscspace_degree 1 -pressure_petscspace_degree 1 \
      -dmts_check .0001 -ts_max_steps 5 -ts_monitor_extreme -dm_view hdf5:./output/sol_3dql.h5 -monitor_solution \
        hdf5:./output/sol_3dtl.h5::append

if [ "$test" = "2dqt" ]; then
  # 2D Quad Trig
  ./ex53  -sol_type quadratic_trig -dm_refine 2 \
      -displacement_petscspace_degree 2 -tracestrain_petscspace_degree 1 -pressure_petscspace_degree 1 \
      -dmts_check .0001 -ts_max_steps 5 -ts_monitor_extreme -dm_view hdf5:./output/sol_2dqt.h5 -monitor_solution \
        hdf5:./output/sol_2dqt.h5::append

elif [ "$test" = "3dqt" ]; then
  # 3D Quad Trig
  ./ex53  -dm_plex_dim 3 -sol_type quadratic_trig -dm_refine 1 \
      -displacement_petscspace_degree 2 -tracestrain_petscspace_degree 1 -pressure_petscspace_degree 1 \
      -dmts_check .0001 -ts_max_steps 5 -ts_monitor_extreme -dm_view hdf5:./output/sol_3dqt.h5 -monitor_solution \
        hdf5:./output/sol_3dqt.h5::append


elif [ "$test" = "terzaghi" ]; then
  # Terzaghi
  ./ex53 -sol_type terzaghi -dm_plex_separate_marker -dm_plex_box_faces 1,8 -simplex 0 -dm_refine 0 \
        -displacement_petscspace_degree 2 -tracestrain_petscspace_degree 1 -pressure_petscspace_degree 1 -niter 16000 \
        -ts_dt 0.0028666667 -ts_max_steps 501 -ts_monitor -dmts_check .0001 -pc_type lu -dm_view hdf5:./output/sol_terzaghi.h5 -monitor_solution \
        hdf5:./output/sol_terzaghi.h5::append

      #  -dm_plex_box_faces 1,64 -ts_max_steps 4 -convest_num_refine 3 gives L_2 convergence rate: [1.1, 1.1, 1.1]
      #  suffix: 2d_terzaghi_tconv
      #  requires: triangle
elif [ "$test" = "terzaghi_tconv" ]; then
   ./ex53 -sol_type terzaghi -dm_plex_separate_marker -dm_plex_box_faces 1,8 -simplex 0 \
        -displacement_petscspace_degree 2 -tracestrain_petscspace_degree 1 -pressure_petscspace_degree 1 -niter 16000 \
        -ts_dt 0.023 -ts_convergence_estimate -pc_type lu -dm_plex_box_faces 1,64 -ts_max_steps 5 -convest_num_refine 7 \
        -dm_view hdf5:./output/sol_terzaghi_tconv.h5 -monitor_solution hdf5:./output/sol_terzaghi_tconv.h5::append -ts_monitor \
        -log_view :./output/run_terzaghi_tconv.py:ascii_info_detail

         #gives L_2 convergence rate: [1.1, 1.1, 1.1]

elif [ "$test" = "terzaghi_sconv" ]; then
  # -dm_plex_box_faces 1,16 -convest_num_refine 4 gives L_2 convergence rate: [1.7, 1.2, 1.1]
  # if we add -displacement_petscspace_degree 3 -tracestrain_petscspace_degree 2 -pressure_petscspace_degree 2, we get [2.1, 1.6, 1.5], so I think we lose an order
  ./ex53 -sol_type terzaghi -dm_plex_separate_marker -dm_plex_box_faces 1,16 -simplex 0 \
    -displacement_petscspace_degree 2 -tracestrain_petscspace_degree 1 -pressure_petscspace_degree 1 -niter 16000 \
    -ts_dt 1e-5 -dt_initial 1e-5 -ts_max_steps 5 -ts_convergence_estimate -ts_convergence_temporal 0 -convest_num_refine 4 -pc_type lu \
    -ts_monitor -convest_monitor -log_view :./output/run_terzaghi_sconv.py:ascii_info_detail


elif [ "$test" = "mandel" ]; then
  ./ex53 -sol_type mandel -dm_plex_separate_marker -simplex 0 -dm_refine 1 \
      -displacement_petscspace_degree 2 -tracestrain_petscspace_degree 1 -pressure_petscspace_degree 1 \
      -ts_dt 0.01 -ts_max_steps 100 -ts_monitor -dmts_check .0001 -pc_type lu -dm_view hdf5:./output/sol_mandel.h5 -monitor_solution hdf5:./output/sol_mandel.h5::append \
      -ts_max_snes_failures -1 -niter 2000
      #-start_in_debugger
    # -dm_plex_box_faces 1,64 -ts_max_steps 4 -convest_num_refine 3 gives L_2 convergence rate: [1.1, 1.1, 1.1]

elif [ "$test" = "mandel_tconv" ]; then
  ./ex53 -sol_type terzaghi -dm_plex_separate_marker -simplex 0 \
      -displacement_petscspace_degree 2 -tracestrain_petscspace_degree 1 -pressure_petscspace_degree 1 -niter 16000 \
      -ts_dt 0.023 -ts_convergence_estimate -pc_type lu -dm_plex_box_faces 1,64 -ts_max_steps 4 -convest_num_refine 3 \
      -dm_view hdf5:./output/sol_mandel_tconv.h5 -monitor_solution hdf5:./output/sol_mandel_tconv.h5::append -convest_monitor -ts_monitor

elif [ "$test" = "cryer" ]; then
  ./ex53 -sol_type cryer \
    -displacement_petscspace_degree 2 -tracestrain_petscspace_degree 1 -pressure_petscspace_degree 1 -dm_refine 1 -bd_dm_refine 1 \
    -ts_dt 0.00028666667 -ts_max_steps 501 -dmts_check .0001 -pc_type lu -pc_factor_shift_type nonzero \
    -dm_view hdf5:./output/sol_cryer.h5 -monitor_solution hdf5:./output/sol_cryer.h5::append -convest_monitor -ts_monitor
# -ts_max_time 0.014333

elif [ "$test" = "cryer_tconv" ]; then
  ./ex53 -sol_type cryer -displacement_petscspace_degree 2 -tracestrain_petscspace_degree 1 -pressure_petscspace_degree 1 -bd_dm_refine 2 \
        -ts_dt 0.023 -ts_max_time 0.092 -ts_convergence_estimate -pc_type lu \
        -bd_dm_refine 3 -ref_limit 0.00666667 -ts_max_steps 5 -convest_num_refine 5 -pc_factor_shift_type nonzero \
        -dm_view hdf5:./output/sol_cryer_tconv.h5 -monitor_solution hdf5:./output/sol_cryer_tconv.h5::append -convest_monitor -ts_monitor
else
  echo "Invalid problem choice."
fi
