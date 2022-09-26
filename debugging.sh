pylith --petsc.fp_trap --petsc.start_in_debugger=noxterm step01_ricker.cfg 


PETSC_OPTIONS="-malloc_debug 0 -fp_trap -start_in_debugger noxterm" ./.libs/test_faultkinporo --test=pylith::mmstests::TestFaultKinPoro2D_MMS_TriP1::testDiscretization

b libsrc/pylith/fekernels/FaultCohesiveKinPoro.cc:506



knepley/feature-hybrid-mass
baagaard/feature-prescribed-slip-dynamic-new

put nanguards in kernel functions


b fekernels/FaultCohesiveKinPoro.cc:1466


p s_x[sOff_x[sOff[1] + 1]]@4
 p s_x[sOff_x[sOff[1]]]@4    
 
 
 p ((Vec_Seq*)x->data)->array@834

 pip install --target=$PWD gmsh
export PYTHONPATH=${PYTHONPATH}:${TOPSRC_DIR}/gmsh


pylith --petsc.start_in_debugger=noxterm step01_slip_cubit_MMS.cfg
