source ~/programs/qp2/quantum_package.rc  

file=$1
m=$2
ezfio=$3

rm -rf $ezfio 

#c=$3
qp create_ezfio -b cc-pvdz $file -m $2 -o $ezfio 
qp set_file ${ezfio}
qp set bsplines_utils bsp_box_size 20. 
qp set bsplines_utils bsp_number 4
qp set bsplines_utils bsp_order 3
qp set bsplines_utils bsp_lmax 3 
qp set bsplines_utils bsp_glp 3

qp run bsplines_ezfio_and_integrals | tee ${ezfio}.change.out
#qp run bsplines_change_basis | tee ${ezfio}.change.out 
#qp run bsplines_check        | tee ${ezfio}.check.out 
