source ~/qp2/quantum_package.rc  

# xyz file
file=$1
# multiplicity 
m=$2
# output ezfio file
ezfio=$3

rm -rf $ezfio 
#c=$3
qp create_ezfio -b cc-pvdz $file -m $2 -o $ezfio 
qp set_file ${ezfio}
qp set bsplines_utils bsp_box_size 10. 
qp set bsplines_utils bsp_number 20
qp set bsplines_utils bsp_order 5
qp set bsplines_utils bsp_lmax 0
qp set bsplines_utils bsp_glp 5

qp run bsplines_change_basis | tee ${ezfio}.change.out 
qp run bsplines_check        | tee ${ezfio}.check.out 
