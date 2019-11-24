source /home/oem/Documentos/UPMC/LCT/LCT_FELIPE/PROJECTS/quantum_package/qp2/quantum_package.rc  


file=$1
m=$2
ezfio=$3

rm -rf $ezfio 
#c=$3
qp create_ezfio -b cc-pvdz $file -m $2 -o $ezfio 
qp set_file ${ezfio}
qp set bsplines_utils bsp_box_size 20. 
qp set bsplines_utils bsp_number 8
qp set bsplines_utils bsp_order 4
qp set bsplines_utils bsp_lmax 4 
qp set bsplines_utils bsp_glp 4

qp run bsplines_two_e_give_integrals  | tee ${ezfio}.change.out
#qp run bsplines_change_basis | tee ${ezfio}.change.out 
#qp run bsplines_check        | tee ${ezfio}.check.out 
