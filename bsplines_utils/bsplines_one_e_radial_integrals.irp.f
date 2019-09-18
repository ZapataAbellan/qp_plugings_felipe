 BEGIN_PROVIDER [double precision,  dr_int  , (bsp_number,bsp_number)]
&BEGIN_PROVIDER [double precision, d2r_int  , (bsp_number,bsp_number)]
&BEGIN_PROVIDER [double precision,   r_int  , (bsp_number,bsp_number)]
&BEGIN_PROVIDER [double precision,  r1_int  , (bsp_number,bsp_number)]
&BEGIN_PROVIDER [double precision,  r2_int  , (bsp_number,bsp_number)]
&BEGIN_PROVIDER [double precision,   o_int  , (bsp_number,bsp_number)]


 implicit none

 !carl de boor's  parameters for b-splines...
 integer::nderiv,left,mflag,leftmk
 real(8),allocatable::work(:,:)
 real(8),allocatable::values(:,:)

 !gauss-legendre points and weights...
 integer::rglp
 real(8)::r1,r2
 real(8),allocatable::ri(:),wi(:)

 !integers for loops...
 integer::i,j,l,v

 allocate(ri(bsp_glp))
 allocate(wi(bsp_glp))
 allocate(work(bsp_order,bsp_order))
 nderiv=3
 allocate(values(bsp_order,nderiv))

 !integration...
 do v=1,bsp_nv
    r1=bsp_grid_knot(bsp_order+v-1)
    r2=bsp_grid_knot(bsp_order+v)
    call gauleg(r1,r2,ri,wi,bsp_glp)
    do rglp=1,bsp_glp
       values(:,:)=0.d0
       call interv(bsp_grid_knot,bsp_number_knot,ri(rglp),left,mflag)
       leftmk = left - bsp_order
       call bsplvd(bsp_grid_knot,bsp_number_knot,bsp_order,ri(rglp),left,work,values,nderiv)
       !do i=1, bsp_order
       !  write(*,*) 'i,values(i,1)',i,values(i,1)
       !enddo
       do j=1,bsp_order
          do i=1,bsp_order
              dr_int(leftmk+i,leftmk+j)=wi(rglp)*values(i,1)*values(j,2)+dr_int(leftmk+i,leftmk+j) 
             d2r_int(leftmk+i,leftmk+j)=wi(rglp)*values(i,1)*values(j,3)+d2r_int(leftmk+i,leftmk+j)
               r_int(leftmk+i,leftmk+j)=wi(rglp)*values(i,1)*values(j,1)*ri(rglp)+r_int(leftmk+i,leftmk+j)
              r1_int(leftmk+i,leftmk+j)=wi(rglp)*values(i,1)*values(j,1)*1.d0/ri(rglp)+r1_int(leftmk+i,leftmk+j)
              r2_int(leftmk+i,leftmk+j)=wi(rglp)*values(i,1)*values(j,1)*1.d0/ri(rglp)**2+r2_int(leftmk+i,leftmk+j)
               o_int(leftmk+i,leftmk+j)=wi(rglp)*values(i,1)*values(j,1)+o_int(leftmk+i,leftmk+j) 
          enddo  
       enddo
    enddo
 enddo
 END_PROVIDER 
