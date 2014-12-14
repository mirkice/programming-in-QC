program hueckel_linear

!use omp_lib

!implicit none

!input form: ./execute num_carbon alpha beta+ beta- distance detail condition

!---------------------delclaration of programming varibles-------------

  integer::i,j,k
!----------------------declaration of parameters-----------------------

  integer,parameter::num_il=5  !number of output in the line
!----------------------declaration of varibles--------------------------

  integer::num_carbon,num_elec,detail,cond,lwork,info
  real*8,allocatable::mat_bp(:,:),mat_bm(:,:),mat_h(:,:),eigs(:),vecs(:,:),mat_hfp(:,:),work(:)  !mat_hfp: matrix h from operator p
  real*8,allocatable::mat_out(:,:)
  real*8,allocatable::mat_r(:,:),mat_p(:,:),mat_u(:,:),mat_v(:,:)         
  real*8::bp,bm,Rsq,R,dist,TMS,TPS,alpha,PS,MS
  real*8,parameter::Pi=3.1415926
  character(len=32)::arg
  
!----------------------construction of the matrix----------------------


  call getarg(1,arg)
  read(arg,*) num_carbon

  call getarg(2,arg)
  read(arg,*) alpha

  call getarg(3,arg)
  read(arg,*) bp             !beta plus

  call getarg(4,arg)
  read(arg,*) bm             !beta minus

  call getarg(5,arg)
  read(arg,*) dist           !distance

  call getarg(6,arg)
  read(arg,*) detail         !detail=1, detailed output. detail=0, simple output, detail=2: simple tps on u,v , detail=3: detailed tps on u,v

  call getarg(7,arg)
  read(arg,*) cond           !cond=1 linear, cond=2 boundary condition


!--------------------------generating hamiltonian--------------------

  num_elec=num_carbon
  allocate(mat_bp(num_elec,num_elec),mat_bm(num_elec,num_elec))
  allocate(mat_h(num_elec,num_elec))
  allocate(eigs(num_elec),vecs(num_elec,num_elec))
  allocate(mat_out(num_elec+3,num_elec))
  mat_bp=0
  mat_bm=0

  do k=1,num_elec
     if( mod(k,2) .ne. 0 ) then
        mat_bp(k,k+1)=1
        mat_bp(k+1,k)=1
     else
        mat_bm(k,k+1)=1
        mat_bm(k+1,k)=1
     endif
  enddo

  if (cond .eq. 2) then
     mat_bm(1,num_elec)=1
     mat_bm(num_elec,1)=1
  endif


  vecs=0

  mat_h=0
  do k=1,num_elec
     mat_h(k,k)=alpha
  enddo
  mat_h=mat_h+mat_bp*bp+mat_bm*bm

  do i=1,num_elec
     do j=i,num_elec
        vecs(i,j)=mat_h(i,j)
     enddo
  enddo

  lwork=240*num_elec

  allocate(work(lwork))

  call dsyev('V','U',num_elec,vecs,num_elec,eigs,work,lwork,info)

  if (info .ne. 0) then
     write(*,*) 'error in calculation'
     stop
  endif
 

!-----------------------output of matrix--------------------
  
  call print_mat(eigs, vecs, num_elec, num_il, mat_out,detail)

!----------------------construction of matrix mat_r----------   

  allocate(mat_r(num_elec,num_elec))
  mat_r=0
  
  do i=1,num_elec
     if (cond .eq. 1) then
        mat_r(i,i)=-dist*(num_elec/2+0.5-i)
     elseif (cond .eq. 2) then
        mat_r(i,i)=num_elec*dist*sin(-2*Pi*(num_elec/2+0.5-i)/(num_elec-0.0))/(2*Pi)
     else
        write(*,*) "condition not specified!"
        stop
     endif
  enddo
!----------------------construction of matrix mat_p----------

  allocate(mat_p(num_elec,num_elec))
  mat_p=0

  do i=1,num_elec-1
     mat_p(i,i+1)=1.0/(2*dist)
     mat_p(i+1,i)=-1.0/(2*dist)
  enddo

  if (cond .eq. 2) then
     mat_p(1,num_elec)=-1.0/(2*dist)
     mat_p(num_elec,1)=1.0/(2*dist)
  endif

  allocate(mat_hfp(num_elec,num_elec))
  mat_hfp=matmul(mat_p,mat_p)/(-2.0)



!----------------------calculation of TPS and TMS------------

select case (detail)
   case (0,1)
      call cal_property(mat_r,mat_r,mat_out,num_elec,TPS)
      call cal_property(mat_p,mat_p,mat_out,num_elec,TMS)

      if (detail .ne. 0) then
         write(*,304) "The result of TPS and TMS with num of carbon:" ,num_elec,TPS,TMS
      else
         write(*,305) num_elec,TPS,TMS
      endif
!--------------------calculation of seperate tensor------------

   case (2,3)
      allocate(mat_u(num_elec,num_elec))
      allocate(mat_v(num_elec,num_elec))

      do i=1,num_elec
         do j=i,num_elec
            call ext_mat(mat_r,num_elec,i,mat_u)
            call ext_mat(mat_r,num_elec,j,mat_v)
            call cal_property(mat_u,mat_v,mat_out,num_elec,PS)
            write(*,306) num_elec,i,j,PS
            if (i .ne. j) then
               write(*,306) num_elec,j,i,PS
            endif
         enddo
      enddo

endselect

306 format(3I10,F20.9)
305 format(I10,2F20.9)
304 format(A50,I10,2F20.9)

303 format(10F20.9)

!-----------------------functions--------------------------






contains

!---------------------output in format---------------
  subroutine print_mat ( eigs, vecs, num_elec, num_dpl, mat_out, detail )  !num_dpl:number of data per line
    integer::num_elec,i,electron(num_elec),j,k,line(num_dpl),orb_count,detail
    real*8,intent(in)::eigs(:),vecs(:,:)
    real*8::eigsr(num_elec),vecsr(num_elec,num_elec)
    real*8,intent(out)::mat_out(num_elec+3,num_elec)
    eigsr=real(eigs)
    vecsr=real(vecs)
    electron=0
    electron(1:num_elec/2)=2
    if (mod(num_elec,2) .eq. 1) then
       electron(num_elec/2+1)=1
    endif
    orb_count=num_elec

!    mat_out(1,:)=(/()/)

    if ( mod(detail,2) .ne. 0 ) then
       write(*,*)   "------------------------------------Here's the information of the orbitals---------------------------------"
       do i=1,num_elec/(num_dpl+1)+1
          line=(/(j,j=num_dpl*(i-1)+1,num_dpl*i)/)
          write(*,300) "Orb_num",line(1:min(num_dpl,orb_count))
          write(*,300) "Electron",electron(line(1:min(num_dpl,orb_count)))
          write(*,301) "Energy ",eigsr(line(1:min(num_dpl,orb_count)))
          write(*,300) "Vectors"
          do k=1,num_elec
             write(*,301) "       ",(vecsr(k,j),j=line(1),line(min(num_dpl,orb_count)))
             mat_out(k+3,line(1:min(num_dpl,orb_count)))=(/(vecsr(k,j),j=line(1),line(min(num_dpl,orb_count)))/)
          enddo
          write(*,*) "----------------------------------------------------------------------------------------------------------"
          mat_out(1,line(1:min(num_dpl,orb_count)))=(/line(1:min(num_dpl,orb_count))/)
          mat_out(2,line(1:min(num_dpl,orb_count)))=(/electron(line(1:min(num_dpl,orb_count)))/)
          mat_out(3,line(1:min(num_dpl,orb_count)))=(/eigsr(line(1:min(num_dpl,orb_count)))/)
 
          orb_count=orb_count-num_dpl
       enddo
    else
       mat_out(2,:)=electron
       do i=1,num_elec
          mat_out(1,i)=i
          mat_out(3,i)=eigsr(i)
          mat_out(4:num_elec+3,i)=vecsr(:,i)
       enddo
    endif


300 format(A15,10I20)
301 format(A15,10F20.9)

    endsubroutine

!-------------------calculate the property ( TPS, TMS ...)-----------------------------------

    subroutine cal_property(mat_oper1,mat_oper2,mat_out,num_elec,output)
      real*8,intent(in)::mat_oper1(num_elec,num_elec),mat_oper2(num_elec,num_elec),mat_out(num_elec+3,num_elec)
      integer,intent(in)::num_elec
      real*8,intent(out)::output
      integer::i,j,k
      real*8::t1,t2
      
      output=0

      !$omp parallel &
      !$omp private(i,j,t1,t2) &
      !$omp shared(num_elec,mat_out,mat_oper1,mat_oper2)
      !$omp do reduction(+:output)
      do i=1,num_elec/2
!         write(*,*) i,OMP_GET_THREAD_NUM()
         do j=num_elec/2+1,num_elec

!               t1=dot_product(matmul(mat_out(4:num_elec+3,i),mat_oper1),mat_out(4:num_elec+3,j))
!               t2=dot_product(matmul(mat_out(4:num_elec+3,i),mat_oper2),mat_out(4:num_elec+3,j))
            call eff_mat_m(mat_oper1,mat_out(4:num_elec+3,i),mat_out(4:num_elec+3,j),num_elec,t1)
            call eff_mat_m(mat_oper2,mat_out(4:num_elec+3,i),mat_out(4:num_elec+3,j),num_elec,t2)
               output=output+t1*t2*2
         enddo
      enddo
      !$omp end do
      !$omp end parallel

    endsubroutine

!------------------extract matrix component-----------------------

    subroutine ext_mat(mat_in,dime,index,mat_out)
      integer,intent(in)::dime,index
      real*8,intent(in)::mat_in(dime,dime)
      real*8,intent(out)::mat_out(dime,dime)
      integer::i,j

      do i=1,dime
         do j=1,dime
            if ( i .eq. index .or. j .eq. index ) then
               mat_out(i,j)=mat_in(i,j)
            else
               mat_out(i,j)=0
            endif
         enddo
      enddo

    endsubroutine
            
        
!-----------------efficient matrix vector multiply------------------

     subroutine eff_mat_m(mat_in,vec_left,vec_right,dime,output)             !vec_left*mat_in*vec_right
       integer,intent(in)::dime
       real*8,intent(in)::mat_in(dime,dime),vec_left(dime),vec_right(dime)
       integer::i,j,k,n_band,lda
       real*8,allocatable::A(:,:)
       real*8::vec_temp(dime)
       real*8,intent(out)::output
       
!       n_band=0
!       do i=2,dime
!          if (mat_in(1,i) .ne. 0) then
!             n_band=n_band+1
!          endif
!       enddo
       
       !lda=2*dime+1
!       allocate(A(lda,dime))
          
!       do 20,j = 1, dime
!          k = n_band + 1 - j
!          do 10,i = MAX( 1, j - n_band ), MIN( dime, j + n_band )
!             A( k + i, j ) = mat_in( i, j )
!10 continue
!20 continue
!       write(*,*) A
 
       write(*,*) mat_in
       stop
       vec_temp=0
       call dgemv('n',dime,dime,1,mat_in,dime,vec_left,1,0,vec_temp,1)
!       write(*,*) vec_temp
       output = dot_product(vec_temp,vec_right)
       
     endsubroutine eff_mat_m


end
