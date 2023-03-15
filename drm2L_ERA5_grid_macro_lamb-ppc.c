#!/bin/csh

@ y1 = 1980
@ y2 = 2020

@ inyr = $y1
while( $inyr <= $y2)

if (($inyr % 4) == 0) then
   set days = 381
else
   set days = 380
endif


cat >! drm2L_grid.f90 << ENDc
module param
! this module is used to define the dimensions of large domain and the time duration of available data
  integer,parameter :: bigx=123  !number of grids in x and y direction
  integer,parameter :: bigy=90
  integer,parameter :: days=${days}
  real(8),parameter    :: dx=75000.00  !unit(m), the distance of each grid
  real(8),parameter    :: dy=75000.00
!__________________________starting dates and ending dates,core region
  integer,parameter :: day1=16, day2=${days} !date starts and ends simulation
  integer,parameter :: domsize=194  !number of grids within the monsoon region
  integer,parameter :: region_num=8932    !number of regions in map file
!______________________________________________________________________
  integer,parameter :: dt=1800  !time interval of each-time back tracing
  integer,parameter :: max_tracing=1000 ! maximum times of back trajactory(>48*15)
  integer,parameter :: max_iteration=100 ! maximum times of iteration
  integer,parameter :: warmup=15 !days that to be excluded at the beginning
  integer,parameter :: velocity_t=3  !time interval of velocity data
  integer,parameter :: valve1=6  !times of tracing per unit velocity time interval (=3600*3/1800)
  integer,parameter :: valve2=48  !times of tracing per day (=3600*24/1800)
  integer,parameter :: NAN_value=-9.99E8
  integer,parameter :: ithyear=${inyr}
end module param
!-------------------------
module typedef
  implicit none
  type start_end
    integer :: start_l,end_l,renum
  end type
end module

!====================================================
program Recycling

! load data first
! data includes daily data and several hourly data(in this case 3 hourly):
! Daily: precipitable water(PW),evaportranspiration(ET),precipitation(PP)
! 3-hourly: u and v
  use param
  use typedef
  implicit none
  real(8), dimension(bigx,bigy,days,2)                 :: grid_1d       ! 1-ET(mm/3h),2-PP(mm/day)
  real(8), dimension(bigx,bigy,days,2)                 :: grid_pw       ! 1-PWU(mm),2-PWL(mm)
  real(8), dimension(bigx,bigy,days*24/velocity_t,2)   :: grid_uv_u     ! 1-U3U(m/s),2-V3U(m/s)
  real(8), dimension(bigx,bigy,days*24/velocity_t,2)   :: grid_uv_l     ! 1-U3L(m/s),2-V3L(m/s)
  real(8), dimension(bigx,bigy,days*24/velocity_t,2)   :: grid_w700     ! 1-W3U(m/s),2-W3D(m/s)
  real(8), dimension(bigx,bigy,days*24/velocity_t  )   :: grid_q700     ! q at 700mb(kg/kg)
  integer, dimension(day1:day2)                     :: ithday        ! read from 1L-DRM
  real(8), dimension(day1:day2)                        :: rho_1l        ! recycling ratio read from 1L-DRM
  integer, dimension(bigx,bigy)                     :: mask_700
  integer, dimension(bigx,bigy)                     :: map
  real(8), dimension(bigx,bigy,days)                   :: grid_qw1d     
  integer, dimension(bigx,bigy,days)                :: qw_mask       ! 0-UL grids, 1-upward qw, 2-downward qw

  real(8)                                              :: EW_RATIO=0.,QWU_RATIO=0.,QWD_RATIO=0.,TPP=0.,EWL_RATIO=0.0
  real(8)                                              :: TU_1=0.,TV_1=0.,TU_2=0.,TV_2=0.
  real(8), dimension(max_tracing,7)                    :: EW_vector_u = NAN_value !matrix used to save the E/W data of each grid on specific day and the region it belongs as well
  real(8), dimension(max_tracing,8)                    :: EW_vector_l = NAN_value
  real(8)                                              :: E_W=0  !matrix to save integrated EW_vector of each grid on specific day
  real(8), dimension(region_num)                       :: sum_upper, sum_upper_pp, sum_lower, sum_lower_pp, sum_upper1l, sum_lower1l  
  integer                                           :: i,j,k,day,time,it,intx,inty,intx1,inty1,flag1,flag2,rn,t
  real(8)                                              :: oldppx,oldppy,diffx,diffy
  real(8), dimension(max_tracing)                      :: ppxu,ppyu,ppxl,ppyl  !variable used to record the position of tracing
  integer                                           :: daycount,hourcount,count1,count2
  real(8)                                              :: temprr,sumrr,sumpp,s, grid_rr
  real(8), dimension(day1:day2)                        :: rr
  integer                                           :: fives, pentads_bef,pentads_aft, file_id    ! these variables are defined for pentads analysis
!-------new variables for computing the contribution of different regions-------
  type(start_end), dimension(max_tracing)           :: reg_count_u
  type(start_end), dimension(max_tracing)           :: reg_count_l
  real(8),            dimension(region_num)            :: daily_rru, daily_rrl, daily_rra, daily_rru1l, daily_rrl1l  !record daily recycling ratio of different regions
  integer,         dimension(domsize,2)             :: domainij
  real(8),            dimension(2)                     :: precip_grid=0  ! 1-upper, 2-lower
  integer                                           :: reg_change,old_reg,new_reg,zeroi,ein,xx,ind,L1,L2
  integer                                           :: reg_change_u, reg_change_l
  real(8)                                              :: s1,pral
  real(8),            dimension(region_num)            :: rho_ave, rho_upper, rho_lower, rho_upper1l, rho_lower1l

! load data into array grid_days and grid_hours
  call openfile(grid_1d,grid_pw,grid_uv_u,grid_uv_l,grid_w700,grid_q700)  !including fixing data <0
  call readmap(map,mask_700)
  xx=map(75,75)
  print *,xx
  xx=mask_700(75,75)
  print *,xx
  call monsoon(domainij)
  call daily_qw(grid_w700,grid_q700,grid_qw1d,qw_mask)

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!############################################################################################
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

! start the daily loop during the chosen duration from day1 to day2
  print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  print *,'!!!  Daily Loop Starts  !!!'
  print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  daily : do day=day1,day2

    print *,'day:', day
    sum_upper = 0.
    sum_lower = 0.
    sum_upper_pp = 0.
    sum_lower_pp = 0.
    sum_upper1l = 0.
    sum_lower1l = 0.


    !------partition precipitation into 2 layers-------
    domain : do k=1,domsize 
      precip_grid(:)=0.
      i=domainij(k,1)
      j=domainij(k,2)
      TPP=grid_1d(i,j,day,2)
      reg_change_u = 0.
      reg_change_l = 0.
      if (mask_700(i,j) == 0) then
         precip_grid(1) = TPP
      else
         precip_grid(1) = TPP*grid_pw(i,j,day,1)/(grid_pw(i,j,day,1)+grid_pw(i,j,day,2))
         precip_grid(2) = TPP*grid_pw(i,j,day,2)/(grid_pw(i,j,day,1)+grid_pw(i,j,day,2))
      end if
!=======================================================
!-------(1) run DRM for the upper level-----------------
!=======================================================
        ppxu = 0.0
        ppyu = 0.0
        ppxu(1)=(i-1)*dx+dx/2
        ppyu(1)=(j-1)*dy+dy/2
        E_W=0
        
        count1=0 
        count2=0
        reg_change=0
        hourcount=day*24/velocity_t   !starting temporal position in velocity files
        daycount=day                  !starting temporal position in other files
     
        ! allocate all to be zero
        do zeroi=1,max_tracing
          do ein=1,7
            EW_vector_u(zeroi,ein)=0
          end do
          reg_count_u(zeroi)%start_l=0
          reg_count_u(zeroi)%end_l=0
          reg_count_u(zeroi)%renum=0
        end do

!====================================================================================================
        ! start the time loop
        timeloopU : do time=1,max_tracing-1

          ! tell ppx and ppy are within large domain or not-----------------------------------
          if ((ppxu(time)<=0).or.(ppxu(time)>=(dx*bigx))) then
           ! print *,'Beyond boundary during timeloop when time is:',time
            exit timeloopU
          else if ((ppyu(time)<=0).or.(ppyu(time)>=(dy*bigy))) then
           ! print *,'Beyond boundary during timeloop when time is:',time
            exit timeloopU
          end if
          !-----------------------------------------------------------------------------------

          !read the Evaporation,precipitable water and wind velocity by subroutine 'finduvEW'
          call finduvEW_U(ppxu(time),ppyu(time),hourcount,daycount,grid_pw(:,:,daycount,:),grid_1d(:,:,daycount,1),&
              grid_uv_u(:,:,hourcount,1),grid_uv_u(:,:,hourcount,2),grid_qw1d(:,:,daycount),mask_700,&
              qw_mask(:,:,daycount),EW_RATIO,QWU_RATIO,TU_1,TV_1,intx1,inty1,flag1)
!        write(500,'(A5,I1)') 'flag=',flag1
        if (flag1==5) then
          exit timeloopU
        else
          EW_vector_u(time,1)=(time-1)*dt
          EW_vector_u(time,2)=EW_RATIO
          EW_vector_u(time,3)=QWU_RATIO
          EW_vector_u(time,4)=map(intx1,inty1)
          EW_vector_u(time,5)=real(flag1)       ! if flag1 == 0, then don't need to be multiplied by rho_lower
          EW_vector_u(time,6)= intx1
          EW_vector_u(time,7)= inty1
          !TET in unit of mm/3h thus TET*8 and EW_vector*8 has the unit of mm/day and day-1 individually
        end if
 
          if (time .eq. 1) then
            reg_count_u(1)%start_l=1
            reg_count_u(1)%renum=map(i,j)
            reg_change=1
          else
            old_reg=EW_vector_u(time-1,4)
            new_reg=EW_vector_u(time,4)
            if (old_reg .ne. new_reg) then
              reg_change=reg_change+1
              reg_count_u(reg_change-1)%end_l=time
              reg_count_u(reg_change)%start_l=time
              reg_count_u(reg_change)%renum=new_reg      
            end if
          end if
  
          ! find the next point using Iterative Technique--------------------------------------
          call finduv(ppxu(time),ppyu(time),hourcount,grid_uv_u(:,:,hourcount,1),& 
               grid_uv_u(:,:,hourcount,2),TU_2,TV_2)
          oldppx=ppxu(time)-(TU_1+TU_2)*dt/2
          oldppy=ppyu(time)-(TV_1+TV_2)*dt/2
          if ((oldppx<=0).or.(oldppx>=(dx*bigx))) exit timeloopU
          if ((oldppy<=0).or.(oldppy>=(dy*bigy))) exit timeloopU
!================================================================================================
          iterationU : do it=1,max_iteration
            !tell if the point still within the large domain
            if ((oldppx<=0).or.(oldppx>=(dx*bigx))) then
            !  print *,'Beyond boundary during iteration when it is',it
              exit iterationU
            else if ((oldppy<=0).or.(oldppy>=(dy*bigy))) then
            !  print *,'Beyond boundary during iteration when it is',it
              exit iterationU
            end if           

            call finduv(oldppx,oldppy,hourcount,grid_uv_u(:,:,hourcount,1),&
                 grid_uv_u(:,:,hourcount,2),TU_2,TV_2,intx,inty,flag2)
          if (flag2==1) then
             exit iterationU
          else
            ppxu(time+1)=ppxu(time)-(TU_1+TU_2)*dt/2
            ppyu(time+1)=ppyu(time)-(TV_1+TV_2)*dt/2
            
            diffx=abs(ppxu(time+1)-oldppx)
            diffy=abs(ppyu(time+1)-oldppy)
            if ((diffx<0.0001).and.(diffy<0.0001)) exit  iterationU
            oldppx=ppxu(time+1)
            oldppy=ppyu(time+1)
          end if
          end do iterationU
!================================================================================================
        ! change the temporal data to use based on back tracing------------------------------
         count1=count1+1
         if (mod(count1,valve1)==0) then
           hourcount=hourcount-1
           if(hourcount<=0) exit timeloopU
         end if
         count2=count2+1
         if (mod(count2,valve2)==0)  then
           daycount=daycount-1
           if(daycount<=0) exit timeloopU
         end if
 
        end do timeloopU
        if (count2>0) then
           reg_count_u(reg_change)%end_l=count2
        end if
        reg_change_u = reg_change

        call DOM_SPTNQU(EW_vector_u,reg_count_u,reg_change_u, rho_upper,rho_upper1l)
!=================================================================
!---------(2) run DRM for the lower level-------------------------
!=================================================================
        ppxl = 0.0
        ppyl = 0.0
        ppxl(1)=(i-1)*dx+dx/2
        ppyl(1)=(j-1)*dy+dy/2

        count1=0
        count2=0
        reg_change=0
        hourcount=day*24/velocity_t   !starting temporal position in velocity files
        daycount=day                  !starting temporal position in other files

        ! allocate all to be zero
        do zeroi=1,max_tracing
          do ein=1,8
            EW_vector_l(zeroi,ein)=0
          end do
          reg_count_l(zeroi)%start_l=0
          reg_count_l(zeroi)%end_l=0
          reg_count_l(zeroi)%renum=0
        end do
        reg_change = 0

!====================================================================================================
        ! start the time loop
        timeloopL : do time=1,max_tracing-1

          ! tell ppx and ppy are within large domain or
          ! not-----------------------------------
          if ((ppxl(time)<=0).or.(ppxl(time)>=(dx*bigx))) then
           ! print *,'Beyond boundary during timeloop when time is:',time
            exit timeloopL
          else if ((ppyl(time)<=0).or.(ppyl(time)>=(dy*bigy))) then
           ! print *,'Beyond boundary during timeloop when time is:',time
            exit timeloopL
          end if
          !-----------------------------------------------------------------------------------

          !read the Evaporation,precipitable water and wind velocity by
          !subroutine 'finduvEW'
          call finduvEW_L(ppxl(time),ppyl(time),hourcount,daycount,grid_pw(:,:,daycount,:),grid_1d(:,:,daycount,1),&
              grid_uv_l(:,:,hourcount,1),grid_uv_l(:,:,hourcount,2),grid_qw1d(:,:,daycount), mask_700,&
              qw_mask(:,:,daycount),EW_RATIO,EWL_RATIO,QWD_RATIO,TU_1,TV_1,intx1,inty1,flag1)
        if (flag1==5) then
          if (reg_change > 0) then
            reg_count_l(reg_change)%end_l=time
          else
            reg_count_l(1)%start_l = 0
            reg_count_l(1)%end_l   = 0
            reg_count_l(1)%renum   = 0
          end if
          exit timeloopL
        else
          EW_vector_l(time,1)=(time-1)*dt
          EW_vector_l(time,2)=EW_RATIO
          EW_vector_l(time,3)=EWL_RATIO
          EW_vector_l(time,4)=QWD_RATIO
          EW_vector_l(time,5)=map(intx1,inty1)
          EW_vector_l(time,6)=real(flag1)       ! if flag1 == 0, then don't need to be multiplied by rho_upper
          EW_vector_l(time,7)=intx1
          EW_vector_l(time,8)=inty1
        end if

          if (time .eq. 1) then
            reg_count_l(1)%start_l=1
            reg_count_l(1)%renum=map(i,j)
            reg_change=1
          else
            old_reg=EW_vector_l(time-1,5)
            new_reg=EW_vector_l(time,5)
            if (old_reg .ne. new_reg) then
              reg_change=reg_change+1
              reg_count_l(reg_change-1)%end_l=time
              reg_count_l(reg_change)%start_l=time
              reg_count_l(reg_change)%renum=new_reg
            end if
          end if

! find the next point using Iterative Technique--------------------------------------
          call finduv(ppxl(time),ppyl(time),hourcount,grid_uv_l(:,:,hourcount,1),&
               grid_uv_l(:,:,hourcount,2),TU_2,TV_2)
!print * ,'after finduv'
          oldppx=ppxl(time)-(TU_1+TU_2)*dt/2
          oldppy=ppyl(time)-(TV_1+TV_2)*dt/2
          if ((oldppx<=0).or.(oldppx>=(dx*bigx))) exit timeloopL
          if ((oldppy<=0).or.(oldppy>=(dy*bigy))) exit timeloopL
!================================================================================================
          iterationL : do it=1,max_iteration
!print *,'in iterationL:'
            !tell if the point still within the large domain
            if ((oldppx<=0).or.(oldppx>=(dx*bigx))) then
            !  print *,'Beyond boundary during iteration when it is',it
              exit iterationL
            else if ((oldppy<=0).or.(oldppy>=(dy*bigy))) then
            !  print *,'Beyond boundary during iteration when it is',it
              exit iterationL
            end if

            call finduv(oldppx,oldppy,hourcount,grid_uv_l(:,:,hourcount,1),&
                 grid_uv_l(:,:,hourcount,2),TU_2,TV_2,intx,inty,flag2)
!print *,'after finduv in iterationL.'
          if (flag2==1) then
             exit iterationL
          else
            ppxl(time+1)=ppxl(time)-(TU_1+TU_2)*dt/2
            ppyl(time+1)=ppyl(time)-(TV_1+TV_2)*dt/2

            diffx=abs(ppxl(time+1)-oldppx)
            diffy=abs(ppyl(time+1)-oldppy)
            if ((diffx<0.0001).and.(diffy<0.0001)) exit  iterationL
            oldppx=ppxl(time+1)
            oldppy=ppyl(time+1)
          end if
          end do iterationL
!================================================================================================
        ! change the temporal data to use based on back
        ! tracing------------------------------
         count1=count1+1
         if (mod(count1,valve1)==0) then
           hourcount=hourcount-1
           if(hourcount<=0) exit timeloopL
         end if
         count2=count2+1
         if (mod(count2,valve2)==0)  then
           daycount=daycount-1
           if(daycount<=0) exit timeloopL
         end if

        end do timeloopL
        if (count2>0) then
           reg_count_l(reg_change)%end_l=count2
        end if
        reg_change_l = reg_change

        call DOM_SPTNQL(EW_vector_l,reg_count_l,reg_change_l, rho_lower, rho_lower1l)

!===================================================================================================
     do rn = 1, region_num
       sum_upper(rn) = sum_upper(rn) + precip_grid(1)*rho_upper(rn)
       sum_lower(rn) = sum_lower(rn) + precip_grid(2)*rho_lower(rn)
       sum_upper1l(rn) = sum_upper1l(rn) + TPP*rho_upper1l(rn)
       sum_lower1l(rn) = sum_lower1l(rn) + TPP*rho_lower1l(rn)
     end do
     sum_upper_pp = sum_upper_pp + precip_grid(1)
     sum_lower_pp = sum_lower_pp + precip_grid(2)

   end do domain 
   daily_rru(:) = sum_upper/sum_upper_pp
   daily_rrl(:) = sum_lower/sum_lower_pp
   daily_rra(:) = (sum_upper+sum_lower)/(sum_upper_pp+sum_lower_pp)
   daily_rru1l(:)  = sum_upper1l/(sum_upper_pp+sum_lower_pp)
   daily_rrl1l(:)  = sum_lower1l/(sum_upper_pp+sum_lower_pp)

   call output_daily_grid(day, daily_rru,   'UPPER2L')
   call output_daily_grid(day, daily_rrl,   'LOWER2L')
   call output_daily_grid(day, daily_rra,   'AVERA2L')
   call output_daily_grid(day, daily_rru1l, 'UPPER1L')
   call output_daily_grid(day, daily_rrl1l, 'LOWER1L')

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   end do daily

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!##################################################################################################
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
end program Recycling


!===========================================================================
subroutine output_daily_grid(dd, outdata, fname)
  use param
  implicit none
  integer, intent(in) :: dd
  real(8), dimension(region_num), intent(in) :: outdata
  character(len=7), intent(in) :: fname
  character(len=90) :: chfile1, path, chfile2, chfile3
  integer :: kk
 
  write(path,'(A45,I4)') '/data/keeling/a/seunguk2/c/2LDRM/pydlmb/GRID/',ithyear
  if (dd < 100) then
    write(chfile3,'(A1,I2)') '0',dd
  else
    write(chfile3,'(I3)') dd
  end if
  write(chfile1,'(A50,A8,A3,A1,A7,A4)') path,'_8MLCpc_',chfile3,'_',fname,'.dat'
  open(unit=700,file=chfile1,status='unknown',form='unformatted')
!  write(700,*)((outdata(kk,rr),rr=1,region_num),kk=1,domsize)
  write(700) real(outdata, kind(1e0))
  close(700)

end subroutine
!--------------------------------------------------------------------------
!subroutine output_ewvector(dd, vector_u, vector_l)
!  use param
!  implicit none
!  integer, intent(in) :: dd
!  real(8), dimension(domsize, max_tracing, 7), intent(in) :: vector_u
!  real(8), dimension(domsize, max_tracing, 8), intent(in) :: vector_l
!  character(len=70) :: chfile1,chfile2,path,chfile3
!  integer :: kk, tt, ii
!
!  write(path,'(A29,I4,A1)') '/data/francina/b/hhu18/2LDRM/',ithyear,'/'
!  if (dd < 100) then
!    write(chfile3,'(A1,I2)') '0',dd
!  else 
!    write(chfile3,'(I3)') dd
!  end if
!  write(chfile1,'(A34,A17,I4,A1,A3,A4)') path,'WRF_OUTPUT_2LBTU_',ithyear,'_',chfile3,'.txt'
!  write(chfile2,'(A34,A17,I4,A1,A3,A4)') path,'WRF_OUTPUT_2LBTL_',ithyear,'_',chfile3,'.txt'
!  open(unit=700,file=chfile1)
!  write(700,*)(((vector_u(kk,tt,ii),ii=1,7),tt=1,max_tracing),kk=1,domsize)
!  close(700)
!  open(unit=800,file=chfile2)
!  write(800,*)(((vector_l(kk,tt,ii),ii=1,8),tt=1,max_tracing),kk=1,domsize)
!  close(800)
!  
!end subroutine
!__________________________________________________________________________
subroutine monsoon(domain)
  use param
  implicit none
  integer,dimension(domsize,2),intent(out) :: domain
  real, dimension(domsize,2) :: domain2
  integer :: ii,error
  character(len=70) :: file_name, file_name2

  file_name='/data/keeling/a/seunguk2/a/2LDRM/pyd075/MW_il.dat'
  file_name2='/data/keeling/a/seunguk2/a/2LDRM/pyd075/MW_jl.dat'
  open(14,file=file_name,status='old',form='unformatted',access='direct',recl=4)
  open(15,file=file_name2,status='old',form='unformatted',access='direct',recl=4)
  do ii = 1,domsize
     read(14,rec=ii) domain2(ii,1)
     read(15,rec=ii) domain2(ii,2)
  enddo
  domain = int(domain2)

  write(*,*) 'Monsoon region has been loaded.'
  close(14)
  close(15)
end subroutine
!___________________________________________________________________________
subroutine openfile(x1d,xpw,xuv1,xuv2,xw,xq)
  use param
  implicit none
  integer :: file_n=5, i,j,k,l,flag
  real(8), dimension(bigx,bigy,days,2),intent(out) :: x1d
  real(8), dimension(bigx,bigy,days,2),intent(out) :: xpw
  real(8), dimension(bigx,bigy,days*24/velocity_t,2),intent(out) :: xuv1
  real(8), dimension(bigx,bigy,days*24/velocity_t,2),intent(out) :: xuv2
  real(8), dimension(bigx,bigy,days*24/velocity_t,2),intent(out) :: xw
  real(8), dimension(bigx,bigy,days*24/velocity_t),  intent(out) :: xq
  real(8), dimension(bigx,bigy,days*24) :: xqi
  real(8), dimension(bigx,bigy,days*24/velocity_t) :: xwi

  character :: n_c
  character(len=70) :: chfile, path, path2
  character(len=90) :: whpath
  character(len=2) :: items(2)=(/'ET','P2'/) ! PP -> P2: precip data corrected 
  character(len=3) :: items2(8)=(/'QWU','QWD','U3U','V3U','U3D','V3D','W3U','W3D'/)
  character(len=4) :: items3(1)=(/'q700'/)
  character(len=5) :: yyyy='_${inyr}'

  write(path,'(A40)') '/data/keeling/a/seunguk2/b/2LDRM/pydlmb/'
  do l=1, 2
    write(chfile,'(A2,A5,A9)') items(l),yyyy,'_lamb.dat'
    write(whpath,'(A40,A16)') path,chfile
    print *, chfile
    open(10,file=whpath,status='old',form='unformatted',access='direct',recl=8*bigx*bigy)
    do k = 1,days
       read(10,rec=k) ((x1d(i,j,k,l),i=1,bigx),j=1,bigy)
    enddo
    write(*,'(A30,A24)') chfile, 'has been opened.'
  end do
  x1d(:,:,:,1) = -x1d(:,:,:,1)
  do i = 1, bigx
    do j = 1, bigy
      do k = 1, days
        x1d(i,j,k,2) = x1d(i,j,k,2) - min(x1d(i,j,k,1),0.0)
        x1d(i,j,k,1) = max(x1d(i,j,k,1),0.0)
      end do
    end do
  end do
  write(*,'(A20)') 'ET/PP correction done'

    write(chfile,'(A4,A5,A9)') items3(1),yyyy,'_lamb.dat'
    write(whpath,'(A40,A18)') path,chfile
    print *, chfile
    open(10,file=whpath,status='old',form='unformatted',access='direct',recl=8*bigx*bigy)
    do k = 1,days*24/velocity_t
       read(10,rec=k) ((xq(i,j,k),i=1,bigx),j=1,bigy)
    enddo
    write(*,'(A30,A20)') chfile, 'has been opened.'

  do l=1, 8
    write(chfile,'(A3,A5,A9)') items2(l),yyyy,'_lamb.dat'
    write(whpath,'(A40,A17)') path,chfile
    print *, chfile
      open(10,file=whpath,status='old',form='unformatted',access='direct',recl=8*bigx*bigy)
    if(l<=2) then
      do k = 1,days
      read(10,rec=k) ((xpw(i,j,k,l),i=1,bigx),j=1,bigy)
      enddo
    else
      if (l<=4) then
        do k = 1,days*24/velocity_t
        read(10,rec=k) ((xuv1(i,j,k,l-2),i=1,bigx),j=1,bigy)
        enddo
      else
        if (l<=6) then
          do k = 1,days*24/velocity_t
          read(10,rec=k) ((xuv2(i,j,k,l-4),i=1,bigx),j=1,bigy)
          enddo
        else
          if (l<=7) then
            do k = 1,days*24/velocity_t
            read(10,rec=k) ((xw(i,j,k,l-6),i=1,bigx),j=1,bigy)
            enddo
          else
            do k = 1,days*24/velocity_t
            read(10,rec=k) ((xwi(i,j,k),i=1,bigx),j=1,bigy)
            enddo
            xw(:,:,:,2) = xwi(:,:,:)
          end if
        end if
      end if
    end if
    write(*,'(A30,A20)') chfile, 'has been opened.'
    close(10)
  end do

end subroutine
!__________________________________________________________________________
subroutine readmap(regions, mask700)
  use param
  implicit none
  integer,dimension(bigx,bigy),intent(out) :: regions
  integer,dimension(bigx,bigy),intent(out) :: mask700
  real,dimension(bigx,bigy) :: regions2
  real,dimension(bigx,bigy) :: mask7002
  character(len=70) :: path,file_name,file_name2
  integer :: i,j

  write(path,'(A40)') '/data/keeling/a/seunguk2/a/2LDRM/pyd075/'
  write(file_name,'(A40,A13)') path, 'ERA5_lcgd.dat'
  write(file_name2,'(A40,A14)') path, '700_masklc.dat'
  open(113,file=file_name2,status='old',form='unformatted',access='direct',recl=4*bigx)
  do j = 1,bigy
    read(113,rec=j) (mask7002(i,j), i=1,bigx)
  enddo
  write(*,*) 'Region Mask has been loaded.'
  close(113)

  open(112,file=file_name,status='old',form='unformatted',access='direct',recl=4*bigx)
  do j = 1,bigy
    read(112,rec=j) (regions2(i,j), i=1,bigx)
  enddo
  close(112)

  write(*,*) '700mb Mask has been loaded.'
  regions = int(regions2)
  mask700 = int(mask7002)
    
end subroutine
!__________________________________________________________________________
subroutine finduvEW(ppx,ppy,hourcount,daycount,W,E,U,V,outW,outE,outU,outV,xint,yint,flag)
  use param
  implicit none
  real(8),intent(in) :: ppx,ppy
  integer,intent(in) :: daycount,hourcount
  real(8),dimension(bigx,bigy),intent(in) ::U,V,E,W
  real(8),intent(out) :: outW,outE,outU,outV
  integer,intent(out) :: xint,yint,flag
  flag=0
  xint=ceiling(ppx/dx)
  yint=ceiling(ppy/dy)
  outW=W(xint,yint)
  outE=E(xint,yint)
  outU=U(xint,yint)
  outV=V(xint,yint)
  if (isnan(outW)) then
    flag=1
  elseif (isnan(outE)) then
    flag=1
  elseif (isnan(outU)) then
    flag=1
  elseif (isnan(outV)) then
    flag=1
  end if

end subroutine
!-------------------------------------------------------------------------
subroutine finduvEW_U(ppx,ppy,hourcount,daycount,W,E,U,V,QW,mask700,maskqw,outEW_ratio,outqwu_ratio,outU,outV,xint,yint,flag)
  use param
  implicit none
  real(8),intent(in) :: ppx,ppy
  integer,intent(in) :: daycount,hourcount
  real(8),dimension(bigx,bigy),intent(in) ::U,V,E,QW
  real(8),dimension(bigx,bigy,2),intent(in) :: W
  integer, dimension(bigx, bigy),intent(in) :: mask700,maskqw
  real(8),intent(out) :: outEW_ratio,outqwu_ratio,outU,outV
  integer,intent(out) :: xint,yint,flag
  integer             :: mask1, mask2
  real(8), dimension(2)  :: W2L
  real(8)                :: ET, W_ratio, W_sum, outW, qwu,outE

  flag=3
  xint=ceiling(ppx/dx)
  yint=ceiling(ppy/dy)
  W2L=W(xint,yint,:)
  outU=U(xint,yint)
  outV=V(xint,yint)
  mask1 = mask700(xint,yint)
  mask2 = maskqw(xint, yint)
  ET=E(xint,yint)
  outE = ET
 
  if (mask1 == 0 .or. mask2 == 0) then  ! only 1-layer
     if (isnan(W2L(2))) then
        outW=W2L(1)
     else
        outW=W2L(1)+W2L(2)
     end if
     W_sum = outW
     qwu  = ET
     outqwu_ratio = qwu/W_sum
     flag = 0    ! only 1-Layer
  else
     outW=W2L(1)
     W_sum = W2L(1)+W2L(2)
     qwu = max(QW(xint,yint),0.0)
     outqwu_ratio = (qwu/outW)*(W_sum/W2L(2))     !(qwu/W_upper)*(W_sum/W_lower)
     flag=mask2     ! upward moisture transport
  end if 
  outEW_ratio = ET/W_sum

  if (isnan(outW)) then
    flag=5
  elseif (isnan(outE)) then
    flag=5
  elseif (isnan(outU)) then
    flag=5
  elseif (isnan(outV)) then
    flag=5
  end if
  
end subroutine
!-------------------------------------------------------------------------
subroutine finduvEW_L(ppx,ppy,hourcount,daycount,W,E,U,V,QW,mask700,maskqw,&
                      outEW_ratio,outEW1_ratio,outQWD_ratio,outU,outV,xint,yint,flag)
  use param
  implicit none
  real(8),intent(in) :: ppx,ppy
  integer,intent(in) :: daycount,hourcount
  real(8),dimension(bigx,bigy),intent(in) ::U,V,E,QW
  real(8),dimension(bigx,bigy,2),intent(in) :: W
  integer, dimension(bigx, bigy),intent(in) :: mask700,maskqw
  real(8),intent(out) :: outEW_ratio,outEW1_ratio,outU,outV,outQWD_ratio
  integer,intent(out) :: xint,yint,flag
  integer             ::  mask1, mask2
  real(8), dimension(2)  :: W2L
  real(8)                ::  outE,ET,Wsum,outW,qwd

  flag=3
  xint=ceiling(ppx/dx)
  yint=ceiling(ppy/dy)
  W2L =W(xint,yint,:)
  outU=U(xint,yint)
  outV=V(xint,yint)
  mask1 = mask700(xint,yint)
  mask2 = maskqw(xint, yint)
  ET    = E(xint,yint)
  outE  = ET

  if (mask1 == 0 .or. mask2 == 0) then
     if (isnan(W2L(2))) then
        outW  = W2L(1)
     else
        outW  = W2L(1)+W2L(2)
     end if
     Wsum  = outW
     outEW1_ratio = NAN_value
     outQWD_ratio  = NAN_value
     flag = 0    ! only 1-Layer
  else
     outW = W2L(1)+W2L(2)
     Wsum = outW
     qwd  = min(QW(xint,yint),0.0)*(-1.0)
     outEW1_ratio = outE/W2L(2)
     outQWD_ratio = (qwd/W2L(2))*(Wsum/W2L(1))
     flag=mask2     ! upward moisture transport
  end if
  outEW_ratio = outE/Wsum

  if (isnan(outW)) then
    flag=5
  elseif (isnan(outE)) then
    flag=5
  elseif (isnan(outU)) then
    flag=5
  elseif (isnan(outV)) then
    flag=5
  end if

end subroutine
!__________________________________________________________________________
subroutine finduv(xp,yp,hourcount,u,v,outu,outv,xint,yint,flag)
  use param
  implicit none
  real(8),intent(in) ::xp,yp
  integer,intent(in) ::hourcount
  real(8),dimension(bigx,bigy),intent(in) ::u,v
  real(8),intent(out)::outu,outv
  integer,intent(out)::xint,yint,flag
  flag=0
  xint=ceiling(xp/dx)
  yint=ceiling(yp/dy)
  outu=u(xint,yint)
  outv=v(xint,yint)
  if (isnan(outu)) then
    flag=1
  else if (isnan(outv)) then
    flag=1
  end if

end subroutine
!____________________________________________________________________________
subroutine SPTNQ(x,y,n,s)
  use param
  implicit none
  integer,intent(in):: n
  real(8),dimension(n),intent(in)::x,y
  real(8),intent(out) :: s
  integer:: i,j,k
  real(8) :: a,b1,b2
  
  s=0
  do i=1,n-1
    a=x(i+1)-x(i)
    b1=y(i)
    b2=y(i+1)
    s=s+a*(b1+b2)/(2*3600*velocity_t)
  !  write(600,*) a, s
  end do
end subroutine
!_____________________________________________________________
subroutine SPTNQ_U(x,y1,y2,z,n,s,rho0)
  use param
  implicit none
  integer,intent(in):: n
  real(8),dimension(n),intent(in)::x,y1,y2,z
  real(8),intent(out) :: s,rho0
  integer:: i,j,k
  real(8) :: a,b1,b2, ewa1,ewa2, ew1,ew2
  real(8) :: ss1, ss2
  real(8) :: bt1, bt2, sst1

  s=0.
  ss1 = 0.0
  ss2 = 0.0
  sst1 = 0.0
  do i=n-1,1,-1
    rho0 = 0.0
    a=x(i+1)-x(i)           ! dt
!-----for total column integration------
    bt1=y1(i)
    bt2=y1(i+1)
    sst1  = sst1 + a*(bt1+bt2)/(2*3600*velocity_t)
    rho0   = 1-exp(-sst1)
!-----for upper layer integration-------
    b1 =y2(i)
    b2 =y2(i+1)
    ewa1=(b1+b2)/(2*3600*velocity_t)
    ss1=ss1+a*ewa1
    ew1=b1*rho0
    ew2=b2*rho0
    if (z(i) == 0.0) then   ! only 1 layer
       ew1=b1
    end if
    if (z(i+1) == 0.0) then
       ew2=b2
    end if
    ewa2=(ew1+ew2)/(2*3600*velocity_t)
    ss2=ss2+exp(ss1)*ewa2*a
  end do
  s = exp(-ss1)*ss2
end subroutine
!-----------------------------------------------------------
subroutine SPTNQ_L(x,y1,y2,y3,z,n,s,rho0)
  use param
  implicit none
  integer,intent(in):: n
  real(8),dimension(n),intent(in)::x,y1,y2,y3,z
  real(8),intent(out) :: s, rho0
  integer:: i,j,k
  real(8) :: a,et1,et2,qw1,qw2
  real(8) :: ewa1,ewa2, ew1,ew2
  real(8) :: ss1, ss2
  real(8) :: sst1, bt1, bt2

  s=0.
  ss1 = 0.0
  ss2 = 0.0
  sst1 = 0.0
  do i=n-1,1,-1
    a=x(i+1)-x(i)           ! dt
!------for total column integration--------
    bt1 = y1(i)
    bt2 = y1(i+1)
    sst1 = sst1 + a*(bt1+bt2)/(2*3600*velocity_t)
    rho0  = 1-exp(-sst1)
!------for lower layer integration-------    
    et1=y2(i)
    et2=y2(i+1)
    qw1=y3(i)
    qw2=y3(i+1)
    if (z(i) == 0.0) then
       et1 = 0.0
       qw1 = 0.0
    end if
    if (z(i+1) == 0.0) then
       et2 = 0.0
       qw2 = 0.0
    end if
    ewa1=(et1+qw1+et2+qw2)/(2*3600*velocity_t)
    ss1=ss1+a*ewa1
    ew1=et1+qw1*rho0
    ew2=et2+qw2*rho0
    ewa2=(ew1+ew2)/(2*3600*velocity_t)
    ss2=ss2+exp(ss1)*ewa2*a
  end do
  s = exp(-ss1)*ss2
end subroutine
!-----------------------------------------------------------
subroutine daily_qw(w1,q1,qw1,mask1)
  use param
  implicit none
  real(8), dimension(bigx, bigy, days*24/velocity_t,2), intent(in) :: w1    !m/s
  real(8), dimension(bigx, bigy, days*24/velocity_t) ,  intent(in) :: q1    !kg/kg
  real(8), dimension(bigx, bigy, days),                 intent(out):: qw1   !convert to mm/3h
  integer, dimension(bigx, bigy, days),               intent(out):: mask1
  integer  :: i, j, k, istep1, istep2, it, flag1
  real(8)     :: sum1, sum2

  do k = 1, days
     istep1 = (k-1)*24/velocity_t+1
     istep2 = k*24/velocity_t
     do j = 1, bigy
        do i = 1, bigx
           sum1 = 0.0
           flag1 = 0
           do it = istep1, istep2
              if (isnan(q1(i,j,it))) then
                  mask1(i,j,k) = 0
                  qw1(i,j,k)   = NAN_value
                  flag1 = 1
                  exit
              else
                sum1 = sum1 + (q1(i,j,it)*w1(i,j,it,1)+q1(i,j,it)*w1(i,j,it,2))*3600*3.0  ! *1000/1000 (m->mm,rho_water)
              end if
           end do
           if (flag1 == 0) then
              if (sum1 >=0) then
                 mask1(i,j,k) = 1
              else
                 mask1(i,j,k) = 2
              end if
              qw1(i,j,k) = sum1/(24/velocity_t)
           end if
        end do
     end do
  end do
  print *, 'wq converted to daily averages with mm/3h unite.'
end subroutine 
!--------------------------------------------------------------------------------
subroutine DOM_SPTNQU(ew, icount, ir, rho2, rho1l)
  use param
  use typedef
  implicit none
  real(8), dimension(max_tracing, 7), intent(in) :: ew
  type(start_end), dimension(max_tracing), intent(in) :: icount
  integer, intent(in)  :: ir     !reg_change
  real(8), dimension(region_num),intent(out)    :: rho2, rho1l
  integer   :: i,j,k, ein, ind, L1, L2
  real(8) ::  s1, s0
  real(8), allocatable :: various_reg(:,:), various_reg0(:,:)
  real(8), dimension(region_num) :: sum_region, sum_region0
  
       allocate(various_reg(ir,3))
       allocate(various_reg0(ir,3))
       various_reg = 0.0
       various_reg0 =0.0
       do ein=1,ir
          L1=icount(1)%start_l
          L2=icount(ein)%end_l
          call SPTNQ_U(ew(L1:L2,1),ew(L1:L2,2),ew(L1:L2,3),ew(L1:L2,5),&
                       L2-L1+1,s1,s0)
          if (ein == 1) then
             various_reg(ein,1)=0.0
             various_reg0(ein,1)=0.0
          else
             various_reg(ein,1)=various_reg(ein-1,2)
             various_reg0(ein,1)=various_reg0(ein-1,2)
          end if
          various_reg(ein,2) = s1
          various_reg(ein,3) = various_reg(ein,2)-various_reg(ein,1)
          various_reg0(ein,2) = s0
          various_reg0(ein,3) = various_reg0(ein,2)-various_reg0(ein,1)
       end do
      ! compute the contribution of each region
      sum_region=0.0
      sum_region0=0.0
      do ein=1,ir
        ind=icount(ein)%renum
        sum_region(ind)=sum_region(ind)+max(various_reg(ein,3),0.0)
        sum_region0(ind)=sum_region0(ind)+max(various_reg0(ein,3),0.0)
      end do
      rho2 = sum_region
      rho1l  = sum_region0

      deallocate(various_reg)
      deallocate(various_reg0)

end subroutine
!-----------------------------------------------------------------
subroutine DOM_SPTNQL(ew, icount, ir, rho2, rho1l)
  use param
  use typedef
  implicit none
  real(8), dimension(max_tracing, 8), intent(in) :: ew
  type(start_end), dimension(max_tracing), intent(in) :: icount
  integer, intent(in)  :: ir
  real(8), dimension(region_num),intent(out)    :: rho2, rho1l
  integer   :: i,j,k, ein, ind, L1, L2
  real(8)              ::  s1, s0
  real(8), allocatable :: various_reg(:,:), various_reg0(:,:)
  real(8), dimension(region_num) :: sum_region, sum_region0

      if (icount(1)%start_l == 0) return
      allocate(various_reg(ir,3))
      allocate(various_reg0(ir,3))
      do ein=1,ir
          L1=icount(1)%start_l
          L2=icount(ein)%end_l
          call SPTNQ_L(ew(L1:L2,1),ew(L1:L2,2),ew(L1:L2,3),ew(L1:L2,4),ew(L1:L2,6),&
                       L2-L1+1,s1,s0)
          if (ein == 1) then
             various_reg(ein,1)=0.0
             various_reg0(ein,1)=0.0
          else
             various_reg(ein,1)=various_reg(ein-1,2)
             various_reg0(ein,1)=various_reg0(ein-1,2)
          end if
          various_reg(ein,2)=s1
          various_reg(ein,3)=various_reg(ein,2)-various_reg(ein,1)
          various_reg0(ein,2)=s0
          various_reg0(ein,3)=various_reg0(ein,2)-various_reg0(ein,1)
      end do

      ! compute the contribution of each region
      sum_region=0.0
      sum_region0=0.0
      do ein=1,ir
        ind=icount(ein)%renum
        sum_region(ind)=sum_region(ind)+max(various_reg(ein,3),0.0)
        sum_region0(ind)=sum_region0(ind)+max(various_reg0(ein,3),0.0)
      end do
      rho2 = sum_region
      rho1l = sum_region0
      deallocate(various_reg)
      deallocate(various_reg0)
   
end subroutine
!------------------------------------------------------------------------
subroutine L2_ave(pp, ru, rl, ra)
  use param
  real(8), dimension(domsize, 2), intent(in) :: pp
  real(8), dimension(domsize, region_num),  intent(in) :: ru, rl
  real(8), dimension(region_num),           intent(out):: ra
  integer                                 :: k, ir
  real(8)                                    :: ppsum, sumpp
  real(8), dimension(region_num)             :: sumrr
   
  sumpp = 0.0
  sumrr = 0.0
  do k = 1, domsize
     ppsum = pp(k,1)+pp(k,2)
     sumpp = sumpp + ppsum 
     do ir = 1, region_num
        sumrr(ir) = sumrr(ir) + pp(k,1)*ru(k,ir)+pp(k,2)*rl(k,ir)
     end do
  end do

  do ir = 1,region_num
     ra(ir) = sumrr(ir)/sumpp
  end do
end subroutine


ENDc

chmod 770 drm2L_grid.f90
gfortran -g -fcheck=all -Wall -mcmodel=medium drm2L_grid.f90
./a.out
rm drm2L_grid.f90

@ inyr++
@ inyr++
@ inyr++
end

