#!/bin/csh

@ y1 = 1999
@ y2 = 2001

@ inyr = $y1
while( $inyr <= $y2)

if (($inyr % 4) == 0) then
   set days = 381
else
   set days = 380
endif

cat >! drm2L_75.f90 << ENDc
module param
! this module is used to define the dimensions of large domain and the time duration of available data
  integer,parameter :: bigx=123          !number of grids in x direction
  integer,parameter :: bigy=90           !number of grids in y direction
  integer,parameter :: days=${days}          !number of days for input fields E, W, U, V and P
  real(8),parameter :: dx=75000.00       !unit(m), the grid-spacing of each grid in x direction
  real(8),parameter :: dy=75000.00       !unit(m), the grid-spacing of each grid in y direction
  integer,parameter :: day1=16           !starting date for drm simulation
  integer,parameter :: day2=${days}          !ending date for drm simulation
  integer,parameter :: domsize=194       !number of grids within the target region
  integer,parameter :: region_num=8      !number of subregions in map file
  integer,parameter :: dt=1800           !unit(second), time interval of each back-tracing step
  integer,parameter :: max_tracing=1000  ! maximum times of back trajactory(>48*15)
  integer,parameter :: max_iteration=100 ! maximum times of iteration
  integer,parameter :: warmup=15         !days that to be excluded at the beginning used for back-tracking (= day1-1)
  integer,parameter :: velocity_t=3      !unit(hr), temporal resolution for velocity U and V fiels
  integer,parameter :: valve1=6          !times of tracing per unit velocity time interval (=3600*3/1800)
  integer,parameter :: valve2=48         !times of tracing per day (=3600*24/1800)
  integer,parameter :: NAN_value=-9.99E8
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
  real(8), dimension(bigx,bigy,days*24/velocity_t,2)   :: grid_uv_u     ! 1-U3U(m/s),2-V3U(m/s); U&V for upper slab
  real(8), dimension(bigx,bigy,days*24/velocity_t,2)   :: grid_uv_l     ! 1-U3L(m/s),2-V3L(m/s); U&V for lower slab
  real(8), dimension(bigx,bigy,days*24/velocity_t,2)   :: grid_w700     ! 1-W3U(m/s),2-W3D(m/s)
  real(8), dimension(bigx,bigy,days*24/velocity_t  )   :: grid_q700     ! q at 700mb(kg/kg)
  integer, dimension(day1:day2)                     :: ithday        ! read from 1L-DRM
  real(8), dimension(day1:day2)                        :: rho_1l        ! recycling ratio read from 1L-DRM
  integer, dimension(bigx,bigy)                     :: mask_700      !
  integer, dimension(bigx,bigy)                     :: map
  real(8), dimension(bigx,bigy,days)                   :: grid_qw1d     
  integer, dimension(bigx,bigy,days)                :: qw_mask       ! 0-UL grids, 1-upward qw, 2-downward qw

  real(8)                                              :: EW_RATIO=0.,QWU_RATIO=0.,QWD_RATIO=0.,TPP=0.,EWL_RATIO=0.0
  real(8)                                              :: TU_1=0.,TV_1=0.,TU_2=0.,TV_2=0.
  real(8), dimension(domsize,max_tracing,7)            :: EW_vector_u = NAN_value !matrix used to save the E/W data of each grid on specific day and the region it belongs as well
  real(8), dimension(domsize,max_tracing,8)            :: EW_vector_l = NAN_value
  real(8)                                              :: E_W=0   !variable to save integrated EW_vector of each grid on a specific day
  real(8), dimension(domsize,day1:day2)                :: rrcalc  !matrix to save recycling ratio of each grid on specific day
  integer                                           :: i,j,k,day,time,it,intx,inty,intx1,inty1,flag1,flag2,rn,t
  real(8)                                              :: oldppx,oldppy,diffx,diffy
  real(8), dimension(max_tracing)                      :: ppxu,ppyu,ppxl,ppyl  !variable used to record the position of tracing
  integer                                           :: daycount,hourcount,count1,count2
  real(8)                                              :: temprr,sumrr,sumpp,s, grid_rr
  real(8), dimension(day1:day2)                        :: rr

  type(start_end), dimension(domsize,max_tracing)   :: reg_count_u
  type(start_end), dimension(domsize,max_tracing)   :: reg_count_l
  real(8),            dimension(day1:day2,region_num)  :: daily_rru, daily_rrl, daily_rra  !record daily recycling ratio of different regions
  integer,         dimension(domsize,2)             :: domainij
  real(8),            dimension(domsize,2)             :: precip_grid=0  ! 1-upper, 2-lower
  integer                                           :: reg_change,old_reg,new_reg,zeroi,ein,xx,ind,L1,L2
  integer,         dimension(domsize)               :: reg_change_u, reg_change_l
  real(8)                                              :: s1,pral
  real(8),            dimension(region_num)            :: rho_ave, rho_upper, rho_lower
  real(8),            dimension(domsize,region_num)    :: rho_upper_k, rho_lower_k

! build files to output data from running the model
  call buildfile

! load data into array grid_days and grid_hours
  call openfile(grid_1d,grid_pw,grid_uv_u,grid_uv_l,grid_w700,grid_q700)  !including fixing data <0
  call readmap(map,mask_700)
  call monsoon(domainij)
  call daily_qw(grid_w700,grid_q700,grid_qw1d,qw_mask)   !read in Q and W at 700mb level (interface between lower and upper slabs)

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!############################################################################################
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

! start the daily loop during the chosen duration from day1 to day2
  print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  print *,'!!!  Daily Loop Starts  !!!'
  print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  !daily : do day=day1,day2
  daily : do day=day1,day2

    print *,'day:', day

    !------partition precipitation into 2 layers-------
    do xx=1,domsize 
      precip_grid(xx,1)=0.
      precip_grid(xx,2)=0.
      i=domainij(xx,1)
      j=domainij(xx,2)
!      TPP=grid_1d(i,j,day,2)
      reg_change_u(xx) = 0.
      reg_change_l(xx) = 0.
      if (mask_700(i,j) == 0) then
         precip_grid(xx,1) = grid_pw(i,j,day,1)
      else
         precip_grid(xx,1) = grid_pw(i,j,day,1)
         precip_grid(xx,2) = grid_pw(i,j,day,2)
      end if
    end do  
    print *, 'precip partitioned into two layers.'
!=======================================================
!-------(1) run DRM for the upper level-----------------
!=======================================================
    print *, 'upper level back tracking...'
    ! start the domain loop for the upper slab
    domainU : do k=1,domsize
    !domainU : do k=1,10
        i=domainij(k,1)
        j=domainij(k,2)
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
     
        ! initialize vars to be zero
        do zeroi=1,max_tracing
          do ein=1,7
            EW_vector_u(k,zeroi,ein)=0
          end do
          reg_count_u(k,zeroi)%start_l=0
          reg_count_u(k,zeroi)%end_l=0
          reg_count_u(k,zeroi)%renum=0
        end do

!====================================================================================================
        ! start the time loop
        timeloopU : do time=1,max_tracing-1
          !print *, ppxu(time)
          ! tell ppx and ppy are within large domain or not-----------------------------------
          if ((ppxu(time)<=0).or.(ppxu(time)>=(dx*bigx))) then
          !  print *,'Beyond boundary during timeloop when time is:',time
            exit timeloopU
          else if ((ppyu(time)<=0).or.(ppyu(time)>=(dy*bigy))) then
          !  print *,'Beyond boundary during timeloop when time is:',time
            exit timeloopU
          end if
          !-----------------------------------------------------------------------------------

          !read the Evaporation,precipitable water and wind velocity by subroutine 'finduvEW_U' for the upper slab
          call finduvEW_U(ppxu(time),ppyu(time),hourcount,daycount,grid_pw(:,:,daycount,:),grid_1d(:,:,daycount,1),&
              grid_uv_u(:,:,hourcount,1),grid_uv_u(:,:,hourcount,2),grid_qw1d(:,:,daycount),mask_700,&
              qw_mask(:,:,daycount),EW_RATIO,QWU_RATIO,TU_1,TV_1,intx1,inty1,flag1)
          if (flag1==5) then
            exit timeloopU
          else
            EW_vector_u(k,time,1) = (time-1)*dt
            EW_vector_u(k,time,2) = EW_RATIO
            EW_vector_u(k,time,3) = QWU_RATIO
            EW_vector_u(k,time,4) = map(intx1,inty1)
            EW_vector_u(k,time,5) = real(flag1)       ! if flag1 == 0, then don't need to be multiplied by rho_lower
            EW_vector_u(k,time,6) = intx1   !** try just use i
            EW_vector_u(k,time,7) = inty1   !** and j instead
            !TET in unit of mm/3h thus TET*8 and EW_vector*8 has the unit of mm/day and day-1 individually
          end if

          !write(*,*) EW_RATIO, QWU_RATIO
 
          if (time .eq. 1) then
            reg_count_u(k,1)%start_l=1
            reg_count_u(k,1)%renum=map(i,j)
            reg_change=1
          else
            old_reg=EW_vector_u(k,time-1,4)
            new_reg=EW_vector_u(k,time,4)
            if (old_reg .ne. new_reg) then
              reg_change=reg_change+1
              reg_count_u(k,reg_change-1)%end_l=time
              reg_count_u(k,reg_change)%start_l=time
              reg_count_u(k,reg_change)%renum=new_reg      
            end if
          end if
  
          ! find the next point using Iterative Technique--------------------------------------
          call finduv(ppxu(time),ppyu(time),hourcount,grid_uv_u(:,:,hourcount,1),& 
               grid_uv_u(:,:,hourcount,2),TU_2,TV_2)
          oldppx=ppxu(time)-(TU_1+TU_2)*dt/2
          oldppy=ppyu(time)-(TV_1+TV_2)*dt/2
          if ((oldppx<=0).or.(oldppx>=(dx*bigx))) exit timeloopU
          if ((oldppy<=0).or.(oldppy>=(dy*bigy))) exit timeloopU
         ! print *, time, ppxu(time), dt, TU_1, TU_2
!================================================================================================
          iterationU : do it=1,max_iteration
            !tell if the point still within the large domain
            if ((oldppx<=0).or.(oldppx>=(dx*bigx))) then
          !    print *,'Beyond boundary during iteration when it is',it,oldppx
              exit iterationU
            else if ((oldppy<=0).or.(oldppy>=(dy*bigy))) then
          !    print *,'Beyond boundary during iteration when it is',it,oldppy
              exit iterationU
            end if           

            call finduv(oldppx,oldppy,hourcount,grid_uv_u(:,:,hourcount,1),&
                 grid_uv_u(:,:,hourcount,2),TU_2,TV_2,intx,inty,flag2)
          if (flag2==1) then
             exit iterationU
          else
            ppxu(time+1)=ppxu(time)-(TU_1+TU_2)*dt/2
            ppyu(time+1)=ppyu(time)-(TV_1+TV_2)*dt/2
            
            !print *, flag2, TV_1, TV_2, ppyu(time)
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
           reg_count_u(k,reg_change)%end_l=count2
        end if
        reg_change_u(k) = reg_change
 
   write(400) real(ppxu, kind(1e0))
   write(400) real(ppyu, kind(1e0))
   end do domainU

   print *,'Doing integral of the upper layer...'
   call DOM_SPTNQU(EW_vector_u,precip_grid(:,1),reg_count_u,reg_change_u, rho_upper,rho_upper_k)
   ! print *, k, shape(rho_upper_k)
   ! print *, 'k, rho_upper_k:', k, rho_upper_k
   write(500) real(rho_upper_k, kind(1e0))
!   write(600) real(EW_vector_u, kind(1e0))

!=================================================================
!---------(2) run DRM for the lower level-------------------------
!=================================================================
    print *, 'lower leve back tracking...'
    ! start the domain loop
    domainL : do k=1,domsize
    !domainL : do k=1,10
        i=domainij(k,1)
        j=domainij(k,2)
        ppxl = 0.0
        ppyl = 0.0
        ppxl(1)=(i-1)*dx+dx/2
        ppyl(1)=(j-1)*dy+dy/2

        count1=0
        count2=0
        reg_change=0
        hourcount=day*24/velocity_t   !starting temporal position in velocity files
        daycount=day                  !starting temporal position in other files

        ! initialized vars to be zero
        do zeroi=1,max_tracing
          do ein=1,8
            EW_vector_l(k,zeroi,ein)=0
          end do
          reg_count_l(k,zeroi)%start_l=0
          reg_count_l(k,zeroi)%end_l=0
          reg_count_l(k,zeroi)%renum=0
        end do
        reg_change = 0

!====================================================================================================
        ! start the time loop
        timeloopL : do time=1,max_tracing-1

          ! tell ppx and ppy are within large domain or
          ! not-----------------------------------
          if ((ppxl(time)<=0).or.(ppxl(time)>=(dx*bigx))) then
           ! print *,'Beyond boundary during timeloop when time is:',ppxl(time)
            exit timeloopL
          else if ((ppyl(time)<=0).or.(ppyl(time)>=(dy*bigy))) then
           ! print *,'Beyond boundary during timeloop when time is:',ppyl(time)
            exit timeloopL
          end if
          !-----------------------------------------------------------------------------------

          !read the Evaporation,precipitable water and wind velocity by
          !subroutine 'finduvEW_L' for the lower slab
          call finduvEW_L(ppxl(time),ppyl(time),hourcount,daycount,grid_pw(:,:,daycount,:),grid_1d(:,:,daycount,1),&
              grid_uv_l(:,:,hourcount,1),grid_uv_l(:,:,hourcount,2),grid_qw1d(:,:,daycount), mask_700,&
              qw_mask(:,:,daycount),EW_RATIO,EWL_RATIO,QWD_RATIO,TU_1,TV_1,intx1,inty1,flag1)
          if (flag1==5) then
            if (reg_change > 0) then
              reg_count_l(k,reg_change)%end_l=time
            else
              reg_count_l(k,1)%start_l = 0
              reg_count_l(k,1)%end_l   = 0
              reg_count_l(k,1)%renum   = 0
            end if
            exit timeloopL
          else
            EW_vector_l(k,time,1)=(time-1)*dt
            EW_vector_l(k,time,2)=EW_RATIO
            EW_vector_l(k,time,3)=EWL_RATIO
            EW_vector_l(k,time,4)=QWD_RATIO
            EW_vector_l(k,time,5)=map(intx1,inty1)
            EW_vector_l(k,time,6)=real(flag1)       ! if flag1 == 0, then don't need to be multiplied by rho_upper
            EW_vector_l(k,time,7)=intx1
            EW_vector_l(k,time,8)=inty1
          end if

          if (time .eq. 1) then
            reg_count_l(k,1)%start_l=1
            reg_count_l(k,1)%renum=map(i,j)
            reg_change=1
          else
            old_reg=EW_vector_l(k,time-1,5)
            new_reg=EW_vector_l(k,time,5)
            if (old_reg .ne. new_reg) then
              reg_change=reg_change+1
              reg_count_l(k,reg_change-1)%end_l=time
              reg_count_l(k,reg_change)%start_l=time
              reg_count_l(k,reg_change)%renum=new_reg
            end if
          end if

! find the next point using Iterative Technique--------------------------------------
          call finduv(ppxl(time),ppyl(time),hourcount,grid_uv_l(:,:,hourcount,1),&
               grid_uv_l(:,:,hourcount,2),TU_2,TV_2)
          oldppx=ppxl(time)-(TU_1+TU_2)*dt/2
          oldppy=ppyl(time)-(TV_1+TV_2)*dt/2
          if ((oldppx<=0).or.(oldppx>=(dx*bigx))) exit timeloopL
          if ((oldppy<=0).or.(oldppy>=(dy*bigy))) exit timeloopL
!================================================================================================
          iterationL : do it=1,max_iteration
            !tell if the point still within the large domain
            if ((oldppx<=0).or.(oldppx>=(dx*bigx))) then
              print *,'Beyond boundary during iteration when it is',it,oldppx
              !print *,ppxl(time),TU_1, TU_2
              exit iterationL
            else if ((oldppy<=0).or.(oldppy>=(dy*bigy))) then
              !print *,'Beyond boundary during iteration when it is',it,oldppy
              exit iterationL
            end if

            call finduv(oldppx,oldppy,hourcount,grid_uv_l(:,:,hourcount,1),&
                 grid_uv_l(:,:,hourcount,2),TU_2,TV_2,intx,inty,flag2)
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
          reg_count_l(k,reg_change)%end_l=count2
        end if
        reg_change_l(k) = reg_change
    write(401) real(ppxl, kind(1e0))
    write(401) real(ppyl, kind(1e0))

    end do domainL

    print *,'Doing integral of the lower layer...'
    call DOM_SPTNQL(EW_vector_l,precip_grid(:,2),reg_count_l,reg_change_l, rho_lower, rho_lower_k)
    !print *, 'k, rho_lower_k:', k, rho_lower_k
    write(501) real(rho_lower_k, kind(1e0))
!    write(601) real(EW_vector_l, kind(1e0))

!===================================================================================================
   
    call L2_ave(precip_grid,rho_upper_k, rho_lower_k, rho_ave)

!    print *, 'it, rho_upper, rho_lower:', it, rho_upper, rho_lower
!    print *, 'it, rho_upper_k         :', it, rho_upper_k
!    print *, 'it, rho_lower_k         :', it, rho_lower_k

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    daily_rru(day,:)=rho_upper
    daily_rrl(day,:)=rho_lower
    daily_rra(day,:)=rho_ave
    write(300,'(I3,8(5X,F6.4))') day,daily_rra(day,:)
    write(301,'(I3,8(5X,F6.4))') day,daily_rru(day,:)
    write(302,'(I3,8(5X,F6.4))') day,daily_rrl(day,:)
    !write(300,'(I3,5X,F5.3,5X,F5.3,5X,F5.3,5X,F5.3,5X,F5.3)') day,daily_rra(day,:)
    !write(301,'(I3,5X,F5.3,5X,F5.3,5X,F5.3,5X,F5.3,5X,F5.3)') day,daily_rru(day,:)
    !write(302,'(I3,5X,F5.3,5X,F5.3,5X,F5.3,5X,F5.3,5X,F5.3)') day,daily_rrl(day,:)
    
  end do daily

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!##################################################################################################
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
  close(300)
  close(301)
  close(302)
  close(400)
  close(401)
  close(500)
  close(501)
!  close(600)
!  close(601)
end program Recycling

!===========================================================================
! subroutine to build files for output
!---------------------------------------------------------------------------
subroutine buildfile
  use param
  implicit none
  integer,parameter :: filen_out=3
  character(len=5) :: items(filen_out)=(/'DAYRA','DAYRU','DAYRL'/) 
  character(len=5) :: items2(2)=(/'trj_U','trj_L'/) 
  character(len=5) :: items3(2)=(/'rat_U','rat_L'/) 
  character(len=5) :: items4(2)=(/'EWv_U','EWv_L'/) 
  character(len=5) :: yyyy='${inyr}_'
  character(len=70) :: chfile1,chfile2,path,chfile3,chfile4,chfile5
  character(len=70) :: chfile6,chfile7,chfile8,chfile9

  !write(path,'(A29,I4,A1)') '/data/francina/b/hhu18/2LDRM/',ithyear,'/'
  !write(path,'(A46)') '/data/keeling/a/seunguk2/a/2LDRM/interp/2LDRM/'
  write(path,'(A40)') '/data/keeling/a/seunguk2/c/2LDRM/pydlmb/'
  write(chfile1,'(A40,A5,A5,A7)') path,yyyy,items(1),'_8M.txt'
  open(unit=300,file=chfile1) 
    
  write(chfile2,'(A40,A5,A5,A7)') path,yyyy,items(2),'_8M.txt'
  open(unit=301,file=chfile2)
  
  write(chfile3,'(A40,A5,A5,A7)') path,yyyy,items(3),'_8M.txt'
  open(unit=302,file=chfile3)

  write(chfile4,'(A40,A5,A5,A7)') path,yyyy,items2(1),'_8M.dat'
  open(unit=400,file=chfile4,status='unknown',form='unformatted')

  write(chfile5,'(A40,A5,A5,A7)') path,yyyy,items2(2),'_8M.dat'
  open(unit=401,file=chfile5,status='unknown',form='unformatted')

  write(chfile6,'(A40,A5,A5,A7)') path,yyyy,items3(1),'_8M.dat'
  open(unit=500,file=chfile6,status='unknown',form='unformatted')

  write(chfile7,'(A40,A5,A5,A7)') path,yyyy,items3(2),'_8M.dat'
  open(unit=501,file=chfile7,status='unknown',form='unformatted')

!  write(chfile8,'(A40,A5,A5,A7)') path,yyyy,items4(1),'_8M.dat'
!  open(unit=600,file=chfile8,status='unknown',form='unformatted')

!  write(chfile9,'(A40,A5,A5,A7)') path,yyyy,items4(2),'_8M.dat'
!  open(unit=601,file=chfile9,status='unknown',form='unformatted')
end subroutine
!==========================================================================
! subroutine to read in the i and j indices for the target region in the domain
!__________________________________________________________________________
subroutine monsoon(domain)
  use param
  implicit none
  integer,dimension(domsize,2),intent(out) :: domain
  real,dimension(domsize,2) :: domain2
  integer :: ii,error
  character(len=70) :: file_name, file_name2
  file_name='/data/keeling/a/seunguk2/a/2LDRM/pyd075/MW_il.dat'
  file_name2='/data/keeling/a/seunguk2/a/2LDRM/pyd075/MW_jl.dat'
  open(14,file=file_name,status='old',form='unformatted',access='direct',recl=4)
  open(15,file=file_name2,status='old',form='unformatted',access='direct',recl=4)
  do ii = 1,domsize
     read(14,rec=ii) domain2(ii,1)
     read(15,rec=ii) domain2(ii,2)
!     write(*,*) domain2(ii,1), domain2(ii,2)
  enddo
  domain = int(domain2)
!  domain(:,2) = bigy+1-domain(:,2)
!  ii=1
!  do while(.true.)
!    read(14,'(2I4)',iostat=error)domain(ii,1),domain(ii,2)
!    if(error/=0) exit
!    ii=ii+1
!  end do
  write(*,*) 'Monsoon region has been loaded.'
  close(14)
  close(15)
end subroutine
!===========================================================================
! subroutine to read in the input fields E, PP
! PW, U and V for upper and lower slabs, W and Q at 700mb level
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
  real(8), dimension(bigx,bigy,days*24) :: xqi
  real(8), dimension(bigx,bigy,days*24/velocity_t),  intent(out) :: xq
  real(8), dimension(bigx,bigy,days*24/velocity_t) :: xwi
  character :: n_c
  character(len=70) :: chfile, path, path2
  character(len=90) :: whpath
  character(len=2) :: items(2)=(/'ET','ET'/)
  character(len=3) :: items2(8)=(/'QWU','QWD','U3U','V3U','U3D','V3D','W3U','W3D'/)
  character(len=4) :: items3(1)=(/'q700'/)
  character(len=5) :: yyyy='_${inyr}'

  write(path,'(A40)') '/data/keeling/a/seunguk2/b/2LDRM/pydlmb/'
!  write(path2,'(A29,I4,A1)') '/data/francina/b/hhu18/2LDRM/',ithyear,'/'

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
  x1d(:,:,:,1) = -x1d(:,:,:,2)
  do i = 1, bigx
    do j = 1, bigy
      do k = 1, days
        x1d(i,j,k,1) = max(x1d(i,j,k,1),0.0)
      end do
    end do
  end do
  
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
            print *, whpath, k
            !read(10,rec=k) ((xw(i,j,k,l-6),i=1,bigx),j=1,bigy)
            read(10,rec=k) ((xw(i,j,k,1),i=1,bigx),j=1,bigy)
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
!==========================================================================
! subroutine to read in masks for each subregion and that for 700 mbdatmask700 is used to tell if a grid only has the upper slab (elevations > 700mb) or both slabs
!__________________________________________________________________________
subroutine readmap(regions, mask700)
  use param
  implicit none
  integer,dimension(bigx,bigy),intent(out) :: regions
  integer,dimension(bigx,bigy),intent(out) :: mask700
  real,dimension(bigx,bigy) :: regions2
  real,dimension(bigx,bigy) :: mask7002
  character(len=70) :: path, file_name, file_name2
  integer :: i,j
      !/data/keeling/a/seunguk2/a/2LDRM/interp/2LDRM/ERA5_6REG.txt

  write(path,'(A40)') '/data/keeling/a/seunguk2/a/2LDRM/pyd075/'
  write(file_name,'(A40,A13)') path, 'ERA5_8rLC.dat'
  write(file_name2,'(A40,A14)') path, '700_masklc.dat'
  open(113,file=file_name2,status='old',form='unformatted',access='direct',recl=4*bigx)
  do j = 1,bigy
    read(113,rec=j) (mask7002(i,j), i=1,bigx)
    !write(*,*) mask7002(200,j)
  enddo
  write(*,*) 'Region Mask has been loaded.'
  open(112,file=file_name,status='old',form='unformatted',access='direct',recl=4*bigx)
  do j = 1,bigy
    read(112,rec=j) (regions2(i,j), i=1,bigx)
    !write(*,*) regions2(200,j)
    ! For UW region setting,
    !regions2(:,j) = amod(regions2(:,j) + 3, 5.0) + 1
    ! For GM region setting,
    ! regions2(:,j) = amod(regions2(:,j) + 3, 5.0) + 1
  enddo
  close(112)
  write(*,*) '700mb Mask has been loaded.'
  close(113)
  regions = int(regions2)
  mask700 = int(mask7002)
end subroutine
!=========================================================================
! subroutine to find U,V,E and W for the upper slab at a certain location and time
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
  !write(*,*) 'ET', ET
  !write(*,*) 'qwu',qwu, outW
  !write(*,*) 'QW',QW(xint,yint), xint, yint
  !write(*,*) 'W2L', W2L

end subroutine
!=========================================================================
! subroutine to find U,V,E and W for the lower slab at a certain location and time
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
     !qwd  = max(QW(xint,yint),0.0)
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

  !write(*,*) 'ET ',outE
  !write(*,*) 'qwd',qwd
  !write(*,*) 'W2L',W2L
end subroutine
!==========================================================================
! subroutine to find U and V for a certain grid at a certain time
! called during back-tracking
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
  !print *, xint, yint, xp, yp, outu, outv, flag
  !if (outu==NAN_value) then
  !  flag=1
  !else if (outv==NAN_value) then
  !  flag=1
  !end if
end subroutine
!=============================================================
! subroutine to do the integral along back-trajectory in upper slab
!_____________________________________________________________
subroutine SPTNQ_U(x,y1,y2,z,n,s)
  use param
  implicit none
  integer,intent(in):: n
  real(8),dimension(n),intent(in)::x,y1,y2,z
  real(8),intent(out) :: s
  integer:: i,j,k
  real(8) :: a,b1,b2, ewa1,ewa2, ew1,ew2
  real(8) :: ss1, ss2
  real(8) :: bt1, bt2, sst1
  real(8) :: rho

  s=0.
  ss1 = 0.0
  ss2 = 0.0
  sst1 = 0.0
  do i=n-1,1,-1
    rho = 0.0
    a=x(i+1)-x(i)           ! dt
!-----for total column integration------
    bt1=y1(i)
    bt2=y1(i+1)
    sst1  = sst1 + a*(bt1+bt2)/(2*3600*velocity_t)
    rho   = 1-exp(-sst1)
!-----for upper layer integration-------
    b1 =y2(i)
    b2 =y2(i+1)
    ewa1=(b1+b2)/(2*3600*velocity_t)
    ss1=ss1+a*ewa1
    ew1=b1*rho
    ew2=b2*rho
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
  !print *, s, ss1, ss2, sst1
end subroutine
!===========================================================
! subroutine to do the integral along back-trajectory in lower slab
!-----------------------------------------------------------
subroutine SPTNQ_L(x,y1,y2,y3,z,n,s)
  use param
  implicit none
  integer,intent(in):: n
  real(8),dimension(n),intent(in)::x,y1,y2,y3,z
  real(8),intent(out) :: s
  integer:: i,j,k
  real(8) :: a,et1,et2,qw1,qw2
  real(8) :: ewa1,ewa2, ew1,ew2
  real(8) :: ss1, ss2
  real(8) :: sst1, bt1, bt2, rho

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
    rho  = 1-exp(-sst1)
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
    ew1=et1+qw1*rho
    ew2=et2+qw2*rho
    ewa2=(ew1+ew2)/(2*3600*velocity_t)
    ss2=ss2+exp(ss1)*ewa2*a
  end do
  s = exp(-ss1)*ss2
end subroutine
!===========================================================
! subroutine to calculate daily values of q*w
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
              !write(*,*) q1(i,j,it), w1(i,j,it,1), w1(i,j,it,2)
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
  print *, 'wq converted to daily averages with mm/3h unit.'
end subroutine
!================================================================================
! subroutine to do the integral for each grid in target region in upper level slab
!--------------------------------------------------------------------------------
subroutine DOM_SPTNQU(ew, pp, icount, ir, rho2, rho2_k)
  use param
  use typedef
  implicit none
  real(8),            dimension(domsize, max_tracing, 7), intent(in) :: ew
  real(8),            dimension(domsize),                 intent(in) :: pp
  type(start_end), dimension(domsize,max_tracing),     intent(in) :: icount
  integer,         dimension(domsize),                 intent(in) :: ir     !reg_change
  real(8),            dimension(region_num),              intent(out):: rho2
  real(8),            dimension(domsize,region_num),      intent(out):: rho2_k
  integer                     :: i,j,k, ein, ind, L1, L2
  real(8)                        :: sumpp, s1, rho
  real(8), dimension(region_num) :: sumrr
  real(8), allocatable           :: various_reg(:,:)
  real(8), dimension(region_num) :: sum_region
  
print *, 'in DOM_SPTNQU:'
   sumrr = 0.0
   sumpp = 0.0
   rho2_k = 0.0
   do k =1, domsize
       allocate(various_reg(ir(k),3))
       various_reg = 0.0
       do ein=1,ir(k)
          L1=icount(k,1)%start_l
          L2=icount(k,ein)%end_l
          call SPTNQ_U(ew(k,L1:L2,1),ew(k,L1:L2,2),ew(k,L1:L2,3),ew(k,L1:L2,5),&
                       L2-L1+1,s1)
          if (ein == 1) then
             various_reg(ein,1)=0.0
          else
             various_reg(ein,1)=various_reg(ein-1,2)
          end if
          various_reg(ein,2) = s1
          various_reg(ein,3) = various_reg(ein,2)-various_reg(ein,1)
       end do
      ! compute the contribution of each region
      sum_region=0.0
      do ein=1,ir(k)
        ind=icount(k,ein)%renum
        sum_region(ind)=sum_region(ind)+max(various_reg(ein,3),0.0)
      end do

      rho2_k(k,:) = sum_region
      do ein = 1, region_num
         sumrr(ein)=sumrr(ein)+rho2_k(k,ein)*pp(k)
      end do
      sumpp=sumpp+pp(k)
      deallocate(various_reg)
      !print *, 'upper rho_k', k, sum_region
   end do

   do ein = 1, region_num
     rho2(ein) = sumrr(ein)/sumpp
   end do
   !print *, sum_region
   print *, 'upper level sumpp', sumpp
   print *, 'upper level sumrr', sumrr
end subroutine
!=================================================================
! subroutine to do the integral for each grid in target region in lower level slab
!-----------------------------------------------------------------
subroutine DOM_SPTNQL(ew, pp, icount, ir, rho2, rho2_k)
  use param
  use typedef
  implicit none
  real(8), dimension(domsize, max_tracing, 8), intent(in) :: ew
  real(8), dimension(domsize), intent(in) :: pp
  type(start_end), dimension(domsize,max_tracing), intent(in) :: icount
  integer,dimension(domsize), intent(in)  :: ir
  real(8), dimension(region_num),intent(out)    :: rho2
  real(8), dimension(domsize,region_num),intent(out)::rho2_k
  integer   :: i,j,k, ein, ind, L1, L2
  real(8), dimension(region_num) :: sumrr
  real(8)              ::  s1, sumpp
  real(8), allocatable :: various_reg(:,:)
  real(8), dimension(region_num) :: sum_region

   sumrr = 0.0
   sumpp = 0.0
   rho2_k = 0.0
   do k =1, domsize
      if (icount(k,1)%start_l == 0) cycle
      allocate(various_reg(ir(k),3))
      do ein=1,ir(k)
          L1=icount(k,1)%start_l
          L2=icount(k,ein)%end_l
          call SPTNQ_L(ew(k,L1:L2,1),ew(k,L1:L2,2),ew(k,L1:L2,3),ew(k,L1:L2,4),ew(k,L1:L2,6),&
                       L2-L1+1,s1)
          if (ein == 1) then
             various_reg(ein,1)=0.0
          else
             various_reg(ein,1)=various_reg(ein-1,2)
          end if
          various_reg(ein,2)=s1
          various_reg(ein,3)=various_reg(ein,2)-various_reg(ein,1)
      end do

      ! compute the contribution of each region
      sum_region=0.0
      do ein=1,ir(k)
        ind=icount(k,ein)%renum
        sum_region(ind)=sum_region(ind)+max(various_reg(ein,3),0.0)
      end do
      rho2_k(k,:) = sum_region
      do ein = 1, region_num
         sumrr(ein)=sumrr(ein)+rho2_k(k,ein)*pp(k)
      end do
      sumpp=sumpp+pp(k)
      deallocate(various_reg)
      !print *, 'lower rho_k', k, sum_region
   end do
   
   do ein = 1,region_num
     rho2(ein) = sumrr(ein)/sumpp
   end do
   print *, 'lower level sumpp', sumpp
   print *, 'lower level sumrr', sumrr
end subroutine
!========================================================================
! subroutine to calculate column-averaged contribution
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
!---------------------------------------------------------------------

ENDc

chmod 770 drm2L_75.f90
gfortran -g -fcheck=all -Wall -mcmodel=medium drm2L_75.f90
./a.out
rm drm2L_75.f90

@ inyr++
end
