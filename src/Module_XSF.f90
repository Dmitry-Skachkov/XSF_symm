

 Module Module_XSF
        implicit none
        real(8), allocatable        :: Pot(:)                    ! excitonic wave function
        real(8), allocatable        :: Pot1(:)
        real(8), allocatable        :: Potnew(:)
        real(8), allocatable        :: Pos(:,:)                  ! grid mesh
        real(8), allocatable        :: Pos1(:,:)
        integer                     :: Nxpot,Nypot,Nzpot         ! mesh size
        integer                     :: Npot                      ! total number of mesh points
        integer                     :: Npot1                     ! the size of the reduced WF
        character(30)               :: name 
        real(8)                     :: Point(3)
        real(8)                     :: Period(3,3)               ! lattice vectors
        integer                     :: Nat                       ! number of atoms
        character(2),allocatable    :: Atoms(:)                  ! atom names
        real(8)     ,allocatable    :: Coord(:,:)                ! positions of the atoms
        real(8)                     :: dx(3),dy(3),dz(3)
        logical                     :: lPot
        real(8)                     :: hole(3)                   ! position of the hole
        real(8)                     :: Rx0                       ! radius of Frenkel exciton
 Contains



  subroutine read_xsf
    integer            :: i,j,k
    integer            :: jj
    print *,'read_xsf: ',name
    lPot = .true.
    open(unit=2,file=trim(adjustl(name)))
    do i=1,46                                   ! xsf created by Yambo has 46 comment lines
     read(2,*)
    enddo
    read(2,*)
    read(2,*) 
    read(2,*) Period(1:3,1)
    read(2,*) Period(1:3,2)
    read(2,*) Period(1:3,3)
    read(2,*) 
    read(2,*) Nat
    allocate(Coord(3,Nat))
    allocate(Atoms(Nat))
    do i=1,Nat
     read(2,*) Atoms(i),Coord(1:3,i)
 !    print 11,i,Atoms(i),Coord(1:3,i)
    enddo
    hole(1:3) = Coord(1:3,1)                      ! first atom = hole
    print *,'hole=',hole
    if(lPot) then
     read(2,*)
     read(2,*)
     read(2,*)
     read(2,*) Nxpot,Nypot,Nzpot
     Npot = Nxpot*Nypot*Nzpot
     Npot1 = Npot
     read(2,*) Point
     Point = -hole                      ! make hole = (0,0,0)
     read(2,*) Period(1:3,1)
     read(2,*) Period(1:3,2)
     read(2,*) Period(1:3,3)
     dx(1:3) = Period(1:3,1)/dfloat(Nxpot-1)
     dy(1:3) = Period(1:3,2)/dfloat(Nypot-1)
     dz(1:3) = Period(1:3,3)/dfloat(Nzpot-1)
     allocate(Pot(Npot))
     allocate(Pot1(Npot))
     allocate(Potnew(Npot))
     allocate(Pos(3,Npot))
     allocate(Pos1(3,Npot))
     print *,'Pot: allocated ',allocated(Pot)
     print *,'Pot1: allocated ',allocated(Pot1)
     print *,'Potnew: allocated ',allocated(Potnew)
     print *,'size(Potnew)=',size(Potnew)
     print *,'Pos: allocated ',allocated(Pos)
     print *,'Pos1: allocated ',allocated(Pos1)
!     read(2,6) (((Pot(i,j,k),i=1,Nxpot),j=1,Nypot),k=1,Nzpot)
     read(2,*) (((Pot(i+(j-1)*Nxpot+(k-1)*Nxpot*Nypot),i=1,Nxpot),j=1,Nypot),k=1,Nzpot)
     Pot1(:) = Pot(:)
     Potnew(:) = Pot(:) 
     print *,'Nxpot=',Nxpot
     print *,'Nypot=',Nypot
     print *,'Nzpot=',Nzpot
     print *,'Npot=',Npot
     dx(1:3) = Period(1:3,1)/dfloat(Nxpot-1)
     dy(1:3) = Period(1:3,2)/dfloat(Nypot-1)
     dz(1:3) = Period(1:3,3)/dfloat(Nzpot-1)
     print *,'dx(1:3)=',dx(1:3)
     print *,'dy(1:3)=',dy(1:3)
     print *,'dz(1:3)=',dz(1:3)
     jj = 0
     do k=1,Nzpot                                          ! calculate the real space mesh
      do j=1,Nypot
       do i=1,Nxpot
        jj = jj + 1
        if(jj > Npot) then
         print *,'read_potential: error jj > Npot'
         stop
        endif
        Pos(1:3,jj) = (i-1)*dx(1:3) + (j-1)*dy(1:3) + (k-1)*dz(1:3) + Point(1:3)
        Pos1(1:3,jj) =  Pos(1:3,jj)
        Pos1(2,jj)   = -Pos(2,jj)                                     ! make V(x,y,z) = V(x,-y,z)
       enddo
      enddo
     enddo
    endif        ! lPot
    close(unit=2)
    print *,'read_xsf:  finish'
6   format(6E14.6)
11  format(I3,2x,A2,3F17.5) 
  end subroutine read_xsf



  subroutine reduce_potential_to_sphere
   integer              :: j
   integer              :: Npotx
   real(8), allocatable :: Posx(:,:),Potx(:) 
   integer              :: jj
   real(8)              :: Rx
   print *,'reduce potential to sphere Rx0'
   print *,'Npot=',Npot
   allocate(Posx(3,Npot))
   allocate(Potx(Npot))
   Npotx = 0
   do j=1,Npot
    Rx = dsqrt(Pos(1,j)**2+Pos(2,j)**2+Pos(3,j)**2)
    if(Rx < Rx0) then                   ! points inside the shpere Rx0  
     Npotx = Npotx + 1
     Posx(1:3,Npotx) = Pos1(1:3,j)
     Potx(Npotx)     = Pot1(j)
    endif 
   enddo
   deallocate(Pos1,Pot1)
   print *,'reduce the size from',Npot,' to',Npotx
   Npot1 = Npotx
   print *,'Npot1=',Npot1
   allocate(Pos1(3,Npot1))
   allocate(Pot1(Npot1))
   do j=1,Npot1
    Pos1(1:3,j) = Posx(1:3,j)                      ! replace Pos1, Pot1 with smaller one
    Pot1(j)     = Potx(j)
   enddo
   deallocate(Posx,Potx)
   print *,'reduce_potential_to_sphere finished'
  end subroutine reduce_potential_to_sphere



     subroutine make_symm
      use qshep
      integer, parameter   :: NR=11                       ! working parameters for qshep3
      integer, parameter   :: NW1=40 
      integer, parameter   :: NQ=17 
      integer              :: LCELL(NR,NR,NR)
      integer, allocatable :: LNEXT(:)
      real(8)              :: xyzmin(3),xyzdel(3)
      real(8)              :: rmax
      real(8), allocatable :: rsq(:)
      real(8), allocatable :: A1(:,:)
      integer              :: ier                       ! working parameters for qshep3
      real(8)              :: dx2,dy2,dz2
      integer              :: i,j,k
      real(8)              :: Px,Py,Pz
      real(8)              :: Rx
      real(8)              :: Potx
      print *
      print *,'make_symm:'
      print *,'Npot1=',Npot1
      allocate(LNEXT(Npot1))
      allocate(rsq(Npot1))
      allocate(A1(9,Npot1))
      print *,'allocated LNEXT ',allocated(LNEXT)
      print *,'allocated rsq ',allocated(rsq)
      print *,'allocated A1 ',allocated(A1)
      print *,'allocated Pos ',allocated(Pos)
      print *,'allocated Pot ',allocated(Pot)
      call qshep3(Npot1,Pos1(1,1:Npot1),Pos1(2,1:Npot1),Pos1(3,1:Npot1),Pot1(1:Npot1),NQ,NW1,NR,LCELL,LNEXT,xyzmin,xyzdel,rmax,rsq,A1,ier)
      if(ier/=0) then
       print *,'make_symm: error =',ier
       stop
      endif

      print *,'calculate symmetric WF'
      do j=1,Npot
       Rx = dsqrt(Pos(1,j)**2+Pos(2,j)**2+Pos(3,j)**2)
       if(Rx < Rx0) then                   ! points inside the sphere Rx0
        print *,'point inside the sphere j=',j
        Px = Pos(1,j)
        Py = Pos(2,j)
        Pz = Pos(3,j)
        Potx = qs3val(Px,Py,Pz,Npot1,Pos1(1,1:Npot1),Pos1(2,1:Npot1),Pos1(3,1:Npot1),Pot1(1:Npot1),NR,LCELL,LNEXT,xyzmin,xyzdel,rmax,rsq,A1)
        Potnew(j) = 0.5d0*(Pot(j) + Potx) 
        print *,'Potnew(j)=',Potnew(j)
       endif
      enddo
      print *,'make_symm: finish'
     end subroutine make_symm





      subroutine write_xsf                
       integer             :: i,j
       print *,'write_xsf'
       open(unit=11,file='new.xsf')
        write(11,1)
        write(11,2)
        write(11,3) Period(1:3,1)
        write(11,3) Period(1:3,2)
        write(11,3) Period(1:3,3)
        write(11,4) Nat
        do i=1,Nat
         write(11,5) Atoms(i),Coord(1:3,i) 
        enddo
        if(lPot) then
         write(11,8)
         write(11,9)
         write(11,10)
         write(11,11) Nxpot,Nypot,Nzpot
         Point = 0.d0
         write(11,12) Point(1:3)
         write(11,3) Period(1:3,1)
         write(11,3) Period(1:3,2)
         write(11,3) Period(1:3,3)
         write(11,13) (Potnew(j),j=1,Npot)
         write(11,14)
         write(11,15)
        endif
       close(unit=11)
       print *,'writing file new.xsf'
1      format('CRYSTAL')
2      format('PRIMVEC')
3      format(3F15.10)
4      format('PRIMCOORD'/I5,'    1')
5      format(1x,A5,3F15.10)    !,2x,3F15.10)
6      format('ANIMSTEPS',I3)
7      format('ATOMS',I3)
8      format('BEGIN_BLOCK_DATAGRID_3D')
9      format('3D_PWSCF')
10     format('BEGIN_DATAGRID_3D_UNKNOWN')
11     format(3I12)
12     format(3F12.6)
13     format(E16.6)
14     format('END_DATAGRID_3D')
15     format('END_BLOCK_DATAGRID_3D')
      end subroutine write_xsf




 end module Module_XSF



