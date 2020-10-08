      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
     4     PERIOD)
C     
      INCLUDE 'ABA_PARAM.INC'

       parameter(zero=0.d0, half=0.5, one=1.d0, two=2.d0, three=3.d0, 
     1 four=4.d0, six=6.d0, eight=8.d0, twelve=12.d0,NGPT=4)
     
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),
     1     SVARS(NSVARS),ENERGY(8),PROPS(*),COORDS(MCRD,NNODE),
     2     U(NDOFEL),DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),
     3     PARAMS(3),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
     4     DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),
     5     JPROPS(*)
        DIMENSION SIG(3),EPS(3)
        DIMENSION B(3,NDOFEL),DB(3,NDOFEL),TB(NDOFEL,3)
        DIMENSION PHI(NNODE),PHIX(NNODE),PHIY(NNODE),PHIC(NNODE)
        DIMENSION PHIE(NNODE)
        DIMENSION GPX(NGPT),GPY(NGPT)
        DIMENSION GWEI(NGPT),GWE(2)
        integer :: numnode, numdof, ndim, nintpt
        parameter (nintpt = 2)
        real :: You, nu
        real :: J_DET
        real*8, dimension(:,:)::W(nintpt,1),Q(nintpt,1)
        real*8,dimension(:,:),allocatable:: E0(:,:), E1(:,:), E2(:,:)
        real*8 :: D(3,3), DINV(3,3)
        character*8 :: type_shape
        parameter (type_shape = 'L2')
        real*8, dimension(:)::shape_fun(2)
        real*8, dimension(:,:)::dNdxi(2,1) ! for 2D as derivative w.r.t xy
        real*8, dimension(:,:)::L1(3,2), L2(3,2)
        integer*8,dimension(:,:),allocatable:: ELEMENT(:,:)
        real*8,dimension(:,:),allocatable:: XX(:,:), YY(:,:), z(:,:)
        INTEGER :: ECON(2)
        integer*8, dimension(:):: gindex(4)
        real*8, dimension(:)::Pp(1,2), qq(1,2), JAC_M(2,2)
        
        real*8, dimension(:,:):: xetaN(2,1), yetaN(2,1)
        real*8, dimension(:,:):: xetadN(2,1), yetadN(2,1)
        real*8, dimension(:,:):: a1N(2,2), a2N(2,2)
        real*8, dimension(:,:):: b1eta(2,2), b2eta(2,2)
        real*8, dimension(:,:):: b1etaT(2,2), b2etaT(2,2)
        real*8, dimension(:,:):: eleezero(4,4), eleeone(4,4)
        real*8, dimension(:,:):: eleetwo(4,4)
        real*8, dimension(:,:):: term1(4,4), term2(4,4)
        real*8, dimension(:,:):: term3(4,4), tmp(16,1)
        real*8,dimension(:,:),allocatable:: E1t(:,:), invE0(:,:)
        real*8,dimension(:,:),allocatable:: term1z(:,:), term2z(:,:)
        real*8,dimension(:,:),allocatable:: term3z(:,:), term4z(:,:)
        real*8,dimension(:,:),allocatable:: Z1t(:,:), kmatt(:,:)
        !real*8, dimension(:,:):: vec(16,8)
        integer*8,dimension(:),allocatable:: marker_v(:)
        real*8,dimension(:),allocatable:: lmbda(:), lmbdaN(:)
        real*8,dimension(:,:),allocatable:: term1k(:,:), term2k(:,:)
        real*8,dimension(:,:),allocatable:: invterm1k(:,:), vec(:,:)
        !real*8, dimension(:,:)::term1K(8,8),term2k(8,8),invterm1k(8,8)
        ! new variables need to ---
        real*8, dimension(:,:):: b1(3,2), b2(3,2)
        real*8, dimension(:,:):: nl(2,4), dnl(2,4)
        real*8, dimension(:,:)::bb1(3,4),bb2(3,4),bb1t(4,3),bb2t(4,3)
        real*8, dimension(:,:):: DB1(3,4), DB2(3,4)
        ! for eigen value calculation
        real*8, dimension(:,:), allocatable:: VL(:,:), VR(:,:)
        real*8, dimension(:), allocatable:: WR(:), WI(:)
        
        INTEGER          NT, LDA, LDVL, LDVR
!        PARAMETER        ( NT = 16 )
        !INTEGER          LDA, LDVL, LDVR
        INTEGER          LWMAX
        PARAMETER        ( LWMAX = 10000 )
!        DOUBLE PRECISION  VL( LDVL, NT ), VR( LDVR, NT )
        DOUBLE PRECISION  WORK( LWMAX )
        INTEGER          INFO, LWORK
        !REAL             ZERO
        !parameter        (zero = 0.0)
        integer          IT,JT
        INTRINSIC        INT, MIN
        EXTERNAL         DGEEV
        !EXTERNAL         PRINT_EIGENVALUES, PRINT_EIGENVECTORS
        !External         PRINT_MATRIX


c
        
        nn=size(coords,2)
        NT=2*ndofel
        LDA = NT
        LDVL = NT 
        LDVR = NT 
        if(.not.allocated(VL))allocate(VL(NT,NT))
        if(.not.allocated(VR))allocate(VR(NT,NT))
        if(.not.allocated(WR))allocate(WR(NT))
        if(.not.allocated(WI))allocate(WI(NT))

        if(.not.allocated(element))allocate(element(nn,2))
        if(.not.allocated(XX))allocate(XX(nn,1))
        if(.not.allocated(YY))allocate(YY(nn,1))
        if(.not.allocated(marker_v))allocate(marker_v(2*ndofel))
        if(.not.allocated(Z))allocate(Z(2*ndofel,2*ndofel))
        if(.not.allocated(E1t))allocate(E1t(ndofel,ndofel))
        if(.not.allocated(invE0))allocate(invE0(ndofel,ndofel))
        if(.not.allocated(kmatt))allocate(kmatt(ndofel,ndofel))
        if(.not.allocated(term1z))allocate(term1z(ndofel,ndofel))
        if(.not.allocated(term2z))allocate(term2z(ndofel,ndofel))
        if(.not.allocated(term3z))allocate(term3z(ndofel,ndofel))
        if(.not.allocated(term4z))allocate(term4z(ndofel,ndofel))
        if(.not.allocated(lmbda))allocate(lmbda(2*ndofel))
        if(.not.allocated(lmbdaN))allocate(lmbdaN(ndofel))
        if(.not.allocated(term1k))allocate(term1k(ndofel,ndofel))
        if(.not.allocated(term2k))allocate(term2k(ndofel,ndofel))
        if(.not.allocated(invterm1k))allocate(invterm1k(ndofel,ndofel))
        if(.not.allocated(vec))allocate(vec(2*ndofel,ndofel))
        if(.not.allocated(z1t))allocate(z1t(2*ndofel,2*ndofel))


        N1=3 
        N2=NDOFEL
        if (size(coords,1).eq. 2) then
            ndim = 2
        else if (size(coords,1).eq.3) then
            ndim = 3
        end if
        numdof = ndim ! Elasticity problem
        numnode = size(coords,2)
        ! define local element connectivity---
	do i=1,nn-1,1
            ELEMENT(i,1)=i
            ELEMENT(i,2)=i+1
             enddo
            element(nn,1) = nn
            element(nn,2) = 1
        !=-------

         if (.not.allocated(E0))allocate(E0(numnode*numdof,
     1    numnode*numdof))  ! De-allocate later on
         if (.not.allocated(E1))allocate(E1(numnode*numdof,
     1      numnode*numdof))   ! De-allocate later on
         if (.not.allocated(E2))allocate(E2(numnode*numdof,
     1        numnode*numdof))   ! De-allocate later on
C==========step(1)=========C
c initialize parameters    C
C==========================C
        You  = PROPS(2)
        nu   = PROPS(3)     
        
  
        RHS(:,NRHS) = 0.0
        AMATRX(:,:) = 0.0
        call quadrature(W,Q,nintpt)
	   
        D=0.0	 
        D1=You/(1-nu*nu)  !***! DEFINE D MATRIX plane stress
        D(1,1)=D1
        D(1,2)=D1*nu
        D(2,1)=D1*nu
        D(2,2)=D1
        D(3,3)=D1*(1.0-nu)*0.5 
        !call inverse(D,Dinv,3)
       ! for planestrain
        !D1 = you/(1+nu)/(1-2*nu)
        !D(1,1) = D1*(1-nu)
        !D(1,2) = D1*nu
        !D(2,1) = D1*nu
        !D(2,2) = D1*(1-nu)
        !D(3,3) = D1*(0.5-nu)


! Define operators
       L1(1,1) = 1; L1(1,2) = 0
       L1(2,1) = 0; L1(2,2) = 0
       L1(3,1) = 0; L1(3,2) = 1

       L2(1,1) = 0; L2(1,2) = 0
       L2(2,1) = 0; L2(2,2) = 1
       L2(3,1) = 1; L2(3,2) = 0
! loop over edges
       
       XX(:,1) = COORDS(1,:)
       YY(:,1) = COORDS(2,:)
       !print *, 'start'
       !print *, coords 
       !print *, XX
       !print *, YY
       !pause

       E0 = 0.0
       E1 = 0.0
       E2 = 0.0
       !print *, numnode
       do iel = 1,numnode
       !print *, iel   
       ECON(1) = ELEMENT(IEL,1)
       ECON(2) = ELEMENT(IEL,2)
       !print *, econ
       PP(1,1) = XX(ECON(1),1)
       pp(1,2) = xx(econ(2),1)
       qq(1,1) = yy(econ(1),1)
       qq(1,2) = yy(econ(2),1)

       do ind = 1,2
       gindex(2*ind-1) = 2*econ(ind)-1
       gindex(2*ind) = 2*econ(ind)
       end do
       !print *, gindex
       !pause

! loop over GAUSS point
       do igp = 1, nintpt
          pt = Q(igp,1)
          call lagrange_basis(type_shape,pt,shape_fun,dNdxi)
!     1    numnode)
        
      dxdc=0.0d0    !**!  DXDC---J11

      x0 = sum(coords(1,:))/size(coords,2)
      y0 = sum(coords(2,:))/size(coords,2)
      !print *, x0, y0
      XS = shape_fun(1)*(pp(1,1)-x0)+shape_fun(2)*(pp(1,2)-x0)
      YS = shape_fun(1)*(qq(1,1)-Y0)+shape_fun(2)*(qq(1,2)-Y0)

    
      XS_DIFF = -HALF*(pp(1,1)-X0)+HALF*(pp(1,2)-X0)
      YS_DIFF = -HALF*(qq(1,1)-Y0)+HALF*(qq(1,2)-Y0)

      J_DET = XS*YS_DIFF - YS*XS_DIFF
      if (J_DET <= 0.0) then 
      print *, J_DET
      end if
      !pause
      do ii = 1,3
      do jj=1,2
       b1(ii,jj) = (1./J_DET)*(YS_DIFF*l1(ii,jj) 
     1                  - l2(ii,jj)*XS_DIFF)
      end do
      end do
      !print *, b1
      do ii = 1,3
      do jj=1,2
       b2(ii,jj) = (1./J_DET)*(-YS*l1(ii,jj) 
     1                  + l2(ii,jj)*XS)
      end do
      end do
      !print *, b2
      !pause 
      nl(1,1) = shape_fun(1)
      nl(1,2) = 0.0
      nl(1,3) = shape_fun(2)
      nl(1,4) = 0.0
      nl(2,1) = 0.0;
      nl(2,2) = shape_fun(1)
      nl(2,3) = 0.0
      nl(2,4) = shape_fun(2)
      dnl(1,1) = dNdxi(1,1)
      dnl(1,2) = 0.0
      dnl(1,3) = dNdxi(2,1)
      dnl(1,4) = 0.0
      dnl(2,1) = 0.0;
      dnl(2,2) = dNdxi(1,1)
      dnl(2,3) = 0.0
      dnl(2,4) = dNdxi(2,1)
      
      bb1 = matmul(b1,nl)
      bb2 = matmul(b2,dnl)
      !print *, bb1
      !print *, bb2
      !pause 
      do jk=1,3
      do Lk=1,4
      bb1t(Lk,jk)=bb1(jk,Lk)
      enddo
      enddo

      do jk=1,3
      do Lk=1,4
      bb2t(Lk,jk)=bb2(jk,Lk)
      enddo
      enddo
      
      DB1  = matmul(D,bb1)
      DB2 = matmul(D,bb2)
      !print *, DB1
      !print *, Db2
      !pause
      eleezero = 0.0
      eleeone = 0.0
      eleetwo = 0.0
      do k1=1,4
      do k2=1,4
      do k3=1,3
      term1 = W(igp,1)*bb1T(k1,k3)*db1(k3,k2)*j_DET
      term2 = W(igp,1)*bb2T(k1,k3)*db1(k3,k2)*J_DET
      term3 = W(igp,1)*bb2T(k1,k3)*db2(k3,k2)*J_DET
      
      eleezero(k1,k2)=eleezero(k1,k2) + term1(k1,k2) 
      eleeone(k1,k2)=eleeone(k1,k2) + term2(k1,k2)
      eleetwo(k1,k2)=eleetwo(k1,k2) + term3(k1,k2)
      end do
      end do
      end do


      
      E0(gindex,gindex) = E0(gindex,gindex) + eleezero
      E1(gindex,gindex) = E1(gindex,gindex) + eleeone
      E2(gindex,gindex) = E2(gindex,gindex) + eleetwo

! correct---  
      !print *,eleezero
      !print *, D
      !CALL PRINT_MATRIX('eleezero',4,4,eleezero,4)
      !pause 
      end do ! of gauss points
       
      end do ! end of edges

      !print *, ndofel
      !CALL PRINT_MATRIX('E0',ndofel,ndofel,E0,ndofel)
      !pause
 
      call inverse(E0, invE0,ndofel)

      do jk=1,ndofel
      do Lk=1,ndofel
      E1t(Lk,jk)=E1(jk,Lk)
      end do
      end do

      term1z = matmul(invE0,E1t)
      term2z = -invE0
      term3z = -E2 + matmul(E1,term1z)
      term4z = -matmul(E1,invE0)
      !CALL PRINT_MATRIX('E0',8,8,term1z,8)
      !pause
      z = 0.0
      z(1:ndofel,1:ndofel) = term1z
      z(1:ndofel,ndofel+1:2*ndofel) = term2z
      z(ndofel+1:2*ndofel,1:ndofel) = term3z
      z(ndofel+1:2*ndofel,ndofel+1:2*ndofel) = term4z
      !CALL PRINT_MATRIX('Z',2*ndofel,2*ndofel,Z,2*ndofel)
      !pause
      VR = 0.0
      LWORK =-1
      CALL DGEEV( 'Vectors', 'Vectors', NT, Z, LDA, WR, WI, VL,
     1 LDVL, VR, LDVR, WORK, LWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      CALL DGEEV( 'Vectors', 'Vectors', NT, Z, LDA, WR, WI, VL, 
     1 LDVL, VR, LDVR, WORK, LWORK, INFO )
      lmbda=0.0
      !call PRINT_EIGENVALUES( 'eigen_val', NT, WR, WI )
      !print *, WR
      !CALL PRINT_MATRIX('VR',16,16,VR,16)
      !pause 
      do Jk = 1,2*ndofel
      lmbda(Jk) = WR(Jk)
      marker_v(Jk) = Jk
      end do

      do ki = 1,2*ndofel
      do kj = 1,2*ndofel
      if ( lmbda(ki) < lmbda(kj) ) then
          t = lmbda(kj)
          lmbda(kj) = lmbda(ki)
          lmbda(ki) = t
          tt = marker_v(kj)
          marker_v(kj) = marker_v(ki)
          marker_v(ki) = tt
    
      endif

      end do
      end do

      !do ki=1,8
      
      !tmp=0.0
      !tmp(:,1) = VR(:,ki)
      !VR(:,ki) = VR(:,marker_v(Ki))
      !VR(:,marker_v(ki))=tmp(:,1)
      
      !end do
      !print *, marker_v
      VR = VR(:,marker_v)
      lmbdaN = lmbda(1:ndofel)
      !print *, lmbdaN
      vec = 0.0
      vec = VR(1:2*ndofel,1:ndofel)
      !CALL PRINT_MATRIX('before: vec',2*ndofel,ndofel,vec,ndofel)
      !print *, 'ndofel', ndofel
      vec(:,ndofel-1)=0.0
      vec(:,ndofel)=0.0

      ! need to check befor enforcing
      ! 0 and 1
      do ii=1,ndofel,2
      vec(ii,ndofel-1)=1.0
      end do
      do ii=1,ndofel,2
      vec(ii+1,ndofel)=1.0
      end do 
      !print *, vec
      !pause 
     
      !print*,'last',vec(:,1)
      !vec(1:2:8,7)=1.0
      !vec((1:2:8)+1,8)=1.0
      !CALL PRINT_MATRIX('vec',2*ndofel,ndofel,vec,ndofel)
      term1K = vec(1:ndofel,1:ndofel)
      term2k = vec(ndofel+1:2*ndofel,1:ndofel)
      call inverse(term1k,invterm1k,ndofel)
      kmatt = 0.0
      kmatt = matmul(term2k,invterm1k)
      AMATRX = kmatt
      rhs(:,1) = -matmul(AMATRX,u)
      !print *, ndofel
      !print *, kmatt
      !CALL PRINT_MATRIX('kmatt',ndofel,ndofel,kmatt,ndofel)
      !CALL PRINT_MATRIX('kmatt',8,8,term2k,8)
      !CALL PRINT_MATRIX('kmatt',8,8,invterm1k,8)
      !pause
      if(allocated(element)) deallocate(element)
      if(allocated(E0)) deallocate(E0)
      if(allocated(E1)) deallocate(E1)
      if(allocated(E2)) deallocate(E2)
      if(allocated(XX)) deallocate(XX)
      if(allocated(YY)) deallocate(YY)
      if(allocated(Z)) deallocate(Z)
      if(allocated(kmatt)) deallocate(kmatt)
      if(allocated(E1t)) deallocate(E1t)
      if(allocated(invE0)) deallocate(invE0)
      if(allocated(term1z)) deallocate(term1z)
      if(allocated(term2z)) deallocate(term2z)
      if(allocated(term3z)) deallocate(term3z)
      if(allocated(term4z)) deallocate(term4z)
      if(allocated(lmbda)) deallocate(lmbda)
      if(allocated(lmbdaN)) deallocate(lmbdaN)
      if(allocated(term1k)) deallocate(term1k)
      if(allocated(term2k)) deallocate(term2k)
      if(allocated(vec)) deallocate(vec)
      if(allocated(invterm1k)) deallocate(invterm1k)
      if(allocated(z1t)) deallocate(z1t)


      if(allocated(marker_v)) deallocate(marker_v) 


        if(allocated(VL))deallocate(VL)
        if(allocated(VR))deallocate(VR)
        if(allocated(WR))deallocate(WR)
        if(allocated(WI))deallocate(WI)
	RETURN
      END

!---  end code -----
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


!  =============================================================================
!
!     Auxiliary routine: printing a matrix.
!
      SUBROUTINE PRINT_MATRIX( DESC, M, N, A, LDA )
      CHARACTER*(*)    DESC
      INTEGER          M, N, LDA
      DOUBLE PRECISION             A( LDA, * )
!
      INTEGER          I, J
!
      WRITE(*,*)
      WRITE(*,*) DESC
      DO I = 1, M
         WRITE(*,9998) ( A( I, J ), J = 1, N )
      END DO
!
 9998 FORMAT( 11(:,1X,F6.4) )
      RETURN
      END

      SUBROUTINE PRINT_EIGENVALUES( DESC, N, WR, WI )
      CHARACTER*(*)    DESC
      INTEGER          N
      DOUBLE PRECISION             WR( * ), WI( * )
!
      DOUBLE PRECISION             ZERO
      PARAMETER        ( ZERO = 0.0 )
      INTEGER          J
!
      WRITE(*,*)
      WRITE(*,*) DESC
      DO J = 1, N
         IF( WI( J ).EQ.ZERO ) THEN
            WRITE(*,9998,ADVANCE='NO') WR( J )
         ELSE
            WRITE(*,9999,ADVANCE='NO') WR( J ), WI( J )
         END IF
      END DO
      WRITE(*,*)
!
 9998 FORMAT( 11(:,1X,F6.2) )
 9999 FORMAT( 11(:,1X,'(',F6.2,',',F6.2,')') )
      RETURN
      END



      SUBROUTINE DER(C,E,PHI,PHIX,PHIY,PHIC,PHIE,
     1 DXDC,DXDE,DYDC,DYDE,AJACOB,COORDS,MCRD,NNODE,NGPT)
      INCLUDE 'aba_param.inc'
      DIMENSION PHI(NNODE),PHIX(NNODE),PHIY(NNODE),PHIC(NNODE),
     1 PHIE(NNODE),COORDS(MCRD,NNODE)
      PARAMETER(ZERO=0.D0,FOURTH=0.25D0,HALF=0.5D0,
	1  ONE=1.D0,TWO=2.D0)
C     INTERPOLATION FUNCTIONS(1) C--zeta, E--eta
C
      PHI(1) = FOURTH*(ONE-C)*(ONE-E)
      PHI(2) = FOURTH*(ONE+C)*(ONE-E)
      PHI(3) = FOURTH*(ONE+C)*(ONE+E)
      PHI(4) = FOURTH*(ONE-C)*(ONE+E)
C
C     DERIVATIVES WRT TO C
C
      PHIC(1) = -FOURTH*(ONE-E)
      PHIC(2) = FOURTH*(ONE-E)
      PHIC(3) = FOURTH*(ONE+E)
      PHIC(4) = -FOURTH*(ONE+E)
C
C     DERIVATIVES WRT TO E
C
      PHIE(1) = -FOURTH*(ONE-C)
      PHIE(2) = -FOURTH*(ONE+C)
      PHIE(3) = FOURTH*(ONE+C)
      PHIE(4) = FOURTH*(ONE-C)
      DXDC=ZERO    !**!  DXDC---J11
      DXDE=ZERO    !**!  DXDE---J21
      DYDC=ZERO    !**!  DYDC---J12
      DYDE=ZERO    !**!  DYDE---J22  


	  DXDC=FOURTH*(-(1-E)*COORDS(1,1)+(1-E)*COORDS(1,2)
	1              +(1+E)*COORDS(1,3)-(1+E)*COORDS(1,4))
	  DYDC=FOURTH*(-(1-E)*COORDS(2,1)+(1-E)*COORDS(2,2)
	1              +(1+E)*COORDS(2,3)-(1+E)*COORDS(2,4))
	  DXDE=FOURTH*(-(1-C)*COORDS(1,1)-(1+C)*COORDS(1,2)
	1              +(1+C)*COORDS(1,3)+(1-C)*COORDS(1,4))
	  DYDE=FOURTH*(-(1-C)*COORDS(2,1)-(1+C)*COORDS(2,2)
	1              +(1+C)*COORDS(2,3)+(1-C)*COORDS(2,4))
C
C     CALCULATION OF JACOBIAN
C
      AJACOB=(DXDC*DYDE-DXDE*DYDC)  !**!AJACOB=|J|

C
C     DERIVATIVES WRT TO X AND Y
C
      DO 5 I=1,NNODE
        PHIX(I)=(DYDE*PHIC(I)-DYDC*PHIE(I))/AJACOB
        PHIY(I)=(-DXDE*PHIC(I)+DXDC*PHIE(I))/AJACOB
    5     CONTINUE
      RETURN
      END







      SUBROUTINE DER1D(C,E,PHI,PHIX,PHIY,PHIC,PHIE,
     1 DXDC,DXDE,DYDC,DYDE,AJACOB,COORDS,MCRD,NNODE,NGPT)
      INCLUDE 'aba_param.inc'
      DIMENSION PHI(NNODE),PHIX(NNODE),PHIY(NNODE),PHIC(NNODE),
     1 PHIE(NNODE),COORDS(MCRD,NNODE)
      PARAMETER(ZERO=0.D0,FOURTH=0.25D0,HALF=0.5D0,
	1  ONE=1.D0,TWO=2.D0)
C     INTERPOLATION FUNCTIONS(1) C--zeta, E--eta
C
   
      REAL :: XS,YS,XS_DIFF,YS_DIFF,J_DET
      REAL, DIMENSION(3,2) :: L1,L2
      REAL, DIMENSION(3,2) :: B_SMALL_1,B_SMALL_2
      REAL, DIMENSION(2,4) :: Nn,Nn_diff
            

      PHI(1) = HALF*(ONE-C)
      PHI(2) = HALF*(ONE+C)
C
C     DERIVATIVES WRT TO C
C
      PHIC(1) = -HALF
      PHIC(2) = HALF
C
C     DERIVATIVES WRT TO E
C
      DXDC=ZERO    !**!  DXDC---J11

      XS = PHI(1)*(COORDS(1,1)-X0)+PHI(2)*(COORDS(1,2)-X0)
      YS = PHI(1)*(COORDS(2,1)-Y0)+PHI(2)*(COORDS(2,2)-Y0)
      
      XS_DIFF = -HALF*(COORDS(1,1)-X0)+HALF*(COORDS(1,2)-X0)
      YS_DIFF = -HALF*(COORDS(2,1)-Y0)+HALF*(COORDS(2,2)-Y0)

      J_DET = XS*YS_DIFF - YS*XS_DIFF

      L1(1,1) = 1; L1(1,2) = 0
      L1(2,1) = 0; L1(2,2) = 0
      L1(3,1) = 0; L1(3,2) = 1

      L2(1,1) = 0; L2(1,2) = 0
      L2(2,1) = 0; L2(2,2) = 1
      L2(3,1) = 1; L2(3,2) = 0

      B_SMALL_1 = (L1*YS_DIFF-L2*XS_DIFF)/J_DET
      B_SMALL_2 = (L1*YS-L2*XS)/J_DET

      Nn = 0.0
      Nn(1,1) = PHI(1); Nn(1,2) = 0.0; Nn(1,3) = PHI(2); Nn(1,4) = 0.0
      Nn(2,1) = 0.0; Nn(2,2) = PHI(1); Nn(2,3) = 0.0; Nn(2,4) = PHI(2)

      Nn_diff(1,1) = -HALF; Nn_diff(1,3) = HALF
      Nn_diff(2,2) = -HALF; Nn_diff(2,4) = HALF
      
!      B1


	  DXDC=FOURTH*(-(1-E)*COORDS(1,1)+(1-E)*COORDS(1,2)
	1              +(1+E)*COORDS(1,3)-(1+E)*COORDS(1,4))
	  DYDC=FOURTH*(-(1-E)*COORDS(2,1)+(1-E)*COORDS(2,2)
	1              +(1+E)*COORDS(2,3)-(1+E)*COORDS(2,4))
	  DXDE=FOURTH*(-(1-C)*COORDS(1,1)-(1+C)*COORDS(1,2)
	1              +(1+C)*COORDS(1,3)+(1-C)*COORDS(1,4))
	  DYDE=FOURTH*(-(1-C)*COORDS(2,1)-(1+C)*COORDS(2,2)
	1              +(1+C)*COORDS(2,3)+(1-C)*COORDS(2,4))
C
C     CALCULATION OF JACOBIAN
C
      AJACOB=(DXDC*DYDE-DXDE*DYDC)  !**!AJACOB=|J|

C
C     DERIVATIVES WRT TO X AND Y
C
      DO 5 I=1,NNODE
        PHIX(I)=(DYDE*PHIC(I)-DYDC*PHIE(I))/AJACOB
        PHIY(I)=(-DXDE*PHIC(I)+DXDC*PHIE(I))/AJACOB
    5     CONTINUE
      RETURN
      END


	subroutine randomsys3(N1)
   	 implicit none
	    real(kind=8), dimension(:),allocatable :: x,b,work
	    real(kind=8), dimension(:,:),allocatable :: a
	    real(kind=8) :: errnorm, xnorm, rcond, anorm, colsum
	    integer :: i, info, lda, ldb, nrhs, N1, j,n
	    integer, dimension(:), allocatable :: ipiv
	    integer, allocatable, dimension(:) :: iwork
	    character, dimension(1) :: norm

    ! initialize random number generator seed
    ! if you remove this, the same numbers will be generated each
    ! time you run this code.
	    call init_random_seed()  

	    print *, "Input n ... "
!  	  read *, N1
           n = N1

	    allocate(a(n,n))
	    allocate(b(n))
	    allocate(x(n))
	    allocate(ipiv(n))
           print*,'allocated'
	    call random_number(a)
	    call random_number(x)

	    b = matmul(a,x)    ! compute RHS
           print*,'b',b

	    anorm = 0.d0
	    do j=1,n
        colsum = 0.d0
        do i=1,n
            colsum = colsum + abs(a(i,j))
            enddo
        anorm = max(anorm, colsum)
        enddo
         print*,anorm
    
	    nrhs = 1 ! number of right hand sides in b
	    lda = n  ! leading dimension of a
	    ldb = n  ! leading dimension of b

	    call dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)

    ! compute 1-norm of error
	    errnorm = 0.d0
	    xnorm = 0.d0
	    do i=1,n
        errnorm = errnorm + abs(x(i)-b(i))
        xnorm = xnorm + abs(x(i))
        enddo

    ! relative error in 1-norm:
	    errnorm = errnorm / xnorm


    ! compute condition number of matrix:
    ! note: uses A returned from dgesv with L,U factors:

	    allocate(work(4*n))
	    allocate(iwork(n))
	    norm = '1'  ! use 1-norm
	    call dgecon(norm,n,a,lda,anorm,rcond,work,iwork,info)

	    if (info /= 0) then
        print *, "*** Error in dgecon: info = ",info
        endif

	    print*, n, 1.d0/rcond, errnorm
!	201 format("For n = ",i4," the approx. condition number is ",&
!            e10.3,/," and the relative error in 1-norm is ",e10.3)    

	    deallocate(a,b,ipiv)
	    deallocate(work,iwork)

	end subroutine
	
       SUBROUTINE init_random_seed()
        INTEGER :: i, n, clock
        INTEGER, DIMENSION(:), ALLOCATABLE :: seed

        CALL RANDOM_SEED(size = n)
        ALLOCATE(seed(n))

        CALL SYSTEM_CLOCK(COUNT=clock)

        seed = clock + 37 * (/ (i - 1, i = 1, n) /)
        CALL RANDOM_SEED(PUT = seed)

        print *, "Using random seed = ", seed
        print *, " "

        DEALLOCATE(seed)
      END SUBROUTINE


        subroutine lagrange_basis(type_shape,pt,shape_fun,dNdxi)
        
        real*8 :: pt
        character*8 :: type_shape
        real*8, dimension(:)::shape_fun(2)
        real*8, dimension(:,:)::dNdxi(2,1)
        real*8 :: xi
        
        xi = pt
     
        if (type_shape .eq. 'L2') then
!	if (.not. allocated(shape_fun)) allocate(shape_fun(2))
!	if (.not. allocated(dNdxi)) allocate(dNdxi(2,1))
            shape_fun(1) = (1.0d0-xi)/2.0d0
            shape_fun(2) = (1.0d0+xi)/2.0d0
            dNdxi(1,1) = -1.0d0/2.0d0
            dNdxi(2,1) = 1.0d0/2.0d0
        endif

        end subroutine 
        
        
        subroutine quadrature(W,Q,nintpt)
        real*8, dimension(:,:)::W(nintpt,1),Q(nintpt,1)

        if (nintpt.eq. 2) then
            W(1,1) = 1.0
            W(2,1) = 1.0
            Q(1,1) = 0.577350269189626
            Q(2,1) =-0.577350269189626
        endif
        end subroutine 


        SUBROUTINE INVERSE(A,a1,n)
        implicit none
        integer n
        real*8 :: A(n,n)
        real*8 :: a1(size(A,1),size(A,2))
        real*8 :: work(size(A,1))            ! work array for LAPACK
        integer*8 :: info,ipiv(size(A,1))     ! pivot indices

      ! Store A in Ainv to prevent it from being overwritten by LAPACK
        a1 = A
      !n = size(A,1)
      ! SGETRF computes an LU factorization of a general M-by-N matrix A
      ! using partial pivoting with row interchanges.
        call DGETRF(n,n,a1,n,ipiv,info)
        !if (info.ne.0) stop 'Matrix is numerically singular!'  ! for real*8 datatype this is commented
      ! SGETRI computes the inverse of a matrix using the LU factorization
      ! computed by SGETRF.
        call DGETRI(n,a1,n,ipiv,work,size(A,1),info)
        !if (info.ne.0) stop 'Matrix inversion failed!'  !  for real*8 datatype this is commented
        end SUBROUTINE INVERSE

! Websites : https://github.com/b-fg/LAPACK_helper/blob/master/lapack_helper.f90
! https://stackoverflow.com/questions/46894223/wrong-inverse-matrix-using-zgetri-in-fortran
! https://stackoverflow.com/questions/20527031/fortran-and-matlab-return-different-eigenvalues-for-same-matrix?rq=1
! https://ww2.odu.edu/~agodunov/computing/programs/book2/Ch06/Inverse.f90 for inverse of the matrix w/o using lapack
! https://stackoverflow.com/questions/26475987/lapack-inversion-routine-strangely-mixes-up-all-variables % Inverse using lapack
! https://www.nag.com/numeric/fl/nagdoc_fl23/examples/source/f07ajfe.f90 Inverse of the matrix
