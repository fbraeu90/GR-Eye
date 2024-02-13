C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
C                             Authors
C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

C                 Fabian A. Braeu & Michael J.A. Girard

C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
C                         Global Variables
C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

      module global_variables
            real(8) :: buffer_r
            integer(8) :: buffer_i

C Note: In case you remesh your geometry, please update the following parameters
            integer, parameter :: nnodes    = 1668
            integer, parameter :: nhexelements = 900
            integer, parameter :: nwedgeelements = 30
            integer, parameter :: nelements = 930

            real(8),    dimension(nnodes,   3) :: n_table
            integer(8), dimension(nelements,8) :: e_table
            integer(8), dimension(8) :: tmp

            real(8),    dimension(nelements,3) :: local_x
            real(8),    dimension(nelements,3) :: local_y
            real(8),    dimension(nelements,3) :: local_z

      end module      

C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
C                         Global Variables
C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     * RPL,DDSDDT,DRPLDE,DRPLDT,
     * STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     * NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     * CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      use global_variables
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     * DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     * STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     * PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
C                         Declarations
C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

      Real(8), parameter :: zero=0.d0, half=0.5d0, one=1.d0, 
     *         two=2.d0, three=3.d0, four=4.d0, five=5.d0,
     *         six=6.d0, seven=7.d0, eight=8.d0, nine=9.d0

      Real(8), dimension(3,3) :: I
      Real(8), dimension(6)   :: Iv

      Real(8), dimension(2,3) :: a0, a, agr
      Real(8), dimension(2)   :: lambda, lambdae, I4e, W4e, W44e, Wle, Wlle
      Real(8)                 :: lambdae_h
      
      Real(8), dimension(3,3)   :: x2k, y2k, z2k, f2k
      Real(8), dimension(6)     :: f2k_v

      Real(8) :: Rho0_t0
      Real(8), dimension(2)   :: Rho0_t0_f, Rho0_f_dot, Rho0_f

      Real(8), dimension(3,3) :: F, Fbe, Fg, Fe, Fg_inverse, Fgr_inverse 
      Real(8), dimension(2,3,3) :: Fr, Fgr, Fef, F_prestretch
      Real(8), dimension(3,3) :: Ce, Cbe
      Real(8), dimension(3,3) :: Be, Bbe
      Real(8), dimension(6)   :: Bbve
      Real(8) :: J,Je

      Real(8) :: I1e, I1be

      Real(8), dimension(6)   :: stress_v, stress_m 
      Real(8), dimension(2,6) :: stress_f

      Real(8), dimension(6,6)   :: stiff_v, stiff_m 
      Real(8), dimension(2,6,6) :: stiff_f
      Real(8), dimension(6,6)   :: DijDkl, S, BbeXI, IXBbe, f4k, cp, Bbe_jau, stiff_f_jau

C     T is the half life time of collagen
      Real(8) :: k_sig, T

C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
C                   General Variables/Definitions
C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

C     Identity matrix (2nd order)
      I = reshape((/ one, zero, zero, 
     *              zero,  one, zero, 
     *              zero, zero, one /), shape(I))   

C     Voigt mapping for identity matrix
      call voigt(I,Iv)

C     Identity matrix (I⨂I - 4th order)
      DijDkl = reshape((/  one,  one,  one, zero, zero, zero, 
     *                     one,  one,  one, zero, zero, zero, 
     *                     one,  one,  one, zero, zero, zero, 
     *                    zero, zero, zero, zero, zero, zero, 
     *                    zero, zero, zero, zero, zero, zero, 
     *                    zero, zero, zero, zero, zero, zero /), shape(DijDkl))  

C     Identity matrix (S - 4th order)
      S = reshape((/  one, zero, zero, zero, zero, zero, 
     *               zero,  one, zero, zero, zero, zero, 
     *               zero, zero,  one, zero, zero, zero, 
     *               zero, zero, zero, half, zero, zero, 
     *               zero, zero, zero, zero, half, zero, 
     *               zero, zero, zero, zero, zero, half /), shape(S))  

C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
C                    Biomechanical Properties
C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

C     Bulk modulus (for volumetric strain energy)
      K = props(1)/props(5)

C     Neo-Hookean shear modulus for ground substance matrix
      c1 = props(2)/props(5)

C     Stress coefficient of collagen fibers
      c3 = props(3)/props(6)

C     Uncrimping rate of collagen fibers
      c4 = props(4)
      
C     Reference mass densities (time 0)
      Rho0_t0_m = props(5)

      if(CMNAME.EQ.'PERIPAP_SCLERA') then
            Rho0_t0_f(1) = 0.9*props(6)
            Rho0_t0_f(2) = 0.1*props(6)
      else if(CMNAME.EQ.'POST_SCLERA') then 
            Rho0_t0_f(1) = 0.5*props(6)
            Rho0_t0_f(2) = 0.5*props(6)
      else if(CMNAME.EQ.'LC') then 
            Rho0_t0_f(1) = 0.1*props(6)
            Rho0_t0_f(2) = 0.9*props(6)
      endif

C     Total tissue reference mass density
      Rho0_t0 = Rho0_t0_m + sum(Rho0_t0_f)

C     Growth parameter [days^-1]
      k_sig = props(7)

C     Growth type: 0 - mass density growth; 1 - transmural volumetric growth
      growth_type = props(8)

C     Half life time of collagen
      T = props(9)

C     Homeostatic stretch 
      lambdae_h = props(10)

C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
C                   Definition of State Variables
C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= 

C     11 -> Mass density (matrix) in the reference configuration
C     12, 22, ... -> Mass density (fiber) in the reference configuration (first digit is the fiber number)
C     13 -> Mass density (matrix) in the current configuration
C     14, 24, ... -> Mass density (fiber) in the current configuration (first digit is the fiber number)
C     15 -> J
C     16, 26, ... -> Remodeling stretch
C     17, 27, ... -> Single fiber stress
C     18, 28, ... -> Homeostatic fiber stress
C     19, 29, ... -> Fiber prestretch 
C     33 -> Shear stiffness (matrix)

C     This is used for display purposes in abaqus cae
      statev(13) = statev(11)/J
      do kf=1,2
            statev(kf*10+4) = statev(kf*10+2)/J
      enddo
      statev(15) = J

C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
C                        Initialization Step
C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= 
      
      if(KSTEP.EQ.1) then
            statev(11) = Rho0_t0_m
            do kf=1,2
                  statev(kf*10+2) = Rho0_t0_f(kf)
                  statev(kf*10+6) = one ! Remodeling stretch
                  statev(kf*10+9) = one ! Fiber prestretch
            enddo

            statev(33) = c1
      endif
      
      if(KSTEP.EQ.3) then
           if(CMNAME.EQ.'PERIPAP_SCLERA' .OR. CMNAME.EQ.'LC') then

C =-=-=-=-=-=-=-=-=-=-=-=-=-=- Stiffness Reduction -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
                  statev(33) = c1 - (TIME(1)+DTIME)/100.d0 * 0.85 * c1

            endif
      endif

C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
C                         Growth Tensor
C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= 
      
C     Kronecker for local fiber coordinate system 
      call kronecker_product(local_x(NOEL,:),local_x(NOEL,:),x2k)
      call kronecker_product(local_y(NOEL,:),local_y(NOEL,:),y2k)
      call kronecker_product(local_z(NOEL,:),local_z(NOEL,:),z2k)

C     Growth tensor isotropic
C      Fg = I
C      Fg = ((statev(12)+statev(22)+statev(11))/Rho0_t0)**(one/three) * x2k +
C     *     ((statev(12)+statev(22)+statev(11))/Rho0_t0)**(one/three) * y2k + 
C     *     ((statev(12)+statev(22)+statev(11))/Rho0_t0)**(one/three) * z2k

C     Growth perpendicular to preferred fiber direction
C     Fg = I
C      Fg =                                                x2k +
C     *     ((statev(12)+statev(22)+statev(11))/Rho0_t0)**(one/two) * y2k + 
C     *     ((statev(12)+statev(22)+statev(11))/Rho0_t0)**(one/two) * z2k


      if(growth_type.EQ.0) then

C           Mass density growth
            Fg = I
            Fg(1,1) = one
            Fg(2,2) = one
            Fg(3,3) = one

      else if(growth_type.EQ.1) then

C           Growth in transmural direction
            Fg = I
            Fg =                                                x2k +
     *                                                          y2k + 
     *           ((statev(12)+statev(22)+statev(11))/Rho0_t0) * z2k

      endif

C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
C                         Remodeling Tensor
C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= 

C     Remodeling tensor
      Fr(1,:,:) = I
      Fr(2,:,:) = I
      Fr(1,:,:) = statev(16)*x2k + one/statev(16)**(one/two)*y2k + one/statev(16)**(one/two)*z2k
      Fr(2,:,:) = statev(26)*y2k + one/statev(26)**(one/two)*x2k + one/statev(26)**(one/two)*z2k

C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
C                         Prestretch Tensor
C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= 

C     Prestretch tensor
      F_prestretch(1,:,:) = I
      F_prestretch(2,:,:) = I
      F_prestretch(1,:,:) = statev(19)*x2k + one/statev(19)**(one/two)*y2k + one/statev(19)**(one/two)*z2k
      F_prestretch(2,:,:) = statev(29)*y2k + one/statev(29)**(one/two)*x2k + one/statev(29)**(one/two)*z2k

C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
C                       Deformation Gradients
C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
C     Deformation gradient tensor
      F = DFGRD1
      call determinant(F,J)

C     Matrix elastic deformation gradient (note: Fr = 1 for the matrix) 
      call MatrixInverse(Fg,Fg_inverse)
      Fe = MATMUL(F,Fg_inverse)

C     Inelastic growth & remodeling deformation gradient 
      do kf=1,2
            Fgr(kf,:,:) = MATMUL(Fg,Fr(kf,:,:))
C     Include prestretch
            Fgr(kf,:,:) = MATMUL(Fgr(kf,:,:),F_prestretch(kf,:,:))
      enddo

C     Fiber elastic deformation gradient
      do kf=1,2
            call MatrixInverse(Fgr(kf,:,:),Fgr_inverse)
            Fef(kf,:,:) = MATMUL(F,Fgr_inverse)
      enddo

C     Determinant of Fe
      call determinant(Fe,Je)

C     Isochoric elastic deformation gradient tensor
      Fbe = Je**(-one/three)*Fe 

C     General and isochoric elastic Right-Cauchy-Green tensor                                    
      Ce  = MATMUL(transpose(Fe ), Fe)
      Cbe = MATMUL(transpose(Fbe),Fbe)

C     General and isochoric elastic Left-Cauchy-Green tensor                                    
      Be   = MATMUL(Fe, transpose(Fe ))
      Bbe  = MATMUL(Fbe,transpose(Fbe))
      call voigt(Bbe,Bbve)

C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
C                    Collagen Fiber Directions
C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

C     Fiber vector/direction a0 (reference configuration)
      a0(1,:) = local_x(NOEL,:)
      a0(2,:) = local_y(NOEL,:)

      do kf=1,2
C           Ensure it is normalized
            a0(kf,:) = a0(kf,:)/norm2(a0(kf,:))

C           Fiber stretch lambda
            lambda(kf) = norm2(MATMUL(F,a0(kf,:)))
C           Fiber vector/direction a (current configuration)
            a(kf,:) = MATMUL(F,a0(kf,:))/lambda(kf)

C           Fiber vector/direction (intermediate configuration)
            agr(kf,:) = MATMUL(Fgr(kf,:,:),a0(kf,:))/norm2(MATMUL(Fgr(kf,:,:),a0(kf,:)))
C           Elastic fiber stretch (intermediate to current configuration)
            lambdae(kf) = norm2(MATMUL(Fef(kf,:,:),agr(kf,:)))
C           
            if(KSTEP.EQ.2) then

C =-=-=-=-=-=-=-=-=-=-=-=-=-=- Prestretching -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
                  lambdae(kf) = one + (TIME(1)+DTIME)/100.d0*(lambdae_h-one)
                  statev(kf*10+9) = lambda(kf)/lambdae(kf)

            endif
      enddo

C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
C                          Invariants
C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

C- - - - - - - - - - - - - - - I1 & I1b - - - - - - - - - - - - - - - - 
      I1e  = Ce(1,1)  + Ce(2,2)  + Ce(3,3)
      I1be = Cbe(1,1) + Cbe(2,2) + Cbe(3,3)

C- - - - - - - - - - - - - - - I4 - - - - - - - - - - - - - - - - - - -  
      do kf=1,2
            I4e(kf) = lambdae(kf)**two      
      enddo

C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
C                         Cauchy Stress
C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=  

C- - - - - - - - - - - V o l u m e t r i c - - - - - - - - - - - - - - -
      stress_v = zero
      do k1 = 1,6
            stress_v(k1) = two/J*K*(Je-one)*Je*Iv(k1)
      enddo        

C- - - - - - - - - - - - - I s o c h o r i c - - - - - - - - - - - - - - - - -
      stress_m = zero
      do k1=1,6
            stress_m(k1) = two/J*statev(33)*( Bbve(k1) - I1be/three*Iv(k1) )
      enddo

C- - - - - - - - - - - - - F i b e r s - - - - - - - - - - - - - - - - -
      stress_f = zero
      do kf=1,2
C           Kronecker product between a and a, and its voigt mapping
            call kronecker_product(a(kf,:),a(kf,:),f2k)
            call voigt(f2k,f2k_v)
            
C     -=-=-=-=-=-= Fung strain energy -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            W4e(kf) = c3*(I4e(kf)-one)*exp(c4*(I4e(kf)-one)**two)
            Wle(kf) = two*lambdae(kf)*W4e(kf)

C     -=-=-=-= No resistance to compression -=-=-=-=-=-=-=-=-=-=-=-=-=-=
            if (lambdae(kf).LT.one) then
                  W4e(kf) = zero
            endif 

            do k1=1,6
                  stress_f(kf,k1) = two/J*I4e(kf)*W4e(kf)*f2k_v(k1)
            enddo

C------------------------------------------------------------------------
C     Local Cauchy stress along the fiber (fiber stress divided by the 
C.    current fiber mass density, hence no J!)
C------------------------------------------------------------------------
            statev(kf*10+7) = two*I4e(kf)*W4e(kf)

C-----------------------------------------------------------------------
C     The homeostatic fiber stress is constant and is saved after the
C     second step in which we apply pre-strech. The pre-stretch is
C     equal to the homeostatic stretch defined in the input file.
C-----------------------------------------------------------------------
            if(KSTEP.EQ.2) then
                  statev(kf*10+8) = statev(kf*10+7)
            endif       
      enddo
      
C- - - - - - - - - - - - - T o t a l - - - - - - - - - - - - - - - - - -
      stress = zero
      do k1=1,6
            stress(k1) = statev(11)*stress_v(k1) + statev(11)*stress_m(k1) 
     *                 + statev(12)*stress_f(1,k1) + statev(22)*stress_f(2,k1)
      enddo

C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
C                        Stiffness Tensor
C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=  
      
C- - - - - - - - - - - V o l u m e t r i c - - - - - - - - - - - - - - -

      stiff_v = zero
C     Jaumann rate of Cauchy stress --> effective elasticity tensor (as defined in Holzapfel's book)
      stiff_v = four*K*Je/J*( (Je-half)*DijDkl)

C- - - - - - - - - - - - - I s o c h o r i c - - - - - - - - - - - - - - - - -
C    Additional 4th order tensor due to Jaumann rate of Cauchy stress used in Abaqus
      BBe_jau = reshape((/two*Bbe(1,1), zero        , zero        , Bbe(1,2)                , Bbe(1,3)     , zero         , 
     *                    zero        , two*Bbe(2,2), zero        , Bbe(1,2)                , zero         , Bbe(2,3)     ,
     *                    zero        , zero        , two*Bbe(3,3), zero                    , Bbe(1,3)     , Bbe(2,3)     , 
     *                    Bbe(1,2)    , Bbe(1,2)    , zero        , half*(Bbe(1,1)+Bbe(2,2)), half*Bbe(2,3), half*Bbe(1,3),
     *                    Bbe(1,3)    , zero        , Bbe(1,3)    , half*Bbe(2,3), half*(Bbe(1,1)+Bbe(3,3)), half*Bbe(1,2),
     *                    zero        , Bbe(2,3)    , Bbe(2,3)    , half*Bbe(1,3), half*Bbe(1,2), half*(Bbe(2,2)+Bbe(3,3))     
     *               /), shape(BBe_jau))
C                       _
C     4th order tensor (Be⨂I)                           
      BbeXI = reshape((/ Bbe(1,1), Bbe(1,1), Bbe(1,1), zero, zero, zero, 
     *                   Bbe(2,2), Bbe(2,2), Bbe(2,2), zero, zero, zero, 
     *                   Bbe(3,3), Bbe(3,3), Bbe(3,3), zero, zero, zero, 
     *                   Bbe(1,2), Bbe(1,2), Bbe(1,2), zero, zero, zero, 
     *                   Bbe(1,3), Bbe(1,3), Bbe(1,3), zero, zero, zero, 
     *                   Bbe(2,3), Bbe(2,3), Bbe(2,3), zero, zero, zero /), shape(BbeXI))  
C                          _
C     4th order tensor (I⨂Be)                           
      IXBbe = reshape((/ Bbe(1,1), Bbe(2,2), Bbe(3,3), Bbe(1,2), Bbe(1,3), Bbe(2,3), 
     *                   Bbe(1,1), Bbe(2,2), Bbe(3,3), Bbe(1,2), Bbe(1,3), Bbe(2,3), 
     *                   Bbe(1,1), Bbe(2,2), Bbe(3,3), Bbe(1,2), Bbe(1,3), Bbe(2,3), 
     *                       zero,     zero,     zero,     zero,     zero,     zero, 
     *                       zero,     zero,     zero,     zero,     zero,     zero, 
     *                       zero,     zero,     zero,     zero,     zero,     zero /), shape(IXBbe))  

      stiff_m = zero
      stiff_m = four/(three*J)*statev(33)*( one/three*I1be*DijDkl - BbeXI - IXBbe + (three/two)*BBe_jau)

C- - - - - - - - - - - - - F i b e r s - - - - - - - - - - - - - - - - -
      stiff_f = zero

      do kf=1,2
C           Kronecker product between a and a
            call kronecker_product(a(kf,:),a(kf,:),f2k)

C     4th order tensor (a⨂a⨂a⨂a)                           
      f4k = reshape((/ f2k(1,1)**two    ,f2k(1,1)*f2k(2,2),f2k(1,1)*f2k(3,3),f2k(1,1)*f2k(1,2),f2k(1,1)*f2k(1,3),f2k(1,1)*f2k(2,3), 
     *                 f2k(2,2)*f2k(1,1),f2k(2,2)**two    ,f2k(2,2)*f2k(3,3),f2k(2,2)*f2k(1,2),f2k(2,2)*f2k(1,3),f2k(2,2)*f2k(2,3),
     *                 f2k(3,3)*f2k(1,1),f2k(3,3)*f2k(2,2),f2k(3,3)**two    ,f2k(3,3)*f2k(1,2),f2k(3,3)*f2k(1,3),f2k(3,3)*f2k(2,3),
     *                 f2k(1,2)*f2k(1,1),f2k(1,2)*f2k(2,2),f2k(1,2)*f2k(3,3),f2k(1,2)**two    ,f2k(1,2)*f2k(1,3),f2k(1,2)*f2k(2,3),
     *                 f2k(1,3)*f2k(1,1),f2k(1,3)*f2k(2,2),f2k(1,3)*f2k(3,3),f2k(1,3)*f2k(1,2),f2k(1,3)**two    ,f2k(1,3)*f2k(2,3),         
     *                 f2k(2,3)*f2k(1,1),f2k(2,3)*f2k(2,2),f2k(2,3)*f2k(3,3),f2k(2,3)*f2k(1,2),f2k(2,3)*f2k(1,3),f2k(2,3)**two     
     *               /), shape(f4k)) ! Need to be here as we can't exceed 132 columns
            
C     Additional 4th order tensor due to Jaumann rate of Cauchy stress used in Abaqus
      stiff_f_jau = reshape((/ two*stress_f(kf,1), zero       , zero, stress_f(kf,4), stress_f(kf,5), zero, 
     *                  zero         , two*stress_f(kf,2) , zero         , stress_f(kf,4), zero, stress_f(kf,6),
     *                  zero         , zero         , two*stress_f(kf,3), zero, stress_f(kf,4), stress_f(kf,6), 
     *                  stress_f(kf,4), stress_f(kf,4), zero, half*(stress_f(kf,1)+stress_f(kf,2)) , half*stress_f(kf,6), 
     *                  half*stress_f(kf,5),
     *                  stress_f(kf,5), zero, stress_f(kf,5), half*stress_f(kf,6), half*(stress_f(kf,1)+stress_f(kf,3)), 
     *                  half*stress_f(kf,4),
     *                  zero, stress_f(kf,6), stress_f(kf,6), half*stress_f(kf,5) , half*stress_f(kf,4), 
     *                  half*(stress_f(kf,2)+stress_f(kf,3))     
     *               /), shape(stiff_f_jau))

C     -=-=-=-=-=-= Fung strain energy derivatives -=-=-=-=-=-=-=-=-=-=-=
            W44e(kf) = c3*(one+two*c4*(I4e(kf)-one)**two)*exp(c4*(I4e(kf)-one)**two)
            Wlle(kf) = two*W4e(kf)+four*I4e(kf)*W44e(kf)

C     -=-=-=-=-=-=- No resistance to compression -=-=-=-=-=-=-=-=-=-=-=-
            if (lambdae(kf).LT.one) then
                  W44e(kf) = zero
            endif

C           Jaumann rate of Cauchy stress --> effective elasticity tensor (as defined in Holzapfel's book)
            stiff_f(kf,:,:) = four/J*W44e(kf)*I4e(kf)**two*f4k + stiff_f_jau
      enddo

C- - - - - - - - - - - - - T o t a l - - - - - - - - - - - - - - - - - -
      do k1=1,6
        do k2=1,6
            ddsdde(k1,k2) = statev(11)*stiff_v(k1,k2)   + statev(11)*stiff_m(k1,k2) 
     *                    + statev(12)*stiff_f(1,k1,k2) + statev(22)*stiff_f(2,k1,k2)
        enddo
      enddo

C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
C                         Growth Evolution
C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= 

C     Growth is activated at step 4 through this evolution equation
      if(KSTEP.EQ.4) then
            
            do kf=1,2

                  if (lambdae(kf).LE.one) then
                        statev(kf*10+2) = statev(kf*10+2)
            
                  else if (lambdae(kf).GT.one) then
                        Rho0_f(kf) = statev(kf*10+2)
                        Rho0_f_dot(kf) = Rho0_f(kf)*k_sig 
     *                                * ( statev(kf*10+7) - statev(kf*10+8) ) / statev(kf*10+8)
                        statev(kf*10+2) = Rho0_f(kf) + DTIME*Rho0_f_dot(kf)
                  endif
            enddo
      endif

C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
C                         Remodeling Evolution
C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= 

C     Remodeling is activated at step 4 through this evolution equation
      if(KSTEP.EQ.4) then
            
            do kf=1,2

                  if (lambdae(kf).LE.one) then
                        statev(kf*10+6) = statev(kf*10+6)
                  
                  else if (lambdae(kf).GT.one) then
                        statev(kf*10+6) = statev(kf*10+6) + DTIME*(
     *                  (Rho0_f_dot(kf)/Rho0_f(kf)+one/T)
     *                  * statev(kf*10+6)/lambdae(kf)
     *                  * (Wle(kf)+lambdae(kf)*Wlle(kf))**(-one)
     *                  * ( statev(kf*10+7) - statev(kf*10+8) )
     *                                          )
                  endif
            enddo
      endif

      RETURN
      END

C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
C                C U S T O M - F U N C T I O N S
C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
C                    To Read External Files
C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

      subroutine uexternaldb(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
C
      use global_variables
      include 'ABA_PARAM.INC'

      character(256) filename
      character(256) jobdir
      character(256) jobname
      character(5) line

      real(8) :: x, y, z, r
      real(8) :: theta, phi, sin_theta, cos_theta, sin_phi, cos_phi
      real(8) :: theta_basis(3), phi_basis(3)


      call getoutdir(jobdir,lenjobdir)
      call getjobname(jobname,lenjobname)
      filename = trim(jobdir) // '/' // trim(jobname)
      filename(lenjobdir+lenjobname+2:lenjobdir+lenjobname+40)='.inp'

      if (lop.EQ.0) then ! Called at the beginning
            open(UNIT=10,file=filename, status='unknown')

C     Read node coordinates table - - - - - - - - - - - - - - - - - - - 
            do while(line.NE.'*Node')
                  read(10,*) line
            enddo

            do i=1,nnodes
                  read(10,*) buffer_r, n_table(i,1), n_table(i,2), n_table(i,3)
            enddo

C     Read element connectivity table - - - - - - - - - - - - - - - - - 
            do while(line.NE.'*Elem')
                  read(10,*) line
            enddo

            do i=1,nhexelements                 
                  read(10,*) buffer_i, tmp(1), tmp(2), tmp(3), tmp(4),
     *                                 tmp(5), tmp(6), tmp(7), tmp(8)
                  e_table(buffer_i,:) = tmp
            enddo

            do while(line.NE.'*Elem')
                  read(10,*) line
            enddo

            read(10,*) line

            do i=1,nwedgeelements
                  read(10,*) buffer_i, tmp(1), tmp(2), tmp(3), tmp(4),
     *                                 tmp(5), tmp(6)
     
C           Dummy coordinates for an approximate center of element computation
                  tmp(7) = tmp(1)
                  tmp(8) = tmp(6)
                  e_table(buffer_i,:) = tmp
            enddo

            close(10)


C           Compute Fiber direction for each element - - - - - - - - - - - - - 
            do i=1,nelements

C                 Calculate element center
                  x = zero
                  y = zero
                  z = zero
                  do j=1,8
                        x = x + n_table(e_table(i,j),1)/8.d0
                        y = y + n_table(e_table(i,j),2)/8.d0
                        z = z + n_table(e_table(i,j),3)/8.d0
                  enddo

C                 Calculate radius of element center
                  r = sqrt(x**2 + y**2 + z**2)

C                 Calculate spherical coordinates
                  phi = acosd(z / r)
                  theta = atand(y / x) + 180.d0
            
C                 Calculate trigonometric functions
                  sin_theta = sind(theta)
                  cos_theta = cosd(theta)
                  sin_phi = sind(phi)
                  cos_phi = cosd(phi)

C                 Calculate basis vectors
                  local_x(i,:) = (/ cos_theta*cos_phi, sin_theta*cos_phi, -sin_phi /)
                  local_y(i,:) = (/ -sin_theta, cos_theta, 0.d0 /)
                  local_z(i,:) = (/ cos_theta*sin_phi, sin_theta*sin_phi, cos_phi /)

            enddo

      endif

      return
      end

C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
C                       Determinant
C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      subroutine determinant(A,det)

      Real(8), dimension(3,3) :: A
      Real(8) :: det
        
      det = A(1, 1)*A(2, 2)*A(3, 3)
     1     +A(1, 2)*A(2, 3)*A(3, 1)
     2     +A(1, 3)*A(3, 2)*A(2, 1)
     3     -A(1, 3)*A(3, 1)*A(2, 2)
     4     -A(2, 3)*A(3, 2)*A(1, 1)
     5     -A(1, 2)*A(2, 1)*A(3, 3)

      return
      end

C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
C                      Cross Product
C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      subroutine cross_product(a,b,c)

      Real(8), dimension(3) :: a, b
      Real(8), dimension(3) :: c
        
      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - a(2)*b(1)

      return
      end

C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
C                       Voigt Mapping
C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      subroutine voigt(A,Av)

      Real(8), dimension(3,3) :: A
      Real(8), dimension(6)   :: Av
        
      Av(1) = A(1,1)
      Av(2) = A(2,2)
      Av(3) = A(3,3)
      Av(4) = A(1,2)
      Av(5) = A(1,3)
      Av(6) = A(2,3)

      return
      end

C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
C                       Kronecker Product
C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      subroutine kronecker_product(a,b,c)
      
      Real(8), dimension(3) :: a, b
      Real(8), dimension(3,3) :: c

      c = 0.D0
      do i = 1,3
            do j = 1,3
                  c(i,j) = a(i)*b(j)
            enddo
      enddo
        
      return
      end   

C================================================================
C Inverse of a 3*3 matrix
C================================================================
      subroutine MatrixInverse(A,B)

      Real(8) :: D
      Real(8), dimension(3,3) :: A, B

      B = 0.D0
C     Determinant of A
      call determinant(A, D)

      if (D.EQ.0.D0) then
            print*,'Error: Determinant is zero for calculation of inverse'
      endif

C     Inverse of A
      B(1,1) = (A(2,2)*A(3,3)-A(2,3)*A(3,2))/D
      B(1,2) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))/D
      B(1,3) = (A(1,2)*A(2,3)-A(1,3)*A(2,2))/D
      B(2,1) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))/D
      B(2,2) = (A(1,1)*A(3,3)-A(1,3)*A(3,1))/D
      B(2,3) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))/D
      B(3,1) = (A(2,1)*A(3,2)-A(2,2)*A(3,1))/D
      B(3,2) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))/D
      B(3,3) = (A(1,1)*A(2,2)-A(1,2)*A(2,1))/D

      return
      end
C==================
