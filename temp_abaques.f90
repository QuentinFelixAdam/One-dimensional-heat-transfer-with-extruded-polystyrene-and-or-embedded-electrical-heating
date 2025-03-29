PROGRAM temp
IMPLICIT NONE
INTEGER :: Nz,stat,i,j,indexl1,indexl2,nlinescase,m
REAL :: dt, dz, t1, t2,alpha,z_R,Delta_R
REAL :: rho, c,k, sigma, eps_s, a,L,hc
REAL:: power,indice_R,rho1
REAL::rho2,L1,L2,rdz2,intercept
REAL :: rk,c_eff_av_1,c_eff_av_2
REAL :: k_eff_2,rho_av_1,rho_av_2
REAL :: k_int_1,k_int_2_m,k_int_2_p,c_int_1,c_int_2_m,c_int_2_p,powerfactor
REAL :: exponentlatent,latentheat
REAL :: c2unfrozen, c2frozen, k2unfrozen, k2frozen
REAL :: theta_uw,beta,temp_ref,theta_0_1
REAL :: d_T_theta_uw,c_eff_1,k_eff_1,c_eff_2
REAL :: c1unfrozen, c1frozen, k1unfrozen, k1frozen,theta_0_2
REAL,DIMENSION(:),ALLOCATABLE::f_s,T_air,w_R,Told,Tnew
REAL,DIMENSION(:),ALLOCATABLE::V_air,fD_loop,fD,n_factor
REAL :: c_1, rho_1,k_1
REAL :: rho_2, c_2_u,c_2_f,k_2_u,k_2_f,theta_tot_2
REAL :: rho_3, c_3_u,c_3_f,k_3_cond,k_3_conv,theta_tot_3
REAL :: rho_4, c_4_u,c_4_f,k_4_u,k_4_f,theta_tot_4
REAL :: rho_5, c_5_u,c_5_f,k_5_u,k_5_f,theta_tot_5
REAL :: rho_6, c_6_u,c_6_f,k_6_u,k_6_f,theta_tot_6
REAL,DIMENSION(:),ALLOCATABLE :: c_eff,k_eff,th_uw,d_T_th_uw
REAL :: rho_avg_1_2,rho_avg_2_3,rho_avg_3_4,rho_avg_4_5,rho_avg_5_6
REAL :: d_1,d_2,d_3,d_4,d_5,d_6,grad_conv,n_factor_winter,n_factor_summer
INTEGER :: Nz_1,Nz_2,Nz_3,Nz_4,Nz_5,Nz_6,countdays
CHARACTER(LEN=21) :: Filename
INTEGER :: m_loop1, m_loop2, m_loop3, m_loop4, Nz_large
REAL :: dz_large,rdz2_large, compteur_prof
INTEGER :: Nz_prof, m_loop5, n_w_to_n_s, n_s_to_n_w
REAL :: kappa_1_u, kappa_1_f, kappa_2_u, kappa_2_f
REAL :: kappa_3_u, kappa_3_f, kappa_4_u, kappa_4_f
REAL :: kappa_5_u, kappa_5_f, kappa_6_u, kappa_6_f
REAL :: k_dry_1, k_dry_2, k_dry_3, k_dry_4, k_dry_5, k_dry_6
REAL :: ksat_1_u, ksat_1_f, ksat_2_u, ksat_2_f
REAL :: ksat_3_u, ksat_3_f, ksat_4_u, ksat_4_f
REAL :: ksat_5_u, ksat_5_f, ksat_6_u, ksat_6_f
REAL :: rho_dry_1, rho_dry_2, rho_dry_3
REAL :: rho_dry_4, rho_dry_5, rho_dry_6
REAL :: beta_1,beta_2,beta_3,beta_4,beta_5,beta_6
INTEGER :: m_loop6, index_FGA, m_loop_spinup
REAL,DIMENSION(:),ALLOCATABLE :: temp_large_depth, temp_shallow_depth
REAL :: Rayleigh_coeff,perma_3,capa_3,condu_3
REAL :: condu_3_rad, air_thermal_exp, air_heat_capa,air_viscosity_k
REAL :: temp_avg_FGA, condu_3_eq, Air_avg, Air_amp
REAL,DIMENSION(:),ALLOCATABLE :: max_depth_FF, Rayleigh_list
REAL,DIMENSION(:),ALLOCATABLE :: k_eff_FGA
INTEGER :: conv_active
REAL,DIMENSION(:),ALLOCATABLE :: list_temp_moy, list_temp_amp
REAL,DIMENSION(:),ALLOCATABLE :: list_FI, list_power
REAL :: d_10_FGA,power_acc

NAMELIST  /para/ dt,Nz,rho,c,k,L,z_R,Delta_R,sigma,eps_s,a,&
               power,rho1,rho2,z_R,L1,L2,c2unfrozen,c2frozen,&
               k2unfrozen,k2frozen,beta,temp_ref,theta_0_1,c1unfrozen,&
               c1frozen,k1unfrozen,k1frozen,latentheat,theta_0_2,&
               c_1, rho_1,k_1,rho_2, c_2_u,c_2_f,k_2_u,k_2_f,theta_tot_2,&
               rho_3, c_3_u,c_3_f,k_3_cond,k_3_conv,theta_tot_3,&
               rho_4, c_4_u,c_4_f,k_4_u,k_4_f,theta_tot_4,&
               rho_5, c_5_u,c_5_f,k_5_u,k_5_f,theta_tot_5,&
               d_1,d_2,d_3,d_4,d_5,grad_conv

OPEN(20,FILE='para.dat')

READ(20,NML=para)

CLOSE(20)

CALL CPU_TIME(t1)

PRINT*,'-----------------------------------'
PRINT*,'------SIMULATION HAS STARTED-------'
PRINT*,'-----------------------------------'

nlinescase = INT(REAL(365.0*86400.0)/REAL(dt))

!nlinescase = 100

ALLOCATE(n_factor(nlinescase))
ALLOCATE(max_depth_FF(nlinescase))
ALLOCATE(Rayleigh_list(nlinescase))
ALLOCATE(k_eff_FGA(nlinescase))


ALLOCATE(list_temp_moy(10))
ALLOCATE(list_temp_amp(10))
ALLOCATE(list_FI(10))

OPEN(36,FILE='Temp_moyenne.dat',STATUS='old')

DO i = 1,10
   READ(36,*) list_temp_moy(i)
END DO

CLOSE(36)

OPEN(37,FILE='Temp_amp.dat',STATUS='old')

DO i = 1,10
   READ(37,*) list_temp_amp(i)
END DO

CLOSE(37)

OPEN(38,FILE='Indices_de_gel.dat',STATUS='old')

DO i = 1,10
   READ(38,*) list_FI(i)
END DO

CLOSE(38)

beta_1 = 0.1
beta_2 = 0.1
beta_3 = 0.1
beta_4 = 0.1
beta_5 = 0.1

! DEFINITION PROPRIETES THERMIQUES COUCHE 2 -- MG20

kappa_2_f = 1.9
kappa_2_u = 4.6

! DEFINITION PROPRIETES THERMIQUES COUCHE 4 -- MG112

kappa_4_f = 0.95
kappa_4_u = 3.55

! DEFINITION PROPRIETES THERMIQUES COUCHE 5 -- ARGILE

kappa_5_f = 0.85
kappa_5_u = 1.9

dz = L/(Nz-1)
alpha = k/(rho*c)
rdz2 = 1.0/(dz**2.0)
rk = 1.0/k
dz_large = 0.1
rdz2_large = 1.0/(dz_large**2.0)
Nz_large = INT(REAL(16.0)/REAL(dz_large)) + 1

ALLOCATE(w_r(Nz))

DO i = 1,Nz

w_R(i) = 0.0

END DO

!indice_R = INT(0.05/dz) + 1
w_R(INT(z_R/dz) + 1) = 1.0

rho_avg_1_2 = (rho_1 + rho_2) * 0.5
rho_avg_2_3 = (rho_2 + rho_3) * 0.5
rho_avg_3_4 = (rho_3 + rho_4) * 0.5
rho_avg_4_5 = (rho_4 + rho_5) * 0.5

alpha = k_1/(rho_1*c_1)

Nz_1 = INT(REAL(d_1)/REAL(dz)) + 1
Nz_2 = INT(REAL(d_1+d_2)/REAL(dz)) + 1
Nz_3 = INT(REAL(d_1+d_2+d_3)/REAL(dz)) + 1
Nz_4= INT(REAL(d_1+d_2+d_3+d_4)/REAL(dz)) + 1
Nz_5= INT(REAL(20.0)/REAL(dz)) + 1

ALLOCATE(Told(Nz))
ALLOCATE(Tnew(Nz))
ALLOCATE(c_eff(Nz))
ALLOCATE(k_eff(Nz))
ALLOCATE(th_uw(Nz))
ALLOCATE(d_T_th_uw(Nz))
ALLOCATE(list_power(15))

DO j = 1,1
   list_power(j) = 0.0
END DO

DO j = 2,15
   list_power(j) = 5.0 + (j-2.0)*2.0
END DO

OPEN(23, IOSTAT=stat, FILE='Surfacetemp.dat', STATUS='old')
IF (stat.EQ.0) THEN
CLOSE(23, STATUS='delete')
END IF

countdays = 0

DO i = 1,500
   !IF((MOD(i-1,0)).EQ.0)THEN
      IF(countdays.LT.10) THEN
         WRITE(FileName,'(A,I1.1,A)') 'Tdepth.',countdays,'.dat'
         OPEN(26, IOSTAT=stat, FILE=FileName, STATUS='old')
         IF (stat.EQ.0) THEN
            CLOSE(26, STATUS='delete')
         END IF
      ELSE IF(countdays.LT.100) THEN
         WRITE(FileName,'(A,I2.2,A)') 'Tdepth.',countdays,'.dat'
         OPEN(26, IOSTAT=stat, FILE=FileName, STATUS='old')
         IF (stat.EQ.0) THEN
            CLOSE(26, STATUS='delete')
         END IF
      ELSE IF(countdays.LT.1000) THEN
         WRITE(FileName,'(A,I3.3,A)') 'Tdepth.',countdays,'.dat'
         OPEN(26, IOSTAT=stat, FILE=FileName, STATUS='old')
         IF (stat.EQ.0) THEN
            CLOSE(26, STATUS='delete')
         END IF
      ELSE IF(countdays.LT.10000) THEN
         WRITE(FileName,'(A,I4.4,A)') 'Tdepth.',countdays,'.dat'
         OPEN(26, IOSTAT=stat, FILE=FileName, STATUS='old')
         IF (stat.EQ.0) THEN
            CLOSE(26, STATUS='delete')
         END IF
      ELSE IF(countdays.LT.100000) THEN
         WRITE(FileName,'(A,I5.5,A)') 'Tdepth.',countdays,'.dat'
         OPEN(26, IOSTAT=stat, FILE=FileName, STATUS='old')
         IF (stat.EQ.0) THEN
            CLOSE(26, STATUS='delete')
         END IF
      END IF
      countdays = countdays + 1
   !END IF
END DO

OPEN(27,FILE='Surfacetemp.dat',POSITION='append')

countdays = 0

! Layer 1 --> Asphalte concrete
! Layer 2 --> MG20
! Layer 3 --> FGA
! Layer 4 --> MG112
! Layer 5 --> Deeper soil

!||||||||| BEGIN: PARAMETERS TO MODIFY ||||||||||||||

rho_1 = 2350.0 ! Bulk density Layer 1 [kg/m**3] Bulk density is equal to dry density here as the volumetric moisture content is equal to 0.

c_1 = 1000.0 ! Heat capacity Layer 1 [J/kg/K]

k_1 = 1.48 ! Thermal conductivity first layer [W/m/K]

d_1 = 0.15 ! Thickness Layer 1 [m]

rho_dry_2 = 2200.0 ! Dry density of Layer 2 [kg/m**3]

theta_tot_2 = 0.088 ! Total volumetric moisture content Layer 2 [--]

d_2 = 0.45 ! Thickness Layer 2 [m]

rho_dry_4 = 1920.0 ! Dry density of Layer 4 [kg/m**3]

theta_tot_4 = 0.1536 ! Total volumetric moisture content Layer 4 [--]

d_4 = 0.30 ! Thickness Layer 1 [m]

rho_dry_5 = 1300.0 ! Dry density of Layer 5 [kg/m**3]

theta_tot_5 = 0.4 ! Total volumetric moisture content Layer 5 [--]

z_R = 0.07

DO i = 1,Nz

w_R(i) = 0.0

END DO

w_R(INT(z_R/dz) + 1) = 1.0

!||||||||| END: PARAMETERS TO MODIFY ||||||||||||||

d_3 = 0.01 - 0.005

d_2 = 0.45

DO m_loop6 = 1,10

Air_amp = list_temp_amp(m_loop6)
Air_avg = list_temp_moy(m_loop6)

d_3 = 0.01 - 0.005

DO m_loop5 = 1,14

d_3 = d_3 + 0.005

DO m_loop4 = 1,1

DO m_loop3 = 1,1

DO m_loop2 = 1,1

DO m_loop1 = 1,1

DO j = 1,Nz
   Tnew(j) = Air_avg
   Told(j) = Air_avg
END DO

IF(countdays.LT.10) THEN
   WRITE(FileName,'(A,I1.1,A)') 'Tdepth.',countdays,'.dat'
ELSE IF(countdays.LT.100) THEN
   WRITE(FileName,'(A,I2.2,A)') 'Tdepth.',countdays,'.dat'
ELSE IF(countdays.LT.1000) THEN
   WRITE(FileName,'(A,I3.3,A)') 'Tdepth.',countdays,'.dat'
ELSE IF(countdays.LT.10000) THEN
   WRITE(FileName,'(A,I4.4,A)') 'Tdepth.',countdays,'.dat'
ELSE IF(countdays.LT.100000) THEN
   WRITE(FileName,'(A,I5.5,A)') 'Tdepth.',countdays,'.dat'
END IF

Nz_1 = INT(REAL(d_1)/REAL(dz)) + 1
Nz_2 = INT(REAL(d_1+d_2)/REAL(dz)) + 1
Nz_3 = INT(REAL(d_1+d_2+d_3)/REAL(dz)) + 1
Nz_4 = INT(REAL(d_1+d_2+d_3+d_4)/REAL(dz)) + 1
Nz_5 = INT(REAL(20.0)/REAL(dz)) + 1

!UPDATE THERMAL PROPERTIES LAYER 2 -- MG20

k_dry_2 = 0.75*10**(-1.2*(1.0-rho_dry_2/2700.0))

ksat_2_f = (2.5**(rho_dry_2/2700.0))*&
(2.24**(1.0-rho_dry_2/2700.0))
ksat_2_u = (2.5**(rho_dry_2/2700.0))*&
(0.6**(1.0-rho_dry_2/2700.0))

k_2_f = ((kappa_2_f*ksat_2_f-k_dry_2)*&
theta_tot_2/(1.0-rho_dry_2/2700.0) + k_dry_2)/&
(1.0 + (kappa_2_f-1.0)*theta_tot_2/(1.0-rho_dry_2/2700.0))

k_2_u = ((kappa_2_u*ksat_2_u-k_dry_2)*&
theta_tot_2/(1.0-rho_dry_2/2700.0) + k_dry_2)/&
(1.0 + (kappa_2_u-1.0)*theta_tot_2/(1.0-rho_dry_2/2700.0))

c_2_f = 800.0*2700.0*(rho_dry_2/2700.0) + 2050.0*900.0*theta_tot_2*1000.0/2700.0
c_2_u = 800.0*2700.0*(rho_dry_2/2700.0) + 4180.0*1000.0*theta_tot_2*1000.0/2700.0

c_2_f = c_2_f/(rho_dry_2)
c_2_u = c_2_u/(rho_dry_2)

rho_2 = 0.5*(rho_dry_2+(rho_dry_2*rho_dry_2+4000.0*theta_tot_2*rho_dry_2)**(0.5))

!UPDATE THERMAL PROPERTIES LAYER 3 -- XPS

rho_dry_3 = 50.0

rho_3 = rho_dry_3

!UPDATE THERMAL PROPERTIES LAYER 4 -- MG1120

k_dry_4 = 0.75*10**(-1.2*(1.0-rho_dry_4/2700.0))

ksat_4_f = (2.5**(rho_dry_4/2700.0))*&
(2.24**(1.0-rho_dry_4/2700.0))
ksat_4_u = (2.5**(rho_dry_4/2700.0))*&
(0.6**(1.0-rho_dry_4/2700.0))

k_4_f = ((kappa_4_f*ksat_4_f-k_dry_4)*&
theta_tot_4/(1.0-rho_dry_4/2700.0) + k_dry_4)/&
(1.0 + (kappa_4_f-1.0)*theta_tot_4/(1.0-rho_dry_4/2700.0))

k_4_u = ((kappa_4_u*ksat_4_u-k_dry_4)*&
theta_tot_4/(1.0-rho_dry_4/2700.0) + k_dry_4)/&
(1.0 + (kappa_4_u-1.0)*theta_tot_4/(1.0-rho_dry_4/2700.0))

c_4_f = 800.0*2700.0*(rho_dry_4/2700.0) + 2050.0*900.0*theta_tot_4*1000.0/2700.0
c_4_u = 800.0*2700.0*(rho_dry_4/2700.0) + 4180.0*1000.0*theta_tot_4*1000.0/2700.0

c_4_f = c_4_f/(rho_dry_4)
c_4_u = c_4_u/(rho_dry_4)

rho_4 = 0.5*(rho_dry_4+(rho_dry_4*rho_dry_4+4000.0*theta_tot_4*rho_dry_4)**(0.5))

!UPDATE THERMAL PROPERTIES LAYER 5 -- ARGILE

k_dry_5 = 0.75*10**(-1.2*(1.0-rho_dry_5/2700.0))

ksat_5_f = (2.5**(rho_dry_5/2700.0))*&
(2.24**(1.0-rho_dry_5/2700.0))
ksat_5_u = (2.5**(rho_dry_5/2700.0))*&
(0.5**(1.0-rho_dry_5/2700.0))

k_5_f = ((kappa_5_f*ksat_5_f-k_dry_5)*&
theta_tot_5/(1.0-rho_dry_5/2700.0) + k_dry_5)/&
(1.0 + (kappa_5_f-1.0)*theta_tot_5/(1.0-rho_dry_5/2700.0))

k_5_u = ((kappa_5_u*ksat_5_u-k_dry_5)*&
theta_tot_5/(1.0-rho_dry_5/2700.0) + k_dry_5)/&
(1.0 + (kappa_5_u-1.0)*theta_tot_5/(1.0-rho_dry_5/2700.0))

c_5_f = 800.0*2700.0*(rho_dry_5/2700.0) + 2050.0*900.0*theta_tot_5*1000.0/2700.0
c_5_u = 800.0*2700.0*(rho_dry_5/2700.0) + 4180.0*1000.0*theta_tot_5*1000.0/2700.0

c_5_f = c_5_f/(rho_dry_5)
c_5_u = c_5_u/(rho_dry_5)

rho_5 = 0.5*(rho_dry_5+(rho_dry_5*rho_dry_5+4000.0*theta_tot_5*rho_dry_5)**(0.5))

PRINT*,'-----------------------------------'
PRINT*,'-------DEBUT CAS OPTI NUMERO = ',countdays
PRINT*,'-----------------------------------'


rho_avg_1_2 = (rho_1 + rho_2) * 0.5
rho_avg_2_3 = (rho_2 + rho_3) * 0.5
rho_avg_3_4 = (rho_3 + rho_4) * 0.5
rho_avg_4_5 = (rho_4 + rho_5) * 0.5

!OPEN(27,FILE='Surfacetemp.dat',POSITION='append')
!WRITE(27,*) countdays,MINVAL(max_depth_FF)
!CLOSE(27)

DO m_loop_spinup = 1,15

IF(m_loop_spinup.EQ.15) THEN
DO j = 1,nlinescase
   max_depth_FF(j) = 0.0
END DO
power_acc = 0.0
END IF

DO i = 1,nlinescase

dz = L/(Nz-1)
powerfactor = 1.0/(c_1*rho_1*dz*Delta_R)

intercept = 0.0

IF(m_loop_spinup.LT.15) THEN
   power = 0.0
END IF

!IF(0.EQ.0) THEN
IF(m_loop_spinup.EQ.15) THEN
IF(((MOD(i-1,INT(86400.0/dt))).EQ.0)) THEN
   DO j = 1,INT(4.0/dz)
      IF(Tnew(j)*Tnew(j+1).LT.0.0) THEN
         intercept = -1.0*(j-1)*dz
         !WRITE(26,*)(m-1)*(nlinescase-1)*dt + (i-1)*dt,intercept
		 max_depth_FF(i) = intercept
      END IF
   END DO
   !WRITE(26,*)(m-1)*(nlinescase-1)*dt + (i-1)*dt,intercept
   max_depth_FF(i) = intercept
END IF
END IF
!END IF

IF(m_loop_spinup.EQ.15) THEN
   IF(i.GE.INT(131040*60.0/dt)) THEN
   IF(((MOD(i-1,INT(86400.0/dt))).EQ.0)) THEN
      IF(intercept.LT.-1.0*(d_1+d_2+d_3)) THEN
	     IF(power.EQ.0.0) THEN
	        power = 40.0
	     ELSE
		    power = 40.0
		 END IF
		 power_acc = power_acc + 1.0
      ELSE
	     IF(power.GT.0.0) THEN
	        power = 0.0
	     ELSE
		    power = 0.0
		 END IF
	  END IF
	  !PRINT*,intercept,power
   END IF
   END IF
END IF

IF(0.EQ.0) THEN
IF(m_loop_spinup.EQ.15) THEN
IF(((MOD(i-1,INT(86400.0/dt))).EQ.0)) THEN
   OPEN(26,FILE=Filename,POSITION='append')
   WRITE(26,*) (i-1)*dt,intercept,&
		 Air_avg,Air_amp,d_3,&
		 list_FI(m_loop6),REAL(40.0)/REAL(Delta_R),&
		 d_1,d_2,d_4,REAL(power_acc), Tnew(1)
    CLOSE(26)
END IF
END IF
END IF

! UPDATE SURFACE TEMPERATURE ON THE ASPHALT

   c_eff(1) = c_1
   k_eff(1) = k_1
   th_uw(1) = 0.0
   d_T_th_uw(1) = 0.0

   IF (Air_avg + Air_amp*(SIN((i-1)*dt*1.99238499086111*10**(-7.0))).LT.0.0) THEN
      Tnew(1) = (Air_avg + Air_amp*(SIN((i-1)*dt*1.99238499086111*10**(-7.0))))*0.9
      Told(1) = (Air_avg + Air_amp*(SIN((i-1)*dt*1.99238499086111*10**(-7.0))))*0.9
   ELSE
      Tnew(1) = (Air_avg + Air_amp*(SIN((i-1)*dt*1.99238499086111*10**(-7.0))))*1.5
      Told(1) = (Air_avg + Air_amp*(SIN((i-1)*dt*1.99238499086111*10**(-7.0))))*1.5
   END IF

!LOOP TO UPDATE TEMPERATURE WITHIN LAYER 1
!FREEZE / THAWING EFFECTS ARE NEGLECTED

   DO j = 2,Nz_1 - 1
   
      c_eff(j) = c_1
      k_eff(j) = k_1
      th_uw(j) = 0.0
      d_T_th_uw(j) = 0.0

      Tnew(j) = Told(j) + alpha*dt*(Told(j-1)+ &
      Told(j+1))*rdz2 -2.0*alpha*dt*rdz2*Told(j) + &
      dt*w_R(j)*power*powerfactor
	  

   END DO

!UPDATE TEMPERATURE OF INTERFACE 1_2 ASPHALT TO MG20

   IF(Told(Nz_1).GE.temp_ref) THEN
         d_T_th_uw(Nz_1) = 0.0
   ELSE
         d_T_th_uw(Nz_1) = ((temp_ref-Told(Nz_1))/(temp_ref+273.15))**beta
         d_T_th_uw(Nz_1) = (th_uw(Nz_1-1)+th_uw(Nz_1+1))*beta*d_T_th_uw(Nz_1)*0.5/&
                        (temp_ref - Told(Nz_1))
   END IF

   c_eff(Nz_1) = ( c_eff(Nz_1-1)+ c_eff(Nz_1+1)) * 0.5
   k_eff(Nz_1) = ( k_eff(Nz_1-1)+ k_eff(Nz_1+1)) * 0.5
   
   Tnew(Nz_1) = Told(Nz_1) + dt*rdz2*(k_eff(Nz_1+1) *Told(Nz_1+1) - &
   (k_eff(Nz_1-1)  + k_eff(Nz_1+1) )*Told(Nz_1) + k_eff(Nz_1-1) *&
   Told(Nz_1-1))/(c_eff(Nz_1)* rho_avg_1_2 + latentheat*d_T_th_uw(Nz_1))


!LOOP TO UPDATE TEMPERATURE WITHIN LAYER 2 -- MG20

   !PRINT*,Tnew(Nz_1+1)

   DO j = Nz_1+1,Nz_2-1

      IF(Told(j).GE.temp_ref) THEN
         th_uw(j) = theta_tot_2
         d_T_th_uw(j) = 0.0
      ELSE
         th_uw(j) = ((temp_ref-Told(j))/(temp_ref+273.15))**beta_2
         th_uw(j) = theta_tot_2*(1.0 - th_uw(j))
         d_T_th_uw(j)  = ((temp_ref-Told(j))/(temp_ref+273.15))**beta_2
         d_T_th_uw(j)  = theta_tot_2*beta_2*d_T_th_uw(j) /&
                        (temp_ref - Told(j))
      END IF
	  
	  IF(d_T_th_uw(j).GT.10.0**(30.0)) THEN
	     d_T_th_uw(j) = 10.0**(30.0)
	  END IF

      exponentlatent = th_uw(j)/theta_tot_2

      c_eff(j) = c_2_f*(1.0-exponentlatent) + &
      c_2_u*exponentlatent
      
      k_eff(j) = k_2_f*(1.0-exponentlatent) + &
      k_2_u*exponentlatent

      Tnew(j) = Told(j) + k_eff(j)*dt*rdz2*(Told(j+1) - 2.0*Told(j) + &
      Told(j-1))/(c_eff(j)*rho_2 + latentheat*d_T_th_uw(j))

   END DO
   
!UPDATE TEMPERATURE OF INTERFACE 2_3 MG20 to XPS

   IF(Told(Nz_2).GE.temp_ref) THEN
         d_T_th_uw(Nz_2) = 0.0
   ELSE
         d_T_th_uw(Nz_2) = ((temp_ref-Told(Nz_2))/(temp_ref+272.15))**(beta_2+beta_3)*0.5
         d_T_th_uw(Nz_2) = (th_uw(Nz_2-1)+th_uw(Nz_2+1))*(beta_2+beta_3)*0.5*d_T_th_uw(Nz_2)*0.5/&
                        (temp_ref - Told(Nz_2))
   END IF

   c_eff(Nz_2) = ( c_eff(Nz_2-1)+ c_eff(Nz_2+1)) * 0.5
   k_eff(Nz_2) = ( k_eff(Nz_2-1)+ k_eff(Nz_2+1)) * 0.5
   
   Tnew(Nz_2) = Told(Nz_2) + dt*rdz2*(k_eff(Nz_2+1) *Told(Nz_2+1) - &
   (k_eff(Nz_2-1)  + k_eff(Nz_2+1) )*Told(Nz_2) + k_eff(Nz_2-1) *&
   Told(Nz_2-1))/(c_eff(Nz_2)* rho_avg_2_3 + latentheat*d_T_th_uw(Nz_2))
 
   
!LOOP TO UPDATE TEMPERATURE WITHIN LAYER 3 -- XPS

   DO j = Nz_2+1,Nz_3-1
	  
	  !IF(Told(j).GE.temp_ref) THEN
      !   th_uw(j) = theta_tot_3
      !   d_T_th_uw(j) = 0.0
      !ELSE
      !   th_uw(j) = ((temp_ref-Told(j))/(temp_ref+273.15))**beta_3
      !   th_uw(j) = theta_tot_3*(1.0 - th_uw(j))
      !   d_T_th_uw(j)  = ((temp_ref-Told(j))/(temp_ref+273.15))**beta_3
      !   d_T_th_uw(j)  = theta_tot_3*beta_3*d_T_th_uw(j) /&
      !                  (temp_ref - Told(j))
      !END IF	

	  !IF(d_T_th_uw(j).GT.10.0**(30.0)) THEN
	  !   d_T_th_uw(j) = 10.0**(30.0)
	  !END IF

	 th_uw(j) = 0.0  
	 d_T_th_uw(j) = 0.0
	 
	 !exponentlatent = th_uw(j)/theta_tot_3
	 !exponentlatent = 0.0

     c_eff(j) = c_3_f
	 
	 k_eff(j) = k_3_cond

     Tnew(j) = Told(j) + k_eff(j)*dt*rdz2*(Told(j+1) - 2.0*Told(j) + &
     Told(j-1))/(c_eff(j)*rho_3 + latentheat*d_T_th_uw(j))
	  

   END DO
   
!UPDATE TEMPERATURE OF INTERFACE 3_4 -- XPS TO MG112

   IF(Told(Nz_3).GE.temp_ref) THEN
         d_T_th_uw(Nz_3) = 0.0
   ELSE
         d_T_th_uw(Nz_3) = ((temp_ref-Told(Nz_3))/(temp_ref+273.15))**(beta_3+beta_4)*0.5
         d_T_th_uw(Nz_3) = (th_uw(Nz_3-1)+th_uw(Nz_3+1))*(beta_3+beta_4)*0.5*d_T_th_uw(Nz_3)*0.5/&
                        (temp_ref - Told(Nz_3))
   END IF

   c_eff(Nz_3) = ( c_eff(Nz_3-1)+ c_eff(Nz_3+1)) * 0.5
   k_eff(Nz_3) = ( k_eff(Nz_3-1)+ k_eff(Nz_3+1)) * 0.5
   
   Tnew(Nz_3) = Told(Nz_3) + dt*rdz2*(k_eff(Nz_3+1) *Told(Nz_3+1) - &
   (k_eff(Nz_3-1)  + k_eff(Nz_3+1) )*Told(Nz_3) + k_eff(Nz_3-1) *&
   Told(Nz_3-1))/(c_eff(Nz_3)* rho_avg_3_4 + latentheat*d_T_th_uw(Nz_3))
   
!LOOP TO UPDATE TEMPERATURE WITHIN LAYER 4 -- MG112

   DO j = Nz_3+1,Nz_4-1

      IF(Told(j).GE.temp_ref) THEN
         th_uw(j) = theta_tot_4
         d_T_th_uw(j) = 0.0
      ELSE
         th_uw(j) = ((temp_ref-Told(j))/(temp_ref+273.15))**beta_4
         th_uw(j) = theta_tot_4*(1.0 - th_uw(j))
         d_T_th_uw(j)  = ((temp_ref-Told(j))/(temp_ref+273.15))**beta_4
         d_T_th_uw(j)  = theta_tot_4*beta_4*d_T_th_uw(j) /&
                        (temp_ref - Told(j))
      END IF
	  
	  IF(d_T_th_uw(j).GT.10.0**(30.0)) THEN
	     d_T_th_uw(j) = 10.0**(30.0)
	  END IF

      exponentlatent = th_uw(j)/theta_tot_4

      c_eff(j) = c_4_f*(1.0-exponentlatent) + &
      c_4_u*exponentlatent
      
      k_eff(j) = k_4_f*(1.0-exponentlatent) + &
      k_4_u*exponentlatent

      Tnew(j) = Told(j) + k_eff(j)*dt*rdz2*(Told(j+1) - 2.0*Told(j) + &
      Told(j-1))/(c_eff(j)*rho_4 + latentheat*d_T_th_uw(j))

   END DO
   
!UPDATE TEMPERATURE OF INTERFACE 4_5 -- MG112 TO ARGILE

   IF(Told(Nz_4).GE.temp_ref) THEN
         d_T_th_uw(Nz_4) = 0.0
   ELSE
         d_T_th_uw(Nz_4) = ((temp_ref-Told(Nz_4))/(temp_ref+273.14))**(beta_4+beta_5)*0.5
         d_T_th_uw(Nz_4) = (th_uw(Nz_4-1)+th_uw(Nz_4+1))*(beta_4+beta_5)*0.5*d_T_th_uw(Nz_4)*0.5/&
                        (temp_ref - Told(Nz_4))
   END IF

   c_eff(Nz_4) = ( c_eff(Nz_4-1)+ c_eff(Nz_4+1)) * 0.5
   k_eff(Nz_4) = ( k_eff(Nz_4-1)+ k_eff(Nz_4+1)) * 0.5
   
   Tnew(Nz_4) = Told(Nz_4) + dt*rdz2*(k_eff(Nz_4+1) *Told(Nz_4+1) - &
   (k_eff(Nz_4-1)  + k_eff(Nz_4+1) )*Told(Nz_4) + k_eff(Nz_4-1) *&
   Told(Nz_4-1))/(c_eff(Nz_4)* rho_avg_4_5 + latentheat*d_T_th_uw(Nz_4))
   
!LOOP TO UPDATE TEMPERATURE WITHIN LAYER 5 -- ARGILE

   DO j = Nz_4+1,Nz_5-1

      IF(Told(j).GE.temp_ref) THEN
         th_uw(j) = theta_tot_5
         d_T_th_uw(j) = 0.0
      ELSE
         th_uw(j) = ((temp_ref-Told(j))/(temp_ref+273.15))**beta_5
         th_uw(j) = theta_tot_5*(1.0 - th_uw(j))
         d_T_th_uw(j)  = ((temp_ref-Told(j))/(temp_ref+273.15))**beta_5
         d_T_th_uw(j)  = theta_tot_5*beta_5*d_T_th_uw(j) /&
                        (temp_ref - Told(j))
      END IF
	  
	  IF(d_T_th_uw(j).GT.10.0**(30.0)) THEN
	     d_T_th_uw(j) = 10.0**(30.0)
	  END IF

      exponentlatent = th_uw(j)/theta_tot_5

      c_eff(j) = c_5_f*(1.0-exponentlatent) + &
      c_5_u*exponentlatent
      
      k_eff(j) = k_5_f*(1.0-exponentlatent) + &
      k_5_u*exponentlatent

      Tnew(j) = Told(j) + k_eff(j)*dt*rdz2*(Told(j+1) - 2.0*Told(j) + &
      Told(j-1))/(c_eff(j)*rho_5 + latentheat*d_T_th_uw(j))

   END DO

!UPDATE TEMPERATURE AT LAST NODE

   Tnew(Nz_5) = Tnew(Nz_5-1)

!SWAP THE TWO VECTORS

   Told = Tnew

END DO ! END FOR nlinescase

PRINT*, "LOOPING YEAR ",m_loop_spinup," WITHIN CASE NUMBER ",countdays

END DO ! END FOR m_loop_spinup

countdays = countdays + 1

END DO ! END FOR m_loop1

END DO ! END FOR m_loop2

END DO ! END FOR m_loop3

END DO ! END FOR m_loop4

END DO ! END FOR m_loop5

END DO ! END FOR m_loop6

CALL CPU_TIME(t2)

CLOSE(27)

END PROGRAM temp
