! This MODULE contains the following FUNCTIONS:
!
! FUNCTION Density(Z_Atomic,r,theta,phi)
! This function is the density of some closed shell atoms I have been using
! It is obtained based on the fittings James, myself and some others have done.
!
! FUNCTION FuncTest3D(r,theta,phi) 
! This function is just a random one chosen to check if the 
! 3D integrating routine works
!
! FUNCTION FuncTest6D(r1,theta1,phi1,r2,theta2,phi2) 
! This function is just a random one chosen to check if the 
! 6D integrating routine works
!


MODULE DensityFunction

CONTAINS


!========================================================================
!========================================================================

! FUNCTION Density(Z_Atomic,r,theta,phi)
! This function is the density of some closed shell atoms I have been using
! It is obtained based on the fittings James, myself and some others have done.
!

FUNCTION Density(Z_Atomic,r,theta,phi) 

USE nrtype		! It uses the module nrtype where the type of used variables is declared


IMPLICIT NONE


INTEGER, INTENT(IN) :: Z_Atomic			!(input) Atomic number
REAL(dp), INTENT(IN) :: r, theta, phi	!(inputs) r, theta, phi coordiantes where the density is
										! evaluated
REAL(dp) :: Coeffs(1:20), Expon(1:20)	! Arrays where the coefficients and exponents of the
										! Gaussian expansion are kept.
REAL(dp) :: density_aux(1:20)			! Auxiliary array to compute the function
REAL(dp) :: density						! The value of the density
INTEGER :: i, j							! Just counters


! The following IF condition is used to choose the coefficients and exponents
! according to the atomic number.
! Since the arrays Coeffs and Expon have been set to be of dimension 20 in some atoms,
! where less than that number of parameters is available, some of them are set to zero.
! A complete set of coefficients and exponents for a number of atoms can be found
! in the file "rho_G_Ayers.txt"

IF (Z_Atomic == 1) THEN

Coeffs(1) = 1.6974473253E-02	 
Expon(1) = 1.1312807919E+02
Coeffs(2) = 5.3025697740E-02	 
Expon(2) = 6.6470125536E+00
Coeffs(3) = 3.0592022212E-02	 
Expon(3) = 5.3772337038E-01
Coeffs(4) = 5.8707537674E-02	 
Expon(4) = 1.2094177622E+00
Coeffs(5) = 3.9342138430E-03	 
Expon(5) = 3.5105124063E+03
Coeffs(6) = 6.1214629017E-03	 
Expon(6) = 9.5506080498E+02
Coeffs(7) = 3.8470423701E-02	 
Expon(7) = 1.6399410871E+01
Coeffs(8) = 4.7015082653E-03	 
Expon(8) = 2.3950580176E-01
Coeffs(9) = 4.2305722453E-03	 
Expon(9) = 4.8487238729E+08
Coeffs(10) = 6.3961685151E-02	 
Expon(10) = 2.7908454796E+00
Coeffs(11) = 1.0719169010E-02	 
Expon(11) = 3.2100441066E+02
Coeffs(12) = 2.6066777824E-02	 
Expon(12) = 4.2095237219E+01
Coeffs(13) = 8.0914384550E-04	 
Expon(13) = 1.4706017501E+03


DO i = 14, 20
	Coeffs(i) = 0.d0
	Expon(i) = 1.d0
END DO


ELSE IF (Z_Atomic == 2) THEN

Coeffs(1) = 3.72587013930776000E-01	
Expon(1) = 3.22049584257601000E+01
Coeffs(2) = 3.03125226640509000E-01	
Expon(2) = 5.92862791831091000E+01
Coeffs(3) = 3.55156193642548000E-01	
Expon(3) = 3.23706355904431000E+00
Coeffs(4) = 6.81872637173492000E-02	
Expon(4) = 1.61835956608224000E+03
Coeffs(5) = 2.62709918168175000E-02	
Expon(5) = 1.35094660425194000E+04
Coeffs(6) = 9.62863544507734000E-02	
Expon(6) = 1.03045514096093000E+00
Coeffs(7) = 3.77920862379295000E-02	
Expon(7) = 8.61567791157559000E+04
Coeffs(8) = 2.22337095698631000E-01	
Expon(8) = 1.83287991704299000E+00
Coeffs(9) = 1.30520673406063000E-01	
Expon(9) = 4.23862904237793000E+02
Coeffs(10) = 2.24296355532833000E-02	
Expon(10) = 5.73160579903299000E-01
Coeffs(11) = 9.37424202617610000E-02	
Expon(11) = 8.26383780782377000E+02
Coeffs(12) = 2.39218920951108000E-01	
Expon(12) = 1.11827115869475000E+02
Coeffs(13) = 1.52543877042631000E-03	
Expon(13) = 3.08911634633236000E-01
Coeffs(14) = 4.33436114721961000E-01	
Expon(14) = 1.78294322243255000E+01
Coeffs(15) = 1.95686993455118000E-02	
Expon(15) = 2.96150976908799000E+04
Coeffs(16) = 4.42597495414955000E-01	
Expon(16) = 5.69232235232830000E+00
Coeffs(17) = 3.60435849926551000E-02	
Expon(17) = 6.49084148432205000E+03
Coeffs(18) = 1.80913727074189000E-01	
Expon(18) = 2.16682704507575000E+02
Coeffs(19) = 4.95652456223779000E-02	
Expon(19) = 3.21124214858286000E+03
Coeffs(20) = 4.64613861668947000E-01	
Expon(20) = 1.00265876351762000E+01


ELSE IF (Z_Atomic == 3) THEN

 Coeffs(1) = 1.73910306538455000E+00	
 Expon(1) = 1.43720065982558000E+01
 Coeffs(2) = 4.68131319934221000E-01	
 Expon(2) = 9.47115407405902000E+02
 Coeffs(3) = 9.32335316863302000E-01	
 Expon(3) = 4.68404489605021000E+00
 Coeffs(4) = 1.06095520157587000E+00	
 Expon(4) = 1.50180919767390000E+02
 Coeffs(5) = 7.48828729260398000E-02	
 Expon(5) = 6.97720209010184000E+04
 Coeffs(6) = 1.34593678495505000E+00	
 Expon(6) = 8.22406324879563000E+01
 Coeffs(7) = 2.59141396717078000E-01	
 Expon(7) = 3.48212688033284000E+03
 Coeffs(8) = 1.44482604658001000E-01	
 Expon(8) = 1.41870295819931000E+04
 Coeffs(9) = 5.56018291363616000E-03	
 Expon(9) = 8.29226753237294000E-01
 Coeffs(10) = 1.93319917802891000E-01	
 Expon(10) = 6.87781643828071000E+03
 Coeffs(11) = 1.04020652065583000E-01	
 Expon(11) = 3.07166124044656000E+04
 Coeffs(12) = 1.77314403758160000E+00	
 Expon(12) = 2.53806527752997000E+01
 Coeffs(13) = 4.08640248217337000E-01	
 Expon(13) = 2.67386831009337000E+00
 Coeffs(14) = 8.11978015398148000E-01	
 Expon(14) = 2.75327500805785000E+02
 Coeffs(15) = 9.18810498055641000E-02	
 Expon(15) = 1.51575486153528000E+00
 Coeffs(16) = 3.50628448503757000E-01	
 Expon(16) = 1.80098086731805000E+03
 Coeffs(17) = 1.38245998390784000E-01	
 Expon(17) = 2.15905017006394000E+05
 Coeffs(18) = 1.44243672851088000E+00	
 Expon(18) = 8.19601785466215000E+00
 Coeffs(19) = 6.16264530266749000E-01	
 Expon(19) = 5.07267465132135000E+02
 Coeffs(20) = 1.61224197560743000E+00	
 Expon(20) = 4.53855482854463000E+01


ELSE IF (Z_Atomic == 4) THEN


Coeffs(1) = 2.16196602500000000E-02	
Expon(1) = 1.71317162000000000E-01
Coeffs(2) = 1.68555951400000000E-02	
Expon(2) = 3.77138325500000000E-01
Coeffs(3) = 7.42880343800000000E-01	
Expon(3) = 4.02346344500000000E+00
Coeffs(4) = 4.00198500800000000E+00	
Expon(4) = 8.85726945700000000E+00
Coeffs(5) = 5.57284805300000000E+00	
Expon(5) = 1.94984304700000000E+01
Coeffs(6) = 7.02282986600000000E+00	
Expon(6) = 4.29239273700000000E+01
Coeffs(7) = 4.89091679300000000E+00	
Expon(7) = 9.44929153800000000E+01
Coeffs(8) = 4.18028750500000000E+00	
Expon(8) = 2.08017103900000000E+02
Coeffs(9) = 3.48295521000000000E+00	
Expon(9) = 4.57929732900000000E+02
Coeffs(10) = 1.77367055400000000E-01	
Expon(10) = 1.00808845200000000E+03
Coeffs(11) = 3.49730837400000000E+00	
Expon(11) = 2.21921018600000000E+03

DO i = 12, 20
	Coeffs(i) = 0.d0
	Expon(i) = 1.d0
END DO


ELSE IF (Z_Atomic == 10) THEN


Coeffs(1) = 8.00273577300000000E-02	
Expon(1) = 5.17334281100000000E-01
Coeffs(2) = 5.42362665100000000E-01	
Expon(2) = 1.13575124700000000E+00
Coeffs(3) = 2.90299108300000000E+00	
Expon(3) = 2.49341855200000000E+00
Coeffs(4) = 8.83031519500000000E-01	
Expon(4) = 5.47402971700000000E+00
Coeffs(5) = 8.39379158500000000E+00	
Expon(5) = 2.63834192900000000E+01
Coeffs(6) = 6.99336489000000000E+01	
Expon(6) = 5.79219325600000000E+01
Coeffs(7) = 1.04383394200000000E+02	
Expon(7) = 1.27161314300000000E+02
Coeffs(8) = 1.17391138800000000E+02	
Expon(8) = 2.79168859400000000E+02
Coeffs(9) = 9.38734165100000000E+01	
Expon(9) = 6.12884921200000000E+02
Coeffs(10) = 6.18959986800000000E+01	
Expon(10) = 1.34552230300000000E+03
Coeffs(11) = 6.90819873000000000E+01	
Expon(11) = 2.95394813100000000E+03
Coeffs(12) = 5.60011458600000000E+01	
Expon(12) = 1.42372723800000000E+04
Coeffs(13) = 5.53432860400000000E+00	
Expon(13) = 3.12563857500000000E+04

DO i = 14, 20
	Coeffs(i) = 0.d0
	Expon(i) = 1.d0
END DO


ELSE IF (Z_Atomic == 12) THEN

Coeffs(1) = 0.1056541802E-01    
Expon(1) = 0.1119931034E+00
Coeffs(2) = 0.1356519906E-01    
Expon(2) = 0.2375019538E+00
Coeffs(3) = 0.1174603920E+00    
Expon(3) = 0.1068117428E+01
Coeffs(4) = 0.2110141931E+01    
Expon(4) = 0.2265139268E+01
Coeffs(5) = 0.7403592664E+01    
Expon(5) = 0.4803644023E+01
Coeffs(6) = 0.1383699028E+00    
Expon(6) = 0.1018700979E+02
Coeffs(7) = 0.2788102492E+02    
Expon(7) = 0.4581403450E+02
Coeffs(8) = 0.1377995932E+03    
Expon(8) = 0.9715707822E+02
Coeffs(9) = 0.1770071493E+03    
Expon(9) = 0.2060394364E+03
Coeffs(10) = 0.2081874665E+03    
Expon(10) = 0.4369444835E+03
Coeffs(11) = 0.1309155764E+03    
Expon(11) = 0.9266210634E+03
Coeffs(12) = 0.1507477824E+03    
Expon(12) = 0.1965070226E+04
Coeffs(13) = 0.3952103802E+02    
Expon(13) = 0.4167292486E+04
Coeffs(14) = 0.1069116883E+03    
Expon(14) = 0.8837509432E+04
Coeffs(15) = 0.4291887253E+02    
Expon(15) = 0.3974492714E+05
Coeffs(16) = 0.3152951047E+02    
Expon(16) = 0.8428642090E+05

DO i = 17, 20
	Coeffs(i) = 0.d0
	Expon(i) = 1.d0
END DO

ELSE IF (Z_Atomic == 18) THEN


Coeffs(1) = 0.5473904910E-02    
Expon(1) = 0.1442207120E+00
Coeffs(2) = 0.6829939177E+00    
Expon(2) = 0.6119681072E+00
Coeffs(3) = 0.1177066792E+00    
Expon(3) = 0.1260605990E+01
Coeffs(4) = 0.1164532302E+00    
Expon(4) = 0.5349097580E+01
Coeffs(5) = 0.4416586452E+02    
Expon(5) = 0.1101871874E+02
Coeffs(6) = 0.1067763036E+02    
Expon(6) = 0.2269769075E+02
Coeffs(7) = 0.4037436211E+02    
Expon(7) = 0.9631253829E+02
Coeffs(8) = 0.4393893376E+03    
Expon(8) = 0.1983962257E+03
Coeffs(9) = 0.5563267702E+03    
Expon(9) = 0.4086805628E+03
Coeffs(10) = 0.7430814518E+03    
Expon(10) = 0.8418496968E+03
Coeffs(11) = 0.4609451848E+03    
Expon(11) = 0.1734143917E+04
Coeffs(12) = 0.5336556378E+03    
Expon(12) = 0.3572199570E+04
Coeffs(13) = 0.2203764729E+03    
Expon(13) = 0.7358449114E+04
Coeffs(14) = 0.2800829655E+03    
Expon(14) = 0.1515782428E+05
Coeffs(15) = 0.1771944629E+03    
Expon(15) = 0.3122392140E+05
Coeffs(16) = 0.2143058729E+03    
Expon(16) = 0.1324916776E+06
Coeffs(17) = 0.7643772860E+01    
Expon(17) = 0.2729223965E+06


DO i = 18, 20
	Coeffs(i) = 0.d0
	Expon(i) = 1.d0
END DO


END IF
! End of the IF condition

! In this DO the product of Coefficient*Exp(-Exponent*r^2)is computed for
! each of the available parameters
DO j = 1,20
	density_aux(j) = Coeffs(j)*Exp(-Expon(j)*(r**2.d0))
END DO

! The value of the density is computed
density = Sum(density_aux)

END FUNCTION Density


!====================================================================================
!====================================================================================

! FUNCTION FuncTest3D(r,theta,phi) 
! This function is just a random one chosen to check if the 
! 3D integrating routine works
! Of course it can be changed to a any desired function.
!

FUNCTION FuncTest3D(r,theta,phi) 

USE nrtype		! It uses the module nrtype where the type of used variables is declared



IMPLICIT NONE

REAL(dp), INTENT(IN) :: r, theta, phi !(inputs) spherical coordinates of the point
									  ! at which the function is being evaluated
REAL(dp) :: FuncTest3D				  ! Value of the function at that point


! This example 
FuncTest3D = Exp(-r)*theta*phi


END FUNCTION FuncTest3D

!========================================================================
!========================================================================


! FUNCTION FuncTest6D(r1,theta1,phi1,r2,theta2,phi2) 
! This function is just a random one chosen to check if the 
! 6D integrating routine works
! Of course it can be changed to a any desired function
!

FUNCTION FuncTest6D(r1,theta1,phi1,r2,theta2,phi2) 

USE nrtype	 ! It uses the module nrtype where the type of used variables is declared

IMPLICIT NONE

REAL(dp), INTENT(IN) :: r1, theta1, phi1	!(inputs) spherical coordinates of the first point
											! at which the function is being evaluated
REAL(dp), INTENT(IN) :: r2, theta2, phi2	!(inputs) spherical coordinates of the second point
											! at which the function is being evaluated
REAL(dp) :: FuncTest6D						! Value of the function at those points


! This example
FuncTest6D = r1*theta1*phi1*Exp(-r2)*theta2*phi2

END FUNCTION FuncTest6D



END MODULE DensityFunction