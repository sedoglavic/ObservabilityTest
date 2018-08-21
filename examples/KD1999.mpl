#-----------------------------------------------------------------------------#
#								     	System 
# Description : A model for a chemical reactor
#
# The exponantial is related to the temperature dependency of the reaction 
# rate (Law of Arrhenius).
#
# Result	: Unidentifiable
#-----------------------------------------------------------------------------#
infolevel[observabilityTest] := 1 :
infolevel[observabilitySymmetries] := 1 :
#-----------------------------------------------------------------------------#
# Bibliography : page 108 in 
#@Book{KumarDaoutidis:1999,
#  author =	 {Aditya Kumar and Prodomos Daoutidis},
#  title =	 {Control of nonlinear differential algebraic equation
#                  systems},
#  publisher =	 {Chapman and Hall / CRC},
#  year =	 1999,
#  number =	 397,
#  series =	 {Research Notes in Mathematics}
#}
#-----------------------------------------------------------------------------#
# We assume that diff(Variables[i],t) = VectorsField[i]

# Arr := k0*exp(-E/(R*T))

dTemp := Fa*(Ta-T)/V - (k0*Arr*Ca*DH + UA*(Tj-T)/V)/(ro*cp):

VectorField:= [

	Fa*(Ca0-Ca)/V - k0*Arr*Ca,
	-Fa*Cb/V + k0*Arr*Ca,
	dTemp,
	Fh*(Th-Tj)/Vh - UA/(roh*cph)*(Tj-T)/Vh,
	E*Arr/(R*T^2)*dTemp

]:     	
#-----------------------------------------------------------------------------#
# We assume that OutputsVariables[i] = OutputSystem[i].
OutputSystem := [

    Cb,
    T

] :
#-----------------------------------------------------------------------------#
Inputs 		:= [Fa,Fh]	:
Parameters 	:= [Ca0,V,Ta,ro,cp,DH,UA,Th,Vh,roh,cph,k0,E,R]	:
#-----------------------------------------------------------------------------#
# The variables have to be ordered as the vector field.
Variables       := [Ca,Cb,T,Tj,Arr]	:
#-----------------------------------------------------------------------------#
# CAUTION read section II of the file INSTALL
libname := cat(currentdir(),"/release"),libname :
#-----------------------------------------------------------------------------#
readlib(observabilityTest) :
NonObservable := observabilityTest(	VectorField	,
					Variables	,
					OutputSystem	,
					Parameters	,
					Inputs			) :
print(%) :					
GroupInfGen := observabilitySymmetries(	VectorField	,
					Variables	,
					OutputSystem	,
					Parameters	,
					Inputs		,
					NonObservable		) :
print(%) :
#-----------------------------------------------------------------------------#
quit :




