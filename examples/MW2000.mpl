#-----------------------------------------------------------------------------#
#								     	System 
# Description 	: Multispecies Transmission Model, 
#		  Mastitis model with control
# Result	:
#-----------------------------------------------------------------------------#
infolevel[observabilityTest] := 1 :
infolevel[observabilitySymmetries] := 1 :
#-----------------------------------------------------------------------------#
# Bibliography : 
#-----------------------------------------------------------------------------#
# We assume that diff(Variables[i],t) = VectorsField[i]

b 	:= mu+c1*(y1+y12)      :
Lambda1 := Beta1*(y1+y12)      :
Lambda2 := Beta2*(y2+y12) + I2 :

VectorField:= [

    (1-Theta1-Theta2)*b - (m1*Lambda1 + m2*Lambda2 + mu)*x12 + 
    (nu1+tau)*y1 + (nu2+tau)*y2 + tau*y12,
    Theta1*b + m1*Lambda1*x12 + nu2*y12 -((1-Pi2)*m2*Lambda2+nu1+mu+c1+tau)*y1,
    Theta2*b + m2*Lambda2*x12 + nu1*y12 -((1-Pi1)*m1*Lambda1+nu2+mu+tau)*y2,
    (1-Pi1)*m1*Lambda1*y2 + (1-Pi2)*m2*Lambda2*y1 - (nu1+nu2+mu+c1+tau)*y12
    
]:     	
# We assume that OutputsVariables[i] = OutputSystem[i].
OutputSystem := [
		    x12+y1+y2+y12,
		    y1+y12,
		    y2+y12
] :
#OutputSystem[1] is an invariant equation.
#-----------------------------------------------------------------------------#
OutputsVariables:= [YO1,YO2,YO3] 					:
Inputs := [] 					:
Parameters 	:= [Theta1,Theta2,mu,nu1,nu2,Pi1,Pi2,Beta1,Beta2,m1,m2,I2,tau,c1] :
# The variables have to be ordered as the vectors field.
Variables 	:= [x12,y1,y2,y12] 			:
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
