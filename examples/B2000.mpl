#-----------------------------------------------------------------------------#
#								     	System 
# Description 	: 
# Result	: 
#-----------------------------------------------------------------------------#
# Bibliography : 
#-----------------------------------------------------------------------------#
infolevel[observabilityTest] := 1 :
infolevel[observabilitySymmetries] := 1 :

# We assume that diff(Variables[i],t) = VectorsField[i]
hx4x3theta4 	:= x4/(1+x3/theta2) 	:
dx2 		:= hx4x3theta4 - 2*x2 	:

VectorField:= [

u - theta1*x1,
dx2,
x1-x3-4*x3*E3,
4*x3*E2-x4,
-theta3*dx2*E3,
-theta2*dx2*E2

]:     	
# We assume that OutputsVariables[i] = OutputSystem[i].
OutputSystem := [

    x4,
    hx4x3theta4

] :
#-----------------------------------------------------------------------------#
OutputsVariables:= [y1,y2]			:
Inputs := [ u ] 					:
Parameters 	:= [theta1,theta2,theta3]	:
# The variables have to be ordered as the vectors field.
Variables 	:= [x1,x2,x3,x4,E3,E2]			:
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




