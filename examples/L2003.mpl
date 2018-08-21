#-----------------------------------------------------------------------------#
#								     	System 
# Description 	: 
# Result	:
#-----------------------------------------------------------------------------#
infolevel[observabilityTest] := 1 :
infolevel[observabilitySymmetries] := 1 :
#-----------------------------------------------------------------------------#
# Bibliography : 
#-----------------------------------------------------------------------------#
# We assume that diff(Variables[i],t) = VectorsField[i]
VectorField:= [

-theta1*(x1/theta2-x2/theta3)/(1+x1/theta2+x2/theta3)+ theta4*(u-x1)/(1
+u/theta5+x1/theta5+x1*u/theta5^2),
theta1*(x1/theta2-x2/theta3)/(1+x1/theta2+x2/theta3)-theta6*(x2/theta7-theta8)/
(1+x2/theta7+theta8)

]:     	
# We assume that OutputsVariables[i] = OutputSystem[i].
OutputSystem := [
	x2

] :
#-----------------------------------------------------------------------------#
OutputsVariables:= [y] 					:
Inputs 		:= [u] 					:
Parameters 	:= [theta1,theta2,theta3,theta4,theta5,theta6,theta7,theta8]					:
# The variables have to be ordered as the vectors field.
Variables 	:= [x1,x2] 					:
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

