#-----------------------------------------------------------------------------#
#								     	System 
# Description 	: 
# Result	:
#-----------------------------------------------------------------------------#
infolevel[observabilityTest] := 1 :
infolevel[observabilitySymmetries] := 1 :
#-----------------------------------------------------------------------------#
# Bibliography : Martin Fransson (private communication)
#-----------------------------------------------------------------------------#
# We assume that diff(Variables[i],t) = VectorsField[i]
VectorField:= [
-(k10+k12+k13)*A1+k21*A2+k31*A3+u1, k12*A1-k21*A2, k13*A1-k31*A3

]:     	
# We assume that OutputsVariables[i] = OutputSystem[i].
OutputSystem := [

(1+B1*u2+B2+B3/(K+A1/V1))*A1/V1

] :
#_____________________________________________________________________________#
OutputsVariables:= [Y] 					:
Inputs 		:= [u1,u2] 					:
Parameters 	:= [k12, k21, k13, k31, k10, V1, B1, B2, B3, K] :
# The variables have to be ordered as the vectors field.
Variables 	:= [A1, A2, A3] :
#_____________________________________________________________________________#
# CAUTION read section II of the file INSTALL
libname := cat(currentdir(),"/release"),libname ;
#_____________________________________________________________________________#
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
#_____________________________________________________________________________#
quit :
