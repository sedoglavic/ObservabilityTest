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

z

]:     	
# We assume that OutputsVariables[i] = OutputSystem[i].
OutputSystem := [

z,
(y-b)^2+a^2-L^2/4,
(x-a)^2+b^2-L^2/4,
x^2+y^2-L^2

] :
#_____________________________________________________________________________#
OutputsVariables:= [z] 					:
Inputs 		:= [] 					:
Parameters 	:= [a,b,x,y,L]					:
# The variables have to be ordered as the vectors field.
Variables 	:= [z] 					:
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
#quit :

