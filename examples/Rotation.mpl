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

z+1

]:     	
# We assume that OutputsVariables[i] = OutputSystem[i].
OutputSystem := [

x^2+y^2-1,
z

] :
#_____________________________________________________________________________#
OutputsVariables:= [] 					:
Inputs 		:= [] 					:
Parameters 	:= [x,y]					:
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

