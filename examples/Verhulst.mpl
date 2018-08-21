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
	
	x*(a-b*x)-c*x

]:     	
# We assume that OutputsVariables[i] = OutputSystem[i].
OutputSystem := [

	d*x

] :
#-----------------------------------------------------------------------------#
OutputsVariables:= [y] 					:
Inputs 		:= [] 					:
Parameters 	:= [a,b,c,d]					:
# The variables have to be ordered as the vectors field.
Variables 	:= [x] 					:
#-----------------------------------------------------------------------------#
# CAUTION read section II of the file INSTALL
libname := cat(currentdir(),"/release"),libname :
#-----------------------------------------------------------------------------#
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
#quit :

