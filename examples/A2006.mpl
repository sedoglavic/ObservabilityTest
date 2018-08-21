#-----------------------------------------------------------------------------#
#								     	System 
# Description 	: 
# Result	:
#
#-----------------------------------------------------------------------------#
infolevel[observabilityTest] := 1 :
infolevel[observabilitySymmetries] := 1 :
#-----------------------------------------------------------------------------#
# Bibliography : Milena Anguelova, private communication 2006
# This examples is observable while previous code say opposite. The bug is know
# fixed. 
#-----------------------------------------------------------------------------#
# We assume that diff(Variables[i],t) = VectorsField[i]
VectorField:= [

	x21 - x11,
	(x10-x11)/r,
	- x21,
	(x20-x21)/r
]:     	
# We assume that OutputsVariables[i] = OutputSystem[i].
OutputSystem := [

	x10
] :
#-----------------------------------------------------------------------------#
OutputsVariables:= [y] 					:
Inputs 		:= [] 					:
Parameters 	:= [r]					:
# The variables have to be ordered as the vectors field.
Variables 	:= [x10,x11,x20,x21] 					:
#-----------------------------------------------------------------------------#
# CAUTION read section II of the file INSTALL
libname := cat(currentdir(),"/release"),libname :
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
