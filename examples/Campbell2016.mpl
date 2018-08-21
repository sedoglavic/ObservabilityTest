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


 -p_3* x_1 - p_4* x_1^2 + p_1 *x_2 + p_2 *x_3 - p_7 *x_1 *u,
p_5* p_6* x_1 - p_3* x_2 - p_4* x_1* x_2 - p_5* x_2 - p_7* x_2* u,
2* p_5* p_6* x_2 - p_3* x_3 - p_4* x_1* x_3 - 2* p_5* x_3 - p_7* x_3 *u

]:     	
# We assume that OutputsVariables[i] = OutputSystem[i].
OutputSystem := [

p_7* x_1* u,
p_7* x_2* u,
p_7* x_3* u



] :

#-----------------------------------------------------------------------------#
OutputsVariables:= [y1,y2,y3] 					:
Inputs 		:= [u] 					:
Parameters 	:= [p_1,p_2,p_3,p_4,p_5,p_6,p_7]					:
# The variables have to be ordered as the vectors field.
Variables 	:= [x_1,x_2,x_3] 					:
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
