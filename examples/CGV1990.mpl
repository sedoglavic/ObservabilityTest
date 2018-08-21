#-----------------------------------------------------------------------------#
#								     	System 
# Description 	: 	A model describing tumor targeting by amtibodies
# Result	:
#-----------------------------------------------------------------------------#
infolevel[observabilityTest] := 1 :
infolevel[observabilitySymmetries] := 1 :
#-----------------------------------------------------------------------------#
# Bibliography : 
#@Article{ChappellGodfreyVajda:1990,
#  author =	 {Michael J. Chappell and Keith R. Godfrey and Sandor
#                  Vajda},
#  title =	 {Global identifiability of the Parameters of
#                  Nonlinear Systems with Specified Inputs: A
#                  Comparison of Methods},
#  journal =	 {Mathematical Biosciences},
#  year =	 1990,
#  volume =	 102,
#  pages =	 {41--73}
#}
#
#-----------------------------------------------------------------------------#
# We assume that diff(Variables[i],t) = VectorsField[i]
VectorField:= [

    k4*q3-(k3+k7)*q1 + u,
    k3*q1-k4*q3 - k5*q3*(R*V3-q35)+k6*q35-k5*q3*(5*V36/V3)*(S*V36-q36)+k6*q36,
    k5*q3*(R*V3-q35)-k6*q35,
    k5*q3*(5*V36/V3)*(S*V36-q36)-k6*q36,
    k7*q1

]:     	
# We assume that OutputsVariables[i] = OutputSystem[i].
OutputSystem := [

q7

] :
#-----------------------------------------------------------------------------#
OutputsVariables:= [y] 					:
Inputs := [u] 					:
Parameters 	:= [k3,k4,k5,k6,k7,R,V3,V36,S]		:
# The variables have to be ordered as the vectors field.
Variables 	:= [q1,q3,q35,q36,q7]			:
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









