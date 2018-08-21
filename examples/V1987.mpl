#-----------------------------------------------------------------------------#
#								     	System 
# Description 	: model of a flow reactor to pyrolyze methane at 
#		T = 1409 K (constant).
# Result	: identifiable
#-----------------------------------------------------------------------------#
infolevel[observabilityTest] := 1 :
infolevel[observabilitySymmetries] := 1 :
#-----------------------------------------------------------------------------#
# Bibliography : 
#@InProceedings{Vajda:1987,
#  author =	 {Sandor Vajda},
#  title =	 {Identifiability of polynomial systems: structural
#                  and numerical aspect},
#  booktitle =	 {Identifiability of parametric models},
#  pages =	 {42-48},
#  year =	 1987,
#  editor =	 {{\'E}ric Walter},
#  publisher =	 {Pergamon Press}
#}
#-----------------------------------------------------------------------------#
# We assume that diff(Variables[i],t) = VectorsField[i]
VectorField:= [

	-x1*(k1+k2*x4) + k5*x3*x4,
	k2*x1*x4 - (k3+k4)*x2,
	k4*x2 - k5*x3*x4,
	x1*(k1+k2*x4) +2*k3*x2 - k5*x3*x4

]:     	
# We assume that OutputsVariables[i] = OutputSystem[i].
OutputSystem := [

    x1,
    x2

] :
#-----------------------------------------------------------------------------#
OutputsVariables:= [y1,y2]				:
Inputs := [] 					:
Parameters 	:= [k1,k2,k3,k4,k5]			:
# The variables have to be ordered as the vectors field.
Variables 	:= [x1,x2,x3,x4]			:
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







