#-----------------------------------------------------------------------------#
#								     	System 
# Description 	: un exemple simple de système plat
# Result	:
#-----------------------------------------------------------------------------#
infolevel[observabilityTest] := 1 :
infolevel[observabilitySymmetries] := 1 :
#-----------------------------------------------------------------------------#
# Bibliography : 
# @Unpublished{MartinPetit:1999,
#  author =	 {Philippe Martin and Nicolas Petit},
#  title =	 {Commande des syst{\`e}mes non-lin{\'e}aires, le
#                  point de vue des syst{\`e}mes plats},
#  note =	 {{\'E}coles Centrale de Paris, notes de cours. Module
#                  th{\'e}matique M'2 num{\'e}ro 21, 1998}
# }
#-----------------------------------------------------------------------------#
# We assume that diff(Variables[i],t) = VectorsField[i]
VectorField:= [
		    x3 - x2*u,
		    u - x2,
		    x2 - x1 + 2*x2*(u-x2)
]:     	
# We assume that OutputsVariables[i] = OutputSystem[i].
OutputSystem := [
		    x1+x2^2/2
] :
#-----------------------------------------------------------------------------#
OutputsVariables:= [y] 					:
Inputs 		:= [u] 					:
Parameters 	:= []					:
# The variables have to be ordered as the vectors field.
Variables 	:= [x1,x2,x3]				:
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
