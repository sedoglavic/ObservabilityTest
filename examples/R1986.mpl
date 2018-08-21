#-----------------------------------------------------------------------------#
#								     	System 
# Description 	:  a pharmacokinetic model
# Result	:
#-----------------------------------------------------------------------------#
infolevel[observabilityTest] := 1 :
infolevel[observabilitySymmetries] := 2 :
#-----------------------------------------------------------------------------#
# Bibliography : 
# @PhdThesis{Raksanyi:1986,
#   author =	 {Atilla Raksanyi},
#   title =	 {Utilisation du calcul formel pour l'{\'e}tude des
#                   syst{\`e}mes d'{\'e}quations polynomiales
#                   (applications en mod{\'e}lisation)},
#   school =	 {Universit{\'e} Paris-Dauphine},
#   year =	 1986
# }
#-----------------------------------------------------------------------------#
# We assume that diff(Variables[i],t) = VectorsField[i]
VectorField:= [

    u - (c1+c2)*x1,
    c1*x1-(c3+c6+c7)*x2 + c5*x4,
    c2*x1 + c3*x2 - c4*x3,
    c6*x2 - c5*x4

]:     	
# We assume that OutputsVariables[i] = OutputSystem[i].
OutputSystem := [
    
   c8*x3,
    c9*x2
		
] :
#-----------------------------------------------------------------------------#
OutputsVariables:= [y1,y2]				:
Inputs 		:= [u] 					:
Parameters 	:= [c1,c2,c3,c4,c5,c6,c7,c8,c9]		:
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
