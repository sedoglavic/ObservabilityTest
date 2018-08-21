#-----------------------------------------------------------------------------#
#								     	System 
# Description 	: 
# Result	:
#-----------------------------------------------------------------------------#
infolevel[observabilityTest] := 1 :
infolevel[observabilitySymmetries] := 1 :
#-----------------------------------------------------------------------------#
# Bibliography : 
#	@Book{Walter1982,
#  		author =	 {Walter, {\'E}ric },
#  		title =	 {Identifiability of state space model},
#  		publisher =	 {Springer},
#  		year =	 1982,
#  		volume =	 46,
#  		series =	 {Lectures notes in biomathematics},
#  		address =	 {New York}
#	}
#-----------------------------------------------------------------------------#
# We assume that diff(Variables[i],t) = VectorsField[i]
infolevel[observabilityTest] := 1 :
infolevel[observabilitySymmetries] := 1 :

VectorField:= [

	-(k01+k21)*x1 + k12*x2 +k13*x1*x3,
	k21*x1 -(k02+k12+k32+k42)*x2-k23*x3*x2+k24*x4,
	k32*x2-(k03+k13*x1+k23*x2-k43*x4)*x3+k34*x4+u,
	k42*x2+k43*x3*x4-(k24+k34)*x4

]:     	

# We assume that OutputsVariables[i] = OutputSystem[i].
OutputSystem := [

	x1

] :
#-----------------------------------------------------------------------------#
OutputsVariables:= [] 					:
Inputs 		:= [u] 					:
Parameters 	:= [k01,k02,k03,k12,k13,k21,k23,k24,k32,k34,k42,k43]:
# The variables have to be ordered as the vectors field.
Variables 	:= [x1,x2,x3,x4] 					:
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



