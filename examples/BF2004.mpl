#-----------------------------------------------------------------------------#
#								     	System 
# Description 	: 
# Result	:
#-----------------------------------------------------------------------------#
infolevel[observabilityTest] := 1 :
infolevel[observabilitySymmetries] := 1 :
#-----------------------------------------------------------------------------#
# Bibliography : 
#  Boubaker, O.    Fourati, A.  
# Inst. Nat. des Sci. Appliquees et de Technol., Tunis, Tunisia;
# This paper appears in: Industrial Technology, 2004. IEEE ICIT
# '04. 2004 IEEE International Conference on
# Publication Date: 8-10 Dec. 2004
# Volume: 3,  On page(s): 1244- 1248 Vol. 3
#-----------------------------------------------------------------------------#
# We assume that diff(Variables[i],t) = VectorsField[i]
VectorField:= [

	mumax*S*X/(S+Ks)-d*X,
	-mumax*S*X/(S+Ks)/Y + d*(Sin - S) 

]:     	
# We assume that OutputsVariables[i] = OutputSystem[i].
OutputSystem := [

S

] :
#-----------------------------------------------------------------------------#
OutputsVariables:= [y] 					:
Inputs 		:= [] 					:
Parameters 	:= [mumax,Ks,d,Y,Sin]					:
# The variables have to be ordered as the vectors field.
Variables 	:= [X,S] 					:
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
quit :

