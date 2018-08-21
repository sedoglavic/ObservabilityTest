#-----------------------------------------------------------------------------#
#								     	System
# Description 	:
# Result	:
#-----------------------------------------------------------------------------#
infolevel[observabilityTest] := 1 :
infolevel[observabilitySymmetries] := 1 :
#-----------------------------------------------------------------------------#
# Bibliography : Agostino Martinelli
#-----------------------------------------------------------------------------#
# We assume that diff(Variables[i],t) = VectorsField[i]
GAMMAdiff:=xi_q*nu-mu*eta_q*nu*(s*phi_c-c*phi_s):
MUdiff:=-mu^2* eta_q*nu*(c*phi_c+s*phi_s):

VectorField:= [
MUdiff,
GAMMAdiff,
-GAMMAdiff*s,
GAMMAdiff*c,
(GAMMAdiff*c/(mu+c)-(MUdiff-GAMMAdiff)*s/((mu+c)^2))/(1+s^2/(mu+c)^2)


]:
# We assume that OutputsVariables[i] = OutputSystem[i].
OutputSystem := [
BETA-psi,
phi_c^2+phi_s^2,
    c^2+s^2
] :
#_____________________________________________________________________________#
OutputsVariables:= [Y] 					:
Inputs 		:= [nu] 					:
Parameters 	:= [phi_c,phi_s,psi,eta_q,xi_q]					:
# The variables have to be ordered as the vectors field.
Variables 	:= [mu,GAMMA,c,s,BETA] 					:
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

