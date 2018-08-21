
#-----------------------------------------------------------------------------#
#								     	System 
# Description : A model for circadian oscillations in the Drosophila
# 		period protein
# Result	: Unidentifiable
#-----------------------------------------------------------------------------#
infolevel[observabilityTest] := 1 :
infolevel[observabilitySymmetries] := 1 :
#-----------------------------------------------------------------------------#
# Bibliography : 
#@Article{Goldbeter:1995,
#  author =	 {A. Goldbeter},
#  title =	 {A model for circadian oscillations in the
#                  {D}rosophila period protein},
#  journal =	 {Proceedings of the Royal Society London},
#  year =	 1995,
#  volume =	 {B},
#  number =	 261,
#  pages =	 {319-324}
#}
#-----------------------------------------------------------------------------#
# We assume that diff(Variables[i],t) = VectorsField[i]
VectorField:= [

    vs*KI^4/(KI^4 + PN^4) - (vm*M/(Km + M)) 				,
    (ks*M - V1*P0/(K1+P0) + V2*P1/(K2+P1)) 				,
    (V1*P0/(K1+P0) - P1*(V2/(K2+P1) + V3/(K3+P1)) + V4*P2/(K4+P2))	,
    (V3*P1/(K3+P1) - P2*(V4/(K4+P2) + k1 + vd/(Kd+P2)) + k2*PN )	,
    k1*P2 - k2*PN 

]:     	
#-----------------------------------------------------------------------------#
# We assume that OutputsVariables[i] = OutputSystem[i].
OutputSystem := [

    PN+P0+P1+P2

] :
#-----------------------------------------------------------------------------#
Inputs 		:= [] 							:
Parameters 	:= [KI,K1,K2,K3,K4,Kd,Km,ks,vd,vs,vm,k2,k1,V1,V2,V3,V4]	:
#-----------------------------------------------------------------------------#
# The variables have to be ordered as the vector field.
Variables       := [M,P0,P1,P2,PN] 					:
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




