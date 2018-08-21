#-----------------------------------------------------------------------------#
#								     	System 
# Description 	: 
# Result	:
#-----------------------------------------------------------------------------#
infolevel[observabilityTest] := 1 :
infolevel[observabilitySymmetries] := 2 :
#-----------------------------------------------------------------------------#
# Bibliography : 
#	@Article{CaumoCobelli1993,
#  		author =	 {Caumo and Cobelli, Claudio},
#  		title =	 {Hepatic glucose production during the labelled
#               	   {IVGTT} estimation by the convolution with a new
#                	  minimal model},
#  		journal =	 {Amer. J. Physiol.},
#  		year =	 1993,
#  		volume =	 264,
#  		pages =	 {829--841}
#	}
#-----------------------------------------------------------------------------#
# We assume that diff(Variables[i],t) = VectorsField[i]
VectorField:= [

-(kp+(F01/V1)/g+k21)*x1 + k12*x2,
k21*x1-(k02+x3+k12)*x2,
-kb*x3+ka*u

]:     	
# We assume that OutputsVariables[i] = OutputSystem[i].
OutputSystem := [

x1/V1

] :
#-----------------------------------------------------------------------------#
OutputsVariables:= [y] 					:
Inputs 		:= [g,u]				:
Parameters 	:= [kp,F01,V1,k21,k12,k02,kb,ka]	:
# The variables have to be ordered as the vectors field.
Variables 	:= [x1,x2,x3] 					:
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
#kernelopts(profile=true);writeto('output');
GroupInfGen := observabilitySymmetries(	VectorField	,
					Variables	,
					OutputSystem	,
					Parameters	,
					Inputs		,
					NonObservable		) :
print(%) :
# kernelopts(profile=false);writeto(terminal);exprofile('output',alpha);

#-----------------------------------------------------------------------------#
quit :
