#-----------------------------------------------------------------------------#
#								     	System 
# Description 	: 
# Result	:
#-----------------------------------------------------------------------------#
infolevel[observabilityTest] := 1 :
infolevel[observabilitySymmetries] := 1 :
#-----------------------------------------------------------------------------#
# Bibliography : Don't remember from where this example is taken :-(
#-----------------------------------------------------------------------------#
# We assume that diff(Variables[i],t) = VectorsField[i]

VectorField:= [

c2*(Theta1*(x6-x1)+u1)/(Theta2*Theta3),
((Theta1-Theta6*Theta7)*(x1-x2)+Theta5*(x3-x2))/((c3*Theta3*(1-Theta2)-Theta8*Theta7))-Theta4*x2,
Theta5*(x2-x3)/(c4*Theta9),
(Theta6*(x1-x4)+Theta10*(x5-x4))/Theta8,
Theta10*(x4-x5)/(Theta11*(1-Theta8)),
(c1*Theta6*Theta7*(x4-x2)+Theta1*(x2-x6))/(Theta2*Theta3)


]:     	
# We assume that OutputsVariables[i] = OutputSystem[i].
OutputSystem := [

55*Theta12*x1/100,
Theta12*(55*Theta8*x4/100+Theta11*(1-Theta8)*x5),
c1,c3,c4
] :
#____________________________________________________________________________#
OutputsVariables:= [y1,y2,y3,y4,y5] 					:
Inputs 		:= [u1] 					:
Parameters 	:= [c1,c2,c3,c4,Theta1,Theta2,Theta3,Theta4,Theta5,Theta6,Theta7,Theta8,Theta9,Theta10,Theta11,Theta12]					:
# The variables have to be ordered as the vectors field.
Variables 	:= [x1,x2,x3,x4,x5,x6] 					:
#____________________________________________________________________________#
# CAUTION read section II of the file INSTALL
libname := cat(currentdir(),"/release"),libname :
#____________________________________________________________________________#
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
