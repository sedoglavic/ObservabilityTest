#-----------------------------------------------------------------------------#
#								     	System 
# Description : 
# 		
# Result      :
#-----------------------------------------------------------------------------#
infolevel[observabilityTest] := 1 :
infolevel[observabilitySymmetries] := 1 :
#-----------------------------------------------------------------------------#
# Bibliography : 
#@InProceedings{LecourtierLamnabhi-LagarrigueWalter:1987,
#  author =	 {Yves Lecourtier and Fran\c{c}oise
#                  Lamnabhi-Lagarrigue and {\'E}ric Walter},
#  title =	 {A method to prove that nonlinear models can be
#                  unidentifiable},
#  booktitle =	 {Proceedings of the 26th Conference on Decision and
#                  Control},
#  pages =	 {2144-2145},
#  year =	 1987,
#  address =	 {Los Angeles},
#  month =	 dec
#}
#-----------------------------------------------------------------------------#
# We assume that diff(Variables[i],t) = VectorsField[i]
VectorField:= [

    - p1*x1 + p2*u,
    - p3*x2 + p4*u,
    - (p1+p3)*x3 + (p4*x1+p2*x2)*u

]:     	
#-----------------------------------------------------------------------------#
# We assume that OutputsVariables[i] = OutputSystem[i].
OutputSystem := [

x3

] :
#-----------------------------------------------------------------------------#
Inputs 		:= [u] 							:
Parameters 	:= [p1,p2,p3,p4]	:
#-----------------------------------------------------------------------------#
# The variables have to be ordered as the vector field.
Variables       := [x1,x2,x3] 					:
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

