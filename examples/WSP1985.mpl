#-----------------------------------------------------------------------------#
#                                                                       System 
# Description   :       Two compartmental model, with a Michaelis-Menten 
#                       elimination pathway from compartment~1 and a linear 
#                       elimination from compartment~2.
# Result        :
#-----------------------------------------------------------------------------#
infolevel[observabilityTest] := 1 :
infolevel[observabilitySymmetries] := 1 :
#-----------------------------------------------------------------------------#
# Bibliography : 
#@Article{WagnerSzunparPerry:1985,
#  author =      {Wagner, J. G.  and Szunpar, G. J.  and Perry, J. J. },
#  title =       {A nonlinear physiologic pharmacokinetic model: 1
#                  steady-state},
#  journal =     {Journal of Pharmacokinetics and Biopharmaceutics},
#  year =        1985,
#  volume =      13,
#  pages =       {73-92}
#}
#-----------------------------------------------------------------------------#
# We assume that diff(Variables[i],t) = VectorsField[i]
VectorField:= [

    - x1*(p1 + p2/(x1+p3)) + p4*x2 + p5*u,
    p1*x1 - (p4+p6)*x2
    
]:      
# We assume that OutputsVariables[i] = OutputSystem[i].
OutputSystem := [

    p7*x1

] :
#-----------------------------------------------------------------------------#
OutputsVariables := [y]                                  :
Inputs		 := [u]                                  :
Parameters       := [p1,p2,p3,p4,p5,p6,p7]               :
# The variables have to be ordered as the vectors field.
Variables       := [x1,x2]                              :
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
