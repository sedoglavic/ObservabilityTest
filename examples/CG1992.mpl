#-----------------------------------------------------------------------------#
#                                                                       System 
# Description   : 
# Result        : Unidentifiable
#-----------------------------------------------------------------------------#
infolevel[observabilityTest] := 1 :
infolevel[observabilitySymmetries] := 1 :
#-----------------------------------------------------------------------------#
# Bibliography : 
#@Article{ChappellGodfrey:1992,
#  author =      {Michael J. Chappell and Keith R. Godfrey},
#  title =       {Structural identifiability of the parameters of a
#                  nonlinear batch reactor model},
#  journal =     {Mathematical Biosciences},
#  year =        1992,
#  number =      108,
#  pages =       {241-251}
#}
#-----------------------------------------------------------------------------#
# We assume that diff(Variables[i],t) = VectorsField[i]
VectorField:= [

    x1*(p1*x2/(p2+x2) - p3) + p4*u,
    - p1*x1*x2/(p5*(p2+x2)) + p6*u

]:      
# We assume that OutputsVariables[i] = OutputSystem[i].
OutputSystem := [

    x1

] :
#-----------------------------------------------------------------------------#
OutputsVariables:= [y]                                  :
Inputs	        := [u]                                  :
Parameters      := [p1,p2,p3,p4,p5,p6]                  :
# The variables have to be ordered as the vectors field.
Variables       := [x1,x2]                              :
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

