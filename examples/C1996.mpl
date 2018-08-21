#-----------------------------------------------------------------------------#
#                                                                       System 
# Description   : 
# Result        :  Unidentifiable 
#-----------------------------------------------------------------------------#
infolevel[observabilityTest] := 1 :
infolevel[observabilitySymmetries] := 1 :
#-----------------------------------------------------------------------------#
# Bibliography : 
#@Article{Chappell:1996,
#  author =      {Michael J. Chappell},
#  title =       {Structural identifiability of models characterizing
#                  saturable binding: Comparison of pseudo-steady-state
#                  and non-pseudo-steady-state model formulations},
#  journal =     {Mathematical Biosciences},
#  year =        1996,
#  number =      133,
#  pages =       {1--20}
#}
#-----------------------------------------------------------------------------#
# We assume that diff(Variables[i],t) = VectorsField[i]
VectorField:= [

    -e*x1 - a*x1*(p-x2) + d*x2,
    a*x1*(p-x2) - d*x2
    
]:      
# We assume that OutputsVariables[i] = OutputSystem[i].
OutputSystem := [

c*(x1+x2)

] :
#-----------------------------------------------------------------------------#
OutputsVariables:= [y]                                  :
Inputs	        := []                                   :
Parameters      := [a,c,d,e,p]                          :
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
