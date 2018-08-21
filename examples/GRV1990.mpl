#-----------------------------------------------------------------------------#
#								     	System 
# Description 	: 
# Result	:
#-----------------------------------------------------------------------------#
infolevel[observabilityTest] := 1 :
infolevel[observabilitySymmetries] := 1 :
#-----------------------------------------------------------------------------#
# Bibliography : 
#@Article{VajdaGodfreyRabitz:1989,
#  author =	 {Sandor Vajda and Keith R. Godfrey and H. Rabitz},
#  title =	 {Similarity transformation approach to
#                  identifiability analysis of non linear
#                  comportemental models },
#  journal =	 {Mathematical Biosciences},
#  year =	 1989,
#  number =	 93,
#  pages =	 {217-248}
#}
#-----------------------------------------------------------------------------#
#\cite{VajdaGodfreyRabitz:1989}:
#  \begin{equation}
#    \label{eq:GRV1990}
#    \left\lbrace
#      \begin{array}{lll}
#        \dot{x}_{1} & = & \theta_{1}x_{1}^{2} + \theta_{2}x_{1}x_{2} + u,\\
#        \dot{x}_{2} & = & \theta_{3}x_{1}^{2} + \theta_{4}x_{1}x_{2},\\
#        y & = & x_{1}.
#      \end{array}
#    \right.  
#\end{equation}
#-----------------------------------------------------------------------------#
# We assume that diff(Variables[i],t) = VectorsField[i]
VectorField:= [

    p1*x1**2 + p2*x1*x2 + u,
    p3*x1**2 + p4*x1*x2

]:     	
# We assume that OutputsVariables[i] = OutputSystem[i].
OutputSystem := [

    x1

] :
#-----------------------------------------------------------------------------#
Inputs := [u] 					:
Parameters 	:= [p1,p2,p3,p4]			:
# The variables have to be ordered as the vectors field.
Variables 	:= [x1,x2] 				:
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




