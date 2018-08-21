#-----------------------------------------------------------------------------#
#								     	System 
# Description 	: 
# Result	:
#-----------------------------------------------------------------------------#
infolevel[observabilityTest] := 1 :
infolevel[observabilitySymmetries] := 1 :
#-----------------------------------------------------------------------------#
# Bibliography : 
#	\begin{eqnarray}
# 		\dot{z} & = & v \\
# 		\dot{v} & = & \frac{k}{m}\left(\frac{i}{c - z} \right)^2 - g
# 	\end{eqnarray}
#	z: position, v: vitesse, i: courant.
#	les paramètres sont k, m, c et g.                   
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# We assume that diff(Variables[i],t) = VectorsField[i]
VectorField:= [

v,
k*(i/(c-z))^2/m - g

]:     	
# We assume that OutputsVariables[i] = OutputSystem[i].
OutputSystem := [

z

] :
#-----------------------------------------------------------------------------#
OutputsVariables:= [y] 					:
Inputs 		:= [i] 					:
Parameters 	:= [k,c,m,g]				:
# The variables have to be ordered as the vectors field.
Variables 	:= [z,v] 					:
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


 
#with(diffalg) ;

#R := differential_ring (derivations = [t], 
#    ranking = [[v,g,k],c,[z,i]]);

#V := [
#k[t],
#c[t],
#g[t],
#z[t]-v[],
#-k[]*i[]^2 + (v[t]+g[])*(c[]-z[])^2
#
#]:     	

#E := Rosenfeld_Groebner(V,R) ;
#map(equations,E) ;

#SOL := map(rewrite_rules,E) ;

#save SOL,"/scratch/sedoglavic/Levitation.m";
quit :
