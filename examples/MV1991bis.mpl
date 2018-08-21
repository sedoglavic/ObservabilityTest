#-----------------------------------------------------------------------------#
#								     	System 
# Description 	: 	model for an induction motor.
#			wa angular speed of the frame of references.
#			m0 load torque    
# Result	:
#-----------------------------------------------------------------------------#
infolevel[observabilityTest] := 1 :
infolevel[observabilitySymmetries] := 1 :
#-----------------------------------------------------------------------------#
# Bibliography : 
#@InProceedings{MarinoValigi:1991,
#  author = 	 {Riccardo Marino and Paolo Valigi},
#  title = 	 {Nonlinear Control of Induction Motors: a Simulation Study},
#  booktitle =	 {Proceedings of First European Control Conference},
#  editor =	 {C. Commault and coll.},
#  pages =	 {1057--1062},
#  year =	 1991,
#  volume =	 1,
#  month =	 jul # {~2--5},
#  address =	 {Grenoble, France},
#  publisher =	 {Herm{\`e}s}
#}
#-----------------------------------------------------------------------------#
# We assume that diff(Variables[i],t) = VectorsField[i]

sigma := Ls - M**2/Lr :
GammaN := (M**2*Rr + Lr**2*Rs)/(sigma*Lr**2) :

VectorField:= [
    
    np*M/(J*Lr)*(Psix*Iy - Psiy*Ix) - Tl/J - frot*omega,

    Rr/Lr * (M*Ix - Psix) - np*omega*Psiy,
    
    np*omega*Psix + Rr/Lr * (M*Iy - Psiy),

   M*Rr/(sigma*Lr**2)*Psix + np*M/(sigma*Lr)*omega*Psiy - GammaN*Ix + Ux/sigma,

   M*Rr/(sigma*Lr**2)*Psiy - np*M/(sigma*Lr)*omega*Psix - GammaN*Iy + Uy/sigma
]:     	

# omega est la vitesse angulaire rototique

# Ix et Iy les courants statoriques

# Psix et Psiy les flux rotoriques

# les  autres parametres  correspondent  aux resistances,  inductance,
# nombre de po^le de la machine


# VectorField[1] est une equation mecanique
# np*M/(J*Lr)*(Psix*Iy - Psiy*Ix) est un couple electromagnetique
# frotk est un coefficient de frottement
# Tl/J represente le couple de charge
# 
# les quatre autres equations sont electriques

# We assume that OutputsVariables[i] = OutputSystem[i].
OutputSystem := [

    M*(Psix*Iy - Psiy*Ix)/Lr,
    omega/np

] :

#les sorties representes le couple electromagnetique
# et la vitesse mecanique

#-----------------------------------------------------------------------------#
OutputsVariables:= [y1,y2] 					:
Inputs 		:= [Ux,Uy]				:
Parameters 	:= [Rs,Rr,Lr,Ls,M,J,Tl,np,frot]		:
# The variables have to be ordered as the vectors field.
Variables 	:= [omega,Psix,Psiy,Ix,Iy]		:
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
