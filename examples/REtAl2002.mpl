#-----------------------------------------------------------------------------#
#								     	System 
# Description 	: 
#			administration par voie intra-musculaire
# Result	:
#-----------------------------------------------------------------------------#
infolevel[observabilityTest] := 1 :
infolevel[observabilitySymmetries] := 2 :
#-----------------------------------------------------------------------------#
# Bibliography : 
# @Article{RougierEtAl2002,
#  author =	 {Rougier, Florent and Claude, Daniel and Maurin,
#                  Michel and Sedoglavic, Alexandre and Ducher, Michel
#                  and Corvaisier, St{\'e}phane and Jelliffe, Roger and
#                  Maire, Pascal},
#  title =	 {Aminoglycoside Nephrotoxicity: Modeling, Simulation
#                  and Control.},
#  journal =	 { Antimicrobial Agents and Chemotherapy},
#  year =	 2003,
#  volume =       47,
#  number =       3,
#  pages =	 {1010--1016},
#  month= mar
#}
#-----------------------------------------------------------------------------#
# We assume that diff(Variables[i],t) = VectorsField[i]

# EQ.8
E := Emax*Qcgamma/(Q50gamma+Qcgamma) ;

# EQ.6
DIFFQc := (-k1)*Qc+V*Qs/(km+Qs) : 

# EQ.10 

CCr := -(CCr0*(E50delta+Edelta)-Edelta*CCrmax)/(-1+eta*(SC*cosphi+CS*sinphi))/(E50delta+Edelta):

# Relations auxilliaires
DIFFE :=  E*gaMMa*DIFFQc*Q50gamma/(Q50gamma+Qcgamma)/Qc :

# le premier ensemble correspond au modele ; 
# le second correspond aux variables auxilliaire permettant d'exprimer
# ce modele sous forme algebrique.

VectorField:= [

-kabs*Qm+u,
-((Ki+(Ks-kreabs)*CCr)+kcp)*Qs + kpc*Qp + kabs*Qm,
-kpc*Qp+kcp*Qs,
(Ki+(Ks-kreabs)*CCr)*Qs + k1*Qc,
DIFFQc,
-CCr/Vol*SCr + k2,

-alpha*V*DIFFE,
Qcgamma*gaMMa*DIFFQc/Qc,
-omega*SC,
omega*CS,
Edelta*delta*DIFFE/E
]:     	
# We assume that OutputsVariables[i] = OutputSystEmax[i].
OutputSystem := [

SCr, u - (Qu+Qs+Qp+Qc)

] :
#-----------------------------------------------------------------------------#
OutputsVariables:= [] 							:
Inputs 		:= [u] 							:
Parameters 	:= [Ki,Ks,kcp,kpc,kreabs,km,k1,Emax,Q50gamma,gaMMa,CCrmax,E50delta,delta,alpha,k2,Vol,CCr0,eta,kabs,omega,cosphi,sinphi]	:
# The variables have to be ordered as the vectors field.
Variables 	:= [Qm,Qs,Qp,Qu,Qc,SCr,V,Qcgamma,CS,SC,Edelta]		:
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
