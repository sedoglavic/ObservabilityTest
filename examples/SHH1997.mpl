#-----------------------------------------------------------------------------#
#								     	System 
# Description : mathematical modelling of a part of the blood
# coagulation mechanism.  The concentration of RVV is supposed to be
# constant (we consider it as a parameter).
# Result	: Non observable, Trans. Degr 1 
#-----------------------------------------------------------------------------#
infolevel[observabilityTest] := 1 :
infolevel[observabilitySymmetries] := 1 :
#-----------------------------------------------------------------------------#
# Bibliography : 
#@TechReport{StortelderHemkerHemker:1997,
#  author =	 {W. J. H. Stortelder and P. W. Hemker and
#                  H. C. Hemker},
#  title =	 {Mathematical modelling in blood coagulation;
#                  simulation and parameter estimation},
#  institution =	 {Centrum voor Wiskunde en Informatica},
#  year =	 1997,
#  number =	 {MAS-R9720},
#  month =	 sep # {~30},
#  address = 	 {available at {\sf http://www.cwi.nl/}}
#}
#-----------------------------------------------------------------------------#
# We assume that diff(Variables[i],t) = VectorsField[i]

RVV := 1001 :

r1 := kcatX*X*RVV/(kmX+X) :
r2 := kiXa*Xa :
r3 := kcatV*V*IIa/(kmV+V) :
r4 := kPT*Va*Xa*PL :
r5 := kPL*PT:
r6 := kcatII*II*PT/(kmII+II):
r7 := kcat2*II*Xa/(km2+II):
r8 := kiIIaa2M*IIa:
r9 := kiIIaATIII*IIa:

VectorField := [

	-r1,
	r1 - r2 -r4 + r5,
	-r3,
	r3-r4+r5,
	-r4+r5,
	r4-r5,
	-r6-r7,
	r6+r7-r8-r9,
	r9
] :

# We assume that OutputsVariables[i] = OutputSystem[i].
OutputSystem := [
       	IIa + 556*IIaa2M/1000
] :
#-----------------------------------------------------------------------------#
Inputs 		:= [] 							:
Parameters := [ kcatX,kmX,kiXa,kPT,kPL,kcatV,kmV,kcatII,kmII,kcat2,km2,
		kiIIaATIII,kiIIaa2M	] :
# The variables have to be ordered as the vectors field.
Variables  := [	X,Xa,V,Va,PL,PT,II,IIa,IIaa2M ] :
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





