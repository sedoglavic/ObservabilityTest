#-----------------------------------------------------------------------------#
# 
# 
# 
#-----------------------------------------------------------------------------#
infolevel[observabilityTest] := 1 :
infolevel[observabilitySymmetries] := 1 :
#-----------------------------------------------------------------------------#
# We assume that diff(Variables[i],t) = VectorsField[i]

# Aiba = exp(-KSMX*(mSM/V))

vX	:= vXm*(mS/V)/((mS/V)+KSX)*Aiba*mX		:
vSM	:= vSMm*(mS/V)/((mS/V)+KM+(mS/V)^2/KI)*mX	:

dmX	:= vX				:
dmS	:= -YSX*vX-YSSM*vSM-YErh*mX+u	:
dmSM	:= vSM				:
dV	:= u				:
dAiba   := KSMX*(mSM*dV/V^2-dmSM/V)*Aiba :

VectorField:= [

    dmX,
    dmS,
    dmSM,
    dV,
    dAiba 

]:     	
# We assume that OutputsVariables[i] = OutputSystem[i].
OutputSystem := [

(mX/V),(mS/V),(mSM/V)

] :
#-----------------------------------------------------------------------------#
OutputsVariables:= [y1,y2,y3] 					:
Inputs := [u] 					:
Parameters 	:= [vXm,KSX,KSMX,vSMm,KM,KI,YSX,YSSM,YErh]		:
# The variables have to be ordered as the vectors field.
Variables 	:= [mX,mS,mSM,V,Aiba]			:
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









