#-----------------------------------------------------------------------------#
# 
# 
# 
#-----------------------------------------------------------------------------#
infolevel[observabilityTest] := 1 :
infolevel[observabilitySymmetries] := 1 :
#-----------------------------------------------------------------------------#
# We assume that diff(Variables[i],t) = VectorsField[i]

vX	:= vXm*(mS/V)/((mS/V)+KSX)*KSMX/((mSM/V)+KSMX)*mX	:
vSM	:= vSMm*(mS/V)/((mS/V)+KM+(mS/V)^2/KI)*mX		:

dmX	:= vX				:
dmS	:= -YSX*vX-YSSM*vSM-YErh*mX+u	:
dmSM	:= vSM				:
dV	:= u				:

VectorField:= [

    dmX,
    dmS,
    dmSM,
    dV

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
Variables 	:= [mX,mS,mSM,V]			:
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









