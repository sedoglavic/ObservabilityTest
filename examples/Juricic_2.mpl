#-----------------------------------------------------------------------------#
#                                                                       System 
# Description   : 
 #Modèle 4C (non identifiable?)
# Result        :
#
#-----------------------------------------------------------------------------#
infolevel[observabilityTest] := 1 :
infolevel[observabilitySymmetries] := 1 :
#-----------------------------------------------------------------------------#
# Bibliography : Milena Anguelova, private communication 2006
# This examples is observable while previous code say opposite. The bug is know
# fixed. 
#-----------------------------------------------------------------------------#
# We assume that diff(Variables[i],t) = VectorsField[i]
VectorField:= [
 -1/(Ci*Ri)*x1 + 1/(Ci*Ri)* x4,
-1/Cwout*(1/Rwout1 +1/Rwout2)*x2 + 1/(Cwout*Rwout1)*x4 + 1/(Cwout*Rwout2)*u2,
-1/Cwn*(1/Rwn1 +1/Rwn2)*x3 + 1/(Cwn*Rwn1)*x4 + 1/(Cwn*Rwn2)*u1,
1/(Ca*Ri)*x1 + 1/(Ca*Rwout1)*x2 + 1/(Ca*Rwn1)*x3 - 1/Ca*(1/Ri + 1/Rwn1 + 1/Rwout1 + zetaD/RD + 1/Rn + zetaW/RW +1/Rout)*x4 + 1/Ca*(zetaD/RD+1/Rn)*u1 + 1/Ca*(zetaW/RW +1/Rout)*u2 + 1/Ca*u3
]:
# We assume that OutputsVariables[i] = OutputSystem[i].
OutputSystem := [

        x4
] :
#-----------------------------------------------------------------------------#
OutputsVariables:= [y]                                  :
Inputs          := [u1,u2,u3,Rn]:
Parameters      := [Ci,Ri,Cwn,Rwn1,Rwn2,Cwout,Rwout1,Rwout2,Ca,RD,Rout,RW,zetaD,zetaW]:
# The variables have to be ordered as the vectors field.
Variables       := [x1,x2,x3,x4]:

#-----------------------------------------------------------------------------#
# CAUTION read section II of the file INSTALL
libname := cat(currentdir(),"/release"),libname :
readlib(observabilityTest) :
NonObservable := observabilityTest(     VectorField     ,
                                        Variables       ,
                                        OutputSystem    ,
                                        Parameters      ,
                                        Inputs                  ) :
print(%) :                                      
GroupInfGen := observabilitySymmetries( VectorField     ,
                                        Variables       ,
                                        OutputSystem    ,
                                        Parameters      ,
                                        Inputs          ,
                                        NonObservable           ) :
print(%) :
#-----------------------------------------------------------------------------#
quit :

 
