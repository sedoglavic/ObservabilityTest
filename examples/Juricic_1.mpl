 
#-----------------------------------------------------------------------------#
#                                                                       System 
# Description   : 
#Modèle 3CV4 (identifiable?)
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
		-1/(Ci*Ri)*x1 + 1/(Ci*Ri)*x3,
		-1/Cw*(1/Rw1 +1/Rw2)*x2 + 1/(Cw*Rw1)*x3 + 1/(Cw*Rw2)*u1,
		1/(Ca*Ri)*x1 + 1/(Ca*Rw1)*x2 - 1/Ca*(1/Ri + 1/Rw1 + zetaD/RD + 1/Rn 
			+ zetaW/RW +1/Rout)*x3 + 1/Ca*(zetaD/RD+1/Rn)*u1 
			+ 1/Ca*(zetaW/RW +1/Rout)*u2 + 1/Ca*u3
]:
# We assume that OutputsVariables[i] = OutputSystem[i].
OutputSystem := [

        x3
] :
#-----------------------------------------------------------------------------#
OutputsVariables:= [y]                                  :
Inputs          := [u1,u2,u3]:
Parameters      := [Ci,Ri,Cw,Rw1,Rw2,Ca,RD,Rout,RW, Rn, zetaD, zetaW]       :
# The variables have to be ordered as the vectors field.
Variables       := [x1,x2,x3]:

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

 
