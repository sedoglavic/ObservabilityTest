infolevel[observabilityTest] := 1 :
infolevel[observabilitySymmetries] := 1 :
M := [-(k10+k12+k13)*A1+k21*A2+k31*A3+u1, k12*A1-k21*A2, k13*A1-k31*A3] :
Y := [(1+B1*u2+B2+B3/(K+A1/V1))*A1/V1] :
X := [A1, A2, A3] :
U := [u1, u2] :
Param := [k12, k21, k13, k31, k10, V1, B1, B2, B3, K] :
# the following lines are related to the installation on my computer
libname := cat(currentdir(),"/release"),libname :
readlib(observabilityTest):
# and should not be used on others.

NonObs := observabilityTest(M, X, Y, Param, U) :


observabilitySymmetries(M,X,Y,Param,U,NonObs) :

quit:
