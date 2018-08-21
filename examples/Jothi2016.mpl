#Current Research Work:
#Observability rank condition on standard vehicle models to calculate the best sensor type and position in order to estimate road and rail roughness from the response of the vehicles.
#I have implemented both Linear and Nonlinear Observability analysis for many small dynamic models in #MATLAB. But I could not implement for the 27 Degree of Freedom system because the Memory is not sufficient for Symbolic Computation. Prof. Manolis from University of Oxford, UK suggested to use your method for Symbolic Computation of Matrices and Rank calculation for Observability analysis.
# kindly request you to help me in implementing the following cases in your code. I am writing both MATLAB #version and Maple version (with my little understanding)

VectorField := [ 
     dyc1,
     dzc1,
     dtc1,
     dsc1,
     dpc1,
     dyb11,
     dzb11,
     dtb11,
     dsb11,
     dpb11,
     dyb12,
     dzb12,
     dtb12,
     dsb12,
     dpb12,
     dyw111,
     dzw111,
     dtw111,
     dyw112,
     dzw112,
     dtw112,
     dyw121,
     dzw121,
     dtw121,
     dyw122,
     dzw122,
     dtw122,
    (1/Mc1)*(-2*ksv*yc1+(ksv)*yb11+(ksv)*yb12-2*Csv*dyc1+(Csv)*dyb11+(Csv)*dyb12),
(1/Mc1)*(-2*ksh*hs*zc1+2*ksh*hs*tc1+(ksh)*zb11+(ksh*hb)*tb11+(ksh)*zb12+(ksh*hb)*tb12-2*Csh*hs*dzc1+2*Csh*hs*dtc1+(Csh)*dzb11+(Csh*hb)*dtb11+(Csh)*dzb12+(Csh*hb)*dtb12),
    (1/Jct1)*(2*ksh*hs*zc1-(2*bs^2*ksv+2*hs^2*ksh)*tc1-(ksh*hs)*zb11-(hs*hb*ksh-bs^2*ksv)*tb11-(ksh*hs)*zb12-(hs*hb*ksh-bs^2*ksv)*tb12+2*Csh*hs*dzc1-(2*bs^2*Csv+2*hs^2*Csh)*dtc1-(Csh*hs)*dzb11-(hs*hb*Csh-bs^2*Csv)*dtb11-(Csh*hs)*dzb12-(hs*hb*Csh-bs^2*Csv)*dtb12),
    (1/Jcs1)*(-2*lb^2*ksh*sc1-(lb*ksh)*zb11-(lb*hb*ksh)*tb11+(lb*ksh)*zb12+(lb*hb*ksh)*tb12-2*lb^2*Csh*dsc1-(lb*Csh)*dzb11-(lb*hb*Csh)*dtb11+(lb*Csh)*dzb12+(lb*hb*Csh)*dtb12),
    (1/Jcp1)*(-2*lb^2*ksv*pc1+(lb*ksv)*yb11+(lb*ksv)*yb12-2*lb^2*Csv*dpc1+(lb*Csv)*dyb11+(lb*Csv)*dyb12),
    (1/Mb11)*(ksv*yc1+lb*ksv*pc1-(ksv+2*kpv)*yb11+kpv*yw111+kpv*yw112+Csv*dyc1+lb*Csv*dpc1-(Csv+2*Cpv)*dyb11+Cpv*dyw111+Cpv*dyw112),
    (1/Mb11)*(ksh*zc1-ksh*hs*tc1-lb*ksh*sc1-(ksh+2*kph)*zb11-(ksh*hb-2*kph*hp)*tb11+kph*zw111+kph*yw112+Csh*dzc1-Csh*hs*dtc1-lb*Csh*dsc1-(Csh+2*Cph)*dzb11-(Csh*hb-2*Cph*hp)*dtb11+Cph*dyw111+Cph*dyw112),
    (1/Jbt11)*(ksh*hb*zc1-(hs*hb*ksh-bs^2*ksv)*tc1-lb*hb*ksh*sc1-(ksh*hb-2*kph*hp)*zb11-(bs^2*ksv+hb^2*ksh+2*bp^2*kpv+2*hp^2*kph)*tb11-kph*hp*zw111+bp^2*kpv*tw111-kph*hp*zw112+bp^2*kpv*tw112+Csh*hb*dzc1-(hs*hb*Csh-bs^2*Csv)*dtc1-lb*hb*Csh*dsc1-(Csh*hb-2*Cph*hp)*dzb11-(bs^2*Csv+hb^2*Csh+2*bp^2*Cpv+2*hp^2*Cph)*dtb11-Cph*hp*dzw111+bp^2*Cpv*dtw111-Cph*hp*dzw112+bp^2*Cpv*dtw112), 
    (1/Jbs11)*(-2*kph*lr^2*sb11-kph*lr*zw111+kph*lr*zw112-2*Cph*lr^2*dsb11-Cph*lr*dzw111+Cph*lr*dzw112),
    (1/Jbp11)*(-2*kpv*lr^2*pb11+lr*kpv*yw111-lr*kpv*yw112-2*Cpv*lr^2*dpb11+lr*Cpv*dyw111-lr*Cpv*dyw112),
    (1/Mb12)*(ksv*yc1+ksv*pc1-(ksv+2*kpv)*yb12+kpv*yw121+kpv*yw122+Csv*dyc1+lb*Csv*dpc1-(Csv+2*Cpv)*dyb12+Cpv*dyw121+Cpv*dyw122),
    (1/Mb12)*(ksh*zc1-ksh*hs*tc1+lb*ksh*sc1-(ksh+2*kph)*zb12-(ksh*hb-2*kph*hp)*tb12+kph*zw121+kph*zw122+Csh*dzc1-Csh*hs*dtc1+lb*Csh*dsc1-(Csh+2*Cph)*dzb12-(Csh*hb-2*Cph*hp)*dtb12+Cph*dzw121+Cph*dzw122),
    (1/Jbt12)*(ksh*hb*zc1-(hs*hb*ksh-bs^2*ksv)*tc1+lb*hb*ksh*sc1-(ksh*hb-2*kph*hp)*zb12-(bs^2*ksv+hb^2*ksh+2*bp^2*kpv+2*hp^2*kph)*tb12-kph*hp*zw121+bp^2*kpv*tw121-kph*hp*zw122+bp^2*kpv*tw122+Csh*hb*dzc1-(hs*hb*Csh-bs^2*Csv)*dtc1+lb*hb*Csh*dsc1-(Csh*hb-2*Cph*hp)*dzb12-(bs^2*Csv+hb^2*Csh+2*bp^2*Cpv+2*hp^2*Cph)*dtb12-Cph*hp*dzw121+bp^2*Cpv*dtw121-Cph*hp*dzw122+bp^2*Cpv*dtw122),
    (1/Jbs12)*(-2*kph*lr^2*sb12-kph*lr*zw121+kph*lr*zw122-2*Cph*lr^2*dsb12-Cph*lr*dzw121+Cph*lr*dzw122),
    (1/Jbp12)*(-2*kpv*lr^2*pb12+lr*kpv*yw121-lr*kpv*yw122-2*Cpv*lr^2*dpb12+lr*Cpv*dyw121-lr*Cpv*dyw122),
    (1/Mw111)*((kpv)*yb11+(lr*kpv)*pb11-kpv*yw111+(Cpv)*dyb11+(lr*Cpv)*dpb11-Cpv*dyw111),
    (1/Mw111)*((kph)*zb11-(kph*hp)*tb11-(kph*lr)*sb11-kph*zw111+(Cph)*dzb11-(Cph*hp)*dtb11-(Cph*lr)*dsb11-Cph*dzw111),
    (1/Jwt111)*((bp^2*kpv)*tb11-bp^2*kpv*tw111+(bp^2*Cpv)*dtb11-bp^2*Cpv*dtw111),
    (1/Mw112)*((kpv)*yb11-(lr*kpv)*pb11-kpv*yw112+(Cpv)*dyb11-(lr*Cpv)*dpb11-Cpv*dyw112),
    (1/Mw112)*((kph)*zb11-(kph*hp)*tb11+(kph*lr)*sb11-kph*zw112+(Cph)*dzb11-(Cph*hp)*dtb11+(Cph*lr)*dsb11-Cph*dzw112),
    (1/Jwt112)*((bp^2*kpv)*tb11-bp^2*kpv*tw112+(bp^2*Cpv)*dtb11-bp^2*Cpv*dtw112),
    (1/Mw121)*((kpv)*yb12+(lr*kpv)*pb12-kpv*yw121+(Cpv)*dyb12+(lr*Cpv)*dpb12-Cpv*dyw121),
    (1/Mw121)*((kph)*zb12-(kph*hp)*tb12-(kph*lr)*sb12-kph*zw121+(Cph)*dzb12-(Cph*hp)*dtb12-(Cph*lr)*dsb12-Cph*dzw121),
    (1/Jwt121)*((bp^2*kpv)*tb12-bp^2*kpv*tw121+(bp^2*Cpv)*dtb12-bp^2*Cpv*dtw121),
    (1/Mw122)*((kpv)*yb12-(lr*kpv)*pb12-kpv*yw122+(Cpv)*dyb12-(lr*Cpv)*dpb12-Cpv*dyw122),
    (1/Mw122)*((kph)*zb12-(kph*hp)*tb12+(kph*lr)*sb12-kph*zw122+(Cph)*dzb12-(Cph*hp)*dtb12+(Cph*lr)*dsb12-Cph*dzw122),
    (1/Jwt122)*((bp^2*kpv)*tb12-bp^2*kpv*tw122+(bp^2*Cpv)*dtb12-bp^2*Cpv*dtw122)]:

#Parameters (mass, stiffness, damping of the dynamic systems) are assumed to be knowm.

#Inputs are already included in the state vector as a state variable.

h1:= (1/Mc1)*(-2*ksv*yc1+(ksv)*yb11+(ksv)*yb12-2*Csv*dyc1+(Csv)*dyb11+(Csv)*dyb12);

h2:= (1/Mc1)*(-2*ksh*hs*zc1+2*ksh*hs*tc1+(ksh)*zb11+(ksh*hb)*tb11+(ksh)*zb12+(ksh*hb)*tb12-2*Csh*hs*dzc1+2*Csh*hs*dtc1+(Csh)*dzb11+(Csh*hb)*dtb11+(Csh)*dzb12+(Csh*hb)*dtb12);
h3:=dtc1;
h4:=dsc1;
h5:=dpc1;


OutputSystem:= [h1,h2,h3,h4,h5]:
Inputs:=[]:
OutputsVariables:=[ o1,o2,o3,o4,o5] : 
Inputs := []:

Parameters := [ Cph, Cpv, Csh, Csv, Jbp11, Jbp12, Jbs11, Jbs12, Jbt11, Jbt12, Jcp1, Jcs1, Jct1, Jwt111, Jwt112, Jwt121, Jwt122, Mb11, Mb12, Mc1, Mw111, Mw112, Mw121, Mw122, bp, bs, hb, hp, hs, kph, kpv, ksh, ksv, lb, lr]:

Variables := [yc1,zc1,tc1,sc1,pc1,yb11,zb11,tb11,sb11,pb11,yb12,zb12,tb12,sb12,pb12,yw111,zw111,tw111,yw112,zw112,tw112,yw121,zw121,tw121,yw122,zw122,tw122,dyc1,dzc1,dtc1,dsc1,dpc1,dyb11,dzb11,dtb11,dsb11,dpb11,dyb12,dzb12,dtb12,dsb12,dpb12,dyw111,dzw111,dtw111,dyw112,dzw112,dtw112,dyw121,dzw121,dtw121,dyw122,dzw122,dtw122]:
libname := cat(currentdir(),"/release"),libname :

readlib(observabilityTest):

NonObservable := observabilityTest(VectorField, Variables, OutputSystem, Parameters, Inputs):

