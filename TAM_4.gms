Set
   t 'time'               / 1*24   /
   g 'generators indices' / g1*g10  /
   k 'cost segments'      / sg1*sg20 /
   bus        / 1*39   /
   slack(bus) / 1     /
   t1(t) /1*9/
   t2(t) /12*24/
   t3(t) /10*12/
   c(bus) /3,39/
   df(bus) /4,8,20/
   sm /sm1*sm5/
   el(bus) /20*29/
;
Alias (t,tt);
Alias (bus,node);

Set GB(bus,g) 'connectivity index of each generating unit to each bus'
/
30.g1
31.g2
32.g3
33.g4
34.g5
35.g6
36.g7
37.g8
38.g9
39.g10
/
;


parameter  gendata(g,*) Information on the Generators ;
$call GDXXRW Udata.xlsx maxDupeErrors=100000000 par=gendata rng=Gen!a1
$GDXIN Udata.gdx
$load gendata
$GDXIN
;

parameter  Xij(bus,node) Inductance of transmission lines  ;
$call GDXXRW Udata.xlsx maxDupeErrors=100000000 par=Xij rng=Xij!a1
$GDXIN Udata.gdx
$load Xij
$GDXIN
;

parameter  Bij(bus,node) Suspedance of transmission lines ;
$call GDXXRW Udata.xlsx maxDupeErrors=100000000 par=Bij rng=Bij!a1
$GDXIN Udata.gdx
$load Bij
$GDXIN
;

parameter  Limit(bus,node) Transmission capacity ;
$call GDXXRW Udata.xlsx maxDupeErrors=100000000 par=Limit rng=Limit!a1
$GDXIN Udata.gdx
$load Limit
$GDXIN
;

parameter  Conex(bus,node) network characteristics ;
$call GDXXRW Udata.xlsx maxDupeErrors=100000000 par=Conex rng=Conex!a1
$GDXIN Udata.gdx
$load Conex
$GDXIN
;


parameter  WD(t,*) Demand proflie ;
$call GDXXRW Udata.xlsx maxDupeErrors=100000000 par=WD rng=Demand!a1
$GDXIN Udata.gdx
$load WD
$GDXIN
;

parameter  BusData(bus,*) Nodal loads ;
$call GDXXRW Udata.xlsx maxDupeErrors=100000000 par=BusData rng=Load!a1
$GDXIN Udata.gdx
$load BusData
$GDXIN
;



Parameter PL(bus,t);
PL(bus,t) = WD(t,'d')*BusData(bus,'Pd');

Scalar
Voll /500/
;

parameter  Price(t,*) Price ;
$call GDXXRW Udata.xlsx maxDupeErrors=100000000 par=Price rng=Price!a1
$GDXIN Udata.gdx
$load Price
$GDXIN
;

parameter  bid(bus,sm) Price2 ;
$call GDXXRW Udata.xlsx maxDupeErrors=100000000 par=bid rng=PC!a1
$GDXIN Udata.gdx
$load bid
$GDXIN
;

parameter  bid2(bus,sm) Price3 ;
$call GDXXRW Udata.xlsx maxDupeErrors=100000000 par=bid2 rng=Pdf!a1
$GDXIN Udata.gdx
$load bid2
$GDXIN
;

parameter  MUU(bus,*) Price4 ;
$call GDXXRW Udata.xlsx maxDupeErrors=100000000 par=MUU rng=MU!a1
$GDXIN Udata.gdx
$load MUU
$GDXIN
;




Variable
OF
delta0(bus,t)
Pij(bus,node,t)
Pc(bus,t)
Pdf(bus,t)
;

Positive Variable
Pg(g,t)
Phs(bus,t)
DL(bus,t)
Del(bus,t)
;
Binary Variable
u(g,t)
;



Equation
obj
eq1
eq2
eq3
eq4
eq5
eq6
eq7

;

obj..                    OF=e=sum(t, sum(bus, VOLL*PL(bus,t)  ) + sum(bus$el(bus),MUU(bus,'MU')*Del(bus,t)) - sum(g,gendata(g,'MC')*Pg(g,t)) ) ;

eq1(g,t)..               Pg(g,t)=l=gendata(g,'Pmax')*u(g,t);

eq2(g,t)..               Pg(g,t)=g=gendata(g,'Pmin')*u(g,t);

eq3(bus,node,t)..        Pij(bus,node,t) =e= Bij(bus,node)*(delta0(bus,t) - delta0(node,t));

eq4(bus,node,t)$Conex(bus,node)..         (Bij(bus,node)*(delta0(bus,t) - delta0(node,t)))=l= Limit(bus,node);

eq5(bus,node,t)$Conex(bus,node)..         (Bij(bus,node)*(delta0(bus,t) - delta0(node,t)))=g=-Limit(bus,node);

eq6(bus,t)..           DL(bus,t) =e= PL(bus,t) + Del(bus,t)$el(bus);

eq7(bus,t)..           sum(g$GB(bus,g), pg(g,t)) - DL(bus,t) =e= sum(node$conex(node,bus),Pij(bus,node,t) );






delta0.up(bus,t)   = pi/2*100;
delta0.lo(bus,t)   =-pi/2*100;
delta0.fx(slack,t) = 0;

model ENGO1 /all/;
option limrow=100;
Option Optcr= 0;
Option Optca= 0.0;
OPTION RESLIM = 100000;
Option solver = CPLEX;
solve ENGO1 using miqcp max OF;

Parameter DL_tot,DLmax,TDpay,TGEN_prof,PLT(t),D_avg,LF,LMP(bus,t),Tcost ,Res_ISO       ;

PLT(t) = sum(bus,PL(bus,t));

DL_tot = sum((bus,t),DL.l(bus,t));

DLmax = 6254.23;

D_avg = sum(t,PLT(t))/card(t);

LF = DLmax/D_avg ;

LMP(bus,t) =  eq7.m(bus,t) ;

Tcost =  sum(t,sum(g,gendata(g,'a')*sqr(Pg.l(g,t))+Pg.l(g,t)*gendata(g,'b')+gendata(g,'a')*u.l(g,t))) ;

TDpay = sum(bus, sum(t$t1(t),LMP(bus,t)*DL.l(bus,t)) );

TGEN_prof = sum(bus, sum(t,LMP(bus,t)*DL.l(bus,t)) - (sum(t,sum(g,gendata(g,'a')*sqr(Pg.l(g,t))+Pg.l(g,t)*gendata(g,'b')+gendata(g,'a')*u.l(g,t)))) );

Res_ISO = TDpay - Tcost;

Display
DL_tot
TDpay
TGEN_prof
LF
Res_ISO
;

execute_unload "Tamrin4.gdx"
Pg.l
Phs.l
DL.l
DL_tot
*DLmax
TDpay
TGEN_prof
OF.l
delta0.l
Pij.l

execute 'gdxxrw.exe Tamrin4.gdx var= Pg          rng=sheet1!a1'
execute 'gdxxrw.exe Tamrin4.gdx var= Phs         rng=sheet2!a1'
execute 'gdxxrw.exe Tamrin4.gdx var= DL          rng=sheet3!a1'
execute 'gdxxrw.exe Tamrin4.gdx par= DL_tot      rng=sheet4!a1'
*execute 'gdxxrw.exe Tamrin4.gdx var= DLmax       rngsheet5 !a1'
execute 'gdxxrw.exe Tamrin4.gdx par= TDpay       rng=sheet6!a1'
execute 'gdxxrw.exe Tamrin4.gdx par= TGEN_prof   rng=sheet7!a1'
execute 'gdxxrw.exe Tamrin4.gdx var= OF          rng=sheet8!a1'
execute 'gdxxrw.exe Tamrin4.gdx var= delta0      rng=sheet9!a1'
execute 'gdxxrw.exe Tamrin4.gdx var= Pij         rng=sheet10!a1'
