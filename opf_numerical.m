function Ccost = opf_numerical(pd,PR)
load Grid.mat Nbus Ngen Nbranch NR PMINgen PMAXgen Igen IR BR_X Cgen_k PTDF_l Sbranch_l losskappa
load('CCarbon.mat','CCarbon');
CCarbon = CCarbon(1:Ngen);
delta = 0.0001;
Cgen_k1 = Cgen_k+delta*CCarbon;

pg = sdpvar(Ngen,1);
pr = sdpvar(NR,1);

obj = Cgen_k1'*pg;
con1 = (1-losskappa)'*(Igen*pg+IR*pr-pd) == 0;
con2 = PTDF_l*(Igen*pg+IR*pr-pd) <= Sbranch_l;
con3 = -PTDF_l*(Igen*pg+IR*pr-pd) <= Sbranch_l;
con4 = PMINgen <= pg <= PMAXgen;
con = [con1,con2,con3,con4,0<=pr<=PR];

optimize(con,obj,sdpsettings('verbose',0));
Ccost = CCarbon'*value(pg);