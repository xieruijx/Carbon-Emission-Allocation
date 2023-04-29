clear
close all
load Grid.mat

tNum = 50;

%% Proposed
tic
[CPrice,SNum] = calCPrice_PTDF_num(PD0(:,tNum),PR(:,tNum));
toc

load('CCarbon.mat','CCarbon');
delta = 0.0001;
Cgen_k1 = Cgen_k+delta*CCarbon;

pg = sdpvar(Ngen,1);
pr = sdpvar(NR,1);

obj_opf = sum(Cgen_k1'*pg*Delta_t);
con1_opf = (1-losskappa)'*(Igen*pg+IR*pr-PD0(:,tNum)) == 0;
con2_opfa = PTDF_l*(Igen*pg+IR*pr-PD0(:,tNum)) <= Sbranch_l;
con2_opfb = -PTDF_l*(Igen*pg+IR*pr-PD0(:,tNum)) <= Sbranch_l;
con3_opf = [PMINgen <= pg <= PMAXgen,0 <= pr <= PR(:,tNum)];
con4_opf = [];
con_opf = [con1_opf,con2_opfa,con2_opfb,con3_opf,con4_opf];
optimize(con_opf,obj_opf,sdpsettings('verbose',0));

PG = value(pg);

(abs((CCarbon')*PG/2-sum(CPrice.*(PD0(:,tNum)))))/((CCarbon')*PG/2)

%% Numerical derivative and integral

Sample = 10;

tic
CPrice = calCPrice_PTDF_numerical_derivative(PD0(:,tNum),PR(:,tNum),Sample);
toc

(abs((CCarbon')*PG/2-sum(CPrice.*(PD0(:,tNum)))))/((CCarbon')*PG/2)


%% Numerical integral

Sample = 10;

tic
CPrice = calCPrice_PTDF_numerical_integral(PD0(:,tNum),PR(:,tNum),Sample);
toc

(abs((CCarbon')*PG/2-sum(CPrice.*(PD0(:,tNum)))))/((CCarbon')*PG/2)