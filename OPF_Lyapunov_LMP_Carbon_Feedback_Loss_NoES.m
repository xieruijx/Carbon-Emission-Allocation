clear
close all



M = 50;
load Grid.mat
load('CCarbon.mat','CCarbon');

delta = 0.0001;
Cgen_k1 = Cgen_k+delta*CCarbon;
kappa_carbon = 0.05;


pgSet = zeros(Ngen,Num_T);
prSet = zeros(NR,Num_T);
lambdaSet = zeros(Nbus,Num_T);
CPriceSet = zeros(Nbus,Num_T);

for t = 1:Num_T
    t

    pg = sdpvar(Ngen,1);
    pr = sdpvar(NR,1);

    obj_opf = sum(Cgen_k1'*pg*Delta_t);
    con1_opf = (1-losskappa)'*(Igen*pg+IR*pr-PD0(:,t)) == 0;
    con2_opfa = PTDF_l*(Igen*pg+IR*pr-PD0(:,t)) <= Sbranch_l;
    con2_opfb = -PTDF_l*(Igen*pg+IR*pr-PD0(:,t)) <= Sbranch_l;
    con3_opf = [PMINgen <= pg <= PMAXgen,0 <= pr <= PR(:,t)];
    con4_opf = [];
    con_opf = [con1_opf,con2_opfa,con2_opfb,con3_opf,con4_opf];
    optimize(con_opf,obj_opf,sdpsettings('verbose',0));
    
    lambdaSet(:,t) = (-dual(con1_opf)*(1-losskappa)-PTDF_l'*dual(con2_opfa)+PTDF_l'*dual(con2_opfb))/Delta_t;
    pgSet(:,t) = value(pg);
    prSet(:,t) = value(pr);

    CPriceSet(:,t) = calCPrice_PTDF(PD0(:,t),PR(:,t));
end

figure;
plot(1:Nbus,CPriceSet(:,50),'LineWidth',2);
xlabel('Bus');
ylabel('Emission price (kgCO_2/kWh)');
set(gca,'FontName','Times New Roman','FontSize',14);

CCostSet_PD = CPriceSet.*PD0(:,1:Num_T);
CCostAccumSys = zeros(Num_T+1,1);
for t = 1:Num_T
    CCostAccumSys(t+1) = CCostAccumSys(t)+2*sum(CCostSet_PD(:,t));
end
ab1 = polyfit((1:Num_T+1)',CCostAccumSys,1);

figure;
plot(0:Num_T,CCostAccumSys,'LineWidth',2);
xlabel('Time (h)');
ylabel('Cumulative system emission (10^5 kgCO_2');
set(gca,'FontName','Times New Roman','FontSize',14);

save('Result_Lyapunov_ES_Feedback_Loss_NoES.mat','PD0','pgSet','lambdaSet','PR','Num_T','Delta_t','CPriceSet','CCostAccumSys');