  function CPrice = calCPrice_PTDF_numerical_integral(PD,PR,Sample)
load Grid.mat Nbus Ngen Nbranch NR PMINgen PMAXgen Igen IR BR_X Cgen_k PTDF_l Sbranch_l losskappa
load('CCarbon.mat','CCarbon');
CCarbon = CCarbon(1:Ngen);
delta = 0.0001;
Cgen_k1 = Cgen_k+delta*CCarbon;
CCarbone = [CCarbon;zeros(NR,1)];
epsilon = 1e-10;
CPriceSet = zeros(Nbus,Sample);
tSet = zeros(Sample,1);
tNum = 0;
for samplenum = 1:Sample
    t = samplenum/Sample
    pd = PD*t;

    pg = sdpvar(Ngen,1);
    pr = sdpvar(NR,1);

    obj = Cgen_k1'*pg;
    con1 = (1-losskappa)'*(Igen*pg+IR*pr-pd) == 0;
    con2 = PTDF_l*(Igen*pg+IR*pr-pd) <= Sbranch_l;
    con3 = -PTDF_l*(Igen*pg+IR*pr-pd) <= Sbranch_l;
    con4 = PMINgen <= pg <= PMAXgen;
    con = [con1,con2,con3,con4,0<=pr<=PR];
%     optimize(con,obj);

    model = export(con,obj);
    Aeq = model.A(model.sense == '=',:);
    Aineq = model.A(model.sense == '<',:);
    Beq = model.rhs(model.sense == '=');
    Bineq = model.rhs(model.sense == '<');
    C = model.obj;
    BineqD = [PTDF_l;-PTDF_l;zeros(2*Ngen+2*NR,Nbus)];
    x = sdpvar(Ngen+NR,1);
    optimize([Aeq*x==Beq,Aineq*x<=Bineq],C'*x,sdpsettings('verbose',0));
    

    IneqEqIndex = abs(Aineq*value(x)-Bineq) < epsilon;
    Aeqt = [Aeq;Aineq(IneqEqIndex,:)];
    CPricet = [(1-losskappa)';BineqD(IneqEqIndex,:)]'*(Aeqt'\CCarbone);
%     x1 = [Aeq;Aineq(IneqEqIndex,:)]\([(1-losskappa)';BineqD(IneqEqIndex,:)]*PD);
% %     assert(sum(abs([Aeq;Aineq(IneqEqIndex,:)]*x1-([(1-losskappa)';BineqD(IneqEqIndex,:)]*PD))) < epsilon);
%     x2 = value(x);
%     
%     Nom = Bineq(~IneqEqIndex)-Aineq(~IneqEqIndex,:)*x2;
%     Denom = Aineq(~IneqEqIndex,:)*x1-BineqD(~IneqEqIndex,:)*PD;
% %     Denom = Aineq(~IneqEqIndex,:)*x1;
%     PIndex = (Denom > 0)&(Nom >= 0);
%     t = t+min(Nom(PIndex)./Denom(PIndex));

%     tNum = tNum+1;
    CPriceSet(:,samplenum) = CPricet(1:Nbus);
    tSet(samplenum) = t;
%     t = t+1/Sample;
end

% CPriceSet = CPriceSet(:,1:tNum);
% tSet = tSet(1:tNum);
% tSet(tNum) = 1;

CPrice = CPriceSet*(tSet-[0;tSet(1:Sample-1)]) / 2;