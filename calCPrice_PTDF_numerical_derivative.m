  function CPrice = calCPrice_PTDF_numerical_derivative(PD,PR,Sample)
load Grid.mat Nbus Ngen Nbranch NR PMINgen PMAXgen Igen IR BR_X Cgen_k PTDF_l Sbranch_l losskappa
load('CCarbon.mat','CCarbon');

CPriceSet = zeros(Nbus,Sample);
tSet = zeros(Sample,1);
for samplenum = 1:Sample
    t = samplenum/Sample
    Ccost0 = opf_numerical(PD*t,PR*t);
    CPricet = zeros(Nbus,1);
    for i = 1:Nbus
        pd = PD*t;
        pd(i) = PD(i)*(t-1/Sample);
        if PD(i) > 0
            CPricet(i) = (Ccost0-opf_numerical(pd,PR*t))/(PD(i)/Sample);
        end
    end

    CPriceSet(:,samplenum) = CPricet;
    tSet(samplenum) = t;
end

CPrice = CPriceSet*(tSet-[0;tSet(1:Sample-1)]) / 2;