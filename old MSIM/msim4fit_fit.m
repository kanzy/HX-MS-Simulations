
function f = msim4fit_fit(x,P1,P2,P3)



START=P1(1);  %start residue number of the peptide
END=P1(2);

T=P1(3);   %HX time (min)

mean=P2;

kc=P3;

Y=0;
for i=(START+2):END   
    Y=Y+1-exp(-kc(i)*T/x);
end

f=Y-mean; %6.4 is the Mean of a subpopulation from double-Gaussian fitting //! should replaced by 6.4/(allD/maxD) for real data