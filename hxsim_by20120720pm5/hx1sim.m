%%%2012-07-20 hx1sim.m: for EX1 HX-MS peptide spectrum (simplified) simulation; refer to hxsim.m


function [obsPeaks, deuPeaks, Distr1,Fract1,Distr2,Fract2]=hx1sim(pepPara,hxPara)

proSeq=pepPara.proSeq; %whole protein sequence
START=pepPara.START; %Start residue# of the peptide
END=pepPara.END; %End residue# of the peptide
Charge=pepPara.Charge; %observing charge state of the peptide in MS spectrum
pepPF=pepPara.pepPF; %protection factors(could be <1) array of the peptide(First two sites and Proline could be input any number except 0)
foldIndex=pepPara.foldIndex; %array of foldon index of each residue(0=uncorrelated; 1,2,3...=diff foldons): residues in the same foldon must open/close together
kOPC=pepPara.kOPC; %opening rate of the correlated residue set
iniD=pepPara.iniD; %D% of every residue of the peptide at the beginning of HX

Temp=hxPara.Temp; %HX temperature (unit: 'C)
pH=hxPara.pH; %HX pH
hxTime=hxPara.hxTime; %HX duration time (unit: sec)
fractionD=hxPara.fractionD; %fraction of D2O in HX buffer
hxDir=hxPara.hxDir; %HX direction (1=H->D, 2=D->H)

%%%size check:
N=END-START+1; %number of residues of the peptide
if size(pepPF,2)~=N || size(foldIndex,2)~=N
    error('Input sizes not match!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%peptide data prep:
kcDH_pro = fbmme_dh(proSeq, pH, Temp, 'poly'); %call fbmme_dh.m
kcHD_pro = fbmme_hd(proSeq, pH, Temp, 'poly', 0); %call fbmme_hd.m
kcDH=kcDH_pro(START:END);
kcHD=kcHD_pro(START:END);

[peptideMass, distND, maxND, maxD]=pepinfo(proSeq(START:END),2); %call pepinfo.m


%%%calculate 'obsDistr':
Distr=1;
ct=0;
for i=(START+2):END
    i
    if foldIndex(i-START+1)==0 %for uncorrelated residues, apply binormial
        if hxDir==1 %H->D
            kex=kcHD(i-START+1)/pepPF(i-START+1);
        else
            kex=kcDH(i-START+1)/pepPF(i-START+1);
        end
        occupyD=(iniD-fractionD)*exp(-kex*hxTime)+fractionD;
        if proSeq(i)=='P'
            occupyD=0; %Proline must be allH
        end
        Distr=conv(Distr, [1-occupyD, occupyD]);
    else
        if proSeq(i)~='P'
            ct=ct+1
        end
    end
end

Distr1=Distr;
Fract1=exp(-kOPC*hxTime);
if hxDir==1 %H->D
    Distr2=[zeros(1,ct), Distr(1:end-ct)];
else %D->H
    Distr2=[Distr(1+ct:end), zeros(1,ct)];
end
Fract2=(1-exp(-kOPC*hxTime));

Distr=Fract1*Distr1 + Fract2*Distr2;

%%%normalization:
if size(Distr,2)>=(1+maxND+maxD)
    Distr=Distr(1:1+maxND+maxD)/sum(Distr(1:1+maxND+maxD));
else
    Distr=[Distr zeros(1, (1+maxND+maxD)-size(Distr,2))]/sum(Distr);
end


%%%get MS observable peaks and plot:
figure

subplot(2,1,1)
obsDistr=conv(distND, Distr);
obsDistr=obsDistr/sum(obsDistr); %normalization
DM=1.00628; %delta mass between a deuteron(2.014102 u) and a proton(1.007825 u)
clear obsPeaks
obsPeaks(1,1)=peptideMass/Charge+1.007276; %m/z of mono; 1.007276 is the mass of proton
obsPeaks(1,2)=obsDistr(1);
for i=2:(maxD+maxND+1)
    obsPeaks(i,1)=obsPeaks(i-1,1)+DM/Charge;
    obsPeaks(i,2)=obsDistr(i);
end
stem(obsPeaks(:,1), obsPeaks(:,2))
v=axis;
axis([obsPeaks(1,1)-0.5/Charge, obsPeaks(end,1)+0.5/Charge, v(3), v(4)])
xlabel('m/z')

subplot(2,1,2)
clear deuPeaks
deuPeaks(1,1)=0;
deuPeaks(1,2)=Distr(1);
for i=2:(maxD+1)
    deuPeaks(i,1)=deuPeaks(i-1,1)+1;
    deuPeaks(i,2)=Distr(i);
end
stem(deuPeaks(:,1), deuPeaks(:,2))
v=axis;
axis([-0.5, maxD+0.5, v(3), v(4)])
xlabel('Deuteron Number')









