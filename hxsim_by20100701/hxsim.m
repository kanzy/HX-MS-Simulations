%%%2010-07-01 revised
%%%2010-04-27 hxsim.m: new basic HDX simulation subroutine for both PHX & NHX;
%%%use ODE to replace Mont Carlo method(the old one:2010-03-08
%%%phxsim_sim.m)

function obsPeaks=hxsim(START,END,Charge,Seq,pepPF,Temp,pH,fractionD,Time,iniD,obsType)

%%%START,END,Charge: peptide ID
%%%Seq:              sequence of the whole protein
%%%pepPF:            protection factors(could be <1) array of the peptide(First two sites and Proline could be input any number except 0)
%%%Temp,pH:          HX experiment conditions
%%%fractionD:        D2O percentage in HX sample
%%%Time:             HX time
%%%iniD:             initial D occupancy for all residues (1=allD, 0=allH //other cases should be writen as another program)

global hxodePara %store parameters for hxode.m use

%%%size check:
if max(size(pepPF))~=END-START+1
    error('START, END and pepPF not match!')
end

%%%data prep:
[peptideMass, distND, maxND, maxD]=pepinfo(Seq(START:END)); %call pepinfo.m
kcDH = fbmme_dh(Seq, pH, Temp, 'poly'); %call fbmme_dh.m
kcHD = fbmme_hd(Seq, pH, Temp, 'poly', 0); %call fbmme_hd.m

%%%calculate 'obsDistr' by convolution of ODE results:
hxodePara.fractionD=fractionD;

switch obsType
    case 1 %with natural isotopes
        obsDistr=distND;
    case 2 %without natural isotopes, just deuteron distribution
        obsDistr=1;
    otherwise
        error('Unknown "obsType"!')
end

j=1;
for i=(START+2):END
    hxodePara.kHD=kcHD(i)/pepPF(j);
    hxodePara.kDH=kcDH(i)/pepPF(j);
    [t,y] = ode15s(@hxode,[0 Time],[1-iniD iniD]); %use hxode.m
    occupyD=y(end,2);
    if Seq(i)=='P'
        occupyD=0; %Proline must be allH
    end
    obsDistr=conv(obsDistr, [1-occupyD, occupyD]);
    j=j+1;
end

%%%normalization:
if size(obsDistr,2)>=(1+maxND+maxD)
    obsDistr=obsDistr(1:1+maxND+maxD)/sum(obsDistr(1:1+maxND+maxD));
else
    obsDistr=[obsDistr zeros(1, (1+maxND+maxD)-size(obsDistr,2))]/sum(obsDistr);
end

%%%get MS observable peaks:
DM=1.00628; %delta mass between a deuteron(2.014102 u) and a proton(1.007825 u)
obsPeaks(1,1)=peptideMass/Charge+1.007276; %m/z of mono; 1.007276 is the mass of proton
obsPeaks(1,2)=obsDistr(1);
for i=2:(maxD+maxND+1)
    obsPeaks(i,1)=obsPeaks(i-1,1)+DM/Charge;
    obsPeaks(i,2)=obsDistr(i);
end


