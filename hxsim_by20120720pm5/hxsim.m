%%%2010-09-20: revised to universal use (both EX1 and EX2)
%%%2010-07-01&04 revised
%%%2010-04-27 hxsim.m: new basic HDX simulation subroutine for both PHX & NHX;
%%%use ODE to replace Mont Carlo method(the old one:2010-03-08
%%%phxsim_sim.m)

function [HDmatrix, Fmatrix, obsPeaks, deuPeaks]=hxsim(pepPara,hxPara)


proSeq=pepPara.proSeq; %whole protein sequence
START=pepPara.START; %Start residue# of the peptide
END=pepPara.END; %End residue# of the peptide
Charge=pepPara.Charge; %observing charge state of the peptide in MS spectrum
kOP=pepPara.kOP; %array of openning rates of each residue of the peptide at HX condition
kCL=pepPara.kCL; %array of closing rates of each residue of the peptide at HX condition
iniD=pepPara.iniD; %array of D% of each residue of the peptide at the beginning of HX
iniF=pepPara.iniF; %array of folding status(1=fully unfolded(open); 2=fully folded(close); 3=equilibrium by HX condition--assuming independent of other residues) of each residue of the peptide at the beginning of HX
foldIndex=pepPara.foldIndex; %array of foldon index of each residue(0=uncorrelated; 1,2,3...=diff foldons): residues in the same foldon must open/close together


Temp=hxPara.Temp; %HX temperature (unit: 'C)
pH=hxPara.pH; %HX pH
hxTime=hxPara.hxTime; %HX duration time (unit: sec)
fractionD=hxPara.fractionD; %fraction of D2O in HX buffer

%%%size check:
N=END-START+1; %number of residues of the peptide
if size(kOP,2)~=N || size(kCL,2)~=N || size(iniD,2)~=N || size(iniF,2)~=N
    error('Input sizes not match!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%peptide data prep:
kcDH_pro = fbmme_dh(proSeq, pH, Temp, 'poly'); %call fbmme_dh.m
kcHD_pro = fbmme_hd(proSeq, pH, Temp, 'poly', 0); %call fbmme_hd.m
kcDH=kcDH_pro(START:END);
kcHD=kcHD_pro(START:END);

[peptideMass, distND, maxND, maxD]=pepinfo(proSeq(START:END),2); %call pepinfo.m


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%SIMULATION PART:

M=1000; %number of simulating molecules

%calculate step time of simulation:
kmax=max(max(kcDH), max(kcHD));
deltaT=1/kmax;

%initialize 'HDmatrix'(0=H; 1=D)
HDmatrix=zeros(M,N);
for i=1:M
    HDmatrix(i,:)=iniD;
end

%initialize 'Fmatrix' (1=open(unfolded); 2=closed(folded)) //this assume no cooperativity among residues, otherwise should use 'Farray' to track diff molecules
Fmatrix=zeros(M,N);
for i=1:M
    for j=1:N
        switch iniF(j)
            case 1
                Fmatrix(i,j)=1;
            case 2
                Fmatrix(i,j)=2;
            case 3
                if rand>kCL(j)/(kCL(j)+kOP(j))
                    Fmatrix(i,j)=1;
                else
                    Fmatrix(i,j)=2;
                end
        end
    end
end

%start simulation:
for time=0:deltaT:hxTime
    ran1=rand(M,N);
    ran2=rand(M,N);
    for i=1:M
        for j=1:N
            %%%consider 'HDmatrix' transition:
            switch HDmatrix(i,j)
                case 0 %H
                    if ran1(i,j)<=(1-exp(-kcHD(j)*deltaT))*fractionD && Fmatrix(i,j)==1
                        HDmatrix(i,j)=1; %H->D
                    end
                case 1 %D
                    if ran1(i,j)<=(1-exp(-kcDH(j)*deltaT))*(1-fractionD) && Fmatrix(i,j)==1
                        HDmatrix(i,j)=0; %D->H
                    end
            end
            
            if foldIndex(j)==0
                switch Fmatrix(i,j)
                    case 1 %open
                        if ran2(i,j)<=(1-exp(-kCL(j)*deltaT))
                            Fmatrix(i,j)=2; %unfolded->folded
                        end
                    case 2 %closed
                        if ran2(i,j)<=(1-exp(-kOP(j)*deltaT))
                            Fmatrix(i,j)=1; %folded->unfolded
                        end
                end
            end
        end
        
        %%%consider 'Fmatrix' transition:
        if max(foldIndex)>0
            for f=1:max(foldIndex)
                x=find(f==foldIndex);
                switch Fmatrix(i,x(1))
                    case 1 %open
                        if ran2(i,x(1))<=(1-exp(-kCL(x(1))*deltaT))
                            for j=1:size(x,2)
                                Fmatrix(i,x(j))=2; %unfolded->folded
                            end
                        end
                    case 2 %closed
                        if ran2(i,x(1))<=(1-exp(-kOP(x(1))*deltaT))
                            for j=1:size(x,2)
                                Fmatrix(i,x(j))=1; %folded->unfolded
                            end
                        end
                end
            end
        end
    end
    disp(time)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%POST-SIMULATION PART:

%%%consider N-term two residues and Prolines
for i=1:N
    if i<3 || proSeq(i+START-1)=='P'
        HDmatrix(:,i)=zeros(M,1); %all must be H
    end
end

%%%get simulation result distribution
Distr=zeros(1,maxD+1);
deltaMass=zeros(1,M);
for i=1:M
    deltaMass(i)=sum(HDmatrix(i,:));
    Distr(deltaMass(i)+1)=Distr(deltaMass(i)+1)+1; %Distr(x) means x-1 units of mass above monoisotopic
end
Distr=Distr/sum(Distr); %normalization


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











