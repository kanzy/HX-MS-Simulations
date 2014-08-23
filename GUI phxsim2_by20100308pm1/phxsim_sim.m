%%%2010-03-08 phxsim_sim.m:

function obsPeaks=phxsim_sim(START,END,Charge,Temp,Time_p,pH_p,QFratio,pf,Seq)

%%%data prep:
[peptideMass, distND, maxND, maxD]=pepinfo(Seq(START:END)); %call pepinfo.m
kcDH = fbmme_dh(Seq, pH_p, Temp, 'poly'); %call fbmme_dh.m
kcHD = fbmme_hd(Seq, pH_p, Temp, 'poly', 0); %call fbmme_hd.m
Dfraction=QFratio(1)/(QFratio(1)+QFratio(2)+QFratio(3));
Hfraction=(QFratio(2)+QFratio(3))/(QFratio(1)+QFratio(2)+QFratio(3));

%%%do simulation:
M=3000; %number of simulating molecules
N=END-START+1;
HDmatrix=ones(M,N); %0=H; 1=D
kmax=max([max(kcDH(START:END)), max(kcHD(START:END))]);
deltaT=0.1/kmax;   %step size of simulation time
for time=0:deltaT:Time_p
    ran1=rand(M,N);
    for i=1:M
        for j=START:END
            switch HDmatrix(i,j-START+1)
                case 0 %H
                    if ran1(i,j-START+1)<=(1-exp(-kcHD(j)*deltaT/pf(j)))*Dfraction
                        HDmatrix(i,j-START+1)=1; %H->D
                    end
                case 1 %D
                    if ran1(i,j-START+1)<=(1-exp(-kcDH(j)*deltaT/pf(j)))*Hfraction
                        HDmatrix(i,j-START+1)=0; %D->H
                    end
            end
        end
    end
    disp(time)
end

%%%consider N-term two residues and Prolines
for i=1:N
    if i<3 || Seq(i+START-1)=='P'
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

%%%do convolution with allH peaks:
obsDistr=conv(distND, Distr);
obsDistr=obsDistr/sum(obsDistr); %normalization

%%%get MS observable peaks:
clear obsPeaks
DM=1.00628; %delta mass between a deuteron(2.014102 u) and a proton(1.007825 u)
obsPeaks(1,1)=peptideMass/Charge+1.007276; %m/z of mono; 1.007276 is the mass of proton
obsPeaks(1,2)=obsDistr(1);
for i=2:(maxD+maxND+1)
    obsPeaks(i,1)=obsPeaks(i-1,1)+DM/Charge;
    obsPeaks(i,2)=obsDistr(i);
end



