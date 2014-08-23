%%%2009-12-20 phxsim.m: simulate pulse-labeling HX-MS experiment with
%%%protein folding models.

global SFfitN
global foldingModel
global k_fm

disp('This program need run phxsim_preload.m first, or import previously saved workspace.')
flag=input('Input: (1=run preload; 2=import workspace; 0=already run or imported)');
switch flag
    case 1
        phxsim_preload
    case 2
        uiimport
    case 0
        %do nothing
    otherwise
        error('Wrong input!')
end

% START=input('Input the START residue number of the simulating peptide: ');
% END=input('Input the END residue number of the simulating peptide: ');
START=77;
END=90;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%start simulation

M=10000; %number of simulating molecules
N=END-START+1;

HDmatrix=ones(M,N); %0=H; 1=D

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%simulation I. Folding phase:
flagSim=1; %to determine whether change 'FoldStates' in phxsim_ppoe1.m, phxsim_iup1.m,...
sampleType=input('Input the sample type(1=real pulse-labeling sample; 2=U-control; 3=N-control): ');
switch sampleType
    case 1
        FoldStates=zeros(M,1); %0=U; 1=I; 2=Ix; 3=N
        kcDH=kcDH_f;
        kcHD=kcHD_f;
    case 2
        FoldStates=zeros(M,1); %0=U; 1=I; 2=Ix; 3=N
        kcDH=kcDH_q; %approximation: unfolding pH = quench pH
        kcHD=kcHD_q;
        flagSim=0; %no change of 'FoldStates' for U-ctrl during folding phase
    case 3
        FoldStates=3*ones(M,1); %0=U; 1=I; 2=Ix; 3=N
        kcDH=kcDH_f;
        kcHD=kcHD_f;
    otherwise
        error('Wrong input of sampleType!')
end

Dfraction=QFratio(1)/(QFratio(1)+QFratio(2));
Hfraction=QFratio(2)/(QFratio(1)+QFratio(2));
simTime=Time_f;

phxsim_sim %call phxsim_sim.m
flagSim=1;

figure
stem(0:N,Distr,'r')
xlabel('delta Mass(above monoisotopic)')
ylabel('Molecules Count')
title({'H/D Distribution at the end of Folding Phase';'(without removing D at N-term & Proline)'})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%simulation II. Pulse-labeling phase:
kcDH=kcDH_p;
kcHD=kcHD_p;
Dfraction=QFratio(1)/(QFratio(1)+QFratio(2)+QFratio(3));
Hfraction=(QFratio(2)+QFratio(3))/(QFratio(1)+QFratio(2)+QFratio(3));
simTime=Time_p;

phxsim_sim %call phxsim_sim.m

figure
stem(0:N,Distr,'r')
xlabel('delta Mass(above monoisotopic)')
ylabel('Molecules Count')
title({'H/D Distribution at the end of Pulse-labeling Phase';'(without removing D at N-term & Proline)'})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%simulation III. Quench until MS detection phase:
clear flag
flag=input('Simulate the phase after pulse-labeling (from quench until MS detection phase)? (1=yes[may take a long time], 0=skip)');
if flag==1
    kcDH=kcDH_q;
    kcHD=kcHD_q;
    Dfraction=0; %approximation
    Hfraction=1;
    simTime=Time_q;
    
    phxsim_sim %call phxsim_sim.m
    
    figure
    stem(0:N,Distr,'r')
    xlabel('delta Mass(above monoisotopic)')
    ylabel('Molecules Count')
    title({'H/D Distribution detected by MS';'(without removing D at N-term & Proline)'})
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%peptide correction: N-term two residues and all Prolines must not be D
maxD=N;
for i=1:N
    if i<3 || currSeq(i+START-1)=='P'
        HDmatrix(:,i)=zeros(M,1); %all must be H
        maxD=maxD-1; %calculate the exchangable hydrogen number of this peptide
    end
end

finalDistr=zeros(1,maxD+1);
deltaMass=zeros(1,M);
for i=1:M
    deltaMass(i)=sum(HDmatrix(i,:));
    finalDistr(deltaMass(i)+1)=finalDistr(deltaMass(i)+1)+1; %finalDistr(x) means x-1 units of mass above monoisotopic
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Convolution with C13 distribution & Components analysis:
clear DistrA carbonDistr obsDistr

DistrA=finalDistr/sum(finalDistr); %normalization

%%%calulate 'carbonDistr':
[peptideMass, peptideCarbonNum]=pepinfo(currSeq(START:END)); %call pepinfo.m
deltaC=binopdf(0:peptideCarbonNum,peptideCarbonNum,0.0111);
for i=3:(N+1)
    if deltaC(i)<1e-3 && deltaC(i-1)<1e-3 && deltaC(i-2)>=1e-3
        k=i-2; break
    end
end
maxC=k-1; %the observable C13 isotope peaks number of this peptide
carbonDistr=deltaC(1:(1+maxC));

obsDistr=conv(carbonDistr, DistrA);
obsDistr=obsDistr/sum(obsDistr); %normalization

%%%plot above results:
figure
stem(0:maxC, carbonDistr,'fill','k'); hold on
stem(0:(maxC+maxD), obsDistr,'b'); hold on
stem(0:maxD, DistrA,'r'); hold on
v=axis; 
axis([-1 maxC+maxD+1 v(3) v(4)])
ylabel('Normalized Intensity')
xlabel('delta Mass (above monoisotopic)')
title([proteinName,' Peptide ',num2str(START),'--',num2str(END)],'FontWeight','bold')










