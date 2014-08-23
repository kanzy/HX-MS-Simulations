%%%2009-12-20 phxsim.m: simulate pulse-labeling HX-MS experiment with
%%%protein folding models.

global SFfitN
global foldingModel
global k_fm
clear SaveResults

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
phxsim_preload_plot %call phxsim_preload_plot.m


% START=input('Input the START residue number of the simulating peptide: ');
% END=input('Input the END residue number of the simulating peptide: ');
START=77;
END=90;


pf_N(77:90)=[1, 1e1, 1e1, 1, 1, 1, 1e1, 1, 1e2, 1e2, 1e2, 1e2, 1e2, 1e2];
pf_I(77:90)=[1, 1e1,   1, 1, 1, 1,   1, 1, 1e1, 1e1, 1e1, 1e1, 1e1, 1e1];
pf_Ix=pf_I; %assume I & Ix have the same protection factors



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%start simulation

% sampleType=input('Input the sample type(1=real pulse-labeling sample; 2=U-control; 3=N-control): ');
for sampleType=1:3
    
    M=10000; %number of simulating molecules
    N=END-START+1;
    HDmatrix=ones(M,N); %0=H; 1=D
    flagSim=1; %to determine whether change 'FoldStates' in phxsim_ppoe1.m, phxsim_iup1.m,...
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%simulation I. Folding phase:
    
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
    
    % figure
    % stem(0:N,Distr,'r')
    % xlabel('delta Mass(above monoisotopic)')
    % ylabel('Molecules Count')
    % title({'D-distribution at the end of Folding Phase';'(without removing D at N-term & Proline)'})
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%simulation II. Pulse-labeling phase:
    kcDH=kcDH_p;
    kcHD=kcHD_p;
    Dfraction=QFratio(1)/(QFratio(1)+QFratio(2)+QFratio(3));
    Hfraction=(QFratio(2)+QFratio(3))/(QFratio(1)+QFratio(2)+QFratio(3));
    simTime=Time_p;
    
    phxsim_sim %call phxsim_sim.m
    
    % figure
    % stem(0:N,Distr,'r')
    % xlabel('delta Mass(above monoisotopic)')
    % ylabel('Molecules Count')
    % title({'D-distribution at the end of Pulse-labeling Phase';'(without removing D at N-term & Proline)'})
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%simulation III. Quench until MS detection phase:
    % clear flag
    % flag=input('Simulate the phase after pulse-labeling (from quench until MS detection phase)? (1=yes[may take a long time], 0=skip)');
    flag=0; %skip this simulation: usually very little additional D loss occurs here 
    if flag==1
        kcDH=kcDH_q;
        kcHD=kcHD_q;
        Dfraction=0; %approximation
        Hfraction=1;
        simTime=Time_q;
        
        phxsim_sim %call phxsim_sim.m
        
        %     figure
        %     stem(0:N,Distr,'r')
        %     xlabel('delta Mass(above monoisotopic)')
        %     ylabel('Molecules Count')
        %     title({'D-distribution detected by MS';'(without removing D at N-term & Proline)'})
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%Convolution with C13 distribution:
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
    % figure
    % stem(0:maxC, carbonDistr,'fill','k'); hold on
    % stem(0:(maxC+maxD), obsDistr,'b'); hold on
    % stem(0:maxD, DistrA,'r'); hold on
    % v=axis;
    % axis([-1 maxC+maxD+1 v(3) v(4)])
    % ylabel('Normalized Intensity')
    % xlabel('delta Mass (above monoisotopic)')
    % title({[proteinName,'_',num2str(START),'-',num2str(END)];'Simulated Final D-distribution(red) & Convoluted with C13(blue)'})
    
    %%%save above results:
    SaveResults.DistrA{sampleType}=DistrA;
    SaveResults.obsDistr{sampleType}=obsDistr;
    SaveResults.carbonDistr=carbonDistr;
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%components fitting & plot: refer to msfit6f.m METHOD C(2009-11-09)
DistrA_x=SaveResults.DistrA{1};
DistrA_U=SaveResults.DistrA{2};
DistrA_N=SaveResults.DistrA{3};

weights=max(DistrA_U, DistrA_N).^3;

ini_fraction_U=0; ini_fraction_N=0;
ini_fractions=[ini_fraction_U; ini_fraction_N];

options=optimset('TolFun', 1e-36, 'TolX', 1e-30);  %seems no difference with different setting here
fitX=lsqnonlin(@phxsim_cafit, ini_fractions, [0;0], [1;1], options, DistrA_U, DistrA_N, DistrA_x, weights);

fraction_U=fitX(1);
fraction_N=fitX(2);
fraction_I=1-(fraction_U + fraction_N);
Distr_I=DistrA_x - (fraction_U*DistrA_U + fraction_N*DistrA_N);

%%%plot above results:
h=figure;

subplot(3,1,1)
stem(0:maxD, DistrA_U,'m');
hold on
title([proteinName,' Peptide ',num2str(START),'--',num2str(END),'  Components Fitting by phxsim.m'],'FontWeight','bold')
axis([-1, maxD+1, 0, max(DistrA_U)+0.02])
text(maxD*0.9, (max(DistrA_U)+0.02)*0.85, 'U-control','BackgroundColor',[.7 .9 .7],'FontAngle','italic');

subplot(3,1,2)
stem(0:maxD, DistrA_x,'k');
hold on
plot(0:maxD, fraction_U*DistrA_U,'m','LineWidth',1.5);
hold on
plot(0:maxD, fraction_N*DistrA_N,'b','LineWidth',1.5);
hold on
plot(0:maxD, Distr_I,'g','LineWidth',1.5);
hold on
ylabel('Normalized Intensity')
axis([-1, maxD+1, min(0, min(Distr_I)), max(DistrA_x)+0.02])
text(maxD*0.9, (max(DistrA_x)+0.02)*0.85, [num2str(Time_f*1000, '%10.1f'), ' msec folding'],'BackgroundColor',[.7 .9 .7],'FontAngle','italic');
text(maxD*0.85,(max(DistrA_x)+0.02)*0.65, ['U=', num2str(fraction_U*100, '%5.1f'),'%'],'FontWeight','bold');
text(maxD*0.85,(max(DistrA_x)+0.02)*0.5, ['I=', num2str(fraction_I*100, '%5.1f'),'%'],'FontWeight','bold');
text(maxD*0.85,(max(DistrA_x)+0.02)*0.35, ['N=', num2str(fraction_N*100, '%5.1f'),'%'],'FontWeight','bold');

subplot(3,1,3)
stem(0:maxD, DistrA_N,'b');
hold on
xlabel('delta Mass (above monoisotopic)')
axis([-1, maxD+1, 0, max(DistrA_N)+0.02])
text(maxD*0.9, (max(DistrA_N)+0.02)*0.85, 'N-control','BackgroundColor',[.7 .9 .7],'FontAngle','italic');

%%%save above figure:
SaveFigureName=[proteinName '_' num2str(START) '-' num2str(END) '_phxsim_ComponentsFitting.fig'];
saveas(figure(h),SaveFigureName)

%%%save above results:
SaveResults.fraction_U=fraction_U;
SaveResults.fraction_N=fraction_N;
SaveResults.fraction_I=fraction_I;
SaveResults.Distr_I=Distr_I;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%final save
SaveResults.proteinName=proteinName;
SaveResults.QFratio=QFratio;
SaveResults.pH=[pH_f, pH_p, pH_q];
SaveResults.Time=[Time_f, Time_p, Time_q];
SaveResults.Temp=[Temp_QF, Temp_MS];
SaveResults.kcDH=[kcDH_f, kcDH_p, kcDH_q];
SaveResults.kcHD=[kcHD_f, kcHD_p, kcHD_q];
SaveResults.foldingModel=foldingModel;
SaveResults.k_fm=k_fm;
if flagSF==1
    SaveResults.SFfitLA=SFfitLA;
end
SaveResults.START=START;
SaveResults.END=END;
SaveResults.maxC=maxC;
SaveResults.maxD=maxD;
switch foldingModel
    case 'iup1'
        SaveResults.pf=[pf_U; pf_N; pf_I];
    case 'ppoe1'
        SaveResults.pf=[pf_U; pf_N; pf_I; pf_Ix];
end

SaveFileName=[proteinName '_' num2str(START) '-' num2str(END) '_phxsim_SaveResults.mat'];
save(SaveFileName,'SaveResults')
disp([SaveFileName, ' has been saved in MATLAB current working directionary!'])














