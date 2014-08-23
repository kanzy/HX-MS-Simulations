%%%2009-12-20 phxsim.m: simulate pulse-labeling HX-MS experiment with
%%%protein folding models.

snase='ATSTKKLHKEPATLIKAIDGDTVKLMYKGQPMTFRLLLVDTPETKHPKKGVEKYGPEASAFTKKMVENAKKIEVEFDKGQRTDKYGRGLAYIYADGKMVNEALVRQGLAKVAYVYKGNNTHEQLLRKSEAQAKKEKLNIWSEDNADSGQ';

% proteinName=input('Input the protein name(snase): ','s');
proteinName='snase';
switch proteinName  %may add more protein cases here
    case 'snase'
        currSeq=snase;
    otherwise
        error('Unknown protein!');
end

%%%quench-flow experiment info:
% QFratio=input('Input the volumes ratio of four buffers [U F P Q]: (e.g. [1 10 2 2])');
% pH_f=input('Input the pH of the folding phase: ');
% pH_p=input('Input the pH of the pulse-labeling phase: ');
% pH_q=input('Input the pH after quench: (should be ~2.5)');
% Time_f=input('Input the time(in sec) of the folding phase: ');
% Time_p=input('Input the time(in sec) of the pulse-labeling phase: ');
% Time_q=input('Input estimated time(in sec) after quench until MS detection: ');
% Temp_QF=input('Input the temperature(Centigrade) in the quench-flow: ');
% Temp_MS=input('Input estimated average temperature(Centigrade) after quench until MS detection: (maybe ~4)');
QFratio=[1 10 2 2];
pH_f=5.3;
pH_p=9.2;
pH_q=2.5;
Time_f=0.235;
Time_p=0.015;
Time_q=500;
Temp_QF=15;
Temp_MS=4;

%%%calculate intrinsic HX rates:
kcDH_f = fbmme_dh(currSeq, pH_f, Temp_QF, 'poly');
kcHD_f = fbmme_hd(currSeq, pH_f, Temp_QF, 'poly', 0);
kcDH_p = fbmme_dh(currSeq, pH_p, Temp_QF, 'poly');
kcHD_p = fbmme_hd(currSeq, pH_p, Temp_QF, 'poly', 0);
kcDH_q = fbmme_dh(currSeq, pH_q, Temp_MS, 'poly');
kcHD_q = fbmme_hd(currSeq, pH_q, Temp_MS, 'poly', 0);

% START=input('Input the START residue number of the simulating peptide: ');
% END=input('Input the END residue number of the simulating peptide: ');
START=77;
END=90;

pf_U=ones(size(currSeq));
pf_N=ones(size(currSeq));
pf_I=ones(size(currSeq));

pf_N(START:END)=[1, 1e3, 1e3, 1, 1, 1, 1e3, 1, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3];
pf_I(START:END)=[1,   1,   1, 1, 1, 1,   1, 1, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%get folding/unfolding rates of a folding model by fitting stop-flow fluor data:
global SFfitN
global foldingModel
global k_fm

SFfitLA=[47, -0.1; 9.7, 0.63; 2.4, 0.17; 0.3, 0.12]; %{lamda,Amp} fitting result of stop-flow fluor data of SNase refolding at zero denaturant, pH5.3, 15'C by Sabrina (2008 PNAS paper)
SFfitN(:,1)=0:0.002:10; %time points
SFfitN(:,2) = 1 - (SFfitLA(1,2)*exp(-SFfitLA(1,1)*SFfitN(:,1)) + SFfitLA(2,2)*exp(-SFfitLA(2,1)*SFfitN(:,1)) ...
    + SFfitLA(3,2)*exp(-SFfitLA(3,1)*SFfitN(:,1)) + SFfitLA(4,2)*exp(-SFfitLA(4,1)*SFfitN(:,1)))...
    /(sum(SFfitLA(:,2))); %fraction of N at each time point

foldingModel=input('Input the folding model to use: (e.g. iup1 or ppoe2)', 's');
switch foldingModel
    case 'iup1'
        %         k_fm_ini=input('Input the initial values for IUP1 Model fitting as [kUN,kNU,kUI,kIU,kIN,kNI]: ');
        k_fm_ini=[1, 0.1, 15, 1, 500, 1];
        lb=zeros(1,6);
        ub=Inf*ones(1,6);
        k_fm=lsqnonlin(@phxsim_kfit,k_fm_ini,lb,ub);
        
        %%%plot fitted results with original experimental curve:
        figure
        semilogx(SFfitN(:,1), SFfitN(:,2),'k','LineWidth',2); hold on %plot stop-flow fluor curve
        [t,y] = ode15s(@fm_iup1,[SFfitN(1,1) SFfitN(end,1)],[1 0 0]);
        
        semilogx(t,y(:,1),'r'); hold on %plot U
        semilogx(t,y(:,2),'g'); hold on %plot I
        semilogx(t,y(:,3),'b'); hold on %plot N
        v=axis;
        axis([1e-4,v(2),-0.01,1.01]);
        xlabel('Time (sec)')
        ylabel('Fraction')
        title('IUP1 Model Fitting (red=U; green=I; blue=N; black=exp fluor)')
        
        %%%assign above folding/unfolding rates for simulation:
        kUN=k_fm(1); kNU=k_fm(2);
        kUI=k_fm(3); kIU=k_fm(4);
        kIN=k_fm(5); kNI=k_fm(6);
        
    case 'ppoe1'
        %         k_fm_ini=input('Input the initial values for IUP1 Model fitting as [kUI,kIU,kIIx,kIxI,kIN,kNI]: ');
        k_fm_ini=[15, 1, 200, 5, 500, 1];
        lb=zeros(1,6);
        ub=Inf*ones(1,6);
        k_fm=lsqnonlin(@phxsim_kfit,k_fm_ini,lb,ub);
        
        %%%plot fitted results with original experimental curve:
        figure
        semilogx(SFfitN(:,1), SFfitN(:,2),'k','LineWidth',2); hold on %plot stop-flow fluor curve
        [t,y] = ode15s(@fm_ppoe1,[SFfitN(1,1) SFfitN(end,1)],[1 0 0 0]);
        semilogx(t,y(:,1),'r'); hold on %plot U
        semilogx(t,y(:,2),'g'); hold on %plot I
        semilogx(t,y(:,3),'y'); hold on %plot Ix
        semilogx(t,y(:,4),'b'); hold on %plot N
        v=axis;
        axis([1e-4,v(2),-0.01,1.01]);
        xlabel('Time (sec)')
        ylabel('Fraction')
        title('PPOE1 Model Fitting (red=U; green=I; yellow=Ix; blue=N; black=exp fluor)')
        
        %%%assign above folding/unfolding rates for simulation:
        kUI=k_fm(1); kIU=k_fm(2);
        kIIx=k_fm(3); kIxI=k_fm(4);
        kIN=k_fm(5); kNI=k_fm(6);
        
        
    otherwise
        error('Unknown folding model!')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%start simulation:

M=1000; %number of simulating molecules

HDmatrix=ones(M,END-START+1); %0=H;1=D
FoldStates=zeros(M,1); %0=U;1=I;2=N

kmax=max([max(kcDH_p), max(kcHD_p), max(k_fm)]);
deltaT=0.01/kmax;   %step size of simulation time







