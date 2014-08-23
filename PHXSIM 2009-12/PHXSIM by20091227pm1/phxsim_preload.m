%%%2009-12-26 phxsim_preload.m: preload variables for phxsim.m use

snase='ATSTKKLHKEPATLIKAIDGDTVKLMYKGQPMTFRLLLVDTPETKHPKKGVEKYGPEASAFTKKMVENAKKIEVEFDKGQRTDKYGRGLAYIYADGKMVNEALVRQGLAKVAYVYKGNNTHEQLLRKSEAQAKKEKLNIWSEDNADSGQ';

proteinName=input('Input the protein name(e.g. snase): ','s');
switch proteinName  %may add more protein cases here
    case 'snase'
        currSeq=snase;
    otherwise
        error('Unknown protein!');
end

%%%quench-flow experiment info:
disp('Now program will collect information of the quench-flow pulse-labeling experiment: ')
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
Time_q=300;
Temp_QF=15;
Temp_MS=4;

%%%calculate intrinsic HX rates:
kcDH_f = fbmme_dh(currSeq, pH_f, Temp_QF, 'poly');
kcHD_f = fbmme_hd(currSeq, pH_f, Temp_QF, 'poly', 0);
kcDH_p = fbmme_dh(currSeq, pH_p, Temp_QF, 'poly');
kcHD_p = fbmme_hd(currSeq, pH_p, Temp_QF, 'poly', 0);
kcDH_q = fbmme_dh(currSeq, pH_q, Temp_MS, 'poly');
kcHD_q = fbmme_hd(currSeq, pH_q, Temp_MS, 'poly', 0);


foldingModel=input('Input the folding model to use: (e.g. iup1 or ppoe1)', 's');


%%%pre-allocate protection factors:
switch foldingModel
    case 'iup1'
        pf_U=ones(size(currSeq));
        pf_N=ones(size(currSeq));
        pf_I=ones(size(currSeq));
        
    case 'ppoe1'
        pf_U=ones(size(currSeq));
        pf_N=ones(size(currSeq));
        pf_I=ones(size(currSeq));
        pf_Ix=ones(size(currSeq));
        
    otherwise
        error('Unknown folding model!')
end


%%%get kinetics rates of the folding model:
flagSF=input('Whether use stop-flow folding (fluorescence) experiment data to fit the kinetics rates?(1=yes,0=no)');
switch flagSF
    case 0
        %%%assign kinetics rates:
        switch foldingModel
            case 'iup1'
                k_fm=input('Input k_fm as [kUN, kNU, kUI, kIU, kIN, kNI]: ');
                kUN=k_fm(1); kNU=k_fm(2);
                kUI=k_fm(3); kIU=k_fm(4);
                kIN=k_fm(5); kNI=k_fm(6);
                
            case 'ppoe1'
                k_fm=input('Input k_fm as [kUI, kIU, kIIx, kIxI, kIN, kNI]: ');
                kUI=k_fm(1); kIU=k_fm(2);
                kIIx=k_fm(3); kIxI=k_fm(4);
                kIN=k_fm(5); kNI=k_fm(6);
                
            otherwise
                error('Unknown folding model!')
        end
        
    case 1
        disp('Get kinetics rates of the folding model by fitting stop-flow fluor data...')
        
        % SFfitLA=input('Input SFfitLA as [{lamda, Amp}pairs]: ');
        SFfitLA=[47, -0.1; 9.7, 0.63; 2.4, 0.17; 0.3, 0.12]; %{lamda,Amp} fitting result of stop-flow fluor data of SNase refolding at zero denaturant, pH5.3, 15'C by Sabrina (2008 PNAS paper)
        
        SFfitN(:,1)=0:0.002:10; %time points
        SFfitN(:,2) = 1 - (SFfitLA(1,2)*exp(-SFfitLA(1,1)*SFfitN(:,1)) + SFfitLA(2,2)*exp(-SFfitLA(2,1)*SFfitN(:,1)) ...
            + SFfitLA(3,2)*exp(-SFfitLA(3,1)*SFfitN(:,1)) + SFfitLA(4,2)*exp(-SFfitLA(4,1)*SFfitN(:,1)))...
            /(sum(SFfitLA(:,2))); %fraction of N at each time point
        
        switch foldingModel
            case 'iup1'
                %         k_fm_ini=input('Input the initial values for IUP1 Model fitting as [kUN,kNU,kUI,kIU,kIN,kNI]: ');
                k_fm_ini=[1, 0.1, 15, 1, 500, 1];
                lb=zeros(1,6);
                ub=Inf*ones(1,6);
                k_fm=lsqnonlin(@phxsim_kfit,k_fm_ini,lb,ub);
                                
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
                
                %%%assign above folding/unfolding rates for simulation:
                kUI=k_fm(1); kIU=k_fm(2);
                kIIx=k_fm(3); kIxI=k_fm(4);
                kIN=k_fm(5); kNI=k_fm(6);
                
            otherwise
                error('Unknown folding model!')
        end
        
    otherwise
        error('Wrong input of flagSF!')
end







