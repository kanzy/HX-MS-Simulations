%%%2010-01-04 phxsim2.m: replace Mont Carlo with Diff Equations simulation
%%%2009-12-20 phxsim.m: simulate pulse-labeling HX-MS experiment with
%%%protein folding models.

global SFfitN
global foldingModel
global k_fm
global kex
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


clear fractionsMatrix
for currResidue=START:END

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%simulation I. Folding phase:
    Dfraction=QFratio(1)/(QFratio(1)+QFratio(2));
    Hfraction=QFratio(2)/(QFratio(1)+QFratio(2));
    kcDH=kcDH_f;
    kcHD=kcHD_f;
    phxsim2_kex %call phxsim2_kex.m
    
    [t,y_f] = ode15s(@fmhx_ppoe1,[0 Time_f],[0 0 0 0 1 0 0 0]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%simulation II. Pulse-labeling phase:
    Dfraction=QFratio(1)/(QFratio(1)+QFratio(2)+QFratio(3));
    Hfraction=(QFratio(2)+QFratio(3))/(QFratio(1)+QFratio(2)+QFratio(3));
    kcDH=kcDH_p;
    kcHD=kcHD_p;
    phxsim2_kex %call phxsim2_kex.m
    
    [t,y_p] = ode15s(@fmhx_ppoe1,[0 Time_p],y_f(end,:));
    
    
    
    fractionsMatrix(currResidue-START+1,:)=y_p(end,:);
end
    
M=10000;
N=END-START+1;
HDmatrix=zeros(M,N);
ran=rand(M,N);
for i=1:M
    for j=1:N
        if ran(i,j)<sum(fractionsMatrix(j,5:8))
            HDmatrix(i,j)=1;
        end
    end
end
    
Distr=zeros(1,N+1);
deltaMass=zeros(1,M);
for i=1:M
    deltaMass(i)=sum(HDmatrix(i,:));
    Distr(deltaMass(i)+1)=Distr(deltaMass(i)+1)+1; %Distr(x) means x-1 units of mass above monoisotopic
end

figure
stem(0:N,Distr,'r')
xlabel('delta Mass(above monoisotopic)')
ylabel('Molecules Count')
    














