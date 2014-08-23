%%%2009-04-21 msim4.m: to investigate 2009Apr19 problem as well as DXMS
%%%software problem

clear

% p=0.011;    %richness of C13
% M=100000; %number of molecules tested
% N=104;  %number of C atom in each molecule (104 for Peptide14-38 of a-Syn)
% 
% for i=1:M
%     mass(i)=0;
%     ran=rand(1,N);  %generate a random number between 0~1
%     for j=1:N
%         if ran(j)<p
%             mass(i)=mass(i)+1;
%         end
%     end
% end
% 
% for i=1:N+1
%     deltaC(i)=0;
%     for j=1:M
%         if mass(j)==i-1
%             deltaC(i)=deltaC(i)+1;  %Notice! deltaC(i) means "plus i-1 units of mass"
%         end
%     end
%     if i>4
%         if deltaC(i)==0 && deltaC(i-1)==0 && deltaC(i-2)==0 && deltaC(i-3)~=0
%             k=i;
%         end
%     end
% end
% deltaC=deltaC(1:k)/M;
% 
% stem(0:(k-1),deltaC) %stem plot of distribution of allH molecule (assume all-C12 molecule weight is 0)


%%below result was calculated out based on M=100,000 & N=104 on 2009-04-21, which can be used as a standard MS peaks distribution caused by C13
%%for a molecule with 104 C atoms --- i.e., Peptide14-38 of a-Syn.
deltaC=[3.1695e-001	3.6614e-001	2.1006e-001	7.7810e-002	2.2870e-002	5.0600e-003	9.3000e-004	1.5000e-004	3.0000e-005];

%%to plot this allH mass distribution:
% sizer=size(deltaC);
% stem(0:sizer(2)-1, deltaC)
% hold on



disp('using kc(/min) of a-Syn Peptide_14-38 at 5"C pD4 H->D, generated by FBMME_HD_aSyn_WT.xls (first two residues kc=0): ')
kc=[0	0	0.07946	0.09268	0.05884	0.05884	0.20721	0.03137	0.02953	0.01725	0.03449	0.10156	0.12306	0.09268	0.20721	0.16848	0.05884	0.10156	0.02671	0.02953	0.01725	0.12199	0.29261	0.12306	0.02459
];

kmax=max(kc);
deltaT=0.01/kmax;   %step-size of time

tf=input('enter HX time(min):');

M=5000;    %number of test molecules
N=25;  %number of exchangable sites of each molecule (25 for Peptide14-38 of a-Syn)

D=zeros(M,N);   %H/D status matrix


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%case I: simple intrinsic(totally unprotected) HX

% for t=0:deltaT:tf
%     for i=1:M    
%         ran=rand(1,N);
%         for j=1:N
%             if 1-exp(-kc(j)*deltaT)>=ran(j)
%                 D(i,j)=1;
%             end
%         end
%     end
%     disp(t) %display the current time
% end
% 
% Distr=zeros(1,N+1);
% for i=1:M
%     deltaMass(i)=0;
%     for j=1:N
%         if D(i,j)==1
%             deltaMass(i)=deltaMass(i)+1;
%         end
%     end
%     Distr(deltaMass(i)+1)=Distr(deltaMass(i)+1)+1; %Notice! Distr(x) means "with additional x-1 units of mass compared with allH"
% end
% 
% stem(0:N,Distr,'r')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%case II: one population, one protection factor

% pf=input('Case II -- enter the Protction Factor: '); 
% 
% for t=0:deltaT:tf
%     for i=1:M    
%         ran=rand(1,N);
%         for j=1:N
%             if 1-exp(-kc(j)*deltaT/pf)>=ran(j)
%                 D(i,j)=1;
%             end
%         end
%     end
%     disp(t) %display the current time
% end
% 
% Distr=zeros(1,N+1);
% for i=1:M
%     deltaMass(i)=0;
%     for j=1:N
%         if D(i,j)==1
%             deltaMass(i)=deltaMass(i)+1;
%         end
%     end
%     Distr(deltaMass(i)+1)=Distr(deltaMass(i)+1)+1; %Notice! Distr(x) means "with additional x-1 units of mass compared with allH"
% end
% 
% stem(0:N,Distr,'r')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%case III: two populations, one same protected, another unprotected

% sub=input('Case III -- enter the percentage of the same protected subpopulation(e.g. 0.5 for half): ');
% pf=input('Case III -- enter Protction Factor for the protected subpopulation: '); 
% 
% M1=sub*M;   %the same protected subpopulation
% 
% for t=0:deltaT:tf
%     for i=1:M1    
%         ran=rand(1,N);  %there was little difference if change ran(j) to rand, but run faster.
%         for j=1:N
%             if 1-exp(-kc(j)*deltaT/pf)>=ran(j)  
%                 D(i,j)=1;
%             end
%         end
%     end
%     for i=(M1+1):M    
%         ran=rand(1,N);
%         for j=1:N
%             if 1-exp(-kc(j)*deltaT)>=ran(j)
%                 D(i,j)=1;
%             end
%         end
%     end
%     disp(t) %display the current time
% end
% 
% Distr=zeros(1,N+1);
% for i=1:M
%     deltaMass(i)=0;
%     for j=1:N
%         if D(i,j)==1
%             deltaMass(i)=deltaMass(i)+1;
%         end
%     end
%     Distr(deltaMass(i)+1)=Distr(deltaMass(i)+1)+1; %Notice! Distr(x) means "with additional x-1 units of mass compared with allH"
% end
% 
% stem(0:N,Distr,'r')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%case IV: one population, two diff protected regions

reg1=input('Case IV -- enter the length of the first region(between 3 and 24, 12 for half): ');
pf1=input('Case IV -- enter the Protction Factor of the first region(1 for unprotected): '); 
disp('the remaining will be the second region.');
pf2=input('Case IV -- enter the Protction Factor of the second region(1 for unprotected): '); 

for t=0:deltaT:tf
    for i=1:M    
        ran=rand(1,N);  %there was little difference if change ran(j) to rand, but run faster.
        for j=1:reg1
            if 1-exp(-kc(j)*deltaT/pf1)>=ran(j)
                D(i,j)=1;
            end
        end
        for j=(reg1+1):N
            if 1-exp(-kc(j)*deltaT/pf2)>=ran(j)
                D(i,j)=1;
            end
        end
    end
    disp(t) %display the current time
end

Distr=zeros(1,N+1);
for i=1:M
    deltaMass(i)=0;
    for j=1:N
        if D(i,j)==1
            deltaMass(i)=deltaMass(i)+1;
        end
    end
    Distr(deltaMass(i)+1)=Distr(deltaMass(i)+1)+1; %Notice! Distr(x) means "with additional x-1 units of mass compared with allH"
end

stem(0:N,Distr,'r')