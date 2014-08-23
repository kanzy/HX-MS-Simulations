%%2008-11-13 msim2.m: to simulate Mass Spec of protein/peptide with
%%HX(H->D) exactly as Konermann 2008 

clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Part 1: first calculate allH peaks distribution due to natural C13 richness

% p=0.011;    %richness of C13
% M=30000; %number of molecules tested
% N=744;  %number of C atom in each molecule (744 for SNase)
% 
% for i=1:M
%     mass(i)=0;
%     ran=rand(1,N);  %generate a random number between 0~1
%     for j=1:N
%         if ran(j)<p
%             mass(i)=mass(i)+1;
%         end
%     end
%     disp(i) %display the current molecule#
% end
% 
% for i=1:N+1
%     deltaC(i)=0;
%     for j=1:M
%         if mass(j)==i-1
%             deltaC(i)=deltaC(i)+1;  %Notice! deltaC(i) means "with additional i-1 units of mass compared with allC-12"
%         end
%     end
%     if i>4
%         if deltaC(i)==0 && deltaC(i-1)==0 && deltaC(i-2)==0 && deltaC(i-3)~=0
%             k=i;
%         end
%     end
% end
% deltaC=deltaC(1:k)/M;
% % stem(0:(k-1),deltaC) %stem plot of distribution of allH molecule (assume all-C12 molecule weight is 0)

%%below result was calculated out based on M=100,000 & N=744 on 2008-11-14
%%12pm, which can be used as a standard MS peaks distribution caused by C13
%%for a molecule with 744 C atoms.
% deltaC=[0.00021	0.0024	0.00888	0.0252	0.05118	0.08577	0.11783	0.137	0.13885	0.12594	0.10096	0.07863	0.05344	0.03437	0.01886	0.01048	0.00555	0.00249	0.00125	0.00049	0.00014	4e-005	3e-005	0	0	1e-005	0 0 0];

%%below was theoretical numbers of C13 distribution from Gaussian fitting
%%based on above deltaC values
%% General model Gauss1:
%%        f(x) =  a1*exp(-((x-b1)/c1)^2)
%% Coefficients (with 95% confidence bounds):
%%        a1 =      0.1402  (0.1361, 0.1443)
%%        b1 =       7.927  (7.83, 8.023)
%%        c1 =       4.021  (3.884, 4.158)
%% Goodness of fit:
%%   SSE: 0.0003535
%%   R-square: 0.9946
%%   Adjusted R-square: 0.9941
%%   RMSE: 0.003687
deltaC=[2.8781e-003	7.2121e-003	1.5970e-002	3.1247e-002	5.4027e-002	8.2544e-002	1.1144e-001	1.3295e-001	1.4015e-001	1.3055e-001	1.0746e-001	7.8165e-002	5.0239e-002	2.8533e-002	1.4320e-002	6.3506e-003	2.4886e-003	8.6177e-004	2.6369e-004	7.1300e-005	1.7035e-005	3.5967e-006	6.7100e-007	1.1062e-007	1.6114e-008	2.0743e-009	2.3595e-010	2.3716e-011	2.1064e-012	1.6532e-013];

sizer=size(deltaC);
stem(0:sizer(2)-1, deltaC)
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Part 2: to simulate a simple two-state folding protein native HX, start from
%%%allH

M=1000; %number of molecules
N=100;  %number of exchangable sites of each molecule (~145 for SNase)

k_ch=300;   %rate of H->D exchange
k_op=2;   %rate of opening(unfolding)
k_cl=3000;    %rate of closing(folding)

kmax=k_ch;
if k_cl>k_ch
    kmax=k_cl;
end

deltaT=0.01/kmax;   %step-size of time

tf=1.4; %folding time

S=zeros(1,M);   %Close/Open status matrix
D=zeros(M,N);   %H/D status matrix

for t=0:deltaT:tf
    for i=1:M
        switch S(i)
            case 0  %status of Close
                if 1-exp(-k_op*deltaT)>=rand
                    S(i)=1;
                end
            case 1  %status of Open
                for j=1:N
                    if 1-exp(-k_ch*deltaT)>=rand
                        D(i,j)=1;
                    end
                end
                if 1-exp(-k_cl*deltaT)>=rand
                    S(i)=0;
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
% stem(0:N,Distr)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Part 3: combine above results to get final MS peak distribution
finalDistr=conv(deltaC, Distr/M);
sizer=size(finalDistr);
stem(0:sizer(2)-1, finalDistr,'r')



