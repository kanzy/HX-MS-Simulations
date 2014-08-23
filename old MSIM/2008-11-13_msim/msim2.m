%%2008-11-13 msim2.m: to simulate Mass Spec of protein/peptide with
%%HX(H->D) exactly as Konermann 2008 

clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Part 1: first calculate allH peaks distribution due to natural C13 richness

p=0.011;    %richness of C13
M=30000; %number of molecules tested
N=744;  %number of C atom in each molecule (744 for SNase)

for i=1:M
    mass(i)=0;
    ran=rand(1,N);  %generate a random number between 0~1
    for j=1:N
        if ran(j)<p
            mass(i)=mass(i)+1;
        end
    end
    disp(i) %display the current molecule#
end

for i=1:N+1
    deltaC(i)=0;
    for j=1:M
        if mass(j)==i-1
            deltaC(i)=deltaC(i)+1;  %Notice! deltaC(i) means "with additional i-1 units of mass compared with allC-12"
        end
    end
    if i>4
        if deltaC(i)==0 && deltaC(i-1)==0 && deltaC(i-2)==0 && deltaC(i-3)~=0
            k=i;
        end
    end
end
deltaC=deltaC(1:k)/M;
% stem(0:(k-1),deltaC) %stem plot of distribution of allH molecule (assume all-C12 molecule weight is 0)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Part 2: to simulate a simple two-state folding protein native HX, start from
%%%allH

M=300; %number of molecules
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
stem(0:sizer(2)-1, finalDistr)



