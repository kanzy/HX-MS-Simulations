%%2008-11-13 msim2.m: 

clear

M=100; %number of molecules
N=100;  %number of exchangable sites of each molecule

k_ch=300;
k_op=0.1;
k_cl=10;

kmax=k_ch;
deltaT=0.01/kmax;   %step-size of time

tf=2; %folding time

S=zeros(M,N);   %status matrix

for t=0:deltaT:tf
    for i=1:M
        for j=1:N
            switch S(i,j)
                case 0  %status of H_cl
                    if 1-exp(-k_op*deltaT)>=rand
                        S(i,j)=1;
                    end
                case 1  %status of H_op
                    if 1-exp(-k_ch*deltaT)>=rand
                        S(i,j)=2;
                    else if 1-exp(-k_cl*deltaT)>=rand
                            S(i,j)=0;
                        end
                    end
                case 2  %status of D
                    %% do nothing
            end
        end
    end
    disp(t) %display the current time
end

Distr=zeros(1,N+1);
for i=1:M
    deltaMass(i)=0;
    for j=1:N
        if S(i,j)==2
            deltaMass(i)=deltaMass(i)+1;
        end
    end
    Distr(deltaMass(i)+1)=Distr(deltaMass(i)+1)+1;
end





        
            