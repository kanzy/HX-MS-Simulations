%%2008-11-07 msim.m: to simulate Mass Spec of protein/peptide with HX(H->D)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Part 1: first calculate allH peaks distribution due to natural C13 richness

p=0.011;    %richness of C13
M=20000; %number of molecules tested
N=744;  %number of C atom in each molecule (744 for SNase)

for i=1:M
    mass(i)=0;
    ran=rand(1,N);  %generate a random number between 0~1
    for j=1:N
        if ran(j)<p
            mass(i)=mass(i)+1;
        end
    end
end

for i=1:N
    deltaC(i)=0;
    for j=1:M
        if mass(j)==i
            deltaC(i)=deltaC(i)+1;
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

N=100;  %number of exchangable H atom in the molecule

t=0:0.1:20;  %folding/HX time from 0~20 sec 

[t,y] = ode15s('nhxf',t,[1 0 0]);   %y(:,3) is the labeled D(%) at different folding/HX time

sizer=size(t);
for i=1:sizer(1)
    for j=1:N
        D(i,j)=(factorial(N)/(factorial(j-1)*factorial(N-j+1)))*(y(i,3)^(j-1))*((1-y(i,3))^(N-j+1));    %D is the result matrix of ( folding/HX time vs. D-distribution ) 
    end
    result(i,:)=conv(deltaC, D(i,:));
end

sizer=size(result);

stem(0:sizer(2)-1,result(1,:))
hold on

stem(0:sizer(2)-1,result(30,:),'r')
hold on

stem(0:sizer(2)-1,result(58,:),'c')
hold on

stem(0:sizer(2)-1,result(201,:),'g')
hold on



