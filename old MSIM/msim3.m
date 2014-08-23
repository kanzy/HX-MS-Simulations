%%%2009-04-24 msim3.m: modified from msim.m, using diff-equations 

clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Part 1: first calculate allH peaks distribution due to natural C13 richness:

% p=0.011;    %richness of C13
% M=100000; %number of molecules tested
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


%below result was calculated out based on the avarage of {M=100,000 & N=744 on 2009-04-24} & {M=100,000 & N=744 on 2008-11-14}, which can be used as a standard MS peaks distribution caused by C13
%for a molecule with 744 C atoms --- i.e., SNase_P117G_H124L.
deltaC=[0.000225	0.002365	0.008655	0.024855	0.05162	0.085585	0.11752	0.13763	0.13893	0.1265	0.10228	0.07759	0.052695	0.03386	0.018935	0.010725	0.00534	0.002605	0.001285	0.000515	0.00017	7e-005	4e-005
];

%to plot this allH mass distribution:
sizer=size(deltaC);
stem(0:sizer(2)-1, deltaC)
hold on



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Part 2: to simulate a simple two-state folding protein native HX, start from
%%%allH

N=100;  %number of exchangable H atom in the molecule
global k_ch k_op k_cl
k_ch=300;
k_op=0.1;
k_cl=10;

t=0:0.1:20;  %folding/HX time from 0~20 sec 

[t,y] = ode15s('nhxf',t,[k_cl/(k_cl+k_op) k_op/(k_cl+k_op) 0 0]);   %y(:,3)+y(:,4) is the labeled D(%) at different folding/HX time

sizer=size(t);
for i=1:sizer(1)
    
    %     for j=1:N         
    %         z=y(i,3)+y(i,4);
    %         if z>0.99
    %             z=0.99;    %to fix a bug upon numerical calculation
    %         end
    %         D(i,j)=(factorial(N)/(factorial(j-1)*factorial(N-j+1)))*(z^(j-1))*((1-z)^(N-j+1));    %D is the result matrix of ( folding/HX time vs. D-distribution ) 
    %     end
    
    M=1000; %number of molecules tested
    for j=1:M
        mass2(j)=0;
        ran2=rand(1,N);  %generate a random number between 0~1
        for k=1:N
            if ran2(k)<y(i,3)+y(i,4)
                mass2(j)=mass2(j)+1;
            end
        end
    end
    
    for w=1:N+1
        deltaD(w)=0;
        for q=1:M
            if mass2(q)==w-1
                deltaD(w)=deltaD(w)+1;
            end
        end
    end
    result(i,:)=conv(deltaC, deltaD/M);
    
end

sizer=size(result);

stem(0:sizer(2)-1,result(1,:))
hold on

stem(0:sizer(2)-1,result(15,:),'r')
hold on

stem(0:sizer(2)-1,result(44,:),'c')
hold on

stem(0:sizer(2)-1,result(201,:),'g')
hold on