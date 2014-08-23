%%mism.m 2008-11-07: Simulation of Mass Spec peaks distribution of a
%%protein or peptide due to isotope C13 (nature richness 0.11%)

clear

p=0.011; %probability of C13 inclusion (natural richness of C13)
n=150;  %number of C atom in protein

m=n;
if(n>20) m=20; end

for i=1:m  %observed m/z peaks
    intensity(i)=(factorial(n)/(factorial(i-1)*factorial(n-i+1)))*(p^(i-1))*((1-p)^(n-i+1));
end

X=0:(m-1);
stem(X,intensity)

hold on

