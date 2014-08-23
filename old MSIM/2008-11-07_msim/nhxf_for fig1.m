function dy=nhxf(t,y)

%% y1:Hcl y2:Hop y3:D
%% k1:kop k2:kcl k3:kch
k1=2;
k2=3000;
k3=300;

dy = zeros(3,1);


dy(1) = -k1*y(1) + k2*y(2);
dy(2) =  k1*y(1) - (k2+k3)*y(2); 
dy(3) =  k3*y(2);

