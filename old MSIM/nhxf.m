function dy=nhxf(t,y)

%% y1:H_cl y2:H_op y3:D_op y4:D_cl
%% k1:k_op k2:k_cl k3:k_ch

global k_ch k_op k_cl

k1=k_op;
k2=k_cl;
k3=k_ch;

dy = zeros(4,1);


dy(1) = -k1*y(1) + k2*y(2);
dy(2) =  k1*y(1) - (k2+k3)*y(2); 
dy(3) =  k3*y(2) + k1*y(4) - k2*y(3);
dy(4) = k2*y(3) - k1*y(4);

