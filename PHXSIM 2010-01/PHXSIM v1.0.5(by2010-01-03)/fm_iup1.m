%%%2009-12-24 fm_iup1.m:
%%%Folding Model IUP 1: U<->N & U<->I & I<->N

function dy=fm_iup1(t,y)

global k_fm

kUN=k_fm(1); kNU=k_fm(2);
kUI=k_fm(3); kIU=k_fm(4);
kIN=k_fm(5); kNI=k_fm(6);

%%% y1:U y2:I y3:N
dy = zeros(3,1);

dy(1) = -(kUN+kUI)*y(1) +       kIU*y(2) +       kNU*y(3);
dy(2) =        kUI*y(1) - (kIN+kIU)*y(2) +       kNI*y(3);
dy(3) =        kUN*y(1) +       kIN*y(2) - (kNU+kNI)*y(3);
