%%%2009-12-24 fm_ppoe1.m:
%%%Folding Model PPOE 1: U<->I & I<->Ix & I<->N

function dy=fm_ppoe1(t,y)

global k_fm

kUI=k_fm(1); kIU=k_fm(2);
kIIx=k_fm(3); kIxI=k_fm(4);
kIN=k_fm(5); kNI=k_fm(6);

%%% y1:U y2:I y3:Ix y4:N
dy = zeros(4,1);

dy(1) = -kUI*y(1) +            kIU*y(2);
dy(2) =  kUI*y(1) - (kIN+kIU+kIIx)*y(2) + kIxI*y(3) + kNI*y(4);
dy(3) =                       kIIx*y(2) - kIxI*y(3);
dy(4) =                        kIN*y(2)             - kNI*y(4);