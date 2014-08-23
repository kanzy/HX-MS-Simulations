%%%2010-04-27 hxode.m: simplest H<->D equation; should be called by hxsim.m

function dy=hxode(t,y)

global hxodePara

kHD=hxodePara.kHD;
kDH=hxodePara.kDH;
fractionD=hxodePara.fractionD;

%%%y1:H, y2:D

dy = zeros(2,1);

dy(1) = -kHD*fractionD*y(1) + kDH*(1-fractionD)*y(2);
dy(2) =  kHD*fractionD*y(1) - kDH*(1-fractionD)*y(2);



