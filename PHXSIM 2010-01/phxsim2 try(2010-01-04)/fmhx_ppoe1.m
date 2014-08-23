%%%2010-01-04 fmhx_ppoe1.m: add hydrogen exchange
%%%2009-12-24 fm_ppoe1.m:
%%%Folding Model PPOE 1: U<->I & I<->Ix & I<->N

function dy=fmhx_ppoe1(t,y)

global k_fm
global kex

kUI=k_fm(1); kIU=k_fm(2);
kIIx=k_fm(3); kIxI=k_fm(4);
kIN=k_fm(5); kNI=k_fm(6);

kexHD_U=kex(1); kexDH_U=kex(2);
kexHD_I=kex(3); kexDH_I=kex(4);
kexHD_Ix=kex(5);kexDH_Ix=kex(6);
kexHD_N=kex(7); kexDH_N=kex(8);

%%% y1:U(H) y2:I(H) y3:Ix(H) y4:N(H)
%%% y5:U(D) y6:I(D) y7:Ix(D) y8:N(D)
dy = zeros(8,1);

dy(1) = -kUI*y(1) +            kIU*y(2)...
        -kexHD_U*y(1) + kexDH_U*y(5);
    
dy(2) =  kUI*y(1) - (kIN+kIU+kIIx)*y(2) + kIxI*y(3) + kNI*y(4)...
        -kexHD_I*y(2) + kexDH_I*y(6);
    
dy(3) =                       kIIx*y(2) - kIxI*y(3)...
        -kexHD_Ix*y(3) + kexDH_Ix*y(7);
        
dy(4) =                        kIN*y(2)             - kNI*y(4)...
        -kexHD_N*y(4) + kexDH_N*y(8);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
dy(5) = -kUI*y(5) +            kIU*y(6)...
        +kexHD_U*y(1) - kexDH_U*y(5);
    
dy(6) =  kUI*y(5) - (kIN+kIU+kIIx)*y(6) + kIxI*y(7) + kNI*y(8)...
        +kexHD_I*y(2) - kexDH_I*y(6);
    
dy(7) =                       kIIx*y(6) - kIxI*y(7)...
        +kexHD_Ix*y(3) - kexDH_Ix*y(7);
        
dy(8) =                        kIN*y(6)             - kNI*y(8)...
        +kexHD_N*y(4) - kexDH_N*y(8);
    
    
    
    
    
    