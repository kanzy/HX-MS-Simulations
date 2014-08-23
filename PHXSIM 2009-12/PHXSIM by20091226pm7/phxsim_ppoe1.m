%%%2009-12-26 phxsim_ppoe1.m: should be called by phxsim.m

ran1=rand(M,N);
ran2=rand(M,1);

for i=1:M
    switch FoldStates(i)
        case 0 %U
            pf=pf_U;
        case 1 %I
            pf=pf_I;
        case 2 %Ix
            pf=pf_Ix;
        case 3 %N
            pf=pf_N;
    end
    
    for j=START:END
        switch HDmatrix(i,j-START+1)
            case 0 %H
                if ran1(i,j-START+1)<=(1-exp(-kcHD(j)*deltaT/pf(j)))*Dfraction
                    HDmatrix(i,j-START+1)=1; %H->D
                end
            case 1 %D
                if ran1(i,j-START+1)<=(1-exp(-kcDH(j)*deltaT/pf(j)))*Hfraction
                    HDmatrix(i,j-START+1)=0; %D->H
                end
        end
    end
    
    if flagSim~=0
        switch FoldStates(i)
            case 0 %U
                if ran2(i)<=1-exp(-kUI*deltaT)
                    FoldStates(i)=1; %U->I
                end
            case 1 %I
                if ran2(i)<=1-exp(-kIU*deltaT)
                    FoldStates(i)=0; %I->U
                end
                if ran2(i)>1-exp(-kIU*deltaT) && ran2(i)<=(1-exp(-kIU*deltaT))+(1-exp(-kIN*deltaT))
                    FoldStates(i)=3; %I->N
                end
                if ran2(i)>(1-exp(-kIU*deltaT))+(1-exp(-kIN*deltaT)) && ...
                        ran2(i)<=(1-exp(-kIU*deltaT))+(1-exp(-kIN*deltaT))+(1-exp(-kIIx*deltaT))
                    FoldStates(i)=2; %I->Ix
                end
            case 2 %Ix
                if ran2(i)<=1-exp(-kIxI*deltaT)
                    FoldStates(i)=1; %Ix->I
                end
            case 3 %N
                if ran2(i)<=1-exp(-kNI*deltaT)
                    FoldStates(i)=1; %N->I
                end
        end
    end
    
end

