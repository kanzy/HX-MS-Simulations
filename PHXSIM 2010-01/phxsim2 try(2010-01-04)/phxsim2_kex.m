%%%2010-01-04 phxsim2_kex.m: prepare 'kex', should be called by phxsim2.m


switch foldingModel
    case 'iup1'
        %add code here
        
    case 'ppoe1'
        kexHD_U = Hfraction*kcHD(currResidue)/pf_U(currResidue);
        kexDH_U = Dfraction*kcDH(currResidue)/pf_U(currResidue);
        kexHD_I = Hfraction*kcHD(currResidue)/pf_I(currResidue);
        kexDH_I = Dfraction*kcDH(currResidue)/pf_I(currResidue);
        kexHD_Ix = Hfraction*kcHD(currResidue)/pf_Ix(currResidue);
        kexDH_Ix = Dfraction*kcDH(currResidue)/pf_Ix(currResidue);
        kexHD_N = Hfraction*kcHD(currResidue)/pf_N(currResidue);
        kexDH_N = Dfraction*kcDH(currResidue)/pf_N(currResidue);
        kex=[kexHD_U,kexDH_U,kexHD_I,kexDH_I,kexHD_Ix,kexDH_Ix,kexHD_N,kexDH_N]; %for fmhx_ppoe1.m use
        
        
        
    otherwise
        error('Unknown folding model!')
end