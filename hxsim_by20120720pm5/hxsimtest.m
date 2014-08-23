%%%2012-07-20 hxsimtest.m: for testing hxsim.m, hx2sim.m, hx1sim.m etc.


clear pepPara hxPara

SNase='ATSTKKLHKEPATLIKAIDGDTVKLMYKGQPMTFRLLLVDTPETKHPKKGVEKYGPEASAFTKKMVENAKKIEVEFDKGQRTDKYGRGLAYIYADGKMVNEALVRQGLAKVAYVYKGNNTHEQLLRKSEAQAKKEKLNIWSEDNADSGQ';

pepPara.proSeq=SNase; %whole protein sequence
pepPara.START=1; %Start residue# of the peptide
pepPara.END=15; %End residue# of the peptide
pepPara.Charge=2; %observing charge state of the peptide in MS spectrum

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sizer=pepPara.END-pepPara.START+1;
pepPara.kOP=[ones(1,5)*1e0, ones(1,sizer-5)*1e-1]; %array of openning rates of each residue of the peptide at HX condition
pepPara.kCL=ones(1,sizer)*1e2; %array of closing rates of each residue of the peptide at HX condition
pepPara.iniD=ones(1,sizer); %array of D% of each residue of the peptide at the beginning of HX
pepPara.iniF=ones(1,sizer)*2; %array of folding status(1=fully unfolded(open); 2=fully folded(close); 3=equilibrium by HX condition--assuming independent of other residues) of each residue of the peptide at the beginning of HX
pepPara.foldIndex=[ones(1,5), zeros(1,sizer-5)]; %array of foldon index of each residue(0=uncorrelated; 1,2,3...=diff foldons): residues in the same foldon must open/close together

hxPara.Temp=20; %HX temperature (unit: 'C)
hxPara.pH=8; %HX pH
hxPara.hxTime=1; %HX duration time (unit: sec)
hxPara.fractionD=0; %fraction of D2O in HX buffer

[HDmatrix, Fmatrix, obsPeaks, deuPeaks]=hxsim(pepPara,hxPara);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sizer=pepPara.END-pepPara.START+1;
% pepPara.kOP=[ones(1,5)*1e0, ones(1,sizer-5)*1e-1]; %array of openning rates of each residue of the peptide at HX condition
% pepPara.kCL=ones(1,sizer)*1e1; %array of closing rates of each residue of the peptide at HX condition
% pepPara.iniD=ones(1,sizer); %array of D% of each residue of the peptide at the beginning of HX
% pepPara.iniF=ones(1,sizer)*2; %array of folding status(1=fully unfolded(open); 2=fully folded(close); 3=equilibrium by HX condition--assuming independent of other residues) of each residue of the peptide at the beginning of HX
pepPara.pepPF=ones(1,sizer)*1e3; %protection factors(could be <1)
pepPara.foldIndex=[ones(1,5), zeros(1,sizer-5)]; %array of foldon index of each residue(0=uncorrelated; 1,2,3...=diff foldons): residues in the same foldon must open/close together
pepPara.kOPC=1e0; %opening rate of the correlated residue set
pepPara.iniD=1; %D% of every residue of the peptide at the beginning of HX

hxPara.Temp=20; %HX temperature (unit: 'C)
hxPara.pH=8; %HX pH
hxPara.hxTime=1; %HX duration time (unit: sec)
hxPara.fractionD=0; %fraction of D2O in HX buffer
hxPara.hxDir=2; %HX direction (1=H->D, 2=D->H)

[obsPeaks, deuPeaks, Distr1,Fract1,Distr2,Fract2]=hx1sim(pepPara,hxPara);









