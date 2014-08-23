%%%2009-12-26 phxsim_sim.m: should be called by phxsim.m

kmax=max([max(kcDH), max(kcHD), max(k_fm)]);
deltaT=0.1/kmax;   %step size of simulation time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%the old implementation(before 2010-01-03):
% for time=0:deltaT:simTime
%     switch foldingModel
%         case 'iup1'
%             %call phxsim_iup1.m
%         case 'ppoe1'
%             phxsim_ppoe1 %call phxsim.ppoe1.m
%     end
%     disp(time) %display current time
% end

%%%2010-01-03 added: new implementation
simSteps=round(simTime/deltaT);
switch foldingModel
    case 'iup1'
        %call hxsim_iup1.mexw32
    case 'ppoe1'
        hxsim_ppoe1(simSteps,deltaT,START,END,Dfraction,pf_U,pf_I,pf_Ix,pf_N,kcHD,kcDH,...
            k_fm,FoldStates,HDmatrix,flagSim) %call hxsim_ppoe1.mexw32
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Distr=zeros(1,N+1);
deltaMass=zeros(1,M);
for i=1:M
    deltaMass(i)=sum(HDmatrix(i,:));
    Distr(deltaMass(i)+1)=Distr(deltaMass(i)+1)+1; %Distr(x) means x-1 units of mass above monoisotopic
end




