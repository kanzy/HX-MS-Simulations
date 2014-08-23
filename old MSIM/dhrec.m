T_folding=1.52e-2;
T_pulse=1.24e-2;

T1 = T_folding + T_pulse;
T2 = 60;
T3 = 300;

% for i=1:142
%     R(i)=exp(-(kc53(i)*T1+kc285_10(i)*T2+kc25(i)*T3));
% end

R=exp(-(kc53*T1+kc285_10*T2+kc25*T3));

recoveryD=sum(R)