%%%2009-12-27 phxsim_cafit.m: modify msfit6_fit4.m for phxsim.m
%%%2009-11-09 msfit6_fit4.m: modify msfit6_fit3.m for msfit6f.m

function y=phxsim_cafit(x,DistrA_U,DistrA_N,DistrA_x,weights)

y=weights.*(DistrA_x - (x(1)*DistrA_U + x(2)*DistrA_N));