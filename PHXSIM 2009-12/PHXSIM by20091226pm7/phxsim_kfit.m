
function f = phxsim_kfit(x)

global SFfitN
global foldingModel
global k_fm

switch foldingModel
    case 'iup1'
        k_fm=x;
        [t,y] = ode15s(@fm_iup1,SFfitN(:,1),[1 0 0]);
        z = y(:,3);
        
    case 'ppoe1'
        k_fm=x;
        [t,y] = ode15s(@fm_ppoe1,SFfitN(:,1),[1 0 0 0]);
        z = y(:,4);
        
    case 'ppoe2'
                
               
    otherwise
        error('Unkown folding model!')
end

f = z - SFfitN(:,2);