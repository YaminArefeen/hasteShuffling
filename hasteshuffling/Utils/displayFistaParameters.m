function [] = displayFistaParameters(p)
%This function prints out all the parameters I use in my current
%implementation of FISTA for sake of visualization while recon is going on
lambda = p.lambda; n_iter = p.n_iter; 
step_size = p.step_size; reg = p.regularizor;

if(reg == 1)
    regstr = 'BART Wavelet';
elseif(reg == 2)
    regstr = 'LLR';
elseif(reg == 3)
    regstr = 'MIKI Wavelet';
elseif(reg == 5 || reg == 0)
    regstr = 'Least Squares';
end

fprintf('Prior: %s\n',regstr)
fprintf('Lambda: %d\n',lambda)
fprintf('Iteratons: %d\n',n_iter)
fprintf('Step Size: %f\n',step_size)
end
