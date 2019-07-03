function history = fista(data,operators,parameters)
%This function attempts to implement FISTA.  Testing f r now,
%add more in depth comments later
%~~~~Unpack all parameters, operators, and data~~~~~~~~~~~~~~~~~
ksp_adj = data.ksp_adj; y = data.y; %assume normalized data.
AhA = operators.AhA; A_for = operators.A_for;
lambda = parameters.lambda; n_iter = parameters.n_iter; 
step_size = parameters.step_size; reg = parameters.regularizor;
res_disp = parameters.residual_disp;

displayFistaParameters(parameters);

%Initialize coefficient dimensions 
[M,N,K] = size(ksp_adj);

%Initializing values to save in history
residuals = ceil(n_iter/res_disp); %residual every res_disp iterations
del_params = ceil(n_iter/res_disp); %change in coeffs every res_disp iterations
all_coeffs = zeros(M,N,K,ceil(n_iter/res_disp)); %coeffs every res_disp iterations
res_ctr = 1; %counter which indexes into above arrays

fprintf('Starting FISTA...\n')

%Begin FISTA
coeffs = ksp_adj; %initialize first guess of coeffs
coeffsold = coeffs; 
z = coeffs; %define our initial auxillary variable 
tk = 1; 

for ii = 1:n_iter
    tic
    %Take our gradient and proximal step
    coeffs = z - step_size*AhA(z) + step_size*ksp_adj; %gradient
    
    %proximal
    if(reg == 1)            
        if(strcmp(class(coeffs),'gpuArray'))
          coeffs = gather(coeffs);
        end

        %Doing wavelet with BART
        w = wavbart(coeffs);
        w = reshape(SoftThresh(w(:),lambda),size(w));
        coeffs = wavbart(w,1,[M,N]);

      if(strcmp(class(coeffs),'gpuArray'))
        coeffs = gpuArrray(coeffs);
      end

    elseif(reg == 2)
        if(strcmp(class(coeffs),'gpuArray'))
          coeffs = gather(coeffs);
        end

        %locally low rank regularization
        coeffs = llr_thresh(coeffs,lambda,[8,8],1);

        if(strcmp(class(coeffs),'gpuArray'))
          coeffs = gpuArrray(coeffs);
        end
    elseif(reg == 3)
        if(strcmp(class(coeffs),'gpuArray'))
          coeffs = gather(coeffs);
        end

        %Doing wavelet with Miki's code
        w = operators.W_for(coeffs);
        w = reshape(SoftThresh(w(:),lambda),size(w));
        coeffs = operators.W_adj(w);        

        if(strcmp(class(coeffs),'gpuArray'))
          coeffs = gpuArrray(coeffs);
        end
	
    end
    
    %Incoorperating the FISTA aspect
    tkp1 = (1 + sqrt(1 + 4 * tk^2))/2;
    z = coeffs + (tk-1)/(tkp1)*(coeffs-coeffsold);
    
    if(mod(ii,res_disp) == 0) 
        del_param = norm(coeffs(:) - coeffsold(:))/norm(coeffsold(:));
    end
    
    coeffsold = coeffs;
    tk = tkp1;
        
    %displaying image each iteration
    display = abs(reshape(coeffs,M,N*K));
    imshow(display,[])
    drawnow
   
    if(mod(ii,res_disp) == 0)
      %Compute residual, comment out if it takes too much time
      forward = A_for(coeffs);
      residual = norm(y(:) - forward(:))/norm(y(:));
      
      residuals(res_ctr) = gather(residual);
      del_params(res_ctr) = gather(del_param);
      all_coeffs(:,:,:,res_ctr) = gather(coeffs);
      res_ctr = res_ctr + 1;
      time = toc;
      fprintf('Iteration: %d || Residual: %f || del_param: %f || Time: %f\n',ii,residual,del_param,time);     
      
    end
end

history.testcoeffs = z - step_size*AhA(z) + step_size*ksp_adj; %gradient
residuals = residuals(1:(res_ctr-1)); 
del_params = del_params(1:(res_ctr-1));
all_coeffs = all_coeffs(:,:,:,1:(res_ctr-1));

history.residuals = residuals;
history.allcoeffs = all_coeffs;
history.del_params = del_params;
end
