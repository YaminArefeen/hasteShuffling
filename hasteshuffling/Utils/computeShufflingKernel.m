function [kernels] = computeShufflingKernel(basis,mask,gpuflag)
%Computes the sampling kernel I^H R I in the shuffling foward adjoint
%operations S_H F_H I_H R I F S x 
%Inputs
%basis ~ T x n_b, temporal basis for shuffling
%mask ~ n_ro x n_pe x n_ch x T, full blown undersampling mask
%gpuflag ~ Indicating whether the operation will be performed in the gpu or not
%Outputs
%kernels ~ n_b x n_b x n_pe

[~,n_pe,~,num_acq] = size(mask);
n_b = size(basis,2);

%For single shot, only need to generate a different sampling kernel for
%each phase encode
kernels = zeros(n_b,n_b,n_pe);
if(gpuflag)
	kernels = gpuArray(kernels);
end

m = squeeze(mask(1,:,1,:)); 

for ii = 1:n_pe
    
    R_tmp = zeros(num_acq,num_acq);
    if(gpuflag)
	R_tmp = gpuArray(R_tmp);
    end

    echo_sampled = find(m(ii,:));
    
    %If phase encode not sampled, sampling kernel becomes 0
    if(~isempty(echo_sampled))
        for ech = 1:length(echo_sampled)
            R_tmp(echo_sampled(ech),echo_sampled(ech)) = 1;
        end
    end
    
    kernels(:,:,ii) = basis'*R_tmp*basis;
end
end

