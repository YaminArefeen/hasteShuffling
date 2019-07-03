function out = kernel_sh(x,kernels,gpuflag)
%Applies the sampling kernel I^H R I in the shuffling foward adjoint
%operations S_H F_H I_H R I F S x 
%Inputs
%kernels ~ n_b x n_b x n_pe
%x ~ n_ro x n_pe x n_ch x n_b, Coefficients taken through fourier
%and coil sensitivity operators
%Outputs
%out ~ n_ro x n_pe x n_ch x n_b

[n_ro,n_pe,n_ch,n_b] = size(x);

out = zeros(n_ro,n_pe,n_ch,n_b);
if(gpuflag)
	out = gpuArray(out);
end

for ii = 1:n_pe
    cur = reshape(x(:,ii,:,:),n_ro*n_ch,n_b).';
    out(:,ii,:,:) = reshape((kernels(:,:,ii)*cur).',n_ro,1,n_ch,n_b);
end
end
