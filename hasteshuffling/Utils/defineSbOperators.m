function opsSb = defineSbOperators(maps,basis,mask,gpuflag)
%This function defines all the single basis operators for shuffling
%that will be used for the recon
%~~~~~~~~~Inputs~~~~~~~~~~~
%maps ~ M x N x C, Coil sensitivity maps
%basis ~ T x n_b, Temporal Basis
%Mask ~ M x N x C x T, undersampling mask for acquired data
%~~~~~~~~~Outputs~~~~~~~~~~
%opsSb, struct containing all of the operators associated with shuffling

if(nargin < 4)
  gpuflag = 0;
end

opsSb.S_for = @(x) bsxfun(@times,maps,permute(x,[1,2,4,3]));
%Forward coil sensitivity operator, should expand object in coil
%dimension
opsSb.S_adj = @(Sx) squeeze(sum(bsxfun(@times,conj(maps),Sx),3));
%Adjoint coil sensitivity operator, should collapse object in coil
%dimension

opsSb.T_for = @(x) temporal_forward(x,basis);
%Forward temporal operator.  Go from coefficient images to temporal images.
opsSb.T_adj = @(x) temporal_adjoint(x,basis);
%Adjoint temporal operator.  Go from temporal images to coefficient images.

opsSb.F_for = @(x) fft2c(x);
opsSb.F_adj = @(x) ifft2c(x);
%Forward and adjoint fourier operator

opsSb.R = @(x) bsxfun(@times,mask,x);
%Undersampling operator
opsSb.A_for = @(x) opsSb.R(opsSb.T_for(opsSb.F_for(opsSb.S_for(x))));
%Applies the forward operator of the shuffling system, R F S T x, which
%takes coefficients x (M x N x n_b) to acquired temporal kspace (M x N x C
%x T).  Note, the function actually applies A_for with the following
%operation order, R T F S x.  (FS) and (T) commute since each voxel sees
%the same basis.  This drastically saves on the number of fourier operators
%which we need to apply.
opsSb.A_adj = @(x) opsSb.S_adj(opsSb.F_adj(opsSb.T_adj(x)));
%Adjoint of the forward operator, x = T_H S_H F_H y.  Again, we apply it as 
%F_H S_H T_H to save computations and the operations commmute for the same
%reason as above.  

kernels = computeShufflingKernel(basis,mask,gpuflag);
opsSb.K = @(x) kernel_sh(x,kernels,gpuflag);
%The forward-adjoint operator (which will need to be computed in the first
%order iterative method used to solve the shuffling problem), will be
%applied as follows, F_H S_H T_H R T F S x.  We can precompute T_H R T as a
%n_b x n_b undersampling kernel applied to each voxel in the image.
%Kernels is this operator and K applies this operator. 

opsSb.AhA = @(x) opsSb.S_adj(opsSb.F_adj(opsSb.K(opsSb.F_for(...
    opsSb.S_for(x)))));
%forward adjoint operator as described above.

W = Wavelet;

opsSb.W_for = @(x) fW(x,W);
opsSb.W_adj = @(w) iW(w,W);
%Define wavelet operator for the sake of speed and see if I can get similar
%results to BART.
end


%forward wavelet transform for a set of images.
function w = fW(x,W)
    K = size(x,3); %Hard coding third dimension as number of images
     
    %Wavelet decomposition
    w = zeros(size(x));
    %wavelet transform of our coefficient images.
    for idx = 1:K
        w(:,:,idx) = W*x(:,:,idx);
    end  
end

%backwards wavelet transform for a set of images.
function x = iW(w,W,K)
    K = size(w,3); %Hard coding third dimension as number of images

    %wavelet reconstruction   
    x = zeros(size(w));
    %wavelet transform of our coefficient images.
    for idx = 1:K
        x(:,:,idx) = W'*w(:,:,idx);
    end         
end
