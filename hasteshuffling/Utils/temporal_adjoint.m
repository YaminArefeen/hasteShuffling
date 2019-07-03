function [alpha] = temporal_adjoint(x,basis)
%This function applies the temporal adjoint operator in a shuffling like
%vein. 
%Inputs
%x ~ (N x N x T) Time series of images 
%basis ~ (T x K) Orthonormal basis for signal decay curves
%Outputs
%alpha ~ (N x N x K) Coefficient images


if(length(size(x))==3 || length(size(x))==2)
    [M,N,T] = size(x);
    %M,N ~ image dimensions, T ~ number of decay curves
    K = size(basis,2);
    %number of coefficients

    x = reshape(x,M*N,T).'; 
    %reshape x so that we can apply adjoint through simple multiplication
    alpha = reshape((basis'*x).',M,N,K);
elseif(length(size(x)) == 4)
    [M,N,N_c,T] = size(x);
    K = size(basis,2);
    
    x = reshape(x,M*N*N_c,T).';
    alpha = reshape((basis'*x).',M,N,N_c,K);
end
end

