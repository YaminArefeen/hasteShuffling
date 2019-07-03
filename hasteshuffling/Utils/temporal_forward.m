function [x] = temporal_forward(alpha,basis)
%This function applies the temporal forward operator in the shuffling like
%vein
%Inputs
%alpha ~ (N x N x K) Coefficient images
%basis ~ (T x K) Orthonormal basis for signal decay curves
%Outputs
%x ~ (N x N x T) Time series of images 

if(length(size(alpha))==3 || length(size(alpha)) == 2)
    [M,N,K] = size(alpha);
    %M,N ~ image dimensions, K ~ number of coefficients
    T = size(basis,1);
%number of images along our decay curve.

    a = reshape(alpha,M*N,K).'; % K x M*N object so that we can project
    %by simple multiplication of basis. 
    x = reshape((basis*a).',M,N,T);
elseif(length(size(alpha))==4)
    [M,N,N_c,K] = size(alpha);
    T = size(basis,1);
    a = reshape(alpha,M*N*N_c,K).';
    x = reshape((basis*a).',M,N,N_c,T);
end
end

