function out = wshflBart(y,maps,basis,mask,n_iter,tolerance,lambda,max_eig,snr)
%In this function, we will perform a first iteration shuffling
%reconstruction using Sidd's Bart code, just to see how fast it runs on gpu
%~~~~~~~~inputs~~~~~~~~~
%y ~ M x N x C x T, Acquired data
%maps ~ M x N x C, Coil sensitivity profiles
%basis ~ T x B, Temporal shuffling basis
%mask ~ M x N x C x T, sampling mask indicating measured samples
%n_iter ~ number of iterations
%tolerance ~ tolerance of residual
%lambda ~ regularization parameter
%max_eig ~ maximum eigen value of AhA if I have it precomputed
%snr ~ Noise level if performing a simulation reconstruction
%~~~~~~~~Outputs~~~~~~~~
%out ~ M x N x B, Shuffling coefficients

if(nargin < 9)
	snr = [];
end

if(nargin < 9)
	max_eig = [];
end

%Create reorder and table
[M,N,C,T] = size(y);
B = size(basis,2);
m = squeeze(mask(1,:,1,:));
n = sum(m(:)); %total number of ky points sampled

reorder = zeros(n,3); %(#ky points sampled by 3)
table = zeros(M,C,n);
ctr = 1; %idx through reorder
for pe = 1:N
    vec = m(pe,:);
    jdx = find(vec(:) == 1);
    %find at which echo times this pe point is sampled
    
    if(isempty(jdx))%if phase encode is not sampled
        continue
    end
    
    for  j = 1:length(jdx)
        reorder(ctr,1) = pe - 1;
        reorder(ctr,3) = jdx(j) - 1;
	table(:,:,ctr) = y(:,pe,:,jdx(j));
        ctr = ctr + 1;
    end
end

if(~isempty(snr))
	table = reshape(awgn(table(:),snr),size(table));
end

%reshape things into appropriate sizes
psf_bart = ones(M,N,1,1,1,1,1);
phi_bart = reshape(basis,1,1,1,1,1,T,B);
maps_bart = reshape(maps,M,N,1,C,1,1,1);

if(isempty(max_eig))
	out = squeeze(bart(sprintf('wshfl -H -f -i %d -t %f -r %f -w',...
	    n_iter,tolerance,lambda),maps_bart,psf_bart,phi_bart,reorder,table));
else
	out = squeeze(bart(sprintf('wshfl -H -f -e %f -i %d -t %f -r %f -w',...
	    max_eig,n_iter,tolerance,lambda),maps_bart,psf_bart,phi_bart,reorder,table));i
end
