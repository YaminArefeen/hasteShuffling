function y_noise = addNoiseShuffle(y,noise_level,fill,gpuflag);
%This function adds white gaussian noise to the shuffling data object y,
%such that the amount of noise added resembles the amount of noise added in
%traditional haste

%Inputs
%y ~ M x N x N_c x L, Shuffling data to which we wish to add noise.
%noise_level ~ 1 x 1, Amount of noise we wish to add.
%fill ~ L x 1, lines acquired in the correct order
%Outputs
%y_noise ~ M x N x N_c x L, Noisy shuffling data 
if(nargin < 4)
  gpuflag = 0;
end

L = size(y,4); %Number of phase encodes acquired (i.e. number of echos) 
n_ro = size(y,1); %number of readoutpoints
n_co = size(y,3); %number of coils

temporary = sum(y(:,fill,:,:),4);
tmp = reshape(awgn(gather(temporary(:)),noise_level,'measured'),n_ro,L,n_co); 
%keep just the data acquired and add noise to it

if(gpuflag)
	y_noise = gpuArray(zeros(size(y))); %Preallocate matrix for y_noise 
	tmp = gpuArray(tmp);
else
	y_noise = zeros(size(y)); %Preallocate matrix for y_noise
end

for ii = 1:size(y,4)
    y_noise(:,fill(ii),:,ii) = tmp(:,ii,:);
end

end

