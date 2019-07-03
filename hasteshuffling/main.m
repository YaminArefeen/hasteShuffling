addpath('Utils')

%Loading Simulation Data
fprintf('Loading Simulation Data...')
load('data/sb256.mat')
load('data/fill.mat') %Phase encode points which will be acquired, not shuffled temporally
maps = sb256.esp;
echs = sb256.echs; %Simulated signal evolution from HASTE acquisition M x N x 140 
ES = sb256.ES;

%Constants
noiseflag = 1; %add white gaussisan noise to our simulated acquired data
bartflag = 0; %Whether we want to reconstruct with BART or my code
snr = 75;

minT2 = .05; %(s)
maxT2 = .4;  %(s)

dictlength = 10;
t2range = linspace(minT2,maxT2,dictlength); %T2 range for dictionary
rfangle = 160; %Flip angle used for simulation
B = 2; %Rank of low dimensional subspace

E = length(fill); %Lenth of echo train
[M,N,C] = size(maps); %Image dimensions

ech = echs(:,:,1:E); %Generating fully sampled dataset. 

%generating undersampled data M x N x C x E for forward model
fillshfl = fill(randperm(length(fill))); %temporally shuffle desired phase encode points, indicating which pe point acquired at each echo time
mask = zeros(M,N,C,E);
for e = 1:E
	mask(:,fillshfl(e),:,e) = 1;
end
fprintf('Done\n')

%Generating Low Dimensional Subspace
fprintf('Generating Low Dimensional Subspace...')
train = ones(E,1)*rfangle;
[Phi,~] = genFSEBasis(train,ES,1,t2range,1,1);
basis = Phi(:,1:B);
fprintf('Done\n')

%Creating Linear Operators
fprintf('Generating Linear Operators...')
op = defineSbOperators(maps,basis,mask);
fprintf('Done\n')

%Generating acquired data and kspadj
fprintf('Computing acquired data...')
y = op.R(op.F_for(op.S_for(ech)));
fprintf('Done\n')

%Defining parameters for iterative method
params.lambda = .5e-3;
params.n_iter = 200;
params.tolerance = 1e-3;
params.regularizor = 3;
params.residual_disp = 20;

if(bartflag) %If using Sidd's waveshuffling code in Bart
	coeffs = wshflBart(y,maps,basis,mask,params.n_iter,params.tolerance,params.lambda,[],snr);
else %Using my Matlab Reconstruction
	fprintf('Computing ksp adj...')
	y = addNoiseShuffle(op.R(op.F_for(op.S_for(ech))),snr,fillshfl);
	kadj = op.A_adj(y);
	fprintf('Done\n')
  
	%Computing Maximum eigenvalue for Stepsize
	fprintf('Computing Stepsize: ')
	stepsize = [];

	if(isempty(stepsize))
	stepsize = computeMaxEig(op.AhA,rand(size(kadj)),0);
	end
	fprintf('%f\n',stepsize)

	%Performing the Iterative Method
	data.y = y;
	data.ksp_adj = kadj;
	params.step_size = stepsize;

	history = fista(data,op,params);

	coeffs = squeeze(history.allcoeffs(:,:,:,end));
end
