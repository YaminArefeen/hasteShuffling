function [Phi,X] = genFSEBasis(train,ES,T1,T2,p,b1)
%Given a range of T1, T2 values and a particular refocusing train, this
%function determines a basis for all of the decay trains characterized by
%the combination of T1/T2 values and the refocusing train.
%
%Inputs
%train (R x 1) ~~ Refocusing train [degrees]
%ES (1 x 1) ~~ Echo spacing [seconds]
%T1 (M x 1) ~~ Range of T1 values [seconds]
%T2 (N x 1) ~~ Range of T2 values [seconds]
%p (S x 1) ~~ Range of proton density values [a.u.]
%
%Outputs
%Phi (R x R) ~~ Basis for the signal evolutions that could be generated from
%the range of T1, T2, and refocusing train values. 
%X (R x (M*N*S)) ~~ Matrix of all the signal evolutions generated, for
%debugging and evaluation purposes.

M = length(T1);
N = length(T2);
S = length(p); %number of T1,T2, and proton density respectively. 
R = length(train); %length of refocusing train.
B = length(b1);

train = train*pi/180;

X = zeros(R,M*N*S); %initial data matrix containing all of the decays

ctr = 1;
for m = 1:M %loop through each T1 value
    for n = 1:N %loop through each T2 value
        for s = 1:S %loop through each proton density value
            for b = 1:B
                X(:,ctr) = p(s)*FSE_signal(train,ES,T1(m),T2(n),b1(b)); 
                %assume b1 one while computing basis for now.
                %calculate the decay for this particular T1, T2, p set.
                ctr = ctr + 1;                
            end
        end
    end
end

[Phi,~,~] = svd(X, 'econ');
%compute the basis for the signal decay curves from x.
end

