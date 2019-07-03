function step_size = computeMaxEig(A,input,verbose)
%Computes the maximum eigenvalue of an operator A, which is assumed to be 
%symmetric, in absolute value.  Then, will compute the maximum step size
%for an iterative method which repeatedly applies A.  Note, for the
%shuffling case, A = G_H G, where G is the forward shuffling operator.
%We utilize the power iteration mtehod.

%Inputs
%A ~ Operator which we want to compute the maximum eigen value of 
%input ~ random input to seed the power iteration algorithm
%verbose ~ print the delta change at each step
%Outputs
%step_size ~ Step_size to be used in iterative algorithm.

epsilon = 1e-2; %Convergence criteria for power iteration method.
delta = inf;
inputold = input/norm(input(:));

while delta > epsilon
    %apply our operator
    input = A(input);
    %normalize our vector
    input = input/norm(input(:));
    
    %Stopping condition check
    delta = norm(input(:) - inputold(:))/norm(inputold(:));
    if(verbose)
    	fprintf('Delta: %f\n',delta)
    end
    inputold = input;
end
max_eig = abs(input(:)'*reshape(A(input),size(input(:)))/(input(:)'*input(:)));
step_size = 1/(2*max_eig);
end
