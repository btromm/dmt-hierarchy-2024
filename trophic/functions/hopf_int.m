function [FC,CV,Cvth,A]=hopf_int(gC,f_diff,sigma)
% Numerical integration of simplified, linearized Hopf model for simulating
% neural dynaic. The model calculates FC and temporal covariances of
% time-series
% gC = current effective connectivity matrix
% f_diff = freq diffs
% sigma = noise parameter
N=size(gC,1); %number of nodes
a=-0.02; %Proximity to bifurcation point
wo = f_diff'*(2*pi); % Angular frequencies for each node (f_diff -> freq diffs; 2*pi -> converts freq diff to angular)

Cvth = zeros(2*N); %Covariance matrix of hopf model

% Jacobian:

s = sum(gC,2); %Row sums of anatomical connectivity matrix Gc
B = diag(s); % diagonal matrix of row sums -> diagonal elements correspond to sum of connections for each node

Axx = a*eye(N) - B + gC; %top-left of Jacobian matrix. Combines damping parameter multiplied by identity matrix, connectivity matrix, and negative connectivity strengths
Ayy = Axx; %sets bottom-right block of Jacobian matrix to same as top-left. This is bcoz Hopf assumes symmetric dynamics for inhibitory nodes
Axy = -diag(wo); % Top-right block of Jacobian, diagonal are equal to negative angular freq (inhibitory self-coupling)
Ayx = diag(wo); % Bottom-left block of Jacobian, diagonal elements are positive angular freq (excitatory self-coupling)

A = [Axx Axy; Ayx Ayy]; % Construct full Jacobian
Qn = (sigma^2)*eye(2*N); % Noise cov matrix, diagonal matrix w/ variance sigma^2. Noise contribution to dynamics

Cvth=sylvester(A,A',-Qn); %Solves sylvester eq to obtain cov. matrix for Hopf model. Equation of form AX+XB=C where X is unknown matrix
FCth=corrcov(Cvth); % Pearson correlations between cov. matrix

FC=FCth(1:N,1:N); % indexed
CV=Cvth(1:N,1:N); %indexed