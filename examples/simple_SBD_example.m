clc; clear;

% Import Manopt and initialize the SBD package
run('../init_sbd');
fprintf('\n\n');

%% I. SIMULATE DATA FOR SBD:
%  =========================
% For the simple example, we will generate a single-slice, random kernel
% with small problem sizes for faster solution.

%% 1. Random kernel and observation
k = [8 8];           	% kernel size

% GENERATE
A0 = proj2oblique(randn(k));

%% 2. Activation map
m = [64 64];            % image size for each slice / observation grid

%   Each pixel has probability theta of being a kernel location
theta = 3e-2;           % activation concentration
eta = 1e-3;             % additive noise variance

% GENERATE
X0_good = false;
while ~X0_good
    X0 = double(rand(m) <= theta);              % activations are on / off
    X0_good = sum(X0(:) ~= 0) > 0;
end

b0 = randn;
Y = convfft2(A0, X0) + b0 +sqrt(eta)*randn(m);     % observation

%% II. Sparse Blind Deconvolution:
%  ===============================
% 1. SETTINGS - refer to documents on details for setting parameters.

% A function for showing updates as RTRM runs
dispfun = @( Y, A, X, k, kplus, idx ) showims(Y,A0,X0,A,X,k,kplus,idx);

params.lambda1 = 1e-1;              % regularization parameter for Phase I

params.phase2 = true;               % whether to do Phase II (refinement)
params.kplus = ceil(0.2 * k);       % padding for sphere lifting
params.lambda2 = 5e-2;              % FINAL reg. param. value for Phase II
params.nrefine = 3;                 % number of refinements

params.signflip = 0.2;              % want entrices of X to be nonnegative
params.xpos     = true;
params.getbias  = true;
params.Xsolve = 'FISTA';            % choose Xsolve: 'FISTA' or 'pdNCG'. 

% RUN SBD
[Aout, Xout, bout, extras] = SBD( Y, k, params, dispfun );
