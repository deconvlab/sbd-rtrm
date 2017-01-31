function [ Aout, Xout, extras ] = SBD( Y, k, params, dispfun )
%SBD Summary of this function goes here
%
%   PARAMS STRUCT:
%   ===============
%   The options struct should include the fields:
%       lambda1,  float > 0  : regularization parameter for Phase I
%       phase2,   bool       : whether to do Phase II (refinement) or not
%       
%   IF phase2 == true, then the following fields should also be included:
%       kplus,    int > 0    : border padding (pixels) for sphere lifting
%       lambda2,  float > 0  : FINAL reg. param. value for Phase II
%       
%       nrefine,  int >= 1   : number of refinements for Phase II.
%           Refinement 1 lifts the sphere and uses lambda1, successive 
%           refinements decrease lambda down to lambda2; 
%           i.e. if nrefine == 1, then no decrease in lambda is made.
%           
%
%   Finally, two optional fields for the struct. These features are
%   automatically disabled if the fields are not included or are empty:
%
%       signflip, float      : attempts to choose the sign of A and X after 
%           each ASolve so the majority of activations in X are positive. 
%           
%           Setting signflip < 0 disables this feature. 
%           
%           Setting signflip >= 0 considers the 'activations' of X to be 
%           entries with abs. value geq signflip * max(X(:)), then 
%           checking forthe potential sign flip.
%
%
%       xpos,     bool    :  Constrain X to have nonnegative entries
%           when running XSolve during Phase II. Should be used in
%           conjunction with signflip to prevent a trivial X from being
%           returned pathologically.
%


%% Process input arguments
starttime = tic;
n = size(Y,3);

if nargin < 4 || isempty(dispfun)
    dispfun = @(Y,A,X,k,kplus,idx) 0;
end

lambda1 = params.lambda1;
if params.phase2
    kplus = params.kplus;
    lambda2 = params.lambda2;
    nrefine = params.nrefine;
end

if ~isfield(params, 'signflip') || isempty(params.signflip)
    signflip = -1;
else
    signflip = params.signflip;
end

if ~isfield(params, 'xpos') || isempty(params.xpos)
    xpos = 0;
else
    xpos = params.xpos;
end

%% PHASE I: First pass at BD
dispfun1 = @(A, X) dispfun(Y, A, X, k, [], 1);

fprintf('PHASE I: \n=========\n');
A = randn([k n]); A = A/norm(A(:));

[A, Xsol, info] = Asolve_Manopt( Y, A, lambda1, [], dispfun1);
extras.phase1.A = A;
extras.phase1.X = Xsol.X;
extras.phase1.info = info;

%% PHASE II: Lift the sphere and do lambda continuation
if params.phase2
    k2 = k + 2*kplus;
    dispfun2 = @(A, X) dispfun(Y, A, X, k2, 0, 1);
    
    A2 = zeros([k2 n]);
    A2(kplus(1)+(1:k(1)), kplus(2)+(1:k(2)), :) = A;
    X2sol = Xsol;
    %X2sol.X = circshift(Xsol.X,-kplus);
    %X2sol.W = circshift(Xsol.W,-kplus);
    % clear A Xsol;

    lambda = lambda1; 
    score = zeros(2*kplus+1);
    fprintf('\n\nPHASE II: \n=========\n');
    lam2fac = (lambda2/lambda1)^(1/nrefine);
    i = 1;
    while i <= nrefine + 1
        fprintf('lambda = %.1e: \n', lambda);    
        [A2, X2sol, info] = Asolve_Manopt( Y, A2, lambda, X2sol, dispfun2 );
        fprintf('\n');

        %Attempt to 'unshift" the a and x by taking the l1-norm over all k-contiguous elements:
        for tau1 = -kplus(1):kplus(1)
            ind1 = tau1+kplus(1)+1;
            for tau2 = -kplus(2):kplus(2)
                ind2 = tau2+kplus(2)+1;
                temp = A2(ind1:(ind1+k(1)-1), ind2:(ind2+k(2))-1,:);
                score(ind1,ind2) = norm(temp(:), 1);
            end
        end
        [temp,ind1] = max(score); [~,ind2] = max(temp);
        tau = [ind1(ind2) ind2]-kplus-1;
        A2 = circshift(A2,-tau);
        X2sol.X = circshift(X2sol.X,tau);
        X2sol.W = circshift(X2sol.W,tau);

        % Save phase 2 extras:
        if i == 1;  idx = 1;    else idx = i;    end
        extras.phase2(idx).A = A2;
        extras.phase2(idx).X = X2sol.X;
        extras.phase2(idx).info = info;
        if i == 1;  extras.phase2 = fliplr(extras.phase2);  end
        
        dispfun2(A2,X2sol.X);
        lambda = lambda*lam2fac;
        i = i+1;
        
    end
end

%% Finished: get the final A, X
Aout = A2(kplus(1)+(1:k(1)), kplus(2)+(1:k(2)), :);
Xout = circshift(X2sol.X,kplus) * norm(Aout(:));
Aout = Aout/norm(Aout(:));

if signflip
    thresh = 0.2*max(abs(Xout(:)));
    sgn = sign(sum(Xout(abs(Xout) >= thresh)));
    Aout = sgn*Aout;
    Xout = sgn*Xout;
end

runtime = toc(starttime);
fprintf('\nDone! Runtime = %.2fs. \n\n', runtime);
end

