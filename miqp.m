function [X,  history] = miqp(Q, r, A, b, Aeq, beq, lb, ub, fmin, varargin)

p = inputParser();
p.addOptional('maxiter', 100)
p.addOptional('diff', 0.01)
p.parse(varargin{:});

N = length(r);

xvars = 1:N;
if fmin > 0
    vvars = N+1:2*N;
    zvar = 2*N+1;
else
    vvars = [];
    zvar = N+1;
end

lb = [colvect(lb); zeros(N+1,1)];
ub = [colvect(ub); ones(N+1,1)];
ub(zvar) = Inf;


fmax = 200;

%Atemp = eye(N);
%Amax = horzcat(Atemp,-Atemp*fmax,zeros(N,1));
%A = [[A, zeros(1,N+1)] ; Amax];
%b = [b;zeros(N,1)];
%Amin = horzcat(-Atemp,Atemp*fmin,zeros(N,1));
%A = [A;Amin];
%b = [b;zeros(N,1)];
Aeq = [Aeq, zeros(nrows(Aeq), N+1)];

Atemp = eye(N);
Amax = horzcat(Atemp,-Atemp*fmax,zeros(N,1));

if isempty(A)
    A = Amax;
else
    A = [ [A, zeros(nrows(A), N+1)]; Amax];
end
b = [b;zeros(N,1)];
Amin = horzcat(-Atemp,Atemp*fmin,zeros(N,1));
A = [A;Amin];
b = [b;zeros(N,1)];

%size(Aeq)
%Aeq = [[Aeq, zeros(nrows(Aeq),N+1)]; zeros(1,2*N+1)]; % Allocate Aeq matrix
%Aeq(end, xvars) = 1;
%beq = [beq; 1];

lambda = 1;

f = [-r,zeros(1, N),lambda];

options = optimoptions(@intlinprog, 'Display', 'none'); % Suppress iterative display
options = optimoptions(options,'LPOptimalityTolerance',1e-6,'MaxTime',200,...
                      'ConstraintTolerance',1e-3,'IntegerTolerance',1e-3);
[xLinInt,fval,exitFlagInt,output] = intlinprog(f,vvars,A,b,Aeq,beq,lb,ub,options);
last_xLinInt = xLinInt;
iter = 1; % iteration counter
assets = xLinInt(xvars); % the x variables
truequadratic = assets'*Q*assets;
zslack = xLinInt(zvar); % slack variable value

history = [];
%A = [A; [zeros(1,N*2), -1]];
%b = [b; -100];
                  
                  
while (abs((zslack - truequadratic)/truequadratic) > p.Results.diff & iter < p.Results.maxiter)% relative error
    newArow = horzcat(2*assets'*Q,zeros(1,N),-1); % Linearized constraint
    rhs = assets'*Q*assets;                       % right hand side of the linearized constraint
    A = [A;newArow];
    b = [b;rhs];
    
    % Solve the problem with the new constraints
    options = optimoptions(options, 'MaxTime', 10);
    tic
    [xLinInt,fval,exitFlagInt,output] = intlinprog(f,vvars,A,b,Aeq,beq,lb,ub,options);
    
    toc
    if exitFlagInt == 0
        warning('Infeasible...')
        xLinInt = last_xLinInt;
    end
    assets = (assets+xLinInt(xvars))/2; % Midway from the previous to the current
    % assets = xLinInt(xvars); % Use the previous line or this one
    
    truequadratic = xLinInt(xvars)'*Q* xLinInt(xvars);
    zslack = xLinInt(zvar);
    history = [history;truequadratic,zslack];
    iter = iter + 1;
    X = xLinInt(xvars);
    tol = abs((zslack - truequadratic)/truequadratic);
    fprintf('Iteration %i: True Quad: %0.5g, tol %0.3g\n',iter,truequadratic, tol)
    last_xLinInt = xLinInt;
end