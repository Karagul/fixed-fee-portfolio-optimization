function Xkm1 = iterative_fee_reweighted_optimization(H,f,A,b,Aeq,beq,lb,ub,varargin)
% https://faculty.washington.edu/mfazel/portfolio-final.pdf
% Portfolio optimization with fixed and linear transaction costs

p = inputParser();
p.addOptional('fixed_fee',0);
p.addOptional('variable_fee',0);
p.addOptional('return_frac',0.50);
p.addOptional('credit',ones(size(f)));
p.addOptional('iter',10);
p.addOptional('eliminatedouplepos', 1);
p.addOptional('display', 'none');
p.addOptional('assert_valid', 0);
p.parse(varargin{:})

% First, solve the intial problem with bounds:
phi_complex_env = ( p.Results.fixed_fee./max(ub) + p.Results.variable_fee );
Xkm1 = triplesolver(H,f + phi_complex_env,A,b,Aeq,beq,lb,ub,p);

for k=1:p.Results.iter
    
    % Solve with the modified objective function with fixed fee costs 
    % "amortized" over the volume of the trade
    mod_objf = f + ( (p.Results.fixed_fee./(Xkm1'+1e-8)) + p.Results.variable_fee);%'
    
    if ~isempty(p.Results.credit)
        mod_objf = mod_objf ./ p.Results.credit;
    end
    
    [Xk, exit_flag] = triplesolver(H,mod_objf,A,b,Aeq,beq,lb,ub,p);
    
    if exit_flag ~= 1
        dispt('FAILED. REVERTING')
        Xk = Xkm1;
    end
    
    % Loop back:
    Xkm1 = Xk;
end



function [Xk, exit_flag] = triplesolver(H,mod_objf,A,b,Aeq,beq,lb,ub, p)
% Run the 3-stage optimization process:


% Find hte maximum return with the modified objective function:
%tic
opts = optimoptions('linprog','Display',p.Results.display, 'Algorithm','dual-simplex');
[~, max_ret] = linprog(mod_objf, A,b,Aeq,beq,lb,ub,opts);

div(sprintf('Max return for solver is %0.4g', max_ret))

% Constrain the next quadratic solver to meet atleast X% of total return:
Aeqq = [Aeq; mod_objf];
beqq = [beq; max_ret + abs(((1-p.Results.return_frac)*max_ret))];

% Quadratically solve the problem:
opts = optimoptions('quadprog','MaxIterations',100, 'Display', p.Results.display, 'ConstraintTolerance',1e-6);
[Xk, ~, exit_flag] = quadprog(H,[],A,b,Aeqq,beqq,lb,ub,[],opts);


% If we de dont have a double allocation problem, leave as is for speed:
%if sum(abs((Xk(1:end/2) - Xk(end/2+1:end)))) > 0.999999 * alloc
%    div('Skipped double allocation reduction optimization step')
%    return
%end

% Resolve so that we prevent any long/short double allocation:
if p.Results.eliminatedouplepos
    Xk = ((Xk(1:end/2) - Xk(end/2+1:end)));
    neg_pos = [Xk'<0, Xk'>0];
    ub(neg_pos) = 0;
    [Xk, ~, exit_flag] = quadprog(H,[],A,b,Aeqq,beqq,lb,ub,[],opts);
end

if p.Results.assert_valid
    assert(exit_flag==1)
end
%t=toc;
%div(sprintf('Elapsed time is %0.5g seconds',t))
					   
					    



