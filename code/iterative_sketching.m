function [x,stats,iter] = iterative_sketching(A, b, d, varargin)
%ITERATIVE_SKETCHING Summary of this function goes here
%   Detailed explanation goes here
    m = size(A, 1);
    n = size(A, 2);

    if length(varargin) >= 1 && ~isempty(varargin{1})
        j = varargin{1};
    else
        j = eps();
    end
    if length(varargin) >= 2 && ~isempty(varargin{2})
        summary = varargin{2};
    else
        summary = [];
    end
    if length(varargin) >= 3 && ~isempty(varargin{3})
        verbose = varargin{3};
    else
        verbose = false;
    end

    stats = [];

    S = sparsesign(d,m,8);
    [Q, R] = qr(full(S*A), 'econ');
    x = R\(Q' * (S*b)); % sketch and solve solution.
    xold = x;
    rold = b - A*x;
    if ~isempty(summary)
        stats(end+1, :) = summary(x);
    end

    Acond = condest(R);
    if j < 1
        z = randn(n, 1);
        for i = 1:ceil(log(n))
            z = z / norm(z); z = R'*z;
            z = z / norm(z); z = R*z;
        end
        Anorm = norm(z);
    end

    iter = 1;
    while true
        xcopy = x;
        x = x + (R\(R'\(A'*rold)));
        xold = xcopy;
        r = b-A*x;
        updatenorm = norm(r-rold);
        rold = r;
        if ~isempty(summary); stats(end+1,:) = summary(x); end
        if verbose; fprintf('Iteration %d\t%e\n',iter,updatenorm); end
        iter = iter + 1;
        if (j >= 1 && iter > j) || (j < 1 &&...
                updatenorm <= j*(Anorm * norm(x) + 0.04*Acond*norm(r)))
            break
        elseif j < 1 && iter >= 100
            warning('Iterative sketching failed to meet tolerance')
            break
        end
    end


end

