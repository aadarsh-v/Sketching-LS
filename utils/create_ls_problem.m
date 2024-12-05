function [A,b,x,r] = create_ls_problem(m,n,cond_num,res_size)
%CREATE_LS_PROBLEM Creates a random LS problem with the specified 
% condition number and residual size. 
    U = randorth(m, n+1); % One extra column for generating b.
    V = randorth(n, n);
    D = diag(logspace(-log10(cond_num),0,n));
    A = U(:, 1:n)*D*V'; % Create A to have same condition number as D.
    
    x = orth(randn(n, 1));
    r = U(:, end) * res_size;
    b = A*x + r; % Has the correct residual as r is orthogonal with wanted norm.
end

