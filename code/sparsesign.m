function S = sparsesign(d, m, zeta)
%SPARSESIGN Create a sparse sign matrix by specifying coordinates and
% values.
    rows = zeros(m*zeta, 1);
    for i = 1:m
        rows((i-1)*zeta+1:i*zeta) = randsample(d, zeta);
    end
    cols = kron((1:m)', ones(zeta, 1));
    signs = (2*randi(2, m*zeta, 1) - 3);
    S = sparse(rows, cols, signs / sqrt(zeta), d, m);
end
