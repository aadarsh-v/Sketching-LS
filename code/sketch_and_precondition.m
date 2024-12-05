function [x, stats] = sketch_and_precondition(A, b, d, j, varargin)
%SKETCH_AND_PRECONDITION Summary of this function goes here
%   Detailed explanation goes here
    m = size(A, 1);
    n = size(A, 2);

    if length(varargin) >= 1 && ~isempty(varargin{1})
        summary = varargin{1};
    else
        summary = [];
    end
    if length(varargin) >= 2 && ~isempty(varargin{2})
        verbose = varargin{2};
    else
        verbose = false;
    end

    S = sparsesign(d,m,8);
    [Q, R] = qr(full(S*A), 'econ');
    y0 = Q' * (S*b); % preconditioned sketch and solve solution.

    if ~isempty(summary)
        summary = @(y) summary(R\y);
    end

    [y, ~, stats] = mylsqr(@(y) A*(R\y), @(y) R'\(A'*y), b, 0, j, ...
        summary, y0, verbose);
    x = R\y;
end

