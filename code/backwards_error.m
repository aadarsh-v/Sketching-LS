function be = backwards_error(A, b, y, theta)
%BACKWARDS_ERROR Summary of this function goes here
%   Detailed explanation goes here
    ynorm = norm(y);
    r = b - A*y;
    rnorm = norm(r);
    if isinf(theta)
        mu = 1;
    else
        mu = theta^2 * ynorm^2;
        mu = mu / (1 + mu);
    end

    phi = sqrt(mu) * rnorm / ynorm;
    be = min(phi, min(svd([A phi*(eye(size(A, 1)) - r*r'/rnorm^2)])));
end

