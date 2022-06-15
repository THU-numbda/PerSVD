function [U, S, V] = rSVD_fp(A, k, p, Omg)
    %s = k/2;
    s = 50;
    l = k+s;
    [Q, ~] = qr(Omg, 0);
    for j = 1:p
        Y = A*Q;
        W = A'*Y;
        [Q, ~] = qr(W, 0);
        %[Q, ~, ~] = svd(W, 'econ');
    end
    Y = A*Q;
    W = A'*Y;
    [Q, S, V] = svd(Y, 'econ');
    B = (W*V*S^-1)';
    [U, S, V] = svd(B, 'econ');
    U = Q*U(:, 1:k);
    S = S(1:k, 1:k);
    V = V(:, 1:k);
end