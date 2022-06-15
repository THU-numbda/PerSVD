function [U, S, V] = PerSVD_once(A, k, p, Omg)
    %s = k/2;
    s = 50;
    l = k+s;
    [Q, ~] = qr(Omg, 0);
    alpha = 0;
    for j = 1:p
        Y = A*Q;
        W = A'*Y;
        [Q, s, ~] = svd(W-alpha*Q, 'econ');
        %if alpha < s(l, l)
        if alpha == 0
            alpha = (alpha + s(l, l))/2;
            if alpha < 1e-16
                alpha = 0;
            end
        end
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