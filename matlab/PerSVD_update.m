function [U, S, V] = PerSVD_update(A, k, p, Omg)
    %s = k/2;
    s = 50;
    l = k+s;
    [Q, ~] = qr(Omg, 0);
    alpha = 0;
    for j = 1:p
        Y = A*Q;
        W = A'*Y;
        D1 = W'*W;
        D2 = Y'*Y;
        I = eye(l);
        alpha1 = alpha;
        alpha2 = 0;
        for u = 1:100
            %[~, s, ~] = svd(A*(A'*Q)-alpha*Q, 'econ');
            M = D1-2*alpha*D2;
            [~, d] = eig(M);
            dl = d(1,1) + alpha*alpha;
            if dl <= 0
                break
            end
            sl = sqrt(dl);
            if alpha1 > sl
                break
            end
            alpha2 = alpha1;
            alpha1 = (alpha1 + sl)/2;
            if (alpha1-alpha2)/alpha1 < 1e-2
                break;
            end
        end
        alpha = alpha1;
        [Q, s, ~] = svd(W-alpha*Q, 'econ');
        
        
        
        if alpha < s(l, l)
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