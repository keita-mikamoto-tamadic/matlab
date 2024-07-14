function P = Ricatti_eq_lqr(A, B, Q, R)
%RICATTI_EQ_LQR Ricatti方程式の解
%   有本ポッター法

    dim_x = size(A,1);
    dim_u = size(B,2);

    % Set Hamilton matrix
    Ham = zeros(2*dim_x, 2*dim_x);
    Ham(1:dim_x, 1:dim_x) = A;
    Ham(1:dim_x, dim_x + 1:end) = -B * inv(R) * B';
    Ham(dim_x + 1:end, 1:dim_x) = -Q;
    Ham(dim_x + 1:end, dim_x + 1:end) = -A';


    % Calculate eigenvalues and eigenvectors
    [eigvec, eigval] = eig(Ham);

    % Extract stable eigenvectors
    stable_eigvec = [];
    for i = 1:2 * dim_x
        if real(eigval(i, i)) < 0
            stable_eigvec = [stable_eigvec, eigvec(:, i)];
        end
    end

    % Calculate P with stable eigenvectors
    Vs_1 = stable_eigvec(1:dim_x, :);
    Vs_2 = stable_eigvec(dim_x + 1:end, :);
    P = real(Vs_2 * inv(Vs_1));

end

