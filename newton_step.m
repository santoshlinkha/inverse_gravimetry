% =========================================================================
% Implement a damped Newton's method to solve the non-linear system
% =========================================================================
function curr_pts = newton_step(init, v)
    max_iter = 1000;
    damping = 0.1;
    h = 1e-6;
    n = numel(init)/4;
    curr_pts = init;

    pts_index = [1:n, n+2:4*n];

    for iter = 1:max_iter
        % disp(vpa([iter, curr_pts(2:end)], 4));
        jac = zeros(4*n, 4*n);

        curr_val = diff_quadrature(curr_pts, v);

        for i = 1:4*n
            t = zeros(1, 4*n);
            if i ~= n+1
                t(i) = h;
            end
            diff_val = (diff_quadrature(curr_pts + t, v) - curr_val) / h;
            jac(:, i) = diff_val.';
        end

        jac = jac(1:end-1, pts_index);  % remove row and column for fixed variable

        if rcond(jac) < 1e-15
            fprintf("Singular Jacobian %d iterations : %f\n", iter, det(jac));
            return;
        end

        step_dir = jac \ curr_val(1:end-1)';   % this may emit a warning

        next_step = curr_pts(pts_index) - damping * (step_dir)';

        cond1 = norm(next_step - curr_pts(pts_index));
        curr_pts(pts_index) = next_step;

        if cond1 < 1e-10
            fprintf("Converged after %d iterations\n", iter);
            break;
        end

    end
    fprintf("Stopped after %d iterations\n", iter);
end

% =========================================================================
% Define the z and c values as a function of x and lambda. 
% =========================================================================
function out = quadrature(d)
    n = numel(d)/4;
    d = d(:).'; % row vector
    x = d(1:n) + 1i*d(n+1:2*n);
    l = d(2*n+1:3*n) + 1i*d(3*n+1:4*n);

    z_val = zeros(1, n);
    c_val = zeros(1, n);

    for k = 1:n
        z_val(k) = sum( l(k) * conj(x) ./ (1 - l(k) * conj(l)) );
        c_val(k) = sum( x(k) * conj(x) ./ (1 - l(k) * conj(l)).^2 );
    end

    out = [real(z_val), imag(z_val), real(c_val), imag(c_val)];
end

% =========================================================================
% Define the difference between the computed and target z and c values
% =========================================================================
function out = diff_quadrature(d, v)
    out = quadrature(d) - v;
end