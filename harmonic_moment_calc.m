clc;
clear;

disp("Section 0 Start")
% =========================================================================
% 1. Define paramters
% 2. Calculate Harmonic Moments
% =========================================================================
a = 3;
b = 2;
n = 5;

harmonic_moments = zeros(3, 2*n);
center = 0 + 0i;
for k=0:(2*n-1)    
    harmonic_moments(1, k+1) = integral2( ...
            @(r, t) (a*r.*cos(t) + 1i*b*r.*sin(t) + center).^k.*r * a*b, ...
            0, 1, 0, 2*pi, 'AbsTol',1e-4,'RelTol',1e-2 ...
        );
end


prony_mat = zeros(n+1, n+1);

for i=1:n
    for j=1:n+1
        prony_mat(i,j) = harmonic_moments(1, i+j-1);
    end
end

disp("Section 1 Complete")

coeffs = zeros(1, n+1);

for i=1:n+1
    coeffs(i) = (-1)^(n+i+1)*det(prony_mat(1:n, [1:i-1, i+1:end]));
end

coeffs = flip(coeffs);


% =========================================================================
% 1. Define Pronys matrix
% 2. Solve for the z_values and c_values
% =========================================================================

zList = roots(coeffs);


c_mat = ones(n, n);

for i=2:n
    c_mat(i, :) = zList.'.^(i-1);
end

cList = c_mat \ harmonic_moments(1, 1:n).';
cList = cList/pi;


for l = 0:2*n - 1
    harmonic_moments(2, l+1) = sum(cList .* zList.^l);
end

harmonic_moments(3, :) =  abs(harmonic_moments(2, :) - harmonic_moments(1, :));


disp("Section 2 Complete");

% =========================================================================
% 1. Use the given equation to solve for the x values and lambda values
% 2. Use fsolve to solve for the non-linear complex equations
% =========================================================================

initGuess = zeros(2*n,1);

w_n = abs(sum(cList));
r_k = sqrt(abs(cList)/pi);
initGuess(1:n) = sqrt(w_n) * r_k / sum(sqrt(abs(cList)));


% Solve with fsolve (suppress output)
options = optimoptions('fsolve', 'Display', 'none', 'MaxFunEvals',1e10,'MaxIter',1e10, 'TolFun', 1e-9);
[sol,fval,exitflag] = fsolve(@(vars) regionPts(vars, zList, cList, n), initGuess, options);

% Extract solutions
xSol = sol(1:n);
lSol = sol(n+1:2*n);

disp("x values");
disp(xSol);

disp("lambda values");
disp(lSol);

disp("Section 3 Complete")
% =========================================================================
% This section of the code is to verify the if the computed values of x and 
% lambda satisfies the equations or not. There might be slight error in the
% solution but as long as it is within the tol range, it is fine.
% =========================================================================

zCheck = zeros(1,n);
for k = 1:n
    sumExpr = 0;
    for j = 1:n
        sumExpr = sumExpr + conj(xSol(j))*lSol(k) / (1 - lSol(k)*conj(lSol(j)));
    end
    zCheck(k) = sumExpr;
end

% c equations
cCheck = zeros(1,n);
for k = 1:n
    sumExpr2 = 0;
    for j = 1:n
        sumExpr2 = sumExpr2 + (conj(xSol(j))*xSol(k)) / (1 - lSol(k)*conj(lSol(j)))^2;
    end
    cCheck(k) = sumExpr2;
end

disp("z values (Computed vs Target)")
disp([zCheck(:), zList(:)]);

disp("c values (Computed vs Target)")
disp([cCheck(:), cList(:)]);

disp("Section 4 Complete")

% ========================================================================
% 1. Define the points on a unit disc.
% 2. Map the ponits on the disc using the conformal map.
% 3. Plot the region on the complex plane.
% ========================================================================

no_of_pts = 1000;
disc_pts = exp(1i*linspace(0, 2*pi, no_of_pts+1));

% figure;
% plot(real(disc_pts), imag(disc_pts), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5); 
% grid on; axis equal;


region_pts = arrayfun(@(z) mapRegion(z, xSol, lSol), disc_pts);


figure;

ellipse_x = a*cos(linspace(0, 2*pi, no_of_pts)) + real(center);
ellipse_y = b*sin(linspace(0, 2*pi, no_of_pts)) + imag(center);

% plot(ellipse_x, ellipse_y, 'b:');
% hold on;
plot(real(region_pts), imag(region_pts), 'r');
grid on; axis equal;

disp("Section 5 Complete")

% ========================================================================
% Define and plot a conformal map based on the computed x and lambda values
% ========================================================================


function G = mapRegion(z, xSol, lSol)
    n = length(xSol);
    G = 0;   % initialize accumulator
    for i = 1:n
        G = G + conj(xSol(i)) * z / (1 - conj(lSol(i)) * z);
    end
end


% =========================================================================
% This calculates the residue for the equation of the system. F solve is 
% used on this to estimate the values of x and lambda. Above code checks it
% =========================================================================
function F = regionPts(vars, zList, cList, n)
    % Split unknowns into complex vectors
    % x = vars(1:n) + 1i*vars(n+1:2*n);
    % l = vars(2*n+1:3*n) + 1i*vars(3*n+1:4*n);
    x = vars(1:n);
    l = vars(n+1:2*n);

    % Initialize residuals
    zEqns = zeros(1,n);
    cEqns = zeros(1,n);

    for k = 1:n
        % z equations
        sumExpr = 0;
        for j = 1:n
            sumExpr = sumExpr + conj(x(j))*l(k) / (1 - l(k)*conj(l(j)));
        end
        zEqns(k) = sumExpr - zList(k);

        % c equations
        sumExpr2 = 0;
        for j = 1:n
            sumExpr2 = sumExpr2 + (conj(x(j))*x(k)) / (1 - l(k)*conj(l(j)))^2;
        end
        cEqns(k) = sumExpr2 - cList(k);
    end
    F = [zEqns cEqns];
end