% Constants start

% Linear mass density (u)
linear_mass_density = 0.01;

% Tension force (T)
tension_force = 100;

% viscous damping coefficent [N*s/m]
b = 0.5;

% String length (L)
string_length = 1;

%Matrix
M = 100;
N = 2000;

% Constants end

beta = b / 2 / linear_mass_density;
x = linspace(0, string_length, M)';
velocity = sqrt(tension_force / linear_mass_density);
dx = x(2) - x(1);
dt = dx / velocity;
p = (velocity * dt / dx)^2;
q = 1 + beta * dt;
u = 1 - beta * dt;

% boundary conditions
% (l, r) Dirichlet condition
% (h, k) Neumann condition

% Neumann - Neumann
s = sin(3 * pi * x + pi / 4);
g = zeros(M, 1);
h = 3 * pi * cos(pi * (1/4 + 3 * 0)) * ones(1, N);
k = 3 * pi * cos(pi * (1/4 + 3 * string_length)) * ones(1, N);

f(:, 1) = s;
f(1, 2) = p * (f(2, 1) - f(1, 1) - dx * h(1)) + f(1, 1) +u * dt * g(1);
f(2:M - 1, 2) = p / 2 * (f(3:M, 1) - 2 * f(2:M - 1, 1) + f(1:M - 2, 1)) + f(2:M - 1, 1) + u * dt * g(2:M - 1);
f(M, 2) = p * (dx * k(2) - f(M, 1) + f(M - 1, 1)) + f(M, 1) + u * dt * g(M);

for n = 2:N - 1
    f(1, n + 1) = 2 * p / q * (f(2, n) - f(1, n) - dx * h(n)) + 2 / q * f(1, n) - u / q * f(1, n - 1);
    f(2:M - 1, n + 1) = (1 / q) * p * (f(3:M, n) - 2 * f(2:M - 1, n) + f(1:M - 2, n)) + 2 * (1 / q) * f(2:M - 1, n) - (u / q) * f(2:M - 1, n - 1);
    f(M, n + 1) = (1 / q) * p * 2 * (dx * k(n) - f(M, n) + f(M - 1, n)) + 2 * (1 / q) * f(M, n) - (u / q) * f(M, n - 1);
end

plot_x = 0:1 / (M - 1):1;
plot_y = 0:1 / (M - 1):1;

figure(1)
surf(f);
zlabel('y');
xlabel('n');
ylabel('m');

% equilibrium is not possible as two neumanns are different

X = 1:M;
fig = figure(3);
axh = axes('Parent', fig);
i = 1;
ph = plot(axh, X, f(:, i));
xlim([0 M]);
ylim([-44 6]);

while i < N
    i = i + 1;
    set(ph, 'YData', f(:, i));
    pause(0.1)
end
