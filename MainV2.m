clear
close all
clc

%% Get Flow and Structure parameters
rho = 1.225;

wh = 50;
wth = 100;
wbeta = 300;
mu = 40;
a = -0.4;
b = 1;
c = 0.6;
m = pi*rho*(b^2)*mu;
xth = 0.2;
xbeta=-0.025;
rth_sq = 0.25;
rbeta_sq = 0.00625;

%% Get Important Parameters
psi1 = 0.165;
psi2 = 0.335;
eps1 = 0.0455;
eps2 = 0.3;

%% Build structural matrices
[Ms, Ks] = GetStructuralMatrices(m, wh, wth, wbeta, a, b, c, xth, xbeta, rth_sq, rbeta_sq);

%% Build unsteady aerodynamic model
[Ma, Ca, Ka, W] = GetUnsteadyFlowForcesV2();

%% Setting Up A matrix

I1 = [1 0   0
      1 0   0
      0 1   0
      0 1   0
      0 0   1
      0 0   1];

Vm = linspace(0, 400, 801);
p = size(Vm);
Re = zeros(p(2), 12);
Im = zeros(p(2), 12);

for i = 1:p(2)
    V = Vm(i);
    M = Ms - 0.5.*rho.*Ma;
    C = 0.5.*(rho*V).*Ca;                  % Multiply by V
    K = Ks - 0.5.*(rho*V*V).*Ka;             % Multiply Aero by V^2
    W_f = (rho*b*b*V*V*V).*W;                % Multiply Aero by V^3
    W1 = [-eps1*V/b 0  0  0  0  0
      0   -eps2*V/b 0 0 0 0 
      0   0   -eps1*V/b 0 0 0
      0   0   0   -eps2*V/b  0   0
      0   0   0   0   -eps1*V/b  0
      0   0   0   0   0       -eps2*V/b];

    A = [-M\C      -M\K         -M\W_f  
         eye(3)             zeros([3 3])               zeros([3 6])
        zeros([6 3])                 I1                W1];
    
    [EigVec, EigVal] = eig(A);
    Re(i, :) = real(diag(EigVal));
    Im(i, :) = imag(diag(EigVal));
    
end

figure(1)
hold on
X = [0 0];
Y = [0 max(Im, [], "all")];
for i = 2:p(2)
    if (Re(i-1, 3) < 0 && Re(i+1, 3) > 0)
        idx = i;
    end
end

plot(X, Y, 'LineWidth', 1, 'Color','#EDB120','HandleVisibility','off')
hold on
for i = 1:p(2)
    plot(Re(i, 3),Im(i, 3),'g.','HandleVisibility','off')
    plot(Re(i, 5),Im(i, 5),'b.','HandleVisibility','off')
    plot(Re(i, 1),Im(i, 1),'r.','HandleVisibility','off')

end
xlim([-10 20])
ylim([0 400])
legend('A', 'B', 'C', 'D')
plot(Re(idx, 3), Im(idx, 3), 'ro', 'DisplayName','Flutter Point')
xlabel('Real (rad/s)')
ylabel('Imaginary (rad/s)')
title('Real vs Imaginary Part of Eigen Values')

figure(2)
hold on
for i = 1:p(2)
    plot(Vm(i),Im(i, 3),'g.','HandleVisibility','off')
    plot(Vm(i),Im(i, 5),'b.','HandleVisibility','off')
    plot(Vm(i),Im(i, 1),'r.','HandleVisibility','off')

end
legend('A', 'B', 'C', 'D')
plot(Vm(idx), Im(idx, 3), 'ro', 'DisplayName','Flutter Point')
xlabel('Velocity (m/s)')
ylabel('Imaginary (rad/s)')
title('Imaginary Part of Eigen Values vs Speed')

figure(3)
hold on
for i = 1:p(2)
    plot(Vm(i),Re(i, 3),'g.','HandleVisibility','off')
    plot(Vm(i),Re(i, 5),'b.','HandleVisibility','off')
    plot(Vm(i),Re(i, 1),'r.','HandleVisibility','off')

end
hold on
X = [0 400];
Y = [0 0];
plot(X, Y, '-k','HandleVisibility','off')
legend('A', 'B', 'C', 'D', 'E')
plot(Vm(idx), Re(idx, 3), 'ro', 'DisplayName','Flutter Point')
xlabel('Velocity (m/s)')
ylabel('Real (rad/s)')
title('Real Part of Eigen Values vs Speed')
