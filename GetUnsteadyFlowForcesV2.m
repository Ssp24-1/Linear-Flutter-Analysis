function [Ma, Ca, Ka, W] = GetUnsteadyFlowForcesV2()

a = -0.4;
b = 1;
c = 0.6;

psi1 = 0.165;
psi2 = 0.335;
eps1 = 0.0455;
eps2 = 0.3;


T1 = -(1/3)*(2+c^2)*((1-c^2)^0.5) + c*acos(c);
T3 = -(1/8)*(1-c^2)*(5*c^2+4) + (1/4)*c*(7+2*c^2)*((1-c^2)^0.5)*acos(c) - (c^2+1/8)*(acos(c))^2;
T4 = -acos(c)+c*((1-c^2)^0.5);
T5 = -(1-c^2)-((acos(c))^2)+(2*c*acos(c)*((1-c^2)^0.5));
T7 = (1/8)*c*(7+2*c^2)*((1-c^2)^0.5) - (c^2+(1/8))*acos(c);
T8 = -(1/3)*((2*c^2)+1)*((1-c^2)^0.5) + c*acos(c);
T9 = (1/2)*((1/3)*((1-c^2)^0.5) + a*T4);
T10 = ((1-c^2)^0.5) + acos(c);
T11 = acos(c)*(1-2*c) + ((1-c^2)^0.5)*(2-c);
T12 = ((1-c^2)^0.5)*(2+c) - acos(c)*(1+2*c);
T13 = -0.5*(T7 + (c-a)*T1);

Wag_0 = 1-psi1-psi2;

%% Aerodynamic Mass Matrix
Ma = (b^2)*[pi    -pi*a*b   -T1*b
      -pi*a*b   pi*(b^2)*((1/8)+(a^2))    -(T7+(c-a)*T1)*(b^2)
      -T1*b     2*T13*b^2           -(1/pi)*T3*(b^2)];

%% Aerodynamic Damping Matrix
Ca_Non_Circ = (b^2)*[0  2*pi     -2*T4
         0   2*pi*(0.5-a)*b     (T1-T8-(c-a)*T4 + (T11/2))*2*b                                % Multiply by (0.5*rho*V^2)/V
         0   ((-2*T9)-T1+T4*(a-0.5))*2*b       -T4*T11*b/pi];

Ca_Circ = [4*pi*b  4*pi*b*b*(0.5-a)     2*b*b*T11
         -4*pi*b*b*(0.5+a)   -4*pi*b*b*b*(0.5-a)*(0.5+a)     -4*b*b*b*(0.5+a)*T11/2           % Multiply by (0.5*rho*V^2)/V
         2*b*b*T12   2*b*b*T12*b*(0.5-a)       2*b*b*T12*b*T11/pi];

Ca = Ca_Non_Circ + Wag_0.*Ca_Circ;

%% Aerodynamic Stiffness Matrix
Ka_Non_Circ = -(b^2)*[0    0    0
       0    0   2*(T4+T10)                                                                   % Multiply by (0.5*rho*V^2)  
       0    0   2*(T5-T4*T10)/pi];

Ka_Circ = [0  -4*pi*b     -4*b*T10
         0   4*pi*b*b*(0.5+a)     4*b*b*(0.5+a)*T10                                           % Multiply by (0.5*rho*V^2) 
         0   -2*b*b*T12       -2*b*b*T12*T10/pi];  

Ka_lag = [4*pi*(-psi1*eps1-psi2*eps2)  2*pi*b*(2*((-psi1*eps1-psi2*eps2)*(0.5-a)))    2*b*((-psi1*eps1-psi2*eps2)*T11)
         -4*pi*b*(0.5+a)*(-psi1*eps1-psi2*eps2)   -4*pi*b*b*(0.5+a)*(-psi1*eps1-psi2*eps2)*(0.5-a)     -2*b*b*(0.5+a)*((-psi1*eps1-psi2*eps2)*T11)         % Multiply by (0.5*rho*V^2) 
         2*b*T12*(-psi1*eps1-psi2*eps2)   b*b*T12*2*(0.5-a)*(-psi1*eps1-psi2*eps2)       b*b*T12*(-psi1*eps1-psi2*eps2)*T11/pi];

Ka = Ka_Non_Circ + Wag_0*Ka_Circ + Ka_lag;
%% Lag State Coeff Matrix

W = [-2*pi*b*psi1*(eps1/b)^2  -2*pi*b*psi2*(eps2/b)^2    2*pi*b*psi1*eps1*(1-eps1*(0.5-a))/b     2*pi*b*psi2*eps2*(1-eps2*(0.5-a))/b    2*pi*b*psi1*eps1*(T10 - eps1*T11/2)/(pi*b)    2*pi*b*psi2*eps2*(T10 - eps2*T11/2)/(pi*b)
    2*pi*b*b*(0.5+a)*psi1*(eps1/b)^2    2*pi*b*b*(0.5+a)*psi2*(eps2/b)^2    -2*pi*b*b*(0.5+a)*psi1*eps1*(1-eps1*(0.5-a))/b  -2*pi*b*b*(0.5+a)*psi2*eps2*(1-eps2*(0.5-a))/b  -2*pi*b*b*(0.5+a)*psi1*eps1*(T10-eps1*T11/2)/(pi*b)  -2*pi*b*b*(0.5+a)*psi2*eps2*(T10-eps2*T11/2)/(pi*b)
         -b*b*T12*psi1*(eps1/b)*(eps1/b)   -b*b*T12*psi2*(eps2/b)*(eps2/b)       b*b*T12*psi1*eps1*(1-eps1*(0.5-a))/b   b*b*T12*psi2*eps2*(1-eps2*(0.5-a))/b    b*b*T12*psi1*eps1*(T10-eps1*T11/2)/(pi*b)   b*b*T12*psi2*eps2*(T10-eps2*T11/2)/(pi*b)];



