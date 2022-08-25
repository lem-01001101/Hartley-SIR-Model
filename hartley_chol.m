beta = 0.01;
gamma = 0.12;
N = 50;
I0 = 1;
b = 0;
%[S I R]

f = @(t,x) [-beta*x(1)*x(2) + b*x(3);
beta*x(1)*x(2) - gamma*x(2);
gamma*x(2) - b*x(3)];

[t,xa] = ode45(f, [0 42], [(N-I0) I0 0]); %6 weeks

plot(t,xa(:,1))
hold on
grid on
plot(t,xa(:,2))
plot(t,xa(:,3))
hold off
legend('S','I','R')
R0 = N*(beta/gamma);
txt = texlabel('R_0 = ');
text(40.5, 31.5, txt)
txt = num2str(R0);
text(42, 32, txt)

%Hartley et al.’s SIR model for Cholera Code:
betaL = 1.1;
betaH = 0.9;
gamma = (1/9);
KL = 10^6;
KH = KL/700;
X = (5/24);
xi = 0.45;
b = (1/12775);
deltaL = (1/30);
N = 500;

I0 = 1;
BL0 = 0;
BH0 = 0;

%[S I R BH BL]
f = @(t,x) [b*N - ( betaL*x(1)*(x(5)/(KL + x(5))) ) - ( betaH*x(1)*(x(4)/(KH + x(4))) ) - b*x(1);
( betaL*x(1)*(x(5)/(KL + x(5))) ) + ( betaH*x(1)*(x(4)/(KH + x(4))) ) - x(2)*(gamma + b);
gamma*x(2) - b*x(3);
xi*x(2) - X*x(4);
X*x(4) - deltaL*x(5)];

[t,xa] = ode45(f, [0 84], [(N-I0) I0 0 BH0 BL0]); %12 weeks

plot(t,xa(:,1))
grid on
hold on
plot(t,xa(:,2))
plot(t,xa(:,3))
hold off
xlabel('Days'); ylabel('Number of individuals');
legend('S','I','R')
R0 = ( (xi*N)/(gamma + b) ) * ( (betaH/(KH*X) ) + (betaL/(KL*deltaL)) );
txt = texlabel('R_0 = ');
text(82, 292, txt)
txt = num2str(R0);
text(85, 300, txt)
