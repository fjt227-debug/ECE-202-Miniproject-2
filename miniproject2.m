clear all; 
close all;

a = 1; %side a is 1m long
b = 2; %side b is 2m long
V1 = 1; %top plate potential of 1V
V2 = -1; %bottom plate potential of -1V
EPS0 = 8.8542*10^(-12);
n = 10; %number of patches per side a
m = 20; %added m = 20 

dp = a/n; %size of individual patch
N = n*m; 
%Patch centers and areas
x = zeros(1,N); 
y = zeros(1,N); 
dS = zeros(1,N);
for i = 1:n
    for j = 1:m 
        x(m*(i-1)+j) = (i-1/2)*dp;
        y(m*(i-1)+j) = (j-1/2)*dp;
        dS(m*(i-1)+j) = dp^2;
    end
end
%Maps charges to potentials
Gself = zeros(N,N);
for i = 1 : N
    for j = 1 : N
        r = sqrt( (x(j)-x(i))^2 + (y(j)-y(i))^2 );
        if (i==j)
            Gself(i,j) = sqrt(dS(j))/(2*sqrt(pi)*EPS0);
        else
            Gself(i,j) = dS(j)/(4*pi*EPS0*r);   
        end
    end
end
%Required seperations
ratios = [0.1 0.5 1 2 10];
Qplus  = zeros(size(ratios));
Ccomp  = zeros(size(ratios));
Canal  = zeros(size(ratios));

%Loop over seperations
for t = 1:numel(ratios)
    d = ratios(t)*a;
    Gcross = zeros(N,N);
    for i = 1:N %Cross plate matrix
        for j = 1:N
            r = sqrt((x(j)-x(i))^2 + (y(j)-y(i))^2 + d^2);
            Gcross(i,j) = dS(j)/(4*pi*EPS0*r);
        end
    end
    A = Gself - Gcross;
    B = V1 * ones(N,1);
    rhos = A \ B; %Surface charge density on the positive plate
    %Totals and Capacitances
    Qtot = dS * rhos;
    Qplus(t) = Qtot;
    Ccomp(t) = Qplus(t) / (V1 - V2); 
    Canal(t) = EPS0 * (a*b) / d;
end    

figure(1);
plot(ratios, Qplus*1e12, 'o-','LineWidth',1.6);
grid on;
xlabel('d / a'); 
ylabel('Q Computational on and Plate [pC]');
title('Total Charge vs Separation');


figure(2);
plot(ratios, Canal./Ccomp, 's-','LineWidth',1.6); 
grid on; 
hold on; 
yline(1,'--');
xlabel('d / a'); 
ylabel('Canalytical / Ccomputational');
title('Capacitance Ratio vs Separation');

disp('Results Miniproject#2:'); %displaying output
for k = 1:numel(ratios)
    fprintf('d/a=%g, Q+=%.4e C, Ccomp=%.4e F, Canal=%.4e F, ratio=%.4f\n',ratios(k), Qplus(k), Ccomp(k), Canal(k), Canal(k)/Ccomp(k));
end
