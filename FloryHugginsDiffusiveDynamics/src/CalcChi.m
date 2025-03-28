function [chipe, chipw, chiew] = CalcChi(pH)
%% return pH-dependent chi-parameters

a = 1003.08;
b = -0.549856;
c = -0.5;

rhop = 0.0738673 - 0.0388822/(1 + exp(-2.19797*(pH -8.786)));
rhoe = 0.021369 + 0.030587/(1 + exp(-2.20309*(pH - 10.1761)));

chipw = -(a*rhop^2 + b);
chiew = -(a*rhoe^2 + b);
chipe = 2*(-a*rhop*rhoe + c) + chipw + chiew; 