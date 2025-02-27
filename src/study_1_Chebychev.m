% =============================================================================
% Project       : wavePoles
% Module name   : study_1_Chebychev
% File name     : study_1_Chebychev.m
% File type     : Matlab script
% Purpose       : decompose a square wave onto Chebychev polynomials
% Author        : QuBi (nitrogenium@outlook.fr)
% Creation date : Sunday, 23 February 2025
% -----------------------------------------------------------------------------
% Best viewed with space indentation (2 spaces)
% =============================================================================

% -----------------------------------------------------------------------------
% DESCRIPTION
% -----------------------------------------------------------------------------
% TODO



close all
clear all
clc

% -----------------------------------------------------------------------------
% SETTINGS
% -----------------------------------------------------------------------------
nTerms = 4;
r = 0.9;      % PWM ratio


% -----------------------------------------------------------------------------
% BASIS FOR CHEBYCHEV POLYNOMIALS
% -----------------------------------------------------------------------------
T = zeros(nTerms+1);
T(1, 1)   = 1;
T(2, 1:2) = [0, 1];
for n = 3:(nTerms+1)
    T(n,:) = 2*[0, T(n-1, 1:(end-1))] - T(n-2, :);
end    

U = zeros(nTerms+1);
U(1, 1)   = 1;
U(2, 1:2) = [0, 2];
for n = 3:(nTerms+1)
    U(n,:) = 2*[0, U(n-1, 1:(end-1))] - U(n-2, :);
end



% -----------------------------------------------------------------------------
% ADDITIVE SYNTHESIS
% -----------------------------------------------------------------------------

% Build the weights for sin() and cos()
n = (1:nTerms)';
s = ones(nTerms,1); s(2:2:end) = -1;
%a = -s.*sin(2*pi*n*r)./(n*pi*sqrt(r*(1-r)));
%b = -s.*(1 - cos(2*pi*n*r))./(n*pi*sqrt(r*(1-r)));
v = sqrt(3); u = -sqrt(3);
a = s.*sqrt(3).*(cos(2*pi*n*r) - 1)./(n.*n*pi*pi*r*(1-r));
b = s.*sqrt(3).*sin(2*pi*n*r)./(n.*n*pi*pi*r*(1-r));

pCos = ([0; a]).' * T; pCos = fliplr(pCos);
pSin = ([b; 0]).' * U; pSin = fliplr(pSin);

% Apply the distorsion polynomial
x = linspace(-1.02, 1.02, 1000);
xCos = cos(2*pi*x);
xSin = sin(2*pi*x);
y = polyval(pCos, xCos) + xSin.*polyval(pSin, xCos);

% -----------------------------------------------------------------------------
% PLOTS: TIME SERIES + ZEROS REPRESENTATION
% -----------------------------------------------------------------------------

% Plot the time series (square wave in time)
figure
plot(x,y)
xlim([-0.75, 0.75])
grid minor

% Alternate representation of the polynomials: zeros plot
figure
rCos = roots(pCos);
rSin = roots(pSin);
plot(real(rCos(2:end)), imag(rCos(2:end)), '+', real(rSin(2:end)), imag(rSin(2:end)), '+')
legend('cos', 'sin')
xlim([-1.3, 1.3])
ylim([-0.5, 0.5])
grid on
title(sprintf('r = %0.3f', r))
axis equal;



