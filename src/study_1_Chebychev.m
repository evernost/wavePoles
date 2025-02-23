% =============================================================================
% Project       : wavePoles
% Module name   : study_1_Chebychev
% File name     : study_1_Chebychev.m
% File type     : Matlab script
% Purpose       : TODO
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

nTerms = 20;

nTries = 5000;
rList = linspace(0.52, 0.38, nTries).';

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

hFig = figure;
set(hFig, 'Visible', 'off');

hAx = axes('Parent', hFig);

writerObj = VideoWriter('roots.avi', 'Uncompressed AVI');
writerObj.FrameRate = 60;
open(writerObj);

lastProgress = 0;

for k = 1:nTries
    r = rList(k);
    
    
    n = (1:nTerms)';
    s = ones(nTerms,1); s(2:2:end) = -1;
    a = -s.*sin(2*pi*n*r)./(n*pi*sqrt(r*(1-r)));
    b = -s.*(1 - cos(2*pi*n*r))./(n*pi*sqrt(r*(1-r)));
    
    pCos = ([0; a]).' * T;
    pSin = ([b; 0]).' * U;
    
    x = linspace(-1.02, 1.02, 1000);
    xCos = cos(2*pi*x);
    xSin = sin(2*pi*x);
    
    rCos = roots(fliplr(pCos));
    rSin = roots(fliplr(pSin));
    plot(real(rCos(2:end)), imag(rCos(2:end)), '+', real(rSin(2:end)), imag(rSin(2:end)), '+')
    xlim([-1.3, 1.3])
    ylim([-0.5, 0.5])
    grid on
    title(sprintf('r = %0.3f', r))

    % Capture current frame
    frame = getframe(gcf);
    
    % Write frame to video file
    writeVideo(writerObj, frame);

    if round(100*k/nTries) >= (lastProgress + 5)
        lastProgress = lastProgress + 5;
        fprintf('%0.1f %%\n', lastProgress);
    end
    %pause(0.001)
end


close(writerObj);

