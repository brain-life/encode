% This function compute the weighted Bundle-based Minimum Distance (wBMD)
% between two bundles. The wBMD is a modified BMD [Garyfallidis et al,
% 2015] that takes into account the fiber weights assigned in LiFE.
function [ distance ] = compute_wBMD(bA, bB, wA, wB, Npoints)
% bA: Bundle A (set of fibers)
% bB: Bundle B (set of fibers)
% wA: weights for Bundle A
% wB: weights for Bundle B
% Npoints: Number of points used to interpolate fibers ([Garyfallidis et al, 2015] uses Npoints=20)
% distance: Bundle-based Minimum Distance (BMD) computed as in eq. (1) in [Garyfallidis et al, 2015]

I = size(bA,1);
J = size(bB,1);

% Interpolate bundles
bA = interpolate(bA,Npoints);
bB = interpolate(bB,Npoints);

vA = zeros(I,1);
vB = zeros(J,1);

% normalize weights
wA = wA./sum(wA);
wB = wB./sum(wB);

% fix one fiber in Bundle A and compute the distances to every fibers in
% Bundle B
for i=1:I
    fA = bA{i};
    parfor j=1:J
        vB(j) = MDF(fA,bB{j});
    end
    vA(i) = min(vB);
end
dA = sum(wA.*vA);

% fix one fiber in Bundle B and compute the distances to every fibers in
% Bundle A
for j=1:J
    fB = bB{j};
    parfor i=1:I
        vA(i) = MDF(bA{i},fB);
    end
    vB(j) = min(vA);
end
dB = sum(wB.*vB);

distance = (dA + dB)/4;

end

% Function that provides an interpolated bundle with a constant N number of
% points each fiber
function [b] = interpolate(b,N)

parfor f=1:size(b,1)
    % Interpolate using Splines Tensor Toolbox
    CS=cscvn(b{f}); % obtain spline model from points
    t_range = linspace(CS.breaks(1),CS.breaks(end),N); % define range of t parameter in order to provide N points
    b{f} = fnval(CS,t_range); % Evaluate curve at N points
end

end

% Function that compute the Minimum average Direct-Flip (MDF) as defined in
% eq. (2) in [Garyfallidis et al, 2015]
function [dist] = MDF(f1,f2)
N = size(f1,2);
d_dir = mean(sqrt(sum((f1'-f2').^2,2)));
d_flip = mean(sqrt(sum((f1'-flip(f2',1)).^2,2)));

dist = min(d_dir,d_flip);

end

