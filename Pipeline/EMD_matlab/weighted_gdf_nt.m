function [E,L2] = weighted_gdf_nt(V1, V2, mw1, mw2, chan_pos, dim_mask, l2_weight, xStep, zStep)
%
% GDF   Ground distance between two vectors
%    [E] = GDF(F1, F2) is the ground distance between two feature vectors.
%
%    Example:
%    -------
%        v1 = [100, 40, 22];
%        v2 = [50, 100, 80];
%        ...
%        [e] = gdf(v1, v2);
%        ...
%
%    This file and its content belong to Ulas Yilmaz.
%    You are welcome to use it for non-commercial purposes, such as
%    student projects, research and personal interest. However,
%    you are not allowed to use it for commercial purposes, without
%    an explicit written and signed license agreement with Ulas Yilmaz.
%    Berlin University of Technology, Germany 2006.
%    http://www.cv.tu-berlin.de/~ulas/RaRF
%
% Modified from Ulas original by adding weights for each dimension and a 
% mask for which dimensions are included. 
% V1, V2 = [ndim,1], doubles, values for the two vectors

% hard coding these in the code for now, to be able to use original emd
% code and only change the distance function.
% dim_mask = [ndim,1], logical, to include or not
% dim_weights = weight for each, generally the standard deviation

% further modified to allow calling l2 distance between 2D waveforms from
% this function. 
% TODO -- change weights from hard coded here to passing in as another
% variable; can then dispense with the mask, just send in 0 weights.

dim_weights(1) = 0.1; % centroid x, 1/um^2
dim_weights(2) = 0.1; % centroid z, 1/um^2
dim_weights(3) = 1; % y from d1d2 data, 1/um^2
dim_weights(4) = 0.0025; % peak to peak amplitude, 1/uV^2
dim_weights(5) = 1000; % duration = peak-to-trough time, 1/msec^2
dim_weights(6) = 20000; % fwhm, 1/msec^2
dim_weights(7) = 7000; % peak to trough amplitude ratio, unitless
dim_weights(8) = 0.25; % pre peak amplitude, 1/uV^2
dim_weights(9) = 0.02; % vertical footprint 1/um^2
dim_weights(10) = l2_weight; % 'L2 norm' (estimated, need to run some bootstrapping)

% V1 = F1(i, 1:a);
% V2 = F2(j, 1:a);
% mw1 = squeeze(mw1(i,:,:));
% mw2 = squeeze(mw2(j,:,:));

% calculate diffZ and diffX before applying mask
diffX = abs(V2(1) - V1(1));
diffZ = abs(V2(2) - V1(2));

V1 = V1(dim_mask(1:9));
V2 = V2(dim_mask(1:9));

w = dim_weights(dim_mask(1:9));
diffsq = (V2 - V1).^2;

z_lim = 500;
x_lim = 100;

if dim_mask(10) ~= 0
    % the L2 calculation compares waveforms centered at peak z for each
    % unit, to allow for moderate drift. Very large drift in z or x is physically
    % unreasonable, so only calculate real L2 for units within limits
    if (diffX < z_lim && diffZ < x_lim)
        % calculate l2 distance between these waveforms
        [~,wave_l2] = calcCenteredCorrL2(mw1,mw2, chan_pos, zStep, 5); % last param = nrow
    else
        wave_l2 = 2; % Largest possible value of normalized L2 diff
    end
    diffsq = [diffsq,wave_l2^2];
    w = [w,dim_weights(10)]; %weights for all parameters
    L2 = wave_l2; % if calculating return L2 for diagnostics
else
    L2 = 0;
end

E = sqrt(dot(diffsq,w));

end