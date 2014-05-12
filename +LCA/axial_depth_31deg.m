function r = axial_depth_31deg(x) 
% 
% Gives axial depth given Sutter X input, assuming Sutter 4th axis setting
% of 31 degrees.
%
% DHO, 5/08.
%

r = sqrt(x.^2 + (1.6667*x).^2); % for 31 deg Sutter setting and pipette angle.