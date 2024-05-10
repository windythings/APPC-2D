function [tempNorm,tempTang] = calcFarFieldInduction(controls,wake)
% calcFarFieldInduction Calculates the induction effects of a far-field
% constant-strength vortex sheet, emanating from an input wake, on a set of
% control points. Implements induction formulation of constant-strength
% vortex sheets from Katz & Plotkin.
%
% Inputs:  controls - a set of control points to calculate induction
%          effects on from a far-field vortex sheet
%          wake - a wake structure that defines the start of the far-field
%          constant-strength vortex sheet
% Outputs: tempNorm - the normal velocities on the control points
%          tempTang - the tangential velocities on the control points

tempNorm = zeros(length(controls.co),1);
tempTang = zeros(length(controls.co),1);
 
% Create two horizontal lines of farfield vorticity starting just aft of 
% the last wake point and ending "very far away," in this case, x/c = 1000.
x1 = wake.endPoints(end,1)+0.0000001;
x2 = wake.endPoints(end,1)+999;
y1 = wake.endPoints(end,2);

gammaInf = wake.gammaInf;

xc = controls.co(:,1);
yc = controls.co(:,2);

theta1 = atan2((yc-y1),(xc-x1));
theta2 = atan2((yc-y1),(xc-x2));

R2 = ((xc-x2).^2)+((yc-y1).^2);
R1 = ((xc-x1).^2)+((yc-y1).^2);

% Solving the potential flow equations for lines of constant vorticty
% yields the following relations for velocities they induce. tempNorm and
% tempTang are velocity components relative to panels.
tmpT = gammaInf./(2.*pi)*(theta2-theta1);
tmpN = gammaInf./(4.*pi)*log(R2./R1);

% Rotate each component vector into the local panel frame
tempTang = cos(controls.theta).*tmpT + sin(controls.theta).*tmpN;
tempNorm = -sin(controls.theta).*tmpT + cos(controls.theta).*tmpN;

end

