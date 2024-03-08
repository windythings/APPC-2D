function [surfaces] = calcAirfoilGamma(surfaces,A,RHS)
% calcAirfoilGamma Calculates the circulation distribution around the
% airfoil elements based on the aerodynamic influence matrix and panelling
% RHS vector.
% 
% Inputs:  surfaces - the structure of surface elements with all details
%          A - the aerodynamic influence coefficient matrix
%          RHS - the V-infinity-normal component vector
% Outputs: surfaces - the structure of surface elements updated with the
%          new values of circulation around the airfoil

% Invert the aerodynamic influence coefficient matrix and multiply by the
% RHS matrix to solve for airfoil circulation distributions.
gammaAirfoil = A\RHS;
    
% gammaAirfoil contains the gammas for all surfaces. This for loop splits
% the vector into parts to be placed in each surfaces's data structure cell
temp = 1;
for i = 1:length(surfaces)
    surfaces(i).gamma = gammaAirfoil(temp:(temp-1)+length(surfaces(i).endPoints(:,1)));
    temp = temp + length(surfaces(i).endPoints(:,1));
end

end

