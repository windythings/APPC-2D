function [A,B] = calcAirfoilMatrices(surfaces)
% calcAirfoilMatrices Calculates the airfoil influence coefficient matrix
% for an arbitrary number of surfaces defined in the input surfaces struct.
% 
% Inputs:  surfaces - a structure containing details of all the panelled
%          airfoil surfaces of interest
% Outputs: A - the influence coefficient matrix used to solve for
%          circulation around the airfoil surfaces
%          B - the influence coefficient matrix used to calculate the
%          induced velocity of a panel at another location in the domain

influence = struct('name', cell(2,2));
totalMatrixSize = 0;

% For loop calculates the influence matrices for each airfoil surface's
% influence on itself and the other(s). Stores all of those influence
% coefficient matrices in the "influence" data structure.
for p = 1:length(surfaces)
    totalMatrixSize = totalMatrixSize + surfaces(p).n;
    
    for q = 1:length(surfaces)
        influence(p,q).name = "On " + num2str(p) + " Due To " + num2str(q);
        [influence(p,q).A,influence(p,q).B] = ...
                                  calcVelMatrices(surfaces(p),surfaces(q));
    end
end

Atemp = [];
Btemp = [];
A = [];
B = [];

% Enforcing Kutta Condition for the appropriate influence matrices (those
% of a surface's influence on itself).
for p = 1:length(surfaces)
    for q = 1:length(surfaces)
        
        influence(p,q).A(end+1,:) = zeros(1,length(influence(p,q).A(end,:)));
        influence(p,q).B(end+1,:) = zeros(1,length(influence(p,q).B(end,:)));
        
        % Enforcing Kutta Condition at trailing edge elements
        if p==q
            influence(p,q).A(end,1) = 1;
            influence(p,q).A(end,end) = 1;
        end
    end
end

% Concatenates all of the influence matrices into a giant A matrix and B
% matrix, plus an RHS vector. This contains all of the influence
% coefficients and can be inverted to find the airfoil gammas.
for p = 1:length(surfaces)
    for q = 1:length(surfaces)
        Atemp = [Atemp influence(p,q).A];
        Btemp = [Btemp influence(p,q).B];
    end
    
    A = [A; Atemp];
    B = [B; Btemp];
    Atemp = [];
    Btemp = [];
end

end