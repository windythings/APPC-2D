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

totalMatrixSize = sum([surfaces.n]);
colOffset = 1;
rowOffset = 1;

A = zeros(totalMatrixSize,totalMatrixSize);
B = zeros(totalMatrixSize,totalMatrixSize);

% For loop calculates the influence matrices for each airfoil surface's
% influence on itself and the other(s). Concatenates these in the 
% appropriate section of the large A and B matrices.
for p = 1:length(surfaces)
    for q = 1:length(surfaces)
        [A(rowOffset:rowOffset+surfaces(p).m-1,colOffset:colOffset+surfaces(q).m),...
         B(rowOffset:rowOffset+surfaces(p).m-1,colOffset:colOffset+surfaces(q).m)] = ...
                            calcVelMatricesFast(surfaces(p),surfaces(q));
        
        % Enforcing Kutta condition
        if p == q
            A(rowOffset+surfaces(p).m,colOffset) = 1;            
            A(rowOffset+surfaces(p).m,colOffset+surfaces(q).m) = 1;
        end
        
        colOffset = colOffset + surfaces(q).n;
    end
    colOffset = 1;
    rowOffset = rowOffset + surfaces(p).n;
end

end