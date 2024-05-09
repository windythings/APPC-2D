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
colOffset = 0;
rowOffset = 0;

A = zeros(totalMatrixSize,totalMatrixSize);
B = zeros(totalMatrixSize,totalMatrixSize);

% For loop calculates the influence matrices for each airfoil surface's
% influence on itself and the other(s). Concatenates these in the 
% appropriate section of the large A and B matrices.
for p = 1:length(surfaces)
    
    indexList = (1:surfaces(p).m)+rowOffset;
    
    for q = 1:length(surfaces)
        
        colIndexList = (1:surfaces(q).n)+colOffset;
        
        [A(indexList,colIndexList),B(indexList,colIndexList)] = ...
                            calcVelMatricesFast(surfaces(p),surfaces(q));
        
        % Enforcing Kutta condition
        if p == q
            A(indexList(end)+1,colIndexList(1)) = 1;            
            A(indexList(end)+1,colIndexList(end)) = 1;
        end
        
        colOffset = colOffset + surfaces(q).n;
    end
    
    colOffset = 0;
    rowOffset = rowOffset + surfaces(p).n;
end

end