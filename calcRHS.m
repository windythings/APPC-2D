function [RHS] = calcRHS(surfaces,wake)
% calcRHS Calculates the RHS vector during the calculation of airfoil
% circulation to account for the freestream velocity and the induced 
% velocities due to the wake and far-field vortex sheets.
% Inputs:  surfaces - structure containing all airfoil surfaces and details
%          wake - structure containing all wake surfaces and details
% Outputs: RHS - augmented RHS vector accounting for freestream and wake 
%                induced velocities on the airfoil


RHS = zeros(sum([surfaces.n]),1);
rowOffset = 0;

for p = 1:length(surfaces)
    indexList = (1:surfaces(p).m)+rowOffset;
    
    % The normalized freestream component normal to a panel is just 
    % sin(theta), where theta is the panel angle relative to the 
    % freestream. This step is standard for most panel codes.
    RHS(indexList) = RHS(indexList) + sin(surfaces(p).theta);
    
    for q = 1:length(wake)-1
        % Calculate the induction on the airfoil surfaces due to the 
        % farfield circulation shed downstream of the system and sum this 
        % with the freestream induction component.
        RHS(indexList) = RHS(indexList) + calcFarFieldInduction(surfaces(p),wake(q));
        
        % Calculate the influence of the wake boundaries themselves on the 
        % RHS of the airfoil surface circulation system. Per definition, 
        % these are subtracted from RHS as they are defined relative to the
        % normal vector inside the airfoil surfaces.
        [A_airfoil_wake,~] = calcVelMatricesFast(surfaces(p),wake(q));
        RHS(indexList) = RHS(indexList) - A_airfoil_wake*wake(q).gamma;
        
        % Augment the Kutta Condition at the trailing edge of the airfoil
        % surfaces such that the airfoil circulation matches the wake
        % circulation, ensuring smooth circulation distributions.
        if p == q
            RHS(indexList(end)+1) = wake(q).gamma(1);
        end
    end
    
    rowOffset = rowOffset + surfaces(p).n;
end

end

