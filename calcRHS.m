function [RHS] = calcRHS(surfaces,wake)
% calcRHS Calculates the RHS vector during the calculation of airfoil
% circulation to account for the freestream velocity and the induced 
% velocities due to the wake and far-field vortex sheets.
% Inputs:  surfaces - structure containing all airfoil surfaces and details
%          wake - structure containing all wake surfaces and details
% Outputs: RHS - augmented RHS vector accounting for freestream and wake 
%                induced velocities on the airfoil

RHS = [];

% The normalized freestream component normal to a panel is just sin(theta),
% where theta is the panel angle relative to the freestream. This step is
% standard for most panel codes.
for p = 1:length(surfaces)
    for i = 1:length(surfaces(p).theta)
         RHStemp(i,1) = sin(surfaces(p).theta(i));
    end
    RHStemp(end+1,1) = 0; % enforce zero Kutta condition for trailing edge elements
    
    RHS = [RHS; RHStemp];
    RHStemp = [];
end

% Calculate the induction on the airfoil surfaces due to the farfield
% circulation shed downstream of the system and sum this with the
% freestream induction component.
vInfNorm_airfoil_farField = [];

for j = 1:length(wake)-1
    for i = 1:length(surfaces)
        [tempNorm,~] = calcFarFieldInduction(surfaces(i),wake(j));
        tempNorm(end+1,1) = 0; % enforce zero Kutta condition for trailing edge elements
        vInfNorm_airfoil_farField = [vInfNorm_airfoil_farField; tempNorm];
    end
    
    RHS = RHS + vInfNorm_airfoil_farField;
    vInfNorm_airfoil_farField = [];
end

% Calculate the influence of the wake boundaries themselves on the RHS of
% the airfoil surface circulation system. Per definition, these are
% subtracted from RHS as they are defined relative to the normal vector
% inside the airfoil surfaces.
temp = [];
vInfNorm_airfoil_wake = [];

for j = 1:length(wake)-1
    for i = 1:length(surfaces)
        % Notation: Matrix Influence_on Airfoil Surface_due to Wake Surface
        [A_airfoil_wake,~] = calcVelMatrices(surfaces(i),wake(j));
        temp = A_airfoil_wake*wake(j).gamma;
        
        % Augment the Kutta Condition at the trailing edge of the airfoil
        % surfaces such that the airfoil circulation matches the wake
        % circulation, ensuring smooth circulation distributions.
        if j == 1
            if i == 1
%                 temp(end+1,1) = -wake(j).gamma(1);
                temp(end+1,1) = 0;
            else
                temp(end+1,1) = 0; 
            end
        elseif j == 2
            if i == 2
%                 temp(end+1,1) = -wake(j).gamma(1);
                temp(end+1,1) = 0;
            else
                temp(end+1,1) = 0;
            end
        else
            temp(end+1,1) = 0;
        end
        vInfNorm_airfoil_wake = [vInfNorm_airfoil_wake; temp];
    end
    
    RHS = RHS - vInfNorm_airfoil_wake;
    temp = [];
    vInfNorm_airfoil_wake = [];
end

end

