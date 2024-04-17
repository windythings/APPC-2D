 function [surfaces,CL] = calcCoefficients(surfaces,wake)
% calcCoefficients Calculates the lift coefficient of the system and the
% pressure coefficients along the airfoil surfaces independently.
% 
% Inputs:  surfaces - the structure containing all the details about the
%          airfoil surfaces
%          wake - the structure containing all the details about the wake
% Outputs: surfaces - the structure containing the airfoil surface details
%          modified with the pressure coefficients
%          CL - the system lift coefficient

CL = 0;

% Structure of this loop is similar to that in calcWakeVels.m, where the
% influence of each potential flow component on each airfoil surface is
% calculated and the tangential induced velocities are summed.
for k = 1:length(surfaces)
    
    % Influence of the airfoil surfaces on each other
    for i = 1:length(surfaces)
        [~,B_surf_surf] = calcVelMatricesFast(surfaces(k),surfaces(i));
        tempTang = B_surf_surf*surfaces(i).gamma;
        
        if i == 1
            vTang = zeros(length(tempTang),1);
        end
        
        vTang = vTang + tempTang;
    end
    
    % Influence of the wake boundaries
    for i = 1:length(wake)-1
        [~,B_surf_wake] = calcVelMatricesFast(surfaces(k),wake(i));
        vTang = vTang + B_surf_wake*wake(i).gamma;
    end
    
    % Influence of the freestream
    vTang = vTang + cos(surfaces(k).theta);
    
    % Influence of the farfield circulation solution
    for i = 1:length(wake)-1
        [~,tempTang] = calcFarFieldInduction(surfaces(k),wake(i));
        vTang = vTang + tempTang;
    end
    
    % Calculate pressure coefficient on each panel using Bernoulli's eq.
    surfaces(k).CP = 1 - (vTang.^2);
    
    % Lift coefficient in 2D is the the panel length multiplied by the
    % change in circulation across the panel.
    surfaces(k).DL = (surfaces(k).pt2(:,1)-surfaces(k).pt1(:,1)).*cos(surfaces(k).theta) + ...
         (surfaces(k).pt2(:,2)-surfaces(k).pt1(:,2)).*sin(surfaces(k).theta);
    temp = surfaces(k).gamma(1:end-1) + surfaces(k).gamma(2:end);
    CL = CL + temp'*surfaces(k).DL;
end

end

