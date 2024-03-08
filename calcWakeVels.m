function [wake] = calcWakeVels(wake,surfaces)
% calcWakeVels Calculates the next iteration of wake average velocities
% bound to the individual wake elements, using induced velocities from the
% airfoil surfaces, wake surfaces, and the freestream velocity.
% 
% Inputs:  wake - the structure containing the wake element details
%          surfaces - the structure containing the airfoil element details
% Outputs: wake - the wake element structure data, with new values of
%          tangent and normal velocity to every wake panel

for j = 1:length(wake)
    
    tempNorm = [];
    tempTang = [];
    vInfNorm = [];
    vInfTang = [];
        
    % For each wake surface, calculate influence coefficients due to the
    % airfoil surfaces, multiply by airfoil circulation to get induced
    % velocities, and add them to the wake structure.
    for i = 1:length(surfaces)
        [A_wake_airfoil,B_wake_airfoil] = calcVelMatrices(wake(j),surfaces(i));
        
        tempNorm = A_wake_airfoil*surfaces(i).gamma;
        tempTang = B_wake_airfoil*surfaces(i).gamma;
        
        if i == 1
            wake(j).vNorm = zeros(length(tempNorm),1);
            wake(j).vTang = zeros(length(tempTang),1);
        end
        
        wake(j).vNorm = wake(j).vNorm + tempNorm;
        wake(j).vTang = wake(j).vTang + tempTang;
    end
    
    % For each wake surface, calculate influence coefficients due to the
    % wake surfaces on themselves and each other. Then multiply by wake
    % circulation to get induced velocities.
    for i = 1:length(wake)-1
        [A_wake_wake,B_wake_wake] = calcVelMatricesWake(wake(j),wake(i));

        tempNorm = A_wake_wake*wake(i).gamma;
        tempTang = B_wake_wake*wake(i).gamma;

        wake(j).vNorm = wake(j).vNorm + tempNorm;
        wake(j).vTang = wake(j).vTang + tempTang;
    end
    
    % Add induced velocities on wake by the freestream velocity. Negative
    % sign on normal component is necessary due to the definition of the
    % panel theta.
    vInfNorm = -sin(wake(j).theta);
    vInfTang = cos(wake(j).theta);
    
    wake(j).vNorm = wake(j).vNorm + vInfNorm;
    wake(j).vTang = wake(j).vTang + vInfTang;
    
    % Calculation of Velocities Induced by Far-Field Constant Vorticity Lines
    for i = 1:length(wake)-1
        [tempNorm,tempTang] = calcFarFieldInduction(wake(j),wake(i));
        
        wake(j).vNorm = wake(j).vNorm + tempNorm;
        wake(j).vTang = wake(j).vTang + tempTang;
    end
end

end

