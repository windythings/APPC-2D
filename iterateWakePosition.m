function [wake] = iterateWakePosition(wake)
% iterateWakePosition Uses the perturbation velocities of airfoil and wake
% elements on the wake to calculate the new y-position of all wake panels.
% 
% Inputs:  wake - the structure containing the wake element details
% Outputs: wake - the wake element structure data, with new endPoint values

for i = 1:length(wake)-1
    % Rotating the tangential and normal velocities induced on the current
    % wake element by the local panel angle to get induced velocities in
    % the inertial frame
    wake(i).vX = wake(i).vTang.*cos(wake(i).theta) - ...
                 wake(i).vNorm.*sin(wake(i).theta);
    wake(i).vY = wake(i).vTang.*sin(wake(i).theta) + ...
                 wake(i).vNorm.*cos(wake(i).theta);
             
    % Shollenberger's algorithm for iterating the wake's y position
    toAdd = 0;
    
    for j = 2:length(wake(i).co)
        for k = 2:j
            toAdd = toAdd + ((wake(i).vY(k)./wake(i).vX(k)).*...
                             (wake(i).co(k,1)-wake(i).co(k-1,1)));
        end
        
        wake(i).co(j,2) = wake(i).co(1,2) + toAdd;
        toAdd = 0;

    end
    
    % Recalculating panel endpoints from control points, extrapolate the
    % last panel endpoint
    co1 = wake(i).co(1:end-1,2);
    co2 = wake(i).co(2:end,2);
    wake(i).endPoints(2:end-1,2) = (co1+co2)./2;
    wake(i).endPoints(end,2) = wake(i).endPoints(end-1,2) + ...
                  (wake(i).endPoints(end-1,2)-wake(i).endPoints(end-2,2));
end

% Recalculate the wake centerline position
wake(i+1).endPoints = (wake(i-1).endPoints+wake(i).endPoints)./2;

end

