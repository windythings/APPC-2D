function [wake] = wakeCirculation(wake,deltaH_ND,chord,chordMultiplier,wakeStart,omega,gammaPrevious,gammaInf)
% wakeCirculation Calculates the wake circulation distribution.
% 
% Inputs:  wake - a structure containing the wake parameters that will be
%          modified to include the gamma distribution for each boundary
%          deltaH_ND - a nondimensional value of the jet pressure jump
%          chord - the airfoil chord length
%          chordMultiplier - the number of chord lengths that the wake
%          should be extended back into the flow for
%          wakeStart - a vector of points where the upper and lower wake
%          should start from, respectively
%          omega - relaxation factor
%          gammaPrevious - previous value of wake circulation
%          gammaInf - freestream gamma
% Outputs: wake - the modified wake struct that includes new circulations

% For both wake boundaries, the wake circulation is calculated as per the
% mathematical formulation and relaxed to alleviate numerical instability.
for i = 1:length(wake)-1
    gamma = ((-1)^(i)).*(deltaH_ND./(wake(i).uBar));
    
    % Wake average velocities are calculated at panel control points,
    % whereas gamma values are required at panel nodes. Therefore, a shift
    % must be done to calculate the gamma values at the first and last wake
    % panel nodes, based on a linear extrapolation.
    gamma = [gamma(1)+(gamma(1)-gamma(2)); 
             gamma; 
             gamma(end)+(gamma(end)-gamma(end-1))];
    gamma1 = gamma(1:end-1);
    gamma2 = gamma(2:end);
    gamma = (gamma1 + gamma2)./2; 
    
    % Wake circulation relaxation
    if i == 1
        gamma = omega.*gamma + (1-omega).*gammaPrevious(1).upper;
    else
        gamma = omega.*gamma + (1-omega).*gammaPrevious(1).lower;
    end
    
    xN = chord*chordMultiplier + wakeStart(i,1);
    xN_Index = find(wake(i).endPoints(:,1) > xN);
    xN_Index = xN_Index(1);
    
    % From xN_Index back to the end of the wake, the circulation is
    % linearly ramped from the final gamma value calculated to the
    % circulation value manifest in the farfield, gammaInf.
    slope = (((((-1)^i).*gammaInf)-gamma(xN_Index))./...
                                        (wake(i).endPoints(end,1)-xN));
    gamma(xN_Index:end,1) = wake(i).endPoints(xN_Index:end,1)*slope + ...
                                            (gamma(xN_Index)-(slope*xN));
 
    % Adding relavant parameters to the wake struct
    wake(i).gamma = gamma;
    wake(i).chord = chord;
    wake(i).chordMultiplier = chordMultiplier;
    wake(i).gammaInf = ((-1).^i).*gammaInf;
end

end

