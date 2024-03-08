function [wake] = createWake(wakeStart,nWake,chord,chordMultiplier,chordsRamp,spacing)
% createWake Assigns the panel-defining points for the wake for both the
% upper and lower wake surfaces, based on the spacing parameter, starting
% at the airfoil chord and ending at a distance chord*chordMultiplier.
% 
% Inputs:  wakeStart - a vector of points where the upper and lower wake
%          should start from, respectively
%          nWake - number of wake endpoints to create
%          chord - the airfoil chord length
%          chordMultiplier - the number of chord lengths that the wake
%          should be extended back into the flow for
%          chordsRamp - the number of chord lengths after the initial wake
%          that the wake should extend for the circulation to be ramped
%          down to the freestream circulation
%          spacing - the panel spacing identifier  (uniform or cosine)
% Outputs: wake - a structure containing the wake parameters for both the
%          upper and lower wake surfaces

% Create wake structure
wake = struct('name', cell(1,2));

% Solves for the x/c distance at which the wake should terminate
wakeEnd = ceil(chord*(chordMultiplier+chordsRamp));

for i = 1:2
    % Uniform panel spacing throughout the wake
    if isequal(spacing,"uniform")
        coordsX = linspace(wakeStart(i,1),wakeEnd,nWake)';
    
    % Concentrates panel density near the surface trailing edges for better
    % wake shape fidelity at the region of highest curvature
    elseif isequal(spacing,"cosine")
        theta = linspace(0,0.5*pi,nWake);
        coordsX = ((wakeEnd-wakeStart(i,1))*...
                        (1-cos(theta))) + wakeStart(i,1);
        coordsX = coordsX';
    end
    
    % The wake is always created parallel to the freestream flow
    coordsY = wakeStart(i,2)*ones(nWake,1);
    
    wake(i).endPoints = [coordsX,coordsY];
    
    if i == 1
        wake(i).name = "Upper";
    else
        wake(i).name = "Lower";
        wake(i+1).endPoints = (wake(i-1).endPoints+wake(i).endPoints)./2;
        wake(i+1).name = "Center";
    end
end

end