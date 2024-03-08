function [wake] = calcUBar(wake)
% calcUBar Calculates the "average" velocity value at the wake boundaries.
% The "average" is treated as simply the speed of the flow at the wake
% boundary itself, which is enough for the algorithm to find a solution.
%
% Inputs:  wake - structure of wake elements
% Outputs: wake - structure of wake elements augmented with uBar cell

for i = 1:length(wake)-1
    uBarCo = sqrt((wake(i).vNorm.^2) + (wake(i).vTang.^2));
    wake(i).uBar = uBarCo;
end

end

