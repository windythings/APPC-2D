function [wake] = panelWake(wake)
% panelWake Assigns the panel-defining points for the wake for both the
% upper and lower wake surfaces, based on the spacing parameter, starting
% at the airfoil chord and ending at a distance chord*chordMultiplier. Also
% assigns the collocation points and calculates local panel angles.
% 
% Inputs:  wake - a structure ready to accept the new wake parameters
%          (panel points, collocation points, and local angles)
% Outputs: wake - a structure containing the wake parameters for both the
%          upper and lower wake surfaces

for i = 1:length(wake)
    % Populate airfoil struct with all details for paneling:
    %         n - number of panel nodes
    %         m - number of panel control points
    %       pt1 - panel nodes set one
    %       pt2 - panel nodes set two, staggered from set one
    %        co - panel control points, midpoint between nodes
    %     theta - panel angle relative to freestream
    wake(i).n = length(wake(i).endPoints(:,1));
    wake(i).m = wake(i).n-1;
    wake(i).pt1 = wake(i).endPoints(1:wake(i).m,:);
    wake(i).pt2 = wake(i).endPoints(2:wake(i).n,:);
    wake(i).co = wake(i).pt1 + ((wake(i).pt2-wake(i).pt1)./2);
    wake(i).theta = atan2((wake(i).pt2(:,2)-wake(i).pt1(:,2)),...
                      (wake(i).pt2(:,1)-wake(i).pt1(:,1)));            
end

end

