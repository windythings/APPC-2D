function [airfoil,chord] = panelAirfoils(surfaceFiles,alpha,mainSurfIndex)
% panelAirfoils Assigns the panel-defining points around the airfoil
% surfaces imported, and calculates the local panel angle to the airfoil
% chord line, collocation point coordinates, and the number of control
% points. Concatenates all of the respective values for each airfoil
% surface requested and returns a structure.
% 
% Inputs:  surfaceFiles - a structure of file names storing coordinates
%          alpha - the angle of attack of the lifting system to rotate by
%          mainSurfIndex - index of main airfoil surface
% Outputs: airfoil - a structure containing the airfoil parameters for all
%          surfaces together in a single location
%          chord - main airfoil surface chord length

% Create airfoil surfaces struct and rotation matrix for AOA
airfoil = struct('name', cell(1,length(surfaceFiles)));
R = [cos(-alpha) -sin(-alpha); sin(-alpha) cos(-alpha)];

for i = 1:length(surfaceFiles)
    coords = importdata(surfaceFiles(i).name);
    
    % Calculate the chord based on the main airfoil
    if i == mainSurfIndex
        chord = abs(max(coords(:,1))); 
    end
    
    % Rotate all airfoil surfaces by alpha relative to the quarter chord
    center = repmat([0.25;0],1,length(coords(:,1)));
    coords = (R*(coords'-center) + center)';
    
    % Populate airfoil struct with all details for paneling:
    % name - airfoil surface name
    % endPoints - panel nodes
    %         n - number of panel nodes
    %         m - number of panel control points
    %       pt1 - panel nodes set one
    %       pt2 - panel nodes set two, staggered from set one
    %        co - panel control points, midpoint between nodes
    %     theta - panel angle relative to freestream
    airfoil(i).name = "Airfoil " + num2str(i);
    airfoil(i).endPoints = coords;
    airfoil(i).n = length(coords);
    airfoil(i).m = airfoil(i).n - 1;
    airfoil(i).pt1 = airfoil(i).endPoints(1:airfoil(i).m,:);
    airfoil(i).pt2 = airfoil(i).endPoints(2:airfoil(i).n,:);
    airfoil(i).co = airfoil(i).pt1 + ((airfoil(i).pt2-airfoil(i).pt1)./2);
    airfoil(i).theta = atan2((airfoil(i).pt2(:,2)-airfoil(i).pt1(:,2)),...
                          (airfoil(i).pt2(:,1)-airfoil(i).pt1(:,1)));
end

end

