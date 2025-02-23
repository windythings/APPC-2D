function [uGrid,vGrid] = createStreamlines(xGrid,yGrid,surfaces,wake,dval)
% createStreamlines Calculates the x and y velocities in the inertial frame
% necessary to implement MATLAB's streamline function in the main solver.
% Uses the Katz & Plotkin adapted algorithm (calcVelMatrices.m) to solve
% for induction effects of all elements with bound circulation on the grid.
%
% Inputs:  xGrid - a set of x values for all points in the inertial grid
%          yGrid - a set of y values for all points in the inertial grid
%          surfaces - structure of all airfoil elements
%          wake - structure of all wake elements
% Outputs: uGrid - x velocities on all the grid points defined
%          vGrid - y velocities on all the grid points defined

if nargin == 4; dval = pi; end

% Create and "panel" a structure containing all the surfaces and wakes to be
% passed into the K & P algorithm for calculations
grid.name = 'Grid';
grid.co = [xGrid(:) yGrid(:)];
%grid.theta = zeros(length(gridCo),1);
grid.theta = [];
grid.m = size(grid.co,1);

% Influence of airfoil surfaces on domain grid
U = 1; V = 0;
for j = 1:length(surfaces)
    [A_grid_surf,B_grid_surf] = calcVelMatricesFast(grid,surfaces(j));

    V = V + A_grid_surf*surfaces(j).gamma;
    U = U + B_grid_surf*surfaces(j).gamma;
end

% Influence of wake  boundaries on domain grid
for j = 1:length(wake)-1
    [A_grid_wake,B_grid_wake] = calcVelMatricesFast(grid,wake(j),dval);
    V = V + A_grid_wake*wake(j).gamma;
    U = U + B_grid_wake*wake(j).gamma;
end

% Calculate induction effects of far-field vortex sheets
for j = 1:length(wake)-1
    [tempNorm,tempTang] = calcFarFieldInduction(grid,wake(j));

    U = U + tempTang;
    V = V + tempNorm;
end

% Package output to same dimensions as input query
uGrid = reshape(U,size(xGrid));
vGrid = reshape(V,size(xGrid));

%% Reshape the velocities calculated at control points on the grid into the
%% meshgrid format for use in MATLAB's streamline function
%uGrid = zeros(length(xGrid(:,1)),length(xGrid(1,:)));
%vGrid = zeros(length(yGrid(:,1)),length(yGrid(1,:)));
%count = 1;
%
%% Add 1 to u component to account for nondimensionalized freestream
%for j = 1:length(xGrid(1,:))
%    for i = 1:length(xGrid(:,1))
%        uGrid(i,j) = U(count)+1;
%        vGrid(i,j) = V(count);
%        count = count + 1;
%    end
%end
%
% for i = 1:length(surfaces)
%     in = inpolygon(xGrid,yGrid,surfaces(i).endPoints(:,1),surfaces(i).endPoints(:,2));
%     uGrid(in) = 0;
%     vGrid(in) = 0;
% end

end

