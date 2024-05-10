% designDriver Runs the main APPC-2D program and plots Cp distributions.

close all; clear; clc;

% File definition - the surface from which the upper wake boundary emanates
% must come first, followed by the surface from which the lower wake
% boundary emanates. All other surfaces must follow these two examples.

vecAng = 0; % degrees

surfaceFiles(1).name = sprintf('./Airfoils/ONR-Coords/nacelleVec%g.dat',vecAng);
surfaceFiles(2).name = sprintf('./Airfoils/ONR-Coords/mainVec%g.dat',vecAng);
% surfaceFiles(3).name = './Airfoils/ONR-Coords/krueger.dat';

% surfaceFiles(1).name = './Airfoils/NACA-643618/nacelle.dat'; % nacelle
% surfaceFiles(2).name = './Airfoils/NACA-643618/airfoil.dat'; % airfoil

% Algorithm call
alpha = 10; % degrees
CT = 10; % non-dimensional thrust coefficient
[surfaces,wake,iter,chord] = main(surfaceFiles,alpha,CT);

%% STREAMLINES/VECTORS
f2 = figure(2);
set(f2,'Position',[100 100 1000 800]);
hold on;
for i = 1:length(surfaces)
    plot(surfaces(i).endPoints(:,1),surfaces(i).endPoints(:,2),'k-');
end

xG = (-chord/2):0.1:(chord*3);
yG = -chord*2:0.05:chord*2;

% Activate starting points for streamlines
xStart = xG(1).*ones(length(yG),1);
yStart = yG';

[xGrid,yGrid] = meshgrid(xG,yG);
[uGrid,vGrid] = createStreamlines(xGrid,yGrid,surfaces,wake);

% Choose whether to plot streamlines or velocity vectors
streamPlot = streamline(xGrid,yGrid,uGrid,vGrid,xStart,yStart);
% quiverPlot = quiver(xGrid,yGrid,uGrid,vGrid,'AutoScaleFactor',0.55,...
%                                            'Color',[0 0.4470 0.7410]);

% Coplot wake boundaries
wakePlot1 = plot(wake(1).endPoints(:,1),wake(1).endPoints(:,2),'-',...
                                      'Color',[0.8500 0.3250 0.0980]);
wakePlot2 = plot(wake(2).endPoints(:,1),wake(2).endPoints(:,2),'-',...
                                      'Color',[0.8500 0.3250 0.0980]);

xlabel('$$x/c$$','Interpreter','latex');
ylabel('$$z/c$$','Interpreter','latex');
axis([-chord/2 chord*3 -chord/2 chord/2]);
axis equal;

%% PRESSURE DISTRIBUTIONS
% Derotate airfoil surfaces to ensure pressure distributions are plotted on
% the correct x/c scale, and then repanel them (see panelAirfoils.m).
alRad = alpha.*pi./180;
R = [cos(alRad) -sin(alRad); sin(alRad) cos(alRad)];

for i = 1:length(surfaces)
    center = repmat([0.25;0],1,length(surfaces(i).endPoints(:,1)));
    surfaces(i).endPoints = (R*(surfaces(i).endPoints'-center) + center)';
    surfaces(i).pt1 = surfaces(i).endPoints(1:surfaces(i).m,:);
    surfaces(i).pt2 = surfaces(i).endPoints(2:surfaces(i).n,:);
    surfaces(i).co = surfaces(i).pt1 + ((surfaces(i).pt2-surfaces(i).pt1)./2);
    surfaces(i).theta = atan2((surfaces(i).pt2(:,2)-surfaces(i).pt1(:,2)),...
                          (surfaces(i).pt2(:,1)-surfaces(i).pt1(:,1)));
end

% Selection of colors for plotting to cycle through
color = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; ...
                              0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880];

f3 = figure(3);
set(f3,'Position',[100 100 1000 800]);
hold on;
set(gca, 'YDir','reverse');

% The following for loop cycles through the surfaces and plots the Cp
% distributions. It also intergrates the Cp distribution to calculate the
% integrated lift and quarter-chord pitching moment coefficients (recall
% the lift coefficient is calculated separately in main.m using
% Kutta-Joukowski theorem and the airfoil circulation distributions).
Cn = 0;
Ca = 0;
Cm = 0;

for i = length(surfaces):-1:1

    % Find the leading edge such that the airfoil upper and lower surfaces
    % can be plotted separately in different styles
    target = find(surfaces(i).co(:,1)==min(surfaces(i).co(:,1)));

    plot((surfaces(i).co(1:target,1)),surfaces(i).CP(1:target),...
                                     'Color',color(i,:),'LineStyle','--');
    plot((surfaces(i).co(target+1:end,1)),surfaces(i).CP(target+1:end,1),...
                                     'Color',color(i,:),'LineStyle','-');

    % Calculation of normal, axial, and moment coefficients from Cp dist.
    for j = 1:length(surfaces(i).CP)
        if j <= target
            Cn = Cn + surfaces(i).CP(j).*abs(surfaces(i).DL(j).*cos(surfaces(i).theta(j)));
            Ca = Ca - surfaces(i).CP(j).*abs(surfaces(i).DL(j).*sin(surfaces(i).theta(j)));
            Cm = Cm + surfaces(i).CP(j).*abs(surfaces(i).DL(j).*cos(surfaces(i).theta(j))).*(0.25-surfaces(i).co(j,1)) - ...
                      surfaces(i).CP(j).*abs(surfaces(i).DL(j).*sin(surfaces(i).theta(j))).*surfaces(i).co(j,2);
        else
            Cn = Cn - surfaces(i).CP(j).*abs(surfaces(i).DL(j).*cos(surfaces(i).theta(j)));
            Ca = Ca + surfaces(i).CP(j).*abs(surfaces(i).DL(j).*sin(surfaces(i).theta(j)));
            Cm = Cm - surfaces(i).CP(j).*abs(surfaces(i).DL(j).*cos(surfaces(i).theta(j))).*(0.25-surfaces(i).co(j,1)) + ...
                      surfaces(i).CP(j).*abs(surfaces(i).DL(j).*sin(surfaces(i).theta(j))).*surfaces(i).co(j,2);
        end
    end
end

% Rotate normal and axial force coefficients to find lift coefficient. Drag
% coefficient can also be calculated but it is useless with the potential
% flow assumptions made in this code.
R = [cos(alRad) -sin(alRad); sin(alRad) cos(alRad)];
coeff = R*[Cn;Ca];
fprintf('\nLift (Cp Dist): Cl = %1.4f\n',coeff(1));
fprintf('1/4 Moment (Cp Dist): Cm_c/4 = %1.4f\n\n',Cm);

axis([-0.2 chord*1.1 -15 2]);
xlabel('$$x/c$$','Interpreter','latex');
ylabel('$$C_p$$','Interpreter','latex');
titleStr = "$\alpha=$" + alpha + "$^\circ$, $C_T=$" + CT;
title(titleStr,'Interpreter','latex');

