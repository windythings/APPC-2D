function [surfaces,wake,iter,chord] = main(surfaceFiles,alphaDeg,CT)
% main Main program of Aero-Propulsive Panel Code (APPC-2D). Assesses a
% given set of 2D aerodynamic surfaces and jet parameters to determine the
% interaction between the jet wake and aerodynamic surfaces. An iterative
% scheme based on a 2D, linear vortex strength panelling method is
% implemented to determine the jet wake shape and strength, as well as
% system lift-coefficient and thrust-coefficient.

% Assumptions: 2D, incompressible, inviscid, irrotational, no body forces
%
% Inputs: surfaceFiles - struct with names of all airfoil coordinate files
%         alphaDeg - system angle of attack, in degrees
%         CT - propulsive system thrust coefficient
%
% Outputs: surfaces - struct containing details of all surface paneling
%          wake - struct containing details of all wake paneling
%          CL -

%% Initialize Constants
% The default values have been found to work well with the example airfoil
% coordinates provided. Change as necessary for specific cases.

CHORD_MULTIPLIER = 5; % no. of chord lengths to extend the jet wake back
CHORDS_RAMP = 2; % no. of chord lengths to ramp circulation to freestream
N_WAKE = 200; % no. of wake endpoints to create
WAKE_SPACING = 'cosine'; % uniform or cosine-type spacing for wake points
RELAX_CIRC = 0.5; % circulation relaxation factor
E_UP = 1E-4; % upper wake shape residual
E_LO = 1E-4; % lower wake shape residual
MAX_ITER = 30; % maximum iterations for wake shape loop calculation

%% Setup Calculations
% Calculate farfield circulation value based on CT and the initial guess
% for the "average" velocity on the wake boundaries.

alpha = alphaDeg * pi/180;
gammaInf = (sqrt(2*CT+1)-1);

%% Airfoil/Nacelle Panelling
% Create initial surface struct containing all paneling details on the
% airfoil surfaces, and calculate influence coefficient matrix.
[surfaces,chord] = panelAirfoils(surfaceFiles,alpha,2);
[A,~] = calcAirfoilMatrices(surfaces);

%% Wake Panelling
% Creates the wake (x,y) endpoints based on the surface trailing edge (x,y)
% points. Index 1 of the wake struct always refers to the "upper" wake, or
% that with defined negative circulation. Index 2 is the "lower" wake and
% index 3 is the (imaginary) wake centerline.

% The wake boundaries always start at the trailing edges of the first two
% surface elements specified.
wakeStart = [surfaces(1).endPoints(end,:);surfaces(2).endPoints(end,:)];
wake = createWake(wakeStart,N_WAKE,chord,CHORD_MULTIPLIER,...
                                                 CHORDS_RAMP,WAKE_SPACING);

% The wake is initialized to have unit velocity.
wake(1).uBar = ones(N_WAKE-1,1);
wake(2).uBar = ones(N_WAKE-1,1);

% The gamma vector is only used to calculate residuals. The current wake
% gamma value is stored in the "wake" date structure.
gamma(1).upper = zeros(length(wake(1).uBar)+1,1);
gamma(1).lower = zeros(length(wake(2).uBar)+1,1);

% Calculates the bound circulation to the wake, along with the decaying
% vorticity down to gammaInf. No relaxation for this iteration.
% In this zeroth iteration, CT=gammaInf and uBar=1 so that CT/uBar=gammaInf
wake = wakeCirculation(wake,gammaInf,chord,CHORD_MULTIPLIER,wakeStart,...
                                                        1,gamma,gammaInf);

% Fill the "previous" (index 1) value of wake gammas for residual
% calculation later
gamma(1).upper = wake(1).gamma;
gamma(1).lower = wake(2).gamma;

iter = 0;
toEnd = 0;

% Set up the figure to plot the wake iteration. Axes are based on the CHORD
% of the airfoil system, but can be adjusted as necessary. This figure is
% closed upon completion of the main algorithm.
f1 = figure(1);
hAx = axes;
set(f1,'Position',[100 100 1000 800]);
hold on;
axis([-chord 6.*chord -2.*chord 2.*chord]);
axis equal;

for i = 1:length(surfaces)
    plot(surfaces(i).endPoints(:,1),surfaces(i).endPoints(:,2),'k-');
end
wakePlot1 = plot(hAx,wake(1).endPoints(:,1),wake(1).endPoints(:,2),'r-');
wakePlot2 = plot(hAx,wake(2).endPoints(:,1),wake(2).endPoints(:,2),'b-');
wakePlot3 = plot(hAx,wake(3).endPoints(:,1),wake(3).endPoints(:,2),'g-');

%% Iterative Scheme

% Main loop used to execute the iterative algorithm. The loop is ended if
% the maximum iteration number is reached or if residuals are within the
% specified constant values.
while (iter <= MAX_ITER) && (~toEnd)
    iter = iter + 1;

    % Wake boundaries are re-plotted with each algorithm iteration.
    set(wakePlot1,'YData',wake(1).endPoints(:,2));
    set(wakePlot2,'YData',wake(2).endPoints(:,2));
    set(wakePlot3,'YData',wake(3).endPoints(:,2));
    drawnow; % Creates "stop-motion animation" of wake shape moving

    % The iterative scheme performs the following steps to calculate the
    % final aero-propulsive solution.
    % 1. Panel the wake boundaries (panelWake.m)
    % 2. Calculate the right-hand-side matrix from freestream and wake
    %    induction effects (calcRHS.m)
    % 3. Invert the aerodynamic influence coefficient matrix and solve for
    %    circulation distribution on airfoil surfaces (calcAirfoilGamma.m)
    % 4. Calculate airfoil and wake induction on the wakes (calcWakeVels.m)
    % 5. Calculate "average" velocity on wake boundaries (calcUBar.m)
    % 6. Recalculate the circulation distribution on the wake boundaries
    %    based on CT and updated uBar (wakeCirculation.m)
    % 7. Change the position of the wake boundaries (iterateWakePosition.m)

    wake = panelWake(wake);
    RHS = calcRHS(surfaces,wake);
    surfaces = calcAirfoilGamma(surfaces,A,RHS);

    wake = calcWakeVels(wake,surfaces);
    wake = calcUBar(wake);
    wake = wakeCirculation(wake,CT,chord,CHORD_MULTIPLIER,wakeStart,...
                                              RELAX_CIRC,gamma,gammaInf);
    wake = iterateWakePosition(wake);

    % Calculation of residuals for loop break condition and resetting of
    % gamma vector with updated wake circulation values.
    gamma(2).upper = wake(1).gamma;
    gamma(2).lower = wake(2).gamma;

    residualUpper = sum(abs(gamma(2).upper-gamma(1).upper))./...
                                                (N_WAKE.*abs(gammaInf)+eps);
    residualLower = sum(abs(gamma(2).lower-gamma(1).lower))./...
                                                (N_WAKE.*abs(gammaInf)+eps);

    gamma(1).upper = gamma(2).upper;
    gamma(1).lower = gamma(2).lower;

    if (residualUpper < E_UP) && (residualLower < E_LO)
        toEnd = 1;
    end
end

% After the iterative algorithm completes, the final wake position will
% induce new velocities on the airfoil surfaces. Therefore, the wake must
% be re-paneled and the circulation distribution on the airfoils must be
% re-calculated.
wake = panelWake(wake);
RHS = calcRHS(surfaces,wake);
surfaces = calcAirfoilGamma(surfaces,A,RHS);

%% CL & CP Calculations
% Calculates 2D inviscid lift and pressure coefficients. Warning: these can
% be misleading as vortex panel codes with no viscous modelling tend to
% overpredict performance, sometimes drastically.

[surfaces,CL] = calcCoefficients(surfaces,wake);
fprintf('Lift (K-J): Cl = %1.4f\n',CL);

% The following is a currently unverified mathematical formulation to
% calculate the contribution of the jet momentum to the system lift
% coefficient. This was not factored into previous calculations.

% Normalized nozzle exit area and applied thrust angle
A9 = sqrt((wake(1).endPoints(1,1)-wake(2).endPoints(1,1)).^2 + ...
          (wake(1).endPoints(1,2)-wake(2).endPoints(1,2)).^2);
thrustAngle = wake(3).theta(1)+pi;

% Calculate this by substituting the thrust equation in for lift and
% normalizing by 1/2*rho*vInf^2*c. Since many parameters are already
% non-dimensional, many things cancel. You must also multiply by the sine
% of the thrust angle.
clJet = 2.*A9.*(wake(3).vTang(1).^2-wake(3).vTang(1)).*sin(thrustAngle);
fprintf('Lift (Jet): Cl = %1.4f\n',clJet);
fprintf('Lift (Total): Cl = %1.4f\n',CL+clJet);
