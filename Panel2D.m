function [Cl,Cd,Cp,xc] = Panel2D(surfaces,alphaDeg,varargin)
% [Cl,Cd,Cp,xc] = Panel2D(surfaces,alphaDeg)
% [Cl,Cd,Cp,xc] = Panel2D(surfaces,alphaDeg,CT)
% [Cl,Cd,Cp,xc] = Panel2D(surfaces,alphaDeg,CT,wakeoptions)
% [Cl,Cd,Cp,xc] = Panel2D(___,Name,Value)
%
% Name-Value Arguments
% Plot: Visualization method, specified as 'none', 'stream', or 'contour'.
% Colormap: Colormap used by contour plots, specified as a method in the Colormap class.
% Mesh: Toggle for showing the mesh used for plotting, specified as 'on' or 'off'.

% Load default options
CT = 0;
opts.plot = 'none';
opts.colormap = 'vik';
opts.mesh = 'off';
opts.internalflow = 'off';
wakeoptions = wakeoptset(Wakesolver);

% Optional argument handling
if nargin > 2
    argflex = 0; % offset for Name-Value arguments
    % Check if a CT was passed
    if isnumeric(varargin{1})
        CT = varargin{1};
        argflex = 1;
        % Check for accompanying wake solver options
        if isa(varargin{min(2,length(varargin))},'Wakesolver')
            wakeoptions = varargin{2};
            argflex = 2;
        end
    end
    for i = 1+argflex:2:length(varargin)
        opts.(lower(varargin{i})) = varargin{i+1};
    end
end

if strcmpi(opts.internalflow,'on')
    [x0(1),k(1)] = min(surfaces{1}(:,1)); % identify LE
    [x0(2),k(2)] = min(surfaces{2}(:,1));
    x1 = max(x0); % set LB to deepest LE
    x2 = min(surfaces{1}(1,1),surfaces{2}(1,1)); % set UB to forwardmost TE
    sz = [50 51];
    x = linspace(x1,x2,sz(2));
    LW = interp1(surfaces{1}(k(1):end,1),surfaces{1}(k(1):end,2),x);
    UW = interp1(surfaces{2}(k(2):-1:1,1),surfaces{2}(k(2):-1:1,2),x);
end

alpha = alphaDeg*pi/180;
CHORD = surfaces{1}(1,1) - min(surfaces{1}(:,1));
nSurf = length(surfaces);

% Rotate surfaces coordinates by the requested angle of attack
for i = 1:nSurf
    surfaces{i} = surfaces{i}*[cos(alpha) -sin(alpha);sin(alpha) cos(alpha)];
    foil.m(i) = size(surfaces{i},1) - 1; % number of panels in each surface
end
M = sum(foil.m); % total number of panels in the airfoil system
foil.xo = zeros(M,1);
foil.yo = zeros(M,1);
foil.dx = zeros(M,1);
foil.dy = zeros(M,1);
AIC = zeros(M+nSurf);
m0 = 0;
for i = 1:nSurf
    foil.xo(m0+(1:foil.m(i))) = surfaces{i}(1:end-1,1);
    foil.yo(m0+(1:foil.m(i))) = surfaces{i}(1:end-1,2);
    foil.dx(m0+(1:foil.m(i))) = diff(surfaces{i}(:,1));
    foil.dy(m0+(1:foil.m(i))) = diff(surfaces{i}(:,2));

    % Kutta condition
    AIC(M+i,m0+i) = 1;
    AIC(M+i,m0+i+foil.m(i)) = 1;

    m0 = m0 + foil.m(i);
end
foil.xc = foil.xo + 0.5*foil.dx;
foil.yc = foil.yo + 0.5*foil.dy;
foil.theta = atan2(foil.dy,foil.dx);

% AIC produces induced velocity components normal to the surface boundaries
% B produces the tangential components
[AIC(1:M,:),B] = influence(foil,foil);
RHS = [sin(foil.theta); zeros(nSurf,1)];

if CT == 0
    foil.G = AIC \ RHS;
    Qtan = B*foil.G + cos(foil.theta);
    wake.m = wakeoptions.NumPanels;
else
    Ginf = sqrt(2*CT + 1) - 1; % far downstream wake circulation
    N = wakeoptions.NumPanels; % number of panels per wake

    % Attach wakes to the interior TE points where the last panel on each wake
    % has constant circulation strength Ginf and extends to downstream infinity
    wakeEnd = surfaces{1}(end,1) + CHORD*wakeoptions.WakeLengthChords;
    wake.m = [N N];
    c_x = (wakeEnd-[surfaces{1}(end,1) surfaces{2}(1,1)]).*(1-cos(linspace(0,pi/2,N).')) + [surfaces{1}(end,1) surfaces{2}(1,1)];
    wake.xo = c_x(:);
    % wake.xo = [linspace(surfaces{1}(end,1),wakeEnd,N) linspace(surfaces{2}(1,1),wakeEnd,N)].';
    wake.dx = [diff(wake.xo(1:N)); 999; diff(wake.xo(N+1:2*N)); 999];
    wake.xc = wake.xo + 0.5*wake.dx;
    % Initial wake shapes are flat sheets
    wake.yo = [surfaces{1}(end,2)+zeros(N,1); surfaces{2}(1,2)+zeros(N,1)];
    wake.dy = zeros(size(wake.xo));
    wake.yc = wake.yo;
    wake.theta = zeros(size(wake.xo));
    % Initial wake circulation
    wake.G = [Ginf+zeros(N+1,1); -Ginf+zeros(N+1,1)];

    % Solve for the wake shape and global circulation solution
    [foil.G,wake,BB] = solvewake(wakeoptions,wake,foil,inv(AIC),RHS,CHORD);

    % Extract results from the solution
    Qtan = B*foil.G + BB*wake.G + cos(foil.theta); % surface tangent velocity
end

Cp = 1 - Qtan.^2;
Cl = -Cp.'*foil.dx;
Cd = Cp.'*foil.dy;

% Package output where each cell contains data for one surface
Cp = mat2cell(Cp,foil.m);
xc = mat2cell(foil.xc*cos(alpha)-foil.yc*sin(alpha),foil.m); % undo AOA rotation

if strcmpi(opts.internalflow,'on')
    XC = repmat(x,sz(1),1);
    Y = (UW - LW).*linspace(0,1,sz(1)+1).' + LW;
    h = diff(Y,1,1);
    YC = Y(1:end-1,:) + 0.5*h;
    tmp = reshape([XC YC],[],2)*[cos(alpha) -sin(alpha);sin(alpha) cos(alpha)];
    duct.xc = tmp(:,1);
    duct.yc = tmp(:,2);
    % set remaining parameters
    duct.m = length(duct.xc);
    duct.xo = duct.xc - 0.5;
    duct.yo = duct.yc;
    duct.dx = ones(duct.m,1);
    duct.dy = zeros(duct.m,1);
    duct.theta = zeros(duct.m,1);

    [A,B] = influence(duct,foil);
    u = B*foil.G + 1;
    v = A*foil.G;
    prj = reshape(u*cos(alpha) - v*sin(alpha), sz);
    mdot = sum(prj.*h,1);

    figure;
    subplot(2,1,1);
    contourf(XC,YC,prj,200,'LineStyle','none');
    if strcmpi(opts.mesh,'on')
        hold on;
        plot(XC,Y,'k');
        plot(XC.',Y.','k');
    end
    daspect([1 1 1]);
    xlim([x1 x2]);
    subplot(2,1,2);
    plot(x,mdot);
    xlim([x1 x2]);
    set(gcf,'Position',[300 200 560 600]);
end

if strcmpi(opts.plot,'none'), return; end
pltC(surfaces,foil,wake,opts);