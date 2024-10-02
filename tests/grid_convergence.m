function grid_convergence
close all

figure
hold on

nmax1 = 10;
nmax2 = 5;

for i = 1:nmax1
    Nx = 2^i + 1;
    surfaces{1} = airfoilCST([0.17312 0.15047 0.16464 0.13569 0.14328 0.13957], ...
                    -[0.17312 0.15047 0.16464 0.13569 0.14328 0.13957], Nx);
    surfaces{2} = 0.5*surfaces{1} + [0.5 0.2];
    for j = 1:nmax2
        fprintf('Running case %d of %d\n', (i-1)*nmax2+j, nmax1*nmax2)
        options = wakeoptset(Wakesolver,'NumPanels',2^(j+3),'FunctionTolerance',1e-4,'Display','off');
        [Cl,~,~,~] = Panel2D(surfaces,10,1,options);
        scatter(i+1,Cl,'magenta','filled','d','MarkerEdgeColor','b')
        drawnow
    end
end
fprintf('Done\n')
end

function coords = airfoilCST(Au,Al,Nx,N1,N2)
% Input:
%     Au: upper surface weights
%     Al: lower surface weights
%     Nx: number of points to discretize the chord (output contains 2*Nx-1 points)
%     N1: (optional) first shape parameter
%     N2: (optional) second shape parameter
%
% Output:
% coords: (2*Nx-1)-by-2 array of airfoil coordinates with x in the first column and y in the second column
%
% Kulfan, B. M., "Universal Parametric Geometry Representation Method," Journal of Aircraft, Vol. 45, No. 1, 2008, pp. 142-158.
% https://doi.org/10.2514/1.29958

if nargin == 3
    N1 = 0.5;
    N2 = 1;
end

n = length(Au)-1; % Set Bernstein polynomial order to no. weights - 1

% Generate x/c grid from TE to LE
x = 0.5*(1 + cos(linspace(0,pi,Nx))).';
% Class function
C = x.^N1 .* (1-x).^N2;
% Component shape functions
i = 0:n;
K = factorial(n)./factorial(i)./factorial(n-i);
S = K.*x.^i .* (1-x).^(n-i);
% Overall shape function
Su = S*Au(:);
Sl = S*Al(:);

yu = C.*Su;
yl = C.*Sl;

coords = [x(1:Nx-1) yl(1:Nx-1);flipud(x) flipud(yu)];
end