function [A,B] = influence(sys1,sys2,D)
% Calculates the influence matrices associated with
% the velocities induced onto sys1 due to sys2.
%
% Pass the optional argument D if the difference in
% segment angles for self-induction is not 180 deg.
%
% The system fields are defined as follows:
%      m: ordered array of the number of panels in each subelement of the system
%   x|yo: left endpoints of the subelements
%   x|yc: control points of the subelements
%   dx|y: spans of panels in dimension x or y
%  theta: panel angles (rad)
if nargin == 2
    D = pi;
end
M1 = sum(sys1.m);
M2 = sum(sys2.m);
nSurf = length(sys2.m);

% Convert control point to local panel coords
xt = sys1.xc - sys2.xo.';
yt = sys1.yc - sys2.yo.';

% Precompute trig expressions as 1-by-M2 arrays
costh = cos(sys2.theta).';
sinth = sin(sys2.theta).';

% Find control point coords in CS of all panels
xp = xt.*costh + yt.*sinth;
yp = -xt.*sinth + yt.*costh;
x2 = sqrt(sys2.dx.^2 + sys2.dy.^2).'; % 1-by-M2

% Find theta1, theta2, r1, r2
theta1 = atan2(yp,xp);
theta2 = atan2(yp,xp-x2);
dtheta = theta2 - theta1; % precompute
if isequal(sys1.id,sys2.id)
    dtheta(logical(eye(M1))) = D; % fix angular difference for self-induction
end

r1 = sqrt(xp.^2 + yp.^2);
r2 = sqrt((xp-x2).^2 + yp.^2);
ln = log(r2./r1); % precompute

% Compute influence coefficients
ap = yp.*ln + xp.*dtheta;
bp = xp.*ln + x2 - yp.*dtheta;
am = -yp.*ln + (x2 - xp).*dtheta;
bm = (x2 - xp).*ln - x2 + yp.*dtheta;
c = 1./(2*pi*x2);

% Velocity decomposition used in Katz and Plotkin Eq. 11.103
ua = c.*(am.*costh - bm.*sinth);
ub = c.*(ap.*costh - bp.*sinth);
va = c.*(am.*sinth + bm.*costh);
vb = c.*(ap.*sinth + bp.*costh);

u = zeros(M1,M2+nSurf);
v = zeros(M1,M2+nSurf);
m0 = 0;
for i = 1:nSurf
    u(:,m0+i+(0:sys2.m(i))) = [ua(:,m0+(1:sys2.m(i))), zeros(M1,1)] + [zeros(M1,1), ub(:,m0+(1:sys2.m(i)))];
    v(:,m0+i+(0:sys2.m(i))) = [va(:,m0+(1:sys2.m(i))), zeros(M1,1)] + [zeros(M1,1), vb(:,m0+(1:sys2.m(i)))];
    m0 = m0 + sys2.m(i);
end

A = -u.*sin(sys1.theta) + v.*cos(sys1.theta);
B = u.*cos(sys1.theta) + v.*sin(sys1.theta);