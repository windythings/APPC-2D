function [A,B] = calcVelMatricesFast(surface1,surface2)

dval = pi;

if strcmp(surface1.name,surface2.name)
    isself = true;
    if ~isempty(regexp(surface1.name,'Wake.*'))
        dval = 0;
    end
else
    isself = false;
end

% convert control points to local panel coordinates
xt = surface1.co(:,1) - surface2.pt1(:,1).';
yt = surface1.co(:,2) - surface2.pt1(:,2).';

% precompute trig expressions as 1-by-surface2.m arrays
costh = cos(surface2.theta).';
sinth = sin(surface2.theta).';

% find control point coordinates in CS of all panels
xp = xt.*costh + yt.*sinth;
yp = -xt.*sinth + yt.*costh;
x2 = sqrt((surface2.pt2(:,1)-surface2.pt1(:,1)).^2 + ...
          (surface2.pt2(:,2)-surface2.pt1(:,2)).^2).';

% find angles and distances between the two systems' panels
theta1 = atan2(yp,xp);
theta2 = atan2(yp,xp-x2);
dtheta = theta2 - theta1; % precomputation

if isself
    dtheta(logical(eye(surface1.m))) = dval;
end

r1 = sqrt(xp.^2 + yp.^2);
r2 = sqrt((xp-x2).^2 + yp.^2);
ln = log(r2./r1); % precomputation

% compute influence coefficients
ap = yp.*ln + xp.*dtheta;
bp = xp.*ln + x2 - yp.*dtheta;
am = -yp.*ln + (x2 - xp).*dtheta;
bm = (x2 - xp).*ln - x2 + yp.*dtheta;

ua = (am.*costh - bm.*sinth)./(2*pi*x2);
ub = (ap.*costh - bp.*sinth)./(2*pi*x2);
va = (am.*sinth + bm.*costh)./(2*pi*x2);
vb = (ap.*sinth + bp.*costh)./(2*pi*x2);

u = [ua, zeros(surface1.m,1)] + [zeros(surface1.m,1), ub];
v = [va, zeros(surface1.m,1)] + [zeros(surface1.m,1), vb];

A = -u.*sin(surface1.theta) + v.*cos(surface1.theta);
B = u.*cos(surface1.theta) + v.*sin(surface1.theta);

end
