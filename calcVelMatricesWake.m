function [A,B] = calcVelMatricesWake(surface1,surface2)
% calcVelMatricesWake Calculates the A and B matrices for a set of
% collocation points defined in the 'surface1' structure, influenced by the
% nodes of the 'surface2' structure. Inspired by the linear-strength,
% vortex panel method given by Katz & Plotkin (2003).
%
% This function is barely modified from calcVelMatrices.m and is only
% called once in calcWakeVels.m when calculating wake self influence. The
% modifications are made between lines 49 and 67 of this file.
%
% Inputs:  surface1 - this surface structure contains the control points on
%          which the influence of other surfaces is calculated upon.
%          surface2 - this surface structure contains the nodes that
%          influence the control points of surface1
% Outputs: A - the matrix that defines the normal component of the velocity
%          induced at the control points on surface1 influenced by the
%          nodes of surface2
%          B - the matrix that defines the horizontal component of the
%          velocity induced at the control points on surface1 influenced by
%          the nodes of surface2

CO = surface1.co;
PT1 = surface2.pt1;
PT2 = surface2.pt2;
THS1 = surface1.theta;
THS2 = surface2.theta;
A = zeros(surface1.m,surface2.n);
B = zeros(surface1.m,surface2.n);

for i = 1:surface1.m
    for j = 1:surface2.m
%       Determine location of collocation point i in terms of panel j
%       coordinates
        XT=CO(i,1)-PT1(j,1);
        ZT=CO(i,2)-PT1(j,2);
        X2T=PT2(j,1)-PT1(j,1);
        Z2T=PT2(j,2)-PT1(j,2);

        X=XT*cos(THS2(j))+ZT*sin(THS2(j));
        Z=-XT*sin(THS2(j))+ZT*cos(THS2(j));
        X2=X2T*cos(THS2(j))+Z2T*sin(THS2(j));

%       Determine radial distance and angle between corner points of jth
%       panel and ith control point
        R1=sqrt(X^2+Z^2);
        R2=sqrt((X-X2)^2+Z^2);
        TH1=atan2(Z,X);
        TH2=atan2(Z,X-X2);

        % In this version of the function, we adjust calculated angles for
        % panel self-induction to ensure the tangential component of
        % self-induction is zero as per intuition. This goes against the
        % Katz & Plotkin code, and has something to do with how MATLAB
        % handles the atan2 function when calculating the panel theta, but
        % it isn't clear why this needs to be handled differently when
        % calculating self-induction from the wake panels.
        if strcmp(surface1.name,surface2.name)
            if i==j
                TH2 = 0;
                TH1 = 0;
            end
        end

        U1L=-(Z*log(R2/R1)+X*(TH2-TH1)-X2*(TH2-TH1))/(2*pi*X2);
        U2L=(Z*log(R2/R1)+X*(TH2-TH1))/(2*pi*X2);
        W1L=-((X2-Z*(TH2-TH1))-X*log(R1/R2)+X2*log(R1/R2))/(2*pi*X2);
        W2L=((X2-Z*(TH2-TH1))-X*log(R1/R2))/(2*pi*X2);

%       Rotate coordinates back from jth panel reference frame to airfoil
%       chord frame
        U1=U1L*cos(-THS2(j))+W1L*sin(-THS2(j));
        U2=U2L*cos(-THS2(j))+W2L*sin(-THS2(j));
        W1=-U1L*sin(-THS2(j))+W1L*cos(-THS2(j));
        W2=-U2L*sin(-THS2(j))+W2L*cos(-THS2(j));

%       A(i,j) is the component of velocity normal to control
%       point i due to panel j
%       B(i,j) is the tangential velocity along control point i due to
%       panel j, used after solving for gammas
        if j==1
            A(i,1)=-U1*sin(THS1(i))+W1*cos(THS1(i));
            HOLDA=-U2*sin(THS1(i))+W2*cos(THS1(i));
            B(i,1)=U1*cos(THS1(i))+W1*sin(THS1(i));
            HOLDB=U2*cos(THS1(i))+W2*sin(THS1(i));
        elseif j==surface2.m
            A(i,surface2.m)=-U1*sin(THS1(i))+W1*cos(THS1(i))+HOLDA;
            A(i,surface2.n)=-U2*sin(THS1(i))+W2*cos(THS1(i));
            B(i,surface2.m)=U1*cos(THS1(i))+W1*sin(THS1(i))+HOLDB;
            B(i,surface2.n)=U2*cos(THS1(i))+W2*sin(THS1(i));
        else
            A(i,j)=-U1*sin(THS1(i))+W1*cos(THS1(i))+HOLDA;
            HOLDA=-U2*sin(THS1(i))+W2*cos(THS1(i));
            B(i,j)=U1*cos(THS1(i))+W1*sin(THS1(i))+HOLDB;
            HOLDB=U2*cos(THS1(i))+W2*sin(THS1(i));
        end
    end
end

end

