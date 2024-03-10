classdef Wakesolver
    properties
        MaxIterations = 50
        FunctionTolerance = 1e-6
        RelaxationFactor = 0.5
        NumPanels = 200
        WakeLengthChords = 7
        GammaDecayChords = 2
        Display = 'iter'
    end
    methods
        function obj = Wakesolver; end
        function obj = wakeoptset(obj,varargin)
            % Updates the solver options
            if isempty(varargin)
                return
            end
            for i = find(ismember(varargin(1:2:length(varargin)),properties(obj)))
                obj.(varargin{2*i-1}) = varargin{2*i};
            end
        end
        function [gamma,wake,B] = solvewake(obj,wake,foil,AIC,RHS,chord)
            % Calculates the global circulation solution according to the Shollenberger algorithm:
            %   1. Guess initial wake shape and bound circulation
            %   2. Calculate induced velocities on the airfoil surfaces due to the wake
            %   3. Solve circulation strengths on the airfoil surfaces such that there is no flow
            %      through the solid boundaries
            %   4. Calculate induced velocities on the wake
            %   5. Move the wake panels to satisfy no flow through the wake boundaries
            %   6. Calculate bound circulation on the wake from the new velocity inside the wake
            %   7. If the new wake circulation is close to the previous solution, exit the program,
            %      else return to step 2
            nSurf = length(foil.m); % number of airfoil surfaces
            N = obj.NumPanels;      % number of panels per wake
            Ginf = abs(wake.G(N+1));
            CT = 0.5*((Ginf + 1)^2 - 1);
            wakeEnd = wake.xo(N);
            xN = wakeEnd - chord*obj.GammaDecayChords;
            idN = [find(wake.xo(1:N)>xN,1) find(wake.xo(N+1:2*N)>xN,1)+N];

            figure
            w1 = plot(wake.xo(1:N),wake.yo(1:N),'b-');
            hold on
            w2 = plot(wake.xo(N+1:2*N),wake.yo(N+1:2*N),'r-');
            m0 = 0;
            for i = 1:nSurf
                plot(foil.xo(m0+[1:foil.m(i),1]),foil.yo(m0+[1:foil.m(i),1]),'k-')
                m0 = m0 + foil.m(i);
            end
            daspect([1 1 1])
            pbaspect([4 3 1])

            iter = 0;
            E = 1;
            while E > obj.FunctionTolerance
                iter = iter + 1;
                if iter > obj.MaxIterations
                    break
                end
                % Solve for the airfoil circulation
                [A,B] = influence(foil,wake);
                RH = [A*wake.G; -wake.G([1 N+2]); zeros(nSurf-2,1)];
                gamma = AIC \ (RHS - RH);

                % Calculate induced velocities on the wake panels
                [A1,B1] = influence(wake,foil);
                [A2,B2] = influence(wake,wake,0);
                up = B1*gamma + B2*wake.G;
                vp = A1*gamma + A2*wake.G;
                u = up.*cos(wake.theta) - vp.*sin(wake.theta) + 1;
                v = up.*sin(wake.theta) + vp.*cos(wake.theta);

                % Move wake position so the wake boundaries are streamlines
                wake.dy = v./u.*wake.dx;
                wake.dy([N 2*N]) = 0;
                wake.yo(2:N) = wake.yo(1) + cumsum(wake.dy(1:N-1));
                wake.yo(N+2:2*N) = wake.yo(N+1) + cumsum(wake.dy(N+1:2*N-1));

                % Handle wakes crossing
                xing = (wake.yo(1:N)-wake.yo(N+1:2*N)) > 0;
                if any(xing) && any(~xing)
                    wake.yo(N+2:2*N) = wake.yo(2:N) - wake.yo(1) + wake.yo(N+1);
                    wake.dy(N+1:2*N) = wake.dy(1:N);
                end

                wake.yc = wake.yo + 0.5*wake.dy;
                wake.theta = atan2(wake.dy,wake.dx);

                % Update wake circulation using velocity at the jet boundary
                Gprev = wake.G;
                G = CT./sqrt(u.^2 + v.^2);
                G(N+1:2*N) = -G(N+1:2*N);
                wake.G(1:N-1) = 0.5*[3*G(1)-G(2); G(1:N-2)+G(2:N-1)]; % interpolate from control points to vertices
                wake.G(N+2:2*N) = 0.5*[3*G(N+1)-G(N+2); G(N+1:2*N-2)+G(N+2:2*N-1)];

                % Decay vorticity to far downstream value
                % GN is the interpolated G at xN
                GN = (wake.G(idN(1))-wake.G(idN(1)-1))/wake.dx(idN(1)-1)*(xN-wake.xo(idN(1))) + wake.G(idN(1));
                wake.G(idN(1):N) = (wake.G(N)-GN)/(wakeEnd-xN).*(wake.xo(idN(1):N)-xN) + GN;
                GN = (wake.G(idN(2)+1)-wake.G(idN(2)))/wake.dx(idN(2)-1)*(xN-wake.xo(idN(2))) + wake.G(idN(2)+1);
                wake.G(idN(2)+1:2*N+1) = (wake.G(2*N+1)-GN)/(wakeEnd-xN).*(wake.xo(idN(2):2*N)-xN) + GN;

                % Calculate residual between current and previous iteration
                E = sum(abs(wake.G - Gprev))/(2*N*Ginf);

                wake.G = obj.RelaxationFactor*wake.G + (1 - obj.RelaxationFactor)*Gprev;

                % Create a movie of the evolving wake shape
                if strcmpi(obj.Display,'iter'), pause(0.01); end
                delete(w1)
                delete(w2)
                w1 = plot(wake.xo(1:N),wake.yo(1:N),'b-');
                w2 = plot(wake.xo(N+1:2*N),wake.yo(N+1:2*N),'r-');
            end
            if strcmpi(obj.Display,'off') || strcmpi(obj.Display,'none')
                close(gcf)
            else
                pause(0.01)
            end
        end
    end
end
