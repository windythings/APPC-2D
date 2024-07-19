function graf = pltC(surfaces,foil,wake,opts)
    graf = figure; hold on;

    nSurf = length(foil.m);
    N = wake.m(1);
    [xmin,xmax] = bounds(foil.xo);
    [ymin,ymax] = bounds(foil.yo);
    ch = xmax - xmin;

    if strcmpi(opts.plot,'contour')
        % Plots the velocity magnitude field on an efficient grid
        filename = mfilename('fullpath');
        filepath = fileparts( filename );
        addpath([filepath,'/mesh2d']);

        initmsh();

        % Define area of interest
        node = [xmax+ch  , ymax+ch/4;
                xmin-ch/2, ymax+ch/4;
                xmin-ch/2, ymin-ch/2;
                xmax+ch  , ymin-ch/2];

        % Place evaluation points barely outside the boundary because
        % exact boundary evaluation is unstable
        bubble = repmat({[]},1,nSurf); % the neighborhood of points near the solid surfaces
        m0 = 0;
        for i = 1:nSurf
            bubble{i} = offsetBoundary(surfaces{i},foil.theta((1:foil.m(i))+m0),5e-5);
            m0 = m0 + foil.m(i);
        end

        if length(wake.m) == 2
            % Get only the wake points relevant to the view area
            i = (wake.xo < xmax+ch/2) & (wake.yo > ymin-ch/4);
            % Set lower wake evaluation points also offset from exact boundary
            j = find(i(1:N));
            j(end+1) = j(end) + 1; % add point that crosses outside domain
            wkl = [wake.xo(j) wake.yo(j)];
            bubble{1} = cat(1,bubble{1}(1,:) - wkl(1,:) + flipud(wkl), ...
                            bubble{1}(2:end-1,:), ...
                            bubble{1}(end,:) - wkl(1,:) + wkl);
            % Repeat for upper wake
            j = find(i(N+1:2*N)) + N;
            j(end+1) = j(end) + 1;
            wku = [wake.xo(j) wake.yo(j)];
            bubble{2} = cat(1,bubble{2}(1,:) - wku(1,:) + flipud(wku), ...
                            bubble{2}(2:end-1,:), ...
                            bubble{2}(end,:) - wku(1,:) + wku);
        end

        node = cat(1,node,bubble{:});

        % Create adjacency matrix to track edges
        M = false(size(node,1));
        M(1:4,1:4) = periodicAdjMatrix(4);
        len = 4;
        for i = 1:nSurf
            r = size(bubble{i},1);
            M(len+(1:r),len+(1:r)) = periodicAdjMatrix(r);
            len = len + r;
        end

        index = repmat((1:len).',1,len);
        edge = reshape(index(M),[2 len]).';

        % Mesh generation
        opts.kind = 'delfront';
        opts.rho2 = 1;
        [vert,etri,tria,tnum] = refine2(node,edge,[],opts);

        % Evaluate velocities at the mesh vertices
        field.xc = vert(:,1);
        field.yc = vert(:,2);
        field.dx = ones(size(field.xc));    % Grid points are the control points
        field.dy = zeros(size(field.xc));   % of the influence calculation.
        field.xo = field.xc - 0.5;          % They are placed at the center of
        field.yo = field.yc;                % flat panels with unit length.
        field.theta = zeros(size(field.xc));
        field.m = length(field.xc);

        [A1,B1] = influence(field,foil);
        u = B1*foil.G + 1;
        v = A1*foil.G;
        if length(wake.m) == 2
            [A2,B2] = influence(field,wake,0);
            u = u + B2*wake.G;
            v = v + A2*wake.G;
        end

        % Make contour plot
        if strcmpi(opts.mesh,'on')
            trisurf(tria,vert(:,1),vert(:,2),sqrt(u.^2+v.^2).','LineStyle','-','FaceColor','interp')
        else
            trisurf(tria,vert(:,1),vert(:,2),sqrt(u.^2+v.^2).','LineStyle','none','FaceColor','interp')
        end
        daspect([1 1 1])
        cm = Colormap;
        colormap(cm.(lower(opts.colormap))(400))
        caxis(max(abs(sqrt(u.^2+v.^2)-1))*[-1 1]+1) % center color data around freestream velocity
        xlim([xmin-ch/2 xmax+ch/2])
        ylim([ymin-ch/4 ymax+ch/4])
    end

    if strcmpi(opts.plot,'stream')
        % Draws streamlines in the flow
        % The grid must be meshgrid-like, but some optimizations were made
        x = [(cos(linspace(pi/2,pi/40,20))-1).'*ch/2+xmin; uniquetol(foil.xo,5.5e-3); (cos(linspace(39*pi/40,pi/2,20))+1).'*ch/2+xmax];
        y = [(cos(linspace(pi/2,pi/20,10))-1).'*ch/4+ymin; uniquetol(foil.yo,2.5e-2); (cos(linspace(19*pi/20,pi/2,10))+1).'*ch/4+ymax];
        [X,Y] = meshgrid(x,y);
        field.xc = X(:);
        field.yc = Y(:);
        field.dx = ones(size(field.xc));
        field.dy = zeros(size(field.xc));
        field.xo = field.xc - 0.5;
        field.yo = field.yc;
        field.theta = zeros(size(field.xc));
        field.m = length(field.xc);
        [A1,B1] = influence(field,foil);
        u = B1*foil.G + 1;
        v = A1*foil.G;
        if length(wake.m) == 2
            [A2,B2] = influence(field,wake,0);
            u = u + B2*wake.G;
            v = v + A2*wake.G;
        end
        for i = 1:nSurf
            in = inpolygon(X,Y,surfaces{i}(:,1),surfaces{i}(:,2));
            u(in) = 0;
            v(in) = 0;
        end
        lineobj = streamline(X, Y, reshape(u,size(X)), reshape(v,size(X)), ...
            [zeros(1,59) (0:19)*ch/40] + x(1), ...
            [(59:-1:1)/60*(y(end)-y(1)) zeros(1,20)] + y(1));
        if length(wake.m) == 2
            cnt = 0;
            for i = 1:length(lineobj)
                % Identify and set style for streamlines in the jet flow
                XYData = get(lineobj(i),{'XData','YData'});
                [in,on] = inpolygon(XYData{1}(end),XYData{2}(end), ...
                    wake.xo([(1:N) (2*N:-1:N+1)]),wake.yo([(1:N) (2*N:-1:N+1)]));
                if in || on
                    set(lineobj(i),'Color',[0.8353 0.3686 0],'LineStyle','--');
                    cnt = cnt + 1;
                    wline = i;
                end
            end
            % Reduce density of streamlines outside of the jet streamtube
            for i = wline-cnt-1:-2:1
                set(lineobj(i),'LineStyle','none');
            end
            for i = wline+2:2:length(lineobj)
                set(lineobj(i),'LineStyle','none');
            end
        end
        % Plot airfoil surfaces on top
        for i = 1:nSurf
            fill(surfaces{i}(:,1),surfaces{i}(:,2),[1 1 1]*0.98)
        end
        if strcmpi(opts.mesh,'on'), plot(X,Y,'r.'); end
        daspect([1 1 1])
        xlim([xmin-ch/2 xmax+ch/2])
        ylim([ymin-ch/4 ymax+ch/4])
    end
end

function newBoundary = offsetBoundary(inputBoundary,angles,offset)
    tnormal = mod(angles+pi/2,2*pi); % angles of the panel normal vectors
    tnormbi = 0.5*(tnormal(1:end-1) + tnormal(2:end)); % normal bisector
    movedir = zeros(length(tnormal)+1,2);
    movedir(2:end-1,:) = [cos(tnormbi) sin(tnormbi)];
    movedir(1,:) = [cos(tnormal(1)) sin(tnormal(1))];
    movedir(end,:) = [cos(tnormal(end)) sin(tnormal(end))];
    newBoundary = inputBoundary + offset*movedir; % expand points outward
end

function M = periodicAdjMatrix(n)
    % PERIODICADJMATRIX helper function for generating the adjacency
    % matrix for a n-segment periodic curve.
    M = logical(eye(n) + diag(ones(1,n-1),1));
    M(n,1) = true;
end