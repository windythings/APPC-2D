function graf = pltC(surfaces,foil,wake,opts)
    graf = figure; hold on;

    nSurf = length(foil.m);
    N = wake.m(1);
    [xmin,xmax] = bounds(foil.xo);
    [ymin,ymax] = bounds(foil.yo);
    ch = xmax - xmin;

    if strcmpi(opts.plot,'contour')
        % Plots the velocity magnitude field on an efficient grid
        bubble = repmat({[]},1,nSurf); % the neighborhood of points near the solid surfaces
        product_info = ver('pde'); % check if the Matlab PDE toolbox is available
        if isempty(product_info)
            % Create a structured mesh using the package-free method
            SW = [xmin ymin];
            NE = [xmax ymax];
            m0 = 0;
            for i = 1:nSurf
                bubble{i} = offsetBoundary(surfaces{i},foil.theta((1:foil.m(i))+m0),0.0125*ch);
                % The resulting boundary can be self-intersecting
                % Generate a concave hull to guarantee a simple polygon
                bubble{i} = ashape(bubble{i},0.1);    % The alpha shape routine is faster than
                % bubble{i} = concavehull(bubble{i}); % k-neighbors concave hull but needs a
                SW = min(SW,min(bubble{i},[],1));     % user specified alpha radius that can
                NE = max(NE,max(bubble{i},[],1));     % differ between geometries
                m0 = m0 + foil.m(i);
            end

            % FINE mesh
            [X,Y] = meshgrid(linspace(SW(1),NE(1),301),linspace(SW(2),NE(2),301));
            inchk = false(size(X));
            for i = 1:nSurf
                in = inpolygon(X,Y,bubble{i}(:,1),bubble{i}(:,2));
                inchk(in) = true; % retain points that are inside the bubble
            end
            for i = 1:nSurf
                in = inpolygon(X(inchk),Y(inchk),surfaces{i}(:,1),surfaces{i}(:,2));
                inchk(inchk) = ~in; % omit points that are inside the solid surfaces
            end
            finex = X(inchk);
            finey = Y(inchk);

            % COARSE mesh
            [X,Y] = meshgrid(linspace(xmin-ch/2,xmax+ch/2,76),linspace(ymin-ch/4,ymax+ch/4,76));
            inchk = true(size(X));
            for i = 1:nSurf
                in = inpolygon(X,Y,bubble{i}(:,1),bubble{i}(:,2));
                inchk(in) = false;
            end
            coarx = X(inchk);
            coary = Y(inchk);

            % WAKE mesh
            wx = [linspace(wake.xo(1)+1e-4,xmax+ch/2,200) linspace(wake.xo(N+1)+1e-4,xmax+ch/2,200)].';
            wy = [interp1(wake.xo(1:N),wake.yo(1:N),wx(1:200))-1e-8; interp1(wake.xo(N+1:2*N),wake.yo(N+1:2*N),wx(201:400))+1e-8];
            wt = [interp1(wake.xo(1:N),wake.theta(1:N),wx(1:200)); interp1(wake.xo(N+1:2*N),wake.theta(N+1:2*N),wx(201:400))];
            i = wy > ymin-ch/4-1e-4;
            ep = (-1:1)*(NE(2)-SW(2))/300;
            wx = reshape(wx(i) + ep.*cos(wt(i) + pi/2),[],1);
            wy = reshape(wy(i) + ep.*sin(wt(i) + pi/2),[],1);

            % TOTAL mesh
            field.xc = [finex; coarx; wx];
            field.yc = [finey; coary; wy];
        else
            % Create unstructured mesh using Matlab PDE Toolbox
            % Place evaluation points barely outside the boundary because
            % exact boundary evaluation is unstable
            m0 = 0;
            for i = 1:nSurf
                bubble{i} = offsetBoundary(surfaces{i},foil.theta((1:foil.m(i))+m0),1e-6);
                m0 = m0 + foil.m(i);
            end

            % Get only the wake points relevant to the view area
            i = (wake.xo < xmax+ch/2) & (wake.yo > ymin-ch/4);

            % Set lower wake evaluation points also offset from exact boundary
            j = find(i(1:N));
            j(end+1) = j(end) + 1; % add point that crosses outside domain
            wkl = [wake.xo(j) wake.yo(j)];
            un = mod(wake.theta(j)+pi/2,2*pi); % angles of wake panel normal vectors
            un = [cos(un) sin(un)]; % convert angles to vectors
            wkl = [wkl-5e-6*un; flipud(wkl+5e-6*un)];
            wkl(1,:) = bubble{1}(1,:); % match the airfoil TE
            wkl(end,:) = bubble{1}(end,:);

            % Repeat for upper wake
            j = find(i(N+1:2*N)) + N;
            j(end+1) = j(end) + 1;
            wku = [wake.xo(j) wake.yo(j)];
            un = mod(wake.theta(j)+pi/2,2*pi);
            un = [cos(un) sin(un)];
            wku = [wku-5e-6*un; flipud(wku+5e-6*un)];
            wku(1,:) = bubble{2}(1,:);
            wku(end,:) = bubble{2}(end,:);

            % Build geometry matrix
            dim1 = max([foil.m+1 length(wkl) length(wku) 10]); % no. of rows containing point data
            gm(:,nSurf+3) = [3 4 xmax+ch/2 xmin([1 1])-ch/2 xmax+ch/2, ...
                          ymin([1 1])-ch/4 ymax([1 1])+ch/4 zeros(1,2*(dim1-4))].'; % domain box
            gm(:,nSurf+2) = [2; length(wku); wku(:); zeros(2*(dim1-length(wku)),1)]; % upper wake
            gm(:,nSurf+1) = [2; length(wkl); wkl(:); zeros(2*(dim1-length(wkl)),1)]; % lower wake
            % Assign W1, W2, and E1 name-spaces to lower and upper wakes
            % and box environment, respectively
            ns(nSurf+(1:3),:) = char('W1','W2','E1');
            m0 = 0;
            for i = 1:nSurf
                % Complete matrix with airfoils
                gm(:,i) = [2; foil.m(i)+1; bubble{i}(:); zeros(2*(dim1-foil.m(i)-1),1)];
                ns(i,:) = ['S' num2str(i)];
                m0 = m0 + foil.m(i);
            end
            sf = ['E1-W1-W2' sprintf('-S%d',1:nSurf)]; % set formula subtracts airfoils from the domain
            g = decsg(gm,sf,ns.');

            % Create mesh using the toolbox functions
            domain = createpde;
            geometryFromEdges(domain,g);
            msh = generateMesh(domain,'Hmin',1e-3); % Hmin may need to adapt to geometry
            field.xc = msh.Nodes(1,:).';
            field.yc = msh.Nodes(2,:).';
        end
        field.dx = ones(size(field.xc));    % Grid points are the control points
        field.dy = zeros(size(field.xc));   % of the influence calculation.
        field.xo = field.xc - 0.5;          % They are placed at the center of
        field.yo = field.yc;                % flat panels with unit length.
        field.theta = zeros(size(field.xc));
        field.m = length(field.xc);
        field.id = 'mesh';

        [A1,B1] = influence(field,foil);
        [A2,B2] = influence(field,wake);

        u = B1*foil.G + B2*wake.G + 1;
        v = A1*foil.G + A2*wake.G;

        % Make contour plot
        DT = delaunay(field.xc,field.yc);
        patch(field.xc(DT).',field.yc(DT).',sqrt(u(DT).^2+v(DT).^2).','LineStyle','none')
        for i = 1:nSurf
            fill(surfaces{i}(:,1),surfaces{i}(:,2),[1 1 1])
        end
        if strcmpi(opts.mesh,'on'), triplot(DT,field.xc,field.yc); end
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
        x = [(cos(linspace(pi/2,pi/40,20))-1).'*ch/2+xmin; uniquetol(foil.xo,5.5e-3); (cos(linspace(pi,21*pi/40,20))+1).'*ch/2+xmax];
        y = [(cos(linspace(pi/2,pi/20,10))-1).'*ch/4+ymin; uniquetol(foil.yo,2.5e-2); (cos(linspace(pi,11*pi/20,10))+1).'*ch/4+ymax];
        [X,Y] = meshgrid(x,y);
        field.xc = X(:);
        field.yc = Y(:);
        field.dx = ones(size(field.xc));
        field.dy = zeros(size(field.xc));
        field.xo = field.xc - 0.5;
        field.yo = field.yc;
        field.theta = zeros(size(field.xc));
        field.m = length(field.xc);
        field.id = 'mesh';
        [A1,B1] = influence(field,foil);
        [A2,B2] = influence(field,wake);
        u = B1*foil.G + B2*wake.G + 1;
        v = A1*foil.G + A2*wake.G;
        for i = 1:nSurf
            in = inpolygon(X,Y,surfaces{i}(:,1),surfaces{i}(:,2));
            u(in) = 0;
            v(in) = 0;
        end
        lineobj = streamline(X, Y, reshape(u,size(X)), reshape(v,size(X)), ...
                             [zeros(1,29) (0:9)*ch/20] + x(1), ...
                             [(1:29)/30*(y(end)-y(1)) zeros(1,10)] + y(1));
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