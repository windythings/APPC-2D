function [VTX,TRI,uGrid,vGrid] = createTriData(surfaces,wake,chord)
% Create path to meshing tools
addpath('./mesh2d'); initmsh();

% Meshing settings
opts.kind = 'delfront';
opts.rho2 = 1;

% Flip orientation of one of the wakes for inner/outer wake evaluation
% consistency
[wake(1).pt1,wake(1).pt2] = deal(wake(1).pt2,wake(1).pt1);
wake(1).theta = wake(1).theta+pi - ceil(wake(1).theta/2/pi)*2*pi;

ductIdxShift = 5; % sample a few points in the duct before mesh interface

% Mesh powered jet region
node = [flipud(wake(1).co); surfaces(1).co(1:ductIdxShift+1,:); ...
        surfaces(2).co(end-ductIdxShift:end,:); wake(2).co];
node(node(:,1)>4*chord,:) = []; % trim domain
edge = (1:size(node,1)).' + [0 1]; edge(end,2) = 1;
[vert{1},~,tria{1},~] = refine2(node,edge,[],opts);
% Evaluate using +pi wake induction
[upow,vpow] = createStreamlines(vert{1}(:,1),vert{1}(:,2),surfaces,wake,pi);

% Mesh outer flow region
node = [flipud(wake(2).co); surfaces(2).co(1:end-ductIdxShift,:); ...
        surfaces(1).co(ductIdxShift+1:end,:); wake(1).co; ...
        4*chord chord;-chord chord;-chord -chord;4*chord -chord];
node(node(:,1)>4*chord,:) = [];
edge = (1:size(node,1)).' + [0 1]; edge(end,2) = 1;
if numel(surfaces) > 2 % add all other surfaces
    for j = 3:numel(surfaces)
        edge = [edge;(1:surfaces(j).m).'+[0 1]+size(node,1)];
        edge(end,2) = size(node,1) + 1;
        node = [node;surfaces(j).co];
    end
end
[vert{2},~,tria{2},~] = refine2(node,edge,[],opts);
% Evaluate using -pi wake induction
[uout,vout] = createStreamlines(vert{2}(:,1),vert{2}(:,2),surfaces,wake,-pi);

% Merge meshes
VTX = [vert{1};vert{2}]; TRI = [tria{1};tria{2}+size(vert{1},1)];
uGrid = [upow;uout]; vGrid = [vpow;vout];
