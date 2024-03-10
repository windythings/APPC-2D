function hull = ashape(points,alpha)
% warning('off','MATLAB:delaunay:DupPtsDelaunayWarnId')
tri = delaunay(points);
p = reshape(points(tri,:),[size(tri,1) 3 2]);
a = vecnorm(diff(p(:,[1 2 3 1],:),1,2),2,3);
s = sum(a,2)/2;
area = sqrt(s.*prod(s - a,2));
shp = sort(tri(prod(a,2)./(4 * area)<alpha,:),2);

edges = reshape(shp(:,[1 2 1 2 3 3]),[3*size(shp,1) 2]);
[C,ia,ic] = unique(edges,'rows');
i = ic==ic.';
i(1:length(ic)+1:numel(i)) = false;
C = C(~any(i(:,ia),1),:);

n = size(C,1);
i = C(:)==C(:).';
i(1:numel(C)+1:numel(i)) = false;
i = reshape(mod(find(i)-1,numel(C)) + 1,n,2);

T = zeros(n+1,1);
curr_pos = i(1);
T(1) = C(curr_pos);
for j = 2:n+1
    curr_pos = i(mod(curr_pos-1,n)+1,(curr_pos<=n)+1);
    T(j) = C(curr_pos);
end

hull = points(T,:);