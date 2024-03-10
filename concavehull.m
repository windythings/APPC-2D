function test_hull = concavehull(points,kstart)
    % Computes the concave hull of a 2D point cloud
    % s = concavehull(points) Returns hull from points, specified as [N,2] scalar
    % s = concavehull(points,kstart) Starts with user-specified k-neighbors
    %
    % Adapted from Python implementation by Joao Paolo Figueira published at
    % https://towardsdatascience.com/the-concave-hull-c649795c0f0f
    %
    % Method for finding line intersections from
    % https://www.mathworks.com/matlabcentral/fileexchange/27205-fast-line-segment-intersection
    %
    % Algorithm by Adriano Moreira and Maribel Yasmina Santos
    % Moreira, A. J. C., and Santos, M. Y., Concave hull: A k-nearest neighbours approach
    % for the computation of the region occupied by a set of points, 2007.
    %
    % Author: Alan S. Hong
    if nargin == 1
        kstart = 3;
    end
    % Remove duplicate points
    [~,w] = unique(points,'stable','rows');
    hull.data_set = [points(w,1) points(w,2)];
    if size(hull.data_set,1) == 3
        return
    end
    % Increment k using a prime table schedule
    prime_k = [3 5 7 11 13 17 21 23 29 31 37 41 43 47 53 59 61 67 71 73 79 83 89 97];
    for k = [kstart prime_k(prime_k > kstart)]
        % Create the initial index
        hull.indices = true(size(hull.data_set,1),1);

        % Make sure no. of neighbors does not exceed available points
        kk = min(k,size(hull.data_set,1));

        [~,first_point] = min(hull.data_set(:,2));
        current_point = first_point;

        % Initialize hull
        hull.hull = hull.data_set(first_point,:);
        test_hull = hull.hull;

        % Remove the first point
        hull.indices(first_point) = false;

        prev_angle = 0;
        step = 2;
        stop = 2 + kk;
        while ((current_point ~= first_point) | (step == 2)) && length(hull.indices(hull.indices)) > 0
            if step == stop
                hull.indices(first_point) = true;
            end

            knn = get_k_nearest(hull,current_point,kk);

            % Calculate the angles between current_point and the knn points
            angles = calculate_angles(hull,current_point,knn,prev_angle);

            % Get descending sorted indices of angles
            [~,candidates] = sort(-angles);

            i = 1;
            invalid_hull = true;
            while invalid_hull && i < length(candidates)
                candidate = candidates(i);

                % Create a test hull to check if there are any self-intersections
                next_point = hull.data_set(knn(candidate),:);
                test_hull = [hull.hull; next_point];

                invalid_hull = ~simple_polygon([test_hull(1:end-2,:) test_hull(2:end-1,:)],[test_hull(end-1,:) test_hull(end,:)]);
                i = i + 1;
            end
            if invalid_hull
                break
            end
            prev_angle = calculate_angles(hull,knn(candidate),current_point,0);
            current_point = knn(candidate);
            hull.hull = test_hull;

            hull.indices(current_point) = false;
            step = step + 1;
        end
        if invalid_hull
            continue
        end
        % Check if all bounds are bounded by the hull
        [in,on] = inpolygon(hull.data_set(:,1),hull.data_set(:,2),hull.hull(:,1),hull.hull(:,2));
        count = length(find(in | on));
        if count == size(hull.data_set,1)
            return
        else
            continue
        end
    end

    function knn = get_k_nearest(self,ix,k)
        % Calculates the k nearest point indices to the point indexed by ix
        ixs = self.indices;
        base_indices = find(ixs);
        distances = sqrt(sum((self.data_set(ix,:) - self.data_set(ixs,:)).^2,2));
        [~,sorted_indices] = sort(distances);
        kk = min(k,length(sorted_indices));
        k_nearest = sorted_indices(1:kk);
        knn = base_indices(k_nearest);
    end

    function angles = calculate_angles(self,ix,ixs,ref_angle)
        % Calculates the angles between from a source point indexed by ix
        % and a set of target points indexed at ixs
        r_ix = self.data_set(ix,:);
        r_ixs = self.data_set(ixs,:);
        angles = atan2(r_ixs(:,2)-r_ix(2),r_ixs(:,1)-r_ix(1)) - ref_angle;
        % Clean up angles to be in the range [0,2*pi)
        angles = mod(angles + 2*pi,2*pi);
        angles(angles<0) = angles(angles<0) + 2*pi;
    end

    function issimple = simple_polygon(XY1,XY2)
        % Checks if adding the new segment to the hull results in a simple polygon
        % Implementation of finding intersections by U. Murat Erdem
        [n_rows_1,n_cols_1] = size(XY1);
        [n_rows_2,n_cols_2] = size(XY2);
        X1 = repmat(XY1(:,1),1,n_rows_2);
        X2 = repmat(XY1(:,3),1,n_rows_2);
        Y1 = repmat(XY1(:,2),1,n_rows_2);
        Y2 = repmat(XY1(:,4),1,n_rows_2);
        XY2 = XY2.';
        X3 = repmat(XY2(1,:),n_rows_1,1);
        X4 = repmat(XY2(3,:),n_rows_1,1);
        Y3 = repmat(XY2(2,:),n_rows_1,1);
        Y4 = repmat(XY2(4,:),n_rows_1,1);
        X4_X3 = (X4-X3);
        Y1_Y3 = (Y1-Y3);
        Y4_Y3 = (Y4-Y3);
        X1_X3 = (X1-X3);
        X2_X1 = (X2-X1);
        Y2_Y1 = (Y2-Y1);
        numerator_a = X4_X3 .* Y1_Y3 - Y4_Y3 .* X1_X3;
        numerator_b = X2_X1 .* Y1_Y3 - Y2_Y1 .* X1_X3;
        denominator = Y4_Y3 .* X2_X1 - X4_X3 .* Y2_Y1;
        u_a = numerator_a ./ denominator;
        u_b = numerator_b ./ denominator;
        INT_B = (u_a > 0) & (u_a < 1) & (u_b > 0) & (u_b < 1);
        issimple = ~any(INT_B(:));
    end
end