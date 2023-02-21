function [num_inliers, av_residual, locations, H] = myRANSAC(matches)

[row, ~] = size(matches);

largest_set_inliers = [];
inliers_avg_residual = 0;

% run ransac for 500 loops looking for homography with most inliers
for j = 1:500

    % getting 4 random matches
    k = randperm(row);
    random_pairs = matches(k(1:4),:);

    x1=random_pairs(1,2); y1=random_pairs(1,1);
    x1p=random_pairs(1,4); yp1=random_pairs(1,3);
    x2=random_pairs(2,2); y2=random_pairs(2,1);
    x2p=random_pairs(2,4); yp2=random_pairs(2,3);
    x3=random_pairs(3,2); y3=random_pairs(3,1);
    x3p=random_pairs(3,4); yp3=random_pairs(3,3);
    x4=random_pairs(4,2); y4=random_pairs(4,1);
    x4p=random_pairs(4,4); yp4=random_pairs(4,3);

    A = [x1 y1 1 0 0 0 -1*x1*x1p -1*y1*x1p -1*x1p;
         0 0 0 x1 y1 1 -1*x1*yp1 -1*y1*yp1 -1*yp1;
         x2 y2 1 0 0 0 -1*x2*x2p -1*y2*x2p -1*x2p;
         0 0 0 x2 y2 1 -1*x2*yp2 -1*y2*yp2 -1*yp2;
         x3 y3 1 0 0 0 -1*x3*x3p -1*y3*x3p -1*x3p;
         0 0 0 x3 y3 1 -1*x3*yp3 -1*y3*yp3 -1*yp3;
         x4 y4 1 0 0 0 -1*x4*x4p -1*y4*x4p -1*x4p;
         0 0 0 x4 y4 1 -1*x4*yp4 -1*y4*yp4 -1*yp4];

    [~,~,V] = svd(A); 
    H=V(:,end);
    H=reshape(H,3,3);

    inliers = []; % store inliers for current homography in this matrix
    sum_residual = 0;

    % loop through all the matches and see if homography maps point in left
    % to point in right
    % if close (distance < epsilon), add to inliers
    for i = 1:row

        % transforming left to right - get left points
        left_x = matches(i,2); % col
        left_y = matches(i,1); % row
   
        projected_point = [left_x left_y 1] * H;
        result_col = projected_point(1,1)/projected_point(1,3); 
        result_row = projected_point(1,2)/projected_point(1,3);
        result_point = [result_row, result_col];

        % find distance between the point we calculated based off our
        % homography * point in left and the paired point in right img
        d = dist2(result_point, [matches(i,3), matches(i,4)]);

        % epsilon set here to 5
        % if distance less than 5 we count as inlier
        if (d < 5)
            inliers = [inliers; matches(i,:)];
            sum_residual = sum_residual + (d^2);
        end

    end

    % if current set of inliers is greater than out set of largest inliers,
    % we have to set the current set to the largest set
    [num_large, ~] = size(largest_set_inliers);
    [num_sample, ~] = size(inliers);
    if num_sample > num_large
        largest_set_inliers = inliers;
        inliers_avg_residual = sum_residual;
    end
  
end

% computing new homography with largest set of inliers
% should be more accurate
[row, ~] = size(largest_set_inliers);
num_inliers = row;

A = [];

for i=1:num_inliers

    xip = largest_set_inliers(i,4);
    yip = largest_set_inliers(i,3);
    xi = largest_set_inliers(i,2);
    yi = largest_set_inliers(i,1);

    A = [A; xi yi 1 0 0 0 (-1*xip*xi) (-1*xip*yi) (-1*xip)];
    A = [A; 0 0 0 xi yi 1 (-1*yip*xi) (-1*yip*yi) (-1*yip)];
    
end

[~, ~, V] = svd(A);
H = V(:,end);
H = reshape(H, 3, 3);

locations = largest_set_inliers; % locations of inlier points
av_residual = inliers_avg_residual/num_inliers; % avg residual of inliers











