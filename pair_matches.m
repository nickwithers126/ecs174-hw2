% top 200 descriptor pairs with smallest pairwise distance
% output array of 200 pairs of points
function out = pair_matches(distances, row_l, col_l, row_r, col_r)

out = [];

for i = 1:200

    min_dist = min(min(distances)); % get smallest distance
    [r, c] = find(distances == min_dist); % find which corners
    next_small = mink(distances(r,:),2); % get 2 smallest for the point in left img
    next_small = next_small(1,2); % getting 2nd smallest

    if (min_dist/next_small) < 0.75 % make sure fetures are unique so check next closest
        out = [out; row_l(r), col_l(r), row_r(c), col_r(c)];
    end

    % set distances to infinity so we don't get them again
    distances(r,:) = inf;
    distances(:,c) = inf;
    
end
