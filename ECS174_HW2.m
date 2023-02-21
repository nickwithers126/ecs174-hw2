% creates panorama given two images of the same scene from different
% angles by estimating a homography

%% --- 1 --- 
% read in images
im_left_color = imread('uttower_left.jpg');
im_right_color = imread('uttower_right.jpg');

% converting to double and grayscale
im_left = rgb2gray(im2double(im_left_color));
im_right = rgb2gray(im2double(im_right_color));


%% --- 2 --- 
% detecting feature points in both images
% Usage:  [cim, r, c] = harris(im, sigma, thresh, radius, disp)
[imc_left, row_l, col_l] = harris(im_left, 3, 0.03, 3, 1);
[imc_right, row_r, col_r] = harris(im_right, 3, 0.03, 3, 1);


%% --- 3 --- 
% creating descripors by extracting neighboring pixels and flattening
% usage: d = get_descriptors(img, row, col, patch_size)
desc_left = get_descriptors(im_left, row_l, col_l, 10); % patch size = 10 -> 21x21 descriptor
desc_right = get_descriptors(im_right, row_r, col_r, 10);


%% --- 4 ---
% normalizing descriptors to have zero mean and unit standard deviation
desc_left = reshape(zscore(desc_left(:)),size(desc_left,1),size(desc_left,2));
desc_right = reshape(zscore(desc_right(:)),size(desc_right,1),size(desc_right,2));

% get matrix of distances between all corners
% distances will be n x m matrix where
% n is # or corners in left img and m is # of corners in right img
distances = dist2(desc_left, desc_right);    


%% --- 5 ---
% finding best matches of our corners according to descriptor distance
matches = pair_matches(distances, row_l, col_l, row_r, col_r);


%% --- 6 ---
% running ransac function to get homography
[num_inliers, av_residual, locations, H] = myRANSAC(matches);
fprintf("number of inliers: %d\n", num_inliers);
fprintf("average residual: %f\n", av_residual);
fprintf("locations of inliers: [displayed -> left_row left_col right_row, right_col]\n");
disp(locations);

% graphing inliers in both images
sideBySide = cat(2, im_left, im_right); %putting images side by side so we can draw a line
figure(3), imshow(sideBySide);
for i=1:size(locations, 1)
   figure(3), line([locations(i,2) locations(i,4)+1024], [locations(i,1) locations(i,3)], 'Color', 'g', 'LineWidth', 1);
end
figure(3),title('RANSAC Inliers');


%% --- 7 ---
% make transformation based on homography computer by ransac 
T = maketform('projective', H);
    
[~, xdata, ydata] = imtransform(im_left_color, T);
[row, col] = size(im_right);
% getting boundaries for large canvas to fit both images when combining
left_border = min(col, xdata(1)); % first col coordinate
right_border = max(col, xdata(2)); % last col coordinate
bottom_border = min(row, ydata(1)); % first row coordinate
top_border = max(row, ydata(2)); % last row coordinate

% must make a transformation for right image so it is same size
% will not change dimensions or warp, only place in correct spot
right_transformation = maketform('affine', [1 0 0; 0 1 0; 0 0 1]);

% transforming both images onto our large canvas with dimensions calculated above
transform_left = imtransform(im_left_color, T, 'XYScale', 1, 'XData', [left_border right_border], 'YData', [bottom_border top_border]);
transform_right = imtransform(im_right_color, right_transformation, 'XYScale', 1, 'XData', [left_border right_border],'YData', [bottom_border top_border]);
 

%% --- 8 --- 
% loop through pixels of the whole canvas and set accordingly
% if left img has pixel and right does not -> left pixel value
% if right img has pixel and left does not -> right pixel value
% if both have pixel values -> average them
result = transform_left;
[m, n, ~] = size(transform_left); % m = rows; n = cols
    
for i = 1:3 % loop through color channels
    for j=1:m % loop through rows
        for k=1:n % loop through cols
            if transform_left(j,k,i) == 0 && transform_right(j,k,i) ~= 0
                result(j,k,i) = transform_right(j,k,i);
            elseif transform_left(j,k,i) ~= 0 && transform_right(j,k,i) == 0
                result(j,k,i) = transform_left(j,k,i);
            else % if overlap, average values from both images
                result(j,k,i) = (transform_left(j,k,i)/2) + (transform_right(j,k,i)/2);
            end
        end  
    end
end

% show result
figure(4),imshow(result);
figure(4),title('Output Panorama');
