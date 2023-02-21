% returns array of descriptors
% rows - numbers of corners
% cols - size of flatten descriptor -> depends on patchsize
function d = get_descriptors(img, row, col, patch_size)

[num_corners,~] = size(col);

d=zeros(num_corners,(2*patch_size+1)^2); % blank descriptor

padded_img = padarray(img, [patch_size patch_size], 'symmetric');

for i=1:num_corners

    neighborhood = padded_img((row(i)+patch_size)-patch_size:(row(i)+patch_size)+patch_size, (col(i)+patch_size)-patch_size:(col(i)+patch_size)+patch_size);
    neighborhood = reshape(neighborhood,1,[]);
    d(i,:) = neighborhood;
    
end