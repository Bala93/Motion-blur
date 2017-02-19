function motionblur()

load cam.csv
img_re = imread('pillars.jpg');
img    = rgb2gray(img_re);
cam1   = cam;

% In plane translation
% cam1 (:,3:end) = 0;
% Out plane rotation
% cam1(:,1:3) = 0;
% cam1(:,end) = 0; 

weightage1 = ones(40,1);
[cam2,weightage2] = get_unique(cam);

% Calculate motion blur
motion_blurred_b1 = uint8(calc_motion_blur(img,cam1,weightage1,1) /40);
motion_blurred_b2 = uint8(calc_motion_blur(img,cam2,weightage2,2)/40);

% Calculating square error
square_error      = min_square_error(motion_blurred_b1,motion_blurred_b2);
disp(square_error);

figure,imshow(motion_blurred_b1);
imwrite(motion_blurred_b1,'motion_blur_1.jpg');

figure,imshow(motion_blurred_b2);
imwrite(motion_blurred_b2,'motion_blur_2.jpg');

end

function [motion_blurred] = calc_motion_blur(img,cam,weightage,choice)
motion_blurred = zeros(size(img));
n = size(cam,1);
% Homography matrix parameters.
% f - focal length
% K - Diagonal matrix with f along the diagonal
% n_ - represents the normal vector coming out of plane
% d - distance between the scene and the camera.
f = 500;
K  = [ f 0 0 ; 0 f 0 ; 0 0 1];
n_ = [0 0 1]';
d = 1000;
file_name = sprintf('motion_blur_%d.gif',choice);
imwrite(uint8(img),file_name,'gif');
for i = 1:n
    tx = cam(i,1);
    ty = cam(i,2);
    tz = cam(i,3);
    rx = cam(i,4);
    ry = cam(i,5);
    rz = cam(i,6);
    % Translational matrix
    t  = [tx;ty;tz];
    % Rotational matrix
    Rx = [ 1 0 0 ; 0 cos(rx) -sin(rx);0 sin(rx) cos(rx)];
    Ry = [ cos(ry) 0 sin(ry) ; 0 1 0 ; -sin(ry) 0 cos(ry)];
    Rz = [ cos(rz) -sin(rz) 0 ; sin(rz) cos(rz) 0 ; 0 0 1];
    % Creating homography
    H = homography(t,Rx,Ry,Rz,K,n_,d);
    % Applying homography to image
    img_after_h = apply_homo(img,H);
    imwrite(uint8(img_after_h),file_name,'gif','DelayTime',0.01,'WriteMode','append');
    motion_blurred = motion_blurred + (weightage(i) * img_after_h);
end
end



function [H] = homography(t,Rx,Ry,Rz,K,n,d)
    R = Rx * Ry * Rz;
    H = K*(R + (t*n')/d)*inv(K);
end

function [img_after_h] = apply_homo(img,H)
  [m,n] = size(img);
  img_after_h   = zeros(m,n);
  H = inv(H);
     for j = 1:m
        for i = 1:n
            % x,y,1 is multiplied with H and is normalized by tmp.
            tmp = H*[i;j;1];
            i1  = tmp(1) / tmp(3);
            j1  = tmp(2) / tmp(3);
            % Since inverse mapping is done bilinearInterp is used.
            v1  = BilinearInterp(i1,j1,img);
            img_after_h(j,i) = v1;
        end
     end
end  

function [intrep] = BilinearInterp(i,j,I)
    % Taking the neigbhorhood into consideration so avoid holes in image. 
    del_x = i - floor(i);
    del_y = j - floor(j);
    i = floor(i);
    j = floor(j);
    [m,n] = size(I);
    if i <= 0 || j <= 0 || i >= n || j >=m
        intrep = 0;
        return
    end
    intrep = (1-del_x)*(1-del_y)*I(j,i) + (1-del_x)*(del_y)*I(j,i+1) + (del_x)*(1-del_y)*I(j+1,i) + (del_x)*(del_y)*I(j+1,i+1);
end

function [unique_rows,weight] = get_unique(A)
    unique_rows = unique(A,'rows');
    weight = zeros(size(unique_rows,1),1);
    for i = 1:size(unique_rows)
        weight(i,1) = sum(ismember(A,unique_rows(i,:),'rows'));
    end
end


function [er] = min_square_error(b1,b2)
    no_pixels = size(b1(:),1);
    squ_diff =  (b1 - b2).^2;
    er       = sum(squ_diff(:))/no_pixels;
end
