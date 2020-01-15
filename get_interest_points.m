% Local Feature Stencil Code
% Written by James Hays for CS 143 @ Brown / CS 4476/6476 @ Georgia Tech

% Returns a set of interest points for the input image

% 'image' can be grayscale or color, your choice.
% 'descriptor_window_image_width', in pixels.
%   This is the local feature descriptor width. It might be useful in this function to 
%   (a) suppress boundary interest points (where a feature wouldn't fit entirely in the image, anyway), or
%   (b) scale the image filters being used. 
% Or you can ignore it.

% 'x' and 'y' are nx1 vectors of x and y coordinates of interest points.
% 'confidence' is an nx1 vector indicating the strength of the interest
%   point. You might use this later or not.
% 'scale' and 'orientation' are nx1 vectors indicating the scale and
%   orientation of each interest point. These are OPTIONAL. By default you
%   do not need to make scale and orientation invariant local features.
function [x, y, confidence, scale, orientation] = get_interest_points(image, descriptor_window_image_width)

% Implement the Harris corner detector (See Szeliski 4.1.1) to start with.
% You can create additional interest point detector functions (e.g. MSER)
% for extra credit.

% If you're finding spurious interest point detections near the boundaries,
% it is safe to simply suppress the gradients / corners near the edges of
% the image.

% The lecture slides and textbook are a bit vague on how to do the
% non-maximum suppression once you've thresholded the cornerness score.
% You are free to experiment. Here are some helpful functions:
%  BWLABEL and the newer BWCONNCOMP will find connected components in 
% thresholded binary image. You could, for instance, take the maximum value
% within each component.
%  COLFILT can be used to run a max() operator on each sliding window. You
% could use this to ensure that every interest point is at a local maximum
% of cornerness.
alpha = 0.06;
threshold = 0.0003;
larger_gaussian = fspecial('gaussian', 9, 1);
gaussian = fspecial('gaussian', 3, 0.5);
image = imfilter(image, gaussian);
Xderivative = [-1,0,1;-1,0,1;-1,0,1];
Yderivative = Xderivative';
Gx = imfilter(image,Xderivative);
Gy = imfilter(image,Yderivative);
Gxy = Gx.*Gy;
Gxy = imfilter(Gxy, larger_gaussian); 
Gxx = Gx.*Gx; 
Gxx = imfilter(Gxx, larger_gaussian); 
Gyy = Gy.*Gy;
Gyy = imfilter(Gyy, larger_gaussian);

harris_value = Gxx.*Gyy-Gxy.^2-alpha.*(Gxx + Gyy).*(Gxx + Gyy);

edges = zeros(size(image));
edges(descriptor_window_image_width+1 : end-descriptor_window_image_width, descriptor_window_image_width+1 : end-descriptor_window_image_width) = 1;
harris_value = harris_value.*edges;
thresholded = harris_value > threshold;
cc = bwconncomp(thresholded);
pixelIdxList = cc.PixelIdxList;
numObjects = cc.NumObjects;
x = zeros(numObjects, 1);
y = zeros(numObjects, 1);
for i = 1:numObjects
    idx = pixelIdxList{i};
    pixels = harris_value(idx);
    [~, maxs_idx] = max(pixels);
    x(i) = floor(idx(maxs_idx)/cc.ImageSize(1));
    y(i) = mod(idx(maxs_idx), cc.ImageSize(1));
end
% After computing interest points, here's roughly how many we return
% For each of the three image pairs
% - Notre Dame: ~1300 and ~1700
% - Mount Rushmore: ~3500 and ~4500
% - Episcopal Gaudi: ~1000 and ~9000

end

