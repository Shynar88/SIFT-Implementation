% Local Feature Stencil Code
% Written by James Hays for CS 143 @ Brown / CS 4476/6476 @ Georgia Tech

% This function is provided for development and debugging but cannot be
% used in the final handin. It 'cheats' by generating interest points from
% known correspondences. It will only work for the three image pairs with
% known correspondences.

% NOTE: These feature points are _subpixel_ precise.
% By default, we round these interest points, but in principle, you can
% interpolate the image to extract descriptors at subpixel locations.

% 'eval_file' is the file path to the list of known correspondences.
% 'scale_factor' is needed to map from the original image coordinates to
%   the resolution being used for the current experiment.

% 'x1' and 'y1' are nx1 vectors of x and y coordinates of interest points
%   in the first image.
% 'x1' and 'y1' are mx1 vectors of x and y coordinates of interest points
%   in the second image. For convenience, n will equal m but don't expect
%   that to be the case when interest points are created independently per
%   image.
function [x1, y1, x2, y2] = cheat_interest_points(eval_file, scale_factor, image1, image2, feature_width)

load(eval_file);

x1 = round( x1 .* scale_factor );
y1 = round( y1 .* scale_factor );
x2 = round( x2 .* scale_factor );
y2 = round( y2 .* scale_factor );

% Check bounds
[m1,n1] = size( image1 );
[m2,n2] = size( image2 );
fw2 = feature_width/2;

ind1 = x1 - fw2 <= 0 | x1 + fw2 > m1 | y1 - fw2 <=0 | y1 + fw2 > n1;
x1( ind1 ) = [];
y1( ind1 ) = [];
x2( ind1 ) = [];
y2( ind1 ) = [];

ind2 = x2 - fw2 <= 0 | x2 + fw2 > m2 | y2 - fw2 <=0 | y2 + fw2 > n2;
x1( ind2 ) = [];
y1( ind2 ) = [];
x2( ind2 ) = [];
y2( ind2 ) = [];

end

