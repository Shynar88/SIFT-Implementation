% Written by James Hays and James Tompkin for CS 143 @ Brown / CS 4476/6476 @ Georgia Tech

% This script 
% - Loads and resizes images
% - Computes correspondence
% - Visualizes the matches
% - Evaluates the matches based on ground truth correspondences

function [accMPEND,accMPEMR,accMPEEG,accMPEAvg] = projLF()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup variables
% Visualize the results
% We will evaluate your code with this set to '[]' (no vis).
visualize = 'dots';

% Amount by which to resize images for speed.
% Feel free to experiment with this for debugging, 
% but we will evaluate your code with this set to 0.5.
scale_factor = 0.5;

% Width and height of the descriptor window around each local feature, in image pixels.
% In SIFT, this is 16 pixels.
% Feel free to experiment with this for debugging or extra credit, 
% but we will evaluate your code with this set to 16.
descriptor_window_image_width = 16;

% Number of points to evaluate for accuracy
% We will evaluate your code on the first 100 matches you return.
maxPtsEval = 100;

% Whether to use the 'cheat' hand-picked interest points
cheatInterestPoints = false; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Notre Dame de Paris
% Easiest
disp( 'Notre Dame de Paris' );
image1 = imread('../data/NotreDame/921919841_a30df938f2_o.jpg');
image2 = imread('../data/NotreDame/4191453057_c86028ce1f_o.jpg');
eval_file = '../data/NotreDame/921919841_a30df938f2_o_to_4191453057_c86028ce1f_o.mat';

image1 = im2single( imresize( rgb2gray(image1), scale_factor, 'bilinear') );
image2 = im2single( imresize( rgb2gray(image2), scale_factor, 'bilinear') );

% Task: implement the following three fuctions
% 1) Find distinctive interest points in each image. Szeliski 4.1.1
if ~cheatInterestPoints
    [x1, y1] = get_interest_points(image1, descriptor_window_image_width);
    [x2, y2] = get_interest_points(image2, descriptor_window_image_width);
else
    % Use cheat_interest_points only for development and debugging!
    [x1, y1, x2, y2] = cheat_interest_points(eval_file, scale_factor, image1, image2, descriptor_window_image_width);
end

% 2) Create feature descriptors at each interest point. Szeliski 4.1.2
[image1_features] = get_descriptors(image1, x1, y1, descriptor_window_image_width);
[image2_features] = get_descriptors(image2, x2, y2, descriptor_window_image_width);

% 3) Match features. Szeliski 4.1.3
[matches, confidences] = match_features(image1_features, image2_features);

% Evaluate matches
[~,~,~,accMPEND] = evaluate_correspondence(image1, image2, eval_file, scale_factor, ... 
                        x1, y1, x2, y2, matches, confidences, ...
                        maxPtsEval, visualize, 'eval_ND.png' );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mount Rushmore
% A little harder than Notre Dame
disp( 'Mount Rushmore' );
image1 = imread('../data/MountRushmore/9021235130_7c2acd9554_o.jpg');
image2 = imread('../data/MountRushmore/9318872612_a255c874fb_o.jpg');
eval_file = '../data/MountRushmore/9021235130_7c2acd9554_o_to_9318872612_a255c874fb_o.mat';

image1 = im2single( imresize( rgb2gray(image1), scale_factor, 'bilinear') );
image2 = im2single( imresize( rgb2gray(image2), scale_factor, 'bilinear') );

% Task: implement the following three fuctions
% 1) Find distinctive interest points in each image. Szeliski 4.1.1
if ~cheatInterestPoints
    [x1, y1] = get_interest_points(image1, descriptor_window_image_width);
    [x2, y2] = get_interest_points(image2, descriptor_window_image_width);
else
    % Use cheat_interest_points only for development and debugging!
    [x1, y1, x2, y2] = cheat_interest_points(eval_file, scale_factor, image1, image2, descriptor_window_image_width);
end

% 2) Create feature descriptors at each interest point. Szeliski 4.1.2
[image1_features] = get_descriptors(image1, x1, y1, descriptor_window_image_width);
[image2_features] = get_descriptors(image2, x2, y2, descriptor_window_image_width);

% 3) Match features. Szeliski 4.1.3
[matches, confidences] = match_features(image1_features, image2_features);

% Evaluate matches
[~,~,~,accMPEMR] = evaluate_correspondence(image1, image2, eval_file, scale_factor, ... 
                        x1, y1, x2, y2, matches, confidences, ...
                        maxPtsEval, visualize, 'eval_MR.png' );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gaudi's Episcopal Palace
% This pair is difficult
disp( 'Gaudi''s Episcopal Palace' );
image1 = imread('../data/EpiscopalGaudi/4386465943_8cf9776378_o.jpg');
image2 = imread('../data/EpiscopalGaudi/3743214471_1b5bbfda98_o.jpg');
eval_file = '../data/EpiscopalGaudi/4386465943_8cf9776378_o_to_3743214471_1b5bbfda98_o.mat';

image1 = im2single( imresize( rgb2gray(image1), scale_factor, 'bilinear') );
image2 = im2single( imresize( rgb2gray(image2), scale_factor, 'bilinear') );

% Task: implement the following three fuctions
% 1) Find distinctive interest points in each image. Szeliski 4.1.1
if ~cheatInterestPoints
    [x1, y1] = get_interest_points(image1, descriptor_window_image_width);
    [x2, y2] = get_interest_points(image2, descriptor_window_image_width);
else
    % Use cheat_interest_points only for development and debugging!
    [x1, y1, x2, y2] = cheat_interest_points(eval_file, scale_factor, image1, image2, descriptor_window_image_width);
end

% 2) Create feature descriptors at each interest point. Szeliski 4.1.2
[image1_features] = get_descriptors(image1, x1, y1, descriptor_window_image_width);
[image2_features] = get_descriptors(image2, x2, y2, descriptor_window_image_width);

% 3) Match feature descriptors. Szeliski 4.1.3
[matches, confidences] = match_features(image1_features, image2_features);

% Evaluate matches
[~,~,~,accMPEEG] = evaluate_correspondence(image1, image2, eval_file, scale_factor, ... 
                        x1, y1, x2, y2, matches, confidences, ...
                        maxPtsEval, visualize, 'eval_EG.png' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute average accuracy
accMPEAvg = (accMPEND + accMPEMR + accMPEEG) / 3;

end