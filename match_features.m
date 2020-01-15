% Local Feature Stencil Code
% Written by James Hays for CS 143 @ Brown / CS 4476/6476 @ Georgia Tech

% Please implement the "nearest neighbor distance ratio test", 
% Equation 4.18 in Section 4.1.3 of Szeliski. 
% For extra credit you can implement spatial verification of matches.

%
% Please assign a confidence, else the evaluation function will not work.
%

% This function does not need to be symmetric (e.g., it can produce
% different numbers of matches depending on the order of the arguments).

% Input:
% 'features1' and 'features2' are the n x feature dimensionality matrices.
% num of features * 128. 1 feature is 128 vertor. sqrt(x1(2)-x1(1))^2 +...+ (x128(2))
% ..... for n dimensions 
%
% Output:
% 'matches' is a k x 2 matrix, where k is the number of matches. The first
%   column is an index in features1, the second column is an index in features2. 
%
% 'confidences' is a k x 1 matrix with a real valued confidence for every match.

function [matches, confidences] = match_features(features1, features2)

% Placeholder random matches and confidences.
% num_features = min(size(features1, 1), size(features2,1));
num_of_f1 = size(features1, 1);
num_of_f2 = size(features2, 1);
matches = [];
confidences = [];
threshold = 0.8;
distances = zeros(num_of_f1, num_of_f2);
for f1 = 1:num_of_f1
    for f2 = 1:num_of_f2
        euclidian_distance = sqrt(sum((features1(f1, :) - features2(f2, :)).^2));
        distances(f1, f2) = euclidian_distance;
    end
    [sorted_distances_for_this_f1, indexes_of_f2_for_corresponding_sorted_f1] = sort(distances(f1,:));
    NNDR = sorted_distances_for_this_f1(1)/sorted_distances_for_this_f1(2);
    if NNDR < threshold
        confidences = [confidences; 1/NNDR];
        matches = [matches; f1, indexes_of_f2_for_corresponding_sorted_f1(1)];
    end
end
    
% Remember that the NNDR test will return a number close to 1 for 
% feature points with similar distances.
% Think about how confidence relates to NNDR.