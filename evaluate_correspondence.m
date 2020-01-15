% Local Feature Stencil Code
% Written by James Hays for CS 143 @ Brown / CS 4476/6476 @ Georgia Tech with Henry Hu <henryhu@gatech.edu>
% Edited by James Tompkin

% Your feature points should be unique within a local region,
% I.E., your detection non-maximal suppression function should work.
% 
% Look at the 'uniqueness test' for how we enforce this.
% It is intentionally simplistic and biased, 
% so make your detector _do the right thing_.

function [numGoodMatches,numBadMatches,accuracyAll,accuracyMaxEval] = evaluate_correspondence(imgA, imgB, ground_truth_correspondence_file, scale_factor, x1i, y1i, x2i, y2i, matches, confidences, maxPtsToEval, vismode, visfilename)

%% Sort matches by confidence
% Sort the matches so that the most confident onces are at the top of the
% list. You should not delete this, so that the evaluation
% functions can be run on the top matches easily.
[~, ind] = sort(confidences, 'descend');
matches = matches(ind,:);

x1_est = x1i(matches(:,1));
y1_est = y1i(matches(:,1));
x2_est = x2i(matches(:,2));
y2_est = y2i(matches(:,2));

x1_est = x1_est ./ scale_factor;
y1_est = y1_est ./ scale_factor;
x2_est = x2_est ./ scale_factor;
y2_est = y2_est ./ scale_factor;

good_matches = zeros(length(x1_est),1); %indicator vector

% Loads `ground truth' positions x1, y1, x2, y2
load(ground_truth_correspondence_file)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uniqueness test
%
x1_est_tmp = x1_est;
y1_est_tmp = y1_est;
x2_est_tmp = x2_est;
y2_est_tmp = y2_est;
uniquenessDist = 5;

% For each ground truth point
numPreMerge = length(x1_est);
for i=1:length(x1)
    % Compute distance of each estimated point to
    % the ground truth point
    x_dists = x1(i) - x1_est_tmp;
    y_dists = y1(i) - y1_est_tmp;
    dists = sqrt( x_dists.^2 + y_dists.^2 );
    toMerge = dists < uniquenessDist;
    
    if sum(toMerge) > 1
        % Do something to remove duplicates. Let's
        % average the coordinates of all points 
        % within 'uniquenessDist' pixels.
        % Also average the corresponded point (!)
        %
        % This part is simplistic, but a real-world
        % computer vision system would not know
        % which correspondences were good.
        avgX1 = mean( x1_est_tmp( toMerge ) );
        avgY1 = mean( y1_est_tmp( toMerge ) );
        avgX2 = mean( x2_est_tmp( toMerge ) );
        avgY2 = mean( y2_est_tmp( toMerge ) );

        x1_est_tmp( toMerge ) = [];
        y1_est_tmp( toMerge ) = [];
        x2_est_tmp( toMerge ) = [];
        y2_est_tmp( toMerge ) = [];

        % Add back point
        x1_est_tmp(end+1) = avgX1;
        y1_est_tmp(end+1) = avgY1;
        x2_est_tmp(end+1) = avgX2;
        y2_est_tmp(end+1) = avgY2;
    end
end
x1_est = x1_est_tmp;
y1_est = y1_est_tmp; 
x2_est = x2_est_tmp;
y2_est = y2_est_tmp;
numPostMerge = length(x1_est);
%
% Uniqueness test end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Distance test
for i = 1:length(x1_est)
    fprintf('( %4.0f, %4.0f) to ( %4.0f, %4.0f)', x1_est(i), y1_est(i), x2_est(i), y2_est(i));

    % For each x1_est, find nearest ground truth point in x1
    x_dists = x1_est(i) - x1;
    y_dists = y1_est(i) - y1;
    dists = sqrt(  x_dists.^2 + y_dists.^2 );
    
    [dists, best_matches] = sort(dists);
    
    current_offset = [x1_est(i) - x2_est(i), y1_est(i) - y2_est(i)];
    most_similar_offset = [x1(best_matches(1)) - x2(best_matches(1)), y1(best_matches(1)) - y2(best_matches(1))];
    
    match_dist = sqrt( sum((current_offset - most_similar_offset).^2));
    
    % A match is bad if there's no ground truth point within 150 pixels 
    % or
    % If nearest ground truth correspondence offset isn't within 40 pixels
    % of the estimated correspondence offset.
    fprintf(' g.t. point %4.0f px. Match error %4.0f px.', dists(1), match_dist);
    
    if(dists(1) > 150 || match_dist > 40)
        good_matches(i) = 0;
        fprintf('  incorrect\n');
    else
        good_matches(i) = 1;
        fprintf('  correct\n');
    end
end


numGoodMatches = sum(good_matches);
numBadMatches = length(x1_est) - sum(good_matches);
disp( ['Uniqueness: Pre-merge:    ' num2str(numPreMerge) '  Post-merge:  ' num2str(numPostMerge)] );
disp( ['Total:      Good matches: ' num2str(numGoodMatches) '  Bad matches: ' num2str(numBadMatches)] );

% For evaluation, we're going to use the word 'accuracy'.
% It's very difficult to count the number of actual correspondences in an
% image, so ideas of recall are tricky to apply.
% The second accuracy measure 'out of maxPtsToEval' captures some of that
% idea: if you return less than maxPtsToEval, your accuracy will drop.
% 
accuracyAll = (sum(good_matches) / length(x1_est)) * 100;
disp( ['Accuracy:  ' num2str(accuracyAll,4) '% (on all ' num2str(length(x1_est)) ' submitted matches)'] );

accuracyMaxEval = (sum(good_matches(1:min(length(good_matches),maxPtsToEval))) / maxPtsToEval) * 100;
disp( ['Accuracy:  ' num2str(accuracyMaxEval,4) '% (on first ' num2str(maxPtsToEval) ' matches sorted by decreasing confidence)'] );

% Visualize the result
if ~isempty(vismode)
    % You may also switch to a different visualization method, by passing
    % 'arrows' into show_correspondence instead of 'dots'.
    show_correspondence(imgA, imgB, ...
                        x1_est * scale_factor, y1_est * scale_factor, ...
                        x2_est * scale_factor, y2_est * scale_factor, ...
                        vismode, visfilename, good_matches);
end