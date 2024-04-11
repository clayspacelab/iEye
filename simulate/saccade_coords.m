function [isacc_coords,fsacc_coords, real_error] = saccade_coords(correct_coords,p)
% Generates saccade
ierr_std = [50, pi/18];
ferr_std = [10, pi/25];
correct_theta = atan2((correct_coords(:,2) - p.yC) , ((correct_coords(:,1) - p.xC)+0.00001));
correct_r = sqrt((correct_coords(:,2) - p.yC).^2 + (correct_coords(:,1) - p.xC).^2);
ierrored_r = correct_r + ierr_std(1) .* randn(p.nTrials, 1);
ferrored_r = correct_r + ferr_std(1) .* randn(p.nTrials, 1);
ierrored_theta = correct_theta + ierr_std(2) .* randn(p.nTrials, 1);
ferrored_theta = correct_theta + ferr_std(2) .* randn(p.nTrials, 1);
ierrored_x = p.xC + ierrored_r .* cos(ierrored_theta);
ierrored_y = p.yC + ierrored_r .* sin(ierrored_theta);
isacc_coords = [ierrored_x, ierrored_y];
ferrored_x = p.xC + ferrored_r .* cos(ferrored_theta);
ferrored_y = p.yC + ferrored_r .* sin(ferrored_theta);
fsacc_coords = [ferrored_x, ferrored_y];
real_error = {};
real_error.isacc_xy = sqrt((isacc_coords(:, 1)-correct_coords(:, 1)).^2 + ...
    (isacc_coords(:, 2)-correct_coords(:, 2)).^2);
real_error.isacc_r = sqrt((ierrored_r-correct_r).^2);
real_error.isacc_theta = sqrt((ierrored_theta-correct_theta).^2);
real_error.fsacc_xy = sqrt((fsacc_coords(:, 1)-correct_coords(:, 1)).^2 + ...
    (fsacc_coords(:, 2)-correct_coords(:, 2)).^2);
real_error.fsacc_r = sqrt((ferrored_r-correct_r).^2);
real_error.fsacc_theta = sqrt((ferrored_theta-correct_theta).^2);
end

