function [X, Y, pup,real_error] = generate_XYPupil(correct_coords,s,sacc_arr,feedback_arr,p)
[isacc_coords,fsacc_coords,real_error] = saccade_coords(correct_coords,p);
dt = 1/p.ifg_freq;
t = dt:dt:s*dt;
X = p.xC;
Y = p.yC;
pup = 1200;

for freq = 0.03:0.03:1
    X = X + 10*exp(-freq) * sin(2*pi*freq*t+10*randn());
    Y = Y + 10*exp(-freq) * sin(2*pi*freq*t+10*randn());
    pup = pup + 50*exp(-freq) * sin(2*pi*freq*t+50*randn());
end

X = X + 2*randn(1, length(t));
Y = Y + 2*randn(1, length(t));
pup = pup + 20*randn(1, length(t));

isacc_smean = sacc_arr(:, 1) + ((sacc_arr(:, 2) - sacc_arr(:, 1))./8);
isacc_emean = sacc_arr(:, 1) + ((sacc_arr(:, 2) - sacc_arr(:, 1))./2);
fsacc_emean = (feedback_arr(:, 1)-10) + ((feedback_arr(:, 2) - (feedback_arr(:, 1)-10))./6);
corrsacc_emean = feedback_arr(:, 2);

isacc_start = round(isacc_smean + random('Exponential', 0, p.nTrials, 1));
isacc_end = round(isacc_emean + random('Exponential', 0, p.nTrials, 1));
fsacc_end = round(fsacc_emean + random('Exponential', 0, p.nTrials, 1));
corrsacc_end = round(corrsacc_emean + random('Exponential', 0, p.nTrials, 1));

for trial = 1:p.nTrials
    isacc_time = isacc_start(trial):isacc_end(trial);
    fsacc_time = isacc_end(trial):fsacc_end(trial);
    corrsacc_time = fsacc_end(trial):corrsacc_end(trial);
    X(isacc_time) = isacc_coords(trial, 1) + randn(length(isacc_time), 1);
    X(fsacc_time) = fsacc_coords(trial, 1) + randn(length(fsacc_time), 1);
    X(corrsacc_time) = correct_coords(trial, 1) + randn(length(corrsacc_time), 1);
    Y(isacc_time) = isacc_coords(trial, 2) + randn(length(isacc_time), 1);
    Y(fsacc_time) = fsacc_coords(trial, 2) + randn(length(fsacc_time), 1);
    Y(corrsacc_time) = correct_coords(trial, 2) + randn(length(corrsacc_time), 1);
end

end

