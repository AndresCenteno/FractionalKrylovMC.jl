% trying stblrnd a lot of times
% this test passes
addpath("stbl")
trials = 1e4;
test_passed = 0;
for i = 1:trials
    alpha = 0.1 + rand()*0.9;
    test_passed = test_passed + (0 < stblrnd(alpha,1,cos(alpha*pi/2)^(1/alpha),0));
    % I don't know if you realize but
    % S(alpha,1,1,0) is equal in distribution for all alpha
end
test_passed