addpath("stbl")
% checking stability of probabilistic representation of
% the exp(-t\lambda**\alpha)
test_passed = 0;
tests = 500;
nsims = 1e5;
errorvec = zeros(tests,1);
for test = 1:tests
    alpha = rand()*0.9 + 0.1;
    skewness = cos(alpha*pi/2)^(1/alpha);
    samples = zeros(nsims,1);
    lambda = rand(); t = rand();
    lambdatalpha = t^(1/alpha)*lambda;
    for sim = 1:nsims
        tau = stblrnd(alpha,1,skewness,0);
        samples(sim) = exp(-tau*lambdatalpha);
    end
    error_vec(test) = abs(mean(samples)-exp(-t*lambda^alpha));
    boolean = abs(mean(samples)-exp(-t*lambda^alpha))<= 1.96*std(samples)/sqrt(nsims);
    test_passed = test_passed + boolean;
end

test_passed
tests
test_passed >= 0.94*tests