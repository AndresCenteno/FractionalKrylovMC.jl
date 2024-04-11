using Random, AlphaStableDistributions, Statistics

Random.seed!(1234)

test_passed = 0;
tests = 500;
nsims = Int(1e5);
error_vec = zeros(tests,1);
for test = 1:tests
    alpha = rand()*0.9 + 0.1;
    skewness = cos(alpha*pi/2)^(1/alpha);
    samples = zeros(nsims,1);
    lambda = rand(); t = rand();
    lambdatalpha = t^(1/alpha)*lambda;
    for sim = 1:nsims
        tau = rand(AlphaStable(alpha,1,skewness,0));
        samples[sim] = exp(-tau*lambdatalpha);
    end
    error_vec[test] = abs(mean(samples)-exp(-t*lambda^alpha));
    boolean = abs(mean(samples)-exp(-t*lambda^alpha))<= 1.96*std(samples)/sqrt(nsims);
    test_passed = test_passed + boolean;
end

@show test_passed
@show tests
@show test_passed >= 0.94*tests