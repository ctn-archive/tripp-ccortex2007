% function lambdae = elevatedRate(lambda, n, N, lambdat, p) returns the
% elevated rate required of n of N neurons, the rest of which fire at rate
% lambda (all Poisson) in order that the grand mean rate rarely drops
% below lambdat. The probability 0<p<1 defines rarely.

function lambdae = elevatedRate(lambda, n, N, lambdat, p)

    z = norminv(1-p, 0, 1);
    
    %Newton's method
    guess = (2*(N*lambda)^.5 + N*(lambdat - lambda))/n + lambda;    
    range = [lambda lambda+5*(guess - lambda)];
    tolerance = .00001;
    while 1
        lambdae = mean(range);
        result = lambda + n*(lambdae-lambda)/N - (z/N)*((N-n)*lambda + n*lambdae)^.5 - lambdat;
        if (abs(result) < tolerance)
            break
        elseif (result > 0)
            range = [range(1) lambdae];
        else 
            range = [lambdae range(2)];
        end
    end
    
    %test(lambda, n, N, lambdat, p, lambdae)
    
function test(lambda, n, N, lambdat, p, lambdae)
    trials = 1000;
    rates =  [lambda + lambda^.5*randn(N-n, trials); lambdae + lambdae^.5*randn(n, trials)];
    m = mean(rates, 1);
    
    sprintf('obtained: %4.4f  ideal: %4.4f', mean(m <= lambdat), p)
    