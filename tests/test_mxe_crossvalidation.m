function test_mxe_crossvalidation

    data = [1 1];
    cvfraction = 0.5;
    numkeepsave = 5;
    maxiter = 10;
    cv_rise_iter = 2;

    cv = mxe_crossvalidation(data, ...
        cvfraction, numkeepsave, @costfun);

    for iter = 0:maxiter
        theta = iter;
        
        last = iter+1;
        stats.iter = iter;
        stats = cv.statsfun(theta, stats);
        info(last) = stats; %#ok<AGROW>
        stopnow = cv.stopfun(theta, info, last);
        if stopnow
            assert(iter >= cv_rise_iter)
            break
        end
    end
    
    assert(cv.hasFinalResult())
    [theta, cost] = cv.finalResult(); %#ok<NASGU>
    iter = theta;
    assert(iter <= cv_rise_iter)
    
    

    function cost = costfun(theta, data) %#ok<INUSD>
    % theta is the actually the iteration number. The returned cost is
    % descending until cv_rise_iter and ascending afterwards.
        cost = abs(theta - cv_rise_iter);
    end
end
