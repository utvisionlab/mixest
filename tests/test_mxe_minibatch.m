function test_mxe_minibatch

    idxAll = 1:5;
    
    options = mxe_options();
    
    % test with minibatchsize = 3;
    options.minibatch.size = 3;
    [mb, idxTrain] = mxe_minibatch(idxAll, options);
    assertEqual(idxTrain, [1 2 3]);
    idxTrain = mb.next();
    assertEqual(idxTrain, [4 5 1]);
    idxTrain = mb.next();
    assertEqual(idxTrain, [2 3 4]);

    % test with minibatchsize > nTrain;
    options.minibatch.size = numel(idxAll)+1;
    [mb, idxTrain] = mxe_minibatch(idxAll, options); %#ok<ASGLU>
    assertEqual(idxTrain, idxAll);
    
    % test with overlap
    options.minibatch.size = 3;
    options.minibatch.overlap = 2;
    [mb, idxTrain] = mxe_minibatch(idxAll, options);
    assertEqual(idxTrain, [1 2 3]);
    idxTrain = mb.next();
    assertEqual(idxTrain, [2 3 4]);
    idxTrain = mb.next();
    assertEqual(idxTrain, [3 4 5]);
    idxTrain = mb.next();
    assertEqual(idxTrain, [4 5 1]);
    idxTrain = mb.next();
    assertEqual(idxTrain, [5 1 2]);
    
end