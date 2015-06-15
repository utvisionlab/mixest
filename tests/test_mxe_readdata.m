function test_mxe_readdata

% matrix input
data = mxe_readdata([10 20 30]);
assertEqual(data.data, [10 20 30]);

data = mxe_readdata([10 20 30], false);
assertEqual(data.data, [10 20 30]);
assertEqual(data.weight, []);
assertEqual(data.index, [1 2 3]);

% struct input
data = mxe_readdata(struct('data',[10 20 30], 'weight',[0.2 0.3], 'index',[2 3]));
assertEqual(data.data, [20 30]);
assertEqual(data.weight, [0.2 0.3]);
assertEqual(data.index, []);

data = mxe_readdata(struct('data',[10 20 30], 'weight',[0.2 0.3], 'index',[2 3]), false);
assertEqual(data.data, [10 20 30]);
assertEqual(data.weight, [0.2 0.3]);
assertEqual(data.index, [2 3]);

data = mxe_readdata(struct('data',[10 20 30], 'weight',[0.1 0.2 0.3]));
assertEqual(data.data, [10 20 30]);
assertEqual(data.weight, [0.1 0.2 0.3]);
assertEqual(data.index, []);

