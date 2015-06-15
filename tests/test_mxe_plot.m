function test_mxe_plot

    h = figure;
    options = mxe_options();
    plot_options = options.plotcost;
    plot_options.axes = gca;
    plot_options.log = false;
    mp = mxe_plot({'sine', 'cosine'}, plot_options);
    x = 0:0.1:10;
    mp.update(x, sin(x))
    mp.update(x, cos(x), 2)
    close(h)
end
