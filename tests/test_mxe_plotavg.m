function test_mxe_plotavg

    plotx = 1:10;
    ploty = 11:20;

    [plotx, ploty] = mxe_plotavg(plotx, ploty, 4);
    assertEqual(plotx, [6 10]);
    assertEqual(ploty, [14.5 18.5]);
end