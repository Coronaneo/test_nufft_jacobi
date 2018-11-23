    %figure('visible','off');
    pic = figure;
    hold on;
    h(1) = plot(vd,log10(timeour1./rank1),'-^r','LineWidth',2);
    h(2) = plot(vd,log10(timeour2./rank2),'-^b','LineWidth',2);

    legend('time-NUFFT','time-FFT','Location','NorthWest');
    if flag > 0
       title('FFT vs NUFFT time, uni_JPT');
    else
        title('FFT vs NUFFT time, non_JPT');
    end
    axis square;
    xlabel('log_2(N)'); ylabel('log_{10}(time)');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);