vd = [7:15];
len = length(vd);
tr = 30;
num = 5;
time = zeros(len,1);
for i = 1:len
    M = randn(tr,2^vd(i))+1i*vd(i)*rand(tr,2^vd(i));
    tic
    for j = 1:num
        [~,R,E]=qr(M,0);
    end
    time(i) = toc/num;
end
    pic1 = figure;
    hold on;
    h(1) = plot(vd,vd,'-r','LineWidth',4);
    h(2) = plot(vd,2*vd,'-b','LineWidth',4);
    h(3) = plot(vd,log2(time),'-^k','LineWidth',2);
    legend('N','N^2','QR','Location','NorthWest');
    axis square;
    xlabel('log_2(N)'); ylabel('log_{2}(time)');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    saveas(pic1,['QRscaling.eps'],'epsc');
    