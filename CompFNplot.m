    flag = 1;
    vd = S(:,1);
rank1 = S(:,2);
rank2 = S(:,3);
rank3 = S(:,4);
timeour1=S(:,5);
timeour2=S(:,6);
timeour3=S(:,7);
errorour1=S(:,8);
errorour2=S(:,9);
errorour3=S(:,10);
timefac1=S(:,12);
timefac2=S(:,13);
timefac3=S(:,14);
    pic = figure;
    hold on;
    ag = (log2(timeour1(1))+log2(timeour2(1))+log2(timeour3(1)))/3;
    h(1) = plot(vd,vd+log2(vd)-vd(1)-log2(vd(1))+ag,'--c','LineWidth',4);
    h(2) = plot(vd,vd+2*log2(vd)-vd(1)-2*log2(vd(1))+ag,'--k','LineWidth',4);
    h(3) = plot(vd,log2(timeour1),'-^r','LineWidth',2);
    h(4) = plot(vd,log2(timeour2),'-^b','LineWidth',2);
    h(5) = plot(vd,log2(timeour3),'-^g','LineWidth',2);
    legend('N log N','N log^2 N','NP2 app','NP1 app','NP0 app','Location','bestoutside');
    %if flag > 0
    %   title('RS FFT vs CHEB NUFFT time, uni JPT');
    %else
    %    title('RS FFT vs CHEB NUFFT time, non JPT');
    %end
    axis tight;
    xlabel('log_2(N)'); ylabel('log_{2}(time)');
    set(gca, 'FontSize', 20);
    b=get(gca);
    set(b.XLabel, 'FontSize', 20);set(b.YLabel, 'FontSize', 20);set(b.ZLabel, 'FontSize', 20);set(b.Title, 'FontSize', 20);
    if flag > 0
       saveas(pic,['Comp_RSFFT_CHEBNF_uni1.eps'],'epsc');
    else 
       saveas(pic,['Comp_RSFFT_CHEBNF_non1.eps'],'epsc');
    end
    hold off;
    pic = figure;
    hold on;
    ag = (log2(timefac1(1))+log2(timefac1(1))+log2(timefac3(1)))/3;
    h(1) = plot(vd,vd+log2(vd)-vd(1)-log2(vd(1))+ag,'--c','LineWidth',4);
    h(2) = plot(vd,vd+2*log2(vd)-vd(1)-2*log2(vd(1))+ag,'--k','LineWidth',4);
    h(3) = plot(vd,log2(timefac1),'-xr','LineWidth',2);
    h(4) = plot(vd,log2(timefac2),'-xb','LineWidth',2);
    h(5) = plot(vd,log2(timefac3),'-xg','LineWidth',2);
    legend('N log N','N log^2 N','NP2 fac','NP1 fac','NP0 fac','Location','bestoutside');
    %if flag > 0
    %   title('RS FFT vs CHEB NUFFT time, uni JPT');
    %else
    %    title('RS FFT vs CHEB NUFFT time, non JPT');
    %end
    axis tight;
    xlabel('log_2(N)'); ylabel('log_{2}(time)');
    set(gca, 'FontSize', 20);
    b=get(gca);
    set(b.XLabel, 'FontSize', 20);set(b.YLabel, 'FontSize', 20);set(b.ZLabel, 'FontSize', 20);set(b.Title, 'FontSize', 20);
    if flag > 0
       saveas(pic,['Comp_RSFFT_CHEBNF_uni2.eps'],'epsc');
    else 
       saveas(pic,['Comp_RSFFT_CHEBNF_non2.eps'],'epsc');
    end
    hold off;
    pic1 = figure;
    hold on;
    h(1) = plot(vd,log10(errorour1),'-^r','LineWidth',2);
    h(2) = plot(vd,log10(errorour2),'-^b','LineWidth',2);
    h(3) = plot(vd,log10(errorour3),'-^g','LineWidth',2);
    legend('NP2 relerr','NP1 relerr','NP0 relerr','Location','bestoutside');
    %if flag > 0
    %   title('relerr, uni JPT');
    %else
    %    title('relerr, non JPT');
    %end
    axis tight;
    xlabel('log_2(N)'); ylabel('log_{10}(relerr)');
    set(gca, 'FontSize', 20);
    b=get(gca);
    set(b.XLabel, 'FontSize', 20);set(b.YLabel, 'FontSize', 20);set(b.ZLabel, 'FontSize', 20);set(b.Title, 'FontSize', 20);
    if flag > 0
       saveas(pic1,['CompFN_relerr_uni.eps'],'epsc');
    else 
       saveas(pic1,['CompFN_relerr_non.eps'],'epsc');
    end
    hold off;
    pic2 = figure;
    hold on;
    h(1) = plot(vd,rank1,'-^r','LineWidth',2);
    h(2) = plot(vd,rank2,'-^b','LineWidth',2);
    h(3) = plot(vd,rank3,'-^g','LineWidth',2);
    legend('NP2 rank','NP1 rank','NP0 rank','Location','bestoutside');
    %if flag > 0
    %   title('rank, uni JPT');
    %else
    %    title('rank, non JPT');
    %end
    axis tight;
    xlabel('log_2(N)'); ylabel('rank');
    set(gca, 'FontSize', 20);
    b=get(gca);
    set(b.XLabel, 'FontSize', 20);set(b.YLabel, 'FontSize', 20);set(b.ZLabel, 'FontSize', 20);set(b.Title, 'FontSize', 20);
    if flag > 0
       saveas(pic2,['CompFN_rank_uni.eps'],'epsc');
    else 
       saveas(pic2,['CompFN_rank_non.eps'],'epsc');
    end