clear
close all
load finalr1.mat

ax = figure;

v = VideoWriter('test.avi');
open(v);

loops = iterations;
F(loops) = struct('cdata',[],'colormap',[]);
for j = 1:loops
    
    ax.NextPlot = 'replace';
    matrix(1,1:2) = cov_matrix(j,1:2);
    matrix(2,1:2) = cov_matrix(j,3:4);
    [ExpSigma,ExpCorrC] = cov2corr(matrix);
    scatter(l(:,1),l(:,2),'k','filled')
    hold on
    scatter(x_kpost(j,1), x_kpost(j,2),10,'r','filled')
    hold on
    scatter(x_true(j,1), y_true(j,1),10,'b','filled')
    grid on
    h = plotEllipses([x_kpost(j,1),x_kpost(j,2)],3*ExpSigma);
    xlabel('X [m]')
    ylabel('Y [m]')
    xlim([-2,10])
    ylim([-3,4])
    drawnow
    F(j) = getframe(gcf);
    writeVideo(v,F(j));
end
close(v);