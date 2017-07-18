function [se, figs, rmse_wVL, rmse_woVL, nFib_tract, nFib_PN, nVoxels] =  feVirtualLesion(fe, fascicleIndices, display)
% 
% Perform a virtual lesion and generate plots of the results.
%
% - Remove a set of fascicles from a whole-brain connectome
% - Find the fascicle's path-neighborhood (the set of fascicles that share the same voxels)
% - Reduce the connectome to only the voxels of the set of fascicels that
%   we are lesioning. These are the voxels that contain the signal and the
%   error useful for the lesion test.
% - Compute a series of statistics on the lesioned connectome.
%
% [se, figure_handls, rmse_wVL, rmse_woVL, nFib_tract, nFib_PN, nVoxels] = ...
%           feVirtualLesion(feNoLesion, fascicleIndices, refitConnectome, display)
%
% Copyright (2016), Franco Pestilli, Indiana University, frakkopesto@gmail.com

 if notDefined('display'),   
    display.tract = 0;
    display.distributions = 0; 
    display.evidence = 0; 
    fig = []; 
end

% First we compute all the variables needed to compute the virtual lesion.
[ rmse_wVL, rmse_woVL, nFib_tract, nFib_PN, nVoxels] = feComputeVirtualLesion(fe, fascicleIndices);

% After that we compute the actual virtual lesion.
se = feComputeEvidence(rmse_wVL,rmse_woVl);

% Make plots of the results.
fig(1).name = sprintf('rmse_distributions_%s',mfilename);
fig(1).h   = nan;
fig(1).type = 'eps';
if display.distributions
    % Raw RMSE distirbutions
    fig(1).name = sprintf('rmse_distributions_%s',mfilename);
    fig(1).h   = figure('name',fig(1).name,'color','w');
    fig(1).type = 'eps';
    set(fig(1).h,'Units','normalized','Position',[0.007 0.55  0.28 0.36]);
    plot(se.lesion.xhist,se.lesion.hist,'-','color', [.95 .45 .1],'linewidth',2); hold on
    plot(se.nolesion.xhist,se.nolesion.hist,'-','linewidth',2, 'color', [.1 .45 .95])
    plot([se.nolesion.rmse.mean,se.nolesion.rmse.mean], [0,0.12],'-','color',[.1 .45 .95] ) 
    plot([se.lesion.rmse.mean,se.lesion.rmse.mean], [0,0.12], '-', 'color',[.95 .45 .1])
    title(sprintf('mean RMSE\nno-lesion %2.3f | lesion %2.2f',se.nolesion.rmse.mean,se.lesion.rmse.mean),'fontsize',16)
    ylabel('Probability', 'fontsize',14);xlabel('RMSE', 'fontsize',14)
    legend({'Lesion','No lesion'},'fontsize',14);
    set(gca,'box','off','xtick',[0 round(se.xrange(2)/2) se.xrange(2)],'ytick',[0 .06 .12],'xlim',[0 se.xrange(2)],'ylim',[0 .125], ...
        'tickdir', 'out', 'ticklength', [0.025 0])

    % Plot the null distribution and the empirical difference
    % We reorganize the variables names hee below just to keep the plotting
    % code more compact.
    ywo_e = se.s.lesioned_e;
    y_e   = se.s.unlesioned_e;
    woxhis = se.s.lesioned.xbins;
    xhis   = se.s.unlesioned.xbins;
    min_x = se.s.min_x;
    max_x = se.s.max_x; 
    fig(2).name = sprintf('virtual_lesion_test_mean_rmse_hist_%s_%s',mfilename,feLesion.name);
    fig(2).h   = figure('name',fig(2).name,'color','w');  
    fig(2).type = 'eps';
    set(fig(2).h,'Units','normalized','Position',[0.007 0.55  0.28 0.36]);
    patch([xhis,xhis],y_e(:),[.1 .45 .95],'FaceColor',[.1 .45 .95],'EdgeColor',[.1 .45 .95]);
    hold on
    patch([woxhis,woxhis],ywo_e(:),[.95 .45 .1],'FaceColor',[.95 .45 .1],'EdgeColor',[.95 .45 .1]);
    set(gca,'tickdir','out', ...
        'box','off', ...
        'ylim',[0 0.25], ...
        'xlim',[min_x,max_x], ...
        'ytick',[0 0.1 0.2], ...
        'xtick',round(linspace(min_x,max_x,4)), ...
        'fontsize',16)
    ylabel('Probability','fontsize',16)
    xlabel('rmse','fontsize',16')
    title(sprintf('Strength of connection evidence %2.3f',(se.s.mean)), 'FontSize',16)
end

fig(3).name = sprintf('tract_%s',mfilename);
fig(3).h   = nan;
fig(3).type = 'jpg';
fig(4).name = sprintf('pathneighborhood_%s',mfilename);
fig(4).h   = nan;
fig(4).type = 'jpg';
fig(5).name = sprintf('pathneighborhood_AND_tract_%s',mfilename);
fig(5).h   = nan;
fig(5).type = 'jpg';
if display.tract
    % Now display the lesion performed.
    % (1) The fascicle to delete.
    fig(3).name = sprintf('tract_%s',mfilename);
    fig(3).h   = figure('name',fig(3).name,'color',[.1 .45 .95]); 
    fig(3).type = 'jpg';
    [fig(3).h, fig(3).light] =  mbaDisplayConnectome(mbaFiberSplitLoops(fasAcpc.fibers),fig(3).h, [.1 .45 .95],'single',[],[],.2);
    view(-23,-23);delete(fig(3).light); fig(3).light = camlight('right');
        
    tmp_fas = feGet(feNoLesion,'fibers acpc');
    fibers2display = randsample(1:length(tmp_fas.fibers),ceil(0.01*length(tmp_fas.fibers)));
    fig(4).name = sprintf('pathneighborhood_%s',mfilename);
    fig(4).h   = figure('name',fig(4).name,'color',[.1 .45 .95]);     
    fig(4).type = 'jpg';
    [fig(4).h, fig(4).light] =  mbaDisplayConnectome(mbaFiberSplitLoops(...
        tmp_fas.fibers(fibers2display)),fig(4).h, [.95 .45 .1],'single',[],.6,.2);
    view(-23,-23);delete(fig(2).light); fig(4).light = camlight('right');

    % (2) The facicle with the pah neighborhood  
    fig(5).name = sprintf('pathneighborhood_AND_tract_%s',mfilename);
    fig(5).h   = figure('name',fig(5).name,'color',[.1 .45 .95]);
    fig(5).type = 'jpg';
    [fig(5).h, fig(5).light] =  mbaDisplayConnectome(mbaFiberSplitLoops(fasAcpc.fibers),fig(5).h, [.1 .45 .95],'single',[],[],.2);
    view(-23,-23);delete(fig(5).light);
    hold on
    [fig(5).h, fig(5).light] =  mbaDisplayConnectome(mbaFiberSplitLoops(...
        tmp_fas.fibers(fibers2display)),fig(5).h, [.95 .45 .1],'single',[],.6,.2);
    view(-23,-23);delete(fig(5).light); fig(5).light = camlight('right');
 
    % Show the voxels coordinates to see if we are in the right spot
    plot_test_fiber_location_debug= false;
    if plot_test_fiber_location_debug
        figure('name','Coordinate check');
        plot3(allCoords(commonCoords,1),allCoords(commonCoords,2), ...
            allCoords(commonCoords,3),'ko','MarkerFaceColor','k','MarkerSize',8);
        hold on;
        plot3(fasCoords(:,1),fasCoords(:,2),fasCoords(:,3),'ro', ...
            'MarkerFaceColor',[.95 .45 .1],'MarkerSize',3);
        axis equal
        view(-23,-23); hold on
        plot3(allCoords(commonCoords,1),allCoords(commonCoords,2), ...
            allCoords(commonCoords,3),'go','MarkerFaceColor','g','MarkerSize',12);
        view(-23,-23);
    end
end

% RMSE distributions
fig(6).name = sprintf('Size_of_effect_of_the_lesion_%s',mfilename);
fig(6).type = 'eps';   
fig(6).h   = nan;
if display.evidence
    fig(6).name = sprintf('Size_of_effect_of_the_lesion_%s',mfilename);
    fig(6).h   = figure('name',fig(6).name,'color','w');
    set(fig(6).h,'Units','normalized','Position',[0.007 0.55  0.28 0.36]);
    subplot(1,4,1)
    plot(1,se.s.mean,'-o','color', [.95 .45 .1],'linewidth',2); hold on
    plot([1,1], [se.s.mean,se.s.mean] + [-se.s.std,se.s.std], '-','color',[.95 .45 .1] )
    ylabel('S (s.d.)', 'fontsize',14);
    set(gca,'box','off','xlim',[0 2], 'ylim',[0 ceil(se.s.mean + se.s.std)], ...
        'tickdir', 'out', 'ticklength', [0.025 0])
    subplot(1,4,2)
    plot(1,se.em.mean,'-o','color', [.95 .45 .1],'linewidth',2); hold on
    ylabel('Earth mover''s distance (raw scanner units)', 'fontsize',14);
    set(gca,'box','off','xlim',[0 2], 'ylim',[0 ceil(se.em.mean)], ...
        'tickdir', 'out', 'ticklength', [0.025 0])
    subplot(1,4,3)
    plot(1,se.kl.mean,'-o','color', [.95 .45 .1],'linewidth',2); hold on
    ylabel('K-L divergence (bits)', 'fontsize',14);
    set(gca,'box','off','xlim',[0 2], 'ylim',[0 ceil(se.kl.mean)], ...
        'tickdir', 'out', 'ticklength', [0.025 0])
    subplot(1,4,4)
    plot(1,se.j.mean,'-o','color', [.95 .45 .1],'linewidth',2); hold on
    ylabel('Jeffrey''s divergence (bits)', 'fontsize',14);
    set(gca,'box','off','xlim',[0 2], 'ylim',[0 ceil(se.j.mean)], ...
        'tickdir', 'out', 'ticklength', [0.025 0])
end

end % Main function
