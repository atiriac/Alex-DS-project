function lookAtCells_redraw(cellNumb, fig_info,...
    fig_dfof, fig_wvfs,fig_tuning,fig_onOffRawDf,fig_onOffAvgDf,thetas, ...
    roiIntDS,roiIntOnOff,DSI,vecSum, vecTheta,wvf_resp_reordered, wvf_resp_mean, sumVar, ...
    rhos_all,rhos1_all,rhos2_all,rhos3_all, rhos_norm_all,...
    wvf_resp_t1,wvf_resp_t2,wvf_resp_t3) 

dfof_llim = -0.3; %Sets the lower limit for dfof in plots below
dfof_ulim = 1; %Sets the upper limit for dfof in plots below

% Write in some stats for the cell PRE
set(0, 'currentfigure', fig_info);
set(gca,'units','normalized','position',[0,0,1,1])
cla
text(0.1,0.6,'DSI')
text(0.1,0.4,num2str(DSI(cellNumb,1)))
text(0.4,0.6,'VS')
text(0.4,0.4,num2str(vecSum(cellNumb,1)))
text(0.7,0.6,'var')
text(0.7,0.4,num2str(sumVar(cellNumb,1)))


% Draw dF/F trace non-ordered
set(0, 'currentfigure', fig_dfof);
plot(roiIntDS(cellNumb,:),'k','LineWidth',1)
ylim([dfof_llim dfof_ulim]);
xlim([0,700]);
xlabel('frame')
set(gca,'units','normalized','position',[.035,.025,.96,.95])


[dummy,wvfLength] = size(wvf_resp_mean);
% Draw wvf trace for each direction (mean and individual trials)
set(0, 'currentfigure', fig_wvfs);
set(fig_wvfs,'Name',['Cell # ', num2str(cellNumb)],'NumberTitle','off')
plot(wvf_resp_mean(cellNumb,:),'r', 'LineWidth',2)
hold on
plot(wvf_resp_t1(cellNumb,:),'k')
plot(wvf_resp_t2(cellNumb,:),'k')
plot(wvf_resp_t3(cellNumb,:),'k')
hold off
ylim([dfof_llim dfof_ulim]);
xlim([0,wvfLength]);
xticks(wvfLength/16:wvfLength/8:wvfLength)
xticklabels({'0','45','90','135','180','225','270','315'})
xlabel('angle')
set(gca,'units','normalized','position',[.05,.12,.945,.86])



% Draw tuning curve on polar plot
set(0, 'currentfigure', fig_tuning);
set(fig_tuning,'Name',['Cell # ', num2str(cellNumb)],'NumberTitle','off')
polarplot(thetas, rhos_all(cellNumb,:),'k','LineWidth',2)
hold on;
polarplot(thetas, rhos1_all(cellNumb,:),'Color', [0.2 0.2 0.2],'LineWidth',0.5)
polarplot(thetas, rhos2_all(cellNumb,:),'Color', [0.2 0.2 0.2],'LineWidth',0.5)
polarplot(thetas, rhos3_all(cellNumb,:),'Color', [0.2 0.2 0.2],'LineWidth',0.5)
polarplot([0, vecTheta(cellNumb,1)], [0, vecSum(cellNumb,1)]*dfof_ulim,'Color', [1,0,0], 'LineWidth',2) %non-normalized vector sum
hold off;
rlim([0 dfof_ulim])

[dummy,lengthOnOffTrial] = size(roiIntOnOff);
% Draw dF/F trace of OnOff trial
set(0, 'currentfigure', fig_onOffRawDf);
plot(roiIntOnOff(cellNumb,:),'k','LineWidth',1)
ylim([dfof_llim dfof_ulim]);
xlim([0,lengthOnOffTrial]);
xlabel('frame')
set(gca,'units','normalized','position',[.035,.025,.96,.95])



         


end




















