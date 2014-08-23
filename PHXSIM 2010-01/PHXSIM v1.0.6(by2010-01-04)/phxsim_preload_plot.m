%%%2009-12-27 phxsim_preload_plot.m: should be called by phxsim.m

h=figure;

switch foldingModel
    case 'iup1'        
        if flagSF==1
            semilogx(SFfitN(:,1), SFfitN(:,2),'k','LineWidth',2); hold on %plot stop-flow fluor curve
            [t,y] = ode15s(@fm_iup1,[SFfitN(1,1) SFfitN(end,1)],[1 0 0]);
        else %flagSF is 0
            [t,y] = ode15s(@fm_iup1,[0 10],[1 0 0]);
        end
        
        semilogx(t,y(:,1),'r'); hold on %plot U
        semilogx(t,y(:,2),'g'); hold on %plot I
        semilogx(t,y(:,3),'b'); hold on %plot N
                
        %%%label U/N/I's fraction at the end of folding phase:
        [t,y] = ode15s(@fm_iup1,0:Time_f:Time_f,[1 0 0]);
        stem(t(end,1),y(end,1),'r','LineStyle','none'); hold on %plot U
        text(t(end,1),y(end,1), ['U=', num2str(y(end,1)*100, '%5.1f'),'%'],'FontWeight','bold');
        stem(t(end,1),y(end,2),'g','LineStyle','none'); hold on %plot I
        text(t(end,1),y(end,2), ['I=', num2str(y(end,2)*100, '%5.1f'),'%'],'FontWeight','bold');
        stem(t(end,1),y(end,4),'b','LineStyle','none'); hold on %plot N
        text(t(end,1),y(end,4), ['N=', num2str(y(end,4)*100, '%5.1f'),'%'],'FontWeight','bold');
        
        v=axis;
        axis([1e-4,v(2),-0.01,1.01]);
        xlabel('Time (sec)')
        ylabel('Fraction')
        if flagSF==1
            title('IUP1 Model Fitting (red=U; green=I; blue=N; black=exp fluor)')
        else %flagSF is 0
            title('IUP1 Model Simulation (red=U; green=I; blue=N)')
        end
        
    case 'ppoe1'
        if flagSF==1
            semilogx(SFfitN(:,1), SFfitN(:,2),'k','LineWidth',2); hold on %plot stop-flow fluor curve
            [t,y] = ode15s(@fm_ppoe1,[SFfitN(1,1) SFfitN(end,1)],[1 0 0 0]);
        else %flagSF is 0
            [t,y] = ode15s(@fm_ppoe1,[0 10],[1 0 0 0]);
        end
        
        semilogx(t,y(:,1),'r'); hold on %plot U
        semilogx(t,y(:,2),'g'); hold on %plot I
        semilogx(t,y(:,3),'y'); hold on %plot Ix
        semilogx(t,y(:,4),'b'); hold on %plot N
                
        %%%label U/N/I's fraction at the end of folding phase:
        [t,y] = ode15s(@fm_ppoe1,0:Time_f:Time_f,[1 0 0 0]);
        stem(t(end,1),y(end,1),'r','LineStyle','none'); hold on %plot U
        text(t(end,1),y(end,1), ['U=', num2str(y(end,1)*100, '%5.1f'),'%'],'FontWeight','bold');
        stem(t(end,1),y(end,2),'g','LineStyle','none'); hold on %plot I
        text(t(end,1),y(end,2), ['I=', num2str(y(end,2)*100, '%5.1f'),'%'],'FontWeight','bold');
        stem(t(end,1),y(end,3),'y','LineStyle','none'); hold on %plot Ix
        text(t(end,1),y(end,3), ['Ix=', num2str(y(end,3)*100, '%5.1f'),'%'],'FontWeight','bold');
        stem(t(end,1),y(end,4),'b','LineStyle','none'); hold on %plot N
        text(t(end,1),y(end,4), ['N=', num2str(y(end,4)*100, '%5.1f'),'%'],'FontWeight','bold');
        
        v=axis;
        axis([1e-4,v(2),-0.01,1.01]);
        xlabel('Time (sec)')
        ylabel('Fraction')
        if flagSF==1
            title('PPOE1 Model Fitting (red=U; green=I; yellow=Ix; blue=N; black=exp fluor)')
        else %flagSF is 0
            title('PPOE1 Model Simulation (red=U; green=I; yellow=Ix; blue=N)')
        end
        
end

%%%save above figure:
SaveFigureName=[proteinName '_phxsim_FoldingModelPlot.fig'];
saveas(figure(h),SaveFigureName)


