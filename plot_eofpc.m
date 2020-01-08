function plot_eofpc(n_mode,time,pc,expvar)


plot(double(time),pc(n_mode,:),'linewidth',2);
datetick('x','yyyy');
xlim([time(1),time(end)]);
title(['PC-',num2str(n_mode),', ',num2str(expvar(n_mode),'%.2f'),'%']);
set(gca,'fontsize',17);
end