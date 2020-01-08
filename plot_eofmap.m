function plot_eofmap(n_mode,lon,lat,eof_map,expvar)


o_contourf(lon,lat,eof_map(:,:,n_mode)');
% plot(time,pcs(n_mode,:),'linewidth',2);
set(gca,'fontsize',17);
title(['EOF mode ',num2str(n_mode),', ',num2str(expvar(n_mode),'%.2f'),'%']);
c = colorbar('fontsize',17);
end
