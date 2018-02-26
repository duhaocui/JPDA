 txtsz=25;
set(gca,'FontSize',txtsz)
screen_size = get(0, 'ScreenSize');
h_plotszz = get(gca, 'title');
kpp_plotszz = get(gca, 'xlabel');
l_plotszz = get(gca, 'ylabel');
m_plotszz = get(gca, 'zlabel');
set(h_plotszz, 'FontName', 'Helvetica', 'FontSize', txtsz*1.5)
set(kpp_plotszz, 'FontName', 'Helvetica', 'FontSize', txtsz)
set(l_plotszz, 'FontName', 'Helvetica', 'FontSize', txtsz)
set(m_plotszz, 'FontName', 'Helvetica', 'FontSize', txtsz)
% set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperType', 'A4');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
% set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);
 set(gcf, 'PaperPosition', [1 1 27 19]);
 set(gcf, 'Position',  0.9*[1 1 screen_size(3) screen_size(4) ] );