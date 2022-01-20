function plotFonts(h)

set(gca,'FontSize',24,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',26,'fontWeight','bold')
h(1).MarkerSize = 8;
h(1).LineWidth = 1;

end

