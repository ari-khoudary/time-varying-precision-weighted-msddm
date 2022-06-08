function makePretty(fontSize)

% makes plots pretty! written by Megan Peters

set(gcf,'color','white')

set(findall(gcf,'-property','FontSize'),'FontSize',fontSize)
set(findall(gcf,'-property','FontName'),'FontName','Arial')
% set(findall(gcf,'-property','MarkerSize'),'MarkerSize',10)
% set(findall(gcf,'-property','LineWidth'),'LineWidth',2)