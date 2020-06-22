% Make scatter plot with regression line
% [rsquare fstat beta] = regScatter(x,y)
% NKD 01-10-2011
function [rsquare fstat beta] = regScatter(x,y)
% [r p]=corrcoef(x,y);
whichstats={'yhat','beta','rsquare','fstat'};
reg=regstats(y,x,'linear',whichstats);
range.x=max(x)-min(x);
range.y=max(y)-min(y);
plotRect=[min(x)-(range.x/10) max(x)+(range.x/10) min(y)-(range.y/10) max(y)+(range.y/10)];
scatter(x,y);hold on;
regressLineY = [reg.beta(1)+reg.beta(2)*min(x) reg.beta(1)+reg.beta(2)*max(x)];
regressLineX = [min(x) max(x)];
plot(regressLineX,regressLineY,'-r'); hold off;
title(sprintf('rsquare = %0.3f p = %0.3f',reg.rsquare,reg.fstat.pval))
axis(plotRect);
rsquare=reg.rsquare;
fstat=reg.fstat;
beta=reg.beta;