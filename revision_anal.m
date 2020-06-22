% New analyses for revision!
% run this only after running scratch and without changing the variable
% space

% anovas looking for order and learning effects
b.xBlock = [b.Geb1,b.Geb2,b.Pan1,b.Pan2];
protocolxblock = [repmat({'Gebuis'},1,108),repmat({'Panamath'},1,108)];
time = [ones(1,54)+0,ones(1,54)+1,ones(1,54)+0,ones(1,54)+1];
partxBlock = [1:54,1:54,1:54,1:54];
[p table stats] = ...
    anovan(b.xBlock(2,:),...
    {time,repmat(firstBlockPart,1,4),protocolxblock,partxBlock},...
    'random',4,'nested',[0,0,0,0;0,0,0,0;0,0,0,0;0,1,0,0],...
    'varnames',{'Time','CounterBalOrder','Protocol','Participant'})
[p table stats] = ...
    anovan(b.xBlock(3,:),...
    {time,repmat(firstBlockPart,1,4),protocolxblock,partxBlock},...
    'random',4,'nested',[0,0,0,0;0,0,0,0;0,0,0,0;0,1,0,0],...
    'varnames',{'Time','CounterBalOrder','Protocol','Participant'})
[p table stats] = ...
    anovan(b.xBlock(4,:),...
    {time,repmat(firstBlockPart,1,4),protocolxblock,partxBlock},...
    'random',4,'nested',[0,0,0,0;0,0,0,0;0,0,0,0;0,1,0,0],...
    'varnames',{'Time','CounterBalOrder','Protocol','Participant'})

gebFirstBetas(:,1) = mean(b.Geb1(:,strcmp(firstBlockPart,'Gebuis')),2);
gebFirstBetas(:,2) = mean(b.Pan1(:,strcmp(firstBlockPart,'Gebuis')),2);
gebFirstBetas(:,3) = mean(b.Geb2(:,strcmp(firstBlockPart,'Gebuis')),2);
gebFirstBetas(:,4) = mean(b.Pan2(:,strcmp(firstBlockPart,'Gebuis')),2);

panFirstBetas(:,1) = mean(b.Pan1(:,strcmp(firstBlockPart,'Panamath')),2);
panFirstBetas(:,2) = mean(b.Geb1(:,strcmp(firstBlockPart,'Panamath')),2);
panFirstBetas(:,3) = mean(b.Pan2(:,strcmp(firstBlockPart,'Panamath')),2);
panFirstBetas(:,4) = mean(b.Geb2(:,strcmp(firstBlockPart,'Panamath')),2);

gebFirstSEMs(:,1) = std(b.Geb1(:,strcmp(firstBlockPart,'Gebuis')),2)/sqrt(sum(strcmp(firstBlockPart,'Gebuis')));
gebFirstSEMs(:,2) = std(b.Pan1(:,strcmp(firstBlockPart,'Gebuis')),2)/sqrt(sum(strcmp(firstBlockPart,'Gebuis')));
gebFirstSEMs(:,3) = std(b.Geb2(:,strcmp(firstBlockPart,'Gebuis')),2)/sqrt(sum(strcmp(firstBlockPart,'Gebuis')));
gebFirstSEMs(:,4) = std(b.Pan2(:,strcmp(firstBlockPart,'Gebuis')),2)/sqrt(sum(strcmp(firstBlockPart,'Gebuis')));

panFirstSEMs(:,1) = std(b.Pan1(:,strcmp(firstBlockPart,'Panamath')),2)/sqrt(sum(strcmp(firstBlockPart,'Panamath')));
panFirstSEMs(:,2) = std(b.Geb1(:,strcmp(firstBlockPart,'Panamath')),2)/sqrt(sum(strcmp(firstBlockPart,'Panamath')));
panFirstSEMs(:,3) = std(b.Pan2(:,strcmp(firstBlockPart,'Panamath')),2)/sqrt(sum(strcmp(firstBlockPart,'Panamath')));
panFirstSEMs(:,4) = std(b.Geb2(:,strcmp(firstBlockPart,'Panamath')),2)/sqrt(sum(strcmp(firstBlockPart,'Panamath')));.



% Figure 4A (beta space)
figureDim = [0 0 9 12]; % left bottom width height?
FontSize = 9;
h(1) = figure('PaperUnits','centimeters','PaperPosition',figureDim,...
    'Units','centimeters','Position',figureDim);
plot([0 0],[-4 8],'color',[.5 .5 .5]);hold on; % number
plot([-4 4],[0 0],'color',[.5 .5 .5]); % size
plot([-4 8]/3,[-4 8],'color',[.5 .5 .5]); % perimeter
plot([-4 8],[-4 8],'color',[.5 .5 .5]); % total surface area
plot([-4 8],[4 -8],'color',[.5 .5 .5]); % item area
% plot means
plot(mean(b.Pan(3,:)),mean(b.Pan(2,:)),'marker','o','markersize',12,...
    'color',[0.75 0 0],'linestyle','none')
plot(mean(b.Geb(3,:)),mean(b.Geb(2,:)),'marker','x','markersize',12,...
    'color',[0 0 0.75],'linestyle','none')
plot(gebFirstBetas(3,:),gebFirstBetas(2,:),'--k')
plot(panFirstBetas(3,:),panFirstBetas(2,:),'--r')
% plot sems
sem.Pan(3) = std(b.Pan(3,:))/sqrt(size(b.Pan,2));
sem.Geb(3) = std(b.Geb(3,:))/sqrt(size(b.Geb,2));
sem.Pan(2) = std(b.Pan(2,:))/sqrt(size(b.Pan,2));
sem.Geb(2) = std(b.Geb(2,:))/sqrt(size(b.Geb,2));
plot([mean(b.Pan(3,:))-sem.Pan(3),mean(b.Pan(3,:))+sem.Pan(3)],...
    [mean(b.Pan(2,:)),mean(b.Pan(2,:))],'-r')
plot([mean(b.Pan(3,:)),mean(b.Pan(3,:))],...
    [mean(b.Pan(2,:))-sem.Pan(2),mean(b.Pan(2,:))+sem.Pan(2)],'-r')
plot([mean(b.Geb(3,:))-sem.Geb(3),mean(b.Geb(3,:))+sem.Geb(3)],...
    [mean(b.Geb(2,:)),mean(b.Geb(2,:))],'-b')
plot([mean(b.Geb(3,:)),mean(b.Geb(3,:))],...
    [mean(b.Geb(2,:))-sem.Geb(2),mean(b.Geb(2,:))+sem.Geb(2)],'-b')
% format figure
set(gca,'DataAspectRatio',[1 1 1]);
xlabel('bSize','FontName','Arial','FontWeight','normal','FontSize',FontSize);
ylabel('bNum','FontName','Arial','FontWeight','normal','FontSize',FontSize);
set(gca,'xtick',-4:4,'ytick',-4:8);
axis([-3 3 -0.25 8])
set(gca,'Units','centimeters','OuterPosition',figureDim,...
    'color','none','Box','off','FontName','Arial',...
    'FontWeight','Normal','FontSize',FontSize,'LineWidth',1,'Clipping','on')

% Figure 4B: Num x Spacing
figureDim = [0 0 9 12]; % left bottom width height?
FontSize = 8;
h(2) = figure('PaperUnits','centimeters','PaperPosition',figureDim,...
    'Units','centimeters','Position',figureDim);
plot([0 0],[-4 8],'color',[.5 .5 .5]);hold on; % number
plot([-4 4],[0 0],'color',[.5 .5 .5]); % spacing
plot([-4 8],[-4 8],'color',[.5 .5 .5]); % Field area
plot([-4 8],[4 -8],'color',[.5 .5 .5]); % Sparsity
% plot means
plot(mean(b.Pan(4,:)),mean(b.Pan(2,:)),'marker','o','markersize',12,...
    'color',[0.75 0 0],'linestyle','none')
plot(mean(b.Geb(4,:)),mean(b.Geb(2,:)),'marker','x','markersize',12,...
    'color',[0 0 0.75],'linestyle','none')
plot(gebFirstBetas(4,:),gebFirstBetas(2,:),'--k')
plot(panFirstBetas(4,:),panFirstBetas(2,:),'--r')
% plot sems
sem.Pan(4) = std(b.Pan(4,:))/sqrt(size(b.Pan,2));
sem.Geb(4) = std(b.Geb(4,:))/sqrt(size(b.Geb,2));
sem.Pan(2) = std(b.Pan(2,:))/sqrt(size(b.Pan,2));
sem.Geb(2) = std(b.Geb(2,:))/sqrt(size(b.Geb,2));
plot([mean(b.Pan(4,:))-sem.Pan(4),mean(b.Pan(4,:))+sem.Pan(4)],...
    [mean(b.Pan(2,:)),mean(b.Pan(2,:))],'-r')
plot([mean(b.Pan(4,:)),mean(b.Pan(4,:))],...
    [mean(b.Pan(2,:))-sem.Pan(2),mean(b.Pan(2,:))+sem.Pan(2)],'-r')
plot([mean(b.Geb(4,:))-sem.Geb(4),mean(b.Geb(4,:))+sem.Geb(4)],...
    [mean(b.Geb(2,:)),mean(b.Geb(2,:))],'-b')
plot([mean(b.Geb(4,:)),mean(b.Geb(4,:))],...
    [mean(b.Geb(2,:))-sem.Geb(2),mean(b.Geb(2,:))+sem.Geb(2)],'-b')
% format figure
set(gca,'DataAspectRatio',[1 1 1]);
xlabel('bSpacing','FontName','Arial','FontWeight','normal','FontSize',FontSize);
ylabel('bNum','FontName','Arial','FontWeight','normal','FontSize',FontSize);
set(gca,'xtick',-4:4,'ytick',-4:8);
axis([-3 3 -0.25 8])
set(gca,'Units','centimeters','OuterPosition',figureDim,...
    'color','none','Box','off','FontName','Arial',...
    'FontWeight','Normal','FontSize',FontSize,'LineWidth',1,'Clipping','on')








