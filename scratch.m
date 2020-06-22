load('ClaytonDataImport')
switches.removeOutliers = 1;

% create size and spacing variables
SizeLeft = pixels1_Left .* DotSize1_Left;
SizeRight = pixels2_Right .* DotSize2_Right;
SparLeft = (ConvexHull1_Left ./ Numerosity1_Left);
SparRight = (ConvexHull2_Right ./ Numerosity2_Right);
SpaceLeft = SparLeft .* ConvexHull1_Left;
SpaceRight = SparRight .* ConvexHull2_Right;

% create right to left ratios and side chosen variables.
bv.Dnum = -log2(LeftRightNumerosityRatio)';
bv.DONSZ = -log2(SizeLeft ./ SizeRight)';
bv.DONSP = -log2(SpaceLeft ./ SpaceRight)';
bv.choice = NaN(numel(Accuracy),1)';
indx = Accuracy & (Numerosity2_Right > Numerosity1_Left);
bv.choice(indx) = 1;
indx = ~Accuracy & ~(Numerosity2_Right > Numerosity1_Left);
bv.choice(indx) = 1;
indx = Accuracy & ~(Numerosity2_Right > Numerosity1_Left);
bv.choice(indx) = 0;
indx = ~Accuracy & (Numerosity2_Right > Numerosity1_Left);
bv.choice(indx) = 0;

% create first 60 index
bv.first60indx = zeros(numel(bv.choice),1);
uniPar = unique(Participant);
for n = 1:numel(uniPar)
    indx = Participant == uniPar(n) & strcmp(Time,'Time1');
    numIndx = find(indx);
    bv.first60indx(numIndx(1:60)) = 1;
    indx = Participant == uniPar(n) & strcmp(Time,'Time2');
    numIndx = find(indx);
    bv.first60indx(numIndx(1:60)) = 1;
end

% create the blockNum and firstBlock variables
GebFirstBlocks = ...
    [zeros(96,1)+1;zeros(60,1)+2;zeros(96,1)+3;zeros(60,1)+4];
PanFirstBlocks = ...
    [zeros(60,1)+1;zeros(96,1)+2;zeros(60,1)+3;zeros(96,1)+4];
bv.BlockNum = NaN(numel(Trialorder),1);
bv.firstBlock = cell(numel(Trialorder),1);
for k = 1:numel(Trialorder)
    indx = Participant==Participant(k) & Trialorder == 1;
    if sum(indx)~=1;keyboard;end;
    if strcmp(Protocol{indx},'Gebuis')
        bv.BlockNum(k) = GebFirstBlocks(Trialorder(k));
        bv.firstBlock{k} = 'Gebuis';
    elseif strcmp(Protocol{indx},'Panamath')
        bv.BlockNum(k) = PanFirstBlocks(Trialorder(k));
        bv.firstBlock{k} = 'Panamath';
    else
        keyboard
    end
end

% loop through participants and calculate beta coefficients for panamath
% and Gebuis. Also calculates the clayton exclusion criteria
ClaytOut = zeros(numel(bv.choice),1);
ClaytOutList = [];
for n = 1:numel(uniPar)
    disp(n);
    
    % which block did they get first? (ABAB or BABA)
    firstBlockPart(n) = unique(bv.firstBlock(Participant == uniPar(n)));
    
    % both blocks combined
    indx = Participant == uniPar(n) & strcmp(Protocol,'Gebuis');
    out = fit_model_data_subset_clayton(bv,indx,'full',1);
    b.Geb(:,n) = out.b;
    p.Geb(n) = out.p;
    stats.Geb{n} = out.stats;
    acc.Geb(n) = mean(Accuracy(indx));
    indx = Participant == uniPar(n) & strcmp(Protocol,'Panamath');
    out = fit_model_data_subset_clayton(bv,indx,'full',1);
    b.Pan(:,n) = out.b;
    p.Pan(n) = out.p;
    stats.Pan{n} = out.stats;
    acc.Pan(n) = mean(Accuracy(indx));
    
    % just first block
    indx = Participant == uniPar(n) & strcmp(Protocol,'Gebuis') & strcmp(Time,'Time1');
    out = fit_model_data_subset_clayton(bv,indx,'full',1);
    b.Geb1(:,n) = out.b; 
    p.Geb1(n) = out.p;
    stats.Geb1{n} = out.stats;
    acc.Geb1(n) = mean(Accuracy(indx));
    indx = Participant == uniPar(n) & strcmp(Protocol,'Gebuis') & bv.first60indx & strcmp(Time,'Time1');
    out = fit_model_data_subset_clayton(bv,indx,'full',1);
    b.Geb601(:,n) = out.b;
    indx = Participant == uniPar(n) & strcmp(Protocol,'Panamath') & strcmp(Time,'Time1');
    out = fit_model_data_subset_clayton(bv,indx,'full',1);
    b.Pan1(:,n) = out.b;
    p.Pan1(n) = out.p;
    stats.Pan1{n} = out.stats;
    acc.Pan1(n) = mean(Accuracy(indx));
    
    % just second block
    indx = Participant == uniPar(n) & strcmp(Protocol,'Gebuis') & strcmp(Time,'Time2');
    out = fit_model_data_subset_clayton(bv,indx,'full',1);
    b.Geb2(:,n) = out.b;
    p.Geb2(n) = out.p;
    stats.Geb2{n} = out.stats;
    acc.Geb2(n) = mean(Accuracy(indx));
    indx = Participant == uniPar(n) & strcmp(Protocol,'Gebuis') & bv.first60indx & strcmp(Time,'Time2');
    out = fit_model_data_subset_clayton(bv,indx,'full',1);
    b.Geb602(:,n) = out.b;
    indx = Participant == uniPar(n) & strcmp(Protocol,'Panamath') & strcmp(Time,'Time2');
    out = fit_model_data_subset_clayton(bv,indx,'full',1);
    b.Pan2(:,n) = out.b;
    p.Pan2(n) = out.p;
    stats.Pan2{n} = out.stats;
    acc.Pan2(n) = mean(Accuracy(indx));
    
    % calculate clayton exclusion criteria
    indxA = Participant == uniPar(n) & strcmp(Protocol,'Gebuis') & strcmp(Time,'Time1');
    indxB = Participant == uniPar(n) & strcmp(Protocol,'Gebuis') & strcmp(Time,'Time2');
    indxC = Participant == uniPar(n) & strcmp(Protocol,'Panamath') & strcmp(Time,'Time1');
    indxD = Participant == uniPar(n) & strcmp(Protocol,'Panamath') & strcmp(Time,'Time2');
    alph.excl = 0.03;
    if (1-cdf('binomial',sum(Accuracy(indxA)),sum(indxA),.5)) > alph.excl |...
            (1-cdf('binomial',sum(Accuracy(indxB)),sum(indxB),.5)) > alph.excl |...
            (1-cdf('binomial',sum(Accuracy(indxC)),sum(indxC),.5)) > alph.excl |...
            (1-cdf('binomial',sum(Accuracy(indxD)),sum(indxD),.5)) > alph.excl
        ClaytOut(Participant == uniPar(n)) = 1;
        ClaytOutList = [ClaytOutList,(n)];
    end
end

% calculate p values for number
for k = 1:57; 
    pNumPan(k) = stats.Pan{k}.p(2); 
    pNumGeb(k) = stats.Pan{k}.p(2);
end;

% remove outliers
% outlier paradigm is to remove anyone that had at least one block where
% the model was not significant (2) or the Beta number was not sig dif from
% 0 for all Geb or all Pan (1)
outIndx = find(p.Geb1 > 0.05 | p.Geb2 > 0.05 | p.Pan1 > 0.05 | p.Pan2 > 0.05 | pNumGeb > 0.05 | pNumGeb > 0.05)
% outIndx = ClaytOutList
if switches.removeOutliers
%    outIndx = b.Geb(2,:) < 1 | b.Pan(2,:) < 1;
    fprintf('%d outliers removed\n',numel(outIndx))
    bFields = fields(b);
    for k = 1:numel(bFields)
        b.(bFields{k})(:,outIndx) = [];
    end
    acc.Pan(:,outIndx) = [];
    acc.Geb(:,outIndx) = [];
    firstBlockPart(outIndx) = [];
end

% inter-test reliability
[rbnum pbnum] = corrcoef([b.Geb(2,:)',b.Pan(2,:)']);
[rbsize pbsize] = corrcoef([b.Geb(3,:)',b.Pan(3,:)']);
[rbspace pbspace] = corrcoef([b.Geb(4,:)',b.Pan(4,:)']);
w.Geb = 1./(sqrt(2)*b.Geb(2,:));
w.Pan = 1./(sqrt(2)*b.Pan(2,:));
[rw pw] = corrcoef(w.Geb,w.Pan); % w
[racc pacc] = corrcoef(acc.Geb,acc.Pan); % accuracy
fprintf('\nIntertest reliability:\nBetaNum,BetaSize,BetaSpace,w,Acc\nr = %0.3f,%0.3f,%0.3f,%0.3f,%0.3f\np = %0.3f,%0.3f,%0.3f,%0.3f,%0.3f\n',...
    rbnum(2),rbsize(2),rbspace(2),rw(2),racc(2),pbnum(2),pbsize(2),pbspace(2),pw(2),pacc(2))

% plot correlation of w (Figure 3)
figureDim = [0 0 9 7];
figure(3);
set(3,'PaperUnits','centimeters','PaperPosition',figureDim,...
    'Units','centimeters','Position',figureDim);
whichstats={'yhat','beta','rsquare','fstat'};
reg=regstats(w.Pan,w.Geb,'linear',whichstats);
range.x=max(w.Geb)-min(w.Geb);
range.y=max(w.Pan)-min(w.Pan);
plotRect=[min(w.Geb)-(range.x/10) max(w.Geb)+(range.x/10) min(w.Pan)-(range.y/10) max(w.Pan)+(range.y/10)];
plot(w.Geb,w.Pan,'xk');hold on;
regressLineY = [reg.beta(1)+reg.beta(2)*min(w.Geb) reg.beta(1)+reg.beta(2)*max(w.Geb)];
regressLineX = [min(w.Geb) max(w.Geb)];
plot(regressLineX,regressLineY,'-k'); 
plot([0 0.7],[0 0.7],':k');hold off;
daspect([1 1 1])
axis([0 0.7 0 0.7])
set(gca,'box','off','xtick',[0.2 0.4 0.6],'ytick',[0.2 0.4 0.6])
xlabel('wGeb');ylabel('wPan')

% test-retest reliability
[rGebNum pGebNum] = corrcoef([b.Geb1(2,:)',b.Geb2(2,:)']);
[rPanNum pPanNum] = corrcoef([b.Pan1(2,:)',b.Pan2(2,:)']);
[rGebSz pGebSz] = corrcoef([b.Geb1(3,:)',b.Geb2(3,:)']);
[rPanSz pPanSz] = corrcoef([b.Pan1(3,:)',b.Pan2(3,:)']);
[rGebSp pGebSp] = corrcoef([b.Geb1(4,:)',b.Geb2(4,:)']);
[rPanSp pPanSp] = corrcoef([b.Pan1(4,:)',b.Pan2(4,:)']);
w.Geb1 = 1./(sqrt(2)*b.Geb1(2,:));
w.Pan1 = 1./(sqrt(2)*b.Pan1(2,:));
w.Geb2 = 1./(sqrt(2)*b.Geb2(2,:));
w.Pan2 = 1./(sqrt(2)*b.Pan2(2,:));
[rGebw pGebw] = corrcoef([w.Geb1',w.Geb2']);
[rPanw pPanw] = corrcoef([w.Pan1',w.Pan2']);
fprintf('\nNum Size Space w (r then p) First row Panamath, Second Gebuis:\n')
Table2 = [rPanNum(2) pPanNum(2) rPanSz(2) pPanSz(2) rPanSp(2) pPanSp(2) rPanw(2) pPanw(2);...
    rGebNum(2) pGebNum(2) rGebSz(2) pGebSz(2) rGebSp(2) pGebSp(2) rGebw(2) pGebw(2)]
[rGeb60 pGeb60] = corrcoef([b.Geb601(2,:)',b.Geb602(2,:)']);

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
% plot participants
plot(b.Pan(3,:),b.Pan(2,:),'marker','o','color',[1 .5 .5],'linestyle','none')
plot(b.Geb(3,:),b.Geb(2,:),'marker','x','color',[.5 .5 1],'linestyle','none')
% plot means
plot(mean(b.Pan(3,:)),mean(b.Pan(2,:)),'marker','o','markersize',12,...
    'color',[0.75 0 0],'linestyle','none')
plot(mean(b.Geb(3,:)),mean(b.Geb(2,:)),'marker','x','markersize',12,...
    'color',[0 0 0.75],'linestyle','none')
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
% plot participants
plot(b.Pan(4,:),b.Pan(2,:),'marker','o','color',[1 .5 .5],'linestyle','none')
plot(b.Geb(4,:),b.Geb(2,:),'marker','x','color',[.5 .5 1],'linestyle','none')
% plot means
plot(mean(b.Pan(4,:)),mean(b.Pan(2,:)),'marker','o','markersize',12,...
    'color',[0.75 0 0],'linestyle','none')
plot(mean(b.Geb(4,:)),mean(b.Geb(2,:)),'marker','x','markersize',12,...
    'color',[0 0 0.75],'linestyle','none')
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

% Figure 2 Model fit to whole population
% model
guessRate = 1;
fl = @(mu)norminv((mu-.5)/guessRate+.5); % link
fd = @(mu)1./(guessRate*normpdf(norminv((mu-.5)/guessRate+.5))); % derivative
fi = @(mu)(guessRate*(normcdf(mu)-.5)+.5); % inverse
% fit model
indx = strcmp(Protocol,'Gebuis');
out = fit_model_data_subset_clayton(bv,indx,'full',1);
b.AllGeb = out.b;
p.AllGeb = out.p;
stats.AllGeb = out.stats;
indx = strcmp(Protocol,'Panamath');
out = fit_model_data_subset_clayton(bv,indx,'full',1);
b.AllPan = out.b;
p.AllPan = out.p;
stats.AllPan = out.stats;
% predict from model fit
indx = strcmp(Protocol,'Panamath');
X = [bv.Dnum(indx)',bv.DONSZ(indx)',bv.DONSP(indx)'];
yhat.Pan = glmval(b.AllPan,X,{fl fd fi});
yhat.All(indx) = yhat.Pan;
indx = strcmp(Protocol,'Gebuis');
X = [bv.Dnum(indx)',bv.DONSZ(indx)',bv.DONSP(indx)'];
yhat.Geb = glmval(b.AllGeb,X,{fl fd fi});
yhat.All(indx) = yhat.Geb;
% get mean acc for each ratio in each protocol
% Panamath
protIndx = strcmp(Protocol,'Panamath');
uniNum.Pan = unique(bv.Dnum(protIndx));
for k = 1:numel(uniNum.Pan)
    indx = protIndx & bv.Dnum' == uniNum.Pan(k);
    meanAcc.Pan(k) = mean(bv.choice(indx));
    meanYhat.Pan(k) = mean(yhat.All(indx));
    count.Pan(k) = sum(indx);
end
% Gebuis
protIndx = strcmp(Protocol,'Gebuis');
uniNum.Geb = unique(bv.Dnum(protIndx));
for k = 1:numel(uniNum.Geb)
    indx = protIndx & bv.Dnum' == uniNum.Geb(k);
    meanAcc.Geb(k) = mean(bv.choice(indx));
    meanYhat.Geb(k) = mean(yhat.All(indx));
    count.Geb(k) = sum(indx);
end
figureDim = [0 0 9 7];
h = figure;
set(h,'PaperUnits','centimeters','PaperPosition',figureDim,...
    'Units','centimeters','Position',figureDim);
plot(uniNum.Pan,meanAcc.Pan,'or'); hold on;
plot(uniNum.Pan,meanYhat.Pan,'-r'); 
plot(uniNum.Geb,meanAcc.Geb,'ob');
plot(uniNum.Geb,meanYhat.Geb,'-b'); 
plot([-1 1],[.5 .5],':k');hold off;
xlabel('Numerical Ratio (Left:Right)');
ylabel('P(choose right)');
set(gca,'box','off','xtick',-1:1,'xticklabel',{'2:1','1:1','1:2'})

% simulation
indx = strcmp(Protocol,'Panamath');
X = [bv.Dnum(indx)',bv.DONSZ(indx)',bv.DONSP(indx)'];
yhat.PanWGebStrat = glmval(mean(b.Pan,2),X,{fl fd fi});
indx = strcmp(Protocol,'Gebuis');
X = [bv.Dnum(indx)',bv.DONSZ(indx)',bv.DONSP(indx)'];
yhat.GebWPanStrat = glmval(mean(b.Pan,2),X,{fl fd fi});

fprintf('mean performance, predicted performance, and Gebuis strat predicted performance\n on Panamath stim set.\n')
[mean(Accuracy(strcmp(Protocol,'Panamath'))),...
    mean( (yhat.Pan < 0.5) == ...
    (bv.Dnum(strcmp(Protocol,'Panamath'))' < 0) ),...
    mean( (yhat.PanWGebStrat < 0.5) ==...
    (bv.Dnum(strcmp(Protocol,'Panamath'))' < 0) ),...
    ]
fprintf('mean performance, predicted performance, and Panamath strat predicted performance\n on Gebuis stim set.\n')
[mean(Accuracy(strcmp(Protocol,'Gebuis'))),...
    mean( (yhat.Geb < 0.5) == ...
    (bv.Dnum(strcmp(Protocol,'Gebuis'))' < 0) ),...
    mean( (yhat.GebWPanStrat < 0.5) ==...
    (bv.Dnum(strcmp(Protocol,'Gebuis'))' < 0) ),...
    ]

