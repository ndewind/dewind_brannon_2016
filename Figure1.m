% Plot Panamath, Gebuis, and DeWind pair ratios in 3D ratio space

%% add the dewind stimuli.
load('559_pairs_info');
load('bv');

% flip side and add features
bv.Dnum(bv.S1side == 0) = -bv.Dnum(bv.S1side == 0);
bv.DONSZ(bv.S1side == 0) = -bv.DONSZ(bv.S1side == 0);
bv.DONSP(bv.S1side == 0) = -bv.DONSP(bv.S1side == 0);
bv.Disa = 0.5*bv.DONSZ - 0.5*bv.Dnum;
bv.Dtsa = 0.5*bv.DONSZ + 0.5*bv.Dnum;
bv.Dspar = 0.5*bv.DONSP - 0.5*bv.Dnum;
bv.Dch = 0.5*bv.DONSP + 0.5*bv.Dnum;

% identify unique pairings
uniRat = [bv.Dnum(1),bv.DONSZ(1),bv.DONSP(1)];
uniIndx = ones(numel(bv.Dnum),1);
for k = 2:numel(bv.Dnum)
    isRepeat = false;
    disp(k)
    for m = 1:size(uniRat,1)
        if bv.Dnum(k) == uniRat(m,1) &...
                bv.DONSZ(k) == uniRat(m,2) &...
                bv.DONSP(k) == uniRat(m,3)
            isRepeat = true;
            break
        end
    end
    if ~isRepeat
        uniRat(end+1,:) = ([bv.Dnum(k),bv.DONSZ(k),bv.DONSP(k)]);
        uniIndx(k) = 1;
    else
        uniIndx(k) = 0;
    end
end

% plot
axisLabel = {'16:1','8:1','4:1','2:1','1:1','1:2','1:4','1:8','1:16'};
figure(1);
indx = logical(uniIndx);
h=subplot(1,2,1);
plot(bv.Disa(indx),bv.Dtsa(indx),'marker','^','color',[0 .75 0],'linestyle','none','markersize',2);hold on;
h=subplot(1,2,2);
plot(bv.Dspar(indx),bv.Dch(indx),'marker','^','color',[0 .75 0],'linestyle','none','markersize',2); hold on;


% correlation matrix for DeWind stimuli
fprintf('Correlation Matrix for DeWind Stimuli\n')
[r.dew p] = ...
    corrcoef([bv.Dnum',bv.Dtsa',bv.Disa',bv.Dch',bv.Dspar',bv.DONSZ',bv.DONSP']);

indx = true(1,numel(bv.correct));
fprintf('Table 1 DeWind:\n')
[r.dew p] = ...
    corrcoef([bv.Dnum(indx)',bv.Dtsa(indx)',bv.Disa(indx)',bv.Dch(indx)',...
    bv.Dspar(indx)',bv.DONSZ(indx)',bv.DONSP(indx)']);
mins.dew = 2.^[min(abs(bv.Dnum(indx))),min(abs(bv.Dtsa(indx))),min(abs(bv.Disa(indx))),min(abs(bv.Dch(indx))),...
    min(abs(bv.Dspar(indx))),min(abs(bv.DONSZ(indx))),min(abs(bv.DONSP(indx)))]
maxs.dew = 2.^[max(abs(bv.Dnum(indx))),max(abs(bv.Dtsa(indx))),max(abs(bv.Disa(indx))),max(abs(bv.Dch(indx))),...
    max(abs(bv.Dspar(indx))),max(abs(bv.DONSZ(indx))),max(abs(bv.DONSP(indx)))]
means.dew = 2.^[mean(abs(bv.Dnum(indx))),mean(abs(bv.Dtsa(indx))),mean(abs(bv.Disa(indx))),mean(abs(bv.Dch(indx))),...
    mean(abs(bv.Dspar(indx))),mean(abs(bv.DONSZ(indx))),mean(abs(bv.DONSP(indx)))]
stds.dew = 2.^[std(abs(bv.Dnum(indx))),std(abs(bv.Dtsa(indx))),std(abs(bv.Disa(indx))),std(abs(bv.Dch(indx))),...
    std(abs(bv.Dspar(indx))),std(abs(bv.DONSZ(indx))),std(abs(bv.DONSP(indx)))]

% keyboard

clear;


%% load gebuis and panamath %%%%%%%%%%%%%%%
load('ClaytonDataImport')

% create size and spacing variables
SizeLeft = pixels1_Left .* DotSize1_Left;
SizeRight = pixels2_Right .* DotSize2_Right;
SparLeft = (ConvexHull1_Left ./ Numerosity1_Left);
SparRight = (ConvexHull2_Right ./ Numerosity2_Right);
SpaceLeft = SparLeft .* ConvexHull1_Left;
SpaceRight = SparRight .* ConvexHull2_Right;

% create right to left ratios and side chosen variables.
bv.Dnum = -log2(LeftRightNumerosityRatio)';
bv.Dsize = -log2(SizeLeft ./ SizeRight)';
bv.Dspace = -log2(SpaceLeft ./ SpaceRight)';
bv.Dtsa = -log2(pixels1_Left ./ pixels2_Right)';
bv.Disa = -log2(DotSize1_Left ./ DotSize2_Right)';
bv.Dspar = -log2(SparLeft ./ SparRight)';
bv.Dch = -log2(ConvexHull1_Left ./ ConvexHull2_Right)';

% identify unique pairings
uniRat = [bv.Dnum(1),bv.Dsize(1),bv.Dspace(1)];
uniProt{1} = Protocol{1};
uniIndx = ones(numel(bv.Dnum),1);
for k = 2:numel(bv.Dnum)
    isRepeat = false;
    disp(k)
    for m = 1:size(uniRat,1)
        if bv.Dnum(k) == uniRat(m,1) &...
                bv.Dsize(k) == uniRat(m,2) &...
                bv.Dspace(k) == uniRat(m,3)
            isRepeat = true;
            break
        end
    end
    if ~isRepeat
        uniRat(end+1,:) = ([bv.Dnum(k),bv.Dsize(k),bv.Dspace(k)]);
        uniProt{end+1} = Protocol{k};
        uniIndx(k) = 1;
    else
        uniIndx(k) = 0;
    end
end

% plot
axisLabel = {'16:1','8:1','4:1','2:1','1:1','1:2','1:4','1:8','1:16'};
figureDim = [0 0 19 10];
figure(1);
set(1,'PaperUnits','centimeters','PaperPosition',figureDim,...
    'Units','centimeters','Position',figureDim);
h=subplot(1,2,1);
plot([-4 1.75],[-4 1.75],':k');hold on;
indx = uniIndx & strcmp(Protocol,'Gebuis');
plot(bv.Disa(indx),bv.Dtsa(indx),'xb','markersize',2); 
indx = uniIndx & strcmp(Protocol,'Panamath');
plot(bv.Disa(indx),bv.Dtsa(indx),'or','markersize',2);
axis([-4 4 -4 4])
daspect([1 1 1])
xlabel('Item area ratio');ylabel('Total area ratio');
set(h,'box','off',...
    'xtick',-4:4,'xticklabel',axisLabel,...
    'ytick',-4:4,'yticklabel',axisLabel);

h=subplot(1,2,2);
plot([-4 1.75],[-4 1.75],':k');hold on;
indx = uniIndx & strcmp(Protocol,'Gebuis');
plot(bv.Dspar(indx),bv.Dch(indx),'xb','markersize',2); 
indx = uniIndx & strcmp(Protocol,'Panamath');
plot(bv.Dspar(indx),bv.Dch(indx),'or','markersize',2);
axis([-4 4 -4 4])
daspect([1 1 1])
xlabel('Sparsity ratio');ylabel('Convex hull ratio');
set(h,'box','off',...
    'xtick',-4:4,'xticklabel',axisLabel,...
    'ytick',-4:4,'yticklabel',axisLabel);

% Table 1
indx = strcmp(Protocol,'Panamath');
fprintf('Table 1 Panamath Stimuli\n')
[r.pan p] = ...
    corrcoef([bv.Dnum(indx)',bv.Dtsa(indx)',bv.Disa(indx)',bv.Dch(indx)',...
    bv.Dspar(indx)',bv.Dsize(indx)',bv.Dspace(indx)']);
mins.pan = 2.^[min(abs(bv.Dnum(indx))),min(abs(bv.Dtsa(indx))),min(abs(bv.Disa(indx))),min(abs(bv.Dch(indx))),...
    min(abs(bv.Dspar(indx))),min(abs(bv.Dsize(indx))),min(abs(bv.Dspace(indx)))]
maxs.pan = 2.^[max(abs(bv.Dnum(indx))),max(abs(bv.Dtsa(indx))),max(abs(bv.Disa(indx))),max(abs(bv.Dch(indx))),...
    max(abs(bv.Dspar(indx))),max(abs(bv.Dsize(indx))),max(abs(bv.Dspace(indx)))]
means.pan = 2.^[mean(abs(bv.Dnum(indx))),mean(abs(bv.Dtsa(indx))),mean(abs(bv.Disa(indx))),mean(abs(bv.Dch(indx))),...
    mean(abs(bv.Dspar(indx))),mean(abs(bv.Dsize(indx))),mean(abs(bv.Dspace(indx)))]
stds.pan = 2.^[std(abs(bv.Dnum(indx))),std(abs(bv.Dtsa(indx))),std(abs(bv.Disa(indx))),std(abs(bv.Dch(indx))),...
    std(abs(bv.Dspar(indx))),std(abs(bv.Dsize(indx))),std(abs(bv.Dspace(indx)))]


indx = strcmp(Protocol,'Gebuis');
fprintf('Correlation Matrix for Gebuis Stimuli\n')
[r.geb p] = ...
    corrcoef([bv.Dnum(indx)',bv.Dtsa(indx)',bv.Disa(indx)',bv.Dch(indx)',...
    bv.Dspar(indx)',bv.Dsize(indx)',bv.Dspace(indx)']);
mins.geb = 2.^[min(abs(bv.Dnum(indx))),min(abs(bv.Dtsa(indx))),min(abs(bv.Disa(indx))),min(abs(bv.Dch(indx))),...
    min(abs(bv.Dspar(indx))),min(abs(bv.Dsize(indx))),min(abs(bv.Dspace(indx)))]
maxs.geb = 2.^[max(abs(bv.Dnum(indx))),max(abs(bv.Dtsa(indx))),max(abs(bv.Disa(indx))),max(abs(bv.Dch(indx))),...
    max(abs(bv.Dspar(indx))),max(abs(bv.Dsize(indx))),max(abs(bv.Dspace(indx)))]
means.geb = 2.^[mean(abs(bv.Dnum(indx))),mean(abs(bv.Dtsa(indx))),mean(abs(bv.Disa(indx))),mean(abs(bv.Dch(indx))),...
    mean(abs(bv.Dspar(indx))),mean(abs(bv.Dsize(indx))),mean(abs(bv.Dspace(indx)))]
stds.geb = 2.^[std(abs(bv.Dnum(indx))),std(abs(bv.Dtsa(indx))),std(abs(bv.Disa(indx))),std(abs(bv.Dch(indx))),...
    std(abs(bv.Dspar(indx))),std(abs(bv.Dsize(indx))),std(abs(bv.Dspace(indx)))]

