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

for n = 1:numel(uniPar)
    rand50Perc = round(rand(numel(bv.Dnum),1));
    
    
end