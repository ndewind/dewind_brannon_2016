load('ClaytonDataImport')

indx = strcmp(Protocol,'Panamath');
uniNum = unique(LeftRightNumerosityRatio(indx));
for k = 1:numel(uniNum)
    indx = strcmp(Protocol,'Panamath') & LeftRightNumerosityRatio == uniNum(k);
    meanAcc(k) = mean(Accuracy(indx));
end
plot(log2(uniNum),meanAcc)