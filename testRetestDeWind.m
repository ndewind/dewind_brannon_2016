load('bv.mat')

% generate overall trial number
bv.subNum = zeros(1,numel(bv.sub));
for k = 1:numel(bv.sub)
    disp(k)
    if strcmp(bv.sub{k},'2d')
        bv.subNum(k) = 2;
    else
        bv.subNum(k) = str2double(bv.sub{k});
    end
end
bv.blockNum = zeros(1,numel(bv.sub));
for k = 1:numel(bv.subNum)
    disp(k)
    bv.blockNum(k) = sum(bv.tnum(1:k) == bv.tnum(k) & ...
        bv.subNum(1:k) == bv.subNum(k));
end
bv.tnumAll = bv.tnum + bv.blockNum*250 - 250;

% prepare binom variables
binom.Dnum(bv.S1side == 0) = -bv.Dnum(bv.S1side == 0);
binom.Dnum(bv.S1side == 1) = bv.Dnum(bv.S1side == 1);
binom.DONSZ(bv.S1side == 0) = -bv.DONSZ(bv.S1side == 0);
binom.DONSZ(bv.S1side == 1) = bv.DONSZ(bv.S1side == 1);
binom.DONSP(bv.S1side == 0) = -bv.DONSP(bv.S1side == 0);
binom.DONSP(bv.S1side == 1) = bv.DONSP(bv.S1side == 1);
binom.choice(bv.respKey == 37) = 0;
binom.choice(bv.respKey == 39) = 1;

uniSub = unique(bv.subNum);
for n = 1:numel(uniSub)
    % first "block"
    indx = bv.tnumAll <= 375 & bv.subNum == uniSub(n);
    out = fit_model_data_subset_clayton(binom,indx,'full',1);
    b.Dew1(:,n) = out.b;
    p.Dew1(n) = out.p;
    stats.Dew1{n} = out.stats;
    % second "block"
    indx = bv.tnumAll >= 376 & bv.tnumAll <= 750 & bv.subNum == uniSub(n);
    out = fit_model_data_subset_clayton(binom,indx,'full',1);
    b.Dew2(:,n) = out.b;
    p.Dew2(n) = out.p;
    stats.Dew2{n} = out.stats;
end
[rNum pNum] = corrcoef(b.Dew1(2,:),b.Dew2(2,:));
[rSz pSz] = corrcoef(b.Dew1(3,:),b.Dew2(3,:));
[rSp pSp] = corrcoef(b.Dew1(4,:),b.Dew2(4,:));
[rNum pNum] = corrcoef(b.Dew1(2,:),b.Dew2(2,:));
w.dew1 = 1./(sqrt(2)*b.Dew1(2,:));
w.dew2 = 1./(sqrt(2)*b.Dew2(2,:));
[rw pw] = corrcoef(w.dew1,w.dew2); % w

fprintf('\nDeWind test-retest correlation (375 trial artificial blocks).\n First row r second row p.\n Bnum Bsize Bspace w\n')
fprintf('%0.3f %0.3f %0.3f %0.3f\n',rNum(2),rSz(2),rSp(2),rw(2))
fprintf('%0.3f %0.3f %0.3f %0.3f\n',pNum(2),pSz(2),pSp(2),pw(2))
