function [out,guessRate] = fit_model_data_subset(bv,indx,model,g)
% fits the full model to a subset of the data.  If a guess rate is supplied
% it will use that, otherwise it will calculate the guess rate

% trim bv down using the indx
binom = bv;
dataFields = fields(binom);
for m = 1:length(dataFields)
    binom.(dataFields{m}) = binom.(dataFields{m})(indx);
end

% if no model is specified the full model with size and spacing is
% calculated.
if ~exist('model','var')
    if ~ischar(model)
        model = 'full';
    end
end

% calculates guessRate or assigns it based on input
if exist('g','var')
    guessRate = g;
else
    error('Error. Provide a ''g''');
%     if strcmp(model,'full')
%         % Calculate guess rate
%         guessRate = calc_guessing_rate_fxn(data);
%     else
%         guessRate = [];
%     end
end

% fit model
switch model
    case 'full'
        fl = @(mu)norminv((mu-.5)/guessRate+.5); % link
        fd = @(mu)1./(guessRate*normpdf(norminv((mu-.5)/guessRate+.5))); % derivative
        fi = @(mu)(guessRate*(normcdf(mu)-.5)+.5); % inverse
        X = [binom.Dnum' binom.DONSZ' binom.DONSP'];
        [out.b,out.dev,out.stats] = glmfit(X,binom.choice','binomial','link',{fl fd fi});
        [b,dev,stats] = glmfit(ones(numel(binom.choice),1),binom.choice','binomial','link',{fl fd fi},'constant','off');
        out.p = 1-chi2cdf(dev-out.dev,3); % whole model significance
    case 'log'
        % piazza fuzzy model
        piazzaModelFuzzy.model = @(w,X)(1/2) * (1 + erf(X/(sqrt(4*w^2))));
        [out.w, r, j] = nlinfit(binom.Dnum',binom.choice',piazzaModelFuzzy.model,.2);
        if strfind(version,'(R14)')
            out.ci = nlparci(out.w,r,j);
        else
            out.ci = nlparci(out.w,r,'jacobian',j);
        end
        
    case 'lin'
        % halberda model
        halberdaModel.model = @(w,X)(1/2) * (1 + erf((X-1) ./ (sqrt(2) * w * sqrt(X.^2+1))));
        [out.w, r, j] = nlinfit(2.^binom.Dnum',binom.choice',halberdaModel.model,.2);
        if strfind(version,'(R14)')
            out.ci = nlparci(out.w,r,j);
        else
            out.ci = nlparci(out.w,r,'jacobian',j);
        end
        
    otherwise
        error('bad model string')
end