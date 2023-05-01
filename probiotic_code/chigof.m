
%%
%
% Get chi-square p-value for the comparison of two categories (can have
% multiple groups to compare)
%
% INPUT:
%   * observed: contigency table with frequencies for each outcome in each
%       group
%
% OUTPUT:
%   * p: p-value
%   * chi2stat: chi-square statistic
%
% Example:
% observed = [100 100;
%     10 100;
%     80 800;
%     10 100];
% 
% [p, chi2stat,df] = chigof(observed);

%%


%%
function [p, chi2stat,df] = chigof(observed)
    N = sum(observed);
    p0 = sum(observed,2)/sum(observed,"all");
    expected = N.*p0;
    
    df = (size(observed,1) - 1)*(size(observed,2)-1);
    chi2stat = sum((reshape(observed,[],1)-reshape(expected,[],1)).^2 ./ reshape(expected,[],1));
    p = 1 - chi2cdf(chi2stat,df);
end

%%
% function [p, chi2stat] = chi2gof(ns,Ns)
%     % Pooled estimate of proportion
%     p0 = sum(ns) / sum(Ns);
%     % Expected counts under H0 (null hypothesis)
%     n0s = NaN(size(ns));
%     observed = NaN(1,length(ns)*2);
%     expected = NaN(1,length(ns)*2);
%     c = 1;
%     for i = 1:length(ns)
%         n0s(i) = Ns(i)*p0;
%         observed(c:c+1) = [ns(i) Ns(i)-ns(i)];
%         expected(c:c+1) = [n0s(i) Ns(i)-n0s(i)];
%         c = c + 2;
%     end
% 
%     % observed/expected
%     df = length(Ns) - 1;
%     chi2stat = sum((observed-expected).^2 ./ expected);
%     p = 1 - chi2cdf(chi2stat,df);
% end
%         