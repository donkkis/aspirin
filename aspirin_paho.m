%   Statistical Methods / TTY Pori 2018
%   Ex. 6-5
%
%   Confidence intervals for aspirin / placebo patients sustaining
%   myocardial infarction
%   P. Aho 8.5.2018
%
%   Originally published in:
%
%   Steering Committee of the Physicians' Health Study Research Group.
%   1989. Final Report on the Aspirin Component of the Ongoing Physicians'
%   Health Study. Available at https://www.nejm.org/doi/full/10.1056/NEJM198907203210301
%
%   

%%
% Setup data

clear; clc;

infarction = [139; 239];
no_infarction = [11037-139; 11034-239];

% Using a 95 % confidence level according to original publication
alfa = 0.05;

data = table(infarction, no_infarction, 'RowNames', {'Aspirin', 'Placebo'});

%%
% Method 1: Binomial fit
% Fit binomial distributions to aspirin and placebo populations
% respectfully
%
% phat_XX = maximum likelihood estimate of a random patient in the 
% relevant population sustaining heart attack
%
% ci_XX = 95 % confidence interval for phat

[phat_aspirin, ci_aspirin_bf] = binofit(...
    data{'Aspirin', 'infarction'},...
    sum(data{'Aspirin', :}));

[phat_placebo, ci_placebo_bf] = binofit(...
    data{'Placebo', 'infarction'},...
    sum(data{'Placebo', :}));

%%
% Method 2: Agresti-Coull
% 
% Agresti, A. and Coull, B. A. 1998. Approximate is better than exact 
% for interval estimation of binomial proportions. The American Statistician, 
% 52(2), 119-126.
%
% Slightly modified version of computing the CI's for binomial distribution's p parameter 
% as proposed originally by Agresti & Coull. Summary of the method can be found here: 
% https://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/agcoulci.htm

%1 - alpha/2 quantile in standard normal distribution
z = norminv(1 - alfa/2);

%successes
X = data{'Aspirin', 'infarction'};
%trials
n = sum(data{'Aspirin', :});
%proportion of successes in observed data
p_ = X/n;

%95 % ci for the aspirin group's p
ci_aspirin_ag = [
    (p_ + z^2/(2*n) - z * sqrt(p_/n*(1-p_) + z^2/(4*n^2))) / (1+z^2/n);...
    (p_ + z^2/(2*n) + z * sqrt(p_/n*(1-p_) + z^2/(4*n^2))) / (1+z^2/n)];

%successes
X = data{'Placebo', 'infarction'};
%trials
n = sum(data{'Placebo', :});
%proportion of successes in observed data
p_ = X/n;

% 95 % ci for the placebo groups p
ci_placebo_ag = [
    (p_ + z^2/(2*n) - z * sqrt(p_/n*(1-p_) + z^2/(4*n^2))) / (1+z^2/n);...
    (p_ + z^2/(2*n) + z * sqrt(p_/n*(1-p_) + z^2/(4*n^2))) / (1+z^2/n)];

%%
% Method 3: Bootstrapping
% Efron, B., & Tibshirani, R. J. 1994. An introduction to the bootstrap. 
% Chapman & hall/crc monographs on statistics & applied probability.

prc_level = [2.5 97.5];
max_rounds = 10000;

na_heart = 139; na_total = 11037;
np_heart = 239; np_total = 11034;

% Encode populations to one-hot vectors: 1 = heart attack, 0 = no heart
% attack
popula_aspirin = [ones(1,na_heart),zeros(1,na_total-na_heart)];
popula_placebo = [ones(1,np_heart),zeros(1,np_total-np_heart)];

% At each bs round, take a random sample with replacement from each population. 
% Compute proportions of heart attack and relative risk ratio between the
% populations.
for i = 1:max_rounds
    na = sum(datasample(popula_aspirin,na_total));
    np = sum(datasample(popula_placebo,np_total));
    p_asp(i) = na / na_total;
    p_pla(i) = np / np_total;
    theta(i) = (na / na_total) / (np / np_total);
end

% compute the CI's by taking 2.5 and 97.5 percentiles from bs result vector
ci_aspirin_bs = prctile(p_asp, prc_level);
ci_placebo_bs = prctile(p_pla, prc_level); 
ratio_ci = prctile(theta,prc_level);

% a point estimate for relative risk ratio
ratio_est = (na_heart / na_total) / (np_heart / np_total);