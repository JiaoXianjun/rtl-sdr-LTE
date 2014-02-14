%  Returns the value x such that chi2cdf(x,k) is p.
function x = chi2cdf_inv(p, k)
x = 50:2000;
y = chi2cdf(x, k);
[~, min_idx] = min(abs(y-p));
x = x(min_idx);
