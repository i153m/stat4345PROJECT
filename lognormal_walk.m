% Function to simulate lognormal random walk model of stock prices
function S = lognormal_walk(S_0, mu, sigma, N)
X = normrnd(mu, sigma, 1, N);
S = S_0 .* exp(cumsum(X));
end

% Function to implement option exercising policy and calculate payoff
function payoff = option_payoff(P, K, m, mu, sigma)
alpha = mu + sigma^2/2;
b = (mmu - log(K/P))/(sigma * sqrt(m));
term1 = P * exp(alpha * (1:m) + sigma * sqrt(1:m) * (b + sigma^2/(2m)));
term2 = K * normcdf(b + sigmasqrt(1:m)) - P * normcdf(sigmasqrt(1:m) + b);
exercise_now = any(term1 > term2);
if (exercise_now)
payoff = P;
else
% Simulate remaining days until option expiration
remaining_days = m - 1;
while (remaining_days > 0)
P_new = P * exp(normrnd(mu, sigma));
if (P_new > K)
payoff = P_new;
break
else
P = P_new;
remaining_days = remaining_days - 1;
end
if (remaining_days == 0)
payoff = K;
end
end
end
payoff = payoff - K;
end

% Set parameters
S_0 = 100; % initial stock price
mu = 0.01; % daily mean return
sigma = 0.2; % daily standard deviation of returns
N = 30; % time period during which option can be exercised
K = 110; % striking price

% Initialize total payoff
payoff_total = 0;

% Simulate many scenarios
for i = 1:10000
% Simulate stock price
S = lognormal_walk(S_0, mu, sigma, N);

% Loop over a number of days remaining until option expiration and determine whether to exercise
m = N;
while m > 0
P = S(N-m+1);
if P > K
if option_payoff(P, K, m, mu, sigma) > 0
payoff_total = payoff_total + P - K;
break
else
m = m - 1;
end
else
m = m - 1;
end
end
end

% Calculate the expected payoff per option
expected_payoff = payoff_total / 10000;