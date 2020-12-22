function out = EffectiveOscillation(t,y)
global lambda kappa eps N gamma D a td b n
hp = @(x) (lambda/kappa)*(1/(x))*(exp(-kappa*(N-x)) - exp(-kappa*(N)));
hm = @(x) (lambda/kappa)*(1/(N - x))*(exp(-kappa*(x)) - exp(-kappa*(N)));
prot1 = @(x) a*(0.5 - 0.5*tanh(n*(x-N/4))) + 0.5 - 0.5*tanh(n*(x-3*N/4)); %a*(1 - heaviside(x - N/4))+ 1 - heaviside(x - 3*N/4);
k1    = @(x) (0.5 + 0.5*tanh(n*(x-N/4))) + b*(0.5 - 0.5*tanh(n*(x-3*N/4))); %heaviside(x-N/4) + b*(heaviside(x - 3*N/4))
prot2 = @(x) (0.5 + 0.5*tanh(n*(x-N/4))) + a*(0.5 + 0.5*tanh(n*(x-3*N/4))); %heaviside(x-N/4) + a*(heaviside(x - 3*N/4))
k2    = @(x) 0.5 - 0.5*tanh(n*(x-3*N/4)) + b*(0.5 - 0.5*tanh(n*(x-N/4))); % 1 - heaviside(x-3*N/4) + b*(1-heaviside(x-N/4))
out = [(hp(y(1))*y(2) - hm(y(1))*y(3));(1/td)*(prot1(y(4))*(1-y(2)) - k1(y(4))*y(2));(1/td)*(prot2(y(4))*(1-y(3)) - k2(y(4))*y(3));(1/td)*(y(1)-y(4))];
