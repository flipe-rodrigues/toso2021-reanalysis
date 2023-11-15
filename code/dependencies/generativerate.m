function x = generativerate(time,gamma,mu,eta,sigma)
    mu_sample = clamp(normrnd(mu,eta),time(1),time(end));
    x_pdf = normpdf(time,mu_sample,sigma);
    x = gamma * x_pdf ./ max(x_pdf);
end