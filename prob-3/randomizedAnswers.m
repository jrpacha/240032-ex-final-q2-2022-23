function fake = randomizedAnswers(x,n)
%      x: correct answers
%      n: number of faked answers you need
%
%  fakeX: fake(1),...,fake(n): faked (randomized?) answers
%

% Lema del día
% "Fake it until you get it"
%                 Elon Musk
% :-)

rng('shuffle');

log10AbsX = log10(abs(x));
exponent = floor(log10(abs(x)));

mantissa = log10AbsX-exponent;

mantissa = 10^mantissa;

fake = (mantissa - 0.5 + rand(n,1))*10^exponent;

if x < 0
    fake = -fake;
end

end