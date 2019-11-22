function e = computeExcitationRaasch(a, vA, tauDeact, tauAct)
% Computes muscle excitation from muscle activation a and derivative of
% activation vA.

% See De Groote et al. (2009) equation (17) (Raasch)

td = (ones(size(a,1),1)*tauDeact);
ta = (ones(size(a,1),1)*tauAct);

e = zeros(size(a));
e(vA<=0) = td(vA<=0) .* vA(vA<=0) + a(vA<=0);

c1 = 1./ta - 1./td;
c2 = 1./td;
D = (c2 + c1 .* a).^2 + 4*c1.*vA;
e(vA>0) = (a(vA>0) .* c1(vA>0) - c2(vA>0) + sqrt(D(vA>0)))./(2*c1(vA>0));

end

