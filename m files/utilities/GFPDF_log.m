function F=GFPDF_log(D2,D1,fcond,G,sigma)
term_num=length(fcond);
F=0;
for j1=1:1:term_num
    coef=fcond(j1)/(2*pi*(log(sigma(j1))).^2)^0.5;
    term=coef.*exp(-0.5*(log(D2./D1./G(j1))).^2./(log(sigma(j1))).^2);
    F=F+term;
end
