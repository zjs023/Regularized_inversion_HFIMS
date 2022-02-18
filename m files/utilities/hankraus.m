function [lambda_opt] = hankraus(Gamma, R, lambda, L, GFbinc)
    for i=1:length(lambda)
        x(i,:) = ncsolve(Gamma, R, lambda(i), L);
        x(i,:) = x(i,:)./(sum(x(i,:))*(GFbinc(2)-GFbinc(1)));
        Axb(i) = norm(Gamma*x(i,:)'./max(Gamma*x(i,:)')-R)^2;
        zeta(i) = norm(Gamma*x(i,:)'-R,2);
        eta(i) = norm(x(i,:),2);
    end

[~,ind_hankeraus] = min(Axb./lambda);
lambda_opt = lambda(ind_hankeraus);
% figure;plot(lambda,Axb./lambda);hold on;plot(lambda_opt,min(Axb./lambda),'ro');
end