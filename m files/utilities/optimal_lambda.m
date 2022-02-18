function [lambda_opt] = optimal_lambda(Gamma, R, lambda, L, GFbinc, c_real)
    for i=1:length(lambda)
        x(i,:) = ncsolve(Gamma, R, lambda(i), L);
        x(i,:) = x(i,:)./(sum(x(i,:))*(GFbinc(2)-GFbinc(1)));
        chi(i) = norm(x(i,:)-c_real);
    end

[~,ind_optimal] = min(chi);
lambda_opt = lambda(ind_optimal);
% figure;plot(lambda,chi);hold on;plot(lambda_opt,min(chi),'ro');
end