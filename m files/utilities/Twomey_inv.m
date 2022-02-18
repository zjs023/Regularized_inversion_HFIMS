function n=Twomey_inv(Gamma,R)

% This function find x that satisfys y=m*x using Twomey method
% m: matrix
% x0: initial guess
% ite: number of iteration

% Define max iteration numbers
max_iter=1000;

if all(R==0)
    n=0;
    return
end

R=R(:);  % Turn inton column vector
len=length(R);

% R_nonzero=R(idx); len_nonzero=length(R_nonzero);
error=sqrt(R);
idx=(error==0);
error(idx)=0.001;
n=lsqnonneg(Gamma,R);
n=smoothing(n,1/3);
% n(1)=0;
% n=0.5*n;

% if mode==1
%     factor=1;
% elseif mode==2;
%     factor=max(max(max(m)),1);
%     m=m/factor;
% end


chi_sqrd=1000;
chi_sqrd_reduction=1;
iter=0;

% original criteria:
while (abs(chi_sqrd_reduction)>0.001)&(chi_sqrd>1)&(iter<=max_iter)
% note that abs(chi_sqrd_reduction)<0.01 is more difficult to reach
% the original criteria has been changed to
% while iter<=max_iter

% to test the influence of iteration on data inversion
%     plot(Dp_array,n);
%     set(gca,'xscale','log')
%     xlabel('Dp (nm)');
%     ylabel('dN/dlogDp');
%     title(['iteration',num2str(iter)]);

    X=R./(Gamma*n);
    
    %for j1=1:len
    for j1=[find(~isnan(X')& ~isinf(X'))]
        n=n.*( (1.0+(X(j1)-1.0)*Gamma(j1,:))' );
    end
    
    R_new=Gamma*n;
    chi_sqrd_new=sum( power((R_new-R)./error, 2) )/...
        len;
    chi_sqrd_reduction=(chi_sqrd-chi_sqrd_new)./chi_sqrd_new;
    
    chi_sqrd=chi_sqrd_new;
    iter=iter+1;
    
end

if iter==max_iter+1
    warning('max iteration reached\n');
end
iter

function n1=smoothing(n,b)
n1=n;  % initilize to the same dimension of n
n1(1)=(1-b)*n(1)+b*n(2);
n1(end)=(1-b)*n(end)+b*n(end-1);
n1(2:end-1)=b*n(1:end-2)+(1-2*b)*n(2:end-1)+b*n(3:end);