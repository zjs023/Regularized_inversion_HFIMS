function Dp_set=Dp2Dp_set(Dp)

len=length(Dp);
Dp_set(1,2:len)=exp(0.5*(log(Dp(1:len-1))+log(Dp(2:len))));
Dp_set(1,1)=exp( log(Dp(1))-0.5*(log(Dp(2))-log(Dp(1))));
Dp_set(3,1:len-1)=Dp_set(1,2:len);
Dp_set(3,len)=exp( log(Dp(len))+0.5*(log(Dp(len))-log(Dp(len-1))));
Dp_set(2,:)=Dp;
