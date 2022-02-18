function F=GFPDF_pie(D2,D1,guess,x)
%%% linspace(a,b,n) generates "n" equally spaced numbers btw a and b
%%% including a and b
%%% the unknown numbers in the GFPDF is n-2, because c at 0.7 and 1.9 are 0
% global range
% x=range;
rangeend=[0 10];
x=cat(2,x,rangeend);
y=zeros(1,length(guess)+4);
for i=1:1:length(guess)
    y(i+1)=guess(i);
end
ratio=D2./D1;
F=interp1(x,y,ratio);