function Dp=Zp2Dp(Zp,T,pressure,chg)

% Dp=Zp2Dp(Zp,T,pressure,chg)
% Calculate Dp from Zp 
% Zp: in SI unit
% pressure: in atm
% T: in K
% chg: unit charge carried by particle

global const
const=get_const;
viscosity=1.8334e-5*power(T/296.15,1.5)*(296.15+110.4)/(T+110.4);

mfp = const.mfp_1atm_298K*(T/296.15)/(pressure/1.0)*(1+110.4/296.15)/(1+110.4/T);
nan_idx=isnan(Zp);
Zp(nan_idx)=1e-6;
for j1=1:length(Zp)
    Dp_high=6.628*mfp/(sqrt(1+4*3.314*mfp*(3*pi*viscosity*Zp(j1)/(const.e*chg)))-1);
    Dp_low=chg*const.e/(3.0*pi*viscosity*Zp(j1));
    
    Dp_Cc=chg*const.e/(3.0*pi*viscosity*Zp(j1));
    
    while 2>0
        Dp_mid=(Dp_high+Dp_low)/2.0;
        Cc=1.0+(2.0*mfp/Dp_mid)*(1.257+0.4*exp(-1.1*Dp_mid/(2.0*mfp)));
        
        if (abs(((Dp_mid/Cc)-Dp_Cc)/Dp_Cc))<0.001
            break;
        elseif (Dp_mid/Cc)>Dp_Cc
            Dp_high=Dp_mid;
        else
            Dp_low=Dp_mid;
        end        
    end
    
    Dp(j1)=Dp_mid;
end
Dp(nan_idx)=NaN;