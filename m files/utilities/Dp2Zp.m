function Zp=Dp2Zp(Dp,T,pressure,chg)
% Calculate Zp from particle diameter
% Zp: si unit
% Dp: meter
% pressure: 1atm
% T: K
% chg: unit charge carried by particle


global const
const=get_const;

viscosity=1.8334e-5*power(T/296.15,1.5)*(296.15+110.4)/(T+110.4);

mfp = const.mfp_1atm_298K*(T/296.15)/(pressure/1.0)*(1+110.4/296.15)/(1+110.4/T);

Cc = 1.0 + (2.0*mfp./Dp).*(1.257 + 0.4*exp(-1.1*Dp/(2.0*mfp)));
Zp = chg*const.e*Cc./(3.0*pi*viscosity*Dp);