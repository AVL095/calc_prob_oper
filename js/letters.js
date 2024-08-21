//X ⁰ ¹ ² ³ ⁴ ⁵ ⁶ ⁷ ⁸ ⁹ ⁺ ⁻ ⁼ ⁽ ⁾ ᵃ ᵇ ᶜ ᵈ ᵉ ᶠ ᶢ ʰ ⁱ ʲ ᵏ ᵐ ⁿ ᵒ ᵖ ʳ ˢ ᵗ ᵘ ᵛ ʷ ˣ ʸ ᶻ ᴬ ᴮ ᴰ ᴱ ᴳ ᴴ ᴵ ᴶ ᴷ ᴸ ᴹ ᴺ ᴼ ᴾ ᴿ ᵀ ᵁ ᵂ
//X ₀ ₁ ₂ ₃ ₄ ₅ ₆ ₇ ₈ ₉ ₊ ₋ ₌ ₍ ₎ ₐ ₑ ᵢ ₒ ᵣ ᵤ ᵥ ₓ
//********************subscripts and superscripts lts and digits**************************
function lts() {
var ij;
lt={
 sb:{a:'0x2090',e:'0x2091',x:'0x2093',r:'0x1D63',
u:'0x1D64',v:'0x1D65',l:'0x2097',m:'0x2098',t:'0x209C',s:'0x209B',
o:'0x2092',p:'0x209A',h:'0x2095',i:'0x1D62',j:'0x2C7C',n:'0x2099',
k:'0x2096',minus:'0x208B',plus:'0x208A',equal:'0x208C',beta:'0x1D66',gamma:'0x1D67',phi:'0x1D69',chi:'0x1D6A'},
  sp: {
l:'0x2E1',u:'0x1D58',o:'0x1D52',w:'0x2B7',p:'0x1D56',i:'0x2071',h:'0x02B0',j:'0x02B2',gamma:'0x02E0',r:'0x02B3',n:'0x207F' ,Yc:'0x0176',Sc:'0x015C',alpha:'0x1D45',beta:'0x1D5D',delta:'0x1D5F',phi:'0x1D60',chi:'0x1D61'},
symbol:{numb:'0x2116',leq:'0x2264',geq:'0x2265',neq:'0x2260',aeq:'0x2248',inf:'0x221E',int:'0x222B',mult:'0x00B7'},
  greek:{ sigma:'0x03C3',delta:'0x03B4',tau:'0x03C4',alpha:'0x03B1',beta:'0x03B2',gamma:'0x03B3',mu:'0x03BC',omega:'0x03C9',epsilon:'0x03B5',ro:'0x03C1',psi:'0x03C8',DELTA:'0x0394',nu:'0x03BD',SIGMA:'0x03A3',chi:'0x03C7'},
 sbd: ['0x2080','0x2081','0x2082','0x2083','0x2084','0x2085','0x2086','0x2087','0x2088','0x2089'],
 spd: ['0x2070','0xB9','0xB2','0xB3','0x2074','0x2075','0x2076','0x77','0x2078','0x2079']
 };
   for(ij in  lt.sb) lt.sb[ij]=String.fromCharCode(lt.sb[ij]);
   for(ij in  lt.sp) lt.sp[ij]=String.fromCharCode(lt.sp[ij]);
   for(ij in  lt.greek) lt.greek[ij]=String.fromCharCode(lt.greek[ij]);
   for(ij in  lt.symbol) lt.symbol[ij]=String.fromCharCode(lt.symbol[ij]);
   for(ij in lt.sbd) {lt.sbd[ij]=String.fromCharCode(lt.sbd[ij]);lt.spd[ij]=String.fromCharCode(lt.spd[ij]);}
}
//******************************************************************************************
function Ltc(col,row) {
   k1=parseInt(col/26);k2=parseInt(col/27);k=parseInt(26*k1);
   if(k2==0)   ltter=String.fromCharCode("A".charCodeAt(0)+(col-1));
   if(k2==1)   ltter=String.fromCharCode("A".charCodeAt(0)+(k1-1))+String.fromCharCode("A".charCodeAt(0)+(col-k-1));
   if(k2>1)   ltter=String.fromCharCode("A".charCodeAt(0)+(k1-1))+String.fromCharCode("A".charCodeAt(0)+(col-k));
   if(ltter=="B@") ltter="AZ";
  return ltter+parseInt(row);
}
//******************************************************************************************
function TextStyle(row,col,otitle,ocolor,talign,fsize,fweight,fstyle,ffamily,txt) {
  obj=elmObject[NumObject(row,col)];
  obj.title=otitle;
  s=obj.style; s.fontStyle=fstyle;s.fontFamily=ffamily;s.color=ocolor;s.fontSize=fsize;
  s.fontWeight=fweight;s.cursor="hand";s.textAlign=talign;
  localStorage[obj.id]=txt;
}
//**************************************************************************