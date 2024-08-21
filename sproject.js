//*********************************************************************************************
function TableData(ind) {
  let cx=[],n=[],variant=[],kx,yd=[],syy=[],pvx=[],sigmax,m0,ep,cp,cv=[],px=[],pv=[],a=[],alpha=[];
switch(ind) {
case 1:
  preg=0.01;
  kx=4;
  cx=[17.0,15.0,13.0,11.0];
  n=[3,10,5,5];
  variant=[
  5.292,5.322,5.681,
  5.538,5.582,5.667,5.732,5.76,5.766,5.839,5.964,6.0,6.09,
  6.082,6.093,6.538,6.587,6.622,
  6.724,7.009,7.111,7.3,-7.344];
  break;
case 2:
  preg=0.1;
  kx=3;
  cx=[18.5,16.5,14.5];
  n=[2,6,5];
  variant=[5.749,5.86,
  5.908,6.079,6.145,6.208,6.589,6.667,
  6.544,6.591,6.752,7,-7.053];
  break;
//###################################Energy################################
case 3:
     sigmax=14.5;
     cp=50.5;
     ep=0.386;
     m0=0.272;
     px=[0.01,0.05,0.1,0.3,0.5,0.7,0.9,0.95,0.99];
     pv=[4.5267857142857135,4.5267857142857135,4.5267857142857135,
             4.5267857142857135,4.5267857142857135,4.5267857142857135,
             4.5267857142857135,4.5267857142857135,4.5267857142857135];
     a=[569.3424278314849,609.1662921734577,631.0603678856719,
          678.190748688179,712.0346046922092,746.8725503950398,
          798.9456476877972,824.6825083888174,874.2872174911809];
     alpha=[2.392265147499767,2.3937375236227156,2.3944913422875898,
                 2.395997471443614,2.3969922651159647,2.397950028627839,
                 2.3992716341320657,2.3998815235976148,2.400985701399871];
     cx=[17.0,15.0,13.0,11.0];
     yd=[5.432,5.794,6.384,7.134];
     syy=[0.216,0.18,0.273,0.297];
     n=[3,10,5,5];
     cv=[1.2,1.32,1.5,1.72,1.8,2.4,2.7,3.,3.6,3.9,4.2,4.44,4.5,4.68,4.8,5.04,5.1,5.28,5.4,5.7,6.];
     pvx=[0.00001,0.005,0.01,0.02,0.025,0.05,0.08,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.925,0.95,0.975,0.99,0.9999];
     break;
//#########################ap#######################
  case 4:
  stat_pvx=4.527;
  stat_AX=712.035;
  stat_ALPHAX=2.397;
  stat_xu=14.5;
  break;
}
//########################################################################
     if(ind==1 || ind==2) {
        TableClear(TableAction1);
        localStorage["AA"]=parseInt(kx);
        nmax=10;//Math.max.apply(Math,n);
        localStorage["AB"]=parseInt(nmax);
        localStorage[str[2][nmax+19]]=parseFloat(preg);
        m2=0;
        for (i=0;i<kx;i++) {
          localStorage[str[i+2][1]]=parseFloat(cx[i]).toFixed(5);
          for (j=0;j<n[i];j++)  localStorage[str[i+2][j+2]]=parseFloat(variant[j+m2]).toFixed(5);
          m2+=n[i];
       }
          TableInit();
          TableAction1();
    }
//#######################################################################
        if(ind==3) {
           TableClear(TableAction2);
            localStorage[str[2][2]]=parseFloat(sigmax);
            localStorage[str[2][3]]=parseFloat(cp);
            localStorage[str[2][4]]=parseFloat(ep);
            localStorage[str[2][5]]=parseFloat(m0);

            for (i=0;i<px.length;i++) {
              localStorage[str[3][i+2]]=parseFloat(px[i]).toFixed(5);
              localStorage[str[4][i+2]]=parseFloat(pv[i]).toFixed(5);
              localStorage[str[5][i+2]]=parseFloat(a[i]).toFixed(5);
              localStorage[str[6][i+2]]=parseFloat(alpha[i]).toFixed(5);
           }

           for (i=0;i<cx.length;i++) {
              localStorage[str[7][i+2]]=parseFloat(cx[i]).toFixed(5);
              localStorage[str[8][i+2]]=parseFloat(yd[i]).toFixed(5);
              localStorage[str[9][i+2]]=parseFloat(syy[i]).toFixed(5);
              localStorage[str[10][i+2]]=parseInt(n[i]).toFixed(5);
           }

           for (i=0;i<cv.length;i++) {
              localStorage[str[11][i+2]]=parseFloat(cv[i]).toFixed(5);
              localStorage[str[12][i+2]]=parseFloat(pvx[i]).toFixed(5);
           }


            TableInit_1();
            TableAction2();
       }
//###########################################################################
       if(ind==4) {
            TableClear(TableAction3);
            localStorage[str[2][2]]=parseFloat(stat_pvx);
            localStorage[str[2][3]]=parseFloat(stat_AX);
            localStorage[str[2][4]]=parseFloat(stat_ALPHAX);
            localStorage[str[2][5]]=parseFloat(stat_xu);
            TableInit_2();
            TableAction3();
        }

        computeAll();
        return;
}
//*******************************************************************************************************
function TableInit() {

     for (i=0;i<colmax;i++)   {
      str[i]=[];
     for (j=0;j<rowmax;j++)  str[i][j]=Ltc(i,j);
   }

 nmax=parseInt(localStorage["AB"]);
elmObject[0].title="Номер_варианта";
for (i=0;i<6;i++) elmObject[i+1].title="Факторы";
for (i=0;i<nmax;i++) localStorage[str[1][i+2]]=(i+1);

TextStyle(1,1,"№","blue","center","12pt","bold","Italic","Times New Roman","№");
TextStyle(1,6,"P","blue","center","12pt","bold","Italic","Times New Roman","P");
TextStyle(1,7,"Нижняя доверительная граница","blue","center","12pt","bold","Italic","Times New Roman","X"+lt.sb.l);
TextStyle(1,8,"Квантиль","blue","center","12pt","bold","Italic","Times New Roman","X"+lt.sb.p);
TextStyle(1,9,"Верхняя доверительная граница","blue","center","12pt","bold","Italic","Times New Roman","X"+lt.sb.u);
TextStyle(1,10,"Нижняя доверительная граница","blue","center","12pt","bold","Italic","Times New Roman","X"+lt.sb.l);
TextStyle(1,11,"Квантиль","blue","center","12pt","bold","Italic","Times New Roman","X"+lt.sb.p);
TextStyle(1,12,"Верхняя доверительная граница","blue","center","12pt","bold","Italic","Times New Roman","X"+lt.sb.u);
TextStyle(1,13,"Нижняя доверительная граница","blue","center","12pt","bold","Italic","Times New Roman","X"+lt.sb.l);
TextStyle(1,14,"Квантиль","blue","center","12pt","bold","Italic","Times New Roman","X"+lt.sb.p);
TextStyle(1,15,"Верхняя доверительная граница","blue","center","12pt","bold","Italic","Times New Roman","X"+lt.sb.u);
TextStyle(1,16,"Нижняя доверительная граница","blue","center","12pt","bold","Italic","Times New Roman","X"+lt.sb.l);
TextStyle(1,17,"Квантиль","blue","center","12pt","bold","Italic","Times New Roman","X"+lt.sb.p);
TextStyle(1,18,"Верхняя доверительная граница","blue","center","12pt","bold","Italic","Times New Roman","X"+lt.sb.u);

TextStyle(nmax+3,1,"Разрушенные объекты","blue","center","12pt","bold","Italic","Times New Roman","r");
TextStyle(nmax+4,1,"Объем испытаний","blue","center","12pt","bold","Italic","Times New Roman","n");
TextStyle(nmax+5,1,"Выборочное среднее","blue","center","12pt","bold","Italic","Times New Roman",lt.greek.mu);
TextStyle(nmax+6,1,"Выборочное ско","blue","center","12pt","bold","Italic","Times New Roman",lt.greek.sigma);
TextStyle(nmax+7,1,"Выборочный коэффициент вариации","blue","center","12pt","bold","Italic","Times New Roman",lt.greek.gamma);
TextStyle(nmax+8,1,"Дисперсия выборочного среднего","blue","center","12pt","bold","Italic","Times New Roman","D{"+lt.greek.mu+"}");
TextStyle(nmax+9,1,"Дисперсия выборочного ско","blue","center","12pt","bold","Italic","Times New Roman","D{"+lt.greek.sigma+"}");
TextStyle(nmax+10,1,"Ковариация оценок","blue","center","12pt","bold","Italic","Times New Roman","D{"+lt.greek.mu+","+lt.greek.sigma+"}");
TextStyle(nmax+11,1,"B","blue","center","12pt","bold","Italic","Times New Roman","B,"+lt.greek.beta);
TextStyle(nmax+12,1,"Расчетная зависимость с.к.о. от среднего","blue","center","12pt","bold","Italic","Times New Roman",lt.greek.sigma+"=B"+lt.symbol.mult+lt.greek.mu+lt.sp.beta);
TextStyle(nmax+13,1,"Статистика критерия нормальности Шапиро-Уилка","blue","center","12pt","bold","Italic","Times New Roman","W"+lt.symbol.geq+"W"+lt.sbd[0]+","+lt.sbd[0]+lt.sbd[5]);
TextStyle(nmax+14,1,"Критическое значение критерия Шапиро-Уилка","blue","center","12pt","bold","Italic","Times New Roman","W"+lt.sbd[0]+","+lt.sbd[0]+lt.sbd[5]);
TextStyle(nmax+15,1,"Критерий Бартлета равенства дисперсий","blue","center","12pt","bold","Italic","Times New Roman",lt.greek.chi+lt.symbol.leq+lt.greek.chi+lt.sbd[0]+","+lt.sbd[0]+lt.sbd[5]);
TextStyle(nmax+16,1,"Критерий Фишера равенства средних","blue","center","12pt","bold","Italic","Times New Roman","F"+lt.symbol.leq+"F"+lt.sbd[0]+","+lt.sbd[0]+lt.sbd[5]);
TextStyle(nmax+17,1,"Критерий Краскелла однородности выборок","blue","center","12pt","bold","Italic","Times New Roman","H"+lt.symbol.leq+"H"+lt.sbd[0]+","+lt.sbd[0]+lt.sbd[5]);
TextStyle(nmax+2,7,"Вероятность","blue","center","12pt","bold","Italic","Times New Roman","P");
TextStyle(nmax+3,7,"Остаточная дисперсия","blue","center","12pt","bold","Italic","Times New Roman","Q");
TextStyle(nmax+4,7,"Показатель степени","blue","center","12pt","bold","Italic","Times New Roman",lt.greek.alpha);
TextStyle(nmax+5,7,"Предел выносливости","blue","center","12pt","bold","Italic","Times New Roman",lt.greek.sigma+lt.sb.minus+lt.sbd[1]);
TextStyle(nmax+6,7,"A=","blue","center","12pt","bold","Italic","Times New Roman","A");
TextStyle(nmax+7,7,"D{a}","blue","center","12pt","bold","Italic","Times New Roman","D{a}");
TextStyle(nmax+8,7,"D{b}=","blue","center","12pt","bold","Italic","Times New Roman","D{b}");
TextStyle(nmax+18,1,"Расчетные значения амплитуд напряжений","blue","center","12pt","bold","Italic","Times New Roman",lt.greek.sigma+lt.sb.a+lt.sb.r);
TextStyle(nmax+19,1,"Вероятность для графика кривой усталости","blue","center","12pt","bold","Italic","Times New Roman","P");
TextStyle(nmax+20,1,"Показатель степени для P","blue","center","12pt","bold","Italic","Times New Roman",lt.greek.alpha);
TextStyle(nmax+21,1,"Предел выносливости для P","blue","center","12pt","bold","Italic","Times New Roman",lt.greek.sigma+lt.sb.minus+lt.sbd[1]);
TextStyle(nmax+22,1,"A для P","blue","center","12pt","bold","Italic","Times New Roman","A");

return;
}

//*******************************************************************************************************
function TableInit_1() {

  for (i=0;i<colmax;i++)   {
      str[i]=[];
     for (j=0;j<rowmax;j++)  str[i][j]=Ltc(i,j);
   } 

  TextStyle(2,1,"Показатель степени","blue","center","12pt","bold","Italic","Times New Roman",lt.greek.sigma+lt.sb.m);
  TextStyle(3,1,"Sigma","blue","center","12pt","bold","Italic","Times New Roman",lt.greek.sigma+lt.sb.s);
  TextStyle(4,1,"Ep","blue","center","12pt","bold","Italic","Times New Roman","e"+lt.sb.p);
  TextStyle(5,1,"M0","blue","center","12pt","bold","Italic","Times New Roman","m"+lt.sb.o);
  TextStyle(1,3,"Вероятности","blue","center","12pt","bold","Italic","Times New Roman","P");
  TextStyle(1,4,"pv","blue","center","12pt","bold","Italic","Times New Roman",lt.greek.sigma+lt.sb.minus+lt.sbd[1]);
  TextStyle(1,5,"a","blue","center","12pt","bold","Italic","Times New Roman","A");
  TextStyle(1,6,"alpha","blue","center","12pt","bold","Italic","Times New Roman",lt.greek.alpha);
  TextStyle(1,7,"cx","blue","center","12pt","bold","Italic","Times New Roman",lt.greek.sigma+lt.sb.a);
  TextStyle(1,8,"lgN","blue","center","12pt","bold","Italic","Times New Roman","lgN");
  TextStyle(1,9,"SlgN","blue","center","12pt","bold","Italic","Times New Roman","SlgN");
  TextStyle(1,10,"n","blue","center","12pt","bold","Italic","Times New Roman","n");
  TextStyle(1,11,"Cv","blue","center","12pt","bold","Italic","Times New Roman","c"+lt.sb.v);
  TextStyle(1,12,"Pvx","blue","center","12pt","bold","Italic","Times New Roman","P");
  TextStyle(1,13,"h[i]","blue","center","12pt","bold","Italic","Times New Roman","h"+lt.sb.i);
  TextStyle(1,14,"cv*smax","blue","center","12pt","bold","Italic","Times New Roman","c"+lt.sb.v+"*"+lt.greek.sigma+lt.sb.m);
  TextStyle(1,15,"alpha","blue","center","12pt","bold","Italic","Times New Roman",lt.greek.alpha+lt.sb.l);
  TextStyle(1,16,"sigma","blue","center","12pt","bold","Italic","Times New Roman",lt.greek.sigma+lt.sb.minus+lt.sbd[1]+lt.sb.l);
  TextStyle(1,17,"a","blue","center","12pt","bold","Italic","Times New Roman","A"+lt.sb.l);
  TextStyle(1,18,"alpha","blue","center","12pt","bold","Italic","Times New Roman","a"+lt.sb.u);
  TextStyle(1,19,"sigma","blue","center","12pt","bold","Italic","Times New Roman",lt.greek.sigma+lt.sb.minus+lt.sbd[1]+lt.sb.u);
  TextStyle(1,20,"a","blue","center","12pt","bold","Italic","Times New Roman","A"+lt.sb.u);

 TextStyle(14,2,"P","blue","center","12pt","bold","Italic","Times New Roman","lgN");
 TextStyle(14,3,"P","blue","center","12pt","bold","Italic","Times New Roman","lgNup");
 TextStyle(14,4,"P","blue","center","12pt","bold","Italic","Times New Roman","lgNlow");
 px=[0.01,0.05,0.1,0.3,0.5,0.7,0.9,0.95,0.99];
 for(i=0;i<9;i++) {
    TextStyle(i+15,1,"P","blue","center","12pt","bold","Italic","Times New Roman",px[i]);
  }

return;
}

//*******************************************************************************************************
function TableInit_2() {

  for (i=0;i<colmax;i++)   {
      str[i]=[];
     for (j=0;j<rowmax;j++)  str[i][j]=Ltc(i,j);
   } 

   TextStyle(2,1,"Предел выносливости","blue","center","12pt","bold","Italic","Times New Roman",lt.greek.sigma+lt.sb.minus+lt.sbd[1]);
   TextStyle(3,1,"A","blue","center","12pt","bold","Italic","Times New Roman","A");
   TextStyle(4,1,"Alpha","blue","center","12pt","bold","Italic","Times New Roman",lt.greek.alpha);
   TextStyle(5,1,"Smax","blue","center","12pt","bold","Italic","Times New Roman",lt.greek.sigma+lt.sb.a);


  return;

  }
