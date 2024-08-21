var igammaepsilon = 0.000000000000001;
var igammabignumber = 4503599627370496.0;     
var igammabignumberinv = 2.22044604925031308085*0.0000000000000001;
var maxrealnumber=Math.pow(10,300);
var minrealnumber=Math.pow(10,-300);
var expm2 = 0.13533528323661269189;
var s2pi = 2.50662827463100050242;
let spi=0.398942280401433;
var logpi = 1.14472988584940017414;
var ls2pi = 0.91893853320467274178;
var machineepsilon=5*Math.pow(10,-16);
//Критические значения критерия Смирнова для уровней значимости alpha=0.01; 0.05; 0.1; и n=3-25;
let uc=[
[1.15,1.15,1.15],[1.49,1.46,1.42], [1.75,1.67,1.60],[1.94,1.82,1.73],[2.10,1.94,1.83],[2.22,2.03,1.91],[2.32,2.11,1.98],[2.41,2.18,2.03],[2.48,2.23,2.09],
[2.55,2.29,2.13],[2.61,2.33,2.17],[2.66,2.37,2.21],[2.70,2.41,2.25],[2.75,2.44,2.28],[2.78,2.48,2.31],[2.80,2.50,2.34],[2.85,2.53,2.36],[2.88,2.56,2.38],
[2.91,2.58,2.41],[2.94,2.60,2.43],[2.96,2.62,2.45],[2.99,2.64,2.47],[3.01,2.66,2.49]
 ];
//Критические значения критерия Шапиро-Уилка для уровней значимости alpha=0.01;0.02;0.05;0.1;0.5 и n=3-50;
//Shapir0.01
wcr001=[0.753,0.687,0.686,0.713,0.730,0.749,0.764,0.781,0.792,0.805, 
0.814,0.825,0.835,0.844,0.851,0.858,0.863,0.868,0.873,0.878, 
0.881,0.884,0.888,0.891,0.894,0.896,0.898,0.900,0.902,0.904, 
0.906,0.908,0.910,0.912,0.914,0.916,0.917,0.919,0.920,0.922, 
0.923,0.924,0.926,0.927,0.928,0.929,0.929,0.930];
//Shapir0.02
wcr02=[0.756,0.707,0.715,0.743,0.760,0.778,0.791,0.806,0.817,0.828, 
0.837,0.846,0.855,0.863,0.869,0.874,0.879,0.884,0.888,0.982, 
0.895,0.898,0.901,0.904,0.906,0.907,0.910,0.912,0.914,0.915, 
0.917,0.919,0.920,0.922,0.924,0.925,0.927,0.928,0.929,0.930, 
0.932,0.933,0.934,0.935,0.936,0.937,0.937,0.938];
//Shapir0.05
wcr005=[0.767,0.748,0.762,0.788,0.803,0.818,0.829,0.842,0.850,0.859, 
0.866,0.874,0.881,0.887,0.892,0.897,0.901,0.905,0.908,0.911, 
0.914,0.916,0.918,0.920,0.923,0.924,0.926,0.927,0.929,0.930, 
0.931,0.933,0.934,0.935,0.936,0.938,0.939,0.940,0.941,0.942, 
0.943,0.944,0.945,0.945,0.946,0.947,0.947,0.947];
//Shapir0.1
wcr01=[0.789,0.792,0.806,0.826,0.838,0.851,0.859,0.869,0.876,0.883, 
0.889,0.895,0.901,0.906,0.910,0.914,0.917,0.920,0.923,0.926, 
0.928,0.930,0.931,0.933,0.935,0.936,0.937,0.939,0.940,0.941, 
0.942,0.943,0.944,0.945,0.946,0.947,0.948,0.949,0.950,0.951, 
0.951,0.952, 0.953,0.953,0.954,0.954,0.955,0.955];
//Shapir0.5
wcr05=[0.959,0.935,0.927,0.927,0.928,0.932,0.935,0.938,0.940,0.943, 
0.945,0.947,0.950,0.952,0.954,0.956,0.957,0.959,0.960,0.961, 
0.962,0.963,0.964,0.965,0.965,0.966,0.966,0.967,0.967,0.967, 
0.968,0.969,0.969,0.970,0.970,0.971,0.971,0.972,0.972,0.972,
0.973,0.973,0.973,0.974,0.974,0.974,0.974,0.974];
//**************************Операции с матрицами****************************
function TransMatrix(A)       //На входе двумерный массив
{
    var m = A.length, n = A[0].length, AT = [];
    for (var i = 0; i < n; i++)
     { AT[i] = [];
       for (var j = 0; j < m; j++) AT[i][j] = A[j][i];
     }
    return AT;
}
//Сложение матриц

function SumMatrix(A,B)       //На входе двумерные массивы одинаковой размерности
{   
    var m = A.length, n = A[0].length, C = [];
    for (var i = 0; i < m; i++)
     { C[i] = [];
       for (var j = 0; j < n; j++) C[i][j] = A[i][j]+B[i][j];
     }
    return C;
}
//Умножение матрицы на число

function multMatrixNumber(a,A)  // a - число, A - матрица (двумерный массив)
{   
    var m = A.length, n = A[0].length, B = [];
    for (var i = 0; i < m; i++)
     { B[i] = [];
       for (var j = 0; j < n; j++) B[i][j] = a*A[i][j];
     }
    return B;
}
//Умножение матриц

function MultiplyMatrix(A,B)
{
    var rowsA = A.length, colsA = A[0].length,
        rowsB = B.length, colsB = B[0].length,
        C = [];
    if (colsA != rowsB) return false;
    for (var i = 0; i < rowsA; i++) C[i] = [];
    for (var k = 0; k < colsB; k++)
     { for (var i = 0; i < rowsA; i++)
        { var t = 0;
          for (var j = 0; j < rowsB; j++) t += A[i][j]*B[j][k];
          C[i][k] = t;
        }
     }
    return C;
}
//Возведение матрицы в степень
//Рекурсивный алгоритм возведения квадратной матрицы A
//A в натуральную степень n.

function MatrixPow(n,A)
{ 
    if (n == 1) return A;     // Функцию MultiplyMatrix см. выше
    else return MultiplyMatrix( A, MatrixPow(n-1,A) );
}
//Определитель матрицы

function Determinant(A)   // Используется алгоритм Барейса, сложность O(n^3)
{
    var N = A.length, B = [], denom = 1, exchanges = 0;
    for (var i = 0; i < N; ++i)
     { B[i] = [];
       for (var j = 0; j < N; ++j) B[i][j] = A[i][j];
     }
    for (var i = 0; i < N-1; ++i)
     { var maxN = i, maxValue = Math.abs(B[i][i]);
       for (var j = i+1; j < N; ++j)
        { var value = Math.abs(B[j][i]);
          if (value > maxValue){ maxN = j; maxValue = value; }
        }
       if (maxN > i)
        { var temp = B[i]; B[i] = B[maxN]; B[maxN] = temp;
          ++exchanges;
        }
       else { if (maxValue == 0) return maxValue; }
       var value1 = B[i][i];
       for (var j = i+1; j < N; ++j)
        { var value2 = B[j][i];
          B[j][i] = 0;
          for (var k = i+1; k < N; ++k) B[j][k] = (B[j][k]*value1-B[i][k]*value2)/denom;
        }
       denom = value1;
     }
    if (exchanges%2) return -B[N-1][N-1];
    else return B[N-1][N-1];
}
//Ранг матрицы

function MatrixRank(A)
{
    var m = A.length, n = A[0].length, k = (m < n ? m : n), r = 1, rank = 0;
    while (r <= k)
     { var B = [];
       for (var i = 0; i < r; i++) B[i] = [];
       for (var a = 0; a < m-r+1; a++)
        { for (var b = 0; b < n-r+1; b++)
           { for (var c = 0; c < r; c++)
              { for (var d = 0; d < r; d++) B[c][d] = A[a+c][b+d]; }
             if (Determinant(B) != 0) rank = r;
           }       // Функцию Determinant см. выше
        }
       r++;
     }
    return rank;
}
//Союзная матрица
//Союзной к матрице A называют матрицу adjA, которая получается из матрицы A, если все ее элементы заменить соответствующими алгебраическими дополнениями Ai,j
//Ai,j и к полученной матрице применить операцию транспонирования.

function AdjugateMatrix(A)   // A - двумерный квадратный массив
{                                        
    var N = A.length, adjA = [];
    for (var i = 0; i < N; i++)
     { adjA[i] = [];
       for (var j = 0; j < N; j++)
        { var B = [], sign = ((i+j)%2==0) ? 1 : -1;
          for (var m = 0; m < j; m++)
           { B[m] = [];
             for (var n = 0; n < i; n++)   B[m][n] = A[m][n];
             for (var n = i+1; n < N; n++) B[m][n-1] = A[m][n];
           }
          for (var m = j+1; m < N; m++)
           { B[m-1] = [];
             for (var n = 0; n < i; n++)   B[m-1][n] = A[m][n];
             for (var n = i+1; n < N; n++) B[m-1][n-1] = A[m][n];
           }
          adjA[i][j] = sign*Determinant(B);   // Функцию Determinant см. выше
        }
     }
    return adjA;
}
//Обратная матрица
//Вычисление обратной матрицы с помощью союзной матрицы (из алгебраических дополнений, см. выше).

function InverseMatrix(A)   // A - двумерный квадратный массив
{   
    var det = Determinant(A);                // Функцию Determinant см. выше
    if (det == 0) return false;
    var N = A.length, A = AdjugateMatrix(A); // Функцию AdjugateMatrix см. выше
    for (var i = 0; i < N; i++)
     { for (var j = 0; j < N; j++) A[i][j] /= det; }
    return A;
}

//*********************Метод наименьших квадратов******************
//с весами V
function mleastsquare(x,y,kmls,V,bmls,dmls) {
   var d=Determinant(V);  
    if (d!=0) {
    V=InverseMatrix(V);
     var dmlsx=[],bmlsx=[];
    for (i=0;i<kmls;i++) {
    dmlsx[i]=[];  
    for (j=0;j<kmls;j++) dmlsx[i][j]=0;
    }
    dmlsx=InverseMatrix(MultiplyMatrix(MultiplyMatrix(TransMatrix(x),V),x));
    bmlsx=MultiplyMatrix(MultiplyMatrix(MultiplyMatrix(dmlsx,TransMatrix(x)),V),y);
    }
   for (i=0;i<kmls;i++) {
      bmls[i]=bmlsx[i];
      for (j=0;j<kmls;j++) dmls[i][j]=dmlsx[i][j];
   }
 }
//********************************************************************
//без весов
 function mleastsquare_line(x,y,kmls,nmls,bmls,dmls,yp,dy) {
    dmlsx=InverseMatrix(MultiplyMatrix(TransMatrix(x),x));
    bmlsx=MultiplyMatrix(MultiplyMatrix(dmlsx,TransMatrix(x)),y);
 for (i=0;i<kmls;i++)  {
      bmls[i]=bmlsx[i];
      for (j=0;j<kmls;j++) dmls[i][j]=dmlsx[i][j];
  }
//**********в матричном виде***********************************************
/*
   var q=[],yp=[],dy=[];
for(i=0;i<nmls;i++) {
     yp[i]=[];q[i]=[],dy[i]=[];
}
   yp=MultiplyMatrix(x,bmls); //расчетное y
   q=MultiplyMatrix(TransMatrix(SumMatrix(yp,multMatrixNumber(-1,y))),SumMatrix(yp,multMatrixNumber(-1,y))); //остаточная дисперсия
    dy=MultiplyMatrix(MultiplyMatrix(x,dmls),TransMatrix(x)); //дисперсия y диагональные элементы dy
    Q=q[0][0];
    //yp[nmls][0];dy[nmls][nmls];
*/
//************в обычном****************************************************
    Q=0;
      for(i=0;i<nmls;i++) {
        s1=0;s2=0;
        for(j=0;j<kmls;j++)   s1+=bmlsx[j]*x[i][j];
         Q+=(s1-y[i][0])*(s1-y[i][0]);
         yp[i]=s1; 
         for(j=0;j<kmls;j++) {
          for(k=0;k<kmls;k++) s2+=dmlsx[j][k]*x[i][j]*x[i][k];
         }
          dy[i]=s2;
        }
        Q=Q/(nmls-kmls);
        return Q;
}
//************Операции с полиномами************************************************

//****************Умножение полиномов***************************************************

function pmpy(x,idimx,y,idimy,z) {
idimz=idimx+idimy-1;
  for(i=0;i<=idimz;i++) z[i]=0;
for(i=1;i<=idimx;i++) {
   for(j=1;j<=idimy;j++) z[i+j-1]=z[i+j-1]+x[i]*y[j];
}
return idimz;
}
//*****************************************************************************
function pnorm(x,idimx,eps) {
  while(1) {
  if(idimx<=0) return;
  z=Math.abs(x[idimx])-eps;
  if(z>0) return;
  idimx--;
  }
}
//********Деление полиномов*********************************************************
function pdiv(x,idimx,y,idimy,z) {
tol=0;
flag="start";
while(1) {
  switch(flag) {
    case "start":
 pnorm(y,idimy,tol);
 if(idimy<=0) {flag="m50";continue;}
 if(idimy>0) {flag="m10";continue;}
 case "m10":
 idimz=idimx-idimy+1;
 if(idimz==0) {flag="m30";continue;}
 if(idimz<0) {flag="m20";continue;}
 if(idimz>0) {flag="m60";continue;}
 case "m20":
 idimz=0;
 case "m30":
 ier=0;
 return idimz;
 case "m50":
 ier=1;
 return idimz;
 case "m60":
 idimx=idimy-1;
 i=idimz;
 case "m70":
 ii=i+idimx;
 z[i]=x[ii]/y[idimy];
 for(k=1;k<=idimx;k++) {
   j=k-1+i;
   x[j]-=z[i]*y[k];
 }
 i--;
 if(i<=0) {flag="m90"; continue;}
 if(i>0) {flag="m70"; continue;}
 case "m90":
 pnorm(x,idimx,tol);
 {flag="m30"; continue;}
 return idimz;
}}
}

//**********Линейная сплайн-интерполяция***************

function linearspline(x,y) {
    n=parseInt(x.length);
    for (i=0;i<n-1;i++) cx[i]=(y[i+1]-y[i])/(x[i+1]-x[i]);
}
function splineinterpolation(xb) {
   if (xb<x[0] || xb>x[n-1]) return 0;
   for (j=0;j<n;j++) {     
        if (xb<=x[j]) break;
   }
    m=j-1;
    xb0=xb-x[m];
    return  y[m]+xb0*cx[m];
   }
//*************Статистические программы********************************

//***Сортировка от 1***************************************************
function sort_1(x,n) {
 var xx; 
for (i=1;i<=n;i++) {
 for (j=i+1;j<=n;j++) {
   if (x[i]>x[j]) {
    xx=x[i];x[i]=x[j];x[j]=xx;
   }
 }
}
}
//****************Сортировка от 0*********************************************
function sort_0(x,n) {
 var xx; 
for (i=0;i<n;i++) {
 for (j=i+1;j<n;j++) {
   if (x[i]>x[j]) {
    xx=x[i];x[i]=x[j];x[j]=xx;
   }
 }
}
}
//********Сортировка связанных массивов*************************************  
function twovar(x,a,n) {
for (i=1;i<=n-1;i++) {
 for (j=i+1;j<=n;j++) {
  if (x[i]>x[j]) {
   xx=x[i];x[i]=x[j];x[j]=xx;
   xt=a[i];a[i]=a[j];a[j]=xt;
  }
 }
 }
}
//*****************Вычисление среднего и ско**********************************
function p24a05(x,n,bmls) {
  s1=0;s2=0;
  for(i=0;i<n;i++) {
    s1=s1+x[i];s2=s2+x[i]*x[i];
  }
  cp=s1/n;
  if(n<=1) {
    cko=0;return;
  }
   cko=Math.sqrt((s2-cp*cp*n)/(n-1));
 bmls[0]=cp;
 bmls[1]=cko;
 } 

//**************Метод наименьших квадратов для порядковых статистик***********************
function mlses(ts,n,l,x,dmls) {
var E= [],V= [],bmls=[],vorder=[];
 kmls=2;
 for (i=0;i<kmls;i++) {
 dmls[i]=[];  
for (j=0;j<kmls;j++) dmls[i][j]=0;
  }

for (i=0;i<l;i++) {
  V[i]= []; 
   for (j=0;j<l;j++) {
    V[i][j]=0;}}

for (i=0;i<l;i++) {
   for (j=i;j<l;j++) {
   if (ts=="normal") vorder=ordern(n, i+1, j+1);
   if (ts=="weibull") vorder=orderw(n, i+1, j+1);
    E[i]=vorder[0];V[j][i]=vorder[1];V[i][j]=vorder[1];
   }
   }

V=InverseMatrix(V);

s1=0;s2=0;s3=0;s4=0;s5=0;s6=0;
for (i=0;i<l;i++) {
ei=E[i];ai=0;ai1=0;bi=0;
for (j=0;j<l;j++) {
Vij=V[j][i];
ai=ai+E[j]*Vij;bi=bi+Vij*x[j];ai1=ai1+Vij;
}
s1=s1+ai;s2=s2+bi;s3=s3+ei*bi;s4=s4+ei*ai;s5=s5+ai1;s6=s6+ai1*ei;
}
d=s5*s4-s1*s1;
bmls[0]=-(s3*s1-s4*s2)/d;
bmls[1]=(s3*s5-s6*s2)/d;
dmls[0][0]=s4/d;dmls[1][1]=s5/d;
dmls[0][1]=-s1/d;dmls[1][0]=-s1/d;
return bmls;
}
//***********Вейбулловкие порядковые статистики**************************************

function orderw(n,r,s) {
pr=(r/(1+n)); qr=1-pr; ps=(s/(1+n)); qs=1-ps;
xr=Math.log(Math.log(1/(1-pr)));
xs=Math.log(Math.log(1/(1-ps)));
xr1=1/(Math.log(1/(1-pr))*(1-pr));
xr2=xr1*(1/(1-pr)-xr1);
xr3=xr2*xr2/xr1+xr1*(1/Math.pow(1-pr,2)-xr2);
xr4=(3*xr1*xr2*xr3-2*Math.pow(xr2,3))/Math.pow(xr1,2)+xr1*(2/Math.pow(1-pr,3)-xr3);
xr55=(-12*xr1*Math.pow(xr2,2)*xr3+3*Math.pow(xr1,2)*Math.pow(xr3,2)+4*Math.pow(xr1,2)*xr2*xr4+6*Math.pow(xr2,4));
xr5=xr55/Math.pow(xr1,3)+xr1*(6/Math.pow(1-pr,4)-xr4);
a1=-12*Math.pow(xr2,3)*xr3-12*xr1*(2*xr2*Math.pow(xr3,2)+Math.pow(xr2,2)*xr4);
b1=6*xr1*xr2*Math.pow(xr3,2)+6*Math.pow(xr1,2)*xr3*xr4;
c1=8*xr1*Math.pow(xr2,2)*xr4+4*Math.pow(xr1,2)*(xr3*xr4+xr2*xr5);
d1=24*Math.pow(xr2,3)*xr3;
xr6=(Math.pow(xr1,3)*(a1+b1+c1+d1)-3*Math.pow(xr1,2)*xr2*xr55)/Math.pow(xr1,6)+xr2*(6/Math.pow(1-pr,4)-xr4)+xr1*(24/Math.pow(1-pr,5)-xr5);
xs1=1/(Math.log(1/(1-ps))*(1-ps));
xs2=xs1*(1/(1-ps)-xs1);
xs3=Math.pow(xs2,2)/xs1+xs1*(1/Math.pow(1-ps,2)-xs2);
xs4=(3*xs1*xs2*xs3-2*Math.pow(xs2,3))/Math.pow(xs1,2)+xs1*(2/Math.pow(1-ps,3)-xs3);
xs5=(-12*xs1*Math.pow(xs2,2)*xs3+3*Math.pow(xs1,2)*Math.pow(xs3,2)+4*Math.pow(xs1,2)*xs2*xs4+6*Math.pow(xs2,4)) /Math.pow(xs1,3)+xs1*(6/Math.pow(1-ps,4)-xs4);
er=xr+pr*qr*xr2/(2*(n + 2))+pr*qr*((qr-pr)*xr3/3+pr*qr*xr4/8)/Math.pow(n+2,2)+pr*qr*(-(qr-pr)*xr3/3+(Math.pow(qr-pr,2)-pr*qr)*xr4/4+qr*pr*(qr-pr)*xr5/6+Math.pow(qr*pr,2)*xr6/48)/Math.pow(n + 2,3);
z1=(qr-pr)*xr2*xs1+(qs-ps)*xr1*xs2+pr*qr*xr3*xs1/2+ps*qs*xr1*xs3/2+pr*qs*xr2*xs2/2;
z1=z1*pr*qs/Math.pow(n+2,2);
z2=-(qr-pr)*xr2*xs1-(qs-ps)*xr1*xs2+(Math.pow(qr-pr,2)-pr*qr)*xr3*xs1;
z3=(Math.pow(qs-ps,2)-ps*qs)*xr1*xs3+(1.5*(qr-pr)*(qs-ps)+0.5*ps*qr-2*pr*qs)*xr2*xs2;
z4=(5/6)*pr*qr*(qr-pr)*xr4*xs1+(5/6)*ps*qs*(qs-ps)*xr1*xs4+(pr*qs*(qr-pr)+0.5*pr*qr*(qs-ps))*xr3*xs2;
z5=(pr*qs*(qs-ps)+0.5*ps*qs*(qr-pr))*xr2*xs3+(1/8)*Math.pow(pr*qr,2)*xr5*xs1+(1/8)*Math.pow(ps*qs,2)*xr1*xs5;
z6=0.25*Math.pow(pr,2)*qr*qs*xr4*xs2+0.25*pr*ps*Math.pow(qs,2)*xr2*xs4+(2*Math.pow(pr*qs,2)+3*pr*qr*ps*qs)*xr3*xs3/12;
z7=z2+z3+z4+z5+z6;
crs=z1+pr*qs*z7/Math.pow(n+2,3)+pr*qs*xr1*xs1/(n+2);
var vorder=[];
vorder[0]=er;vorder[1]=crs;
return vorder;
}

//********Нормальные порядковые статистики****************************************

function ordern(n,r,s) {
 pr=r/(1+n);qr=1-pr; ps=s/(1+n);qs=1-ps;p=1;k=n; 
xr=invnormaldistribution(pr);
xs=invnormaldistribution(ps);
dr=Math.exp(-xr*xr/2.)/s2pi;
ds=Math.exp(-xs*xs/2.)/s2pi;
xr1=p/dr;
xr2=xr*Math.pow(p/dr,2);
xr3=(2*xr*xr+1)*Math.pow(p/dr,3);
xr4=(6*xr*xr*xr+7*xr)*Math.pow(p/dr,4);
xr5=(24*Math.pow(xr,4)+46*xr*xr+7)*Math.pow(p/dr,5);
xr6=(120*Math.pow(xr,5)+326*Math.pow(xr,3)+127*xr)*Math.pow(p/dr,6);
xs1=p/ds;
xs2=xs*Math.pow(p/ds,2);
xs3=(2*xs*xs + 1)*Math.pow(p/ds,3);
xs4=(6*Math.pow(xs,3)+7*xs)*Math.pow(p/ds,4);
xs5=(24*Math.pow(xs,4)+46*Math.pow(xs,2)+7)*Math.pow(p/ds,5);
xs6=(120*Math.pow(xs,5)+326*Math.pow(xs,3)+127*xs)*Math.pow(p/ds,6);
er=xr+pr*qr*xr2/(2*(k+2))+pr*qr*((qr-pr)*xr3/3+pr*qr*xr4/8)/Math.pow(k+2,2)+pr*qr*(-(qr-pr)*xr3/3+(Math.pow(qr-pr,2)-pr*qr)*xr4/4+qr*pr*(qr-pr)*xr5/6+Math.pow(qr*pr,2)*xr6/48)/Math.pow(k+2,3);

z1=(qr-pr)*xr2*xs1+(qs-ps)*xr1*xs2+pr*qr*xr3*xs1/2+ps*qs*xr1*xs3/2+pr*qs*xr2*xs2/2;
z1=z1*pr*qs/Math.pow(k+2,2);
z2=-(qr-pr)*xr2*xs1-(qs-ps)*xr1*xs2+(Math.pow(qr-pr,2)-pr*qr)*xr3*xs1;
z3=(Math.pow(qs-ps,2)-ps*qs)*xr1*xs3+(1.5*(qr-pr)*(qs-ps)+0.5*ps*qr-2*pr*qs)*xr2*xs2;
z4=(5/6)*pr*qr*(qr-pr)*xr4*xs1+(5/6)*ps*qs*(qs-ps)*xr1*xs4+(pr*qs*(qr-pr)+0.5*pr*qr*(qs-ps))*xr3*xs2;
z5=(pr*qs*(qs-ps)+0.5*ps*qs*(qr-pr))*xr2*xs3 + (1 /8)*Math.pow(pr*qr,2)*xr5*xs1+(1/8)*Math.pow(ps*qs,2)*xr1*xs5;
z6=0.25*Math.pow(pr,2)*qr*qs*xr4*xs2+0.25*pr*ps*Math.pow(qs,2)*xr2*xs4+(2*Math.pow(pr*qs,2)+3*pr*qr*ps*qs)*xr3*xs3/12;
z7=z2+z3+z4+z5+z6;
crs=z1+pr*qs*z7/Math.pow(k+2,3)+pr*qs*xr1*xs1/(k+2);
var vorder=[];
vorder[0]=er;vorder[1]=crs;
return vorder;
}

//*****************Приближенные доверительные граница для квантиля**************************
function lmtaprn(beta,n,t1,t2,t12,zp,xpborder){
zb=invnormaldistribution(beta);
//zp=invnormaldistribution(p);
f1x=t2/(n-1);
f2x=2*t12/Math.sqrt(n);
f4x=1-f1x/2;
e3x=f4x*f4x-zb*zb*f1x;
f3x=zp*Math.sqrt(n);
e11=f4x*f3x+zb*zb*f2x/2;
e2x=f3x*f3x-zb*zb*t1;
e44=Math.sqrt(Math.abs(e11*e11-e2x*e3x));
xpborder[0]=(e11+e44)/e3x;
xpborder[1]=(e11-e44)/e3x;
}
//*************Моделирование случайного цензурирования***************************
function random_cens(n,xc,x) {

var r=0;
xmax=0;
var k=0;
var z=0;
var yc=[];

for (j=0;j<n;j++) {
z=xc[j];
if (z<0) {
 yc[r]=-z;r++;
 }
 else {
 x[k]=z;k++;
 if (z>xmax) xmax=z;
 }
}

 j=0;m=0;
 while(j<r) {
  if (yc[j]<xmax) {
     jj=Math.round((k-1)*Math.random());
     if (x[jj]<yc[j]) {
      j--;
     }
    else {
    yc[j] = x[jj];
    }
  }
  else {
  m++;
  }
 j++; 
 }

  r=r-m;
  for (j=0;j<r;j++) x[j+k]=yc[j];
 return m;
}

//************Нецентральное распределение Стьюдента*******************

  function dovq(beta,delta,f) {
  ki=501;
  tp=invnontap(beta,f,delta);
  tl=tp/1.2;tu=tp*1.2;step=(tu-tl)/(ki-1);
   for (j=0;j<ki;j++) {
     xx=tl+j*step;
     yy=sf54r(xx,delta,f);
     if (xx>=0 && yy>=beta) break;
     if (xx<0 && yy<=beta) break;
     xx0=xx;yy0=yy;
     }
    cx=(yy-yy0)/step
    return xx0+(beta-yy0)/cx;
}
//***********************************************************************
  function dovq1(beta,delta,f) {
       ki = 5001;
       tp=invnontap(beta,f,delta);
       tl=tp/1.4;tu=tp*1.4;step=(tu-tl)/(ki-1);
     for (j=1;j<=ki;j++) {
       x=tl+(j-1)*step;
       bet=sf54r(x,delta,f);
       if (x>=0 && bet>=beta) return x;
       if (x<0 && bet<=beta) return x;
     } 
  return x;
}
//************************************************************************
function invnontap(beta,f,d) {
 zb=invnormaldistribution(beta);
 f4x=1-1/(4*f);
 z=(f4x*d+zb*Math.sqrt(f4x*f4x-zb*zb/(2.*f)+d*d/(2.*f)))/(f4x*f4x-zb*zb/(2.*f));
 return z;
}
//*********************************************************************
 function sf35r(x) {
   t = Math.abs(x);s=0;
   if (t<6.5) {
s1 =Math.exp(-t*t)*((((((0.56419*t+6.802899)*t+38.71143)*t+131.1266)*t+278.5978)*t+355.969)*t+224.1828);
s2=(((((((t+12.05784)*t+69.11384)*t+238.4503)*t+527.5538)*t+741.5214)*t+608.9322)*t+224.1828);
      s=s1/s2;
   }
      if (x>=0) return s;
      if (x<0) return (2-s);
 }
//***********************************************************************
function sf49r(x) {
 z=0.5*sf35r(-0.7071067*x);
 return z;
}
//*********************************************************************
function sf53r(y,z,eps) {
  sys076 = 88.72283;sys017 = 0.1591549;expov = 88.72283;c = 0.1591549;ep1 = eps;
 if (eps==0) ep1=0.000001;
 t=0;b=Math.abs(y);a=Math.abs(z);
 
 if (a==0) return t;
 ta =Math.atan(a);
 
 if (a*b>4) {
 t = sf49r(b);
 t =c*(ta+Math.atan(1/a))-0.5*(t-0.5);
 if (z<0) t=-t;
 return t;
 }

 hsqb=0.5*b*b;
  if (hsqb>expov) return t;
   bexp=Math.exp(-hsqb);asq=a*a;a4=asq*asq;b4=hsqb*hsqb;a4b4=a4*b4;ahsqb=a*hsqb;ab4=a*b4*0.5;
    f=1;sum=0;g=3;

//m2:
 while(1<0) {
  g1=g;ber=0;ter=ab4;

//m3:
while (1<0) {
 ber=ber+ter;
 if (ter<=ber*ep1) break; //break m3;
 ter=ter*hsqb/g1;g1=g1+1;
 }

  d1=(ber+ahsqb)/f;d2=ber*asq/(f+2);d=d1-d2;sum=sum+d;t=ta-sum*bexp;aeps=ep1*t;
  ahsqb=ahsqb*a4b4/((g-1)*g);ab4=ab4*a4b4/((g+1)*g);f=f+4;g=g+2;
   if (d2*bexp<aeps) break; //break m2;
 }
   t=t*c;
 
  if (z<0) t=-t;
  return t;
} 
  
//*************************************************************************  
function sf54r(x,d,idf) {

sys059 = maxrealnumber;
c = 0.1591549;a2 = 0.1591549;a3 = 2.506628;a1 = 0.3989423;sys089 = 0.00001;sys029 = 0.7071068;

if (idf<=0) {ierr=65;return sys059;}
df=idf;tval=x;
i1=idf-2*Math.round(idf/2);  //??????????????? was int
a=tval/Math.sqrt(df);b=df/(df+tval*tval);
sb =Math.sqrt(b);da=d*a;dsb=d*sb;dasb=a*dsb;
p1 = sf49r(dasb);
f2=a*sb*Math.exp(-0.5*dsb*dsb)*p1*a1;
f1=b*(da*f2+a*a2*Math.exp(-0.5*d*d));
sum = 0;
if (idf!=1) {
 if (i1>0) {
  sum = f1;
 } 
else {
  sum = f2;
}
if (idf>=4) {
idfm2 = idf - 2;az = 1;fz = 2;
 for (l=2;l<=idfm2;l+=2) {
     fkm1 = fz - 1; f2 = b * (da * az * f1 + f2) * fkm1 / fz;
      az = 1/(az * fkm1);f1 = b * (da * az * f2 + f1) * fz / (fz + 1);
      if (i1<=0) {
       sum = sum + f2;
      }
      else {
        sum = sum + f1;
      }
      az = 1/(az*fz);fz = fz + 2;
}
}
}
       if (i1>0) {
         p1=0.5 * sf35r(sys029*dsb);
         p=sf53r(dsb,a,sys089);
         z=p1+2*(p+sum);
        }
        else {
         p1=0.5*sf35r(sys029*d);
         z=p1+sum*a3;
       }
       
       if (z<0) return 0;
       return z;
} 
//**************************************************************************
//********Критерии проверки статистических гипотез****************************
//#############################################################################
function omega1(kx,m2,xz,w,p) {
    var i, j, km, i1, j1, x, k, mx,kck,m3;
    var c=[];
    var pr=[];
    var s,nnew, nn;
    var a=[],anit=[],h=[],m=[];				
    
    n=0;
    for(i=0;i<kx;i++) n+=m2[i];
    for (i=1;i<=n;i++) a[i] = i;

    km = 0;
   for (i = 1; i <= kx; i++) {
      for (j = 1; j <= m2[i - 1]; j++) anit[j + km] = i;
    m[i]=m2[i-1];
    km += m2[i-1];
   }
//########################################################################
    kck = 0;
    for (i = 1; i <= n; i++) {
        c[i] = 1;pr[i]=true;
    }
    c[n] = 0;
    h[kck] = kruskalstatistic(a,kx,n,m);
    kck++;
    i = 1;

    while (i < n) {
        i = 1; x = 0;
        while (c[i]==(n - i + 1)) {
            pr[i]=!pr[i];
            c[i]=1;
            if (pr[i]) x++;
            i++;
        }
        if (i<n) {
            if (pr[i]) {
                k=c[i]+x;
            }
            else {
                k=n-i+1-c[i]+x;
            }
            mx=a[k];a[k]=a[k+1];a[k+1]=mx;
            m3=0;
            for (i1 = 1; i1 <= n; i1++) {
                for (j1 = i1; j1 <= n; j1++) {
                    if ((a[i1] > a[j1]) && (anit[a[i1]] == anit[a[j1]])) {m3=1;break;}
                }
                  if(m3==1) break;
            }
            if(m3!=1) {h[kck] = kruskalstatistic(a,kx,n,m);kck++;}
            c[i]+= 1;
        }
    }
//#######################################################################
   h.sort(function(a,b){return a-b});
   nnew=0;nn=0; 
     
  for (i=0;i<kck;i++) {
  if(parseFloat(h[i]).toFixed(znaki)==parseFloat(h[i+1]).toFixed(znaki)) {
    nnew++;
   }
   else {
   w[nn]=nnew+1;xz[nn]=h[i];nnew=0;nn++;
   }
   }
     s=0;kxt=0;
    for (i=0;i<nn;i++) {
      s+=parseFloat(w[i]/kck);
      p[i]=s;
      //if (p[i]<alfa/2) kxt=i;
    }

    return kck;
}
//##########################################################################
//********Перебор всех перестановок элементов множества**********************
  function omega(crit,kx,m,a,anit,n,mx,h) {
     var i,j,k,s,ss,test,r,r1,r2,msmal,km;
    test = 0;
    if (mx==1) {
     for (i=1;i<=n;i++) {
        for (j=i;j<=n;j++) {
           if (a[i]>a[j] && anit[a[i]]==anit[a[j]])  {test=1;break;}
      }
        if (test==1) break;
     }
  
//**********************************************************************
    if (test==0) {
     switch(crit) {
     case "Leman":
      h[kk]=lemanstatistic(a,m[1],m[2]);kk++;
      break;
//***************************************************************************
   case "Wilcoxon":
   case "STukey":
    h[kk]=wilcoxonstatistic(a,m[1],m[2]);kk++;
    break;
 //******************************************************************************
    case "Kruskal":
        h[kk]=kruskalstatistic(a,kx,n,m); kk++;
        break;
//******************************************************************************
  case "Mood":
        h[kk]=moodstatistic(a,m[1],m[2]); kk++;
        break;
     }
     }
//***********************************************************************
  }
     else {
      for (i=1;i<=mx;i++) {
       ss=a[1];
       for (j=1;j<=mx-1;j++) a[j]=a[j+1];
       a[mx]=ss;
       omega(crit,kx,m,a,anit,n,mx-1,h);
      }
     }
     return kk;
    } 

//************cns2=n!/n1!*n2!...nkx!;mm[kx]=n1+n2+...+nkx******************
function cns2(kx,n,mm) {
  var s1,s2,s3,i,j,w2;
  s1=0;
  for (i=1;i<=n;i++) s1=s1+Math.log(i);
   s3=0;
  for (j=1;j<=kx;j++) {
   s2=0;
    for (i=1;i<=mm[j];i++) s2=s2+Math.log(i);
   s3=s3+s2;
  }

return Math.exp(s1-s3);
}

//***********cnm=n!/m!*(n-m)!********************************
function cnm(n,m) {
  var s1,s2,i;
  s1=0; s2=0;
  for (i=m+1;i<=n;i++) s1=s1+Math.log(i);
    for (i=1;i<=n-m;i++)   s2=s2+Math.log(i);
return Math.exp(s1-s2);
}
//**************************************************
function moodstatistic(awa,m1,m2) {
      msmal=Math.min(m1,m2);
      n=m1+m2;
         r=0;
         for (i=1;i<=n;i++) {
          for (j=1;j<=msmal;j++) if (awa[i]==j) r=r+(i-(n+1)/2)**2;
         }
         return r;
}
//**************************************************
function lemanstatistic(ala,m1,m2) {
 //var r1,r2,i,j,k;
   r1=0;r2=0;n=m1+m2;
     for (k=1;k<=m1;k++) {
       for (i=1;i<=n;i++) if (ala[i]==k) r1=r1+(i-k)*(i-k);
      }
     for (k=1;k<=m2;k++) {
       for (i=1;i<=n;i++) if (ala[i]==k+m1) r2=r2+(i-k)*(i-k);
      }
      zleman=(r1*m1+r2*m2+m1*m2/6)/(m1*m1*m2*m2)-2/3;
      return zleman;
}
//**************************************************
 function kruskalstatistic(aa,kx,n,m) {
   //var s,km,i,j,k,r;
       s=0;km=0;
        for (j=1;j<=kx;j++) {
          r=0;
          for (i=1;i<=n;i++) {
            for (k=1;k<=m[j];k++) if (aa[i]==k+km) r=r+i;
          }
            s=s+r*r/m[j];km=km+m[j];
        }
        return 12.00*s/(n*(n+1))-3*(n+1);
}
//**************************************************
 function wilcoxonstatistic(awa,m1,m2) {
     msmal=Math.min(m1,m2);
     n=m1+m2;
         r=0;
         for (i=1;i<=n;i++) {
          for (j=1;j<=msmal;j++) if (awa[i]==j) r=r+i;
         }
         return r;
}
//**************************************************
 function seriesstatistic(n,asa,col) {
  ks=1;ksr=0;
  for (i=1;i<=n;i++) {
   if (asa[i]==asa[i+1]) { 
     ks++;
   }
   else {
    ksr++;
    col[ksr]=ks;ks=1;
   } 
  }
  return ksr;
}
//*****************************************************************************
function mann(n,l,x,a) {
let e=[],i,s1,s2,k,k1,vorder=[],wmann;

 for(i=0;i<n;i++) {
    vorder=ordern(n,i+1,i+1);
    e[i]=vorder[0];
}

  k=l-1;k1=parseInt(l/2)+1;

  s1=0;s2=0;
 for(i=0;i<k;i++) {
   a[i]=e[i+1]-e[i];
   z=x[i+1]-x[i];
   s1=s1+z/a[i];
   if((i+1)>=k1) s2=s2+z/a[i];
 }
  wmann=s2/s1;
  return wmann;
}
//*****************************************************************************
function manw(n,l,x,a) {
let e=[],i,s1,s2,k,k1,vorder=[],wmanw;

 for(i=0;i<n;i++) {
    vorder=orderw(n,i+1,i+1);
    e[i]=vorder[0];
}

  k=l-1;k1=parseInt(l/2)+1;
  s1=0;s2=0;
 for(i=0;i<k;i++) {
   a[i]=e[i+1]-e[i];
   z=Math.log(x[i+1])-Math.log(x[i]);
   s1=s1+z/a[i];
   if((i+1)>=k1) s2=s2+z/a[i];
 }
  wmanw=s2/s1;
  return wmanw;
}
//*****************************************************************************
function shapir(n,x) {
let e=[],v=[],z=[],i,j,ck,s1,s2,c,vorder=[],wshapir,a=[];
 for (i=0;i<n;i++) {
 v[i]=[];e[i]=0;  
for (j=0;j<n;j++) v[i][j]=0;
  }
for(i=0;i<n;i++) {
  for(j=i;j<n;j++) {
    vorder=ordern(n,i+1,j+1);
    e[i]=vorder[0];v[i][j]=vorder[1];v[j][i]=v[i][j];
  }
}
v=InverseMatrix(v); 
for(i=0;i<n;i++) {
  s1=0;
  for(j=0;j<n;j++) s1+=v[j][i]*e[j];
  a[i]=s1;
}
s3=0;
for(i=0;i<n;i++) s3+=a[i]*a[i];
c=Math.sqrt(s3);
s1=0;s2=0;
for(i=0;i<n;i++) {
 s1+=x[i]*x[i];
 s2+=x[i];
}
ck=s1-(s2*s2/n);
s2=0;
for(i=0;i<n;i++) {
  a[i]=a[i]/c;s2+=a[i]*x[i];
}
wshapir=s2*s2/ck;
return wshapir;
}
//*****************************************************************************
function shapir_francia(n,x,a) {
let e=[],i,ck,s1,s2,k,wshapir_fr,z;

s1=0;
for(i=0;i<n;i++) {
    z=(i+1-3./8.)/(n+0.25);
    a[i]=invnormaldistribution(z);
    s1=s1+a[i]*a[i];
  }
 s1=Math.sqrt(s1);

 for(i=0;i<n;i++) a[i]=a[i]/s1;
 s1=0;s2=0;
 for(i=0;i<n;i++) {
   s1+=x[i]*x[i];
   s2+=x[i];
 }
 ck=s1-(s2*s2/n);
 s2=0;
 for(i=0;i<n/2;i++) s2+=a[n-i-1]*(x[n-i-1]-x[i]);
 wshapir_fr=s2*s2/ck;
 return wshapir_fr;
}

//*********************Моделирование критического значения**wm****************
function shapir_simulation(ll,a,n,wshapir,wcrit) {

 let s1,s2,s3,i,j,x=[],ck,wm=[],zx,pvalue;
// ll=45050 число повторений

for(j=0;j<ll;j++) {
  for(i=0;i<n;i++) {
       zx=parseFloat(Math.random());
       x[i]=invnormaldistribution(zx);
  }
  s1=0;s2=0;s3=0;
  x.sort(function(a,b){return a-b});
  for(i=0;i<n;i++) {
    s1=s1+x[i]*x[i];s2=s2+x[i];
  }
    ck=s1-(s2*s2/n);

    for(i=0;i<n;i++) s3+=a[i]*x[i];
    wm[j]=s3*s3/ck;
}
    wm.sort(function(a,b){return a-b});
    for(i=0;i<ll;i++) {
         if(wm[i]>=wshapir) {
          pvalue=parseFloat((i+1.-3./8.)/(ll+0.25));
          break;
         }
    }

    wcrit[0]=wm[0.01*ll];
    wcrit[1]=wm[0.02*ll];
    wcrit[2]=wm[0.05*ll];
    wcrit[3]=wm[0.1*ll];
    wcrit[4]=wm[0.5*ll];

  return pvalue;
}
//*********************Моделирование критического значения критерия Мана****************
function mann_simulation(critname,mm,a,n,l,wmann,wcrit) {

 let s1,s2,i,j,x=[],wm=[],zx,k,k1;
//mm=45050 число повторений

  k=l-1;
  k1=parseInt(l/2)+1;

for(j=0;j<mm;j++) {
  for(i=0;i<n;i++) {
       zx=parseFloat(Math.random());
       if(critname=="Mann") x[i]=invnormaldistribution(zx);
       if(critname=="Manw") x[i]=Math.log(Math.log(1./(1.-zx)));
  }
   x.sort(function(a,b){return a-b});

   s1=0;s2=0;
  for(i=0;i<k;i++) {
    z=x[i+1]-x[i];
    s1=s1+z/a[i];
    if((i+1)>=k1) s2=s2+z/a[i];
  }
   wm[j]=s2/s1;
}
    wm.sort(function(a,b){return a-b});
    for(i=0;i<mm;i++) {
         if(wm[i]>=wmann) {
          pvalue=1.-parseFloat((i+0.5)/mm);
          break;
         }
    }

    wcrit[0]=wm[0.99*mm];
    wcrit[1]=wm[0.98*mm];
    wcrit[2]=wm[0.95*mm];
    wcrit[3]=wm[0.90*mm];
    wcrit[4]=wm[0.50*mm];

  return pvalue;
}
//*********************Моделирование функции распределения и критических значений*критериев согласия****************
function gfit_simulation(mm,n,critname,critstat,wcrit) {

 let s1,s2,i,j,x=[],cko,cp,zx,pand,z,and,wm=[],pvalue,crit,ks,wemp,pks;

//******************************************************************************
for(j=0;j<mm;j++) {
  for(i=0;i<n;i++) {
       zx=parseFloat(Math.random());
       x[i]=invnormaldistribution(zx);
  }
  s1=0;s2=0;
  x.sort(function(a,b){return a-b});
  for(i=0;i<n;i++) {
    s1=s1+x[i]*x[i];s2=s2+x[i];
  }
    cp=s2/n;
    cko=Math.sqrt((s1-cp*cp*n)/(n-1.));

 //  Select critname
    if(critname=="Anderson") {
     and=0;
      for(i=0;i<n;i++) {
        z=(x[i]-cp)/cko;
        pand=normaldistribution(z);
        and=and+((i+0.5)/n)*Math.log(pand)+(1.-(i+0.5)/n)*Math.log(1.-pand);
      }
        and=-n-2.*and;
        crit=(and-0.7/n)*(1.+3.6/n-8./(n*n));
    }

    if(critname=="KS") {
     ksm=0;
      for(i=0;i<n;i++) {
        z=(x[i]-cp)/cko;
        pks=normaldistribution(z);
        wemp=(i+0.5)/n;
        ksm=ksm+(wemp-pks)**2;
      }
        crit=(1./(12.*n)+ksm)*(1.+0.5/n);
    }
  wm[j]=crit;

 }
//*********************************************************************************
    wm.sort(function(a,b){return a-b});

    for(i=0;i<mm;i++) {
         if(wm[i]>=critstat) {
          pvalue=parseFloat((i+0.5)/mm);
          break;
         }
    }

    wcrit[0]=wm[0.99*mm];
    wcrit[1]=wm[0.95*mm];
    wcrit[2]=wm[0.90*mm];
    wcrit[3]=wm[0.85*mm];
    wcrit[4]=wm[0.50*mm];

  return pvalue;
}


//***********Стандартные статистические функции************
/*
minint
maxint
gammastirf
gamma
lngamma
incompletebetapseta
incompletebetafe
incompletebetafe2
incompletebeta
fdistribution
chisquaredistribution
studenttdistr
incompletegamma
incompletegammac
normaldistribution

invnormaldistribution
invincompletegammac
invchisquaredistribution
invfdistribution
invincompletebeta
invstudenttdistr
*/
//*****************************************************************************
function incompletebetaps(a,b,x,maxgam) {
  sg=0;
  ai=1/a;
  u=(1-b)*x;
  v=u/(a+1);
  t1=v;
  t=u;
  n=2;
  s=0;
  z=machineepsilon * ai;
  while(Math.abs(v)>z) {
  u=(n-b)*x/n;t=t*u;v=t/(a+n);v=v/Math.pow(10,25);s=s+v;n++;v=v*Math.pow(10,25);
  }
  s=s*Math.pow(10,25);
  //alert(s);
  //return;
  s=s+t1;
  s=s+ai;
  u=a*Math.log(x);
  if (((a+b)<maxgam) && (Math.abs(u)<Math.log(maxrealnumber))) {
    t=gamma(a+b)/(gamma(a)*gamma(b));
    s=s*t*Math.pow(x,a);
   } 
    else {
    t=lngamma(a+b,sg)-lngamma(a,sg)-lngamma(b,sg)+u+Math.log(s);
    if(t<Math.log(minrealnumber)) {
      s=0;
    } else {
      s=Math.exp(t);
    }
  }
  return s;
}
//*****************************************************************************
function incompletebetafe(a,b,x,big,biginv) {

    k1 = a;
    k2 = a + b;
    k3 = a;
    k4 = a + 1;
    k5 = 1;
    k6 = b - 1;
    k7 = k4;
    k8 = a + 2;
    pkm2 = 0;
    qkm2 = 1;
    pkm1 = 1;
    qkm1 = 1;
    ans = 1;
    r = 1;
    n = 0;
    thresh = 3 * machineepsilon;
    do {
        xk = -(x * k1 * k2 / (k3 * k4));
        pk = pkm1 + pkm2 * xk;
        qk = qkm1 + qkm2 * xk;
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;
        xk = x * k5 * k6 / (k7 * k8);
        pk = pkm1 + pkm2 * xk;
        qk = qkm1 + qkm2 * xk;
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;
        if(qk!=0) {
            r = pk / qk;
        }
        if(r!=0) {
            t = Math.abs((ans - r) / r);
            ans = r;
        } else {
            t = 1;
        }
        if(t<thresh) {
            break;
        }
        k1 = k1 + 1;
        k2 = k2 + 1;
        k3 = k3 + 2;
        k4 = k4 + 2;
        k5 = k5 + 1;
        k6 = k6 - 1;
        k7 = k7 + 2;
        k8 = k8 + 2;
        if((Math.abs(qk) + Math.abs(pk))>big) {
            pkm2*=biginv;
            pkm1*=biginv;
            qkm2*=biginv;
            qkm1*=biginv;
        }
        if((Math.abs(qk)<biginv) || (Math.abs(pk)<biginv)) {
            pkm2*=big;
            pkm1*=big;
            qkm2*=big;
            qkm1*=big;
        }
        n++;
    }
    while(n<300);
    result = ans;
    return result;
}
//*****************************************************************************

function incompletebetafe2(a,b,x,big,biginv) {


    k1 = a;
    k2 = b - 1;
    k3 = a;
    k4 = a + 1;
    k5 = 1;
    k6 = a + b;
    k7 = a + 1;
    k8 = a + 2;
    pkm2 = 0;
    qkm2 = 1;
    pkm1 = 1;
    qkm1 = 1;
    z = x / (1 - x);
    ans = 1;
    r = 1;
    n = 0;
    thresh = 3*machineepsilon;
    do {
        xk = -(z * k1 * k2 / (k3 * k4));
        pk = pkm1 + pkm2 * xk;
        qk = qkm1 + qkm2 * xk;
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;
        xk = z * k5 * k6 / (k7 * k8);
        pk = pkm1 + pkm2 * xk;
        qk = qkm1 + qkm2 * xk;
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;
        if(qk!=0) {
            r = pk / qk;
        }
        if(r!=0) {
            t = Math.abs((ans - r) / r);
            ans = r;
        } else {
            t = 1;
        }
        if(t<thresh) {
            break;
        }
        k1++;
        k2--;
        k3 = k3 + 2;
        k4 = k4 + 2;
        k5++;
        k6++;
        k7 = k7 + 2;
        k8 = k8 + 2;
        if((Math.abs(qk)+Math.abs(pk))>big) {
            pkm2*=biginv;
            pkm1*=biginv;
            qkm2*=biginv;
            qkm1*=biginv;
        }
        if((Math.abs(qk)<biginv) || (Math.abs(pk)<biginv)) {
            pkm2*=big;
            pkm1*=big;
            qkm2*=big;
            qkm1*=big;
        }
        n++;
    }
    while(n<300);
    result = ans;
    return result;
}
//*****************************************************************************

function incompletebeta(a,b,x) {

    big = 4.5035996273705*Math.pow(10,15);
    biginv = 2.22044604925031*Math.pow(10,-16);
    maxgam = 171.624376956303;
    minlog = Math.log(minrealnumber);
    maxlog = Math.log(maxrealnumber);

    if(x==0) {
      result=0;
      return result;
    }
    if(x==1) {
      result=1;
      return result;
    }
    flag=0;
    if((b*x<=1) && (x<=0.95)) {
      result=incompletebetaps(a, b, x, maxgam);
      return result;
    }
    w=1-x;

    if(x>(a/(a+b))) {
      flag=1;
      t=a;
      a=b;
      b=t;
      xc=x;
      x=w;
    } else {
      xc=w;
    }

    if((flag==1) && ((b*x<=1)) && (x<=0.95)) {
      t=incompletebetaps(a, b, x, maxgam);
      if(t<=machineepsilon) {
        result=1-machineepsilon;
      } else {
        result=1-t;
      }
      return result;
    }

    y=x*(a+b-2)-(a-1);
    if(y<0) {
      w = incompletebetafe(a, b, x, big, biginv);
    } else {
      w = incompletebetafe2(a, b, x, big, biginv) / xc;
    }
    y=a*Math.log(x);
    t=b*Math.log(xc);
    if(((a+b)<maxgam) && (Math.abs(y)<maxlog) && (Math.abs(t)<maxlog)) {
      t=Math.pow(xc,b);
      t*=Math.pow(x,a);
      t/=a;
      t*=w;
      t*=(gamma(a+b)/(gamma(a)*gamma(b)));
      if(flag==1) {
        if(t<=machineepsilon) {
          result=1-machineepsilon;
        } else {
          result=1-t;
        }
      } else {
        result=t;
      }
      return result;
    }
    y+=t+lngamma(a+b,sg)-lngamma(a,sg)-lngamma(b,sg);
    y+=Math.log(w/a);
    if(y<minlog) {
      t=0;
    } else {
      t=Math.exp(y);
    }
    if(flag==1) {
      if(t<=machineepsilon) {
        t=1-machineepsilon;
      } else {
        t=1-t;
      }
    }
    result=t;
    return result;
}



//*****************************************************************************
function fdistribution(a,b,x) {
    //w=(a*x)/(b+a*x);
    result=incompletebeta(0.5*a,0.5*b,a*x/(b+a*x));
  return result;
}
//*****************************************************************************
function studenttdistr(k,t) {


  if(t==0) {
    result=0.5;
    return result;
  }
  /*
  if(t<-2) {
    rk=k;
    z=rk/(rk+t*t);
    result=0.5*incompletebeta(0.5*rk,0.5,z);
    return result;
  }
  */
  
  if(t<0) {
    x=-t;
  } 
  else 
  {
   x=t;
  }
  rk=k;
  z=1+x*x/rk;
  if((k%2)!=0) {
    xsqk=x/Math.sqrt(rk);
    p=Math.atan(xsqk);
    if(k>1) {
      f=1;
      tz=1;
      j=3;
      while((j<=(k-2)) && ((tz/f)>machineepsilon)) {
        tz*=(j-1)/(z*j);
        f+=tz;
        j+=2;
      }
      p+=f*xsqk/z;
    }
    p*=2/Math.PI;
  } else {
    f=1;tz=1;j=2;
    while(j<=(k-2) && ((tz/f)>machineepsilon)) {
      tz*=(j-1)/(z*j);
      f+=tz;
      j+=2;
    }
    p=f*x/Math.sqrt(z*rk);
  }

  if(t<0) {
   result=0.5-0.5*p;
   return result;
  }
   
  result=0.5+0.5*p;
  return result;
 
}

//*****************************************************************************
function invincompletebeta(a,b,y)   {
             result = 0;
             aaa = 0;
             bbb = 0;
             y0 = 0;
             d = 0;
             yyy = 0;
             x = 0;
             x0 = 0;
             x1 = 0;
             lgm = 0;
             yp = 0;
             di = 0;
             dithresh = 0;
             yl = 0;
             yh = 0;
             xt = 0;
             i = 0;
             rflg = 0;
             dir = 0;
             nflg = 0;
             s = 0;
             mainlooppos = 0;
             ihalve = 0;
             ihalvecycle = 0;
             newt = 0;
             newtcycle = 0;
             breaknewtcycle = 0;
             breakihalvecycle = 0;

            i = 0;
            
            
            if(y==0)
            {
                result = 0;
                return result;
            }
            if(y==1.0)
            {
                result = 1;
                return result;
            }
            
           
            dithresh = 0;
            rflg = 0;
            aaa = 0;
            bbb = 0;
            y0 = 0;
            x = 0;
            yyy = 0;
            lgm = 0;
            dir = 0;
            di = 0;
            x0 = 0;
            yl = 0;
            x1 = 1.0;
            yh = 1.0;
            nflg = 0;
            mainlooppos = 0;
            ihalve = 1;
            ihalvecycle = 2;
            newt = 3;
            newtcycle = 4;
            breaknewtcycle = 5;
            breakihalvecycle = 6;
            
          
            while(true)  {
                
                if( mainlooppos==0 )
                {
                    if(a<=1.0 || b<=1.0)
                    {
                        dithresh = 0.000001;
                        rflg = 0;
                        aaa = a;
                        bbb = b;
                        y0 = y;
                        x = aaa/(aaa+bbb);
                        yyy = incompletebeta(aaa, bbb, x);
                        mainlooppos = ihalve;
                        continue;
                    }
                    else
                    {
                        dithresh = 0.0001;
                    }
                    yp = -invnormaldistribution(y);
                    if(y>0.5)
                    {
                        rflg = 1;
                        aaa = b;
                        bbb = a;
                        y0 = 1.0-y;
                        yp = -yp;
                    }
                    else
                    {
                        rflg = 0;
                        aaa = a;
                        bbb = b;
                        y0 = y;
                    }
                    lgm = (yp*yp-3.0)/6.0;
                    x = 2.0/(1.0/(2.0*aaa-1.0)+1.0/(2.0*bbb-1.0));
                    d = yp*Math.sqrt(x+lgm)/x-(1.0/(2.0*bbb-1.0)-1.0/(2.0*aaa-1.0))*(lgm+5.0/6.0-2.0/(3.0*x));
                    d = 2.0*d;
                    if(d<Math.log(minrealnumber))
                    {
                        x = 0;
                        break;
                    }
                    x = aaa/(aaa+bbb*Math.exp(d));
                    yyy = incompletebeta(aaa,bbb,x);
                    yp = (yyy-y0)/y0;
                    if(Math.abs(yp)<0.2)
                    {
                        mainlooppos = newt;
                        continue;
                    }
                    mainlooppos = ihalve;
                    continue;
                }
                
                
                if( mainlooppos==ihalve )
                {
                    dir = 0;
                    di = 0.5;
                    i = 0;
                    mainlooppos = ihalvecycle;
                    continue;
                }
                
               
                if( mainlooppos==ihalvecycle )
                {
                    if( i<=99 )
                    {
                        if( i!=0 )
                        {
                            x = x0+di*(x1-x0);
                            if(x==1.0)
                            {
                                x = 1.0-machineepsilon;
                            }
                            if(x==0)
                            {
                                di = 0.5;
                                x = x0+di*(x1-x0);
                                if(x==0)
                                {
                                    break;
                                }
                            }
                            yyy = incompletebeta(aaa, bbb, x);
                            yp = (x1-x0)/(x1+x0);
                            if(Math.abs(yp)<dithresh)
                            {
                                mainlooppos = newt;
                                continue;
                            }
                            yp = (yyy-y0)/y0;
                            if(Math.abs(yp)<dithresh)
                            {
                                mainlooppos = newt;
                                continue;
                            }
                        }
                        if(yyy<y0)
                        {
                            x0 = x;
                            yl = yyy;
                            if( dir<0 )
                            {
                                dir = 0;
                                di = 0.5;
                            }
                            else
                            {
                                if( dir>3 )
                                {
                                    di = 1.0-(1.0-di)*(1.0-di);
                                }
                                else
                                {
                                    if( dir>1 )
                                    {
                                        di = 0.5*di+0.5;
                                    }
                                    else
                                    {
                                        di = (y0-yyy)/(yh-yl);
                                    }
                                }
                            }
                            dir = dir+1;
                            if(x0>0.75)
                            {
                                if( rflg==1 )
                                {
                                    rflg = 0;
                                    aaa = a;
                                    bbb = b;
                                    y0 = y;
                                }
                                else
                                {
                                    rflg = 1;
                                    aaa = b;
                                    bbb = a;
                                    y0 = 1.0-y;
                                }
                                x = 1.0-x;
                                yyy = incompletebeta(aaa, bbb, x);
                                x0 = 0;
                                yl = 0;
                                x1 = 1.0;
                                yh = 1.0;
                                mainlooppos = ihalve;
                                continue;
                            }
                        }
                        else
                        {
                            x1 = x;
                            if( rflg==1 && x1<machineepsilon)
                            {
                                x = 0;
                                break;
                            }
                            yh = yyy;
                            if( dir>0 )
                            {
                                dir = 0;
                                di = 0.5;
                            }
                            else
                            {
                                if( dir<-3 )
                                {
                                    di = di*di;
                                }
                                else
                                {
                                    if( dir<-1 )
                                    {
                                        di = 0.5*di;
                                    }
                                    else
                                    {
                                        di = (yyy-y0)/(yh-yl);
                                    }
                                }
                            }
                            dir = dir-1;
                        }
                        i = i+1;
                        mainlooppos = ihalvecycle;
                        continue;
                    }
                    else
                    {
                        mainlooppos = breakihalvecycle;
                        continue;
                    }
                }
                
                
                if( mainlooppos==breakihalvecycle )
                {
                    if( x0>=1.0)
                    {
                        x = 1.0-machineepsilon;
                        break;
                    }
                    if(x<=0)
                    {
                        x = 0;
                        break;
                    }
                    mainlooppos = newt;
                    continue;
                }
                
               
                if( mainlooppos==newt )
                {
                    if( nflg!=0 )
                    {
                        break;
                    }
                    nflg = 1;
                    lgm = lngamma(aaa+bbb,s)-lngamma(aaa,s)-lngamma(bbb,s);
                    i = 0;
                    mainlooppos = newtcycle;
                    continue;
                }
                
               
                if( mainlooppos==newtcycle )
                {
                    if( i<=7 )
                    {
                        if( i!=0 )
                        {
                            yyy = incompletebeta(aaa, bbb, x);
                        }
                        if(yyy<yl)
                        {
                            x = x0;
                            yyy = yl;
                        }
                        else
                        {
                            if(yyy>yh)
                            {
                                x = x1;
                                yyy = yh;
                            }
                            else
                            {
                                if(yyy<y0)
                                {
                                    x0 = x;
                                    yl = yyy;
                                }
                                else
                                {
                                    x1 = x;
                                    yh = yyy;
                                }
                            }
                        }
                        if(x==1.0 || x==0)
                        {
                            mainlooppos = breaknewtcycle;
                            continue;
                        }
                        d = (aaa-1.0)*Math.log(x)+(bbb-1.0)*Math.log(1.0-x)+lgm;
                        if(d<Math.log(minrealnumber))
                        {
                            break;
                        }
                        if(d>Math.log(maxrealnumber))
                        {
                            mainlooppos = breaknewtcycle;
                            continue;
                        }
                        d = Math.exp(d);
                        d = (yyy-y0)/d;
                        xt = x-d;
                        if(xt<=x0)
                        {
                            yyy = (x-x0)/(x1-x0);
                            xt = x0+0.5*yyy*(x-x0);
                            if(xt<=0)
                            {
                                mainlooppos = breaknewtcycle;
                                continue;
                            }
                        }
                        if(xt>=x1)
                        {
                            yyy = (x1-x)/(x1-x0);
                            xt = x1-0.5*yyy*(x1-x);
                            if(xt>=1)
                            {
                                mainlooppos = breaknewtcycle;
                                continue;
                            }
                        }
                        x = xt;
                        if(Math.abs(d/x)<128.0*machineepsilon)
                        {
                            break;
                        }
                        i = i+1;
                        mainlooppos = newtcycle;
                        continue;
                    }
                    else
                    {
                        mainlooppos = breaknewtcycle;
                        continue;
                    }
                }
                
                
                if( mainlooppos==breaknewtcycle )
                {
                    dithresh = 256.0*machineepsilon;
                    mainlooppos = ihalve;
                    continue;
                }
            }
            
           
            if( rflg!=0 )
            {
                if(x<=machineepsilon)
                {
                    x = 1.0-machineepsilon;
                }
                else
                {
                    x = 1.0-x;
                }
            }
            result = x;
            return result;
        }

//*************************************************************************
function invstudenttdistr(k,p) {

if (p==0.5) {t=0;return t};

//if((p>0.25) && (p<0.75)) {
  z=1-2*p;
  z=invincompletebeta(0.5,0.5*k,Math.abs(z));
  t=Math.sqrt(k*z/(1-z));
  if (p<0.5) t=-t;
  return t;
//}


/*
rflg=-1;
if(p>0.5) {
  p=1-pt;
  rflg=1;
}
z=invincompletebeta(0.5*k,0.5,2*p);

if((maxrealnumber*z)<k) {
  result=rflg*maxrealnumber;
  return result;
}

t=Math.sqrt(k/z-k);
result=rflg*t;
return result;
*/
}


//*****************************************************************************
function minint(m1,m2) {
  if(m1<m2) {
    return m1;
  } else {
    return m2;
  }
}
//*****************************************************************************
function maxint(m1,m2) {
  if(m1>m2) {
    return m1;
  } else {
    return m2;
  }
}
//*****************************************************************************


function invnormaldistribution(y0)  {
            var result = 0;
            var x = 0;
            var y = 0;
            var z = 0;
            var y2 = 0;
            var x0 = 0;
            var x1 = 0;
            var code = 0;
            var p0 = 0;
            var q0 = 0;
            var p1 = 0;
            var q1 = 0;
            var p2 = 0;
            var q2 = 0;

            
            if(y0<=0) {
                result = -maxrealnumber;
                return result;
            }
            if(y0>=1) {
                result = maxrealnumber;
                return result;
            }
            code = 1;
            y = y0;
            if(y>1.0-expm2) {
                y = 1.0-y;
                code = 0;
                }
            if(y>expm2) {
                y = y-0.5;
                y2 = y*y;
                p0 = -59.9633501014107895267;
                p0 = 98.0010754185999661536+y2*p0;
                p0 = -56.6762857469070293439+y2*p0;
                p0 = 13.9312609387279679503+y2*p0;
                p0 = -1.23916583867381258016+y2*p0;
                q0 = 1;
                q0 = 1.95448858338141759834+y2*q0;
                q0 = 4.67627912898881538453+y2*q0;
                q0 = 86.3602421390890590575+y2*q0;
                q0 = -225.462687854119370527+y2*q0;
                q0 = 200.260212380060660359+y2*q0;
                q0 = -82.0372256168333339912+y2*q0;
                q0 = 15.9056225126211695515+y2*q0;
                q0 = -1.18331621121330003142+y2*q0;
                x = y+y*y2*p0/q0;
                x = x*s2pi;
                result = x;
                return result;
            }
            
            zz=Math.log(y);
            x = Math.sqrt(-(2.0*zz));
            x0 = x-Math.log(x)/x;
            z = 1.0/x;
            
            if(x<8.0) {
                p1 = 4.05544892305962419923;
                p1 = 31.5251094599893866154+z*p1;
                p1 = 57.1628192246421288162+z*p1;
                p1 = 44.0805073893200834700+z*p1;
                p1 = 14.6849561928858024014+z*p1;
                p1 = 2.18663306850790267539+z*p1;
                p1 = -(1.40256079171354495875*0.1)+z*p1;
                p1 = -(3.50424626827848203418*0.01)+z*p1;
                p1 = -(8.57456785154685413611*0.0001)+z*p1;
                q1 = 1;
                q1 = 15.7799883256466749731+z*q1;
                q1 = 45.3907635128879210584+z*q1;
                q1 = 41.3172038254672030440+z*q1;
                q1 = 15.0425385692907503408+z*q1;
                q1 = 2.50464946208309415979+z*q1;
                q1 = -(1.42182922854787788574*0.1)+z*q1;
                q1 = -(3.80806407691578277194*0.01)+z*q1;
                q1 = -(9.33259480895457427372*0.0001)+z*q1;
                x1 = z*p1/q1;
            }
            else  {
                p2 = 3.23774891776946035970;
                p2 = 6.91522889068984211695+z*p2;
                p2 = 3.93881025292474443415+z*p2;
                p2 = 1.33303460815807542389+z*p2;
                p2 = 2.01485389549179081538*0.1+z*p2;
                p2 = 1.23716634817820021358*0.01+z*p2;
                p2 = 3.01581553508235416007*0.0001+z*p2;
                p2 = 2.65806974686737550832*0.000001+z*p2;
                p2 = 6.23974539184983293730*0.000000001+z*p2;
                q2 = 1;
                q2 = 6.02427039364742014255+z*q2;
                q2 = 3.67983563856160859403+z*q2;
                q2 = 1.37702099489081330271+z*q2;
                q2 = 2.16236993594496635890*0.1+z*q2;
                q2 = 1.34204006088543189037*0.01+z*q2;
                q2 = 3.28014464682127739104*0.0001+z*q2;
                q2 = 2.89247864745380683936*0.000001+z*q2;
                q2 = 6.79019408009981274425*0.000000001+z*q2;
                x1 = z*p2/q2;
            }
            x = x0-x1;
            if( code!=0 ) x = -x;
            result = x;
            return result;
          }


//************************************************************************
 function errorfunctionc(x) {
            var result = 0;
            var p = 0;
            var q = 0;

            if(x<0) {
                result = 2-errorfunctionc(-x);
                return result;
            }
            if(x<0.5) {
                result = 1.0-errorfunction(x);
                return result;
            }
            if(x>=10) {
                result = 0;
                return result;
            }
            p = 0.0;
            p = 0.5641877825507397413087057563+x*p;
            p = 9.675807882987265400604202961+x*p;
            p = 77.08161730368428609781633646+x*p;
            p = 368.5196154710010637133875746+x*p;
            p = 1143.262070703886173606073338+x*p;
            p = 2320.439590251635247384768711+x*p;
            p = 2898.0293292167655611275846+x*p;
            p = 1826.3348842295112592168999+x*p;
            q = 1.0;
            q = 17.14980943627607849376131193+x*q;
            q = 137.1255960500622202878443578+x*q;
            q = 661.7361207107653469211984771+x*q;
            q = 2094.384367789539593790281779+x*q;
            q = 4429.612803883682726711528526+x*q;
            q = 6089.5424232724435504633068+x*q;
            q = 4958.82756472114071495438422+x*q;
            q = 1826.3348842295112595576438+x*q;
            result = Math.exp(-x*x)*p/q;
            return result;
        }

//*****************************************************************
 function errorfunction(x) {
            var result = 0;
            var xsq = 0;
            var s = 0;
            var p = 0;
            var q = 0;

            s = Math.sign(x);
            x = Math.abs(x);
            if(x<0.5) {
                xsq = x*x;
                p = 0.007547728033418631287834;
                p = -0.288805137207594084924010+xsq*p;
                p = 14.3383842191748205576712+xsq*p;
                p = 38.0140318123903008244444+xsq*p;
                p = 3017.82788536507577809226+xsq*p;
                p = 7404.07142710151470082064+xsq*p;
                p = 80437.3630960840172832162+xsq*p;
                q = 0.0;
                q = 1.0+xsq*q;
                q = 38.0190713951939403753468+xsq*q;
                q = 658.070155459240506326937+xsq*q;
                q = 6379.60017324428279487120+xsq*q;
                q = 34216.5257924628539769006+xsq*q;
                q = 80437.3630960840172826266+xsq*q;
                result = s*1.1283791670955125738961589031*x*p/q;
                return result;
            }
            if(x>=10) {
                result = s;
                return result;
            }
            result = s*(1-errorfunctionc(x));
            return result;
        }

//**************************************************************************

 function normaldistribution(x) {
            result = 0;
            result = 0.5*(errorfunction(x/1.41421356237309504880)+1);
            return result;
        }

//*************************************************************************
function invfdistribution(a,b,y) {
        result = 0;
        w = 0;
        w = incompletebeta(0.5*b, 0.5*a, 0.5);
        if((w>y) || (y<0.001))
            {
                w = invincompletebeta(0.5*b, 0.5*a, y);
                result = (b-b*w)/(a*w);
            }
            else
            {
                w = invincompletebeta(0.5*a, 0.5*b, 1.0-y);
                result = b*w/(a*(1.0-w));
            }
            return result;
        }

//Gamma class
//*******************************************************************
function gammastirf(x) {
 
var stir,w,y,v;   
w = 1./x;
stir = 7.87311395793093628397*0.0001;
stir = -2.29549961613378126380*0.0001+w*stir;
stir = -2.68132617805781232825*0.001+w*stir;
stir = 3.47222221605458667310*0.001+w*stir;
stir = 0.0833333333333482257126+w*stir;
w = 1.+w*stir;
y = Math.exp(x);
if (x>143.01608)  {
   v = Math.pow(x, 0.5*x-0.25);
   y = v*v/y;
   }
 else  {
   y = Math.pow(x, x-0.5)/y;
   }
   result =2.50662827463100050242*y*w;
   return result;
}

//*****************************************************************************
function gamma(x) {
sgngam=1;
q=Math.abs(x);
if(q>33) {
  if(x<0) {
    p=Math.floor(q);
    i=Math.round(p);
    if(i%2==0) sgngam=-1;
    z=q-p;
    if(z>0.5) {
      p++;
      z=q-p;
    }
    z=q*Math.sin(Math.PI*z);
    z=Math.abs(z);
    z=Math.PI/(z*gammastirf(q));
  } else {
    z=gammastirf(x);
  }
  result=sgngam*z;
  return result;
}

z=1;
while(x>=3) {
  x--;
  z*=x;
}
while(x<0) {
  if(x>-0.000000001) {
    result=z/((1+0.577215664901533*x)*x);
    return result;
  }
  z/=x;
  x++;
}
while(x<2) {
  if(x<0.000000001) {
    result = z/((1 + 0.577215664901533*x)*x);
    return result;
  }
  z/=x;
  x++;
}

if(x==2) {
  result=z;
  return result;
}
    x = x - 2;
    pp = 1.60119522476752*0.0001;
    pp = 1.19135147006586*0.001 + x * pp;
    pp = 1.04213797561762*0.01 + x * pp;
    pp = 4.76367800457137*0.01 + x * pp;
    pp = 0.207448227648436 + x * pp;
    pp = 0.494214826801497 + x * pp;
    pp = 1 + x*pp;
    qq = -2.3158187332412*0.00001;
    qq = 5.39605580493303*0.0001 + x*qq;
    qq = -4.45641913851797*0.001 + x*qq;
    qq = 0.011813978522206 + x*qq;
    qq = 3.58236398605499*0.01 + x*qq;
    qq = -0.234591795718243 + x*qq;
    qq = 7.14304917030273*0.01 + x*qq;
    qq = 1 + x*qq;
    result = z * pp / qq;
    return result;
}
//*****************************************************************************

function invincompletegammac(a,y0) {

            tmp=0;     
            x0 = igammabignumber;
            yl = 0;
            x1 = 0;
            yh = 1;
            dithresh = 5*igammaepsilon;
            d = 1/(9*a);
            y = 1-d-invnormaldistribution(y0)*Math.sqrt(d);
            x = a*y*y*y;
            lgm = lngamma(a,tmp);
            i = 0;
            while( i<10 )
            {
                if(x>x0 || x<x)
                {
                    d = 0.0625;
                    break;
                }
                y = incompletegammac(a, x);
                if(y<yl || y>yh)
                {
                    d = 0.0625;
                    break;
                }
                if(y<y0)
                {
                    x0 = x;
                    yl = y;
                }
                else
                {
                    x1 = x;
                    yh = y;
                }
                d = (a-1)*Math.log(x)-x-lgm;
                if (d<-709.78271289338399)
                {
                    d = 0.0625;
                    break;
                }
                d = -Math.exp(d);
                d = (y-y0)/d;
                if(Math.abs(d/x)<igammaepsilon)
                {
                    result = x;
                    return result;
                }
                x = x-d;
                i = i+1;
            }
            if (x0==igammabignumber) {
                if (x<=0) x = 1;
                 while(x0==igammabignumber) {
                    x = (1+d)*x;
                    y = incompletegammac(a, x);
                    if (y<y0)   {
                        x0 = x;
                        yl = y;
                        break;
                    }
                    d = d+d;
                }
            }
            d = 0.5;
            dir = 0;
            i = 0;
            while (i<400) {
                x = x1+d*(x0-x1);
                y = incompletegammac(a, x);
                lgm = (x0-x1)/(x1+x0);
                if (Math.abs(lgm)<dithresh) break;
                lgm = (y-y0)/y0;
                if (Math.abs(lgm)<dithresh) break;
                if (x<=0)  break;
                if(y>=y0) {
                    x1 = x;
                    yh = y;
                    if (dir<0) {
                        dir = 0;d = 0.5;
                    }
                    else
                    {
                        if (dir>1) {
                            d = 0.5*d+0.5;
                        }
                        else  {
                            d = (y0-yl)/(yh-yl);
                        }
                    }
                    dir = dir+1;
                }
                else
                {
                    x0 = x;
                    yl = y;
                    if (dir>0)
                    {
                        dir = 0;
                        d = 0.5;
                    }
                    else
                    {
                        if (dir<-1)
                        {
                            d = 0.5*d;
                        }
                        else
                        {
                            d = (y0-yl)/(yh-yl);
                        }
                    }
                    dir = dir-1;
                }
                i = i+1;
            }
            result = x;
            return result;
        }

//*****************************************************************************
function invchisquaredistribution(v, y) {
result = 2*invincompletegammac(0.5*v,y);
return result;
}


//**************************************************************************

function incompletegamma(a,x) {
tmp=0;
if(x<=0 || a<=0) {
  result=0;
  return result;
}
if(x>1 && x>a) {
   result=1-incompletegammac(a, x);
   return result;
}
ax=a*Math.log(x)-x-lngamma(a,tmp);
if(ax<-709.782712893384) {
  result=0;
  return result;
}
ax=Math.exp(ax);
r=a;c=1;ans=1;
do {
  r++;c*=x/r;ans+=c;
} while((c/ans)>igammaepsilon);
result=ans*ax/a;
return result;
}


//*****************************************************************************

 function incompletegammac(a,x) {
        tmp=0;
        if(x<=0 || a<=0) {
                result = 1;
                return result;
            }
            if( x<1 || x<a)
            {
                result = 1-incompletegamma(a, x);
                return result;
            }
            ax = a*Math.log(x)-x-lngamma(a,tmp);
            if(ax<-709.78271289338399)
            {
                result = 0;
                return result;
            }
            ax = Math.exp(ax);
            y = 1-a;
            z = x+y+1;
            c = 0;
            pkm2 = 1;
            qkm2 = x;
            pkm1 = x+1;
            qkm1 = z*x;
            ans = pkm1/qkm1;
            do
            {
                c = c+1;
                y = y+1;
                z = z+2;
                yc = y*c;
                pk = pkm1*z-pkm2*yc;
                qk = qkm1*z-qkm2*yc;
                if (qk!=0)
                {
                    r = pk/qk;
                    t = Math.abs((ans-r)/r);
                    ans = r;
                }
                else
                {
                    t = 1;
                }
                pkm2 = pkm1;
                pkm1 = pk;
                qkm2 = qkm1;
                qkm1 = qk;
                if(Math.abs(pk)>igammabignumber)
                {
                    pkm2 = pkm2*igammabignumberinv;
                    pkm1 = pkm1*igammabignumberinv;
                    qkm2 = qkm2*igammabignumberinv;
                    qkm1 = qkm1*igammabignumberinv;
                }
            } while(t>igammaepsilon);
            result = ans*ax;
            return result;
        }


function lngamma(x,sgngam) {

    sgngam = 1;
    
    if(x<-34) {
      q=-x;
      w=lngamma(q,tmp);
      p=Math.floor(q);
      i=Math.round(p);
      if(i%2==0) {
        sgngam=-1;
      } else {
        sgngam=1;
      }
      z=q-p;
      if(z>0.5) {
        p++;
        z=p-q;
      }
      z=q*Math.sin(Math.PI*z);
       result=loqpi-Math.log(z)-w;
       return result;
    }
    if(x<13) {
      z=1;
      p=0;
      u=x;
      while(u>=3) {
        p--;
        u=x+p;
        z*=u;
      }
      while(u<2) {
        z/=u;
        p++;
        u=x+p;
      }
      if(z<0) {
        sgngam=-1;
        z=-z;
      } else {
        sgngam=1;
      }
      if(u==2) {
        result=Math.log(z);
        return result;
      }
        p = p - 2;
        x = x + p;
        b = -1378.25152569121;
        b = -38801.6315134638 + x * b;
        b = -331612.992738871 + x * b;
        b = -1162370.97492762 + x * b;
        b = -1721737.0082084 + x * b;
        b = -853555.664245765 + x * b;
        c = 1;
        c = -351.815701436523 + x * c;
        c = -17064.2106651881 + x * c;
        c = -220528.590553854 + x * c;
        c = -1139334.44367983 + x * c;
        c = -2532523.07177583 + x * c;
        c = -2018891.41433533 + x * c;
        p = x * b / c;
        result = Math.log(z) + p;
        return result;
    }
    q=(x-0.5)*Math.log(x)-x+ls2pi;
    if(x>100000000) {
      result=q;
      return result;
    }
    p=1/(x*x);
    if(x>=1000) {
      q+=((7.93650793650794 * 0.0001 * p - 2.77777777777778 * 0.001) * p + 0.0833333333333333) / x;
    } else {
        a = 8.11614167470508 * 0.0001;
        a = -(5.95061904284301 * 0.0001) + p * a;
        a = 7.93650340457717 * 0.0001 + p * a;
        a = -(2.777777777301 * 0.001) + p * a;
        a = 8.33333333333332 * 0.01 + p * a;
        q = q + a / x;
    }
    result=q;
    return result;
}
//*************************************************************
  function chisquaredistribution(v,x) {
            result = incompletegamma(v/2.0, x/2.0);
            return result;
        }
//*********************Weibull********************************************
function weibulldistribution(x,x0,b,c) {
   return 1-Math.exp(-Math.pow((x-x0)/c,b));
}
//**********************
function weibulldensity(x,x0,b,c) {
  return (b/c)*Math.pow((x-x0)/c,b-1)*Math.exp(-Math.pow((x-x0)/c,b));
}
//***********************
function invweibulldistribution_1(p,x0,b,c) {
 if(p<=0 || p>=1 || x<=x0 || c<=0) return 0;
 zw=Math.log(Math.log(1/(1-p)));aw=Math.log(c);sw=1/b;
 x=aw+zw*sw;x=Math.exp(x)+x0;
 return x;
}
//*********************
function normaldensity(x) {
  return Math.exp(-(x*x/2))/Math.sqrt(2*Math.PI);
}

function invweibulldistribution(p) {
 x=Math.log(Math.log(1/(1-p)));
 return x;
}
//************************Cumulative Distribution*****************************
function cum(xcum,fcum,ycum) {
 n=xcum.length;
    for (i = 0;i<n-1;i++) {
      for (j = i+1;j<n;j++) {
          if (Math.abs(xcum[i])>Math.abs(xcum[j])) {
           xr = xcum[i];xcum[i]=xcum[j];xcum[j]=xr;
          }
   }}
//*******************************************************
/*
//m-number of bases; nx[m]-number of censorized values on the base
 var nx=[],xb=[];
 m=0;ni=0;
 for (i=0;i<n;i++) {
  if  (xcum[i]<0) {
    ni++;
    }
   else {
    nx[m]=ni;ni=0;    
    xb[m]=xcum[i];m++;
   }
  }
*/
//*************************************************************
   //fcum[n]-cumulative distribution function; ycum[j] - failures
   f0=0.5/(n+1);Neqv=n;j=0;ni=0;
  for (i=0;i<n;i++) {
         if(xcum[i]<0) {
            ni++;
            }
            else {
           Neqv=Neqv-ni/(1-f0);ds=0;
           if(Neqv>0) ds=1/(Neqv+1);
           fcum[j]=f0+ds;f0=fcum[j];ycum[j]=xcum[i];j++;ni=0;
       }
   }  
 }
//***************************************************
function OrderNew(ts,p,n,outorder) {
   r=p*(n+1);
   if (ts=="normal") {zp=invnormaldistribution(p);dp=normaldensity(zp);dpdz=-dp*zp;}
   if (ts=="weibull") {zp=Math.log(Math.log(1./(1.-p)));dp=Math.exp(zp-Math.exp(zp));dpdz=dp*(1-Math.exp(zp));}
   d1=r*(n-r+1)/((n+1)*(n+1)*(n+2));
   d2=-0.5*d1*dpdz/(dp*dp*dp); 
   er=zp+d2; 
   dr=d1/(dp*dp)-d2*d2;
   outorder[0]=er;outorder[1]=dr;outorder[2]=zp;outorder[3]=dp;outorder[4]=dpdz;
}
//******************************************************************************
function Mnk_regress(kx,wreg,xreg,yreg,xparam) {
      a=0;xcp=0;znam=0;s1=0;s2=0;
    for(i=0;i<kx;i++) {
       a+=yreg[i]*wreg[i]; 
       xcp+=xreg[i]*wreg[i];
       znam+=wreg[i];
   }
       a=a/znam;xcp=xcp/znam;
       for(i=0;i<kx;i++) {
        s1+=yreg[i]*(xreg[i]-xcp)*wreg[i];
        s2+=(xreg[i]-xcp)*(xreg[i]-xcp)*wreg[i]; 
      }
      b=s1/s2;
      xparam[0]=a; xparam[1]=b; xparam[2]=xcp;
      Q=0;
      for(i=0;i<kx;i++) {
         Q+=(yreg[i]-a-b*(xreg[i]-xcp))*(yreg[i]-a-b*(xreg[i]-xcp))*wreg[i];
      }
        xparam[3]=Q/(kx-2);
        xparam[4]=1/znam; //D{a}
        xparam[5]=1/s2; //D{b}
   }
//*****************************************************
 function CritByName(crit,a,kx,n,m) {
    let i,j,k,msmal,s,km,r1,r2,r,hcrit,er;
    let vorder=[];

    msmal=Math.min(m[1],m[2]);

   if (crit=="Kruskal") {
    s=0;km=0;
        for (j=1;j<=kx;j++) {
          r=0;
          for (i=1;i<=n;i++) {
            for (k=1;k<=m[j];k++) if (a[i]==k+km) r=r+i;
          }
            s=s+r*r/m[j];km=km+m[j];
        }
        hcrit=12*s/(n*(n+1))-3*(n+1);
     }
  if (crit=="Leman") {
   r1=0;r2=0;
     for (k=1;k<=m[1];k++) {
       for (i=1;i<=n;i++) if (a[i]==k) r1=r1+(i-k)*(i-k);
      }
     for (k=1;k<=m[2];k++) {
       for (i=1;i<=n;i++) if (a[i]==k+m[1]) r2=r2+(i-k)*(i-k);
      }
      hcrit=(r1*m[1]+r2*m[2]+m[1]*m[2]/6)/(m[1]*m[1]*m[2]*m[2])-2/3;
   }
  if (crit=="Mood") {
      r=0;
      for (i=1;i<=n;i++) {
          for (j=1;j<=msmal;j++) if(a[i]==j) r=r+(i-0.5*(n+1))**2;
         }
         hcrit=r;
  }
  if (crit=="Wilcoxon") {
        r=0;
         for (i=1;i<=n;i++) {
          for (j=1;j<=msmal;j++)    if (a[i]==j) r+=i;
         }
         hcrit=r;
  }
  if (crit=="Ansari") {
      r=0;
      for (i=1;i<=n;i++) {
          for (j=1;j<=msmal;j++) if (a[i]==j) r+=((n+1)/2-Math.abs(i-(n+1)/2));
        }
         hcrit=r;
  }
  if(crit=="Fisher") {
      r=0;
       for (i=1;i<=n;i++) { 
           vorder=ordern(n,i,i);
           er=vorder[0];
         for (j=1;j<=msmal;j++) if (a[i]==j) r=r+er;
        }
         hcrit=r; 
  }
 if(crit=="VanDerVarden") {
      r=0;
       for (i=1;i<=n;i++) { 
           er=invnormaldistribution(i/(m[1]+m[2]+1));
         for (j=1;j<=msmal;j++) if (a[i]==j) r=r+er;
        }
         hcrit=r; 
  }
 if(crit=="Capon") {
      r=0;
       for (i=1;i<=n;i++) { 
        vorder=ordern(n,i,i);
        er=vorder[0];
        cr=vorder[1];
        er=cr+er*er;
         for (j=1;j<=msmal;j++) if (a[i]==j) r=r+er;
        }
         hcrit=r; 
  }
if(crit=="Klotz") {
      r=0;
       for (i=1;i<=n;i++) { 
           er=invnormaldistribution(i/(m[1]+m[2]+1));
         for (j=1;j<=msmal;j++) if (a[i]==j) r=r+er*er;
        }
         hcrit=r; 
  }
   return hcrit;
}
//*********************************************
 function omega2(kx,m2,crit,xz,w,p) {
        var m3=[],h=[],a=[],m=[];
         n=0;
        for(i=0;i<kx;i++) n+=m2[i];
       kk=0;
      for(j=1;j<=kx;j++) {
        m[j]=m2[j-1];
        for(j1=1;j1<=m[j];j1++) m3[j1+kk]=j;
        kk+=m[j];
        }
//**********************************************************
      kk=0;
      for (i1=1;i1<=n;i1++) {
      a[1]=i1;

     if(n<2) return;

 M2: for (i2=1;i2<=n;i2++) { 
       a[2]=i2;
      for (j=1;j<2;j++)  if( (i2<=a[j]) && (m3[i2]==m3[a[j]])) continue M2;
       if(n==2) {h[kk]=CritByName(crit,a,kx,n,m);kk++;break;}

     if(n<3) return;

 M3: for (i3=1;i3<=n;i3++) { 
      a[3]=i3;
      for (j=1;j<3;j++)   if((i3<=a[j]) && (m3[i3]==m3[a[j]])) continue M3;
       if(n==3) {h[kk]= CritByName(crit,a,kx,n,m);kk++;break;}

  if(n<4) return;

 M4: for (i4=1;i4<=n;i4++) { 
       a[4]=i4;
      for (j=1;j<4;j++)    if((i4<=a[j])&& (m3[i4]==m3[a[j]])) continue M4;
       if(n==4) {h[kk]= CritByName(crit,a,kx,n,m);kk++;break;}

  if(n<5) return;

 M5: for (i5=1;i5<=n;i5++) { 
       a[5]=i5;
      for (j=1;j<5;j++)   if((i5<=a[j]) && (m3[i5]==m3[a[j]])) continue M5;
       if(n==5) {h[kk]= CritByName(crit,a,kx,n,m);kk++;break;}

       if(n<6) return;
M6: for (i6=1;i6<=n;i6++) { 
       a[6]=i6;
      for (j=1;j<6;j++)   if((i6<=a[j])&& (m3[i6]==m3[a[j]])) continue M6;
      if(n==6) {h[kk]= CritByName(crit,a,kx,n,m);kk++;break;}

       if(n<7) return;
M7: for (i7=1;i7<=n;i7++) { 
      a[7]=i7;
      for (j=1;j<7;j++)   if((i7<=a[j])&& (m3[i7]==m3[a[j]])) continue M7;
      if(n==7) {h[kk]= CritByName(crit,a,kx,n,m);kk++;break;}

       if(n<8) return;
M8: for (i8=1;i8<=n;i8++) { 
      a[8]=i8;
      for (j=1;j<8;j++)    if((i8<=a[j])&& (m3[i8]==m3[a[j]])) continue M8;
      if(n==8) {h[kk]= CritByName(crit,a,kx,n,m);kk++;break;}

 if(n<9) return;

M9: for (i9=1;i9<=n;i9++) { 
      a[9]=i9;
      for (j=1;j<9;j++)   if((i9<=a[j])&& (m3[i9]==m3[a[j]])) continue M9;
      if(n==9) {h[kk]= CritByName(crit,a,kx,n,m);kk++;break;}

   if(n<10) return;

M10: for (i10=1;i10<=n;i10++) { 
      a[10]=i10;
      for (j=1;j<10;j++)   if((i10<=a[j]) && (m3[i10]==m3[a[j]])) continue M10;
      if(n==10) {h[kk]= CritByName(crit,a,kx,n,m);kk++;break;}

 if(n<11) return; 

M11: for (i11=1;i11<=n;i11++) {
      a[11]=i11;
      for (j=1;j<11;j++)   if((i11<=a[j]) && (m3[i11]==m3[a[j]])) continue M11;
      if(n==11) {h[kk]= CritByName(crit,a,kx,n,m);kk++;break;}

  if(n<12) return;

 M12: for (i12=1;i12<=n;i12++) { 
       a[12]=i12;
      for (j=1;j<12;j++)   if((i12<=a[j]) && (m3[i12]==m3[a[j]])) continue M12;
      if(n==12) {h[kk]= CritByName(crit,a,kx,n,m);kk++;break;}

 if(n<13) return;

 M13: for (i13=1;i13<=n;i13++) { 
      a[13]=i13;
      for (j=1;j<13;j++)  if((i13<=a[j]) && (m3[i13]==m3[a[j]])) continue M13;
      if(n==13) {h[kk]= CritByName(crit,a,kx,n,m);kk++;break;}

 if(n<14) return;

 M14: for (i14=1;i14<=n;i14++) { 
      
       a[14]=i14;
      for (j=1;j<14;j++)   if((i14<=a[j])&& (m3[i14]==m3[a[j]])) continue M14;
      if(n==14) {h[kk]= CritByName(crit,a,kx,n,m);kk++;break;}

  if(n<15) return;

 M15: for (i15=1;i15<=n;i15++) { 
    
       a[15]=i15;
      for (j=1;j<15;j++)   if((i15<=a[j]) && (m3[i15]==m3[a[j]])) continue M15;
      if(n==15) {h[kk]= CritByName(crit,a,kx,n,m);kk++;break;}

     if(n<16) return;

M16: for (i16=1;i16<=n;i16++) { 
       a[16]=i16;
       for (j=1;j<16;j++)   if((i16<=a[j])&& (m3[i16]==m3[a[j]])) continue M16;
       if(n==16) {h[kk]= CritByName(crit,a,kx,n,m);kk++;break;}

      if(n<17) return;

M17: for (i17=1;i17<=n;i17++) { 
       a[17]=i17;
       for (j=1;j<17;j++)   if((i17<=a[j])&& (m3[i17]==m3[a[j]])) continue M17;
       if(n==17) {h[kk]= CritByName(crit,a,kx,n,m);kk++;break;}

        if(n<18) return;

M18: for (i18=1;i18<=n;i18++) { 
        a[18]=i18;
      for (j=1;j<18;j++) if((i18<=a[j])&& (m3[i18]==m3[a[j]])) continue M18;
      if(n==18) {h[kk]= CritByName(crit,a,kx,n,m);kk++;break;}

      if(n<19) return;

M19: for (i19=1;i19<=n;i19++) { 
       a[19]=i19;
      for (j=1;j<19;j++) if((i19<=a[j])&& (m3[i19]==m3[a[j]])) continue M19;
      if(n==19) {h[kk]= CritByName(crit,a,kx,n,m);kk++;break;}

      if(n<20) return;

M20: for (i20=1;i20<=n;i20++) { 
       a[20]=i20;
      for (j=1;j<20;j++) if((i20<=a[j])&& (m3[i20]==m3[a[j]])) continue M20;
      if(n==20) {h[kk]= CritByName(crit,a,kx,n,m);kk++;break;}

      if(n<21) return;

M21: for (i21=1;i21<=n;i21++) { 
       a[21]=i21;
      for (j=1;j<21;j++) if((i21<=a[j])&& (m3[i21]==m3[a[j]])) continue M21;
      if(n==21) {h[kk]= CritByName(crit,a,kx,n,m);kk++;break;}
    }}}}}}}}}}}
    }}}}}}}}}} 
//********************************************************************
  h.sort(function(a,b){return a-b});
  nnew=0;nn=0; 
  for (i=0;i<kk;i++) {
  if(parseFloat(h[i]).toFixed(znaki)==parseFloat(h[i+1]).toFixed(znaki)) {
    nnew++;
   }
   else {
   w[nn]=nnew+1;xz[nn]=h[i];nnew=0;nn++;
   }
   }
     s=0;kxt=0;
    for (i=0;i<nn;i++) {
      s+=w[i]/kk;p[i]=s;
    }
    return kk;
 } // end function

//***********Точное распределение критерия Ансари-Брэдли*********
function data_ansari(m,n,astat,w,pw) {
  let min_val,max_val,i,ir,nrows,sum,a0,s;
  
  min_val=m;
  max_val=n;
 
  nrows=gscale(min_val,max_val,w);
  //nrows=P.length;

  a0=Math.floor((min_val+1)/2)*(1+Math.floor(min_val/2));

   sum=0;
   for(i=1;i<=nrows;i++) {
     astat[i]=a0+i-1;
     sum=sum+w[i];
    }
    
    s=0;
    for(i=1;i<=nrows;i++) {
       s=s+w[i];
       pw[i]=s/sum;
    }
    return nrows;
}
//**********************************************************************
  function gscale(test,other,pw) {

/*
  algorithm as 93 appl. statist. (1976) vol.25, no.1
  from the sizes of two samples the distribution of the
  ansari-bradley test for scale is generated in array pw.
*/
	let ai,one,i,symm,m,lres,mm1,nm1,nm2,ier,mnow,ks,j,z,ndo,l1out;
        let n,ln1,ln2,nc,l2out,n2b1,n2b2,ln3,kk;
        let fadd=[],fimply=[],a2=[],a3=[];

//pw.length=1+Math.floor((m*n)/2);  

          one = 1;
          m = Math.min(test, other);
          if(m < 0) return(0);
          n=Math.max(test,other);
          lres = 1 + Math.floor((m * n)/2);
          symm=false;
          z=(m+n)%2;
          if(z==0) symm=true;
          mm1 = m - 1;
//*****************************************************         
       if(m<=2) {
          if(mm1<0) {
            pw[1]=1;return(lres);
          }
          
          if(mm1==0) ln1=start1(n, pw);
          if(mm1>0)  ln1=start2(n, pw);
            if(symm || (other > test)) return(lres);
            j = lres;
            ndo =Math.floor(lres/2);
          for(i=1;i<=ndo;i++) {
              ai = pw[i];
              pw[i]=pw[j];
              pw[j]=ai;
              j = j - 1;
          }
            return(lres);
       }
//***********************************************************
          nm1=n-1;nm2=n-2;mnow=3;nc=3;ier=0;
//**************************************************************
   while(true) {
    if(ier==0) {
        if((n % 2)!=1) {
          n2b1 = 3;
          n2b2 = 2;
          ln1=start2(n,pw);
          ln3=start2(nm2,a3);
          ln2=start1(nm1,a2);
//***********************************************************************          
          fadd=frqadd(a2,a3,ln2,l2out,ln3,n2b2);
          l2out=fadd[1];n2b2=fadd[2];
          ln2 = ln2 + nm1;
          fimply=imply(a2,a3,l2out,ln2,j,nc);
          ln2=fimply[1];j=fimply[2];nc=fimply[3];
//***************************************************************************
          nc = nc + 1;
          if(mnow==m) break;
          mnow = mnow + 1;
         } else {
          n2b1 = 2;
          n2b2 = 3;
          ln1=start1(n,pw);
          ln2=start2(nm1,a2);
        }
     }
//*******************************************************************
          fadd=frqadd(pw,a2,ln1, l1out,ln2, n2b1);
          l1out=fadd[1];n2b1=fadd[2];
          ln1 = ln1 + n;
          fimply=imply(pw,a3,l1out, ln1,ln3, nc);
          ln1=fimply[1];ln3=fimply[2];nc=fimply[3];
          nc = nc + 1;
          if(mnow==m) break;
          mnow = mnow + 1;
          fadd=frqadd(a2,a3, ln2, l2out,ln3, n2b2);
          l2out=fadd[1];n2b2=fadd[2];
          ln2 = ln2 + nm1;
          fimply=imply(a2,a3,l2out, ln2,j, nc);
          ln2=fimply[1];j=fimply[2];nc=fimply[3];
//***************************************************************************
          nc = nc + 1;
          if(mnow==m) break;
          mnow = mnow + 1;
          ier = 1;
  }
//********************************************************************
          if(symm) return(lres);
          ks = Math.floor((m + 3)/ 2);
          j = 1;
          for(i = ks;i<=lres;i++) {
           if(i>ln1) {
            pw[i]=a2[j];
           } else {
            pw[i]=pw[i]+a2[j];
          }
             j = j + 1;
          }
          if(other < test) return(lres);
          j = lres;
          ndo = Math.floor(lres/2);
          for(i=1;i<=ndo;i++) {
              ai = pw[i];
              pw[i]=pw[j];
              pw[j]=ai;
              j = j - 1;
          }
 return(lres);
	}
//********************************************************************
function start1(n,f) {
/*
	  algorithm as 93.1 appl. statist. (1976) vol.25, no.1
	  generates a 1,n ansari-bradley distribution in f.
*/
	let i,lout;
	lout=Math.floor(1+n/2);
	for(i=1;i<=lout;i++) f[i]=2;
	if ((n % 2)==0) f[lout]=1;
        return(lout);
  }
//************************************************************
  function start2(n,f) {

/*
	  algorithm as 93.2 appl. statist. (1976) vol.25, no.1
	  generates a 2,n ansari-bradley distribution in f.
*/
	let one,two,three,four,i,j,a,b,lt1,ndo,nu,lout;

	one=1;two=2;three=3;four=4;
	nu=n-n % 2;
	j=nu+1;lout=j;lt1=lout+1;
        ndo=Math.floor(lt1/2);
        a=one;b=three;
      for (i=1;i<=ndo;i++) {
	f[i]=a;f[j]=a;j=j-1;a=a+b;b=four-b;
      }
	if(nu==n) return(lout);
	nu=ndo+1;
	for(i=nu;i<=lout;i++) f[i]=f[i]+two;
	f[lt1]=two;lout=lt1;
	return(lout);
  }
//************************************************************
  function frqadd(f1,f2,l1in,l1out,l2,nstart) {
         

/*
	  algorithm as 93.3 appl. statist. (1976) vol.25, no.1
	  array f1 has twice the contents of array f2 added into it
	  starting with elements nstart and 1 in f1 and f2 respectively.
*/

	let i1,i2,nxt,fadd=[];

        i2=1;
          for(i1=nstart;i1<=l1in;i1++) {
              f1[i1]=f1[i1]+2*f2[i2];
              i2=i2+1;
          }
          nxt=l1in+1;
          l1out=l2+nstart-1;
          for(i1=nxt;i1<=l1out;i1++) {
              f1[i1]=2*f2[i2];
              i2=i2+1;
          }
          nstart=nstart+1;
          fadd[1]=l1out;fadd[2]=nstart;
          return(fadd);
 }
//*******************************************************************
  function imply(f1,f2,l1in,l1out,l2,noff) {
/*
	  algorithm as 93.4 appl. statist. (1976) vol.25, no.1
	  given l1in elements of an array f1, a symmetrical
	  array f2 is derived and added onto f1, leaving the
	  first noff elements of f1 unchanged and giving a
	  symmetrical result of l1out elements in f1.
*/
	let sum,diff,i2,i1,j2,j1,j2min,ndo,fimply=[];

    i2=1-noff;j1=l1out;j2=l1out-noff;l2=j2;
    j2min=Math.floor((j2 + 1)/2);
    ndo=Math.floor((l1out+1)/2);

          for(i1=1;i1<=ndo;i1++) {
              if(i2>0) {
                  sum=f1[i1]+f2[i2];
                  f1[i1]=sum;
              }
              else {
                  sum=f1[i1];
              }
              i2=i2+1;
              if(j2>=j2min) {
                  if(j1<=l1in) {
                      diff=sum-f1[j1];
                   }
                  else {
                     diff=sum;
                   }

              f2[i1]=diff;f2[j2]=diff;j2=j2-1;
             }
             f1[j1]=sum;j1=j1-1;
          }
            fimply[1]=l1out;fimply[2]=l2;fimply[3]=noff;
  return(fimply);
}
//***********Точное распределение критерия Уилкоксона******************************************
function wilc_Rec(m,n,range,w,pw) {

/*
  AS 62 generates the frequencies for the Mann-Whitney U-statistic.
  Users are much more likely to need the distribution function.
  Code to return the distribution function has been added at the end
  of AS 62 by Alan Miller
*/


 let work=[],minmn,maxmn,inx,i;
 let n1,l,k,j,sum,mn,s;
  
  mn=m*n+1;
  maxmn=Math.max(m,n);
  minmn=Math.min(m,n);
  n1=maxmn+1;
 
     for(i=1;i<=n1;i++) w[i]=1;
      n1=n1+1;
      for(i=n1;i<=mn;i++) w[i]=0;
  
     work[1]=0;inx=maxmn;

    for(i=2;i<=minmn;i++) {
        work[i]=0;inx=inx+maxmn;n1=inx+2;l=1+inx/2;k=i;

       for(j=1;j<=l;j++) {
          k=k+1;n1=n1-1;
          sum=w[j]+work[j];
          w[j]=sum;work[k]=sum-w[n1];w[n1]=sum;
       }
    }
  
  
  //w[i] - Frequencies
   //for(i=1;i<=mn;i++) pw[i]=w[i];
 
sum=0;
for(i=1;i<=mn;i++) sum=sum+w[i];

    s=0;
     for(i=1;i<=mn;i++) {
        s=s+w[i];
        range[i]=i-1+minmn*(minmn+1)/2; //Wilcoxon
        //range[i]=i-1; //Mann-Whitney
        pw[i]=s/sum;
     }
 
    return(mn);
 }
//*****************Точное распределение критерия Уилкоксона по рекурентной формуле*********************
function wilc_Rec_1(m,n,P) { 
       var W=[];
        w=m*n;s=0;

       for(i=1;i<=m;i++) {
           W[i]=[];
           for(j=1;j<=n;j++) {
             W[i][j]=[];
             for(k=0;k<=w;k++)  W[i][j][k]=0;
              }}

       for(k=0;k<=w;k++) {
         for(i=1;i<=m;i++) {
             for(j=1;j<=n;j++) {
                 h1=test_Rec(i,j-1,k-i);
                 if(h1==-1)  h1=W[i][j-1][k-i];
                 h2=test_Rec(i-1,j,k);
                 if(h2==-1)  h2=W[i-1][j][k];
                 W[i][j][k]=(h1*j+h2*i)/(i+j);
           }}
               s+=W[m][n][k];
              P[k]=s;
        }
 }
//*****************************************************
function test_Rec(a,b,c) {
 if(c<0) return 0;
 if(a==0 || b==0) {
   if(c==0) return 1;
   return 0;
 }
 return -1;
}
//*****Точное распределение критерия знаковых рангов Уилкоксона************************************************
 function wsigne_Rec(n,P) { 
   var a=[],aa=[];
   k1=3;k3=4;
  for(i=1;i<=4;i++) {aa[i]=1;a[i]=1;}
  for(j=3;j<=n;j++) {
      k1+=1;k2=k3; k3+=j;
      for(i=k1;i<=k3;i++) {
        s1=0;
        if(i<=k2) s1=a[i];
        a[i]=aa[i-k1+1]+s1;
      }
     for(i=k1;i<=k3;i++) aa[i]=a[i];
  }
  s=0;s1=Math.pow(2,n);
 for(i=1;i<=k3;i++) {
    s+=a[i]/s1;aa[i]=i-1;P[i]=s;
   }
}
//**Точное распределение критерия знаковых рангов Уилкоксона по рекурентной формуле(в книгах ошибка, смотри test_signe)***
function signe_Rec(n,P) { 
      var i,k,s;
       var W=[];
       var w=n*(n+1)/2+1;
      var  z=Math.pow(2,n);

     for(i=1;i<=n;i++) {
           W[i]=[];
             for(k=0;k<=w;k++)  W[i][k]=0;
             }
       s=0;
     
        for(k=0;k<=w;k++) {
             for(i=1;i<=n;i++) {
                h1=test_signe(i-1,k);
                if(h1==-1)  h1=W[i-1][k];
                h2=test_signe(i-1,k-i);
                if(h2==-1)  h2=W[i-1][k-i];
                W[i][k]=h1+h2; 
              }
               s+=W[n][k];
               P[k]=s/z;
        }
 }
//*****************************************************
function test_signe(a,b) {
   if(b<0) return 0;
   if(a==0) {
     if(b==0) return 1;
     return 0;
   }
 return -1;
}
//************************************************************
function make_order(ts,nx1,nx,kmls,x,bmls,dmls) {
   var V= [],xmls=[],yud=[],vorder=[];
      for (i=0;i<nx1;i++) {
        V[i]= [];xmls[i]=[]; yud[i]=[];yud[i][0]=x[i];
        for (j=0;j<kmls;j++) xmls[i][j]=1;
        for (j=0;j<nx1;j++) V[i][j]=0;
       }
       for (i=0;i<nx1;i++) {
          for (j=i;j<nx1;j++) {
           if (ts=="normal") vorder=ordern(nx, i+1, j+1);
           if (ts=="weibull") vorder=orderw(nx, i+1, j+1);
           xmls[i][1]=vorder[0];V[j][i]=vorder[1];V[i][j]=vorder[1];
       }
    }
       mleastsquare(xmls,yud,kmls,V,bmls,dmls);
}
//************************************************************************
function regress_matrix(kx,kmls,wreg,yreg,xreg,yr,Q,emls,dmls) {
      var V=[],y=[],x=[];
       for(i=0;i<kx;i++) {
        V[i]=[];
          for(j=0;j<kx;j++) V[i][j]=1;
      }
      for(i=0;i<kx;i++) V[i][i]=wreg[i];
       for(i=0;i<kx;i++) {
        y[i]=[];y[i][0]=yreg[i];
       }
   for(i=0;i<kx;i++) {
         x[i]=[];
           for(j=0;j<kmls;j++) x[i][j]=Math.pow(xreg[i],j);
   }
//**********************************************************
     mleastsquare(x,y,kmls,V,emls,dmls);
    s2=0;znam=0;
 for(i=0;i<kx;i++) {
      znam=znam+1/wreg[i];s1=0;
      for (j=0;j<kmls;j++)  s1+=emls[j]*x[i][j];
       yr[i]=s1;
      s2+=(yreg[i]-s1)*(yreg[i]-s1)/wreg[i];
   }
      Q[0]=s2/znam;
 }

//#######################Nalder-Mead################################################
function simpl(xsimpl,step,eps,lim,nx,fn) {
    let x1=[], sum=[], xnx, step1, step2, difer, dif;
    let  alfa, beta, gama, suml, sums, sumh, sum2;
    let istep, k1, k2, k3, k4, vn, i, j, k, l, l1, l2, l3, kount;
    let ik, inx, index;

    k = 0; q = 0;
    let dop=[];

    for (i=0;i<500;i++) x1[i]=0;
    for (i=0;i<25; i++) {dop[i]=0;sum[i]=0;}
        
    alfa = 1.0; beta = 0.45; gama = 2.8; istep = 0; ier = 0; k1 = nx + 1;
    k2 = nx + 2; k3 = nx + 3; k4 = nx + 4;
    vn = nx; xnx = 1. / (1.0*vn);
    step1=step/ (1.0*vn * Math.sqrt(2.)) * (Math.sqrt(vn + 1.) + vn - 1.);
    step2 = step / (1.0*vn * Math.sqrt(2.)) * (Math.sqrt(vn + 1.) - 1.);
    for (i = 2; i <= k1; i++) {
        l = (i - 1) * nx; l1 = l + i - 1;
        for (j = 1; j <= nx; j++) {
            l2 = l + j; x1[l2] = xsimpl[j-1] + step2;
        }
        x1[l1] = xsimpl[i-2]+step1;
    }
    for (j=1;j <= nx;j++) x1[j]=xsimpl[j-1];
    for (i=1;i<=k1;i++) {
        l=(i-1)*nx+1;
        for (ik=1;ik<=nx;ik++) {
            dop[ik-1]=x1[l];l++;
        }
        sum[i] = (fn(dop));
    }

  flag="a28";
while(1>0) {
  switch(flag) {
 case "a28":
    sumh = sum[1]; index = 1;
    for (i = 2;i<=k1;i++) {
        if (sum[i] >sumh) {
            sumh=sum[i];index=i;
        }
    }
    suml = sum[1]; kount = 1;
    for (i = 2; i <= k1; i++) {
        if (suml > sum[i]) {
            suml = sum[i]; kount = i;
        }
    }
    istep++; difer = 0.;
    for (i = 1; i <= k1; i++) {
        dif = sum[i] - sum[kount]; difer = difer + dif * dif;
    }
    difer=xnx*Math.sqrt(difer);
    if (suml<=eps) {flag="a32";break;}
    if (suml>eps) {flag="a31";break;}
case "a32":
    if (difer<= eps) {flag="a30";break;}
    if (difer>eps) {flag="a26";break;}
case "a31":
    if ((difer / Math.abs(suml)) <= eps) {flag="a30";break;}
    if ((difer / Math.abs(suml)) > eps) {flag="a26";break;}
case "a30":
    difer = 0.; l1 = (kount - 1) * nx;
    for (i = 1; i <= k1; i++) {
        if (kount != i) {
            l = (i - 1) * nx;
            for (j = 1; j <= nx; j++) {
                l2 = l + j; l3 = l1 + j;
                if (x1[l3] <= eps) dif = x1[l2] - x1[l3];
                if (x1[l3] > eps)  dif = (x1[l2] - x1[l3]) / x1[l3];
                difer = difer + dif * dif;
            }
        }
    }
    difer = xnx *Math.sqrt(difer);
    if (difer <= eps) {flag="a23";break;}
case "a26":
    if (istep >= lim) {flag="a38";break;}
    if (k == 1) {flag="a17";break;}
    for (j = 1; j <= nx; j++) {
        sum2 = 0.;
        for (i = 1; i <= k1; i++) {
            l = (i - 1) * nx + j; sum2 = sum2 + x1[l];
        }
        l2 = (index - 1) * nx + j; l1 = k1 * nx + j; x1[l1] = (sum2 - x1[l2]) * xnx;
        l3 = l1 + nx; x1[l3] = (1. + alfa) * x1[l1] - alfa * x1[l2];
    }
    inx = k2 * nx + 1;
    for (ik = 1; ik <= nx; ik++) {
        dop[ik-1] = x1[inx]; inx++;
    }
    sum[k3] =(fn(dop));
    if (sum[k3] < suml) {flag="a11";break;}
    sums = suml;

//??????????????????????????????????
    for (i = 1; i <= k1; i++) {
        if (index == i) continue;
        if (sum[i] <= sums)  sums = sum[i];
    }
    if (sum[k3] > sums) {flag="a13";break;}
    flag="a14";break;
case "a11":
    for (j = 1; j <= nx; j++) {
        l = k1 * nx + j; l1 = l + nx; l2 = l1 + nx; x1[l2] = (1. - gama) * x1[l] + gama * x1[l1];
    }
    inx = k3 * nx + 1;
    for (ik = 1; ik <= nx; ik++) {
        dop[ik-1] = x1[inx]; inx++;
    }
    sum[k4] =(fn(dop));
    if (sum[k4] < suml) {flag="a16";break;}
    flag="a14";break;
case "a13":
    if (sum[k3] > sumh) {flag="a17";break;}
     k = 1; flag="a14";break;
case "a17":
    for (j = 1; j <= nx; j++) {
        l = k3 * nx + j; l1 = (index - 1) * nx + j; l2 = l - nx - nx;
        x1[l] = beta * x1[l1] + (1. - beta) * x1[l2];
    }
    k = 0; inx = k3 * nx + 1;
    for (ik = 1; ik <= nx; ik++) {
        dop[ik-1] = x1[inx]; inx++;
    }
    sum[k4] =(fn(dop));
    if (sumh > sum[k4]) {flag="a16";break;}
    for (j = 1; j <= nx; j++) {
        l1 = (kount - 1) * nx + j;
        for (i = 1; i <= k1; i++) {

//??????????????????????????????
         if (i != kount) l = (i - 1) * nx + j;
          x1[l] = 0.5 * (x1[l] + x1[l1]);
        }
    }
    for (i = 1; i <= k1; i++) {
        if (i != kount) {
        l = (i - 1) * nx + 1;
         for (ik = 1; ik <= nx; ik++) {
             dop[ik-1] = x1[l]; l++;
         }
        sum[i]=(fn(dop));
       }
    }
    flag="a28";break;
case "a16":
    for (j = 1; j <= nx; j++) {
        l = (index - 1) * nx + j; l1 = k3 * nx + j; x1[l] = x1[l1];
        sum[index] = sum[k4];
    }
    flag="a28";break;
case "a14":
    for (j = 1; j <= nx; j++) {
        l = (index - 1) * nx + j; l1 = k2 * nx + j; x1[l] = x1[l1];
        sum[index] = sum[k3];
    }
    flag="a28";break;
case "a38":
    ier=1;
case "a23":
    for (j = 1; j <= nx; j++) {
        l = (kount - 1) * nx + j; xsimpl[j-1] = x1[l];
    }
    sum[1]=suml;lim=istep;q = sum[1];
    return lim;

} //end switch
}  //end while


}

//******************************************************************
   function simplex(p,y,pr,n,eps,itmax,iregim) {
    var prr=[];
    var pbar=[];  
     one=1;two=2;mpts=n+1;alphas=1;betas=0.5;gammas=2;zero=0;one=1;half=0.5;two=2;mpts=n+1;
     iflag=0;iter=0;
    if(iregim==0) {
      s1=(Math.sqrt(one+n)+n-one)/(n*Math.sqrt(two));
      s2=(Math.sqrt(one+n)-one)/(n*Math.sqrt(two));
    for (i=2;i<=mpts;i++){
     for (j=1;j<=n;j++) {
      p[i][j]=s2+p[1][j];
      L = i-1;
      p[i][L]=s1+p[1][j];
     }
    }
  }
 for (i=1;i<=mpts;i++){
  for (j=1;j<=n;j++)  pr[j]=p[i][j];
   y[i]=funx(pr);
 }
//***************************************************
  flag="50";
 while(1>0) {
 switch(flag) {
  case "50":
    ilo=1;
    if(y[1]>y[2]) {
      ihi=1;inhi=2;
      }
      else {
      ihi=2;inhi=1;
    }
    for (i=1;i<=mpts;i++) {
      if(y[i]<y[ilo]) ilo=i;
     if(y[i]>y[ihi]) {
      inhi=ihi;ihi=i;
      }
      else {
      if(y[i]>y[inhi]) {
        if(i != ihi) inhi=i;
      }
    }
   }
   rtol=2*Math.abs(y[ihi]-y[ilo])/(Math.abs(y[ihi])+Math.abs(y[ilo]));
   if(rtol<Math.abs(eps)) {
    iflag=1;return iter;
   }
   if(iter==itmax) return iter;
   iter++;
    for(j=1;j<=n;j++) pbar[j]=zero;
      for(i=1;i<=mpts;i++) {
      if(i!=ihi) {
       for(j=1;j<=n;j++) pbar[j]=pbar[j]+p[i][j];
      }
     }
    for(j=1;j<=n;j++) {
       pbar[j]=pbar[j]/n;pr[j]=(one+alphas)*pbar[j]-alphas*p[ihi][j];
    }
       ypr=funx(pr);
   if(ypr<=y[ilo]){
        for(j=1;j<=n;j++)  prr[j]=gammas*pr[j]+(one-gammas)*pbar[j];
        yprr=funx(prr);
      if(yprr>y[ilo]){
        for(j=1;j<=n;j++)  p[ihi][j]=prr[j];
       y[ihi]=yprr;
       flag="50";break;
      }
     for(j=1;j<=n;j++) p[ihi][j]=pr[j];
      y[ihi]=ypr;
        flag="50";break;
   }   
   if(ypr<y[inhi]){
        for(j=1;j<=n;j++)  p[ihi][j]=pr[j];
        y[ihi]=ypr;
         flag="50";break;
    }
   if(ypr<y[ihi]){
       for(j=1;j<=n;j++)  p[ihi][j]=pr[j];
       y[ihi]=ypr;
  }
    for(j=1;j<=n;j++)  prr[j]=betas*p[ihi][j]+(one-betas)*pbar[j];
    yprr=funx(prr);
      if(yprr<y[ihi]){
       for(j=1;j<=n;j++)  p[ihi][j]=prr[j];
       y[ihi]=yprr;
         flag="50";break;
      }
      for(i=1;i<=mpts;i++) {
      if(i != ilo) {
       for(j=1;j<=n;j++) {
        pr[j]=half*(p[i][j]+p[ilo][j]);
        p[i][j]=pr[j];
      }
      y[i]=funx(pr);
      }
      }
        flag="50";break;
      }//end switch
      } //while
 return iter;
 }
//*************************************************************************************
function mlses_cum(ts,nn,fcum,x,bmls,dmls) {
var E= [],V= [];
l=fcum.length;
for (i=0;i<l;i++) {
  V[i]= []; 
   for (j=0;j<l;j++) {
    V[i][j]=0;
 }
}

 kmls=2;
 for (i=0;i<kmls;i++) {
 dmls[i]=[];  
for (j=0;j<kmls;j++) dmls[i][j]=0;
  }

for (i=0;i<l;i++) {
   for (j=i;j<l;j++) {
   if (ts=="normal") vorder=ordern_cum(nn, fcum[i],fcum[j]);
   if (ts=="weibull") vorder=orderw_cum(nn, fcum[i], fcum[j]);
    E[i]=vorder[0];
    V[j][i]=parseFloat(vorder[1]);
    V[i][j]=parseFloat(vorder[1]);
   }
   }

V=InverseMatrix(V);
s1=0;s2=0;s3=0;s4=0;s5=0;s6=0;
for (i=0;i<l;i++) {
ei=E[i];ai=0;ai1=0;bi=0;
for (j=0;j<l;j++) {
Vij=V[j][i];
ai=ai+E[j]*Vij;bi=bi+Vij*x[j];ai1=ai1+Vij;
}
s1=s1+ai;s2=s2+bi;s3=s3+ei*bi;s4=s4+ei*ai;s5=s5+ai1;s6=s6+ai1*ei;
}
d=s5*s4-s1*s1;
bmls[0]=-(s3*s1-s4*s2)/d;
bmls[1]=(s3*s5-s6*s2)/d;
dmls[0][0]=s4/d;dmls[1][1]=s5/d;
dmls[0][1]=-s1/d;dmls[1][0]=-s1/d;
}
//***********Вейбулловкие порядковые статистики*fcum*************************************

function orderw_cum(n,pr,ps) {
 qr=1-pr; qs=1-ps;
xr=Math.log(Math.log(1/(1-pr)));
xs=Math.log(Math.log(1/(1-ps)));
xr1=1/(Math.log(1/(1-pr))*(1-pr));
xr2=xr1*(1/(1-pr)-xr1);
xr3=xr2*xr2/xr1+xr1*(1/Math.pow(1-pr,2)-xr2);
xr4=(3*xr1*xr2*xr3-2*Math.pow(xr2,3))/Math.pow(xr1,2)+xr1*(2/Math.pow(1-pr,3)-xr3);
xr55=(-12*xr1*Math.pow(xr2,2)*xr3+3*Math.pow(xr1,2)*Math.pow(xr3,2)+4*Math.pow(xr1,2)*xr2*xr4+6*Math.pow(xr2,4));
xr5=xr55/Math.pow(xr1,3)+xr1*(6/Math.pow(1-pr,4)-xr4);
a1=-12*Math.pow(xr2,3)*xr3-12*xr1*(2*xr2*Math.pow(xr3,2)+Math.pow(xr2,2)*xr4);
b1=6*xr1*xr2*Math.pow(xr3,2)+6*Math.pow(xr1,2)*xr3*xr4;
c1=8*xr1*Math.pow(xr2,2)*xr4+4*Math.pow(xr1,2)*(xr3*xr4+xr2*xr5);
d1=24*Math.pow(xr2,3)*xr3;
xr6=(Math.pow(xr1,3)*(a1+b1+c1+d1)-3*Math.pow(xr1,2)*xr2*xr55)/Math.pow(xr1,6)+xr2*(6/Math.pow(1-pr,4)-xr4)+xr1*(24/Math.pow(1-pr,5)-xr5);
xs1=1/(Math.log(1/(1-ps))*(1-ps));
xs2=xs1*(1/(1-ps)-xs1);
xs3=Math.pow(xs2,2)/xs1+xs1*(1/Math.pow(1-ps,2)-xs2);
xs4=(3*xs1*xs2*xs3-2*Math.pow(xs2,3))/Math.pow(xs1,2)+xs1*(2/Math.pow(1-ps,3)-xs3);
xs5=(-12*xs1*Math.pow(xs2,2)*xs3+3*Math.pow(xs1,2)*Math.pow(xs3,2)+4*Math.pow(xs1,2)*xs2*xs4+6*Math.pow(xs2,4)) /Math.pow(xs1,3)+xs1*(6/Math.pow(1-ps,4)-xs4);
er=xr+pr*qr*xr2/(2*(n + 2))+pr*qr*((qr-pr)*xr3/3+pr*qr*xr4/8)/Math.pow(n+2,2)+pr*qr*(-(qr-pr)*xr3/3+(Math.pow(qr-pr,2)-pr*qr)*xr4/4+qr*pr*(qr-pr)*xr5/6+Math.pow(qr*pr,2)*xr6/48)/Math.pow(n + 2,3);
z1=(qr-pr)*xr2*xs1+(qs-ps)*xr1*xs2+pr*qr*xr3*xs1/2+ps*qs*xr1*xs3/2+pr*qs*xr2*xs2/2;
z1=z1*pr*qs/Math.pow(n+2,2);
z2=-(qr-pr)*xr2*xs1-(qs-ps)*xr1*xs2+(Math.pow(qr-pr,2)-pr*qr)*xr3*xs1;
z3=(Math.pow(qs-ps,2)-ps*qs)*xr1*xs3+(1.5*(qr-pr)*(qs-ps)+0.5*ps*qr-2*pr*qs)*xr2*xs2;
z4=(5/6)*pr*qr*(qr-pr)*xr4*xs1+(5/6)*ps*qs*(qs-ps)*xr1*xs4+(pr*qs*(qr-pr)+0.5*pr*qr*(qs-ps))*xr3*xs2;
z5=(pr*qs*(qs-ps)+0.5*ps*qs*(qr-pr))*xr2*xs3+(1/8)*Math.pow(pr*qr,2)*xr5*xs1+(1/8)*Math.pow(ps*qs,2)*xr1*xs5;
z6=0.25*Math.pow(pr,2)*qr*qs*xr4*xs2+0.25*pr*ps*Math.pow(qs,2)*xr2*xs4+(2*Math.pow(pr*qs,2)+3*pr*qr*ps*qs)*xr3*xs3/12;
z7=z2+z3+z4+z5+z6;
crs=z1+pr*qs*z7/Math.pow(n+2,3)+pr*qs*xr1*xs1/(n+2);
var vorder=[];
vorder[0]=er;vorder[1]=crs;
return vorder;
}

//********Нормальные порядковые статистики**fcum**************************************

function ordern_cum(n,pr,ps) {
 qr=1-pr;qs=1-ps;p=1;k=n; 
xr=invnormaldistribution(pr);
xs=invnormaldistribution(ps);
dr=Math.exp(-xr*xr/2.)/s2pi;
ds=Math.exp(-xs*xs/2.)/s2pi;
xr1=p/dr;
xr2=xr*Math.pow(p/dr,2);
xr3=(2*xr*xr+1)*Math.pow(p/dr,3);
xr4=(6*xr*xr*xr+7*xr)*Math.pow(p/dr,4);
xr5=(24*Math.pow(xr,4)+46*xr*xr+7)*Math.pow(p/dr,5);
xr6=(120*Math.pow(xr,5)+326*Math.pow(xr,3)+127*xr)*Math.pow(p/dr,6);
xs1=p/ds;
xs2=xs*Math.pow(p/ds,2);
xs3=(2*xs*xs + 1)*Math.pow(p/ds,3);
xs4=(6*Math.pow(xs,3)+7*xs)*Math.pow(p/ds,4);
xs5=(24*Math.pow(xs,4)+46*Math.pow(xs,2)+7)*Math.pow(p/ds,5);
xs6=(120*Math.pow(xs,5)+326*Math.pow(xs,3)+127*xs)*Math.pow(p/ds,6);
er=xr+pr*qr*xr2/(2*(k+2))+pr*qr*((qr-pr)*xr3/3+pr*qr*xr4/8)/Math.pow(k+2,2)+pr*qr*(-(qr-pr)*xr3/3+(Math.pow(qr-pr,2)-pr*qr)*xr4/4+qr*pr*(qr-pr)*xr5/6+Math.pow(qr*pr,2)*xr6/48)/Math.pow(k+2,3);

z1=(qr-pr)*xr2*xs1+(qs-ps)*xr1*xs2+pr*qr*xr3*xs1/2+ps*qs*xr1*xs3/2+pr*qs*xr2*xs2/2;
z1=z1*pr*qs/Math.pow(k+2,2);
z2=-(qr-pr)*xr2*xs1-(qs-ps)*xr1*xs2+(Math.pow(qr-pr,2)-pr*qr)*xr3*xs1;
z3=(Math.pow(qs-ps,2)-ps*qs)*xr1*xs3+(1.5*(qr-pr)*(qs-ps)+0.5*ps*qr-2*pr*qs)*xr2*xs2;
z4=(5/6)*pr*qr*(qr-pr)*xr4*xs1+(5/6)*ps*qs*(qs-ps)*xr1*xs4+(pr*qs*(qr-pr)+0.5*pr*qr*(qs-ps))*xr3*xs2;
z5=(pr*qs*(qs-ps)+0.5*ps*qs*(qr-pr))*xr2*xs3 + (1 /8)*Math.pow(pr*qr,2)*xr5*xs1+(1/8)*Math.pow(ps*qs,2)*xr1*xs5;
z6=0.25*Math.pow(pr,2)*qr*qs*xr4*xs2+0.25*pr*ps*Math.pow(qs,2)*xr2*xs4+(2*Math.pow(pr*qs,2)+3*pr*qr*ps*qs)*xr3*xs3/12;
z7=z2+z3+z4+z5+z6;
crs=z1+pr*qs*z7/Math.pow(k+2,3)+pr*qs*xr1*xs1/(k+2);
var vorder=[];
vorder[0]=er;vorder[1]=crs;
return vorder;
}


//######################Доверительные граница для квантиля#####################
function nctlimit(beta,f,t1,t2,t12,d){

/*
 input:
 beta-доверительная вероятность
 f-число степеней свободы
 t1,t22-диагональные элементы ковариационной матрица оценок параметров
 t12,t21-недиагональные элементы ковариационной матрица оценок параметров
 output:
 t- квантили (процентные точки) нецентрального распределения Стьюдента
*/
 
  let zb,p,B,C,m,s,b,c,z,z1,x1,x2,x3,x4,t;


 if(beta==0.8) {
  x1=0.5546390;x2=1.5745702;x3=2.3510674;x4=0.8906575;
 }

 if(beta==0.9) {
  x1=1.4713227;x2=0.9015831;x3=1.3143355;x4=0.8148939;
 }
 if(beta==0.95) {
   x1=3.7717452;x2=-1.4945123;x3=-0.5361564;x4=0.3030390;
 }
if(beta==0.975) {
  x1=2.7859170;x2=-5.2013471;x3=0.7889855;x4=0.8579880;
}
 if(beta==0.99) {
  x1=6.2277134;x2=-10.0225250;x3=-3.5862127;x4=2.1059900;
}

    zb=invnormaldistribution(beta);
    p=normaldistribution(d/Math.sqrt(f+1.));
    
    if(p>=0.2 && p<=0.5) {   
      b=zb+(Math.pow(zb,3)+zb)/(4.*f)+(5.*Math.pow(zb,5)+16.*Math.pow(zb,3)+3.*zb)/(96*f*f);
      c=1.+(2.*zb*zb+1.)/(4.*f)+zb*d/(4.*f)+(4.*Math.pow(zb,4)+12.*zb*zb+1.)/(32.*f*f)+(Math.pow(zb,3)+4.*zb)*d/(16.*f*f)-(zb*zb-1)*d*d/(24.*f*f)-zb*Math.pow(d,3)/(32.*f*f);
      t=b+d*c;
    }
    else {
     C=d**2-zb**2;   
     z1=1.-2.*t2*Math.sqrt(2./f)*Math.exp(lngamma((f+1.)/2.,0)-lngamma(f/2.,0));
     m=1.-z1*(x1+x2*p+x3/beta+x4/f)/f;
     z=Math.exp(lngamma(0.5*(f+1.),0)-lngamma(0.5*f,0));
     s=(1.-(2./f)*(z**2))/m;
     B=m*(1.-s*(zb**2));
     t=(d+zb*Math.sqrt(t1+s*C))/B;
  } 
  
   return t;
}

