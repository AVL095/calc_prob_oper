  let sm_points=[],s_type=[],s_color=[],s_marker=[],s_x=[],s_y=[],s_markersize=[],s_thick=[],s_dash=[];
  let ychart=["chart1","chart2","chart3","chart4","chart5"];
  let str=[];
  let ch_check;
  let rsimpl=[];
  let xdata=[];

  let xprob=[0.001,0.005,0.01,0.025,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.975,0.99,0.995,0.999];
  let stat_cv=[1.2,1.32,1.5,1.72,1.8,2.4,2.7,3.,3.6,3.9,4.2,4.44,4.5,4.68,4.8,5.04,5.1,5.28,5.4,5.7,6.];
  let stat_cvm,stat_cvstd,stat_pvx,stat_AX,stat_ALPHAX,stat_sy,stat_xu,stat_zp;
//*************************************************************************************
function All() {
  for(i=0;i<ychart.length;i++)  $(document.getElementById(ychart[i])).hide();
  setupMessage("Graphics","","left","D");
}


function D1() {
  if(ch_check==4 || ch_check==5) return;
  for(i=0;i<ychart.length;i++)  $(document.getElementById(ychart[i])).hide();
  Distribution(ychart[0]);
  $(document.getElementById(ychart[0])).show();

}

function D2() {
  if(ch_check==4 || ch_check==5) return;
  for(i=0;i<ychart.length;i++)  $(document.getElementById(ychart[i])).hide();
  Sy(ychart[1]);
 $(document.getElementById(ychart[1])).show();
}

function D3() {
  if(ch_check==4 || ch_check==5) return;
  for(i=0;i<ychart.length;i++)  $(document.getElementById(ychart[i])).hide();
  Regress(ychart[2]);
  $(document.getElementById(ychart[2])).show();
}

function D4() {
  if(ch_check==1 || ch_check==2 || ch_check==3 || ch_check==5) return;
  for(i=0;i<ychart.length;i++)  $(document.getElementById(ychart[i])).hide();
  Resurs(ychart[3]);
  $(document.getElementById(ychart[3])).show();
}

function D5() {
  if(ch_check==1 || ch_check==2 || ch_check==3 || ch_check==4) return;
  for(i=0;i<ychart.length;i++)  $(document.getElementById(ychart[i])).hide();
  Resurs_ap(ychart[4]);
  $(document.getElementById(ychart[4])).show();

}

//***********************************************************************************
function getGraphOptions(outg) {
 r='';  a='';
 r+='<button id="b1" class="btn btn-success" onClick=D1(ychart[0]); title="F1"><span class="glyphicon glyphicon-pencil"></span>Distribution</button> ';
 r+='<button id="b2" class="btn btn-success" onClick=D2(ychart[1]); title="F2"><span class="glyphicon glyphicon-pencil"></span>Sy=f(y)</button>'; 
 r+='<button id="b3" class="btn btn-success" onClick=D3(ychart[2]); title="F3"><span class="glyphicon glyphicon-pencil"></span>Regress</button> ';
 r+='<button id="b4" class="btn btn-success" onClick=D4(ychart[3]); title="F4"><span class="glyphicon glyphicon-pencil"></span>Resurs</button> ';
 r+='<button id="b5" class="btn btn-success" onClick=D5(ychart[4]); title="F5"><span class="glyphicon glyphicon-pencil"></span>Resurs_ap</button> ';
 for(i=0;i<ychart.length;i++)  a+='<div id='+'"'+ychart[i]+'"'+' style="width:90%;  height: 80%;  display: inline-block"></div>';
 outg[0]=r;outg[1]=a;
  }
//***********************************************************************
function TableSaveData() {
 var n=[];
 kx=parseInt(localStorage["AA"]);
  message="";message+="case 0:";message+="<br>";message+="var kx="+kx+";";message+="<br>";
  message+="var cx=[";
   for (i=0;i<kx;i++) {
     z=parseFloat(localStorage[str[i+2][1]]);
      message+=z;
      if(i<(kx-1)) message+=",";
  }
     message+="];<br>";
    message+="var n=[";
    for (i=0;i<kx;i++) {
      n[i]=parseInt(localStorage[str[i+2][nmax+4]]);
      message+=n[i];
      if(i<(kx-1)) message+=",";
  }
     message+="];<br>";
     message+="var variant=[";
       for (i=0;i<kx;i++) {
            for (j=0;j<n[i];j++) {
               z=parseFloat(localStorage[str[i+2][j+2]]);
               message+=z;
               if(i<(kx-1) || j<(n[i]-1)) message+=",";
            }
        if(i<(kx-1)) message+="<br>";
      }
      message+="];";message+="<br>";message+="break;";
   setupMessage("Save",message,"left","");
  }
//**********************************************************************************************************
function TableAction1() {
  let n=[],r=[],dmls=[],bmls=[],zpp=[],xp=[],xpl=[],xpu=[],xpborder=[],dv=[],dvsko=[],slgn=[],cx=[],rx=[],xycum=[];
  let xparam=[],xreg=[],yreg=[],wreg=[],dv_mp=[],wregp=[],femp=[],yr=[],yrp=[],emls=[],emlsp=[],Q=[],dmlsf=[];
  let dv_t=[],dvsko_t=[],slgn_t=[],bmls_t=[],t11reg=[],t22reg=[],t12reg=[],dv_tp=[],wcrit=[];
  let w=[]; 
  let xw=[0.01,0.05,0.1,0.3,0.5,0.7,0.9,0.95,0.99];
  ch_check=1;
   breg=0.95;
   nmax=parseInt(localStorage["AB"]);
   kx=parseInt(localStorage["AA"]);
    ts="normal";kmls=2;
    kmlsf=parseInt(localStorage["H"+(nmax+7)]);
    if(kmlsf>kx) {kmlsf=kx;localStorage["H"+(nmax+7)]=kmlsf;}
    if(kmlsf<2) {kmlsf=2;localStorage["H"+(nmax+7)]=kmlsf;}
   for (i=0;i<kmlsf;i++) {
      dmlsf[i]=[];  
    for (j=0;j<kmlsf;j++) dmlsf[i][j]=0;
    }
    for (ix=0;ix<kmls;ix++) {
    dmls[ix]=[];bmls[ix]=0;  
    for (jx=0;jx<kmls;jx++) dmls[ix][jx]=0;
    }
//**********************Data******************************************
    m2=0;
   preg=parseFloat(localStorage[str[2][nmax+19]]);

   for(ix=0;ix<kx;ix++) {
     nt=0;
     var xcum=[],fcum=[],ycum=[];
    for(jx=0;jx<nmax;jx++) {
     z=localStorage[str[ix+2][jx+2]];
      if(z==undefined || z=="" || z==0) break;
      nt++;
     }
     n[ix]=parseInt(nt);
     localStorage[str[ix+2][nmax+4]]= n[ix];
     cx[ix]=parseFloat(localStorage[str[ix+2][1]]);
     w[ix]=n[ix]; 
      
     for(jx=0;jx<n[ix];jx++) {
      xcum[jx]=parseFloat(localStorage[str[ix+2][jx+2]]); 
      }

// ############################Цензурирование##########################
        s1=0;s2=0;
        let xsimpl=[],rsimpl=[];
        kr=0;
        for (jx=0;jx<n[ix];jx++) {
           rsimpl[jx]=0;
           if(xcum[jx]<0) {rsimpl[jx]=1;kr=kr+1;xcum[jx]=-xcum[jx];}
           xdata[jx]=xcum[jx];
           s1=s1+(1-rsimpl[jx])*xcum[jx]; 
           s2=s2+(1-rsimpl[jx])*xcum[jx]*xcum[jx];
        }
        krr=n[ix]-kr;
        localStorage[str[ix+2][nmax+3]]=krr;
        dv[ix]=s1/krr;
        dvsko[ix]=Math.sqrt((s2-dv[ix]*dv[ix]*krr)/(krr-1)); 
        
         if(kr!=0) {
             xsimpl[0]=dv[ix];
             xsimpl[1]=dvsko[ix];
             stepx=0.01;eps=0.000001;lim=1000;
             limx=simpl(xsimpl,stepx,eps,lim,2,normal_minimize);
             dv[ix]=xsimpl[0];dvsko[ix]=xsimpl[1];
        }
           r[ix]=krr;
           rx[ix]=r[ix];
//######################################################
       t11reg[ix]=dvsko[ix]**2/n[ix];
       t22reg[ix]=0.5*t11reg[ix];
       t12reg[ix]=0;
       localStorage[str[ix+2][nmax+5]]= dv[ix].toFixed(4);
       localStorage[str[ix+2][nmax+6]]= dvsko[ix].toFixed(4);
       localStorage[str[ix+2][nmax+7]]= (dvsko[ix]/dv[ix]).toFixed(4);
       localStorage[str[ix+2][nmax+8]]=t11reg[ix].toFixed(4);
       localStorage[str[ix+2][nmax+9]]=t22reg[ix].toFixed(4);
       localStorage[str[ix+2][nmax+10]]=t12reg[ix].toFixed(4);

     for (jx=0;jx<xw.length;jx++) {
        z=invnormaldistribution(xw[jx]);
        xp[jx]= dv[ix]+dvsko[ix]*z;
        lmtaprn(breg,n[ix],t11reg[ix],t22reg[ix],t12reg[ix],z,xpborder);
        xpl[jx]=(dv[ix]+dvsko[ix]*xpborder[1]/Math.sqrt(n[ix]));
        xpu[jx]=(dv[ix]+dvsko[ix]*xpborder[0]/Math.sqrt(n[ix]));
        localStorage[str[6][jx+2]]=xw[jx].toFixed(4);
        localStorage[str[7+ix*3][jx+2]]=xpl[jx].toFixed(4);
        localStorage[str[8+ix*3][jx+2]]=xp[jx].toFixed(4);
        localStorage[str[9+ix*3][jx+2]]=xpu[jx].toFixed(4);
     }
//****************Критерий нормальности Шапиро-Уилка***********************
       wshapir=shapir(r[ix],xcum);
       wc=parseFloat(wcr005[n[ix]-3]);
       localStorage[str[ix+2][nmax+13]]=wshapir.toFixed(znaki);
       if(wshapir>=wc) localStorage[str[ix+2][nmax+14]]="H"+lt.sbd[0]+"+ "+wc;
       if(wshapir<wc) localStorage[str[ix+2][nmax+14]]="H"+lt.sbd[0]+"- "+wc;
       m2+=r[ix];
  }

//*************************************************************************
 if(kx>1) {
localStorage["D"+(nmax+16)]="";localStorage["D"+(nmax+16)]="";
//****************************Bartlet****************************************
     s1=0;s2=0;s3=0;nx=m2;xcp=0;
  for (i=0;i<kx;i++) {
    s1+=1/(r[i]-1);
    s2+=(r[i]-1)*Math.log(dvsko[i]*dvsko[i]);
    s3+=dvsko[i]*dvsko[i]*(r[i]-1);
    xcp+=r[i]*dv[i];
}
    df2=nx-kx;df1=kx-1;s3=s3/df2;cko=Math.sqrt(s3);c=1+(s1-1/df2)/(3*df1);
    xistat=(df2*Math.log(s3)-s2)/c;
    xicrit=invchisquaredistribution(df1,0.05);
    localStorage["B"+(nmax+15)]=lt.greek.chi+"="+xistat.toFixed(znaki);
    localStorage["C"+(nmax+15)]=lt.greek.chi+lt.sbd[0]+","+lt.sbd[0]+lt.sbd[5]+"="+xicrit.toFixed(znaki);
    if(xistat<=xicrit)   {
        localStorage["D"+(nmax+15)]="s="+cko.toFixed(znaki);
        localStorage["D"+(nmax+16)]="H"+lt.sbd[0]+"+";
    }
    if(xistat>xicrit)   localStorage["D"+(nmax+15)]="H"+lt.sbd[0]+"-";
//********************Fisher***************************************************
   cko1=0;xcp=xcp/nx;
   for (i=0;i<kx;i++) cko1+=r[i]*(dv[i]-xcp)*(dv[i]-xcp);
   cko1=cko1/df1;fstat=cko1/s3;
   fcrit=invfdistribution(df1,df2,0.05);
   localStorage["B"+(nmax+16)]="F="+fstat.toFixed(znaki);
   localStorage["C"+(nmax+16)]="F"+lt.sbd[0]+","+lt.sbd[0]+lt.sbd[5]+"="+fcrit.toFixed(znaki);
   if(fstat<=fcrit)   {
       localStorage["D"+(nmax+16)]="a="+xcp.toFixed(znaki);
       localStorage["D"+(nmax+17)]="H"+lt.sbd[0]+"+";
  }
   if(fstat>fcrit)     localStorage["D"+(nmax+16)]="H"+lt.sbd[0]+"-";
//********************Kruskal*************************************************************
 var outh=[];
 kruscount(0.05,r,xycum,outh);
 hstat=outh[1];hcrit=outh[4];
 localStorage["B"+(nmax+17)]="H="+hstat.toFixed(4);
 localStorage["C"+(nmax+17)]="H"+lt.sbd[0]+","+lt.sbd[0]+lt.sbd[5]+"="+hcrit.toFixed(4);
 if(hstat<=hcrit)   localStorage["D"+(nmax+17)]="H"+lt.sbd[0]+"+";
 if(hstat>hcrit)     localStorage["D"+(nmax+17)]="H"+lt.sbd[0]+"-";
}

  //************************Sy=A*(y)^alfa*******************************************************************
 if(kx>2) {
        for(ix=0;ix<kx;ix++) {
          yreg[ix]=Math.log(dvsko[ix]); 
          xreg[ix]=Math.log(dv[ix]);
          wreg[ix]=1;
       }

      Mnk_regress(kx,wreg,xreg,yreg,xparam);
 
      alpha_slgn=xparam[1];
      a_slgn=Math.exp(xparam[0]-alpha_slgn*xparam[2]);
      localStorage[str[2][nmax+11]]=a_slgn.toFixed(4);
      localStorage[str[3][nmax+11]]=alpha_slgn.toFixed(4);
      for(i=0;i<kx;i++) {
        slgn[i]=a_slgn*Math.pow(dv[i],alpha_slgn);
      }
}
//**********************************Regress*****************************************************************
  if(kx<=2) return;

  var letr=[];

 a_slgn=parseFloat(localStorage[str[2][nmax+11]]);
 alpha_slgn=parseFloat(localStorage[str[3][nmax+11]]);

for (i=0;i<kx;i++) {
   wreg[i]=slgn[i]/dv[i];
   wreg[i]=n[i]/wreg[i]**2;
   xreg[i]=cx[i];
   localStorage[str[i+2][nmax+12]]=slgn[i]; 
 }
  
  emls[1]=0;emls[2]=0;
  fatiq(kx,wreg,cx,dv,emls);
  for (i=0;i<kx;i++) {
   yr[i]=emls[2]+emls[3]*dv[i]**(-emls[1]);
   localStorage[str[i+2][nmax+18]]=yr[i]; 
  }

  for(i=0;i<11;i++)  {
       localStorage[str[i+8][nmax+3]]=""; 
       localStorage[str[i+8][nmax+4]]=""; 
       localStorage[str[i+8][nmax+5]]=""; 
       localStorage[str[i+8][nmax+6]]="";
       localStorage[str[i+9][nmax+3]]=""; 
       localStorage[str[i+9][nmax+4]]=""; 
       localStorage[str[i+9][nmax+5]]=""; 
       localStorage[str[i+9][nmax+6]]="";
 }

    kp=xw.length;
    for(ip=0;ip<kp;ip++) {
       preg=xw[ip];
       localStorage[str[ip+8][nmax+2]]=preg;
       zpreg=invnormaldistribution(preg);
       emls[1]=0;emls[2]=0;
 
       for (jp=0;jp<kx;jp++) {
         dv_mp[jp]=dv[jp]+zpreg*slgn[jp];
         wregp[jp]=wreg[jp]/(1.+zpreg**2/2.);
      }

      fatiq(kx,wregp,cx,dv_mp,emls);

      localStorage[str[ip+8][nmax+3]]=emls[0]; //qp
      localStorage[str[ip+8][nmax+4]]=emls[1];  //alphap
      localStorage[str[ip+8][nmax+5]]=emls[2];  //pvp
      localStorage[str[ip+8][nmax+6]]=emls[3];  //anp
      localStorage[str[ip+8][nmax+7]]=emls[4];  //dap
      localStorage[str[ip+8][nmax+8]]=emls[5];  //dbp   
  }

     preg=parseFloat(localStorage[str[2][nmax+19]]);
     zpreg=invnormaldistribution(preg);
     for (jp=0;jp<kx;jp++) {
         dv_mp[jp]=dv[jp]+zpreg*slgn[jp];
         wregp[jp]=wreg[jp]/(1.+zpreg**2/2.);
      }
      emls[1]=0;emls[2]=0;
      fatiq(kx,wregp,cx,dv_mp,emls);
      alphap=emls[1];pvp=emls[2];anp=emls[3];
      localStorage[str[2][nmax+20]]=alphap;
      localStorage[str[2][nmax+21]]=pvp;
      localStorage[str[2][nmax+22]]=anp;
return;
  }
//**********************************************************************************************************
function Distribution(gchart) {

  var xw=[0.01,0.05,0.1,0.3,0.5,0.7,0.9,0.95,0.99];
  var dv=[],ds=[];
  var xycum=[],femp=[];
 
  kx=parseInt(localStorage["AA"]);nmax=parseInt(localStorage["AB"]);
  sm=4*kx;
   xtext="Y";ytext="Zp";
  ts="normal";StatTitle=ts+" Distribution";graphtext=StatTitle;
 var pzx=[];
 pzx[0]=xw[0];pzx[1]=xw[xw.length-1];
 xmax=-100000;xmin=100000;

 for (ix=0;ix<kx;ix++) {
   s_x[ix]=[];s_y[ix]=[];
   dv[ix]=parseFloat(localStorage[str[ix+2][nmax+5]]);
   ds[ix]=parseFloat(localStorage[str[ix+2][nmax+6]]);
 
    //xp (по двум точкам)
   for (j=0;j<2;j++) {
    if(ts=="normal") zxx=invnormaldistribution(pzx[j]);
    if(ts=="weibull") zxx= Math.log(Math.log(1/(1-pzx[j])));
    s_x[ix][j]=dv[ix]+zxx*ds[ix];
    s_y[ix][j]=zxx;
    s_type[ix]="line";s_color[ix]="black";s_marker[ix]="none";s_markersize[ix]="9";s_thick[ix]="2";s_dash[ix]="solid";
 }

 //low
  s_x[ix+kx]=[];s_y[ix+kx]=[];
 for (j=0;j<xw.length;j++) {
      if(ts=="normal") zxx=invnormaldistribution(xw[j]);
      if(ts=="weibull") zxx= Math.log(Math.log(1/(1-xw[j])));
      s_x[ix+kx][j]=parseFloat(localStorage[str[ix*3+7][j+2]]);
      if(s_x[ix+kx][j]<xmin) xmin=s_x[ix+kx][j];
      s_y[ix+kx][j]=zxx;
      s_type[ix+kx]="spline";s_color[ix+kx]="green";s_marker[ix+kx]="none";s_markersize[ix+kx]="9";s_thick[ix+kx]="1";s_dash[ix+kx]="solid";
 }

   //upper
     s_x[ix+2*kx]=[];s_y[ix+2*kx]=[];
   for (j=0;j<xw.length;j++) {
     if(ts=="normal") zxx=invnormaldistribution(xw[j]);
     if(ts=="weibull") zxx= Math.log(Math.log(1/(1-xw[j])));
     s_x[ix+2*kx][j]=parseFloat(localStorage[str[ix*3+9][j+2]]);
     if(s_x[ix+2*kx][j]>xmax) xmax=s_x[ix+2*kx][j];
     s_y[ix+2*kx][j]=zxx;
     s_type[ix+2*kx]="spline";s_color[ix+2*kx]="green";s_marker[ix+2*kx]="none";s_markersize[ix+2*kx]="9";s_thick[ix+2*kx]="1";s_dash[ix+2*kx]="solid";
  }

     //data point

     var xcum=[],fcum=[],ycum=[];
      var rx=[],n=[];
      n[ix]=parseInt(localStorage[str[ix+2][nmax+4]]);

     for(jx=0;jx<n[ix];jx++)  xcum[jx]=parseFloat(localStorage[str[ix+2][jx+2]]); 
     cum(xcum,fcum,ycum);
     rx[ix]=fcum.length;
       
      s_x[ix+3*kx]=[];s_y[ix+3*kx]=[];

    for (j=0;j<rx[ix];j++) {
      pz=fcum[j];
      if(ts=="normal") zx=invnormaldistribution(pz);
      if(ts=="weibull") zx= Math.log(Math.log(1/(1-pz)));
      s_x[ix+3*kx][j]=ycum[j];
      s_y[ix+3*kx][j]=zx;
   }
    sm_points[ix]=2;
    sm_points[ix+kx]=xw.length;
    sm_points[ix+2*kx]=xw.length;
    sm_points[ix+3*kx]=rx[ix];
    s_type[ix+3*kx]="scatter";s_color[ix+3*kx]="red";s_marker[ix+3*kx]="cross";s_markersize[ix+3*kx]="9";s_thick[ix+3*kx]="2";s_dash[ix+3*kx]="solid";

}

 xmax=Math.round(xmax)+0.5;xmin=Math.round(xmin)-0.5;xinterval=0.5;ymin=0;ymax=0;yinterval=0.5;
 draw_graph(gchart,sm,sm_points,s_x,s_y,s_type,s_color,s_marker,s_markersize,s_thick,s_dash,graphtext,xtext,ytext,xmin,xmax,xinterval,ymin,ymax,yinterval,ts);
}
//***************************************************************************************************************
function Sy(gchart) {
  sm=2;xtext="Y";ytext="Sy";graphtext="Sy=f(y)";
  sm_points[0]=kx;s_x[0]=[];s_y[0]=[]; sm_points[1]=kx;s_x[1]=[];s_y[1]=[];
  var dv=[],ds=[];
  kx=parseInt(localStorage["AA"]);nmax=parseInt(localStorage["AB"]);
    for(i=0;i<kx;i++) {
      dv[i]=parseFloat(localStorage[str[i+2][nmax+5]]);
      ds[i]=parseFloat(localStorage[str[i+2][nmax+6]]);
      slgn=parseFloat(localStorage[str[i+2][nmax+12]]);
      s_x[0][i]=dv[i]; s_y[0][i]=slgn;
      s_x[1][i]=dv[i]; s_y[1][i]=ds[i];
  }
   s_type[0]="spline";s_color[0]="black";s_marker[0]="none";s_markersize[0]="9";s_thick[0]="2";s_dash[0]="solid";
   s_type[1]="scatter";s_color[1]="red";s_marker[1]="circle";s_markersize[1]="9";s_thick[1]="2";s_dash[1]="solid";
    ymax=Math.max.apply(Math,ds)+0.1;ymin=Math.min.apply(Math,ds)-0.1;yinterval=0.05;
    xmax=Math.max.apply(Math,dv);xmin=Math.min.apply(Math,dv);
    xmax=Math.round(xmax)+0.5;xmin=Math.round(xmin)-0.5;xinterval=1;
draw_graph(gchart,sm,sm_points,s_x,s_y,s_type,s_color,s_marker,s_markersize,s_thick,s_dash,graphtext,xtext,ytext,xmin,xmax,xinterval,ymin,ymax,yinterval,"");
}
//**********************************Regress***********************************************************************************
function Regress(gchart) {
  kx=parseInt(localStorage["AA"]);nmax=parseInt(localStorage["AB"]);
  let yreg=[];
  if(kx<3) return;
  
   //0.5
   alpha=parseFloat(localStorage[str[12][nmax+4]]);
   pv=parseFloat(localStorage[str[12][nmax+5]]);
   an=parseFloat(localStorage[str[12][nmax+6]]);
   //preg
   alphap=parseFloat(localStorage[str[2][nmax+20]]);
   pvp=parseFloat(localStorage[str[2][nmax+21]]);
   anp=parseFloat(localStorage[str[2][nmax+22]]);

   s_x[0]=[];s_y[0]=[],s_x[1]=[];s_y[1]=[],s_x[2]=[];s_y[2]=[];

   for(i=0;i<kx;i++) {
      yreg[i]=parseFloat(localStorage[str[i+2][1]]);
      s_y[0][i]=yreg[i];
      s_x[0][i]=parseFloat(localStorage[str[i+2][nmax+5]]);
   }

  sm=3;xtext="lgN";ytext="Sa";graphtext="Fatique Curve"; kint=20;
  sm_points[0]=kx;sm_points[1]=kint;sm_points[2]=kint;
    
    yinterval=10;ymax=yreg[0]+50;ymin=yreg[kx-1]-20;
    xmin=4;xmax=8;xinterval=0.5;

  for (j=0;j<kint;j++) {
    xcur=xmin+(xmax-xmin)*j/(kint-1);
    s_x[1][j]=xcur;
    s_y[1][j]=pv+an*Math.pow(xcur,-alpha);
    s_x[2][j]=xcur;
    s_y[2][j]=pvp+anp*Math.pow(xcur,-alphap);
  }
  
   s_type[0]="scatter";s_color[0]="red";s_marker[0]="circle";s_markersize[0]="9";s_thick[0]="2";s_dash[0]="solid";
   s_type[1]="spline";s_color[1]="black";s_marker[1]="none";s_markersize[1]="9";s_thick[1]="2";s_dash[1]="solid";
   s_type[2]="spline";s_color[2]="green";s_marker[2]="none";s_markersize[2]="9";s_thick[2]="2";s_dash[2]="solid";
draw_graph(gchart,sm,sm_points,s_x,s_y,s_type,s_color,s_marker,s_markersize,s_thick,s_dash,graphtext,xtext,ytext,xmin,xmax,xinterval,ymin,ymax,yinterval,"");
}

//**********************************Resurs***********************************************************************************
function Resurs(gchart) {
    kp=9;
    sm=3;
   s_x[0]=[];s_y[0]=[]; s_x[1]=[];s_y[1]=[]; s_x[2]=[];s_y[2]=[];
   xtext="lgN";ytext="P";graphtext="Resurs Distribution";
   sm_points[0]=kp;sm_points[1]=kp;sm_points[2]=kp;

    for(i=0;i<kp;i++) {
         s_y[0][i] =parseFloat(localStorage[str[3][i+2]]);
         s_x[0][i]=parseFloat(localStorage[str[2][i+15]]);

         s_y[1][i] =parseFloat(localStorage[str[3][i+2]]);
         s_x[1][i]=parseFloat(localStorage[str[3][i+15]]);

         s_y[2][i] =parseFloat(localStorage[str[3][i+2]]);
         s_x[2][i]=parseFloat(localStorage[str[4][i+15]]);
    }
    xmin=5;xmax=8;xinterval=0.2;
    ymin=0;ymax=1;yinterval=0.1;
   s_type[0]="scatter";s_color[0]="red";s_marker[0]="circle";s_markersize[0]="9";s_thick[0]="2";s_dash[0]="solid";
   s_type[1]="scatter";s_color[1]="black";s_marker[1]="circle";s_markersize[1]="9";s_thick[1]="2";s_dash[1]="solid";
   s_type[2]="scatter";s_color[2]="green";s_marker[2]="circle";s_markersize[2]="9";s_thick[2]="2";s_dash[2]="solid";
   draw_graph(gchart,sm,sm_points,s_x,s_y,s_type,s_color,s_marker,s_markersize,s_thick,s_dash,graphtext,xtext,ytext,xmin,xmax,xinterval,ymin,ymax,yinterval,"");
}

//***************************************************************************
function kruscount(alfa,m2x,x,outh) {
  var metka=[];
   kx=m2x.length;
   xmetka(x,m2x,metka);
       n=0;s=0;
   for (i=0;i<kx;i++) n+=m2x[i];
     for (i=0;i<kx;i++){
        r=0;
        for (j=0;j<n;j++)  if (metka[j]==i) r+=(j+1);
         s+=r*r/m2x[i];
    }
 
     H=12*s/(n*(n+1))-3*(n+1); H1=0.5*H*(1+(n-kx)/(n-1-H));
     nx=n;
     falfa=invfdistribution(kx-1,nx-kx,alfa);chialfa=invchisquaredistribution(kx-1,alfa);Halfa=0.5*(falfa*(kx-1)+chialfa);
     n=nx;
     outh[0]=H; //H-statistic
     outh[1]=H1; //Hcor-statistic
     outh[2]=falfa; //F(kx-1,n-kx,alfa)
     outh[3]=chialfa; //Chi(kx-1,alfa)
     outh[4]=Halfa; //Hcr 
     outh[5]=H1+"<="+Halfa; 
}
//***********************************************************************************
function xmetka(x,m2x,metka) {
    n=0;
   kx=m2x.length;
   for (i=0;i<kx;i++) {
    for (j=0;j<m2x[i];j++) metka[j+n]=i;
     n+=m2x[i];
  }  
   for(i=0;i<n;i++) x[i]=parseFloat(x[i]);
     for (i=0;i<n-1;i++) {
       for (j=i+1;j<n;j++) {
        if (x[i]>x[j]) {
            xx=x[i];x[i]=x[j];x[j]=xx;
            xt=metka[i];metka[i]=metka[j];metka[j]=xt;
       }
     }
   }
}

//###########################################################
function TableAction2() {
     ch_check=4;
     sigmax=parseFloat(localStorage[str[2][2]]);
     cp=parseFloat(localStorage[str[2][3]]);
     ep=parseFloat(localStorage[str[2][4]]);
     m0=parseFloat(localStorage[str[2][5]]);
  
      let emls=[];

      beta = 0.95;
      kpolinom = 8;
      kp = 9;
      kpt = 4;
      ry = 1;
      nk = 21;
       
      let px=[],pv=[],a=[],alpha=[],ccx=[],yd=[],syy=[],hx=[],cv=[],pvx=[];
      let zp=[],n=[],w=[];

      for(i=0;i<kp;i++) {
        px[i] =parseFloat(localStorage[str[3][i+2]]);
        pv[i] =parseFloat(localStorage[str[4][i+2]]);
        a[i] =parseFloat(localStorage[str[5][i+2]]);
        alpha[i] =parseFloat(localStorage[str[6][i+2]]);
      }
      

      for(i=0;i<kpt;i++) {
        ccx[i] =parseFloat(localStorage[str[7][i+2]]);
        yd[i] =parseFloat(localStorage[str[8][i+2]]);
        syy[i] =parseFloat(localStorage[str[9][i+2]]);
        n[i] =parseInt(localStorage[str[10][i+2]]);
        syy[i]=syy[i]/yd[i];
     }
     for(i=0;i<nk;i++) {
       cv[i] =parseFloat(localStorage[str[11][i+2]]);
       pvx[i] =parseFloat(localStorage[str[12][i+2]]);
       localStorage[str[14][i+2]]=sigmax*cv[i];
      }

      //Расчет плотности распределения эксплуатационной нагруженности h[i]

      //1 Вариант. Аппроксимация на основе плотности нормального закона
      h=[];
      s1=0.;s2=0;
      for(i=0;i<nk;i++) {
         s1=s1+cv[i];s2=s2+(cv[i])**2;
      }
      cvm=s1/nk;
      cvstd=Math.sqrt((s2-cvm**2*nk)/(nk-1.)); 
      for(i=0;i<nk;i++) {
         z=(cv[i]-cvm)/cvstd;
         h[i]=(Math.exp(-0.5*z*z))/(s2pi*cvstd);
         localStorage[str[13][i+2]]=h[i];
      }
      k1=nk;  
       
       //2 Вариант. Полиномиальная аппроксимация спектра нагруженности
       //let h=[0.016742374,0.043890843,0.170680793,0.303101719,0.18763615,0.071089133,
              //0.061550515,0.074894518,0.039857085,0.004677445,0.01584167];
       //k1=h.length;
 
      //3 Вариант. Аппроксимация по эмпирической кривой распределения амплитуд cv=f(pvx)
      //
      //h=[];k1=cv.length;
      //for(j=0;j<nk;j++) h[j]=(pvx[j+1]-pvx[j])/(cv[j+1]-cv[j]);

      let c=[],yb=[],alphab=[],alphah=[],pvb=[],pvh=[],ab=[],ah=[],ypb=[],yph=[],emlsx=[]; 
       //Quantile curves
    for(jj=0;jj<kp;jj++) {
      for(j=0;j<7;j++) {
         yb[j]=5.+j*0.5;
         c[j]=pv[jj]+a[jj]*yb[j]**(-alpha[jj]);
         zpx=invnormaldistribution(px[jj]);
         yb[j]=Math.log(yb[j]);
   
        Fatiqp(zpx,kpt,ccx,syy,n,beta,a[jj],pv[jj],alpha[jj],yb[j],emlsx);
   
        ypb[j]=Math.exp(emlsx[0]);
        yph[j]=Math.exp(emlsx[1]);
        w[j]=1.;
      }

        emls[0]=0;emls[1]=0;emls[2]=0;emls[3]=0;emls[4]=0;emls[5]=0;

        fatiq(7,w,c,ypb,emls);

        alphab[jj]=emls[1];
        pvb[jj]=emls[2];
        ab[jj]=emls[3];

        localStorage[str[15][jj+2]]=alphab[jj];
        localStorage[str[16][jj+2]]=pvb[jj];
        localStorage[str[17][jj+2]]=ab[jj];

        emls[0]=0;emls[1]=0;emls[2]=0;emls[3]=0;emls[4]=0;emls[5]=0;

        fatiq(7,w,c,yph,emls);

        alphah[jj]=emls[1];
        pvh[jj]=emls[2];
        ah[jj]=emls[3];
        localStorage[str[18][jj+2]]=alphah[jj];
        localStorage[str[19][jj+2]]=pvh[jj];
        localStorage[str[20][jj+2]]=ah[jj];
   }
 
 
//##################Resurs###########################################
  gammax = (m0 + 1.) / m0;
  delta = 0.1;
  rx=0;
  for(jk = 0;jk<3;jk++) { //Цикл по доверительным границам
      for(jj=0;jj<kp;jj++) { //Цикл по вероятностям

//********Долговечности по исходной кривой усталости****************************************
         if(jk==0) rx = ((sigmax - pv[jj]) / a[jj])** (-1 / alpha[jj]); //0,5
         if(jk==1) rx = ((sigmax - pvb[jj]) / ab[jj])** (-1 / alphab[jj]); //up
         if(jk==2) rx = ((sigmax - pvh[jj]) / ah[jj]) ** (-1 / alphah[jj]); //low
         k = 0;
//***********Долговечности по вторичной кривой усталости*********************
kc=0;
while(kc==0) {
      ch = sigmax; cp1 = cp; ck = ch - ry; ck1 = ck - ry;
       if(jk==0) {
         rnh = ((ch - pv[jj]) / a[jj]) ** (-1 / alpha[jj]);
         rnk = ((ck - pv[jj]) / a[jj]) ** (-1 / alpha[jj]);
         rnk1 = ((ck1 - pv[jj]) / a[jj])** (-1 / alpha[jj]);
      }
      if (jk== 1) {
        rnh = ((ch - pvb[jj]) / ab[jj])** (-1 / alphab[jj]);
        rnk = ((ck - pvb[jj]) / ab[jj]) ** (-1 / alphab[jj]);
        rnk1 = ((ck1 - pvb[jj]) / ab[jj]) ** (-1 / alphab[jj]);
      }
       if(jk==2) {
          rnh = ((ch - pvh[jj]) / ah[jj])**(-1/alphah[jj]);
          rnk = ((ck - pvh[jj]) / ah[jj]) **(-1/alphah[jj]);
          rnk1 = ((ck1 - pvh[jj]) / ah[jj])**(-1/ alphah[jj]);
      }
     //***************Снижение статических свойств, цикл по спектру****************
  for(j=0;j<k1;j++) {
         x = h[j]*10**(rx-rnh);
         //***********Проверка условия разрушения по долговечности********************
         if(x >= 1) break;
         rkx = (ch / cp1) ** gammax;
         cpx = cp1 * (1 - (1 - rkx) * x) ** (1. / gammax);
         
         //***********Проверка условия разрушения по пределу прочности****************
         if(cpx <= ch) break;
         epx = ep * (cpx / cp1) ** (1./ (m0+1.));
         rnx1= Second(x, rnh, rnk, ch, ck, ep, epx, cp1, cpx);
         rnx2= Second(x, rnh, rnk1, ch, ck1, ep, epx, cp1, cpx);
         cp1 = cpx; rnh = rnx1; rnk = rnx2; ch = ch - ry; ck = ck - ry; ck1 = ck - ry;
         b =Math.log(ck / ch) /Math.log(rnx1 / rnx2);
         a2 =Math.log(ch) + b *Math.log((rnx1));
         rnk1 = Math.exp((a2 - Math.log((ck1))) / b);
}
        kc=0;
  if(x >= 1 || cpx <= ch) {
        if(k==0) {
           rx = rx - delta;
        }
       else {
          kc=1;
        }
   }
   else {
        k = 1;
        rx = rx + delta;
  }
} //end while

     localStorage[str[2+jk][jj+15]]=rx-delta/2;
  }  // Next jj
  }  //Next jk

//###############################################################
}
//####################################################
function fatiq(k,w,cx,y,emls) {
	
                let bx,alpha,pv;
                 alpha=emls[1];pv=emls[2];bx=1;xx=0;s1x=1;
                
	q = 100000.;
	s2 = 0; s3 = 0;
	for (i = 0; i < k; i++) {
		s2 += w[i];
		s3 += w[i] * Math.log(y[i]);
	}
	a = s3 / s2;
                 step = 0.1;
                 xmin=0;
                 xmax=cx[0]-step;
                 kn=parseInt((xmax-xmin)/step);
	if (pv != 0) kn = 1;
                 if (pv != 0) pvx = pv;

	for (j = 0; j <=kn; j++) {
		pvx =xmin+ j*(xmax-xmin)/kn;
		xcp = 0.;
		for (i = 0; i < k; i++)  xcp += w[i] *Math.log(cx[i] - pvx);
		xcp = xcp / s2;
		s1 = 0.; s3 = 0.;
		for (i = 0; i < k; i++) {
			s1 += w[i] * Math.pow(Math.log(cx[i] - pvx) - xcp, 2);
			s3 += w[i] * (Math.log(cx[i] - pvx) - xcp) *Math.log(y[i]);
		}
		if (alpha != 0) b = -1. / alpha;
		if (alpha == 0) b = s3 / s1;
		s3 = 0;
		for (i = 0; i < k; i++) {
			s3 += w[i] * Math.pow((Math.log(y[i]) - a - b * (Math.log(cx[i] - pvx) - xcp)), 2);
		}
		if (s3 < q) {
			q = s3; pv = pvx; bx = b; xx = xcp; s1x = s1;
		}
	}
                 emls[0]=q;
	emls[1]= -1. / bx; //alpha
	emls[2]=pv; //pv
	emls[3]=Math.exp(xx - a / bx); //an
	emls[4] = q / s2; //da
	emls[5]= q / s1x; //db
    return;
}
//##############################################
function Fatiqp(zp,k,ccx,syy,n,beta,a,pv,alpha,yb,emlsx) {

 zb=invnormaldistribution(beta);
 s1 = 0; s2 = 0; s3 = 0;
 for(i=0;i<k;i++) {
    w=(syy[i])**2*(1+zp**2/2.)/n[i];
    w=1/w;
    xx=Math.log(ccx[i]-pv);
    s1=s1+w*xx;
    s2=s2+w;
 }
    xcp = s1 / s2;
    da=1./s2;
    for(i=0;i<k;i++) {
      xx=Math.log(ccx[i]-pv);
      w=(syy[i])**2*(1+zp**2/2.)/n[i];
      w=1/w;
      s3=s3+w*(xx-xcp)**2;
    }
     db=1./s3;
     xp=Math.log(a)-alpha*yb;
     sp=Math.sqrt(da+db*(xp-xcp)**2);
     ypb=yb+zb*sp;
     yph=yb-zb*sp;
     emlsx[0]=ypb;
     emlsx[1]=yph;
}
//################################################################

function Second(x, rnh, rnk, ch, ck, ep, epx,cp, cpx) {
       teta = 1. / (0.085 + 0.37 * ep);
       z = Math.log10(1.- x);
       z1 = 1. - 0.675 *Math.log10(ch / cpx);
       z2 = 0.675 * Math.log10(cpx / cp);
       tetax = 1./ (0.085 + 0.37 * epx);
       rnx =Math.log10(ck / ch) * (tetax - teta) + rnk * (1 - 0.675 * Math.log10(ck / cp)) + z1 * z + rnh * z2;
       rnx = Math.abs(rnx) / (1. - 0.675 * Math.log10(ck / cpx));
       return(rnx);
}

//#############################################################
  function TableAction3() {
     ch_check=5;
     let lgN_r=[];
 
     pvx=parseFloat(localStorage[str[2][2]]);
     AX=parseFloat(localStorage[str[2][3]]);
     ALPHAX=parseFloat(localStorage[str[2][4]]);
     xu=parseFloat(localStorage[str[2][5]]);
  
     FF=1;

  
      ncv=stat_cv.length;
      s1=0.;s2=0;
      for(i=0;i<ncv;i++) {
         s1=s1+stat_cv[i];s2=s2+(stat_cv[i])**2;
      }
      cvm=s1/ncv;
      cvstd=Math.sqrt((s2-cvm**2*ncv)/(ncv-1.)); 

     stat_pvx=pvx;
     stat_AX=AX;
     stat_ALPHAX=ALPHAX;
     stat_xu=xu;
     stat_cvm=cvm;
     stat_cvstd=cvstd;

     u=0.5*pvx;
     z1=(xu-cvm)/cvstd;
     z2=(u-cvm)/cvstd;

     ap2=normaldistribution(z1)-normaldistribution(z2);
     ap1=simpson(fresurs1,u,xu,1000);
     ap=(ap1/ap2-u)/(xu-u);

     for (ix=0;ix<xprob.length;ix++) {
        stat_zp=invnormaldistribution(xprob[ix]);
        s3=simpson(fresurs2,pvx,xu,1000);
        lgN=Math.log10(ap/s3);
        lgN_r[ix]=lgN;
        tsum=ap/(s3*FF*3600);
        localStorage[str[3][ix+2]]=xprob[ix];
        localStorage[str[4][ix+2]]=lgN;
      }
   
  
  } //end function

//*****************************************************************************
    simpson=function(fun,xl,xu,nk) {
      sm=[1,4,1];
      sum=0;dm=(xu-xl)/nk;
      for(i=0;i<nk;i++)   for(j=0;j<3;j++)  sum+=sm[j]*fun(xl+i*dm+j*dm/2);
      return  sum*dm/6;
}
//**********************Integrated function resurs***********************************************
     function fresurs2(x) {
       z=(x/stat_xu-stat_cvm)/Math.sqrt(s2pi*stat_cvstd);
       phi=Math.exp(-0.5*z*z);
       z1=(x-stat_pvx)/stat_AX;
       if(z1<=0) return 0; 
       lgN=(z1)**(-1./stat_ALPHAX); 
       slgN=0.0384*(lgN)**0.979;
       lgNp=lgN+stat_zp*slgN;
       return phi/Math.pow(10,lgNp);
     }
//**********************************************************************************
   function fresurs1(x) {
       z=(x/stat_xu-stat_cvm)/Math.sqrt(s2pi*stat_cvstd);
       phi=Math.exp(-0.5*z*z);
       return x*phi;
       }

//**********************************Resurs_ap***********************************************************************************
function Resurs_ap(gchart) {
    kp=xprob.length;
    sm=1;
    s_x[0]=[];s_y[0]=[];
    xtext="lgN";ytext="P";graphtext="Resurs Distribution";
    sm_points[0]=kp;

    for(i=0;i<kp;i++) {
         s_y[0][i] =parseFloat(localStorage[str[3][i+2]]);
         s_x[0][i]=parseFloat(localStorage[str[4][i+2]]);
    }
    xmin=5;xmax=8;xinterval=0.2;
    ymin=0;ymax=1;yinterval=0.05;
   s_type[0]="scatter";s_color[0]="red";s_marker[0]="circle";s_markersize[0]="9";s_thick[0]="2";s_dash[0]="solid";
   draw_graph(gchart,sm,sm_points,s_x,s_y,s_type,s_color,s_marker,s_markersize,s_thick,s_dash,graphtext,xtext,ytext,xmin,xmax,xinterval,ymin,ymax,yinterval,"");
 
}


//############MLE Normal Minimized Function############################

function normal_minimize(xsimpl) {
    let i,kx,s1,s2,s3,s4,cp,cko,nsim,z,d,p,psi,c1,c2;
    s1 = 0; s2 = 0; s3 = 0; s4 = 0; kx = 0;
    cp=xsimpl[0];cko=xsimpl[1];
    nsim=xdata.length;
    if(cko<=0) return 10000;
    for (i=0;i<nsim;i++) {
            z=(xdata[i]-cp)/cko;
            d=spi * Math.exp(-z*z/2.);
            p=normaldistribution(z);
            psi=d/(1.-p);
            s1 +=(1.-rsimpl[i])*(xdata[i]-cp);
            s2 += (1.-rsimpl[i])*(xdata[i]-cp)**2;
            s3 += rsimpl[i]*psi;
            s4 += rsimpl[i]*psi*z;
            kx+=1-rsimpl[i];
    }
    c1=s1+cko*s3;
    c2=s2+(s4-kx)*(cko)**2;
    zsim=c1**2+c2**2;
    return zsim;
}