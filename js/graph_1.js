//************************************************************************************************************
 function draw_graph(chartCont,sm,smpoints,sx,sy,stype,scolor,smarker,smarkersize,sthick,sdash,graphtext,xtext,ytext,xmin,xmax,xinterval,ymin,ymax,yinterval,ts) {
  
  //входные параметры
   //sm-количество графиков
   //smpoints[sm]-вектор количества точек в каждом графике размерности sm
   //sx[sm][smpoints[sm]]-двумерный массив значений x размерности sm*smpoints[sm]
   //sy[sm][smpoints[sm]]-двумерный массив значений y размерности sm*smpoints[sm]
   //stype[sm]-вектор типов графиков размерности sm
   //scolor[sm]-вектор цветов графиков размерности sm
  //smarker[sm]-вектор маркеров графиков размерности sm
  //graphtext - имя графика
  //xtext - имя оси х
  //ytext - имя оси y
  //xmin,xmax,xinterval-минимальное x, максимальное x и интервал по х
  //ymin,ymax,yinterval-минимальное y, максимальное y и интервал по y

  // document.getElementById(chartCont).style.border="4px double blue";
   var datagraph=[];
   if(ts=="normal")   var p=[0.99,0.95,0.9,0.7,0.5,0.3,0.1,0.05,0.01];
   if(ts=="weibull")   var p=[0.99,0.9,0.7,0.5,0.3,0.2,0.1,0.05,0.025,0.01];
      k=0;sgrid=1;
      if(ts!="") k=p.length;
   for (ix=0;ix<sm+k;ix++) {
     if(ix>=sm) {
       sx[ix]=[];
       sy[ix]=[];
       stype[ix]="line";scolor[ix]="blue";smarker[ix]="none";smarkersize[ix]="0";sthick[ix]="1";sdash[ix]="solid"; smpoints[ix]=2;
       sx[ix][0]=xmin;sx[ix][1]=xmax;
       if(ts=="normal") {sy[ix][0]=invnormaldistribution(p[ix-sm]);  ymin=-2.33;ymax=2.33;yinterval=0.5;}
       if(ts=="weibull") {sy[ix][0]=invweibulldistribution(p[ix-sm]);  ymin=-5.0;ymax=2.0;yinterval=0.5;}
       sy[ix][1]=sy[ix][0];
       sgrid=0;
    }
     var dataSeries = {type:stype[ix],color:scolor[ix],markerType:smarker[ix],markerSize:smarkersize[ix],lineThickness:sthick[ix],lineDashType:sdash[ix],markerBorderThickness:"1",click:onClick};
     var dataPoints = [];
    for(jx=0;jx<smpoints[ix];jx++) {
        if(ix<sm) dataPoints.push({x:sx[ix][jx],y:sy[ix][jx]});
        if(ix>=sm && jx==0) dataPoints.push({x:sx[ix][jx],y:sy[ix][jx]});
        if(ix>=sm && jx==1) dataPoints.push({x:sx[ix][jx],y:sy[ix][jx],indexLabel: '"'+p[ix-sm]+'"',indexLabelFontSize:"8",indexLabelFontColor:"blue",indexLabelFontFamily:"Times New Roman"});
    }
     dataSeries.dataPoints = dataPoints;
     datagraph.push(dataSeries);   
  } 
     var chart = new CanvasJS.Chart(chartCont,   {
        exportFileName: "Chart",
        exportEnabled: true,
         title:{
        text: graphtext,
        fontColor: "red",
        fontSize:"14",
        fontFamily:"Times New Roman",
        fontWeight:"bold",
        fontStyle:"italic",
        verticalAlign: "top",
        horizontalAlign: "center",
        size: 4
         },
        //exportFileName: "Chart",
        //exportEnabled: true,
        zoomEnabled: true,
        animationEnabled: true,
      axisX:{
        title:xtext,
        titleFontColor:"red",
        titleFontSize:"14",
        titleFontFamily:"Times New Roman",
        titleFontWeight:"bold",
        titleFontStyle:"italic",
        labelFontColor:"blue",
        labelFontSize:"12",
        labelFontFamily:"Times New Roman",
        labelFontWeight:"bold",
        labelFontStyle:"normal",
        lineColor: "blue",
        tickColor: "red",
        interval: xinterval,
        valueFormatString: "#0.0",
        gridThickness: 1,
        gridColor: "blue",
        minimum: xmin,
        maximum: xmax      
      },
       axisY:{
        title:ytext,
        titleFontColor:"red",
        titleFontSize:"14",
        titleFontFamily:"Times New Roman",
        titleFontWeight:"bold",
        titleFontStyle:"italic",
        labelFontColor:"blue",
        labelFontSize:"12",
        labelFontFamily:"Times New Roman",
        labelFontWeight:"bold",
        labelFontStyle:"normal",
        lineColor: "blue",
        tickColor: "red",
        interval: yinterval,
        valueFormatString: "#0.00",
        gridThickness: sgrid,
        gridColor: "blue",
        minimum: ymin,
        maximum: ymax      
      },
      data: datagraph
      });

     function onClick(e) {
  //     MsgClick(e);
        }
     chart.render();
 }
//**********************************************************************************
function createElem(chartCont) {
    var parent = document.getElementsByTagName('BODY')[0];
    var div = document.createElement('div');
    //div.className = 'elemClass';
    div.id =chartCont;
    div.style="width: 50%; height: 300px;border:4px double blue;display: inline-block";
    //div.style.background = '#c5c5c5';
    //div.style.color = '#000';
    //var text = "текст для вставки";
    //var textNode = document.createTextNode(text);
   //div.appendChild(textNode);
    parent.appendChild(div);
}