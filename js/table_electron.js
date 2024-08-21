   var elmObject=[];
//***********************************************************************
function TableSaveDefault() {
         m1="";
         for (i=0;i<elmObject.length;i++) {
            elmid=elmObject[i].id;
            if(localStorage[elmid] !=undefined && localStorage[elmid] !="") m1=m1+elmid+':'+'"'+localStorage[elmid]+'",';
            }
  setupMessage("Save",m1,"center","");
 }
//***********************************************************************
function TablePrint() {
       mg="";
       mg+="<table border='2' align='center' bordercolor='black' width='auto' bgcolor='black'>";
       mg+="<th width='10%' bgcolor='white'>"+"</th>";
    for (var i=0;i<colmax-1; i++)   mg+="<th width='auto' bgcolor='white'>"+String.fromCharCode("A".charCodeAt(0)+i)+"</th>";
       mg+="</tr>";
         for (i=0; i<rowmax-1; i++) {
            for (j=0;j<colmax; j++) {
              ltter = String.fromCharCode("A".charCodeAt(0)+(j-1));
               mg+="<td  align='center' bgcolor='white'>";
               z=localStorage[ltter+(i+1)];
                    if(z==undefined) {
                    if(j==0) mg+=i+1;
                    if(j!=0) mg+=String.fromCharCode(" ".charCodeAt(0));
                }
                if(z!=undefined)  mg+=z;
             }
              mg+="</tr>";
         }
  mg+="</tr>";
  mg+="<table>";
  setupMessage("Result",mg,"center","");
 }
//************************************************************************
function TableCreate(fun) {
  for (var i=0; i<rowmax; i++) {
    var row = document.querySelector("table").insertRow(-1);
     for (var j=0;j<colmax;j++) {

           k1=parseInt(j/26);
           k2=parseInt(j/27);
           k=parseInt(26*k1);
           if(k2==0)   ltter=String.fromCharCode("A".charCodeAt(0)+(j-1));
           if(k2==1)   ltter=String.fromCharCode("A".charCodeAt(0)+(k1-1))+String.fromCharCode("A".charCodeAt(0)+(j-k-1));
           if(k2>1)   ltter=String.fromCharCode("A".charCodeAt(0)+(k1-1))+String.fromCharCode("A".charCodeAt(0)+(j-k));
           if(ltter=="B@") ltter="AZ";
        
        if(i==0 && j==0) ltter=lt.symbol.numb;
          letinput="<input id='"+ ltter+i+"'/>";
          row.insertCell(-1).innerHTML = i&&j ? letinput:i||ltter;
 }
}
          elmObject=[].slice.call(document.querySelectorAll("input"));
           TableInit();
           TableAction(fun);
   return;        
}
//***********************************************************************
function TableClear(fun) {
      localStorage.clear();TableAction(fun);
}
//****************************************************************************
function TableAction(fun) {
var CurrentData={};
elmObject.forEach(function(elm) {
       elm.onfocus = function(e) {
        e.target.style.backgroundColor ="white";
        e.target.value=localStorage[e.target.id] || "";
     };

elm.onchange = function(e) {
  localStorage[e.target.id] =e.target.value;
         funobject(fun,e.target);
         computeAll();
}
     var getter = function() {
        var value =localStorage[elm.id] || "";
           if (value.charAt(0) == "=") {
           with (CurrentData) return eval(value.substring(1));
        } else { return isNaN(parseFloat(value)) ? value : parseFloat(value);}
    };
Object.defineProperty(CurrentData, elm.id, {get:getter});Object.defineProperty(CurrentData, elm.id.toLowerCase(), {get:getter});
});
(window.computeAll = function() {
elmObject.forEach(function(elm) {
try {
      elm.value=CurrentData[elm.id];
     //if(localStorage[elm.id] !=undefined) elm.title=localStorage[elm.id];
   } 
catch(e) {} });
})();
}

//***********************************
var funobject=function(fun,etarget) {
      return fun(etarget);
}
function NumObject(row,col) {
   return (row-1)*(colmax-1)+col-1;
}
