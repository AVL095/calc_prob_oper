function setupMessage(title,body,style,graph) {
  container = document.createElement('div');
  tst='';
  if(style=='left') tst='style="text-align:left"';
  if(style=='center') tst='style="text-align:center"';
  if(style=='right') tst='style="text-align:right"';
  t1='<div id="msg" class="my-message"><div class="my-message-title">'+title+'</div>';
  t2='<div id="body" class="my-message-body"'+tst+'>'+body+'</div>';
  t8='<button class="btn btn-success" onClick=CloseOnClick(); title="Exit"><span class="glyphicon glyphicon-pencil"></span>Exit</button>';
  t9='<button class="btn btn-success" onClick=CopyOnClick("body"); title="Copy"><span class="glyphicon glyphicon-pencil"></span>Copy</button>';
   if(graph!="" ) {
      var outg=[];
      getGraphOptions(outg);
      t=t1+t8+outg[0]+t2+outg[1];
     }
     else {
     t=t1+t8+t9+t2;
    }
    container.innerHTML=t;
    messageElem=container.firstChild;
    positionMessage(messageElem);
    document.body.appendChild(messageElem);
 }
//***************************************************
function positionMessage(elem) {
  elem.style.position = 'absolute';
  var scroll = document.documentElement.scrollTop || document.body.scrollTop;
  elem.style.top = scroll;
  elem.style.left = Math.floor(document.body.clientWidth/10) - 1 + 'px';
}
//**************************************************
function CloseOnClick() {
    messageElem.parentNode.removeChild(messageElem);
 }
//****************************************************

 function CopyOnClick(t) {
   window.getSelection().removeAllRanges();
   var ta = document.getElementById(t);
   var range = document.createRange();
   range.selectNode(ta);
   window.getSelection().addRange(range);
   document.execCommand('copy');
   window.getSelection().removeAllRanges();
}

 