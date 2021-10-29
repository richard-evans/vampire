/** X3DOM Runtime, http://www.x3dom.org/ 1.7.2 - 61a235203deb34329fe615cbbf21314db6ebf49f - Mon Dec 19 19:17:05 2016 +0100 */
if(!Array.forEach){Array.forEach=function(array,fun,thisp){var len=array.length;for(var i=0;i<len;i++){if(i in array){fun.call(thisp,array[i],i,array);}}};}
if(!Array.map){Array.map=function(array,fun,thisp){var len=array.length;var res=[];for(var i=0;i<len;i++){if(i in array){res[i]=fun.call(thisp,array[i],i,array);}}
return res;};}
if(!Array.filter){Array.filter=function(array,fun,thisp){var len=array.length;var res=[];for(var i=0;i<len;i++){if(i in array){var val=array[i];if(fun.call(thisp,val,i,array)){res.push(val);}}}
return res;};}
var x3dom={canvases:[],x3dNS:'http://www.web3d.org/specifications/x3d-namespace',x3dextNS:'http://philip.html5.org/x3d/ext',xsltNS:'http://www.w3.org/1999/XSL/x3dom.Transform',xhtmlNS:'http://www.w3.org/1999/xhtml'};x3dom.nodeTypes={};x3dom.nodeTypesLC={};x3dom.components={};x3dom.geoCache=[];x3dom.caps={PLATFORM:navigator.platform,AGENT:navigator.userAgent,RENDERMODE:"HARDWARE"};x3dom.registerNodeType=function(nodeTypeName,componentName,nodeDef){if(x3dom.components[componentName]===undefined){x3dom.components[componentName]={};}
nodeDef._typeName=nodeTypeName;nodeDef._compName=componentName;x3dom.components[componentName][nodeTypeName]=nodeDef;x3dom.nodeTypes[nodeTypeName]=nodeDef;x3dom.nodeTypesLC[nodeTypeName.toLowerCase()]=nodeDef;};x3dom.isX3DElement=function(node){var name=(node.nodeType===Node.ELEMENT_NODE&&node.localName)?node.localName.toLowerCase():null;return(name&&(x3dom.nodeTypes[node.localName]||x3dom.nodeTypesLC[name]||name=="x3d"||name=="websg"||name=="route"));};x3dom.extend=function(f){function G(){}
G.prototype=f.prototype||f;return new G();};x3dom.getStyle=function(oElm,strCssRule){var strValue="";var style=document.defaultView.getComputedStyle?document.defaultView.getComputedStyle(oElm,null):null;if(style){strValue=style.getPropertyValue(strCssRule);}
else if(oElm.currentStyle){strCssRule=strCssRule.replace(/\-(\w)/g,function(strMatch,p1){return p1.toUpperCase();});strValue=oElm.currentStyle[strCssRule];}
return strValue;};function defineClass(parent,ctor,methods){if(parent){function Inheritance(){}
Inheritance.prototype=parent.prototype;ctor.prototype=new Inheritance();ctor.prototype.constructor=ctor;ctor.superClass=parent;}
if(methods){for(var m in methods){ctor.prototype[m]=methods[m];}}
return ctor;}
x3dom.isa=function(object,clazz){return(object instanceof clazz);};x3dom.getGlobal=function(){return(function(){return this;}).call(null);};x3dom.loadJS=function(src,path_prefix,blocking){blocking=(blocking===false)?blocking:true;if(blocking){var url=(path_prefix)?path_prefix.trim()+src:src;var req=new XMLHttpRequest();if(req){req.open("GET",url,false);req.send(null);eval(req.responseText);}}else{var head=document.getElementsByTagName('HEAD').item(0);var script=document.createElement("script");var loadpath=(path_prefix)?path_prefix.trim()+src:src;if(head){x3dom.debug.logError("Trying to load external JS file: "+loadpath);script.type="text/javascript";script.src=loadpath;head.appendChild(script);}else{alert("No document object found. Can't load components!");}}};function array_to_object(a){var o={};for(var i=0;i<a.length;i++){o[a[i]]='';}
return o;}
window.requestAnimFrame=(function(){return window.requestAnimationFrame||window.webkitRequestAnimationFrame||window.mozRequestAnimationFrame||window.oRequestAnimationFrame||window.msRequestAnimationFrame||function(callback,element){window.setTimeout(callback,16);};})();x3dom.toggleFullScreen=function(){if(document.fullScreen||document.mozFullScreen||document.webkitIsFullScreen){if(document.cancelFullScreen){document.cancelFullScreen();}
else if(document.mozCancelFullScreen){document.mozCancelFullScreen();}
else if(document.webkitCancelFullScreen){document.webkitCancelFullScreen();}}
else{var docElem=document.documentElement;if(docElem.requestFullScreen){docElem.requestFullScreen();}
else if(docElem.mozRequestFullScreen){docElem.mozRequestFullScreen();}
else if(docElem.webkitRequestFullScreen){docElem.webkitRequestFullScreen();}}};x3dom.debug={INFO:"INFO",WARNING:"WARNING",ERROR:"ERROR",EXCEPTION:"EXCEPTION",isActive:false,isFirebugAvailable:false,isSetup:false,isAppend:false,numLinesLogged:0,maxLinesToLog:10000,logContainer:null,setup:function(){if(x3dom.debug.isSetup){return;}
try{if(window.console.firebug!==undefined){x3dom.debug.isFirebugAvailable=true;}}
catch(err){x3dom.debug.isFirebugAvailable=false;}
x3dom.debug.setupLogContainer();x3dom.debug.isSetup=true;},activate:function(visible){x3dom.debug.isActive=true;x3dom.debug.logContainer.style.display=(visible)?"block":"none";if(!x3dom.debug.isAppend){if(navigator.appName=="Microsoft Internet Explorer"){x3dom.debug.logContainer.style.marginLeft="8px";document.documentElement.appendChild(x3dom.debug.logContainer);}else{document.body.appendChild(x3dom.debug.logContainer);}
x3dom.debug.isAppend=true;}},setupLogContainer:function(){x3dom.debug.logContainer=document.createElement("div");x3dom.debug.logContainer.id="x3dom_logdiv";x3dom.debug.logContainer.setAttribute("class","x3dom-logContainer");x3dom.debug.logContainer.style.clear="both";},doLog:function(msg,logType){if(!x3dom.debug.isActive){return;}
if(x3dom.debug.numLinesLogged===x3dom.debug.maxLinesToLog){msg="Maximum number of log lines (="+x3dom.debug.maxLinesToLog+") reached. Deactivating logging...";}
if(x3dom.debug.numLinesLogged>x3dom.debug.maxLinesToLog){return;}
var node=document.createElement("p");node.style.margin=0;switch(logType){case x3dom.debug.INFO:node.style.color="#00ff00";break;case x3dom.debug.WARNING:node.style.color="#cd853f";break;case x3dom.debug.ERROR:node.style.color="#ff4500";break;case x3dom.debug.EXCEPTION:node.style.color="#ffff00";break;default:node.style.color="#00ff00";break;}
try{node.innerHTML=logType+": "+msg;x3dom.debug.logContainer.insertBefore(node,x3dom.debug.logContainer.firstChild);}catch(err){if(window.console.firebug!==undefined){window.console.warn(msg);}}
if(x3dom.debug.isFirebugAvailable){switch(logType){case x3dom.debug.INFO:window.console.info(msg);break;case x3dom.debug.WARNING:window.console.warn(msg);break;case x3dom.debug.ERROR:window.console.error(msg);break;case x3dom.debug.EXCEPTION:window.console.debug(msg);break;default:break;}}
x3dom.debug.numLinesLogged++;},logInfo:function(msg){x3dom.debug.doLog(msg,x3dom.debug.INFO);},logWarning:function(msg){x3dom.debug.doLog(msg,x3dom.debug.WARNING);},logError:function(msg){x3dom.debug.doLog(msg,x3dom.debug.ERROR);},logException:function(msg){x3dom.debug.doLog(msg,x3dom.debug.EXCEPTION);},assert:function(c,msg){if(!c){x3dom.debug.doLog("Assertion failed in "+
x3dom.debug.assert.caller.name+': '+
msg,x3dom.debug.ERROR);}},typeOf:function(obj){var type=typeof obj;return type==="object"&&!obj?"null":type;},exists:function(obj,name,type){type=type||"function";return(obj?this.typeOf(obj[name]):"null")===type;},dumpFields:function(node){var str="";for(var fName in node){str+=(fName+", ");}
str+='\n';x3dom.debug.logInfo(str);return str;}};x3dom.debug.setup();x3dom.arc={};x3dom.arc.instance=null;x3dom.arc.Limits=function(min,max,initial)
{this._min=min;this._max=max;this.getValue=function(value)
{value=this._min+(this._max-this._min)*value;return this._max>=value?(this._min<=value?value:this._min):this._max;};};x3dom.arc.ARF=function(name,min,max,dirFac,factorGetterFunc,factorSetterFunc,getterFunc,setterFunc)
{this._name=name;this._stateValue=[0.5,0.5];this._limits=new x3dom.arc.Limits(min,max);this._factorGetterFunc=factorGetterFunc;this._factorSetterFunc=factorSetterFunc;this._setterFunc=setterFunc;this._getterFunc=getterFunc;this._dirFac=dirFac;this.getFactor=function()
{return this._factorGetterFunc();};this.update=function(state,step)
{var stateVal=this._stateValue[state]+step*this._dirFac;this._stateValue[state]=0<=stateVal?(1>=stateVal?stateVal:1):0;this._setterFunc(this._limits.getValue(this._stateValue[state]));};this.reset=function()
{this._stateValue[0]=0.5;this._stateValue[1]=0.5;};};x3dom.arc.AdaptiveRenderControl=defineClass(null,function(scene)
{x3dom.arc.instance=this;this._scene=scene;this._targetFrameRate=[];this._targetFrameRate[0]=this._scene._vf.minFrameRate;this._targetFrameRate[1]=this._scene._vf.maxFrameRate;this._currentState=0;var that=this;var environment=that._scene.getEnvironment();this._arfs=[];this._arfs.push(new x3dom.arc.ARF("smallFeatureCulling",0,10,-1,function()
{return environment._vf.smallFeatureFactor;},function(value)
{environment._vf.smallFeatureFactor=value;},function()
{return environment._vf.smallFeatureThreshold;},function(value)
{environment._vf.smallFeatureThreshold=value;}));this._arfs.push(new x3dom.arc.ARF("lowPriorityCulling",0,100,1,function()
{return environment._vf.lowPriorityFactor;},function(value)
{environment._vf.lowPriorityFactor=value;},function()
{return environment._vf.lowPriorityThreshold*100;},function(value)
{environment._vf.lowPriorityThreshold=value/100;}));this._arfs.push(new x3dom.arc.ARF("tessellationDetailCulling",1,12,-1,function()
{return environment._vf.tessellationErrorFactor;},function(value)
{environment._vf.tessellationErrorFactor=value;},function()
{return environment.tessellationErrorThreshold;},function(value)
{environment.tessellationErrorThreshold=value;}));this._stepWidth=0.1;},{update:function(state,fps)
{this._currentState=state;var delta=fps-this._targetFrameRate[state];this._stepWidth=Math.abs(delta)>10?0.1:0.01;var factorSum=0;var normFactors=[];var i,n=this._arfs.length;for(i=0;i<n;++i)
{normFactors[i]=this._arfs[i].getFactor();if(normFactors[i]>0)
factorSum+=normFactors[i];}
var dirFac=delta<0?-1:1;for(i=0;i<n;++i)
{if(normFactors[i]>0)
{normFactors[i]/=factorSum;this._arfs[i].update(state,this._stepWidth*normFactors[i]*dirFac);}}},reset:function()
{for(var i=0,n=this._arfs.length;i<n;++i)
{this._arfs[i].reset();}}});x3dom.Request=function(url,onloadCallback,priority){this.url=url;this.priority=priority;this.xhr=new XMLHttpRequest();this.onloadCallbacks=[onloadCallback];var self=this;this.xhr.onload=function(){if(x3dom.DownloadManager.debugOutput){x3dom.debug.logInfo('Download manager received data for URL \''+self.url+'\'.');}
--x3dom.DownloadManager.activeDownloads;if((x3dom.DownloadManager.stallToKeepOrder===false)||(x3dom.DownloadManager.resultGetsStalled(self.priority)===false)){var i;for(i=0;i<self.onloadCallbacks.length;++i){self.onloadCallbacks[i](self.xhr.response);}
x3dom.DownloadManager.removeDownload(self);x3dom.DownloadManager.updateStalledResults();}
else if(x3dom.DownloadManager.debugOutput){x3dom.debug.logInfo('Download manager stalled downloaded result for URL \''+self.url+'\'.');}
x3dom.DownloadManager.tryNextDownload();};};x3dom.Request.prototype.send=function(){this.xhr.open('GET',encodeURI(this.url),true);this.xhr.responseType='arraybuffer';this.xhr.send(null);if(x3dom.DownloadManager.debugOutput){x3dom.debug.logInfo('Download manager posted XHR for URL \''+this.url+'\'.');}};x3dom.DownloadManager={requests:[],maxDownloads:6,activeDownloads:0,debugOutput:false,stallToKeepOrder:false,toggleDebugOutput:function(flag){this.debugOutput=flag;},toggleStrictReturnOrder:function(flag){this.stallToKeepOrder=false;},removeDownload:function(req){var i,j;var done=false;for(i=0;i<this.requests.length&&!done;++i){if(this.requests[i]){for(j=0;j<this.requests[i].length;++j){if(this.requests[i][j]===req){this.requests[i].splice(j,1);done=true;break;}}}}},tryNextDownload:function(){var firstRequest;var i,j;if(this.activeDownloads<this.maxDownloads){for(i=0;i<this.requests.length&&!firstRequest;++i){if(this.requests[i]){for(j=0;j<this.requests[i].length;++j){if(this.requests[i][j].xhr.readyState===XMLHttpRequest.UNSENT){firstRequest=this.requests[i][j];break;}}}}
if(firstRequest){firstRequest.send();++this.activeDownloads;}}},resultGetsStalled:function(priority){var i;for(i=0;i<priority;++i){if(this.requests[i]&&this.requests[i].length){return true;}}
return false;},updateStalledResults:function(){if(x3dom.DownloadManager.stallToKeepOrder){var i,j,k;var req,pendingRequestFound=false;for(i=0;i<this.requests.length&&!pendingRequestFound;++i){if(this.requests[i]){for(j=0;j<this.requests[i].length;++j){req=this.requests[i][j];if(req.xhr.readyState===XMLHttpRequest.DONE){if(x3dom.DownloadManager.debugOutput){x3dom.debug.logInfo('Download manager releases stalled result for URL \''+req.url+'\'.');}
for(k=0;k<req.onloadCallbacks.length;++k){req.onloadCallbacks[k](req.xhr.response);}
this.requests[i].splice(j,1);}
else{pendingRequestFound=true;}}}}}},get:function(urls,onloadCallbacks,priorities){var i,j,k,r;var found=false;var url,onloadCallback,priority;if(urls.length!==onloadCallbacks.length||urls.length!==priorities.length)
{x3dom.debug.logError('DownloadManager: The number of given urls, onload callbacks and priorities is not equal. Ignoring requests.');return;}
for(k=0;k<urls.length;++k){if(!onloadCallbacks[k]===undefined||!priorities[k]===undefined){x3dom.debug.logError('DownloadManager: No onload callback and / or priority specified. Ignoring request for \"'+url+'\"');continue;}
else{url=urls[k];onloadCallback=onloadCallbacks[k];priority=priorities[k];for(i=0;i<this.requests.length&&!found;++i){if(this.requests[i]){for(j=0;j<this.requests[i].length;++j){if(this.requests[i][j].url===url){this.requests[i][j].onloadCallbacks.push(onloadCallback);if(x3dom.DownloadManager.debugOutput){x3dom.debug.logInfo('Download manager appended onload callback for URL \''+url+'\' to a registered request using the same URL.');}
found=true;break;}}}}
if(!found){r=new x3dom.Request(url,onloadCallback,priority);if(this.requests[priority]!=undefined){this.requests[priority].push(r);}
else{this.requests[priority]=[r];}}}}
for(i=0;i<urls.length&&this.activeDownloads<this.maxDownloads;++i){this.tryNextDownload();}},abortAllDownloads:function()
{var request;for(var i=0;i<this.requests.length;i++)
{if(this.requests[i]!=undefined)
{for(var j=0;j<this.requests[i].length;j++)
{request=this.requests[i][j];request.xhr.abort();this.removeDownload(request);}}}}};x3dom.RequestManager={};x3dom.RequestManager.requests=[];x3dom.RequestManager.maxParallelRequests=40;x3dom.RequestManager.failedRequests=0;x3dom.RequestManager.loadedRequests=0;x3dom.RequestManager.totalRequests=0;x3dom.RequestManager.activeRequests=[];x3dom.RequestManager.requestHeaders=[];x3dom.RequestManager.withCredentials=false;x3dom.RequestManager.addRequestHeader=function(header,value)
{this.requestHeaders.push({header:header,value:value});};x3dom.RequestManager._sendRequest=function()
{if(this.activeRequests.length>this.maxParallelRequests)
{return;}
var request=this.requests.pop();if(request)
{this.activeRequests.push(request);request.send(null);this._sendRequest();}};x3dom.RequestManager.addRequest=function(request)
{if(!(request instanceof XMLHttpRequest))
{return;}
this.totalRequests++;request.withCredentials=this.withCredentials;for(var i=0;i<this.requestHeaders.length;i++)
{var header=this.requestHeaders[i].header;var value=this.requestHeaders[i].value;request.setRequestHeader(header,value);}
request.addEventListener("load",this._onLoadHandler.bind(this));request.addEventListener("error",this._onErrorHandler.bind(this));this.requests.push(request);this._sendRequest();};x3dom.RequestManager.abortAllRequests=function()
{for(var i=0;i<this.activeRequests.length;i++)
{this.activeRequests[i].abort();}
this.requests=this.activeRequests=[];};x3dom.RequestManager._removeActiveRequest=function(request)
{var idx=this.activeRequests.indexOf(request);return this.activeRequests.splice(idx,1);};x3dom.RequestManager._onLoadHandler=function(e)
{this._removeActiveRequest(e.target);this.loadedRequests++;this._sendRequest();};x3dom.RequestManager._onErrorHandler=function(e)
{this._removeActiveRequest(e.target);this.failedRequests++;this._sendRequest();};x3dom.MultiMaterial=function(params)
{this._origAmbientIntensity=params.ambientIntensity;this._origDiffuseColor=params.diffuseColor;this._origEmissiveColor=params.emissiveColor;this._origShininess=params.shininess;this._origSpeclarColor=params.specularColor;this._origTransparency=params.transparency;this._origBackAmbientIntensity=params.backAmbientIntensity;this._origBackDiffuseColor=params.backDiffuseColor;this._origBackEmissiveColor=params.backEmissiveColor;this._origBackShininess=params.backShininess;this._origBackSpecularColor=params.backSpecularColor;this._origBackTransparency=params.backTransparency;this._ambientIntensity=params.ambientIntensity;this._diffuseColor=params.diffuseColor;this._emissiveColor=params.emissiveColor;this._shininess=params.shininess;this._specularColor=params.specularColor;this._transparency=params.transparency;this._backAmbientIntensity=params.backAmbientIntensity;this._backDiffuseColor=params.backDiffuseColor;this._backEmissiveColor=params.backEmissiveColor;this._backShininess=params.backShininess;this._backSpecularColor=params.backSpecularColor;this._backTransparency=params.backTransparency;this._highlighted=false;this.reset=function(){this._ambientIntensity=this._origAmbientIntensity;this._diffuseColor=this._origDiffuseColor;this._emissiveColor=this._origEmissiveColor;this._shininess=this._origShininess;this._specularColor=this._origSpeclarColor;this._transparency=this._origTransparency;this._backAmbientIntensity=this._origBackAmbientIntensity;this._backDiffuseColor=this._origBackDiffuseColor;this._backEmissiveColor=this._origBackEmissiveColor;this._backShininess=this._origBackShininess;this._backSpecularColor=this._origBackSpecularColor;this._backTransparency=this._origBackTransparency;};};x3dom.Parts=function(multiPart,ids,colorMap,emissiveMap,specularMap,visibilityMap)
{var parts=this;this.multiPart=multiPart;this.ids=ids;this.colorMap=colorMap;this.emissiveMap=emissiveMap;this.specularMap=specularMap;this.visibilityMap=visibilityMap;this.width=parts.colorMap.getWidth();this.widthTwo=this.width*this.width;this.setDiffuseColor=function(color,side)
{var i,partID,pixelIDFront,pixelIDBack;if(side==undefined&&side!="front"&&side!="back"&&side!="both"){side="both";}
color=x3dom.fields.SFColor.parse(color);if(ids.length&&ids.length>1)
{var pixels=parts.colorMap.getPixels();for(i=0;i<parts.ids.length;i++)
{partID=parts.ids[i];pixelIDFront=partID;pixelIDBack=partID+this.widthTwo;if(side=="front")
{this.multiPart._materials[partID]._diffuseColor=color;}
else if(side=="back")
{this.multiPart._materials[partID]._backDiffuseColor=color;}
else if(side=="both")
{this.multiPart._materials[partID]._diffuseColor=color;this.multiPart._materials[partID]._backDiffuseColor=color;}
if(!this.multiPart._materials[partID]._highlighted)
{if(side=="front"){pixels[pixelIDFront].r=color.r;pixels[pixelIDFront].g=color.g;pixels[pixelIDFront].b=color.b;}else if(side=="back"){pixels[pixelIDBack].r=color.r;pixels[pixelIDBack].g=color.g;pixels[pixelIDBack].b=color.b;}else if(side=="both"){pixels[pixelIDFront].r=color.r;pixels[pixelIDFront].g=color.g;pixels[pixelIDFront].b=color.b;pixels[pixelIDBack].r=color.r;pixels[pixelIDBack].g=color.g;pixels[pixelIDBack].b=color.b;}}}
parts.colorMap.setPixels(pixels);}
else
{var xFront,yFront,xBack,yBack,pixelFront,pixelBack;partID=parts.ids[0];pixelIDFront=partID;pixelIDBack=partID+this.widthTwo;if(side=="front")
{xFront=pixelIDFront%this.width;yFront=Math.floor(pixelIDFront/this.width);pixelFront=parts.colorMap.getPixel(xFront,yFront);this.multiPart._materials[partID]._diffuseColor=color;}
else if(side=="back")
{xBack=pixelIDBack%this.width;yBack=Math.floor(pixelIDBack/this.width);pixelBack=parts.colorMap.getPixel(xBack,yBack);this.multiPart._materials[partID]._backDiffuseColor=color;}
else if(side=="both")
{xFront=pixelIDFront%this.width;yFront=Math.floor(pixelIDFront/this.width);xBack=pixelIDBack%this.width;yBack=Math.floor(pixelIDBack/this.width);pixelFront=parts.colorMap.getPixel(xFront,yFront);pixelBack=parts.colorMap.getPixel(xBack,yBack);this.multiPart._materials[partID]._diffuseColor=color;this.multiPart._materials[partID]._backDiffuseColor=color;}
if(!this.multiPart._materials[partID]._highlighted)
{if(side=="front")
{pixelFront.r=color.r;pixelFront.g=color.g;pixelFront.b=color.b;parts.colorMap.setPixel(xFront,yFront,pixelFront);}
else if(side=="back")
{pixelBack.r=color.r;pixelBack.g=color.g;pixelBack.b=color.b;parts.colorMap.setPixel(xBack,yBack,pixelBack);}
else if(side=="both")
{pixelFront.r=color.r;pixelFront.g=color.g;pixelFront.b=color.b;pixelBack.r=color.r;pixelBack.g=color.g;pixelBack.b=color.b;parts.colorMap.setPixel(xFront,yFront,pixelFront);parts.colorMap.setPixel(xBack,yBack,pixelBack);}}}};this.getDiffuseColor=function(side)
{var i,partID;if(side==undefined&&side!="front"&&side!="back"){side="front";}
if(ids.length&&ids.length>1)
{var diffuseColors=[];for(i=0;i<parts.ids.length;i++)
{partID=parts.ids[i];if(side=="front")
{diffuseColors.push(this.multiPart._materials[partID]._diffuseColor);}
else if(side=="back")
{diffuseColors.push(this.multiPart._materials[partID]._backDiffuseColor);}}
return diffuseColors;}
else
{partID=parts.ids[0];if(side=="front")
{return this.multiPart._materials[partID]._diffuseColor;}
else if(side=="back")
{return this.multiPart._materials[partID]._backDiffuseColor;}}};this.setEmissiveColor=function(color,side)
{var i,partID,pixelIDFront,pixelIDBack;if(side==undefined&&side!="front"&&side!="back"&&side!="both"){side="both";}
color=x3dom.fields.SFColor.parse(color);if(ids.length&&ids.length>1)
{var pixels=parts.emissiveMap.getPixels();for(i=0;i<parts.ids.length;i++)
{partID=parts.ids[i];pixelIDFront=partID;pixelIDBack=partID+this.widthTwo;if(side=="front")
{this.multiPart._materials[partID]._emissiveColor=color;}
else if(side=="back")
{this.multiPart._materials[partID]._backEmissiveColor=color;}
else if(side=="both")
{this.multiPart._materials[partID]._emissiveColor=color;this.multiPart._materials[partID]._backEmissiveColor=color;}
if(!this.multiPart._materials[partID]._highlighted)
{if(side=="front"){pixels[pixelIDFront].r=color.r;pixels[pixelIDFront].g=color.g;pixels[pixelIDFront].b=color.b;}else if(side=="back"){pixels[pixelIDBack].r=color.r;pixels[pixelIDBack].g=color.g;pixels[pixelIDBack].b=color.b;}else if(side=="both"){pixels[pixelIDFront].r=color.r;pixels[pixelIDFront].g=color.g;pixels[pixelIDFront].b=color.b;pixels[pixelIDBack].r=color.r;pixels[pixelIDBack].g=color.g;pixels[pixelIDBack].b=color.b;}}}
parts.emissiveMap.setPixels(pixels);}
else
{var xFront,yFront,xBack,yBack,pixelFront,pixelBack;partID=parts.ids[0];pixelIDFront=partID;pixelIDBack=partID+this.widthTwo;if(side=="front")
{xFront=pixelIDFront%this.width;yFront=Math.floor(pixelIDFront/this.width);pixelFront=parts.emissiveMap.getPixel(xFront,yFront);this.multiPart._materials[partID]._emissiveColor=color;}
else if(side=="back")
{xBack=pixelIDBack%this.width;yBack=Math.floor(pixelIDBack/this.width);pixelBack=parts.emissiveMap.getPixel(xBack,yBack);this.multiPart._materials[partID]._backEmissiveColor=color;}
else if(side=="both")
{xFront=pixelIDFront%this.width;yFront=Math.floor(pixelIDFront/this.width);xBack=pixelIDBack%this.width;yBack=Math.floor(pixelIDBack/this.width);pixelFront=parts.emissiveMap.getPixel(xFront,yFront);pixelBack=parts.emissiveMap.getPixel(xBack,yBack);this.multiPart._materials[partID]._emissiveColor=color;this.multiPart._materials[partID]._backEmissiveColor=color;}
if(!this.multiPart._materials[partID]._highlighted)
{if(side=="front")
{pixelFront.r=color.r;pixelFront.g=color.g;pixelFront.b=color.b;parts.emissiveMap.setPixel(xFront,yFront,pixelFront);}
else if(side=="back")
{pixelBack.r=color.r;pixelBack.g=color.g;pixelBack.b=color.b;parts.emissiveMap.setPixel(xBack,yBack,pixelBack);}
else if(side=="both")
{pixelFront.r=color.r;pixelFront.g=color.g;pixelFront.b=color.b;pixelBack.r=color.r;pixelBack.g=color.g;pixelBack.b=color.b;parts.emissiveMap.setPixel(xFront,yFront,pixelFront);parts.emissiveMap.setPixel(xBack,yBack,pixelBack);}}}};this.getEmissiveColor=function(side)
{var i,partID;if(side==undefined&&side!="front"&&side!="back"){side="front";}
if(ids.length&&ids.length>1)
{var emissiveColors=[];for(i=0;i<parts.ids.length;i++)
{partID=parts.ids[i];if(side=="front")
{emissiveColors.push(this.multiPart._materials[partID]._emissiveColor);}
else if(side=="back")
{emissiveColors.push(this.multiPart._materials[partID]._backEmissiveColor);}}
return emissiveColors;}
else
{partID=parts.ids[0];if(side=="front")
{return this.multiPart._materials[partID]._emissiveColor;}
else if(side=="back")
{return this.multiPart._materials[partID]._backEmissiveColor;}}};this.setSpecularColor=function(color,side)
{var i,partID,pixelIDFront,pixelIDBack;if(side==undefined&&side!="front"&&side!="back"&&side!="both"){side="both";}
color=x3dom.fields.SFColor.parse(color);if(ids.length&&ids.length>1)
{var pixels=parts.specularMap.getPixels();for(i=0;i<parts.ids.length;i++)
{partID=parts.ids[i];pixelIDFront=partID;pixelIDBack=partID+this.widthTwo;if(side=="front")
{this.multiPart._materials[partID]._specularColor=color;}
else if(side=="back")
{this.multiPart._materials[partID]._backSpecularColor=color;}
else if(side=="both")
{this.multiPart._materials[partID]._specularColor=color;this.multiPart._materials[partID]._backSpecularColor=color;}
if(!this.multiPart._materials[partID]._highlighted)
{if(side=="front"){pixels[pixelIDFront].r=color.r;pixels[pixelIDFront].g=color.g;pixels[pixelIDFront].b=color.b;}else if(side=="back"){pixels[pixelIDBack].r=color.r;pixels[pixelIDBack].g=color.g;pixels[pixelIDBack].b=color.b;}else if(side=="both"){pixels[pixelIDFront].r=color.r;pixels[pixelIDFront].g=color.g;pixels[pixelIDFront].b=color.b;pixels[pixelIDBack].r=color.r;pixels[pixelIDBack].g=color.g;pixels[pixelIDBack].b=color.b;}}}
parts.specularMap.setPixels(pixels);}
else
{var xFront,yFront,xBack,yBack,pixelFront,pixelBack;partID=parts.ids[0];pixelIDFront=partID;pixelIDBack=partID+this.widthTwo;if(side=="front")
{xFront=pixelIDFront%this.width;yFront=Math.floor(pixelIDFront/this.width);pixelFront=parts.specularMap.getPixel(xFront,yFront);this.multiPart._materials[partID]._specularColor=color;}
else if(side=="back")
{xBack=pixelIDBack%this.width;yBack=Math.floor(pixelIDBack/this.width);pixelBack=parts.specularMap.getPixel(xBack,yBack);this.multiPart._materials[partID]._backSpecularColor=color;}
else if(side=="both")
{xFront=pixelIDFront%this.width;yFront=Math.floor(pixelIDFront/this.width);xBack=pixelIDBack%this.width;yBack=Math.floor(pixelIDBack/this.width);pixelFront=parts.specularMap.getPixel(xFront,yFront);pixelBack=parts.specularMap.getPixel(xBack,yBack);this.multiPart._materials[partID]._specularColor=color;this.multiPart._materials[partID]._backSpecularColor=color;}
if(!this.multiPart._materials[partID]._highlighted)
{if(side=="front")
{pixelFront.r=color.r;pixelFront.g=color.g;pixelFront.b=color.b;parts.specularMap.setPixel(xFront,yFront,pixelFront);}
else if(side=="back")
{pixelBack.r=color.r;pixelBack.g=color.g;pixelBack.b=color.b;parts.specularMap.setPixel(xBack,yBack,pixelBack);}
else if(side=="both")
{pixelFront.r=color.r;pixelFront.g=color.g;pixelFront.b=color.b;pixelBack.r=color.r;pixelBack.g=color.g;pixelBack.b=color.b;parts.specularMap.setPixel(xFront,yFront,pixelFront);parts.specularMap.setPixel(xBack,yBack,pixelBack);}}}};this.getSpecularColor=function(side)
{var i,partID;if(side==undefined&&side!="front"&&side!="back"){side="front";}
if(ids.length&&ids.length>1)
{var specularColors=[];for(i=0;i<parts.ids.length;i++)
{partID=parts.ids[i];if(side=="front")
{specularColors.push(this.multiPart._materials[partID]._specularColor);}
else if(side=="back")
{specularColors.push(this.multiPart._materials[partID]._backSpecularColor);}}
return specularColors;}
else
{partID=parts.ids[0];if(side=="front")
{return this.multiPart._materials[partID]._specularColor;}
else if(side=="back")
{return this.multiPart._materials[partID]._backSpecularColor;}}};this.setTransparency=function(transparency,side)
{var i,partID,pixelIDFront,pixelIDBack;if(side==undefined&&side!="front"&&side!="back"&&side!="both"){side="both";}
if(ids.length&&ids.length>1)
{var pixels=parts.colorMap.getPixels();for(i=0;i<parts.ids.length;i++)
{partID=parts.ids[i];pixelIDFront=partID;pixelIDBack=partID+this.widthTwo;if(side=="front")
{this.multiPart._materials[partID]._transparency=transparency;}
else if(side=="back")
{this.multiPart._materials[partID]._backTransparency=transparency;}
else if(side=="both")
{this.multiPart._materials[partID]._transparency=transparency;this.multiPart._materials[partID]._backTransparency=transparency;}
if(!this.multiPart._materials[partID]._highlighted)
{if(side=="front"){pixels[pixelIDFront].a=1.0-transparency;}else if(side=="back"){pixels[pixelIDBack].a=1.0-transparency;}else if(side=="both"){pixels[pixelIDFront].a=1.0-transparency;pixels[pixelIDBack].a=1.0-transparency;}}}
parts.colorMap.setPixels(pixels);}
else
{var xFront,yFront,xBack,yBack,pixelFront,pixelBack;partID=parts.ids[0];pixelIDFront=partID;pixelIDBack=partID+this.widthTwo;if(side=="front")
{xFront=pixelIDFront%this.width;yFront=Math.floor(pixelIDFront/this.width);pixelFront=parts.colorMap.getPixel(xFront,yFront);this.multiPart._materials[partID]._transparency=transparency;}
else if(side=="back")
{xBack=pixelIDBack%this.width;yBack=Math.floor(pixelIDBack/this.width);pixelBack=parts.colorMap.getPixel(xBack,yBack);this.multiPart._materials[partID]._backTransparency=transparency;}
else if(side=="both")
{xFront=pixelIDFront%this.width;yFront=Math.floor(pixelIDFront/this.width);xBack=pixelIDBack%this.width;yBack=Math.floor(pixelIDBack/this.width);pixelFront=parts.colorMap.getPixel(xFront,yFront);pixelBack=parts.colorMap.getPixel(xBack,yBack);this.multiPart._materials[partID]._transparency=transparency;this.multiPart._materials[partID]._backTransparency=transparency;}
if(!this.multiPart._materials[partID]._highlighted)
{if(side=="front")
{pixelFront.a=1.0-transparency;parts.colorMap.setPixel(xFront,yFront,pixelFront);}
else if(side=="back")
{pixelBack.a=1.0-transparency;parts.colorMap.setPixel(xBack,yBack,pixelBack);}
else if(side=="both")
{pixelFront.a=1.0-transparency;pixelBack.a=1.0-transparency;parts.colorMap.setPixel(xFront,yFront,pixelFront);parts.colorMap.setPixel(xBack,yBack,pixelBack);}}}};this.getTransparency=function(side)
{var i,partID;if(side==undefined&&side!="front"&&side!="back"){side="front";}
if(ids.length&&ids.length>1)
{var transparencies=[];for(i=0;i<parts.ids.length;i++)
{partID=parts.ids[i];if(side=="front")
{transparencies.push(this.multiPart._materials[partID]._transparency);}
else if(side=="back")
{transparencies.push(this.multiPart._materials[partID]._backTransparency);}}
return transparencies;}
else
{partID=parts.ids[0];if(side=="front")
{return this.multiPart._materials[partID]._transparency;}
else if(side=="back")
{return this.multiPart._materials[partID]._backTransparency;}}};this.setShininess=function(shininess,side)
{var i,partID,pixelIDFront,pixelIDBack;if(side==undefined&&side!="front"&&side!="back"&&side!="both"){side="both";}
if(ids.length&&ids.length>1)
{var pixels=parts.specularMap.getPixels();for(i=0;i<parts.ids.length;i++)
{partID=parts.ids[i];pixelIDFront=partID;pixelIDBack=partID+this.widthTwo;if(side=="front")
{this.multiPart._materials[partID]._shininess=shininess;}
else if(side=="back")
{this.multiPart._materials[partID]._backShininess=shininess;}
else if(side=="both")
{this.multiPart._materials[partID]._shininess=shininess;this.multiPart._materials[partID]._backShininess=shininess;}
if(!this.multiPart._materials[partID]._highlighted)
{if(side=="front"){pixels[pixelIDFront].a=shininess;}else if(side=="back"){pixels[pixelIDBack].a=shininess;}else if(side=="both"){pixels[pixelIDFront].a=shininess;pixels[pixelIDBack].a=shininess;}}}
parts.specularMap.setPixels(pixels);}
else
{var xFront,yFront,xBack,yBack,pixelFront,pixelBack;partID=parts.ids[0];pixelIDFront=partID;pixelIDBack=partID+this.widthTwo;if(side=="front")
{xFront=pixelIDFront%this.width;yFront=Math.floor(pixelIDFront/this.width);pixelFront=parts.specularMap.getPixel(xFront,yFront);this.multiPart._materials[partID]._shininess=shininess;}
else if(side=="back")
{xBack=pixelIDBack%this.width;yBack=Math.floor(pixelIDBack/this.width);pixelBack=parts.specularMap.getPixel(xBack,yBack);this.multiPart._materials[partID]._backShininess=shininess;}
else if(side=="both")
{xFront=pixelIDFront%this.width;yFront=Math.floor(pixelIDFront/this.width);xBack=pixelIDBack%this.width;yBack=Math.floor(pixelIDBack/this.width);pixelFront=parts.specularMap.getPixel(xFront,yFront);pixelBack=parts.specularMap.getPixel(xBack,yBack);this.multiPart._materials[partID]._shininess=shininess;this.multiPart._materials[partID]._backShininess=shininess;}
if(!this.multiPart._materials[partID]._highlighted)
{if(side=="front")
{pixelFront.a=shininess;parts.specularMap.setPixel(xFront,yFront,pixelFront);}
else if(side=="back")
{pixelBack.a=shininess;parts.specularMap.setPixel(xBack,yBack,pixelBack);}
else if(side=="both")
{pixelFront.a=shininess;pixelBack.a=shininess;parts.specularMap.setPixel(xFront,yFront,pixelFront);parts.specularMap.setPixel(xBack,yBack,pixelBack);}}}};this.getShininess=function(side)
{var i,partID;if(side==undefined&&side!="front"&&side!="back"){side="front";}
if(ids.length&&ids.length>1)
{var shininesses=[];for(i=0;i<parts.ids.length;i++)
{partID=parts.ids[i];if(side=="front")
{shininesses.push(this.multiPart._materials[partID]._shininess);}
else if(side=="back")
{shininesses.push(this.multiPart._materials[partID]._backShininess);}}
return shininesses;}
else
{partID=parts.ids[0];if(side=="front")
{return this.multiPart._materials[partID]._shininess;}
else if(side=="back")
{return this.multiPart._materials[partID]._backShininess;}}};this.setAmbientIntensity=function(ambientIntensity,side)
{var i,partID,pixelIDFront,pixelIDBack;if(side==undefined&&side!="front"&&side!="back"&&side!="both"){side="both";}
if(ids.length&&ids.length>1)
{var pixels=parts.emissiveMap.getPixels();for(i=0;i<parts.ids.length;i++)
{partID=parts.ids[i];pixelIDFront=partID;pixelIDBack=partID+this.widthTwo;if(side=="front")
{this.multiPart._materials[partID]._ambientIntensity=ambientIntensity;}
else if(side=="back")
{this.multiPart._materials[partID]._backAmbientIntensity=ambientIntensity;}
else if(side=="both")
{this.multiPart._materials[partID]._ambientIntensity=ambientIntensity;this.multiPart._materials[partID]._backAmbientIntensity=ambientIntensity;}
if(!this.multiPart._materials[partID]._highlighted)
{if(side=="front"){pixels[pixelIDFront].a=ambientIntensity;}else if(side=="back"){pixels[pixelIDBack].a=ambientIntensity;}else if(side=="both"){pixels[pixelIDFront].a=ambientIntensity;pixels[pixelIDBack].a=ambientIntensity;}}}
parts.emissiveMap.setPixels(pixels);}
else
{var xFront,yFront,xBack,yBack,pixelFront,pixelBack;partID=parts.ids[0];pixelIDFront=partID;pixelIDBack=partID+this.widthTwo;if(side=="front")
{xFront=pixelIDFront%this.width;yFront=Math.floor(pixelIDFront/this.width);pixelFront=parts.emissiveMap.getPixel(xFront,yFront);this.multiPart._materials[partID]._ambientIntensity=ambientIntensity;}
else if(side=="back")
{xBack=pixelIDBack%this.width;yBack=Math.floor(pixelIDBack/this.width);pixelBack=parts.emissiveMap.getPixel(xBack,yBack);this.multiPart._materials[partID]._backAmbientIntensity=ambientIntensity;}
else if(side=="both")
{xFront=pixelIDFront%this.width;yFront=Math.floor(pixelIDFront/this.width);xBack=pixelIDBack%this.width;yBack=Math.floor(pixelIDBack/this.width);pixelFront=parts.emissiveMap.getPixel(xFront,yFront);pixelBack=parts.emissiveMap.getPixel(xBack,yBack);this.multiPart._materials[partID]._ambientIntensity=ambientIntensity;this.multiPart._materials[partID]._backAmbientIntensity=ambientIntensity;}
if(!this.multiPart._materials[partID]._highlighted)
{if(side=="front")
{pixelFront.a=ambientIntensity;parts.emissiveMap.setPixel(xFront,yFront,pixelFront);}
else if(side=="back")
{pixelBack.a=ambientIntensity;parts.emissiveMap.setPixel(xBack,yBack,pixelBack);}
else if(side=="both")
{pixelFront.a=ambientIntensity;pixelBack.a=ambientIntensity;parts.emissiveMap.setPixel(xFront,yFront,pixelFront);parts.emissiveMap.setPixel(xBack,yBack,pixelBack);}}}};this.getAmbientIntensity=function(side)
{var i,partID;if(side==undefined&&side!="front"&&side!="back"){side="front";}
if(ids.length&&ids.length>1)
{var ambientIntensities=[];for(i=0;i<parts.ids.length;i++)
{partID=parts.ids[i];if(side=="front")
{ambientIntensities.push(this.multiPart._materials[partID]._ambientIntensity);}
else if(side=="back")
{ambientIntensities.push(this.multiPart._materials[partID]._backAmbientIntensity);}}
return ambientIntensities;}
else
{partID=parts.ids[0];if(side=="front")
{return this.multiPart._materials[partID]._ambientIntensity;}
else if(side=="back")
{return this.multiPart._materials[partID]._backAmbientIntensity;}}};this.highlight=function(color)
{var i,partID,pixelIDFront,pixelIDBack,dtColor,eaColor,ssColor;color=x3dom.fields.SFColor.parse(color);if(ids.length&&ids.length>1)
{var dtPixels=parts.colorMap.getPixels();var eaPixels=parts.emissiveMap.getPixels();var ssPixels=parts.specularMap.getPixels();dtColor=new x3dom.fields.SFColorRGBA(0,0,0,1.0);eaColor=new x3dom.fields.SFColorRGBA(color.r,color.g,color.b,0);ssColor=new x3dom.fields.SFColorRGBA(0,0,0,0);for(i=0;i<parts.ids.length;i++){partID=parts.ids[i];pixelIDFront=partID;pixelIDBack=(parseInt(partID)+parseInt(this.widthTwo)).toString();if(!this.multiPart._materials[partID]._highlighted)
{this.multiPart._materials[partID]._highlighted=true;dtPixels[pixelIDFront]=dtColor;eaPixels[pixelIDFront]=eaColor;ssPixels[pixelIDFront]=ssColor;dtPixels[pixelIDBack]=dtColor;eaPixels[pixelIDBack]=eaColor;ssPixels[pixelIDBack]=ssColor;}}
this.colorMap.setPixels(dtPixels,false);this.emissiveMap.setPixels(eaPixels,false);this.specularMap.setPixels(ssPixels,true);}
else
{partID=parts.ids[0];pixelIDFront=partID;pixelIDBack=partID+this.widthTwo;var xFront=pixelIDFront%this.width;var yFront=Math.floor(pixelIDFront/this.width);var xBack=pixelIDBack%this.width;var yBack=Math.floor(pixelIDBack/this.width);if(!this.multiPart._materials[partID]._highlighted)
{this.multiPart._materials[partID]._highlighted=true;dtColor=new x3dom.fields.SFColorRGBA(0,0,0,1);eaColor=new x3dom.fields.SFColorRGBA(color.r,color.g,color.b,0);ssColor=new x3dom.fields.SFColorRGBA(0,0,0,0);this.colorMap.setPixel(xFront,yFront,dtColor,false);this.emissiveMap.setPixel(xFront,yFront,eaColor,false);this.specularMap.setPixel(xFront,yFront,ssColor,false);this.colorMap.setPixel(xBack,yBack,dtColor,false);this.emissiveMap.setPixel(xBack,yBack,eaColor,false);this.specularMap.setPixel(xBack,yBack,ssColor,true);}}};this.unhighlight=function(){var i,partID,pixelIDFront,pixelIDBack,material;var dtColorFront,eaColorFront,ssColorFront;var dtColorBack,eaColorBack,ssColorBack;if(ids.length&&ids.length>1)
{var dtPixels=parts.colorMap.getPixels();var eaPixels=parts.emissiveMap.getPixels();var ssPixels=parts.specularMap.getPixels();for(i=0;i<parts.ids.length;i++){partID=parts.ids[i];pixelIDFront=partID;pixelIDBack=partID+this.widthTwo;material=this.multiPart._materials[partID];if(material._highlighted)
{material._highlighted=false;dtPixels[pixelIDFront]=new x3dom.fields.SFColorRGBA(material._diffuseColor.r,material._diffuseColor.g,material._diffuseColor.b,1.0-material._transparency);eaPixels[pixelIDFront]=new x3dom.fields.SFColorRGBA(material._emissiveColor.r,material._emissiveColor.g,material._emissiveColor.b,material._ambientIntensity);ssPixels[pixelIDFront]=new x3dom.fields.SFColorRGBA(material._specularColor.r,material._specularColor.g,material._specularColor.b,material._shininess);dtPixels[pixelIDBack]=new x3dom.fields.SFColorRGBA(material._backDiffuseColor.r,material._backDiffuseColor.g,material._backDiffuseColor.b,1.0-material._backTransparency);eaPixels[pixelIDBack]=new x3dom.fields.SFColorRGBA(material._backEmissiveColor.r,material._backEmissiveColor.g,material._backEmissiveColor.b,material._backAmbientIntensity);ssPixels[pixelIDBack]=new x3dom.fields.SFColorRGBA(material._backSpecularColor.r,material._backSpecularColor.g,material._backSpecularColor.b,material._backShininess);}}
this.colorMap.setPixels(dtPixels,false);this.emissiveMap.setPixels(eaPixels,false);this.specularMap.setPixels(ssPixels,true);}
else
{partID=parts.ids[0];pixelIDFront=partID;pixelIDBack=partID+this.widthTwo;var xFront=pixelIDFront%this.width;var yFront=Math.floor(pixelIDFront/this.width);var xBack=pixelIDBack%this.width;var yBack=Math.floor(pixelIDBack/this.width);material=this.multiPart._materials[partID];if(material._highlighted)
{material._highlighted=false;dtColorFront=new x3dom.fields.SFColorRGBA(material._diffuseColor.r,material._diffuseColor.g,material._diffuseColor.b,1.0-material._transparency);eaColorFront=new x3dom.fields.SFColorRGBA(material._emissiveColor.r,material._emissiveColor.g,material._emissiveColor.b,material._ambientIntensity);ssColorFront=new x3dom.fields.SFColorRGBA(material._specularColor.r,material._specularColor.g,material._specularColor.b,material._shininess);dtColorBack=new x3dom.fields.SFColorRGBA(material._backDiffuseColor.r,material._backDiffuseColor.g,material._backDiffuseColor.b,1.0-material._backTransparency);eaColorBack=new x3dom.fields.SFColorRGBA(material._backEmissiveColor.r,material._backEmissiveColor.g,material._backEmissiveColor.b,material._backAmbientIntensity);ssColorBack=new x3dom.fields.SFColorRGBA(material._backSpecularColor.r,material._backSpecularColor.g,material._backSpecularColor.b,material._backShininess);this.colorMap.setPixel(xFront,yFront,dtColorFront,false);this.emissiveMap.setPixel(xFront,yFront,eaColorFront,false);this.specularMap.setPixel(xFront,yFront,ssColorFront,false);this.colorMap.setPixel(xBack,yBack,dtColorBack,false);this.emissiveMap.setPixel(xBack,yBack,eaColorBack,false);this.specularMap.setPixel(xBack,yBack,ssColorBack,true);}}};this.toggleHighlight=function(color){for(var i=0;i<parts.ids.length;i++){if(this.multiPart._materials[parts.ids[i]]._highlighted){this.unhighlight();}else{this.highlight(color);}}};this.setColor=function(color,side){this.setDiffuseColor(color,side);};this.getColorRGB=function(){var str=this.getColorRGBA();var values=str.split(" ");return values[0]+" "+values[1]+" "+values[2];};this.getColorRGBA=function(){var x,y;var colorRGBA=this.multiPart._originalColor[parts.ids[0]];if(this.multiPart._highlightedParts[parts.ids[0]]){colorRGBA=this.multiPart._highlightedParts[parts.ids[0]];}else{x=parts.ids[0]%parts.colorMap.getWidth();y=Math.floor(parts.ids[0]/parts.colorMap.getWidth());colorRGBA=parts.colorMap.getPixel(x,y);}
return colorRGBA.toString();};this.resetColor=function(){var i,partID,pixelIDFront,pixelIDBack,material;var dtColorFront,eaColorFront,ssColorFront;var dtColorBack,eaColorBack,ssColorBack;if(ids.length&&ids.length>1)
{var dtPixels=parts.colorMap.getPixels();var eaPixels=parts.emissiveMap.getPixels();var ssPixels=parts.specularMap.getPixels();for(i=0;i<parts.ids.length;i++){partID=parts.ids[i];pixelIDFront=partID;pixelIDBack=partID+this.widthTwo;material=this.multiPart._materials[partID];material.reset();if(!material._highlighted)
{dtPixels[pixelIDFront]=new x3dom.fields.SFColorRGBA(material._diffuseColor.r,material._diffuseColor.g,material._diffuseColor.b,1.0-material._transparency);eaPixels[pixelIDFront]=new x3dom.fields.SFColorRGBA(material._emissiveColor.r,material._emissiveColor.g,material._emissiveColor.b,material._ambientIntensity);ssPixels[pixelIDFront]=new x3dom.fields.SFColorRGBA(material._specularColor.r,material._specularColor.g,material._specularColor.b,material._shininess);dtPixels[pixelIDBack]=new x3dom.fields.SFColorRGBA(material._backDiffuseColor.r,material._backDiffuseColor.g,material._backDiffuseColor.b,1.0-material._backTransparency);eaPixels[pixelIDBack]=new x3dom.fields.SFColorRGBA(material._backEmissiveColor.r,material._backEmissiveColor.g,material._backEmissiveColor.b,material._backAmbientIntensity);ssPixels[pixelIDBack]=new x3dom.fields.SFColorRGBA(material._backSpecularColor.r,material._backSpecularColor.g,material._backSpecularColor.b,material._backShininess);}}
this.colorMap.setPixels(dtPixels,false);this.emissiveMap.setPixels(eaPixels,false);this.specularMap.setPixels(ssPixels,true);}
else
{partID=parts.ids[0];pixelIDFront=partID;pixelIDBack=partID+this.widthTwo;var xFront=pixelIDFront%this.width;var yFront=Math.floor(pixelIDFront/this.width);var xBack=pixelIDBack%this.width;var yBack=Math.floor(pixelIDBack/this.width);material=this.multiPart._materials[partID];material.reset();if(!material._highlighted)
{dtColorFront=new x3dom.fields.SFColorRGBA(material._diffuseColor.r,material._diffuseColor.g,material._diffuseColor.b,1.0-material._transparency);eaColorFront=new x3dom.fields.SFColorRGBA(material._emissiveColor.r,material._emissiveColor.g,material._emissiveColor.b,material._ambientIntensity);ssColorFront=new x3dom.fields.SFColorRGBA(material._specularColor.r,material._specularColor.g,material._specularColor.b,material._shininess);dtColorBack=new x3dom.fields.SFColorRGBA(material._backDiffuseColor.r,material._backDiffuseColor.g,material._backDiffuseColor.b,1.0-material._backTransparency);eaColorBack=new x3dom.fields.SFColorRGBA(material._backEmissiveColor.r,material._backEmissiveColor.g,material._backEmissiveColor.b,material._backAmbientIntensity);ssColorBack=new x3dom.fields.SFColorRGBA(material._backSpecularColor.r,material._backSpecularColor.g,material._backSpecularColor.b,material._backShininess);this.colorMap.setPixel(xFront,yFront,dtColorFront,false);this.emissiveMap.setPixel(xFront,yFront,eaColorFront,false);this.specularMap.setPixel(xFront,yFront,ssColorFront,false);this.colorMap.setPixel(xBack,yBack,dtColorBack,false);this.emissiveMap.setPixel(xBack,yBack,eaColorBack,false);this.specularMap.setPixel(xBack,yBack,ssColorBack,true);}}};this.setVisibility=function(visibility){var i,j,x,y,usage,visibleCount,visibilityAsInt;if(!(ids.length&&ids.length>1)){x=parts.ids[0]%parts.colorMap.getWidth();y=Math.floor(parts.ids[0]/parts.colorMap.getWidth());var pixel=parts.visibilityMap.getPixel(x,y);visibilityAsInt=(visibility)?1:0;if(pixel.r!=visibilityAsInt){pixel.r=visibilityAsInt;this.multiPart._partVisibility[parts.ids[0]]=visibility;usage=this.multiPart._idMap.mapping[parts.ids[0]].usage;for(j=0;j<usage.length;j++){visibleCount=this.multiPart._visiblePartsPerShape[usage[j]];if(visibility&&visibleCount.val<visibleCount.max){visibleCount.val++;}else if(!visibility&&visibleCount.val>0){visibleCount.val--;}
if(visibleCount.val){this.multiPart._inlineNamespace.defMap[usage[j]]._vf.render=true;}else{this.multiPart._inlineNamespace.defMap[usage[j]]._vf.render=false;}}}
parts.visibilityMap.setPixel(x,y,pixel);this.multiPart.invalidateVolume();}
else
{var pixels=parts.visibilityMap.getPixels();for(i=0;i<parts.ids.length;i++){visibilityAsInt=(visibility)?1:0;if(pixels[parts.ids[i]].r!=visibilityAsInt){pixels[parts.ids[i]].r=visibilityAsInt;this.multiPart._partVisibility[parts.ids[i]]=visibility;usage=this.multiPart._idMap.mapping[parts.ids[i]].usage;for(j=0;j<usage.length;j++){visibleCount=this.multiPart._visiblePartsPerShape[usage[j]];if(visibility&&visibleCount.val<visibleCount.max){visibleCount.val++;}else if(!visibility&&visibleCount.val>0){visibleCount.val--;}
if(visibleCount.val){this.multiPart._inlineNamespace.defMap[usage[j]]._vf.render=true;}else{this.multiPart._inlineNamespace.defMap[usage[j]]._vf.render=false;}}}}
parts.visibilityMap.setPixels(pixels);this.multiPart.invalidateVolume();}};this.getVolume=function(){var volume;var transmat=this.multiPart.getCurrentTransform();if(ids.length&&ids.length>1)
{volume=new x3dom.fields.BoxVolume();for(var i=0;i<parts.ids.length;i++){volume.extendBounds(this.multiPart._partVolume[parts.ids[i]].min,this.multiPart._partVolume[parts.ids[i]].max);}
volume.transform(transmat);return volume;}
else
{volume=x3dom.fields.BoxVolume.copy(this.multiPart._partVolume[parts.ids[0]]);volume.transform(transmat);return volume;}};this.fit=function(updateCenterOfRotation){var volume=this.getVolume();this.multiPart._nameSpace.doc._viewarea.fit(volume.min,volume.max,updateCenterOfRotation);};};x3dom.Properties=function(){this.properties={};};x3dom.Properties.prototype.setProperty=function(name,value){x3dom.debug.logInfo("Properties: Setting property '"+name+"' to value '"+value+"'");this.properties[name]=value;};x3dom.Properties.prototype.getProperty=function(name,def){if(this.properties[name]){return this.properties[name]}else{return def;}};x3dom.Properties.prototype.merge=function(other){for(var attrname in other.properties){this.properties[attrname]=other.properties[attrname];}};x3dom.Properties.prototype.toString=function(){var str="";for(var name in this.properties){str+="Name: "+name+" Value: "+this.properties[name]+"\n";}
return str;};x3dom.DoublyLinkedList=function(){this.length=0;this.first=null;this.last=null;};x3dom.DoublyLinkedList.ListNode=function(point,point_index,normals,colors,texCoords){this.point=point;this.point_index=point_index;this.normals=normals;this.colors=colors;this.texCoords=texCoords;this.next=null;this.prev=null;};x3dom.DoublyLinkedList.prototype.appendNode=function(node){if(this.first===null){node.prev=node;node.next=node;this.first=node;this.last=node;}else{node.prev=this.last;node.next=this.first;this.first.prev=node;this.last.next=node;this.last=node;}
this.length++;};x3dom.DoublyLinkedList.prototype.insertAfterNode=function(node,newNode){newNode.prev=node;newNode.next=node.next;node.next.prev=newNode;node.next=newNode;if(newNode.prev==this.last){this.last=newNode;}
this.length++;};x3dom.DoublyLinkedList.prototype.deleteNode=function(node){if(this.length>1){node.prev.next=node.next;node.next.prev=node.prev;if(node==this.first){this.first=node.next;}
if(node==this.last){this.last=node.prev;}}else{this.first=null;this.last=null;}
node.prev=null;node.next=null;this.length--;};x3dom.DoublyLinkedList.prototype.getNode=function(index){var node=null;if(index>this.length){return node;}
for(var i=0;i<this.length;i++){if(i==0){node=this.first;}else{node=node.next;}
if(i==index){return node;}}
return null;};x3dom.DoublyLinkedList.prototype.invert=function(){var tmp=null;var node=this.first;for(var i=0;i<this.length;i++){tmp=node.prev;node.prev=node.next;node.next=tmp;node=node.prev;}
tmp=this.first;this.first=this.last;this.last=tmp;};x3dom.EarClipping={getIndexes:function(linklist){var node=linklist.first.next;var plane=this.identifyPlane(node.prev.point,node.point,node.next.point);var i,points,x,y;points=[];point_indexes=[];for(i=0;i<linklist.length;i++){node=linklist.getNode(i);switch(plane){case"XY":{x=node.point.x;y=node.point.y;break;}
case"XZ":{x=node.point.z;y=node.point.x;break;}
default:{x=node.point.y;y=node.point.z;}}
points.push(y);points.push(x);point_indexes.push(node.point_index);}
var triangles=x3dom.EarCut.triangulate(points,null,2);triangles=triangles.map(function(m){return point_indexes[m];});return triangles;},getMultiIndexes:function(linklist){var node=linklist.first.next;var plane=this.identifyPlane(node.prev.point,node.point,node.next.point);var data={};data.indices=[];data.point=[];data.normals=[];data.colors=[];data.texCoords=[];var mapped={};mapped.indices=[];mapped.point=[];mapped.normals=[];mapped.colors=[];mapped.texCoords=[];points=[];for(i=0;i<linklist.length;i++){node=linklist.getNode(i);switch(plane){case"XY":{x=node.point.x;y=node.point.y;break;}
case"XZ":{x=node.point.z;y=node.point.x;break;}
default:{x=node.point.y;y=node.point.z;}}
points.push(y);points.push(x);mapped.indices.push(node.point_index);mapped.point.push(node.point);if(node.normals)mapped.normals.push(node.normals);if(node.colors)mapped.colors.push(node.colors);if(node.texCoords)mapped.texCoords.push(node.texCoords);}
var triangles=x3dom.EarCut.triangulate(points,null,2);data.indices=triangles.map(function(m){return mapped.indices[m];});data.point=triangles.map(function(m){return mapped.point[m];});if(node.normals)data.normals=triangles.map(function(m){return mapped.normals[m];});if(node.colors)data.colors=triangles.map(function(m){return mapped.colors[m];});if(node.texCoords)data.texCoords=triangles.map(function(m){return mapped.texCoords[m];});return data;},identifyPlane:function(p1,p2,p3){var v1x,v1y,v1z;var v2x,v2y,v2z;var v3x,v3y,v3z;v1x=p2.x-p1.x;v1y=p2.y-p1.y;v1z=p2.z-p1.z;v2x=p3.x-p1.x;v2y=p3.y-p1.y;v2z=p3.z-p1.z;v3x=Math.abs(v1y*v2z-v1z*v2y);v3y=Math.abs(v1z*v2x-v1x*v2z);v3z=Math.abs(v1x*v2y-v1y*v2x);var angle=Math.max(v3x,v3y,v3z);if(angle==v3x){return'YZ';}else if(angle==v3y){return'XZ';}else if(angle==v3z){return'XY';}else{return'XZ';}}};x3dom.EarCut={triangulate:function mapEarcut(data,holes,dim){return earcut(data,holes,dim);function earcut(data,holeIndices,dim){dim=dim||2;var hasHoles=holeIndices&&holeIndices.length,outerLen=hasHoles?holeIndices[0]*dim:data.length,clockwise=windingOrder(data,0,outerLen,dim),outerNode=linkedList(data,0,outerLen,dim,true,clockwise),triangles=[];if(!outerNode)return triangles;var minX,minY,maxX,maxY,x,y,size;if(hasHoles)outerNode=eliminateHoles(data,holeIndices,outerNode,dim);if(data.length>80*dim){minX=maxX=data[0];minY=maxY=data[1];for(var i=dim;i<outerLen;i+=dim){x=data[i];y=data[i+1];if(x<minX)minX=x;if(y<minY)minY=y;if(x>maxX)maxX=x;if(y>maxY)maxY=y;}
size=Math.max(maxX-minX,maxY-minY);}
earcutLinked(outerNode,triangles,dim,minX,minY,size);if(clockwise===false){triangles.reverse();}
return triangles;}
function windingOrder(data,start,end,dim){var sum=0;for(i=start,j=end-dim;i<end;i+=dim){sum+=(data[j]-data[i])*(data[i+1]+data[j+1]);j=i;}
return sum>0;}
function linkedList(data,start,end,dim,clockwise,oclockwise){var i,j,last;if(clockwise===oclockwise){for(i=start;i<end;i+=dim)last=insertNode(i,data[i],data[i+1],last);}else{for(i=end-dim;i>=start;i-=dim)last=insertNode(i,data[i],data[i+1],last);}
return last;}
function filterPoints(start,end){if(!start)return start;if(!end)end=start;var p=start,again;do{again=false;if(!p.steiner&&(equals(p,p.next)||area(p.prev,p,p.next)===0)){removeNode(p);p=end=p.prev;if(p===p.next)return null;again=true;}else{p=p.next;}}while(again||p!==end);return end;}
function earcutLinked(ear,triangles,dim,minX,minY,size,pass){if(!ear)return;if(!pass&&size)indexCurve(ear,minX,minY,size);var stop=ear,prev,next;while(ear.prev!==ear.next){prev=ear.prev;next=ear.next;if(size?isEarHashed(ear,minX,minY,size):isEar(ear)){triangles.push(prev.i/dim);triangles.push(ear.i/dim);triangles.push(next.i/dim);removeNode(ear);ear=next.next;stop=next.next;continue;}
ear=next;if(ear===stop){if(!pass){earcutLinked(filterPoints(ear),triangles,dim,minX,minY,size,1);}else if(pass===1){ear=cureLocalIntersections(ear,triangles,dim);earcutLinked(ear,triangles,dim,minX,minY,size,2);}else if(pass===2){splitEarcut(ear,triangles,dim,minX,minY,size);}
break;}}}
function isEar(ear){var a=ear.prev,b=ear,c=ear.next;if(area(a,b,c)>=0)return false;var p=ear.next.next;while(p!==ear.prev){if(pointInTriangle(a.x,a.y,b.x,b.y,c.x,c.y,p.x,p.y)&&area(p.prev,p,p.next)>=0)return false;p=p.next;}
return true;}
function isEarHashed(ear,minX,minY,size){var a=ear.prev,b=ear,c=ear.next;if(area(a,b,c)>=0)return false;var minTX=a.x<b.x?(a.x<c.x?a.x:c.x):(b.x<c.x?b.x:c.x),minTY=a.y<b.y?(a.y<c.y?a.y:c.y):(b.y<c.y?b.y:c.y),maxTX=a.x>b.x?(a.x>c.x?a.x:c.x):(b.x>c.x?b.x:c.x),maxTY=a.y>b.y?(a.y>c.y?a.y:c.y):(b.y>c.y?b.y:c.y);var minZ=zOrder(minTX,minTY,minX,minY,size),maxZ=zOrder(maxTX,maxTY,minX,minY,size);var p=ear.nextZ;while(p&&p.z<=maxZ){if(p!==ear.prev&&p!==ear.next&&pointInTriangle(a.x,a.y,b.x,b.y,c.x,c.y,p.x,p.y)&&area(p.prev,p,p.next)>=0)return false;p=p.nextZ;}
p=ear.prevZ;while(p&&p.z>=minZ){if(p!==ear.prev&&p!==ear.next&&pointInTriangle(a.x,a.y,b.x,b.y,c.x,c.y,p.x,p.y)&&area(p.prev,p,p.next)>=0)return false;p=p.prevZ;}
return true;}
function cureLocalIntersections(start,triangles,dim){var p=start;do{var a=p.prev,b=p.next.next;if(intersects(a,p,p.next,b)&&locallyInside(a,b)&&locallyInside(b,a)){triangles.push(a.i/dim);triangles.push(p.i/dim);triangles.push(b.i/dim);removeNode(p);removeNode(p.next);p=start=b;}
p=p.next;}while(p!==start);return p;}
function splitEarcut(start,triangles,dim,minX,minY,size){var a=start;do{var b=a.next.next;while(b!==a.prev){if(a.i!==b.i&&isValidDiagonal(a,b)){var c=splitPolygon(a,b);a=filterPoints(a,a.next);c=filterPoints(c,c.next);earcutLinked(a,triangles,dim,minX,minY,size);earcutLinked(c,triangles,dim,minX,minY,size);return;}
b=b.next;}
a=a.next;}while(a!==start);}
function eliminateHoles(data,holeIndices,outerNode,dim){var queue=[],i,len,start,end,list;for(i=0,len=holeIndices.length;i<len;i++){start=holeIndices[i]*dim;end=i<len-1?holeIndices[i+1]*dim:data.length;list=linkedList(data,start,end,dim,false);if(list===list.next)list.steiner=true;queue.push(getLeftmost(list));}
queue.sort(compareX);for(i=0;i<queue.length;i++){eliminateHole(queue[i],outerNode);outerNode=filterPoints(outerNode,outerNode.next);}
return outerNode;}
function compareX(a,b){return a.x-b.x;}
function eliminateHole(hole,outerNode){outerNode=findHoleBridge(hole,outerNode);if(outerNode){var b=splitPolygon(outerNode,hole);filterPoints(b,b.next);}}
function findHoleBridge(hole,outerNode){var p=outerNode,hx=hole.x,hy=hole.y,qx=-Infinity,m;do{if(hy<=p.y&&hy>=p.next.y){var x=p.x+(hy-p.y)*(p.next.x-p.x)/(p.next.y-p.y);if(x<=hx&&x>qx){qx=x;m=p.x<p.next.x?p:p.next;}}
p=p.next;}while(p!==outerNode);if(!m)return null;var stop=m,tanMin=Infinity,tan;p=m.next;while(p!==stop){if(hx>=p.x&&p.x>=m.x&&pointInTriangle(hy<m.y?hx:qx,hy,m.x,m.y,hy<m.y?qx:hx,hy,p.x,p.y)){tan=Math.abs(hy-p.y)/(hx-p.x);if((tan<tanMin||(tan===tanMin&&p.x>m.x))&&locallyInside(p,hole)){m=p;tanMin=tan;}}
p=p.next;}
return m;}
function indexCurve(start,minX,minY,size){var p=start;do{if(p.z===null)p.z=zOrder(p.x,p.y,minX,minY,size);p.prevZ=p.prev;p.nextZ=p.next;p=p.next;}while(p!==start);p.prevZ.nextZ=null;p.prevZ=null;sortLinked(p);}
function sortLinked(list){var i,p,q,e,tail,numMerges,pSize,qSize,inSize=1;do{p=list;list=null;tail=null;numMerges=0;while(p){numMerges++;q=p;pSize=0;for(i=0;i<inSize;i++){pSize++;q=q.nextZ;if(!q)break;}
qSize=inSize;while(pSize>0||(qSize>0&&q)){if(pSize===0){e=q;q=q.nextZ;qSize--;}else if(qSize===0||!q){e=p;p=p.nextZ;pSize--;}else if(p.z<=q.z){e=p;p=p.nextZ;pSize--;}else{e=q;q=q.nextZ;qSize--;}
if(tail)tail.nextZ=e;else list=e;e.prevZ=tail;tail=e;}
p=q;}
tail.nextZ=null;inSize*=2;}while(numMerges>1);return list;}
function zOrder(x,y,minX,minY,size){x=32767*(x-minX)/size;y=32767*(y-minY)/size;x=(x|(x<<8))&0x00FF00FF;x=(x|(x<<4))&0x0F0F0F0F;x=(x|(x<<2))&0x33333333;x=(x|(x<<1))&0x55555555;y=(y|(y<<8))&0x00FF00FF;y=(y|(y<<4))&0x0F0F0F0F;y=(y|(y<<2))&0x33333333;y=(y|(y<<1))&0x55555555;return x|(y<<1);}
function getLeftmost(start){var p=start,leftmost=start;do{if(p.x<leftmost.x)leftmost=p;p=p.next;}while(p!==start);return leftmost;}
function pointInTriangle(ax,ay,bx,by,cx,cy,px,py){return(cx-px)*(ay-py)-(ax-px)*(cy-py)>=0&&(ax-px)*(by-py)-(bx-px)*(ay-py)>=0&&(bx-px)*(cy-py)-(cx-px)*(by-py)>=0;}
function isValidDiagonal(a,b){return equals(a,b)||a.next.i!==b.i&&a.prev.i!==b.i&&!intersectsPolygon(a,b)&&locallyInside(a,b)&&locallyInside(b,a)&&middleInside(a,b);}
function area(p,q,r){return(q.y-p.y)*(r.x-q.x)-(q.x-p.x)*(r.y-q.y);}
function equals(p1,p2){return p1.x===p2.x&&p1.y===p2.y;}
function intersects(p1,q1,p2,q2){return area(p1,q1,p2)>0!==area(p1,q1,q2)>0&&area(p2,q2,p1)>0!==area(p2,q2,q1)>0;}
function intersectsPolygon(a,b){var p=a;do{if(p.i!==a.i&&p.next.i!==a.i&&p.i!==b.i&&p.next.i!==b.i&&intersects(p,p.next,a,b))return true;p=p.next;}while(p!==a);return false;}
function locallyInside(a,b){return area(a.prev,a,a.next)<0?area(a,b,a.next)>=0&&area(a,a.prev,b)>=0:area(a,b,a.prev)<0||area(a,a.next,b)<0;}
function middleInside(a,b){var p=a,inside=false,px=(a.x+b.x)/2,py=(a.y+b.y)/2;do{if(((p.y>py)!==(p.next.y>py))&&(px<(p.next.x-p.x)*(py-p.y)/(p.next.y-p.y)+p.x))
inside=!inside;p=p.next;}while(p!==a);return inside;}
function splitPolygon(a,b){var a2=new Node(a.i,a.x,a.y),b2=new Node(b.i,b.x,b.y),an=a.next,bp=b.prev;a.next=b;b.prev=a;a2.next=an;an.prev=a2;b2.next=a2;a2.prev=b2;bp.next=b2;b2.prev=bp;return b2;}
function insertNode(i,x,y,last){var p=new Node(i,x,y);if(!last){p.prev=p;p.next=p;}else{p.next=last.next;p.prev=last;last.next.prev=p;last.next=p;}
return p;}
function removeNode(p){p.next.prev=p.prev;p.prev.next=p.next;if(p.prevZ)p.prevZ.nextZ=p.nextZ;if(p.nextZ)p.nextZ.prevZ=p.prevZ;}
function Node(i,x,y){this.i=i;this.x=x;this.y=y;this.prev=null;this.next=null;this.z=null;this.prevZ=null;this.nextZ=null;this.steiner=false;}}}
x3dom.FieldInterpolator=function(beginTime,endTime,beginValue,endValue)
{this.beginTime=beginTime||0;this.endTime=endTime||1;this.beginValue=beginValue||0;this.endValue=endValue||0;this.isInterpolating=false;};x3dom.FieldInterpolator.prototype.isActive=function()
{return(this.beginTime>0);};x3dom.FieldInterpolator.prototype.calcFraction=function(time)
{var fraction=(time-this.beginTime)/(this.endTime-this.beginTime);return(Math.sin((fraction*Math.PI)-(Math.PI/2))+1)/2.0;};x3dom.FieldInterpolator.prototype.reset=function()
{this.isInterpolating=false;this.beginTime=0;this.endTime=1;this.beginValue=0;this.endValue=0;};x3dom.FieldInterpolator.prototype.interpolate=function(time)
{if(time<this.beginTime)
{return this.beginValue;}
else if(time>=this.endTime)
{var endValue=this.endValue;this.reset();return endValue;}
else
{this.isInterpolating=true;return this.beginValue+(this.endValue-this.beginValue)*this.calcFraction(time);}};x3dom.Utils={};x3dom.Utils.maxIndexableCoords=65535;x3dom.Utils.needLineWidth=false;x3dom.Utils.measurements=[];window.performance=window.performance||{};performance.now=(function(){return performance.now||performance.mozNow||performance.msNow||performance.oNow||performance.webkitNow||function(){return new Date().getTime();};})();x3dom.Utils.startMeasure=function(name){var uname=name.toUpperCase();if(!x3dom.Utils.measurements[uname]){if(performance&&performance.now){x3dom.Utils.measurements[uname]=performance.now();}else{x3dom.Utils.measurements[uname]=new Date().getTime();}}};x3dom.Utils.stopMeasure=function(name){var uname=name.toUpperCase();if(x3dom.Utils.measurements[uname]){var startTime=x3dom.Utils.measurements[uname];delete x3dom.Utils.measurements[uname];if(performance&&performance.now){return performance.now()-startTime;}else{return new Date().getTime()-startTime;}}
return 0;};x3dom.Utils.isNumber=function(n){return!isNaN(parseFloat(n))&&isFinite(n);};x3dom.Utils.createTexture2D=function(gl,doc,src,bgnd,crossOrigin,scale,genMipMaps)
{var texture=gl.createTexture();var data=new Uint8Array([0,0,0,255,0,0,0,255,0,0,0,255,0,0,0,255]);gl.bindTexture(gl.TEXTURE_2D,texture);gl.texImage2D(gl.TEXTURE_2D,0,gl.RGBA,2,2,0,gl.RGBA,gl.UNSIGNED_BYTE,data);if(genMipMaps){gl.generateMipmap(gl.TEXTURE_2D);}
gl.bindTexture(gl.TEXTURE_2D,null);texture.ready=false;if(src==null||src=='')
return texture;var image=new Image();switch(crossOrigin.toLowerCase()){case'anonymous':{image.crossOrigin='anonymous';}break;case'use-credentials':{image.crossOrigin='use-credentials'}break;case'none':{}break;default:{if(x3dom.Utils.forbiddenBySOP(src)){image.crossOrigin='anonymous';}}}
image.src=src;doc.downloadCount++;image.onload=function(){texture.originalWidth=image.width;texture.originalHeight=image.height;if(scale)
image=x3dom.Utils.scaleImage(image);if(bgnd==true){gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL,true);}
gl.bindTexture(gl.TEXTURE_2D,texture);gl.texImage2D(gl.TEXTURE_2D,0,gl.RGBA,gl.RGBA,gl.UNSIGNED_BYTE,image);if(genMipMaps){gl.generateMipmap(gl.TEXTURE_2D);}
gl.bindTexture(gl.TEXTURE_2D,null);if(bgnd==true){gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL,false);}
texture.width=image.width;texture.height=image.height;texture.ready=true;doc.downloadCount--;doc.needRender=true;};image.onerror=function(error){if(x3dom.caps.EXTENSIONS.indexOf('WEBGL_compressed_texture_s3tc')!==-1){x3dom.Utils.tryCompressedTexture2D(texture,gl,doc,src,bgnd,crossOrigin,genMipMaps,function(success){if(success){}else{x3dom.debug.logError("[Utils|createTexture2D] Can't load Image: "+src);}
doc.downloadCount--;});}else{x3dom.debug.logError("[Utils|createTexture2D] Can't load Image: "+src);doc.downloadCount--;}};return texture;};x3dom.Utils.createCompressedTexture2D=function(gl,doc,src,bgnd,crossOrigin,genMipMaps)
{var texture=gl.createTexture();var data=new Uint8Array([0,0,0,255,0,0,0,255,0,0,0,255,0,0,0,255]);gl.bindTexture(gl.TEXTURE_2D,texture);gl.texImage2D(gl.TEXTURE_2D,0,gl.RGBA,2,2,0,gl.RGBA,gl.UNSIGNED_BYTE,data);if(genMipMaps){gl.generateMipmap(gl.TEXTURE_2D);}
gl.bindTexture(gl.TEXTURE_2D,null);texture.ready=false;if(src==null||src=='')
return texture;ddsXhr=new XMLHttpRequest();var ext=gl.getExtension('WEBGL_compressed_texture_s3tc');ddsXhr.open('GET',src,true);ddsXhr.responseType="arraybuffer";ddsXhr.onload=function(){gl.bindTexture(gl.TEXTURE_2D,texture);var mipmaps=uploadDDSLevels(gl,ext,this.response);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MAG_FILTER,gl.LINEAR);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MIN_FILTER,mipmaps>1?gl.LINEAR_MIPMAP_LINEAR:gl.LINEAR);texture.ready=true;doc.downloadCount--;doc.needRender=true;};doc.downloadCount++;x3dom.RequestManager.addRequest(ddsXhr);return texture;};x3dom.Utils.tryCompressedTexture2D=function(texture,gl,doc,src,bgnd,crossOrigin,genMipMaps,cb)
{ddsXhr=new XMLHttpRequest();var ext=gl.getExtension('WEBGL_compressed_texture_s3tc');ddsXhr.open('GET',src,true);ddsXhr.responseType="arraybuffer";ddsXhr.onload=function(){gl.bindTexture(gl.TEXTURE_2D,texture);var mipmaps=uploadDDSLevels(gl,ext,this.response);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MAG_FILTER,gl.LINEAR);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MIN_FILTER,mipmaps>1?gl.LINEAR_MIPMAP_LINEAR:gl.LINEAR);texture.ready=true;doc.needRender=true;cb(true);};ddsXhr.onerror=function(){cb(false);};x3dom.RequestManager.addRequest(ddsXhr);};function uploadDDSLevels(gl,ext,arrayBuffer,loadMipmaps){var DDS_MAGIC=0x20534444;var DDSD_CAPS=0x1,DDSD_HEIGHT=0x2,DDSD_WIDTH=0x4,DDSD_PITCH=0x8,DDSD_PIXELFORMAT=0x1000,DDSD_MIPMAPCOUNT=0x20000,DDSD_LINEARSIZE=0x80000,DDSD_DEPTH=0x800000;var DDSCAPS_COMPLEX=0x8,DDSCAPS_MIPMAP=0x400000,DDSCAPS_TEXTURE=0x1000;var DDSCAPS2_CUBEMAP=0x200,DDSCAPS2_CUBEMAP_POSITIVEX=0x400,DDSCAPS2_CUBEMAP_NEGATIVEX=0x800,DDSCAPS2_CUBEMAP_POSITIVEY=0x1000,DDSCAPS2_CUBEMAP_NEGATIVEY=0x2000,DDSCAPS2_CUBEMAP_POSITIVEZ=0x4000,DDSCAPS2_CUBEMAP_NEGATIVEZ=0x8000,DDSCAPS2_VOLUME=0x200000;var DDPF_ALPHAPIXELS=0x1,DDPF_ALPHA=0x2,DDPF_FOURCC=0x4,DDPF_RGB=0x40,DDPF_YUV=0x200,DDPF_LUMINANCE=0x20000;function FourCCToInt32(value){return value.charCodeAt(0)+
(value.charCodeAt(1)<<8)+
(value.charCodeAt(2)<<16)+
(value.charCodeAt(3)<<24);}
function Int32ToFourCC(value){return String.fromCharCode(value&0xff,(value>>8)&0xff,(value>>16)&0xff,(value>>24)&0xff);}
var FOURCC_DXT1=FourCCToInt32("DXT1");var FOURCC_DXT5=FourCCToInt32("DXT5");var headerLengthInt=31;var off_magic=0;var off_size=1;var off_flags=2;var off_height=3;var off_width=4;var off_mipmapCount=7;var off_pfFlags=20;var off_pfFourCC=21;var header=new Int32Array(arrayBuffer,0,headerLengthInt),fourCC,blockBytes,internalFormat,width,height,dataLength,dataOffset,byteArray,mipmapCount,i;if(header[off_magic]!=DDS_MAGIC){console.error("Invalid magic number in DDS header");return 0;}
if(!header[off_pfFlags]&DDPF_FOURCC){console.error("Unsupported format, must contain a FourCC code");return 0;}
fourCC=header[off_pfFourCC];switch(fourCC){case FOURCC_DXT1:blockBytes=8;internalFormat=ext.COMPRESSED_RGBA_S3TC_DXT1_EXT;break;case FOURCC_DXT5:blockBytes=16;internalFormat=ext.COMPRESSED_RGBA_S3TC_DXT5_EXT;break;default:console.error("Unsupported FourCC code:",Int32ToFourCC(fourCC));return null;}
mipmapCount=1;if(header[off_flags]&DDSD_MIPMAPCOUNT&&loadMipmaps!==false){mipmapCount=Math.max(1,header[off_mipmapCount]);}
width=header[off_width];height=header[off_height];dataOffset=header[off_size]+4;for(i=0;i<mipmapCount;++i){dataLength=Math.max(4,width)/4*Math.max(4,height)/4*blockBytes;byteArray=new Uint8Array(arrayBuffer,dataOffset,dataLength);gl.compressedTexImage2D(gl.TEXTURE_2D,i,internalFormat,width,height,0,byteArray);dataOffset+=dataLength;width*=0.5;height*=0.5;}
return mipmapCount;};x3dom.Utils.createTextureCube=function(gl,doc,src,bgnd,crossOrigin,scale,genMipMaps)
{var texture=gl.createTexture();var faces;if(bgnd){faces=[gl.TEXTURE_CUBE_MAP_POSITIVE_Z,gl.TEXTURE_CUBE_MAP_NEGATIVE_Z,gl.TEXTURE_CUBE_MAP_POSITIVE_Y,gl.TEXTURE_CUBE_MAP_NEGATIVE_Y,gl.TEXTURE_CUBE_MAP_POSITIVE_X,gl.TEXTURE_CUBE_MAP_NEGATIVE_X];}
else
{faces=[gl.TEXTURE_CUBE_MAP_NEGATIVE_Z,gl.TEXTURE_CUBE_MAP_POSITIVE_Z,gl.TEXTURE_CUBE_MAP_NEGATIVE_Y,gl.TEXTURE_CUBE_MAP_POSITIVE_Y,gl.TEXTURE_CUBE_MAP_NEGATIVE_X,gl.TEXTURE_CUBE_MAP_POSITIVE_X];}
texture.ready=false;texture.pendingTextureLoads=-1;texture.textureCubeReady=false;var width=0,height=0;for(var i=0;i<faces.length;i++){var face=faces[i];var image=new Image();switch(crossOrigin.toLowerCase()){case'anonymous':{image.crossOrigin='anonymous';}break;case'use-credentials':{image.crossOrigin='use-credentials'}break;case'none':{}break;default:{if(x3dom.Utils.forbiddenBySOP(src[i])){image.crossOrigin='anonymous';}}}
texture.pendingTextureLoads++;doc.downloadCount++;image.onload=(function(texture,face,image,swap){return function(){if(width==0&&height==0){width=image.width;height=image.height;}
else if(scale&&(width!=image.width||height!=image.height)){x3dom.debug.logWarning("[Utils|createTextureCube] Rescaling CubeMap images, which are of different size!");image=x3dom.Utils.rescaleImage(image,width,height);}
gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL,swap);gl.bindTexture(gl.TEXTURE_CUBE_MAP,texture);gl.texImage2D(face,0,gl.RGBA,gl.RGBA,gl.UNSIGNED_BYTE,image);gl.bindTexture(gl.TEXTURE_CUBE_MAP,null);gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL,false);texture.pendingTextureLoads--;doc.downloadCount--;if(texture.pendingTextureLoads<0){texture.width=width;texture.height=height;texture.textureCubeReady=true;if(genMipMaps){gl.bindTexture(gl.TEXTURE_CUBE_MAP,texture);gl.generateMipmap(gl.TEXTURE_CUBE_MAP);gl.bindTexture(gl.TEXTURE_CUBE_MAP,null);}
x3dom.debug.logInfo("[Utils|createTextureCube] Loading CubeMap finished...");doc.needRender=true;}};})(texture,face,image,bgnd);image.onerror=function()
{doc.downloadCount--;x3dom.debug.logError("[Utils|createTextureCube] Can't load CubeMap!");};image.src=src[i];}
return texture;};x3dom.Utils.initFBO=function(gl,w,h,type,mipMap,needDepthBuf,numMrt){var tex=gl.createTexture();tex.width=w;tex.height=h;gl.bindTexture(gl.TEXTURE_2D,tex);gl.texImage2D(gl.TEXTURE_2D,0,gl.RGBA,w,h,0,gl.RGBA,type,null);if(mipMap)
gl.generateMipmap(gl.TEXTURE_2D);gl.bindTexture(gl.TEXTURE_2D,null);var i,mrts=null;if(x3dom.caps.DRAW_BUFFERS&&numMrt!==undefined){mrts=[tex];for(i=1;i<numMrt;i++){mrts[i]=gl.createTexture();mrts[i].width=w;mrts[i].height=h;gl.bindTexture(gl.TEXTURE_2D,mrts[i]);gl.texImage2D(gl.TEXTURE_2D,0,gl.RGBA,w,h,0,gl.RGBA,type,null);if(mipMap)
gl.generateMipmap(gl.TEXTURE_2D);gl.bindTexture(gl.TEXTURE_2D,null);}}
var fbo=gl.createFramebuffer();var dtex=null;var rb=null;if(needDepthBuf){if(x3dom.caps.DEPTH_TEXTURE!==null){dtex=gl.createTexture();gl.bindTexture(gl.TEXTURE_2D,dtex);gl.texImage2D(gl.TEXTURE_2D,0,gl.DEPTH_COMPONENT,w,h,0,gl.DEPTH_COMPONENT,gl.UNSIGNED_SHORT,null);if(mipMap)
gl.generateMipmap(gl.TEXTURE_2D);gl.bindTexture(gl.TEXTURE_2D,null);dtex.width=w;dtex.height=h;}
else{rb=gl.createRenderbuffer();gl.bindRenderbuffer(gl.RENDERBUFFER,rb);gl.renderbufferStorage(gl.RENDERBUFFER,gl.DEPTH_COMPONENT16,w,h);gl.bindRenderbuffer(gl.RENDERBUFFER,null);}}
gl.bindFramebuffer(gl.FRAMEBUFFER,fbo);gl.framebufferTexture2D(gl.FRAMEBUFFER,gl.COLOR_ATTACHMENT0,gl.TEXTURE_2D,tex,0);if(x3dom.caps.DRAW_BUFFERS&&numMrt!==undefined){for(i=1;i<numMrt;i++){gl.framebufferTexture2D(gl.FRAMEBUFFER,gl.COLOR_ATTACHMENT0+i,gl.TEXTURE_2D,mrts[i],0);}}
if(needDepthBuf&&x3dom.caps.DEPTH_TEXTURE!==null){gl.framebufferTexture2D(gl.FRAMEBUFFER,gl.DEPTH_ATTACHMENT,gl.TEXTURE_2D,dtex,0);}
else{gl.framebufferRenderbuffer(gl.FRAMEBUFFER,gl.DEPTH_ATTACHMENT,gl.RENDERBUFFER,rb);}
var status=gl.checkFramebufferStatus(gl.FRAMEBUFFER);if(status!=gl.FRAMEBUFFER_COMPLETE){x3dom.debug.logWarning("[Utils|InitFBO] FBO-Status: "+status);}
gl.bindFramebuffer(gl.FRAMEBUFFER,null);return{fbo:fbo,dtex:dtex,rbo:rb,tex:tex,texTargets:mrts,width:w,height:h,type:type,mipMap:mipMap};};x3dom.Utils.getFileName=function(url)
{var filename;if(url.lastIndexOf("/")>-1){filename=url.substr(url.lastIndexOf("/")+1);}
else if(url.lastIndexOf("\\")>-1){filename=url.substr(url.lastIndexOf("\\")+1);}
else{filename=url;}
return filename;};x3dom.Utils.isWebGL2Enabled=function()
{var canvas=document.createElement("canvas");var webgl2=canvas.getContext("webgl2")||canvas.getContext("experimental-webgl2");return(webgl2)?true:false;};x3dom.Utils.findTextureByName=function(texture,name)
{for(var i=0;i<texture.length;++i)
{if(name==texture[i].samplerName)
return texture[i];}
return false;};x3dom.Utils.rescaleImage=function(image,width,height)
{var canvas=document.createElement("canvas");canvas.width=width;canvas.height=height;canvas.getContext("2d").drawImage(image,0,0,image.width,image.height,0,0,canvas.width,canvas.height);return canvas;};x3dom.Utils.scaleImage=function(image)
{if(!x3dom.Utils.isPowerOfTwo(image.width)||!x3dom.Utils.isPowerOfTwo(image.height)){var canvas=document.createElement("canvas");canvas.width=x3dom.Utils.nextHighestPowerOfTwo(image.width);canvas.height=x3dom.Utils.nextHighestPowerOfTwo(image.height);var ctx=canvas.getContext("2d");ctx.drawImage(image,0,0,image.width,image.height,0,0,canvas.width,canvas.height);image=canvas;}
return image;};x3dom.Utils.isPowerOfTwo=function(x)
{return((x&(x-1))===0);};x3dom.Utils.nextHighestPowerOfTwo=function(x)
{--x;for(var i=1;i<32;i<<=1){x=x|x>>i;}
return(x+1);};x3dom.Utils.nextBestPowerOfTwo=function(x)
{var log2x=Math.log(x)/0.693147180559945;return Math.pow(2,Math.round(log2x));};x3dom.Utils.getDataTypeSize=function(type)
{switch(type)
{case"Int8":case"Uint8":return 1;case"Int16":case"Uint16":return 2;case"Int32":case"Uint32":case"Float32":return 4;case"Float64":default:return 8;}};x3dom.Utils.getOffsetMultiplier=function(indexType,gl)
{switch(indexType)
{case gl.UNSIGNED_SHORT:return 1;case gl.UNSIGNED_INT:return 2;case gl.UNSIGNED_BYTE:return 0.5;default:return 1;}};x3dom.Utils.getByteAwareOffset=function(offset,indexType,gl)
{switch(indexType)
{case gl.UNSIGNED_SHORT:return 2*offset;case gl.UNSIGNED_INT:return 4*offset;case gl.UNSIGNED_BYTE:return offset;default:return 2*offset;}};x3dom.Utils.getVertexAttribType=function(type,gl)
{var dataType=gl.NONE;switch(type)
{case"Int8":dataType=gl.BYTE;break;case"Uint8":dataType=gl.UNSIGNED_BYTE;break;case"Int16":dataType=gl.SHORT;break;case"Uint16":dataType=gl.UNSIGNED_SHORT;break;case"Int32":dataType=gl.INT;break;case"Uint32":dataType=gl.UNSIGNED_INT;break;case"Float32":dataType=gl.FLOAT;break;case"Float64":default:x3dom.debug.logError("Can't find this.gl data type for "+type+", getting FLOAT...");dataType=gl.FLOAT;break;}
return dataType;};x3dom.Utils.getArrayBufferView=function(type,buffer)
{var array=null;switch(type)
{case"Int8":array=new Int8Array(buffer);break;case"Uint8":array=new Uint8Array(buffer);break;case"Int16":array=new Int16Array(buffer);break;case"Uint16":array=new Uint16Array(buffer);break;case"Int32":array=new Int32Array(buffer);break;case"Uint32":array=new Uint32Array(buffer);break;case"Float32":array=new Float32Array(buffer);break;case"Float64":array=new Float64Array(buffer);break;default:x3dom.debug.logError("Can't create typed array view of type "+type+", trying Float32...");array=new Float32Array(buffer);break;}
return array;};x3dom.Utils.isUnsignedType=function(str)
{return(str=="Uint8"||str=="Uint16"||str=="Uint16"||str=="Uint32");};x3dom.Utils.checkDirtyLighting=function(viewarea)
{return(viewarea.getLights().length+viewarea._scene.getNavigationInfo()._vf.headlight);};x3dom.Utils.checkDirtyEnvironment=function(viewarea,shaderProperties)
{var environment=viewarea._scene.getEnvironment();return(shaderProperties.GAMMACORRECTION!=environment._vf.gammaCorrectionDefault);};x3dom.Utils.minFilterDic=function(gl,minFilter)
{switch(minFilter.toUpperCase())
{case"NEAREST":return gl.NEAREST;case"LINEAR":return gl.LINEAR;case"NEAREST_MIPMAP_NEAREST":return gl.NEAREST_MIPMAP_NEAREST;case"NEAREST_MIPMAP_LINEAR":return gl.NEAREST_MIPMAP_LINEAR;case"LINEAR_MIPMAP_NEAREST":return gl.LINEAR_MIPMAP_NEAREST;case"LINEAR_MIPMAP_LINEAR":return gl.LINEAR_MIPMAP_LINEAR;case"AVG_PIXEL":return gl.LINEAR;case"AVG_PIXEL_AVG_MIPMAP":return gl.LINEAR_MIPMAP_LINEAR;case"AVG_PIXEL_NEAREST_MIPMAP":return gl.LINEAR_MIPMAP_NEAREST;case"DEFAULT":return gl.LINEAR_MIPMAP_LINEAR;case"FASTEST":return gl.NEAREST;case"NEAREST_PIXEL":return gl.NEAREST;case"NEAREST_PIXEL_AVG_MIPMAP":return gl.NEAREST_MIPMAP_LINEAR;case"NEAREST_PIXEL_NEAREST_MIPMAP":return gl.NEAREST_MIPMAP_NEAREST;case"NICEST":return gl.LINEAR_MIPMAP_LINEAR;default:return gl.LINEAR;}};x3dom.Utils.magFilterDic=function(gl,magFilter)
{switch(magFilter.toUpperCase())
{case"NEAREST":return gl.NEAREST;case"LINEAR":return gl.LINEAR;case"AVG_PIXEL":return gl.LINEAR;case"DEFAULT":return gl.LINEAR;case"FASTEST":return gl.NEAREST;case"NEAREST_PIXEL":return gl.NEAREST;case"NICEST":return gl.LINEAR;default:return gl.LINEAR;}};x3dom.Utils.boundaryModesDic=function(gl,mode)
{switch(mode.toUpperCase())
{case"CLAMP":return gl.CLAMP_TO_EDGE;case"CLAMP_TO_EDGE":return gl.CLAMP_TO_EDGE;case"CLAMP_TO_BOUNDARY":return gl.CLAMP_TO_EDGE;case"MIRRORED_REPEAT":return gl.MIRRORED_REPEAT;case"REPEAT":return gl.REPEAT;default:return gl.REPEAT;}};x3dom.Utils.primTypeDic=function(gl,type)
{switch(type.toUpperCase())
{case"POINTS":return gl.POINTS;case"LINES":return gl.LINES;case"LINELOOP":return gl.LINE_LOOP;case"LINESTRIP":return gl.LINE_STRIP;case"TRIANGLES":return gl.TRIANGLES;case"TRIANGLESTRIP":return gl.TRIANGLE_STRIP;case"TRIANGLEFAN":return gl.TRIANGLE_FAN;default:return gl.TRIANGLES;}};x3dom.Utils.depthFunc=function(gl,func)
{switch(func.toUpperCase())
{case"NEVER":return gl.NEVER;case"ALWAYS":return gl.ALWAYS;case"LESS":return gl.LESS;case"EQUAL":return gl.EQUAL;case"LEQUAL":return gl.LEQUAL;case"GREATER":return gl.GREATER;case"GEQUAL":return gl.GEQUAL;case"NOTEQUAL":return gl.NOTEQUAL;default:return gl.LEQUAL;}};x3dom.Utils.blendFunc=function(gl,func)
{switch(func.toLowerCase())
{case"zero":return gl.ZERO;case"one":return gl.ONE;case"dst_color":return gl.DST_COLOR;case"dst_alpha":return gl.DST_ALPHA;case"src_color":return gl.SRC_COLOR;case"src_alpha":return gl.SRC_ALPHA;case"one_minus_dst_color":return gl.ONE_MINUS_DST_COLOR;case"one_minus_dst_alpha":return gl.ONE_MINUS_DST_ALPHA;case"one_minus_src_color":return gl.ONE_MINUS_SRC_COLOR;case"one_minus_src_alpha":return gl.ONE_MINUS_SRC_ALPHA;case"src_alpha_saturate":return gl.SRC_ALPHA_SATURATE;case"constant_color":return gl.CONSTANT_COLOR;case"constant_alpha":return gl.CONSTANT_ALPHA;case"one_minus_constant_color":return gl.ONE_MINUS_CONSTANT_COLOR;case"one_minus_constant_alpha":return gl.ONE_MINUS_CONSTANT_ALPHA;default:return 0;}};x3dom.Utils.blendEquation=function(gl,func)
{switch(func.toLowerCase())
{case"func_add":return gl.FUNC_ADD;case"func_subtract":return gl.FUNC_SUBTRACT;case"func_reverse_subtract":return gl.FUNC_REVERSE_SUBTRACT;case"min":return 0;case"max":return 0;case"logic_op":return 0;default:return 0;}};x3dom.Utils.gunzip=function(arraybuffer)
{var byteArray=new Uint8Array(arraybuffer);try{arraybuffer=new Zlib.Gunzip(byteArray).decompress().buffer;}catch(e){}
return arraybuffer;};x3dom.Utils.generateProperties=function(viewarea,shape)
{var property={};var geometry=shape._cf.geometry.node;var appearance=shape._cf.appearance.node;var texture=appearance?appearance._cf.texture.node:null;var material=appearance?appearance._cf.material.node:null;var environment=viewarea._scene.getEnvironment();if(appearance&&appearance._shader&&x3dom.isa(appearance._shader,x3dom.nodeTypes.ComposedShader)){property.CSHADER=appearance._shader._id;}
else if(geometry){property.CSHADER=-1;property.SOLID=(shape.isSolid())?1:0;property.TEXT=(x3dom.isa(geometry,x3dom.nodeTypes.Text))?1:0;property.POPGEOMETRY=(x3dom.isa(geometry,x3dom.nodeTypes.PopGeometry))?1:0;property.IMAGEGEOMETRY=(x3dom.isa(geometry,x3dom.nodeTypes.ImageGeometry))?1:0;property.BINARYGEOMETRY=(x3dom.isa(geometry,x3dom.nodeTypes.BinaryGeometry))?1:0;property.EXTERNALGEOMETRY=(x3dom.isa(geometry,x3dom.nodeTypes.ExternalGeometry))?1:0;property.IG_PRECISION=(property.IMAGEGEOMETRY)?geometry.numCoordinateTextures():0;property.IG_INDEXED=(property.IMAGEGEOMETRY&&geometry.getIndexTexture()!=null)?1:0;property.POINTLINE2D=!geometry.needLighting()?1:0;property.VERTEXID=((property.BINARYGEOMETRY||property.EXTERNALGEOMETRY)&&geometry._vf.idsPerVertex)?1:0;property.IS_PARTICLE=(x3dom.isa(geometry,x3dom.nodeTypes.ParticleSet))?1:0;property.TWOSIDEDMAT=(property.APPMAT&&x3dom.isa(material,x3dom.nodeTypes.TwoSidedMaterial))?1:0;property.SEPARATEBACKMAT=(property.TWOSIDEDMAT&&material._vf.separateBackColor)?1:0;property.SHADOW=(viewarea.getLightsShadow())?1:0;property.FOG=(viewarea._scene.getFog()._vf.visibilityRange>0)?1:0;property.CSSHADER=(appearance&&appearance._shader&&x3dom.isa(appearance._shader,x3dom.nodeTypes.CommonSurfaceShader))?1:0;property.APPMAT=(appearance&&(material||property.CSSHADER))?1:0;property.LIGHTS=(!property.POINTLINE2D&&appearance&&shape.isLit()&&(material||property.CSSHADER))?viewarea.getLights().length+(viewarea._scene.getNavigationInfo()._vf.headlight):0;property.TEXTURED=(texture||property.TEXT||(property.CSSHADER&&appearance._shader.needTexcoords()))?1:0;property.CUBEMAP=(texture&&x3dom.isa(texture,x3dom.nodeTypes.X3DEnvironmentTextureNode))||(property.CSSHADER&&appearance._shader.getEnvironmentMap())?1:0;property.PIXELTEX=(texture&&x3dom.isa(texture,x3dom.nodeTypes.PixelTexture))?1:0;property.TEXTRAFO=(appearance&&appearance._cf.textureTransform.node)?1:0;property.DIFFUSEMAP=(texture&&!x3dom.isa(texture,x3dom.nodeTypes.X3DEnvironmentTextureNode))||(property.CSSHADER&&appearance._shader.getDiffuseMap())?1:0;property.NORMALMAP=(property.CSSHADER&&appearance._shader.getNormalMap())?1:0;property.NORMALSPACE=(property.NORMALMAP)?appearance._shader._vf.normalSpace.toUpperCase():"";property.SPECMAP=(property.CSSHADER&&appearance._shader.getSpecularMap())?1:0;property.SHINMAP=(property.CSSHADER&&appearance._shader.getShininessMap())?1:0;property.DISPLACEMENTMAP=(property.CSSHADER&&appearance._shader.getDisplacementMap())?1:0;property.DIFFPLACEMENTMAP=(property.CSSHADER&&appearance._shader.getDiffuseDisplacementMap())?1:0;property.MULTIDIFFALPMAP=(property.VERTEXID&&property.CSSHADER&&appearance._shader.getMultiDiffuseAlphaMap())?1:0;property.MULTIEMIAMBMAP=(property.VERTEXID&&property.CSSHADER&&appearance._shader.getMultiEmissiveAmbientMap())?1:0;property.MULTISPECSHINMAP=(property.VERTEXID&&property.CSSHADER&&appearance._shader.getMultiSpecularShininessMap())?1:0;property.MULTIVISMAP=(property.VERTEXID&&property.CSSHADER&&appearance._shader.getMultiVisibilityMap())?1:0;property.BLENDING=(property.TEXT||property.CUBEMAP||property.CSSHADER||(texture&&texture._blending))?1:0;property.REQUIREBBOX=(geometry._vf.coordType!==undefined&&geometry._vf.coordType!="Float32")?1:0;property.REQUIREBBOXNOR=(geometry._vf.normalType!==undefined&&geometry._vf.normalType!="Float32")?1:0;property.REQUIREBBOXCOL=(geometry._vf.colorType!==undefined&&geometry._vf.colorType!="Float32")?1:0;property.REQUIREBBOXTEX=(geometry._vf.texCoordType!==undefined&&geometry._vf.texCoordType!="Float32")?1:0;property.COLCOMPONENTS=geometry._mesh._numColComponents;property.NORCOMPONENTS=geometry._mesh._numNormComponents;property.POSCOMPONENTS=geometry._mesh._numPosComponents;property.SPHEREMAPPING=(geometry._cf.texCoord!==undefined&&geometry._cf.texCoord.node!==null&&geometry._cf.texCoord.node._vf.mode&&geometry._cf.texCoord.node._vf.mode.toLowerCase()=="sphere")?1:0;property.VERTEXCOLOR=(geometry._mesh._colors[0].length>0||(property.IMAGEGEOMETRY&&geometry.getColorTexture())||(property.POPGEOMETRY&&geometry.hasColor())||(geometry._vf.color!==undefined&&geometry._vf.color.length>0))?1:0;property.CLIPPLANES=shape._clipPlanes.length;property.ALPHATHRESHOLD=(appearance)?appearance._vf.alphaClipThreshold.toFixed(2):0.1;property.GAMMACORRECTION=environment._vf.gammaCorrectionDefault;property.KHR_MATERIAL_COMMONS=0;}
property.toIdentifier=function(){delete this.id;var id="";for(var p in this){if(this[p]!=this.toIdentifier&&this[p]!=this.toString){id+=this[p];}}
this.id=id;return id;};property.toString=function(){var str="";for(var p in this){if(this[p]!=this.toIdentifier&&this[p]!=this.toString){str+=p+": "+this[p]+", ";}}
return str;};property.toIdentifier();return property;};x3dom.Utils.wrapProgram=function(gl,program,shaderID)
{var shader={shaderID:shaderID,program:program};shader.bind=function(){gl.useProgram(program);};var loc=null;var obj=null;var i,glErr;var numUniforms=gl.getProgramParameter(program,gl.ACTIVE_UNIFORMS);for(i=0;i<numUniforms;++i){try{obj=gl.getActiveUniform(program,i);}
catch(eu){if(!obj)continue;}
glErr=gl.getError();if(glErr){x3dom.debug.logError("GL-Error (on searching uniforms): "+glErr);}
loc=gl.getUniformLocation(program,obj.name);switch(obj.type){case gl.SAMPLER_2D:shader.__defineSetter__(obj.name,(function(loc){return function(val){gl.uniform1i(loc,val);};})(loc));break;case gl.SAMPLER_CUBE:shader.__defineSetter__(obj.name,(function(loc){return function(val){gl.uniform1i(loc,val);};})(loc));break;case gl.BOOL:shader.__defineSetter__(obj.name,(function(loc){return function(val){gl.uniform1i(loc,val);};})(loc));break;case gl.FLOAT:if(obj.name.indexOf("[0]")!=-1)
shader.__defineSetter__(obj.name.substring(0,obj.name.length-3),(function(loc){return function(val){gl.uniform1fv(loc,new Float32Array(val));};})(loc));else
shader.__defineSetter__(obj.name,(function(loc){return function(val){gl.uniform1f(loc,val);};})(loc));break;case gl.FLOAT_VEC2:shader.__defineSetter__(obj.name,(function(loc){return function(val){gl.uniform2f(loc,val[0],val[1]);};})(loc));break;case gl.FLOAT_VEC3:if(obj.name.indexOf("[0]")!=-1)
shader.__defineSetter__(obj.name.substring(0,obj.name.length-3),(function(loc){return function(val){gl.uniform3fv(loc,new Float32Array(val));};})(loc));else
shader.__defineSetter__(obj.name,(function(loc){return function(val){gl.uniform3f(loc,val[0],val[1],val[2]);};})(loc));break;case gl.FLOAT_VEC4:shader.__defineSetter__(obj.name,(function(loc){return function(val){gl.uniform4f(loc,val[0],val[1],val[2],val[3]);};})(loc));break;case gl.FLOAT_MAT2:shader.__defineSetter__(obj.name,(function(loc){return function(val){gl.uniformMatrix2fv(loc,false,new Float32Array(val));};})(loc));break;case gl.FLOAT_MAT3:shader.__defineSetter__(obj.name,(function(loc){return function(val){gl.uniformMatrix3fv(loc,false,new Float32Array(val));};})(loc));break;case gl.FLOAT_MAT4:shader.__defineSetter__(obj.name,(function(loc){return function(val){gl.uniformMatrix4fv(loc,false,new Float32Array(val));};})(loc));break;case gl.INT:shader.__defineSetter__(obj.name,(function(loc){return function(val){gl.uniform1i(loc,val);};})(loc));break;default:x3dom.debug.logWarning('GLSL program variable '+obj.name+' has unknown type '+obj.type);}}
var numAttribs=gl.getProgramParameter(program,gl.ACTIVE_ATTRIBUTES);for(i=0;i<numAttribs;++i){try{obj=gl.getActiveAttrib(program,i);}
catch(ea){if(!obj)continue;}
glErr=gl.getError();if(glErr){x3dom.debug.logError("GL-Error (on searching attributes): "+glErr);}
loc=gl.getAttribLocation(program,obj.name);shader[obj.name]=loc;}
return shader;};x3dom.Utils.forbiddenBySOP=function(uri_string){uri_string=uri_string.toLowerCase();var Scheme_AuthorityPQF=uri_string.split('//');var Scheme;var AuthorityPQF;var Authority;var UserInfo_HostPort;var HostPort;var Host_Port;var Port;var Host;var originPort=document.location.port===""?"80":document.location.port;if(Scheme_AuthorityPQF.length===2){Scheme=Scheme_AuthorityPQF[0];AuthorityPQF=Scheme_AuthorityPQF[1];Authority=AuthorityPQF.split('/')[0].split('?')[0].split('#')[0];UserInfo_HostPort=Authority.split('@');if(UserInfo_HostPort.length===1){HostPort=UserInfo_HostPort[0];}else{HostPort=UserInfo_HostPort[1];}
Host_Port=HostPort.split(':');Host=Host_Port[0];Port=Host_Port[1];}
Port=Port||"80";Host=Host||document.location.host;Scheme=Scheme||document.location.protocol;return!(Port===originPort&&Host===document.location.host&&Scheme===document.location.protocol);};x3dom.States=function(x3dElem){var that=this;this.active=false;this.viewer=document.createElement('div');this.viewer.id='x3dom-state-viewer';var title=document.createElement('div');title.className='x3dom-states-head';title.appendChild(document.createTextNode('x3dom'));var subTitle=document.createElement('span');subTitle.className='x3dom-states-head2';subTitle.appendChild(document.createTextNode('stats'));title.appendChild(subTitle);this.renderMode=document.createElement('div');this.renderMode.className='x3dom-states-rendermode-hardware';this.measureList=document.createElement('ul');this.measureList.className='x3dom-states-list';this.infoList=document.createElement('ul');this.infoList.className='x3dom-states-list';this.requestList=document.createElement('ul');this.requestList.className='x3dom-states-list';this.viewer.appendChild(this.renderMode);this.viewer.appendChild(this.measureList);this.viewer.appendChild(this.infoList);this.viewer.appendChild(this.requestList);this.disableContextMenu=function(e){e.preventDefault();e.stopPropagation();e.returnValue=false;return false;};this.thousandSeperator=function(value){return value.toString().replace(/\B(?=(\d{3})+(?!\d))/g,",");};this.toFixed=function(value){var fixed=(value<1)?2:(value<10)?2:2;return value.toFixed(fixed);};this.addItem=function(list,key,value){var item=document.createElement('li');item.className='x3dom-states-item';var keyDiv=document.createElement('div');keyDiv.className='x3dom-states-item-title';keyDiv.appendChild(document.createTextNode(key));var valueDiv=document.createElement('div');valueDiv.className='x3dom-states-item-value';valueDiv.appendChild(document.createTextNode(value));item.appendChild(keyDiv);item.appendChild(valueDiv);list.appendChild(item);};this.update=function(){if(!x3dElem.runtime&&this.updateMethodID!==undefined){clearInterval(this.updateMethodID);return;}
var infos=x3dElem.runtime.states.infos;var measurements=x3dElem.runtime.states.measurements;var renderMode=x3dom.caps.RENDERMODE;if(renderMode=="HARDWARE"){this.renderMode.innerHTML="Hardware-Rendering";this.renderMode.className='x3dom-states-rendermode-hardware';}else if(renderMode=="SOFTWARE"){this.renderMode.innerHTML="Software-Rendering";this.renderMode.className='x3dom-states-rendermode-software';}
this.measureList.innerHTML="";for(var m in measurements)
{if(measurements.hasOwnProperty(m))
{this.addItem(this.measureList,m,this.toFixed(measurements[m]));}}
this.infoList.innerHTML="";for(var i in infos)
{if(infos.hasOwnProperty(i))
{this.addItem(this.infoList,i,this.thousandSeperator(infos[i]));}}
this.requestList.innerHTML="";this.addItem(this.requestList,"#ACTIVE",x3dom.RequestManager.activeRequests.length);this.addItem(this.requestList,"#TOTAL",x3dom.RequestManager.totalRequests);this.addItem(this.requestList,"#LOADED",x3dom.RequestManager.loadedRequests);this.addItem(this.requestList,"#FAILED",x3dom.RequestManager.failedRequests);};this.updateMethodID=window.setInterval(function(){that.update();},1000);this.viewer.addEventListener("contextmenu",that.disableContextMenu);};x3dom.States.prototype.display=function(value){this.active=(value!==undefined)?value:!this.active;this.viewer.style.display=(this.active)?"block":"none";};x3dom.StateManager=function(ctx3d)
{this.gl=ctx3d;this.states=[];this.initStates();};x3dom.StateManager.prototype.initStates=function()
{this.states['shaderID']=null;this.states['colorMask']={red:null,green:null,blue:null,alpha:null};this.states['depthMask']=null;this.states['stencilMask']=null;this.states['cullFace']=null;this.states['frontFace']=null;this.states['lineWidth']=null;this.states['blendColor']={red:null,green:null,blue:null,alpha:null};this.states['blendEquation']=null;this.states['blendEquationSeparate']={modeRGB:null,modeAlpha:null};this.states['blendFunc']={sfactor:null,dfactor:null};this.states['blendFuncSeparate']={srcRGB:null,dstRGB:null,srcAlpha:null,dstAlpha:null};this.states['depthFunc']=null;this.states['viewport']={x:null,y:null,width:null,height:null};this.states['depthRange']={zNear:null,zFar:null};};x3dom.StateManager.prototype.useProgram=function(shader)
{if(this.states['shaderID']!=shader.shaderID)
{this.gl.useProgram(shader.program);this.states['shaderID']=shader.shaderID;return true;}
return false;};x3dom.StateManager.prototype.unsetProgram=function()
{this.states['shaderID']=null;};x3dom.StateManager.prototype.enable=function(cap)
{if(this.states[cap]!==true)
{this.gl.enable(cap);this.states[cap]=true;}};x3dom.StateManager.prototype.disable=function(cap)
{if(this.states[cap]!==false)
{this.gl.disable(cap);this.states[cap]=false;}};x3dom.StateManager.prototype.colorMask=function(red,green,blue,alpha)
{if(this.states['colorMask'].red!=red||this.states['colorMask'].green!=green||this.states['colorMask'].blue!=blue||this.states['colorMask'].alpha!=alpha)
{this.gl.colorMask(red,green,blue,alpha);this.states['colorMask'].red=red;this.states['colorMask'].green=green;this.states['colorMask'].blue=blue;this.states['colorMask'].alpha=alpha;}};x3dom.StateManager.prototype.depthMask=function(flag)
{if(this.states['depthMask']!=flag)
{this.gl.depthMask(flag);this.states['depthMask']=flag;}};x3dom.StateManager.prototype.stencilMask=function(mask)
{if(this.states['stencilMask']!=mask)
{this.gl.stencilMask(mask);this.states['stencilMask']=mask;}};x3dom.StateManager.prototype.cullFace=function(mode)
{if(this.states['cullFace']!=mode)
{this.gl.cullFace(mode);this.states['cullFace']=mode;}};x3dom.StateManager.prototype.frontFace=function(mode)
{if(this.states['frontFace']!=mode)
{this.gl.frontFace(mode);this.states['frontFace']=mode;}};x3dom.StateManager.prototype.lineWidth=function(width)
{width=(width<=1)?1:width;if(this.states['lineWidth']!=width)
{this.gl.lineWidth(width);this.states['lineWidth']=width;}};x3dom.StateManager.prototype.blendColor=function(red,green,blue,alpha)
{if(this.states['blendColor'].red!=red||this.states['blendColor'].green!=green||this.states['blendColor'].blue!=blue||this.states['blendColor'].alpha!=alpha)
{this.gl.blendColor(red,green,blue,alpha);this.states['blendColor'].red=red;this.states['blendColor'].green=green;this.states['blendColor'].blue=blue;this.states['blendColor'].alpha=alpha;}};x3dom.StateManager.prototype.blendEquation=function(mode)
{if(mode&&this.states['blendEquation']!=mode)
{this.gl.blendEquation(mode);this.states['blendEquation']=mode;}};x3dom.StateManager.prototype.blendEquationSeparate=function(modeRGB,modeAlpha)
{if(this.states['blendEquationSeparate'].modeRGB!=modeRGB||this.states['blendEquationSeparate'].modeAlpha!=modeAlpha)
{this.gl.blendEquationSeparate(modeRGB,modeAlpha);this.states['blendEquationSeparate'].modeRGB=modeRGB;this.states['blendEquationSeparate'].modeAlpha=modeAlpha;}};x3dom.StateManager.prototype.blendFunc=function(sfactor,dfactor)
{if(this.states['blendFunc'].sfactor!=sfactor||this.states['blendFunc'].dfactor!=dfactor)
{this.gl.blendFunc(sfactor,dfactor);this.states['blendFunc'].sfactor=sfactor;this.states['blendFunc'].dfactor=dfactor;}};x3dom.StateManager.prototype.blendFuncSeparate=function(srcRGB,dstRGB,srcAlpha,dstAlpha)
{if(this.states['blendFuncSeparate'].srcRGB!=srcRGB||this.states['blendFuncSeparate'].dstRGB!=dstRGB||this.states['blendFuncSeparate'].srcAlpha!=srcAlpha||this.states['blendFuncSeparate'].dstAlpha!=dstAlpha)
{this.gl.blendFuncSeparate(srcRGB,dstRGB,srcAlpha,dstAlpha);this.states['blendFuncSeparate'].srcRGB=srcRGB;this.states['blendFuncSeparate'].dstRGB=dstRGB;this.states['blendFuncSeparate'].srcAlpha=srcAlpha;this.states['blendFuncSeparate'].dstAlpha=dstAlpha;}};x3dom.StateManager.prototype.depthFunc=function(func)
{if(this.states['depthFunc']!=func)
{this.gl.depthFunc(func);this.states['depthFunc']=func;}};x3dom.StateManager.prototype.depthRange=function(zNear,zFar)
{if(zNear<0||zFar<0||zNear>zFar)
{return;}
zNear=(zNear>1)?1:zNear;zFar=(zFar>1)?1:zFar;if(this.states['depthRange'].zNear!=zNear||this.states['depthRange'].zFar!=zFar)
{this.gl.depthRange(zNear,zFar);this.states['depthRange'].zNear=zNear;this.states['depthRange'].zFar=zFar;}};x3dom.StateManager.prototype.viewport=function(x,y,width,height)
{if(this.states['viewport'].x!=x||this.states['viewport'].y!=y||this.states['viewport'].width!=width||this.states['viewport'].height!=height)
{this.gl.viewport(x,y,width,height);this.states['viewport'].x=x;this.states['viewport'].y=y;this.states['viewport'].width=width;this.states['viewport'].height=height;}};x3dom.StateManager.prototype.bindFramebuffer=function(target,framebuffer)
{this.gl.bindFramebuffer(target,framebuffer);this.initStates();};x3dom.BinaryContainerLoader={outOfMemory:false,checkError:function(gl){var glErr=gl.getError();if(glErr){if(glErr==gl.OUT_OF_MEMORY){this.outOfMemory=true;x3dom.debug.logError("GL-Error "+glErr+" on loading binary container (out of memory).");console.error("WebGL: OUT_OF_MEMORY");}
else{x3dom.debug.logError("GL-Error "+glErr+" on loading binary container.");}}}};x3dom.BinaryContainerLoader.setupBinGeo=function(shape,sp,gl,viewarea,currContext)
{if(this.outOfMemory){return;}
var t00=new Date().getTime();var that=this;var binGeo=shape._cf.geometry.node;shape._webgl.binaryGeometry=-1;shape._webgl.internalDownloadCount=((binGeo._vf.index.length>0)?1:0)+
((binGeo._hasStrideOffset&&binGeo._vf.coord.length>0)?1:0)+
((!binGeo._hasStrideOffset&&binGeo._vf.coord.length>0)?1:0)+
((!binGeo._hasStrideOffset&&binGeo._vf.normal.length>0)?1:0)+
((!binGeo._hasStrideOffset&&binGeo._vf.texCoord.length>0)?1:0)+
((!binGeo._hasStrideOffset&&binGeo._vf.color.length>0)?1:0);var createTriangleSoup=(binGeo._vf.normalPerVertex==false)||((binGeo._vf.index.length>0)&&(binGeo._vf.indexType=="Int32"||(binGeo._vf.indexType=="Uint32"&&!x3dom.caps.INDEX_UINT)));shape._webgl.makeSeparateTris={index:null,coord:null,normal:null,texCoord:null,color:null,pushBuffer:function(name,buf){this[name]=buf;if(--shape._webgl.internalDownloadCount==0){if(this.coord)
this.createMesh();shape._nameSpace.doc.needRender=true;}
if(--shape._nameSpace.doc.downloadCount==0)
shape._nameSpace.doc.needRender=true;},createMesh:function(){var geoNode=binGeo;if(geoNode._hasStrideOffset){x3dom.debug.logError(geoNode._vf.indexType+" index type and per-face normals not supported for interleaved arrays.");return;}
for(var k=0;k<shape._webgl.primType.length;k++){if(shape._webgl.primType[k]==gl.TRIANGLE_STRIP){x3dom.debug.logError("makeSeparateTris: triangle strips not yet supported for per-face normals.");return;}}
var attribTypeStr=geoNode._vf.coordType;shape._webgl.coordType=x3dom.Utils.getVertexAttribType(attribTypeStr,gl);var bgCenter,bgSize,bgPrecisionMax;if(shape._webgl.coordType!=gl.FLOAT)
{if(geoNode._mesh._numPosComponents==4&&x3dom.Utils.isUnsignedType(geoNode._vf.coordType))
bgCenter=x3dom.fields.SFVec3f.copy(geoNode.getMin());else
bgCenter=x3dom.fields.SFVec3f.copy(geoNode._vf.position);bgSize=x3dom.fields.SFVec3f.copy(geoNode._vf.size);bgPrecisionMax=geoNode.getPrecisionMax('coordType');}
else
{bgCenter=new x3dom.fields.SFVec3f(0,0,0);bgSize=new x3dom.fields.SFVec3f(1,1,1);bgPrecisionMax=1.0;}
var dataLen=shape._coordStrideOffset[0]/x3dom.Utils.getDataTypeSize(geoNode._vf.coordType);dataLen=(dataLen==0)?3:dataLen;x3dom.debug.logWarning("makeSeparateTris.createMesh called with coord length "+dataLen);if(this.color&&dataLen!=shape._colorStrideOffset[0]/x3dom.Utils.getDataTypeSize(geoNode._vf.colorType))
{this.color=null;x3dom.debug.logWarning("Color format not supported.");}
var texDataLen=this.texCoord?(shape._texCoordStrideOffset[0]/x3dom.Utils.getDataTypeSize(geoNode._vf.texCoordType)):0;geoNode._vf.normalType="Float32";shape._webgl.normalType=gl.FLOAT;geoNode._mesh._numNormComponents=3;shape._normalStrideOffset=[0,0];var posBuf=[],normBuf=[],texcBuf=[],colBuf=[];var i,j,l,n=this.index?(this.index.length-2):(this.coord.length/3-2);for(i=0;i<n;i+=3)
{j=dataLen*(this.index?this.index[i]:i);var p0=new x3dom.fields.SFVec3f(bgSize.x*this.coord[j]/bgPrecisionMax,bgSize.y*this.coord[j+1]/bgPrecisionMax,bgSize.z*this.coord[j+2]/bgPrecisionMax);posBuf.push(this.coord[j]);posBuf.push(this.coord[j+1]);posBuf.push(this.coord[j+2]);if(dataLen>3)posBuf.push(this.coord[j+3]);if(this.color){colBuf.push(this.color[j]);colBuf.push(this.color[j+1]);colBuf.push(this.color[j+2]);if(dataLen>3)colBuf.push(this.color[j+3]);}
if(this.texCoord){l=texDataLen*(this.index?this.index[i]:i);texcBuf.push(this.texCoord[l]);texcBuf.push(this.texCoord[l+1]);if(texDataLen>3){texcBuf.push(this.texCoord[l+2]);texcBuf.push(this.texCoord[l+3]);}}
j=dataLen*(this.index?this.index[i+1]:i+1);var p1=new x3dom.fields.SFVec3f(bgSize.x*this.coord[j]/bgPrecisionMax,bgSize.y*this.coord[j+1]/bgPrecisionMax,bgSize.z*this.coord[j+2]/bgPrecisionMax);posBuf.push(this.coord[j]);posBuf.push(this.coord[j+1]);posBuf.push(this.coord[j+2]);if(dataLen>3)posBuf.push(this.coord[j+3]);if(this.color){colBuf.push(this.color[j]);colBuf.push(this.color[j+1]);colBuf.push(this.color[j+2]);if(dataLen>3)colBuf.push(this.color[j+3]);}
if(this.texCoord){l=texDataLen*(this.index?this.index[i+1]:i+1);texcBuf.push(this.texCoord[l]);texcBuf.push(this.texCoord[l+1]);if(texDataLen>3){texcBuf.push(this.texCoord[l+2]);texcBuf.push(this.texCoord[l+3]);}}
j=dataLen*(this.index?this.index[i+2]:i+2);var p2=new x3dom.fields.SFVec3f(bgSize.x*this.coord[j]/bgPrecisionMax,bgSize.y*this.coord[j+1]/bgPrecisionMax,bgSize.z*this.coord[j+2]/bgPrecisionMax);posBuf.push(this.coord[j]);posBuf.push(this.coord[j+1]);posBuf.push(this.coord[j+2]);if(dataLen>3)posBuf.push(this.coord[j+3]);if(this.color){colBuf.push(this.color[j]);colBuf.push(this.color[j+1]);colBuf.push(this.color[j+2]);if(dataLen>3)colBuf.push(this.color[j+3]);}
if(this.texCoord){l=texDataLen*(this.index?this.index[i+2]:i+2);texcBuf.push(this.texCoord[l]);texcBuf.push(this.texCoord[l+1]);if(texDataLen>3){texcBuf.push(this.texCoord[l+2]);texcBuf.push(this.texCoord[l+3]);}}
var a=p0.subtract(p1);var b=p1.subtract(p2);var norm=a.cross(b).normalize();for(j=0;j<3;j++){normBuf.push(norm.x);normBuf.push(norm.y);normBuf.push(norm.z);}}
var buffer=gl.createBuffer();shape._webgl.buffers[1]=buffer;gl.bindBuffer(gl.ARRAY_BUFFER,buffer);gl.bufferData(gl.ARRAY_BUFFER,x3dom.Utils.getArrayBufferView(geoNode._vf.coordType,posBuf),gl.STATIC_DRAW);gl.vertexAttribPointer(sp.position,geoNode._mesh._numPosComponents,shape._webgl.coordType,false,shape._coordStrideOffset[0],shape._coordStrideOffset[1]);gl.enableVertexAttribArray(sp.position);buffer=gl.createBuffer();shape._webgl.buffers[2]=buffer;gl.bindBuffer(gl.ARRAY_BUFFER,buffer);gl.bufferData(gl.ARRAY_BUFFER,new Float32Array(normBuf),gl.STATIC_DRAW);gl.vertexAttribPointer(sp.normal,geoNode._mesh._numNormComponents,shape._webgl.normalType,false,shape._normalStrideOffset[0],shape._normalStrideOffset[1]);gl.enableVertexAttribArray(sp.normal);if(this.texCoord)
{buffer=gl.createBuffer();shape._webgl.buffers[3]=buffer;gl.bindBuffer(gl.ARRAY_BUFFER,buffer);gl.bufferData(gl.ARRAY_BUFFER,x3dom.Utils.getArrayBufferView(geoNode._vf.texCoordType,texcBuf),gl.STATIC_DRAW);gl.vertexAttribPointer(sp.texcoord,geoNode._mesh._numTexComponents,shape._webgl.texCoordType,false,shape._texCoordStrideOffset[0],shape._texCoordStrideOffset[1]);gl.enableVertexAttribArray(sp.texcoord);}
if(this.color)
{buffer=gl.createBuffer();shape._webgl.buffers[4]=buffer;gl.bindBuffer(gl.ARRAY_BUFFER,buffer);gl.bufferData(gl.ARRAY_BUFFER,x3dom.Utils.getArrayBufferView(geoNode._vf.colorType,colBuf),gl.STATIC_DRAW);gl.vertexAttribPointer(sp.color,geoNode._mesh._numColComponents,shape._webgl.colorType,false,shape._colorStrideOffset[0],shape._colorStrideOffset[1]);gl.enableVertexAttribArray(sp.color);}
geoNode._vf.vertexCount=[];geoNode._vf.vertexCount[0]=posBuf.length/dataLen;geoNode._mesh._numCoords=geoNode._vf.vertexCount[0];geoNode._mesh._numFaces=geoNode._vf.vertexCount[0]/3;shape._webgl.primType=[];shape._webgl.primType[0]=gl.TRIANGLES;posBuf=null;normBuf=null;texcBuf=null;colBuf=null;this.index=null;this.coord=null;this.normal=null;this.texCoord=null;this.color=null;that.checkError(gl);delete shape._webgl.shader;shape._webgl.shader=currContext.cache.getDynamicShader(gl,viewarea,shape);}};if(binGeo._vf.index.length>0)
{shape._webgl.binaryGeometry=1;var xmlhttp0=new XMLHttpRequest();xmlhttp0.open("GET",shape._nameSpace.getURL(binGeo._vf.index),true);xmlhttp0.responseType="arraybuffer";shape._nameSpace.doc.downloadCount+=1;x3dom.RequestManager.addRequest(xmlhttp0);xmlhttp0.onload=function()
{shape._nameSpace.doc.downloadCount-=1;shape._webgl.internalDownloadCount-=1;if(xmlhttp0.status!=200){x3dom.debug.logError("XHR1/ index load failed with status: "+xmlhttp0.status);return;}
if(!shape._webgl)
return;var XHR_buffer=binGeo._vf.compressed==true?x3dom.Utils.gunzip(xmlhttp0.response):xmlhttp0.response;var geoNode=binGeo;var attribTypeStr=geoNode._vf.indexType;var indexArray=x3dom.Utils.getArrayBufferView(attribTypeStr,XHR_buffer);if(createTriangleSoup){shape._webgl.makeSeparateTris.pushBuffer("index",indexArray);return;}
var indicesBuffer=gl.createBuffer();if(x3dom.caps.INDEX_UINT&&attribTypeStr=="Uint32"){shape._webgl.indexType=gl.UNSIGNED_INT;}
else{shape._webgl.indexType=gl.UNSIGNED_SHORT;}
gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER,indicesBuffer);gl.bufferData(gl.ELEMENT_ARRAY_BUFFER,indexArray,gl.STATIC_DRAW);gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER,null);if(geoNode._vf.vertexCount[0]==0)
geoNode._vf.vertexCount[0]=indexArray.length;geoNode._mesh._numFaces=0;for(var i=0;i<geoNode._vf.vertexCount.length;i++){if(shape._webgl.primType[i]==gl.TRIANGLE_STRIP)
geoNode._mesh._numFaces+=geoNode._vf.vertexCount[i]-2;else
geoNode._mesh._numFaces+=geoNode._vf.vertexCount[i]/3;}
indexArray=null;if(shape._webgl.internalDownloadCount==0)
{shape._nameSpace.doc.needRender=true;}
that.checkError(gl);var t11=new Date().getTime()-t00;x3dom.debug.logInfo("XHR0/ index load time: "+t11+" ms");shape._webgl.buffers[0]=indicesBuffer;};}
if(binGeo._hasStrideOffset&&binGeo._vf.coord.length>0)
{var xmlhttp=new XMLHttpRequest();xmlhttp.open("GET",shape._nameSpace.getURL(binGeo._vf.coord),true);xmlhttp.responseType="arraybuffer";shape._nameSpace.doc.downloadCount+=1;x3dom.RequestManager.addRequest(xmlhttp);xmlhttp.onload=function()
{shape._nameSpace.doc.downloadCount-=1;shape._webgl.internalDownloadCount-=1;if(xmlhttp.status!=200){x3dom.debug.logError("XHR1/ interleaved array load failed with status: "+xmlhttp.status);return;}
if(!shape._webgl)
return;var XHR_buffer=binGeo._vf.compressed==true?x3dom.Utils.gunzip(xmlhttp.response):xmlhttp.response;var geoNode=binGeo;var attribTypeStr=geoNode._vf.coordType;shape._webgl.coordType=x3dom.Utils.getVertexAttribType(attribTypeStr,gl);shape._webgl.normalType=shape._webgl.coordType;shape._webgl.texCoordType=shape._webgl.coordType;shape._webgl.colorType=shape._webgl.coordType;var attributes=x3dom.Utils.getArrayBufferView(attribTypeStr,XHR_buffer);var dataLen=shape._coordStrideOffset[0]/x3dom.Utils.getDataTypeSize(attribTypeStr);if(dataLen)
geoNode._mesh._numCoords=attributes.length/dataLen;if(geoNode._vf.index.length==0){for(var i=0;i<geoNode._vf.vertexCount.length;i++){if(shape._webgl.primType[i]==gl.TRIANGLE_STRIP)
geoNode._mesh._numFaces+=geoNode._vf.vertexCount[i]-2;else
geoNode._mesh._numFaces+=geoNode._vf.vertexCount[i]/3;}}
var buffer=gl.createBuffer();gl.bindBuffer(gl.ARRAY_BUFFER,buffer);gl.bufferData(gl.ARRAY_BUFFER,attributes,gl.STATIC_DRAW);gl.vertexAttribPointer(sp.position,geoNode._mesh._numPosComponents,shape._webgl.coordType,false,shape._coordStrideOffset[0],shape._coordStrideOffset[1]);gl.enableVertexAttribArray(sp.position);if(geoNode._vf.normal.length>0)
{shape._webgl.buffers[2]=buffer;gl.bindBuffer(gl.ARRAY_BUFFER,buffer);gl.bufferData(gl.ARRAY_BUFFER,attributes,gl.STATIC_DRAW);gl.vertexAttribPointer(sp.normal,geoNode._mesh._numNormComponents,shape._webgl.normalType,false,shape._normalStrideOffset[0],shape._normalStrideOffset[1]);gl.enableVertexAttribArray(sp.normal);}
if(geoNode._vf.texCoord.length>0)
{shape._webgl.buffers[3]=buffer;gl.bindBuffer(gl.ARRAY_BUFFER,buffer);gl.bufferData(gl.ARRAY_BUFFER,attributes,gl.STATIC_DRAW);gl.vertexAttribPointer(sp.texcoord,geoNode._mesh._numTexComponents,shape._webgl.texCoordType,false,shape._texCoordStrideOffset[0],shape._texCoordStrideOffset[1]);gl.enableVertexAttribArray(sp.texcoord);}
if(geoNode._vf.color.length>0)
{shape._webgl.buffers[4]=buffer;gl.bindBuffer(gl.ARRAY_BUFFER,buffer);gl.bufferData(gl.ARRAY_BUFFER,attributes,gl.STATIC_DRAW);gl.vertexAttribPointer(sp.color,geoNode._mesh._numColComponents,shape._webgl.colorType,false,shape._colorStrideOffset[0],shape._colorStrideOffset[1]);gl.enableVertexAttribArray(sp.color);}
attributes=null;if(shape._webgl.internalDownloadCount==0)
{shape._nameSpace.doc.needRender=true;}
that.checkError(gl);var t11=new Date().getTime()-t00;x3dom.debug.logInfo("XHR/ interleaved array load time: "+t11+" ms");shape._webgl.buffers[1]=buffer;};}
if(!binGeo._hasStrideOffset&&binGeo._vf.coord.length>0)
{var xmlhttp1=new XMLHttpRequest();xmlhttp1.open("GET",shape._nameSpace.getURL(binGeo._vf.coord),true);xmlhttp1.responseType="arraybuffer";shape._nameSpace.doc.downloadCount+=1;x3dom.RequestManager.addRequest(xmlhttp1);xmlhttp1.onload=function()
{shape._nameSpace.doc.downloadCount-=1;shape._webgl.internalDownloadCount-=1;if(xmlhttp1.status!=200){x3dom.debug.logError("XHR1/ coord load failed with status: "+xmlhttp1.status);return;}
if(!shape._webgl)
return;var XHR_buffer=binGeo._vf.compressed==true?x3dom.Utils.gunzip(xmlhttp1.response):xmlhttp1.response;var geoNode=binGeo;var i=0;var attribTypeStr=geoNode._vf.coordType;shape._webgl.coordType=x3dom.Utils.getVertexAttribType(attribTypeStr,gl);var vertices=x3dom.Utils.getArrayBufferView(attribTypeStr,XHR_buffer);if(createTriangleSoup){shape._webgl.makeSeparateTris.pushBuffer("coord",vertices);return;}
gl.bindAttribLocation(sp.program,0,"position");var positionBuffer=gl.createBuffer();gl.bindBuffer(gl.ARRAY_BUFFER,positionBuffer);gl.bufferData(gl.ARRAY_BUFFER,vertices,gl.STATIC_DRAW);gl.bindBuffer(gl.ARRAY_BUFFER,null);geoNode._mesh._numCoords=vertices.length/geoNode._mesh._numPosComponents;if(geoNode._vf.index.length==0){for(i=0;i<geoNode._vf.vertexCount.length;i++){if(shape._webgl.primType[i]==gl.TRIANGLE_STRIP)
geoNode._mesh._numFaces+=geoNode._vf.vertexCount[i]-2;else
geoNode._mesh._numFaces+=geoNode._vf.vertexCount[i]/3;}}
if((attribTypeStr=="Float32")&&(shape._vf.bboxSize.x<0||shape._vf.bboxSize.y<0||shape._vf.bboxSize.z<0))
{var min=new x3dom.fields.SFVec3f(vertices[0],vertices[1],vertices[2]);var max=new x3dom.fields.SFVec3f(vertices[0],vertices[1],vertices[2]);for(i=3;i<vertices.length;i+=3)
{if(min.x>vertices[i+0]){min.x=vertices[i+0];}
if(min.y>vertices[i+1]){min.y=vertices[i+1];}
if(min.z>vertices[i+2]){min.z=vertices[i+2];}
if(max.x<vertices[i+0]){max.x=vertices[i+0];}
if(max.y<vertices[i+1]){max.y=vertices[i+1];}
if(max.z<vertices[i+2]){max.z=vertices[i+2];}}
shape._vf.bboxCenter.setValues(min.add(max).multiply(0.5));shape._vf.bboxSize.setValues(max.subtract(min));}
vertices=null;if(shape._webgl.internalDownloadCount==0)
{shape._nameSpace.doc.needRender=true;}
that.checkError(gl);var t11=new Date().getTime()-t00;x3dom.debug.logInfo("XHR1/ coord load time: "+t11+" ms");shape._webgl.buffers[1]=positionBuffer;};}
if(!binGeo._hasStrideOffset&&binGeo._vf.normal.length>0)
{var xmlhttp2=new XMLHttpRequest();xmlhttp2.open("GET",shape._nameSpace.getURL(binGeo._vf.normal),true);xmlhttp2.responseType="arraybuffer";shape._nameSpace.doc.downloadCount+=1;x3dom.RequestManager.addRequest(xmlhttp2);xmlhttp2.onload=function()
{shape._nameSpace.doc.downloadCount-=1;shape._webgl.internalDownloadCount-=1;if(xmlhttp2.status!=200){x3dom.debug.logError("XHR2/ normal load failed with status: "+xmlhttp2.status);return;}
if(!shape._webgl)
return;var XHR_buffer=binGeo._vf.compressed==true?x3dom.Utils.gunzip(xmlhttp2.response):xmlhttp2.response;var attribTypeStr=binGeo._vf.normalType;shape._webgl.normalType=x3dom.Utils.getVertexAttribType(attribTypeStr,gl);var normals=x3dom.Utils.getArrayBufferView(attribTypeStr,XHR_buffer);if(createTriangleSoup){shape._webgl.makeSeparateTris.pushBuffer("normal",normals);return;}
var normalBuffer=gl.createBuffer();gl.bindBuffer(gl.ARRAY_BUFFER,normalBuffer);gl.bufferData(gl.ARRAY_BUFFER,normals,gl.STATIC_DRAW);gl.bindBuffer(gl.ARRAY_BUFFER,null);normals=null;if(shape._webgl.internalDownloadCount==0)
{shape._nameSpace.doc.needRender=true;}
that.checkError(gl);var t11=new Date().getTime()-t00;x3dom.debug.logInfo("XHR2/ normal load time: "+t11+" ms");shape._webgl.buffers[2]=normalBuffer;};}
if(!binGeo._hasStrideOffset&&binGeo._vf.texCoord.length>0)
{var xmlhttp3=new XMLHttpRequest();xmlhttp3.open("GET",shape._nameSpace.getURL(binGeo._vf.texCoord),true);xmlhttp3.responseType="arraybuffer";shape._nameSpace.doc.downloadCount+=1;x3dom.RequestManager.addRequest(xmlhttp3);xmlhttp3.onload=function()
{var i,j;var tmp;shape._nameSpace.doc.downloadCount-=1;shape._webgl.internalDownloadCount-=1;if(xmlhttp3.status!=200){x3dom.debug.logError("XHR3/ texcoord load failed with status: "+xmlhttp3.status);return;}
if(!shape._webgl)
return;var XHR_buffer=binGeo._vf.compressed==true?x3dom.Utils.gunzip(xmlhttp3.response):xmlhttp3.response;var attribTypeStr=binGeo._vf.texCoordType;shape._webgl.texCoordType=x3dom.Utils.getVertexAttribType(attribTypeStr,gl);var texCoords=x3dom.Utils.getArrayBufferView(attribTypeStr,XHR_buffer);if(createTriangleSoup){shape._webgl.makeSeparateTris.pushBuffer("texCoord",texCoords);return;}
if(binGeo._vf["idsPerVertex"])
{var idBuffer=gl.createBuffer();gl.bindBuffer(gl.ARRAY_BUFFER,idBuffer);var ids=x3dom.Utils.getArrayBufferView("Float32",texCoords.length/2);for(i=0,j=0;i<texCoords.length;i+=2,j++)
{ids[j]=texCoords[i+1]*65536+texCoords[i];}
gl.bufferData(gl.ARRAY_BUFFER,ids,gl.STATIC_DRAW);gl.bindBuffer(gl.ARRAY_BUFFER,null);shape._webgl.buffers[5]=idBuffer;}
else
{var texcBuffer=gl.createBuffer();gl.bindBuffer(gl.ARRAY_BUFFER,texcBuffer);gl.bufferData(gl.ARRAY_BUFFER,texCoords,gl.STATIC_DRAW);gl.bindBuffer(gl.ARRAY_BUFFER,null);shape._webgl.buffers[3]=texcBuffer;}
texCoords=null;if(shape._webgl.internalDownloadCount==0)
{shape._nameSpace.doc.needRender=true;}
that.checkError(gl);var t11=new Date().getTime()-t00;x3dom.debug.logInfo("XHR3/ texCoord load time: "+t11+" ms");};}
if(!binGeo._hasStrideOffset&&binGeo._vf.color.length>0)
{var xmlhttp4=new XMLHttpRequest();xmlhttp4.open("GET",shape._nameSpace.getURL(binGeo._vf.color),true);xmlhttp4.responseType="arraybuffer";shape._nameSpace.doc.downloadCount+=1;x3dom.RequestManager.addRequest(xmlhttp4);xmlhttp4.onload=function()
{shape._nameSpace.doc.downloadCount-=1;shape._webgl.internalDownloadCount-=1;if(xmlhttp4.status!=200){x3dom.debug.logError("XHR4/ color load failed with status: "+xmlhttp4.status);return;}
if(!shape._webgl)
return;var XHR_buffer=binGeo._vf.compressed==true?x3dom.Utils.gunzip(xmlhttp4.response):xmlhttp4.response;var attribTypeStr=binGeo._vf.colorType;shape._webgl.colorType=x3dom.Utils.getVertexAttribType(attribTypeStr,gl);var colors=x3dom.Utils.getArrayBufferView(attribTypeStr,XHR_buffer);if(createTriangleSoup){shape._webgl.makeSeparateTris.pushBuffer("color",colors);return;}
var colorBuffer=gl.createBuffer();gl.bindBuffer(gl.ARRAY_BUFFER,colorBuffer);gl.bufferData(gl.ARRAY_BUFFER,colors,gl.STATIC_DRAW);gl.bindBuffer(gl.ARRAY_BUFFER,null);colors=null;if(shape._webgl.internalDownloadCount==0)
{shape._nameSpace.doc.needRender=true;}
that.checkError(gl);var t11=new Date().getTime()-t00;x3dom.debug.logInfo("XHR4/ color load time: "+t11+" ms");shape._webgl.buffers[4]=colorBuffer;};}};x3dom.BinaryContainerLoader.setupPopGeo=function(shape,sp,gl,viewarea,currContext)
{if(this.outOfMemory){return;}
var popGeo=shape._cf.geometry.node;if(popGeo.hasIndex()){shape._webgl.popGeometry=1;shape._webgl.buffers[0]=gl.createBuffer();gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER,shape._webgl.buffers[0]);gl.bufferData(gl.ELEMENT_ARRAY_BUFFER,popGeo.getTotalNumberOfIndices()*2,gl.STATIC_DRAW);shape._webgl.buffers[5]=gl.createBuffer();var idBuffer=new Float32Array(popGeo._vf.vertexBufferSize);(function(){for(var i=0;i<idBuffer.length;++i)idBuffer[i]=i;})();gl.bindBuffer(gl.ARRAY_BUFFER,shape._webgl.buffers[5]);gl.bufferData(gl.ARRAY_BUFFER,idBuffer,gl.STATIC_DRAW);}
else{shape._webgl.popGeometry=-1;}
shape._webgl.buffers[1]=gl.createBuffer();gl.bindBuffer(gl.ARRAY_BUFFER,shape._webgl.buffers[1]);gl.bufferData(gl.ARRAY_BUFFER,(popGeo._vf.attributeStride*popGeo._vf.vertexBufferSize),gl.STATIC_DRAW);var attribTypeStr=popGeo._vf.coordType;shape._webgl.coordType=x3dom.Utils.getVertexAttribType(attribTypeStr,gl);shape._coordStrideOffset[0]=popGeo.getAttributeStride();shape._coordStrideOffset[1]=popGeo.getPositionOffset();gl.vertexAttribPointer(sp.position,shape._cf.geometry.node._mesh._numPosComponents,shape._webgl.coordType,false,shape._coordStrideOffset[0],shape._coordStrideOffset[1]);gl.enableVertexAttribArray(sp.position);if(popGeo.hasNormal()){attribTypeStr=popGeo._vf.normalType;shape._webgl.normalType=x3dom.Utils.getVertexAttribType(attribTypeStr,gl);shape._normalStrideOffset[0]=popGeo.getAttributeStride();shape._normalStrideOffset[1]=popGeo.getNormalOffset();shape._webgl.buffers[2]=shape._webgl.buffers[1];gl.vertexAttribPointer(sp.normal,shape._cf.geometry.node._mesh._numNormComponents,shape._webgl.normalType,false,shape._normalStrideOffset[0],shape._normalStrideOffset[1]);gl.enableVertexAttribArray(sp.normal);}
if(popGeo.hasTexCoord()){attribTypeStr=popGeo._vf.texCoordType;shape._webgl.texCoordType=x3dom.Utils.getVertexAttribType(attribTypeStr,gl);shape._webgl.buffers[3]=shape._webgl.buffers[1];shape._texCoordStrideOffset[0]=popGeo.getAttributeStride();shape._texCoordStrideOffset[1]=popGeo.getTexCoordOffset();gl.vertexAttribPointer(sp.texcoord,shape._cf.geometry.node._mesh._numTexComponents,shape._webgl.texCoordType,false,shape._texCoordStrideOffset[0],shape._texCoordStrideOffset[1]);gl.enableVertexAttribArray(sp.texcoord);}
if(popGeo.hasColor()){attribTypeStr=popGeo._vf.colorType;shape._webgl.colorType=x3dom.Utils.getVertexAttribType(attribTypeStr,gl);shape._webgl.buffers[4]=shape._webgl.buffers[1];shape._colorStrideOffset[0]=popGeo.getAttributeStride();shape._colorStrideOffset[1]=popGeo.getColorOffset();gl.vertexAttribPointer(sp.color,shape._cf.geometry.node._mesh._numColComponents,shape._webgl.colorType,false,shape._colorStrideOffset[0],shape._colorStrideOffset[1]);gl.enableVertexAttribArray(sp.color);}
shape._webgl.currentNumIndices=0;shape._webgl.currentNumVertices=0;shape._webgl.numVerticesAtLevel=[];shape._webgl.levelsAvailable=0;this.checkError(gl);shape._webgl.levelLoaded=[];(function(){for(var i=0;i<popGeo.getNumLevels();++i)
shape._webgl.levelLoaded.push(false);})();var uploadDataToGPU=function(data,lvl){shape._webgl.levelLoaded[lvl]=true;shape._webgl.numVerticesAtLevel[lvl]=0;if(data){var indexDataLengthInBytes=0;var redrawNeeded=false;if(popGeo.hasIndex()){indexDataLengthInBytes=popGeo.getNumIndicesByLevel(lvl)*2;if(indexDataLengthInBytes>0){redrawNeeded=true;var indexDataView=new Uint8Array(data,0,indexDataLengthInBytes);gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER,shape._webgl.buffers[0]);(function(){var indexDataOffset=0;for(var i=0;i<lvl;++i){indexDataOffset+=popGeo.getNumIndicesByLevel(i);}
gl.bufferSubData(gl.ELEMENT_ARRAY_BUFFER,indexDataOffset*2,indexDataView);})();}}
var vertexDataLengthInBytes=data.byteLength-indexDataLengthInBytes;if(vertexDataLengthInBytes>0){redrawNeeded=true;var attributeDataView=new Uint8Array(data,indexDataLengthInBytes,vertexDataLengthInBytes);gl.bindBuffer(gl.ARRAY_BUFFER,shape._webgl.buffers[1]);if(!popGeo.hasIndex()){gl.bufferSubData(gl.ARRAY_BUFFER,shape._webgl.currentNumVertices*popGeo.getAttributeStride(),attributeDataView);}
else{gl.bufferSubData(gl.ARRAY_BUFFER,popGeo.getVertexDataBufferOffset(lvl)*popGeo.getAttributeStride(),attributeDataView);}
shape._webgl.numVerticesAtLevel[lvl]=vertexDataLengthInBytes/popGeo.getAttributeStride();shape._webgl.currentNumVertices+=shape._webgl.numVerticesAtLevel[lvl];}
(function(){var numValidIndices=0;for(var i=shape._webgl.levelsAvailable;i<popGeo.getNumLevels();++i){if(shape._webgl.levelLoaded[i]===false){break;}
else{numValidIndices+=popGeo.getNumIndicesByLevel(i);++shape._webgl.levelsAvailable;}}
shape._webgl.currentNumIndices=numValidIndices;})();popGeo._mesh._numCoords=shape._webgl.currentNumVertices;popGeo._mesh._numFaces=(popGeo.hasIndex()?shape._webgl.currentNumIndices:shape._webgl.currentNumVertices)/3;popGeo.adaptVertexCount(popGeo.hasIndex()?popGeo._mesh._numFaces*3:popGeo._mesh._numCoords);if(redrawNeeded){shape._nameSpace.doc.needRender=true;}}};var dataURLs=popGeo.getDataURLs();var downloadCallbacks=[];var priorities=[];shape._webgl.downloadStartTimer=new Date().getTime();for(var i=0;i<dataURLs.length;++i){shape._nameSpace.doc.downloadCount+=1;(function(idx){downloadCallbacks.push(function(data){shape._nameSpace.doc.downloadCount-=1;return uploadDataToGPU(data,idx);});})(i);priorities.push(i);}
x3dom.DownloadManager.get(dataURLs,downloadCallbacks,priorities);};x3dom.BinaryContainerLoader.setupImgGeo=function(shape,sp,gl,viewarea,currContext)
{if(this.outOfMemory){return;}
var imageGeometry=shape._cf.geometry.node;if(imageGeometry.getIndexTexture()){shape._webgl.imageGeometry=1;}else{shape._webgl.imageGeometry=-1;}
imageGeometry.unsetGeoDirty();if(currContext.IG_PositionBuffer==null){currContext.IG_PositionBuffer=gl.createBuffer();}
shape._webgl.buffers[1]=currContext.IG_PositionBuffer;gl.bindBuffer(gl.ARRAY_BUFFER,currContext.IG_PositionBuffer);var vertices=new Float32Array(shape._webgl.positions[0]);gl.bufferData(gl.ARRAY_BUFFER,vertices,gl.STATIC_DRAW);gl.bindBuffer(gl.ARRAY_BUFFER,currContext.IG_PositionBuffer);gl.vertexAttribPointer(sp.position,imageGeometry._mesh._numPosComponents,shape._webgl.coordType,false,shape._coordStrideOffset[0],shape._coordStrideOffset[1]);gl.enableVertexAttribArray(sp.position);vertices=null;this.checkError(gl);};x3dom.DrawableCollection=function(drawableCollectionConfig){this.collection=[];this.viewMatrix=drawableCollectionConfig.viewMatrix;this.projMatrix=drawableCollectionConfig.projMatrix;this.sceneMatrix=drawableCollectionConfig.sceneMatrix;this.viewarea=drawableCollectionConfig.viewArea;var scene=this.viewarea._scene;var env=scene.getEnvironment();var viewpoint=scene.getViewpoint();this.near=viewpoint.getNear();this.pixelHeightAtDistOne=viewpoint.getImgPlaneHeightAtDistOne()/this.viewarea._height;this.context=drawableCollectionConfig.context;this.gl=drawableCollectionConfig.gl;this.viewFrustum=this.viewarea.getViewfrustum(this.sceneMatrix);this.worldVol=new x3dom.fields.BoxVolume();this.frustumCulling=drawableCollectionConfig.frustumCulling&&(this.viewFrustum!=null);this.smallFeatureThreshold=drawableCollectionConfig.smallFeatureThreshold;this.sortOpaque=(this.smallFeatureThreshold>0&&env._lowPriorityThreshold<1);this.sortTrans=drawableCollectionConfig.sortTrans;this.prioLevels=10;this.maxTreshold=100;this.sortBySortKey=false;this.sortByPriority=false;this.numberOfNodes=0;this.length=0;};x3dom.DrawableCollection.prototype.cull=function(transform,graphState,singlePath,planeMask){var node=graphState.boundedNode;if(!node||!node._vf.render){return-1;}
var volume=node.getVolume();var MASK_SET=63;if(this.frustumCulling&&graphState.needCulling){var wvol;if(singlePath&&!graphState.worldVolume.isValid()){graphState.worldVolume.transformFrom(transform,volume);wvol=graphState.worldVolume;}
else if(planeMask<MASK_SET){this.worldVol.transformFrom(transform,volume);wvol=this.worldVol;}
if(planeMask<MASK_SET)
planeMask=this.viewFrustum.intersect(wvol,planeMask);if(planeMask==-1)
{return-1;}}
else{planeMask=MASK_SET;}
graphState.coverage=-1;if(this.smallFeatureThreshold>0||node.forceUpdateCoverage()){var modelViewMat=this.viewMatrix.mult(transform);graphState.center=modelViewMat.multMatrixPnt(volume.getCenter());var rVec=modelViewMat.multMatrixVec(volume.getRadialVec());var r=rVec.length();var dist=Math.max(-graphState.center.z-r,this.near);var projPixelLength=dist*this.pixelHeightAtDistOne;graphState.coverage=(r*2.0)/projPixelLength;if(this.smallFeatureThreshold>0&&graphState.coverage<this.smallFeatureThreshold&&graphState.needCulling){return-1;}}
this.numberOfNodes++;return planeMask;};x3dom.DrawableCollection.prototype.addShape=function(shape,transform,graphState){var drawable={};drawable.shape=shape;drawable.transform=transform;drawable.localTransform=graphState.localMatrix;drawable.localVolume=graphState.volume;drawable.worldVolume=x3dom.fields.BoxVolume.copy(graphState.worldVolume);drawable.priority=Math.max(0,graphState.coverage);drawable.shaderID=shape.getShaderProperties(this.viewarea).id;var appearance=shape._cf.appearance.node;drawable.sortType=appearance?appearance._vf.sortType.toLowerCase():"opaque";drawable.sortKey=appearance?appearance._vf.sortKey:0;if(drawable.sortType=='transparent'){if(this.smallFeatureThreshold>0){drawable.zPos=graphState.center.z;}
else{var center=transform.multMatrixPnt(shape.getCenter());center=this.viewMatrix.multMatrixPnt(center);drawable.zPos=center.z;}}
if(!this.sortBySortKey&&drawable.sortKey!=0){this.sortBySortKey=true;}
if(this.collection[drawable.sortType]===undefined){this.collection[drawable.sortType]=[];}
this.collection[drawable.sortType].push(drawable);this.length++;if(this.context&&this.gl){this.context.setupShape(this.gl,drawable,this.viewarea);}};x3dom.DrawableCollection.prototype.addDrawable=function(drawable){drawable.shaderID=drawable.shape.getShaderProperties(this.viewarea).id;var appearance=drawable.shape._cf.appearance.node;drawable.sortType=appearance?appearance._vf.sortType.toLowerCase():"opaque";drawable.sortKey=appearance?appearance._vf.sortKey:0;if(drawable.sortType=='transparent'){var center=drawable.transform.multMatrixPnt(drawable.shape.getCenter());center=this.viewMatrix.multMatrixPnt(center);drawable.zPos=center.z;}
if(!this.sortBySortKey&&drawable.sortKey!=0){this.sortBySortKey=true;}
if(this.collection[drawable.sortType]===undefined){this.collection[drawable.sortType]=[];}
this.collection[drawable.sortType].push(drawable);this.length++;if(this.context&&this.gl){this.context.setupShape(this.gl,drawable,this.viewarea);}};x3dom.DrawableCollection.prototype.calculatePriority=function(graphState){var priority=Math.max(0,graphState.coverage);var pl=this.prioLevels-1;priority=Math.min(Math.round(priority/(this.maxTreshold/pl)),pl);return priority;};x3dom.DrawableCollection.prototype.concat=function(){var opaque=(this.collection['opaque']!==undefined)?this.collection['opaque']:[];var transparent=(this.collection['transparent']!==undefined)?this.collection['transparent']:[];this.collection=opaque.concat(transparent);};x3dom.DrawableCollection.prototype.get=function(idx){return this.collection[idx];};x3dom.DrawableCollection.prototype.sort=function(){var opaque=[];var transparent=[];var that=this;if(this.collection['opaque']!==undefined){if(this.sortOpaque){this.collection['opaque'].sort(function(a,b){if(a.sortKey==b.sortKey||!that.sortBySortKey){return b.priority-a.priority;}
return a.sortKey-b.sortKey;});}
opaque=this.collection['opaque'];}
if(this.collection['transparent']!==undefined){if(this.sortTrans){this.collection['transparent'].sort(function(a,b){if(a.sortKey==b.sortKey||!that.sortBySortKey){if(a.priority==b.priority||!that.sortByPriority){return a.zPos-b.zPos;}
return b.priority-a.priority;}
return a.sortKey-b.sortKey;});}
transparent=this.collection['transparent'];}
this.collection=opaque.concat(transparent);};x3dom.DrawableCollection.prototype.forEach=function(fnc,maxPriority){maxPriority=(maxPriority!==undefined)?Math.min(maxPriority,this.prioLevels):this.prioLevels;var sortKey,priority,shaderID,drawable;for(sortKey=0;sortKey<this.collection['opaque'].length;++sortKey)
{if(this.collection['opaque'][sortKey]!==undefined)
{for(priority=this.collection['opaque'][sortKey].length;priority>0;--priority)
{if(this.collection['opaque'][sortKey][priority]!==undefined)
{for(shaderID in this.collection['opaque'][sortKey][priority])
{for(drawable=0;drawable<this.collection['opaque'][sortKey][priority][shaderID].length;++drawable)
{fnc(this.collection['opaque'][sortKey][priority][shaderID][drawable]);}}}}}}
for(sortKey=0;sortKey<this.collection['transparent'].length;++sortKey)
{if(this.collection['transparent'][sortKey]!==undefined)
{for(priority=this.collection['transparent'][sortKey].length;priority>0;--priority)
{if(this.collection['transparent'][sortKey][priority]!==undefined)
{for(var shaderId in this.collection['transparent'][sortKey][priority])
{this.collection['transparent'][sortKey][priority][shaderId].sort(function(a,b){return a.zPos-b.zPos});for(drawable=0;drawable<this.collection['transparent'][sortKey][priority][shaderId].length;++drawable)
{fnc(this.collection['transparent'][sortKey][priority][shaderId][drawable]);}}}}}}};x3dom.Moveable=function(x3domElem,boundedObj,callback,gridSize,mode){this._x3domRoot=x3domElem;this._runtime=x3domElem.runtime;this._callback=callback;this._gridSize=gridSize?gridSize:0;this._moveable=boundedObj;this._drag=false;this._w=0;this._h=0;this._uPlane=null;this._vPlane=null;this._pPlane=null;this._isect=null;this._translationOffset=null;this._rotationOffset=null;this._scaleOffset=null;this._lastX=0;this._lastY=0;this._buttonState=0;this._mode=(mode&&mode.length)?mode.toLowerCase():"translation";this._firstRay=null;this._matrixTrafo=null;this._navType="examine";this.attachHandlers();};x3dom.Moveable.prototype.setGridSize=function(gridSize){this._gridSize=gridSize;};x3dom.Moveable.prototype.setMode=function(mode){this._mode=mode.toLowerCase();};x3dom.Moveable.prototype.attachHandlers=function(){this._moveable._iMove=this;if(!this._x3domRoot._iMove)
this._x3domRoot._iMove=[];this._x3domRoot._iMove.push(this);this._moveable.addEventListener('mousedown',this.start,false);this._moveable.addEventListener('mouseover',this.over,false);this._moveable.addEventListener('mouseout',this.out,false);if(this._x3domRoot._iMove.length==1){this._x3domRoot.addEventListener('mouseup',this.stop,false);this._x3domRoot.addEventListener('mouseout',this.stop,false);this._x3domRoot.addEventListener('mousemove',this.move,true);if(!this._runtime.canvas.disableTouch){this._x3domRoot.addEventListener('MozTouchDown',this.touchStartHandlerMoz,false);this._x3domRoot.addEventListener('MozTouchMove',this.touchMoveHandlerMoz,true);this._x3domRoot.addEventListener('MozTouchUp',this.touchEndHandlerMoz,false);this._x3domRoot.addEventListener('touchstart',this.touchStartHandler,false);this._x3domRoot.addEventListener('touchmove',this.touchMoveHandler,true);this._x3domRoot.addEventListener('touchend',this.touchEndHandler,false);}}};x3dom.Moveable.prototype.detachHandlers=function(){var iMove=this._x3domRoot._iMove;if(iMove){for(var i=0,n=iMove.length;i<n;i++){if(iMove[i]==this){iMove.splice(i,1);break;}}}
this._moveable.removeEventListener('mousedown',this.start,false);this._moveable.removeEventListener('mouseover',this.over,false);this._moveable.removeEventListener('mouseout',this.out,false);if(iMove.length==0){this._x3domRoot.removeEventListener('mouseup',this.stop,false);this._x3domRoot.removeEventListener('mouseout',this.stop,false);this._x3domRoot.removeEventListener('mousemove',this.move,true);if(!this._runtime.canvas.disableTouch){this._x3domRoot.removeEventListener('MozTouchDown',this.touchStartHandlerMoz,false);this._x3domRoot.removeEventListener('MozTouchMove',this.touchMoveHandlerMoz,true);this._x3domRoot.removeEventListener('MozTouchUp',this.touchEndHandlerMoz,false);this._x3domRoot.removeEventListener('touchstart',this.touchStartHandler,false);this._x3domRoot.removeEventListener('touchmove',this.touchMoveHandler,true);this._x3domRoot.removeEventListener('touchend',this.touchEndHandler,false);}}
if(this._moveable._iMove)
delete this._moveable._iMove;};x3dom.Moveable.prototype.calcViewPlane=function(origin){this._w=this._runtime.getWidth();this._h=this._runtime.getHeight();var ray=this._runtime.getViewingRay(0,this._h-1);var r=ray.pos.add(ray.dir);ray=this._runtime.getViewingRay(this._w-1,this._h-1);var s=ray.pos.add(ray.dir);ray=this._runtime.getViewingRay(0,0);var t=ray.pos.add(ray.dir);this._uPlane=s.subtract(r).normalize();this._vPlane=t.subtract(r).normalize();if(arguments.length===0)
this._pPlane=r;else
this._pPlane=x3dom.fields.SFVec3f.copy(origin);};x3dom.Moveable.prototype.det=function(mat){return mat[0][0]*mat[1][1]*mat[2][2]+mat[0][1]*mat[1][2]*mat[2][0]+
mat[0][2]*mat[2][1]*mat[1][0]-mat[2][0]*mat[1][1]*mat[0][2]-
mat[0][0]*mat[2][1]*mat[1][2]-mat[1][0]*mat[0][1]*mat[2][2];};x3dom.Moveable.prototype.translateXY=function(l){var track=null;var z=[],n=[];for(var i=0;i<3;i++){z[i]=[];n[i]=[];z[i][0]=this._uPlane.at(i);n[i][0]=z[i][0];z[i][1]=this._vPlane.at(i);n[i][1]=z[i][1];z[i][2]=(l.pos.subtract(this._pPlane)).at(i);n[i][2]=-l.dir.at(i);}
var s=this.det(n);if(s!==0){var t=this.det(z)/s;track=l.pos.addScaled(l.dir,t);}
if(track){if(this._isect){track=track.subtract(this._isect);}
track=track.add(this._translationOffset);}
return track;};x3dom.Moveable.prototype.translateZ=function(l,currY){var vol=this._runtime.getSceneBBox();var sign=(currY<this._lastY)?1:-1;var fact=sign*(vol.max.subtract(vol.min)).length()/100;this._translationOffset=this._translationOffset.addScaled(l.dir,fact);return this._translationOffset;};x3dom.Moveable.prototype.rotate=function(posX,posY){var twoPi=2*Math.PI;var alpha=((posY-this._lastY)*twoPi)/this._w;var beta=((posX-this._lastX)*twoPi)/this._h;var q=x3dom.fields.Quaternion.axisAngle(this._uPlane,alpha);var h=q.toMatrix();this._rotationOffset=h.mult(this._rotationOffset);q=x3dom.fields.Quaternion.axisAngle(this._vPlane,beta);h=q.toMatrix();this._rotationOffset=h.mult(this._rotationOffset);var mat=this._rotationOffset.mult(x3dom.fields.SFMatrix4f.scale(this._scaleOffset));var rot=new x3dom.fields.Quaternion(0,0,1,0);rot.setValue(mat);return rot;};x3dom.Moveable.prototype.over=function(event){var that=this._iMove;that._runtime.getCanvas().style.cursor="crosshair";};x3dom.Moveable.prototype.out=function(event){var that=this._iMove;if(!that._drag)
that._runtime.getCanvas().style.cursor="pointer";};x3dom.Moveable.prototype.start=function(event){var that=this._iMove;switch(that._mode){case"translation":that._buttonState=(event.button==4)?1:(event.button&3);break;case"rotation":that._buttonState=4;break;case"all":default:that._buttonState=event.button;break;}
if(!that._drag&&that._buttonState){that._lastX=event.layerX;that._lastY=event.layerY;that._drag=true;that._navType=that._runtime.navigationType();that._runtime.noNav();that._isect=new x3dom.fields.SFVec3f(event.worldX,event.worldY,event.worldZ);that.calcViewPlane(that._isect);that._firstRay=that._runtime.getViewingRay(event.layerX,event.layerY);var mTrans=that._moveable.getAttribute("translation");that._matrixTrafo=null;if(mTrans){that._translationOffset=x3dom.fields.SFVec3f.parse(mTrans);var mRot=that._moveable.getAttribute("rotation");mRot=mRot?x3dom.fields.Quaternion.parseAxisAngle(mRot):new x3dom.fields.Quaternion(0,0,1,0);that._rotationOffset=mRot.toMatrix();var mScal=that._moveable.getAttribute("scale");that._scaleOffset=mScal?x3dom.fields.SFVec3f.parse(mScal):new x3dom.fields.SFVec3f(1,1,1);}
else{mTrans=that._moveable.getAttribute("matrix");if(mTrans){that._matrixTrafo=x3dom.fields.SFMatrix4f.parse(mTrans).transpose();var translation=new x3dom.fields.SFVec3f(0,0,0),scaleFactor=new x3dom.fields.SFVec3f(1,1,1);var rotation=new x3dom.fields.Quaternion(0,0,1,0),scaleOrientation=new x3dom.fields.Quaternion(0,0,1,0);that._matrixTrafo.getTransform(translation,rotation,scaleFactor,scaleOrientation);that._translationOffset=translation;that._rotationOffset=rotation.toMatrix();that._scaleOffset=scaleFactor;}
else{that._translationOffset=new x3dom.fields.SFVec3f(0,0,0);that._rotationOffset=new x3dom.fields.SFMatrix4f();that._scaleOffset=new x3dom.fields.SFVec3f(1,1,1);}}
that._runtime.getCanvas().style.cursor="crosshair";}};x3dom.Moveable.prototype.move=function(event){for(var i=0,n=this._iMove.length;i<n;i++){var that=this._iMove[i];if(that._drag){var pos=that._runtime.mousePosition(event);var ray=that._runtime.getViewingRay(pos[0],pos[1]);var track=null;if(that._buttonState==2)
track=that.translateZ(that._firstRay,pos[1]);else if(that._buttonState==1)
track=that.translateXY(ray);else
track=that.rotate(pos[0],pos[1]);if(track){if(that._gridSize>0&&that._buttonState!=4){var x=that._gridSize*Math.round(track.x/that._gridSize);var y=that._gridSize*Math.round(track.y/that._gridSize);var z=that._gridSize*Math.round(track.z/that._gridSize);track=new x3dom.fields.SFVec3f(x,y,z);}
if(!that._matrixTrafo){if(that._buttonState==4){that._moveable.setAttribute("rotation",track.toAxisAngle().toString());}
else{that._moveable.setAttribute("translation",track.toString());}}
else{if(that._buttonState==4){that._matrixTrafo.setRotate(track);}
else{that._matrixTrafo.setTranslate(track);}
that._moveable.setAttribute("matrix",that._matrixTrafo.toGL().toString());}
if(that._callback){that._callback(that._moveable,track);}}
that._lastX=pos[0];that._lastY=pos[1];}}};x3dom.Moveable.prototype.stop=function(event){for(var i=0,n=this._iMove.length;i<n;i++){var that=this._iMove[i];if(that._drag){that._lastX=event.layerX;that._lastY=event.layerY;that._isect=null;that._drag=false;var navi=that._runtime.canvas.doc._scene.getNavigationInfo();navi.setType(that._navType);that._runtime.getCanvas().style.cursor="pointer";}}};x3dom.Moveable.prototype.touchStartHandler=function(evt){evt.preventDefault();};x3dom.Moveable.prototype.touchStartHandlerMoz=function(evt){evt.preventDefault();};x3dom.Moveable.prototype.touchMoveHandler=function(evt){evt.preventDefault();};x3dom.Moveable.prototype.touchMoveHandlerMoz=function(evt){evt.preventDefault();};x3dom.Moveable.prototype.touchEndHandler=function(evt){if(this._iMove.length){var that=this._iMove[0];that.stop.apply(that._x3domRoot,[evt]);}
evt.preventDefault();};x3dom.Moveable.prototype.touchEndHandlerMoz=function(evt){if(this._iMove.length){var that=this._iMove[0];that.stop.apply(that._x3domRoot,[evt]);}
evt.preventDefault();};(function(){'use strict';function q(b){throw b;}var t=void 0,u=!0,aa=this;function A(b,a){var c=b.split("."),d=aa;!(c[0]in d)&&d.execScript&&d.execScript("var "+c[0]);for(var e;c.length&&(e=c.shift());)!c.length&&a!==t?d[e]=a:d=d[e]?d[e]:d[e]={}};var B="undefined"!==typeof Uint8Array&&"undefined"!==typeof Uint16Array&&"undefined"!==typeof Uint32Array&&"undefined"!==typeof DataView;function F(b,a){this.index="number"===typeof a?a:0;this.m=0;this.buffer=b instanceof(B?Uint8Array:Array)?b:new(B?Uint8Array:Array)(32768);2*this.buffer.length<=this.index&&q(Error("invalid index"));this.buffer.length<=this.index&&this.f()}F.prototype.f=function(){var b=this.buffer,a,c=b.length,d=new(B?Uint8Array:Array)(c<<1);if(B)d.set(b);else for(a=0;a<c;++a)d[a]=b[a];return this.buffer=d};F.prototype.d=function(b,a,c){var d=this.buffer,e=this.index,f=this.m,g=d[e],k;c&&1<a&&(b=8<a?(H[b&255]<<24|H[b>>>8&255]<<16|H[b>>>16&255]<<8|H[b>>>24&255])>>32-a:H[b]>>8-a);if(8>a+f)g=g<<a|b,f+=a;else for(k=0;k<a;++k)g=g<<1|b>>a-k-1&1,8===++f&&(f=0,d[e++]=H[g],g=0,e===d.length&&(d=this.f()));d[e]=g;this.buffer=d;this.m=f;this.index=e};F.prototype.finish=function(){var b=this.buffer,a=this.index,c;0<this.m&&(b[a]<<=8-this.m,b[a]=H[b[a]],a++);B?c=b.subarray(0,a):(b.length=a,c=b);return c};var ba=new(B?Uint8Array:Array)(256),ca;for(ca=0;256>ca;++ca){for(var K=ca,da=K,ea=7,K=K>>>1;K;K>>>=1)da<<=1,da|=K&1,--ea;ba[ca]=(da<<ea&255)>>>0}var H=ba;function ja(b,a,c){var d,e="number"===typeof a?a:a=0,f="number"===typeof c?c:b.length;d=-1;for(e=f&7;e--;++a)d=d>>>8^O[(d^b[a])&255];for(e=f>>3;e--;a+=8)d=d>>>8^O[(d^b[a])&255],d=d>>>8^O[(d^b[a+1])&255],d=d>>>8^O[(d^b[a+2])&255],d=d>>>8^O[(d^b[a+3])&255],d=d>>>8^O[(d^b[a+4])&255],d=d>>>8^O[(d^b[a+5])&255],d=d>>>8^O[(d^b[a+6])&255],d=d>>>8^O[(d^b[a+7])&255];return(d^4294967295)>>>0}
var ka=[0,1996959894,3993919788,2567524794,124634137,1886057615,3915621685,2657392035,249268274,2044508324,3772115230,2547177864,162941995,2125561021,3887607047,2428444049,498536548,1789927666,4089016648,2227061214,450548861,1843258603,4107580753,2211677639,325883990,1684777152,4251122042,2321926636,335633487,1661365465,4195302755,2366115317,997073096,1281953886,3579855332,2724688242,1006888145,1258607687,3524101629,2768942443,901097722,1119000684,3686517206,2898065728,853044451,1172266101,3705015759,2882616665,651767980,1373503546,3369554304,3218104598,565507253,1454621731,3485111705,3099436303,671266974,1594198024,3322730930,2970347812,795835527,1483230225,3244367275,3060149565,1994146192,31158534,2563907772,4023717930,1907459465,112637215,2680153253,3904427059,2013776290,251722036,2517215374,3775830040,2137656763,141376813,2439277719,3865271297,1802195444,476864866,2238001368,4066508878,1812370925,453092731,2181625025,4111451223,1706088902,314042704,2344532202,4240017532,1658658271,366619977,2362670323,4224994405,1303535960,984961486,2747007092,3569037538,1256170817,1037604311,2765210733,3554079995,1131014506,879679996,2909243462,3663771856,1141124467,855842277,2852801631,3708648649,1342533948,654459306,3188396048,3373015174,1466479909,544179635,3110523913,3462522015,1591671054,702138776,2966460450,3352799412,1504918807,783551873,3082640443,3233442989,3988292384,2596254646,62317068,1957810842,3939845945,2647816111,81470997,1943803523,3814918930,2489596804,225274430,2053790376,3826175755,2466906013,167816743,2097651377,4027552580,2265490386,503444072,1762050814,4150417245,2154129355,426522225,1852507879,4275313526,2312317920,282753626,1742555852,4189708143,2394877945,397917763,1622183637,3604390888,2714866558,953729732,1340076626,3518719985,2797360999,1068828381,1219638859,3624741850,2936675148,906185462,1090812512,3747672003,2825379669,829329135,1181335161,3412177804,3160834842,628085408,1382605366,3423369109,3138078467,570562233,1426400815,3317316542,2998733608,733239954,1555261956,3268935591,3050360625,752459403,1541320221,2607071920,3965973030,1969922972,40735498,2617837225,3943577151,1913087877,83908371,2512341634,3803740692,2075208622,213261112,2463272603,3855990285,2094854071,198958881,2262029012,4057260610,1759359992,534414190,2176718541,4139329115,1873836001,414664567,2282248934,4279200368,1711684554,285281116,2405801727,4167216745,1634467795,376229701,2685067896,3608007406,1308918612,956543938,2808555105,3495958263,1231636301,1047427035,2932959818,3654703836,1088359270,936918E3,2847714899,3736837829,1202900863,817233897,3183342108,3401237130,1404277552,615818150,3134207493,3453421203,1423857449,601450431,3009837614,3294710456,1567103746,711928724,3020668471,3272380065,1510334235,755167117],O=B?new Uint32Array(ka):ka;function P(){}P.prototype.getName=function(){return this.name};P.prototype.getData=function(){return this.data};P.prototype.Y=function(){return this.Z};A("Zlib.GunzipMember",P);A("Zlib.GunzipMember.prototype.getName",P.prototype.getName);A("Zlib.GunzipMember.prototype.getData",P.prototype.getData);A("Zlib.GunzipMember.prototype.getMtime",P.prototype.Y);function la(b){this.buffer=new(B?Uint16Array:Array)(2*b);this.length=0}la.prototype.getParent=function(b){return 2*((b-2)/4|0)};la.prototype.push=function(b,a){var c,d,e=this.buffer,f;c=this.length;e[this.length++]=a;for(e[this.length++]=b;0<c;)if(d=this.getParent(c),e[c]>e[d])f=e[c],e[c]=e[d],e[d]=f,f=e[c+1],e[c+1]=e[d+1],e[d+1]=f,c=d;else break;return this.length};la.prototype.pop=function(){var b,a,c=this.buffer,d,e,f;a=c[0];b=c[1];this.length-=2;c[0]=c[this.length];c[1]=c[this.length+1];for(f=0;;){e=2*f+2;if(e>=this.length)break;e+2<this.length&&c[e+2]>c[e]&&(e+=2);if(c[e]>c[f])d=c[f],c[f]=c[e],c[e]=d,d=c[f+1],c[f+1]=c[e+1],c[e+1]=d;else break;f=e}return{index:b,value:a,length:this.length}};function ma(b){var a=b.length,c=0,d=Number.POSITIVE_INFINITY,e,f,g,k,h,l,s,p,m,n;for(p=0;p<a;++p)b[p]>c&&(c=b[p]),b[p]<d&&(d=b[p]);e=1<<c;f=new(B?Uint32Array:Array)(e);g=1;k=0;for(h=2;g<=c;){for(p=0;p<a;++p)if(b[p]===g){l=0;s=k;for(m=0;m<g;++m)l=l<<1|s&1,s>>=1;n=g<<16|p;for(m=l;m<e;m+=h)f[m]=n;++k}++g;k<<=1;h<<=1}return[f,c,d]};function na(b,a){this.k=qa;this.I=0;this.input=B&&b instanceof Array?new Uint8Array(b):b;this.b=0;a&&(a.lazy&&(this.I=a.lazy),"number"===typeof a.compressionType&&(this.k=a.compressionType),a.outputBuffer&&(this.a=B&&a.outputBuffer instanceof Array?new Uint8Array(a.outputBuffer):a.outputBuffer),"number"===typeof a.outputIndex&&(this.b=a.outputIndex));this.a||(this.a=new(B?Uint8Array:Array)(32768))}var qa=2,ra={NONE:0,v:1,o:qa,ba:3},sa=[],S;for(S=0;288>S;S++)switch(u){case 143>=S:sa.push([S+48,8]);break;case 255>=S:sa.push([S-144+400,9]);break;case 279>=S:sa.push([S-256+0,7]);break;case 287>=S:sa.push([S-280+192,8]);break;default:q("invalid literal: "+S)}
na.prototype.g=function(){var b,a,c,d,e=this.input;switch(this.k){case 0:c=0;for(d=e.length;c<d;){a=B?e.subarray(c,c+65535):e.slice(c,c+65535);c+=a.length;var f=a,g=c===d,k=t,h=t,l=t,s=t,p=t,m=this.a,n=this.b;if(B){for(m=new Uint8Array(this.a.buffer);m.length<=n+f.length+5;)m=new Uint8Array(m.length<<1);m.set(this.a)}k=g?1:0;m[n++]=k|0;h=f.length;l=~h+65536&65535;m[n++]=h&255;m[n++]=h>>>8&255;m[n++]=l&255;m[n++]=l>>>8&255;if(B)m.set(f,n),n+=f.length,m=m.subarray(0,n);else{s=0;for(p=f.length;s<p;++s)m[n++]=f[s];m.length=n}this.b=n;this.a=m}break;case 1:var r=new F(B?new Uint8Array(this.a.buffer):this.a,this.b);r.d(1,1,u);r.d(1,2,u);var v=ta(this,e),x,Q,y;x=0;for(Q=v.length;x<Q;x++)if(y=v[x],F.prototype.d.apply(r,sa[y]),256<y)r.d(v[++x],v[++x],u),r.d(v[++x],5),r.d(v[++x],v[++x],u);else if(256===y)break;this.a=r.finish();this.b=this.a.length;break;case qa:var E=new F(B?new Uint8Array(this.a.buffer):this.a,this.b),Ka,R,X,Y,Z,pb=[16,17,18,0,8,7,9,6,10,5,11,4,12,3,13,2,14,1,15],fa,La,ga,Ma,oa,wa=Array(19),Na,$,pa,C,Oa;Ka=qa;E.d(1,1,u);E.d(Ka,2,u);R=ta(this,e);fa=ua(this.W,15);La=va(fa);ga=ua(this.V,7);Ma=va(ga);for(X=286;257<X&&0===fa[X-1];X--);for(Y=30;1<Y&&0===ga[Y-1];Y--);var Pa=X,Qa=Y,J=new(B?Uint32Array:Array)(Pa+Qa),w,L,z,ha,I=new(B?Uint32Array:Array)(316),G,D,M=new(B?Uint8Array:Array)(19);for(w=L=0;w<Pa;w++)J[L++]=fa[w];for(w=0;w<Qa;w++)J[L++]=ga[w];if(!B){w=0;for(ha=M.length;w<ha;++w)M[w]=0}w=G=0;for(ha=J.length;w<ha;w+=L){for(L=1;w+L<ha&&J[w+L]===J[w];++L);z=L;if(0===J[w])if(3>z)for(;0<z--;)I[G++]=0,M[0]++;else for(;0<z;)D=138>z?z:138,D>z-3&&D<z&&(D=z-3),10>=D?(I[G++]=17,I[G++]=D-3,M[17]++):(I[G++]=18,I[G++]=D-11,M[18]++),z-=D;else if(I[G++]=J[w],M[J[w]]++,z--,3>z)for(;0<z--;)I[G++]=J[w],M[J[w]]++;else for(;0<z;)D=6>z?z:6,D>z-3&&D<z&&(D=z-3),I[G++]=16,I[G++]=D-3,M[16]++,z-=D}b=B?I.subarray(0,G):I.slice(0,G);oa=ua(M,7);for(C=0;19>C;C++)wa[C]=oa[pb[C]];for(Z=19;4<Z&&0===wa[Z-1];Z--);Na=va(oa);E.d(X-257,5,u);E.d(Y-1,5,u);E.d(Z-4,4,u);for(C=0;C<Z;C++)E.d(wa[C],3,u);C=0;for(Oa=b.length;C<Oa;C++)if($=b[C],E.d(Na[$],oa[$],u),16<=$){C++;switch($){case 16:pa=2;break;case 17:pa=3;break;case 18:pa=7;break;default:q("invalid code: "+$)}E.d(b[C],pa,u)}var Ra=[La,fa],Sa=[Ma,ga],N,Ta,ia,za,Ua,Va,Wa,Xa;Ua=Ra[0];Va=Ra[1];Wa=Sa[0];Xa=Sa[1];N=0;for(Ta=R.length;N<Ta;++N)if(ia=R[N],E.d(Ua[ia],Va[ia],u),256<ia)E.d(R[++N],R[++N],u),za=R[++N],E.d(Wa[za],Xa[za],u),E.d(R[++N],R[++N],u);else if(256===ia)break;this.a=E.finish();this.b=this.a.length;break;default:q("invalid compression type")}return this.a};function xa(b,a){this.length=b;this.Q=a}
var ya=function(){function b(a){switch(u){case 3===a:return[257,a-3,0];case 4===a:return[258,a-4,0];case 5===a:return[259,a-5,0];case 6===a:return[260,a-6,0];case 7===a:return[261,a-7,0];case 8===a:return[262,a-8,0];case 9===a:return[263,a-9,0];case 10===a:return[264,a-10,0];case 12>=a:return[265,a-11,1];case 14>=a:return[266,a-13,1];case 16>=a:return[267,a-15,1];case 18>=a:return[268,a-17,1];case 22>=a:return[269,a-19,2];case 26>=a:return[270,a-23,2];case 30>=a:return[271,a-27,2];case 34>=a:return[272,a-31,2];case 42>=a:return[273,a-35,3];case 50>=a:return[274,a-43,3];case 58>=a:return[275,a-51,3];case 66>=a:return[276,a-59,3];case 82>=a:return[277,a-67,4];case 98>=a:return[278,a-83,4];case 114>=a:return[279,a-99,4];case 130>=a:return[280,a-115,4];case 162>=a:return[281,a-131,5];case 194>=a:return[282,a-163,5];case 226>=a:return[283,a-195,5];case 257>=a:return[284,a-227,5];case 258===a:return[285,a-258,0];default:q("invalid length: "+a)}}var a=[],c,d;for(c=3;258>=c;c++)d=b(c),a[c]=d[2]<<24|d[1]<<16|d[0];return a}(),Aa=B?new Uint32Array(ya):ya;function ta(b,a){function c(a,c){var b=a.Q,d=[],e=0,f;f=Aa[a.length];d[e++]=f&65535;d[e++]=f>>16&255;d[e++]=f>>24;var g;switch(u){case 1===b:g=[0,b-1,0];break;case 2===b:g=[1,b-2,0];break;case 3===b:g=[2,b-3,0];break;case 4===b:g=[3,b-4,0];break;case 6>=b:g=[4,b-5,1];break;case 8>=b:g=[5,b-7,1];break;case 12>=b:g=[6,b-9,2];break;case 16>=b:g=[7,b-13,2];break;case 24>=b:g=[8,b-17,3];break;case 32>=b:g=[9,b-25,3];break;case 48>=b:g=[10,b-33,4];break;case 64>=b:g=[11,b-49,4];break;case 96>=b:g=[12,b-
65,5];break;case 128>=b:g=[13,b-97,5];break;case 192>=b:g=[14,b-129,6];break;case 256>=b:g=[15,b-193,6];break;case 384>=b:g=[16,b-257,7];break;case 512>=b:g=[17,b-385,7];break;case 768>=b:g=[18,b-513,8];break;case 1024>=b:g=[19,b-769,8];break;case 1536>=b:g=[20,b-1025,9];break;case 2048>=b:g=[21,b-1537,9];break;case 3072>=b:g=[22,b-2049,10];break;case 4096>=b:g=[23,b-3073,10];break;case 6144>=b:g=[24,b-4097,11];break;case 8192>=b:g=[25,b-6145,11];break;case 12288>=b:g=[26,b-8193,12];break;case 16384>=b:g=[27,b-12289,12];break;case 24576>=b:g=[28,b-16385,13];break;case 32768>=b:g=[29,b-24577,13];break;default:q("invalid distance")}f=g;d[e++]=f[0];d[e++]=f[1];d[e++]=f[2];var h,k;h=0;for(k=d.length;h<k;++h)m[n++]=d[h];v[d[0]]++;x[d[3]]++;r=a.length+c-1;p=null}var d,e,f,g,k,h={},l,s,p,m=B?new Uint16Array(2*a.length):[],n=0,r=0,v=new(B?Uint32Array:Array)(286),x=new(B?Uint32Array:Array)(30),Q=b.I,y;if(!B){for(f=0;285>=f;)v[f++]=0;for(f=0;29>=f;)x[f++]=0}v[256]=1;d=0;for(e=a.length;d<e;++d){f=k=0;for(g=3;f<g&&d+f!==e;++f)k=k<<8|a[d+f];h[k]===t&&(h[k]=[]);l=h[k];if(!(0<r--)){for(;0<l.length&&32768<d-l[0];)l.shift();if(d+3>=e){p&&c(p,-1);f=0;for(g=e-d;f<g;++f)y=a[d+f],m[n++]=y,++v[y];break}0<l.length?(s=Ba(a,d,l),p?p.length<s.length?(y=a[d-1],m[n++]=y,++v[y],c(s,0)):c(p,-1):s.length<Q?p=s:c(s,0)):p?c(p,-1):(y=a[d],m[n++]=y,++v[y])}l.push(d)}m[n++]=256;v[256]++;b.W=v;b.V=x;return B?m.subarray(0,n):m}
function Ba(b,a,c){var d,e,f=0,g,k,h,l,s=b.length;k=0;l=c.length;a:for(;k<l;k++){d=c[l-k-1];g=3;if(3<f){for(h=f;3<h;h--)if(b[d+h-1]!==b[a+h-1])continue a;g=f}for(;258>g&&a+g<s&&b[d+g]===b[a+g];)++g;g>f&&(e=d,f=g);if(258===g)break}return new xa(f,a-e)}
function ua(b,a){var c=b.length,d=new la(572),e=new(B?Uint8Array:Array)(c),f,g,k,h,l;if(!B)for(h=0;h<c;h++)e[h]=0;for(h=0;h<c;++h)0<b[h]&&d.push(h,b[h]);f=Array(d.length/2);g=new(B?Uint32Array:Array)(d.length/2);if(1===f.length)return e[d.pop().index]=1,e;h=0;for(l=d.length/2;h<l;++h)f[h]=d.pop(),g[h]=f[h].value;k=Ca(g,g.length,a);h=0;for(l=f.length;h<l;++h)e[f[h].index]=k[h];return e}
function Ca(b,a,c){function d(b){var c=h[b][l[b]];c===a?(d(b+1),d(b+1)):--g[c];++l[b]}var e=new(B?Uint16Array:Array)(c),f=new(B?Uint8Array:Array)(c),g=new(B?Uint8Array:Array)(a),k=Array(c),h=Array(c),l=Array(c),s=(1<<c)-a,p=1<<c-1,m,n,r,v,x;e[c-1]=a;for(n=0;n<c;++n)s<p?f[n]=0:(f[n]=1,s-=p),s<<=1,e[c-2-n]=(e[c-1-n]/2|0)+a;e[0]=f[0];k[0]=Array(e[0]);h[0]=Array(e[0]);for(n=1;n<c;++n)e[n]>2*e[n-1]+f[n]&&(e[n]=2*e[n-1]+f[n]),k[n]=Array(e[n]),h[n]=Array(e[n]);for(m=0;m<a;++m)g[m]=c;for(r=0;r<e[c-1];++r)k[c-
1][r]=b[r],h[c-1][r]=r;for(m=0;m<c;++m)l[m]=0;1===f[c-1]&&(--g[0],++l[c-1]);for(n=c-2;0<=n;--n){v=m=0;x=l[n+1];for(r=0;r<e[n];r++)v=k[n+1][x]+k[n+1][x+1],v>b[m]?(k[n][r]=v,h[n][r]=a,x+=2):(k[n][r]=b[m],h[n][r]=m,++m);l[n]=0;1===f[n]&&d(n)}return g}
function va(b){var a=new(B?Uint16Array:Array)(b.length),c=[],d=[],e=0,f,g,k,h;f=0;for(g=b.length;f<g;f++)c[b[f]]=(c[b[f]]|0)+1;f=1;for(g=16;f<=g;f++)d[f]=e,e+=c[f]|0,e<<=1;f=0;for(g=b.length;f<g;f++){e=d[b[f]];d[b[f]]+=1;k=a[f]=0;for(h=b[f];k<h;k++)a[f]=a[f]<<1|e&1,e>>>=1}return a};function Da(b,a){this.input=b;this.b=this.c=0;this.i={};a&&(a.flags&&(this.i=a.flags),"string"===typeof a.filename&&(this.filename=a.filename),"string"===typeof a.comment&&(this.A=a.comment),a.deflateOptions&&(this.l=a.deflateOptions));this.l||(this.l={})}
Da.prototype.g=function(){var b,a,c,d,e,f,g,k,h=new(B?Uint8Array:Array)(32768),l=0,s=this.input,p=this.c,m=this.filename,n=this.A;h[l++]=31;h[l++]=139;h[l++]=8;b=0;this.i.fname&&(b|=Ea);this.i.fcomment&&(b|=Fa);this.i.fhcrc&&(b|=Ga);h[l++]=b;a=(Date.now?Date.now():+new Date)/1E3|0;h[l++]=a&255;h[l++]=a>>>8&255;h[l++]=a>>>16&255;h[l++]=a>>>24&255;h[l++]=0;h[l++]=Ha;if(this.i.fname!==t){g=0;for(k=m.length;g<k;++g)f=m.charCodeAt(g),255<f&&(h[l++]=f>>>8&255),h[l++]=f&255;h[l++]=0}if(this.i.comment){g=0;for(k=n.length;g<k;++g)f=n.charCodeAt(g),255<f&&(h[l++]=f>>>8&255),h[l++]=f&255;h[l++]=0}this.i.fhcrc&&(c=ja(h,0,l)&65535,h[l++]=c&255,h[l++]=c>>>8&255);this.l.outputBuffer=h;this.l.outputIndex=l;e=new na(s,this.l);h=e.g();l=e.b;B&&(l+8>h.buffer.byteLength?(this.a=new Uint8Array(l+8),this.a.set(new Uint8Array(h.buffer)),h=this.a):h=new Uint8Array(h.buffer));d=ja(s,t,t);h[l++]=d&255;h[l++]=d>>>8&255;h[l++]=d>>>16&255;h[l++]=d>>>24&255;k=s.length;h[l++]=k&255;h[l++]=k>>>8&255;h[l++]=k>>>16&255;h[l++]=k>>>24&255;this.c=p;B&&l<h.length&&(this.a=h=h.subarray(0,l));return h};var Ha=255,Ga=2,Ea=8,Fa=16;A("Zlib.Gzip",Da);A("Zlib.Gzip.prototype.compress",Da.prototype.g);function T(b,a){this.p=[];this.q=32768;this.e=this.j=this.c=this.u=0;this.input=B?new Uint8Array(b):b;this.w=!1;this.r=Ia;this.M=!1;if(a||!(a={}))a.index&&(this.c=a.index),a.bufferSize&&(this.q=a.bufferSize),a.bufferType&&(this.r=a.bufferType),a.resize&&(this.M=a.resize);switch(this.r){case Ja:this.b=32768;this.a=new(B?Uint8Array:Array)(32768+this.q+258);break;case Ia:this.b=0;this.a=new(B?Uint8Array:Array)(this.q);this.f=this.U;this.B=this.R;this.s=this.T;break;default:q(Error("invalid inflate mode"))}}
var Ja=0,Ia=1,Ya={O:Ja,N:Ia};T.prototype.h=function(){for(;!this.w;){var b=U(this,3);b&1&&(this.w=u);b>>>=1;switch(b){case 0:var a=this.input,c=this.c,d=this.a,e=this.b,f=a.length,g=t,k=t,h=d.length,l=t;this.e=this.j=0;c+1>=f&&q(Error("invalid uncompressed block header: LEN"));g=a[c++]|a[c++]<<8;c+1>=f&&q(Error("invalid uncompressed block header: NLEN"));k=a[c++]|a[c++]<<8;g===~k&&q(Error("invalid uncompressed block header: length verify"));c+g>a.length&&q(Error("input buffer is broken"));switch(this.r){case Ja:for(;e+g>d.length;){l=h-e;g-=l;if(B)d.set(a.subarray(c,c+l),e),e+=l,c+=l;else for(;l--;)d[e++]=a[c++];this.b=e;d=this.f();e=this.b}break;case Ia:for(;e+g>d.length;)d=this.f({F:2});break;default:q(Error("invalid inflate mode"))}if(B)d.set(a.subarray(c,c+g),e),e+=g,c+=g;else for(;g--;)d[e++]=a[c++];this.c=c;this.b=e;this.a=d;break;case 1:this.s(Za,$a);break;case 2:ab(this);break;default:q(Error("unknown BTYPE: "+b))}}return this.B()};var bb=[16,17,18,0,8,7,9,6,10,5,11,4,12,3,13,2,14,1,15],cb=B?new Uint16Array(bb):bb,db=[3,4,5,6,7,8,9,10,11,13,15,17,19,23,27,31,35,43,51,59,67,83,99,115,131,163,195,227,258,258,258],eb=B?new Uint16Array(db):db,fb=[0,0,0,0,0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,0,0,0],gb=B?new Uint8Array(fb):fb,hb=[1,2,3,4,5,7,9,13,17,25,33,49,65,97,129,193,257,385,513,769,1025,1537,2049,3073,4097,6145,8193,12289,16385,24577],ib=B?new Uint16Array(hb):hb,jb=[0,0,0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13],kb=B?new Uint8Array(jb):jb,lb=new(B?Uint8Array:Array)(288),V,mb;V=0;for(mb=lb.length;V<mb;++V)lb[V]=143>=V?8:255>=V?9:279>=V?7:8;var Za=ma(lb),nb=new(B?Uint8Array:Array)(30),ob,qb;ob=0;for(qb=nb.length;ob<qb;++ob)nb[ob]=5;var $a=ma(nb);function U(b,a){for(var c=b.j,d=b.e,e=b.input,f=b.c,g=e.length,k;d<a;)f>=g&&q(Error("input buffer is broken")),c|=e[f++]<<d,d+=8;k=c&(1<<a)-1;b.j=c>>>a;b.e=d-a;b.c=f;return k}
function rb(b,a){for(var c=b.j,d=b.e,e=b.input,f=b.c,g=e.length,k=a[0],h=a[1],l,s;d<h&&!(f>=g);)c|=e[f++]<<d,d+=8;l=k[c&(1<<h)-1];s=l>>>16;b.j=c>>s;b.e=d-s;b.c=f;return l&65535}
function ab(b){function a(a,b,c){var d,e=this.J,f,g;for(g=0;g<a;)switch(d=rb(this,b),d){case 16:for(f=3+U(this,2);f--;)c[g++]=e;break;case 17:for(f=3+U(this,3);f--;)c[g++]=0;e=0;break;case 18:for(f=11+U(this,7);f--;)c[g++]=0;e=0;break;default:e=c[g++]=d}this.J=e;return c}var c=U(b,5)+257,d=U(b,5)+1,e=U(b,4)+4,f=new(B?Uint8Array:Array)(cb.length),g,k,h,l;for(l=0;l<e;++l)f[cb[l]]=U(b,3);if(!B){l=e;for(e=f.length;l<e;++l)f[cb[l]]=0}g=ma(f);k=new(B?Uint8Array:Array)(c);h=new(B?Uint8Array:Array)(d);b.J=0;b.s(ma(a.call(b,c,g,k)),ma(a.call(b,d,g,h)))}T.prototype.s=function(b,a){var c=this.a,d=this.b;this.C=b;for(var e=c.length-258,f,g,k,h;256!==(f=rb(this,b));)if(256>f)d>=e&&(this.b=d,c=this.f(),d=this.b),c[d++]=f;else{g=f-257;h=eb[g];0<gb[g]&&(h+=U(this,gb[g]));f=rb(this,a);k=ib[f];0<kb[f]&&(k+=U(this,kb[f]));d>=e&&(this.b=d,c=this.f(),d=this.b);for(;h--;)c[d]=c[d++-k]}for(;8<=this.e;)this.e-=8,this.c--;this.b=d};T.prototype.T=function(b,a){var c=this.a,d=this.b;this.C=b;for(var e=c.length,f,g,k,h;256!==(f=rb(this,b));)if(256>f)d>=e&&(c=this.f(),e=c.length),c[d++]=f;else{g=f-257;h=eb[g];0<gb[g]&&(h+=U(this,gb[g]));f=rb(this,a);k=ib[f];0<kb[f]&&(k+=U(this,kb[f]));d+h>e&&(c=this.f(),e=c.length);for(;h--;)c[d]=c[d++-k]}for(;8<=this.e;)this.e-=8,this.c--;this.b=d};T.prototype.f=function(){var b=new(B?Uint8Array:Array)(this.b-32768),a=this.b-32768,c,d,e=this.a;if(B)b.set(e.subarray(32768,b.length));else{c=0;for(d=b.length;c<d;++c)b[c]=e[c+32768]}this.p.push(b);this.u+=b.length;if(B)e.set(e.subarray(a,a+32768));else for(c=0;32768>c;++c)e[c]=e[a+c];this.b=32768;return e};T.prototype.U=function(b){var a,c=this.input.length/this.c+1|0,d,e,f,g=this.input,k=this.a;b&&("number"===typeof b.F&&(c=b.F),"number"===typeof b.P&&(c+=b.P));2>c?(d=(g.length-this.c)/this.C[2],f=258*(d/2)|0,e=f<k.length?k.length+f:k.length<<1):e=k.length*c;B?(a=new Uint8Array(e),a.set(k)):a=k;return this.a=a};T.prototype.B=function(){var b=0,a=this.a,c=this.p,d,e=new(B?Uint8Array:Array)(this.u+(this.b-32768)),f,g,k,h;if(0===c.length)return B?this.a.subarray(32768,this.b):this.a.slice(32768,this.b);f=0;for(g=c.length;f<g;++f){d=c[f];k=0;for(h=d.length;k<h;++k)e[b++]=d[k]}f=32768;for(g=this.b;f<g;++f)e[b++]=a[f];this.p=[];return this.buffer=e};T.prototype.R=function(){var b,a=this.b;B?this.M?(b=new Uint8Array(a),b.set(this.a.subarray(0,a))):b=this.a.subarray(0,a):(this.a.length>a&&(this.a.length=a),b=this.a);return this.buffer=b};function sb(b){this.input=b;this.c=0;this.t=[];this.D=!1}sb.prototype.X=function(){this.D||this.h();return this.t.slice()};sb.prototype.h=function(){for(var b=this.input.length;this.c<b;){var a=new P,c=t,d=t,e=t,f=t,g=t,k=t,h=t,l=t,s=t,p=this.input,m=this.c;a.G=p[m++];a.H=p[m++];(31!==a.G||139!==a.H)&&q(Error("invalid file signature:"+a.G+","+a.H));a.z=p[m++];switch(a.z){case 8:break;default:q(Error("unknown compression method: "+a.z))}a.n=p[m++];l=p[m++]|p[m++]<<8|p[m++]<<16|p[m++]<<24;a.Z=new Date(1E3*l);a.fa=p[m++];a.ea=p[m++];0<(a.n&4)&&(a.aa=p[m++]|p[m++]<<8,m+=a.aa);if(0<(a.n&Ea)){h=[];for(k=0;0<(g=p[m++]);)h[k++]=String.fromCharCode(g);a.name=h.join("")}if(0<(a.n&Fa)){h=[];for(k=0;0<(g=p[m++]);)h[k++]=String.fromCharCode(g);a.A=h.join("")}0<(a.n&Ga)&&(a.S=ja(p,0,m)&65535,a.S!==(p[m++]|p[m++]<<8)&&q(Error("invalid header crc16")));c=p[p.length-4]|p[p.length-3]<<8|p[p.length-2]<<16|p[p.length-1]<<24;p.length-m-4-4<512*c&&(f=c);d=new T(p,{index:m,bufferSize:f});a.data=e=d.h();m=d.c;a.ca=s=(p[m++]|p[m++]<<8|p[m++]<<16|p[m++]<<24)>>>0;ja(e,t,t)!==s&&q(Error("invalid CRC-32 checksum: 0x"+ja(e,t,t).toString(16)+" / 0x"+s.toString(16)));a.da=c=(p[m++]|p[m++]<<8|p[m++]<<16|p[m++]<<24)>>>0;(e.length&4294967295)!==c&&q(Error("invalid input size: "+(e.length&4294967295)+" / "+c));this.t.push(a);this.c=m}this.D=u;var n=this.t,r,v,x=0,Q=0,y;r=0;for(v=n.length;r<v;++r)Q+=n[r].data.length;if(B){y=new Uint8Array(Q);for(r=0;r<v;++r)y.set(n[r].data,x),x+=n[r].data.length}else{y=[];for(r=0;r<v;++r)y[r]=n[r].data;y=Array.prototype.concat.apply([],y)}return y};A("Zlib.Gunzip",sb);A("Zlib.Gunzip.prototype.decompress",sb.prototype.h);A("Zlib.Gunzip.prototype.getMembers",sb.prototype.X);function tb(b){if("string"===typeof b){var a=b.split(""),c,d;c=0;for(d=a.length;c<d;c++)a[c]=(a[c].charCodeAt(0)&255)>>>0;b=a}for(var e=1,f=0,g=b.length,k,h=0;0<g;){k=1024<g?1024:g;g-=k;do e+=b[h++],f+=e;while(--k);e%=65521;f%=65521}return(f<<16|e)>>>0};function ub(b,a){var c,d;this.input=b;this.c=0;if(a||!(a={}))a.index&&(this.c=a.index),a.verify&&(this.$=a.verify);c=b[this.c++];d=b[this.c++];switch(c&15){case vb:this.method=vb;break;default:q(Error("unsupported compression method"))}0!==((c<<8)+d)%31&&q(Error("invalid fcheck flag:"+((c<<8)+d)%31));d&32&&q(Error("fdict flag is not supported"));this.L=new T(b,{index:this.c,bufferSize:a.bufferSize,bufferType:a.bufferType,resize:a.resize})}
ub.prototype.h=function(){var b=this.input,a,c;a=this.L.h();this.c=this.L.c;this.$&&(c=(b[this.c++]<<24|b[this.c++]<<16|b[this.c++]<<8|b[this.c++])>>>0,c!==tb(a)&&q(Error("invalid adler-32 checksum")));return a};var vb=8;function wb(b,a){this.input=b;this.a=new(B?Uint8Array:Array)(32768);this.k=W.o;var c={},d;if((a||!(a={}))&&"number"===typeof a.compressionType)this.k=a.compressionType;for(d in a)c[d]=a[d];c.outputBuffer=this.a;this.K=new na(this.input,c)}var W=ra;wb.prototype.g=function(){var b,a,c,d,e,f,g,k=0;g=this.a;b=vb;switch(b){case vb:a=Math.LOG2E*Math.log(32768)-8;break;default:q(Error("invalid compression method"))}c=a<<4|b;g[k++]=c;switch(b){case vb:switch(this.k){case W.NONE:e=0;break;case W.v:e=1;break;case W.o:e=2;break;default:q(Error("unsupported compression type"))}break;default:q(Error("invalid compression method"))}d=e<<6|0;g[k++]=d|31-(256*c+d)%31;f=tb(this.input);this.K.b=k;g=this.K.g();k=g.length;B&&(g=new Uint8Array(g.buffer),g.length<=k+4&&(this.a=new Uint8Array(g.length+4),this.a.set(g),g=this.a),g=g.subarray(0,k+4));g[k++]=f>>24&255;g[k++]=f>>16&255;g[k++]=f>>8&255;g[k++]=f&255;return g};function xb(b,a){var c,d,e,f;if(Object.keys)c=Object.keys(a);else for(d in c=[],e=0,a)c[e++]=d;e=0;for(f=c.length;e<f;++e)d=c[e],A(b+"."+d,a[d])};A("Zlib.Inflate",ub);A("Zlib.Inflate.prototype.decompress",ub.prototype.h);xb("Zlib.Inflate.BufferType",{ADAPTIVE:Ya.N,BLOCK:Ya.O});A("Zlib.Deflate",wb);A("Zlib.Deflate.compress",function(b,a){return(new wb(b,a)).g()});A("Zlib.Deflate.prototype.compress",wb.prototype.g);xb("Zlib.Deflate.CompressionType",{NONE:W.NONE,FIXED:W.v,DYNAMIC:W.o});}).call(this);if(x3dom.glTF==null)
x3dom.glTF={};x3dom.glTF.glTFLoader=function(response,meshOnly)
{this.meshOnly=meshOnly;this.header=this.readHeader(response);if(this.header.sceneLength>0){this.scene=this.readScene(response,this.header);this.body=this.readBody(response,this.header);}
this._mesh={};};x3dom.glTF.glTFLoader.prototype.getScene=function(shape,shaderProgram,gl,sceneName)
{this.reset(shape,gl);if(sceneName==null)
{sceneName=this.scene["scene"];}
var scene=this.scene.scenes[sceneName];this.updateScene(shape,shaderProgram,gl,scene);};x3dom.glTF.glTFLoader.prototype.getMesh=function(shape,shaderProgram,gl,meshName)
{this.reset(shape,gl);var mesh;if(meshName==null)
{mesh=Object.keys(this.scene.meshes)[0];}else
{for(var key in this.scene.meshes){if(this.scene.meshes.hasOwnProperty(key)&&key==meshName)
{mesh=this.scene.meshes[key];break;}}}
this.updateMesh(shape,shaderProgram,gl,mesh);};x3dom.glTF.glTFLoader.prototype.reset=function(shape,gl)
{this._mesh._numCoords=0;this._mesh._numFaces=0;shape._webgl.externalGeometry=-1;if(this.loaded.bufferViews==null)
this.loaded.bufferViews=this.loadBufferViews(shape,gl);};x3dom.glTF.glTFLoader.prototype.updateScene=function(shape,shaderProgram,gl,scene)
{var nodes=scene["nodes"];for(var i=0;i<nodes.length;++i)
{var nodeID=nodes[i];this.traverseNode(shape,shaderProgram,gl,this.scene.nodes[nodeID]);}};x3dom.glTF.glTFLoader.prototype.traverseNode=function(shape,shaderProgram,gl,node)
{var children=node["children"];if(children!=null)
for(var i=0;i<children.length;++i)
{var childID=children[i];this.traverseNode(shape,shaderProgram,gl,this.scene.nodes[childID]);}
var meshes=node["meshes"];if(meshes!=null&&meshes.length>0)
for(var i=0;i<meshes.length;++i){var meshID=meshes[i];if(this.loaded.meshes[meshID]==null){this.updateMesh(shape,shaderProgram,gl,this.scene.meshes[meshID]);this.loaded.meshes[meshID]=1;}}};x3dom.glTF.glTFLoader.prototype.updateMesh=function(shape,shaderProgram,gl,mesh)
{var primitives=mesh["primitives"];for(var i=0;i<primitives.length;++i){this.loadglTFMesh(shape,shaderProgram,gl,primitives[i]);}};x3dom.glTF.glTFLoader.prototype.loadPrimitive=function(shape,shaderProgram,gl,primitive)
{var INDEX_BUFFER_IDX=0;var POSITION_BUFFER_IDX=1;var NORMAL_BUFFER_IDX=2;var TEXCOORD_BUFFER_IDX=3;var COLOR_BUFFER_IDX=4;var x3domTypeID,x3domShortTypeID;var meshIdx=this.loaded.meshCount;var bufferOffset=meshIdx*6;shape._webgl.primType[meshIdx]=primitive["mode"];var indexed=(primitive.indices!=null&&primitive.indices!="");if(indexed==true){var indicesAccessor=this.scene.accessors[primitive.indices];shape._webgl.indexOffset[meshIdx]=indicesAccessor["byteOffset"];shape._webgl.drawCount[meshIdx]=indicesAccessor["count"];shape._webgl.buffers[INDEX_BUFFER_IDX+bufferOffset]=this.loaded.bufferViews[indicesAccessor["bufferView"]];this._mesh._numFaces+=indicesAccessor["count"]/3;}
var attributes=primitive["attributes"];for(var attributeID in attributes)
{var accessorName=attributes[attributeID];var accessor=this.scene.accessors[accessorName];switch(attributeID)
{case"POSITION":x3domTypeID="coord";x3domShortTypeID="Pos";shape._webgl.buffers[POSITION_BUFFER_IDX+bufferOffset]=this.loaded.bufferViews[accessor["bufferView"]];if(indexed==false)
{shape._webgl.drawCount[meshIdx]=accessor["count"];this._mesh._numFaces+=accessor["count"]/3;}
this._mesh._numCoords+=accessor["count"];break;case"NORMAL":x3domTypeID="normal";x3domShortTypeID="Norm";shape._webgl.buffers[NORMAL_BUFFER_IDX+bufferOffset]=this.loaded.bufferViews[accessor["bufferView"]];break;case"TEXCOORD_0":x3domTypeID="texCoord";x3domShortTypeID="Tex";shape._webgl.buffers[TEXCOORD_BUFFER_IDX+bufferOffset]=this.loaded.bufferViews[accessor["bufferView"]];break;case"COLOR":x3domTypeID="color";x3domShortTypeID="Col";shape._webgl.buffers[COLOR_BUFFER_IDX+bufferOffset]=this.loaded.bufferViews[accessor["bufferView"]];break;}
if(x3domTypeID!=null){shape["_"+x3domTypeID+"StrideOffset"][meshIdx]=[];shape["_"+x3domTypeID+"StrideOffset"][meshIdx][0]=accessor["byteStride"];shape["_"+x3domTypeID+"StrideOffset"][meshIdx][1]=accessor["byteOffset"];shape._webgl[x3domTypeID+"Type"]=accessor["componentType"];this._mesh["_num"+x3domShortTypeID+"Components"]=this.getNumComponentsForType(accessor["type"]);}}
this.loaded.meshCount+=1;shape._dirty.shader=true;shape._nameSpace.doc.needRender=true;x3dom.BinaryContainerLoader.checkError(gl);};x3dom.glTF.glTFLoader.prototype.loadglTFMesh=function(shape,shaderProgram,gl,primitive)
{"use strict";var mesh=new x3dom.glTF.glTFMesh();mesh.primitiveType=primitive["mode"];var indexed=(primitive.indices!=null&&primitive.indices!="");if(indexed==true){var indicesAccessor=this.scene.accessors[primitive.indices];mesh.buffers[glTF_BUFFER_IDX.INDEX]={};mesh.buffers[glTF_BUFFER_IDX.INDEX].offset=indicesAccessor["byteOffset"];mesh.buffers[glTF_BUFFER_IDX.INDEX].type=indicesAccessor["componentType"];mesh.buffers[glTF_BUFFER_IDX.INDEX].idx=this.loaded.bufferViews[indicesAccessor["bufferView"]];mesh.drawCount=indicesAccessor["count"];this._mesh._numFaces+=indicesAccessor["count"]/3;}
var attributes=primitive["attributes"];for(var attributeID in attributes)
{var accessorName=attributes[attributeID];var accessor=this.scene.accessors[accessorName];var idx=null;switch(attributeID)
{case"POSITION":idx=glTF_BUFFER_IDX.POSITION;if(indexed==false)
{mesh.drawCount=accessor["count"];this._mesh._numFaces+=indicesAccessor["count"]/3;}
this._mesh.numCoords+=accessor["count"];break;case"NORMAL":idx=glTF_BUFFER_IDX.NORMAL;break;case"TEXCOORD_0":idx=glTF_BUFFER_IDX.TEXCOORD;break;case"COLOR":idx=glTF_BUFFER_IDX.COLOR;break;}
if(idx!=null){mesh.buffers[idx]={};mesh.buffers[idx].idx=this.loaded.bufferViews[accessor["bufferView"]];mesh.buffers[idx].offset=accessor["byteOffset"];mesh.buffers[idx].stride=accessor["byteStride"];mesh.buffers[idx].type=accessor["componentType"];mesh.buffers[idx].numComponents=this.getNumComponentsForType(accessor["type"]);}}
this.loaded.meshCount+=1;shape._dirty.shader=true;shape._nameSpace.doc.needRender=true;x3dom.BinaryContainerLoader.checkError(gl);if(primitive.material!=null&&!this.meshOnly)
mesh.material=this.loadMaterial(gl,this.scene.materials[primitive.material]);if(shape.meshes==null)
shape.meshes=[];shape.meshes.push(mesh);};x3dom.glTF.glTFLoader.prototype.loadBufferViews=function(shape,gl)
{var buffers={};var bufferViews=this.scene.bufferViews;for(var bufferId in bufferViews)
{if(!bufferViews.hasOwnProperty(bufferId))continue;var bufferView=bufferViews[bufferId];if(bufferView.target==null&&bufferView.target!=gl.ARRAY_BUFFER&&bufferView.target!=gl.ELEMENT_ARRAY_BUFFER)
continue;if(bufferView.target==gl.ELEMENT_ARRAY_BUFFER)
shape._webgl.externalGeometry=1;var data=new Uint8Array(this.body.buffer,this.header.bodyOffset+bufferView["byteOffset"],bufferView["byteLength"]);var newBuffer=gl.createBuffer();gl.bindBuffer(bufferView["target"],newBuffer);gl.bufferData(bufferView["target"],data,gl.STATIC_DRAW);buffers[bufferId]=newBuffer;}
return buffers;};x3dom.glTF.glTFLoader.prototype.readHeader=function(response)
{var header={};var magicBytes=new Uint8Array(response,0,4);var versionBytes=new Uint32Array(response,4,1);var lengthBytes=new Uint32Array(response,8,1);var sceneLengthBytes=new Uint32Array(response,12,1);var sceneFormatBytes=new Uint32Array(response,16,1);header.magic=new TextDecoder("ascii").decode(magicBytes);if(versionBytes[0]==1)
header.version="Version 1";header.length=lengthBytes[0];header.sceneLength=sceneLengthBytes[0];if(sceneFormatBytes[0]==0)
header.sceneFormat="JSON";header.bodyOffset=header.sceneLength+20;return header;};x3dom.glTF.glTFLoader.prototype.readScene=function(response,header)
{var sceneBytes=new Uint8Array(response,20,header.sceneLength);var json=JSON.parse(new TextDecoder("utf-8").decode(sceneBytes));return json;};x3dom.glTF.glTFLoader.prototype.readBody=function(response,header)
{var offset=header.sceneLength+20;var body=new Uint8Array(response,offset,header.length-offset);return body;};x3dom.glTF.glTFLoader.prototype.getNumComponentsForType=function(type)
{switch(type)
{case"SCALAR":return 1;case"VEC2":return 2;case"VEC3":return 3;case"VEC4":return 4;default:return 0;}};x3dom.glTF.glTFLoader.prototype.loadImage=function(imageNodeName,mimeType)
{if(this.loaded.images==null)
this.loaded.images={};if(this.loaded.images[imageNodeName]!=null)
return this.loaded.images[imageNodeName];var imageNode=this.scene.images[imageNodeName];if(imageNode.extensions!=null&&imageNode.extensions.KHR_binary_glTF!=null)
{var ext=imageNode.extensions.KHR_binary_glTF;var bufferView=this.scene.bufferViews[ext.bufferView];var uint8Array=new Uint8Array(this.body.buffer,this.header.bodyOffset+bufferView.byteOffset,bufferView.byteLength);var blob=new Blob([uint8Array],{type:ext.mimeType});var blobUrl=window.URL.createObjectURL(blob);var image=new Image();image.src=blobUrl;this.loaded.images[imageNodeName]=image;return image;}
return null;};x3dom.glTF.glTFLoader.prototype.loadTexture=function(gl,textureNode)
{var format=textureNode.format;var internalFormat=textureNode.internalFormat;var sampler={};var samplerNode=this.scene.samplers[textureNode.sampler];if(samplerNode!=null)
{for(var key in samplerNode){if(samplerNode.hasOwnProperty(key))
sampler[key]=samplerNode[key];}}
var image=this.loadImage(textureNode.source);var target=textureNode.target;var type=textureNode.type;var glTFTexture=new x3dom.glTF.glTFTexture(gl,format,internalFormat,sampler,target,type,image);return glTFTexture;};x3dom.glTF.glTFLoader.prototype.loadMaterial=function(gl,materialNode)
{if(materialNode.extensions!=null&&materialNode.extensions.KHR_materials_common!=null)
{materialNode=materialNode.extensions.KHR_materials_common;var material=new x3dom.glTF.glTFKHRMaterialCommons();material.technique=glTF_KHR_MATERIAL_COMMON_TECHNIQUE[materialNode.technique];material.doubleSided=materialNode.doubleSided;for(var key in materialNode.values)
if(materialNode.values.hasOwnProperty(key))
{var value=materialNode.values[key];if(typeof value==='string')
{var textureNode=this.scene.textures[value];material[key+"Tex"]=this.loadTexture(gl,textureNode);}
else
{material[key]=value;}}
return material;}else
{var technique=this.scene.techniques[materialNode.technique];var program=this.loadShaderProgram(gl,technique.program);var material=new x3dom.glTF.glTFMaterial(technique);material.program=program;for(var key in materialNode.values)
if(materialNode.values.hasOwnProperty(key))
{var value=materialNode.values[key];if(typeof value==='string')
{var textureNode=this.scene.textures[value];material.textures[key]=this.loadTexture(gl,textureNode);}
else
{material.values[key]=value;}}
return material;}
return new x3dom.glTF.glTFKHRMaterialCommons();};x3dom.glTF.glTFLoader.prototype.loadShaderProgram=function(gl,shaderProgramName)
{if(this.loaded.programs==null)
this.loaded.programs={};if(this.loaded.programs[shaderProgramName]!=null)
return this.loaded.programs[shaderProgramName];var shaderProgramNode=this.scene.programs[shaderProgramName];var vertexShaderNode=this.scene.shaders[shaderProgramNode.vertexShader];var vertexShaderSrc=this._loadShaderSource(vertexShaderNode);var fragmentShaderNode=this.scene.shaders[shaderProgramNode.fragmentShader];var fragmentShaderSrc=this._loadShaderSource(fragmentShaderNode);var program=gl.createProgram();var vertexShader=gl.createShader(gl.VERTEX_SHADER);gl.shaderSource(vertexShader,vertexShaderSrc);gl.compileShader(vertexShader);if(!gl.getShaderParameter(vertexShader,gl.COMPILE_STATUS)){x3dom.debug.logError("[glTF binary] VertexShader "+gl.getShaderInfoLog(vertexShader));}
var fragmentShader=gl.createShader(gl.FRAGMENT_SHADER);gl.shaderSource(fragmentShader,fragmentShaderSrc);gl.compileShader(fragmentShader);if(!gl.getShaderParameter(fragmentShader,gl.COMPILE_STATUS)){x3dom.debug.logError("[glTF binary] FragmentShader "+gl.getShaderInfoLog(fragmentShader));}
gl.attachShader(program,vertexShader);gl.attachShader(program,fragmentShader);gl.bindAttribLocation(program,0,"position");gl.linkProgram(program);var program=x3dom.Utils.wrapProgram(gl,program);this.loaded.programs[shaderProgramName]=program;return program;};x3dom.glTF.glTFLoader.prototype._loadShaderSource=function(shaderNode)
{var bufferView=this.scene.bufferViews[shaderNode.extensions.KHR_binary_glTF.bufferView];var shaderBytes=new Uint8Array(this.body.buffer,this.header.bodyOffset+bufferView.byteOffset,bufferView.byteLength);var src=new TextDecoder("ascii").decode(shaderBytes);return src;};if(x3dom.glTF==null)
x3dom.glTF={};glTF_BUFFER_IDX={INDEX:0,POSITION:1,NORMAL:2,TEXCOORD:3,COLOR:4};glTF_KHR_MATERIAL_COMMON_TECHNIQUE={BLINN:0,PHONG:1,LAMBERT:2,CONSTANT:3};x3dom.glTF.glTFMesh=function()
{this.indexOffset=0;this.drawCount=0;this.numFaces=0;this.primitiveType=0;this.numCoords=0;this.buffers={};this.material=null;};x3dom.glTF.glTFMesh.prototype.bindVertexAttribPointer=function(gl,shaderProgram)
{if(this.buffers[glTF_BUFFER_IDX.INDEX]){gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER,this.buffers[glTF_BUFFER_IDX.INDEX].idx);}
if(this.material!=null&&this.material.attributeMapping!=null)
{var mapping=this.material.attributeMapping;this._bindVertexAttribPointer(gl,shaderProgram[mapping[glTF_BUFFER_IDX.POSITION]],this.buffers[glTF_BUFFER_IDX.POSITION]);this._bindVertexAttribPointer(gl,shaderProgram[mapping[glTF_BUFFER_IDX.NORMAL]],this.buffers[glTF_BUFFER_IDX.NORMAL]);this._bindVertexAttribPointer(gl,shaderProgram[mapping[glTF_BUFFER_IDX.TEXCOORD]],this.buffers[glTF_BUFFER_IDX.TEXCOORD]);this._bindVertexAttribPointer(gl,shaderProgram[mapping[glTF_BUFFER_IDX.COLOR]],this.buffers[glTF_BUFFER_IDX.COLOR]);}
else
{this._bindVertexAttribPointer(gl,shaderProgram.position,this.buffers[glTF_BUFFER_IDX.POSITION]);this._bindVertexAttribPointer(gl,shaderProgram.normal,this.buffers[glTF_BUFFER_IDX.NORMAL]);this._bindVertexAttribPointer(gl,shaderProgram.texcoord,this.buffers[glTF_BUFFER_IDX.TEXCOORD]);this._bindVertexAttribPointer(gl,shaderProgram.color,this.buffers[glTF_BUFFER_IDX.COLOR]);}};x3dom.glTF.glTFMesh.prototype.bindVertexAttribPointerPosition=function(gl,shaderProgram,useMaterial)
{if(this.buffers[glTF_BUFFER_IDX.INDEX]){gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER,this.buffers[glTF_BUFFER_IDX.INDEX].idx);}
if(useMaterial==true&&this.material!=null&&this.material.attributeMapping!=null)
{var mapping=this.material.attributeMapping;this._bindVertexAttribPointer(gl,shaderProgram[mapping[glTF_BUFFER_IDX.POSITION]],this.buffers[glTF_BUFFER_IDX.POSITION]);}
else
{this._bindVertexAttribPointer(gl,shaderProgram.position,this.buffers[glTF_BUFFER_IDX.POSITION]);}};x3dom.glTF.glTFMesh.prototype._bindVertexAttribPointer=function(gl,shaderPosition,buffer)
{if(shaderPosition!=null&&buffer!=null)
{gl.bindBuffer(gl.ARRAY_BUFFER,buffer.idx);gl.vertexAttribPointer(shaderPosition,buffer.numComponents,buffer.type,false,buffer.stride,buffer.offset);gl.enableVertexAttribArray(shaderPosition);}};x3dom.glTF.glTFMesh.prototype.render=function(gl,polyMode)
{if(this.material!=null&&!this.material.created())
return;if(polyMode==null)
polyMode=this.primitiveType;if(this.buffers[glTF_BUFFER_IDX.INDEX])
gl.drawElements(polyMode,this.drawCount,this.buffers[glTF_BUFFER_IDX.INDEX].type,this.buffers[glTF_BUFFER_IDX.INDEX].offset);else
gl.drawArrays(polyMode,0,this.drawCount);};x3dom.glTF.glTFTexture=function(gl,format,internalFormat,sampler,target,type,image)
{this.format=format;this.internalFormat=internalFormat;this.sampler=sampler;this.target=target;this.type=type;this.image=image;this.created=false;this.create(gl);};x3dom.glTF.glTFTexture.prototype.isPowerOfTwo=function(x)
{var powerOfTwo=!(x==0)&&!(x&(x-1));return powerOfTwo;};x3dom.glTF.glTFTexture.prototype.create=function(gl)
{if(this.image.complete==false)
return;this.glTexture=gl.createTexture();gl.bindTexture(gl.TEXTURE_2D,this.glTexture);gl.texImage2D(gl.TEXTURE_2D,0,this.internalFormat,this.format,this.type,this.image);if(this.sampler.magFilter!=null)
gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MAG_FILTER,this.sampler.magFilter);if(this.sampler.minFilter!=null)
gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MIN_FILTER,this.sampler.minFilter);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MIN_FILTER,gl.LINEAR);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_S,gl.CLAMP_TO_EDGE);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_T,gl.CLAMP_TO_EDGE);gl.bindTexture(gl.TEXTURE_2D,null);this.created=true;};x3dom.glTF.glTFTexture.prototype.bind=function(gl,textureUnit,shaderProgram,uniformName)
{if(!this.created)
this.create(gl);gl.activeTexture(gl.TEXTURE0+textureUnit);gl.bindTexture(gl.TEXTURE_2D,this.glTexture);gl.uniform1i(gl.getUniformLocation(shaderProgram,uniformName),textureUnit);};x3dom.glTF.glTFKHRMaterialCommons=function()
{this.diffuse=[0.3,0.1,0.1,1];this.diffuseTex=null;this.emission=[0.0,0.0,0.0,1];this.emissionTex=null;this.specular=[0.8,0.8,0.8,1];this.specularTex=null;this.ambient=[0,0,0,1];this.shininess=2;this.transparency=0.0;this.globalAmbient=[0,0,0,1];this.lightVector=[1,0,0,1];this.doubleSided=false;this.technique=glTF_KHR_MATERIAL_COMMON_TECHNIQUE.BLINN;};x3dom.glTF.glTFKHRMaterialCommons.prototype.created=function()
{if(this.diffuseTex!=null&&this.diffuseTex.created!=true)
return false;if(this.emissionTex!=null&&this.emissionTex.created!=true)
return false;if(this.specularTex!=null&&this.specularTex.created!=true)
return false;return true;};x3dom.glTF.glTFKHRMaterialCommons.prototype.setShader=function(gl,cache,shape,properties)
{properties.EMPTY_SHADER=0;properties.KHR_MATERIAL_COMMONS=1;if(this.diffuseTex!=null)
properties.USE_DIFFUSE_TEX=1;else
properties.USE_DIFFUSE_TEX=0;if(this.emissionTex!=null)
properties.USE_SPECULAR_TEX=1;else
properties.USE_SPECULAR_TEX=0;if(this.specularTex!=null)
properties.USE_EMISSION_TEX=1;else
properties.USE_EMISSION_TEX=0;properties.toIdentifier();this.program=cache.getShaderByProperties(gl,shape,properties);};x3dom.glTF.glTFKHRMaterialCommons.prototype.bind=function(gl,shaderProgram)
{this.program.bind();for(var key in shaderProgram){if(!shaderProgram.hasOwnProperty(key))
continue;if(this.program.hasOwnProperty(key))
this.program[key]=shaderProgram[key];}
if(this.diffuseTex!=null)
this.diffuseTex.bind(gl,0,this.program.program,"diffuseTex");else
this.program.diffuse=this.diffuse;if(this.emissionTex!=null)
this.emissionTex.bind(gl,0,this.program.program,"emissionTex");else
this.program.emission=this.emission;if(this.specularTex!=null)
this.specularTex.bind(gl,0,this.program.program,"specularTex");else
this.program.specular=this.specular;this.program.shininess=this.shininess;this.program.transparency=this.transparency;this.program.globalAmbient=this.globalAmbient;this.program.lightVector=this.lightVector;this.program.technique=this.technique;};x3dom.glTF.glTFMaterial=function(technique)
{this.technique=technique;this.values={};this.semanticMapping={};this.attributeMapping={};this.textures={};for(var key in this.technique.uniforms)
{if(this.technique.uniforms.hasOwnProperty(key))
{var parameter=this.technique.parameters[this.technique.uniforms[key]];if(parameter.semantic!=null)
switch(parameter.semantic)
{case"MODELVIEW":this.semanticMapping["modelViewMatrix"]=key;break;case"MODELVIEWINVERSETRANSPOSE":this.semanticMapping["modelViewInverseTransposeMatrix"]=key;break;case"PROJECTION":this.semanticMapping["projectionMatrix"]=key;break;case"MODEL":this.semanticMapping["modelMatrix"]=key;break;case"MODELVIEWPROJECTION":this.semanticMapping["modelViewProjectionMatrix"]=key;break;case"VIEW":this.semanticMapping["viewMatrix"]=key;break;case"MODELVIEWINVERSE":this.semanticMapping["modelViewInverseMatrix"]=key;break;default:break;}}}
for(var key in this.technique.attributes){if(this.technique.attributes.hasOwnProperty(key)){var parameter=this.technique.parameters[this.technique.attributes[key]];if(parameter.semantic!=null)
switch(parameter.semantic){case"POSITION":this.attributeMapping[glTF_BUFFER_IDX.POSITION]=key;break;case"NORMAL":this.attributeMapping[glTF_BUFFER_IDX.NORMAL]=key;break;case"TEXCOORD_0":this.attributeMapping[glTF_BUFFER_IDX.TEXCOORD]=key;break;case"COLOR":this.attributeMapping[glTF_BUFFER_IDX.COLOR]=key;break;default:break;}}}};x3dom.glTF.glTFMaterial.prototype.created=function()
{for(var key in this.textures){if(!this.textures.hasOwnProperty(key))continue;if(this.textures[key].created!=true)
return false;}
return true;};x3dom.glTF.glTFMaterial.prototype.bind=function(gl,shaderParameter)
{if(this.program!=null)
this.program.bind();this.updateTransforms(shaderParameter);for(var key in this.technique.uniforms)
if(this.technique.uniforms.hasOwnProperty(key))
{var uniformName=this.technique.uniforms[key];if(this.textures[uniformName]!=null){var texture=this.textures[uniformName];texture.bind(gl,0,this.program.program,key);}
else if(this.values[uniformName]!=null)
this.program[key]=this.values[uniformName];}};x3dom.glTF.glTFMaterial.prototype.updateTransforms=function(shaderParameter)
{if(this.program!=null)
{this.program.bind();if(this.semanticMapping["modelViewMatrix"]!=null)
this.program[this.semanticMapping["modelViewMatrix"]]=shaderParameter.modelViewMatrix;if(this.semanticMapping["viewMatrix"]!=null)
this.program[this.semanticMapping["viewMatrix"]]=shaderParameter.viewMatrix;if(this.semanticMapping["modelViewInverseTransposeMatrix"]!=null){var mat=shaderParameter.normalMatrix;var model_view_inv_gl=[mat[0],mat[1],mat[2],mat[4],mat[5],mat[6],mat[8],mat[9],mat[10]];this.program[this.semanticMapping["modelViewInverseTransposeMatrix"]]=model_view_inv_gl;}
if(this.semanticMapping["modelViewInverseMatrix"]!=null)
this.program[this.semanticMapping["modelViewInverseMatrix"]]=shaderParameter.modelViewMatrixInverse;if(this.semanticMapping["modelViewProjectionMatrix"]!=null)
this.program[this.semanticMapping["modelViewProjectionMatrix"]]=shaderParameter.modelViewProjectionMatrix;if(this.semanticMapping["modelMatrix"]!=null)
this.program[this.semanticMapping["modelMatrix"]]=shaderParameter.model;if(this.semanticMapping["projectionMatrix"]!=null)
this.program[this.semanticMapping["projectionMatrix"]]=shaderParameter.projectionMatrix;}};x3dom.X3DCanvas=function(x3dElem,canvasIdx)
{var that=this;this._canvasIdx=canvasIdx;this.x3dElem=x3dElem;this._current_dim=[0,0];this.fps_t0=new Date().getTime();this.lastTimeFPSWasTaken=0;this.framesSinceLastTime=0;this._totalTime=0;this._elapsedTime=0;this.doc=null;this.devicePixelRatio=window.devicePixelRatio||1;this.lastMousePos={x:0,y:0};x3dom.caps.DOMNodeInsertedEvent_perSubtree=!(navigator.userAgent.indexOf('MSIE')!=-1||navigator.userAgent.indexOf('Trident')!=-1);x3dElem.__setAttribute=x3dElem.setAttribute;x3dElem.setAttribute=function(attrName,newVal)
{this.__setAttribute(attrName,newVal);newVal=parseInt(newVal)*that.devicePixelRatio;switch(attrName){case"width":that.canvas.setAttribute("width",newVal);if(that.doc&&that.doc._viewarea){that.doc._viewarea._width=parseInt(that.canvas.getAttribute("width"),0);that.doc.needRender=true;}
break;case"height":that.canvas.setAttribute("height",newVal);if(that.doc&&that.doc._viewarea){that.doc._viewarea._height=parseInt(that.canvas.getAttribute("height"),0);that.doc.needRender=true;}
break;default:break;}};x3dom.caps.MOBILE=(navigator.appVersion.indexOf("Mobile")>-1);this.backend=this.x3dElem.getAttribute('backend');this.backend=(this.backend)?this.backend.toLowerCase():'none';this.canvas=this._createHTMLCanvas(x3dElem);this.canvas.parent=this;this.gl=this._initContext(this.canvas,(this.backend.search("desktop")>=0),(this.backend.search("mobile")>=0),(this.backend.search("flashie")>=0),(this.backend.search("webgl2")>=0));this.backend='webgl';if(this.gl==null)
{this.hasRuntime=false;this._createInitFailedDiv(x3dElem);return;}
x3dom.caps.BACKEND=this.backend;var runtimeEnabled=x3dElem.getAttribute("runtimeEnabled");if(runtimeEnabled!==null)
{this.hasRuntime=(runtimeEnabled.toLowerCase()=="true");}
else
{this.hasRuntime=x3dElem.hasRuntime;}
this.showStat=x3dElem.getAttribute("showStat");this.stateViewer=new x3dom.States(x3dElem);if(this.showStat!==null&&this.showStat=="true")
{this.stateViewer.display(true);}
this.x3dElem.appendChild(this.stateViewer.viewer);this.showProgress=x3dElem.getAttribute("showProgress");this.progressDiv=this._createProgressDiv();this.progressDiv.style.display=(this.showProgress!==null&&this.showProgress=="true")?"inline":"none";this.x3dElem.appendChild(this.progressDiv);this.showTouchpoints=x3dElem.getAttribute("showTouchpoints");this.showTouchpoints=this.showTouchpoints?this.showTouchpoints:false;this.disableTouch=x3dElem.getAttribute("disableTouch");this.disableTouch=this.disableTouch?(this.disableTouch.toLowerCase()=="true"):false;this.disableKeys=x3dElem.getAttribute("keysEnabled");this.disableKeys=this.disableKeys?(this.disableKeys.toLowerCase()=="true"):false;this.disableRightDrag=x3dElem.getAttribute("disableRightDrag");this.disableRightDrag=this.disableRightDrag?(this.disableRightDrag.toLowerCase()=="true"):false;this.disableLeftDrag=x3dElem.getAttribute("disableLeftDrag");this.disableLeftDrag=this.disableLeftDrag?(this.disableLeftDrag.toLowerCase()=="true"):false;this.disableMiddleDrag=x3dElem.getAttribute("disableMiddleDrag");this.disableMiddleDrag=this.disableMiddleDrag?(this.disableMiddleDrag.toLowerCase()=="true"):false;this.bindEventListeners();};x3dom.X3DCanvas.prototype.bindEventListeners=function(){var that=this;this.onMouseDown=function(evt){if(!this.isMulti){this.focus();this.classList.add('x3dom-canvas-mousedown');switch(evt.button){case 0:this.mouse_button=1;break;case 1:this.mouse_button=4;break;case 2:this.mouse_button=2;break;default:this.mouse_button=0;break;}
if(evt.shiftKey){this.mouse_button=1;}
if(evt.ctrlKey){this.mouse_button=4;}
if(evt.altKey){this.mouse_button=2;}
var pos=this.parent.mousePosition(evt);this.mouse_drag_x=pos.x;this.mouse_drag_y=pos.y;this.mouse_dragging=true;this.parent.doc.onMousePress(that.gl,this.mouse_drag_x,this.mouse_drag_y,this.mouse_button);this.parent.doc.needRender=true;}}
this.onMouseUp=function(evt){if(!this.isMulti){var prev_mouse_button=this.mouse_button;this.classList.remove('x3dom-canvas-mousedown');this.mouse_button=0;this.mouse_dragging=false;this.parent.doc.onMouseRelease(that.gl,this.mouse_drag_x,this.mouse_drag_y,this.mouse_button,prev_mouse_button);this.parent.doc.needRender=true;}}
this.onMouseOver=function(evt){if(!this.isMulti){this.mouse_button=0;this.mouse_dragging=false;this.parent.doc.onMouseOver(that.gl,this.mouse_drag_x,this.mouse_drag_y,this.mouse_button);this.parent.doc.needRender=true;}}
this.onMouseAlt=function(evt){if(!this.isMulti){this.mouse_button=0;this.mouse_dragging=false;this.classList.remove('x3dom-canvas-mousedown');this.parent.doc.onMouseOut(that.gl,this.mouse_drag_x,this.mouse_drag_y,this.mouse_button);this.parent.doc.needRender=true;}}
this.onDoubleClick=function(evt){if(!this.isMulti){this.mouse_button=0;var pos=this.parent.mousePosition(evt);this.mouse_drag_x=pos.x;this.mouse_drag_y=pos.y;this.mouse_dragging=false;this.parent.doc.onDoubleClick(that.gl,this.mouse_drag_x,this.mouse_drag_y);this.parent.doc.needRender=true;}}
this.onMouseMove=function(evt){if(!this.isMulti){var pos=this.parent.mousePosition(evt);if(pos.x!=that.lastMousePos.x||pos.y!=that.lastMousePos.y){that.lastMousePos=pos;if(evt.shiftKey){this.mouse_button=1;}
if(evt.ctrlKey){this.mouse_button=4;}
if(evt.altKey){this.mouse_button=2;}
this.mouse_drag_x=pos.x;this.mouse_drag_y=pos.y;if(this.mouse_dragging){if(this.mouse_button==1&&!this.parent.disableLeftDrag||this.mouse_button==2&&!this.parent.disableRightDrag||this.mouse_button==4&&!this.parent.disableMiddleDrag)
{this.parent.doc.onDrag(that.gl,this.mouse_drag_x,this.mouse_drag_y,this.mouse_button);}}
else{this.parent.doc.onMove(that.gl,this.mouse_drag_x,this.mouse_drag_y,this.mouse_button);}
this.parent.doc.needRender=true;evt.preventDefault();evt.stopPropagation();}}}
this.onDOMMouseScroll=function(evt){if(!this.isMulti){this.focus();var originalY=this.parent.mousePosition(evt).y;this.mouse_drag_y+=2*evt.detail;this.parent.doc.onWheel(that.gl,this.mouse_drag_x,this.mouse_drag_y,originalY);this.parent.doc.needRender=true;evt.preventDefault();evt.stopPropagation();}}
this.onKeyPress=function(evt){if(!this.parent.disableKeys){this.parent.doc.onKeyPress(evt.charCode);}
this.parent.doc.needRender=true;}
this.onMouseWheel=function(evt){if(!this.isMulti){this.focus();var originalY=this.parent.mousePosition(evt).y;this.mouse_drag_y-=0.1*evt.wheelDelta;this.parent.doc.onWheel(that.gl,this.mouse_drag_x,this.mouse_drag_y,originalY);this.parent.doc.needRender=true;evt.preventDefault();evt.stopPropagation();}}
this.onKeyUp=function(evt){if(!this.parent.disableKeys){this.parent.doc.onKeyUp(evt.keyCode);}
this.parent.doc.needRender=true;}
this.onKeyDown=function(evt){if(!this.parent.disableKeys){this.parent.doc.onKeyDown(evt.keyCode);}
this.parent.doc.needRender=true;}
if(this.canvas!==null&&this.gl!==null&&this.hasRuntime){this.canvas.mouse_dragging=false;this.canvas.mouse_button=0;this.canvas.mouse_drag_x=0;this.canvas.mouse_drag_y=0;this.canvas.isMulti=false;this.canvas.oncontextmenu=function(evt){evt.preventDefault();evt.stopPropagation();return false;};this.canvas.addEventListener("webglcontextlost",function(event){x3dom.debug.logError("WebGL context lost");event.preventDefault();},false);this.canvas.addEventListener("webglcontextrestored",function(event){x3dom.debug.logError("recover WebGL state and resources on context lost NYI");event.preventDefault();},false);this.canvas.addEventListener('mousedown',this.onMouseDown,false);this.canvas.addEventListener('mouseup',this.onMouseUp,false);this.canvas.addEventListener('mouseover',this.onMouseOver,false);this.canvas.addEventListener('mouseout',this.onMouseOut,false);this.canvas.addEventListener('dblclick',this.onDoubleClick,false);this.canvas.addEventListener('mousemove',this.onMouseMove,false);this.canvas.addEventListener('DOMMouseScroll',this.onDOMMouseScroll,false);this.canvas.addEventListener('mousewheel',this.onMouseWheel,false);this.canvas.addEventListener('keypress',this.onKeyPress,true);this.canvas.addEventListener('keyup',this.onKeyUp,true);this.canvas.addEventListener('keydown',this.onKeyDown,true);var touches={numTouches:0,firstTouchTime:new Date().getTime(),firstTouchPoint:new x3dom.fields.SFVec2f(0,0),lastPos:new x3dom.fields.SFVec2f(),lastDrag:new x3dom.fields.SFVec2f(),lastMiddle:new x3dom.fields.SFVec2f(),lastSquareDistance:0,lastAngle:0,lastLayer:[],examineNavType:1,calcAngle:function(vector)
{var rotation=vector.normalize().dot(new x3dom.fields.SFVec2f(1,0));rotation=Math.acos(rotation);if(vector.y<0)
rotation=Math.PI+(Math.PI-rotation);return rotation;},disableTouch:this.disableTouch,visMarker:this.showTouchpoints,visMarkerBag:[],visualizeTouches:function(evt)
{if(!this.visMarker)
return;var touchBag=[];var marker=null;for(var i=0;i<evt.touches.length;i++){var id=evt.touches[i].identifier||evt.touches[i].streamId;if(!id)id=0;var index=this.visMarkerBag.indexOf(id);if(index>=0){marker=document.getElementById("visMarker"+id);marker.style.left=(evt.touches[i].pageX)+"px";marker.style.top=(evt.touches[i].pageY)+"px";}
else{marker=document.createElement("div");marker.appendChild(document.createTextNode("#"+id));marker.id="visMarker"+id;marker.className="x3dom-touch-marker";document.body.appendChild(marker);index=this.visMarkerBag.length;this.visMarkerBag[index]=id;}
touchBag.push(id);}
for(var j=this.visMarkerBag.length-1;j>=0;j--){var oldId=this.visMarkerBag[j];if(touchBag.indexOf(oldId)<0){this.visMarkerBag.splice(j,1);marker=document.getElementById("visMarker"+oldId);document.body.removeChild(marker);}}}};var touchStartHandler=function(evt,doc)
{this.isMulti=true;evt.preventDefault();touches.visualizeTouches(evt);this.focus();if(doc==null)
doc=this.parent.doc;var navi=doc._scene.getNavigationInfo();switch(navi.getType()){case"examine":touches.examineNavType=1;break;case"turntable":touches.examineNavType=2;break;default:touches.examineNavType=0;break;}
touches.lastLayer=[];var i,pos;for(i=0;i<evt.touches.length;i++){pos=this.parent.mousePosition(evt.touches[i]);touches.lastLayer.push([evt.touches[i].identifier,new x3dom.fields.SFVec2f(pos.x,pos.y)]);}
if(touches.numTouches<1&&evt.touches.length==1){touches.numTouches=1;touches.lastDrag=new x3dom.fields.SFVec2f(evt.touches[0].screenX,evt.touches[0].screenY);}
else if(touches.numTouches<2&&evt.touches.length>=2){touches.numTouches=2;var touch0=new x3dom.fields.SFVec2f(evt.touches[0].screenX,evt.touches[0].screenY);var touch1=new x3dom.fields.SFVec2f(evt.touches[1].screenX,evt.touches[1].screenY);var distance=touch1.subtract(touch0);var middle=distance.multiply(0.5).add(touch0);var squareDistance=distance.dot(distance);touches.lastMiddle=middle;touches.lastSquareDistance=squareDistance;touches.lastAngle=touches.calcAngle(distance);touches.lastPos=this.parent.mousePosition(evt.touches[0]);}
doc._scene.updateVolume();if(touches.examineNavType==1){for(i=0;i<evt.touches.length;i++){pos=this.parent.mousePosition(evt.touches[i]);doc.onPick(that.gl,pos.x,pos.y);doc._viewarea.prepareEvents(pos.x,pos.y,1,"onmousedown");doc._viewarea._pickingInfo.lastClickObj=doc._viewarea._pickingInfo.pickObj;}}
else if(evt.touches.length){pos=this.parent.mousePosition(evt.touches[0]);doc.onMousePress(that.gl,pos.x,pos.y,1);}
doc.needRender=true;};var touchMoveHandler=function(evt,doc)
{evt.preventDefault();touches.visualizeTouches(evt);if(doc==null)
doc=this.parent.doc;var pos=null;var rotMatrix=null;var touch0,touch1,distance,middle,squareDistance,deltaMiddle,deltaZoom,deltaMove;if(touches.examineNavType==1){if(evt.touches.length==1){var currentDrag=new x3dom.fields.SFVec2f(evt.touches[0].screenX,evt.touches[0].screenY);var deltaDrag=currentDrag.subtract(touches.lastDrag);touches.lastDrag=currentDrag;var mx=x3dom.fields.SFMatrix4f.rotationY(deltaDrag.x/100);var my=x3dom.fields.SFMatrix4f.rotationX(deltaDrag.y/100);rotMatrix=mx.mult(my);doc.onMoveView(that.gl,evt,touches,null,rotMatrix);}
else if(evt.touches.length>=2){touch0=new x3dom.fields.SFVec2f(evt.touches[0].screenX,evt.touches[0].screenY);touch1=new x3dom.fields.SFVec2f(evt.touches[1].screenX,evt.touches[1].screenY);distance=touch1.subtract(touch0);middle=distance.multiply(0.5).add(touch0);squareDistance=distance.dot(distance);deltaMiddle=middle.subtract(touches.lastMiddle);deltaZoom=squareDistance-touches.lastSquareDistance;deltaMove=new x3dom.fields.SFVec3f(deltaMiddle.x/screen.width,-deltaMiddle.y/screen.height,deltaZoom/(screen.width*screen.height*0.2));var rotation=touches.calcAngle(distance);var angleDelta=touches.lastAngle-rotation;touches.lastAngle=rotation;rotMatrix=x3dom.fields.SFMatrix4f.rotationZ(angleDelta);touches.lastMiddle=middle;touches.lastSquareDistance=squareDistance;doc.onMoveView(that.gl,evt,touches,deltaMove,rotMatrix);}}
else if(evt.touches.length){if(touches.examineNavType==2&&evt.touches.length>=2){touch0=new x3dom.fields.SFVec2f(evt.touches[0].screenX,evt.touches[0].screenY);touch1=new x3dom.fields.SFVec2f(evt.touches[1].screenX,evt.touches[1].screenY);distance=touch1.subtract(touch0);squareDistance=distance.dot(distance);deltaZoom=(squareDistance-touches.lastSquareDistance)/(0.1*(screen.width+screen.height));touches.lastPos.y+=deltaZoom;touches.lastSquareDistance=squareDistance;doc.onDrag(that.gl,touches.lastPos.x,touches.lastPos.y,2);}
else{pos=this.parent.mousePosition(evt.touches[0]);doc.onDrag(that.gl,pos.x,pos.y,1);}}
doc.needRender=true;};var touchEndHandler=function(evt,doc)
{this.isMulti=false;evt.preventDefault();touches.visualizeTouches(evt);if(doc==null)
doc=this.parent.doc;doc._viewarea._isMoving=false;if(touches.numTouches==2&&evt.touches.length==1)
touches.lastDrag=new x3dom.fields.SFVec2f(evt.touches[0].screenX,evt.touches[0].screenY);var dblClick=false;if(evt.touches.length<2){if(touches.numTouches==1)
dblClick=true;touches.numTouches=evt.touches.length;}
if(touches.examineNavType==1){for(var i=0;i<touches.lastLayer.length;i++){var pos=touches.lastLayer[i][1];doc.onPick(that.gl,pos.x,pos.y);if(doc._scene._vf.pickMode.toLowerCase()!=="box"){doc._viewarea.prepareEvents(pos.x,pos.y,1,"onmouseup");doc._viewarea._pickingInfo.lastClickObj=doc._viewarea._pickingInfo.pickObj;if(doc._viewarea._pickingInfo.pickObj&&doc._viewarea._pickingInfo.pickObj===doc._viewarea._pickingInfo.lastClickObj){doc._viewarea.prepareEvents(pos.x,pos.y,1,"onclick");}}
else{var line=doc._viewarea.calcViewRay(pos.x,pos.y);var isect=doc._scene.doIntersect(line);var obj=line.hitObject;if(isect&&obj){doc._viewarea._pick.setValues(line.hitPoint);doc._viewarea.checkEvents(obj,pos.x,pos.y,1,"onclick");x3dom.debug.logInfo("Hit '"+obj._xmlNode.localName+"/ "+
obj._DEF+"' at pos "+doc._viewarea._pick);}}}
if(dblClick){var now=new Date().getTime();var dist=touches.firstTouchPoint.subtract(touches.lastDrag).length();if(dist<18&&now-touches.firstTouchTime<180)
doc.onDoubleClick(that.gl,0,0);touches.firstTouchTime=now;touches.firstTouchPoint=touches.lastDrag;}}
else if(touches.lastLayer.length){pos=touches.lastLayer[0][1];doc.onMouseRelease(that.gl,pos.x,pos.y,0,1);}
doc.needRender=true;};if(!this.disableTouch)
{this.canvas.addEventListener('touchstart',touchStartHandler,true);this.canvas.addEventListener('touchmove',touchMoveHandler,true);this.canvas.addEventListener('touchend',touchEndHandler,true);}}}
x3dom.X3DCanvas.prototype._initContext=function(canvas,forbidMobileShaders,forceMobileShaders,tryWebGL2)
{x3dom.debug.logInfo("Initializing X3DCanvas for ["+canvas.id+"]");var gl=x3dom.gfx_webgl(canvas,forbidMobileShaders,forceMobileShaders,tryWebGL2,this.x3dElem);if(!gl)
{x3dom.debug.logError("No 3D context found...");this.x3dElem.removeChild(canvas);return null;}
else
{var webglVersion=parseFloat(x3dom.caps.VERSION.match(/\d+\.\d+/)[0]);if(webglVersion<1.0){x3dom.debug.logError("WebGL version "+x3dom.caps.VERSION+" lacks important WebGL/GLSL features needed for shadows, special vertex attribute types, etc.!");}}
return gl;};x3dom.X3DCanvas.prototype.appendParam=function(node,name,value){var param=document.createElement('param');param.setAttribute('name',name);param.setAttribute('value',value);node.appendChild(param);};x3dom.X3DCanvas.prototype._createInitFailedDiv=function(x3dElem){var div=document.createElement('div');div.setAttribute("id","x3dom-create-init-failed");div.style.width=x3dElem.getAttribute("width");div.style.height=x3dElem.getAttribute("height");div.style.backgroundColor="#C00";div.style.color="#FFF";div.style.fontSize="20px";div.style.fontWidth="bold";div.style.padding="10px 10px 10px 10px";div.style.display="inline-block";div.style.fontFamily="Helvetica";div.style.textAlign="center";div.appendChild(document.createTextNode('Your Browser does not support X3DOM'));div.appendChild(document.createElement('br'));div.appendChild(document.createTextNode('Read more about Browser support on:'));div.appendChild(document.createElement('br'));var link=document.createElement('a');link.setAttribute('href','http://www.x3dom.org/?page_id=9');link.appendChild(document.createTextNode('X3DOM | Browser Support'));div.appendChild(link);var altImg=x3dElem.getAttribute("altImg")||null;if(altImg){var altImgObj=new Image();altImgObj.src=altImg;div.style.backgroundImage="url("+altImg+")";div.style.backgroundRepeat="no-repeat";div.style.backgroundPosition="50% 50%";}
x3dElem.appendChild(div);x3dom.debug.logError("Your Browser does not support X3DOM!");};x3dom.X3DCanvas.prototype._createHTMLCanvas=function(x3dElem)
{x3dom.debug.logInfo("Creating canvas for (X)3D element...");var canvas=document.createElement('canvas');canvas.setAttribute("class","x3dom-canvas");var userStyle=x3dElem.getAttribute("style");if(userStyle){x3dom.debug.logInfo("Inline X3D styles detected");}
var evtArr=["onmousedown","onmousemove","onmouseout","onmouseover","onmouseup","onclick","ondblclick","onkeydown","onkeypress","onkeyup","ontouchstart","ontouchmove","ontouchend","ontouchcancel","ontouchleave","ontouchenter","ondragstart","ondrop","ondragover"];for(var i=0;i<evtArr.length;i++)
{var evtName=evtArr[i];var userEvt=x3dElem.getAttribute(evtName);if(userEvt){x3dom.debug.logInfo(evtName+", "+userEvt);canvas.setAttribute(evtName,userEvt);x3dElem.removeAttribute(evtName);}}
var userProp=x3dElem.getAttribute("draggable");if(userProp){x3dom.debug.logInfo("draggable="+userProp);canvas.setAttribute("draggable",userProp);}
if(!x3dElem.__addEventListener&&!x3dElem.__removeEventListener)
{x3dElem.__addEventListener=x3dElem.addEventListener;x3dElem.__removeEventListener=x3dElem.removeEventListener;x3dElem.addEventListener=function(type,func,phase){var j,found=false;for(j=0;j<evtArr.length&&!found;j++){if(evtArr[j]===type){found=true;}}
if(found){x3dom.debug.logInfo('addEventListener for div.on'+type);canvas.addEventListener(type,func,phase);}else{x3dom.debug.logInfo('addEventListener for X3D.on'+type);this.__addEventListener(type,func,phase);}};x3dElem.removeEventListener=function(type,func,phase){var j,found=false;for(j=0;j<evtArr.length&&!found;j++){if(evtArr[j]===type){found=true;}}
if(found){x3dom.debug.logInfo('removeEventListener for div.on'+type);canvas.removeEventListener(type,func,phase);}else{x3dom.debug.logInfo('removeEventListener for X3D.on'+type);this.__removeEventListener(type,func,phase);}};}
if(x3dElem.hasAttribute("ondownloadsfinished"))
{x3dElem.addEventListener("downloadsfinished",function()
{var eventObject={target:x3dElem,type:"downloadsfinished"};var funcStr=x3dElem.getAttribute("ondownloadsfinished");var func=new Function('event',funcStr);func.call(x3dElem,eventObject);},true);}
x3dElem.appendChild(canvas);var id=x3dElem.getAttribute("id");if(id!==null){canvas.id="x3dom-"+id+"-canvas";}else{var index=new Date().getTime();canvas.id="x3dom-"+index+"-canvas";}
var w,h;if((w=x3dElem.getAttribute("width"))!==null){if(w.indexOf("%")>=0){x3dom.debug.logWarning("The width attribute is to be specified in pixels not in percent.");}
canvas.style.width=w;canvas.setAttribute("width",w);}
if((h=x3dElem.getAttribute("height"))!==null){if(h.indexOf("%")>=0){x3dom.debug.logWarning("The height attribute is to be specified in pixels not in percent.");}
canvas.style.height=h;canvas.setAttribute("height",h);}
canvas.setAttribute("tabindex","0");return canvas;};x3dom.X3DCanvas.prototype._watchForResize=function(){var new_dim=[parseInt(x3dom.getStyle(this.canvas,"width")),parseInt(x3dom.getStyle(this.canvas,"height"))];if((this._current_dim[0]!=new_dim[0])||(this._current_dim[1]!=new_dim[1])){this._current_dim=new_dim;this.x3dElem.setAttribute("width",new_dim[0]+"px");this.x3dElem.setAttribute("height",new_dim[1]+"px");}};x3dom.X3DCanvas.prototype._createProgressDiv=function(){var progressDiv=document.createElement('div');progressDiv.setAttribute("class","x3dom-progress");var _text=document.createElement('strong');_text.appendChild(document.createTextNode('Loading...'));progressDiv.appendChild(_text);var _inner=document.createElement('span');_inner.setAttribute('style',"width: 25%;");_inner.appendChild(document.createTextNode(' '));progressDiv.appendChild(_inner);progressDiv.oncontextmenu=progressDiv.onmousedown=function(evt){evt.preventDefault();evt.stopPropagation();return false;};return progressDiv;};x3dom.X3DCanvas.prototype.mousePosition=function(evt)
{var rect=evt.target.getBoundingClientRect();var offsetX=Math.round(evt.clientX-rect.left)*this.devicePixelRatio;var offsetY=Math.round(evt.clientY-rect.top)*this.devicePixelRatio;return new x3dom.fields.SFVec2f(offsetX,offsetY);};x3dom.X3DCanvas.prototype.tick=function(timestamp)
{var that=this;this._elapsedTime=(this._totalTime)?timestamp-this._totalTime:0;this._totalTime=timestamp;var runtime=this.x3dElem.runtime;var d=new Date().getTime();var diff=d-this.lastTimeFPSWasTaken;var fps=1000.0/(d-this.fps_t0);this.fps_t0=d;this.doc.advanceTime(d/1000.0);var animD=new Date().getTime()-d;if(this.doc.needRender){if(diff>=1000){runtime.fps=this.framesSinceLastTime/(diff/1000.0);runtime.addMeasurement('FPS',runtime.fps);this.framesSinceLastTime=0;this.lastTimeFPSWasTaken=d;}
this.framesSinceLastTime++;runtime.addMeasurement('ANIM',animD);if(runtime.isReady==false){runtime.ready();runtime.isReady=true;}
runtime.enterFrame({"total":this._totalTime,"elapsed":this._elapsedTime});this.doc.needRender=false;this.doc.render(this.gl);if(!this.doc._scene._vf.doPickPass)
{runtime.removeMeasurement('PICKING');}
runtime.exitFrame({"total":this._totalTime,"elapsed":this._elapsedTime});}
if(this.progressDiv){if(this.doc.downloadCount>0){runtime.addInfo("#LOADS:",this.doc.downloadCount);}else{runtime.removeInfo("#LOADS:");}
if(this.doc.properties.getProperty("showProgress")!=='false'){if(this.progressDiv){this.progressDiv.childNodes[0].textContent='Loading: '+(+this.doc.downloadCount);if(this.doc.downloadCount>0){this.progressDiv.style.display='inline';}else{this.progressDiv.style.display='none';}}}else{this.progressDiv.style.display='none';}}
if(this.doc.downloadCount==0&&this.doc.previousDownloadCount>0)
{var evt;if(document.createEvent){evt=document.createEvent("Events");evt.initEvent("downloadsfinished",true,true);that.x3dElem.dispatchEvent(evt);}else if(document.createEventObject){evt=document.createEventObject();that.x3dElem.fireEvent("ondownloadsfinished",evt);}}
this.doc.previousDownloadCount=this.doc.downloadCount;};x3dom.X3DCanvas.prototype.load=function(uri,sceneElemPos,settings){this.doc=new x3dom.X3DDocument(this.canvas,this.gl,settings);var x3dCanvas=this;this.doc.onload=function(){if(x3dCanvas.hasRuntime){(function mainloop(timestamp){if(x3dCanvas.doc&&x3dCanvas.x3dElem.runtime){x3dCanvas._watchForResize();x3dCanvas.tick(timestamp);window.requestAnimFrame(mainloop,x3dCanvas);}})();}else{x3dCanvas.tick();}};this.x3dElem.render=function(){if(x3dCanvas.hasRuntime){x3dCanvas.doc.needRender=true;}else{x3dCanvas.doc.render(x3dCanvas.gl);}};this.x3dElem.context=x3dCanvas.gl.ctx3d;this.doc.onerror=function(){alert('Failed to load X3D document');};this.doc.load(uri,sceneElemPos);};x3dom.runtime={};x3dom.Runtime=function(doc,canvas){this.doc=doc;this.canvas=canvas;this.config={};this.isReady=false;this.fps=0;this.states={measurements:[],infos:[]};};x3dom.Runtime.prototype.addMeasurement=function(title,value){this.states.measurements[title]=value;};x3dom.Runtime.prototype.removeMeasurement=function(title){if(this.states.measurements[title]){delete this.states.measurements[title];}};x3dom.Runtime.prototype.addInfo=function(title,value){this.states.infos[title]=value;};x3dom.Runtime.prototype.removeInfo=function(title){delete this.states.infos[title];};x3dom.Runtime.prototype.initialize=function(doc,canvas){this.doc=doc;this.canvas=canvas;this.config={};this.isReady=false;this.fps=0;};x3dom.Runtime.prototype.noBackendFound=function(){x3dom.debug.logInfo('No backend found. Unable to render.');};x3dom.Runtime.prototype.ready=function(){x3dom.debug.logInfo('System ready.');};x3dom.Runtime.prototype.enterFrame=function(){};x3dom.Runtime.prototype.exitFrame=function(){};x3dom.Runtime.prototype.triggerRedraw=function(){this.canvas.doc.needRender=true;};x3dom.Runtime.prototype.getActiveBindable=function(typeName){var stacks;var i,current,result;var type;stacks=this.canvas.doc._bindableBag._stacks;result=[];type=x3dom.nodeTypesLC[typeName.toLowerCase()];if(!type){x3dom.debug.logError('No node of type "'+typeName+'" found.');return null;}
for(i=0;i<stacks.length;i++){current=stacks[i].getActive();if(current._xmlNode!==undefined&&x3dom.isa(current,type)){result.push(current);}}
return result[0]?result[0]._xmlNode:null;};x3dom.Runtime.prototype.nextView=function(){var stack=this.canvas.doc._scene.getViewpoint()._stack;if(stack){stack.switchTo('next');}else{x3dom.debug.logError('No valid ViewBindable stack.');}};x3dom.Runtime.prototype.prevView=function(){var stack=this.canvas.doc._scene.getViewpoint()._stack;if(stack){stack.switchTo('prev');}else{x3dom.debug.logError('No valid ViewBindable stack.');}};x3dom.Runtime.prototype.viewpoint=function(){return this.canvas.doc._scene.getViewpoint();};x3dom.Runtime.prototype.viewMatrix=function(){return this.canvas.doc._viewarea.getViewMatrix();};x3dom.Runtime.prototype.projectionMatrix=function(){return this.canvas.doc._viewarea.getProjectionMatrix();};x3dom.Runtime.prototype.getWorldToCameraCoordinatesMatrix=function(){return this.canvas.doc._viewarea.getWCtoCCMatrix();};x3dom.Runtime.prototype.getCameraToWorldCoordinatesMatrix=function(){return this.canvas.doc._viewarea.getCCtoWCMatrix();};x3dom.Runtime.prototype.getViewingRay=function(x,y){return this.canvas.doc._viewarea.calcViewRay(x,y);};x3dom.Runtime.prototype.shootRay=function(x,y){var doc=this.canvas.doc;var info=doc._viewarea._pickingInfo;doc.onPick(this.canvas.gl,x,y);return{pickPosition:info.pickObj?info.pickPos:null,pickNormal:info.pickObj?info.pickNorm:null,pickObject:info.pickObj?info.pickObj._xmlNode:null};};x3dom.Runtime.prototype.getWidth=function(){return this.canvas.doc._viewarea._width;};x3dom.Runtime.prototype.getHeight=function(){return this.canvas.doc._viewarea._height;};x3dom.Runtime.prototype.mousePosition=function(event){var pos=this.canvas.mousePosition(event);return[pos.x,pos.y];};x3dom.Runtime.prototype.calcCanvasPos=function(wx,wy,wz){var pnt=new x3dom.fields.SFVec3f(wx,wy,wz);var mat=this.canvas.doc._viewarea.getWCtoCCMatrix();var pos=mat.multFullMatrixPnt(pnt);var w=this.canvas.doc._viewarea._width;var h=this.canvas.doc._viewarea._height;var x=Math.round((pos.x+1)*(w-1)/2);var y=Math.round((h-1)*(1-pos.y)/2);return[x,y];};x3dom.Runtime.prototype.calcPagePos=function(wx,wy,wz){var elem=this.canvas.canvas.offsetParent;if(!elem){x3dom.debug.logError("Can't calc page pos without offsetParent.");return[0,0];}
var canvasPos=elem.getBoundingClientRect();var mousePos=this.calcCanvasPos(wx,wy,wz);var scrollLeft=window.pageXOffset||document.body.scrollLeft;var scrollTop=window.pageYOffset||document.body.scrollTop;var compStyle=document.defaultView.getComputedStyle(elem,null);var paddingLeft=parseFloat(compStyle.getPropertyValue('padding-left'));var borderLeftWidth=parseFloat(compStyle.getPropertyValue('border-left-width'));var paddingTop=parseFloat(compStyle.getPropertyValue('padding-top'));var borderTopWidth=parseFloat(compStyle.getPropertyValue('border-top-width'));var x=canvasPos.left+paddingLeft+borderLeftWidth+scrollLeft+mousePos[0];var y=canvasPos.top+paddingTop+borderTopWidth+scrollTop+mousePos[1];return[x,y];};x3dom.Runtime.prototype.calcClientPos=function(wx,wy,wz){var elem=this.canvas.canvas.offsetParent;if(!elem){x3dom.debug.logError("Can't calc client pos without offsetParent.");return[0,0];}
var canvasPos=elem.getBoundingClientRect();var mousePos=this.calcCanvasPos(wx,wy,wz);var compStyle=document.defaultView.getComputedStyle(elem,null);var paddingLeft=parseFloat(compStyle.getPropertyValue('padding-left'));var borderLeftWidth=parseFloat(compStyle.getPropertyValue('border-left-width'));var paddingTop=parseFloat(compStyle.getPropertyValue('padding-top'));var borderTopWidth=parseFloat(compStyle.getPropertyValue('border-top-width'));var x=canvasPos.left+paddingLeft+borderLeftWidth+mousePos[0];var y=canvasPos.top+paddingTop+borderTopWidth+mousePos[1];return[x,y];};x3dom.Runtime.prototype.getScreenshot=function(){var url="";var backend=this.canvas.backend;var canvas=this.canvas.canvas;if(canvas){if(backend=="flash"){url=canvas.getScreenshot();}
else{var canvas2d=document.createElement("canvas");canvas2d.width=canvas.width;canvas2d.height=canvas.height;var ctx=canvas2d.getContext("2d");ctx.drawImage(canvas,0,0,canvas.width,canvas.height);ctx.scale(1,-1);ctx.translate(0,-canvas.height);url=canvas2d.toDataURL();}}
return url;};x3dom.Runtime.prototype.getCanvas=function(){return this.canvas.canvas;};x3dom.Runtime.prototype.lightMatrix=function(){this.canvas.doc._viewarea.getLightMatrix();};x3dom.Runtime.prototype.resetView=function(){this.canvas.doc._viewarea.resetView();};x3dom.Runtime.prototype.lightView=function(){if(this.canvas.doc._nodeBag.lights.length>0){this.canvas.doc._viewarea.animateTo(this.canvas.doc._viewarea.getLightMatrix()[0],this.canvas.doc._scene.getViewpoint());return true;}else{x3dom.debug.logInfo("No lights to navigate to.");return false;}};x3dom.Runtime.prototype.uprightView=function(){this.canvas.doc._viewarea.uprightView();};x3dom.Runtime.prototype.fitAll=function(updateCenterOfRotation)
{if(updateCenterOfRotation===undefined){updateCenterOfRotation=true;}
var scene=this.canvas.doc._scene;scene.updateVolume();this.canvas.doc._viewarea.fit(scene._lastMin,scene._lastMax,updateCenterOfRotation);};x3dom.Runtime.prototype.fitObject=function(obj,updateCenterOfRotation)
{if(obj&&obj._x3domNode)
{if(updateCenterOfRotation===undefined){updateCenterOfRotation=true;}
var min=x3dom.fields.SFVec3f.MAX();var max=x3dom.fields.SFVec3f.MIN();var vol=obj._x3domNode.getVolume();vol.getBounds(min,max);var mat=obj._x3domNode.getCurrentTransform();min=mat.multMatrixPnt(min);max=mat.multMatrixPnt(max);if(x3dom.isa(obj._x3domNode,x3dom.nodeTypes.X3DTransformNode))
{var invMat=obj._x3domNode._trafo.inverse();min=invMat.multMatrixPnt(min);max=invMat.multMatrixPnt(max);}
this.canvas.doc._viewarea.fit(min,max,updateCenterOfRotation);}};x3dom.Runtime.prototype.showAll=function(axis,updateCenterOfRotation){this.canvas.doc._viewarea.showAll(axis,updateCenterOfRotation);};x3dom.Runtime.prototype.showObject=function(obj,axis)
{if(obj&&obj._x3domNode)
{var min=x3dom.fields.SFVec3f.MAX();var max=x3dom.fields.SFVec3f.MIN();var vol=obj._x3domNode.getVolume();vol.getBounds(min,max);var mat=obj._x3domNode.getCurrentTransform();min=mat.multMatrixPnt(min);max=mat.multMatrixPnt(max);var viewarea=this.canvas.doc._viewarea;var focalLen=(viewarea._width<viewarea._height)?viewarea._width:viewarea._height;var n0;switch(axis)
{case"posX":n0=new x3dom.fields.SFVec3f(1,0,0);break;case"negX":n0=new x3dom.fields.SFVec3f(-1,0,0);break;case"posY":n0=new x3dom.fields.SFVec3f(0,1,0);break;case"negY":n0=new x3dom.fields.SFVec3f(1,-1,0);break;case"posZ":n0=new x3dom.fields.SFVec3f(0,0,1);break;case"negZ":n0=new x3dom.fields.SFVec3f(0,0,-1);break;}
var viewpoint=this.canvas.doc._scene.getViewpoint();var fov=viewpoint.getFieldOfView()/2.0;var ta=Math.tan(fov);if(Math.abs(ta)>x3dom.fields.Eps){focalLen/=ta;}
var w=viewarea._width-1;var h=viewarea._height-1;var frame=0.25;var minScreenPos=new x3dom.fields.SFVec2f(frame*w,frame*h);frame=0.75;var maxScreenPos=new x3dom.fields.SFVec2f(frame*w,frame*h);var dia2=max.subtract(min).multiply(0.5);var rw=dia2.length();var pc=min.add(dia2);var vc=maxScreenPos.subtract(minScreenPos).multiply(0.5);var rs=1.5*vc.length();vc=vc.add(minScreenPos);var dist=1.0;if(rs>x3dom.fields.Eps){dist=(rw/rs)*Math.sqrt(vc.x*vc.x+vc.y*vc.y+focalLen*focalLen);}
n0=mat.multMatrixVec(n0).normalize();n0=n0.multiply(dist);var p0=pc.add(n0);var qDir=x3dom.fields.Quaternion.rotateFromTo(new x3dom.fields.SFVec3f(0,0,1),n0);var R=qDir.toMatrix();var T=x3dom.fields.SFMatrix4f.translation(p0.negate());var M=x3dom.fields.SFMatrix4f.translation(p0);M=M.mult(R).mult(T).mult(M);var viewmat=M.inverse();viewarea.animateTo(viewmat,viewpoint);}};x3dom.Runtime.prototype.getCenter=function(domNode){if(domNode&&domNode._x3domNode&&(this.isA(domNode,"X3DShapeNode")||this.isA(domNode,"X3DGeometryNode")))
{return domNode._x3domNode.getCenter();}
return null;};x3dom.Runtime.prototype.getCurrentTransform=function(domNode){if(domNode&&domNode._x3domNode)
{return domNode._x3domNode.getCurrentTransform();}
return null;};x3dom.Runtime.prototype.getBBox=function(domNode){if(domNode&&domNode._x3domNode&&this.isA(domNode,"X3DBoundedObject"))
{var vol=domNode._x3domNode.getVolume();return{min:x3dom.fields.SFVec3f.copy(vol.min),max:x3dom.fields.SFVec3f.copy(vol.max)}}
return null;};x3dom.Runtime.prototype.getSceneBBox=function(){var scene=this.canvas.doc._scene;scene.updateVolume();return{min:x3dom.fields.SFVec3f.copy(scene._lastMin),max:x3dom.fields.SFVec3f.copy(scene._lastMax)}};x3dom.Runtime.prototype.debug=function(show){var doc=this.canvas.doc;if(doc._viewarea._visDbgBuf===undefined)
doc._viewarea._visDbgBuf=(doc._x3dElem.getAttribute("showLog")==='true');if(arguments.length>0){if(show===true){doc._viewarea._visDbgBuf=true;x3dom.debug.logContainer.style.display="block";}
else{doc._viewarea._visDbgBuf=false;x3dom.debug.logContainer.style.display="none";}}
else{doc._viewarea._visDbgBuf=!doc._viewarea._visDbgBuf;x3dom.debug.logContainer.style.display=(doc._viewarea._visDbgBuf==true)?"block":"none";}
doc.needRender=true;return doc._viewarea._visDbgBuf;};x3dom.Runtime.prototype.navigationType=function(){return this.canvas.doc._scene.getNavigationInfo().getType();};x3dom.Runtime.prototype.noNav=function(){this.canvas.doc._scene.getNavigationInfo().setType("none");};x3dom.Runtime.prototype.examine=function(){this.canvas.doc._scene.getNavigationInfo().setType("examine");};x3dom.Runtime.prototype.turnTable=function(){this.canvas.doc._scene.getNavigationInfo().setType("turntable");};x3dom.Runtime.prototype.fly=function(){this.canvas.doc._scene.getNavigationInfo().setType("fly");};x3dom.Runtime.prototype.freeFly=function(){this.canvas.doc._scene.getNavigationInfo().setType("freefly");};x3dom.Runtime.prototype.lookAt=function(){this.canvas.doc._scene.getNavigationInfo().setType("lookat");};x3dom.Runtime.prototype.lookAround=function(){this.canvas.doc._scene.getNavigationInfo().setType("lookaround");};x3dom.Runtime.prototype.walk=function(){this.canvas.doc._scene.getNavigationInfo().setType("walk");};x3dom.Runtime.prototype.game=function(){this.canvas.doc._scene.getNavigationInfo().setType("game");};x3dom.Runtime.prototype.helicopter=function(){this.canvas.doc._scene.getNavigationInfo().setType("helicopter");};x3dom.Runtime.prototype.resetExamin=function(){var viewarea=this.canvas.doc._viewarea;viewarea._rotMat=x3dom.fields.SFMatrix4f.identity();viewarea._transMat=x3dom.fields.SFMatrix4f.identity();viewarea._movement=new x3dom.fields.SFVec3f(0,0,0);viewarea._needNavigationMatrixUpdate=true;this.canvas.doc.needRender=true;};x3dom.Runtime.prototype.disableKeys=function(){this.canvas.disableKeys=true;};x3dom.Runtime.prototype.enableKeys=function(){this.canvas.disableKeys=false;};x3dom.Runtime.prototype.disableLeftDrag=function(){this.canvas.disableLeftDrag=true;};x3dom.Runtime.prototype.enableLeftDrag=function(){this.canvas.disableLeftDrag=false;};x3dom.Runtime.prototype.disableRightDrag=function(){this.canvas.disableRightDrag=true;};x3dom.Runtime.prototype.enableRightDrag=function(){this.canvas.disableRightDrag=false;};x3dom.Runtime.prototype.disableMiddleDrag=function(){this.canvas.disableMiddleDrag=true;};x3dom.Runtime.prototype.enableMiddleDrag=function(){this.canvas.disableMiddleDrag=false;};x3dom.Runtime.prototype.togglePoints=function(lines){var doc=this.canvas.doc;var mod=(lines===true)?3:2;doc._viewarea._points=++doc._viewarea._points%mod;doc.needRender=true;return doc._viewarea._points;};x3dom.Runtime.prototype.pickRect=function(x1,y1,x2,y2){return this.canvas.doc.onPickRect(this.canvas.gl,x1,y1,x2,y2);};x3dom.Runtime.prototype.pickMode=function(options){if(options&&options.internal===true){return this.canvas.doc._scene._vf.pickMode;}
return this.canvas.doc._scene._vf.pickMode.toLowerCase();};x3dom.Runtime.prototype.changePickMode=function(type){type=type.toLowerCase();switch(type){case'idbuf':type='idBuf';break;case'idbuf24':type='idBuf24';break;case'idbufid':type='idBufId';break;case'texcoord':type='texCoord';break;case'color':type='color';break;case'box':type='box';break;default:x3dom.debug.logWarning("Switch pickMode to "+type+' unknown intersect type');type=undefined;}
if(type!==undefined){this.canvas.doc._scene._vf.pickMode=type;x3dom.debug.logInfo("Switched pickMode to '"+type+"'.");return true;}
return false;};x3dom.Runtime.prototype.speed=function(newSpeed){var navi=this.canvas.doc._scene.getNavigationInfo();if(newSpeed){navi._vf.speed=newSpeed;x3dom.debug.logInfo("Changed navigation speed to "+navi._vf.speed);}
return navi._vf.speed;};x3dom.Runtime.prototype.zoom=function(zoomAmount){this.canvas.doc._viewarea.zoom(zoomAmount);this.canvas.doc.needRender=true;};x3dom.Runtime.prototype.statistics=function(mode){var states=this.canvas.stateViewer;if(states){this.canvas.doc.needRender=true;if(mode===true){states.display(mode);return true;}
else if(mode===false){states.display(mode);return false;}
else{states.display(!states.active);return states.active;}}
return false;};x3dom.Runtime.prototype.processIndicator=function(mode){var processDiv=this.canvas.progressDiv;if(processDiv){if(mode===true){processDiv.style.display='inline';return true;}
else if(mode===false){processDiv.style.display='none';return false;}
return processDiv.style.display!='none'}
return false;};x3dom.Runtime.prototype.properties=function(){return this.canvas.doc.properties;};x3dom.Runtime.prototype.backendName=function(){return this.canvas.backend;};x3dom.Runtime.prototype.getFPS=function(){return this.fps;};x3dom.Runtime.prototype.isA=function(domNode,nodeType){var inherits=false;if(nodeType&&domNode&&domNode._x3domNode){if(nodeType===""){nodeType="X3DNode";}
inherits=x3dom.isa(domNode._x3domNode,x3dom.nodeTypesLC[nodeType.toLowerCase()]);}
return inherits;};x3dom.Runtime.prototype.getPixelScale=function(){var vp=this.viewpoint();if(!x3dom.isa(vp,x3dom.nodeTypes.OrthoViewpoint)){x3dom.debug.logError("getPixelScale is only implemented for orthographic Viewpoints");return null;}
var zoomLevel=vp.getZoom();var left=zoomLevel[0];var bottom=zoomLevel[1];var right=zoomLevel[2];var top=zoomLevel[3];var x=right-left;var y=top-bottom;var pixelScaleX=x/this.getWidth();var pixelScaleY=y/this.getHeight();return new x3dom.fields.SFVec3f(pixelScaleX,pixelScaleY,0.0);};x3dom.Runtime.prototype.toggleProjection=function(perspViewID,orthoViewID)
{var dist;var factor=2.2;var runtime=document.getElementById("x3d").runtime;var navInfo=runtime.canvas.doc._scene.getNavigationInfo();var speed=navInfo._vf.transitionTime;var persp=document.getElementById(perspViewID)._x3domNode;var ortho=document.getElementById(orthoViewID)._x3domNode;navInfo._vf.transitionTime=0;ortho._bindAnimation=false;persp._bindAnimation=false;if(persp._vf.isActive){ortho._viewMatrix=persp._viewMatrix;document.getElementById(orthoViewID).setAttribute("set_bind","true");dist=persp._viewMatrix.e3().length()/factor;ortho.setZoom(dist);}
else if(ortho._vf.isActive){persp._viewMatrix=ortho._viewMatrix;document.getElementById(perspViewID).setAttribute("set_bind","true");dist=ortho._fieldOfView[2]*factor;var translation=ortho._viewMatrix.e3().normalize().multiply(dist);persp._viewMatrix.setTranslate(translation);}
navInfo._vf.transitionTime=speed;ortho._bindAnimation=true;persp._bindAnimation=true;return(persp._vf.isActive)?0:1;};x3dom.userAgentFeature={supportsDOMAttrModified:false};(function loadX3DOM(){"use strict";var onload=function(){var i,j;var x3ds_unfiltered=document.getElementsByTagName('X3D');var x3ds=[];for(i=0;i<x3ds_unfiltered.length;i++){if(x3ds_unfiltered[i].hasRuntime===undefined)
x3ds.push(x3ds_unfiltered[i]);}
var params;var settings=new x3dom.Properties();var validParams=array_to_object(['showLog','showStat','showProgress','PrimitiveQuality','components','loadpath','disableDoubleClick','backend','altImg','flashrenderer','swfpath','runtimeEnabled','keysEnabled','showTouchpoints','disableTouch','maxActiveDownloads']);var components,prefix;var showLoggingConsole=false;for(i=0;i<x3ds.length;i++){settings.setProperty("showLog",x3ds[i].getAttribute("showLog")||'false');settings.setProperty("showStat",x3ds[i].getAttribute("showStat")||'false');settings.setProperty("showProgress",x3ds[i].getAttribute("showProgress")||'true');settings.setProperty("PrimitiveQuality",x3ds[i].getAttribute("PrimitiveQuality")||'High');params=x3ds[i].getElementsByTagName('PARAM');for(j=0;j<params.length;j++){if(params[j].getAttribute('name')in validParams){settings.setProperty(params[j].getAttribute('name'),params[j].getAttribute('value'));}else{}}
if(settings.getProperty('showLog')==='true'){showLoggingConsole=true;}
if(typeof X3DOM_SECURITY_OFF!='undefined'&&X3DOM_SECURITY_OFF===true){components=settings.getProperty('components',x3ds[i].getAttribute("components"));if(components){prefix=settings.getProperty('loadpath',x3ds[i].getAttribute("loadpath"));components=components.trim().split(',');for(j=0;j<components.length;j++){x3dom.loadJS(components[j]+".js",prefix);}}
if(x3ds[i].getAttribute("src")){var _scene=document.createElement("scene");var _inl=document.createElement("Inline");_inl.setAttribute("url",x3ds[i].getAttribute("src"));_scene.appendChild(_inl);x3ds[i].appendChild(_scene);}}}
if(showLoggingConsole==true){x3dom.debug.activate(true);}else{x3dom.debug.activate(false);}
x3ds=Array.map(x3ds,function(n){n.hasRuntime=true;return n;});if(x3dom.versionInfo!==undefined){x3dom.debug.logInfo("X3DOM version "+x3dom.versionInfo.version+", "+"Revison <a href='https://github.com/x3dom/x3dom/tree/"+x3dom.versionInfo.revision+"'>"
+x3dom.versionInfo.revision+"</a>, "+"Date "+x3dom.versionInfo.date);}
x3dom.debug.logInfo("Found "+x3ds.length+" X3D and nodes...");var x3d_element;var x3dcanvas;var altDiv,altP,aLnk,altImg;var t0,t1;for(i=0;i<x3ds.length;i++)
{x3d_element=x3ds[i];x3dcanvas=new x3dom.X3DCanvas(x3d_element,x3dom.canvases.length);x3dom.canvases.push(x3dcanvas);if(x3dcanvas.gl===null){altDiv=document.createElement("div");altDiv.setAttribute("class","x3dom-nox3d");altDiv.setAttribute("id","x3dom-nox3d");altP=document.createElement("p");altP.appendChild(document.createTextNode("WebGL is not yet supported in your browser. "));aLnk=document.createElement("a");aLnk.setAttribute("href","http://www.x3dom.org/?page_id=9");aLnk.appendChild(document.createTextNode("Follow link for a list of supported browsers... "));altDiv.appendChild(altP);altDiv.appendChild(aLnk);x3dcanvas.x3dElem.appendChild(altDiv);if(x3dcanvas.stateViewer){x3d_element.removeChild(x3dcanvas.stateViewer.viewer);}
continue;}
t0=new Date().getTime();x3ds[i].runtime=new x3dom.Runtime(x3ds[i],x3dcanvas);x3ds[i].runtime.initialize(x3ds[i],x3dcanvas);if(x3dom.runtime.ready){x3ds[i].runtime.ready=x3dom.runtime.ready;}
if(x3dcanvas.backend==''){x3dom.runtime.noBackendFound();}
x3dcanvas.load(x3ds[i],i,settings);if(settings.getProperty('showStat')==='true'){x3ds[i].runtime.statistics(true);}else{x3ds[i].runtime.statistics(false);}
if(settings.getProperty('showProgress')==='true'){if(settings.getProperty('showProgress')==='bar'){x3dcanvas.progressDiv.setAttribute("class","x3dom-progress bar");}
x3ds[i].runtime.processIndicator(true);}else{x3ds[i].runtime.processIndicator(false);}
t1=new Date().getTime()-t0;x3dom.debug.logInfo("Time for setup and init of GL element no. "+i+": "+t1+" ms.");}
var ready=(function(eventType){var evt=null;if(document.createEvent){evt=document.createEvent("Events");evt.initEvent(eventType,true,true);document.dispatchEvent(evt);}else if(document.createEventObject){evt=document.createEventObject();document.body.fireEvent('on'+eventType,evt);}})('load');};var onunload=function(){if(x3dom.canvases){for(var i=0;i<x3dom.canvases.length;i++){x3dom.canvases[i].doc.shutdown(x3dom.canvases[i].gl);}
x3dom.canvases=[];}};x3dom.reload=function(){onload();};if(window.addEventListener){window.addEventListener('load',onload,false);window.addEventListener('unload',onunload,false);window.addEventListener('reload',onunload,false);}else if(window.attachEvent){window.attachEvent('onload',onload);window.attachEvent('onunload',onunload);window.attachEvent('onreload',onunload);}
if(document.readyState==="complete"){window.setTimeout(function(){onload();},20);}})();x3dom.Cache=function(){this.textures=[];this.shaders=[];};x3dom.Cache.prototype.getTexture2D=function(gl,doc,url,bgnd,crossOrigin,scale,genMipMaps){var textureIdentifier=url;if(this.textures[textureIdentifier]===undefined){this.textures[textureIdentifier]=x3dom.Utils.createTexture2D(gl,doc,url,bgnd,crossOrigin,scale,genMipMaps);}
return this.textures[textureIdentifier];};x3dom.Cache.prototype.getTexture2DByDEF=function(gl,nameSpace,def){var textureIdentifier=nameSpace.name+"_"+def;if(this.textures[textureIdentifier]===undefined){this.textures[textureIdentifier]=gl.createTexture();}
return this.textures[textureIdentifier];};x3dom.Cache.prototype.getTextureCube=function(gl,doc,url,bgnd,crossOrigin,scale,genMipMaps){var textureIdentifier="";for(var i=0;i<url.length;++i){textureIdentifier+=url[i]+"|";}
if(this.textures[textureIdentifier]===undefined){this.textures[textureIdentifier]=x3dom.Utils.createTextureCube(gl,doc,url,bgnd,crossOrigin,scale,genMipMaps);}
return this.textures[textureIdentifier];};x3dom.Cache.prototype.getShader=function(gl,shaderIdentifier){var program=null;if(this.shaders[shaderIdentifier]===undefined){switch(shaderIdentifier){case x3dom.shader.PICKING:program=new x3dom.shader.PickingShader(gl);break;case x3dom.shader.PICKING_24:program=new x3dom.shader.Picking24Shader(gl);break;case x3dom.shader.PICKING_ID:program=new x3dom.shader.PickingIdShader(gl);break;case x3dom.shader.PICKING_COLOR:program=new x3dom.shader.PickingColorShader(gl);break;case x3dom.shader.PICKING_TEXCOORD:program=new x3dom.shader.PickingTexcoordShader(gl);break;case x3dom.shader.FRONTGROUND_TEXTURE:program=new x3dom.shader.FrontgroundTextureShader(gl);break;case x3dom.shader.BACKGROUND_TEXTURE:program=new x3dom.shader.BackgroundTextureShader(gl);break;case x3dom.shader.BACKGROUND_SKYTEXTURE:program=new x3dom.shader.BackgroundSkyTextureShader(gl);break;case x3dom.shader.BACKGROUND_CUBETEXTURE:program=new x3dom.shader.BackgroundCubeTextureShader(gl);break;case x3dom.shader.SHADOW:program=new x3dom.shader.ShadowShader(gl);break;case x3dom.shader.BLUR:program=new x3dom.shader.BlurShader(gl);break;case x3dom.shader.DEPTH:break;case x3dom.shader.NORMAL:program=new x3dom.shader.NormalShader(gl);break;case x3dom.shader.TEXTURE_REFINEMENT:program=new x3dom.shader.TextureRefinementShader(gl);break;default:break;}
if(program)
this.shaders[shaderIdentifier]=x3dom.Utils.wrapProgram(gl,program,shaderIdentifier);else
x3dom.debug.logError("Couldn't create shader: "+shaderIdentifier);}
return this.shaders[shaderIdentifier];};x3dom.Cache.prototype.getDynamicShader=function(gl,viewarea,shape){var properties=x3dom.Utils.generateProperties(viewarea,shape);var shaderID=properties.id;if(this.shaders[shaderID]===undefined){var program=null;if(properties.CSHADER!=-1){program=new x3dom.shader.ComposedShader(gl,shape);}else{program=(x3dom.caps.MOBILE&&!properties.CSSHADER)?new x3dom.shader.DynamicMobileShader(gl,properties):new x3dom.shader.DynamicShader(gl,properties);}
this.shaders[shaderID]=x3dom.Utils.wrapProgram(gl,program,shaderID);}
return this.shaders[shaderID];};x3dom.Cache.prototype.getShaderByProperties=function(gl,shape,properties,pickMode,shadows){var shaderID=properties.id;if(pickMode!==undefined&&pickMode!==null){shaderID+=pickMode;}
if(shadows!==undefined&&shadows!==null){shaderID+="S";}
if(this.shaders[shaderID]===undefined)
{var program=null;if(pickMode!==undefined&&pickMode!==null){program=new x3dom.shader.DynamicShaderPicking(gl,properties,pickMode);}
else if(shadows!==undefined&&shadows!==null){program=new x3dom.shader.DynamicShadowShader(gl,properties);}
else if(properties.CSHADER!=-1)
program=new x3dom.shader.ComposedShader(gl,shape);else if(properties.KHR_MATERIAL_COMMONS!=null&&properties.KHR_MATERIAL_COMMONS!=0)
program=new x3dom.shader.KHRMaterialCommonsShader(gl,properties);else if(properties.EMPTY_SHADER!=null&&properties.EMPTY_SHADER!=0)
return{"shaderID":shaderID};else{program=(x3dom.caps.MOBILE&&!properties.CSSHADER)?new x3dom.shader.DynamicMobileShader(gl,properties):new x3dom.shader.DynamicShader(gl,properties);}
this.shaders[shaderID]=x3dom.Utils.wrapProgram(gl,program,shaderID);}
return this.shaders[shaderID];};x3dom.Cache.prototype.getShadowRenderingShader=function(gl,shadowedLights){var ID="shadow";for(var i=0;i<shadowedLights.length;i++){if(x3dom.isa(shadowedLights[i],x3dom.nodeTypes.SpotLight))
ID+="S";else if(x3dom.isa(shadowedLights[i],x3dom.nodeTypes.PointLight))
ID+="P";else
ID+="D";}
if(this.shaders[ID]===undefined){var program=new x3dom.shader.ShadowRenderingShader(gl,shadowedLights);this.shaders[ID]=x3dom.Utils.wrapProgram(gl,program,ID);}
return this.shaders[ID];};x3dom.Cache.prototype.Release=function(gl){for(var texture in this.textures){gl.deleteTexture(this.textures[texture]);}
this.textures=[];for(var shaderId in this.shaders){var shader=this.shaders[shaderId];var glShaders=gl.getAttachedShaders(shader.program);for(var i=0;i<glShaders.length;++i){gl.detachShader(shader.program,glShaders[i]);gl.deleteShader(glShaders[i]);}
gl.deleteProgram(shader.program)}
this.shaders=[];};function startDashVideo(recurl,texturediv){var vars=function(){var vars={};var parts=window.location.href.replace(/[?&]+([^=&]+)=([^&]*)/gi,function(m,key,value){vars[key]=value;});return vars;},url=recurl,video,context,player;if(vars&&vars.hasOwnProperty("url")){url=vars.url;}
video=document.querySelector(texturediv);context=new Dash.di.DashContext();player=new MediaPlayer(context);player.startup();player.attachView(video);player.setAutoPlay(false);player.attachSource(url);}
x3dom.Texture=function(gl,doc,cache,node){this.gl=gl;this.doc=doc;this.cache=cache;this.node=node;this.samplerName="diffuseMap";this.type=gl.TEXTURE_2D;this.format=gl.RGBA;this.magFilter=gl.LINEAR;this.minFilter=gl.LINEAR;this.wrapS=gl.REPEAT;this.wrapT=gl.REPEAT;this.genMipMaps=false;this.texture=null;this.ready=false;this.dashtexture=false;var tex=this.node;var suffix="mpd";this.node._x3domTexture=this;if(x3dom.isa(tex,x3dom.nodeTypes.MovieTexture)){if(tex._vf.url[0].indexOf(suffix,tex._vf.url[0].length-suffix.length)!==-1){this.dashtexture=true;var js=document.getElementById("AdditionalDashVideoScript");if(!js){js=document.createElement("script");js.setAttribute("type","text/javascript");js.setAttribute("src",x3dom.Texture.dashVideoScriptFile);js.setAttribute("id","AdditionalDashVideoScript");js.onload=function(){var texObj;while((texObj=x3dom.Texture.loadDashVideos.pop())){x3dom.Texture.textNum++;texObj.update();}
js.ready=true;};document.getElementsByTagName('head')[0].appendChild(js);}
if(js.ready===true){x3dom.Texture.textNum++;this.update();}
else{x3dom.Texture.loadDashVideos.push(this);}}}
if(!this.dashtexture){this.update();}};x3dom.Texture.dashVideoScriptFile="dash.all.js";x3dom.Texture.loadDashVideos=[];x3dom.Texture.textNum=0;x3dom.Texture.clampFontSize=false;x3dom.Texture.minFontQuality=0.5;x3dom.Texture.maxFontQuality=10;x3dom.Texture.prototype.update=function()
{if(x3dom.isa(this.node,x3dom.nodeTypes.Text))
{this.updateText();}
else
{this.updateTexture();}};x3dom.Texture.prototype.setPixel=function(x,y,pixel,update)
{var gl=this.gl;var pixels=new Uint8Array(pixel);gl.bindTexture(this.type,this.texture);gl.pixelStorei(gl.UNPACK_ALIGNMENT,1);gl.texSubImage2D(this.type,0,x,y,1,1,this.format,gl.UNSIGNED_BYTE,pixels);gl.bindTexture(this.type,null);if(update){this.doc.needRender=true;}};x3dom.Texture.prototype.updateTexture=function()
{var gl=this.gl;var doc=this.doc;var tex=this.node;this.samplerName=tex._type;if(x3dom.isa(tex,x3dom.nodeTypes.X3DEnvironmentTextureNode)){this.type=gl.TEXTURE_CUBE_MAP;}else{this.type=gl.TEXTURE_2D;}
if(x3dom.isa(tex,x3dom.nodeTypes.PixelTexture)){switch(tex._vf.image.comp)
{case 1:this.format=gl.LUMINANCE;break;case 2:this.format=gl.LUMINANCE_ALPHA;break;case 3:this.format=gl.RGB;break;case 4:this.format=gl.RGBA;break;}}else{this.format=gl.RGBA;}
if(tex._cf.textureProperties.node!==null){var texProp=tex._cf.textureProperties.node;this.wrapS=x3dom.Utils.boundaryModesDic(gl,texProp._vf.boundaryModeS);this.wrapT=x3dom.Utils.boundaryModesDic(gl,texProp._vf.boundaryModeT);this.minFilter=x3dom.Utils.minFilterDic(gl,texProp._vf.minificationFilter);this.magFilter=x3dom.Utils.magFilterDic(gl,texProp._vf.magnificationFilter);if(texProp._vf.generateMipMaps===true){this.genMipMaps=true;if(this.minFilter==gl.NEAREST){this.minFilter=gl.NEAREST_MIPMAP_NEAREST;}else if(this.minFilter==gl.LINEAR){this.minFilter=gl.LINEAR_MIPMAP_LINEAR;}
if(this.texture&&(this.texture.ready||this.texture.textureCubeReady)){gl.bindTexture(this.type,this.texture);gl.generateMipmap(this.type);gl.bindTexture(this.type,null);}}else{this.genMipMaps=false;if((this.minFilter==gl.LINEAR_MIPMAP_LINEAR)||(this.minFilter==gl.LINEAR_MIPMAP_NEAREST)){this.minFilter=gl.LINEAR;}else if((this.minFilter==gl.NEAREST_MIPMAP_LINEAR)||(this.minFilter==gl.NEAREST_MIPMAP_NEAREST)){this.minFilter=gl.NEAREST;}}}else{if(tex._vf.repeatS==false){this.wrapS=gl.CLAMP_TO_EDGE;}
else
{this.wrapS=gl.REPEAT;}
if(tex._vf.repeatT==false){this.wrapT=gl.CLAMP_TO_EDGE;}
else
{this.wrapT=gl.REPEAT;}
if(this.samplerName=="displacementMap"||this.samplerName=="multiDiffuseAlphaMap"||this.samplerName=="multiVisibilityMap"||this.samplerName=="multiEmissiveAmbientMap"||this.samplerName=="multiSpecularShininessMap")
{this.wrapS=gl.CLAMP_TO_EDGE;this.wrapT=gl.CLAMP_TO_EDGE;this.minFilter=gl.NEAREST;this.magFilter=gl.NEAREST;}}
var childTex=(tex._video&&tex._needPerFrameUpdate===true);if(tex._isCanvas&&tex._canvas)
{if(this.texture==null){this.texture=gl.createTexture()}
this.texture.width=tex._canvas.width;this.texture.height=tex._canvas.height;this.texture.ready=true;gl.bindTexture(this.type,this.texture);gl.texImage2D(this.type,0,this.format,this.format,gl.UNSIGNED_BYTE,tex._canvas);if(this.genMipMaps){gl.generateMipmap(this.type);}
gl.bindTexture(this.type,null);}
else if(x3dom.isa(tex,x3dom.nodeTypes.RenderedTexture))
{if(tex._webgl&&tex._webgl.fbo){if(tex._webgl.fbo.dtex&&tex._vf.depthMap)
this.texture=tex._webgl.fbo.dtex;else
this.texture=tex._webgl.fbo.tex;}
else{this.texture=null;x3dom.debug.logError("Try updating RenderedTexture without FBO initialized!");}
if(this.texture){this.texture.ready=true;}}
else if(x3dom.isa(tex,x3dom.nodeTypes.PixelTexture))
{if(this.texture==null){if(this.node._DEF){this.texture=this.cache.getTexture2DByDEF(gl,this.node._nameSpace,this.node._DEF);}else{this.texture=gl.createTexture();}}
this.texture.width=tex._vf.image.width;this.texture.height=tex._vf.image.height;this.texture.ready=true;var pixelArr=tex._vf.image.array;var pixelArrfont_size=tex._vf.image.width*tex._vf.image.height*tex._vf.image.comp;if(pixelArr.length<pixelArrfont_size)
{pixelArr=tex._vf.image.toGL();while(pixelArr.length<pixelArrfont_size){pixelArr.push(0);}}
var pixels=new Uint8Array(pixelArr);gl.bindTexture(this.type,this.texture);gl.pixelStorei(gl.UNPACK_ALIGNMENT,1);gl.texImage2D(this.type,0,this.format,tex._vf.image.width,tex._vf.image.height,0,this.format,gl.UNSIGNED_BYTE,pixels);if(this.genMipMaps){gl.generateMipmap(this.type);}
gl.bindTexture(this.type,null);}
else if(x3dom.isa(tex,x3dom.nodeTypes.MovieTexture)||childTex)
{var that=this;var p=document.getElementsByTagName('body')[0];if(this.texture==null){this.texture=gl.createTexture();}
if(this.dashtexture){var element_vid=document.createElement('div');element_vid.setAttribute('class','dash-video-player'+x3dom.Texture.textNum);tex._video=document.createElement('video');tex._video.setAttribute('preload','auto');tex._video.setAttribute('muted','muted');var scriptToRun=document.createElement('script');scriptToRun.setAttribute('type','text/javascript');scriptToRun.innerHTML='startDashVideo("'+tex._vf.url[0]+'",".dash-video-player'+x3dom.Texture.textNum+' video")';element_vid.appendChild(scriptToRun);element_vid.appendChild(tex._video);p.appendChild(element_vid);tex._video.style.visibility="hidden";tex._video.style.display="none";}
else{if(!childTex){tex._video=document.createElement('video');tex._video.setAttribute('preload','auto');tex._video.setAttribute('muted','muted');p.appendChild(tex._video);tex._video.style.visibility="hidden";tex._video.style.display="none";}
for(var i=0;i<tex._vf.url.length;i++){var videoUrl=tex._nameSpace.getURL(tex._vf.url[i]);x3dom.debug.logInfo('Adding video file: '+videoUrl);var src=document.createElement('source');src.setAttribute('src',videoUrl);tex._video.appendChild(src);}}
var updateMovie=function()
{gl.bindTexture(that.type,that.texture);gl.texImage2D(that.type,0,that.format,that.format,gl.UNSIGNED_BYTE,tex._video);if(that.genMipMaps){gl.generateMipmap(that.type);}
gl.bindTexture(that.type,null);that.texture.ready=true;doc.needRender=true;};var startVideo=function()
{tex._video.play();tex._intervalID=setInterval(updateMovie,16);};var videoDone=function()
{clearInterval(tex._intervalID);if(tex._vf.loop===true)
{tex._video.play();tex._intervalID=setInterval(updateMovie,16);}};tex._video.addEventListener("canplaythrough",startVideo,true);tex._video.addEventListener("ended",videoDone,true);}
else if(x3dom.isa(tex,x3dom.nodeTypes.X3DEnvironmentTextureNode))
{this.texture=this.cache.getTextureCube(gl,doc,tex.getTexUrl(),false,tex._vf.crossOrigin,tex._vf.scale,this.genMipMaps);}
else
{this.texture=this.cache.getTexture2D(gl,doc,tex._nameSpace.getURL(tex._vf.url[0]),false,tex._vf.crossOrigin,tex._vf.scale,this.genMipMaps);}};x3dom.Texture.prototype.updateText=function()
{var gl=this.gl;this.wrapS=gl.CLAMP_TO_EDGE;this.wrapT=gl.CLAMP_TO_EDGE;this.type=gl.TEXTURE_2D;this.format=gl.RGBA;this.magFilter=gl.LINEAR;this.minFilter=gl.LINEAR;var fontStyleNode=this.node._cf.fontStyle.node;var font_family='serif';var font_style='normal';var font_justify='left';var font_size=1.0;var font_spacing=1.0;var font_horizontal=true;var font_language="";var oversample=2.0;var minor_alignment='FIRST';if(fontStyleNode!==null)
{var fonts=fontStyleNode._vf.family.toString();fonts=fonts.trim().replace(/\'/g,'').replace(/\,/,' ');fonts=fonts.split(" ");font_family=Array.map(fonts,function(s){if(s=='SANS'){return'sans-serif';}
else if(s=='SERIF'){return'serif';}
else if(s=='TYPEWRITER'){return'monospace';}
else{return''+s+'';}}).join(",");font_style=fontStyleNode._vf.style.toString().replace(/\'/g,'');switch(font_style.toUpperCase()){case'PLAIN':font_style='normal';break;case'BOLD':font_style='bold';break;case'ITALIC':font_style='italic';break;case'BOLDITALIC':font_style='italic bold';break;default:font_style='normal';}
var leftToRight=fontStyleNode._vf.leftToRight?'ltr':'rtl';var topToBottom=fontStyleNode._vf.topToBottom;font_justify=fontStyleNode._vf.justify[0].toString().replace(/\'/g,'');switch(font_justify.toUpperCase()){case'BEGIN':font_justify='left';break;case'END':font_justify='right';break;case'FIRST':font_justify='left';break;case'MIDDLE':font_justify='center';break;default:font_justify='left';break;}
if(fontStyleNode._vf.justify[1]===undefined){minor_alignment='FIRST';}
else{minor_alignment=fontStyleNode._vf.justify[1].toString().replace(/\'/g,'');switch(minor_alignment.toUpperCase()){case'BEGIN':minor_alignment='BEGIN';break;case'FIRST':minor_alignment='FIRST';break;case'MIDDLE':minor_alignment='MIDDLE';break;case'END':minor_alignment='END';break;default:minor_alignment='FIRST';break;}}
font_size=fontStyleNode._vf.size;font_spacing=fontStyleNode._vf.spacing;font_horizontal=fontStyleNode._vf.horizontal;font_language=fontStyleNode._vf.language;oversample=fontStyleNode._vf.quality;oversample=Math.max(x3dom.Texture.minFontQuality,oversample);oversample=Math.min(x3dom.Texture.maxFontQuality,oversample);if(font_size<0.1)font_size=0.1;if(x3dom.Texture.clampFontSize&&font_size>2.3)
{font_size=2.3;}}
var textX,textY;var paragraph=this.node._vf.string;var maxExtent=this.node._vf.maxExtent;var lengths=[];var text_canvas=document.createElement('canvas');text_canvas.dir=leftToRight;var x3dToPx=42;var textHeight=font_size*x3dToPx;var textAlignment=font_justify;document.body.appendChild(text_canvas);var text_ctx=text_canvas.getContext('2d');text_ctx.font=font_style+" "+textHeight+"px "+font_family;var maxWidth=0,pWidth,pLength;var i,j;for(i=0;i<paragraph.length;i++){pWidth=text_ctx.measureText(paragraph[i]).width;if(pWidth>maxWidth){maxWidth=pWidth;}
pLength=this.node._vf.length[i]|0;if(maxExtent>0&&(pLength>maxExtent||pLength==0)){pLength=maxExtent;}
lengths[i]=pLength<=0?pWidth:pLength*x3dToPx;}
var canvas_extra=0.1*textHeight;var txtW=maxWidth;var txtH=textHeight*font_spacing*paragraph.length+canvas_extra;textX=0;textY=0;var x_offset=0,y_offset=0,baseLine='top';switch(font_justify){case"center":x_offset=-txtW/2;textX=txtW/2;break;case"left":x_offset=leftToRight=='ltr'?0:-txtW;textX=0;break;case"right":x_offset=leftToRight=='ltr'?-txtW:0;textX=txtW;break;}
switch(minor_alignment){case"MIDDLE":y_offset=txtH/2;break;case"BEGIN":y_offset=topToBottom?0:txtH-canvas_extra;baseLine=topToBottom?'top':'bottom';textY=topToBottom?0:textHeight;break;case"FIRST":y_offset=topToBottom?textHeight:txtH-canvas_extra;baseLine=topToBottom?'alphabetic':'bottom';textY=topToBottom?textHeight:textHeight;break;case"END":y_offset=topToBottom?txtH-canvas_extra:0;baseLine=topToBottom?'bottom':'top';textY=topToBottom?textHeight:0;break;}
var pxToX3d=1/42.0;var w=txtW*pxToX3d;var h=txtH*pxToX3d;x_offset*=pxToX3d;y_offset*=pxToX3d;text_canvas.width=txtW*oversample;text_canvas.height=txtH*oversample;text_canvas.dir=leftToRight;text_ctx.scale(oversample,oversample);text_ctx.fillStyle='rgba(0,0,0,0)';text_ctx.fillRect(0,0,text_ctx.canvas.width,text_ctx.canvas.height);text_ctx.fillStyle='white';text_ctx.textBaseline=baseLine;text_ctx.font=font_style+" "+textHeight+"px "+font_family;text_ctx.textAlign=textAlignment;for(i=0;i<paragraph.length;i++){j=topToBottom?i:paragraph.length-1-i;text_ctx.fillText(paragraph[j],textX,textY,lengths[j]);textY+=textHeight*font_spacing;}
if(this.texture===null)
{this.texture=gl.createTexture();}
gl.bindTexture(this.type,this.texture);gl.texImage2D(this.type,0,this.format,this.format,gl.UNSIGNED_BYTE,text_canvas);gl.bindTexture(this.type,null);document.body.removeChild(text_canvas);this.node._mesh._positions[0]=[0+x_offset,-h+y_offset,0,w+x_offset,-h+y_offset,0,w+x_offset,0+y_offset,0,0+x_offset,0+y_offset,0];this.node.invalidateVolume();Array.forEach(this.node._parentNodes,function(node){node.setAllDirty();});};x3dom.X3DDocument=function(canvas,ctx,settings){this.canvas=canvas;this.ctx=ctx;this.properties=settings;this.needRender=true;this._x3dElem=null;this._scene=null;this._viewarea=null;this.downloadCount=0;this.previousDownloadCount=0;this._nodeBag={timer:[],lights:[],clipPlanes:[],followers:[],trans:[],renderTextures:[],viewarea:[],affectedPointingSensors:[]};this.onload=function(){};this.onerror=function(){};};x3dom.X3DDocument.prototype.load=function(uri,sceneElemPos){var uri_docs={};var queued_uris=[uri];var doc=this;function next_step(){if(queued_uris.length===0){doc._setup(uri_docs[uri],uri_docs,sceneElemPos);doc.onload();return;}
var next_uri=queued_uris.shift();if(x3dom.isX3DElement(next_uri)&&(next_uri.localName.toLowerCase()==='x3d'||next_uri.localName.toLowerCase()==='websg'))
{uri_docs[next_uri]=next_uri;doc._x3dElem=next_uri;next_step();}}
next_step();};x3dom.findScene=function(x3dElem){var sceneElems=[];for(var i=0;i<x3dElem.childNodes.length;i++){var sceneElem=x3dElem.childNodes[i];if(sceneElem&&sceneElem.localName&&sceneElem.localName.toLowerCase()==="scene"){sceneElems.push(sceneElem);}}
if(sceneElems.length>1){x3dom.debug.logError("X3D element has more than one Scene child (has "+
x3dElem.childNodes.length+").");}
else{return sceneElems[0];}
return null;};x3dom.X3DDocument.prototype._setup=function(sceneDoc,uriDocs,sceneElemPos){var doc=this;function cleanNodeBag(bag,node){for(var i=0,n=bag.length;i<n;i++){if(bag[i]===node){bag.splice(i,1);break;}}}
function removeX3DOMBackendGraph(domNode){var children=domNode.childNodes;for(var i=0,n=children.length;i<n;i++){removeX3DOMBackendGraph(children[i]);}
if(domNode._x3domNode){var node=domNode._x3domNode;var nameSpace=node._nameSpace;if(x3dom.isa(node,x3dom.nodeTypes.X3DShapeNode)){if(node._cleanupGLObjects){node._cleanupGLObjects(true);}
if(x3dom.nodeTypes.Shape.idMap.nodeID[node._objectID]){delete x3dom.nodeTypes.Shape.idMap.nodeID[node._objectID];}}
else if(x3dom.isa(node,x3dom.nodeTypes.TimeSensor)){cleanNodeBag(doc._nodeBag.timer,node);}
else if(x3dom.isa(node,x3dom.nodeTypes.X3DLightNode)){cleanNodeBag(doc._nodeBag.lights,node);}
else if(x3dom.isa(node,x3dom.nodeTypes.X3DFollowerNode)){cleanNodeBag(doc._nodeBag.followers,node);}
else if(x3dom.isa(node,x3dom.nodeTypes.X3DTransformNode)){cleanNodeBag(doc._nodeBag.trans,node);}
else if(x3dom.isa(node,x3dom.nodeTypes.RenderedTexture)){cleanNodeBag(doc._nodeBag.renderTextures,node);if(node._cleanupGLObjects){node._cleanupGLObjects();}}
else if(x3dom.isa(node,x3dom.nodeTypes.X3DPointingDeviceSensorNode)){cleanNodeBag(doc._nodeBag.affectedPointingSensors,node);}
else if(x3dom.isa(node,x3dom.nodeTypes.Texture)){node.shutdown();}
else if(x3dom.isa(node,x3dom.nodeTypes.AudioClip)){node.shutdown();}
else if(x3dom.isa(node,x3dom.nodeTypes.X3DBindableNode)){var stack=node._stack;if(stack){node.bind(false);cleanNodeBag(stack._bindBag,node);}
if(node._cleanupGLObjects){node._cleanupGLObjects();}}
else if(x3dom.isa(node,x3dom.nodeTypes.Scene)){if(node._webgl){node._webgl=null;}}
if(nameSpace&&!(domNode.getAttribute('use')||domNode.getAttribute('USE')))
{nameSpace.removeNode(node._DEF);}
node._xmlNode=null;delete domNode._x3domNode;}}
var domEventListener={onAttrModified:function(e){if('_x3domNode'in e.target){var attrToString={1:"MODIFICATION",2:"ADDITION",3:"REMOVAL"};e.target._x3domNode.updateField(e.attrName,e.newValue);doc.needRender=true;}},onNodeRemoved:function(e){var domNode=e.target;if(!domNode)
return;if('_x3domNode'in domNode.parentNode&&'_x3domNode'in domNode){var parent=domNode.parentNode._x3domNode;var child=domNode._x3domNode;if(parent&&child){parent.removeChild(child);parent.nodeChanged();removeX3DOMBackendGraph(domNode);if(doc._viewarea&&doc._viewarea._scene){doc._viewarea._scene.nodeChanged();doc._viewarea._scene.updateVolume();doc.needRender=true;}}}
else if(domNode.localName&&domNode.localName.toUpperCase()=="ROUTE"&&domNode._nodeNameSpace){var fromNode=domNode._nodeNameSpace.defMap[domNode.getAttribute('fromNode')];var toNode=domNode._nodeNameSpace.defMap[domNode.getAttribute('toNode')];if(fromNode&&toNode){fromNode.removeRoute(domNode.getAttribute('fromField'),toNode,domNode.getAttribute('toField'));}}
else if(domNode.localName&&domNode.localName.toUpperCase()=="X3D"){var runtime=domNode.runtime;if(runtime&&runtime.canvas&&runtime.canvas.doc&&runtime.canvas.doc._scene){var sceneNode=runtime.canvas.doc._scene._xmlNode;removeX3DOMBackendGraph(sceneNode);for(var i=0;i<x3dom.canvases.length;i++){if(x3dom.canvases[i]===runtime.canvas){x3dom.canvases[i].doc.shutdown(x3dom.canvases[i].gl);x3dom.canvases.splice(i,1);break;}}
runtime.canvas.doc._scene=null;runtime.canvas.doc._viewarea=null;runtime.canvas.doc=null;runtime.canvas=null;runtime=null;domNode.context=null;domNode.runtime=null;}}},onNodeInserted:function(e){var child=e.target;var parentNode=child.parentNode;if('_x3domNode'in parentNode){if(parentNode.tagName&&parentNode.tagName.toLowerCase()=='inline'||parentNode.tagName.toLowerCase()=='multipart'){}
else{var parent=parentNode._x3domNode;if(parent&&parent._nameSpace&&(child instanceof Element)){if(x3dom.caps.DOMNodeInsertedEvent_perSubtree)
{removeX3DOMBackendGraph(child);}
var newNode=parent._nameSpace.setupTree(child);parent.addChild(newNode,child.getAttribute("containerField"));parent.nodeChanged();var grandParentNode=parentNode.parentNode;if(grandParentNode&&grandParentNode._x3domNode)
grandParentNode._x3domNode.nodeChanged();if(doc._viewarea&&doc._viewarea._scene){doc._viewarea._scene.nodeChanged();doc._viewarea._scene.updateVolume();doc.needRender=true;}}
else{x3dom.debug.logWarning("No _nameSpace in onNodeInserted");}}}}};sceneDoc.addEventListener('DOMNodeRemoved',domEventListener.onNodeRemoved,true);sceneDoc.addEventListener('DOMNodeInserted',domEventListener.onNodeInserted,true);if((x3dom.userAgentFeature.supportsDOMAttrModified===true)){sceneDoc.addEventListener('DOMAttrModified',domEventListener.onAttrModified,true);}
var sceneElem=x3dom.findScene(sceneDoc);this._bindableBag=new x3dom.BindableBag(this);var nameSpace=new x3dom.NodeNameSpace("scene",doc);var scene=nameSpace.setupTree(sceneElem);this._scene=scene;this._bindableBag.setRefNode(scene);this._viewarea=new x3dom.Viewarea(this,scene);this._viewarea._width=this.canvas.width;this._viewarea._height=this.canvas.height;};x3dom.X3DDocument.prototype.advanceTime=function(t){var i=0;if(this._nodeBag.timer.length){for(i=0;i<this._nodeBag.timer.length;i++)
{this.needRender|=this._nodeBag.timer[i].tick(t);}}
if(this._nodeBag.followers.length){for(i=0;i<this._nodeBag.followers.length;i++)
{this.needRender|=this._nodeBag.followers[i].tick(t);}}
if(this._nodeBag.trans.length){for(i=0;i<this._nodeBag.trans.length;i++)
{this.needRender|=this._nodeBag.trans[i].tick(t);}}
if(this._nodeBag.viewarea.length){for(i=0;i<this._nodeBag.viewarea.length;i++)
{this.needRender|=this._nodeBag.viewarea[i].tick(t);}}};x3dom.X3DDocument.prototype.render=function(ctx){if(!ctx||!this._viewarea){return;}
ctx.renderScene(this._viewarea);};x3dom.X3DDocument.prototype.onPick=function(ctx,x,y){if(!ctx||!this._viewarea){return;}
ctx.pickValue(this._viewarea,x,y,1);};x3dom.X3DDocument.prototype.onPickRect=function(ctx,x1,y1,x2,y2){if(!ctx||!this._viewarea){return[];}
return ctx.pickRect(this._viewarea,x1,y1,x2,y2);};x3dom.X3DDocument.prototype.onMove=function(ctx,x,y,buttonState){if(!ctx||!this._viewarea){return;}
if(this._viewarea._scene._vf.doPickPass)
ctx.pickValue(this._viewarea,x,y,buttonState);this._viewarea.onMove(x,y,buttonState);};x3dom.X3DDocument.prototype.onMoveView=function(ctx,evt,touches,translation,rotation){if(!ctx||!this._viewarea){return;}
this._scene.getNavigationInfo()._impl.onTouchDrag(this._viewarea,evt,touches,translation,rotation);};x3dom.X3DDocument.prototype.onDrag=function(ctx,x,y,buttonState){if(!ctx||!this._viewarea){return;}
if(this._viewarea._scene._vf.doPickPass)
ctx.pickValue(this._viewarea,x,y,buttonState);this._viewarea.onDrag(x,y,buttonState);};x3dom.X3DDocument.prototype.onWheel=function(ctx,x,y,originalY){if(!ctx||!this._viewarea){return;}
if(this._viewarea._scene._vf.doPickPass)
ctx.pickValue(this._viewarea,x,originalY,0);this._viewarea.onDrag(x,y,2);};x3dom.X3DDocument.prototype.onMousePress=function(ctx,x,y,buttonState){if(!ctx||!this._viewarea){return;}
this._viewarea._scene.updateVolume();ctx.pickValue(this._viewarea,x,y,buttonState);this._viewarea.onMousePress(x,y,buttonState);};x3dom.X3DDocument.prototype.onMouseRelease=function(ctx,x,y,buttonState,prevButton){if(!ctx||!this._viewarea){return;}
var button=(prevButton<<8)|buttonState;ctx.pickValue(this._viewarea,x,y,button);this._viewarea.onMouseRelease(x,y,buttonState,prevButton);};x3dom.X3DDocument.prototype.onMouseOver=function(ctx,x,y,buttonState){if(!ctx||!this._viewarea){return;}
ctx.pickValue(this._viewarea,x,y,buttonState);this._viewarea.onMouseOver(x,y,buttonState);};x3dom.X3DDocument.prototype.onMouseOut=function(ctx,x,y,buttonState){if(!ctx||!this._viewarea){return;}
ctx.pickValue(this._viewarea,x,y,buttonState);this._viewarea.onMouseOut(x,y,buttonState);};x3dom.X3DDocument.prototype.onDoubleClick=function(ctx,x,y){if(!ctx||!this._viewarea){return;}
this._viewarea.onDoubleClick(x,y);};x3dom.X3DDocument.prototype.onKeyDown=function(keyCode)
{switch(keyCode){case 37:this._viewarea.strafeLeft();break;case 38:this._viewarea.moveFwd();break;case 39:this._viewarea.strafeRight();break;case 40:this._viewarea.moveBwd();break;default:}};x3dom.X3DDocument.prototype.onKeyUp=function(keyCode)
{var stack=null;switch(keyCode){case 13:x3dom.toggleFullScreen();break;case 33:stack=this._scene.getViewpoint()._stack;if(stack){stack.switchTo('next');}
else{x3dom.debug.logError('No valid ViewBindable stack.');}
break;case 34:stack=this._scene.getViewpoint()._stack;if(stack){stack.switchTo('prev');}
else{x3dom.debug.logError('No valid ViewBindable stack.');}
break;case 37:break;case 38:break;case 39:break;case 40:break;default:}};x3dom.X3DDocument.prototype.onKeyPress=function(charCode)
{var nav=this._scene.getNavigationInfo();var env=this._scene.getEnvironment();switch(charCode)
{case 32:var states=this.canvas.parent.stateViewer;if(states){states.display();}
x3dom.debug.logInfo("a: show all | d: show helper buffers | s: small feature culling | t: light view | "+"m: toggle render mode | c: frustum culling | p: intersect type | r: reset view | \n"+"e: examine mode | f: fly mode | y: freefly mode | w: walk mode | h: helicopter mode | "+"l: lookAt mode | o: lookaround | g: game mode | n: turntable | u: upright position | \n"+"v: print viewpoint info | pageUp: next view | pageDown: prev. view | "+"+: increase speed | -: decrease speed ");break;case 43:nav._vf.speed=2*nav._vf.speed;x3dom.debug.logInfo("Changed navigation speed to "+nav._vf.speed);break;case 45:nav._vf.speed=0.5*nav._vf.speed;x3dom.debug.logInfo("Changed navigation speed to "+nav._vf.speed);break;case 51:x3dom.nodeTypes.PopGeometry.ErrorToleranceFactor+=0.5;x3dom.debug.logInfo("Changed POP error tolerance to "+x3dom.nodeTypes.PopGeometry.ErrorToleranceFactor);break;case 52:x3dom.nodeTypes.PopGeometry.ErrorToleranceFactor-=0.5;x3dom.debug.logInfo("Changed POP error tolerance to "+x3dom.nodeTypes.PopGeometry.ErrorToleranceFactor);break;case 54:nav._vf.typeParams[1]+=1.0;nav._heliUpdated=false;x3dom.debug.logInfo("Changed helicopter height to "+nav._vf.typeParams[1]);break;case 55:nav._vf.typeParams[1]-=1.0;nav._heliUpdated=false;x3dom.debug.logInfo("Changed helicopter height to "+nav._vf.typeParams[1]);break;case 56:nav._vf.typeParams[0]-=0.02;nav._heliUpdated=false;x3dom.debug.logInfo("Changed helicopter angle to "+nav._vf.typeParams[0]);break;case 57:nav._vf.typeParams[0]+=0.02;nav._heliUpdated=false;x3dom.debug.logInfo("Changed helicopter angle to "+nav._vf.typeParams[0]);break;case 97:this._viewarea.showAll();break;case 99:env._vf.frustumCulling=!env._vf.frustumCulling;x3dom.debug.logInfo("Viewfrustum culling "+(env._vf.frustumCulling?"on":"off"));break;case 100:if(this._viewarea._visDbgBuf===undefined){this._viewarea._visDbgBuf=(this._x3dElem.getAttribute("showLog")==='true');}
this._viewarea._visDbgBuf=!this._viewarea._visDbgBuf;x3dom.debug.logContainer.style.display=(this._viewarea._visDbgBuf==true)?"block":"none";break;case 101:nav.setType("examine",this._viewarea);break;case 102:nav.setType("fly",this._viewarea);break;case 103:nav.setType("game",this._viewarea);break;case 104:nav.setType("helicopter",this._viewarea);break;case 105:this._viewarea.fit(this._scene._lastMin,this._scene._lastMax);break;case 108:nav.setType("lookat",this._viewarea);break;case 109:this._viewarea._points=++this._viewarea._points%2;break;case 110:nav.setType("turntable",this._viewarea);break;case 111:nav.setType("lookaround",this._viewarea);break;case 112:switch(this._scene._vf.pickMode.toLowerCase())
{case"idbuf":this._scene._vf.pickMode="color";break;case"color":this._scene._vf.pickMode="texCoord";break;case"texcoord":this._scene._vf.pickMode="box";break;default:this._scene._vf.pickMode="idBuf";break;}
x3dom.debug.logInfo("Switch pickMode to '"+this._scene._vf.pickMode+"'.");break;case 114:this._viewarea.resetView();break;case 115:env._vf.smallFeatureCulling=!env._vf.smallFeatureCulling;x3dom.debug.logInfo("Small feature culling "+(env._vf.smallFeatureCulling?"on":"off"));break;case 116:if(this._nodeBag.lights.length>0){this._viewarea.animateTo(this._viewarea.getLightMatrix()[0],this._scene.getViewpoint());}
break;case 117:this._viewarea.uprightView();break;case 118:var that=this;(function(){var viewpoint=that._viewarea._scene.getViewpoint();var mat_view=that._viewarea.getViewMatrix().inverse();var rotation=new x3dom.fields.Quaternion(0,0,1,0);rotation.setValue(mat_view);var rot=rotation.toAxisAngle();var translation=mat_view.e3();x3dom.debug.logInfo('\n&lt;Viewpoint position="'+translation.x.toFixed(5)+' '
+translation.y.toFixed(5)+' '+translation.z.toFixed(5)+'" '+'orientation="'+rot[0].x.toFixed(5)+' '+rot[0].y.toFixed(5)+' '
+rot[0].z.toFixed(5)+' '+rot[1].toFixed(5)+'" \n\t'+'zNear="'+viewpoint.getNear().toFixed(5)+'" '+'zFar="'+viewpoint.getFar().toFixed(5)+'" '+'description="'+viewpoint._vf.description+'"&gt;'+'&lt;/Viewpoint&gt;');})();break;case 119:nav.setType("walk",this._viewarea);break;case 121:nav.setType("freefly",this._viewarea);break;default:}};x3dom.X3DDocument.prototype.shutdown=function(ctx)
{if(!ctx){return;}
ctx.shutdown(this._viewarea);};x3dom.MatrixMixer=function(beginTime,endTime)
{this.beginTime=beginTime||0;this.endTime=endTime||1;this.isMixing=false;this._beginMat=x3dom.fields.SFMatrix4f.identity();this._beginInvMat=x3dom.fields.SFMatrix4f.identity();this._beginLogMat=x3dom.fields.SFMatrix4f.identity();this._endMat=x3dom.fields.SFMatrix4f.identity();this._endLogMat=x3dom.fields.SFMatrix4f.identity();this._beginRot=new x3dom.fields.Quaternion();this._endRot=new x3dom.fields.Quaternion();this._beginPos=new x3dom.fields.SFVec3f();this._endPos=new x3dom.fields.SFVec3f();this._result=x3dom.fields.SFMatrix4f.identity();this._useQuaternion=false;};x3dom.MatrixMixer.prototype._calcFraction=function(time)
{var fraction=(time-this.beginTime)/(this.endTime-this.beginTime);return(Math.sin((fraction*Math.PI)-(Math.PI/2))+1)/2.0;};x3dom.MatrixMixer.prototype._isValid=function()
{var angles=this._beginMat.inverse().mult(this._endMat).getEulerAngles();return(Math.abs(angles[0])!=Math.PI&&Math.abs(angles[1])!=Math.PI&&Math.abs(angles[2])!=Math.PI);};x3dom.MatrixMixer.prototype._prepareQuaternionAnimation=function()
{this._beginRot.setValue(this._beginMat);this._endRot.setValue(this._endMat);this._beginPos=this._beginMat.e3();this._endPos=this._endMat.e3();this._useQuaternion=true;};x3dom.MatrixMixer.prototype._reset=function()
{this.beginTime=0;this.endTime=0;this._useQuaternion=false;this.isMixing=false;};x3dom.MatrixMixer.prototype.isActive=function()
{return(this.beginTime>0);};x3dom.MatrixMixer.prototype.setBeginMatrix=function(mat)
{this._beginMat.setValues(mat);this._beginInvMat=mat.inverse();this._beginLogMat=x3dom.fields.SFMatrix4f.zeroMatrix();};x3dom.MatrixMixer.prototype.setEndMatrix=function(mat)
{this._endMat.setValues(mat);if(!this._isValid())
{this._prepareQuaternionAnimation();}
this._endLogMat=this._endMat.mult(this._beginInvMat).log();this._logDiffMat=this._endLogMat.addScaled(this._beginLogMat,-1);};x3dom.MatrixMixer.prototype._mixQuaternion=function(fraction)
{var rotation=this._beginRot.slerp(this._endRot,fraction);var translation=this._beginPos.addScaled(this._endPos.subtract(this._beginPos),fraction);this._result.setRotate(rotation);this._result.setTranslate(translation);return this._result.copy();};x3dom.MatrixMixer.prototype._mixMatrix=function(fraction)
{return this._logDiffMat.multiply(fraction).add(this._beginLogMat).exp().mult(this._beginMat);};x3dom.MatrixMixer.prototype.mix=function(time)
{if(time<=this.beginTime)
{return this._beginMat;}
else if(time>=this.endTime)
{this._reset();return this._endMat;}
else
{this.isMixing=true;var fraction=this._calcFraction(time);if(this._useQuaternion)
{return this._mixQuaternion(fraction);}
else
{return this._mixMatrix(fraction);}}};x3dom.InputTypes={NAVIGATION:1,INTERACTION:2};x3dom.Viewarea=function(document,scene){this._doc=document;this._scene=scene;document._nodeBag.viewarea.push(this);this._pickingInfo={pickPos:new x3dom.fields.SFVec3f(0,0,0),pickNorm:new x3dom.fields.SFVec3f(0,0,1),pickObj:null,firstObj:null,lastObj:null,lastClickObj:null,shadowObjectId:-1};this._currentInputType=x3dom.InputTypes.NAVIGATION;this._rotMat=x3dom.fields.SFMatrix4f.identity();this._transMat=x3dom.fields.SFMatrix4f.identity();this._movement=new x3dom.fields.SFVec3f(0,0,0);this._needNavigationMatrixUpdate=true;this._deltaT=0;this._flyMat=null;this._pitch=0;this._yaw=0;this._eyePos=new x3dom.fields.SFVec3f(0,0,0);this._width=400;this._height=300;this._dx=0;this._dy=0;this._lastX=-1;this._lastY=-1;this._pressX=-1;this._pressY=-1;this._lastButton=0;this._points=0;this._numRenderedNodes=0;this._pick=new x3dom.fields.SFVec3f(0,0,0);this._pickNorm=new x3dom.fields.SFVec3f(0,0,1);this._isAnimating=false;this._isMoving=false;this._lastTS=0;this._mixer=new x3dom.MatrixMixer();this._interpolator=new x3dom.FieldInterpolator();this.arc=null;};x3dom.Viewarea.prototype.tick=function(timeStamp)
{var needMixAnim=false;var env=this._scene.getEnvironment();if(env._vf.enableARC&&this.arc==null)
{this.arc=new x3dom.arc.AdaptiveRenderControl(this._scene);}
if(this._mixer.isActive())
{var mat=this._mixer.mix(timeStamp);this._scene.getViewpoint().setView(mat);}
if(this._interpolator.isActive())
{var value=this._interpolator.interpolate(timeStamp);this._scene.getViewpoint().setZoom(value);}
var needNavAnim=this.navigateTo(timeStamp);var lastIsAnimating=this._isAnimating;this._lastTS=timeStamp;this._isAnimating=(this._mixer.isMixing||this._interpolator.isInterpolating||needNavAnim);if(this.arc!=null)
{this.arc.update(this.isMovingOrAnimating()?1:0,this._doc._x3dElem.runtime.getFPS());}
return(this._isAnimating||lastIsAnimating);};x3dom.Viewarea.prototype.isMoving=function()
{return this._isMoving;};x3dom.Viewarea.prototype.isAnimating=function()
{return this._isAnimating;};x3dom.Viewarea.prototype.isMovingOrAnimating=function()
{return(this._isMoving||this._isAnimating);};x3dom.Viewarea.prototype.navigateTo=function(timeStamp)
{var navi=this._scene.getNavigationInfo();return navi._impl.navigateTo(this,timeStamp);};x3dom.Viewarea.prototype.moveFwd=function()
{var navi=this._scene.getNavigationInfo();navi._impl.moveForward(this);};x3dom.Viewarea.prototype.moveBwd=function()
{var navi=this._scene.getNavigationInfo();navi._impl.moveBackwards(this);};x3dom.Viewarea.prototype.strafeRight=function()
{var navi=this._scene.getNavigationInfo();navi._impl.strafeRight(this);};x3dom.Viewarea.prototype.strafeLeft=function()
{var navi=this._scene.getNavigationInfo();navi._impl.strafeLeft(this);};x3dom.Viewarea.prototype.animateTo=function(target,prev,dur)
{var navi=this._scene.getNavigationInfo();navi._impl.animateTo(this,target,prev,dur);};x3dom.Viewarea.prototype.orthoAnimateTo=function(target,prev,duration)
{var navi=this._scene.getNavigationInfo();navi._impl.orthoAnimateTo(this,target,prev,duration);};x3dom.Viewarea.prototype.zoom=function(zoomAmount)
{var navi=this._scene.getNavigationInfo();navi._impl.zoom(this,zoomAmount);};x3dom.Viewarea.prototype.getLights=function(){var enabledLights=[];for(var i=0;i<this._doc._nodeBag.lights.length;i++)
{if(this._doc._nodeBag.lights[i]._vf.on==true)
{enabledLights.push(this._doc._nodeBag.lights[i]);}}
return enabledLights;};x3dom.Viewarea.prototype.getLightsShadow=function(){var lights=this._doc._nodeBag.lights;for(var l=0;l<lights.length;l++){if(lights[l]._vf.shadowIntensity>0.0){return true;}}
return false;};x3dom.Viewarea.prototype.updateSpecialNavigation=function(viewpoint,mat_viewpoint){var navi=this._scene.getNavigationInfo();var navType=navi.getType();if(navType=="helicopter"&&!navi._heliUpdated)
{var typeParams=navi.getTypeParams();var theta=typeParams[0];var currViewMat=viewpoint.getViewMatrix().mult(mat_viewpoint.inverse()).inverse();this._from=currViewMat.e3();this._at=this._from.subtract(currViewMat.e2());this._up=new x3dom.fields.SFVec3f(0,1,0);this._from.y=typeParams[1];this._at.y=this._from.y;var sv=currViewMat.e0();var q=x3dom.fields.Quaternion.axisAngle(sv,theta);var temp=q.toMatrix();var fin=x3dom.fields.SFMatrix4f.translation(this._from);fin=fin.mult(temp);temp=x3dom.fields.SFMatrix4f.translation(this._from.negate());fin=fin.mult(temp);this._at=fin.multMatrixPnt(this._at);this._flyMat=x3dom.fields.SFMatrix4f.lookAt(this._from,this._at,this._up);this._scene.getViewpoint().setView(this._flyMat.inverse());navi._heliUpdated=true;}};x3dom.Viewarea.prototype.getViewpointMatrix=function()
{var viewpoint=this._scene.getViewpoint();var mat_viewpoint=viewpoint.getCurrentTransform();this.updateSpecialNavigation(viewpoint,mat_viewpoint);return viewpoint.getViewMatrix().mult(mat_viewpoint.inverse());};x3dom.Viewarea.prototype.getViewMatrix=function()
{return this.getViewpointMatrix().mult(this._transMat).mult(this._rotMat);};x3dom.Viewarea.prototype.getLightMatrix=function()
{var lights=this._doc._nodeBag.lights;var i,n=lights.length;if(n>0)
{var vol=this._scene.getVolume();if(vol.isValid())
{var min=x3dom.fields.SFVec3f.MAX();var max=x3dom.fields.SFVec3f.MIN();vol.getBounds(min,max);var l_arr=[];var viewpoint=this._scene.getViewpoint();var fov=viewpoint.getFieldOfView();var dia=max.subtract(min);var dist1=(dia.y/2.0)/Math.tan(fov/2.0)+(dia.z/2.0);var dist2=(dia.x/2.0)/Math.tan(fov/2.0)+(dia.z/2.0);dia=min.add(dia.multiply(0.5));for(i=0;i<n;i++)
{if(x3dom.isa(lights[i],x3dom.nodeTypes.PointLight)){var wcLoc=lights[i].getCurrentTransform().multMatrixPnt(lights[i]._vf.location);dia=dia.subtract(wcLoc).normalize();}
else{var dir=lights[i].getCurrentTransform().multMatrixVec(lights[i]._vf.direction);dir=dir.normalize().negate();dia=dia.add(dir.multiply(1.2*(dist1>dist2?dist1:dist2)));}
l_arr[i]=lights[i].getViewMatrix(dia);}
return l_arr;}}
return[this.getViewMatrix()];};x3dom.Viewarea.prototype.getWCtoLCMatrix=function(lMat)
{var proj=this.getProjectionMatrix();var view;if(arguments.length===0){view=this.getLightMatrix()[0];}
else{view=lMat;}
return proj.mult(view);};x3dom.Viewarea.prototype.getWCtoLCMatricesPointLight=function(view,lightNode,mat_proj)
{var zNear=lightNode._vf.zNear;var zFar=lightNode._vf.zFar;var proj=this.getLightProjectionMatrix(view,zNear,zFar,false,mat_proj);proj._00=1;proj._11=1;var matrices=[];matrices[0]=proj.mult(view);var rotationMatrix;for(var i=1;i<=3;i++){rotationMatrix=x3dom.fields.SFMatrix4f.rotationY(i*Math.PI/2);matrices[i]=proj.mult(rotationMatrix.mult(view));}
rotationMatrix=x3dom.fields.SFMatrix4f.rotationX(Math.PI/2);matrices[4]=proj.mult(rotationMatrix.mult(view));rotationMatrix=x3dom.fields.SFMatrix4f.rotationX(3*Math.PI/2);matrices[5]=proj.mult(rotationMatrix.mult(view));return matrices;};x3dom.Viewarea.prototype.getWCtoLCMatricesCascaded=function(view,lightNode,mat_proj)
{var numCascades=Math.max(1,Math.min(lightNode._vf.shadowCascades,6));var splitFactor=Math.max(0,Math.min(lightNode._vf.shadowSplitFactor,1));var splitOffset=Math.max(0,Math.min(lightNode._vf.shadowSplitOffset,1));var isSpotLight=x3dom.isa(lightNode,x3dom.nodeTypes.SpotLight);var zNear=lightNode._vf.zNear;var zFar=lightNode._vf.zFar;var proj=this.getLightProjectionMatrix(view,zNear,zFar,true,mat_proj);if(isSpotLight){proj._00=1;proj._11=1;}
var viewProj=proj.mult(view);var matrices=[];if(numCascades==1){matrices[0]=viewProj;return matrices;}
var cascadeSplits=this.getShadowSplitDepths(numCascades,splitFactor,splitOffset,true,mat_proj);for(var i=0;i<numCascades;i++){var fittingMat=this.getLightFittingMatrix(viewProj,cascadeSplits[i],cascadeSplits[i+1],mat_proj);matrices[i]=fittingMat.mult(viewProj);}
return matrices;};x3dom.Viewarea.prototype.getLightProjectionMatrix=function(lMat,zNear,zFar,highPrecision,mat_proj)
{var proj=x3dom.fields.SFMatrix4f.copy(mat_proj);if(!highPrecision||zNear>0||zFar>0){var lightPos=lMat.inverse().e3();var nearScale=0.8;var farScale=1.2;var min=x3dom.fields.SFVec3f.copy(this._scene._lastMin);var max=x3dom.fields.SFVec3f.copy(this._scene._lastMax);var dia=max.subtract(min);var sRad=dia.length()/2;var sCenter=min.add(dia.multiply(0.5));var vDist=(lightPos.subtract(sCenter)).length();var near,far;if(sRad){if(vDist>sRad)
near=(vDist-sRad)*nearScale;else
near=1;far=(vDist+sRad)*farScale;}
if(zNear>0)near=zNear;if(zFar>0)far=zFar;proj._22=-(far+near)/(far-near);proj._23=-2.0*far*near/(far-near);return proj;}
else{var cropMatrix=this.getLightCropMatrix(proj.mult(lMat));return cropMatrix.mult(proj);}};x3dom.Viewarea.prototype.getProjectionMatrix=function()
{var viewpoint=this._scene.getViewpoint();return viewpoint.getProjectionMatrix(this._width/this._height);};x3dom.Viewarea.prototype.getViewfrustum=function(clipMat)
{var env=this._scene.getEnvironment();if(env._vf.frustumCulling==true)
{if(arguments.length==0){var proj=this.getProjectionMatrix();var view=this.getViewMatrix();return new x3dom.fields.FrustumVolume(proj.mult(view));}
else{return new x3dom.fields.FrustumVolume(clipMat);}}
return null;};x3dom.Viewarea.prototype.getWCtoCCMatrix=function()
{var view=this.getViewMatrix();var proj=this.getProjectionMatrix();return proj.mult(view);};x3dom.Viewarea.prototype.getCCtoWCMatrix=function()
{var mat=this.getWCtoCCMatrix();return mat.inverse();};x3dom.Viewarea.prototype.calcViewRay=function(x,y,mat)
{var cctowc=mat?mat:this.getCCtoWCMatrix();var rx=x/(this._width-1.0)*2.0-1.0;var ry=(this._height-1.0-y)/(this._height-1.0)*2.0-1.0;var from=cctowc.multFullMatrixPnt(new x3dom.fields.SFVec3f(rx,ry,-1));var at=cctowc.multFullMatrixPnt(new x3dom.fields.SFVec3f(rx,ry,1));var dir=at.subtract(from);return new x3dom.fields.Ray(from,dir);};x3dom.Viewarea.prototype.showAll=function(axis,updateCenterOfRotation)
{if(axis===undefined)
axis="negZ";if(updateCenterOfRotation===undefined){updateCenterOfRotation=false;}
var scene=this._scene;scene.updateVolume();var min=x3dom.fields.SFVec3f.copy(scene._lastMin);var max=x3dom.fields.SFVec3f.copy(scene._lastMax);var x="x",y="y",z="z";var sign=1;var to,from=new x3dom.fields.SFVec3f(0,0,-1);switch(axis){case"posX":sign=-1;case"negX":z="x";x="y";y="z";to=new x3dom.fields.SFVec3f(sign,0,0);break;case"posY":sign=-1;case"negY":z="y";x="z";y="x";to=new x3dom.fields.SFVec3f(0,sign,0);break;case"posZ":sign=-1;case"negZ":default:to=new x3dom.fields.SFVec3f(0,0,-sign);break;}
var viewpoint=scene.getViewpoint();var fov=viewpoint.getFieldOfView();var isOrtho=x3dom.isa(viewpoint,x3dom.nodeTypes.OrthoViewpoint);var dia=max.subtract(min);var dia2=dia.multiply(0.5);var center=min.add(dia2);if(updateCenterOfRotation){viewpoint.setCenterOfRotation(center);}
var diaz2=dia[z]/2.0,tanfov2=Math.tan(fov/2.0);var dist1=(dia[y]/2.0)/tanfov2+diaz2;var dist2=(dia[x]/2.0)/tanfov2+diaz2;dia=min.add(dia.multiply(0.5));if(isOrtho)
{dia[z]+=sign*(dist1>dist2?dist1:dist2)*3.01;}
else
{dia[z]+=sign*(dist1>dist2?dist1:dist2)*1.01;}
var quat=x3dom.fields.Quaternion.rotateFromTo(from,to);var viewmat=quat.toMatrix();viewmat=viewmat.mult(x3dom.fields.SFMatrix4f.translation(dia.negate()));if(isOrtho)
{this.orthoAnimateTo(dist1,Math.abs(viewpoint._fieldOfView[0]));this.animateTo(viewmat,viewpoint);}
else
{this.animateTo(viewmat,viewpoint);}};x3dom.Viewarea.prototype.fit=function(min,max,updateCenterOfRotation)
{if(updateCenterOfRotation===undefined){updateCenterOfRotation=true;}
var dia2=max.subtract(min).multiply(0.5);var center=min.add(dia2);var bsr=dia2.length();var viewpoint=this._scene.getViewpoint();var fov=viewpoint.getFieldOfView();var viewmat=x3dom.fields.SFMatrix4f.copy(this.getViewMatrix());var rightDir=new x3dom.fields.SFVec3f(viewmat._00,viewmat._01,viewmat._02);var upDir=new x3dom.fields.SFVec3f(viewmat._10,viewmat._11,viewmat._12);var viewDir=new x3dom.fields.SFVec3f(viewmat._20,viewmat._21,viewmat._22);var tanfov2=Math.tan(fov/2.0);var dist=bsr/tanfov2;var eyePos=center.add(viewDir.multiply(dist));viewmat._03=-rightDir.dot(eyePos);viewmat._13=-upDir.dot(eyePos);viewmat._23=-viewDir.dot(eyePos);if(updateCenterOfRotation){viewpoint.setCenterOfRotation(center);}
if(x3dom.isa(viewpoint,x3dom.nodeTypes.OrthoViewpoint))
{this.orthoAnimateTo(dist/2.01,Math.abs(viewpoint._fieldOfView[0]));this.animateTo(viewmat,viewpoint);}
else
{this.animateTo(viewmat,viewpoint);}};x3dom.Viewarea.prototype.resetView=function()
{var navi=this._scene.getNavigationInfo();navi._impl.resetView(this);};x3dom.Viewarea.prototype.resetNavHelpers=function()
{this._rotMat=x3dom.fields.SFMatrix4f.identity();this._transMat=x3dom.fields.SFMatrix4f.identity();this._movement=new x3dom.fields.SFVec3f(0,0,0);this._needNavigationMatrixUpdate=true;};x3dom.Viewarea.prototype.uprightView=function()
{var mat=this.getViewMatrix().inverse();var from=mat.e3();var at=from.subtract(mat.e2());var up=new x3dom.fields.SFVec3f(0,1,0);var s=mat.e2().cross(up).normalize();var v=s.cross(up).normalize();at=from.add(v);mat=x3dom.fields.SFMatrix4f.lookAt(from,at,up);mat=mat.inverse();this.animateTo(mat,this._scene.getViewpoint());};x3dom.Viewarea.prototype.callEvtHandler=function(node,eventType,event)
{if(!node||!node._xmlNode)
return null;try{var attrib=node._xmlNode[eventType];if(typeof(attrib)==="function"){attrib.call(node._xmlNode,event);}
else{var funcStr=node._xmlNode.getAttribute(eventType);var func=new Function('event',funcStr);func.call(node._xmlNode,event);}
var list=node._listeners[event.type];if(list){for(var it=0;it<list.length;it++){list[it].call(node._xmlNode,event);}}}
catch(e){x3dom.debug.logException(e);}
return event.cancelBubble;};x3dom.Viewarea.prototype.checkEvents=function(obj,x,y,buttonState,eventType)
{var that=this;var needRecurse=true;var childNode;var i,n;var target=(obj&&obj._xmlNode)?obj._xmlNode:{};var affectedPointingSensorsList=this._doc._nodeBag.affectedPointingSensors;var event={viewarea:that,target:target,type:eventType.substr(2,eventType.length-2),button:buttonState,layerX:x,layerY:y,worldX:that._pick.x,worldY:that._pick.y,worldZ:that._pick.z,normalX:that._pickNorm.x,normalY:that._pickNorm.y,normalZ:that._pickNorm.z,hitPnt:that._pick.toGL(),hitObject:target,shadowObjectId:that._pickingInfo.shadowObjectId,cancelBubble:false,stopPropagation:function(){this.cancelBubble=true;},preventDefault:function(){this.cancelBubble=true;}};try{var anObj=obj;if(anObj&&anObj._xmlNode&&anObj._cf.geometry&&!anObj._xmlNode[eventType]&&!anObj._xmlNode.hasAttribute(eventType)&&!anObj._listeners[event.type]){anObj=anObj._cf.geometry.node;}
if(anObj&&that.callEvtHandler(anObj,eventType,event)===true){needRecurse=false;}}
catch(e){x3dom.debug.logException(e);}
var recurse=function(obj){Array.forEach(obj._parentNodes,function(node){if(node._xmlNode&&(node._xmlNode[eventType]||node._xmlNode.hasAttribute(eventType)||node._listeners[event.type]))
{if(that.callEvtHandler(node,eventType,event)===true){needRecurse=false;}}
if(buttonState==0&&affectedPointingSensorsList.length==0&&(eventType=='onmousemove'||eventType=='onmouseover'||eventType=='onmouseout'))
{n=node._childNodes.length;for(i=0;i<n;++i)
{childNode=node._childNodes[i];if(x3dom.isa(childNode,x3dom.nodeTypes.X3DPointingDeviceSensorNode)&&childNode._vf.enabled)
{affectedPointingSensorsList.push(childNode);}}}
if(x3dom.isa(node,x3dom.nodeTypes.Anchor)&&eventType==='onclick'){node.handleTouch();needRecurse=false;}
else if(needRecurse){recurse(node);}});};if(needRecurse&&obj){recurse(obj);}
return needRecurse;};x3dom.Viewarea.prototype._notifyAffectedPointingSensors=function(event)
{var funcDict={"mousedown":"pointerPressedOverSibling","mousemove":"pointerMoved","mouseover":"pointerMovedOver","mouseout":"pointerMovedOut"};var func=funcDict[event.type];var affectedPointingSensorsList=this._doc._nodeBag.affectedPointingSensors;var i,n=affectedPointingSensorsList.length;if(n>0&&func!==undefined)
{for(i=0;i<n;i++)
affectedPointingSensorsList[i][func](event);}};x3dom.Viewarea.prototype.initMouseState=function()
{this._deltaT=0;this._dx=0;this._dy=0;this._lastX=-1;this._lastY=-1;this._pressX=-1;this._pressY=-1;this._lastButton=0;this._isMoving=false;this._needNavigationMatrixUpdate=true;};x3dom.Viewarea.prototype.onMousePress=function(x,y,buttonState)
{this._needNavigationMatrixUpdate=true;this.prepareEvents(x,y,buttonState,"onmousedown");this._pickingInfo.lastClickObj=this._pickingInfo.pickObj;this._pickingInfo.firstObj=this._pickingInfo.pickObj;this._dx=0;this._dy=0;this._lastX=x;this._lastY=y;this._pressX=x;this._pressY=y;this._lastButton=buttonState;this._isMoving=false;if(this._currentInputType==x3dom.InputTypes.NAVIGATION)
{var navi=this._scene.getNavigationInfo();navi._impl.onMousePress(this,x,y,buttonState);}};x3dom.Viewarea.prototype.onMouseRelease=function(x,y,buttonState,prevButton)
{var i;var affectedPointingSensorsList=this._doc._nodeBag.affectedPointingSensors;for(i=0;i<affectedPointingSensorsList.length;++i)
{affectedPointingSensorsList[i].pointerReleased();}
this._doc._nodeBag.affectedPointingSensors=[];var tDist=3.0;var dir;var navi=this._scene.getNavigationInfo();var navType=navi.getType();if(this._scene._vf.pickMode.toLowerCase()!=="box"){this.prepareEvents(x,y,prevButton,"onmouseup");if(this._pickingInfo.pickObj&&this._pickingInfo.pickObj===this._pickingInfo.lastClickObj)
{this.prepareEvents(x,y,prevButton,"onclick");}
else if(!this._pickingInfo.pickObj&&!this._pickingInfo.lastClickObj&&!this._pickingInfo.firstObj)
{var eventType="backgroundClicked";try{if(this._scene._xmlNode&&(this._scene._xmlNode["on"+eventType]||this._scene._xmlNode.hasAttribute("on"+eventType)||this._scene._listeners[eventType])){var event={target:this._scene._xmlNode,type:eventType,button:prevButton,layerX:x,layerY:y,cancelBubble:false,stopPropagation:function(){this.cancelBubble=true;},preventDefault:function(){this.cancelBubble=true;}};this._scene.callEvtHandler(("on"+eventType),event);}}
catch(e){x3dom.debug.logException("backgroundClicked: "+e);}}}
else{var t0=new Date().getTime();var line=this.calcViewRay(x,y);var isect=this._scene.doIntersect(line);var obj=line.hitObject;if(isect&&obj)
{this._pick.setValues(line.hitPoint);this.checkEvents(obj,x,y,buttonState,"onclick");x3dom.debug.logInfo("Hit '"+obj._xmlNode.localName+"/ "+
obj._DEF+"' at dist="+line.dist.toFixed(4));x3dom.debug.logInfo("Ray hit at position "+this._pick);}
var t1=new Date().getTime()-t0;x3dom.debug.logInfo("Picking time (box): "+t1+"ms");if(!isect){dir=this.getViewMatrix().e2().negate();var u=dir.dot(line.pos.negate())/dir.dot(line.dir);this._pick=line.pos.add(line.dir.multiply(u));}}
this._pickingInfo.firstObj=null;if(this._currentInputType==x3dom.InputTypes.NAVIGATION&&(this._pickingInfo.pickObj||this._pickingInfo.shadowObjectId>=0)&&navType==="lookat"&&this._pressX===x&&this._pressY===y)
{var step=(this._lastButton&2)?-1:1;var dist=this._pickingInfo.pickPos.subtract(this._from).length()/tDist;var laMat=new x3dom.fields.SFMatrix4f();laMat.setValues(this.getViewMatrix());laMat=laMat.inverse();var from=laMat.e3();var at=from.subtract(laMat.e2());var up=laMat.e1();dir=this._pickingInfo.pickPos.subtract(from);var len=dir.length();dir=dir.normalize();var newAt=from.addScaled(dir,len);var s=dir.cross(up).normalize();dir=s.cross(up).normalize();if(step<0){dist=(0.5+len+dist)*2;}
var newFrom=newAt.addScaled(dir,dist);laMat=x3dom.fields.SFMatrix4f.lookAt(newFrom,newAt,up);laMat=laMat.inverse();dist=newFrom.subtract(from).length();var dur=Math.max(0.5,Math.log((1+dist)/navi._vf.speed));this.animateTo(laMat,this._scene.getViewpoint(),dur);}
this._dx=0;this._dy=0;this._lastX=x;this._lastY=y;this._lastButton=buttonState;this._isMoving=false;};x3dom.Viewarea.prototype.onMouseOver=function(x,y,buttonState)
{this._dx=0;this._dy=0;this._lastButton=0;this._isMoving=false;this._lastX=x;this._lastY=y;this._deltaT=0;};x3dom.Viewarea.prototype.onMouseOut=function(x,y,buttonState)
{this._dx=0;this._dy=0;this._lastButton=0;this._isMoving=false;this._lastX=x;this._lastY=y;this._deltaT=0;var i;var affectedPointingSensorsList=this._doc._nodeBag.affectedPointingSensors;for(i=0;i<affectedPointingSensorsList.length;++i)
{affectedPointingSensorsList[i].pointerReleased();}
this._doc._nodeBag.affectedPointingSensors=[];};x3dom.Viewarea.prototype.onDoubleClick=function(x,y)
{if(this._doc._x3dElem.hasAttribute('disableDoubleClick')&&this._doc._x3dElem.getAttribute('disableDoubleClick')==='true'){return;}
var navi=this._scene.getNavigationInfo();navi._impl.onDoubleClick(this,x,y);};x3dom.Viewarea.prototype.handleMoveEvt=function(x,y,buttonState)
{if(buttonState==0)
{this._doc._nodeBag.affectedPointingSensors=[];}
this.prepareEvents(x,y,buttonState,"onmousemove");if(this._pickingInfo.pickObj!==this._pickingInfo.lastObj)
{if(this._pickingInfo.lastObj){var obj=this._pickingInfo.pickObj;this._pickingInfo.pickObj=this._pickingInfo.lastObj;this.prepareEvents(x,y,buttonState,"onmouseout");this._pickingInfo.pickObj=obj;}
if(this._pickingInfo.pickObj){this.prepareEvents(x,y,buttonState,"onmouseover");}
this._pickingInfo.lastObj=this._pickingInfo.pickObj;}};x3dom.Viewarea.prototype.onMove=function(x,y,buttonState)
{this.handleMoveEvt(x,y,buttonState);if(this._lastX<0||this._lastY<0){this._lastX=x;this._lastY=y;}
this._dx=x-this._lastX;this._dy=y-this._lastY;this._lastX=x;this._lastY=y;};x3dom.Viewarea.prototype.onMoveView=function(translation,rotation)
{if(this._currentInputType==x3dom.InputTypes.NAVIGATION)
{var navi=this._scene.getNavigationInfo();var viewpoint=this._scene.getViewpoint();if(navi.getType()==="examine")
{if(translation)
{var distance=(this._scene._lastMax.subtract(this._scene._lastMin)).length();distance=((distance<x3dom.fields.Eps)?1:distance)*navi._vf.speed;translation=translation.multiply(distance);this._movement=this._movement.add(translation);this._transMat=viewpoint.getViewMatrix().inverse().mult(x3dom.fields.SFMatrix4f.translation(this._movement)).mult(viewpoint.getViewMatrix());}
if(rotation)
{var center=viewpoint.getCenterOfRotation();var mat=this.getViewMatrix();mat.setTranslate(new x3dom.fields.SFVec3f(0,0,0));this._rotMat=this._rotMat.mult(x3dom.fields.SFMatrix4f.translation(center)).mult(mat.inverse()).mult(rotation).mult(mat).mult(x3dom.fields.SFMatrix4f.translation(center.negate()));}
this._isMoving=true;}}};x3dom.Viewarea.prototype.onDrag=function(x,y,buttonState)
{this.handleMoveEvt(x,y,buttonState);if(this._currentInputType==x3dom.InputTypes.NAVIGATION)
{this._scene.getNavigationInfo()._impl.onDrag(this,x,y,buttonState);}};x3dom.Viewarea.prototype.prepareEvents=function(x,y,buttonState,eventType)
{var affectedPointingSensorsList=this._doc._nodeBag.affectedPointingSensors;var pickMode=this._scene._vf.pickMode.toLowerCase();var avoidTraversal=(pickMode.indexOf("idbuf")==0||pickMode=="color"||pickMode=="texcoord");var obj=null;if(avoidTraversal){obj=this._pickingInfo.pickObj;if(obj){this._pick.setValues(this._pickingInfo.pickPos);this._pickNorm.setValues(this._pickingInfo.pickNorm);this.checkEvents(obj,x,y,buttonState,eventType);if(eventType==="onclick"){if(obj._xmlNode)
x3dom.debug.logInfo("Hit \""+obj._xmlNode.localName+"/ "+obj._DEF+"\"");x3dom.debug.logInfo("Ray hit at position "+this._pick);}}}
var event={viewarea:this,target:{},type:eventType.substr(2,eventType.length-2),button:buttonState,layerX:x,layerY:y,worldX:this._pick.x,worldY:this._pick.y,worldZ:this._pick.z,normalX:this._pickNorm.x,normalY:this._pickNorm.y,normalZ:this._pickNorm.z,hitPnt:this._pick.toGL(),hitObject:(obj&&obj._xmlNode)?obj._xmlNode:null,shadowObjectId:this._pickingInfo.shadowObjectId,cancelBubble:false,stopPropagation:function(){this.cancelBubble=true;},preventDefault:function(){this.cancelBubble=true;}};this._notifyAffectedPointingSensors(event);if(affectedPointingSensorsList.length>0)
{this._currentInputType=x3dom.InputTypes.INTERACTION;}
else
{this._currentInputType=x3dom.InputTypes.NAVIGATION;}};x3dom.Viewarea.prototype.getRenderMode=function()
{return this._points;};x3dom.Viewarea.prototype.getShadowedLights=function()
{var shadowedLights=[];var shadowIndex=0;var slights=this.getLights();for(var i=0;i<slights.length;i++){if(slights[i]._vf.shadowIntensity>0.0){shadowedLights[shadowIndex]=slights[i];shadowIndex++;}}
return shadowedLights;};x3dom.Viewarea.prototype.getShadowSplitDepths=function(numCascades,splitFactor,splitOffset,postProject,mat_proj)
{var logSplit;var practSplit=[];var viewPoint=this._scene.getViewpoint();var zNear=viewPoint.getNear();var zFar=viewPoint.getFar();practSplit[0]=zNear;zNear=zNear+splitOffset*(zFar-zNear)/10;for(var i=1;i<numCascades;i++){logSplit=zNear*Math.pow((zFar/zNear),i/numCascades);practSplit[i]=splitFactor*logSplit+(1-splitFactor)*(zNear+i/(numCascades*(zNear-zFar)));}
practSplit[numCascades]=zFar;if(!postProject)
return practSplit;var postProj=[];for(var j=0;j<=numCascades;j++){postProj[j]=mat_proj.multFullMatrixPnt(new x3dom.fields.SFVec3f(0,0,-practSplit[j])).z;}
return postProj;};x3dom.Viewarea.prototype.getLightCropMatrix=function(WCToLCMatrix)
{var sceneMin=x3dom.fields.SFVec3f.copy(this._scene._lastMin);var sceneMax=x3dom.fields.SFVec3f.copy(this._scene._lastMax);var sceneCorners=[];sceneCorners[0]=new x3dom.fields.SFVec3f(sceneMin.x,sceneMin.y,sceneMin.z);sceneCorners[1]=new x3dom.fields.SFVec3f(sceneMin.x,sceneMin.y,sceneMax.z);sceneCorners[2]=new x3dom.fields.SFVec3f(sceneMin.x,sceneMax.y,sceneMin.z);sceneCorners[3]=new x3dom.fields.SFVec3f(sceneMin.x,sceneMax.y,sceneMax.z);sceneCorners[4]=new x3dom.fields.SFVec3f(sceneMax.x,sceneMin.y,sceneMin.z);sceneCorners[5]=new x3dom.fields.SFVec3f(sceneMax.x,sceneMin.y,sceneMax.z);sceneCorners[6]=new x3dom.fields.SFVec3f(sceneMax.x,sceneMax.y,sceneMin.z);sceneCorners[7]=new x3dom.fields.SFVec3f(sceneMax.x,sceneMax.y,sceneMax.z);var i;for(i=0;i<8;i++){sceneCorners[i]=WCToLCMatrix.multFullMatrixPnt(sceneCorners[i]);}
var minScene=x3dom.fields.SFVec3f.copy(sceneCorners[0]);var maxScene=x3dom.fields.SFVec3f.copy(sceneCorners[0]);for(i=1;i<8;i++){minScene.z=Math.min(sceneCorners[i].z,minScene.z);maxScene.z=Math.max(sceneCorners[i].z,maxScene.z);}
var scaleZ=2.0/(maxScene.z-minScene.z);var offsetZ=-(scaleZ*(maxScene.z+minScene.z))/2.0;var cropMatrix=x3dom.fields.SFMatrix4f.identity();cropMatrix._22=scaleZ;cropMatrix._23=offsetZ;return cropMatrix;};x3dom.Viewarea.prototype.getLightFittingMatrix=function(WCToLCMatrix,zNear,zFar,mat_proj)
{var mat_view=this.getViewMatrix();var mat_view_proj=mat_proj.mult(mat_view);var mat_view_proj_inverse=mat_view_proj.inverse();var frustumCorners=[];frustumCorners[0]=new x3dom.fields.SFVec3f(-1,-1,zFar);frustumCorners[1]=new x3dom.fields.SFVec3f(-1,-1,zNear);frustumCorners[2]=new x3dom.fields.SFVec3f(-1,1,zFar);frustumCorners[3]=new x3dom.fields.SFVec3f(-1,1,zNear);frustumCorners[4]=new x3dom.fields.SFVec3f(1,-1,zFar);frustumCorners[5]=new x3dom.fields.SFVec3f(1,-1,zNear);frustumCorners[6]=new x3dom.fields.SFVec3f(1,1,zFar);frustumCorners[7]=new x3dom.fields.SFVec3f(1,1,zNear);var i;for(i=0;i<8;i++){frustumCorners[i]=mat_view_proj_inverse.multFullMatrixPnt(frustumCorners[i]);frustumCorners[i]=WCToLCMatrix.multFullMatrixPnt(frustumCorners[i]);}
var minFrustum=x3dom.fields.SFVec3f.copy(frustumCorners[0]);var maxFrustum=x3dom.fields.SFVec3f.copy(frustumCorners[0]);for(i=1;i<8;i++){minFrustum.x=Math.min(frustumCorners[i].x,minFrustum.x);minFrustum.y=Math.min(frustumCorners[i].y,minFrustum.y);minFrustum.z=Math.min(frustumCorners[i].z,minFrustum.z);maxFrustum.x=Math.max(frustumCorners[i].x,maxFrustum.x);maxFrustum.y=Math.max(frustumCorners[i].y,maxFrustum.y);maxFrustum.z=Math.max(frustumCorners[i].z,maxFrustum.z);}
function clip(min,max)
{var xMin=min.x;var yMin=min.y;var zMin=min.z;var xMax=max.x;var yMax=max.y;var zMax=max.z;if(xMin>1.0||xMax<-1.0){xMin=-1.0;xMax=1.0;}else{xMin=Math.max(xMin,-1.0);xMax=Math.min(xMax,1.0);}
if(yMin>1.0||yMax<-1.0){yMin=-1.0;yMax=1.0;}else{yMin=Math.max(yMin,-1.0);yMax=Math.min(yMax,1.0);}
if(zMin>1.0||zMax<-1.0){zMin=-1.0;zMax=1.0;}else{zMin=Math.max(zMin,-1.0);zMax=Math.min(zMax,1.0);}
var minValues=new x3dom.fields.SFVec3f(xMin,yMin,zMin);var maxValues=new x3dom.fields.SFVec3f(xMax,yMax,zMax);return new x3dom.fields.BoxVolume(minValues,maxValues);}
var frustumBB=clip(minFrustum,maxFrustum);var scaleX=2.0/(frustumBB.max.x-frustumBB.min.x);var scaleY=2.0/(frustumBB.max.y-frustumBB.min.y);var offsetX=-(scaleX*(frustumBB.max.x+frustumBB.min.x))/2.0;var offsetY=-(scaleY*(frustumBB.max.y+frustumBB.min.y))/2.0;var fittingMatrix=x3dom.fields.SFMatrix4f.identity();fittingMatrix._00=scaleX;fittingMatrix._11=scaleY;fittingMatrix._03=offsetX;fittingMatrix._13=offsetY;return fittingMatrix;};x3dom.Mesh=function(parent)
{this._parent=parent;this._vol=new x3dom.fields.BoxVolume();this._invalidate=true;this._numFaces=0;this._numCoords=0;this._primType='TRIANGLES';this._positions=[];this._normals=[];this._texCoords=[];this._colors=[];this._indices=[];this._positions[0]=[];this._normals[0]=[];this._texCoords[0]=[];this._colors[0]=[];this._indices[0]=[];};x3dom.Mesh.prototype._dynamicFields={};x3dom.Mesh.prototype._numPosComponents=3;x3dom.Mesh.prototype._numTexComponents=2;x3dom.Mesh.prototype._numColComponents=3;x3dom.Mesh.prototype._numNormComponents=3;x3dom.Mesh.prototype._lit=true;x3dom.Mesh.prototype._vol=null;x3dom.Mesh.prototype._invalidate=true;x3dom.Mesh.prototype._numFaces=0;x3dom.Mesh.prototype._numCoords=0;x3dom.Mesh.prototype.setMeshData=function(positions,normals,texCoords,colors,indices)
{this._positions[0]=positions;this._normals[0]=normals;this._texCoords[0]=texCoords;this._colors[0]=colors;this._indices[0]=indices;this._invalidate=true;this._numFaces=this._indices[0].length/3;this._numCoords=this._positions[0].length/3;};x3dom.Mesh.prototype.getVolume=function()
{if(this._invalidate==true&&!this._vol.isValid())
{var coords=this._positions[0];var n=coords.length;if(n>3)
{var initVal=new x3dom.fields.SFVec3f(coords[0],coords[1],coords[2]);this._vol.setBounds(initVal,initVal);for(var i=3;i<n;i+=3)
{if(this._vol.min.x>coords[i]){this._vol.min.x=coords[i];}
if(this._vol.min.y>coords[i+1]){this._vol.min.y=coords[i+1];}
if(this._vol.min.z>coords[i+2]){this._vol.min.z=coords[i+2];}
if(this._vol.max.x<coords[i]){this._vol.max.x=coords[i];}
if(this._vol.max.y<coords[i+1]){this._vol.max.y=coords[i+1];}
if(this._vol.max.z<coords[i+2]){this._vol.max.z=coords[i+2];}}
this._invalidate=false;}}
return this._vol;};x3dom.Mesh.prototype.invalidate=function()
{this._invalidate=true;this._vol.invalidate();};x3dom.Mesh.prototype.isValid=function()
{return this._vol.isValid();};x3dom.Mesh.prototype.getCenter=function()
{return this.getVolume().getCenter();};x3dom.Mesh.prototype.getDiameter=function()
{return this.getVolume().getDiameter();};x3dom.Mesh.prototype.doIntersect=function(line)
{var vol=this.getVolume();var isect=line.intersect(vol.min,vol.max);if(isect&&line.enter<line.dist)
{line.dist=line.enter;line.hitObject=this._parent;line.hitPoint=line.pos.add(line.dir.multiply(line.enter));}
return isect;};x3dom.Mesh.prototype.calcNormals=function(creaseAngle,ccw)
{if(ccw===undefined)
ccw=true;var multInd=this._multiIndIndices&&this._multiIndIndices.length;var idxs=multInd?this._multiIndIndices:this._indices[0];var coords=this._positions[0];var vertNormals=[];var vertFaceNormals=[];var i,j,m=coords.length;var a,b,n=null;var num=(this._posSize!==undefined&&this._posSize>m)?this._posSize/3:m/3;num=3*((num-Math.floor(num)>0)?Math.floor(num+1):num);for(i=0;i<num;++i){vertFaceNormals[i]=[];}
num=idxs.length;for(i=0;i<num;i+=3){var ind_i0,ind_i1,ind_i2;var t;if(!multInd){ind_i0=idxs[i]*3;ind_i1=idxs[i+1]*3;ind_i2=idxs[i+2]*3;t=new x3dom.fields.SFVec3f(coords[ind_i1],coords[ind_i1+1],coords[ind_i1+2]);a=new x3dom.fields.SFVec3f(coords[ind_i0],coords[ind_i0+1],coords[ind_i0+2]).subtract(t);b=t.subtract(new x3dom.fields.SFVec3f(coords[ind_i2],coords[ind_i2+1],coords[ind_i2+2]));ind_i0=i*3;ind_i1=(i+1)*3;ind_i2=(i+2)*3;}
else{ind_i0=i*3;ind_i1=(i+1)*3;ind_i2=(i+2)*3;t=new x3dom.fields.SFVec3f(coords[ind_i1],coords[ind_i1+1],coords[ind_i1+2]);a=new x3dom.fields.SFVec3f(coords[ind_i0],coords[ind_i0+1],coords[ind_i0+2]).subtract(t);b=t.subtract(new x3dom.fields.SFVec3f(coords[ind_i2],coords[ind_i2+1],coords[ind_i2+2]));}
n=a.cross(b).normalize();if(!ccw)
n=n.negate();if(creaseAngle<=x3dom.fields.Eps){vertNormals[ind_i0]=vertNormals[ind_i1]=vertNormals[ind_i2]=n.x;vertNormals[ind_i0+1]=vertNormals[ind_i1+1]=vertNormals[ind_i2+1]=n.y;vertNormals[ind_i0+2]=vertNormals[ind_i1+2]=vertNormals[ind_i2+2]=n.z;}
else{vertFaceNormals[idxs[i]].push(n);vertFaceNormals[idxs[i+1]].push(n);vertFaceNormals[idxs[i+2]].push(n);}}
if(creaseAngle>x3dom.fields.Eps)
{for(i=0;i<m;i+=3){var iThird=i/3;var arr;if(!multInd){arr=vertFaceNormals[iThird];}
else{arr=vertFaceNormals[idxs[iThird]];}
num=arr.length;n=new x3dom.fields.SFVec3f(0,0,0);for(j=0;j<num;++j){n=n.add(arr[j]);}
n=n.normalize();vertNormals[i]=n.x;vertNormals[i+1]=n.y;vertNormals[i+2]=n.z;}}
this._normals[0]=vertNormals;};x3dom.Mesh.prototype.splitMesh=function(primStride,checkMultiIndIndices)
{var pStride;var isMultiInd;if(typeof primStride===undefined){pStride=3;}else{pStride=primStride;}
if(typeof checkMultiIndIndices===undefined){checkMultiIndIndices=false;}
var MAX=x3dom.Utils.maxIndexableCoords;MAX=Math.floor(MAX/pStride)*pStride;if(this._positions[0].length/3<=MAX&&!checkMultiIndIndices){return;}
if(checkMultiIndIndices){isMultiInd=this._multiIndIndices&&this._multiIndIndices.length;}else{isMultiInd=false;}
var positions=this._positions[0];var normals=this._normals[0];var texCoords=this._texCoords[0];var colors=this._colors[0];var indices=isMultiInd?this._multiIndIndices:this._indices[0];var i=0;do
{this._positions[i]=[];this._normals[i]=[];this._texCoords[i]=[];this._colors[i]=[];this._indices[i]=[];var k=(indices.length-((i+1)*MAX)>=0);if(k){this._indices[i]=indices.slice(i*MAX,(i+1)*MAX);}else{this._indices[i]=indices.slice(i*MAX);}
if(!isMultiInd){if(i){var m=i*MAX;for(var j=0,l=this._indices[i].length;j<l;j++){this._indices[i][j]-=m;}}}else{for(var j=0,l=this._indices[i].length;j<l;j++){this._indices[i][j]=j;}}
if(k){this._positions[i]=positions.slice(i*MAX*3,3*(i+1)*MAX);}else{this._positions[i]=positions.slice(i*MAX*3);}
if(normals.length){if(k){this._normals[i]=normals.slice(i*MAX*3,3*(i+1)*MAX);}else{this._normals[i]=normals.slice(i*MAX*3);}}
if(texCoords.length){if(k){this._texCoords[i]=texCoords.slice(i*MAX*this._numTexComponents,this._numTexComponents*(i+1)*MAX);}else{this._texCoords[i]=texCoords.slice(i*MAX*this._numTexComponents);}}
if(colors.length){if(k){this._colors[i]=colors.slice(i*MAX*this._numColComponents,this._numColComponents*(i+1)*MAX);}else{this._colors[i]=colors.slice(i*MAX*this._numColComponents);}}}
while(positions.length>++i*MAX*3);};x3dom.Mesh.prototype.calcTexCoords=function(mode)
{this._texCoords[0]=[];if(mode.toLowerCase()==="sphere-local")
{for(var i=0,j=0,n=this._normals[0].length;i<n;i+=3)
{this._texCoords[0][j++]=0.5+this._normals[0][i]/2.0;this._texCoords[0][j++]=0.5+this._normals[0][i+1]/2.0;}}
else
{var min=new x3dom.fields.SFVec3f(0,0,0),max=new x3dom.fields.SFVec3f(0,0,0);var vol=this.getVolume();vol.getBounds(min,max);var dia=max.subtract(min);var S=0,T=1;if(dia.x>=dia.y)
{if(dia.x>=dia.z)
{S=0;T=dia.y>=dia.z?1:2;}
else
{S=2;T=0;}}
else
{if(dia.y>=dia.z)
{S=1;T=dia.x>=dia.z?0:2;}
else
{S=2;T=1;}}
var sDenom=1,tDenom=1;var sMin=0,tMin=0;switch(S){case 0:sDenom=dia.x;sMin=min.x;break;case 1:sDenom=dia.y;sMin=min.y;break;case 2:sDenom=dia.z;sMin=min.z;break;}
switch(T){case 0:tDenom=dia.x;tMin=min.x;break;case 1:tDenom=dia.y;tMin=min.y;break;case 2:tDenom=dia.z;tMin=min.z;break;}
for(var k=0,l=0,m=this._positions[0].length;k<m;k+=3)
{this._texCoords[0][l++]=(this._positions[0][k+S]-sMin)/sDenom;this._texCoords[0][l++]=(this._positions[0][k+T]-tMin)/tDenom;}}};if(typeof x3dom==="undefined")
{x3dom={extend:function(f){function G(){}
G.prototype=f.prototype||f;return new G();},debug:{logInfo:function(msg){console.log(msg);},logWarning:function(msg){console.warn(msg);},logError:function(msg){console.error(msg);}}};if(!Array.map){Array.map=function(array,fun,thisp){var len=array.length;var res=[];for(var i=0;i<len;i++){if(i in array){res[i]=fun.call(thisp,array[i],i,array);}}
return res;};}
console.log("Using x3dom fields.js as standalone math and/or base types library.");}
x3dom.fields={};var VecMath=x3dom.fields;x3dom.fields.Eps=0.000001;x3dom.fields.SFMatrix4f=function(_00,_01,_02,_03,_10,_11,_12,_13,_20,_21,_22,_23,_30,_31,_32,_33)
{if(arguments.length===0){this._00=1;this._01=0;this._02=0;this._03=0;this._10=0;this._11=1;this._12=0;this._13=0;this._20=0;this._21=0;this._22=1;this._23=0;this._30=0;this._31=0;this._32=0;this._33=1;}
else{this._00=_00;this._01=_01;this._02=_02;this._03=_03;this._10=_10;this._11=_11;this._12=_12;this._13=_13;this._20=_20;this._21=_21;this._22=_22;this._23=_23;this._30=_30;this._31=_31;this._32=_32;this._33=_33;}};x3dom.fields.SFMatrix4f.prototype.e0=function(){var baseVec=new x3dom.fields.SFVec3f(this._00,this._10,this._20);return baseVec.normalize();};x3dom.fields.SFMatrix4f.prototype.e1=function(){var baseVec=new x3dom.fields.SFVec3f(this._01,this._11,this._21);return baseVec.normalize();};x3dom.fields.SFMatrix4f.prototype.e2=function(){var baseVec=new x3dom.fields.SFVec3f(this._02,this._12,this._22);return baseVec.normalize();};x3dom.fields.SFMatrix4f.prototype.e3=function(){return new x3dom.fields.SFVec3f(this._03,this._13,this._23);};x3dom.fields.SFMatrix4f.copy=function(that){return new x3dom.fields.SFMatrix4f(that._00,that._01,that._02,that._03,that._10,that._11,that._12,that._13,that._20,that._21,that._22,that._23,that._30,that._31,that._32,that._33);};x3dom.fields.SFMatrix4f.prototype.copy=function(){return x3dom.fields.SFMatrix4f.copy(this);};x3dom.fields.SFMatrix4f.identity=function(){return new x3dom.fields.SFMatrix4f(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1);};x3dom.fields.SFMatrix4f.zeroMatrix=function(){return new x3dom.fields.SFMatrix4f(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);};x3dom.fields.SFMatrix4f.translation=function(vec){return new x3dom.fields.SFMatrix4f(1,0,0,vec.x,0,1,0,vec.y,0,0,1,vec.z,0,0,0,1);};x3dom.fields.SFMatrix4f.rotationX=function(a){var c=Math.cos(a);var s=Math.sin(a);return new x3dom.fields.SFMatrix4f(1,0,0,0,0,c,-s,0,0,s,c,0,0,0,0,1);};x3dom.fields.SFMatrix4f.rotationY=function(a){var c=Math.cos(a);var s=Math.sin(a);return new x3dom.fields.SFMatrix4f(c,0,s,0,0,1,0,0,-s,0,c,0,0,0,0,1);};x3dom.fields.SFMatrix4f.rotationZ=function(a){var c=Math.cos(a);var s=Math.sin(a);return new x3dom.fields.SFMatrix4f(c,-s,0,0,s,c,0,0,0,0,1,0,0,0,0,1);};x3dom.fields.SFMatrix4f.scale=function(vec){return new x3dom.fields.SFMatrix4f(vec.x,0,0,0,0,vec.y,0,0,0,0,vec.z,0,0,0,0,1);};x3dom.fields.SFMatrix4f.lookAt=function(from,at,up)
{var view=from.subtract(at).normalize();var right=up.normalize().cross(view).normalize();if(right.dot(right)<x3dom.fields.Eps){x3dom.debug.logWarning("View matrix is linearly dependent.");return x3dom.fields.SFMatrix4f.translation(from);}
var newUp=view.cross(right).normalize();var tmp=x3dom.fields.SFMatrix4f.identity();tmp.setValue(right,newUp,view,from);return tmp;};x3dom.fields.SFMatrix4f.perspectiveFrustum=function(left,right,bottom,top,near,far)
{return new x3dom.fields.SFMatrix4f(2*near/(right-left),0,(right+left)/(right-left),0,0,2*near/(top-bottom),(top+bottom)/(top-bottom),0,0,0,-(far+near)/(far-near),-2*far*near/(far-near),0,0,-1,0);};x3dom.fields.SFMatrix4f.perspective=function(fov,aspect,near,far)
{var f=1/Math.tan(fov/2);return new x3dom.fields.SFMatrix4f(f/aspect,0,0,0,0,f,0,0,0,0,(near+far)/(near-far),2*near*far/(near-far),0,0,-1,0);};x3dom.fields.SFMatrix4f.ortho=function(left,right,bottom,top,near,far,aspect)
{var rl=(right-left)/2;var tb=(top-bottom)/2;var fn=far-near;if(aspect===undefined)
aspect=1.0;if(aspect<(rl/tb))
tb=rl/aspect;else
rl=tb*aspect;left=-rl;right=rl;bottom=-tb;top=tb;rl*=2;tb*=2;return new x3dom.fields.SFMatrix4f(2/rl,0,0,-(right+left)/rl,0,2/tb,0,-(top+bottom)/tb,0,0,-2/fn,-(far+near)/fn,0,0,0,1);};x3dom.fields.SFMatrix4f.prototype.setTranslate=function(vec){this._03=vec.x;this._13=vec.y;this._23=vec.z;};x3dom.fields.SFMatrix4f.prototype.setScale=function(vec){this._00=vec.x;this._11=vec.y;this._22=vec.z;};x3dom.fields.SFMatrix4f.prototype.setRotate=function(quat){var xx=quat.x*quat.x;var xy=quat.x*quat.y;var xz=quat.x*quat.z;var yy=quat.y*quat.y;var yz=quat.y*quat.z;var zz=quat.z*quat.z;var wx=quat.w*quat.x;var wy=quat.w*quat.y;var wz=quat.w*quat.z;this._00=1-2*(yy+zz);this._01=2*(xy-wz);this._02=2*(xz+wy);this._10=2*(xy+wz);this._11=1-2*(xx+zz);this._12=2*(yz-wx);this._20=2*(xz-wy);this._21=2*(yz+wx);this._22=1-2*(xx+yy);};x3dom.fields.SFMatrix4f.parseRotation=function(str){var m=/^([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)\s*,?\s*([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)\s*,?\s*([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)\s*,?\s*([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)$/.exec(str);var x=+m[1],y=+m[2],z=+m[3],a=+m[4];var d=Math.sqrt(x*x+y*y+z*z);if(d===0){x=1;y=z=0;}else{x/=d;y/=d;z/=d;}
var c=Math.cos(a);var s=Math.sin(a);var t=1-c;return new x3dom.fields.SFMatrix4f(t*x*x+c,t*x*y+s*z,t*x*z-s*y,0,t*x*y-s*z,t*y*y+c,t*y*z+s*x,0,t*x*z+s*y,t*y*z-s*x,t*z*z+c,0,0,0,0,1).transpose();};x3dom.fields.SFMatrix4f.parse=function(str){var needTranspose=false;var val=/matrix.*\((.+)\)/;if(val.exec(str)){str=RegExp.$1;needTranspose=true;}
var arr=Array.map(str.split(/[,\s]+/),function(n){return+n;});if(arr.length>=16)
{if(!needTranspose){return new x3dom.fields.SFMatrix4f(arr[0],arr[1],arr[2],arr[3],arr[4],arr[5],arr[6],arr[7],arr[8],arr[9],arr[10],arr[11],arr[12],arr[13],arr[14],arr[15]);}
else{return new x3dom.fields.SFMatrix4f(arr[0],arr[4],arr[8],arr[12],arr[1],arr[5],arr[9],arr[13],arr[2],arr[6],arr[10],arr[14],arr[3],arr[7],arr[11],arr[15]);}}
else if(arr.length===6){return new x3dom.fields.SFMatrix4f(arr[0],arr[1],0,arr[4],arr[2],arr[3],0,arr[5],0,0,1,0,0,0,0,1);}
else{x3dom.debug.logWarning("SFMatrix4f - can't parse string: "+str);return x3dom.fields.SFMatrix4f.identity();}};x3dom.fields.SFMatrix4f.prototype.mult=function(that){return new x3dom.fields.SFMatrix4f(this._00*that._00+this._01*that._10+this._02*that._20+this._03*that._30,this._00*that._01+this._01*that._11+this._02*that._21+this._03*that._31,this._00*that._02+this._01*that._12+this._02*that._22+this._03*that._32,this._00*that._03+this._01*that._13+this._02*that._23+this._03*that._33,this._10*that._00+this._11*that._10+this._12*that._20+this._13*that._30,this._10*that._01+this._11*that._11+this._12*that._21+this._13*that._31,this._10*that._02+this._11*that._12+this._12*that._22+this._13*that._32,this._10*that._03+this._11*that._13+this._12*that._23+this._13*that._33,this._20*that._00+this._21*that._10+this._22*that._20+this._23*that._30,this._20*that._01+this._21*that._11+this._22*that._21+this._23*that._31,this._20*that._02+this._21*that._12+this._22*that._22+this._23*that._32,this._20*that._03+this._21*that._13+this._22*that._23+this._23*that._33,this._30*that._00+this._31*that._10+this._32*that._20+this._33*that._30,this._30*that._01+this._31*that._11+this._32*that._21+this._33*that._31,this._30*that._02+this._31*that._12+this._32*that._22+this._33*that._32,this._30*that._03+this._31*that._13+this._32*that._23+this._33*that._33);};x3dom.fields.SFMatrix4f.prototype.multMatrixPnt=function(vec){return new x3dom.fields.SFVec3f(this._00*vec.x+this._01*vec.y+this._02*vec.z+this._03,this._10*vec.x+this._11*vec.y+this._12*vec.z+this._13,this._20*vec.x+this._21*vec.y+this._22*vec.z+this._23);};x3dom.fields.SFMatrix4f.prototype.multMatrixVec=function(vec){return new x3dom.fields.SFVec3f(this._00*vec.x+this._01*vec.y+this._02*vec.z,this._10*vec.x+this._11*vec.y+this._12*vec.z,this._20*vec.x+this._21*vec.y+this._22*vec.z);};x3dom.fields.SFMatrix4f.prototype.multFullMatrixPnt=function(vec){var w=this._30*vec.x+this._31*vec.y+this._32*vec.z+this._33;if(w){w=1.0/w;}
return new x3dom.fields.SFVec3f((this._00*vec.x+this._01*vec.y+this._02*vec.z+this._03)*w,(this._10*vec.x+this._11*vec.y+this._12*vec.z+this._13)*w,(this._20*vec.x+this._21*vec.y+this._22*vec.z+this._23)*w);};x3dom.fields.SFMatrix4f.prototype.multMatrixPlane=function(plane){var normal=new x3dom.fields.SFVec3f(plane.x,plane.y,plane.z);var memberPnt=normal.multiply(-plane.w);memberPnt=this.multMatrixPnt(memberPnt);var invTranspose=this.inverse().transpose();normal=invTranspose.multMatrixVec(normal);var d=-normal.dot(memberPnt);return new x3dom.fields.SFVec4f(normal.x,normal.y,normal.z,d);};x3dom.fields.SFMatrix4f.prototype.transpose=function(){return new x3dom.fields.SFMatrix4f(this._00,this._10,this._20,this._30,this._01,this._11,this._21,this._31,this._02,this._12,this._22,this._32,this._03,this._13,this._23,this._33);};x3dom.fields.SFMatrix4f.prototype.negate=function(){return new x3dom.fields.SFMatrix4f(-this._00,-this._01,-this._02,-this._03,-this._10,-this._11,-this._12,-this._13,-this._20,-this._21,-this._22,-this._23,-this._30,-this._31,-this._32,-this._33);};x3dom.fields.SFMatrix4f.prototype.multiply=function(s){return new x3dom.fields.SFMatrix4f(s*this._00,s*this._01,s*this._02,s*this._03,s*this._10,s*this._11,s*this._12,s*this._13,s*this._20,s*this._21,s*this._22,s*this._23,s*this._30,s*this._31,s*this._32,s*this._33);};x3dom.fields.SFMatrix4f.prototype.add=function(that){return new x3dom.fields.SFMatrix4f(this._00+that._00,this._01+that._01,this._02+that._02,this._03+that._03,this._10+that._10,this._11+that._11,this._12+that._12,this._13+that._13,this._20+that._20,this._21+that._21,this._22+that._22,this._23+that._23,this._30+that._30,this._31+that._31,this._32+that._32,this._33+that._33);};x3dom.fields.SFMatrix4f.prototype.addScaled=function(that,s){return new x3dom.fields.SFMatrix4f(this._00+s*that._00,this._01+s*that._01,this._02+s*that._02,this._03+s*that._03,this._10+s*that._10,this._11+s*that._11,this._12+s*that._12,this._13+s*that._13,this._20+s*that._20,this._21+s*that._21,this._22+s*that._22,this._23+s*that._23,this._30+s*that._30,this._31+s*that._31,this._32+s*that._32,this._33+s*that._33);};x3dom.fields.SFMatrix4f.prototype.setValues=function(that){this._00=that._00;this._01=that._01;this._02=that._02;this._03=that._03;this._10=that._10;this._11=that._11;this._12=that._12;this._13=that._13;this._20=that._20;this._21=that._21;this._22=that._22;this._23=that._23;this._30=that._30;this._31=that._31;this._32=that._32;this._33=that._33;};x3dom.fields.SFMatrix4f.prototype.setValue=function(v1,v2,v3,v4){this._00=v1.x;this._01=v2.x;this._02=v3.x;this._10=v1.y;this._11=v2.y;this._12=v3.y;this._20=v1.z;this._21=v2.z;this._22=v3.z;this._30=0;this._31=0;this._32=0;if(arguments.length>3){this._03=v4.x;this._13=v4.y;this._23=v4.z;this._33=1;}};x3dom.fields.SFMatrix4f.prototype.setFromArray=function(a){this._00=a[0];this._01=a[4];this._02=a[8];this._03=a[12];this._10=a[1];this._11=a[5];this._12=a[9];this._13=a[13];this._20=a[2];this._21=a[6];this._22=a[10];this._23=a[14];this._30=a[3];this._31=a[7];this._32=a[11];this._33=a[15];};x3dom.fields.SFMatrix4f.prototype.toGL=function(){return[this._00,this._10,this._20,this._30,this._01,this._11,this._21,this._31,this._02,this._12,this._22,this._32,this._03,this._13,this._23,this._33];};x3dom.fields.SFMatrix4f.prototype.at=function(i,j){var field="_"+i+j;return this[field];};x3dom.fields.SFMatrix4f.prototype.sqrt=function(){var Y=x3dom.fields.SFMatrix4f.identity();var result=x3dom.fields.SFMatrix4f.copy(this);for(var i=0;i<6;i++)
{var iX=result.inverse();var iY=(i==0)?x3dom.fields.SFMatrix4f.identity():Y.inverse();var rd=result.det(),yd=Y.det();var g=Math.abs(Math.pow(rd*yd,-0.125));var ig=1.0/g;result=result.multiply(g);result=result.addScaled(iY,ig);result=result.multiply(0.5);Y=Y.multiply(g);Y=Y.addScaled(iX,ig);Y=Y.multiply(0.5);}
return result;};x3dom.fields.SFMatrix4f.prototype.normInfinity=function(){var t=0,m=0;if((t=Math.abs(this._00))>m){m=t;}
if((t=Math.abs(this._01))>m){m=t;}
if((t=Math.abs(this._02))>m){m=t;}
if((t=Math.abs(this._03))>m){m=t;}
if((t=Math.abs(this._10))>m){m=t;}
if((t=Math.abs(this._11))>m){m=t;}
if((t=Math.abs(this._12))>m){m=t;}
if((t=Math.abs(this._13))>m){m=t;}
if((t=Math.abs(this._20))>m){m=t;}
if((t=Math.abs(this._21))>m){m=t;}
if((t=Math.abs(this._22))>m){m=t;}
if((t=Math.abs(this._23))>m){m=t;}
if((t=Math.abs(this._30))>m){m=t;}
if((t=Math.abs(this._31))>m){m=t;}
if((t=Math.abs(this._32))>m){m=t;}
if((t=Math.abs(this._33))>m){m=t;}
return m;};x3dom.fields.SFMatrix4f.prototype.norm1_3x3=function(){var max=Math.abs(this._00)+
Math.abs(this._10)+
Math.abs(this._20);var t=0;if((t=Math.abs(this._01)+
Math.abs(this._11)+
Math.abs(this._21))>max){max=t;}
if((t=Math.abs(this._02)+
Math.abs(this._12)+
Math.abs(this._22))>max){max=t;}
return max;};x3dom.fields.SFMatrix4f.prototype.normInf_3x3=function(){var max=Math.abs(this._00)+
Math.abs(this._01)+
Math.abs(this._02);var t=0;if((t=Math.abs(this._10)+
Math.abs(this._11)+
Math.abs(this._12))>max){max=t;}
if((t=Math.abs(this._20)+
Math.abs(this._21)+
Math.abs(this._22))>max){max=t;}
return max;};x3dom.fields.SFMatrix4f.prototype.adjointT_3x3=function(){var result=x3dom.fields.SFMatrix4f.identity();result._00=this._11*this._22-this._12*this._21;result._01=this._12*this._20-this._10*this._22;result._02=this._10*this._21-this._11*this._20;result._10=this._21*this._02-this._22*this._01;result._11=this._22*this._00-this._20*this._02;result._12=this._20*this._01-this._21*this._00;result._20=this._01*this._12-this._02*this._11;result._21=this._02*this._10-this._00*this._12;result._22=this._00*this._11-this._01*this._10;return result;};x3dom.fields.SFMatrix4f.prototype.equals=function(that){var eps=0.000000000001;return Math.abs(this._00-that._00)<eps&&Math.abs(this._01-that._01)<eps&&Math.abs(this._02-that._02)<eps&&Math.abs(this._03-that._03)<eps&&Math.abs(this._10-that._10)<eps&&Math.abs(this._11-that._11)<eps&&Math.abs(this._12-that._12)<eps&&Math.abs(this._13-that._13)<eps&&Math.abs(this._20-that._20)<eps&&Math.abs(this._21-that._21)<eps&&Math.abs(this._22-that._22)<eps&&Math.abs(this._23-that._23)<eps&&Math.abs(this._30-that._30)<eps&&Math.abs(this._31-that._31)<eps&&Math.abs(this._32-that._32)<eps&&Math.abs(this._33-that._33)<eps;};x3dom.fields.SFMatrix4f.prototype.getTransform=function(translation,rotation,scaleFactor,scaleOrientation,center)
{var m=null;if(arguments.length>4){m=x3dom.fields.SFMatrix4f.translation(center.negate());m=m.mult(this);var c=x3dom.fields.SFMatrix4f.translation(center);m=m.mult(c);}
else{m=x3dom.fields.SFMatrix4f.copy(this);}
var flip=m.decompose(translation,rotation,scaleFactor,scaleOrientation);scaleFactor.setValues(scaleFactor.multiply(flip));};x3dom.fields.SFMatrix4f.prototype.decompose=function(t,r,s,so)
{var A=x3dom.fields.SFMatrix4f.copy(this);var Q=x3dom.fields.SFMatrix4f.identity(),S=x3dom.fields.SFMatrix4f.identity(),SO=x3dom.fields.SFMatrix4f.identity();t.x=A._03;t.y=A._13;t.z=A._23;A._03=0.0;A._13=0.0;A._23=0.0;A._30=0.0;A._31=0.0;A._32=0.0;var det=A.polarDecompose(Q,S);var f=1.0;if(det<0.0){Q=Q.negate();f=-1.0;}
r.setValue(Q);S.spectralDecompose(SO,s);so.setValue(SO);return f;};x3dom.fields.SFMatrix4f.prototype.polarDecompose=function(Q,S)
{var TOL=0.000000000001;var Mk=this.transpose();var Ek=x3dom.fields.SFMatrix4f.identity();var Mk_one=Mk.norm1_3x3();var Mk_inf=Mk.normInf_3x3();var MkAdjT;var MkAdjT_one,MkAdjT_inf;var Ek_one,Mk_det;do
{MkAdjT=Mk.adjointT_3x3();Mk_det=Mk._00*MkAdjT._00+
Mk._01*MkAdjT._01+
Mk._02*MkAdjT._02;if(Mk_det==0.0)
{x3dom.debug.logWarning("polarDecompose: Mk_det == 0.0");break;}
MkAdjT_one=MkAdjT.norm1_3x3();MkAdjT_inf=MkAdjT.normInf_3x3();var gamma=Math.sqrt(Math.sqrt((MkAdjT_one*MkAdjT_inf)/(Mk_one*Mk_inf))/Math.abs(Mk_det));var g1=0.5*gamma;var g2=0.5/(gamma*Mk_det);Ek.setValues(Mk);Mk=Mk.multiply(g1);Mk=Mk.addScaled(MkAdjT,g2);Ek=Ek.addScaled(Mk,-1.0);Ek_one=Ek.norm1_3x3();Mk_one=Mk.norm1_3x3();Mk_inf=Mk.normInf_3x3();}while(Ek_one>(Mk_one*TOL));Q.setValues(Mk.transpose());S.setValues(Mk.mult(this));for(var i=0;i<3;++i)
{for(var j=i;j<3;++j)
{S['_'+j+i]=0.5*(S['_'+j+i]+S['_'+i+j]);S['_'+i+j]=0.5*(S['_'+j+i]+S['_'+i+j]);}}
return Mk_det;};x3dom.fields.SFMatrix4f.prototype.spectralDecompose=function(SO,k)
{var next=[1,2,0];var maxIterations=20;var diag=[this._00,this._11,this._22];var offDiag=[this._12,this._20,this._01];for(var iter=0;iter<maxIterations;++iter)
{var sm=Math.abs(offDiag[0])+Math.abs(offDiag[1])+Math.abs(offDiag[2]);if(sm==0){break;}
for(var i=2;i>=0;--i)
{var p=next[i];var q=next[p];var absOffDiag=Math.abs(offDiag[i]);var g=100.0*absOffDiag;if(absOffDiag>0.0)
{var t=0,h=diag[q]-diag[p];var absh=Math.abs(h);if(absh+g==absh)
{t=offDiag[i]/h;}
else
{var theta=0.5*h/offDiag[i];t=1.0/(Math.abs(theta)+Math.sqrt(theta*theta+1.0));t=theta<0.0?-t:t;}
var c=1.0/Math.sqrt(t*t+1.0);var s=t*c;var tau=s/(c+1.0);var ta=t*offDiag[i];offDiag[i]=0.0;diag[p]-=ta;diag[q]+=ta;var offDiagq=offDiag[q];offDiag[q]-=s*(offDiag[p]+tau*offDiagq);offDiag[p]+=s*(offDiagq-tau*offDiag[p]);for(var j=2;j>=0;--j)
{var a=SO['_'+j+p];var b=SO['_'+j+q];SO['_'+j+p]-=s*(b+tau*a);SO['_'+j+q]+=s*(a-tau*b);}}}}
k.x=diag[0];k.y=diag[1];k.z=diag[2];};x3dom.fields.SFMatrix4f.prototype.log=function(){var maxiter=12;var eps=1e-12;var A=x3dom.fields.SFMatrix4f.copy(this),Z=x3dom.fields.SFMatrix4f.copy(this);Z._00-=1;Z._11-=1;Z._22-=1;Z._33-=1;var k=0;while(Z.normInfinity()>0.5)
{A=A.sqrt();Z.setValues(A);Z._00-=1;Z._11-=1;Z._22-=1;Z._33-=1;k++;}
A._00-=1;A._11-=1;A._22-=1;A._33-=1;A=A.negate();Z.setValues(A);var result=x3dom.fields.SFMatrix4f.copy(A);var i=1;while(Z.normInfinity()>eps&&i<maxiter)
{Z=Z.mult(A);i++;result=result.addScaled(Z,1.0/i);}
return result.multiply(-(1<<k));};x3dom.fields.SFMatrix4f.prototype.exp=function(){var q=6;var A=x3dom.fields.SFMatrix4f.copy(this),D=x3dom.fields.SFMatrix4f.identity(),N=x3dom.fields.SFMatrix4f.identity(),result=x3dom.fields.SFMatrix4f.identity();var k=0,c=1.0;var j=1.0+parseInt(Math.log(A.normInfinity()/0.693));if(j<0){j=0;}
A=A.multiply(1.0/(1<<j));for(k=1;k<=q;k++)
{c*=(q-k+1)/(k*(2*q-k+1));result=A.mult(result);N=N.addScaled(result,c);if(k%2){D=D.addScaled(result,-c);}
else{D=D.addScaled(result,c);}}
result=D.inverse().mult(N);for(k=0;k<j;k++)
{result=result.mult(result);}
return result;};x3dom.fields.SFMatrix4f.prototype.det3=function(a1,a2,a3,b1,b2,b3,c1,c2,c3){return((a1*b2*c3)+(a2*b3*c1)+(a3*b1*c2)-
(a1*b3*c2)-(a2*b1*c3)-(a3*b2*c1));};x3dom.fields.SFMatrix4f.prototype.det=function(){var a1=this._00;var b1=this._10;var c1=this._20;var d1=this._30;var a2=this._01;var b2=this._11;var c2=this._21;var d2=this._31;var a3=this._02;var b3=this._12;var c3=this._22;var d3=this._32;var a4=this._03;var b4=this._13;var c4=this._23;var d4=this._33;return(a1*this.det3(b2,b3,b4,c2,c3,c4,d2,d3,d4)-
b1*this.det3(a2,a3,a4,c2,c3,c4,d2,d3,d4)+
c1*this.det3(a2,a3,a4,b2,b3,b4,d2,d3,d4)-
d1*this.det3(a2,a3,a4,b2,b3,b4,c2,c3,c4));};x3dom.fields.SFMatrix4f.prototype.inverse=function(){var a1=this._00;var b1=this._10;var c1=this._20;var d1=this._30;var a2=this._01;var b2=this._11;var c2=this._21;var d2=this._31;var a3=this._02;var b3=this._12;var c3=this._22;var d3=this._32;var a4=this._03;var b4=this._13;var c4=this._23;var d4=this._33;var rDet=this.det();if(rDet==0)
{x3dom.debug.logWarning("Invert matrix: singular matrix, no inverse!");return x3dom.fields.SFMatrix4f.identity();}
rDet=1.0/rDet;return new x3dom.fields.SFMatrix4f(+this.det3(b2,b3,b4,c2,c3,c4,d2,d3,d4)*rDet,-this.det3(a2,a3,a4,c2,c3,c4,d2,d3,d4)*rDet,+this.det3(a2,a3,a4,b2,b3,b4,d2,d3,d4)*rDet,-this.det3(a2,a3,a4,b2,b3,b4,c2,c3,c4)*rDet,-this.det3(b1,b3,b4,c1,c3,c4,d1,d3,d4)*rDet,+this.det3(a1,a3,a4,c1,c3,c4,d1,d3,d4)*rDet,-this.det3(a1,a3,a4,b1,b3,b4,d1,d3,d4)*rDet,+this.det3(a1,a3,a4,b1,b3,b4,c1,c3,c4)*rDet,+this.det3(b1,b2,b4,c1,c2,c4,d1,d2,d4)*rDet,-this.det3(a1,a2,a4,c1,c2,c4,d1,d2,d4)*rDet,+this.det3(a1,a2,a4,b1,b2,b4,d1,d2,d4)*rDet,-this.det3(a1,a2,a4,b1,b2,b4,c1,c2,c4)*rDet,-this.det3(b1,b2,b3,c1,c2,c3,d1,d2,d3)*rDet,+this.det3(a1,a2,a3,c1,c2,c3,d1,d2,d3)*rDet,-this.det3(a1,a2,a3,b1,b2,b3,d1,d2,d3)*rDet,+this.det3(a1,a2,a3,b1,b2,b3,c1,c2,c3)*rDet);};x3dom.fields.SFMatrix4f.prototype.getEulerAngles=function(){var theta_1,theta_2,theta;var phi_1,phi_2,phi;var psi_1,psi_2,psi;var cos_theta_1,cos_theta_2;if(Math.abs((Math.abs(this._20)-1.0))>0.0001){theta_1=-Math.asin(this._20);theta_2=Math.PI-theta_1;cos_theta_1=Math.cos(theta_1);cos_theta_2=Math.cos(theta_2);psi_1=Math.atan2(this._21/cos_theta_1,this._22/cos_theta_1);psi_2=Math.atan2(this._21/cos_theta_2,this._22/cos_theta_2);phi_1=Math.atan2(this._10/cos_theta_1,this._00/cos_theta_1);phi_2=Math.atan2(this._10/cos_theta_2,this._00/cos_theta_2);return[psi_1,theta_1,phi_1,psi_2,theta_2,phi_2];}
else{phi=0;if(this._20==-1.0){theta=Math.PI/2.0;psi=phi+Math.atan2(this._01,this._02);}
else{theta=-(Math.PI/2.0);psi=-phi+Math.atan2(-this._01,-this._02);}
return[psi,theta,phi,psi,theta,phi];}};x3dom.fields.SFMatrix4f.prototype.toString=function(){return'\n'+
this._00.toFixed(6)+', '+this._01.toFixed(6)+', '+
this._02.toFixed(6)+', '+this._03.toFixed(6)+', \n'+
this._10.toFixed(6)+', '+this._11.toFixed(6)+', '+
this._12.toFixed(6)+', '+this._13.toFixed(6)+', \n'+
this._20.toFixed(6)+', '+this._21.toFixed(6)+', '+
this._22.toFixed(6)+', '+this._23.toFixed(6)+', \n'+
this._30.toFixed(6)+', '+this._31.toFixed(6)+', '+
this._32.toFixed(6)+', '+this._33.toFixed(6);};x3dom.fields.SFMatrix4f.prototype.setValueByStr=function(str){var needTranspose=false;var val=/matrix.*\((.+)\)/;if(val.exec(str)){str=RegExp.$1;needTranspose=true;}
var arr=Array.map(str.split(/[,\s]+/),function(n){return+n;});if(arr.length>=16)
{if(!needTranspose){this._00=arr[0];this._01=arr[1];this._02=arr[2];this._03=arr[3];this._10=arr[4];this._11=arr[5];this._12=arr[6];this._13=arr[7];this._20=arr[8];this._21=arr[9];this._22=arr[10];this._23=arr[11];this._30=arr[12];this._31=arr[13];this._32=arr[14];this._33=arr[15];}
else{this._00=arr[0];this._01=arr[4];this._02=arr[8];this._03=arr[12];this._10=arr[1];this._11=arr[5];this._12=arr[9];this._13=arr[13];this._20=arr[2];this._21=arr[6];this._22=arr[10];this._23=arr[14];this._30=arr[3];this._31=arr[7];this._32=arr[11];this._33=arr[15];}}
else if(arr.length===6){this._00=arr[0];this._01=arr[1];this._02=0;this._03=arr[4];this._10=arr[2];this._11=arr[3];this._12=0;this._13=arr[5];this._20=0;this._21=0;this._22=1;this._23=0;this._30=0;this._31=0;this._32=0;this._33=1;}
else{x3dom.debug.logWarning("SFMatrix4f - can't parse string: "+str);}
return this;};x3dom.fields.SFVec2f=function(x,y){if(arguments.length===0){this.x=0;this.y=0;}
else{this.x=x;this.y=y;}};x3dom.fields.SFVec2f.copy=function(v){return new x3dom.fields.SFVec2f(v.x,v.y);};x3dom.fields.SFVec2f.parse=function(str){var m=/^\s*([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)\s*,?\s*([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)\s*$/.exec(str);return new x3dom.fields.SFVec2f(+m[1],+m[2]);};x3dom.fields.SFVec2f.prototype.copy=function(){return x3dom.fields.SFVec2f.copy(this);};x3dom.fields.SFVec2f.prototype.setValues=function(that){this.x=that.x;this.y=that.y;};x3dom.fields.SFVec2f.prototype.at=function(i){switch(i){case 0:return this.x;case 1:return this.y;default:return this.x;}};x3dom.fields.SFVec2f.prototype.add=function(that){return new x3dom.fields.SFVec2f(this.x+that.x,this.y+that.y);};x3dom.fields.SFVec2f.prototype.subtract=function(that){return new x3dom.fields.SFVec2f(this.x-that.x,this.y-that.y);};x3dom.fields.SFVec2f.prototype.negate=function(){return new x3dom.fields.SFVec2f(-this.x,-this.y);};x3dom.fields.SFVec2f.prototype.dot=function(that){return this.x*that.x+this.y*that.y;};x3dom.fields.SFVec2f.prototype.reflect=function(n){var d2=this.dot(n)*2;return new x3dom.fields.SFVec2f(this.x-d2*n.x,this.y-d2*n.y);};x3dom.fields.SFVec2f.prototype.normalize=function(){var n=this.length();if(n){n=1.0/n;}
return new x3dom.fields.SFVec2f(this.x*n,this.y*n);};x3dom.fields.SFVec2f.prototype.multComponents=function(that){return new x3dom.fields.SFVec2f(this.x*that.x,this.y*that.y);};x3dom.fields.SFVec2f.prototype.multiply=function(n){return new x3dom.fields.SFVec2f(this.x*n,this.y*n);};x3dom.fields.SFVec2f.prototype.divideComponents=function(that){return new x3dom.fields.SFVec2f(this.x/that.x,this.y/that.y);};x3dom.fields.SFVec2f.prototype.divide=function(n){var denom=n?(1.0/n):1.0;return new x3dom.fields.SFVec2f(this.x*denom,this.y*denom);};x3dom.fields.SFVec2f.prototype.equals=function(that,eps){return Math.abs(this.x-that.x)<eps&&Math.abs(this.y-that.y)<eps;};x3dom.fields.SFVec2f.prototype.length=function(){return Math.sqrt((this.x*this.x)+(this.y*this.y));};x3dom.fields.SFVec2f.prototype.toGL=function(){return[this.x,this.y];};x3dom.fields.SFVec2f.prototype.toString=function(){return this.x+" "+this.y;};x3dom.fields.SFVec2f.prototype.setValueByStr=function(str){var m=/^\s*([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)\s*,?\s*([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)\s*$/.exec(str);this.x=+m[1];this.y=+m[2];return this;};x3dom.fields.SFVec3f=function(x,y,z){if(arguments.length===0){this.x=0;this.y=0;this.z=0;}
else{this.x=x;this.y=y;this.z=z;}};x3dom.fields.SFVec3f.NullVector=new x3dom.fields.SFVec3f(0,0,0);x3dom.fields.SFVec3f.OneVector=new x3dom.fields.SFVec3f(1,1,1);x3dom.fields.SFVec3f.copy=function(v){return new x3dom.fields.SFVec3f(v.x,v.y,v.z);};x3dom.fields.SFVec3f.MIN=function(){return new x3dom.fields.SFVec3f(-Number.MAX_VALUE,-Number.MAX_VALUE,-Number.MAX_VALUE);};x3dom.fields.SFVec3f.MAX=function(){return new x3dom.fields.SFVec3f(Number.MAX_VALUE,Number.MAX_VALUE,Number.MAX_VALUE);};x3dom.fields.SFVec3f.parse=function(str){try{var m=/^\s*([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)\s*,?\s*([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)\s*,?\s*([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)\s*$/.exec(str);return new x3dom.fields.SFVec3f(+m[1],+m[2],+m[3]);}
catch(e){var c=x3dom.fields.SFColor.colorParse(str);return new x3dom.fields.SFVec3f(c.r,c.g,c.b);}};x3dom.fields.SFVec3f.prototype.copy=function(){return x3dom.fields.SFVec3f.copy(this);};x3dom.fields.SFVec3f.prototype.setValues=function(that){this.x=that.x;this.y=that.y;this.z=that.z;};x3dom.fields.SFVec3f.prototype.at=function(i){switch(i){case 0:return this.x;case 1:return this.y;case 2:return this.z;default:return this.x;}};x3dom.fields.SFVec3f.prototype.add=function(that){return new x3dom.fields.SFVec3f(this.x+that.x,this.y+that.y,this.z+that.z);};x3dom.fields.SFVec3f.prototype.addScaled=function(that,s){return new x3dom.fields.SFVec3f(this.x+s*that.x,this.y+s*that.y,this.z+s*that.z);};x3dom.fields.SFVec3f.prototype.subtract=function(that){return new x3dom.fields.SFVec3f(this.x-that.x,this.y-that.y,this.z-that.z);};x3dom.fields.SFVec3f.prototype.negate=function(){return new x3dom.fields.SFVec3f(-this.x,-this.y,-this.z);};x3dom.fields.SFVec3f.prototype.dot=function(that){return(this.x*that.x+this.y*that.y+this.z*that.z);};x3dom.fields.SFVec3f.prototype.cross=function(that){return new x3dom.fields.SFVec3f(this.y*that.z-this.z*that.y,this.z*that.x-this.x*that.z,this.x*that.y-this.y*that.x);};x3dom.fields.SFVec3f.prototype.reflect=function(n){var d2=this.dot(n)*2;return new x3dom.fields.SFVec3f(this.x-d2*n.x,this.y-d2*n.y,this.z-d2*n.z);};x3dom.fields.SFVec3f.prototype.length=function(){return Math.sqrt((this.x*this.x)+(this.y*this.y)+(this.z*this.z));};x3dom.fields.SFVec3f.prototype.normalize=function(){var n=this.length();if(n){n=1.0/n;}
return new x3dom.fields.SFVec3f(this.x*n,this.y*n,this.z*n);};x3dom.fields.SFVec3f.prototype.multComponents=function(that){return new x3dom.fields.SFVec3f(this.x*that.x,this.y*that.y,this.z*that.z);};x3dom.fields.SFVec3f.prototype.multiply=function(n){return new x3dom.fields.SFVec3f(this.x*n,this.y*n,this.z*n);};x3dom.fields.SFVec3f.prototype.divide=function(n){var denom=n?(1.0/n):1.0;return new x3dom.fields.SFVec3f(this.x*denom,this.y*denom,this.z*denom);};x3dom.fields.SFVec3f.prototype.equals=function(that,eps){return Math.abs(this.x-that.x)<eps&&Math.abs(this.y-that.y)<eps&&Math.abs(this.z-that.z)<eps;};x3dom.fields.SFVec3f.prototype.toGL=function(){return[this.x,this.y,this.z];};x3dom.fields.SFVec3f.prototype.toString=function(){return this.x+" "+this.y+" "+this.z;};x3dom.fields.SFVec3f.prototype.setValueByStr=function(str){try{var m=/^\s*([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)\s*,?\s*([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)\s*,?\s*([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)\s*$/.exec(str);this.x=+m[1];this.y=+m[2];this.z=+m[3];}
catch(e){var c=x3dom.fields.SFColor.colorParse(str);this.x=c.r;this.y=c.g;this.z=c.b;}
return this;};x3dom.fields.SFVec4f=function(x,y,z,w){if(arguments.length===0){this.x=0;this.y=0;this.z=0;this.w=0;}
else{this.x=x;this.y=y;this.z=z;this.w=w;}};x3dom.fields.SFVec4f.copy=function(v){return new x3dom.fields.SFVec4f(v.x,v.y,v.z,v.w);};x3dom.fields.SFVec4f.prototype.copy=function(){return x3dom.fields.SFVec4f(this);};x3dom.fields.SFVec4f.parse=function(str){var m=/^\s*([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)\s*,?\s*([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)\s*,?\s*([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)\s*,?\s*([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)\s*$/.exec(str);return new x3dom.fields.SFVec4f(+m[1],+m[2],+m[3],+m[4]);};x3dom.fields.SFVec4f.prototype.setValueByStr=function(str){var m=/^\s*([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)\s*,?\s*([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)\s*,?\s*([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)\s*,?\s*([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)\s*$/.exec(str);this.x=+m[1];this.y=+m[2];this.z=+m[3];this.w=+m[4];return this;};x3dom.fields.SFVec4f.prototype.toGL=function(){return[this.x,this.y,this.z,this.w];};x3dom.fields.SFVec4f.prototype.toString=function(){return this.x+" "+this.y+" "+this.z+" "+this.w;};x3dom.fields.Quaternion=function(x,y,z,w){if(arguments.length===0){this.x=0;this.y=0;this.z=0;this.w=1;}
else{this.x=x;this.y=y;this.z=z;this.w=w;}};x3dom.fields.Quaternion.copy=function(v){return new x3dom.fields.Quaternion(v.x,v.y,v.z,v.w);};x3dom.fields.Quaternion.prototype.multiply=function(that){return new x3dom.fields.Quaternion(this.w*that.x+this.x*that.w+this.y*that.z-this.z*that.y,this.w*that.y+this.y*that.w+this.z*that.x-this.x*that.z,this.w*that.z+this.z*that.w+this.x*that.y-this.y*that.x,this.w*that.w-this.x*that.x-this.y*that.y-this.z*that.z);};x3dom.fields.Quaternion.parseAxisAngle=function(str){var m=/^\s*([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)\s*,?\s*([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)\s*,?\s*([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)\s*,?\s*([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)\s*$/.exec(str);return x3dom.fields.Quaternion.axisAngle(new x3dom.fields.SFVec3f(+m[1],+m[2],+m[3]),+m[4]);};x3dom.fields.Quaternion.axisAngle=function(axis,a){var t=axis.length();if(t>x3dom.fields.Eps)
{var s=Math.sin(a/2)/t;var c=Math.cos(a/2);return new x3dom.fields.Quaternion(axis.x*s,axis.y*s,axis.z*s,c);}
else
{return new x3dom.fields.Quaternion(0,0,0,1);}};x3dom.fields.Quaternion.prototype.copy=function(){return x3dom.fields.Quaternion.copy(this);};x3dom.fields.Quaternion.prototype.toMatrix=function(){var xx=this.x*this.x;var xy=this.x*this.y;var xz=this.x*this.z;var yy=this.y*this.y;var yz=this.y*this.z;var zz=this.z*this.z;var wx=this.w*this.x;var wy=this.w*this.y;var wz=this.w*this.z;return new x3dom.fields.SFMatrix4f(1-2*(yy+zz),2*(xy-wz),2*(xz+wy),0,2*(xy+wz),1-2*(xx+zz),2*(yz-wx),0,2*(xz-wy),2*(yz+wx),1-2*(xx+yy),0,0,0,0,1);};x3dom.fields.Quaternion.prototype.toAxisAngle=function()
{var x=0,y=0,z=0;var s=0,a=0;var that=this;if(this.w>1)
{that=x3dom.fields.Quaternion.normalize(this);}
a=2*Math.acos(that.w);s=Math.sqrt(1-that.w*that.w);if(s==0)
{x=that.x;y=that.y;z=that.z;}
else
{x=that.x/s;y=that.y/s;z=that.z/s;}
return[new x3dom.fields.SFVec3f(x,y,z),a];};x3dom.fields.Quaternion.prototype.angle=function()
{return 2*Math.acos(this.w);};x3dom.fields.Quaternion.prototype.setValue=function(matrix)
{var tr,s=1;var qt=[0,0,0];var i=0,j=0,k=0;var nxt=[1,2,0];tr=matrix._00+matrix._11+matrix._22;if(tr>0.0)
{s=Math.sqrt(tr+1.0);this.w=s*0.5;s=0.5/s;this.x=(matrix._21-matrix._12)*s;this.y=(matrix._02-matrix._20)*s;this.z=(matrix._10-matrix._01)*s;}
else
{if(matrix._11>matrix._00){i=1;}
else{i=0;}
if(matrix._22>matrix.at(i,i)){i=2;}
j=nxt[i];k=nxt[j];s=Math.sqrt(matrix.at(i,i)-(matrix.at(j,j)+matrix.at(k,k))+1.0);qt[i]=s*0.5;s=0.5/s;this.w=(matrix.at(k,j)-matrix.at(j,k))*s;qt[j]=(matrix.at(j,i)+matrix.at(i,j))*s;qt[k]=(matrix.at(k,i)+matrix.at(i,k))*s;this.x=qt[0];this.y=qt[1];this.z=qt[2];}
if(this.w>1.0||this.w<-1.0)
{var errThreshold=1+(x3dom.fields.Eps*100);if(this.w>errThreshold||this.w<-errThreshold)
{x3dom.debug.logInfo("MatToQuat: BUG: |quat[4]| ("+this.w+") >> 1.0 !");}
if(this.w>1.0){this.w=1.0;}
else{this.w=-1.0;}}};x3dom.fields.Quaternion.prototype.setFromEuler=function(alpha,beta,gamma){var sx=Math.sin(alpha*0.5);var cx=Math.cos(alpha*0.5);var sy=Math.sin(beta*0.5);var cy=Math.cos(beta*0.5);var sz=Math.sin(gamma*0.5);var cz=Math.cos(gamma*0.5);this.x=(sx*cy*cz)-(cx*sy*sz);this.y=(cx*sy*cz)+(sx*cy*sz);this.z=(cx*cy*sz)-(sx*sy*cz);this.w=(cx*cy*cz)+(sx*sy*sz);};x3dom.fields.Quaternion.prototype.dot=function(that){return this.x*that.x+this.y*that.y+this.z*that.z+this.w*that.w;};x3dom.fields.Quaternion.prototype.add=function(that){return new x3dom.fields.Quaternion(this.x+that.x,this.y+that.y,this.z+that.z,this.w+that.w);};x3dom.fields.Quaternion.prototype.subtract=function(that){return new x3dom.fields.Quaternion(this.x-that.x,this.y-that.y,this.z-that.z,this.w-that.w);};x3dom.fields.Quaternion.prototype.setValues=function(that){this.x=that.x;this.y=that.y;this.z=that.z;this.w=that.w;};x3dom.fields.Quaternion.prototype.equals=function(that,eps){return(this.dot(that)>=1.0-eps);};x3dom.fields.Quaternion.prototype.multScalar=function(s){return new x3dom.fields.Quaternion(this.x*s,this.y*s,this.z*s,this.w*s);};x3dom.fields.Quaternion.prototype.normalize=function(that){var d2=this.dot(that);var id=1.0;if(d2){id=1.0/Math.sqrt(d2);}
return new x3dom.fields.Quaternion(this.x*id,this.y*id,this.z*id,this.w*id);};x3dom.fields.Quaternion.prototype.negate=function(){return new x3dom.fields.Quaternion(-this.x,-this.y,-this.z,-this.w);};x3dom.fields.Quaternion.prototype.inverse=function(){return new x3dom.fields.Quaternion(-this.x,-this.y,-this.z,this.w);};x3dom.fields.Quaternion.prototype.slerp=function(that,t){var cosom=this.dot(that);var rot1;if(cosom<0.0)
{cosom=-cosom;rot1=that.negate();}
else
{rot1=new x3dom.fields.Quaternion(that.x,that.y,that.z,that.w);}
var scalerot0,scalerot1;if((1.0-cosom)>0.00001)
{var omega=Math.acos(cosom);var sinom=Math.sin(omega);scalerot0=Math.sin((1.0-t)*omega)/sinom;scalerot1=Math.sin(t*omega)/sinom;}
else
{scalerot0=1.0-t;scalerot1=t;}
return this.multScalar(scalerot0).add(rot1.multScalar(scalerot1));};x3dom.fields.Quaternion.rotateFromTo=function(fromVec,toVec){var from=fromVec.normalize();var to=toVec.normalize();var cost=from.dot(to);if(cost>0.99999)
{return new x3dom.fields.Quaternion(0,0,0,1);}
else if(cost<-0.99999)
{var cAxis=new x3dom.fields.SFVec3f(1,0,0);var tmp=from.cross(cAxis);if(tmp.length()<0.00001)
{cAxis.x=0;cAxis.y=1;cAxis.z=0;tmp=from.cross(cAxis);}
tmp=tmp.normalize();return x3dom.fields.Quaternion.axisAngle(tmp,Math.PI);}
var axis=fromVec.cross(toVec);axis=axis.normalize();var s=Math.sqrt(0.5*(1.0-cost));axis=axis.multiply(s);s=Math.sqrt(0.5*(1.0+cost));return new x3dom.fields.Quaternion(axis.x,axis.y,axis.z,s);};x3dom.fields.Quaternion.prototype.toGL=function(){var val=this.toAxisAngle();return[val[0].x,val[0].y,val[0].z,val[1]];};x3dom.fields.Quaternion.prototype.toString=function(){return this.x+" "+this.y+" "+this.z+", "+this.w;};x3dom.fields.Quaternion.prototype.setValueByStr=function(str){var m=/^\s*([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)\s*,?\s*([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)\s*,?\s*([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)\s*,?\s*([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)\s*$/.exec(str);var quat=x3dom.fields.Quaternion.axisAngle(new x3dom.fields.SFVec3f(+m[1],+m[2],+m[3]),+m[4]);this.x=quat.x;this.y=quat.y;this.z=quat.z;this.w=quat.w;return this;};x3dom.fields.SFColor=function(r,g,b){if(arguments.length===0){this.r=0;this.g=0;this.b=0;}
else{this.r=r;this.g=g;this.b=b;}};x3dom.fields.SFColor.parse=function(str){try{var m=/^\s*([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)\s*,?\s*([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)\s*,?\s*([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)\s*$/.exec(str);return new x3dom.fields.SFColor(+m[1],+m[2],+m[3]);}
catch(e){return x3dom.fields.SFColor.colorParse(str);}};x3dom.fields.SFColor.copy=function(that){return new x3dom.fields.SFColor(that.r,that.g,that.b);};x3dom.fields.SFColor.prototype.copy=function(){return x3dom.fields.SFColor.copy(this);};x3dom.fields.SFColor.prototype.setHSV=function(h,s,v){x3dom.debug.logWarning("SFColor.setHSV() NYI");};x3dom.fields.SFColor.prototype.getHSV=function(){var h=0,s=0,v=0;x3dom.debug.logWarning("SFColor.getHSV() NYI");return[h,s,v];};x3dom.fields.SFColor.prototype.setValues=function(color){this.r=color.r;this.g=color.g;this.b=color.b;};x3dom.fields.SFColor.prototype.equals=function(that,eps){return Math.abs(this.r-that.r)<eps&&Math.abs(this.g-that.g)<eps&&Math.abs(this.b-that.b)<eps;};x3dom.fields.SFColor.prototype.add=function(that){return new x3dom.fields.SFColor(this.r+that.r,this.g+that.g,this.b+that.b);};x3dom.fields.SFColor.prototype.subtract=function(that){return new x3dom.fields.SFColor(this.r-that.r,this.g-that.g,this.b-that.b);};x3dom.fields.SFColor.prototype.multiply=function(n){return new x3dom.fields.SFColor(this.r*n,this.g*n,this.b*n);};x3dom.fields.SFColor.prototype.toGL=function(){return[this.r,this.g,this.b];};x3dom.fields.SFColor.prototype.toString=function(){return this.r+" "+this.g+" "+this.b;};x3dom.fields.SFColor.prototype.setValueByStr=function(str){try{var m=/^\s*([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)\s*,?\s*([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)\s*,?\s*([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)\s*$/.exec(str);this.r=+m[1];this.g=+m[2];this.b=+m[3];}
catch(e){var c=x3dom.fields.SFColor.colorParse(str);this.r=c.r;this.g=c.g;this.b=c.b;}
return this;};x3dom.fields.SFColor.colorParse=function(color){var red=0,green=0,blue=0;var color_names={aliceblue:'f0f8ff',antiquewhite:'faebd7',aqua:'00ffff',aquamarine:'7fffd4',azure:'f0ffff',beige:'f5f5dc',bisque:'ffe4c4',black:'000000',blanchedalmond:'ffebcd',blue:'0000ff',blueviolet:'8a2be2',brown:'a52a2a',burlywood:'deb887',cadetblue:'5f9ea0',chartreuse:'7fff00',chocolate:'d2691e',coral:'ff7f50',cornflowerblue:'6495ed',cornsilk:'fff8dc',crimson:'dc143c',cyan:'00ffff',darkblue:'00008b',darkcyan:'008b8b',darkgoldenrod:'b8860b',darkgray:'a9a9a9',darkgreen:'006400',darkkhaki:'bdb76b',darkmagenta:'8b008b',darkolivegreen:'556b2f',darkorange:'ff8c00',darkorchid:'9932cc',darkred:'8b0000',darksalmon:'e9967a',darkseagreen:'8fbc8f',darkslateblue:'483d8b',darkslategray:'2f4f4f',darkturquoise:'00ced1',darkviolet:'9400d3',deeppink:'ff1493',deepskyblue:'00bfff',dimgray:'696969',dodgerblue:'1e90ff',feldspar:'d19275',firebrick:'b22222',floralwhite:'fffaf0',forestgreen:'228b22',fuchsia:'ff00ff',gainsboro:'dcdcdc',ghostwhite:'f8f8ff',gold:'ffd700',goldenrod:'daa520',gray:'808080',green:'008000',greenyellow:'adff2f',honeydew:'f0fff0',hotpink:'ff69b4',indianred:'cd5c5c',indigo:'4b0082',ivory:'fffff0',khaki:'f0e68c',lavender:'e6e6fa',lavenderblush:'fff0f5',lawngreen:'7cfc00',lemonchiffon:'fffacd',lightblue:'add8e6',lightcoral:'f08080',lightcyan:'e0ffff',lightgoldenrodyellow:'fafad2',lightgrey:'d3d3d3',lightgreen:'90ee90',lightpink:'ffb6c1',lightsalmon:'ffa07a',lightseagreen:'20b2aa',lightskyblue:'87cefa',lightslateblue:'8470ff',lightslategray:'778899',lightsteelblue:'b0c4de',lightyellow:'ffffe0',lime:'00ff00',limegreen:'32cd32',linen:'faf0e6',magenta:'ff00ff',maroon:'800000',mediumaquamarine:'66cdaa',mediumblue:'0000cd',mediumorchid:'ba55d3',mediumpurple:'9370d8',mediumseagreen:'3cb371',mediumslateblue:'7b68ee',mediumspringgreen:'00fa9a',mediumturquoise:'48d1cc',mediumvioletred:'c71585',midnightblue:'191970',mintcream:'f5fffa',mistyrose:'ffe4e1',moccasin:'ffe4b5',navajowhite:'ffdead',navy:'000080',oldlace:'fdf5e6',olive:'808000',olivedrab:'6b8e23',orange:'ffa500',orangered:'ff4500',orchid:'da70d6',palegoldenrod:'eee8aa',palegreen:'98fb98',paleturquoise:'afeeee',palevioletred:'d87093',papayawhip:'ffefd5',peachpuff:'ffdab9',peru:'cd853f',pink:'ffc0cb',plum:'dda0dd',powderblue:'b0e0e6',purple:'800080',red:'ff0000',rosybrown:'bc8f8f',royalblue:'4169e1',saddlebrown:'8b4513',salmon:'fa8072',sandybrown:'f4a460',seagreen:'2e8b57',seashell:'fff5ee',sienna:'a0522d',silver:'c0c0c0',skyblue:'87ceeb',slateblue:'6a5acd',slategray:'708090',snow:'fffafa',springgreen:'00ff7f',steelblue:'4682b4',tan:'d2b48c',teal:'008080',thistle:'d8bfd8',tomato:'ff6347',turquoise:'40e0d0',violet:'ee82ee',violetred:'d02090',wheat:'f5deb3',white:'ffffff',whitesmoke:'f5f5f5',yellow:'ffff00',yellowgreen:'9acd32'};if(color_names[color]){color="#"+color_names[color];}
if(color.substr&&color.substr(0,1)==="#"){color=color.substr(1);var len=color.length;if(len===6){red=parseInt("0x"+color.substr(0,2),16)/255.0;green=parseInt("0x"+color.substr(2,2),16)/255.0;blue=parseInt("0x"+color.substr(4,2),16)/255.0;}
else if(len===3){red=parseInt("0x"+color.substr(0,1),16)/15.0;green=parseInt("0x"+color.substr(1,1),16)/15.0;blue=parseInt("0x"+color.substr(2,1),16)/15.0;}}
return new x3dom.fields.SFColor(red,green,blue);};x3dom.fields.SFColorRGBA=function(r,g,b,a){if(arguments.length===0){this.r=0;this.g=0;this.b=0;this.a=1;}
else{this.r=r;this.g=g;this.b=b;this.a=a;}};x3dom.fields.SFColorRGBA.parse=function(str){try{var m=/^([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)\s*,?\s*([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)\s*,?\s*([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)\s*,?\s*([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)$/.exec(str);return new x3dom.fields.SFColorRGBA(+m[1],+m[2],+m[3],+m[4]);}
catch(e){return x3dom.fields.SFColorRGBA.colorParse(str);}};x3dom.fields.SFColorRGBA.copy=function(that){return new x3dom.fields.SFColorRGBA(that.r,that.g,that.b,that.a);};x3dom.fields.SFColorRGBA.prototype.copy=function(){return x3dom.fields.SFColorRGBA.copy(this);};x3dom.fields.SFColorRGBA.prototype.setValues=function(color){this.r=color.r;this.g=color.g;this.b=color.b;this.a=color.a;};x3dom.fields.SFColorRGBA.prototype.equals=function(that,eps){return Math.abs(this.r-that.r)<eps&&Math.abs(this.g-that.g)<eps&&Math.abs(this.b-that.b)<eps&&Math.abs(this.a-that.a)<eps;};x3dom.fields.SFColorRGBA.prototype.toGL=function(){return[this.r,this.g,this.b,this.a];};x3dom.fields.SFColorRGBA.prototype.toString=function(){return this.r+" "+this.g+" "+this.b+" "+this.a;};x3dom.fields.SFColorRGBA.prototype.setValueByStr=function(str){try{var m=/^([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)\s*,?\s*([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)\s*,?\s*([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)\s*,?\s*([+\-]?\d*\.*\d*[eE]?[+\-]?\d*?)$/.exec(str);this.r=+m[1];this.g=+m[2];this.b=+m[3];this.a=+m[4];}
catch(e){var c=x3dom.fields.SFColorRGBA.colorParse(str);this.r=c.r;this.g=c.g;this.b=c.b;this.a=c.a;}
return this;};x3dom.fields.SFColorRGBA.prototype.toUint=function(){return((Math.round(this.r*255)<<24)|(Math.round(this.g*255)<<16)|(Math.round(this.b*255)<<8)|Math.round(this.a*255))>>>0;};x3dom.fields.SFColorRGBA.colorParse=function(color){var red=0,green=0,blue=0,alpha=0;var color_names={aliceblue:'f0f8ff',antiquewhite:'faebd7',aqua:'00ffff',aquamarine:'7fffd4',azure:'f0ffff',beige:'f5f5dc',bisque:'ffe4c4',black:'000000',blanchedalmond:'ffebcd',blue:'0000ff',blueviolet:'8a2be2',brown:'a52a2a',burlywood:'deb887',cadetblue:'5f9ea0',chartreuse:'7fff00',chocolate:'d2691e',coral:'ff7f50',cornflowerblue:'6495ed',cornsilk:'fff8dc',crimson:'dc143c',cyan:'00ffff',darkblue:'00008b',darkcyan:'008b8b',darkgoldenrod:'b8860b',darkgray:'a9a9a9',darkgreen:'006400',darkkhaki:'bdb76b',darkmagenta:'8b008b',darkolivegreen:'556b2f',darkorange:'ff8c00',darkorchid:'9932cc',darkred:'8b0000',darksalmon:'e9967a',darkseagreen:'8fbc8f',darkslateblue:'483d8b',darkslategray:'2f4f4f',darkturquoise:'00ced1',darkviolet:'9400d3',deeppink:'ff1493',deepskyblue:'00bfff',dimgray:'696969',dodgerblue:'1e90ff',feldspar:'d19275',firebrick:'b22222',floralwhite:'fffaf0',forestgreen:'228b22',fuchsia:'ff00ff',gainsboro:'dcdcdc',ghostwhite:'f8f8ff',gold:'ffd700',goldenrod:'daa520',gray:'808080',green:'008000',greenyellow:'adff2f',honeydew:'f0fff0',hotpink:'ff69b4',indianred:'cd5c5c',indigo:'4b0082',ivory:'fffff0',khaki:'f0e68c',lavender:'e6e6fa',lavenderblush:'fff0f5',lawngreen:'7cfc00',lemonchiffon:'fffacd',lightblue:'add8e6',lightcoral:'f08080',lightcyan:'e0ffff',lightgoldenrodyellow:'fafad2',lightgrey:'d3d3d3',lightgreen:'90ee90',lightpink:'ffb6c1',lightsalmon:'ffa07a',lightseagreen:'20b2aa',lightskyblue:'87cefa',lightslateblue:'8470ff',lightslategray:'778899',lightsteelblue:'b0c4de',lightyellow:'ffffe0',lime:'00ff00',limegreen:'32cd32',linen:'faf0e6',magenta:'ff00ff',maroon:'800000',mediumaquamarine:'66cdaa',mediumblue:'0000cd',mediumorchid:'ba55d3',mediumpurple:'9370d8',mediumseagreen:'3cb371',mediumslateblue:'7b68ee',mediumspringgreen:'00fa9a',mediumturquoise:'48d1cc',mediumvioletred:'c71585',midnightblue:'191970',mintcream:'f5fffa',mistyrose:'ffe4e1',moccasin:'ffe4b5',navajowhite:'ffdead',navy:'000080',oldlace:'fdf5e6',olive:'808000',olivedrab:'6b8e23',orange:'ffa500',orangered:'ff4500',orchid:'da70d6',palegoldenrod:'eee8aa',palegreen:'98fb98',paleturquoise:'afeeee',palevioletred:'d87093',papayawhip:'ffefd5',peachpuff:'ffdab9',peru:'cd853f',pink:'ffc0cb',plum:'dda0dd',powderblue:'b0e0e6',purple:'800080',red:'ff0000',rosybrown:'bc8f8f',royalblue:'4169e1',saddlebrown:'8b4513',salmon:'fa8072',sandybrown:'f4a460',seagreen:'2e8b57',seashell:'fff5ee',sienna:'a0522d',silver:'c0c0c0',skyblue:'87ceeb',slateblue:'6a5acd',slategray:'708090',snow:'fffafa',springgreen:'00ff7f',steelblue:'4682b4',tan:'d2b48c',teal:'008080',thistle:'d8bfd8',tomato:'ff6347',turquoise:'40e0d0',violet:'ee82ee',violetred:'d02090',wheat:'f5deb3',white:'ffffff',whitesmoke:'f5f5f5',yellow:'ffff00',yellowgreen:'9acd32'};if(color_names[color]){color="#"+color_names[color]+"ff";}
if(color.substr&&color.substr(0,1)==="#"){color=color.substr(1);var len=color.length;if(len===8){red=parseInt("0x"+color.substr(0,2),16)/255.0;green=parseInt("0x"+color.substr(2,2),16)/255.0;blue=parseInt("0x"+color.substr(4,2),16)/255.0;alpha=parseInt("0x"+color.substr(6,2),16)/255.0;}
else if(len===6){red=parseInt("0x"+color.substr(0,2),16)/255.0;green=parseInt("0x"+color.substr(2,2),16)/255.0;blue=parseInt("0x"+color.substr(4,2),16)/255.0;alpha=1.0;}
else if(len===4){red=parseInt("0x"+color.substr(0,1),16)/15.0;green=parseInt("0x"+color.substr(1,1),16)/15.0;blue=parseInt("0x"+color.substr(2,1),16)/15.0;alpha=parseInt("0x"+color.substr(3,1),16)/15.0;}
else if(len===3){red=parseInt("0x"+color.substr(0,1),16)/15.0;green=parseInt("0x"+color.substr(1,1),16)/15.0;blue=parseInt("0x"+color.substr(2,1),16)/15.0;alpha=1.0;}}
return new x3dom.fields.SFColorRGBA(red,green,blue,alpha);};x3dom.fields.SFImage=function(w,h,c,arr){if(arguments.length===0||!(arr&&arr.map)){this.width=0;this.height=0;this.comp=0;this.array=[];}
else{this.width=w;this.height=h;this.comp=c;var that=this.array;arr.map(function(v){that.push(v);},this.array);}};x3dom.fields.SFImage.parse=function(str){var img=new x3dom.fields.SFImage();img.setValueByStr(str);return img;};x3dom.fields.SFImage.copy=function(that){var destination=new x3dom.fields.SFImage();destination.width=that.width;destination.height=that.height;destination.comp=that.comp;destination.setPixels(that.getPixels());return destination;};x3dom.fields.SFImage.prototype.copy=function(){return x3dom.fields.SFImage.copy(this);};x3dom.fields.SFImage.prototype.setValueByStr=function(str){var mc=str.match(/(\w+)/g);var n=mc.length;var c2=0;var hex="0123456789ABCDEF";this.array=[];if(n>2){this.width=+mc[0];this.height=+mc[1];this.comp=+mc[2];c2=2*this.comp;}else{this.width=0;this.height=0;this.comp=0;return;}
var len,i;for(i=3;i<n;i++){var r,g,b,a;if(!mc[i].substr){continue;}
if(mc[i].substr(1,1).toLowerCase()!=="x"){var inp=parseInt(mc[i],10);if(this.comp===1){r=inp;this.array.push(r);}
else if(this.comp===2){r=inp>>8&255;g=inp&255;this.array.push(r,g);}
else if(this.comp===3){r=inp>>16&255;g=inp>>8&255;b=inp&255;this.array.push(r,g,b);}
else if(this.comp===4){r=inp>>24&255;g=inp>>16&255;b=inp>>8&255;a=inp&255;this.array.push(r,g,b,a);}}
else if(mc[i].substr(1,1).toLowerCase()==="x"){mc[i]=mc[i].substr(2);len=mc[i].length;if(len===c2){if(this.comp===1){r=parseInt("0x"+mc[i].substr(0,2),16);this.array.push(r);}
else if(this.comp===2){r=parseInt("0x"+mc[i].substr(0,2),16);g=parseInt("0x"+mc[i].substr(2,2),16);this.array.push(r,g);}
else if(this.comp===3){r=parseInt("0x"+mc[i].substr(0,2),16);g=parseInt("0x"+mc[i].substr(2,2),16);b=parseInt("0x"+mc[i].substr(4,2),16);this.array.push(r,g,b);}
else if(this.comp===4){r=parseInt("0x"+mc[i].substr(0,2),16);g=parseInt("0x"+mc[i].substr(2,2),16);b=parseInt("0x"+mc[i].substr(4,2),16);a=parseInt("0x"+mc[i].substr(6,2),16);this.array.push(r,g,b,a);}}}}};x3dom.fields.SFImage.prototype.setPixel=function(x,y,color){var startIdx=(y*this.width+x)*this.comp;if(this.comp===1&&startIdx<this.array.length){this.array[startIdx]=color.r*255;}
else if(this.comp===2&&(startIdx+1)<this.array.length){this.array[startIdx]=color.r*255;this.array[startIdx+1]=color.g*255;}
else if(this.comp===3&&(startIdx+2)<this.array.length){this.array[startIdx]=color.r*255;this.array[startIdx+1]=color.g*255;this.array[startIdx+2]=color.b*255;}
else if(this.comp===4&&(startIdx+3)<this.array.length){this.array[startIdx]=color.r*255;this.array[startIdx+1]=color.g*255;this.array[startIdx+2]=color.b*255;this.array[startIdx+3]=color.a*255;}};x3dom.fields.SFImage.prototype.getPixel=function(x,y){var startIdx=(y*this.width+x)*this.comp;if(this.comp===1&&startIdx<this.array.length){return new x3dom.fields.SFColorRGBA(this.array[startIdx]/255,0,0,1);}
else if(this.comp===2&&(startIdx+1)<this.array.length){return new x3dom.fields.SFColorRGBA(this.array[startIdx]/255,this.array[startIdx+1]/255,0,1);}
else if(this.comp===3&&(startIdx+2)<this.array.length){return new x3dom.fields.SFColorRGBA(this.array[startIdx]/255,this.array[startIdx+1]/255,this.array[startIdx+2]/255,1);}
else if(this.comp===4&&(startIdx+3)<this.array.length){return new x3dom.fields.SFColorRGBA(this.array[startIdx]/255,this.array[startIdx+1]/255,this.array[startIdx+2]/255,this.array[startIdx+3]/255);}};x3dom.fields.SFImage.prototype.setPixels=function(pixels){var i,idx=0;if(this.comp===1){for(i=0;i<pixels.length;i++){this.array[idx++]=pixels[i].r*255;}}
else if(this.comp===2){for(i=0;i<pixels.length;i++){this.array[idx++]=pixels[i].r*255;this.array[idx++]=pixels[i].g*255;}}
else if(this.comp===3){for(i=0;i<pixels.length;i++){this.array[idx++]=pixels[i].r*255;this.array[idx++]=pixels[i].g*255;this.array[idx++]=pixels[i].b*255;}}
else if(this.comp===4){for(i=0;i<pixels.length;i++){this.array[idx++]=pixels[i].r*255;this.array[idx++]=pixels[i].g*255;this.array[idx++]=pixels[i].b*255;this.array[idx++]=pixels[i].a*255;}}};x3dom.fields.SFImage.prototype.getPixels=function(){var i;var pixels=[];if(this.comp===1){for(i=0;i<this.array.length;i+=this.comp){pixels.push(new x3dom.fields.SFColorRGBA(this.array[i]/255,0,0,1));}}
else if(this.comp===2){for(i=0;i<this.array.length;i+=this.comp){pixels.push(new x3dom.fields.SFColorRGBA(this.array[i]/255,this.array[i+1]/255,0,1));}}
else if(this.comp===3){for(i=0;i<this.array.length;i+=this.comp){pixels.push(new x3dom.fields.SFColorRGBA(this.array[i]/255,this.array[i+1]/255,this.array[i+2]/255,1));}}
else if(this.comp===4){for(i=0;i<this.array.length;i+=this.comp){pixels.push(new x3dom.fields.SFColorRGBA(this.array[i]/255,this.array[i+1]/255,this.array[i+2]/255,this.array[i+3]/255));}}
return pixels;};x3dom.fields.SFImage.prototype.toGL=function(){var a=[];Array.map(this.array,function(c){a.push(c);});return a;};x3dom.fields.MFColor=function(colorArray){if(colorArray){var that=this;colorArray.map(function(c){that.push(c);},this);}};x3dom.fields.MFColor.copy=function(colorArray){var destination=new x3dom.fields.MFColor();colorArray.map(function(v){destination.push(v.copy());},this);return destination;};x3dom.fields.MFColor.prototype=x3dom.extend([]);x3dom.fields.MFColor.parse=function(str){var mc=str.match(/([+\-0-9eE\.]+)/g);var colors=[];for(var i=0,n=mc?mc.length:0;i<n;i+=3){colors.push(new x3dom.fields.SFColor(+mc[i+0],+mc[i+1],+mc[i+2]));}
return new x3dom.fields.MFColor(colors);};x3dom.fields.MFColor.prototype.copy=function(){return x3dom.fields.MFColor.copy(this);};x3dom.fields.MFColor.prototype.setValueByStr=function(str){this.length=0;var mc=str.match(/([+\-0-9eE\.]+)/g);for(var i=0,n=mc?mc.length:0;i<n;i+=3){this.push(new x3dom.fields.SFColor(+mc[i+0],+mc[i+1],+mc[i+2]));}};x3dom.fields.MFColor.prototype.toGL=function(){var a=[];Array.map(this,function(c){a.push(c.r);a.push(c.g);a.push(c.b);});return a;};x3dom.fields.MFColorRGBA=function(colorArray){if(colorArray){var that=this;colorArray.map(function(c){that.push(c);},this);}};x3dom.fields.MFColorRGBA.copy=function(colorArray){var destination=new x3dom.fields.MFColorRGBA();colorArray.map(function(v){destination.push(v.copy());},this);return destination;};x3dom.fields.MFColorRGBA.prototype=x3dom.extend([]);x3dom.fields.MFColorRGBA.parse=function(str){var mc=str.match(/([+\-0-9eE\.]+)/g);var colors=[];for(var i=0,n=mc?mc.length:0;i<n;i+=4){colors.push(new x3dom.fields.SFColorRGBA(+mc[i+0],+mc[i+1],+mc[i+2],+mc[i+3]));}
return new x3dom.fields.MFColorRGBA(colors);};x3dom.fields.MFColorRGBA.prototype.copy=function(){return x3dom.fields.MFColorRGBA.copy(this);};x3dom.fields.MFColorRGBA.prototype.setValueByStr=function(str){this.length=0;var mc=str.match(/([+\-0-9eE\.]+)/g);for(var i=0,n=mc?mc.length:0;i<n;i+=4){this.push(new x3dom.fields.SFColorRGBA(+mc[i+0],+mc[i+1],+mc[i+2],+mc[i+3]));}};x3dom.fields.MFColorRGBA.prototype.toGL=function(){var a=[];Array.map(this,function(c){a.push(c.r);a.push(c.g);a.push(c.b);a.push(c.a);});return a;};x3dom.fields.MFRotation=function(rotArray){if(rotArray){var that=this;rotArray.map(function(v){that.push(v);},this);}};x3dom.fields.MFRotation.prototype=x3dom.extend([]);x3dom.fields.MFRotation.copy=function(rotationArray){var destination=new x3dom.fields.MFRotation();rotationArray.map(function(v){destination.push(v.copy());},this);return destination;};x3dom.fields.MFRotation.prototype.copy=function(){return x3dom.fields.MFRotation.copy(this);};x3dom.fields.MFRotation.parse=function(str){var mc=str.match(/([+\-0-9eE\.]+)/g);var vecs=[];for(var i=0,n=mc?mc.length:0;i<n;i+=4){vecs.push(x3dom.fields.Quaternion.axisAngle(new x3dom.fields.SFVec3f(+mc[i+0],+mc[i+1],+mc[i+2]),+mc[i+3]));}
return new x3dom.fields.MFRotation(vecs);};x3dom.fields.MFRotation.prototype.setValueByStr=function(str){this.length=0;var mc=str.match(/([+\-0-9eE\.]+)/g);for(var i=0,n=mc?mc.length:0;i<n;i+=4){this.push(x3dom.fields.Quaternion.axisAngle(new x3dom.fields.SFVec3f(+mc[i+0],+mc[i+1],+mc[i+2]),+mc[i+3]));}};x3dom.fields.MFRotation.prototype.toGL=function(){var a=[];Array.map(this,function(c){var val=c.toAxisAngle();a.push(val[0].x);a.push(val[0].y);a.push(val[0].z);a.push(val[1]);});return a;};x3dom.fields.MFVec3f=function(vec3Array){if(vec3Array){var that=this;vec3Array.map(function(v){that.push(v);},this);}};x3dom.fields.MFVec3f.prototype=x3dom.extend(Array);x3dom.fields.MFVec3f.copy=function(vec3Array){var destination=new x3dom.fields.MFVec3f();vec3Array.map(function(v){destination.push(v.copy());},this);return destination;};x3dom.fields.MFVec3f.parse=function(str){var mc=str.match(/([+\-0-9eE\.]+)/g);var vecs=[];for(var i=0,n=mc?mc.length:0;i<n;i+=3){vecs.push(new x3dom.fields.SFVec3f(+mc[i+0],+mc[i+1],+mc[i+2]));}
return new x3dom.fields.MFVec3f(vecs);};x3dom.fields.MFVec3f.prototype.copy=function()
{x3dom.fields.MFVec3f.copy(this);};x3dom.fields.MFVec3f.prototype.setValueByStr=function(str){this.length=0;var mc=str.match(/([+\-0-9eE\.]+)/g);for(var i=0,n=mc?mc.length:0;i<n;i+=3){this.push(new x3dom.fields.SFVec3f(+mc[i+0],+mc[i+1],+mc[i+2]));}};x3dom.fields.MFVec3f.prototype.toGL=function(){var a=[];Array.map(this,function(c){a.push(c.x);a.push(c.y);a.push(c.z);});return a;};x3dom.fields.MFVec2f=function(vec2Array){if(vec2Array){var that=this;vec2Array.map(function(v){that.push(v);},this);}};x3dom.fields.MFVec2f.prototype=x3dom.extend([]);x3dom.fields.MFVec2f.copy=function(vec2Array){var destination=new x3dom.fields.MFVec2f();vec2Array.map(function(v){destination.push(v.copy());},this);return destination;};x3dom.fields.MFVec2f.parse=function(str){var mc=str.match(/([+\-0-9eE\.]+)/g);var vecs=[];for(var i=0,n=mc?mc.length:0;i<n;i+=2){vecs.push(new x3dom.fields.SFVec2f(+mc[i+0],+mc[i+1]));}
return new x3dom.fields.MFVec2f(vecs);};x3dom.fields.MFVec2f.prototype.copy=function(){return x3dom.fields.MFVec2f.copy(this);};x3dom.fields.MFVec2f.prototype.setValueByStr=function(str){this.length=0;var mc=str.match(/([+\-0-9eE\.]+)/g);for(var i=0,n=mc?mc.length:0;i<n;i+=2){this.push(new x3dom.fields.SFVec2f(+mc[i+0],+mc[i+1]));}};x3dom.fields.MFVec2f.prototype.toGL=function(){var a=[];Array.map(this,function(v){a.push(v.x);a.push(v.y);});return a;};x3dom.fields.MFInt32=function(array){if(array){var that=this;array.map(function(v){that.push(v);},this);}};x3dom.fields.MFInt32.prototype=x3dom.extend([]);x3dom.fields.MFInt32.copy=function(intArray){var destination=new x3dom.fields.MFInt32();intArray.map(function(v){destination.push(v);},this);return destination;};x3dom.fields.MFInt32.parse=function(str){var mc=str.match(/([+\-]?\d+\s*){1},?\s*/g);var vals=[];for(var i=0,n=mc?mc.length:0;i<n;++i){vals.push(parseInt(mc[i],10));}
return new x3dom.fields.MFInt32(vals);};x3dom.fields.MFInt32.prototype.copy=function(){return x3dom.fields.MFInt32.copy(this);};x3dom.fields.MFInt32.prototype.setValueByStr=function(str){this.length=0;var mc=str.match(/([+\-]?\d+\s*){1},?\s*/g);for(var i=0,n=mc?mc.length:0;i<n;++i){this.push(parseInt(mc[i],10));}};x3dom.fields.MFInt32.prototype.toGL=function(){var a=[];Array.map(this,function(v){a.push(v);});return a;};x3dom.fields.MFFloat=function(array){if(array){var that=this;array.map(function(v){that.push(v);},this);}};x3dom.fields.MFFloat.prototype=x3dom.extend([]);x3dom.fields.MFFloat.copy=function(floatArray){var destination=new x3dom.fields.MFFloat();floatArray.map(function(v){destination.push(v);},this);return destination;};x3dom.fields.MFFloat.parse=function(str){var mc=str.match(/([+\-0-9eE\.]+)/g);var vals=[];for(var i=0,n=mc?mc.length:0;i<n;i++){vals.push(+mc[i]);}
return new x3dom.fields.MFFloat(vals);};x3dom.fields.MFFloat.prototype.copy=function(){return x3dom.fields.MFFloat.copy(this);};x3dom.fields.MFFloat.prototype.setValueByStr=function(str){this.length=0;var mc=str.match(/([+\-0-9eE\.]+)/g);for(var i=0,n=mc?mc.length:0;i<n;i++){this.push(+mc[i]);}};x3dom.fields.MFFloat.prototype.toGL=function(){var a=[];Array.map(this,function(v){a.push(v);});return a;};x3dom.fields.MFBoolean=function(array){if(array){var that=this;array.map(function(v){that.push(v);},this);}};x3dom.fields.MFBoolean.prototype=x3dom.extend([]);x3dom.fields.MFBoolean.copy=function(boolArray){var destination=new x3dom.fields.MFBoolean();boolArray.map(function(v){destination.push(v);},this);return destination;};x3dom.fields.MFBoolean.parse=function(str){var mc=str.match(/(true|false|1|0)/ig);var vals=[];for(var i=0,n=mc?mc.length:0;i<n;i++){vals.push((mc[i]=='1'||mc[i].toLowerCase()=='true'));}
return new x3dom.fields.MFBoolean(vals);};x3dom.fields.MFBoolean.prototype.copy=function(){return x3dom.fields.MFBoolean.copy(this);};x3dom.fields.MFBoolean.prototype.setValueByStr=function(str){this.length=0;var mc=str.match(/(true|false|1|0)/ig);for(var i=0,n=mc?mc.length:0;i<n;i++){this.push((mc[i]=='1'||mc[i].toLowerCase()=='true'));}};x3dom.fields.MFBoolean.prototype.toGL=function(){var a=[];Array.map(this,function(v){a.push(v?1:0);});return a;};x3dom.fields.MFString=function(strArray){if(strArray&&strArray.map){var that=this;strArray.map(function(v){that.push(v);},this);}};x3dom.fields.MFString.prototype=x3dom.extend([]);x3dom.fields.MFString.copy=function(stringArray){var destination=new x3dom.fields.MFString();stringArray.map(function(v){destination.push(v);},this);return destination;};x3dom.fields.MFString.parse=function(str){var arr=[];if(str.length&&str[0]=='"'){var m,re=/"((?:[^\\"]|\\\\|\\")*)"/g;while((m=re.exec(str))){var s=m[1].replace(/\\([\\"])/g,"$1");if(s!==undefined){arr.push(s);}}}
else{arr.push(str);}
return new x3dom.fields.MFString(arr);};x3dom.fields.MFString.prototype.copy=function(){return x3dom.fields.MFString.copy(this);};x3dom.fields.MFString.prototype.setValueByStr=function(str){this.length=0;if(str.length&&str[0]=='"'){var m,re=/"((?:[^\\"]|\\\\|\\")*)"/g;while((m=re.exec(str))){var s=m[1].replace(/\\([\\"])/,"$1");if(s!==undefined){this.push(s);}}}
else{this.push(str);}
return this;};x3dom.fields.MFString.prototype.toString=function(){var str="";for(var i=0,n=this.length;i<n;i++){str=str+this[i]+" ";}
return str;};x3dom.fields.SFNode=function(type){this.type=type;this.node=null;};x3dom.fields.SFNode.prototype.hasLink=function(node){return(node?(this.node===node):this.node);};x3dom.fields.SFNode.prototype.addLink=function(node){this.node=node;return true;};x3dom.fields.SFNode.prototype.rmLink=function(node){if(this.node===node){this.node=null;return true;}
else{return false;}};x3dom.fields.MFNode=function(type){this.type=type;this.nodes=[];};x3dom.fields.MFNode.prototype.hasLink=function(node){if(node){for(var i=0,n=this.nodes.length;i<n;i++){if(this.nodes[i]===node){return true;}}}
else{return(this.length>0);}
return false;};x3dom.fields.MFNode.prototype.addLink=function(node){this.nodes.push(node);return true;};x3dom.fields.MFNode.prototype.rmLink=function(node){for(var i=0,n=this.nodes.length;i<n;i++){if(this.nodes[i]===node){this.nodes.splice(i,1);return true;}}
return false;};x3dom.fields.MFNode.prototype.length=function(){return this.nodes.length;};x3dom.fields.Line=function(pos,dir)
{if(arguments.length===0)
{this.pos=new x3dom.fields.SFVec3f(0,0,0);this.dir=new x3dom.fields.SFVec3f(0,0,1);}
this.pos=x3dom.fields.SFVec3f.copy(pos);this.dir=x3dom.fields.SFVec3f.copy(dir);};x3dom.fields.Line.prototype.closestPoint=function(p)
{var distVec=p.subtract(this.pos);var projDist=distVec.dot(this.dir);return this.pos.add(this.dir.multiply(projDist));};x3dom.fields.Line.prototype.shortestDistance=function(p)
{var distVec=p.subtract(this.pos);var projDist=distVec.dot(this.dir);return distVec.subtract(this.dir.multiply(projDist)).length();};x3dom.fields.Ray=function(pos,dir)
{if(arguments.length===0)
{this.pos=new x3dom.fields.SFVec3f(0,0,0);this.dir=new x3dom.fields.SFVec3f(0,0,1);}
else
{this.pos=new x3dom.fields.SFVec3f(pos.x,pos.y,pos.z);var n=dir.length();if(n){n=1.0/n;}
this.dir=new x3dom.fields.SFVec3f(dir.x*n,dir.y*n,dir.z*n);}
this.enter=0;this.exit=0;this.hitObject=null;this.hitPoint={};this.dist=Number.MAX_VALUE;};x3dom.fields.Ray.prototype.toString=function(){return'Ray: ['+this.pos.toString()+'; '+this.dir.toString()+']';};x3dom.fields.Ray.prototype.intersectPlane=function(p,n)
{var result=null;var alpha;var nDotDir=n.dot(this.dir);if(nDotDir<0.0)
{alpha=(p.dot(n)-this.pos.dot(n))/nDotDir;result=this.pos.addScaled(this.dir,alpha);}
return result;};x3dom.fields.Ray.prototype.intersect=function(low,high)
{var isect=0.0;var out=Number.MAX_VALUE;var r,te,tl;if(this.dir.x>x3dom.fields.Eps)
{r=1.0/this.dir.x;te=(low.x-this.pos.x)*r;tl=(high.x-this.pos.x)*r;if(tl<out){out=tl;}
if(te>isect){isect=te;}}
else if(this.dir.x<-x3dom.fields.Eps)
{r=1.0/this.dir.x;te=(high.x-this.pos.x)*r;tl=(low.x-this.pos.x)*r;if(tl<out){out=tl;}
if(te>isect){isect=te;}}
else if(this.pos.x<low.x||this.pos.x>high.x)
{return false;}
if(this.dir.y>x3dom.fields.Eps)
{r=1.0/this.dir.y;te=(low.y-this.pos.y)*r;tl=(high.y-this.pos.y)*r;if(tl<out){out=tl;}
if(te>isect){isect=te;}
if(isect-out>=x3dom.fields.Eps){return false;}}
else if(this.dir.y<-x3dom.fields.Eps)
{r=1.0/this.dir.y;te=(high.y-this.pos.y)*r;tl=(low.y-this.pos.y)*r;if(tl<out){out=tl;}
if(te>isect){isect=te;}
if(isect-out>=x3dom.fields.Eps){return false;}}
else if(this.pos.y<low.y||this.pos.y>high.y)
{return false;}
if(this.dir.z>x3dom.fields.Eps)
{r=1.0/this.dir.z;te=(low.z-this.pos.z)*r;tl=(high.z-this.pos.z)*r;if(tl<out){out=tl;}
if(te>isect){isect=te;}}
else if(this.dir.z<-x3dom.fields.Eps)
{r=1.0/this.dir.z;te=(high.z-this.pos.z)*r;tl=(low.z-this.pos.z)*r;if(tl<out){out=tl;}
if(te>isect){isect=te;}}
else if(this.pos.z<low.z||this.pos.z>high.z)
{return false;}
this.enter=isect;this.exit=out;return(isect-out<x3dom.fields.Eps);};x3dom.fields.BoxVolume=function(min,max)
{if(arguments.length<2){this.min=new x3dom.fields.SFVec3f(0,0,0);this.max=new x3dom.fields.SFVec3f(0,0,0);this.valid=false;}
else{this.min=x3dom.fields.SFVec3f.copy(min);this.max=x3dom.fields.SFVec3f.copy(max);this.valid=true;}
this.updateInternals();};x3dom.fields.BoxVolume.prototype.getScalarValue=function()
{var extent=this.max.subtract(this.min);return(extent.x*extent.y*extent.z);};x3dom.fields.BoxVolume.copy=function(other)
{var volume=new x3dom.fields.BoxVolume(other.min,other.max);volume.valid=other.valid;return volume;};x3dom.fields.BoxVolume.prototype.equals=function(other)
{return(this.min.equals(other.min,0.000000000001)&&this.max.equals(other.max,0.000000000001));};x3dom.fields.BoxVolume.prototype.updateInternals=function()
{this.radialVec=this.max.subtract(this.min).multiply(0.5);this.center=this.min.add(this.radialVec);this.diameter=2*this.radialVec.length();};x3dom.fields.BoxVolume.prototype.setBounds=function(min,max)
{this.min.setValues(min);this.max.setValues(max);this.updateInternals();this.valid=true;};x3dom.fields.BoxVolume.prototype.setBoundsByCenterSize=function(center,size)
{var halfSize=size.multiply(0.5);this.min=center.subtract(halfSize);this.max=center.add(halfSize);this.updateInternals();this.valid=true;};x3dom.fields.BoxVolume.prototype.extendBounds=function(min,max)
{if(this.valid)
{if(this.min.x>min.x){this.min.x=min.x;}
if(this.min.y>min.y){this.min.y=min.y;}
if(this.min.z>min.z){this.min.z=min.z;}
if(this.max.x<max.x){this.max.x=max.x;}
if(this.max.y<max.y){this.max.y=max.y;}
if(this.max.z<max.z){this.max.z=max.z;}
this.updateInternals();}
else
{this.setBounds(min,max);}};x3dom.fields.BoxVolume.prototype.getBounds=function(min,max)
{min.setValues(this.min);max.setValues(this.max);};x3dom.fields.BoxVolume.prototype.getRadialVec=function()
{return this.radialVec;};x3dom.fields.BoxVolume.prototype.invalidate=function()
{this.valid=false;this.min=new x3dom.fields.SFVec3f(0,0,0);this.max=new x3dom.fields.SFVec3f(0,0,0);this.updateInternals();};x3dom.fields.BoxVolume.prototype.isValid=function()
{return this.valid;};x3dom.fields.BoxVolume.prototype.getCenter=function()
{return this.center;};x3dom.fields.BoxVolume.prototype.getDiameter=function()
{return this.diameter;};x3dom.fields.BoxVolume.prototype.transform=function(m)
{var xmin,ymin,zmin;var xmax,ymax,zmax;xmin=xmax=m._03;ymin=ymax=m._13;zmin=zmax=m._23;var a=this.max.x*m._00;var b=this.min.x*m._00;if(a>=b){xmax+=a;xmin+=b;}
else{xmax+=b;xmin+=a;}
a=this.max.y*m._01;b=this.min.y*m._01;if(a>=b){xmax+=a;xmin+=b;}
else{xmax+=b;xmin+=a;}
a=this.max.z*m._02;b=this.min.z*m._02;if(a>=b){xmax+=a;xmin+=b;}
else{xmax+=b;xmin+=a;}
a=this.max.x*m._10;b=this.min.x*m._10;if(a>=b){ymax+=a;ymin+=b;}
else{ymax+=b;ymin+=a;}
a=this.max.y*m._11;b=this.min.y*m._11;if(a>=b){ymax+=a;ymin+=b;}
else{ymax+=b;ymin+=a;}
a=this.max.z*m._12;b=this.min.z*m._12;if(a>=b){ymax+=a;ymin+=b;}
else{ymax+=b;ymin+=a;}
a=this.max.x*m._20;b=this.min.x*m._20;if(a>=b){zmax+=a;zmin+=b;}
else{zmax+=b;zmin+=a;}
a=this.max.y*m._21;b=this.min.y*m._21;if(a>=b){zmax+=a;zmin+=b;}
else{zmax+=b;zmin+=a;}
a=this.max.z*m._22;b=this.min.z*m._22;if(a>=b){zmax+=a;zmin+=b;}
else{zmax+=b;zmin+=a;}
this.min.x=xmin;this.min.y=ymin;this.min.z=zmin;this.max.x=xmax;this.max.y=ymax;this.max.z=zmax;this.updateInternals();};x3dom.fields.BoxVolume.prototype.transformFrom=function(m,other)
{var xmin,ymin,zmin;var xmax,ymax,zmax;xmin=xmax=m._03;ymin=ymax=m._13;zmin=zmax=m._23;var a=other.max.x*m._00;var b=other.min.x*m._00;if(a>=b){xmax+=a;xmin+=b;}
else{xmax+=b;xmin+=a;}
a=other.max.y*m._01;b=other.min.y*m._01;if(a>=b){xmax+=a;xmin+=b;}
else{xmax+=b;xmin+=a;}
a=other.max.z*m._02;b=other.min.z*m._02;if(a>=b){xmax+=a;xmin+=b;}
else{xmax+=b;xmin+=a;}
a=other.max.x*m._10;b=other.min.x*m._10;if(a>=b){ymax+=a;ymin+=b;}
else{ymax+=b;ymin+=a;}
a=other.max.y*m._11;b=other.min.y*m._11;if(a>=b){ymax+=a;ymin+=b;}
else{ymax+=b;ymin+=a;}
a=other.max.z*m._12;b=other.min.z*m._12;if(a>=b){ymax+=a;ymin+=b;}
else{ymax+=b;ymin+=a;}
a=other.max.x*m._20;b=other.min.x*m._20;if(a>=b){zmax+=a;zmin+=b;}
else{zmax+=b;zmin+=a;}
a=other.max.y*m._21;b=other.min.y*m._21;if(a>=b){zmax+=a;zmin+=b;}
else{zmax+=b;zmin+=a;}
a=other.max.z*m._22;b=other.min.z*m._22;if(a>=b){zmax+=a;zmin+=b;}
else{zmax+=b;zmin+=a;}
this.min.x=xmin;this.min.y=ymin;this.min.z=zmin;this.max.x=xmax;this.max.y=ymax;this.max.z=zmax;this.updateInternals();this.valid=true;};x3dom.fields.FrustumVolume=function(clipMat)
{this.planeNormals=[];this.planeDistances=[];this.directionIndex=[];if(arguments.length===0){return;}
var planeEquation=[];for(var i=0;i<6;i++){this.planeNormals[i]=new x3dom.fields.SFVec3f(0,0,0);this.planeDistances[i]=0;this.directionIndex[i]=0;planeEquation[i]=new x3dom.fields.SFVec4f(0,0,0,0);}
planeEquation[0].x=clipMat._30-clipMat._00;planeEquation[0].y=clipMat._31-clipMat._01;planeEquation[0].z=clipMat._32-clipMat._02;planeEquation[0].w=clipMat._33-clipMat._03;planeEquation[1].x=clipMat._30+clipMat._00;planeEquation[1].y=clipMat._31+clipMat._01;planeEquation[1].z=clipMat._32+clipMat._02;planeEquation[1].w=clipMat._33+clipMat._03;planeEquation[2].x=clipMat._30+clipMat._10;planeEquation[2].y=clipMat._31+clipMat._11;planeEquation[2].z=clipMat._32+clipMat._12;planeEquation[2].w=clipMat._33+clipMat._13;planeEquation[3].x=clipMat._30-clipMat._10;planeEquation[3].y=clipMat._31-clipMat._11;planeEquation[3].z=clipMat._32-clipMat._12;planeEquation[3].w=clipMat._33-clipMat._13;planeEquation[4].x=clipMat._30+clipMat._20;planeEquation[4].y=clipMat._31+clipMat._21;planeEquation[4].z=clipMat._32+clipMat._22;planeEquation[4].w=clipMat._33+clipMat._23;planeEquation[5].x=clipMat._30-clipMat._20;planeEquation[5].y=clipMat._31-clipMat._21;planeEquation[5].z=clipMat._32-clipMat._22;planeEquation[5].w=clipMat._33-clipMat._23;for(i=0;i<6;i++){var vectorLength=Math.sqrt(planeEquation[i].x*planeEquation[i].x+
planeEquation[i].y*planeEquation[i].y+
planeEquation[i].z*planeEquation[i].z);planeEquation[i].x/=vectorLength;planeEquation[i].y/=vectorLength;planeEquation[i].z/=vectorLength;planeEquation[i].w/=-vectorLength;}
var updateDirectionIndex=function(normalVec){var ind=0;if(normalVec.x>0)ind|=1;if(normalVec.y>0)ind|=2;if(normalVec.z>0)ind|=4;return ind;};this.planeNormals[3].setValues(planeEquation[0]);this.planeDistances[3]=planeEquation[0].w;this.directionIndex[3]=updateDirectionIndex(this.planeNormals[3]);this.planeNormals[2].setValues(planeEquation[1]);this.planeDistances[2]=planeEquation[1].w;this.directionIndex[2]=updateDirectionIndex(this.planeNormals[2]);this.planeNormals[5].setValues(planeEquation[2]);this.planeDistances[5]=planeEquation[2].w;this.directionIndex[5]=updateDirectionIndex(this.planeNormals[5]);this.planeNormals[4].setValues(planeEquation[3]);this.planeDistances[4]=planeEquation[3].w;this.directionIndex[4]=updateDirectionIndex(this.planeNormals[4]);this.planeNormals[0].setValues(planeEquation[4]);this.planeDistances[0]=planeEquation[4].w;this.directionIndex[0]=updateDirectionIndex(this.planeNormals[0]);this.planeNormals[1].setValues(planeEquation[5]);this.planeDistances[1]=planeEquation[5].w;this.directionIndex[1]=updateDirectionIndex(this.planeNormals[1]);};x3dom.fields.FrustumVolume.prototype.intersect=function(vol,planeMask)
{if(this.planeNormals.length<6){x3dom.debug.logWarning("FrustumVolume not initialized!");return false;}
var that=this;var min=vol.min,max=vol.max;var setDirectionIndexPoint=function(index){var pnt=new x3dom.fields.SFVec3f(0,0,0);if(index&1){pnt.x=min.x;}
else{pnt.x=max.x;}
if(index&2){pnt.y=min.y;}
else{pnt.y=max.y;}
if(index&4){pnt.z=min.z;}
else{pnt.z=max.z;}
return pnt;};var pntIsInHalfSpace=function(i,pnt){var s=that.planeNormals[i].dot(pnt)-that.planeDistances[i];return(s>=0);};var isInHalfSpace=function(i){var p=setDirectionIndexPoint(that.directionIndex[i]);return pntIsInHalfSpace(i,p);};var isOutHalfSpace=function(i){var p=setDirectionIndexPoint(that.directionIndex[i]^7);return!pntIsInHalfSpace(i,p);};var mask=1;if(planeMask<0)planeMask=0;for(var i=0;i<6;i++,mask<<=1){if((planeMask&mask)!=0)
continue;if(isOutHalfSpace(i))
return-1;if(isInHalfSpace(i))
planeMask|=mask;}
return planeMask;};x3dom.docs={};x3dom.docs.specURLMap={CADGeometry:"CADGeometry.html",Core:"core.html",DIS:"dis.html",CubeMapTexturing:"env_texture.html",EnvironmentalEffects:"enveffects.html",EnvironmentalSensor:"envsensor.html",Followers:"followers.html",Geospatial:"geodata.html",Geometry2D:"geometry2D.html",Geometry3D:"geometry3D.html",Grouping:"group.html","H-Anim":"hanim.html",Interpolation:"interp.html",KeyDeviceSensor:"keyboard.html",Layering:"layering.html",Layout:"layout.html",Lighting:"lighting.html",Navigation:"navigation.html",Networking:"networking.html",NURBS:"nurbs.html",ParticleSystems:"particle_systems.html",Picking:"picking.html",PointingDeviceSensor:"pointingsensor.html",Rendering:"rendering.html",RigidBodyPhysics:"rigid_physics.html",Scripting:"scripting.html",Shaders:"shaders.html",Shape:"shape.html",Sound:"sound.html",Text:"text.html",Texturing3D:"texture3D.html",Texturing:"texturing.html",Time:"time.html",EventUtilities:"utils.html",VolumeRendering:"volume.html"};x3dom.docs.specBaseURL="http://www.web3d.org/x3d/specifications/ISO-IEC-19775-1.2-X3D-AbstractSpecification/Part01/components/";x3dom.docs.getNodeTreeInfo=function(){var tn,t;var types="";var objInArray=function(array,obj){for(var i=0;i<array.length;i++){if(array[i]===obj){return true;}}
return false;};var dump=function(t,indent){for(var i=0;i<indent;i++){types+="&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;";}
types+="<a href='"+
x3dom.docs.specBaseURL+x3dom.docs.specURLMap[x3dom.nodeTypes[t]._compName]+"#"+t+"' style='color:black; text-decoration:none; font-weight:bold;'>"+
t+"</a> &nbsp; <a href='"+
x3dom.docs.specBaseURL+x3dom.docs.specURLMap[x3dom.nodeTypes[t]._compName]+"' style='color:black; text-decoration:none; font-style:italic;'>"+
x3dom.nodeTypes[t]._compName+"</a><br/>";for(var i in x3dom.nodeTypes[t].childTypes[t]){dump(x3dom.nodeTypes[t].childTypes[t][i],indent+1);}};for(tn in x3dom.nodeTypes){var t=x3dom.nodeTypes[tn];if(t.childTypes===undefined){t.childTypes={};}
while(t.superClass){if(t.superClass.childTypes[t.superClass._typeName]===undefined){t.superClass.childTypes[t.superClass._typeName]=[];}
if(!objInArray(t.superClass.childTypes[t.superClass._typeName],t._typeName)){t.superClass.childTypes[t.superClass._typeName].push(t._typeName);}
t=t.superClass;}}
dump("X3DNode",0);return"<div class='x3dom-doc-nodes-tree'>"+types+"</div>";};x3dom.docs.getComponentInfo=function(){var components=[];var component;var result="";var c,cn;for(c in x3dom.components){components.push(c);}
components.sort();for(cn in components){c=components[cn];component=x3dom.components[c];result+="<h2><a href='"+
x3dom.docs.specBaseURL+x3dom.docs.specURLMap[c]+"' style='color:black; text-decoration:none; font-style:italic;'>"+
c+"</a></h2>";result+="<ul style='list-style-type:circle;'>";for(var t in component){result+="<li><a href='"+
x3dom.docs.specBaseURL+x3dom.docs.specURLMap[c]+"#"+t+"' style='color:black; text-decoration:none; font-weight:bold;'>"+
t+"</a></li>";}
result+="</ul>";}
return result;};x3dom.shader={};x3dom.shader.PICKING="picking";x3dom.shader.PICKING_24="picking24";x3dom.shader.PICKING_ID="pickingId";x3dom.shader.PICKING_COLOR="pickingColor";x3dom.shader.PICKING_TEXCOORD="pickingTexCoord";x3dom.shader.FRONTGROUND_TEXTURE="frontgroundTexture";x3dom.shader.BACKGROUND_TEXTURE="backgroundTexture";x3dom.shader.BACKGROUND_SKYTEXTURE="backgroundSkyTexture";x3dom.shader.BACKGROUND_CUBETEXTURE="backgroundCubeTexture";x3dom.shader.BLUR="blur";x3dom.shader.DEPTH="depth";x3dom.shader.NORMAL="normal";x3dom.shader.TEXTURE_REFINEMENT="textureRefinement";x3dom.shader.SSAO="ssao";x3dom.shader.material=function(){var shaderPart="uniform vec3  diffuseColor;\n"+"uniform vec3  specularColor;\n"+"uniform vec3  emissiveColor;\n"+"uniform float shininess;\n"+"uniform float transparency;\n"+"uniform float ambientIntensity;\n";return shaderPart;};x3dom.shader.twoSidedMaterial=function(){var shaderPart="uniform vec3  backDiffuseColor;\n"+"uniform vec3  backSpecularColor;\n"+"uniform vec3  backEmissiveColor;\n"+"uniform float backShininess;\n"+"uniform float backTransparency;\n"+"uniform float backAmbientIntensity;\n";return shaderPart;};x3dom.shader.fog=function(){var shaderPart="uniform vec3  fogColor;\n"+"uniform float fogType;\n"+"uniform float fogRange;\n"+"varying vec3 fragEyePosition;\n"+"float calcFog(in vec3 eye) {\n"+"   float f0 = 0.0;\n"+"   if(fogType == 0.0) {\n"+"       if(length(eye) < fogRange){\n"+"           f0 = (fogRange-length(eye)) / fogRange;\n"+"       }\n"+"   }else{\n"+"       if(length(eye) < fogRange){\n"+"           f0 = exp(-length(eye) / (fogRange-length(eye) ) );\n"+"       }\n"+"   }\n"+"   f0 = clamp(f0, 0.0, 1.0);\n"+"   return f0;\n"+"}\n";return shaderPart;};x3dom.shader.clipPlanes=function(numClipPlanes){var shaderPart="",c;for(c=0;c<numClipPlanes;c++){shaderPart+="uniform vec4 clipPlane"+c+"_Plane;\n";shaderPart+="uniform float clipPlane"+c+"_CappingStrength;\n";shaderPart+="uniform vec3 clipPlane"+c+"_CappingColor;\n";}
shaderPart+="vec3 calculateClipPlanes() {\n";for(c=0;c<numClipPlanes;c++){shaderPart+="    vec4 clipPlane"+c+" = clipPlane"+c+"_Plane * viewMatrixInverse;\n";shaderPart+="    float dist"+c+" = dot(fragPosition, clipPlane"+c+");\n";}
shaderPart+="    if( ";for(c=0;c<numClipPlanes;c++){if(c!=0){shaderPart+=" || ";}
shaderPart+="dist"+c+" < 0.0";}
shaderPart+=" ) ";shaderPart+="{ discard; }\n";for(c=0;c<numClipPlanes;c++){shaderPart+="    if( abs(dist"+c+") < clipPlane"+c+"_CappingStrength ) ";shaderPart+="{ return clipPlane"+c+"_CappingColor; }\n";}
shaderPart+="    return vec3(-1.0, -1.0, -1.0);\n";shaderPart+="}\n";return shaderPart;};x3dom.shader.gammaCorrectionDecl=function(properties){var shaderPart="";if(properties.GAMMACORRECTION==="none"){}else if(properties.GAMMACORRECTION==="fastlinear"){shaderPart+="vec4 gammaEncode(vec4 color){\n"+"  vec4 tmp = sqrt(color);\n"+"  return vec4(tmp.rgb, color.a);\n"+"}\n";shaderPart+="vec4 gammaDecode(vec4 color){\n"+"  vec4 tmp = color * color;\n"+"  return vec4(tmp.rgb, color.a);\n"+"}\n";shaderPart+="vec3 gammaEncode(vec3 color){\n"+"  return sqrt(color);\n"+"}\n";shaderPart+="vec3 gammaDecode(vec3 color){\n"+"  return (color * color);\n"+"}\n";}else{shaderPart+="const vec4 gammaEncode4Vector = vec4(0.4545454545454545, 0.4545454545454545, 0.4545454545454545, 1.0);\n";shaderPart+="const vec4 gammaDecode4Vector = vec4(2.2, 2.2, 2.2, 1.0);\n";shaderPart+="vec4 gammaEncode(vec4 color){\n"+"    return pow(color, gammaEncode4Vector);\n"+"}\n";shaderPart+="vec4 gammaDecode(vec4 color){\n"+"    return pow(color, gammaDecode4Vector);\n"+"}\n";shaderPart+="const vec3 gammaEncode3Vector = vec3(0.4545454545454545, 0.4545454545454545, 0.4545454545454545);\n";shaderPart+="const vec3 gammaDecode3Vector = vec3(2.2, 2.2, 2.2);\n";shaderPart+="vec3 gammaEncode(vec3 color){\n"+"    return pow(color, gammaEncode3Vector);\n"+"}\n";shaderPart+="vec3 gammaDecode(vec3 color){\n"+"    return pow(color, gammaDecode3Vector);\n"+"}\n";}
return shaderPart;};x3dom.shader.encodeGamma=function(properties,expr){if(properties.GAMMACORRECTION==="none"){return expr;}else{return"gammaEncode ("+expr+")";}};x3dom.shader.decodeGamma=function(properties,expr){if(properties.GAMMACORRECTION==="none"){return expr;}else{return"gammaDecode ("+expr+")";}};x3dom.shader.rgbaPacking=function(){var shaderPart="";shaderPart+="vec4 packDepth(float depth){\n"+" depth = (depth + 1.0)*0.5;\n"+" vec4 outVal = vec4(1.0, 255.0, 65025.0, 160581375.0) * depth;\n"+" outVal = fract(outVal);\n"+"   outVal -= outVal.yzww * vec4(1.0/255.0, 1.0/255.0, 1.0/255.0, 0.0);\n"+"   return outVal;\n"+"}\n";shaderPart+="float unpackDepth(vec4 color){\n"+" float depth = dot(color, vec4(1.0, 1.0/255.0, 1.0/65025.0, 1.0/160581375.0));\n"+" return (2.0*depth - 1.0);\n"+"}\n";return shaderPart;};x3dom.shader.shadowRendering=function(){var shaderPart="";shaderPart+="float getLightInfluence(float lType, float lShadowIntensity, float lOn, vec3 lLocation, vec3 lDirection, "+"float lCutOffAngle, float lBeamWidth, vec3 lAttenuation, float lRadius, vec3 eyeCoords) {\n"+" if (lOn == 0.0 || lShadowIntensity == 0.0){ return 0.0;\n"+" } else if (lType == 0.0) {\n"+"  return 1.0;\n"+" } else {\n"+"    float attenuation = 0.0;\n"+"    vec3 lightVec = (lLocation - (eyeCoords));\n"+"    float distance = length(lightVec);\n"+"  lightVec = normalize(lightVec);\n"+"  eyeCoords = normalize(-eyeCoords);\n"+"    if(lRadius == 0.0 || distance <= lRadius) {\n"+"        attenuation = 1.0 / max(lAttenuation.x + lAttenuation.y * distance + lAttenuation.z * (distance * distance), 1.0);\n"+"  }\n"+"   if (lType == 1.0) return attenuation;\n"+"    float spotAngle = acos(max(0.0, dot(-lightVec, normalize(lDirection))));\n"+"    if(spotAngle >= lCutOffAngle) return 0.0;\n"+"    else if(spotAngle <= lBeamWidth) return attenuation;\n"+"    else return attenuation * (spotAngle - lCutOffAngle) / (lBeamWidth - lCutOffAngle);\n"+" }\n"+"}\n";shaderPart+="void getShadowValues(inout vec4 shadowMapValues, inout float viewSampleDepth, in mat4 lightMatrix, in vec4 worldCoords, in sampler2D shadowMap){\n"+" vec4 lightSpaceCoords = lightMatrix*worldCoords;\n"+" vec3 lightSpaceCoordsCart = lightSpaceCoords.xyz / lightSpaceCoords.w;\n"+" vec2 textureCoords = (lightSpaceCoordsCart.xy + 1.0)*0.5;\n"+" viewSampleDepth = lightSpaceCoordsCart.z;\n"+" shadowMapValues = texture2D(shadowMap, textureCoords);\n";if(!x3dom.caps.FP_TEXTURES||x3dom.caps.MOBILE)
shaderPart+=" shadowMapValues = vec4(1.0,1.0,unpackDepth(shadowMapValues),1.0);\n";shaderPart+="}\n";shaderPart+="void getShadowValuesPointLight(inout vec4 shadowMapValues, inout float viewSampleDepth, in vec3 lLocation, in vec4 worldCoords, in mat4 lightViewMatrix,"+"in mat4 lMatrix_0, in mat4 lMatrix_1, in mat4 lMatrix_2, in mat4 lMatrix_3, in mat4 lMatrix_4, in mat4 lMatrix_5,"+"in sampler2D shadowMap_0, in sampler2D shadowMap_1, in sampler2D shadowMap_2, in sampler2D shadowMap_3,"+"in sampler2D shadowMap_4, in sampler2D shadowMap_5){\n"+" vec4 transformed = lightViewMatrix * worldCoords;\n"+" vec3 lightVec = normalize(transformed.xyz/transformed.w);\n"+" vec3 lightVecAbs = abs(lightVec);\n"+" float maximum = max(max(lightVecAbs.x, lightVecAbs.y),lightVecAbs.z);\n"+" if (lightVecAbs.x == maximum) {\n"+"  if (lightVec.x < 0.0) getShadowValues(shadowMapValues, viewSampleDepth, lMatrix_3,worldCoords,shadowMap_3);\n"+"  else getShadowValues(shadowMapValues, viewSampleDepth, lMatrix_1,worldCoords,shadowMap_1);\n"+" }\n"+" else if (lightVecAbs.y == maximum) {\n"+"  if (lightVec.y < 0.0) getShadowValues(shadowMapValues, viewSampleDepth, lMatrix_4,worldCoords,shadowMap_4);\n"+"  else getShadowValues(shadowMapValues, viewSampleDepth, lMatrix_5,worldCoords,shadowMap_5);\n"+" }\n"+" else if (lightVec.z < 0.0) getShadowValues(shadowMapValues, viewSampleDepth, lMatrix_0,worldCoords,shadowMap_0);\n"+" else getShadowValues(shadowMapValues, viewSampleDepth, lMatrix_2,worldCoords,shadowMap_2);\n"+"}\n";shaderPart+="void getShadowValuesCascaded(inout vec4 shadowMapValues, inout float viewSampleDepth, in vec4 worldCoords, in float eyeDepth, in mat4 lMatrix_0, in mat4 lMatrix_1, in mat4 lMatrix_2,"+"in mat4 lMatrix_3, in mat4 lMatrix_4, in mat4 lMatrix_5, in sampler2D shadowMap_0, in sampler2D shadowMap_1, in sampler2D shadowMap_2,"+"in sampler2D shadowMap_3, in sampler2D shadowMap_4, in sampler2D shadowMap_5, in float split_0, in float split_1, in float split_2, in float split_3, in float split_4){\n"+" if (eyeDepth < split_0) getShadowValues(shadowMapValues, viewSampleDepth, lMatrix_0, worldCoords, shadowMap_0);\n"+" else if (eyeDepth < split_1) getShadowValues(shadowMapValues, viewSampleDepth, lMatrix_1, worldCoords, shadowMap_1);\n"+" else if (eyeDepth < split_2) getShadowValues(shadowMapValues, viewSampleDepth, lMatrix_2, worldCoords, shadowMap_2);\n"+" else if (eyeDepth < split_3) getShadowValues(shadowMapValues, viewSampleDepth, lMatrix_3, worldCoords, shadowMap_3);\n"+" else if (eyeDepth < split_4) getShadowValues(shadowMapValues, viewSampleDepth, lMatrix_4, worldCoords, shadowMap_4);\n"+" else getShadowValues(shadowMapValues, viewSampleDepth, lMatrix_5, worldCoords, shadowMap_5);\n"+"}\n";shaderPart+="float ESM(float shadowMapDepth, float viewSampleDepth, float offset){\n";if(!x3dom.caps.FP_TEXTURES||x3dom.caps.MOBILE)
shaderPart+=" return exp(-80.0*(1.0-offset)*(viewSampleDepth - shadowMapDepth));\n";else shaderPart+=" return shadowMapDepth * exp(-80.0*(1.0-offset)*viewSampleDepth);\n";shaderPart+="}\n";shaderPart+="float VSM(vec2 moments, float viewSampleDepth, float offset){\n"+" viewSampleDepth = (viewSampleDepth + 1.0) * 0.5;\n"+" if (viewSampleDepth <= moments.x) return 1.0;\n"+" float variance = moments.y - moments.x * moments.x;\n"+" variance = max(variance, 0.00002 + offset*0.01);\n"+" float d = viewSampleDepth - moments.x;\n"+" return variance/(variance + d*d);\n"+"}\n";return shaderPart;};x3dom.shader.light=function(numLights){var shaderPart="";for(var l=0;l<numLights;l++){shaderPart+="uniform float light"+l+"_On;\n"+"uniform float light"+l+"_Type;\n"+"uniform vec3  light"+l+"_Location;\n"+"uniform vec3  light"+l+"_Direction;\n"+"uniform vec3  light"+l+"_Color;\n"+"uniform vec3  light"+l+"_Attenuation;\n"+"uniform float light"+l+"_Radius;\n"+"uniform float light"+l+"_Intensity;\n"+"uniform float light"+l+"_AmbientIntensity;\n"+"uniform float light"+l+"_BeamWidth;\n"+"uniform float light"+l+"_CutOffAngle;\n"+"uniform float light"+l+"_ShadowIntensity;\n";}
shaderPart+="vec3 lighting(in float lType, in vec3 lLocation, in vec3 lDirection, in vec3 lColor, in vec3 lAttenuation, "+"in float lRadius, in float lIntensity, in float lAmbientIntensity, in float lBeamWidth, "+"in float lCutOffAngle, in vec3 N, in vec3 V, float shin, float ambIntensity)\n"+"{\n"+"   vec3 L;\n"+"   float spot = 1.0, attentuation = 0.0;\n"+"   if(lType == 0.0) {\n"+"       L = -normalize(lDirection);\n"+"  V = normalize(V);\n"+"  attentuation = 1.0;\n"+"   } else{\n"+"       L = (lLocation - (-V));\n"+"       float d = length(L);\n"+"  L = normalize(L);\n"+"  V = normalize(V);\n"+"       if(lRadius == 0.0 || d <= lRadius) {\n"+"        attentuation = 1.0 / max(lAttenuation.x + lAttenuation.y * d + lAttenuation.z * (d * d), 1.0);\n"+"  }\n"+"       if(lType == 2.0) {\n"+"           float spotAngle = acos(max(0.0, dot(-L, normalize(lDirection))));\n"+"           if(spotAngle >= lCutOffAngle) spot = 0.0;\n"+"           else if(spotAngle <= lBeamWidth) spot = 1.0;\n"+"           else spot = (spotAngle - lCutOffAngle ) / (lBeamWidth - lCutOffAngle);\n"+"       }\n"+"   }\n"+"   vec3  H = normalize( L + V );\n"+"   float NdotL = clamp(dot(L, N), 0.0, 1.0);\n"+"   float NdotH = clamp(dot(H, N), 0.0, 1.0);\n"+"   float ambientFactor  = lAmbientIntensity * ambIntensity;\n"+"   float diffuseFactor  = lIntensity * NdotL;\n"+"   float specularFactor = lIntensity * pow(NdotH, shin*128.0);\n"+"   return vec3(ambientFactor, diffuseFactor, specularFactor) * attentuation * spot;\n"+"}\n";return shaderPart;};x3dom.shader.TBNCalculation=function(){var shaderPart="";shaderPart+="mat3 cotangent_frame(vec3 N, vec3 p, vec2 uv)\n"+"{\n"+"    // get edge vectors of the pixel triangle\n"+"    vec3 dp1 = dFdx( p );\n"+"    vec3 dp2 = dFdy( p );\n"+"    vec2 duv1 = dFdx( uv );\n"+"    vec2 duv2 = dFdy( uv );\n"+"\n"+"    // solve the linear system\n"+"    vec3 dp2perp = cross( dp2, N );\n"+"    vec3 dp1perp = cross( N, dp1 );\n"+"    vec3 T = dp2perp * duv1.x + dp1perp * duv2.x;\n"+"    vec3 B = dp2perp * duv1.y + dp1perp * duv2.y;\n"+"\n"+"    // construct a scale-invariant frame\n"+"    float invmax = inversesqrt( max( dot(T,T), dot(B,B) ) );\n"+"    return mat3( T * invmax, B * invmax, N );\n"+"}\n\n";shaderPart+="vec3 perturb_normal( vec3 N, vec3 V, vec2 texcoord )\n"+"{\n"+"    // assume N, the interpolated vertex normal and\n"+"    // V, the view vector (vertex to eye)\n"+"    vec3 map = texture2D(normalMap, texcoord ).xyz;\n"+"    map = 2.0 * map - 1.0;\n"+"    mat3 TBN = cotangent_frame(N, -V, texcoord);\n"+"    return normalize(TBN * map);\n"+"}\n\n";return shaderPart;};x3dom.shader.DynamicShader=function(gl,properties)
{this.program=gl.createProgram();var vertexShader=this.generateVertexShader(gl,properties);var fragmentShader=this.generateFragmentShader(gl,properties);gl.attachShader(this.program,vertexShader);gl.attachShader(this.program,fragmentShader);gl.bindAttribLocation(this.program,0,"position");gl.linkProgram(this.program);return this.program;};x3dom.shader.DynamicShader.prototype.generateVertexShader=function(gl,properties)
{var shader="";shader+="uniform mat4 modelViewMatrix;\n";shader+="uniform mat4 modelViewProjectionMatrix;\n";if(properties.POSCOMPONENTS==3){shader+="attribute vec3 position;\n";}else if(properties.POSCOMPONENTS==4){shader+="attribute vec4 position;\n";}
if(properties.IMAGEGEOMETRY){shader+="uniform vec3 IG_bboxMin;\n";shader+="uniform vec3 IG_bboxMax;\n";shader+="uniform float IG_coordTextureWidth;\n";shader+="uniform float IG_coordTextureHeight;\n";shader+="uniform vec2 IG_implicitMeshSize;\n";for(var i=0;i<properties.IG_PRECISION;i++){shader+="uniform sampler2D IG_coords"+i+"\n;";}
if(properties.IG_INDEXED){shader+="uniform sampler2D IG_index;\n";shader+="uniform float IG_indexTextureWidth;\n";shader+="uniform float IG_indexTextureHeight;\n";}}
if(properties.POPGEOMETRY){shader+="uniform float PG_precisionLevel;\n";shader+="uniform float PG_powPrecision;\n";shader+="uniform vec3 PG_maxBBSize;\n";shader+="uniform vec3 PG_bbMin;\n";shader+="uniform vec3 PG_bbMaxModF;\n";shader+="uniform vec3 PG_bboxShiftVec;\n";shader+="uniform float PG_numAnchorVertices;\n";shader+="attribute float PG_vertexID;\n";}
if(properties.LIGHTS){if(properties.NORMALMAP&&properties.NORMALSPACE=="OBJECT"){}else{shader+="varying vec3 fragNormal;\n";shader+="uniform mat4 normalMatrix;\n";if(properties.IMAGEGEOMETRY){shader+="uniform sampler2D IG_normals;\n";}else{if(properties.NORCOMPONENTS==2){if(properties.POSCOMPONENTS!=4){shader+="attribute vec2 normal;\n";}}else if(properties.NORCOMPONENTS==3){shader+="attribute vec3 normal;\n";}}}}
if(properties.VERTEXCOLOR){if(properties.IMAGEGEOMETRY){shader+="uniform sampler2D IG_colors;\n";if(properties.COLCOMPONENTS==3){shader+="varying vec3 fragColor;\n";}else if(properties.COLCOMPONENTS==4){shader+="varying vec4 fragColor;\n";}}else{if(properties.COLCOMPONENTS==3){shader+="attribute vec3 color;\n";shader+="varying vec3 fragColor;\n";}else if(properties.COLCOMPONENTS==4){shader+="attribute vec4 color;\n";shader+="varying vec4 fragColor;\n";}}}
if(properties.TEXTURED){shader+="varying vec2 fragTexcoord;\n";if(!properties.SPHEREMAPPING){if(properties.IMAGEGEOMETRY){shader+="uniform sampler2D IG_texCoords;\n";}else if(!properties.IS_PARTICLE){shader+="attribute vec2 texcoord;\n";}}
if(properties.TEXTRAFO){shader+="uniform mat4 texTrafoMatrix;\n";}
if(properties.NORMALMAP&&properties.NORMALSPACE=="TANGENT"&&!x3dom.caps.STD_DERIVATIVES){x3dom.debug.logWarning("Your System doesn't support the 'OES_STANDARD_DERIVATIVES' Extension. "+"You must set tangents and binormals manually via the FloatVertexAttribute-Node "+"to use normal maps");shader+="attribute vec3 tangent;\n";shader+="attribute vec3 binormal;\n";shader+="varying vec3 fragTangent;\n";shader+="varying vec3 fragBinormal;\n";}
if(properties.CUBEMAP){shader+="varying vec3 fragViewDir;\n";shader+="uniform mat4 viewMatrix;\n";}
if(properties.DISPLACEMENTMAP){shader+="uniform sampler2D displacementMap;\n";shader+="uniform float displacementFactor;\n";shader+="uniform float displacementWidth;\n";shader+="uniform float displacementHeight;\n";shader+="uniform float displacementAxis;\n";}
if(properties.DIFFPLACEMENTMAP){shader+="uniform sampler2D diffuseDisplacementMap;\n";shader+="uniform float displacementFactor;\n";shader+="uniform float displacementWidth;\n";shader+="uniform float displacementHeight;\n";shader+="uniform float displacementAxis;\n";}}
if(properties.VERTEXID){shader+="attribute float id;\n";shader+="varying float fragID;\n";}
if(properties.IS_PARTICLE){shader+="attribute vec3 particleSize;\n";}
if(properties.LIGHTS||properties.FOG||properties.CLIPPLANES){shader+="uniform vec3 eyePosition;\n";shader+="varying vec4 fragPosition;\n";if(properties.FOG){shader+="varying vec3 fragEyePosition;\n";}}
if(properties.REQUIREBBOX){shader+="uniform vec3 bgCenter;\n";shader+="uniform vec3 bgSize;\n";shader+="uniform float bgPrecisionMax;\n";}
if(properties.REQUIREBBOXNOR){shader+="uniform float bgPrecisionNorMax;\n";}
if(properties.REQUIREBBOXCOL){shader+="uniform float bgPrecisionColMax;\n";}
if(properties.REQUIREBBOXTEX){shader+="uniform float bgPrecisionTexMax;\n";}
shader+="void main(void) {\n";if(properties.IMAGEGEOMETRY){if(properties.IG_INDEXED){shader+="vec2 halfPixel = vec2(0.5/IG_indexTextureWidth,0.5/IG_indexTextureHeight);\n";shader+="vec2 IG_texCoord = vec2(position.x*(IG_implicitMeshSize.x/IG_indexTextureWidth), position.y*(IG_implicitMeshSize.y/IG_indexTextureHeight)) + halfPixel;\n";shader+="vec2 IG_indices = texture2D( IG_index, IG_texCoord ).rg;\n";shader+="halfPixel = vec2(0.5/IG_coordTextureWidth,0.5/IG_coordTextureHeight);\n";shader+="IG_texCoord = (IG_indices * 0.996108948) + halfPixel;\n";}else{shader+="vec2 halfPixel = vec2(0.5/IG_coordTextureWidth, 0.5/IG_coordTextureHeight);\n";shader+="vec2 IG_texCoord = vec2(position.x*(IG_implicitMeshSize.x/IG_coordTextureWidth), position.y*(IG_implicitMeshSize.y/IG_coordTextureHeight)) + halfPixel;\n";}
shader+="vec3 temp = vec3(0.0, 0.0, 0.0);\n";shader+="vec3 vertPosition = vec3(0.0, 0.0, 0.0);\n";for(var i=0;i<properties.IG_PRECISION;i++){shader+="temp = 255.0 * texture2D( IG_coords"+i+", IG_texCoord ).rgb;\n";shader+="vertPosition *= 256.0;\n";shader+="vertPosition += temp;\n";}
shader+="vertPosition /= (pow(2.0, 8.0 * "+properties.IG_PRECISION+".0) - 1.0);\n";shader+="vertPosition = vertPosition * (IG_bboxMax - IG_bboxMin) + IG_bboxMin;\n";if(properties.LIGHTS){shader+="vec3 vertNormal = texture2D( IG_normals, IG_texCoord ).rgb;\n";shader+="vertNormal = vertNormal * 2.0 - 1.0;\n";}
if(properties.VERTEXCOLOR){if(properties.COLCOMPONENTS==3){shader+="fragColor = texture2D( IG_colors, IG_texCoord ).rgb;\n";}else if(properties.COLCOMPONENTS==4){shader+="fragColor = texture2D( IG_colors, IG_texCoord ).rgba;\n";}}
if(properties.TEXTURED){shader+="vec4 IG_doubleTexCoords = texture2D( IG_texCoords, IG_texCoord );\n";shader+="vec2 vertTexCoord;";shader+="vertTexCoord.r = (IG_doubleTexCoords.r * 0.996108948) + (IG_doubleTexCoords.b * 0.003891051);\n";shader+="vertTexCoord.g = (IG_doubleTexCoords.g * 0.996108948) + (IG_doubleTexCoords.a * 0.003891051);\n";}}else{shader+="vec3 vertPosition = position.xyz;\n";if(properties.POPGEOMETRY){shader+="vec3 offsetVec = step(vertPosition / bgPrecisionMax, PG_bbMaxModF) * PG_bboxShiftVec;\n";shader+="if ((PG_precisionLevel <= 2.0) || PG_vertexID >= PG_numAnchorVertices) {\n";shader+="   vertPosition = floor(vertPosition / PG_powPrecision) * PG_powPrecision;\n";shader+="   vertPosition /= (65536.0 - PG_powPrecision);\n";shader+="}\n";shader+="else {\n";shader+="   vertPosition /= bgPrecisionMax;\n";shader+="}\n";shader+="vertPosition = (vertPosition + offsetVec + PG_bbMin) * PG_maxBBSize;\n";}
else if(properties.REQUIREBBOX){shader+="vertPosition = bgCenter + bgSize * vertPosition / bgPrecisionMax;\n";}
if(properties.LIGHTS){if(properties.NORCOMPONENTS==2){if(properties.POSCOMPONENTS==4){shader+="vec3 vertNormal = vec3(position.w / 256.0); \n";shader+="vertNormal.x = floor(vertNormal.x) / 255.0; \n";shader+="vertNormal.y = fract(vertNormal.y) * 1.00392156862745; \n";}
else if(properties.REQUIREBBOXNOR){shader+="vec3 vertNormal = vec3(normal.xy, 0.0) / bgPrecisionNorMax;\n";}
shader+="vec2 thetaPhi = 3.14159265358979 * vec2(vertNormal.x, vertNormal.y*2.0-1.0); \n";shader+="vec4 sinCosThetaPhi = sin( vec4(thetaPhi, thetaPhi + 1.5707963267949) ); \n";shader+="vertNormal.x = sinCosThetaPhi.x * sinCosThetaPhi.w; \n";shader+="vertNormal.y = sinCosThetaPhi.x * sinCosThetaPhi.y; \n";shader+="vertNormal.z = sinCosThetaPhi.z; \n";}else{if(properties.NORMALMAP&&properties.NORMALSPACE=="OBJECT"){}else{shader+="vec3 vertNormal = normal;\n";if(properties.REQUIREBBOXNOR){shader+="vertNormal = vertNormal / bgPrecisionNorMax;\n";}
if(properties.POPGEOMETRY){shader+="vertNormal = 2.0*vertNormal - 1.0;\n";}}}}
if(properties.VERTEXCOLOR){shader+="fragColor = color;\n";if(properties.REQUIREBBOXCOL){shader+="fragColor = fragColor / bgPrecisionColMax;\n";}}
if((properties.TEXTURED)&&!properties.SPHEREMAPPING){if(properties.IS_PARTICLE){shader+="vec2 vertTexCoord = vec2(0.0);\n";}
else{shader+="vec2 vertTexCoord = texcoord;\n";if(properties.REQUIREBBOXTEX){shader+="vertTexCoord = vertTexCoord / bgPrecisionTexMax;\n";}}}}
if(properties.LIGHTS){if((properties.DISPLACEMENTMAP||properties.DIFFPLACEMENTMAP)&&!properties.NORMALMAP){shader+="float dx = 1.0 / displacementWidth;\n";shader+="float dy = 1.0 / displacementHeight;\n";if(properties.DISPLACEMENTMAP)
{shader+="float s1 = texture2D(displacementMap, vec2(vertTexCoord.x - dx, 1.0 - vertTexCoord.y)).r;\n";shader+="float s2 = texture2D(displacementMap, vec2(vertTexCoord.x, 1.0 - vertTexCoord.y - dy)).r;\n";shader+="float s3 = texture2D(displacementMap, vec2(vertTexCoord.x + dx, 1.0 - vertTexCoord.y)).r;\n";shader+="float s4 = texture2D(displacementMap, vec2(vertTexCoord.x, 1.0 - vertTexCoord.y + dy)).r;\n";}
else if(properties.DIFFPLACEMENTMAP)
{shader+="float s1 = texture2D(diffuseDisplacementMap, vec2(vertTexCoord.x - dx, 1.0 - vertTexCoord.y)).a;\n";shader+="float s2 = texture2D(diffuseDisplacementMap, vec2(vertTexCoord.x, 1.0 - vertTexCoord.y - dy)).a;\n";shader+="float s3 = texture2D(diffuseDisplacementMap, vec2(vertTexCoord.x + dx, 1.0 - vertTexCoord.y)).a;\n";shader+="float s4 = texture2D(diffuseDisplacementMap, vec2(vertTexCoord.x, 1.0 - vertTexCoord.y + dy)).a;\n";}
shader+="float coef = displacementFactor;\n";shader+="vec3 calcNormal;\n";shader+="if (displacementAxis == 0.0) {\n";shader+="calcNormal = vec3((s1 - s3) * coef, -5.0, (s2 - s4) * coef);\n";shader+="} else if(displacementAxis == 1.0) {\n";shader+="calcNormal = vec3((s1 - s3) * coef, -5.0, (s2 - s4) * coef);\n";shader+="} else {\n";shader+="calcNormal = vec3((s1 - s3) * coef, -(s2 - s4) * coef, 5.0);\n";shader+="}\n";shader+="calcNormal = normalize(calcNormal);\n";shader+="fragNormal = (normalMatrix * vec4(calcNormal, 0.0)).xyz;\n";}
else if(properties.NORMALMAP&&properties.NORMALSPACE=="OBJECT"){}
else
{shader+="fragNormal = (normalMatrix * vec4(vertNormal, 0.0)).xyz;\n";}}
if(properties.TEXTURED){if(properties.CUBEMAP){shader+="fragViewDir = (viewMatrix[3].xyz);\n";}
if(properties.SPHEREMAPPING){shader+=" fragTexcoord = 0.5 + fragNormal.xy / 2.0;\n";}else if(properties.TEXTRAFO){shader+=" fragTexcoord = (texTrafoMatrix * vec4(vertTexCoord, 1.0, 1.0)).xy;\n";}else{shader+=" fragTexcoord = vertTexCoord;\n";if(properties.POPGEOMETRY&&x3dom.debug.usePrecisionLevelAsTexCoord===true)
shader+="fragTexcoord = vec2(0.03125 + 0.9375 * (PG_precisionLevel / 16.0), 1.0);";}
if(properties.NORMALMAP&&properties.NORMALSPACE=="TANGENT"&&!x3dom.caps.STD_DERIVATIVES){shader+="fragTangent  = (normalMatrix * vec4(tangent, 0.0)).xyz;\n";shader+="fragBinormal = (normalMatrix * vec4(binormal, 0.0)).xyz;\n";}}
if(properties.LIGHTS||properties.FOG||properties.CLIPPLANES){shader+="fragPosition = (modelViewMatrix * vec4(vertPosition, 1.0));\n";if(properties.FOG){shader+="fragEyePosition = eyePosition - fragPosition.xyz;\n";}}
if(properties.MULTIDIFFALPMAP){shader+="fragID = id;\n";}
if(properties.DISPLACEMENTMAP){shader+="vertPosition += normalize(vertNormal) * texture2D(displacementMap, vec2(fragTexcoord.x, 1.0-fragTexcoord.y)).r * displacementFactor;\n";}
else if(properties.DIFFPLACEMENTMAP)
{shader+="vertPosition += normalize(vertNormal) * texture2D(diffuseDisplacementMap, vec2(fragTexcoord.x, 1.0-fragTexcoord.y)).a * displacementFactor;\n";}
shader+="gl_Position = modelViewProjectionMatrix * vec4(vertPosition, 1.0);\n";if(properties.IS_PARTICLE){shader+="float spriteDist = (gl_Position.w > 0.000001) ? gl_Position.w : 0.000001;\n";shader+="float pointSize = floor(length(particleSize) * 256.0 / spriteDist + 0.5);\n";shader+="gl_PointSize = clamp(pointSize, 2.0, 256.0);\n";}
else{shader+="gl_PointSize = 2.0;\n";}
shader+="}\n";var vertexShader=gl.createShader(gl.VERTEX_SHADER);gl.shaderSource(vertexShader,shader);gl.compileShader(vertexShader);if(!gl.getShaderParameter(vertexShader,gl.COMPILE_STATUS)){x3dom.debug.logInfo("VERTEX:\n"+shader);x3dom.debug.logError("VertexShader "+gl.getShaderInfoLog(vertexShader));}
return vertexShader;};x3dom.shader.DynamicShader.prototype.generateFragmentShader=function(gl,properties)
{var shader="#ifdef GL_FRAGMENT_PRECISION_HIGH\n";shader+=" precision highp float;\n";shader+="#else\n";shader+=" precision mediump float;\n";shader+="#endif\n\n";shader+="uniform mat4 modelMatrix;\n";shader+="uniform mat4 modelViewMatrix;\n";shader+="uniform mat4 viewMatrixInverse;\n";shader+=x3dom.shader.material();if(properties.TWOSIDEDMAT){shader+=x3dom.shader.twoSidedMaterial();}
if(properties.VERTEXCOLOR){if(properties.COLCOMPONENTS==3){shader+="varying vec3 fragColor;  \n";}else if(properties.COLCOMPONENTS==4){shader+="varying vec4 fragColor;  \n";}}
if(properties.CUBEMAP||properties.CLIPPLANES)
{shader+="uniform mat4 modelViewMatrixInverse;\n";}
if(properties.VERTEXID){shader+="varying float fragID;\n";if(properties.MULTIDIFFALPMAP){shader+="uniform sampler2D multiDiffuseAlphaMap;\n";shader+="uniform float multiDiffuseAlphaWidth;\n";shader+="uniform float multiDiffuseAlphaHeight;\n";}
if(properties.MULTIEMIAMBMAP){shader+="uniform sampler2D multiEmissiveAmbientMap;\n";shader+="uniform float multiEmissiveAmbientWidth;\n";shader+="uniform float multiEmissiveAmbientHeight;\n";}
if(properties.MULTISPECSHINMAP){shader+="uniform sampler2D multiSpecularShininessMap;\n";shader+="uniform float multiSpecularShininessWidth;\n";shader+="uniform float multiSpecularShininessHeight;\n";}
if(properties.MULTIVISMAP){shader+="uniform sampler2D multiVisibilityMap;\n";shader+="uniform float multiVisibilityWidth;\n";shader+="uniform float multiVisibilityHeight;\n";}}
if(properties.TEXTURED){shader+="varying vec2 fragTexcoord;\n";if((properties.TEXTURED||properties.DIFFUSEMAP)){shader+="uniform sampler2D diffuseMap;\n";}
if(properties.CUBEMAP){shader+="uniform samplerCube environmentMap;\n";shader+="varying vec3 fragViewDir;\n";shader+="uniform float environmentFactor;\n";}
if(properties.SPECMAP){shader+="uniform sampler2D specularMap;\n";}
if(properties.SHINMAP){shader+="uniform sampler2D shininessMap;\n";}
if(properties.DISPLACEMENTMAP){shader+="uniform sampler2D displacementMap;\n";shader+="uniform float displacementWidth;\n";shader+="uniform float displacementHeight;\n";}
if(properties.DIFFPLACEMENTMAP){shader+="uniform sampler2D diffuseDisplacementMap;\n";shader+="uniform float displacementWidth;\n";shader+="uniform float displacementHeight;\n";}
if(properties.NORMALMAP){shader+="uniform sampler2D normalMap;\n";if(properties.NORMALSPACE=="TANGENT"){if(x3dom.caps.STD_DERIVATIVES){shader+="#extension GL_OES_standard_derivatives:enable\n";shader+=x3dom.shader.TBNCalculation();}else{shader+="varying vec3 fragTangent;\n";shader+="varying vec3 fragBinormal;\n";}}else if(properties.NORMALSPACE=="OBJECT"){shader+="uniform mat4 normalMatrix;\n";}}}
if(properties.FOG){shader+=x3dom.shader.fog();}
if(properties.LIGHTS||properties.CLIPPLANES)
{shader+="varying vec4 fragPosition;\n";shader+="uniform float isOrthoView;\n";}
if(properties.LIGHTS){if(properties.NORMALMAP&&properties.NORMALSPACE=="OBJECT"){}else{shader+="varying vec3 fragNormal;\n";}
shader+=x3dom.shader.light(properties.LIGHTS);}
if(properties.CLIPPLANES){shader+=x3dom.shader.clipPlanes(properties.CLIPPLANES);}
shader+=x3dom.shader.gammaCorrectionDecl(properties);shader+="void main(void) {\n";if(properties.CLIPPLANES)
{shader+="vec3 cappingColor = calculateClipPlanes();\n";}
shader+="vec4 color;\n";shader+="color.rgb = "+x3dom.shader.decodeGamma(properties,"diffuseColor")+";\n";shader+="color.a = 1.0 - transparency;\n";shader+="vec3 _emissiveColor     = emissiveColor;\n";shader+="float _shininess        = shininess;\n";shader+="vec3 _specularColor     = specularColor;\n";shader+="float _ambientIntensity = ambientIntensity;\n";shader+="float _transparency     = transparency;\n";if(properties.MULTIVISMAP||properties.MULTIDIFFALPMAP||properties.MULTISPECSHINMAP||properties.MULTIEMIAMBMAP){shader+="vec2 idCoord;\n";shader+="float roundedIDVisibility = floor(fragID+0.5);\n";shader+="float roundedIDMaterial = floor(fragID+0.5);\n";shader+="if(!gl_FrontFacing) {\n";shader+="    roundedIDMaterial = floor(fragID + (multiDiffuseAlphaWidth*multiDiffuseAlphaWidth) + 0.5);\n";shader+="}\n";}
if(properties.MULTIVISMAP){shader+="idCoord.x = mod(roundedIDVisibility, multiVisibilityWidth) * (1.0 / multiVisibilityWidth) + (0.5 / multiVisibilityWidth);\n";shader+="idCoord.y = floor(roundedIDVisibility / multiVisibilityWidth) * (1.0 / multiVisibilityHeight) + (0.5 / multiVisibilityHeight);\n";shader+="vec4 visibility = texture2D( multiVisibilityMap, idCoord );\n";shader+="if (visibility.r < 1.0) discard; \n";}
if(properties.MULTIDIFFALPMAP){shader+="idCoord.x = mod(roundedIDMaterial, multiDiffuseAlphaWidth) * (1.0 / multiDiffuseAlphaWidth) + (0.5 / multiDiffuseAlphaWidth);\n";shader+="idCoord.y = floor(roundedIDMaterial / multiDiffuseAlphaWidth) * (1.0 / multiDiffuseAlphaHeight) + (0.5 / multiDiffuseAlphaHeight);\n";shader+="vec4 diffAlpha = texture2D( multiDiffuseAlphaMap, idCoord );\n";shader+="color.rgb = "+x3dom.shader.decodeGamma(properties,"diffAlpha.rgb")+";\n";shader+="_transparency = 1.0 - diffAlpha.a;\n";shader+="color.a = diffAlpha.a;\n";}
if(properties.MULTIEMIAMBMAP){shader+="idCoord.x = mod(roundedIDMaterial, multiDiffuseAlphaWidth) * (1.0 / multiDiffuseAlphaWidth) + (0.5 / multiDiffuseAlphaWidth);\n";shader+="idCoord.y = floor(roundedIDMaterial / multiDiffuseAlphaWidth) * (1.0 / multiDiffuseAlphaHeight) + (0.5 / multiDiffuseAlphaHeight);\n";shader+="vec4 emiAmb = texture2D( multiEmissiveAmbientMap, idCoord );\n";shader+="_emissiveColor = emiAmb.rgb;\n";shader+="_ambientIntensity = emiAmb.a;\n";}
if(properties.VERTEXCOLOR){if(properties.COLCOMPONENTS===3){shader+="color.rgb = "+x3dom.shader.decodeGamma(properties,"fragColor")+";\n";}else if(properties.COLCOMPONENTS===4){shader+="color = "+x3dom.shader.decodeGamma(properties,"fragColor")+";\n";}}
if(properties.LIGHTS){shader+="vec3 ambient   = vec3(0.0, 0.0, 0.0);\n";shader+="vec3 diffuse   = vec3(0.0, 0.0, 0.0);\n";shader+="vec3 specular  = vec3(0.0, 0.0, 0.0);\n";shader+="vec3 eye;\n";shader+="if ( isOrthoView > 0.0 ) {\n";shader+="    eye = vec3(0.0, 0.0, 1.0);\n";shader+="} else {\n";shader+="    eye = -fragPosition.xyz;\n";shader+="}\n";if(properties.NORMALMAP&&properties.NORMALSPACE=="OBJECT"){shader+="vec3 normal  = vec3(0.0, 0.0, 0.0);\n";}else{shader+="vec3 normal    = normalize(fragNormal);\n";}
if(properties.NORMALMAP){if(properties.NORMALSPACE=="TANGENT"){shader+="vec3 n = normal;\n";if(x3dom.caps.STD_DERIVATIVES){shader+="normal = perturb_normal( n, fragPosition.xyz, vec2(fragTexcoord.x, 1.0-fragTexcoord.y) );\n";}else{shader+="vec3 t = normalize( fragTangent );\n";shader+="vec3 b = normalize( fragBinormal );\n";shader+="mat3 tangentToWorld = mat3(t, b, n);\n";shader+="normal = texture2D( normalMap, vec2(fragTexcoord.x, 1.0-fragTexcoord.y) ).rgb;\n";shader+="normal = 2.0 * normal - 1.0;\n";shader+="normal = normalize( normal * tangentToWorld );\n";shader+="normal.y = -normal.y;\n";shader+="normal.x = -normal.x;\n";}}else if(properties.NORMALSPACE=="OBJECT"){shader+="normal = texture2D( normalMap, vec2(fragTexcoord.x, 1.0-fragTexcoord.y) ).rgb;\n";shader+="normal = 2.0 * normal - 1.0;\n";shader+="normal = (normalMatrix * vec4(normal, 0.0)).xyz;\n";shader+="normal = normalize(normal);\n";}}
if(properties.SHINMAP){shader+="_shininess = texture2D( shininessMap, vec2(fragTexcoord.x, 1.0-fragTexcoord.y) ).r;\n";}
if(properties.SPECMAP){shader+="_specularColor = "+x3dom.shader.decodeGamma(properties,"texture2D(specularMap, vec2(fragTexcoord.x, 1.0-fragTexcoord.y)).rgb")+";\n";}
if(properties.MULTISPECSHINMAP){shader+="idCoord.x = mod(roundedIDMaterial, multiSpecularShininessWidth) * (1.0 / multiSpecularShininessWidth) + (0.5 / multiSpecularShininessWidth);\n";shader+="idCoord.y = floor(roundedIDMaterial / multiSpecularShininessWidth) * (1.0 / multiSpecularShininessHeight) + (0.5 / multiSpecularShininessHeight);\n";shader+="vec4 specShin = texture2D( multiSpecularShininessMap, idCoord );\n";shader+="_specularColor = specShin.rgb;\n";shader+="_shininess = specShin.a;\n";}
if(!properties.SOLID||properties.TWOSIDEDMAT){shader+="if (dot(normal, eye) < 0.0) {\n";shader+="  normal *= -1.0;\n";shader+="}\n";}
if(properties.SEPARATEBACKMAT){shader+="  if(!gl_FrontFacing) {\n";shader+="    color.rgb = "+x3dom.shader.decodeGamma(properties,"backDiffuseColor")+";\n";shader+="    color.a = 1.0 - backTransparency;\n";shader+="    _transparency = 1.0 - backTransparency;\n";shader+="    _shininess = backShininess;\n";shader+="    _emissiveColor = backEmissiveColor;\n";shader+="    _specularColor = backSpecularColor;\n";shader+="    _ambientIntensity = backAmbientIntensity;\n";shader+="  }\n";}
if(properties.LIGHTS){shader+="vec3 ads;\n";for(var l=0;l<properties.LIGHTS;l++){var lightCol="light"+l+"_Color";shader+="ads = lighting(light"+l+"_Type, "+"light"+l+"_Location, "+"light"+l+"_Direction, "+
lightCol+", "+"light"+l+"_Attenuation, "+"light"+l+"_Radius, "+"light"+l+"_Intensity, "+"light"+l+"_AmbientIntensity, "+"light"+l+"_BeamWidth, "+"light"+l+"_CutOffAngle, "+"normal, eye, _shininess, _ambientIntensity);\n";shader+="   ambient  += "+lightCol+" * ads.r;\n"+"   diffuse  += "+lightCol+" * ads.g;\n"+"   specular += "+lightCol+" * ads.b;\n";}
shader+="ambient = max(ambient, 0.0);\n";shader+="diffuse = max(diffuse, 0.0);\n";shader+="specular = max(specular, 0.0);\n";}
if(properties.TEXTURED&&(properties.DIFFUSEMAP||properties.DIFFPLACEMENTMAP||properties.TEXT||properties.CUBEMAP)){if(properties.CUBEMAP){shader+="vec3 viewDir = normalize(fragViewDir);\n";shader+="vec3 reflected = reflect(-eye, normal);\n";shader+="reflected = (modelViewMatrixInverse * vec4(reflected, 0.0)).xyz;\n";shader+="vec4 envColor = "+x3dom.shader.decodeGamma(properties,"textureCube(environmentMap, reflected)")+";\n";shader+="color.a *= envColor.a;\n";}
if(properties.DIFFPLACEMENTMAP)
{shader+="vec2 texCoord = vec2(fragTexcoord.x, 1.0-fragTexcoord.y);\n";shader+="vec4 texColor = texture2D(diffuseDisplacementMap, texCoord);\n";}
else if(properties.DIFFUSEMAP||properties.TEXT)
{if(properties.PIXELTEX){shader+="vec2 texCoord = fragTexcoord;\n";}else{shader+="vec2 texCoord = vec2(fragTexcoord.x, 1.0-fragTexcoord.y);\n";}
shader+="vec4 texColor = "+x3dom.shader.decodeGamma(properties,"texture2D(diffuseMap, texCoord)")+";\n";shader+="color.a *= texColor.a;\n";}
if(properties.BLENDING){shader+="color.rgb = (_emissiveColor + max(ambient + diffuse, 0.0) * color.rgb + specular*_specularColor);\n";if(properties.DIFFUSEMAP||properties.TEXT){shader+="color.rgb *= texColor.rgb;\n";}
if(properties.CUBEMAP){shader+="color.rgb *= mix(vec3(1.0,1.0,1.0), envColor.rgb, environmentFactor);\n";}}else{shader+="color.rgb = (_emissiveColor + max(ambient + diffuse, 0.0) * texColor.rgb + specular*_specularColor);\n";}}else{shader+="color.rgb = (_emissiveColor + max(ambient + diffuse, 0.0) * color.rgb + specular*_specularColor);\n";}}else{if(properties.APPMAT&&!properties.VERTEXCOLOR&&!properties.TEXTURED){shader+="color = vec4(0.0, 0.0, 0.0, 1.0 - _transparency);\n";}
if(properties.TEXTURED&&(properties.DIFFUSEMAP||properties.DIFFPLACEMENTMAP||properties.TEXT)){if(properties.PIXELTEX){shader+="vec2 texCoord = fragTexcoord;\n";}else{shader+="vec2 texCoord = vec2(fragTexcoord.x, 1.0-fragTexcoord.y);\n";}
if(properties.IS_PARTICLE){shader+="texCoord = clamp(gl_PointCoord, 0.01, 0.99);\n";}
shader+="vec4 texColor = "+x3dom.shader.decodeGamma(properties,"texture2D(diffuseMap, texCoord)")+";\n";shader+="color.a = texColor.a;\n";if(properties.BLENDING||properties.IS_PARTICLE){shader+="color.rgb += _emissiveColor.rgb;\n";shader+="color.rgb *= texColor.rgb;\n";}else{shader+="color = texColor;\n";}}else if(!properties.VERTEXCOLOR&&!properties.POINTLINE2D){shader+="color.rgb += _emissiveColor;\n";}else if(!properties.VERTEXCOLOR&&properties.POINTLINE2D){shader+="color.rgb = _emissiveColor;\n";if(properties.IS_PARTICLE){shader+="float pAlpha = 1.0 - clamp(length((gl_PointCoord - 0.5) * 2.0), 0.0, 1.0);\n";shader+="color.rgb *= vec3(pAlpha);\n";shader+="color.a = pAlpha;\n";}}else if(properties.IS_PARTICLE){shader+="float pAlpha = 1.0 - clamp(length((gl_PointCoord - 0.5) * 2.0), 0.0, 1.0);\n";shader+="color.rgb *= vec3(pAlpha);\n";shader+="color.a = pAlpha;\n";}}
if(properties.CLIPPLANES)
{shader+="if (cappingColor.r != -1.0) {\n";shader+="    color.rgb = cappingColor;\n";shader+="}\n";}
if(properties.TEXT){shader+="if (color.a <= 0.5) discard;\n";}else{shader+="if (color.a <= "+properties.ALPHATHRESHOLD+") discard;\n";}
shader+="color = clamp(color, 0.0, 1.0);\n";shader+="color = "+x3dom.shader.encodeGamma(properties,"color")+";\n";if(properties.FOG){shader+="float f0 = calcFog(fragEyePosition);\n";shader+="color.rgb = fogColor * (1.0-f0) + f0 * (color.rgb);\n";}
shader+="gl_FragColor = color;\n";shader+="}\n";var fragmentShader=gl.createShader(gl.FRAGMENT_SHADER);gl.shaderSource(fragmentShader,shader);gl.compileShader(fragmentShader);if(!gl.getShaderParameter(fragmentShader,gl.COMPILE_STATUS)){x3dom.debug.logInfo("FRAGMENT:\n"+shader);x3dom.debug.logError("FragmentShader "+gl.getShaderInfoLog(fragmentShader));}
return fragmentShader;};x3dom.shader.DynamicMobileShader=function(gl,properties)
{this.program=gl.createProgram();var vertexShader=this.generateVertexShader(gl,properties);var fragmentShader=this.generateFragmentShader(gl,properties);gl.attachShader(this.program,vertexShader);gl.attachShader(this.program,fragmentShader);gl.bindAttribLocation(this.program,0,"position");gl.linkProgram(this.program);return this.program;};x3dom.shader.DynamicMobileShader.prototype.generateVertexShader=function(gl,properties)
{var shader="";shader+=x3dom.shader.material();if(properties.TWOSIDEDMAT){shader+=x3dom.shader.twoSidedMaterial();}
shader+="uniform mat4 normalMatrix;\n";shader+="uniform mat4 modelViewMatrix;\n";shader+="uniform mat4 modelViewProjectionMatrix;\n";if(properties.POSCOMPONENTS==3){shader+="attribute vec3 position;\n";}else if(properties.POSCOMPONENTS==4){shader+="attribute vec4 position;\n";}
if(properties.IMAGEGEOMETRY){shader+="uniform vec3 IG_bboxMin;\n";shader+="uniform vec3 IG_bboxMax;\n";shader+="uniform float IG_coordTextureWidth;\n";shader+="uniform float IG_coordTextureHeight;\n";shader+="uniform vec2 IG_implicitMeshSize;\n";for(var i=0;i<properties.IG_PRECISION;i++){shader+="uniform sampler2D IG_coords"+i+"\n;";}
if(properties.IG_INDEXED){shader+="uniform sampler2D IG_index;\n";shader+="uniform float IG_indexTextureWidth;\n";shader+="uniform float IG_indexTextureHeight;\n";}}
if(properties.POPGEOMETRY){shader+="uniform float PG_precisionLevel;\n";shader+="uniform float PG_powPrecision;\n";shader+="uniform vec3 PG_maxBBSize;\n";shader+="uniform vec3 PG_bbMin;\n";shader+="uniform vec3 PG_bbMaxModF;\n";shader+="uniform vec3 PG_bboxShiftVec;\n";shader+="uniform float PG_numAnchorVertices;\n";shader+="attribute float PG_vertexID;\n";}
if(!properties.POINTLINE2D){if(properties.IMAGEGEOMETRY){shader+="uniform sampler2D IG_normals;\n";}else{if(properties.NORCOMPONENTS==2){if(properties.POSCOMPONENTS!=4){shader+="attribute vec2 normal;\n";}}else if(properties.NORCOMPONENTS==3){shader+="attribute vec3 normal;\n";}}}
shader+="varying vec4 fragColor;\n";if(properties.VERTEXCOLOR){if(properties.IMAGEGEOMETRY){shader+="uniform sampler2D IG_colors;";}else{if(properties.COLCOMPONENTS==3){shader+="attribute vec3 color;";}else if(properties.COLCOMPONENTS==4){shader+="attribute vec4 color;";}}}
if(properties.TEXTURED){shader+="varying vec2 fragTexcoord;\n";if(properties.IMAGEGEOMETRY){shader+="uniform sampler2D IG_texCoords;";}else{shader+="attribute vec2 texcoord;\n";}
if(properties.TEXTRAFO){shader+="uniform mat4 texTrafoMatrix;\n";}
if(!properties.BLENDING){shader+="varying vec3 fragAmbient;\n";shader+="varying vec3 fragDiffuse;\n";}
if(properties.CUBEMAP){shader+="varying vec3 fragViewDir;\n";shader+="varying vec3 fragNormal;\n";shader+="uniform mat4 viewMatrix;\n";}}
if(properties.FOG){shader+=x3dom.shader.fog();}
if(properties.LIGHTS){shader+=x3dom.shader.light(properties.LIGHTS);}
if(properties.REQUIREBBOX){shader+="uniform vec3 bgCenter;\n";shader+="uniform vec3 bgSize;\n";shader+="uniform float bgPrecisionMax;\n";}
if(properties.REQUIREBBOXNOR){shader+="uniform float bgPrecisionNorMax;\n";}
if(properties.REQUIREBBOXCOL){shader+="uniform float bgPrecisionColMax;\n";}
if(properties.REQUIREBBOXTEX){shader+="uniform float bgPrecisionTexMax;\n";}
shader+="void main(void) {\n";shader+="gl_PointSize = 2.0;\n";if(properties.IMAGEGEOMETRY){if(properties.IG_INDEXED){shader+="vec2 halfPixel = vec2(0.5/IG_indexTextureWidth,0.5/IG_indexTextureHeight);\n";shader+="vec2 IG_texCoord = vec2(position.x*(IG_implicitMeshSize.x/IG_indexTextureWidth), position.y*(IG_implicitMeshSize.y/IG_indexTextureHeight)) + halfPixel;\n";shader+="vec2 IG_indices = texture2D( IG_index, IG_texCoord ).rg;\n";shader+="halfPixel = vec2(0.5/IG_coordTextureWidth,0.5/IG_coordTextureHeight);\n";shader+="IG_texCoord = (IG_indices * 0.996108948) + halfPixel;\n";}else{shader+="vec2 halfPixel = vec2(0.5/IG_coordTextureWidth, 0.5/IG_coordTextureHeight);\n";shader+="vec2 IG_texCoord = vec2(position.x*(IG_implicitMeshSize.x/IG_coordTextureWidth), position.y*(IG_implicitMeshSize.y/IG_coordTextureHeight)) + halfPixel;\n";}
shader+="vec3 temp = vec3(0.0, 0.0, 0.0);\n";shader+="vec3 vertPosition = vec3(0.0, 0.0, 0.0);\n";for(var i=0;i<properties.IG_PRECISION;i++){shader+="temp = 255.0 * texture2D( IG_coords"+i+", IG_texCoord ).rgb;\n";shader+="vertPosition *= 256.0;\n";shader+="vertPosition += temp;\n";}
shader+="vertPosition /= (pow(2.0, 8.0 * "+properties.IG_PRECISION+".0) - 1.0);\n";shader+="vertPosition = vertPosition * (IG_bboxMax - IG_bboxMin) + IG_bboxMin;\n";if(!properties.POINTLINE2D){shader+="vec3 vertNormal = texture2D( IG_normals, IG_texCoord ).rgb;\n";shader+="vertNormal = vertNormal * 2.0 - 1.0;\n";}
if(properties.VERTEXCOLOR){if(properties.COLCOMPONENTS==3){shader+="vec3 vertColor = texture2D( IG_colors, IG_texCoord ).rgb;";}else if(properties.COLCOMPONENTS==4){shader+="vec4 vertColor = texture2D( IG_colors, IG_texCoord ).rgba;";}}
if(properties.TEXTURED){shader+="vec4 IG_doubleTexCoords = texture2D( IG_texCoords, IG_texCoord );\n";shader+="vec2 vertTexCoord;";shader+="vertTexCoord.r = (IG_doubleTexCoords.r * 0.996108948) + (IG_doubleTexCoords.b * 0.003891051);\n";shader+="vertTexCoord.g = (IG_doubleTexCoords.g * 0.996108948) + (IG_doubleTexCoords.a * 0.003891051);\n";}}else{shader+="vec3 vertPosition = position.xyz;\n";if(properties.POPGEOMETRY){shader+="vec3 offsetVec = step(vertPosition / bgPrecisionMax, PG_bbMaxModF) * PG_bboxShiftVec;\n";shader+="if ((PG_precisionLevel <= 2.0) || PG_vertexID >= PG_numAnchorVertices) {\n";shader+="   vertPosition = floor(vertPosition / PG_powPrecision) * PG_powPrecision;\n";shader+="   vertPosition /= (65536.0 - PG_powPrecision);\n";shader+="}\n";shader+="else {\n";shader+="   vertPosition /= bgPrecisionMax;\n";shader+="}\n";shader+="vertPosition = (vertPosition + offsetVec + PG_bbMin) * PG_maxBBSize;\n";}
else if(properties.REQUIREBBOX){shader+="vertPosition = bgCenter + bgSize * vertPosition / bgPrecisionMax;\n";}
if(!properties.POINTLINE2D){if(properties.NORCOMPONENTS==2){if(properties.POSCOMPONENTS==4){shader+="vec3 vertNormal = vec3(position.w / 256.0); \n";shader+="vertNormal.x = floor(vertNormal.x) / 255.0; \n";shader+="vertNormal.y = fract(vertNormal.y) * 1.00392156862745; \n";}else if(properties.REQUIREBBOXNOR){shader+="vec3 vertNormal = vec3(normal.xy, 0.0) / bgPrecisionNorMax;\n";}else{shader+="vec3 vertNormal = vec3(normal.xy, 0.0);\n";}
shader+="vec2 thetaPhi = 3.14159265358979 * vec2(vertNormal.x, vertNormal.y*2.0-1.0); \n";shader+="vec4 sinCosThetaPhi = vec4(thetaPhi, thetaPhi + 1.5707963267949); \n";shader+="vec4 thetaPhiPow2 = sinCosThetaPhi * sinCosThetaPhi; \n";shader+="vec4 thetaPhiPow3 =  thetaPhiPow2  * sinCosThetaPhi; \n";shader+="vec4 thetaPhiPow5 =  thetaPhiPow3  * thetaPhiPow2; \n";shader+="vec4 thetaPhiPow7 =  thetaPhiPow5  * thetaPhiPow2; \n";shader+="vec4 thetaPhiPow9 =  thetaPhiPow7  * thetaPhiPow2; \n";shader+="sinCosThetaPhi +=  -0.16666666667   * thetaPhiPow3; \n";shader+="sinCosThetaPhi +=   0.00833333333   * thetaPhiPow5; \n";shader+="sinCosThetaPhi +=  -0.000198412698  * thetaPhiPow7; \n";shader+="sinCosThetaPhi +=   0.0000027557319 * thetaPhiPow9; \n";shader+="vertNormal.x = sinCosThetaPhi.x * sinCosThetaPhi.w; \n";shader+="vertNormal.y = sinCosThetaPhi.x * sinCosThetaPhi.y; \n";shader+="vertNormal.z = sinCosThetaPhi.z; \n";}else{shader+="vec3 vertNormal = normal;\n";if(properties.REQUIREBBOXNOR){shader+="vertNormal = vertNormal / bgPrecisionNorMax;\n";}
if(properties.POPGEOMETRY){shader+="vertNormal = 2.0*vertNormal - 1.0;\n";}}}
if(properties.VERTEXCOLOR){if(properties.COLCOMPONENTS==3){shader+="vec3 vertColor = color;";}else if(properties.COLCOMPONENTS==4){shader+="vec4 vertColor = color;";}
if(properties.REQUIREBBOXCOL){shader+="vertColor = vertColor / bgPrecisionColMax;\n";}}
if(properties.TEXTURED){shader+="vec2 vertTexCoord = texcoord;\n";if(properties.REQUIREBBOXTEX){shader+="vertTexCoord = vertTexCoord / bgPrecisionTexMax;\n";}}}
shader+="vec3 positionMV = (modelViewMatrix * vec4(vertPosition, 1.0)).xyz;\n";if(!properties.POINTLINE2D){shader+="vec3 normalMV = normalize( (normalMatrix * vec4(vertNormal, 0.0)).xyz );\n";}
shader+="vec3 eye = -positionMV;\n";if(properties.VERTEXCOLOR){shader+="vec3 rgb = vertColor.rgb;\n";if(properties.COLCOMPONENTS==4){shader+="float alpha = vertColor.a;\n";}else if(properties.COLCOMPONENTS==3){shader+="float alpha = 1.0 - transparency;\n";}}else{shader+="vec3 rgb = diffuseColor;\n";shader+="float alpha = 1.0 - transparency;\n";}
if(properties.TEXTURED){if(properties.CUBEMAP){shader+="fragViewDir = viewMatrix[3].xyz;\n";shader+="fragNormal = normalMV;\n";}else if(properties.SPHEREMAPPING){shader+=" fragTexcoord = 0.5 + normalMV.xy / 2.0;\n";}else if(properties.TEXTRAFO){shader+=" fragTexcoord = (texTrafoMatrix * vec4(vertTexCoord, 1.0, 1.0)).xy;\n";}else{shader+=" fragTexcoord = vertTexCoord;\n";if(properties.POPGEOMETRY&&x3dom.debug.usePrecisionLevelAsTexCoord===true)
shader+="fragTexcoord = vec2(0.03125 + 0.9375 * (PG_precisionLevel / 16.0), 1.0);";}}
if(properties.LIGHTS){shader+="vec3 ambient   = vec3(0.0, 0.0, 0.0);\n";shader+="vec3 diffuse   = vec3(0.0, 0.0, 0.0);\n";shader+="vec3 specular  = vec3(0.0, 0.0, 0.0);\n";shader+="float _shininess     = shininess;\n";shader+="vec3 _specularColor  = specularColor;\n";shader+="vec3 _emissiveColor  = emissiveColor;\n";shader+="float _ambientIntensity = ambientIntensity;\n";if(!properties.SOLID||properties.TWOSIDEDMAT){shader+="if (dot(normalMV, eye) < 0.0) {\n";shader+="  normalMV *= -1.0;\n";if(properties.SEPARATEBACKMAT){shader+="    rgb = backDiffuseColor;\n";shader+="    alpha = 1.0 - backTransparency;\n";shader+="    _shininess = backShininess;\n";shader+="    _emissiveColor = backEmissiveColor;\n";shader+="    _specularColor = backSpecularColor;\n";shader+="    _ambientIntensity = backAmbientIntensity;\n";}
shader+="  }\n";}
if(properties.LIGHTS){shader+="vec3 ads;\n";for(var l=0;l<properties.LIGHTS;l++){var lightCol="light"+l+"_Color";shader+="ads = lighting(light"+l+"_Type, "+"light"+l+"_Location, "+"light"+l+"_Direction, "+
lightCol+", "+"light"+l+"_Attenuation, "+"light"+l+"_Radius, "+"light"+l+"_Intensity, "+"light"+l+"_AmbientIntensity, "+"light"+l+"_BeamWidth, "+"light"+l+"_CutOffAngle, "+"normalMV, eye, _shininess, _ambientIntensity);\n";shader+="   ambient  += "+lightCol+" * ads.r;\n"+"   diffuse  += "+lightCol+" * ads.g;\n"+"   specular += "+lightCol+" * ads.b;\n";}
shader+="ambient = max(ambient, 0.0);\n";shader+="diffuse = max(diffuse, 0.0);\n";shader+="specular = max(specular, 0.0);\n";}
if(properties.TEXTURED&&!properties.BLENDING){shader+="fragAmbient = ambient;\n";shader+="fragDiffuse = diffuse;\n";shader+="fragColor.rgb = (_emissiveColor + specular*_specularColor);\n";shader+="fragColor.a = alpha;\n";}else{shader+="fragColor.rgb = (_emissiveColor + max(ambient + diffuse, 0.0) * rgb + specular*_specularColor);\n";shader+="fragColor.a = alpha;\n";}}else{if(properties.APPMAT&&!properties.VERTEXCOLOR){shader+="rgb = vec3(0.0, 0.0, 0.0);\n";}
if(properties.TEXTURED&&!properties.BLENDING){shader+="fragAmbient = vec3(0.0);\n";shader+="fragDiffuse = vec3(1.0);\n";shader+="fragColor.rgb = vec3(0.0);\n";shader+="fragColor.a = alpha;\n";}else if(!properties.VERTEXCOLOR&&properties.POINTLINE2D){shader+="fragColor.rgb = emissiveColor;\n";shader+="fragColor.a = alpha;\n";}else{shader+="fragColor.rgb = rgb + emissiveColor;\n";shader+="fragColor.a = alpha;\n";}}
if(properties.FOG){shader+="float f0 = calcFog(-positionMV);\n";shader+="fragColor.rgb = fogColor * (1.0-f0) + f0 * (fragColor.rgb);\n";}
shader+="gl_Position = modelViewProjectionMatrix * vec4(vertPosition, 1.0);\n";shader+="}\n";var vertexShader=gl.createShader(gl.VERTEX_SHADER);gl.shaderSource(vertexShader,shader);gl.compileShader(vertexShader);if(!gl.getShaderParameter(vertexShader,gl.COMPILE_STATUS)){x3dom.debug.logError("[DynamicMobileShader] VertexShader "+gl.getShaderInfoLog(vertexShader));}
return vertexShader;};x3dom.shader.DynamicMobileShader.prototype.generateFragmentShader=function(gl,properties)
{var shader="#ifdef GL_FRAGMENT_PRECISION_HIGH\n";shader+="precision highp float;\n";shader+="#else\n";shader+=" precision mediump float;\n";shader+="#endif\n\n";shader+="varying vec4 fragColor;\n";if(properties.TEXTURED){if(properties.CUBEMAP){shader+="uniform samplerCube cubeMap;\n";shader+="varying vec3 fragViewDir;\n";shader+="varying vec3 fragNormal;\n";shader+="uniform mat4 modelViewMatrixInverse;\n";}else{shader+="uniform sampler2D diffuseMap;           \n";shader+="varying vec2 fragTexcoord;       \n";}
if(!properties.BLENDING){shader+="varying vec3 fragAmbient;\n";shader+="varying vec3 fragDiffuse;\n";}}
shader+="void main(void) {\n";shader+="vec4 color = fragColor;\n";if(properties.TEXTURED){if(properties.CUBEMAP){shader+="vec3 normal = normalize(fragNormal);\n";shader+="vec3 viewDir = normalize(fragViewDir);\n";shader+="vec3 reflected = reflect(viewDir, normal);\n";shader+="reflected = (modelViewMatrixInverse * vec4(reflected,0.0)).xyz;\n";shader+="vec4 texColor = textureCube(cubeMap, reflected);\n";}else{shader+="vec4 texColor = texture2D(diffuseMap, vec2(fragTexcoord.s, 1.0-fragTexcoord.t));\n";}
if(properties.BLENDING){if(properties.CUBEMAP){shader+="color.rgb = mix(color.rgb, texColor.rgb, vec3(0.75));\n";shader+="color.a = texColor.a;\n";}else{shader+="color.rgb *= texColor.rgb;\n";shader+="color.a *= texColor.a;\n";}}else{shader+="color.rgb += max(fragAmbient + fragDiffuse, 0.0) * texColor.rgb;\n";shader+="color.a *= texColor.a;\n";}}
if(properties.TEXT){shader+="if (color.a <= 0.5) discard;\n";}else{shader+="if (color.a <= 0.1) discard;\n";}
shader+="gl_FragColor = clamp(color, 0.0, 1.0);\n";shader+="}\n";var fragmentShader=gl.createShader(gl.FRAGMENT_SHADER);gl.shaderSource(fragmentShader,shader);gl.compileShader(fragmentShader);if(!gl.getShaderParameter(fragmentShader,gl.COMPILE_STATUS)){x3dom.debug.logError("[DynamicMobileShader] FragmentShader "+gl.getShaderInfoLog(fragmentShader));}
return fragmentShader;};x3dom.shader.DynamicShaderPicking=function(gl,properties,pickMode)
{this.program=gl.createProgram();var vertexShader=this.generateVertexShader(gl,properties,pickMode);var fragmentShader=this.generateFragmentShader(gl,properties,pickMode);gl.attachShader(this.program,vertexShader);gl.attachShader(this.program,fragmentShader);gl.bindAttribLocation(this.program,0,"position");gl.linkProgram(this.program);return this.program;};x3dom.shader.DynamicShaderPicking.prototype.generateVertexShader=function(gl,properties,pickMode)
{var shader="";shader+="uniform mat4 modelMatrix;\n";shader+="uniform mat4 modelViewProjectionMatrix;\n";shader+="attribute vec3 position;\n";shader+="uniform vec3 from;\n";shader+="varying vec3 worldCoord;\n";if(pickMode==1){shader+="attribute vec3 color;\n";shader+="varying vec3 fragColor;\n";}else if(pickMode==2){shader+="attribute vec2 texcoord;\n";shader+="varying vec3 fragColor;\n";}
if(properties.REQUIREBBOX){shader+="uniform vec3 bgCenter;\n";shader+="uniform vec3 bgSize;\n";shader+="uniform float bgPrecisionMax;\n";}
if(properties.REQUIREBBOXCOL){shader+="uniform float bgPrecisionColMax;\n";}
if(properties.REQUIREBBOXTEX){shader+="uniform float bgPrecisionTexMax;\n";}
if(properties.VERTEXID){shader+="uniform float shadowIDs;\n";if(pickMode==3){shader+="varying vec3 idCoord;\n";}else{shader+="varying vec2 idCoord;\n";}
shader+="varying float fragID;\n";shader+="attribute float id;\n";}
if(properties.IMAGEGEOMETRY){shader+="uniform vec3 IG_bboxMin;\n";shader+="uniform vec3 IG_bboxMax;\n";shader+="uniform float IG_coordTextureWidth;\n";shader+="uniform float IG_coordTextureHeight;\n";shader+="uniform vec2 IG_implicitMeshSize;\n";for(var i=0;i<properties.IG_PRECISION;i++){shader+="uniform sampler2D IG_coords"+i+"\n;";}
if(properties.IG_INDEXED){shader+="uniform sampler2D IG_index;\n";shader+="uniform float IG_indexTextureWidth;\n";shader+="uniform float IG_indexTextureHeight;\n";}}
if(properties.POPGEOMETRY){shader+="uniform float PG_precisionLevel;\n";shader+="uniform float PG_powPrecision;\n";shader+="uniform vec3 PG_maxBBSize;\n";shader+="uniform vec3 PG_bbMin;\n";shader+="uniform vec3 PG_bbMaxModF;\n";shader+="uniform vec3 PG_bboxShiftVec;\n";shader+="uniform float PG_numAnchorVertices;\n";shader+="attribute float PG_vertexID;\n";}
if(properties.CLIPPLANES){shader+="uniform mat4 modelViewMatrix;\n";shader+="varying vec4 fragPosition;\n";}
shader+="void main(void) {\n";shader+="gl_PointSize = 2.0;\n";shader+="vec3 pos = position;\n";if(properties.VERTEXID){if(pickMode==0){shader+="idCoord = vec2((id + shadowIDs) / 256.0);\n";shader+="idCoord.x = floor(idCoord.x) / 255.0;\n";shader+="idCoord.y = fract(idCoord.y) * 1.00392156862745;\n";shader+="fragID = id;\n";}else if(pickMode==3){shader+="float ID = id + shadowIDs;\n";shader+="float h = floor(ID / 256.0);\n";shader+="idCoord.x = ID - (h * 256.0);\n";shader+="idCoord.z = floor(h / 256.0);\n";shader+="idCoord.y = h - (idCoord.z * 256.0);\n";shader+="idCoord = idCoord.zyx / 255.0;\n";shader+="fragID = id;\n";}else if(pickMode==4){shader+="idCoord = vec2((id + shadowIDs) / 256.0);\n";shader+="idCoord.x = floor(idCoord.x) / 255.0;\n";shader+="idCoord.y = fract(idCoord.y) * 1.00392156862745;\n";shader+="fragID = id;\n";}}
if(properties.IMAGEGEOMETRY){if(properties.IG_INDEXED){shader+="vec2 halfPixel = vec2(0.5/IG_indexTextureWidth,0.5/IG_indexTextureHeight);\n";shader+="vec2 IG_texCoord = vec2(position.x*(IG_implicitMeshSize.x/IG_indexTextureWidth), position.y*(IG_implicitMeshSize.y/IG_indexTextureHeight)) + halfPixel;\n";shader+="vec2 IG_indices = texture2D( IG_index, IG_texCoord ).rg;\n";shader+="halfPixel = vec2(0.5/IG_coordTextureWidth,0.5/IG_coordTextureHeight);\n";shader+="IG_texCoord = (IG_indices * 0.996108948) + halfPixel;\n";}else{shader+="vec2 halfPixel = vec2(0.5/IG_coordTextureWidth, 0.5/IG_coordTextureHeight);\n";shader+="vec2 IG_texCoord = vec2(position.x*(IG_implicitMeshSize.x/IG_coordTextureWidth), position.y*(IG_implicitMeshSize.y/IG_coordTextureHeight)) + halfPixel;\n";}
shader+="pos = texture2D( IG_coordinateTexture, IG_texCoord ).rgb;\n";shader+="pos = pos * (IG_bboxMax - IG_bboxMin) + IG_bboxMin;\n";}else if(properties.POPGEOMETRY){shader+="vec3 offsetVec = step(pos / bgPrecisionMax, PG_bbMaxModF) * PG_bboxShiftVec;\n";shader+="if (PG_precisionLevel <= 2.0) {\n";shader+="pos = floor(pos / PG_powPrecision) * PG_powPrecision;\n";shader+="pos /= (65536.0 - PG_powPrecision);\n";shader+="}\n";shader+="else {\n";shader+="pos /= bgPrecisionMax;\n";shader+="}\n";shader+="pos = (pos + offsetVec + PG_bbMin) * PG_maxBBSize;\n";}else{if(properties.REQUIREBBOX){shader+="pos = bgCenter + bgSize * pos / bgPrecisionMax;\n";}
if(pickMode==1&&!properties.REQUIREBBOXCOL){shader+="fragColor = color;\n";}else if(pickMode==1&&properties.REQUIREBBOXCOL){shader+="fragColor = color / bgPrecisionColMax;\n";}else if(pickMode==2&&!properties.REQUIREBBOXTEX){shader+="fragColor = vec3(abs(texcoord.x), abs(texcoord.y), 0.0);\n";}else if(pickMode==2&&properties.REQUIREBBOXTEX){shader+="vec2 texCoord = texcoord / bgPrecisionTexMax;\n";shader+="fragColor = vec3(abs(texCoord.x), abs(texCoord.y), 0.0);\n";}}
if(properties.CLIPPLANES){shader+="fragPosition = (modelViewMatrix * vec4(pos, 1.0));\n";}
shader+="worldCoord = (modelMatrix * vec4(pos, 1.0)).xyz - from;\n";shader+="gl_Position = modelViewProjectionMatrix * vec4(pos, 1.0);\n";shader+="}\n";var vertexShader=gl.createShader(gl.VERTEX_SHADER);gl.shaderSource(vertexShader,shader);gl.compileShader(vertexShader);if(!gl.getShaderParameter(vertexShader,gl.COMPILE_STATUS)){x3dom.debug.logInfo("VERTEX:\n"+shader);x3dom.debug.logError("VertexShader "+gl.getShaderInfoLog(vertexShader));}
return vertexShader;};x3dom.shader.DynamicShaderPicking.prototype.generateFragmentShader=function(gl,properties,pickMode)
{var shader="#ifdef GL_FRAGMENT_PRECISION_HIGH\n";shader+=" precision highp float;\n";shader+="#else\n";shader+=" precision mediump float;\n";shader+="#endif\n\n";shader+="uniform float highBit;\n";shader+="uniform float lowBit;\n";shader+="uniform float sceneSize;\n";shader+="varying vec3 worldCoord;\n";if(pickMode==1||pickMode==2){shader+="varying vec3 fragColor;\n";}
if(properties.VERTEXID){if(pickMode==3){shader+="varying vec3 idCoord;\n";}else{shader+="varying vec2 idCoord;\n";}
shader+="varying float fragID;\n";}
if(properties.CLIPPLANES){shader+="uniform mat4 viewMatrixInverse;\n";shader+="varying vec4 fragPosition;\n";}
if(properties.MULTIVISMAP){shader+="uniform sampler2D multiVisibilityMap;\n";shader+="uniform float multiVisibilityWidth;\n";shader+="uniform float multiVisibilityHeight;\n";}
if(properties.CLIPPLANES){shader+=x3dom.shader.clipPlanes(properties.CLIPPLANES);}
shader+="void main(void) {\n";if(properties.CLIPPLANES)
{shader+="calculateClipPlanes();\n";}
if(pickMode==1||pickMode==2){shader+="vec4 color = vec4(fragColor, lowBit);\n";}else if(pickMode==4){shader+="vec4 color = vec4(highBit, lowBit, 0.0, 0.0);\n";}else{shader+="vec4 color = vec4(0.0, 0.0, highBit, lowBit);\n";}
if(properties.VERTEXID){if(pickMode==0||pickMode==4){shader+="color.ba = idCoord;\n";}else if(pickMode==3){shader+="color.gba = idCoord;\n";}
if(properties.MULTIVISMAP){shader+="vec2 idTexCoord;\n";shader+="float roundedID = floor(fragID+0.5);\n";shader+="idTexCoord.x = (mod(roundedID, multiVisibilityWidth)) * (1.0 / multiVisibilityWidth) + (0.5 / multiVisibilityWidth);\n";shader+="idTexCoord.y = (floor(roundedID / multiVisibilityHeight)) * (1.0 / multiVisibilityHeight) + (0.5 / multiVisibilityHeight);\n";shader+="vec4 visibility = texture2D( multiVisibilityMap, idTexCoord );\n";shader+="if (visibility.r < 1.0) discard; \n";}}
if(pickMode!=1&&pickMode!=2){shader+="float d = length(worldCoord) / sceneSize;\n";}
if(pickMode==0){shader+="vec2 comp = fract(d * vec2(256.0, 1.0));\n";shader+="color.rg = comp - (comp.rr * vec2(0.0, 1.0/256.0));\n";}else if(pickMode==3){shader+="color.r = d;\n";}
shader+="gl_FragColor = color;\n";shader+="}\n";var fragmentShader=gl.createShader(gl.FRAGMENT_SHADER);gl.shaderSource(fragmentShader,shader);gl.compileShader(fragmentShader);if(!gl.getShaderParameter(fragmentShader,gl.COMPILE_STATUS)){x3dom.debug.logInfo("FRAGMENT:\n"+shader);x3dom.debug.logError("FragmentShader "+gl.getShaderInfoLog(fragmentShader));}
return fragmentShader;};x3dom.shader.DynamicShadowShader=function(gl,properties)
{this.program=gl.createProgram();var vertexShader=this.generateVertexShader(gl,properties);var fragmentShader=this.generateFragmentShader(gl,properties);gl.attachShader(this.program,vertexShader);gl.attachShader(this.program,fragmentShader);gl.bindAttribLocation(this.program,0,"position");gl.linkProgram(this.program);return this.program;};x3dom.shader.DynamicShadowShader.prototype.generateVertexShader=function(gl,properties)
{var shader="";shader+="attribute vec3 position;\n";shader+="uniform mat4 modelViewProjectionMatrix;\n";shader+="varying vec4 projCoords;\n";if(properties.VERTEXID){shader+="varying float fragID;\n";shader+="attribute float id;\n";}
if(properties.REQUIREBBOX){shader+="uniform vec3 bgCenter;\n";shader+="uniform vec3 bgSize;\n";shader+="uniform float bgPrecisionMax;\n";}
if(properties.IMAGEGEOMETRY){shader+="uniform vec3 IG_bboxMin;\n";shader+="uniform vec3 IG_bboxMax;\n";shader+="uniform float IG_coordTextureWidth;\n";shader+="uniform float IG_coordTextureHeight;\n";shader+="uniform vec2 IG_implicitMeshSize;\n";for(var i=0;i<properties.IG_PRECISION;i++){shader+="uniform sampler2D IG_coords"+i+"\n;";}
if(properties.IG_INDEXED){shader+="uniform sampler2D IG_index;\n";shader+="uniform float IG_indexTextureWidth;\n";shader+="uniform float IG_indexTextureHeight;\n";}}
if(properties.POPGEOMETRY){shader+="uniform float PG_precisionLevel;\n";shader+="uniform float PG_powPrecision;\n";shader+="uniform vec3 PG_maxBBSize;\n";shader+="uniform vec3 PG_bbMin;\n";shader+="uniform vec3 PG_bbMaxModF;\n";shader+="uniform vec3 PG_bboxShiftVec;\n";shader+="uniform float PG_numAnchorVertices;\n";shader+="attribute float PG_vertexID;\n";}
if(properties.CLIPPLANES){shader+="uniform mat4 modelViewMatrix;\n";shader+="varying vec4 fragPosition;\n";}
shader+="void main(void) {\n";shader+="    vec3 pos = position;\n";if(properties.IMAGEGEOMETRY){if(properties.IG_INDEXED){shader+="    vec2 halfPixel = vec2(0.5/IG_indexTextureWidth,0.5/IG_indexTextureHeight);\n";shader+="    vec2 IG_texCoord = vec2(position.x*(IG_implicitMeshSize.x/IG_indexTextureWidth), position.y*(IG_implicitMeshSize.y/IG_indexTextureHeight)) + halfPixel;\n";shader+="    vec2 IG_indices = texture2D( IG_index, IG_texCoord ).rg;\n";shader+="    halfPixel = vec2(0.5/IG_coordTextureWidth,0.5/IG_coordTextureHeight);\n";shader+="    IG_texCoord = (IG_indices * 0.996108948) + halfPixel;\n";}else{shader+="    vec2 halfPixel = vec2(0.5/IG_coordTextureWidth, 0.5/IG_coordTextureHeight);\n";shader+="    vec2 IG_texCoord = vec2(position.x*(IG_implicitMeshSize.x/IG_coordTextureWidth), position.y*(IG_implicitMeshSize.y/IG_coordTextureHeight)) + halfPixel;\n";}
shader+="    pos = texture2D( IG_coordinateTexture, IG_texCoord ).rgb;\n";shader+="    pos = pos * (IG_bboxMax - IG_bboxMin) + IG_bboxMin;\n";}else if(properties.POPGEOMETRY){shader+="    vec3 offsetVec = step(pos / bgPrecisionMax, PG_bbMaxModF) * PG_bboxShiftVec;\n";shader+="    if (PG_precisionLevel <= 2.0) {\n";shader+="        pos = floor(pos / PG_powPrecision) * PG_powPrecision;\n";shader+="        pos /= (65536.0 - PG_powPrecision);\n";shader+="    }\n";shader+="    else {\n";shader+="        pos /= bgPrecisionMax;\n";shader+="    }\n";shader+="    pos = (pos + offsetVec + PG_bbMin) * PG_maxBBSize;\n";}else{if(properties.REQUIREBBOX){shader+="    pos = bgCenter + bgSize * pos / bgPrecisionMax;\n";}}
if(properties.VERTEXID){shader+="    fragID = id;\n";}
if(properties.CLIPPLANES){shader+="    fragPosition = (modelViewMatrix * vec4(pos, 1.0));\n";}
shader+="    projCoords = modelViewProjectionMatrix * vec4(pos, 1.0);\n";shader+="    gl_Position = projCoords;\n";shader+="}\n";var vertexShader=gl.createShader(gl.VERTEX_SHADER);gl.shaderSource(vertexShader,shader);gl.compileShader(vertexShader);if(!gl.getShaderParameter(vertexShader,gl.COMPILE_STATUS)){x3dom.debug.logError("[ShadowShader] VertexShader "+gl.getShaderInfoLog(vertexShader));}
return vertexShader;};x3dom.shader.DynamicShadowShader.prototype.generateFragmentShader=function(gl,properties)
{var shader="";shader+="#ifdef GL_FRAGMENT_PRECISION_HIGH\n";shader+="    precision highp float;\n";shader+="#else\n";shader+="    precision mediump float;\n";shader+="#endif\n\n";shader+="varying vec4 projCoords;\n";shader+="uniform float offset;\n";shader+="uniform bool cameraView;\n";if(properties.VERTEXID){shader+="varying float fragID;\n";}
if(properties.MULTIVISMAP){shader+="uniform sampler2D multiVisibilityMap;\n";shader+="uniform float multiVisibilityWidth;\n";shader+="uniform float multiVisibilityHeight;\n";}
if(properties.CLIPPLANES){shader+="uniform mat4 viewMatrixInverse;\n";shader+="varying vec4 fragPosition;\n";shader+=x3dom.shader.clipPlanes(properties.CLIPPLANES);}
if(!x3dom.caps.FP_TEXTURES){shader+=x3dom.shader.rgbaPacking();}
shader+="void main(void) {\n";if(properties.CLIPPLANES)
{shader+="calculateClipPlanes();\n";}
if(properties.MULTIVISMAP){shader+="    vec2 idTexCoord;\n";shader+="    float roundedID = floor(fragID+0.5);\n";shader+="    idTexCoord.x = (mod(roundedID, multiVisibilityWidth)) * (1.0 / multiVisibilityWidth) + (0.5 / multiVisibilityWidth);\n";shader+="    idTexCoord.y = (floor(roundedID / multiVisibilityHeight)) * (1.0 / multiVisibilityHeight) + (0.5 / multiVisibilityHeight);\n";shader+="    vec4 visibility = texture2D( multiVisibilityMap, idTexCoord );\n";shader+="    if (visibility.r < 1.0) discard; \n";}
shader+="    vec3 proj = (projCoords.xyz / projCoords.w);\n";if(!x3dom.caps.FP_TEXTURES){shader+="    gl_FragColor = packDepth(proj.z);\n";}else{shader+="    if (!cameraView){\n";shader+="     proj.z = (proj.z + 1.0)*0.5;\n";shader+="     proj.y = proj.z * proj.z;\n";shader+="    }\n";shader+="    gl_FragColor = vec4(proj, 1.0);\n";}
shader+="}\n";var fragmentShader=gl.createShader(gl.FRAGMENT_SHADER);gl.shaderSource(fragmentShader,shader);gl.compileShader(fragmentShader);if(!gl.getShaderParameter(fragmentShader,gl.COMPILE_STATUS)){x3dom.debug.logError("[ShadowShader] FragmentShader "+gl.getShaderInfoLog(fragmentShader));}
return fragmentShader;};x3dom.shader.ComposedShader=function(gl,shape)
{this.program=gl.createProgram();var vertexShader=this.generateVertexShader(gl,shape);var fragmentShader=this.generateFragmentShader(gl,shape);gl.attachShader(this.program,vertexShader);gl.attachShader(this.program,fragmentShader);gl.bindAttribLocation(this.program,0,"position");gl.linkProgram(this.program);return this.program;};x3dom.shader.ComposedShader.prototype.generateVertexShader=function(gl,shape)
{var shader=shape._cf.appearance.node._shader._vertex._vf.url[0];var vertexShader=gl.createShader(gl.VERTEX_SHADER);gl.shaderSource(vertexShader,shader);gl.compileShader(vertexShader);if(!gl.getShaderParameter(vertexShader,gl.COMPILE_STATUS)){x3dom.debug.logError("[ComposedShader] VertexShader "+gl.getShaderInfoLog(vertexShader));}
return vertexShader;};x3dom.shader.ComposedShader.prototype.generateFragmentShader=function(gl,shape)
{var shader=shape._cf.appearance.node._shader._fragment._vf.url[0];var fragmentShader=gl.createShader(gl.FRAGMENT_SHADER);gl.shaderSource(fragmentShader,shader);gl.compileShader(fragmentShader);if(!gl.getShaderParameter(fragmentShader,gl.COMPILE_STATUS)){x3dom.debug.logError("[ComposedShader] FragmentShader "+gl.getShaderInfoLog(fragmentShader));}
return fragmentShader;};x3dom.shader.NormalShader=function(gl)
{this.program=gl.createProgram();var vertexShader=this.generateVertexShader(gl);var fragmentShader=this.generateFragmentShader(gl);gl.attachShader(this.program,vertexShader);gl.attachShader(this.program,fragmentShader);gl.bindAttribLocation(this.program,0,"position");gl.linkProgram(this.program);return this.program;};x3dom.shader.NormalShader.prototype.generateVertexShader=function(gl)
{var shader="attribute vec3 position;\n"+"attribute vec3 normal;\n"+"uniform vec3 bgCenter;\n"+"uniform vec3 bgSize;\n"+"uniform float bgPrecisionMax;\n"+"uniform float bgPrecisionNorMax;\n"+"uniform mat4 normalMatrix;\n"+"uniform mat4 modelViewProjectionMatrix;\n"+"varying vec3 fragNormal;\n"+"void main(void) {\n"+"    vec3 pos = bgCenter + bgSize * position / bgPrecisionMax;\n"+"    fragNormal = (normalMatrix * vec4(normal / bgPrecisionNorMax, 0.0)).xyz;\n"+"    gl_Position = modelViewProjectionMatrix * vec4(pos, 1.0);\n"+"}\n";var vertexShader=gl.createShader(gl.VERTEX_SHADER);gl.shaderSource(vertexShader,shader);gl.compileShader(vertexShader);if(!gl.getShaderParameter(vertexShader,gl.COMPILE_STATUS)){x3dom.debug.logError("[NormalShader] VertexShader "+gl.getShaderInfoLog(vertexShader));}
return vertexShader;};x3dom.shader.NormalShader.prototype.generateFragmentShader=function(gl)
{var shader="#ifdef GL_FRAGMENT_PRECISION_HIGH\n";shader+="precision highp float;\n";shader+="#else\n";shader+=" precision mediump float;\n";shader+="#endif\n\n";shader+="varying vec3 fragNormal;\n"+"void main(void) {\n"+"    gl_FragColor = vec4(normalize(fragNormal) / 2.0 + 0.5, 1.0);\n"+"}\n";var fragmentShader=gl.createShader(gl.FRAGMENT_SHADER);gl.shaderSource(fragmentShader,shader);gl.compileShader(fragmentShader);if(!gl.getShaderParameter(fragmentShader,gl.COMPILE_STATUS)){x3dom.debug.logError("[NormalShader] FragmentShader "+gl.getShaderInfoLog(fragmentShader));}
return fragmentShader;};x3dom.shader.FrontgroundTextureShader=function(gl)
{this.program=gl.createProgram();var vertexShader=this.generateVertexShader(gl);var fragmentShader=this.generateFragmentShader(gl);gl.attachShader(this.program,vertexShader);gl.attachShader(this.program,fragmentShader);gl.bindAttribLocation(this.program,0,"position");gl.linkProgram(this.program);return this.program;};x3dom.shader.FrontgroundTextureShader.prototype.generateVertexShader=function(gl)
{var shader="attribute vec3 position;\n"+"varying vec2 fragTexCoord;\n"+"\n"+"void main(void) {\n"+"    vec2 texCoord = (position.xy + 1.0) * 0.5;\n"+"    fragTexCoord = texCoord;\n"+"    gl_Position = vec4(position.xy, 0.0, 1.0);\n"+"}\n";var vertexShader=gl.createShader(gl.VERTEX_SHADER);gl.shaderSource(vertexShader,shader);gl.compileShader(vertexShader);if(!gl.getShaderParameter(vertexShader,gl.COMPILE_STATUS)){x3dom.debug.logError("[FrontgroundTextureShader] VertexShader "+gl.getShaderInfoLog(vertexShader));}
return vertexShader;};x3dom.shader.FrontgroundTextureShader.prototype.generateFragmentShader=function(gl)
{var shader="#ifdef GL_FRAGMENT_PRECISION_HIGH\n";shader+="precision highp float;\n";shader+="#else\n";shader+=" precision mediump float;\n";shader+="#endif\n\n";shader+="uniform sampler2D tex;\n"+"varying vec2 fragTexCoord;\n"+"\n"+"void main(void) {\n"+"    vec4 col = texture2D(tex, fragTexCoord);\n"+"    gl_FragColor = vec4(col.rgb, 1.0);\n"+"}\n";var fragmentShader=gl.createShader(gl.FRAGMENT_SHADER);gl.shaderSource(fragmentShader,shader);gl.compileShader(fragmentShader);if(!gl.getShaderParameter(fragmentShader,gl.COMPILE_STATUS)){x3dom.debug.logError("[FrontgroundTextureShader] FragmentShader "+gl.getShaderInfoLog(fragmentShader));}
return fragmentShader;};x3dom.shader.BackgroundTextureShader=function(gl)
{this.program=gl.createProgram();var vertexShader=this.generateVertexShader(gl);var fragmentShader=this.generateFragmentShader(gl);gl.attachShader(this.program,vertexShader);gl.attachShader(this.program,fragmentShader);gl.bindAttribLocation(this.program,0,"position");gl.linkProgram(this.program);return this.program;};x3dom.shader.BackgroundTextureShader.prototype.generateVertexShader=function(gl)
{var shader="attribute vec3 position;\n"+"varying vec2 fragTexCoord;\n"+"uniform vec2 scale;\n"+"uniform vec2 translation;\n"+"\n"+"void main(void) {\n"+"    vec2 texCoord = (position.xy + 1.0) * 0.5;\n"+"    fragTexCoord = texCoord * scale + translation;\n"+"    gl_Position = vec4(position.xy, 0.0, 1.0);\n"+"}\n";var vertexShader=gl.createShader(gl.VERTEX_SHADER);gl.shaderSource(vertexShader,shader);gl.compileShader(vertexShader);if(!gl.getShaderParameter(vertexShader,gl.COMPILE_STATUS)){x3dom.debug.logError("[BackgroundTextureShader] VertexShader "+gl.getShaderInfoLog(vertexShader));}
return vertexShader;};x3dom.shader.BackgroundTextureShader.prototype.generateFragmentShader=function(gl)
{var shader="#ifdef GL_FRAGMENT_PRECISION_HIGH\n";shader+="precision highp float;\n";shader+="#else\n";shader+=" precision mediump float;\n";shader+="#endif\n\n";shader+="uniform sampler2D tex;\n"+"varying vec2 fragTexCoord;\n"+"\n"+"void main(void) {\n"+"    gl_FragColor = texture2D(tex, fragTexCoord);\n"+"}";var fragmentShader=gl.createShader(gl.FRAGMENT_SHADER);gl.shaderSource(fragmentShader,shader);gl.compileShader(fragmentShader);if(!gl.getShaderParameter(fragmentShader,gl.COMPILE_STATUS)){x3dom.debug.logError("[BackgroundTextureShader] FragmentShader "+gl.getShaderInfoLog(fragmentShader));}
return fragmentShader;};x3dom.shader.BackgroundSkyTextureShader=function(gl)
{this.program=gl.createProgram();var vertexShader=this.generateVertexShader(gl);var fragmentShader=this.generateFragmentShader(gl);gl.attachShader(this.program,vertexShader);gl.attachShader(this.program,fragmentShader);gl.bindAttribLocation(this.program,0,"position");gl.linkProgram(this.program);return this.program;};x3dom.shader.BackgroundSkyTextureShader.prototype.generateVertexShader=function(gl)
{var shader="attribute vec3 position;\n"+"attribute vec2 texcoord;\n"+"uniform mat4 modelViewProjectionMatrix;\n"+"varying vec2 fragTexCoord;\n"+"\n"+"void main(void) {\n"+"    fragTexCoord = texcoord;\n"+"    gl_Position = modelViewProjectionMatrix * vec4(position, 1.0);\n"+"}\n";var vertexShader=gl.createShader(gl.VERTEX_SHADER);gl.shaderSource(vertexShader,shader);gl.compileShader(vertexShader);if(!gl.getShaderParameter(vertexShader,gl.COMPILE_STATUS)){x3dom.debug.logError("[BackgroundSkyTextureShader] VertexShader "+gl.getShaderInfoLog(vertexShader));}
return vertexShader;};x3dom.shader.BackgroundSkyTextureShader.prototype.generateFragmentShader=function(gl)
{var shader="#ifdef GL_FRAGMENT_PRECISION_HIGH\n";shader+="precision highp float;\n";shader+="#else\n";shader+=" precision mediump float;\n";shader+="#endif\n\n";shader+="uniform sampler2D tex;\n"+"varying vec2 fragTexCoord;\n"+"\n"+"void main(void) {\n"+"    gl_FragColor = texture2D(tex, fragTexCoord);\n"+"}\n";var fragmentShader=gl.createShader(gl.FRAGMENT_SHADER);gl.shaderSource(fragmentShader,shader);gl.compileShader(fragmentShader);if(!gl.getShaderParameter(fragmentShader,gl.COMPILE_STATUS)){x3dom.debug.logError("[BackgroundSkyTextureShader] FragmentShader "+gl.getShaderInfoLog(fragmentShader));}
return fragmentShader;};x3dom.shader.BackgroundCubeTextureShader=function(gl)
{this.program=gl.createProgram();var vertexShader=this.generateVertexShader(gl);var fragmentShader=this.generateFragmentShader(gl);gl.attachShader(this.program,vertexShader);gl.attachShader(this.program,fragmentShader);gl.bindAttribLocation(this.program,0,"position");gl.linkProgram(this.program);return this.program;};x3dom.shader.BackgroundCubeTextureShader.prototype.generateVertexShader=function(gl)
{var shader="attribute vec3 position;\n"+"uniform mat4 modelViewProjectionMatrix;\n"+"varying vec3 fragNormal;\n"+"\n"+"void main(void) {\n"+"    fragNormal = normalize(position);\n"+"    gl_Position = modelViewProjectionMatrix * vec4(position, 1.0);\n"+"}\n";var vertexShader=gl.createShader(gl.VERTEX_SHADER);gl.shaderSource(vertexShader,shader);gl.compileShader(vertexShader);if(!gl.getShaderParameter(vertexShader,gl.COMPILE_STATUS)){x3dom.debug.logError("[BackgroundCubeTextureShader] VertexShader "+gl.getShaderInfoLog(vertexShader));}
return vertexShader;};x3dom.shader.BackgroundCubeTextureShader.prototype.generateFragmentShader=function(gl)
{var shader="#ifdef GL_FRAGMENT_PRECISION_HIGH\n";shader+="precision highp float;\n";shader+="#else\n";shader+=" precision mediump float;\n";shader+="#endif\n\n";shader+="uniform samplerCube tex;\n"+"varying vec3 fragNormal;\n"+"\n"+"float magn(float val) {\n"+"    return ((val >= 0.0) ? val : -1.0 * val);\n"+"}"+"\n"+"void main(void) {\n"+"    vec3 normal = -reflect(normalize(fragNormal), vec3(0.0,0.0,1.0));\n"+"    if (magn(normal.y) >= magn(normal.x) && magn(normal.y) >= magn(normal.z))\n"+"        normal.xz = -normal.xz;\n"+"    gl_FragColor = textureCube(tex, normal);\n"+"}\n";var fragmentShader=gl.createShader(gl.FRAGMENT_SHADER);gl.shaderSource(fragmentShader,shader);gl.compileShader(fragmentShader);if(!gl.getShaderParameter(fragmentShader,gl.COMPILE_STATUS)){x3dom.debug.logError("[BackgroundCubeTextureShader] FragmentShader "+gl.getShaderInfoLog(fragmentShader));}
return fragmentShader;};x3dom.shader.ShadowRenderingShader=function(gl,shadowedLights)
{this.program=gl.createProgram();var vertexShader=this.generateVertexShader(gl);var fragmentShader=this.generateFragmentShader(gl,shadowedLights);gl.attachShader(this.program,vertexShader);gl.attachShader(this.program,fragmentShader);gl.bindAttribLocation(this.program,0,"position");gl.linkProgram(this.program);return this.program;};x3dom.shader.ShadowRenderingShader.prototype.generateVertexShader=function(gl)
{var shader="";shader+="attribute vec2 position;\n";shader+="varying vec2 vPosition;\n";shader+="void main(void) {\n";shader+=" vPosition = position;\n";shader+=" gl_Position = vec4(position, -1.0, 1.0);\n";shader+="}\n";var vertexShader=gl.createShader(gl.VERTEX_SHADER);gl.shaderSource(vertexShader,shader);gl.compileShader(vertexShader);if(!gl.getShaderParameter(vertexShader,gl.COMPILE_STATUS)){x3dom.debug.logError("[ShadowRendering] VertexShader "+gl.getShaderInfoLog(vertexShader));}
return vertexShader;};x3dom.shader.ShadowRenderingShader.prototype.generateFragmentShader=function(gl,shadowedLights)
{var shader="#ifdef GL_FRAGMENT_PRECISION_HIGH\n";shader+="precision highp float;\n";shader+="#else\n";shader+=" precision mediump float;\n";shader+="#endif\n\n";shader+="uniform mat4 inverseViewProj;\n";shader+="uniform mat4 inverseProj;\n";shader+="varying vec2 vPosition;\n";shader+="uniform sampler2D sceneMap;\n";for(var i=0;i<5;i++)
shader+="uniform float cascade"+i+"_Depth;\n";for(var l=0;l<shadowedLights.length;l++){shader+="uniform float light"+l+"_On;\n"+"uniform float light"+l+"_Type;\n"+"uniform vec3  light"+l+"_Location;\n"+"uniform vec3  light"+l+"_Direction;\n"+"uniform vec3  light"+l+"_Attenuation;\n"+"uniform float light"+l+"_Radius;\n"+"uniform float light"+l+"_BeamWidth;\n"+"uniform float light"+l+"_CutOffAngle;\n"+"uniform float light"+l+"_ShadowIntensity;\n"+"uniform float light"+l+"_ShadowOffset;\n"+"uniform mat4 light"+l+"_ViewMatrix;\n";for(var j=0;j<6;j++){shader+="uniform mat4 light"+l+"_"+j+"_Matrix;\n";shader+="uniform sampler2D light"+l+"_"+j+"_ShadowMap;\n";}
for(var j=0;j<5;j++)
shader+="uniform float light"+l+"_"+j+"_Split;\n";}
if(!x3dom.caps.FP_TEXTURES||x3dom.caps.MOBILE)
shader+=x3dom.shader.rgbaPacking();shader+=x3dom.shader.shadowRendering();shader+=x3dom.shader.gammaCorrectionDecl({});shader+="void main(void) {\n"+" float shadowValue = 1.0;\n"+" vec2 texCoordsSceneMap = (vPosition + 1.0)*0.5;\n"+" vec4 projCoords = texture2D(sceneMap, texCoordsSceneMap);\n"+" if (projCoords != vec4(1.0,1.0,1.0,0.0)){\n";if(!x3dom.caps.FP_TEXTURES||x3dom.caps.MOBILE){shader+=" projCoords.z = unpackDepth(projCoords);\n"+" projCoords.w = 1.0;\n";}
shader+=" projCoords = projCoords / projCoords.w;\n"+" projCoords.xy = vPosition;\n"+" vec4 eyeCoords = inverseProj*projCoords;\n"+" vec4 worldCoords = inverseViewProj*projCoords;\n"+" float lightInfluence = 0.0;\n";for(var l=0;l<shadowedLights.length;l++){shader+=" lightInfluence = getLightInfluence(light"+l+"_Type, light"+l+"_ShadowIntensity, light"+l+"_On, light"+l+"_Location, light"+l+"_Direction, "+"light"+l+"_CutOffAngle, light"+l+"_BeamWidth, light"+l+"_Attenuation, light"+l+"_Radius, eyeCoords.xyz/eyeCoords.w);\n"+" if (lightInfluence != 0.0){\n"+"  vec4 shadowMapValues;\n"+"  float viewSampleDepth;\n";if(!x3dom.isa(shadowedLights[l],x3dom.nodeTypes.PointLight)){shader+="  getShadowValuesCascaded(shadowMapValues, viewSampleDepth, worldCoords, -eyeCoords.z/eyeCoords.w,"+"light"+l+"_0_Matrix,light"+l+"_1_Matrix,light"+l+"_2_Matrix,light"+l+"_3_Matrix,light"+l+"_4_Matrix,light"+l+"_5_Matrix,"+"light"+l+"_0_ShadowMap,light"+l+"_1_ShadowMap,light"+l+"_2_ShadowMap,light"+l+"_3_ShadowMap,"+"light"+l+"_4_ShadowMap,light"+l+"_5_ShadowMap, light"+l+"_0_Split, light"+l+"_1_Split, light"+l+"_2_Split, light"+l+"_3_Split, \n"+"light"+l+"_4_Split);\n";}else{shader+="  getShadowValuesPointLight(shadowMapValues, viewSampleDepth, light"+l+"_Location, worldCoords, light"+l+"_ViewMatrix, "+"light"+l+"_0_Matrix,light"+l+"_1_Matrix,light"+l+"_2_Matrix,light"+l+"_3_Matrix,light"+l+"_4_Matrix,light"+l+"_5_Matrix,"+"light"+l+"_0_ShadowMap,light"+l+"_1_ShadowMap,light"+l+"_2_ShadowMap,light"+l+"_3_ShadowMap,"+"light"+l+"_4_ShadowMap,light"+l+"_5_ShadowMap);\n";}
if(!x3dom.caps.FP_TEXTURES||x3dom.caps.MOBILE)
shader+=" shadowValue *= clamp(ESM(shadowMapValues.z, viewSampleDepth, light"+l+"_ShadowOffset), "+"    1.0 - light"+l+"_ShadowIntensity*lightInfluence, 1.0);\n";else
shader+="  shadowValue *= clamp(VSM(shadowMapValues.zy, viewSampleDepth, light"+l+"_ShadowOffset), "+"    1.0 - light"+l+"_ShadowIntensity*lightInfluence, 1.0);\n";shader+=" }\n";}
shader+="}\n"+" gl_FragColor = "+x3dom.shader.encodeGamma({},"vec4(shadowValue, shadowValue, shadowValue, 1.0)")+";\n"+"}\n";var fragmentShader=gl.createShader(gl.FRAGMENT_SHADER);gl.shaderSource(fragmentShader,shader);gl.compileShader(fragmentShader);if(!gl.getShaderParameter(fragmentShader,gl.COMPILE_STATUS)){x3dom.debug.logError("[ShadowRendering] FragmentShader "+gl.getShaderInfoLog(fragmentShader));}
return fragmentShader;};x3dom.shader.TextureRefinementShader=function(gl){this.program=gl.createProgram();var vertexShader=this.generateVertexShader(gl);var fragmentShader=this.generateFragmentShader(gl);gl.attachShader(this.program,vertexShader);gl.attachShader(this.program,fragmentShader);gl.bindAttribLocation(this.program,0,"position");gl.linkProgram(this.program);return this.program;};x3dom.shader.TextureRefinementShader.prototype.generateVertexShader=function(gl){var shader="attribute vec2 position;\n"+"varying vec2 fragTexCoord;\n"+"\n"+"void main(void) {\n"+"    fragTexCoord = (position.xy + 1.0) / 2.0;\n"+"    gl_Position = vec4(position, -1.0, 1.0);\n"+"}\n";var vertexShader=gl.createShader(gl.VERTEX_SHADER);gl.shaderSource(vertexShader,shader);gl.compileShader(vertexShader);if(!gl.getShaderParameter(vertexShader,gl.COMPILE_STATUS)){x3dom.debug.logError("[TextureRefinementShader] VertexShader "+gl.getShaderInfoLog(vertexShader));}
return vertexShader;};x3dom.shader.TextureRefinementShader.prototype.generateFragmentShader=function(gl){var shader="#ifdef GL_FRAGMENT_PRECISION_HIGH\n"+" precision highp float;\n"+"#else\n"+" precision mediump float;\n"+"#endif\n\n";shader+="uniform sampler2D stamp;\n"+"uniform sampler2D lastTex;\n"+"uniform sampler2D curTex;\n"+"uniform int mode;\n"+"uniform vec2 repeat;\n"+"varying vec2 fragTexCoord;\n"+"\n"+"void init(void);\n"+"void refine(void);\n"+"\n"+"void main(void) {\n"+"    if (mode == 0) { init(); }\n"+"    else { refine(); }\n"+"}\n"+"\n"+"void init(void) {\n"+"    gl_FragColor = texture2D(curTex, fragTexCoord);\n"+"}\n"+"\n"+"void refine(void) {\n"+"    vec3 red = texture2D(stamp, repeat * fragTexCoord).rgb;\n"+"    vec3 v1  = texture2D(lastTex, fragTexCoord).rgb;\n"+"    vec3 v2  = texture2D(curTex, fragTexCoord).rgb;\n"+"    if (red.r <= 0.5) {\n"+"        gl_FragColor = vec4(v1, 1.0);\n"+"    }\n"+"    else {\n"+"        gl_FragColor = vec4(v2, 1.0);\n"+"    }\n"+"}\n";var fragmentShader=gl.createShader(gl.FRAGMENT_SHADER);gl.shaderSource(fragmentShader,shader);gl.compileShader(fragmentShader);if(!gl.getShaderParameter(fragmentShader,gl.COMPILE_STATUS)){x3dom.debug.logError("[TextureRefinementShader] FragmentShader "+gl.getShaderInfoLog(fragmentShader));}
return fragmentShader;};x3dom.shader.BlurShader=function(gl)
{this.program=gl.createProgram();var vertexShader=this.generateVertexShader(gl);var fragmentShader=this.generateFragmentShader(gl);gl.attachShader(this.program,vertexShader);gl.attachShader(this.program,fragmentShader);gl.bindAttribLocation(this.program,0,"position");gl.linkProgram(this.program);return this.program;};x3dom.shader.BlurShader.prototype.generateVertexShader=function(gl)
{var shader="";shader+="attribute vec2 position;\n";shader+="varying vec2 vPosition;\n";shader+="void main(void) {\n";shader+=" vPosition = position;\n";shader+=" gl_Position = vec4(position, -1.0, 1.0);\n";shader+="}\n";var vertexShader=gl.createShader(gl.VERTEX_SHADER);gl.shaderSource(vertexShader,shader);gl.compileShader(vertexShader);if(!gl.getShaderParameter(vertexShader,gl.COMPILE_STATUS)){x3dom.debug.logError("[BlurShader] VertexShader "+gl.getShaderInfoLog(vertexShader));}
return vertexShader;};x3dom.shader.BlurShader.prototype.generateFragmentShader=function(gl)
{var shader="#ifdef GL_FRAGMENT_PRECISION_HIGH\n";shader+="precision highp float;\n";shader+="#else\n";shader+=" precision mediump float;\n";shader+="#endif\n\n";shader+="varying vec2 vPosition;\n"+"uniform sampler2D texture;\n"+"uniform bool horizontal;\n"+"uniform float pixelSizeHor;\n"+"uniform float pixelSizeVert;\n"+"uniform int filterSize;\n";if(!x3dom.caps.FP_TEXTURES||x3dom.caps.MOBILE){shader+=x3dom.shader.rgbaPacking()+"void main(void) {\n"+" vec2 texCoords = (vPosition + 1.0)*0.5;\n"+" vec2 offset;\n"+" if (horizontal) offset = vec2(pixelSizeHor, 0.0);\n"+" else offset = vec2(0.0, pixelSizeVert);\n"+" float depth = unpackDepth(texture2D(texture, texCoords));\n"+" if (filterSize == 3){\n"+"  depth = depth * 0.3844;\n"+"  depth += 0.3078*unpackDepth(texture2D(texture, texCoords-offset));\n"+"  depth += 0.3078*unpackDepth(texture2D(texture, texCoords+offset));\n"+" } else if (filterSize == 5){\n"+"  depth = depth * 0.2921;\n"+"  depth += 0.2339*unpackDepth(texture2D(texture, texCoords-offset));\n"+"  depth += 0.2339*unpackDepth(texture2D(texture, texCoords+offset));\n"+"  depth += 0.1201*unpackDepth(texture2D(texture, texCoords-2.0*offset));\n"+"  depth += 0.1201*unpackDepth(texture2D(texture, texCoords+2.0*offset));\n"+" } else if (filterSize == 7){\n"+"  depth = depth * 0.2161;\n"+"  depth += 0.1907*unpackDepth(texture2D(texture, texCoords-offset));\n"+"  depth += 0.1907*unpackDepth(texture2D(texture, texCoords+offset));\n"+"  depth += 0.1311*unpackDepth(texture2D(texture, texCoords-2.0*offset));\n"+"  depth += 0.1311*unpackDepth(texture2D(texture, texCoords+2.0*offset));\n"+"  depth += 0.0702*unpackDepth(texture2D(texture, texCoords-3.0*offset));\n"+"  depth += 0.0702*unpackDepth(texture2D(texture, texCoords+3.0*offset));\n"+" }\n"+" gl_FragColor = packDepth(depth);\n"+"}\n";}else{shader+="void main(void) {\n"+" vec2 texCoords = (vPosition + 1.0)*0.5;\n"+" vec2 offset;\n"+" if (horizontal) offset = vec2(pixelSizeHor, 0.0);\n"+" else offset = vec2(0.0, pixelSizeVert);\n"+" vec4 color = texture2D(texture, texCoords);\n"+" if (filterSize == 3){\n"+"  color = color * 0.3844;\n"+"  color += 0.3078*texture2D(texture, texCoords-offset);\n"+"  color += 0.3078*texture2D(texture, texCoords+offset);\n"+" } else if (filterSize == 5){\n"+"  color = color * 0.2921;\n"+"  color += 0.2339*texture2D(texture, texCoords-offset);\n"+"  color += 0.2339*texture2D(texture, texCoords+offset);\n"+"  color += 0.1201*texture2D(texture, texCoords-2.0*offset);\n"+"  color += 0.1201*texture2D(texture, texCoords+2.0*offset);\n"+" } else if (filterSize == 7){\n"+"  color = color * 0.2161;\n"+"  color += 0.1907*texture2D(texture, texCoords-offset);\n"+"  color += 0.1907*texture2D(texture, texCoords+offset);\n"+"  color += 0.1311*texture2D(texture, texCoords-2.0*offset);\n"+"  color += 0.1311*texture2D(texture, texCoords+2.0*offset);\n"+"  color += 0.0702*texture2D(texture, texCoords-3.0*offset);\n"+"  color += 0.0702*texture2D(texture, texCoords+3.0*offset);\n"+" }\n"+" gl_FragColor = color;\n"+"}\n";}
var fragmentShader=gl.createShader(gl.FRAGMENT_SHADER);gl.shaderSource(fragmentShader,shader);gl.compileShader(fragmentShader);if(!gl.getShaderParameter(fragmentShader,gl.COMPILE_STATUS)){x3dom.debug.logError("[BlurShader] FragmentShader "+gl.getShaderInfoLog(fragmentShader));}
return fragmentShader;};x3dom.shader.KHRMaterialCommonsShader=function(gl,properties)
{this.program=gl.createProgram();var vertexShader=this.generateVertexShader(gl);var fragmentShader=this.generateFragmentShader(gl,properties);gl.attachShader(this.program,vertexShader);gl.attachShader(this.program,fragmentShader);gl.bindAttribLocation(this.program,0,"position");gl.linkProgram(this.program);return this.program;};x3dom.shader.KHRMaterialCommonsShader.prototype.generateVertexShader=function(gl)
{var shader="precision highp float;\n"+"attribute vec3 position;"+"attribute vec3 normal;"+"attribute vec3 texcoord;"+"varying vec3 v_eye;"+"varying vec3 v_normal;"+"varying vec3 v_texcoord;"+"uniform mat4 modelViewProjectionMatrix;"+"uniform mat4 modelViewMatrix;"+"uniform mat4 normalMatrix;"+"void main (void)"+"{"+"    vec4 pos = modelViewProjectionMatrix * vec4(position, 1.0);"+"    v_eye = (modelViewMatrix * vec4(position, 1.0)).xyz;"+"    v_normal = (normalMatrix * vec4(normal,1.0)).xyz;"+"    v_texcoord = texcoord;"+"    gl_Position = pos;"+"}";var vertexShader=gl.createShader(gl.VERTEX_SHADER);gl.shaderSource(vertexShader,shader);gl.compileShader(vertexShader);if(!gl.getShaderParameter(vertexShader,gl.COMPILE_STATUS)){x3dom.debug.logError("[KHRMaterialCommonsShader] VertexShader "+gl.getShaderInfoLog(vertexShader));}
return vertexShader;};x3dom.shader.KHRMaterialCommonsShader.prototype.generateFragmentShader=function(gl,properties)
{var shader="precision highp float;\n"+"varying vec3 v_eye;\n"+"varying vec3 v_normal;\n"+"varying vec3 v_texcoord;\n"+"uniform vec4 lightVector;\n"+"uniform vec4 ambient;\n";if(properties.LIGHTS||properties.CLIPPLANES)
{shader+="varying vec4 fragPosition;\n";shader+="uniform float isOrthoView;\n";}
if(properties.LIGHTS){if(properties.NORMALMAP&&properties.NORMALSPACE=="OBJECT"){}else{shader+="varying vec3 fragNormal;\n";}
shader+=x3dom.shader.light(properties.LIGHTS);}
if(properties.USE_DIFFUSE_TEX==0)
shader+="uniform vec4 diffuse;\n";else
shader+="uniform sampler2D diffuseTex;\n";if(properties.USE_EMISSION_TEX==0)
shader+="uniform vec4 emission;\n";else
shader+="uniform sampler2D emissionTex;\n";if(properties.USE_SPECULAR_TEX==0)
shader+="uniform vec4 specular;\n";else
shader+="uniform sampler2D specularTex;\n";shader+="uniform float shininess;\n"+"uniform float transparency;\n"+"uniform float ambientIntensity;\n"+"uniform vec4 ambientLight;\n"+"uniform int technique;\n"+"void main(void)\n"+"{\n"+"vec4 I = -vec4(normalize(v_eye),1.0);\n"+"vec4 N = vec4(normalize(v_normal),1.0);\n"+"vec4 al = ambientLight;\n"+"vec4 L = normalize(lightVector-vec4(v_eye,1.0));\n";if(properties.USE_DIFFUSE_TEX==0)
shader+="vec4 _diffuse = diffuse;\n";else
shader+="vec4 _diffuse = texture2D(diffuseTex, v_texcoord.xy);\n";if(properties.USE_SPECULAR_TEX==0)
shader+="vec4 _specularColor = specular;\n";else
shader+="vec4 _specularColor = texture2D(specularTex, v_texcoord.xy);\n";if(properties.USE_EMISSION_TEX==0)
shader+="vec4 _emission = emission;\n";else
shader+="vec4 _emission = texture2D(emissionTex, v_texcoord.xy);\n";shader+="vec4 color;\n"+"if(technique == 0) // BLINN\n"+"{\n"+"vec4 H = normalize(I+L);\n"+"color = _emission + ambient * al + _diffuse * max(dot(N,L),0.0) + _specularColor * pow(max(dot(H,N),0.0),shininess);\n"+"}\n"+"else if(technique==1) // PHONG\n"+"{\n"+"vec4 R = -reflect(L,N);\n"+"color = _emission + ambient * al + _diffuse * max(dot(N,L),0.0) + _specularColor * pow(max(dot(R,I),0.0),shininess);\n"+"}\n"+"else if(technique==2) // LAMBERT\n"+"{\n"+"color = _emission + ambient * al + _diffuse * max(dot(N,L), 0.0);\n"+"}\n"+"else if(technique==3) // CONSTANT\n"+"{\n"+"color = _emission + ambient * al;\n"+"}\n";if(properties.LIGHTS){shader+="vec3 ambient   = vec3(0.0, 0.0, 0.0);\n";shader+="vec3 diffuse   = vec3(0.0, 0.0, 0.0);\n";shader+="vec3 specular  = vec3(0.0, 0.0, 0.0);\n";shader+="vec3 eye;\n";shader+="if ( isOrthoView > 0.0 ) {\n";shader+="    eye = vec3(0.0, 0.0, 1.0);\n";shader+="} else {\n";shader+="    eye = -v_eye.xyz;\n";shader+="}\n";shader+="vec3 ads;\n";for(var l=0;l<properties.LIGHTS;l++){var lightCol="light"+l+"_Color";shader+="ads = lighting(light"+l+"_Type, "+"light"+l+"_Location, "+"light"+l+"_Direction, "+
lightCol+", "+"light"+l+"_Attenuation, "+"light"+l+"_Radius, "+"light"+l+"_Intensity, "+"light"+l+"_AmbientIntensity, "+"light"+l+"_BeamWidth, "+"light"+l+"_CutOffAngle, "+"v_normal, eye, shininess, ambientIntensity);\n";shader+="ambient  += "+lightCol+" * ads.r;\n"+"diffuse  += "+lightCol+" * ads.g;\n"+"specular += "+lightCol+" * ads.b;\n";}
shader+="ambient = max(ambient, 0.0);\n";shader+="diffuse = max(diffuse, 0.0);\n";shader+="specular = max(specular, 0.0);\n";shader+="color.rgb = (_emission.rgb + max(ambient + diffuse, 0.0) * color.rgb + specular*_specularColor.rgb);\n";}
shader+="gl_FragColor = vec4(color.rgb, 1.0-transparency);\n"+"}";var fragmentShader=gl.createShader(gl.FRAGMENT_SHADER);gl.shaderSource(fragmentShader,shader);gl.compileShader(fragmentShader);if(!gl.getShaderParameter(fragmentShader,gl.COMPILE_STATUS)){x3dom.debug.logError("[KHRMaterialCommonsShader] FragmentShader "+gl.getShaderInfoLog(fragmentShader));}
return fragmentShader;};x3dom.shader.SSAOShader=function(gl)
{this.program=gl.createProgram();var vertexShader=this.generateVertexShader(gl);var fragmentShader=this.generateFragmentShader(gl);gl.attachShader(this.program,vertexShader);gl.attachShader(this.program,fragmentShader);gl.bindAttribLocation(this.program,0,"position");gl.linkProgram(this.program);return this.program;};x3dom.shader.SSAOShader.prototype.generateVertexShader=function(gl)
{var shader="attribute vec3 position;\n"+"varying vec2 depthTexCoord;\n"+"varying vec2 randomTexCoord;\n"+"uniform vec2 randomTextureTilingFactor;\n"+"\n"+"void main(void) {\n"+"    vec2 texCoord = (position.xy + 1.0) * 0.5;\n"+"    depthTexCoord = texCoord;\n"+"  randomTexCoord = randomTextureTilingFactor*texCoord;\n"+"    gl_Position = vec4(position.xy, 0.0, 1.0);\n"+"}\n";var vertexShader=gl.createShader(gl.VERTEX_SHADER);gl.shaderSource(vertexShader,shader);gl.compileShader(vertexShader);if(!gl.getShaderParameter(vertexShader,gl.COMPILE_STATUS)){x3dom.debug.logError("[SSAOShader] VertexShader "+gl.getShaderInfoLog(vertexShader));}
return vertexShader;};x3dom.shader.SSAOShader.depthReconsructionFunctionCode=function()
{var code="uniform float depthReconstructionConstantA;\n"+"uniform float depthReconstructionConstantB;\n";if(!x3dom.caps.FP_TEXTURES||x3dom.caps.MOBILE)
code+=x3dom.shader.rgbaPacking();code+="float getDepth(vec2 depthTexCoord) {\n"+"    vec4 col = texture2D(depthTexture, depthTexCoord);\n"+"    float d;\n";if(!x3dom.caps.FP_TEXTURES||x3dom.caps.MOBILE){code+="    d = unpackDepth(col);\n";}else{code+="    d = col.b;\n"}
code+="    return depthReconstructionConstantB/(depthReconstructionConstantA+d);\n";code+="}\n";return code;}
x3dom.shader.SSAOShader.prototype.generateFragmentShader=function(gl)
{var shader="#ifdef GL_FRAGMENT_PRECISION_HIGH\n";shader+="precision highp float;\n";shader+="#else\n";shader+=" precision mediump float;\n";shader+="#endif\n\n";shader+="uniform sampler2D depthTexture;\n"+"uniform sampler2D randomTexture;\n"+"uniform float nearPlane;\n"+"uniform float farPlane;\n"+"uniform float radius;\n"+"uniform float depthBufferEpsilon;\n"+"uniform vec3 samples[16];\n"+"varying vec2 depthTexCoord;\n"+"varying vec2 randomTexCoord;\n";shader+=x3dom.shader.SSAOShader.depthReconsructionFunctionCode();shader+="void main(void) {\n"+"    float referenceDepth = getDepth(depthTexCoord);\n"+"    if(referenceDepth == 1.0)\n"+"    {\n"+"        gl_FragColor = vec4(1.0,1.0,1.0, 1.0);\n"+"        return;\n"+"    }\n"+"    int numOcclusions = 0;\n"+"    for(int i = 0; i<16; ++i){\n"+"        float scale  = 1.0/referenceDepth;\n"+"        vec3 samplepos = reflect(samples[i],texture2D(randomTexture,randomTexCoord).xyz*2.0-vec3(1.0,1.0,1.0));\n"+"        float sampleDepth = getDepth(depthTexCoord+samplepos.xy*scale*radius);\n"+"        //if(abs(sampleDepth-referenceDepth)<=radius*(1.0/nearPlane))\n"+"        if( sampleDepth < referenceDepth-depthBufferEpsilon) {\n"+"            ++numOcclusions;\n"+"        }\n"+"    }\n"+"    float r = 1.0-float(numOcclusions)/16.0;\n"+"    r*=2.0;\n"+"    gl_FragColor = vec4(r,r,r, 1.0);\n"+"}\n";var fragmentShader=gl.createShader(gl.FRAGMENT_SHADER);gl.shaderSource(fragmentShader,shader);gl.compileShader(fragmentShader);if(!gl.getShaderParameter(fragmentShader,gl.COMPILE_STATUS)){x3dom.debug.logError("[SSAOhader] FragmentShader "+gl.getShaderInfoLog(fragmentShader));}
return fragmentShader;};x3dom.shader.SSAOBlurShader=function(gl)
{this.program=gl.createProgram();var vertexShader=this.generateVertexShader(gl);var fragmentShader=this.generateFragmentShader(gl);gl.attachShader(this.program,vertexShader);gl.attachShader(this.program,fragmentShader);gl.bindAttribLocation(this.program,0,"position");gl.linkProgram(this.program);return this.program;};x3dom.shader.SSAOBlurShader.prototype.generateVertexShader=function(gl)
{var shader="attribute vec3 position;\n"+"varying vec2 fragTexCoord;\n"+"\n"+"void main(void) {\n"+"    vec2 texCoord = (position.xy + 1.0) * 0.5;\n"+"    fragTexCoord = texCoord;\n"+"    gl_Position = vec4(position.xy, 0.0, 1.0);\n"+"}\n";var vertexShader=gl.createShader(gl.VERTEX_SHADER);gl.shaderSource(vertexShader,shader);gl.compileShader(vertexShader);if(!gl.getShaderParameter(vertexShader,gl.COMPILE_STATUS)){x3dom.debug.logError("[SSAOShader] VertexShader "+gl.getShaderInfoLog(vertexShader));}
return vertexShader;};x3dom.shader.SSAOBlurShader.prototype.generateFragmentShader=function(gl)
{var shader="#ifdef GL_FRAGMENT_PRECISION_HIGH\n";shader+="precision highp float;\n";shader+="#else\n";shader+=" precision mediump float;\n";shader+="#endif\n\n";shader+="uniform sampler2D SSAOTexture;\n"+"uniform sampler2D depthTexture;\n"+"uniform float nearPlane;\n"+"uniform float farPlane;\n"+"uniform float amount;\n"+"uniform vec2 pixelSize;\n"+"uniform float depthThreshold;\n"+"varying vec2 fragTexCoord;\n";shader+=x3dom.shader.SSAOShader.depthReconsructionFunctionCode();shader+="void main(void) {\n"+"    float sum = 0.0;\n"+"    float numSamples = 0.0;\n"+"    float referenceDepth = getDepth(fragTexCoord);\n"+"    for(int i = -2; i<2;i++){\n"+"        for(int j = -2; j<2;j++){\n"+"            vec2 sampleTexCoord = fragTexCoord+vec2(pixelSize.x*float(i),pixelSize.y*float(j));\n"+"            if(abs(referenceDepth - getDepth(sampleTexCoord))<depthThreshold){\n"+"                sum+= texture2D(SSAOTexture,sampleTexCoord).r;\n"+"                numSamples++;\n"+"    }}}\n"+"    float intensity = mix(1.0,sum/numSamples,amount);\n"+"    gl_FragColor = vec4(intensity,intensity,intensity,1.0);\n"+"}\n";var fragmentShader=gl.createShader(gl.FRAGMENT_SHADER);gl.shaderSource(fragmentShader,shader);gl.compileShader(fragmentShader);if(!gl.getShaderParameter(fragmentShader,gl.COMPILE_STATUS)){x3dom.debug.logError("[SSAOhader] FragmentShader "+gl.getShaderInfoLog(fragmentShader));}
return fragmentShader;};x3dom.SSAO={};x3dom.SSAO.isEnabled=function(scene){return scene.getEnvironment()._vf.SSAO};x3dom.SSAO.reinitializeShadersIfNecessary=function(gl){if(x3dom.SSAO.shaderProgram===undefined){x3dom.SSAO.shaderProgram=x3dom.Utils.wrapProgram(gl,new x3dom.shader.SSAOShader(gl),"ssao");}
if(x3dom.SSAO.blurShaderProgram===undefined){x3dom.SSAO.blurShaderProgram=x3dom.Utils.wrapProgram(gl,new x3dom.shader.SSAOBlurShader(gl),"ssao-blur");}};x3dom.SSAO.reinitializeRandomTextureIfNecessary=function(gl,scene){var sizeHasChanged=scene.getEnvironment()._vf.SSAOrandomTextureSize!=x3dom.SSAO.currentRandomTextureSize;if(x3dom.SSAO.randomTexture===undefined){x3dom.SSAO.randomTexture=gl.createTexture();}
if(x3dom.SSAO.randomTexture===undefined||sizeHasChanged){gl.bindTexture(gl.TEXTURE_2D,x3dom.SSAO.randomTexture);var rTexSize=x3dom.SSAO.currentRandomTextureSize=scene.getEnvironment()._vf.SSAOrandomTextureSize;var randomImageBuffer=new ArrayBuffer(rTexSize*rTexSize*4);var randomImageView=new Uint8Array(randomImageBuffer);for(var i=0;i<rTexSize*rTexSize;++i){var x=Math.random()*2.0-1.0;var y=Math.random()*2.0-1.0;var z=0;var length=Math.sqrt(x*x+y*y+z*z);x/=length;y/=length;randomImageView[4*i]=(x+1.0)*0.5*255.0;randomImageView[4*i+1]=(y+1.0)*0.5*255.0;randomImageView[4*i+2]=(z+1.0)*0.5*255.0;randomImageView[4*i+3]=255;}
gl.texImage2D(gl.TEXTURE_2D,0,gl.RGBA,rTexSize,rTexSize,0,gl.RGBA,gl.UNSIGNED_BYTE,randomImageView);gl.bindTexture(gl.TEXTURE_2D,null);}};x3dom.SSAO.reinitializeFBOIfNecessary=function(gl,canvas){var dimensionsHaveChanged=x3dom.SSAO.currentFBOWidth!=canvas.width||x3dom.SSAO.currentFBOHeight!=canvas.height;if(x3dom.SSAO.fbo===undefined||dimensionsHaveChanged)
{x3dom.SSAO.currentFBOWidth=canvas.width;x3dom.SSAO.currentFBOHeight=canvas.height;var oldfbo=gl.getParameter(gl.FRAMEBUFFER_BINDING);if(x3dom.SSAO.fbo===undefined){x3dom.SSAO.fbo=gl.createFramebuffer();}
gl.bindFramebuffer(gl.FRAMEBUFFER,x3dom.SSAO.fbo);if(x3dom.SSAO.fbotex===undefined){x3dom.SSAO.fbotex=gl.createTexture();}
gl.bindTexture(gl.TEXTURE_2D,x3dom.SSAO.fbotex);gl.texImage2D(gl.TEXTURE_2D,0,gl.RGBA,x3dom.SSAO.currentFBOWidth,x3dom.SSAO.currentFBOHeight,0,gl.RGBA,gl.UNSIGNED_BYTE,null);gl.bindTexture(gl.TEXTURE_2D,null);gl.framebufferTexture2D(gl.FRAMEBUFFER,gl.COLOR_ATTACHMENT0,gl.TEXTURE_2D,x3dom.SSAO.fbotex,0);gl.bindFramebuffer(gl.FRAMEBUFFER,oldfbo);}};x3dom.SSAO.render=function(stateManager,gl,scene,tex,canvas,fbo){var oldfbo=gl.getParameter(gl.FRAMEBUFFER_BINDING);if(fbo!=null){gl.bindFramebuffer(gl.FRAMEBUFFER,fbo);}
stateManager.frontFace(gl.CCW);stateManager.disable(gl.CULL_FACE);stateManager.disable(gl.DEPTH_TEST);var sp=x3dom.SSAO.shaderProgram;stateManager.useProgram(sp);sp.depthTexture=0;sp.randomTexture=1;sp.radius=scene.getEnvironment()._vf.SSAOradius;sp.randomTextureTilingFactor=[canvas.width/x3dom.SSAO.currentRandomTextureSize,canvas.height/x3dom.SSAO.currentRandomTextureSize];var viewpoint=scene.getViewpoint();var nearPlane=viewpoint.getNear();var farPlane=viewpoint.getFar();sp.nearPlane=nearPlane;sp.farPlane=farPlane;sp.depthReconstructionConstantA=(farPlane+nearPlane)/(nearPlane-farPlane);sp.depthReconstructionConstantB=(2.0*farPlane*nearPlane)/(nearPlane-farPlane);sp.depthBufferEpsilon=0.0001*(farPlane-nearPlane);sp.samples=[0.03800223814729654,0.10441029119843426,-0.04479934806797181,-0.03800223814729654,-0.10441029119843426,0.04479934806797181,-0.17023209847088397,0.1428416910414532,0.6154407640895228,0.17023209847088397,-0.1428416910414532,-0.6154407640895228,-0.288675134594813,-0.16666666666666646,-0.3774214123135722,0.288675134594813,0.16666666666666646,0.3774214123135722,0.07717696785196887,-0.43769233467209245,-0.5201284112706428,-0.07717696785196887,0.43769233467209245,0.5201284112706428,0.5471154183401156,-0.09647120981496134,-0.15886420745887797,-0.5471154183401156,0.09647120981496134,0.15886420745887797,0.3333333333333342,0.5773502691896253,-0.8012446019636266,-0.3333333333333342,-0.5773502691896253,0.8012446019636266,-0.49994591864508653,0.5958123446480936,-0.15385106176844343,0.49994591864508653,-0.5958123446480936,0.15385106176844343,-0.8352823295874743,-0.3040179051783715,0.7825440557226517,0.8352823295874743,0.3040179051783715,-0.7825440557226517];if(!sp.tex){sp.tex=0;}
gl.activeTexture(gl.TEXTURE0);gl.bindTexture(gl.TEXTURE_2D,tex);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MAG_FILTER,gl.LINEAR);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MIN_FILTER,gl.LINEAR);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_S,gl.CLAMP_TO_EDGE);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_T,gl.CLAMP_TO_EDGE);gl.activeTexture(gl.TEXTURE1);gl.bindTexture(gl.TEXTURE_2D,x3dom.SSAO.randomTexture);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MAG_FILTER,gl.NEAREST);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MIN_FILTER,gl.NEAREST);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_S,gl.REPEAT);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_T,gl.REPEAT);gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER,scene._fgnd._webgl.buffers[0]);gl.bindBuffer(gl.ARRAY_BUFFER,scene._fgnd._webgl.buffers[1]);gl.vertexAttribPointer(sp.position,3,gl.FLOAT,false,0,0);gl.enableVertexAttribArray(sp.position);gl.drawElements(scene._fgnd._webgl.primType,scene._fgnd._webgl.indexes.length,gl.UNSIGNED_SHORT,0);gl.disableVertexAttribArray(sp.position);gl.activeTexture(gl.TEXTURE0);gl.bindTexture(gl.TEXTURE_2D,null);gl.activeTexture(gl.TEXTURE1);gl.bindTexture(gl.TEXTURE_2D,null);if(fbo!=null){gl.bindFramebuffer(gl.FRAMEBUFFER,oldfbo);}};x3dom.SSAO.blur=function(stateManager,gl,scene,ssaoTexture,depthTexture,canvas,fbo){var oldfbo=gl.getParameter(gl.FRAMEBUFFER_BINDING);if(fbo!=null){gl.bindFramebuffer(gl.FRAMEBUFFER,fbo);}
stateManager.frontFace(gl.CCW);stateManager.disable(gl.CULL_FACE);stateManager.disable(gl.DEPTH_TEST);var sp=x3dom.SSAO.blurShaderProgram;stateManager.useProgram(sp);sp.SSAOTexture=0;sp.depthTexture=1;sp.depthThreshold=scene.getEnvironment()._vf.SSAOblurDepthTreshold;var viewpoint=scene.getViewpoint();var nearPlane=viewpoint.getNear();var farPlane=viewpoint.getFar();sp.nearPlane=nearPlane;sp.farPlane=farPlane;sp.depthReconstructionConstantA=(farPlane+nearPlane)/(nearPlane-farPlane);sp.depthReconstructionConstantB=(2.0*farPlane*nearPlane)/(nearPlane-farPlane);sp.pixelSize=[1.0/canvas.width,1.0/canvas.height];sp.amount=scene.getEnvironment()._vf.SSAOamount;gl.activeTexture(gl.TEXTURE0);gl.bindTexture(gl.TEXTURE_2D,ssaoTexture);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MAG_FILTER,gl.LINEAR);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MIN_FILTER,gl.LINEAR);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_S,gl.CLAMP_TO_EDGE);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_T,gl.CLAMP_TO_EDGE);gl.activeTexture(gl.TEXTURE1);gl.bindTexture(gl.TEXTURE_2D,depthTexture);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MAG_FILTER,gl.NEAREST);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MIN_FILTER,gl.NEAREST);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_S,gl.CLAMP_TO_EDGE);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_T,gl.CLAMP_TO_EDGE);gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER,scene._fgnd._webgl.buffers[0]);gl.bindBuffer(gl.ARRAY_BUFFER,scene._fgnd._webgl.buffers[1]);gl.vertexAttribPointer(sp.position,3,gl.FLOAT,false,0,0);gl.enableVertexAttribArray(sp.position);gl.drawElements(scene._fgnd._webgl.primType,scene._fgnd._webgl.indexes.length,gl.UNSIGNED_SHORT,0);gl.disableVertexAttribArray(sp.position);gl.activeTexture(gl.TEXTURE0);gl.bindTexture(gl.TEXTURE_2D,null);gl.activeTexture(gl.TEXTURE1);gl.bindTexture(gl.TEXTURE_2D,null);if(fbo!=null){gl.bindFramebuffer(gl.FRAMEBUFFER,oldfbo);}};x3dom.SSAO.renderSSAO=function(stateManager,gl,scene,canvas){this.reinitializeShadersIfNecessary(gl);this.reinitializeRandomTextureIfNecessary(gl,scene);this.reinitializeFBOIfNecessary(gl,canvas);stateManager.viewport(0,0,canvas.width,canvas.height);this.render(stateManager,gl,scene,scene._webgl.fboScene.tex,canvas,x3dom.SSAO.fbo);gl.enable(gl.BLEND);gl.blendFunc(gl.DST_COLOR,gl.ONE_MINUS_SRC_ALPHA);this.blur(stateManager,gl,scene,x3dom.SSAO.fbotex,scene._webgl.fboScene.tex,canvas,null);gl.disable(gl.BLEND);};x3dom.gfx_webgl=(function(){"use strict";function Context(ctx3d,canvas,name,x3dElem){this.ctx3d=ctx3d;this.canvas=canvas;this.name=name;this.x3dElem=x3dElem;this.IG_PositionBuffer=null;this.cache=new x3dom.Cache();this.stateManager=new x3dom.StateManager(ctx3d);}
Context.prototype.getName=function(){return this.name;};function setupContext(canvas,forbidMobileShaders,forceMobileShaders,tryWebGL2,x3dElem){var validContextNames=['webgl','experimental-webgl','moz-webgl','webkit-3d'];if(tryWebGL2){validContextNames=['experimental-webgl2'].concat(validContextNames);}
var ctx=null;var envNodes=x3dElem.getElementsByTagName("Environment");var ssaoEnabled=(envNodes&&envNodes[0]&&envNodes[0].hasAttribute("SSAO")&&envNodes[0].getAttribute("SSAO").toLowerCase()==='true')?true:false;var ctxAttribs={alpha:true,depth:true,stencil:true,antialias:!ssaoEnabled,premultipliedAlpha:false,preserveDrawingBuffer:true,failIfMajorPerformanceCaveat:true};for(var i=0;i<validContextNames.length;i++){try{ctx=canvas.getContext(validContextNames[i],ctxAttribs);if(!ctx){x3dom.caps.RENDERMODE="SOFTWARE";ctxAttribs.failIfMajorPerformanceCaveat=false;ctx=canvas.getContext(validContextNames[i],ctxAttribs);}
if(ctx){var newCtx=new Context(ctx,canvas,'webgl',x3dElem);try{x3dom.caps.VENDOR=ctx.getParameter(ctx.VENDOR);x3dom.caps.VERSION=ctx.getParameter(ctx.VERSION);x3dom.caps.RENDERER=ctx.getParameter(ctx.RENDERER);x3dom.caps.SHADING_LANGUAGE_VERSION=ctx.getParameter(ctx.SHADING_LANGUAGE_VERSION);x3dom.caps.RED_BITS=ctx.getParameter(ctx.RED_BITS);x3dom.caps.GREEN_BITS=ctx.getParameter(ctx.GREEN_BITS);x3dom.caps.BLUE_BITS=ctx.getParameter(ctx.BLUE_BITS);x3dom.caps.ALPHA_BITS=ctx.getParameter(ctx.ALPHA_BITS);x3dom.caps.DEPTH_BITS=ctx.getParameter(ctx.DEPTH_BITS);x3dom.caps.MAX_VERTEX_ATTRIBS=ctx.getParameter(ctx.MAX_VERTEX_ATTRIBS);x3dom.caps.MAX_VERTEX_TEXTURE_IMAGE_UNITS=ctx.getParameter(ctx.MAX_VERTEX_TEXTURE_IMAGE_UNITS);x3dom.caps.MAX_VARYING_VECTORS=ctx.getParameter(ctx.MAX_VARYING_VECTORS);x3dom.caps.MAX_VERTEX_UNIFORM_VECTORS=ctx.getParameter(ctx.MAX_VERTEX_UNIFORM_VECTORS);x3dom.caps.MAX_COMBINED_TEXTURE_IMAGE_UNITS=ctx.getParameter(ctx.MAX_COMBINED_TEXTURE_IMAGE_UNITS);x3dom.caps.MAX_TEXTURE_SIZE=ctx.getParameter(ctx.MAX_TEXTURE_SIZE);x3dom.caps.MAX_TEXTURE_IMAGE_UNITS=ctx.getParameter(ctx.MAX_TEXTURE_IMAGE_UNITS);x3dom.caps.MAX_CUBE_MAP_TEXTURE_SIZE=ctx.getParameter(ctx.MAX_CUBE_MAP_TEXTURE_SIZE);x3dom.caps.COMPRESSED_TEXTURE_FORMATS=ctx.getParameter(ctx.COMPRESSED_TEXTURE_FORMATS);x3dom.caps.MAX_RENDERBUFFER_SIZE=ctx.getParameter(ctx.MAX_RENDERBUFFER_SIZE);x3dom.caps.MAX_VIEWPORT_DIMS=ctx.getParameter(ctx.MAX_VIEWPORT_DIMS);x3dom.caps.ALIASED_LINE_WIDTH_RANGE=ctx.getParameter(ctx.ALIASED_LINE_WIDTH_RANGE);x3dom.caps.ALIASED_POINT_SIZE_RANGE=ctx.getParameter(ctx.ALIASED_POINT_SIZE_RANGE);x3dom.caps.SAMPLES=ctx.getParameter(ctx.SAMPLES);x3dom.caps.INDEX_UINT=ctx.getExtension("OES_element_index_uint");x3dom.caps.FP_TEXTURES=ctx.getExtension("OES_texture_float");x3dom.caps.FPL_TEXTURES=ctx.getExtension("OES_texture_float_linear");x3dom.caps.STD_DERIVATIVES=ctx.getExtension("OES_standard_derivatives");x3dom.caps.DRAW_BUFFERS=ctx.getExtension("WEBGL_draw_buffers");x3dom.caps.DEPTH_TEXTURE=ctx.getExtension("WEBGL_depth_texture");x3dom.caps.DEBUGRENDERINFO=ctx.getExtension("WEBGL_debug_renderer_info");x3dom.caps.EXTENSIONS=ctx.getSupportedExtensions();if(x3dom.Utils.isWebGL2Enabled())
{x3dom.caps.DEPTH_TEXTURE=null;}
if(x3dom.caps.DEBUGRENDERINFO){x3dom.caps.UNMASKED_RENDERER_WEBGL=ctx.getParameter(x3dom.caps.DEBUGRENDERINFO.UNMASKED_RENDERER_WEBGL);x3dom.caps.UNMASKED_VENDOR_WEBGL=ctx.getParameter(x3dom.caps.DEBUGRENDERINFO.UNMASKED_VENDOR_WEBGL);}else{x3dom.caps.UNMASKED_RENDERER_WEBGL="";x3dom.caps.UNMASKED_VENDOR_WEBGL="";}
var extString=x3dom.caps.EXTENSIONS.toString().replace(/,/g,", ");x3dom.debug.logInfo(validContextNames[i]+" context found\nVendor: "+x3dom.caps.VENDOR+" "+x3dom.caps.UNMASKED_VENDOR_WEBGL+", Renderer: "+x3dom.caps.RENDERER+" "+x3dom.caps.UNMASKED_RENDERER_WEBGL+", "+"Version: "+x3dom.caps.VERSION+", "+"ShadingLangV.: "+x3dom.caps.SHADING_LANGUAGE_VERSION
+", MSAA samples: "+x3dom.caps.SAMPLES+"\nExtensions: "+extString);if(x3dom.caps.INDEX_UINT){x3dom.Utils.maxIndexableCoords=4294967295;}
x3dom.caps.MOBILE=(function(a){return(/android.+mobile|avantgo|bada\/|blackberry|blazer|compal|elaine|fennec|hiptop|iemobile|ip(hone|od)|iris|kindle|lge |maemo|midp|mmp|opera m(ob|in)i|palm( os)?|phone|p(ixi|re)\/|plucker|pocket|psp|symbian|treo|up\.(browser|link)|vodafone|wap|windows (ce|phone)|xda|xiino/i.test(a)||/1207|6310|6590|3gso|4thp|50[1-6]i|770s|802s|a wa|abac|ac(er|oo|s\-)|ai(ko|rn)|al(av|ca|co)|amoi|an(ex|ny|yw)|aptu|ar(ch|go)|as(te|us)|attw|au(di|\-m|r |s )|avan|be(ck|ll|nq)|bi(lb|rd)|bl(ac|az)|br(e|v)w|bumb|bw\-(n|u)|c55\/|capi|ccwa|cdm\-|cell|chtm|cldc|cmd\-|co(mp|nd)|craw|da(it|ll|ng)|dbte|dc\-s|devi|dica|dmob|do(c|p)o|ds(12|\-d)|el(49|ai)|em(l2|ul)|er(ic|k0)|esl8|ez([4-7]0|os|wa|ze)|fetc|fly(\-|_)|g1 u|g560|gene|gf\-5|g\-mo|go(\.w|od)|gr(ad|un)|haie|hcit|hd\-(m|p|t)|hei\-|hi(pt|ta)|hp( i|ip)|hs\-c|ht(c(\-| |_|a|g|p|s|t)|tp)|hu(aw|tc)|i\-(20|go|ma)|i230|iac( |\-|\/)|ibro|idea|ig01|ikom|im1k|inno|ipaq|iris|ja(t|v)a|jbro|jemu|jigs|kddi|keji|kgt( |\/)|klon|kpt |kwc\-|kyo(c|k)|le(no|xi)|lg( g|\/(k|l|u)|50|54|e\-|e\/|\-[a-w])|libw|lynx|m1\-w|m3ga|m50\/|ma(te|ui|xo)|mc(01|21|ca)|m\-cr|me(di|rc|ri)|mi(o8|oa|ts)|mmef|mo(01|02|bi|de|do|t(\-| |o|v)|zz)|mt(50|p1|v )|mwbp|mywa|n10[0-2]|n20[2-3]|n30(0|2)|n50(0|2|5)|n7(0(0|1)|10)|ne((c|m)\-|on|tf|wf|wg|wt)|nok(6|i)|nzph|o2im|op(ti|wv)|oran|owg1|p800|pan(a|d|t)|pdxg|pg(13|\-([1-8]|c))|phil|pire|pl(ay|uc)|pn\-2|po(ck|rt|se)|prox|psio|pt\-g|qa\-a|qc(07|12|21|32|60|\-[2-7]|i\-)|qtek|r380|r600|raks|rim9|ro(ve|zo)|s55\/|sa(ge|ma|mm|ms|ny|va)|sc(01|h\-|oo|p\-)|sdk\/|se(c(\-|0|1)|47|mc|nd|ri)|sgh\-|shar|sie(\-|m)|sk\-0|sl(45|id)|sm(al|ar|b3|it|t5)|so(ft|ny)|sp(01|h\-|v\-|v )|sy(01|mb)|t2(18|50)|t6(00|10|18)|ta(gt|lk)|tcl\-|tdg\-|tel(i|m)|tim\-|t\-mo|to(pl|sh)|ts(70|m\-|m3|m5)|tx\-9|up(\.b|g1|si)|utst|v400|v750|veri|vi(rg|te)|vk(40|5[0-3]|\-v)|vm40|voda|vulc|vx(52|53|60|61|70|80|81|83|85|98)|w3c(\-| )|webc|whit|wi(g |nc|nw)|wmlb|wonu|x700|xda(\-|2|g)|yas\-|your|zeto|zte\-/i.test(a.substr(0,4)))})(navigator.userAgent||navigator.vendor||window.opera);if(x3dom.caps.RENDERER.indexOf("PowerVR")>=0||navigator.appVersion.indexOf("Mobile")>-1||x3dom.caps.MAX_VARYING_VECTORS<=8||x3dom.caps.MAX_VERTEX_TEXTURE_IMAGE_UNITS<2){x3dom.caps.MOBILE=true;}
if(x3dom.caps.MOBILE){if(forbidMobileShaders){x3dom.caps.MOBILE=false;x3dom.debug.logWarning("Detected mobile graphics card! "+"But being forced to desktop shaders which might not work!");}
else{x3dom.debug.logWarning("Detected mobile graphics card! "+"Using low quality shaders without ImageGeometry support!");}}
else{if(forceMobileShaders){x3dom.caps.MOBILE=true;x3dom.debug.logWarning("Detected desktop graphics card! "+"But being forced to mobile shaders with lower quality!");}}}
catch(ex){x3dom.debug.logWarning("Your browser probably supports an older WebGL version. "+"Please try the old mobile runtime instead:\n"+"http://www.x3dom.org/x3dom/src_mobile/x3dom.js");newCtx=null;}
return newCtx;}}
catch(e){x3dom.debug.logWarning(e);}}
return null;}
Context.prototype.setupShape=function(gl,drawable,viewarea){var q=0,q6;var textures,t;var vertices,positionBuffer;var texCoordBuffer,normalBuffer,colorBuffer;var indicesBuffer,indexArray;var shape=drawable.shape;var geoNode=shape._cf.geometry.node;if(shape._webgl!==undefined){var needFullReInit=false;if(shape._dirty.colors===true&&shape._webgl.shader.color===undefined&&geoNode._mesh._colors[0].length){needFullReInit=true;}
if(needFullReInit&&shape._cleanupGLObjects){shape._cleanupGLObjects(true,false);}
if(shape._dirty.texture===true){if(shape._webgl.texture.length!=shape.getTextures().length){for(t=0;t<shape._webgl.texture.length;++t){shape._webgl.texture.pop();}
textures=shape.getTextures();for(t=0;t<textures.length;++t){shape._webgl.texture.push(new x3dom.Texture(gl,shape._nameSpace.doc,this.cache,textures[t]));}
shape._dirty.shader=true;if(shape._webgl.shader.texcoord===undefined)
shape._dirty.texcoords=true;}
else{textures=shape.getTextures();for(t=0;t<textures.length;++t){if(textures[t]===shape._webgl.texture[t].node){shape._webgl.texture[t].update();}
else{shape._webgl.texture[t].texture=null;shape._webgl.texture[t].node=textures[t];shape._webgl.texture[t].update();}}}
shape._dirty.texture=false;}
shape._webgl.shader=this.cache.getShaderByProperties(gl,shape,shape.getShaderProperties(viewarea));if(!needFullReInit&&shape._webgl.binaryGeometry==0)
{for(q=0;q<shape._webgl.positions.length;q++)
{q6=6*q;if(shape._dirty.positions==true||shape._dirty.indexes==true){if(shape._webgl.shader.position!==undefined){shape._webgl.indexes[q]=geoNode._mesh._indices[q];gl.deleteBuffer(shape._webgl.buffers[q6]);indicesBuffer=gl.createBuffer();shape._webgl.buffers[q6]=indicesBuffer;if(x3dom.caps.INDEX_UINT&&(geoNode._mesh._positions[0].length/3>65535)){indexArray=new Uint32Array(shape._webgl.indexes[q]);shape._webgl.indexType=gl.UNSIGNED_INT;}
else{indexArray=new Uint16Array(shape._webgl.indexes[q]);shape._webgl.indexType=gl.UNSIGNED_SHORT;}
gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER,indicesBuffer);gl.bufferData(gl.ELEMENT_ARRAY_BUFFER,indexArray,gl.STATIC_DRAW);indexArray=null;shape._webgl.positions[q]=geoNode._mesh._positions[q];gl.deleteBuffer(shape._webgl.buffers[q6+1]);positionBuffer=gl.createBuffer();shape._webgl.buffers[q6+1]=positionBuffer;gl.bindBuffer(gl.ARRAY_BUFFER,positionBuffer);gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER,shape._webgl.buffers[q6]);vertices=new Float32Array(shape._webgl.positions[q]);gl.bufferData(gl.ARRAY_BUFFER,vertices,gl.STATIC_DRAW);gl.bindBuffer(gl.ARRAY_BUFFER,positionBuffer);gl.vertexAttribPointer(shape._webgl.shader.position,geoNode._mesh._numPosComponents,shape._webgl.coordType,false,shape._coordStrideOffset[0],shape._coordStrideOffset[1]);vertices=null;}
shape._dirty.positions=false;shape._dirty.indexes=false;}
if(shape._dirty.colors==true){if(shape._webgl.shader.color!==undefined){shape._webgl.colors[q]=geoNode._mesh._colors[q];gl.deleteBuffer(shape._webgl.buffers[q6+4]);colorBuffer=gl.createBuffer();shape._webgl.buffers[q6+4]=colorBuffer;colors=new Float32Array(shape._webgl.colors[q]);gl.bindBuffer(gl.ARRAY_BUFFER,colorBuffer);gl.bufferData(gl.ARRAY_BUFFER,colors,gl.STATIC_DRAW);gl.vertexAttribPointer(shape._webgl.shader.color,geoNode._mesh._numColComponents,shape._webgl.colorType,false,shape._colorStrideOffset[0],shape._colorStrideOffset[1]);colors=null;}
shape._dirty.colors=false;}
if(shape._dirty.normals==true){if(shape._webgl.shader.normal!==undefined){shape._webgl.normals[q]=geoNode._mesh._normals[q];gl.deleteBuffer(shape._webgl.buffers[q6+2]);normalBuffer=gl.createBuffer();shape._webgl.buffers[q6+2]=normalBuffer;normals=new Float32Array(shape._webgl.normals[q]);gl.bindBuffer(gl.ARRAY_BUFFER,normalBuffer);gl.bufferData(gl.ARRAY_BUFFER,normals,gl.STATIC_DRAW);gl.vertexAttribPointer(shape._webgl.shader.normal,geoNode._mesh._numNormComponents,shape._webgl.normalType,false,shape._normalStrideOffset[0],shape._normalStrideOffset[1]);normals=null;}
shape._dirty.normals=false;}
if(shape._dirty.texcoords==true){if(shape._webgl.shader.texcoord!==undefined){shape._webgl.texcoords[q]=geoNode._mesh._texCoords[q];gl.deleteBuffer(shape._webgl.buffers[q6+3]);texCoordBuffer=gl.createBuffer();shape._webgl.buffers[q6+3]=texCoordBuffer;texCoords=new Float32Array(shape._webgl.texcoords[q]);gl.bindBuffer(gl.ARRAY_BUFFER,texCoordBuffer);gl.bufferData(gl.ARRAY_BUFFER,texCoords,gl.STATIC_DRAW);gl.vertexAttribPointer(shape._webgl.shader.texCoord,geoNode._mesh._numTexComponents,shape._webgl.texCoordType,false,shape._texCoordStrideOffset[0],shape._texCoordStrideOffset[1]);texCoords=null;}
shape._dirty.texcoords=false;}
if(shape._dirty.specialAttribs==true){if(shape._webgl.shader.particleSize!==undefined){var szArr=geoNode._vf.size.toGL();if(szArr.length){gl.deleteBuffer(shape._webgl.buffers[q6+5]);shape._webgl.buffers[q6+5]=gl.createBuffer();gl.bindBuffer(gl.ARRAY_BUFFER,shape._webgl.buffers[q6+5]);gl.bufferData(gl.ARRAY_BUFFER,new Float32Array(szArr),gl.STATIC_DRAW);}
shape._dirty.specialAttribs=false;}}}}
else
{}
if(shape._webgl.imageGeometry!=0){for(t=0;t<shape._webgl.texture.length;++t){shape._webgl.texture[t].updateTexture();}
geoNode.unsetGeoDirty();shape.unsetGeoDirty();}
if(!needFullReInit){return;}}
else if(!(x3dom.isa(geoNode,x3dom.nodeTypes.Text)||x3dom.isa(geoNode,x3dom.nodeTypes.BinaryGeometry)||x3dom.isa(geoNode,x3dom.nodeTypes.PopGeometry)||x3dom.isa(geoNode,x3dom.nodeTypes.ExternalGeometry)||x3dom.isa(shape,x3dom.nodeTypes.ExternalShape))&&(!geoNode||geoNode._mesh._positions[0].length<1))
{if(x3dom.caps.MAX_VERTEX_TEXTURE_IMAGE_UNITS<2&&x3dom.isa(geoNode,x3dom.nodeTypes.ImageGeometry)){x3dom.debug.logError("Can't render ImageGeometry nodes with only "+
x3dom.caps.MAX_VERTEX_TEXTURE_IMAGE_UNITS+" vertex texture units. Please upgrade your GPU!");}
else{x3dom.debug.logError("NO VALID MESH OR NO VERTEX POSITIONS SET!");}
return;}
shape.unsetDirty();if(!shape._cleanupGLObjects)
{shape._cleanupGLObjects=function(force,delGL)
{if(this._webgl&&((arguments.length>0&&force)||this._parentNodes.length==0))
{var sp=this._webgl.shader;for(var q=0;q<this._webgl.positions.length;q++){var q6=6*q;if(sp.position!==undefined){gl.deleteBuffer(this._webgl.buffers[q6+1]);gl.deleteBuffer(this._webgl.buffers[q6]);}
if(sp.normal!==undefined){gl.deleteBuffer(this._webgl.buffers[q6+2]);}
if(sp.texcoord!==undefined){gl.deleteBuffer(this._webgl.buffers[q6+3]);}
if(sp.color!==undefined){gl.deleteBuffer(this._webgl.buffers[q6+4]);}
if(sp.id!==undefined){gl.deleteBuffer(this._webgl.buffers[q6+5]);}}
for(var df=0;df<this._webgl.dynamicFields.length;df++){var attrib=this._webgl.dynamicFields[df];if(sp[attrib.name]!==undefined){gl.deleteBuffer(attrib.buf);}}
if(delGL===undefined)
delGL=true;if(delGL){delete this._webgl;x3dom.BinaryContainerLoader.outOfMemory=false;}}};}
shape._webgl={positions:geoNode._mesh._positions,normals:geoNode._mesh._normals,texcoords:geoNode._mesh._texCoords,colors:geoNode._mesh._colors,indexes:geoNode._mesh._indices,indexType:gl.UNSIGNED_SHORT,coordType:gl.FLOAT,normalType:gl.FLOAT,texCoordType:gl.FLOAT,colorType:gl.FLOAT,texture:[],dirtyLighting:x3dom.Utils.checkDirtyLighting(viewarea),imageGeometry:0,binaryGeometry:0,popGeometry:0,externalGeometry:0};textures=shape.getTextures();for(t=0;t<textures.length;++t){shape._webgl.texture.push(new x3dom.Texture(gl,shape._nameSpace.doc,this.cache,textures[t]));}
shape._webgl.shader=this.cache.getShaderByProperties(gl,shape,shape.getShaderProperties(viewarea));var sp=shape._webgl.shader;var currAttribs=0;shape._webgl.buffers=[];shape._webgl.dynamicFields=[];if(x3dom.isa(geoNode,x3dom.nodeTypes.X3DBinaryContainerGeometryNode))
{shape._webgl.primType=[];for(var primCnt=0;primCnt<geoNode._vf.primType.length;++primCnt)
{shape._webgl.primType.push(x3dom.Utils.primTypeDic(gl,geoNode._vf.primType[primCnt]));}}
else
{shape._webgl.primType=x3dom.Utils.primTypeDic(gl,geoNode._mesh._primType);}
if(x3dom.isa(geoNode,x3dom.nodeTypes.ExternalGeometry))
{geoNode.update(shape,sp,gl,viewarea,this);}
else if(x3dom.isa(shape,x3dom.nodeTypes.ExternalShape))
{shape.update(shape,sp,gl,viewarea,this);}
else if(x3dom.isa(geoNode,x3dom.nodeTypes.BinaryGeometry))
{x3dom.BinaryContainerLoader.setupBinGeo(shape,sp,gl,viewarea,this);}
else if(x3dom.isa(geoNode,x3dom.nodeTypes.PopGeometry))
{x3dom.BinaryContainerLoader.setupPopGeo(shape,sp,gl,viewarea,this);}
else if(x3dom.isa(geoNode,x3dom.nodeTypes.ImageGeometry))
{x3dom.BinaryContainerLoader.setupImgGeo(shape,sp,gl,viewarea,this);}
else
{for(q=0;q<shape._webgl.positions.length;q++)
{q6=6*q;if(sp.position!==undefined){indicesBuffer=gl.createBuffer();shape._webgl.buffers[q6]=indicesBuffer;if(x3dom.caps.INDEX_UINT&&(shape._webgl.positions[0].length/3>65535)){indexArray=new Uint32Array(shape._webgl.indexes[q]);shape._webgl.indexType=gl.UNSIGNED_INT;}
else{indexArray=new Uint16Array(shape._webgl.indexes[q]);shape._webgl.indexType=gl.UNSIGNED_SHORT;}
gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER,indicesBuffer);gl.bufferData(gl.ELEMENT_ARRAY_BUFFER,indexArray,gl.STATIC_DRAW);indexArray=null;positionBuffer=gl.createBuffer();shape._webgl.buffers[q6+1]=positionBuffer;gl.bindBuffer(gl.ARRAY_BUFFER,positionBuffer);vertices=new Float32Array(shape._webgl.positions[q]);gl.bufferData(gl.ARRAY_BUFFER,vertices,gl.STATIC_DRAW);gl.bindBuffer(gl.ARRAY_BUFFER,positionBuffer);gl.vertexAttribPointer(sp.position,geoNode._mesh._numPosComponents,shape._webgl.coordType,false,shape._coordStrideOffset[0],shape._coordStrideOffset[1]);gl.enableVertexAttribArray(sp.position);vertices=null;}
if(sp.normal!==undefined||shape._webgl.normals[q]){normalBuffer=gl.createBuffer();shape._webgl.buffers[q6+2]=normalBuffer;var normals=new Float32Array(shape._webgl.normals[q]);gl.bindBuffer(gl.ARRAY_BUFFER,normalBuffer);gl.bufferData(gl.ARRAY_BUFFER,normals,gl.STATIC_DRAW);gl.vertexAttribPointer(sp.normal,geoNode._mesh._numNormComponents,shape._webgl.normalType,false,shape._normalStrideOffset[0],shape._normalStrideOffset[1]);gl.enableVertexAttribArray(sp.normal);normals=null;}
if(sp.texcoord!==undefined){var texcBuffer=gl.createBuffer();shape._webgl.buffers[q6+3]=texcBuffer;var texCoords=new Float32Array(shape._webgl.texcoords[q]);gl.bindBuffer(gl.ARRAY_BUFFER,texcBuffer);gl.bufferData(gl.ARRAY_BUFFER,texCoords,gl.STATIC_DRAW);gl.vertexAttribPointer(sp.texcoord,geoNode._mesh._numTexComponents,shape._webgl.texCoordType,false,shape._texCoordStrideOffset[0],shape._texCoordStrideOffset[1]);gl.enableVertexAttribArray(sp.texcoord);texCoords=null;}
if(sp.color!==undefined){colorBuffer=gl.createBuffer();shape._webgl.buffers[q6+4]=colorBuffer;var colors=new Float32Array(shape._webgl.colors[q]);gl.bindBuffer(gl.ARRAY_BUFFER,colorBuffer);gl.bufferData(gl.ARRAY_BUFFER,colors,gl.STATIC_DRAW);gl.vertexAttribPointer(sp.color,geoNode._mesh._numColComponents,shape._webgl.colorType,false,shape._colorStrideOffset[0],shape._colorStrideOffset[1]);gl.enableVertexAttribArray(sp.color);colors=null;}
if(sp.particleSize!==undefined){var sizeArr=geoNode._vf.size.toGL();if(sizeArr.length){var sizeBuffer=gl.createBuffer();shape._webgl.buffers[q6+5]=sizeBuffer;gl.bindBuffer(gl.ARRAY_BUFFER,sizeBuffer);gl.bufferData(gl.ARRAY_BUFFER,new Float32Array(sizeArr),gl.STATIC_DRAW);}}}
for(var df in geoNode._mesh._dynamicFields)
{if(!geoNode._mesh._dynamicFields.hasOwnProperty(df))
continue;var attrib=geoNode._mesh._dynamicFields[df];shape._webgl.dynamicFields[currAttribs]={buf:{},name:df,numComponents:attrib.numComponents};if(sp[df]!==undefined){var attribBuffer=gl.createBuffer();shape._webgl.dynamicFields[currAttribs++].buf=attribBuffer;var attribs=new Float32Array(attrib.value);gl.bindBuffer(gl.ARRAY_BUFFER,attribBuffer);gl.bufferData(gl.ARRAY_BUFFER,attribs,gl.STATIC_DRAW);gl.vertexAttribPointer(sp[df],attrib.numComponents,gl.FLOAT,false,0,0);attribs=null;}}}};Context.prototype.setupScene=function(gl,bgnd){var sphere=null;var texture=null;var that=this;if(bgnd._webgl!==undefined){if(!bgnd._dirty){return;}
if(bgnd._webgl.texture!==undefined&&bgnd._webgl.texture){gl.deleteTexture(bgnd._webgl.texture);}
if(bgnd._cleanupGLObjects){bgnd._cleanupGLObjects();}
bgnd._webgl={};}
bgnd._dirty=false;var url=bgnd.getTexUrl();var i=0;var w=1,h=1;if(url.length>0&&url[0].length>0){if(url.length>=6&&url[1].length>0&&url[2].length>0&&url[3].length>0&&url[4].length>0&&url[5].length>0){sphere=new x3dom.nodeTypes.Sphere();bgnd._webgl={positions:sphere._mesh._positions[0],indexes:sphere._mesh._indices[0],buffers:[{},{}]};bgnd._webgl.primType=gl.TRIANGLES;bgnd._webgl.shader=this.cache.getShader(gl,x3dom.shader.BACKGROUND_CUBETEXTURE);bgnd._webgl.texture=x3dom.Utils.createTextureCube(gl,bgnd._nameSpace.doc,url,true,bgnd._vf.crossOrigin,true,false);}
else{bgnd._webgl={positions:[-w,-h,0,-w,h,0,w,-h,0,w,h,0],indexes:[0,1,2,3],buffers:[{},{}]};url=bgnd._nameSpace.getURL(url[0]);bgnd._webgl.texture=x3dom.Utils.createTexture2D(gl,bgnd._nameSpace.doc,url,true,bgnd._vf.crossOrigin,false,false);bgnd._webgl.primType=gl.TRIANGLE_STRIP;bgnd._webgl.shader=this.cache.getShader(gl,x3dom.shader.BACKGROUND_TEXTURE);}}
else{if(bgnd.getSkyColor().length>1||bgnd.getGroundColor().length){sphere=new x3dom.nodeTypes.Sphere();texture=gl.createTexture();bgnd._webgl={positions:sphere._mesh._positions[0],texcoords:sphere._mesh._texCoords[0],indexes:sphere._mesh._indices[0],buffers:[{},{},{}],texture:texture,primType:gl.TRIANGLES};var N=x3dom.Utils.nextHighestPowerOfTwo(bgnd.getSkyColor().length+bgnd.getGroundColor().length+2);N=(N<512)?512:N;var n=bgnd._vf.groundAngle.length;var tmp=[],arr=[];var colors=[],sky=[0];for(i=0;i<bgnd._vf.skyColor.length;i++){colors[i]=bgnd._vf.skyColor[i];}
for(i=0;i<bgnd._vf.skyAngle.length;i++){sky[i+1]=bgnd._vf.skyAngle[i];}
if(n>0||bgnd._vf.groundColor.length==1){if(sky[sky.length-1]<Math.PI/2){sky[sky.length]=Math.PI/2-x3dom.fields.Eps;colors[colors.length]=colors[colors.length-1];}
for(i=n-1;i>=0;i--){if((i==n-1)&&(Math.PI-bgnd._vf.groundAngle[i]<=Math.PI/2)){sky[sky.length]=Math.PI/2;colors[colors.length]=bgnd._vf.groundColor[bgnd._vf.groundColor.length-1];}
sky[sky.length]=Math.PI-bgnd._vf.groundAngle[i];colors[colors.length]=bgnd._vf.groundColor[i+1];}
if(n==0&&bgnd._vf.groundColor.length==1){sky[sky.length]=Math.PI/2;colors[colors.length]=bgnd._vf.groundColor[0];}
sky[sky.length]=Math.PI;colors[colors.length]=bgnd._vf.groundColor[0];}
else{if(sky[sky.length-1]<Math.PI){sky[sky.length]=Math.PI;colors[colors.length]=colors[colors.length-1];}}
for(i=0;i<sky.length;i++){sky[i]/=Math.PI;}
if(sky.length!=colors.length){x3dom.debug.logError("Number of background colors and corresponding angles are different!");var minArrayLength=(sky.length<colors.length)?sky.length:colors.length;sky.length=minArrayLength;colors.length=minArrayLength;}
var interp=new x3dom.nodeTypes.ColorInterpolator();interp._vf.key=new x3dom.fields.MFFloat(sky);interp._vf.keyValue=new x3dom.fields.MFColor(colors);for(i=0;i<N;i++){interp._vf.set_fraction=i/(N-1.0);interp.fieldChanged("set_fraction");tmp[i]=interp._vf.value_changed;}
tmp.reverse();var alpha=Math.floor((1.0-bgnd.getTransparency())*255);for(i=0;i<tmp.length;i++){arr.push(Math.floor(tmp[i].r*255),Math.floor(tmp[i].g*255),Math.floor(tmp[i].b*255),alpha);}
var pixels=new Uint8Array(arr);var format=gl.RGBA;N=pixels.length/4;gl.bindTexture(gl.TEXTURE_2D,texture);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MIN_FILTER,gl.LINEAR);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MAG_FILTER,gl.LINEAR);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_S,gl.CLAMP_TO_EDGE);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_T,gl.CLAMP_TO_EDGE);gl.pixelStorei(gl.UNPACK_ALIGNMENT,1);gl.texImage2D(gl.TEXTURE_2D,0,format,1,N,0,format,gl.UNSIGNED_BYTE,pixels);gl.bindTexture(gl.TEXTURE_2D,null);bgnd._webgl.shader=this.cache.getShader(gl,x3dom.shader.BACKGROUND_SKYTEXTURE);}
else{bgnd._webgl={};}}
if(bgnd._webgl.shader){var sp=bgnd._webgl.shader;var positionBuffer=gl.createBuffer();bgnd._webgl.buffers[1]=positionBuffer;gl.bindBuffer(gl.ARRAY_BUFFER,positionBuffer);var vertices=new Float32Array(bgnd._webgl.positions);gl.bufferData(gl.ARRAY_BUFFER,vertices,gl.STATIC_DRAW);gl.bindBuffer(gl.ARRAY_BUFFER,positionBuffer);gl.vertexAttribPointer(sp.position,3,gl.FLOAT,false,0,0);gl.enableVertexAttribArray(sp.position);var indicesBuffer=gl.createBuffer();bgnd._webgl.buffers[0]=indicesBuffer;var indexArray=new Uint16Array(bgnd._webgl.indexes);gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER,indicesBuffer);gl.bufferData(gl.ELEMENT_ARRAY_BUFFER,indexArray,gl.STATIC_DRAW);vertices=null;indexArray=null;if(sp.texcoord!==undefined){var texcBuffer=gl.createBuffer();bgnd._webgl.buffers[2]=texcBuffer;var texcoords=new Float32Array(bgnd._webgl.texcoords);gl.bindBuffer(gl.ARRAY_BUFFER,texcBuffer);gl.bufferData(gl.ARRAY_BUFFER,texcoords,gl.STATIC_DRAW);gl.vertexAttribPointer(sp.texcoord,2,gl.FLOAT,false,0,0);gl.enableVertexAttribArray(sp.texcoord);texcoords=null;}
bgnd._cleanupGLObjects=function(){var sp=this._webgl.shader;if(sp.position!==undefined){gl.deleteBuffer(this._webgl.buffers[0]);gl.deleteBuffer(this._webgl.buffers[1]);}
if(sp.texcoord!==undefined){gl.deleteBuffer(this._webgl.buffers[2]);}};}
bgnd._webgl.render=function(gl,mat_view,mat_proj)
{var sp=bgnd._webgl.shader;var alpha=1.0-bgnd.getTransparency();var mat_scene=null;var projMatrix_22=mat_proj._22,projMatrix_23=mat_proj._23;var camPos=mat_view.e3();if((sp!==undefined&&sp!==null)&&(sp.texcoord!==undefined&&sp.texcoord!==null)&&(bgnd._webgl.texture!==undefined&&bgnd._webgl.texture!==null)){gl.clearColor(0,0,0,alpha);gl.clearDepth(1.0);gl.clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER_BIT|gl.STENCIL_BUFFER_BIT);that.stateManager.frontFace(gl.CCW);that.stateManager.disable(gl.CULL_FACE);that.stateManager.disable(gl.DEPTH_TEST);that.stateManager.disable(gl.BLEND);that.stateManager.useProgram(sp);if(!sp.tex){sp.tex=0;}
mat_proj._22=100001/99999;mat_proj._23=200000/99999;mat_view._03=0;mat_view._13=0;mat_view._23=0;mat_scene=mat_proj.mult(mat_view);sp.modelViewProjectionMatrix=mat_scene.toGL();mat_view._03=camPos.x;mat_view._13=camPos.y;mat_view._23=camPos.z;mat_proj._22=projMatrix_22;mat_proj._23=projMatrix_23;gl.activeTexture(gl.TEXTURE0);gl.bindTexture(gl.TEXTURE_2D,bgnd._webgl.texture);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MIN_FILTER,gl.LINEAR);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MAG_FILTER,gl.LINEAR);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_S,gl.CLAMP_TO_EDGE);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_T,gl.CLAMP_TO_EDGE);gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER,bgnd._webgl.buffers[0]);gl.bindBuffer(gl.ARRAY_BUFFER,bgnd._webgl.buffers[1]);gl.vertexAttribPointer(sp.position,3,gl.FLOAT,false,0,0);gl.enableVertexAttribArray(sp.position);gl.bindBuffer(gl.ARRAY_BUFFER,bgnd._webgl.buffers[2]);gl.vertexAttribPointer(sp.texcoord,2,gl.FLOAT,false,0,0);gl.enableVertexAttribArray(sp.texcoord);gl.drawElements(bgnd._webgl.primType,bgnd._webgl.indexes.length,gl.UNSIGNED_SHORT,0);gl.activeTexture(gl.TEXTURE0);gl.bindTexture(gl.TEXTURE_2D,null);gl.disableVertexAttribArray(sp.position);gl.disableVertexAttribArray(sp.texcoord);gl.clear(gl.DEPTH_BUFFER_BIT);}
else if(!sp||!bgnd._webgl.texture||(bgnd._webgl.texture.textureCubeReady!==undefined&&bgnd._webgl.texture.textureCubeReady!==true)){var bgCol=bgnd.getSkyColor().toGL();gl.clearColor(bgCol[0],bgCol[1],bgCol[2],alpha);gl.clearDepth(1.0);gl.clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER_BIT|gl.STENCIL_BUFFER_BIT);}
else{gl.clearColor(0,0,0,alpha);gl.clearDepth(1.0);gl.clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER_BIT|gl.STENCIL_BUFFER_BIT);that.stateManager.frontFace(gl.CCW);that.stateManager.disable(gl.CULL_FACE);that.stateManager.disable(gl.DEPTH_TEST);that.stateManager.disable(gl.BLEND);that.stateManager.useProgram(sp);if(!sp.tex){sp.tex=0;}
if(bgnd._webgl.texture.textureCubeReady){mat_proj._22=100001/99999;mat_proj._23=200000/99999;mat_view._03=0;mat_view._13=0;mat_view._23=0;mat_scene=mat_proj.mult(mat_view);sp.modelViewProjectionMatrix=mat_scene.toGL();mat_view._03=camPos.x;mat_view._13=camPos.y;mat_view._23=camPos.z;mat_proj._22=projMatrix_22;mat_proj._23=projMatrix_23;gl.activeTexture(gl.TEXTURE0);gl.bindTexture(gl.TEXTURE_CUBE_MAP,bgnd._webgl.texture);gl.texParameteri(gl.TEXTURE_CUBE_MAP,gl.TEXTURE_MIN_FILTER,gl.LINEAR);gl.texParameteri(gl.TEXTURE_CUBE_MAP,gl.TEXTURE_MAG_FILTER,gl.LINEAR);gl.texParameteri(gl.TEXTURE_CUBE_MAP,gl.TEXTURE_WRAP_S,gl.CLAMP_TO_EDGE);gl.texParameteri(gl.TEXTURE_CUBE_MAP,gl.TEXTURE_WRAP_T,gl.CLAMP_TO_EDGE);}
else{gl.activeTexture(gl.TEXTURE0);gl.bindTexture(gl.TEXTURE_2D,bgnd._webgl.texture);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MIN_FILTER,gl.LINEAR);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MAG_FILTER,gl.LINEAR);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_S,gl.CLAMP_TO_EDGE);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_T,gl.CLAMP_TO_EDGE);if(bgnd._vf.scaling&&bgnd._webgl.texture.ready)
{var ratio=1.0;var viewport=new x3dom.fields.SFVec2f(that.canvas.width,that.canvas.height);var texture=new x3dom.fields.SFVec2f(bgnd._webgl.texture.width,bgnd._webgl.texture.height);if(viewport.x>viewport.y)
{ratio=viewport.x/texture.x
texture.x=viewport.x;texture.y=texture.y*ratio;}
else
{ratio=viewport.y/texture.y
texture.y=viewport.y;texture.x=texture.x*ratio;}
var scale=viewport.divideComponents(texture);var translation=texture.subtract(viewport).multiply(0.5).divideComponents(texture);}
else
{var scale=new x3dom.fields.SFVec2f(1.0,1.0);var translation=new x3dom.fields.SFVec2f(0.0,0.0);}
sp.scale=scale.toGL();sp.translation=translation.toGL();}
gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER,bgnd._webgl.buffers[0]);gl.bindBuffer(gl.ARRAY_BUFFER,bgnd._webgl.buffers[1]);gl.vertexAttribPointer(sp.position,3,gl.FLOAT,false,0,0);gl.enableVertexAttribArray(sp.position);gl.drawElements(bgnd._webgl.primType,bgnd._webgl.indexes.length,gl.UNSIGNED_SHORT,0);gl.disableVertexAttribArray(sp.position);gl.activeTexture(gl.TEXTURE0);if(bgnd._webgl.texture.textureCubeReady){gl.bindTexture(gl.TEXTURE_CUBE_MAP,null);}
else{gl.bindTexture(gl.TEXTURE_2D,null);}
gl.clear(gl.DEPTH_BUFFER_BIT);}};};Context.prototype.setupFgnds=function(gl,scene){if(scene._fgnd!==undefined){return;}
var that=this;var w=1,h=1;scene._fgnd={};scene._fgnd._webgl={positions:[-w,-h,0,-w,h,0,w,-h,0,w,h,0],indexes:[0,1,2,3],buffers:[{},{}]};scene._fgnd._webgl.primType=gl.TRIANGLE_STRIP;scene._fgnd._webgl.shader=this.cache.getShader(gl,x3dom.shader.FRONTGROUND_TEXTURE);var sp=scene._fgnd._webgl.shader;var positionBuffer=gl.createBuffer();scene._fgnd._webgl.buffers[1]=positionBuffer;gl.bindBuffer(gl.ARRAY_BUFFER,positionBuffer);var vertices=new Float32Array(scene._fgnd._webgl.positions);gl.bufferData(gl.ARRAY_BUFFER,vertices,gl.STATIC_DRAW);gl.bindBuffer(gl.ARRAY_BUFFER,positionBuffer);gl.vertexAttribPointer(sp.position,3,gl.FLOAT,false,0,0);var indicesBuffer=gl.createBuffer();scene._fgnd._webgl.buffers[0]=indicesBuffer;var indexArray=new Uint16Array(scene._fgnd._webgl.indexes);gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER,indicesBuffer);gl.bufferData(gl.ELEMENT_ARRAY_BUFFER,indexArray,gl.STATIC_DRAW);vertices=null;indexArray=null;scene._fgnd._webgl.render=function(gl,tex){scene._fgnd._webgl.texture=tex;that.stateManager.frontFace(gl.CCW);that.stateManager.disable(gl.CULL_FACE);that.stateManager.disable(gl.DEPTH_TEST);that.stateManager.useProgram(sp);if(!sp.tex){sp.tex=0;}
gl.activeTexture(gl.TEXTURE0);gl.bindTexture(gl.TEXTURE_2D,scene._fgnd._webgl.texture);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MAG_FILTER,gl.NEAREST);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MIN_FILTER,gl.NEAREST);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_S,gl.CLAMP_TO_EDGE);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_T,gl.CLAMP_TO_EDGE);gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER,scene._fgnd._webgl.buffers[0]);gl.bindBuffer(gl.ARRAY_BUFFER,scene._fgnd._webgl.buffers[1]);gl.vertexAttribPointer(sp.position,3,gl.FLOAT,false,0,0);gl.enableVertexAttribArray(sp.position);gl.drawElements(scene._fgnd._webgl.primType,scene._fgnd._webgl.indexes.length,gl.UNSIGNED_SHORT,0);gl.disableVertexAttribArray(sp.position);gl.activeTexture(gl.TEXTURE0);gl.bindTexture(gl.TEXTURE_2D,null);};};Context.prototype.renderShadowPass=function(gl,viewarea,mat_scene,mat_view,targetFbo,camOffset,isCameraView)
{var scene=viewarea._scene;var indicesReady=false;this.stateManager.bindFramebuffer(gl.FRAMEBUFFER,targetFbo.fbo);this.stateManager.viewport(0,0,targetFbo.width,targetFbo.height);gl.clearColor(1.0,1.0,1.0,0.0);gl.clearDepth(1.0);gl.clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER_BIT);this.stateManager.depthFunc(gl.LEQUAL);this.stateManager.enable(gl.DEPTH_TEST);this.stateManager.enable(gl.CULL_FACE);this.stateManager.disable(gl.BLEND);var bgCenter=x3dom.fields.SFVec3f.NullVector.toGL();var bgSize=x3dom.fields.SFVec3f.OneVector.toGL();var env=scene.getEnvironment();var excludeTrans=env._vf.shadowExcludeTransparentObjects;var i,n=scene.drawableCollection.length;for(i=0;i<n;i++)
{var drawable=scene.drawableCollection.get(i);var trafo=drawable.transform;var shape=drawable.shape;var s_gl=shape._webgl;if(!s_gl||(excludeTrans&&drawable.sortType=='transparent')){continue;}
var s_geo=shape._cf.geometry.node;var s_app=shape._cf.appearance.node;var s_msh=s_geo._mesh;var properties=shape.getShaderProperties(viewarea);var sp=this.cache.getShaderByProperties(gl,shape,properties,null,true);if(!sp){return;}
this.stateManager.useProgram(sp);sp.cameraView=isCameraView;sp.offset=camOffset;sp.modelViewProjectionMatrix=mat_scene.mult(trafo).toGL();if(s_gl.coordType!=gl.FLOAT){if(!s_gl.popGeometry&&(x3dom.Utils.isUnsignedType(s_geo._vf.coordType))){sp.bgCenter=s_geo.getMin().toGL();}
else{sp.bgCenter=s_geo._vf.position.toGL();}
sp.bgSize=s_geo._vf.size.toGL();sp.bgPrecisionMax=s_geo.getPrecisionMax('coordType');}
if(shape._clipPlanes){sp.modelViewMatrix=mat_view.mult(trafo).toGL();sp.viewMatrixInverse=mat_view.inverse().toGL();for(var cp=0;cp<shape._clipPlanes.length;cp++){var clip_plane=shape._clipPlanes[cp].plane;var clip_trafo=shape._clipPlanes[cp].trafo;sp['clipPlane'+cp+'_Plane']=clip_trafo.multMatrixPlane(clip_plane._vf.plane).toGL();sp['clipPlane'+cp+'_CappingStrength']=clip_plane._vf.cappingStrength;sp['clipPlane'+cp+'_CappingColor']=clip_plane._vf.cappingColor.toGL();}}
if(s_gl.imageGeometry!=0&&!x3dom.caps.MOBILE)
{sp.IG_bboxMin=s_geo.getMin().toGL();sp.IG_bboxMax=s_geo.getMax().toGL();sp.IG_implicitMeshSize=s_geo._vf.implicitMeshSize.toGL();var coordTex=x3dom.Utils.findTextureByName(s_gl.texture,"IG_coords0");if(coordTex){sp.IG_coordTextureWidth=coordTex.texture.width;sp.IG_coordTextureHeight=coordTex.texture.height;}
if(s_gl.imageGeometry==1){var indexTex=x3dom.Utils.findTextureByName(s_gl.texture,"IG_index");if(indexTex){sp.IG_indexTextureWidth=indexTex.texture.width;sp.IG_indexTextureHeight=indexTex.texture.height;}
gl.activeTexture(gl.TEXTURE0);gl.bindTexture(gl.TEXTURE_2D,indexTex.texture);gl.activeTexture(gl.TEXTURE1);gl.bindTexture(gl.TEXTURE_2D,coordTex.texture);}
else{gl.activeTexture(gl.TEXTURE0);gl.bindTexture(gl.TEXTURE_2D,coordTex.texture);}
gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_S,gl.CLAMP_TO_EDGE);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_T,gl.CLAMP_TO_EDGE);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MAG_FILTER,gl.NEAREST);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MIN_FILTER,gl.NEAREST);var texUnit=0;if(s_geo.getIndexTexture()){if(!sp.IG_indexTexture){sp.IG_indexTexture=texUnit++;}}
if(s_geo.getCoordinateTexture(0)){if(!sp.IG_coordinateTexture){sp.IG_coordinateTexture=texUnit++;}}}
else if((s_gl.binaryGeometry!=0||s_gl.externalGeometry!=0)&&s_geo._vf["idsPerVertex"]==true){var shader=s_app._shader;if(shader&&x3dom.isa(s_app._shader,x3dom.nodeTypes.CommonSurfaceShader)){if(shader.getMultiVisibilityMap()){sp.multiVisibilityMap=0;var visTex=x3dom.Utils.findTextureByName(s_gl.texture,"multiVisibilityMap");sp.multiVisibilityWidth=visTex.texture.width;sp.multiVisibilityHeight=visTex.texture.height;gl.activeTexture(gl.TEXTURE0);gl.bindTexture(gl.TEXTURE_2D,visTex.texture);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_S,gl.CLAMP_TO_EDGE);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_T,gl.CLAMP_TO_EDGE);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MAG_FILTER,gl.NEAREST);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MIN_FILTER,gl.NEAREST);}}}
if(shape.isSolid()){this.stateManager.enable(gl.CULL_FACE);if(shape.isCCW()){this.stateManager.frontFace(gl.CCW);}
else{this.stateManager.frontFace(gl.CW);}}
else{this.stateManager.disable(gl.CULL_FACE);}
var depthMode=s_app?s_app._cf.depthMode.node:null;if(depthMode)
{if(depthMode._vf.enableDepthTest)
{this.stateManager.enable(gl.DEPTH_TEST);this.stateManager.depthMask(!depthMode._vf.readOnly);this.stateManager.depthFunc(x3dom.Utils.depthFunc(gl,depthMode._vf.depthFunc));this.stateManager.depthRange(depthMode._vf.zNearRange,depthMode._vf.zFarRange);}
else
{this.stateManager.disable(gl.DEPTH_TEST);}}
else
{this.stateManager.enable(gl.DEPTH_TEST);this.stateManager.depthMask(true);this.stateManager.depthFunc(gl.LEQUAL);}
if(s_gl.popGeometry){var model_view=mat_view.mult(trafo);this.updatePopState(drawable,s_geo,sp,s_gl,scene,model_view,viewarea,this.x3dElem.runtime.fps);}
var q_n;if(s_gl.externalGeometry!=0)
{q_n=shape.meshes.length;}
else
{q_n=s_gl.positions.length;}
for(var q=0;q<q_n;q++){var q6=6*q;var v,v_n,offset;if(s_gl.externalGeometry!=0){var mesh=shape.meshes[q];mesh.bindVertexAttribPointerPosition(gl,sp,false);mesh.render(gl,null);}
else
if(!(sp.position!==undefined&&s_gl.buffers[q6+1]&&s_gl.indexes[q]))
continue;indicesReady=false;if(s_gl.externalGeometry==0){if(s_gl.buffers[q6]){gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER,s_gl.buffers[q6]);indicesReady=true;}
this.setVertexAttribPointerPosition(gl,shape,q6,q);if(sp.id!==undefined&&s_gl.buffers[q6+5]){gl.bindBuffer(gl.ARRAY_BUFFER,s_gl.buffers[q6+5]);if((s_gl.binaryGeometry!=0||s_gl.externalGeometry!=0)&&s_geo._vf["idsPerVertex"]==true){gl.vertexAttribPointer(sp.id,1,gl.FLOAT,false,4,0);gl.enableVertexAttribArray(sp.id);}
else{}}
if(indicesReady&&(s_gl.binaryGeometry>0||s_gl.popGeometry>0)){for(v=0,offset=0,v_n=s_geo._vf.vertexCount.length;v<v_n;v++){gl.drawElements(s_gl.primType[v],s_geo._vf.vertexCount[v],s_gl.indexType,x3dom.Utils.getByteAwareOffset(offset,s_gl.indexType,gl));offset+=s_geo._vf.vertexCount[v];}}
else if(s_gl.binaryGeometry<0||s_gl.popGeometry<0||s_gl.imageGeometry){for(v=0,offset=0,v_n=s_geo._vf.vertexCount.length;v<v_n;v++){gl.drawArrays(s_gl.primType[v],offset,s_geo._vf.vertexCount[v]);offset+=s_geo._vf.vertexCount[v];}}
else if(s_geo.hasIndexOffset()){var indOff=shape.tessellationProperties();for(v=0,v_n=indOff.length;v<v_n;v++){gl.drawElements(s_gl.primType,indOff[v].count,s_gl.indexType,indOff[v].offset*x3dom.Utils.getOffsetMultiplier(s_gl.indexType,gl));}}
else if(s_gl.indexes[q].length==0){gl.drawArrays(s_gl.primType,0,s_gl.positions[q].length/3);}
else{gl.drawElements(s_gl.primType,s_gl.indexes[q].length,s_gl.indexType,0);}}
gl.disableVertexAttribArray(sp.position);if(sp.texcoord!==undefined&&s_gl.buffers[q6+3]){gl.disableVertexAttribArray(sp.texcoord);}
if(sp.color!==undefined&&s_gl.buffers[q6+4]){gl.disableVertexAttribArray(sp.color);}
if(sp.id!==undefined&&s_gl.buffers[q6+5]){gl.disableVertexAttribArray(sp.id);}}
if(s_gl.imageGeometry!=0&&!x3dom.caps.MOBILE){gl.activeTexture(gl.TEXTURE0);gl.bindTexture(gl.TEXTURE_2D,null);if(s_gl.imageGeometry==1){gl.activeTexture(gl.TEXTURE1);gl.bindTexture(gl.TEXTURE_2D,null);}}}
if(x3dom.Utils.needLineWidth){this.stateManager.lineWidth(1);}
if(depthMode){this.stateManager.enable(gl.DEPTH_TEST);this.stateManager.depthMask(true);this.stateManager.depthFunc(gl.LEQUAL);this.stateManager.depthRange(0,1);}
gl.flush();this.stateManager.bindFramebuffer(gl.FRAMEBUFFER,null);};Context.prototype.renderPickingPass=function(gl,scene,mat_view,mat_scene,from,sceneSize,pickMode,lastX,lastY,width,height)
{var ps=scene._webgl.pickScale;var bufHeight=scene._webgl.fboPick.height;var x=lastX*ps;var y=(bufHeight-1)-lastY*ps;var indicesReady=false;this.stateManager.bindFramebuffer(gl.FRAMEBUFFER,scene._webgl.fboPick.fbo);this.stateManager.viewport(0,0,scene._webgl.fboPick.width,bufHeight);gl.clearColor(0.0,0.0,0.0,0.0);gl.clearDepth(1.0);gl.clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER_BIT);var viewarea=scene.drawableCollection.viewarea;var env=scene.getEnvironment();var n=scene.drawableCollection.length;if(env._vf.smallFeatureCulling&&env._lowPriorityThreshold<1&&viewarea.isMovingOrAnimating()){n=Math.floor(n*env._lowPriorityThreshold);if(!n&&scene.drawableCollection.length)
n=1;}
var bgCenter=x3dom.fields.SFVec3f.NullVector.toGL();var bgSize=x3dom.fields.SFVec3f.OneVector.toGL();this.stateManager.depthFunc(gl.LEQUAL);this.stateManager.enable(gl.DEPTH_TEST);this.stateManager.enable(gl.CULL_FACE);this.stateManager.disable(gl.BLEND);if(x3dom.Utils.needLineWidth){this.stateManager.lineWidth(2);}
for(var i=0;i<n;i++)
{var drawable=scene.drawableCollection.get(i);var trafo=drawable.transform;var shape=drawable.shape;var s_gl=shape._webgl;if(!s_gl||shape._objectID<1||!shape._vf.isPickable){continue;}
var s_geo=shape._cf.geometry.node;var s_app=shape._cf.appearance.node;var s_msh=s_geo._mesh;var properties=shape.getShaderProperties(viewarea);var sp=this.cache.getShaderByProperties(gl,shape,properties,pickMode);if(!sp){return;}
this.stateManager.useProgram(sp);sp.modelMatrix=trafo.toGL();sp.modelViewProjectionMatrix=mat_scene.mult(trafo).toGL();sp.lowBit=(shape._objectID&255)/255.0;sp.highBit=(shape._objectID>>>8)/255.0;sp.from=from.toGL();sp.sceneSize=sceneSize;if((s_gl.binaryGeometry!=0||s_gl.externalGeometry!=0)&&s_geo._vf["idsPerVertex"]==true){sp.shadowIDs=(shape._vf.idOffset+x3dom.nodeTypes.Shape.objectID+2);}
if(s_gl.coordType!=gl.FLOAT){if(!s_gl.popGeometry&&(x3dom.Utils.isUnsignedType(s_geo._vf.coordType))){sp.bgCenter=s_geo.getMin().toGL();}
else{sp.bgCenter=s_geo._vf.position.toGL();}
sp.bgSize=s_geo._vf.size.toGL();sp.bgPrecisionMax=s_geo.getPrecisionMax('coordType');}
if(pickMode==1&&s_gl.colorType!=gl.FLOAT){sp.bgPrecisionColMax=s_geo.getPrecisionMax('colorType');}
if(pickMode==2&&s_gl.texCoordType!=gl.FLOAT){sp.bgPrecisionTexMax=s_geo.getPrecisionMax('texCoordType');}
if(shape._clipPlanes){sp.modelViewMatrix=mat_view.mult(trafo).toGL();sp.viewMatrixInverse=mat_view.inverse().toGL();for(var cp=0;cp<shape._clipPlanes.length;cp++){var clip_plane=shape._clipPlanes[cp].plane;var clip_trafo=shape._clipPlanes[cp].trafo;sp['clipPlane'+cp+'_Plane']=clip_trafo.multMatrixPlane(clip_plane._vf.plane).toGL();sp['clipPlane'+cp+'_CappingStrength']=clip_plane._vf.cappingStrength;sp['clipPlane'+cp+'_CappingColor']=clip_plane._vf.cappingColor.toGL();}}
if(s_gl.imageGeometry!=0&&!x3dom.caps.MOBILE)
{sp.IG_bboxMin=s_geo.getMin().toGL();sp.IG_bboxMax=s_geo.getMax().toGL();sp.IG_implicitMeshSize=s_geo._vf.implicitMeshSize.toGL();var coordTex=x3dom.Utils.findTextureByName(s_gl.texture,"IG_coords0");if(coordTex){sp.IG_coordTextureWidth=coordTex.texture.width;sp.IG_coordTextureHeight=coordTex.texture.height;}
if(s_gl.imageGeometry==1){var indexTex=x3dom.Utils.findTextureByName(s_gl.texture,"IG_index");if(indexTex){sp.IG_indexTextureWidth=indexTex.texture.width;sp.IG_indexTextureHeight=indexTex.texture.height;}
gl.activeTexture(gl.TEXTURE0);gl.bindTexture(gl.TEXTURE_2D,indexTex.texture);gl.activeTexture(gl.TEXTURE1);gl.bindTexture(gl.TEXTURE_2D,coordTex.texture);}
else{gl.activeTexture(gl.TEXTURE0);gl.bindTexture(gl.TEXTURE_2D,coordTex.texture);}
gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_S,gl.CLAMP_TO_EDGE);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_T,gl.CLAMP_TO_EDGE);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MAG_FILTER,gl.NEAREST);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MIN_FILTER,gl.NEAREST);var texUnit=0;if(s_geo.getIndexTexture()){if(!sp.IG_indexTexture){sp.IG_indexTexture=texUnit++;}}
if(s_geo.getCoordinateTexture(0)){if(!sp.IG_coordinateTexture){sp.IG_coordinateTexture=texUnit++;}}}
else if((s_gl.binaryGeometry!=0||s_gl.externalGeometry!=0)&&s_geo._vf["idsPerVertex"]==true){var shader=s_app._shader;if(shader&&x3dom.isa(s_app._shader,x3dom.nodeTypes.CommonSurfaceShader)){if(shader.getMultiVisibilityMap()){sp.multiVisibilityMap=0;var visTex=x3dom.Utils.findTextureByName(s_gl.texture,"multiVisibilityMap");sp.multiVisibilityWidth=visTex.texture.width;sp.multiVisibilityHeight=visTex.texture.height;gl.activeTexture(gl.TEXTURE0);gl.bindTexture(gl.TEXTURE_2D,visTex.texture);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_S,gl.CLAMP_TO_EDGE);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_T,gl.CLAMP_TO_EDGE);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MAG_FILTER,gl.NEAREST);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MIN_FILTER,gl.NEAREST);}}}
if(shape.isSolid()){this.stateManager.enable(gl.CULL_FACE);if(shape.isCCW()){this.stateManager.frontFace(gl.CCW);}
else{this.stateManager.frontFace(gl.CW);}}
else{this.stateManager.disable(gl.CULL_FACE);}
var depthMode=s_app?s_app._cf.depthMode.node:null;if(depthMode)
{if(depthMode._vf.enableDepthTest)
{this.stateManager.enable(gl.DEPTH_TEST);this.stateManager.depthMask(!depthMode._vf.readOnly);this.stateManager.depthFunc(x3dom.Utils.depthFunc(gl,depthMode._vf.depthFunc));this.stateManager.depthRange(depthMode._vf.zNearRange,depthMode._vf.zFarRange);}
else
{this.stateManager.disable(gl.DEPTH_TEST);}}
else
{this.stateManager.enable(gl.DEPTH_TEST);this.stateManager.depthMask(true);this.stateManager.depthFunc(gl.LEQUAL);}
if(s_gl.popGeometry){var model_view=mat_view.mult(trafo);this.updatePopState(drawable,s_geo,sp,s_gl,scene,model_view,viewarea,this.x3dElem.runtime.fps);}
var q_n;if(s_gl.externalGeometry!=0)
{q_n=shape.meshes.length;}
else
{q_n=s_gl.positions.length;}
for(var q=0;q<q_n;q++){var q6=6*q;var v,v_n,offset;if(s_gl.externalGeometry!=0){var mesh=shape.meshes[q];mesh.bindVertexAttribPointerPosition(gl,sp,false);mesh.render(gl,null);}
else
if(!(sp.position!==undefined&&s_gl.buffers[q6+1]&&s_gl.indexes[q]))
continue;indicesReady=false;if(s_gl.externalGeometry==0){if(s_gl.buffers[q6]){gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER,s_gl.buffers[q6]);indicesReady=true;}
this.setVertexAttribPointerPosition(gl,shape,q6,q);if(pickMode==1){this.setVertexAttribPointerColor(gl,shape,q6,q);}
if(pickMode==2&&sp.texcoord!==undefined&&s_gl.buffers[q6+3]){this.setVertexAttribPointerTexCoord(gl,shape,q6,q);}
if(sp.id!==undefined&&s_gl.buffers[q6+5]){gl.bindBuffer(gl.ARRAY_BUFFER,s_gl.buffers[q6+5]);if((s_gl.binaryGeometry!=0||s_gl.externalGeometry!=0)&&s_geo._vf["idsPerVertex"]==true){gl.vertexAttribPointer(sp.id,1,gl.FLOAT,false,4,0);gl.enableVertexAttribArray(sp.id);}
else{}}
if(indicesReady&&(s_gl.binaryGeometry>0||s_gl.popGeometry>0)){for(v=0,offset=0,v_n=s_geo._vf.vertexCount.length;v<v_n;v++){gl.drawElements(s_gl.primType[v],s_geo._vf.vertexCount[v],s_gl.indexType,x3dom.Utils.getByteAwareOffset(offset,s_gl.indexType,gl));offset+=s_geo._vf.vertexCount[v];}}
else if(s_gl.binaryGeometry<0||s_gl.popGeometry<0||s_gl.imageGeometry){for(v=0,offset=0,v_n=s_geo._vf.vertexCount.length;v<v_n;v++){gl.drawArrays(s_gl.primType[v],offset,s_geo._vf.vertexCount[v]);offset+=s_geo._vf.vertexCount[v];}}
else if(s_geo.hasIndexOffset()){var indOff=shape.tessellationProperties();for(v=0,v_n=indOff.length;v<v_n;v++){gl.drawElements(s_gl.primType,indOff[v].count,s_gl.indexType,indOff[v].offset*x3dom.Utils.getOffsetMultiplier(s_gl.indexType,gl));}}
else if(s_gl.indexes[q].length==0){gl.drawArrays(s_gl.primType,0,s_gl.positions[q].length/3);}
else{gl.drawElements(s_gl.primType,s_gl.indexes[q].length,s_gl.indexType,0);}}
gl.disableVertexAttribArray(sp.position);if(sp.texcoord!==undefined&&s_gl.buffers[q6+3]){gl.disableVertexAttribArray(sp.texcoord);}
if(sp.color!==undefined&&s_gl.buffers[q6+4]){gl.disableVertexAttribArray(sp.color);}
if(sp.id!==undefined&&s_gl.buffers[q6+5]){gl.disableVertexAttribArray(sp.id);}}
if(s_gl.imageGeometry!=0&&!x3dom.caps.MOBILE){gl.activeTexture(gl.TEXTURE0);gl.bindTexture(gl.TEXTURE_2D,null);if(s_gl.imageGeometry==1){gl.activeTexture(gl.TEXTURE1);gl.bindTexture(gl.TEXTURE_2D,null);}}}
if(x3dom.Utils.needLineWidth){this.stateManager.lineWidth(1);}
if(depthMode){this.stateManager.enable(gl.DEPTH_TEST);this.stateManager.depthMask(true);this.stateManager.depthFunc(gl.LEQUAL);this.stateManager.depthRange(0,1);}
gl.flush();try{var data=new Uint8Array(4*width*height);gl.readPixels(x,y,width,height,gl.RGBA,gl.UNSIGNED_BYTE,data);scene._webgl.fboPick.pixelData=data;}
catch(se){scene._webgl.fboPick.pixelData=[];x3dom.debug.logException(se+" (cannot pick)");}
this.stateManager.bindFramebuffer(gl.FRAMEBUFFER,null);};Context.prototype.renderShape=function(drawable,viewarea,slights,numLights,mat_view,mat_scene,mat_light,mat_proj,gl)
{var indicesReady=false;var shape=drawable.shape;var transform=drawable.transform;if(!shape||!shape._webgl||!transform){x3dom.debug.logError("[Context|RenderShape] No valid Shape!");return;}
var s_gl=shape._webgl;var sp=s_gl.shader;if(!sp){x3dom.debug.logError("[Context|RenderShape] No Shader is set!");return;}
var changed=this.stateManager.useProgram(sp);var s_app=shape._cf.appearance.node;var s_geo=shape._cf.geometry.node;var s_msh=s_geo._mesh;var scene=viewarea._scene;var tex=null;if(s_gl.coordType!=gl.FLOAT){if(!s_gl.popGeometry&&(x3dom.Utils.isUnsignedType(s_geo._vf.coordType))){sp.bgCenter=s_geo.getMin().toGL();}
else{sp.bgCenter=s_geo._vf.position.toGL();}
sp.bgSize=s_geo._vf.size.toGL();sp.bgPrecisionMax=s_geo.getPrecisionMax('coordType');}
else{sp.bgCenter=[0,0,0];sp.bgSize=[1,1,1];sp.bgPrecisionMax=1;}
if(s_gl.colorType!=gl.FLOAT){sp.bgPrecisionColMax=s_geo.getPrecisionMax('colorType');}
else{sp.bgPrecisionColMax=1;}
if(s_gl.texCoordType!=gl.FLOAT){sp.bgPrecisionTexMax=s_geo.getPrecisionMax('texCoordType');}
else{sp.bgPrecisionTexMax=1;}
if(s_gl.normalType!=gl.FLOAT){sp.bgPrecisionNorMax=s_geo.getPrecisionMax('normalType');}
else{sp.bgPrecisionNorMax=1;}
if(s_gl.imageGeometry!=0){sp.IG_bboxMin=s_geo.getMin().toGL();sp.IG_bboxMax=s_geo.getMax().toGL();sp.IG_implicitMeshSize=s_geo._vf.implicitMeshSize.toGL();tex=x3dom.Utils.findTextureByName(s_gl.texture,"IG_coords0");if(tex){sp.IG_coordTextureWidth=tex.texture.width;sp.IG_coordTextureHeight=tex.texture.height;}
if(s_gl.imageGeometry==1){tex=x3dom.Utils.findTextureByName(s_gl.texture,"IG_index");if(tex){sp.IG_indexTextureWidth=tex.texture.width;sp.IG_indexTextureHeight=tex.texture.height;}}
tex=null;}
var fog=scene.getFog();if(fog&&changed){sp.fogColor=fog._vf.color.toGL();sp.fogRange=fog._vf.visibilityRange;sp.fogType=(fog._vf.fogType=="LINEAR")?0.0:1.0;}
var mat=s_app?s_app._cf.material.node:null;var shader=s_app?s_app._shader:null;var twoSidedMat=false;var isUserDefinedShader=shader&&x3dom.isa(shader,x3dom.nodeTypes.ComposedShader);if(s_gl.csshader){sp.diffuseColor=shader._vf.diffuseFactor.toGL();sp.specularColor=shader._vf.specularFactor.toGL();sp.emissiveColor=shader._vf.emissiveFactor.toGL();sp.shininess=shader._vf.shininessFactor;sp.ambientIntensity=(shader._vf.ambientFactor.x+
shader._vf.ambientFactor.y+
shader._vf.ambientFactor.z)/3;sp.transparency=1.0-shader._vf.alphaFactor;sp.environmentFactor=shader._vf.environmentFactor.x;if(shader.getDisplacementMap()){tex=x3dom.Utils.findTextureByName(s_gl.texture,"displacementMap");sp.displacementWidth=tex.texture.width;sp.displacementHeight=tex.texture.height;sp.displacementFactor=shader._vf.displacementFactor;sp.displacementAxis=(shader._vf.displacementAxis=="x")?0.0:(shader._vf.displacementAxis=="y")?1.0:2.0;}
else if(shader.getDiffuseDisplacementMap()){tex=x3dom.Utils.findTextureByName(s_gl.texture,"diffuseDisplacementMap");sp.displacementWidth=tex.texture.width;sp.displacementHeight=tex.texture.height;sp.displacementFactor=shader._vf.displacementFactor;sp.displacementAxis=(shader._vf.displacementAxis=="x")?0.0:(shader._vf.displacementAxis=="y")?1.0:2.0;}
if(shader.getMultiDiffuseAlphaMap()){tex=x3dom.Utils.findTextureByName(s_gl.texture,"multiDiffuseAlphaMap");sp.multiDiffuseAlphaWidth=tex.texture.width;sp.multiDiffuseAlphaHeight=tex.texture.height;}
if(shader.getMultiEmissiveAmbientMap()){tex=x3dom.Utils.findTextureByName(s_gl.texture,"multiEmissiveAmbientMap");sp.multiEmissiveAmbientWidth=tex.texture.width;sp.multiEmissiveAmbientHeight=tex.texture.height;}
if(shader.getMultiSpecularShininessMap()){tex=x3dom.Utils.findTextureByName(s_gl.texture,"multiSpecularShininessMap");sp.multiSpecularShininessWidth=tex.texture.width;sp.multiSpecularShininessHeight=tex.texture.height;}
if(shader.getMultiVisibilityMap()){tex=x3dom.Utils.findTextureByName(s_gl.texture,"multiVisibilityMap");sp.multiVisibilityWidth=tex.texture.width;sp.multiVisibilityHeight=tex.texture.height;}}
else if(mat){sp.diffuseColor=mat._vf.diffuseColor.toGL();sp.specularColor=mat._vf.specularColor.toGL();sp.emissiveColor=mat._vf.emissiveColor.toGL();sp.shininess=mat._vf.shininess;sp.ambientIntensity=mat._vf.ambientIntensity;sp.transparency=mat._vf.transparency;sp.environmentFactor=0.0;if(x3dom.isa(mat,x3dom.nodeTypes.TwoSidedMaterial)){twoSidedMat=true;sp.backDiffuseColor=mat._vf.backDiffuseColor.toGL();sp.backSpecularColor=mat._vf.backSpecularColor.toGL();sp.backEmissiveColor=mat._vf.backEmissiveColor.toGL();sp.backShininess=mat._vf.backShininess;sp.backAmbientIntensity=mat._vf.backAmbientIntensity;sp.backTransparency=mat._vf.backTransparency;}}
else{sp.diffuseColor=[1.0,1.0,1.0];sp.specularColor=[0.0,0.0,0.0];sp.emissiveColor=[0.0,0.0,0.0];sp.shininess=0.0;sp.ambientIntensity=1.0;sp.transparency=0.0;}
if(shader){if(isUserDefinedShader){for(var fName in shader._vf){if(shader._vf.hasOwnProperty(fName)&&fName!=='language'){var field=shader._vf[fName];if(field!==undefined&&field!==null){if(field.toGL){sp[fName]=field.toGL();}
else{sp[fName]=field;}}}}}
else if(x3dom.isa(shader,x3dom.nodeTypes.CommonSurfaceShader)){s_gl.csshader=shader;}}
for(var p=0;p<numLights&&changed;p++){var light_transform=mat_view.mult(slights[p].getCurrentTransform());if(x3dom.isa(slights[p],x3dom.nodeTypes.DirectionalLight)){sp['light'+p+'_Type']=0.0;sp['light'+p+'_On']=(slights[p]._vf.on)?1.0:0.0;sp['light'+p+'_Color']=slights[p]._vf.color.toGL();sp['light'+p+'_Intensity']=slights[p]._vf.intensity;sp['light'+p+'_AmbientIntensity']=slights[p]._vf.ambientIntensity;sp['light'+p+'_Direction']=light_transform.multMatrixVec(slights[p]._vf.direction).toGL();sp['light'+p+'_Attenuation']=[1.0,1.0,1.0];sp['light'+p+'_Location']=[1.0,1.0,1.0];sp['light'+p+'_Radius']=0.0;sp['light'+p+'_BeamWidth']=0.0;sp['light'+p+'_CutOffAngle']=0.0;sp['light'+p+'_ShadowIntensity']=slights[p]._vf.shadowIntensity;}
else if(x3dom.isa(slights[p],x3dom.nodeTypes.PointLight)){sp['light'+p+'_Type']=1.0;sp['light'+p+'_On']=(slights[p]._vf.on)?1.0:0.0;sp['light'+p+'_Color']=slights[p]._vf.color.toGL();sp['light'+p+'_Intensity']=slights[p]._vf.intensity;sp['light'+p+'_AmbientIntensity']=slights[p]._vf.ambientIntensity;sp['light'+p+'_Direction']=[1.0,1.0,1.0];sp['light'+p+'_Attenuation']=slights[p]._vf.attenuation.toGL();sp['light'+p+'_Location']=light_transform.multMatrixPnt(slights[p]._vf.location).toGL();sp['light'+p+'_Radius']=slights[p]._vf.radius;sp['light'+p+'_BeamWidth']=0.0;sp['light'+p+'_CutOffAngle']=0.0;sp['light'+p+'_ShadowIntensity']=slights[p]._vf.shadowIntensity;}
else if(x3dom.isa(slights[p],x3dom.nodeTypes.SpotLight)){sp['light'+p+'_Type']=2.0;sp['light'+p+'_On']=(slights[p]._vf.on)?1.0:0.0;sp['light'+p+'_Color']=slights[p]._vf.color.toGL();sp['light'+p+'_Intensity']=slights[p]._vf.intensity;sp['light'+p+'_AmbientIntensity']=slights[p]._vf.ambientIntensity;sp['light'+p+'_Direction']=light_transform.multMatrixVec(slights[p]._vf.direction).toGL();sp['light'+p+'_Attenuation']=slights[p]._vf.attenuation.toGL();sp['light'+p+'_Location']=light_transform.multMatrixPnt(slights[p]._vf.location).toGL();sp['light'+p+'_Radius']=slights[p]._vf.radius;sp['light'+p+'_BeamWidth']=slights[p]._vf.beamWidth;sp['light'+p+'_CutOffAngle']=slights[p]._vf.cutOffAngle;sp['light'+p+'_ShadowIntensity']=slights[p]._vf.shadowIntensity;}}
var nav=scene.getNavigationInfo();if(nav._vf.headlight&&changed){numLights=(numLights)?numLights:0;sp['light'+numLights+'_Type']=0.0;sp['light'+numLights+'_On']=1.0;sp['light'+numLights+'_Color']=[1.0,1.0,1.0];sp['light'+numLights+'_Intensity']=1.0;sp['light'+numLights+'_AmbientIntensity']=0.0;sp['light'+numLights+'_Direction']=[0.0,0.0,-1.0];sp['light'+numLights+'_Attenuation']=[1.0,1.0,1.0];sp['light'+numLights+'_Location']=[1.0,1.0,1.0];sp['light'+numLights+'_Radius']=0.0;sp['light'+numLights+'_BeamWidth']=0.0;sp['light'+numLights+'_CutOffAngle']=0.0;sp['light'+numLights+'_ShadowIntensity']=0.0;}
if(shape._clipPlanes){for(var cp=0;cp<shape._clipPlanes.length;cp++){var clip_plane=shape._clipPlanes[cp].plane;var clip_trafo=shape._clipPlanes[cp].trafo;sp['clipPlane'+cp+'_Plane']=clip_trafo.multMatrixPlane(clip_plane._vf.plane).toGL();sp['clipPlane'+cp+'_CappingStrength']=clip_plane._vf.cappingStrength;sp['clipPlane'+cp+'_CappingColor']=clip_plane._vf.cappingColor.toGL();}}
var depthMode=s_app?s_app._cf.depthMode.node:null;if(depthMode)
{if(depthMode._vf.enableDepthTest)
{this.stateManager.enable(gl.DEPTH_TEST);this.stateManager.depthMask(!depthMode._vf.readOnly);this.stateManager.depthFunc(x3dom.Utils.depthFunc(gl,depthMode._vf.depthFunc));this.stateManager.depthRange(depthMode._vf.zNearRange,depthMode._vf.zFarRange);}
else
{this.stateManager.disable(gl.DEPTH_TEST);}}
else
{this.stateManager.enable(gl.DEPTH_TEST);this.stateManager.depthMask(true);this.stateManager.depthFunc(gl.LEQUAL);}
var blendMode=s_app?s_app._cf.blendMode.node:null;if(blendMode)
{var srcFactor=x3dom.Utils.blendFunc(gl,blendMode._vf.srcFactor);var destFactor=x3dom.Utils.blendFunc(gl,blendMode._vf.destFactor);if(srcFactor&&destFactor)
{this.stateManager.enable(gl.BLEND);this.stateManager.blendFuncSeparate(srcFactor,destFactor,gl.ONE,gl.ONE);this.stateManager.blendColor(blendMode._vf.color.r,blendMode._vf.color.g,blendMode._vf.color.b,1.0-blendMode._vf.colorTransparency);this.stateManager.blendEquation(x3dom.Utils.blendEquation(gl,blendMode._vf.equation));}
else
{this.stateManager.disable(gl.BLEND);}}
else
{this.stateManager.enable(gl.BLEND);this.stateManager.blendFuncSeparate(gl.SRC_ALPHA,gl.ONE_MINUS_SRC_ALPHA,gl.ONE,gl.ONE);}
var colorMaskMode=s_app?s_app._cf.colorMaskMode.node:null;if(colorMaskMode)
{this.stateManager.colorMask(colorMaskMode._vf.maskR,colorMaskMode._vf.maskG,colorMaskMode._vf.maskB,colorMaskMode._vf.maskA);}
else
{this.stateManager.colorMask(true,true,true,true);}
var lineProperties=s_app?s_app._cf.lineProperties.node:null;if(lineProperties)
{this.stateManager.lineWidth(lineProperties._vf.linewidthScaleFactor);}
else if(x3dom.Utils.needLineWidth)
{this.stateManager.lineWidth(1);}
if(shape.isSolid()&&!twoSidedMat){this.stateManager.enable(gl.CULL_FACE);if(shape.isCCW()){this.stateManager.frontFace(gl.CCW);}
else{this.stateManager.frontFace(gl.CW);}}
else{this.stateManager.disable(gl.CULL_FACE);}
var model_view=mat_view.mult(transform);var model_view_inv=model_view.inverse();sp.isOrthoView=(mat_proj._33==1)?1.0:0.0;sp.modelViewMatrix=model_view.toGL();sp.viewMatrix=mat_view.toGL();sp.normalMatrix=model_view_inv.transpose().toGL();sp.modelViewMatrixInverse=model_view_inv.toGL();sp.modelViewProjectionMatrix=mat_scene.mult(transform).toGL();if(isUserDefinedShader||shape._clipPlanes&&shape._clipPlanes.length)
{sp.viewMatrixInverse=mat_view.inverse().toGL();}
if(isUserDefinedShader||s_gl.externalGeometry!=0){sp.model=transform.toGL();sp.projectionMatrix=mat_proj.toGL();sp.worldMatrix=transform.toGL();sp.worldInverseTranspose=transform.inverse().transpose().toGL();}
if(s_gl.popGeometry){this.updatePopState(drawable,s_geo,sp,s_gl,scene,model_view,viewarea,this.x3dElem.runtime.fps);}
for(var cnt=0,cnt_n=s_gl.texture.length;cnt<cnt_n;cnt++){tex=s_gl.texture[cnt];gl.activeTexture(gl.TEXTURE0+cnt);gl.bindTexture(tex.type,tex.texture);gl.texParameteri(tex.type,gl.TEXTURE_WRAP_S,tex.wrapS);gl.texParameteri(tex.type,gl.TEXTURE_WRAP_T,tex.wrapT);gl.texParameteri(tex.type,gl.TEXTURE_MAG_FILTER,tex.magFilter);gl.texParameteri(tex.type,gl.TEXTURE_MIN_FILTER,tex.minFilter);if(!shader||!isUserDefinedShader){if(!sp[tex.samplerName])
sp[tex.samplerName]=cnt;}}
if(s_app&&s_app._cf.textureTransform.node){var texTrafo=s_app.texTransformMatrix();sp.texTrafoMatrix=texTrafo.toGL();}
var attrib=null;var df,df_n=s_gl.dynamicFields.length;for(df=0;df<df_n;df++){attrib=s_gl.dynamicFields[df];if(sp[attrib.name]!==undefined){gl.bindBuffer(gl.ARRAY_BUFFER,attrib.buf);gl.vertexAttribPointer(sp[attrib.name],attrib.numComponents,gl.FLOAT,false,0,0);gl.enableVertexAttribArray(sp[attrib.name]);}}
var v,v_n,offset,q_n;var isParticleSet=false;if(x3dom.isa(s_geo,x3dom.nodeTypes.ParticleSet)){isParticleSet=true;}
if(s_gl.externalGeometry!=0){q_n=shape.meshes.length;}
else
{q_n=s_gl.positions.length;}
for(var q=0;q<q_n;q++){var q6=6*q;if(s_gl.externalGeometry!=0){var mesh=shape.meshes[q];var exGeomShaderProgram=sp;if(mesh.material!=null){if(mesh.material.program!=null){exGeomShaderProgram=mesh.material.program;}
if(mesh.material.setShader!=null)
mesh.material.setShader(gl,this.cache,shape,shape.getShaderProperties(viewarea));mesh.material.bind(gl,sp,this.cache,shape.getShaderProperties(viewarea));}
mesh.bindVertexAttribPointer(gl,exGeomShaderProgram);var renderMode=viewarea.getRenderMode();var polyMode=null;if(renderMode>0)
polyMode=(renderMode==1)?gl.POINTS:gl.LINES;mesh.render(gl,polyMode);}
else
if(!(sp.position!==undefined&&s_gl.buffers[q6+1]&&s_gl.indexes[q]))
continue;indicesReady=false;if(s_gl.externalGeometry==0){if(!(sp.position!==undefined&&s_gl.buffers[q6+1]&&(s_gl.indexes[q])))
continue;if(s_gl.buffers[q6]){if(isParticleSet&&s_geo.drawOrder()!="any"){var indexArray,zPos=[];var pnts=s_geo._cf.coord.node.getPoints();var pn=(pnts.length==s_gl.indexes[q].length)?s_gl.indexes[q].length:0;for(var i=0;i<pn;i++){var center=model_view.multMatrixPnt(pnts[i]);zPos.push([i,center.z]);}
if(s_geo.drawOrder()=="backtofront")
zPos.sort(function(a,b){return a[1]-b[1];});else
zPos.sort(function(b,a){return a[1]-b[1];});for(i=0;i<pn;i++){shape._webgl.indexes[q][i]=zPos[i][0];}
if(x3dom.caps.INDEX_UINT&&(pn>65535)){indexArray=new Uint32Array(shape._webgl.indexes[q]);shape._webgl.indexType=gl.UNSIGNED_INT;}
else{indexArray=new Uint16Array(shape._webgl.indexes[q]);shape._webgl.indexType=gl.UNSIGNED_SHORT;}
gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER,s_gl.buffers[q6]);gl.bufferData(gl.ELEMENT_ARRAY_BUFFER,indexArray,gl.DYNAMIC_DRAW);indexArray=null;}
gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER,s_gl.buffers[q6]);indicesReady=true;}
this.setVertexAttribPointerPosition(gl,shape,q6,q);this.setVertexAttribPointerNormal(gl,shape,q6,q);this.setVertexAttribPointerTexCoord(gl,shape,q6,q);this.setVertexAttribPointerColor(gl,shape,q6,q);if((sp.id!==undefined||sp.particleSize!==undefined)&&shape._webgl.buffers[q6+5]){gl.bindBuffer(gl.ARRAY_BUFFER,shape._webgl.buffers[q6+5]);if((s_gl.binaryGeometry!=0||s_gl.externalGeometry!=0)&&s_geo._vf["idsPerVertex"]==true){gl.vertexAttribPointer(sp.id,1,gl.FLOAT,false,4,0);gl.enableVertexAttribArray(sp.id);}
else if(isParticleSet){gl.vertexAttribPointer(sp.particleSize,3,gl.FLOAT,false,0,0);gl.enableVertexAttribArray(sp.particleSize);}}
if(s_gl.popGeometry!=0&&s_gl.buffers[q6+5]){gl.bindBuffer(gl.ARRAY_BUFFER,s_gl.buffers[q6+5]);gl.vertexAttribPointer(sp.PG_vertexID,1,gl.FLOAT,false,4,0);gl.enableVertexAttribArray(sp.PG_vertexID);}
var indOff,renderMode=viewarea.getRenderMode();if(renderMode>0){var polyMode=(renderMode==1)?gl.POINTS:gl.LINES;if(indicesReady&&(s_gl.binaryGeometry>0||s_gl.popGeometry>0)){for(v=0,offset=0,v_n=s_geo._vf.vertexCount.length;v<v_n;v++){gl.drawElements(polyMode,s_geo._vf.vertexCount[v],s_gl.indexType,x3dom.Utils.getByteAwareOffset(offset,s_gl.indexType,gl));offset+=s_geo._vf.vertexCount[v];}}
else if(s_gl.binaryGeometry<0||s_gl.popGeometry<0||s_gl.imageGeometry){for(v=0,offset=0,v_n=s_geo._vf.vertexCount.length;v<v_n;v++){gl.drawArrays(polyMode,offset,s_geo._vf.vertexCount[v]);offset+=s_geo._vf.vertexCount[v];}}
else if(s_geo.hasIndexOffset()){indOff=shape.tessellationProperties();for(v=0,v_n=indOff.length;v<v_n;v++){gl.drawElements(polyMode,indOff[v].count,s_gl.indexType,indOff[v].offset*x3dom.Utils.getOffsetMultiplier(s_gl.indexType,gl));}}
else if(s_gl.indexes[q].length==0){gl.drawArrays(polyMode,0,s_gl.positions[q].length/3);}
else{gl.drawElements(polyMode,s_gl.indexes[q].length,s_gl.indexType,0);}}
else{if(indicesReady&&(s_gl.binaryGeometry>0||s_gl.popGeometry>0)){for(v=0,offset=0,v_n=s_geo._vf.vertexCount.length;v<v_n;v++){gl.drawElements(s_gl.primType[v],s_geo._vf.vertexCount[v],s_gl.indexType,x3dom.Utils.getByteAwareOffset(offset,s_gl.indexType,gl));offset+=s_geo._vf.vertexCount[v];}}
else if(s_gl.binaryGeometry<0||s_gl.popGeometry<0||s_gl.imageGeometry){for(v=0,offset=0,v_n=s_geo._vf.vertexCount.length;v<v_n;v++){gl.drawArrays(s_gl.primType[v],offset,s_geo._vf.vertexCount[v]);offset+=s_geo._vf.vertexCount[v];}}
else if(s_geo.hasIndexOffset()){indOff=shape.tessellationProperties();for(v=0,v_n=indOff.length;v<v_n;v++){gl.drawElements(s_gl.primType,indOff[v].count,s_gl.indexType,indOff[v].offset*x3dom.Utils.getOffsetMultiplier(s_gl.indexType,gl));}}
else if(s_gl.indexes[q].length==0){gl.drawArrays(s_gl.primType,0,s_gl.positions[q].length/3);}
else{gl.drawElements(s_gl.primType,s_gl.indexes[q].length,s_gl.indexType,0);}}}
gl.disableVertexAttribArray(sp.position);if(sp.normal!==undefined){gl.disableVertexAttribArray(sp.normal);}
if(sp.texcoord!==undefined){gl.disableVertexAttribArray(sp.texcoord);}
if(sp.color!==undefined){gl.disableVertexAttribArray(sp.color);}
if(s_gl.buffers[q6+5]){if(sp.id!==undefined)
gl.disableVertexAttribArray(sp.id);else if(sp.particleSize!==undefined)
gl.disableVertexAttribArray(sp.particleSize);}
if(s_gl.popGeometry!=0&&sp.PG_vertexID!==undefined){gl.disableVertexAttribArray(sp.PG_vertexID);}}
for(df=0;df<df_n;df++){attrib=s_gl.dynamicFields[df];if(sp[attrib.name]!==undefined){gl.disableVertexAttribArray(sp[attrib.name]);}}
if(s_gl.imageGeometry){v_n=s_geo._vf.vertexCount.length;this.numDrawCalls+=v_n;for(v=0;v<v_n;v++){if(s_gl.primType[v]==gl.TRIANGLE_STRIP)
this.numFaces+=(s_geo._vf.vertexCount[v]-2);else
this.numFaces+=(s_geo._vf.vertexCount[v]/3);this.numCoords+=s_geo._vf.vertexCount[v];}}
else{this.numCoords+=s_msh._numCoords;this.numFaces+=s_msh._numFaces;if(s_gl.binaryGeometry||s_gl.popGeometry){this.numDrawCalls+=s_geo._vf.vertexCount.length;}
else if(s_geo.hasIndexOffset()){this.numDrawCalls+=shape.tessellationProperties().length;}
else{this.numDrawCalls+=q_n;}}
if(depthMode){this.stateManager.enable(gl.DEPTH_TEST);this.stateManager.depthMask(true);this.stateManager.depthFunc(gl.LEQUAL);this.stateManager.depthRange(0,1);}
if(blendMode){this.stateManager.enable(gl.BLEND);this.stateManager.blendFuncSeparate(gl.SRC_ALPHA,gl.ONE_MINUS_SRC_ALPHA,gl.ONE,gl.ONE);this.stateManager.blendColor(1,1,1,1);this.stateManager.blendEquation(gl.FUNC_ADD);}
if(colorMaskMode){this.stateManager.colorMask(true,true,true,true);}
if(lineProperties){this.stateManager.lineWidth(1);}
var s_gl_tex=s_gl.texture;cnt_n=s_gl_tex?s_gl_tex.length:0;for(cnt=0;cnt<cnt_n;cnt++){if(!s_gl_tex[cnt])
continue;if(s_app&&s_app._cf.texture.node){tex=s_app._cf.texture.node.getTexture(cnt);gl.activeTexture(gl.TEXTURE0+cnt);if(x3dom.isa(tex,x3dom.nodeTypes.X3DEnvironmentTextureNode)){gl.bindTexture(gl.TEXTURE_CUBE_MAP,null);}
else{gl.bindTexture(gl.TEXTURE_2D,null);}}}};Context.prototype.updatePopState=function(drawable,popGeo,sp,s_gl,scene,model_view,viewarea,currFps)
{var tol=x3dom.nodeTypes.PopGeometry.ErrorToleranceFactor*popGeo._vf.precisionFactor;if(currFps<=1||viewarea.isMovingOrAnimating()){tol*=x3dom.nodeTypes.PopGeometry.PrecisionFactorOnMove;}
var currentLOD=16;if(tol>0){var viewpoint=scene.getViewpoint();var imgPlaneHeightAtDistOne=viewpoint.getImgPlaneHeightAtDistOne();var near=viewpoint.getNear();var center=model_view.multMatrixPnt(popGeo._vf.position);var tightRad=model_view.multMatrixVec(popGeo._vf.size).length()*0.5;var largestRad=model_view.multMatrixVec(popGeo._vf.maxBBSize).length()*0.5;var dist=Math.max(-center.z-tightRad,near);var projPixelLength=dist*(imgPlaneHeightAtDistOne/viewarea._height);var arg=(2*largestRad)/(tol*projPixelLength);currentLOD=Math.ceil(Math.log(arg)/0.693147180559945);currentLOD=(currentLOD<1)?1:((currentLOD>16)?16:currentLOD);}
var minPrec=popGeo._vf.minPrecisionLevel,maxPrec=popGeo._vf.maxPrecisionLevel;currentLOD=(minPrec!=-1&&currentLOD<minPrec)?minPrec:currentLOD;currentLOD=(maxPrec!=-1&&currentLOD>maxPrec)?maxPrec:currentLOD;var currentLOD_min=(s_gl.levelsAvailable<currentLOD)?s_gl.levelsAvailable:currentLOD;currentLOD=currentLOD_min;if(tol<=1)
currentLOD=(currentLOD==popGeo.getNumLevels())?16:currentLOD;var hasIndex=popGeo._vf.indexedRendering;var p_msh=popGeo._mesh;p_msh._numCoords=0;p_msh._numFaces=0;for(var i=0;i<currentLOD_min;++i){var numVerticesAtLevel_i=s_gl.numVerticesAtLevel[i];p_msh._numCoords+=numVerticesAtLevel_i;p_msh._numFaces+=(hasIndex?popGeo.getNumIndicesByLevel(i):numVerticesAtLevel_i)/3;}
x3dom.nodeTypes.PopGeometry.numRenderedVerts+=p_msh._numCoords;x3dom.nodeTypes.PopGeometry.numRenderedTris+=p_msh._numFaces;p_msh.currentLOD=currentLOD;popGeo.adaptVertexCount(hasIndex?p_msh._numFaces*3:p_msh._numCoords);sp.PG_maxBBSize=popGeo._vf.maxBBSize.toGL();sp.PG_bbMin=popGeo._bbMinBySize;sp.PG_numAnchorVertices=popGeo._vf.numAnchorVertices;sp.PG_bbMaxModF=popGeo._vf.bbMaxModF.toGL();sp.PG_bboxShiftVec=popGeo._vf.bbShiftVec.toGL();sp.PG_precisionLevel=currentLOD;sp.PG_powPrecision=x3dom.nodeTypes.PopGeometry.powLUT[currentLOD-1];};Context.prototype.pickValue=function(viewarea,x,y,buttonState,viewMat,sceneMat)
{x3dom.Utils.startMeasure("picking");var scene=viewarea._scene;var gl=this.ctx3d;if(!gl||!scene||!scene._webgl||!scene.drawableCollection){return false;}
var pm=scene._vf.pickMode.toLowerCase();var pickMode=0;switch(pm){case"box":return false;case"idbuf":pickMode=0;break;case"idbuf24":pickMode=3;break;case"idbufid":pickMode=4;break;case"color":pickMode=1;break;case"texcoord":pickMode=2;break;}
var mat_view,mat_scene;if(arguments.length>4){mat_view=viewMat;mat_scene=sceneMat;}
else{mat_view=viewarea._last_mat_view;mat_scene=viewarea._last_mat_scene;}
var min=x3dom.fields.SFVec3f.copy(scene._lastMin);var max=x3dom.fields.SFVec3f.copy(scene._lastMax);var from=mat_view.inverse().e3();var _min=x3dom.fields.SFVec3f.copy(from);var _max=x3dom.fields.SFVec3f.copy(from);if(_min.x>min.x){_min.x=min.x;}
if(_min.y>min.y){_min.y=min.y;}
if(_min.z>min.z){_min.z=min.z;}
if(_max.x<max.x){_max.x=max.x;}
if(_max.y<max.y){_max.y=max.y;}
if(_max.z<max.z){_max.z=max.z;}
scene._lastMin.setValues(_min);scene._lastMax.setValues(_max);var sceneSize=scene._lastMax.subtract(scene._lastMin).length();var cctowc=viewarea.getCCtoWCMatrix();scene._lastMin.setValues(min);scene._lastMax.setValues(max);var baseID=x3dom.nodeTypes.Shape.objectID+2;this.renderPickingPass(gl,scene,mat_view,mat_scene,from,sceneSize,pickMode,x,y,2,2);var pixelData=scene._webgl.fboPick.pixelData;if(pixelData&&pixelData.length)
{var pickPos=new x3dom.fields.SFVec3f(0,0,0);var pickNorm=new x3dom.fields.SFVec3f(0,0,1);var index=0;var objId=pixelData[index+3],shapeId;var pixelOffset=1.0/scene._webgl.pickScale;var denom=1.0/256.0;var dist,line,lineoff,right,up;if(pickMode==0){objId+=256*pixelData[index+2];dist=(pixelData[index]/255.0)*denom+
(pixelData[index+1]/255.0);line=viewarea.calcViewRay(x,y,cctowc);pickPos=line.pos.add(line.dir.multiply(dist*sceneSize));index=4;dist=(pixelData[index]/255.0)*denom+
(pixelData[index+1]/255.0);lineoff=viewarea.calcViewRay(x+pixelOffset,y,cctowc);right=lineoff.pos.add(lineoff.dir.multiply(dist*sceneSize));right=right.subtract(pickPos).normalize();index=8;dist=(pixelData[index]/255.0)*denom+
(pixelData[index+1]/255.0);lineoff=viewarea.calcViewRay(x,y-pixelOffset,cctowc);up=lineoff.pos.add(lineoff.dir.multiply(dist*sceneSize));up=up.subtract(pickPos).normalize();pickNorm=right.cross(up).normalize();}
else if(pickMode==3){objId+=256*pixelData[index+2]+
65536*pixelData[index+1];dist=pixelData[index]/255.0;line=viewarea.calcViewRay(x,y,cctowc);pickPos=line.pos.add(line.dir.multiply(dist*sceneSize));index=4;dist=pixelData[index]/255.0;lineoff=viewarea.calcViewRay(x+pixelOffset,y,cctowc);right=lineoff.pos.add(lineoff.dir.multiply(dist*sceneSize));right=right.subtract(pickPos).normalize();index=8;dist=pixelData[index]/255.0;lineoff=viewarea.calcViewRay(x,y-pixelOffset,cctowc);up=lineoff.pos.add(lineoff.dir.multiply(dist*sceneSize));up=up.subtract(pickPos).normalize();pickNorm=right.cross(up).normalize();}
else if(pickMode==4){objId+=256*pixelData[index+2];shapeId=pixelData[index+1];shapeId+=256*pixelData[index];if(objId==0&&(shapeId>0&&shapeId<baseID)){objId=shapeId;}}
else{pickPos.x=pixelData[index];pickPos.y=pixelData[index+1];pickPos.z=pixelData[index+2];}
var eventType="shadowObjectIdChanged";var shadowObjectIdChanged,event;var button=Math.max(buttonState>>>8,buttonState&255);if(objId>=baseID){objId-=baseID;var hitObject;if(pickMode!=4){viewarea._pickingInfo.pickPos=pickPos;viewarea._pick.setValues(pickPos);viewarea._pickingInfo.pickNorm=pickNorm;viewarea._pickNorm.setValues(pickNorm);viewarea._pickingInfo.pickObj=null;viewarea._pickingInfo.lastClickObj=null;hitObject=scene._xmlNode;}
else{viewarea._pickingInfo.pickObj=x3dom.nodeTypes.Shape.idMap.nodeID[shapeId];hitObject=viewarea._pickingInfo.pickObj._xmlNode;}
if(scene._multiPartMap){var mp,multiPart;for(mp=0;mp<scene._multiPartMap.multiParts.length;mp++)
{multiPart=scene._multiPartMap.multiParts[mp];if(objId>=multiPart._minId&&objId<=multiPart._maxId)
{hitObject=multiPart._xmlNode;event={target:multiPart._xmlNode,button:button,mouseup:((buttonState>>>8)>0),layerX:x,layerY:y,pickedId:objId,worldX:pickPos.x,worldY:pickPos.y,worldZ:pickPos.z,normalX:pickNorm.x,normalY:pickNorm.y,normalZ:pickNorm.z,hitPnt:pickPos.toGL(),hitObject:hitObject,cancelBubble:false,stopPropagation:function(){this.cancelBubble=true;},preventDefault:function(){this.cancelBubble=true;}};multiPart.handleEvents(event);}
else
{event={target:multiPart._xmlNode,button:button,mouseup:((buttonState>>>8)>0),layerX:x,layerY:y,pickedId:-1,cancelBubble:false,stopPropagation:function(){this.cancelBubble=true;},preventDefault:function(){this.cancelBubble=true;}};multiPart.handleEvents(event);}}}
shadowObjectIdChanged=(viewarea._pickingInfo.shadowObjectId!=objId);viewarea._pickingInfo.lastShadowObjectId=viewarea._pickingInfo.shadowObjectId;viewarea._pickingInfo.shadowObjectId=objId;if((shadowObjectIdChanged||button)&&scene._xmlNode&&(scene._xmlNode["on"+eventType]||scene._xmlNode.hasAttribute("on"+eventType)||scene._listeners[eventType]))
{event={target:scene._xmlNode,type:eventType,button:button,mouseup:((buttonState>>>8)>0),layerX:x,layerY:y,shadowObjectId:objId,worldX:pickPos.x,worldY:pickPos.y,worldZ:pickPos.z,normalX:pickNorm.x,normalY:pickNorm.y,normalZ:pickNorm.z,hitPnt:pickPos.toGL(),hitObject:hitObject,cancelBubble:false,stopPropagation:function(){this.cancelBubble=true;},preventDefault:function(){this.cancelBubble=true;}};scene.callEvtHandler(("on"+eventType),event);}
if(scene._shadowIdMap&&scene._shadowIdMap.mapping&&objId<scene._shadowIdMap.mapping.length){var shIds=scene._shadowIdMap.mapping[objId].usage;var n,c,shObj;if(!line){line=viewarea.calcViewRay(x,y,cctowc);}
for(c=0;c<shIds.length;c++){shObj=scene._nameSpace.defMap[shIds[c]];if(shObj&&shObj.doIntersect(line)){viewarea._pickingInfo.pickObj=shObj;break;}}
for(n=0;n<scene._nameSpace.childSpaces.length;n++)
{for(c=0;c<shIds.length;c++){shObj=scene._nameSpace.childSpaces[n].defMap[shIds[c]];if(shObj&&shObj.doIntersect(line)){viewarea._pickingInfo.pickObj=shObj;break;}}}}}
else{if(scene._multiPartMap){for(mp=0;mp<scene._multiPartMap.multiParts.length;mp++)
{multiPart=scene._multiPartMap.multiParts[mp];event={target:multiPart._xmlNode,button:button,mouseup:((buttonState>>>8)>0),layerX:x,layerY:y,pickedId:-1,cancelBubble:false,stopPropagation:function(){this.cancelBubble=true;},preventDefault:function(){this.cancelBubble=true;}};multiPart.handleEvents(event);}}
shadowObjectIdChanged=(viewarea._pickingInfo.shadowObjectId!=-1);viewarea._pickingInfo.shadowObjectId=-1;if(shadowObjectIdChanged&&scene._xmlNode&&(scene._xmlNode["on"+eventType]||scene._xmlNode.hasAttribute("on"+eventType)||scene._listeners[eventType]))
{event={target:scene._xmlNode,type:eventType,button:button,mouseup:((buttonState>>>8)>0),layerX:x,layerY:y,shadowObjectId:viewarea._pickingInfo.shadowObjectId,cancelBubble:false,stopPropagation:function(){this.cancelBubble=true;},preventDefault:function(){this.cancelBubble=true;}};scene.callEvtHandler(("on"+eventType),event);}
if(objId>0){viewarea._pickingInfo.pickPos=pickPos;viewarea._pickingInfo.pickNorm=pickNorm;viewarea._pickingInfo.pickObj=x3dom.nodeTypes.Shape.idMap.nodeID[objId];}
else{viewarea._pickingInfo.pickObj=null;viewarea._pickingInfo.lastClickObj=null;}}}
var pickTime=x3dom.Utils.stopMeasure("picking");this.x3dElem.runtime.addMeasurement('PICKING',pickTime);return true;};Context.prototype.pickRect=function(viewarea,x1,y1,x2,y2)
{var gl=this.ctx3d;var scene=viewarea?viewarea._scene:null;if(!gl||!scene||!scene._webgl||!scene.drawableCollection)
return false;var from=viewarea._last_mat_view.inverse().e3();var sceneSize=scene._lastMax.subtract(scene._lastMin).length();var x=(x1<=x2)?x1:x2;var y=(y1>=y2)?y1:y2;var width=(1+Math.abs(x2-x1))*scene._webgl.pickScale;var height=(1+Math.abs(y2-y1))*scene._webgl.pickScale;this.renderPickingPass(gl,scene,viewarea._last_mat_view,viewarea._last_mat_scene,from,sceneSize,0,x,y,(width<1)?1:width,(height<1)?1:height);var index;var pickedObjects=[];for(index=0;scene._webgl.fboPick.pixelData&&index<scene._webgl.fboPick.pixelData.length;index+=4){var objId=scene._webgl.fboPick.pixelData[index+3]+
scene._webgl.fboPick.pixelData[index+2]*256;if(objId>0)
pickedObjects.push(objId);}
pickedObjects.sort();var pickedObjectsTemp=(function(arr){var a=[],l=arr.length;for(var i=0;i<l;i++){for(var j=i+1;j<l;j++){if(arr[i]===arr[j])
j=++i;}
a.push(arr[i]);}
return a;})(pickedObjects);pickedObjects=pickedObjectsTemp;var pickedNode,pickedNodes=[];var hitObject;var baseID=x3dom.nodeTypes.Shape.objectID+2;for(index=0;index<pickedObjects.length;index++){objId=pickedObjects[index];if(objId>=baseID)
{objId-=baseID;if(scene._multiPartMap){var mp,multiPart,colorMap,emissiveMap,specularMap,visibilityMap,partID;for(mp=0;mp<scene._multiPartMap.multiParts.length;mp++){multiPart=scene._multiPartMap.multiParts[mp];colorMap=multiPart._inlineNamespace.defMap["MultiMaterial_ColorMap"];emissiveMap=multiPart._inlineNamespace.defMap["MultiMaterial_EmissiveMap"];specularMap=multiPart._inlineNamespace.defMap["MultiMaterial_SpecularMap"];visibilityMap=multiPart._inlineNamespace.defMap["MultiMaterial_VisibilityMap"];if(objId>=multiPart._minId&&objId<=multiPart._maxId){partID=multiPart._idMap.mapping[objId-multiPart._minId].name;hitObject=new x3dom.Parts(multiPart,[objId],colorMap,emissiveMap,specularMap,visibilityMap);pickedNode={"partID":partID,"part":hitObject};pickedNodes.push(pickedNode);}}}}
else
{hitObject=x3dom.nodeTypes.Shape.idMap.nodeID[objId];hitObject=(hitObject&&hitObject._xmlNode)?hitObject._xmlNode:null;if(hitObject)
pickedNodes.push(hitObject);}}
return pickedNodes;};Context.prototype.renderScene=function(viewarea)
{var gl=this.ctx3d;var scene=viewarea._scene;if(gl===null||scene===null){return;}
var rentex=viewarea._doc._nodeBag.renderTextures;var rt_tex,rtl_i,rtl_n=rentex.length;var texProp=null;var type=gl.UNSIGNED_BYTE;var shadowType=gl.UNSIGNED_BYTE;var nearestFilt=false;if(x3dom.caps.FP_TEXTURES&&!x3dom.caps.MOBILE){type=gl.FLOAT;shadowType=gl.FLOAT;if(!x3dom.caps.FPL_TEXTURES){nearestFilt=true;}}
var shadowedLights,numShadowMaps;var i,j,n,size,sizeAvailable;var texType,refinementPos;var vertices=[-1,-1,1,-1,-1,1,-1,1,1,-1,1,1];scene.updateVolume();if(!scene._webgl)
{scene._webgl={};this.setupFgnds(gl,scene);scene._webgl.pickScale=0.5;scene._webgl._currFboWidth=Math.round(this.canvas.width*scene._webgl.pickScale);scene._webgl._currFboHeight=Math.round(this.canvas.height*scene._webgl.pickScale);scene._webgl.fboPick=x3dom.Utils.initFBO(gl,scene._webgl._currFboWidth,scene._webgl._currFboHeight,gl.UNSIGNED_BYTE,false,true);scene._webgl.fboPick.pixelData=null;scene._webgl.normalShader=this.cache.getShader(gl,x3dom.shader.NORMAL);scene._webgl.fboShadow=[];shadowedLights=viewarea.getShadowedLights();n=shadowedLights.length;for(i=0;i<n;i++)
{size=shadowedLights[i]._vf.shadowMapSize;if(!x3dom.isa(shadowedLights[i],x3dom.nodeTypes.PointLight))
numShadowMaps=Math.max(1,Math.min(shadowedLights[i]._vf.shadowCascades,6));else
numShadowMaps=6;scene._webgl.fboShadow[i]=[];for(j=0;j<numShadowMaps;j++)
scene._webgl.fboShadow[i][j]=x3dom.Utils.initFBO(gl,size,size,shadowType,false,true);}
if(scene._webgl.fboShadow.length>0||x3dom.SSAO.isEnabled(scene))
scene._webgl.fboScene=x3dom.Utils.initFBO(gl,this.canvas.width,this.canvas.height,shadowType,false,true);scene._webgl.fboBlur=[];for(i=0;i<n;i++)
{size=scene._webgl.fboShadow[i][0].height;sizeAvailable=false;for(j=0;j<scene._webgl.fboBlur.length;j++){if(size==scene._webgl.fboBlur[j].height)
sizeAvailable=true;}
if(!sizeAvailable)
scene._webgl.fboBlur[scene._webgl.fboBlur.length]=x3dom.Utils.initFBO(gl,size,size,shadowType,false,true);}
scene._webgl.ppBuffer=gl.createBuffer();gl.bindBuffer(gl.ARRAY_BUFFER,scene._webgl.ppBuffer);gl.bufferData(gl.ARRAY_BUFFER,new Float32Array(vertices),gl.STATIC_DRAW);scene._webgl.refinement={stamps:new Array(2),positionBuffer:gl.createBuffer()};gl.bindBuffer(gl.ARRAY_BUFFER,scene._webgl.refinement.positionBuffer);gl.bufferData(gl.ARRAY_BUFFER,new Float32Array(vertices),gl.STATIC_DRAW);for(rtl_i=0;rtl_i<rtl_n;rtl_i++){rt_tex=rentex[rtl_i];texProp=rt_tex._cf.textureProperties.node;texType=rt_tex.requirePingPong()?gl.UNSIGNED_BYTE:type;rt_tex._webgl={};rt_tex._webgl.fbo=x3dom.Utils.initFBO(gl,rt_tex._vf.dimensions[0],rt_tex._vf.dimensions[1],texType,(texProp&&texProp._vf.generateMipMaps),rt_tex._vf.depthMap||!rt_tex.requirePingPong());rt_tex._cleanupGLObjects=function(retainTex){if(!retainTex)
gl.deleteTexture(this._webgl.fbo.tex);if(this._webgl.fbo.dtex)
gl.deleteTexture(this._webgl.fbo.dtex);if(this._webgl.fbo.rbo)
gl.deleteFramebuffer(this._webgl.fbo.rbo);gl.bindFramebuffer(gl.FRAMEBUFFER,null);gl.deleteFramebuffer(this._webgl.fbo.fbo);this._webgl.fbo.rbo=null;this._webgl.fbo.fbo=null;};if(rt_tex.requirePingPong()){refinementPos=rt_tex._vf.dimensions[0]+"x"+rt_tex._vf.dimensions[1];if(scene._webgl.refinement[refinementPos]===undefined){scene._webgl.refinement[refinementPos]=x3dom.Utils.initFBO(gl,rt_tex._vf.dimensions[0],rt_tex._vf.dimensions[1],texType,false,false);}
rt_tex._webgl.texture=null;}}
viewarea._last_mat_view=x3dom.fields.SFMatrix4f.identity();viewarea._last_mat_proj=x3dom.fields.SFMatrix4f.identity();viewarea._last_mat_scene=x3dom.fields.SFMatrix4f.identity();this._calledViewpointChangedHandler=false;}
else
{var fboWidth=Math.round(this.canvas.width*scene._webgl.pickScale);var fboHeight=Math.round(this.canvas.height*scene._webgl.pickScale);if(scene._webgl._currFboWidth!==fboWidth||scene._webgl._currFboHeight!==fboHeight){scene._webgl._currFboWidth=fboWidth;scene._webgl._currFboHeight=fboHeight;scene._webgl.fboPick=x3dom.Utils.initFBO(gl,fboWidth,fboHeight,scene._webgl.fboPick.type,false,true);scene._webgl.fboPick.pixelData=null;x3dom.debug.logInfo("Refreshed picking FBO to size ("+fboWidth+", "+fboHeight+")");}
for(rtl_i=0;rtl_i<rtl_n;rtl_i++){rt_tex=rentex[rtl_i];if(rt_tex._webgl&&rt_tex._webgl.fbo&&rt_tex._webgl.fbo.width==rt_tex._vf.dimensions[0]&&rt_tex._webgl.fbo.height==rt_tex._vf.dimensions[1])
continue;rt_tex.invalidateGLObject();if(rt_tex._cleanupGLObjects)
rt_tex._cleanupGLObjects();else
rt_tex._cleanupGLObjects=function(retainTex){if(!retainTex)
gl.deleteTexture(this._webgl.fbo.tex);if(this._webgl.fbo.dtex)
gl.deleteTexture(this._webgl.fbo.dtex);if(this._webgl.fbo.rbo)
gl.deleteRenderbuffer(this._webgl.fbo.rbo);gl.bindFramebuffer(gl.FRAMEBUFFER,null);gl.deleteFramebuffer(this._webgl.fbo.fbo);this._webgl.fbo.rbo=null;this._webgl.fbo.fbo=null;};texProp=rt_tex._cf.textureProperties.node;texType=rt_tex.requirePingPong()?gl.UNSIGNED_BYTE:type;rt_tex._webgl={};rt_tex._webgl.fbo=x3dom.Utils.initFBO(gl,rt_tex._vf.dimensions[0],rt_tex._vf.dimensions[1],texType,(texProp&&texProp._vf.generateMipMaps),rt_tex._vf.depthMap||!rt_tex.requirePingPong());if(rt_tex.requirePingPong()){refinementPos=rt_tex._vf.dimensions[0]+"x"+rt_tex._vf.dimensions[1];if(scene._webgl.refinement[refinementPos]===undefined){scene._webgl.refinement[refinementPos]=x3dom.Utils.initFBO(gl,rt_tex._vf.dimensions[0],rt_tex._vf.dimensions[1],texType,false,false);}
rt_tex._webgl.texture=null;}
x3dom.debug.logInfo("Init/resize RenderedTexture_"+rtl_i+" to size "+
rt_tex._vf.dimensions[0]+" x "+rt_tex._vf.dimensions[1]);}
shadowedLights=viewarea.getShadowedLights();n=shadowedLights.length;for(i=0;i<n;i++){size=shadowedLights[i]._vf.shadowMapSize;if(!x3dom.isa(shadowedLights[i],x3dom.nodeTypes.PointLight))
numShadowMaps=Math.max(1,Math.min(shadowedLights[i]._vf.shadowCascades,6));else
numShadowMaps=6;if(typeof scene._webgl.fboShadow[i]==="undefined"||scene._webgl.fboShadow[i].length!=numShadowMaps||scene._webgl.fboShadow[i][0].height!=size){scene._webgl.fboShadow[i]=[];for(j=0;j<numShadowMaps;j++){scene._webgl.fboShadow[i][j]=x3dom.Utils.initFBO(gl,size,size,shadowType,false,true);}}}
for(i=0;i<n;i++){size=scene._webgl.fboShadow[i][0].height;sizeAvailable=false;for(j=0;j<scene._webgl.fboBlur.length;j++){if(size==scene._webgl.fboBlur[j].height)
sizeAvailable=true;}
if(!sizeAvailable)
scene._webgl.fboBlur[scene._webgl.fboBlur.length]=x3dom.Utils.initFBO(gl,size,size,shadowType,false,true);}
if((x3dom.SSAO.isEnabled(scene)||scene._webgl.fboShadow.length>0)&&typeof scene._webgl.fboScene=="undefined"||scene._webgl.fboScene&&(this.canvas.width!=scene._webgl.fboScene.width||this.canvas.height!=scene._webgl.fboScene.height)){scene._webgl.fboScene=x3dom.Utils.initFBO(gl,this.canvas.width,this.canvas.height,shadowType,false,true);}}
var env=scene.getEnvironment();env.checkSanity();var bgnd=scene.getBackground();this.setupScene(gl,bgnd);this.numFaces=0;this.numCoords=0;this.numDrawCalls=0;var mat_proj=viewarea.getProjectionMatrix();var mat_view=viewarea.getViewMatrix();if(!this._calledViewpointChangedHandler||!viewarea._last_mat_view.equals(mat_view)){var e_viewpoint=scene.getViewpoint();var e_eventType="viewpointChanged";try{if(e_viewpoint._xmlNode&&(e_viewpoint._xmlNode["on"+e_eventType]||e_viewpoint._xmlNode.hasAttribute("on"+e_eventType)||e_viewpoint._listeners[e_eventType])){var e_viewtrafo=e_viewpoint.getCurrentTransform();e_viewtrafo=e_viewtrafo.inverse().mult(mat_view);var e_mat=e_viewtrafo.inverse();var e_rotation=new x3dom.fields.Quaternion(0,0,1,0);e_rotation.setValue(e_mat);var e_translation=e_mat.e3();var e_event={target:e_viewpoint._xmlNode,type:e_eventType,matrix:e_viewtrafo,position:e_translation,orientation:e_rotation.toAxisAngle(),cancelBubble:false,stopPropagation:function(){this.cancelBubble=true;},preventDefault:function(){this.cancelBubble=true;}};e_viewpoint.callEvtHandler(("on"+e_eventType),e_event);this._calledViewpointChangedHandler=true;}}
catch(e_e){x3dom.debug.logException(e_e);}}
viewarea._last_mat_view=mat_view;viewarea._last_mat_proj=mat_proj;var mat_scene=mat_proj.mult(mat_view);viewarea._last_mat_scene=mat_scene;scene.drawableCollection=null;if(!scene.drawableCollection)
{var drawableCollectionConfig={viewArea:viewarea,sortTrans:env._vf.sortTrans,viewMatrix:mat_view,projMatrix:mat_proj,sceneMatrix:mat_scene,frustumCulling:true,smallFeatureThreshold:env._smallFeatureThreshold,context:this,gl:gl};scene.drawableCollection=new x3dom.DrawableCollection(drawableCollectionConfig);x3dom.Utils.startMeasure('traverse');scene.collectDrawableObjects(x3dom.fields.SFMatrix4f.identity(),scene.drawableCollection,true,false,0,[]);var traverseTime=x3dom.Utils.stopMeasure('traverse');this.x3dElem.runtime.addMeasurement('TRAVERSE',traverseTime);}
x3dom.Utils.startMeasure('sorting');scene.drawableCollection.sort();var sortTime=x3dom.Utils.stopMeasure('sorting');this.x3dElem.runtime.addMeasurement('SORT',sortTime);var slights=viewarea.getLights();var numLights=slights.length;var mat_light;var WCToLCMatrices=[];var lMatrices=[];var shadowCount=0;x3dom.Utils.startMeasure('shadow');for(var p=0;p<numLights;p++){if(slights[p]._vf.shadowIntensity>0.0){var lightMatrix=viewarea.getLightMatrix()[p];shadowMaps=scene._webgl.fboShadow[shadowCount];var offset=Math.max(0.0,Math.min(1.0,slights[p]._vf.shadowOffset));if(!x3dom.isa(slights[p],x3dom.nodeTypes.PointLight)){var numCascades=Math.max(1,Math.min(slights[p]._vf.shadowCascades,6));mat_light=viewarea.getWCtoLCMatricesCascaded(lightMatrix,slights[p],mat_proj);for(i=0;i<numCascades;i++){this.renderShadowPass(gl,viewarea,mat_light[i],mat_view,shadowMaps[i],offset,false);}}
else{mat_light=viewarea.getWCtoLCMatricesPointLight(lightMatrix,slights[p],mat_proj);for(i=0;i<6;i++){this.renderShadowPass(gl,viewarea,mat_light[i],mat_view,shadowMaps[i],offset,false);}}
shadowCount++;WCToLCMatrices[WCToLCMatrices.length]=mat_light;lMatrices[lMatrices.length]=lightMatrix;}}
if(shadowCount>0||x3dom.SSAO.isEnabled(scene)){this.renderShadowPass(gl,viewarea,mat_scene,mat_view,scene._webgl.fboScene,0.0,true);var shadowTime=x3dom.Utils.stopMeasure('shadow');this.x3dElem.runtime.addMeasurement('SHADOW',shadowTime);}
else{this.x3dElem.runtime.removeMeasurement('SHADOW');}
mat_light=viewarea.getWCtoLCMatrix(viewarea.getLightMatrix()[0]);for(rtl_i=0;rtl_i<rtl_n;rtl_i++){this.renderRTPass(gl,viewarea,rentex[rtl_i]);}
x3dom.Utils.startMeasure('render');this.stateManager.viewport(0,0,this.canvas.width,this.canvas.height);bgnd._webgl.render(gl,mat_view,mat_proj);x3dom.nodeTypes.PopGeometry.numRenderedVerts=0;x3dom.nodeTypes.PopGeometry.numRenderedTris=0;n=scene.drawableCollection.length;if(env._vf.smallFeatureCulling&&env._lowPriorityThreshold<1&&viewarea.isMovingOrAnimating()){n=Math.floor(n*env._lowPriorityThreshold);if(!n&&scene.drawableCollection.length)
n=1;}
this.stateManager.unsetProgram();for(i=0;i<n;i++){var drawable=scene.drawableCollection.get(i);this.renderShape(drawable,viewarea,slights,numLights,mat_view,mat_scene,mat_light,mat_proj,gl);}
if(shadowCount>0)
this.renderShadows(gl,viewarea,shadowedLights,WCToLCMatrices,lMatrices,mat_view,mat_proj,mat_scene);this.stateManager.disable(gl.BLEND);this.stateManager.disable(gl.DEPTH_TEST);viewarea._numRenderedNodes=n;if(x3dom.SSAO.isEnabled(scene))
x3dom.SSAO.renderSSAO(this.stateManager,gl,scene,this.canvas);if(viewarea._visDbgBuf!==undefined&&viewarea._visDbgBuf)
{var pm=scene._vf.pickMode.toLowerCase();if(pm.indexOf("idbuf")==0||pm=="color"||pm=="texcoord"){this.stateManager.viewport(0,3*this.canvas.height/4,this.canvas.width/4,this.canvas.height/4);scene._fgnd._webgl.render(gl,scene._webgl.fboPick.tex);}
if(shadowCount>0||x3dom.SSAO.isEnabled(scene)){this.stateManager.viewport(this.canvas.width/4,3*this.canvas.height/4,this.canvas.width/4,this.canvas.height/4);scene._fgnd._webgl.render(gl,scene._webgl.fboScene.tex);}
var row=3,col=2;for(i=0;i<shadowCount;i++){var shadowMaps=scene._webgl.fboShadow[i];for(j=0;j<shadowMaps.length;j++){this.stateManager.viewport(col*this.canvas.width/4,row*this.canvas.height/4,this.canvas.width/4,this.canvas.height/4);scene._fgnd._webgl.render(gl,shadowMaps[j].tex);if(col<2){col++;}else{col=0;row--;}}}
for(rtl_i=0;rtl_i<rtl_n;rtl_i++){rt_tex=rentex[rtl_i];if(!rt_tex._webgl.fbo.fbo)
continue;this.stateManager.viewport(rtl_i*this.canvas.width/8,5*this.canvas.height/8,this.canvas.width/8,this.canvas.height/8);scene._fgnd._webgl.render(gl,rt_tex._webgl.fbo.tex);}}
gl.finish();var renderTime=x3dom.Utils.stopMeasure('render');this.x3dElem.runtime.addMeasurement('RENDER',renderTime);this.x3dElem.runtime.addMeasurement('DRAW',(n?renderTime/n:0));this.x3dElem.runtime.addInfo('#NODES:',scene.drawableCollection.numberOfNodes);this.x3dElem.runtime.addInfo('#SHAPES:',viewarea._numRenderedNodes);this.x3dElem.runtime.addInfo("#DRAWS:",this.numDrawCalls);this.x3dElem.runtime.addInfo("#POINTS:",this.numCoords);this.x3dElem.runtime.addInfo("#TRIS:",this.numFaces);};Context.prototype.renderPingPongPass=function(gl,viewarea,rt){var scene=viewarea._scene;var refinementPos=rt._vf.dimensions[0]+"x"+rt._vf.dimensions[1];var refinementFbo=scene._webgl.refinement[refinementPos];if(rt._currLoadLevel==0&&(!scene._webgl.refinement.stamps[0]||!scene._webgl.refinement.stamps[1])){scene._webgl.refinement.stamps[0]=this.cache.getTexture2D(gl,rt._nameSpace.doc,rt._nameSpace.getURL(rt._vf.stamp0),false,false,false,false);scene._webgl.refinement.stamps[1]=this.cache.getTexture2D(gl,rt._nameSpace.doc,rt._nameSpace.getURL(rt._vf.stamp1),false,false,false,false);}
if(rt._currLoadLevel<rt._loadLevel){rt._currLoadLevel++;if(rt._webgl.texture)
gl.deleteTexture(rt._webgl.texture);var filename=rt._vf.url[0]+"/"+rt._currLoadLevel+"."+rt._vf.format;rt._webgl.texture=x3dom.Utils.createTexture2D(gl,rt._nameSpace.doc,rt._nameSpace.getURL(filename),false,false,false,false);if(rt._vf.iterations%2===0)
(rt._currLoadLevel%2!==0)?rt._repeat.x*=2.0:rt._repeat.y*=2.0;else
(rt._currLoadLevel%2===0)?rt._repeat.x*=2.0:rt._repeat.y*=2.0;}
if(!rt._webgl.texture.ready||!scene._webgl.refinement.stamps[0].ready||!scene._webgl.refinement.stamps[1].ready)
return;this.stateManager.bindFramebuffer(gl.FRAMEBUFFER,refinementFbo.fbo);this.stateManager.viewport(0,0,refinementFbo.width,refinementFbo.height);this.stateManager.disable(gl.BLEND);this.stateManager.disable(gl.CULL_FACE);this.stateManager.disable(gl.DEPTH_TEST);gl.clearColor(0,0,0,1);gl.clearDepth(1);gl.clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER_BIT);var sp=this.cache.getShader(gl,x3dom.shader.TEXTURE_REFINEMENT);this.stateManager.useProgram(sp);gl.bindBuffer(gl.ARRAY_BUFFER,scene._webgl.refinement.positionBuffer);gl.vertexAttribPointer(sp.position,2,gl.FLOAT,false,0,0);gl.enableVertexAttribArray(sp.position);sp.stamp=0;gl.activeTexture(gl.TEXTURE0);gl.bindTexture(gl.TEXTURE_2D,scene._webgl.refinement.stamps[(rt._currLoadLevel+1)%2]);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MAG_FILTER,gl.NEAREST);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MIN_FILTER,gl.NEAREST);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_S,gl.REPEAT);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_T,gl.REPEAT);if(rt._currLoadLevel>1){sp.lastTex=1;gl.activeTexture(gl.TEXTURE1);gl.bindTexture(gl.TEXTURE_2D,rt._webgl.fbo.tex);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MAG_FILTER,gl.NEAREST);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MIN_FILTER,gl.NEAREST);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_S,gl.CLAMP_TO_EDGE);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_T,gl.CLAMP_TO_EDGE);}
sp.curTex=2;gl.activeTexture(gl.TEXTURE2);gl.bindTexture(gl.TEXTURE_2D,rt._webgl.texture);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MAG_FILTER,gl.NEAREST);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MIN_FILTER,gl.NEAREST);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_S,gl.CLAMP_TO_EDGE);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_T,gl.CLAMP_TO_EDGE);sp.mode=rt._currLoadLevel-1;sp.repeat=rt._repeat.toGL();gl.drawArrays(gl.TRIANGLES,0,6);this.stateManager.bindFramebuffer(gl.FRAMEBUFFER,rt._webgl.fbo.fbo);gl.clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER_BIT);sp.mode=0;sp.curTex=2;gl.activeTexture(gl.TEXTURE2);gl.bindTexture(gl.TEXTURE_2D,refinementFbo.tex);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MAG_FILTER,gl.NEAREST);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MIN_FILTER,gl.NEAREST);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_S,gl.CLAMP_TO_EDGE);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_T,gl.CLAMP_TO_EDGE);gl.drawArrays(gl.TRIANGLES,0,6);gl.activeTexture(gl.TEXTURE0);gl.bindTexture(gl.TEXTURE_2D,null);gl.disableVertexAttribArray(sp.position);this.stateManager.bindFramebuffer(gl.FRAMEBUFFER,null);this.stateManager.viewport(0,0,this.canvas.width,this.canvas.height);if(rt._vf.autoRefinement)
rt.nextLevel();if(rt._currLoadLevel==rt._vf.maxLevel)
rt._currLoadLevel++;if(rt._webgl.fbo.mipMap){gl.bindTexture(gl.TEXTURE_2D,rt._webgl.fbo.tex);gl.generateMipmap(gl.TEXTURE_2D);gl.bindTexture(gl.TEXTURE_2D,null);}
if(!rt.requirePingPong()){gl.deleteTexture(rt._webgl.texture);delete rt._webgl.texture;rt._cleanupGLObjects(true);}
rt._renderedImage++;};Context.prototype.renderRTPass=function(gl,viewarea,rt)
{if(x3dom.isa(rt,x3dom.nodeTypes.RefinementTexture)){if(rt.requirePingPong()){this.renderPingPongPass(gl,viewarea,rt);}
return;}
switch(rt._vf.update.toUpperCase()){case"NONE":return;case"NEXT_FRAME_ONLY":if(!rt._needRenderUpdate){return;}
rt._needRenderUpdate=false;break;case"ALWAYS":default:break;}
var scene=viewarea._scene;var bgnd=null;var mat_view=rt.getViewMatrix();var mat_proj=rt.getProjectionMatrix();var mat_scene=mat_proj.mult(mat_view);var lightMatrix=viewarea.getLightMatrix()[0];var mat_light=viewarea.getWCtoLCMatrix(lightMatrix);var i,n,m=rt._cf.excludeNodes.nodes.length;var arr=new Array(m);for(i=0;i<m;i++){var render=rt._cf.excludeNodes.nodes[i]._vf.render;if(render===undefined){arr[i]=-1;}
else{if(render===true){arr[i]=1;}else{arr[i]=0;}}
rt._cf.excludeNodes.nodes[i]._vf.render=false;}
this.stateManager.bindFramebuffer(gl.FRAMEBUFFER,rt._webgl.fbo.fbo);this.stateManager.viewport(0,0,rt._webgl.fbo.width,rt._webgl.fbo.height);if(rt._cf.background.node===null){gl.clearColor(0,0,0,1);gl.clearDepth(1.0);gl.clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER_BIT|gl.STENCIL_BUFFER_BIT);}
else if(rt._cf.background.node===scene.getBackground()){bgnd=scene.getBackground();bgnd._webgl.render(gl,mat_view,mat_proj);}
else{bgnd=rt._cf.background.node;this.setupScene(gl,bgnd);bgnd._webgl.render(gl,mat_view,mat_proj);}
this.stateManager.depthFunc(gl.LEQUAL);this.stateManager.enable(gl.DEPTH_TEST);this.stateManager.enable(gl.CULL_FACE);this.stateManager.blendFuncSeparate(gl.SRC_ALPHA,gl.ONE_MINUS_SRC_ALPHA,gl.ONE,gl.ONE);this.stateManager.enable(gl.BLEND);var slights=viewarea.getLights();var numLights=slights.length;var transform,shape,drawable;var locScene=rt._cf.scene.node;if(!locScene||locScene===scene){n=scene.drawableCollection.length;if(rt._vf.showNormals){this.renderNormals(gl,scene,scene._webgl.normalShader,mat_view,mat_scene);}
else{this.stateManager.unsetProgram();for(i=0;i<n;i++){drawable=scene.drawableCollection.get(i);this.renderShape(drawable,viewarea,slights,numLights,mat_view,mat_scene,mat_light,mat_proj,gl);}}}
else{var env=scene.getEnvironment();var drawableCollectionConfig={viewArea:viewarea,sortTrans:env._vf.sortTrans,viewMatrix:mat_view,projMatrix:mat_proj,sceneMatrix:mat_scene,frustumCulling:false,smallFeatureThreshold:1,context:this,gl:gl};locScene.numberOfNodes=0;locScene.drawableCollection=new x3dom.DrawableCollection(drawableCollectionConfig);locScene.collectDrawableObjects(x3dom.fields.SFMatrix4f.identity(),locScene.drawableCollection,true,false,0,[]);locScene.drawableCollection.sort();n=locScene.drawableCollection.length;if(rt._vf.showNormals){this.renderNormals(gl,locScene,scene._webgl.normalShader,mat_view,mat_scene);}
else{this.stateManager.unsetProgram();for(i=0;i<n;i++){drawable=locScene.drawableCollection.get(i);if(!drawable.shape._vf.render){continue;}
this.renderShape(drawable,viewarea,slights,numLights,mat_view,mat_scene,mat_light,mat_proj,gl);}}}
this.stateManager.disable(gl.BLEND);this.stateManager.disable(gl.DEPTH_TEST);gl.flush();this.stateManager.bindFramebuffer(gl.FRAMEBUFFER,null);if(rt._webgl.fbo.mipMap){gl.bindTexture(gl.TEXTURE_2D,rt._webgl.fbo.tex);gl.generateMipmap(gl.TEXTURE_2D);gl.bindTexture(gl.TEXTURE_2D,null);}
for(i=0;i<m;i++){if(arr[i]!==0){rt._cf.excludeNodes.nodes[i]._vf.render=true;}}};Context.prototype.renderNormals=function(gl,scene,sp,mat_view,mat_scene)
{if(!sp||!scene){return;}
this.stateManager.depthFunc(gl.LEQUAL);this.stateManager.enable(gl.DEPTH_TEST);this.stateManager.enable(gl.CULL_FACE);this.stateManager.disable(gl.BLEND);this.stateManager.useProgram(sp);var bgCenter=x3dom.fields.SFVec3f.NullVector.toGL();var bgSize=x3dom.fields.SFVec3f.OneVector.toGL();for(var i=0,n=scene.drawableCollection.length;i<n;i++)
{var drawable=scene.drawableCollection.get(i);var trafo=drawable.transform;var shape=drawable.shape;var s_gl=shape._webgl;if(!s_gl||!shape||!shape._vf.render){continue;}
var s_geo=shape._cf.geometry.node;var s_msh=s_geo._mesh;var model_view_inv=mat_view.mult(trafo).inverse();sp.normalMatrix=model_view_inv.transpose().toGL();sp.modelViewProjectionMatrix=mat_scene.mult(trafo).toGL();sp.imageGeometry=s_gl.imageGeometry;if(s_gl.coordType!=gl.FLOAT){if(s_gl.popGeometry!=0||(s_msh._numPosComponents==4&&x3dom.Utils.isUnsignedType(s_geo._vf.coordType)))
sp.bgCenter=s_geo.getMin().toGL();else
sp.bgCenter=s_geo._vf.position.toGL();sp.bgSize=s_geo._vf.size.toGL();sp.bgPrecisionMax=s_geo.getPrecisionMax('coordType');}
else{sp.bgCenter=bgCenter;sp.bgSize=bgSize;sp.bgPrecisionMax=1;}
if(s_gl.normalType!=gl.FLOAT){sp.bgPrecisionNorMax=s_geo.getPrecisionMax('normalType');}
else{sp.bgPrecisionNorMax=1;}
if(shape.isSolid()){this.stateManager.enable(gl.CULL_FACE);if(shape.isCCW()){this.stateManager.frontFace(gl.CCW);}
else{this.stateManager.frontFace(gl.CW);}}
else{this.stateManager.disable(gl.CULL_FACE);}
for(var q=0,q_n=s_gl.positions.length;q<q_n;q++){var q6=6*q;var v,v_n,offset;if(s_gl.externalGeometry!=0){var mesh=shape.meshes[q];mesh.bindVertexAttribPointer(gl,sp);mesh.render(gl);}
else
if(!(sp.position!==undefined&&s_gl.buffers[q6+1]&&s_gl.indexes[q]))
continue;if(s_gl.buffers[q6]){gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER,s_gl.buffers[q6]);}
this.setVertexAttribPointerPosition(gl,shape,q6,q);this.setVertexAttribPointerNormal(gl,shape,q6,q);if(s_gl.binaryGeometry>0||s_gl.popGeometry>0){for(v=0,offset=0,v_n=s_geo._vf.vertexCount.length;v<v_n;v++){gl.drawElements(s_gl.primType[v],s_geo._vf.vertexCount[v],s_gl.indexType,x3dom.Utils.getByteAwareOffset(offset,s_gl.indexType,gl));offset+=s_geo._vf.vertexCount[v];}}
else if(s_gl.binaryGeometry<0||s_gl.popGeometry<0||s_gl.imageGeometry){for(v=0,offset=0,v_n=s_geo._vf.vertexCount.length;v<v_n;v++){gl.drawArrays(s_gl.primType[v],offset,s_geo._vf.vertexCount[v]);offset+=s_geo._vf.vertexCount[v];}}
else if(s_geo.hasIndexOffset()){var indOff=shape.tessellationProperties();for(v=0,v_n=indOff.length;v<v_n;v++){gl.drawElements(s_gl.primType,indOff[v].count,s_gl.indexType,indOff[v].offset*x3dom.Utils.getOffsetMultiplier(s_gl.indexType,gl));}}
else if(s_gl.indexes[q].length==0){gl.drawArrays(s_gl.primType,0,s_gl.positions[q].length/3);}
else{gl.drawElements(s_gl.primType,s_gl.indexes[q].length,s_gl.indexType,0);}
gl.disableVertexAttribArray(sp.position);if(sp.normal!==undefined){gl.disableVertexAttribArray(sp.normal);}}}};Context.prototype.shutdown=function(viewarea){var gl=this.ctx3d;var scene=viewarea._scene;if(gl==null||!scene){return;}
var bgnd=scene.getBackground();if(bgnd._webgl.position!==undefined){gl.deleteBuffer(bgnd._webgl.buffers[1]);gl.deleteBuffer(bgnd._webgl.buffers[0]);}
var fgnd=scene._fgnd;if(fgnd._webgl.position!==undefined){gl.deleteBuffer(fgnd._webgl.buffers[1]);gl.deleteBuffer(fgnd._webgl.buffers[0]);}
var n=scene.drawableCollection?scene.drawableCollection.length:0;for(var i=0;i<n;i++){var shape=scene.drawableCollection.get(i).shape;if(shape._cleanupGLObjects)
shape._cleanupGLObjects(true);}
this.cache.Release(gl);};Context.prototype.renderShadows=function(gl,viewarea,shadowedLights,wctolc,lMatrices,mat_view,mat_proj,mat_scene)
{var scene=viewarea._scene;var texLimit=x3dom.caps.MAX_TEXTURE_IMAGE_UNITS;if(texLimit<7)
return;var texUnits=1;var renderSplit=[0];var shadowMaps,numShadowMaps;var i,j,k;for(i=0;i<shadowedLights.length;i++)
{var filterSize=shadowedLights[i]._vf.shadowFilterSize;shadowMaps=scene._webgl.fboShadow[i];numShadowMaps=shadowMaps.length;for(j=0;j<numShadowMaps;j++){this.blurTex(gl,scene,shadowMaps[j],filterSize);}
texUnits+=6;if(texUnits>texLimit){renderSplit[renderSplit.length]=i;texUnits=7;}}
renderSplit[renderSplit.length]=shadowedLights.length;var n=renderSplit.length-1;var mat_proj_inv=mat_proj.inverse();var mat_scene_inv=mat_scene.inverse();this.stateManager.enable(gl.BLEND);this.stateManager.blendFunc(gl.DST_COLOR,gl.ZERO);for(var s=0;s<n;s++)
{var startIndex=renderSplit[s];var endIndex=renderSplit[s+1];var currentLights=[];for(k=startIndex;k<endIndex;k++)
currentLights[currentLights.length]=shadowedLights[k];var sp=this.cache.getShadowRenderingShader(gl,currentLights);this.stateManager.useProgram(sp);gl.bindBuffer(gl.ARRAY_BUFFER,scene._webgl.ppBuffer);gl.vertexAttribPointer(sp.position,2,gl.FLOAT,false,0,0);gl.enableVertexAttribArray(sp.position);sp.sceneMap=0;gl.activeTexture(gl.TEXTURE0);gl.bindTexture(gl.TEXTURE_2D,scene._webgl.fboScene.tex);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MAG_FILTER,gl.LINEAR);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MIN_FILTER,gl.LINEAR);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_S,gl.CLAMP_TO_EDGE);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_T,gl.CLAMP_TO_EDGE);sp.inverseProj=mat_proj_inv.toGL();sp.inverseViewProj=mat_scene_inv.toGL();var mat_light;var lightMatrix;var shadowIndex=0;for(var p=0,pn=currentLights.length;p<pn;p++){lightMatrix=lMatrices[p+startIndex];mat_light=wctolc[p+startIndex];shadowMaps=scene._webgl.fboShadow[p+startIndex];numShadowMaps=mat_light.length;for(i=0;i<numShadowMaps;i++){gl.activeTexture(gl.TEXTURE1+shadowIndex);gl.bindTexture(gl.TEXTURE_2D,shadowMaps[i].tex);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MAG_FILTER,gl.LINEAR);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MIN_FILTER,gl.LINEAR);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_S,gl.CLAMP_TO_EDGE);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_T,gl.CLAMP_TO_EDGE);sp['light'+p+'_'+i+'_ShadowMap']=shadowIndex+1;sp['light'+p+'_'+i+'_Matrix']=mat_light[i].toGL();shadowIndex++;}
sp['light'+p+'_ViewMatrix']=lightMatrix.toGL();if(!x3dom.isa(currentLights[p],x3dom.nodeTypes.PointLight)){for(j=0;j<numShadowMaps;j++){var numCascades=Math.max(1,Math.min(currentLights[p]._vf.shadowCascades,6));var splitFactor=Math.max(0,Math.min(currentLights[p]._vf.shadowSplitFactor,1));var splitOffset=Math.max(0,Math.min(currentLights[p]._vf.shadowSplitOffset,1));var splitDepths=viewarea.getShadowSplitDepths(numCascades,splitFactor,splitOffset,false,mat_proj);sp['light'+p+'_'+j+'_Split']=splitDepths[j+1];}}
var light_transform=mat_view.mult(currentLights[p].getCurrentTransform());if(x3dom.isa(currentLights[p],x3dom.nodeTypes.DirectionalLight))
{sp['light'+p+'_Type']=0.0;sp['light'+p+'_On']=(currentLights[p]._vf.on)?1.0:0.0;sp['light'+p+'_Direction']=light_transform.multMatrixVec(currentLights[p]._vf.direction).toGL();sp['light'+p+'_Attenuation']=[1.0,1.0,1.0];sp['light'+p+'_Location']=[1.0,1.0,1.0];sp['light'+p+'_Radius']=0.0;sp['light'+p+'_BeamWidth']=0.0;sp['light'+p+'_CutOffAngle']=0.0;sp['light'+p+'_ShadowIntensity']=currentLights[p]._vf.shadowIntensity;sp['light'+p+'_ShadowCascades']=currentLights[p]._vf.shadowCascades;sp['light'+p+'_ShadowOffset']=Math.max(0.0,Math.min(1.0,currentLights[p]._vf.shadowOffset));}
else if(x3dom.isa(currentLights[p],x3dom.nodeTypes.PointLight))
{sp['light'+p+'_Type']=1.0;sp['light'+p+'_On']=(currentLights[p]._vf.on)?1.0:0.0;sp['light'+p+'_Direction']=[1.0,1.0,1.0];sp['light'+p+'_Attenuation']=currentLights[p]._vf.attenuation.toGL();sp['light'+p+'_Location']=light_transform.multMatrixPnt(currentLights[p]._vf.location).toGL();sp['light'+p+'_Radius']=currentLights[p]._vf.radius;sp['light'+p+'_BeamWidth']=0.0;sp['light'+p+'_CutOffAngle']=0.0;sp['light'+p+'_ShadowIntensity']=currentLights[p]._vf.shadowIntensity;sp['light'+p+'_ShadowOffset']=Math.max(0.0,Math.min(1.0,currentLights[p]._vf.shadowOffset));}
else if(x3dom.isa(currentLights[p],x3dom.nodeTypes.SpotLight))
{sp['light'+p+'_Type']=2.0;sp['light'+p+'_On']=(currentLights[p]._vf.on)?1.0:0.0;sp['light'+p+'_Direction']=light_transform.multMatrixVec(currentLights[p]._vf.direction).toGL();sp['light'+p+'_Attenuation']=currentLights[p]._vf.attenuation.toGL();sp['light'+p+'_Location']=light_transform.multMatrixPnt(currentLights[p]._vf.location).toGL();sp['light'+p+'_Radius']=currentLights[p]._vf.radius;sp['light'+p+'_BeamWidth']=currentLights[p]._vf.beamWidth;sp['light'+p+'_CutOffAngle']=currentLights[p]._vf.cutOffAngle;sp['light'+p+'_ShadowIntensity']=currentLights[p]._vf.shadowIntensity;sp['light'+p+'_ShadowCascades']=currentLights[p]._vf.shadowCascades;sp['light'+p+'_ShadowOffset']=Math.max(0.0,Math.min(1.0,currentLights[p]._vf.shadowOffset));}}
gl.drawArrays(gl.TRIANGLES,0,6);var nk=shadowIndex+1;for(k=0;k<nk;k++){gl.activeTexture(gl.TEXTURE0+k);gl.bindTexture(gl.TEXTURE_2D,null);}
gl.disableVertexAttribArray(sp.position);}
this.stateManager.blendFunc(gl.SRC_ALPHA,gl.ONE_MINUS_SRC_ALPHA);};Context.prototype.blurTex=function(gl,scene,targetFbo,filterSize)
{if(filterSize<=0)
return;else if(filterSize<5)
filterSize=3;else if(filterSize<7)
filterSize=5;else
filterSize=7;var width=targetFbo.width;var height=targetFbo.height;var fboBlur=null;for(var i=0,n=scene._webgl.fboBlur.length;i<n;i++)
if(height==scene._webgl.fboBlur[i].height){fboBlur=scene._webgl.fboBlur[i];break;}
this.stateManager.bindFramebuffer(gl.FRAMEBUFFER,fboBlur.fbo);this.stateManager.viewport(0,0,width,height);this.stateManager.enable(gl.BLEND);this.stateManager.blendFunc(gl.ONE,gl.ZERO);this.stateManager.disable(gl.CULL_FACE);this.stateManager.disable(gl.DEPTH_TEST);gl.clearColor(1.0,1.0,1.0,0.0);gl.clearDepth(1.0);gl.clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER_BIT);var sp=this.cache.getShader(gl,x3dom.shader.BLUR);this.stateManager.useProgram(sp);gl.bindBuffer(gl.ARRAY_BUFFER,scene._webgl.ppBuffer);gl.vertexAttribPointer(sp.position,2,gl.FLOAT,false,0,0);gl.enableVertexAttribArray(sp.position);sp.pixelSizeHor=1.0/width;sp.pixelSizeVert=1.0/height;sp.filterSize=filterSize;sp.horizontal=true;sp.texture=0;gl.activeTexture(gl.TEXTURE0);gl.bindTexture(gl.TEXTURE_2D,targetFbo.tex);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MAG_FILTER,gl.LINEAR);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MIN_FILTER,gl.LINEAR);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_S,gl.CLAMP_TO_EDGE);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_T,gl.CLAMP_TO_EDGE);gl.drawArrays(gl.TRIANGLES,0,6);this.stateManager.bindFramebuffer(gl.FRAMEBUFFER,targetFbo.fbo);gl.clearColor(1.0,1.0,1.0,0.0);gl.clearDepth(1.0);gl.clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER_BIT);sp.horizontal=false;gl.activeTexture(gl.TEXTURE0);gl.bindTexture(gl.TEXTURE_2D,fboBlur.tex);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MAG_FILTER,gl.LINEAR);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MIN_FILTER,gl.LINEAR);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_S,gl.CLAMP_TO_EDGE);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_T,gl.CLAMP_TO_EDGE);gl.drawArrays(gl.TRIANGLES,0,6);gl.activeTexture(gl.TEXTURE0);gl.bindTexture(gl.TEXTURE_2D,null);gl.disableVertexAttribArray(sp.position);gl.flush();this.stateManager.blendFunc(gl.SRC_ALPHA,gl.ONE_MINUS_SRC_ALPHA);this.stateManager.bindFramebuffer(gl.FRAMEBUFFER,null);this.stateManager.viewport(0,0,this.canvas.width,this.canvas.height);};Context.prototype.setVertexAttribPointerPosition=function(gl,shape,q6,q)
{var sp=shape._webgl.shader;if(sp.position!==undefined&&shape._webgl.buffers[q6+1])
{var s_geo=shape._cf.geometry.node;gl.bindBuffer(gl.ARRAY_BUFFER,shape._webgl.buffers[q6+1]);gl.vertexAttribPointer(sp.position,s_geo._mesh._numPosComponents,shape._webgl.coordType,false,shape._coordStrideOffset[0],shape._coordStrideOffset[1]);gl.enableVertexAttribArray(sp.position);}};Context.prototype.setVertexAttribPointerNormal=function(gl,shape,q6,q)
{var sp=shape._webgl.shader;if(sp.normal!==undefined&&shape._webgl.buffers[q6+2])
{var s_geo=shape._cf.geometry.node;gl.bindBuffer(gl.ARRAY_BUFFER,shape._webgl.buffers[q6+2]);gl.vertexAttribPointer(sp.normal,s_geo._mesh._numNormComponents,shape._webgl.normalType,false,shape._normalStrideOffset[0],shape._normalStrideOffset[1]);gl.enableVertexAttribArray(sp.normal);}};Context.prototype.setVertexAttribPointerTexCoord=function(gl,shape,q6,q)
{var sp=shape._webgl.shader;if(sp.texcoord!==undefined&&shape._webgl.buffers[q6+3])
{var s_geo=shape._cf.geometry.node;gl.bindBuffer(gl.ARRAY_BUFFER,shape._webgl.buffers[q6+3]);gl.vertexAttribPointer(sp.texcoord,s_geo._mesh._numTexComponents,shape._webgl.texCoordType,false,shape._texCoordStrideOffset[0],shape._texCoordStrideOffset[1]);gl.enableVertexAttribArray(sp.texcoord);}};Context.prototype.setVertexAttribPointerColor=function(gl,shape,q6,q)
{var sp=shape._webgl.shader;if(sp.color!==undefined&&shape._webgl.buffers[q6+4])
{var s_geo=shape._cf.geometry.node;gl.bindBuffer(gl.ARRAY_BUFFER,shape._webgl.buffers[q6+4]);gl.vertexAttribPointer(sp.color,s_geo._mesh._numColComponents,shape._webgl.colorType,false,shape._colorStrideOffset[0],shape._colorStrideOffset[1]);gl.enableVertexAttribArray(sp.color);}};return setupContext;})();x3dom.NodeNameSpace=function(name,document){this.name=name;this.doc=document;this.baseURL="";this.defMap={};this.parent=null;this.childSpaces=[];};x3dom.NodeNameSpace.prototype.addNode=function(node,name){this.defMap[name]=node;node._nameSpace=this;};x3dom.NodeNameSpace.prototype.removeNode=function(name){var node=name?this.defMap[name]:null;if(node){delete this.defMap[name];node._nameSpace=null;}};x3dom.NodeNameSpace.prototype.getNamedNode=function(name){return this.defMap[name];};x3dom.NodeNameSpace.prototype.getNamedElement=function(name){var node=this.defMap[name];return(node?node._xmlNode:null);};x3dom.NodeNameSpace.prototype.addSpace=function(space){this.childSpaces.push(space);space.parent=this;};x3dom.NodeNameSpace.prototype.removeSpace=function(space){space.parent=null;for(var it=0;it<this.childSpaces.length;it++){if(this.childSpaces[it]==space){this.childSpaces.splice(it,1);}}};x3dom.NodeNameSpace.prototype.setBaseURL=function(url){var i=url.lastIndexOf("/");this.baseURL=(i>=0)?url.substr(0,i+1):"";x3dom.debug.logInfo("setBaseURL: "+this.baseURL);};x3dom.NodeNameSpace.prototype.getURL=function(url){if(url===undefined||!url.length){return"";}
else{return((url[0]==='/')||(url.indexOf(":")>=0))?url:(this.baseURL+url);}};x3dom.hasElementAttribute=function(attrName)
{var ok=this.__hasAttribute(attrName);if(!ok&&attrName){ok=this.__hasAttribute(attrName.toLowerCase());}
return ok;};x3dom.getElementAttribute=function(attrName)
{var attrib=this.__getAttribute(attrName);if(!attrib&&attrib!=""&&attrName){attrib=this.__getAttribute(attrName.toLowerCase());}
if(attrib||!this._x3domNode){return attrib;}
else{return this._x3domNode._vf[attrName];}};x3dom.setElementAttribute=function(attrName,newVal)
{this.__setAttribute(attrName,newVal);var x3dNode=this._x3domNode;if(x3dNode){x3dNode.updateField(attrName,newVal);x3dNode._nameSpace.doc.needRender=true;}};x3dom.getFieldValue=function(fieldName)
{var x3dNode=this._x3domNode;if(x3dNode&&(x3dNode._vf[fieldName]!==undefined)){var fieldValue=x3dNode._vf[fieldName];if(fieldValue instanceof Object&&'copy'in fieldValue)
{return x3dNode._vf[fieldName].copy();}
else
{return x3dNode._vf[fieldName];}}
return null;};x3dom.setFieldValue=function(fieldName,fieldvalue){var x3dNode=this._x3domNode;if(x3dNode&&(x3dNode._vf[fieldName]!==undefined)){if(fieldvalue instanceof Object&&'copy'in fieldvalue)
{x3dNode._vf[fieldName]=fieldvalue.copy();}
else
x3dNode._vf[fieldName]=fieldvalue;x3dNode.fieldChanged(fieldName);x3dNode._nameSpace.doc.needRender=true;}};x3dom.requestFieldRef=function(fieldName)
{var x3dNode=this._x3domNode;if(x3dNode&&x3dNode._vf[fieldName])
{return x3dNode._vf[fieldName];}
return null;};x3dom.releaseFieldRef=function(fieldName)
{var x3dNode=this._x3domNode;if(x3dNode&&x3dNode._vf[fieldName])
{x3dNode.fieldChanged(fieldName);x3dNode._nameSpace.doc.needRender=true;}};x3dom.NodeNameSpace.prototype.setupTree=function(domNode,parent){var n=null;parent=parent||null;if(x3dom.isX3DElement(domNode)){if(domNode._x3domNode){x3dom.debug.logWarning('Tree is already initialized');return null;}
if((domNode.tagName!==undefined)&&(!domNode.__addEventListener)&&(!domNode.__removeEventListener))
{domNode.__addEventListener=domNode.addEventListener;domNode.addEventListener=function(type,func,phase){if(!this._x3domNode._listeners[type]){this._x3domNode._listeners[type]=[];}
this._x3domNode._listeners[type].push(func);this.__addEventListener(type,func,phase);};domNode.__removeEventListener=domNode.removeEventListener;domNode.removeEventListener=function(type,func,phase){var list=this._x3domNode._listeners[type];if(list){for(var it=0;it<list.length;it++){if(list[it]==func){list.splice(it,1);}}}
this.__removeEventListener(type,func,phase);};}
if(domNode.hasAttribute('USE')||domNode.hasAttribute('use'))
{if(!domNode.hasAttribute('USE')){domNode.setAttribute('USE',domNode.getAttribute('use'));}
n=this.defMap[domNode.getAttribute('USE')];if(!n){var nsName=domNode.getAttribute('USE').split('__');if(nsName.length>=2){var otherNS=this;while(otherNS){if(otherNS.name==nsName[0])
n=otherNS.defMap[nsName[1]];if(n)
otherNS=null;else
otherNS=otherNS.parent;}
if(!n){n=null;x3dom.debug.logWarning('Could not USE: '+domNode.getAttribute('USE'));}}}
if(n){domNode._x3domNode=n;}
return n;}
else{if(domNode.localName.toLowerCase()==='route'){var route=domNode;var fnAtt=route.getAttribute('fromNode')||route.getAttribute('fromnode');var tnAtt=route.getAttribute('toNode')||route.getAttribute('tonode');var fromNode=this.defMap[fnAtt];var toNode=this.defMap[tnAtt];if(!(fromNode&&toNode)){x3dom.debug.logWarning("Broken route - can't find all DEFs for "+fnAtt+" -> "+tnAtt);}
else{fnAtt=route.getAttribute('fromField')||route.getAttribute('fromfield');tnAtt=route.getAttribute('toField')||route.getAttribute('tofield');fromNode.setupRoute(fnAtt,toNode,tnAtt);route._nodeNameSpace=this;}
return null;}
domNode.requestFieldRef=x3dom.requestFieldRef;domNode.releaseFieldRef=x3dom.releaseFieldRef;domNode.getFieldValue=x3dom.getFieldValue;domNode.setFieldValue=x3dom.setFieldValue;var nodeType=x3dom.nodeTypesLC[domNode.localName.toLowerCase()];if(nodeType===undefined){x3dom.debug.logWarning("Unrecognised X3D element &lt;"+domNode.localName+"&gt;.");}
else{if((x3dom.userAgentFeature.supportsDOMAttrModified===false)&&(domNode instanceof Element)){if(domNode.setAttribute&&!domNode.__setAttribute){domNode.__setAttribute=domNode.setAttribute;domNode.setAttribute=x3dom.setElementAttribute;}
if(domNode.getAttribute&&!domNode.__getAttribute){domNode.__getAttribute=domNode.getAttribute;domNode.getAttribute=x3dom.getElementAttribute;}
if(domNode.hasAttribute&&!domNode.__hasAttribute){domNode.__hasAttribute=domNode.hasAttribute;domNode.hasAttribute=x3dom.hasElementAttribute;}}
var ctx={doc:this.doc,xmlNode:domNode,nameSpace:this};n=new nodeType(ctx);if(domNode.hasAttribute('DEF')){n._DEF=domNode.getAttribute('DEF');this.defMap[n._DEF]=n;}
else{if(domNode.hasAttribute('id')){n._DEF=domNode.getAttribute('id');this.defMap[n._DEF]=n;}}
if(domNode.highlight===undefined)
{domNode.highlight=function(enable,colorStr){var color=x3dom.fields.SFColor.parse(colorStr);this._x3domNode.highlight(enable,color);this._x3domNode._nameSpace.doc.needRender=true;};}
n._xmlNode=domNode;domNode._x3domNode=n;var that=this;Array.forEach(domNode.childNodes,function(childDomNode){var c=that.setupTree(childDomNode,n);if(c){n.addChild(c,childDomNode.getAttribute("containerField"));}});n.nodeChanged();return n;}}}
else if(domNode.localName)
{if(parent&&domNode.localName.toLowerCase()=="x3dommetagroup")
{Array.forEach(domNode.childNodes,function(childDomNode){var c=this.setupTree(childDomNode,parent);if(c){parent.addChild(c,childDomNode.getAttribute("containerField"));}}.bind(this));}
else
{x3dom.debug.logWarning("Unrecognised X3D element &lt;"+domNode.localName+"&gt;.");n=null;}}
return n;};x3dom.registerNodeType("X3DNode","Core",defineClass(null,function(ctx){this._xmlNode=null;this._DEF=null;this._nameSpace=(ctx&&ctx.nameSpace)?ctx.nameSpace:null;this._vf={};this._vfFieldTypes={};this._cf={};this._cfFieldTypes={};this._fieldWatchers={};this._routes={};this._listeners={};this._parentNodes=[];this._childNodes=[];this.addField_SFNode('metadata',x3dom.nodeTypes.X3DMetadataObject);},{type:function(){return this.constructor;},typeName:function(){return this.constructor._typeName;},addChild:function(node,containerFieldName){if(node){var field=null;if(containerFieldName){field=this._cf[containerFieldName];}
else{for(var fieldName in this._cf){if(this._cf.hasOwnProperty(fieldName)){var testField=this._cf[fieldName];if(x3dom.isa(node,testField.type)){field=testField;break;}}}}
if(field&&field.addLink(node)){node._parentNodes.push(this);this._childNodes.push(node);node.parentAdded(this);return true;}}
return false;},removeChild:function(node){if(node){for(var fieldName in this._cf){if(this._cf.hasOwnProperty(fieldName)){var field=this._cf[fieldName];if(field.rmLink(node)){for(var i=node._parentNodes.length-1;i>=0;i--){if(node._parentNodes[i]===this){node._parentNodes.splice(i,1);node.parentRemoved(this);}}
for(var j=this._childNodes.length-1;j>=0;j--){if(this._childNodes[j]===node){node.onRemove();this._childNodes.splice(j,1);return true;}}}}}}
return false;},onRemove:function(){},parentAdded:function(parent){},parentRemoved:function(parent){for(var i=0,n=this._childNodes.length;i<n;i++){if(this._childNodes[i]){this._childNodes[i].parentRemoved(this);}}},getCurrentTransform:function(){if(this._parentNodes.length>=1){return this.transformMatrix(this._parentNodes[0].getCurrentTransform());}
else{return x3dom.fields.SFMatrix4f.identity();}},transformMatrix:function(transform){return transform;},getVolume:function(){return null;},invalidateVolume:function(){},invalidateCache:function(){},volumeValid:function(){return false;},collectDrawableObjects:function(transform,drawableCollection,singlePath,invalidateCache,planeMask,clipPlanes){},highlight:function(enable,color)
{if(this._vf.hasOwnProperty("diffuseColor"))
{if(enable){if(this._actDiffuseColor===undefined){this._actDiffuseColor=new x3dom.fields.SFColor();this._highlightOn=false;}
if(!this._highlightOn){this._actDiffuseColor.setValues(this._vf.diffuseColor);this._highlightOn=true;}
this._vf.diffuseColor.setValues(color);}
else{if(this._actDiffuseColor!==undefined){this._vf.diffuseColor.setValues(this._actDiffuseColor);this._highlightOn=false;delete this._actDiffuseColor;}}}
for(var i=0,n=this._childNodes.length;i<n;i++)
{if(this._childNodes[i])
this._childNodes[i].highlight(enable,color);}},findX3DDoc:function(){return this._nameSpace.doc;},doIntersect:function(line){var isect=false;for(var i=0;i<this._childNodes.length;i++){if(this._childNodes[i]){isect=this._childNodes[i].doIntersect(line)||isect;}}
return isect;},postMessage:function(field,msg){this._vf[field]=msg;var listeners=this._fieldWatchers[field];var that=this;if(listeners){Array.forEach(listeners,function(l){l.call(that,msg);});}
var eventObject={target:that._xmlNode,type:"outputchange",fieldName:field,value:msg};this.callEvtHandler("onoutputchange",eventObject);},updateField:function(field,msg){var f=this._vf[field];if(f===undefined){for(var key in this._vf){if(key.toLowerCase()==field){field=key;f=this._vf[field];break;}}
var pre="set_";if(f===undefined&&field.indexOf(pre)==0){var fieldName=field.substr(pre.length,field.length-1);if(this._vf[fieldName]!==undefined){field=fieldName;f=this._vf[field];}}
if(f===undefined){f=null;this._vf[field]=f;}}
if(f!==null){try{this._vf[field].setValueByStr(msg);}
catch(exc1){try{switch((typeof(this._vf[field])).toString()){case"number":if(typeof(msg)=="number")
this._vf[field]=msg;else
this._vf[field]=+msg;break;case"boolean":if(typeof(msg)=="boolean")
this._vf[field]=msg;else
this._vf[field]=(msg.toLowerCase()=="true");break;case"string":this._vf[field]=msg;break;}}
catch(exc2){x3dom.debug.logError("updateField: setValueByStr() NYI for "+typeof(f));}}
this.fieldChanged(field);}},setupRoute:function(fromField,toNode,toField){var pos;var fieldName;var pre="set_",post="_changed";if(!this._vf[fromField]){pos=fromField.indexOf(pre);if(pos===0){fieldName=fromField.substr(pre.length,fromField.length-1);if(this._vf[fieldName]){fromField=fieldName;}}else{pos=fromField.indexOf(post);if(pos>0){fieldName=fromField.substr(0,fromField.length-post.length);if(this._vf[fieldName]){fromField=fieldName;}}}}
if(!toNode._vf[toField]){pos=toField.indexOf(pre);if(pos===0){fieldName=toField.substr(pre.length,toField.length-1);if(toNode._vf[fieldName]){toField=fieldName;}}
else{pos=toField.indexOf(post);if(pos>0){fieldName=toField.substr(0,toField.length-post.length);if(toNode._vf[fieldName]){toField=fieldName;}}}}
var where=this._DEF+"&"+fromField+"&"+toNode._DEF+"&"+toField;if(!this._routes[where]){if(!this._fieldWatchers[fromField]){this._fieldWatchers[fromField]=[];}
this._fieldWatchers[fromField].push(function(msg){toNode.postMessage(toField,msg);});if(!toNode._fieldWatchers[toField]){toNode._fieldWatchers[toField]=[];}
toNode._fieldWatchers[toField].push(function(msg){toNode._vf[toField]=msg;toNode.fieldChanged(toField);});this._routes[where]={from:this._fieldWatchers[fromField].length-1,to:toNode._fieldWatchers[toField].length-1};}},removeRoute:function(fromField,toNode,toField){var pos;var fieldName;var pre="set_",post="_changed";if(!this._vf[fromField]){pos=fromField.indexOf(pre);if(pos===0){fieldName=fromField.substr(pre.length,fromField.length-1);if(this._vf[fieldName]){fromField=fieldName;}}else{pos=fromField.indexOf(post);if(pos>0){fieldName=fromField.substr(0,fromField.length-post.length);if(this._vf[fieldName]){fromField=fieldName;}}}}
if(!toNode._vf[toField]){pos=toField.indexOf(pre);if(pos===0){fieldName=toField.substr(pre.length,toField.length-1);if(toNode._vf[fieldName]){toField=fieldName;}}
else{pos=toField.indexOf(post);if(pos>0){fieldName=toField.substr(0,toField.length-post.length);if(toNode._vf[fieldName]){toField=fieldName;}}}}
var where=this._DEF+"&"+fromField+"&"+toNode._DEF+"&"+toField;if(this._routes[where]){this._fieldWatchers[fromField].splice(this._routes[where].from,1);toNode._fieldWatchers[toField].splice(this._routes[where].to,1);delete this._routes[where];}},fieldChanged:function(fieldName){},nodeChanged:function(){},callEvtHandler:function(eventType,event){var node=this;if(!node._xmlNode){return event.cancelBubble;}
try{var attrib=node._xmlNode[eventType];event.target=node._xmlNode;if(typeof(attrib)==="function"){attrib.call(node._xmlNode,event);}
else{var funcStr=node._xmlNode.getAttribute(eventType);var func=new Function('event',funcStr);func.call(node._xmlNode,event);}
var list=node._listeners[event.type];if(list){for(var it=0;it<list.length;it++){list[it].call(node._xmlNode,event);}}}
catch(ex){x3dom.debug.logException(ex);}
return event.cancelBubble;},initSetter:function(xmlNode,name){if(!xmlNode||!name)
return;var nameLC=name.toLowerCase();if(xmlNode.__defineSetter__&&xmlNode.__defineGetter__){xmlNode.__defineSetter__(name,function(value){xmlNode.setAttribute(name,value);});xmlNode.__defineGetter__(name,function(){return xmlNode.getAttribute(name);});if(nameLC!=name){xmlNode.__defineSetter__(nameLC,function(value){xmlNode.setAttribute(name,value);});xmlNode.__defineGetter__(nameLC,function(){return xmlNode.getAttribute(name);});}}
else{Object.defineProperty(xmlNode,name,{set:function(value){xmlNode.setAttribute(name,value);},get:function(){return xmlNode.getAttribute(name);},configurable:true,enumerable:true});}
if(this._vf[name]&&!xmlNode.attributes[name]&&!xmlNode.attributes[name.toLowerCase()]){var str="";try{if(this._vf[name].toGL)
str=this._vf[name].toGL().toString();else
str=this._vf[name].toString();}
catch(e){str=this._vf[name].toString();}
if(!str){str="";}
xmlNode.setAttribute(name,str);}},addField_SFInt32:function(ctx,name,n){this._vf[name]=ctx&&ctx.xmlNode&&ctx.xmlNode.hasAttribute(name)?parseInt(ctx.xmlNode.getAttribute(name),10):n;if(ctx&&ctx.xmlNode){this.initSetter(ctx.xmlNode,name);}
this._vfFieldTypes[name]="SFInt32";},addField_SFFloat:function(ctx,name,n){this._vf[name]=ctx&&ctx.xmlNode&&ctx.xmlNode.hasAttribute(name)?+ctx.xmlNode.getAttribute(name):n;if(ctx&&ctx.xmlNode){this.initSetter(ctx.xmlNode,name);}
this._vfFieldTypes[name]="SFFloat";},addField_SFDouble:function(ctx,name,n){this._vf[name]=ctx&&ctx.xmlNode&&ctx.xmlNode.hasAttribute(name)?+ctx.xmlNode.getAttribute(name):n;if(ctx&&ctx.xmlNode){this.initSetter(ctx.xmlNode,name);}
this._vfFieldTypes[name]="SFDouble";},addField_SFTime:function(ctx,name,n){this._vf[name]=ctx&&ctx.xmlNode&&ctx.xmlNode.hasAttribute(name)?+ctx.xmlNode.getAttribute(name):n;if(ctx&&ctx.xmlNode){this.initSetter(ctx.xmlNode,name);}
this._vfFieldTypes[name]="SFTime";},addField_SFBool:function(ctx,name,n){this._vf[name]=ctx&&ctx.xmlNode&&ctx.xmlNode.hasAttribute(name)?ctx.xmlNode.getAttribute(name).toLowerCase()==="true":n;if(ctx&&ctx.xmlNode){this.initSetter(ctx.xmlNode,name);}
this._vfFieldTypes[name]="SFBool";},addField_SFString:function(ctx,name,n){this._vf[name]=ctx&&ctx.xmlNode&&ctx.xmlNode.hasAttribute(name)?ctx.xmlNode.getAttribute(name):n;if(ctx&&ctx.xmlNode){this.initSetter(ctx.xmlNode,name);}
this._vfFieldTypes[name]="SFString";},addField_SFColor:function(ctx,name,r,g,b){this._vf[name]=ctx&&ctx.xmlNode&&ctx.xmlNode.hasAttribute(name)?x3dom.fields.SFColor.parse(ctx.xmlNode.getAttribute(name)):new x3dom.fields.SFColor(r,g,b);if(ctx&&ctx.xmlNode){this.initSetter(ctx.xmlNode,name);}
this._vfFieldTypes[name]="SFColor";},addField_SFColorRGBA:function(ctx,name,r,g,b,a){this._vf[name]=ctx&&ctx.xmlNode&&ctx.xmlNode.hasAttribute(name)?x3dom.fields.SFColorRGBA.parse(ctx.xmlNode.getAttribute(name)):new x3dom.fields.SFColorRGBA(r,g,b,a);if(ctx&&ctx.xmlNode){this.initSetter(ctx.xmlNode,name);}
this._vfFieldTypes[name]="SFColorRGBA";},addField_SFVec2f:function(ctx,name,x,y){this._vf[name]=ctx&&ctx.xmlNode&&ctx.xmlNode.hasAttribute(name)?x3dom.fields.SFVec2f.parse(ctx.xmlNode.getAttribute(name)):new x3dom.fields.SFVec2f(x,y);if(ctx&&ctx.xmlNode){this.initSetter(ctx.xmlNode,name);}
this._vfFieldTypes[name]="SFVec2f";},addField_SFVec3f:function(ctx,name,x,y,z){this._vf[name]=ctx&&ctx.xmlNode&&ctx.xmlNode.hasAttribute(name)?x3dom.fields.SFVec3f.parse(ctx.xmlNode.getAttribute(name)):new x3dom.fields.SFVec3f(x,y,z);if(ctx&&ctx.xmlNode){this.initSetter(ctx.xmlNode,name);}
this._vfFieldTypes[name]="SFVec3f";},addField_SFVec4f:function(ctx,name,x,y,z,w){this._vf[name]=ctx&&ctx.xmlNode&&ctx.xmlNode.hasAttribute(name)?x3dom.fields.SFVec4f.parse(ctx.xmlNode.getAttribute(name)):new x3dom.fields.SFVec4f(x,y,z,w);if(ctx&&ctx.xmlNode){this.initSetter(ctx.xmlNode,name);}
this._vfFieldTypes[name]="SFVec4f";},addField_SFVec3d:function(ctx,name,x,y,z){this.addField_SFVec3f(ctx,name,x,y,z);this._vfFieldTypes[name]="SFVec3d";},addField_SFRotation:function(ctx,name,x,y,z,a){this._vf[name]=ctx&&ctx.xmlNode&&ctx.xmlNode.hasAttribute(name)?x3dom.fields.Quaternion.parseAxisAngle(ctx.xmlNode.getAttribute(name)):x3dom.fields.Quaternion.axisAngle(new x3dom.fields.SFVec3f(x,y,z),a);if(ctx&&ctx.xmlNode){this.initSetter(ctx.xmlNode,name);}
this._vfFieldTypes[name]="SFRotation";},addField_SFMatrix4f:function(ctx,name,_00,_01,_02,_03,_10,_11,_12,_13,_20,_21,_22,_23,_30,_31,_32,_33){this._vf[name]=ctx&&ctx.xmlNode&&ctx.xmlNode.hasAttribute(name)?x3dom.fields.SFMatrix4f.parse(ctx.xmlNode.getAttribute(name)):new x3dom.fields.SFMatrix4f(_00,_01,_02,_03,_10,_11,_12,_13,_20,_21,_22,_23,_30,_31,_32,_33);if(ctx&&ctx.xmlNode){this.initSetter(ctx.xmlNode,name);}
this._vfFieldTypes[name]="SFMatrix4f";},addField_SFImage:function(ctx,name,def){this._vf[name]=ctx&&ctx.xmlNode&&ctx.xmlNode.hasAttribute(name)?x3dom.fields.SFImage.parse(ctx.xmlNode.getAttribute(name)):new x3dom.fields.SFImage(def);if(ctx&&ctx.xmlNode){this.initSetter(ctx.xmlNode,name);}
this._vfFieldTypes[name]="SFImage";},addField_MFString:function(ctx,name,def){this._vf[name]=ctx&&ctx.xmlNode&&ctx.xmlNode.hasAttribute(name)?x3dom.fields.MFString.parse(ctx.xmlNode.getAttribute(name)):new x3dom.fields.MFString(def);if(ctx&&ctx.xmlNode){this.initSetter(ctx.xmlNode,name);}
this._vfFieldTypes[name]="MFString";},addField_MFBoolean:function(ctx,name,def){this._vf[name]=ctx&&ctx.xmlNode&&ctx.xmlNode.hasAttribute(name)?x3dom.fields.MFBoolean.parse(ctx.xmlNode.getAttribute(name)):new x3dom.fields.MFBoolean(def);if(ctx&&ctx.xmlNode){this.initSetter(ctx.xmlNode,name);}
this._vfFieldTypes[name]="MFBoolean";},addField_MFInt32:function(ctx,name,def){this._vf[name]=ctx&&ctx.xmlNode&&ctx.xmlNode.hasAttribute(name)?x3dom.fields.MFInt32.parse(ctx.xmlNode.getAttribute(name)):new x3dom.fields.MFInt32(def);if(ctx&&ctx.xmlNode){this.initSetter(ctx.xmlNode,name);}
this._vfFieldTypes[name]="MFInt32";},addField_MFFloat:function(ctx,name,def){this._vf[name]=ctx&&ctx.xmlNode&&ctx.xmlNode.hasAttribute(name)?x3dom.fields.MFFloat.parse(ctx.xmlNode.getAttribute(name)):new x3dom.fields.MFFloat(def);if(ctx&&ctx.xmlNode){this.initSetter(ctx.xmlNode,name);}
this._vfFieldTypes[name]="MFFloat";},addField_MFDouble:function(ctx,name,def){this._vf[name]=ctx&&ctx.xmlNode&&ctx.xmlNode.hasAttribute(name)?x3dom.fields.MFFloat.parse(ctx.xmlNode.getAttribute(name)):new x3dom.fields.MFFloat(def);if(ctx&&ctx.xmlNode){this.initSetter(ctx.xmlNode,name);}
this._vfFieldTypes[name]="MFDouble";},addField_MFColor:function(ctx,name,def){this._vf[name]=ctx&&ctx.xmlNode&&ctx.xmlNode.hasAttribute(name)?x3dom.fields.MFColor.parse(ctx.xmlNode.getAttribute(name)):new x3dom.fields.MFColor(def);if(ctx&&ctx.xmlNode){this.initSetter(ctx.xmlNode,name);}
this._vfFieldTypes[name]="MFColor";},addField_MFColorRGBA:function(ctx,name,def){this._vf[name]=ctx&&ctx.xmlNode&&ctx.xmlNode.hasAttribute(name)?x3dom.fields.MFColorRGBA.parse(ctx.xmlNode.getAttribute(name)):new x3dom.fields.MFColorRGBA(def);if(ctx&&ctx.xmlNode){this.initSetter(ctx.xmlNode,name);}
this._vfFieldTypes[name]="MFColorRGBA";},addField_MFVec2f:function(ctx,name,def){this._vf[name]=ctx&&ctx.xmlNode&&ctx.xmlNode.hasAttribute(name)?x3dom.fields.MFVec2f.parse(ctx.xmlNode.getAttribute(name)):new x3dom.fields.MFVec2f(def);if(ctx&&ctx.xmlNode){this.initSetter(ctx.xmlNode,name);}
this._vfFieldTypes[name]="MFVec2f";},addField_MFVec3f:function(ctx,name,def){this._vf[name]=ctx&&ctx.xmlNode&&ctx.xmlNode.hasAttribute(name)?x3dom.fields.MFVec3f.parse(ctx.xmlNode.getAttribute(name)):new x3dom.fields.MFVec3f(def);if(ctx&&ctx.xmlNode){this.initSetter(ctx.xmlNode,name);}
this._vfFieldTypes[name]="MFVec3f";},addField_MFVec3d:function(ctx,name,def){this.addField_MFVec3f(ctx,name,def);this._vfFieldTypes[name]="MFVec3d";},addField_MFRotation:function(ctx,name,def){this._vf[name]=ctx&&ctx.xmlNode&&ctx.xmlNode.hasAttribute(name)?x3dom.fields.MFRotation.parse(ctx.xmlNode.getAttribute(name)):new x3dom.fields.MFRotation(def);if(ctx&&ctx.xmlNode){this.initSetter(ctx.xmlNode,name);}
this._vfFieldTypes[name]="MFRotation";},addField_SFNode:function(name,type){this._cf[name]=new x3dom.fields.SFNode(type);this._cfFieldTypes[name]="SFNode";},addField_MFNode:function(name,type){this._cf[name]=new x3dom.fields.MFNode(type);this._cfFieldTypes[name]="MFNode";}}));x3dom.registerNodeType("X3DMetadataObject","Core",defineClass(x3dom.nodeTypes.X3DNode,function(ctx){x3dom.nodeTypes.X3DMetadataObject.superClass.call(this,ctx);this.addField_SFString(ctx,'name',"");this.addField_SFString(ctx,'reference',"");}));x3dom.registerNodeType("MetadataBoolean","Core",defineClass(x3dom.nodeTypes.X3DMetadataObject,function(ctx){x3dom.nodeTypes.MetadataBoolean.superClass.call(this,ctx);this.addField_MFBoolean(ctx,'value',[]);}));x3dom.registerNodeType("MetadataDouble","Core",defineClass(x3dom.nodeTypes.X3DMetadataObject,function(ctx){x3dom.nodeTypes.MetadataDouble.superClass.call(this,ctx);this.addField_MFDouble(ctx,'value',[]);}));x3dom.registerNodeType("MetadataFloat","Core",defineClass(x3dom.nodeTypes.X3DMetadataObject,function(ctx){x3dom.nodeTypes.MetadataFloat.superClass.call(this,ctx);this.addField_MFFloat(ctx,'value',[]);}));x3dom.registerNodeType("MetadataInteger","Core",defineClass(x3dom.nodeTypes.X3DMetadataObject,function(ctx){x3dom.nodeTypes.MetadataInteger.superClass.call(this,ctx);this.addField_MFInt32(ctx,'value',[]);}));x3dom.registerNodeType("MetadataSet","Core",defineClass(x3dom.nodeTypes.X3DMetadataObject,function(ctx){x3dom.nodeTypes.MetadataSet.superClass.call(this,ctx);this.addField_MFNode('value',x3dom.nodeTypes.X3DMetadataObject);}));x3dom.registerNodeType("MetadataString","Core",defineClass(x3dom.nodeTypes.X3DMetadataObject,function(ctx){x3dom.nodeTypes.MetadataString.superClass.call(this,ctx);this.addField_MFString(ctx,'value',[]);}));x3dom.registerNodeType("Field","Core",defineClass(x3dom.nodeTypes.X3DNode,function(ctx){x3dom.nodeTypes.Field.superClass.call(this,ctx);this.addField_SFString(ctx,'name',"");this.addField_SFString(ctx,'type',"");this.addField_SFString(ctx,'value',"");},{fieldChanged:function(fieldName){var that=this;if(fieldName==='value'){Array.forEach(this._parentNodes,function(node){node.fieldChanged(that._vf.name);});}}}));x3dom.registerNodeType("X3DChildNode","Core",defineClass(x3dom.nodeTypes.X3DNode,function(ctx){x3dom.nodeTypes.X3DChildNode.superClass.call(this,ctx);}));x3dom.registerNodeType("X3DBindableNode","Core",defineClass(x3dom.nodeTypes.X3DChildNode,function(ctx){x3dom.nodeTypes.X3DBindableNode.superClass.call(this,ctx);this.addField_SFBool(ctx,'bind',false);this.addField_SFString(ctx,'description',"");this.addField_SFBool(ctx,'isActive',false);this._autoGen=(ctx&&ctx.autoGen?true:false);if(this._autoGen)
this._vf.description="default"+this.constructor.superClass._typeName;this._stack=null;this._bindAnimation=true;},{bind:function(value){if(this._stack){if(value){this._stack.push(this);}
else{this._stack.pop(this);}}
else{x3dom.debug.logError('No BindStack in '+this.typeName()+'Bindable');}},activate:function(prev){this.postMessage('isActive',true);x3dom.debug.logInfo('activate '+this.typeName()+'Bindable '+
this._DEF+'/'+this._vf.description);},deactivate:function(prev){this.postMessage('isActive',false);x3dom.debug.logInfo('deactivate '+this.typeName()+'Bindable '+
this._DEF+'/'+this._vf.description);},fieldChanged:function(fieldName){if(fieldName.indexOf("bind")>=0){this.bind(this._vf.bind);}},nodeChanged:function(){this._stack=this._nameSpace.doc._bindableBag.addBindable(this);}}));x3dom.registerNodeType("X3DInfoNode","Core",defineClass(x3dom.nodeTypes.X3DChildNode,function(ctx){x3dom.nodeTypes.X3DInfoNode.superClass.call(this,ctx);}));x3dom.registerNodeType("WorldInfo","Core",defineClass(x3dom.nodeTypes.X3DInfoNode,function(ctx){x3dom.nodeTypes.WorldInfo.superClass.call(this,ctx);this.addField_MFString(ctx,'info',[]);this.addField_SFString(ctx,'title',"");x3dom.debug.logInfo(this._vf.info);x3dom.debug.logInfo(this._vf.title);}));x3dom.registerNodeType("X3DSensorNode","Core",defineClass(x3dom.nodeTypes.X3DChildNode,function(ctx){x3dom.nodeTypes.X3DSensorNode.superClass.call(this,ctx);this.addField_SFBool(ctx,'enabled',true);}));x3dom.registerNodeType("Param","Core",defineClass(x3dom.nodeTypes.X3DNode,function(ctx){x3dom.nodeTypes.Param.superClass.call(this,ctx);x3dom.debug.logWarning('DEPRECATED: Param element needs to be child of X3D element '
+'[<a href="http://x3dom.org/docs/latest/configuration.html">DOCS</a>]');}));x3dom.registerNodeType("X3DBoundedObject","Grouping",defineClass(x3dom.nodeTypes.X3DChildNode,function(ctx){x3dom.nodeTypes.X3DBoundedObject.superClass.call(this,ctx);this.addField_SFBool(ctx,'render',true);this.addField_SFVec3f(ctx,'bboxCenter',0,0,0);this.addField_SFVec3f(ctx,'bboxSize',-1,-1,-1);this._graph={boundedNode:this,localMatrix:x3dom.fields.SFMatrix4f.identity(),globalMatrix:null,volume:new x3dom.fields.BoxVolume(),lastVolume:new x3dom.fields.BoxVolume(),worldVolume:new x3dom.fields.BoxVolume(),center:new x3dom.fields.SFVec3f(0,0,0),coverage:-1,needCulling:true};},{fieldChanged:function(fieldName){if(this._vf.hasOwnProperty(fieldName)){this.invalidateVolume();}},nodeChanged:function(){this.invalidateVolume();},parentAdded:function(parent){this.invalidateVolume();},getVolume:function()
{var vol=this._graph.volume;if(!this.volumeValid()&&this._vf.render)
{for(var i=0,n=this._childNodes.length;i<n;i++)
{var child=this._childNodes[i];if(!child||child._vf.render!==true)
continue;var childVol=child.getVolume();if(childVol&&childVol.isValid())
vol.extendBounds(childVol.min,childVol.max);}}
if(!vol.equals(this._graph.lastVolume))
{this._graph.lastVolume=x3dom.fields.BoxVolume.copy(vol);var event={target:this._xmlNode,type:"volumechanged",volume:x3dom.fields.BoxVolume.copy(vol)};this.callEvtHandler("onvolumechanged",event);}
return vol;},invalidateVolume:function()
{var graph=this._graph;graph.volume.invalidate();graph.worldVolume.invalidate();graph.globalMatrix=null;for(var i=0,n=this._parentNodes.length;i<n;i++){var node=this._parentNodes[i];if(node)
node.invalidateVolume();}},invalidateCache:function()
{var graph=this._graph;graph.worldVolume.invalidate();graph.globalMatrix=null;},cacheInvalid:function()
{return(this._graph.globalMatrix==null||!this._graph.worldVolume.isValid());},volumeValid:function()
{return this._graph.volume.isValid();},graphState:function()
{return this._graph;},forceUpdateCoverage:function()
{return false;}}));x3dom.registerNodeType("X3DGroupingNode","Grouping",defineClass(x3dom.nodeTypes.X3DBoundedObject,function(ctx){x3dom.nodeTypes.X3DGroupingNode.superClass.call(this,ctx);this.addField_MFNode('children',x3dom.nodeTypes.X3DChildNode);},{collectDrawableObjects:function(transform,drawableCollection,singlePath,invalidateCache,planeMask,clipPlanes)
{if(singlePath&&(this._parentNodes.length>1))
singlePath=false;if(singlePath&&(invalidateCache=invalidateCache||this.cacheInvalid()))
this.invalidateCache();planeMask=drawableCollection.cull(transform,this.graphState(),singlePath,planeMask);if(planeMask<0){return;}
var cnode,childTransform;if(singlePath){if(!this._graph.globalMatrix){this._graph.globalMatrix=this.transformMatrix(transform);}
childTransform=this._graph.globalMatrix;}
else{childTransform=this.transformMatrix(transform);}
var n=this._childNodes.length;if(x3dom.nodeTypes.ClipPlane.count>0){var localClipPlanes=[];for(var j=0;j<n;j++){if((cnode=this._childNodes[j])){if(x3dom.isa(cnode,x3dom.nodeTypes.ClipPlane)&&cnode._vf.on&&cnode._vf.enabled){localClipPlanes.push({plane:cnode,trafo:childTransform});}}}
clipPlanes=localClipPlanes.concat(clipPlanes);}
for(var i=0;i<n;i++){if((cnode=this._childNodes[i])){cnode.collectDrawableObjects(childTransform,drawableCollection,singlePath,invalidateCache,planeMask,clipPlanes);}}}}));x3dom.registerNodeType("Switch","Grouping",defineClass(x3dom.nodeTypes.X3DGroupingNode,function(ctx){x3dom.nodeTypes.Switch.superClass.call(this,ctx);this.addField_SFInt32(ctx,'whichChoice',-1);},{fieldChanged:function(fieldName){if(fieldName=="whichChoice"){this.invalidateVolume();}},getVolume:function()
{var vol=this._graph.volume;if(!this.volumeValid()&&this._vf.render)
{if(this._vf.whichChoice>=0&&this._vf.whichChoice<this._childNodes.length)
{var child=this._childNodes[this._vf.whichChoice];var childVol=(child&&child._vf.render===true)?child.getVolume():null;if(childVol&&childVol.isValid())
vol.extendBounds(childVol.min,childVol.max);}}
return vol;},collectDrawableObjects:function(transform,drawableCollection,singlePath,invalidateCache,planeMask,clipPlanes)
{if(singlePath&&(this._parentNodes.length>1))
singlePath=false;if(singlePath&&(invalidateCache=invalidateCache||this.cacheInvalid()))
this.invalidateCache();if(this._vf.whichChoice<0||this._vf.whichChoice>=this._childNodes.length||(planeMask=drawableCollection.cull(transform,this.graphState(),singlePath,planeMask))<0){return;}
var cnode,childTransform;if(singlePath){if(!this._graph.globalMatrix){this._graph.globalMatrix=this.transformMatrix(transform);}
childTransform=this._graph.globalMatrix;}
else{childTransform=this.transformMatrix(transform);}
if((cnode=this._childNodes[this._vf.whichChoice])){cnode.collectDrawableObjects(childTransform,drawableCollection,singlePath,invalidateCache,planeMask,clipPlanes);}},doIntersect:function(line)
{if(this._vf.whichChoice<0||this._vf.whichChoice>=this._childNodes.length){return false;}
var child=this._childNodes[this._vf.whichChoice];if(child){return child.doIntersect(line);}
return false;}}));x3dom.registerNodeType("X3DTransformNode","Grouping",defineClass(x3dom.nodeTypes.X3DGroupingNode,function(ctx){x3dom.nodeTypes.X3DTransformNode.superClass.call(this,ctx);if(ctx)
ctx.doc._nodeBag.trans.push(this);else
x3dom.debug.logWarning("X3DTransformNode: No runtime context found!");this._trafo=null;this._needCssStyleUpdates=true;},{tick:function(t)
{var dom=this._xmlNode;if(dom&&(dom['ontransform']||dom.hasAttribute('ontransform')||this._listeners['transform'])){var transMatrix=this.getCurrentTransform();var event={target:dom,type:'transform',worldX:transMatrix._03,worldY:transMatrix._13,worldZ:transMatrix._23,cancelBubble:false,stopPropagation:function(){this.cancelBubble=true;}};this.callEvtHandler("ontransform",event);}
if(this._needCssStyleUpdates&&dom){var trans=x3dom.getStyle(dom,"-webkit-transform")||x3dom.getStyle(dom,"-moz-transform")||x3dom.getStyle(dom,"-ms-transform")||x3dom.getStyle(dom,"transform");if(trans&&(trans!='none')){this._trafo.setValueByStr(trans);this.invalidateVolume();return true;}
this._needCssStyleUpdates=false;}
return false;},transformMatrix:function(transform){return transform.mult(this._trafo);},getVolume:function()
{var vol=this._graph.volume;if(!this.volumeValid()&&this._vf.render)
{this._graph.localMatrix=this._trafo;for(var i=0,n=this._childNodes.length;i<n;i++)
{var child=this._childNodes[i];if(!child||child._vf.render!==true)
continue;var childVol=child.getVolume();if(childVol&&childVol.isValid())
vol.extendBounds(childVol.min,childVol.max);}
if(vol.isValid())
vol.transform(this._trafo);}
return vol;},doIntersect:function(line)
{var isect=false;var mat=this._trafo.inverse();var tmpPos=new x3dom.fields.SFVec3f(line.pos.x,line.pos.y,line.pos.z);var tmpDir=new x3dom.fields.SFVec3f(line.dir.x,line.dir.y,line.dir.z);line.pos=mat.multMatrixPnt(line.pos);line.dir=mat.multMatrixVec(line.dir);if(line.hitObject){line.dist*=line.dir.length();}
for(var i=0;i<this._childNodes.length;i++)
{if(this._childNodes[i]){isect=this._childNodes[i].doIntersect(line)||isect;}}
line.pos.setValues(tmpPos);line.dir.setValues(tmpDir);if(isect){line.hitPoint=this._trafo.multMatrixPnt(line.hitPoint);line.dist*=line.dir.length();}
return isect;},parentRemoved:function(parent)
{var i,n;if(this._parentNodes.length==0){var doc=this.findX3DDoc();for(i=0,n=doc._nodeBag.trans.length;i<n;i++){if(doc._nodeBag.trans[i]===this){doc._nodeBag.trans.splice(i,1);}}}
for(i=0,n=this._childNodes.length;i<n;i++){if(this._childNodes[i]){this._childNodes[i].parentRemoved(this);}}}}));x3dom.registerNodeType("Transform","Grouping",defineClass(x3dom.nodeTypes.X3DTransformNode,function(ctx){x3dom.nodeTypes.Transform.superClass.call(this,ctx);this.addField_SFVec3f(ctx,'center',0,0,0);this.addField_SFVec3f(ctx,'translation',0,0,0);this.addField_SFRotation(ctx,'rotation',0,0,1,0);this.addField_SFVec3f(ctx,'scale',1,1,1);this.addField_SFRotation(ctx,'scaleOrientation',0,0,1,0);this._trafo=x3dom.fields.SFMatrix4f.translation(this._vf.translation.add(this._vf.center)).mult(this._vf.rotation.toMatrix()).mult(this._vf.scaleOrientation.toMatrix()).mult(x3dom.fields.SFMatrix4f.scale(this._vf.scale)).mult(this._vf.scaleOrientation.toMatrix().inverse()).mult(x3dom.fields.SFMatrix4f.translation(this._vf.center.negate()));},{fieldChanged:function(fieldName)
{if(fieldName=="center"||fieldName=="translation"||fieldName=="rotation"||fieldName=="scale"||fieldName=="scaleOrientation")
{this._trafo=x3dom.fields.SFMatrix4f.translation(this._vf.translation.add(this._vf.center)).mult(this._vf.rotation.toMatrix()).mult(this._vf.scaleOrientation.toMatrix()).mult(x3dom.fields.SFMatrix4f.scale(this._vf.scale)).mult(this._vf.scaleOrientation.toMatrix().inverse()).mult(x3dom.fields.SFMatrix4f.translation(this._vf.center.negate()));this.invalidateVolume();}
else if(fieldName=="render"){this.invalidateVolume();}}}));x3dom.registerNodeType("MatrixTransform","Grouping",defineClass(x3dom.nodeTypes.X3DTransformNode,function(ctx){x3dom.nodeTypes.MatrixTransform.superClass.call(this,ctx);this.addField_SFMatrix4f(ctx,'matrix',1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1);this._trafo=this._vf.matrix.transpose();},{fieldChanged:function(fieldName){if(fieldName=="matrix"){this._trafo=this._vf.matrix.transpose();this.invalidateVolume();}
else if(fieldName=="render"){this.invalidateVolume();}}}));x3dom.registerNodeType("Group","Grouping",defineClass(x3dom.nodeTypes.X3DGroupingNode,function(ctx){x3dom.nodeTypes.Group.superClass.call(this,ctx);}));x3dom.registerNodeType("Block","Grouping",defineClass(x3dom.nodeTypes.X3DGroupingNode,function(ctx){x3dom.nodeTypes.Block.superClass.call(this,ctx);this.addField_MFString(ctx,'nameSpaceName',[]);}));x3dom.registerNodeType("StaticGroup","Grouping",defineClass(x3dom.nodeTypes.X3DGroupingNode,function(ctx){x3dom.nodeTypes.StaticGroup.superClass.call(this,ctx);x3dom.debug.logWarning("StaticGroup erroneously also bakes parent transforms, if happens use Group node!");this.addField_SFBool(ctx,'debug',false);this.addField_SFBool(ctx,'showDebugBoxVolumes',false);this.addField_SFString(ctx,'bvhType','jsBIH');this.addField_SFInt32(ctx,'maxObjectsPerNode',1);this.addField_SFInt32(ctx,'maxDepth',-1);this.addField_SFFloat(ctx,'minRelativeBBoxSize',0.01);this.needBvhRebuild=true;this.drawableCollection=null;this.bvh=null;},{getMaxDepth:function()
{if(this._vf.maxDepth==-1)
{return(this._vf.bvhType==('jsBIH'||'BIH'))?50:4;}
return this._vf.maxDepth;},collectDrawableObjects:function(transform,drawableCollection,singlePath,invalidateCache,planeMask,clipPlanes)
{if(singlePath&&(this._parentNodes.length>1))
singlePath=false;if(singlePath&&(invalidateCache=invalidateCache||this.cacheInvalid()))
this.invalidateCache();planeMask=drawableCollection.cull(transform,this.graphState(),singlePath,planeMask);if(planeMask<0){return;}
var cnode,childTransform;if(singlePath){if(!this._graph.globalMatrix){this._graph.globalMatrix=this.transformMatrix(transform);}
childTransform=this._graph.globalMatrix;}
else{childTransform=this.transformMatrix(transform);}
if(this.needBvhRebuild)
{var drawableCollectionConfig={viewArea:drawableCollection.viewarea,sortTrans:drawableCollection.sortTrans,viewMatrix:drawableCollection.viewMatrix,projMatrix:drawableCollection.projMatrix,sceneMatrix:drawableCollection.sceneMatrix,frustumCulling:false,smallFeatureThreshold:0,context:drawableCollection.context};this.drawableCollection=new x3dom.DrawableCollection(drawableCollectionConfig);var i,n=this._childNodes.length;for(i=0;i<n;i++){if((cnode=this._childNodes[i])){cnode.collectDrawableObjects(childTransform,this.drawableCollection,singlePath,invalidateCache,planeMask,clipPlanes);}}
this.drawableCollection.concat();var scene=this._nameSpace.doc._scene;var bvhSettings=new x3dom.bvh.Settings(this._vf.debug,this._vf.showDebugBoxVolumes,this._vf.bvhType,this._vf.maxObjectsPerNode,this.getMaxDepth(),this._vf.minRelativeBBoxSize);this.bvh=(this._vf.bvhType=='jsBIH')?new x3dom.bvh.BIH(scene,bvhSettings):new x3dom.bvh.Culler(this.drawableCollection,scene,bvhSettings);if(this._vf.debug||this._vf.showDebugBoxVolumes)
this.bvh=new x3dom.bvh.DebugDecorator(this.bvh,scene,bvhSettings);n=this.drawableCollection.length;for(i=0;i<n;i++)
{this.bvh.addDrawable(this.drawableCollection.get(i))}
this.bvh.compile();if(this._vf.debug)
this.bvh.showCompileStats();this.needBvhRebuild=false;}
x3dom.Utils.startMeasure('bvhTraverse');this.bvh.collectDrawables(drawableCollection);var dt=x3dom.Utils.stopMeasure('bvhTraverse');this._nameSpace.doc.ctx.x3dElem.runtime.addMeasurement('BVH',dt);this.bvh.showTraverseStats(this._nameSpace.doc.ctx.x3dElem.runtime);}}));x3dom.registerNodeType("RemoteSelectionGroup","Grouping",defineClass(x3dom.nodeTypes.X3DGroupingNode,function(ctx){x3dom.nodeTypes.RemoteSelectionGroup.superClass.call(this,ctx);this.addField_MFString(ctx,'url',["ws://localhost:35668/cstreams/0"]);this.addField_MFString(ctx,'label',[]);this.addField_SFInt32(ctx,'maxRenderedIds',-1);this.addField_SFBool(ctx,'reconnect',true);this.addField_SFFloat(ctx,'scaleRenderedIdsOnMove',1.0);this.addField_SFBool(ctx,'enableCulling',true);this.addField_MFString(ctx,'invisibleNodes',[]);this._idList=[];this._websocket=null;this._nameObjMap={};this._createTime=[];this._visibleList=[];if(ctx)
this.initializeSocket();else
x3dom.debug.logWarning("RemoteSelectionGroup: No runtime context found!");},{initializeSocket:function()
{var that=this;if("WebSocket"in window)
{var wsUrl="ws://localhost:35668/cstreams/0";if(this._vf.url.length&&this._vf.url[0].length)
wsUrl=this._vf.url[0];this._websocket=new WebSocket(wsUrl);this._websocket._lastMsg=null;this._websocket._lastData="";this._websocket.onopen=function(evt)
{x3dom.debug.logInfo("WS Connected");var view=that._nameSpace.doc._viewarea.getViewMatrix();this._lastMsg=view.toGL().toString();view=that._nameSpace.doc._viewarea.getProjectionMatrix();this._lastMsg+=(","+view.toGL().toString());this.send(this._lastMsg);x3dom.debug.logInfo("WS Sent: "+this._lastMsg);this._lastMsg="";this._lastData="";};this._websocket.onclose=function(evt)
{x3dom.debug.logInfo("WS Disconnected");if(that._vf.reconnect)
{window.setTimeout(function(){that.initializeSocket();},2000);}};this._websocket.onmessage=function(evt)
{if(that._vf.maxRenderedIds<0)
{that._idList=x3dom.fields.MFString.parse(evt.data);}
else if(that._vf.maxRenderedIds>0)
{that._idList=[];var arr=x3dom.fields.MFString.parse(evt.data);var n=Math.min(arr.length,Math.abs(that._vf.maxRenderedIds));for(var i=0;i<n;++i){that._idList[i]=arr[i];}}
if(that._vf.maxRenderedIds!=0&&this._lastData!=evt.data)
{this._lastData=evt.data;that._nameSpace.doc.needRender=true;that.invalidateVolume();}};this._websocket.onerror=function(evt)
{x3dom.debug.logError(evt.data);};this._websocket.updateCamera=function()
{var view=that._nameSpace.doc._viewarea.getViewMatrix();var message=view.toGL().toString();view=that._nameSpace.doc._viewarea.getProjectionMatrix();message+=(","+view.toGL().toString());if(this._lastMsg!=null&&this._lastMsg!=message)
{this._lastMsg=message;this.send(message);}};}
else
{x3dom.debug.logError("Browser has no WebSocket support!");}},nodeChanged:function()
{var n=this._vf.label.length;this._nameObjMap={};this._createTime=new Array(n);this._visibleList=new Array(n);for(var i=0;i<n;++i)
{var shape=this._childNodes[i];if(shape&&x3dom.isa(shape,x3dom.nodeTypes.X3DShapeNode))
{this._nameObjMap[this._vf.label[i]]={shape:shape,pos:i};this._visibleList[i]=true;}
else{this._visibleList[i]=false;x3dom.debug.logError("Invalid children: "+this._vf.label[i]);}
this._createTime[i]=0;}
this.invalidateVolume();x3dom.debug.logInfo("RemoteSelectionGroup has "+n+" entries.");},fieldChanged:function(fieldName)
{if(fieldName=="url")
{if(this._websocket){this._websocket.close();this._websocket=null;}
this.initializeSocket();}
else if(fieldName=="invisibleNodes")
{for(var i=0,n=this._vf.label.length;i<n;++i)
{var shape=this._childNodes[i];if(shape&&x3dom.isa(shape,x3dom.nodeTypes.X3DShapeNode))
{this._visibleList[i]=true;for(var j=0,numInvis=this._vf.invisibleNodes.length;j<numInvis;++j)
{var nodeName=this._vf.invisibleNodes[j];var starInd=nodeName.lastIndexOf('*');var matchNameBegin=false;if(starInd>0){nodeName=nodeName.substring(0,starInd);matchNameBegin=true;}
if(nodeName.length<=1)
continue;if((matchNameBegin&&this._vf.label[i].indexOf(nodeName)==0)||this._vf.label[i]==nodeName){this._visibleList[i]=false;break;}}}
else{this._visibleList[i]=false;}}
this.invalidateVolume();}
else if(fieldName=="render"){this.invalidateVolume();}},getNumRenderedObjects:function(len,isMoving)
{var n=len;if(this._vf.maxRenderedIds>0)
{var num=Math.max(this._vf.maxRenderedIds,16);var scale=1;if(isMoving)
scale=Math.min(this._vf.scaleRenderedIdsOnMove,1);num=Math.max(Math.round(scale*num),0);n=Math.min(n,num);}
return n;},collectDrawableObjects:function(transform,drawableCollection,singlePath,invalidateCache,planeMask,clipPlanes)
{if(singlePath&&(this._parentNodes.length>1))
singlePath=false;if(singlePath&&(invalidateCache=invalidateCache||this.cacheInvalid()))
this.invalidateCache();planeMask=drawableCollection.cull(transform,this.graphState(),singlePath,planeMask);if(planeMask<=0){return;}
var viewarea=this._nameSpace.doc._viewarea;var isMoving=viewarea.isMovingOrAnimating();var ts=new Date().getTime();var maxLiveTime=10000;var i,n,numChild=this._childNodes.length;if(!this._vf.enableCulling)
{n=this.getNumRenderedObjects(numChild,isMoving);var cnt=0;for(i=0;i<numChild;i++)
{var shape=this._childNodes[i];if(shape)
{var needCleanup=true;if(this._visibleList[i]&&cnt<n&&shape.collectDrawableObjects(transform,drawableCollection,singlePath,invalidateCache,planeMask,clipPlanes))
{this._createTime[i]=ts;cnt++;needCleanup=false;}
if(needCleanup&&!isMoving&&this._createTime[i]>0&&ts-this._createTime[i]>maxLiveTime&&shape._cleanupGLObjects)
{shape._cleanupGLObjects(true);this._createTime[i]=0;}}}
return;}
if(this._websocket)
this._websocket.updateCamera();if(this._vf.label.length)
{n=this.getNumRenderedObjects(this._idList.length,isMoving);for(i=0;i<n;i++)
{var obj=this._nameObjMap[this._idList[i]];if(obj&&obj.shape){obj.shape.collectDrawableObjects(transform,drawableCollection,singlePath,invalidateCache,planeMask,clipPlanes);this._createTime[obj.pos]=ts;}
else
x3dom.debug.logError("Invalid label: "+this._idList[i]);}
for(i=0;i<this._childNodes.length;i++)
{if(this._childNodes[i]&&!isMoving&&this._createTime[i]>0&&ts-this._createTime[i]>maxLiveTime&&this._childNodes[i]._cleanupGLObjects)
{this._childNodes[i]._cleanupGLObjects(true);this._createTime[i]=0;}}}}}));x3dom.registerNodeType("Scene","Grouping",defineClass(x3dom.nodeTypes.X3DGroupingNode,function(ctx){x3dom.nodeTypes.Scene.superClass.call(this,ctx);this.addField_SFString(ctx,'pickMode',"idBuf");this.addField_SFBool(ctx,'doPickPass',true);this.addField_SFString(ctx,'shadowObjectIdMapping',"");this._lastMin=new x3dom.fields.SFVec3f(0,0,0);this._lastMax=new x3dom.fields.SFVec3f(1,1,1);this._shadowIdMap=null;this.loadMapping();this._multiPartMap=null;},{fieldChanged:function(fieldName)
{if(fieldName=="shadowObjectIdMapping")
{this.loadMapping();}},updateVolume:function()
{var vol=this.getVolume();if(vol.isValid())
{this._lastMin=x3dom.fields.SFVec3f.copy(vol.min);this._lastMax=x3dom.fields.SFVec3f.copy(vol.max);}},loadMapping:function()
{this._shadowIdMap=null;if(this._vf.shadowObjectIdMapping.length==0){return;}
var that=this;var xhr=new XMLHttpRequest();xhr.open("GET",this._nameSpace.getURL(this._vf.shadowObjectIdMapping),true);x3dom.RequestManager.addRequest(xhr);xhr.onload=function()
{that._shadowIdMap=eval("("+xhr.response+")");if(!that._shadowIdMap||!that._shadowIdMap.mapping){x3dom.debug.logWarning("Invalid ID map: "+that._vf.shadowObjectIdMapping);}
else{x3dom.debug.assert(that._shadowIdMap.maxID<=that._shadowIdMap.mapping.length,"Too few ID map entries in "+that._vf.shadowObjectIdMapping+", "+"length of mapping array is only "+that._shadowIdMap.mapping.length+" instead of "+that._shadowIdMap.ids.length+"!");}};}}));x3dom.BindableStack=function(doc,type,defaultType,getter){this._doc=doc;this._type=type;this._defaultType=defaultType;this._defaultRoot=null;this._getter=getter;this._bindBag=[];this._bindStack=[];};x3dom.BindableStack.prototype.top=function(){return((this._bindStack.length>0)?this._bindStack[this._bindStack.length-1]:null);};x3dom.BindableStack.prototype.push=function(bindable){var top=this.top();if(top===bindable){return;}
if(top){top.deactivate();}
this._bindStack.push(bindable);bindable.activate(top);};x3dom.BindableStack.prototype.replaceTop=function(bindable){var top=this.top();if(top===bindable){return;}
if(top){top.deactivate();this._bindStack[this._bindStack.length-1]=bindable;bindable.activate(top);}};x3dom.BindableStack.prototype.pop=function(bindable){var top;if(bindable){top=this.top();if(bindable!==top){return null;}}
top=this._bindStack.pop();if(top){top.deactivate();}
return top;};x3dom.BindableStack.prototype.switchTo=function(target){var last=this.getActive();var n=this._bindBag.length;var toBind=0;var i=0,lastIndex=-1;if(n<=1){return;}
switch(target)
{case'first':toBind=this._bindBag[0];break;case'last':toBind=this._bindBag[n-1];break;default:for(i=0;i<n;i++){if(this._bindBag[i]==last){lastIndex=i;break;}}
if(lastIndex>=0){i=lastIndex;while(!toBind){if(target=='next'){i=(i<(n-1))?(i+1):0;}else{i=(i>0)?(i-1):(n-1);}
if(i==lastIndex){break;}
if(this._bindBag[i]._vf.description.length>=0){toBind=this._bindBag[i];}}}
break;}
if(toBind){this.replaceTop(toBind);}else{x3dom.debug.logWarning('Cannot switch bindable; no other bindable with description found.');}};x3dom.BindableStack.prototype.getActive=function(){if(this._bindStack.length===0){if(this._bindBag.length===0){if(this._defaultRoot){x3dom.debug.logInfo('create new '+this._defaultType._typeName+' for '+this._type._typeName+'-stack');var obj=new this._defaultType({doc:this._doc,nameSpace:this._defaultRoot._nameSpace,autoGen:true});this._defaultRoot.addChild(obj);obj.nodeChanged();}
else{x3dom.debug.logError('stack without defaultRoot');}}
else{x3dom.debug.logInfo('activate first '+this._type._typeName+' for '+this._type._typeName+'-stack');}
this._bindStack.push(this._bindBag[0]);this._bindBag[0].activate();}
return this._bindStack[this._bindStack.length-1];};x3dom.BindableBag=function(doc){this._stacks=[];this.addType("X3DViewpointNode","Viewpoint","getViewpoint",doc);this.addType("X3DNavigationInfoNode","NavigationInfo","getNavigationInfo",doc);this.addType("X3DBackgroundNode","Background","getBackground",doc);this.addType("X3DFogNode","Fog","getFog",doc);this.addType("X3DEnvironmentNode","Environment","getEnvironment",doc);};x3dom.BindableBag.prototype.addType=function(typeName,defaultTypeName,getter,doc){var type=x3dom.nodeTypes[typeName];var defaultType=x3dom.nodeTypes[defaultTypeName];if(type&&defaultType){var stack=new x3dom.BindableStack(doc,type,defaultType,getter);this._stacks.push(stack);}
else{x3dom.debug.logWarning('Invalid Bindable type/defaultType: '+
typeName+'/'+defaultType);}};x3dom.BindableBag.prototype.setRefNode=function(node){Array.forEach(this._stacks,function(stack){stack._defaultRoot=node;node[stack._getter]=function(){return stack.getActive();};});};x3dom.BindableBag.prototype.addBindable=function(node){for(var i=0,n=this._stacks.length;i<n;i++){var stack=this._stacks[i];if(x3dom.isa(node,stack._type)){x3dom.debug.logInfo('register '+node.typeName()+'Bindable '+
node._DEF+'/'+node._vf.description);stack._bindBag.push(node);var top=stack.top();if(top&&top._autoGen){stack.replaceTop(node);for(var j=0,m=stack._bindBag.length;j<m;j++){if(stack._bindBag[j]===top){stack._bindBag.splice(j,1);break;}}
stack._defaultRoot.removeChild(top);}
return stack;}}
x3dom.debug.logError(node.typeName()+' is not a valid bindable');return null;};x3dom.registerNodeType("X3DGeometryNode","Rendering",defineClass(x3dom.nodeTypes.X3DNode,function(ctx){x3dom.nodeTypes.X3DGeometryNode.superClass.call(this,ctx);this.addField_SFBool(ctx,'solid',true);this.addField_SFBool(ctx,'ccw',true);this.addField_SFBool(ctx,'useGeoCache',true);this.addField_SFBool(ctx,'lit',true);this._mesh=new x3dom.Mesh(this);},{getVolume:function(){return this._mesh.getVolume();},invalidateVolume:function(){this._mesh.invalidate();},getCenter:function(){return this._mesh.getCenter();},getDiameter:function(){return this._mesh.getDiameter();},doIntersect:function(line){return this._mesh.doIntersect(line);},forceUpdateCoverage:function(){return false;},hasIndexOffset:function(){return false;},getColorTexture:function(){return null;},getColorTextureURL:function(){return null;},parentAdded:function(parent){if(x3dom.isa(parent,x3dom.nodeTypes.X3DShapeNode)){if(parent._cleanupGLObjects){parent._cleanupGLObjects(true);}
parent.setAllDirty();parent.invalidateVolume();}},needLighting:function(){var hasTris=this._mesh._primType.indexOf("TRIANGLE")>=0;return(this._vf.lit&&hasTris);}}));x3dom.registerNodeType("Mesh","Rendering",defineClass(x3dom.nodeTypes.X3DGeometryNode,function(ctx){x3dom.nodeTypes.Mesh.superClass.call(this,ctx);this.addField_SFString(ctx,'primType',"triangle");this.addField_MFInt32(ctx,'index',[]);this.addField_MFNode('vertexAttributes',x3dom.nodeTypes.X3DVertexAttributeNode);},{nodeChanged:function()
{var time0=new Date().getTime();var i,n=this._cf.vertexAttributes.nodes.length;for(i=0;i<n;i++)
{var name=this._cf.vertexAttributes.nodes[i]._vf.name;switch(name.toLowerCase())
{case"position":this._mesh._positions[0]=this._cf.vertexAttributes.nodes[i]._vf.value.toGL();break;case"normal":this._mesh._normals[0]=this._cf.vertexAttributes.nodes[i]._vf.value.toGL();break;case"texcoord":this._mesh._texCoords[0]=this._cf.vertexAttributes.nodes[i]._vf.value.toGL();break;case"color":this._mesh._colors[0]=this._cf.vertexAttributes.nodes[i]._vf.value.toGL();break;default:this._mesh._dynamicFields[name]={};this._mesh._dynamicFields[name].numComponents=this._cf.vertexAttributes.nodes[i]._vf.numComponents;this._mesh._dynamicFields[name].value=this._cf.vertexAttributes.nodes[i]._vf.value.toGL();break;}}
this._mesh._indices[0]=this._vf.index.toGL();this.invalidateVolume();this._mesh._numFaces=this._mesh._indices[0].length/3;this._mesh._numCoords=this._mesh._positions[0].length/3;var time1=new Date().getTime()-time0;x3dom.debug.logWarning("Mesh load time: "+time1+" ms");}}));x3dom.registerNodeType("PointSet","Rendering",defineClass(x3dom.nodeTypes.X3DGeometryNode,function(ctx){x3dom.nodeTypes.PointSet.superClass.call(this,ctx);this.addField_SFNode('coord',x3dom.nodeTypes.X3DCoordinateNode);this.addField_SFNode('color',x3dom.nodeTypes.X3DColorNode);this._mesh._primType='POINTS';},{nodeChanged:function()
{var time0=new Date().getTime();var coordNode=this._cf.coord.node;x3dom.debug.assert(coordNode,"PointSet without coord node!");var positions=coordNode.getPoints();var numColComponents=3;var colorNode=this._cf.color.node;var colors=new x3dom.fields.MFColor();if(colorNode){colors=colorNode._vf.color;x3dom.debug.assert(positions.length==colors.length,"Size of color and coord array differs!");if(x3dom.isa(colorNode,x3dom.nodeTypes.ColorRGBA)){numColComponents=4;}}
this._mesh._numColComponents=numColComponents;this._mesh._lit=false;this._mesh._indices[0]=[];this._mesh._positions[0]=positions.toGL();this._mesh._colors[0]=colors.toGL();this._mesh._normals[0]=[];this._mesh._texCoords[0]=[];this.invalidateVolume();this._mesh._numCoords=this._mesh._positions[0].length/3;var time1=new Date().getTime()-time0;},fieldChanged:function(fieldName)
{var pnts=null;if(fieldName=="coord")
{pnts=this._cf.coord.node.getPoints();this._mesh._positions[0]=pnts.toGL();this.invalidateVolume();Array.forEach(this._parentNodes,function(node){node._dirty.positions=true;node.invalidateVolume();});}
else if(fieldName=="color")
{pnts=this._cf.color.node._vf.color;this._mesh._colors[0]=pnts.toGL();Array.forEach(this._parentNodes,function(node){node._dirty.colors=true;});}}}));x3dom.registerNodeType("X3DComposedGeometryNode","Rendering",defineClass(x3dom.nodeTypes.X3DGeometryNode,function(ctx){x3dom.nodeTypes.X3DComposedGeometryNode.superClass.call(this,ctx);this.addField_SFBool(ctx,'colorPerVertex',true);this.addField_SFBool(ctx,'normalPerVertex',true);this.addField_SFString(ctx,'normalUpdateMode','fast');this.addField_MFNode('attrib',x3dom.nodeTypes.X3DVertexAttributeNode);this.addField_SFNode('coord',x3dom.nodeTypes.X3DCoordinateNode);this.addField_SFNode('normal',x3dom.nodeTypes.Normal);this.addField_SFNode('color',x3dom.nodeTypes.X3DColorNode);this.addField_SFNode('texCoord',x3dom.nodeTypes.X3DTextureCoordinateNode);},{handleAttribs:function()
{var i,n=this._cf.attrib.nodes.length;for(i=0;i<n;i++)
{var name=this._cf.attrib.nodes[i]._vf.name;switch(name.toLowerCase())
{case"position":this._mesh._positions[0]=this._cf.attrib.nodes[i]._vf.value.toGL();break;case"normal":this._mesh._normals[0]=this._cf.attrib.nodes[i]._vf.value.toGL();break;case"texcoord":this._mesh._texCoords[0]=this._cf.attrib.nodes[i]._vf.value.toGL();break;case"color":this._mesh._colors[0]=this._cf.attrib.nodes[i]._vf.value.toGL();break;default:this._mesh._dynamicFields[name]={};this._mesh._dynamicFields[name].numComponents=this._cf.attrib.nodes[i]._vf.numComponents;this._mesh._dynamicFields[name].value=this._cf.attrib.nodes[i]._vf.value.toGL();break;}}}}));x3dom.registerNodeType("LineSet","Rendering",defineClass(x3dom.nodeTypes.X3DGeometryNode,function(ctx){x3dom.nodeTypes.LineSet.superClass.call(this,ctx);this.addField_MFInt32(ctx,'vertexCount',[]);this.addField_MFNode('attrib',x3dom.nodeTypes.X3DVertexAttributeNode);this.addField_SFNode('coord',x3dom.nodeTypes.X3DCoordinateNode);this.addField_SFNode('color',x3dom.nodeTypes.X3DColorNode);this._mesh._primType="LINES";x3dom.Utils.needLineWidth=true;},{nodeChanged:function(){var coordNode=this._cf.coord.node;x3dom.debug.assert(coordNode);var positions=coordNode.getPoints();this._mesh._positions[0]=positions.toGL();var colorNode=this._cf.color.node;if(colorNode){var colors=colorNode._vf.color;this._mesh._colors[0]=colors.toGL();this._mesh._numColComponents=3;if(x3dom.isa(colorNode,x3dom.nodeTypes.ColorRGBA)){this._mesh._numColComponents=4;}}
var cnt=0;this._mesh._indices[0]=[];for(var i=0,n=this._vf.vertexCount.length;i<n;i++){var vc=this._vf.vertexCount[i];if(vc<2){x3dom.debug.logError("LineSet.vertexCount must not be smaller than 2!");break;}
for(var j=vc-2;j>=0;j--){this._mesh._indices[0].push(cnt++,cnt);if(j==0)cnt++;}}},fieldChanged:function(fieldName){if(fieldName=="coord"){var pnts=this._cf.coord.node.getPoints();this._mesh._positions[0]=pnts.toGL();this.invalidateVolume();Array.forEach(this._parentNodes,function(node){node._dirty.positions=true;node.invalidateVolume();});}
else if(fieldName=="color"){var cols=this._cf.color.node._vf.color;this._mesh._colors[0]=cols.toGL();Array.forEach(this._parentNodes,function(node){node._dirty.colors=true;});}}}));x3dom.registerNodeType("IndexedLineSet","Rendering",defineClass(x3dom.nodeTypes.X3DGeometryNode,function(ctx){x3dom.nodeTypes.IndexedLineSet.superClass.call(this,ctx);this.addField_SFBool(ctx,'colorPerVertex',true);this.addField_MFNode('attrib',x3dom.nodeTypes.X3DVertexAttributeNode);this.addField_SFNode('coord',x3dom.nodeTypes.X3DCoordinateNode);this.addField_SFNode('color',x3dom.nodeTypes.X3DColorNode);this.addField_MFInt32(ctx,'coordIndex',[]);this.addField_MFInt32(ctx,'colorIndex',[]);this._mesh._primType='LINES';x3dom.Utils.needLineWidth=true;},{_buildGeometry:function()
{var time0=new Date().getTime();var indexes=this._vf.coordIndex;var colorInd=this._vf.colorIndex;var hasColor=false,hasColorInd=false;var colPerVert=this._vf.colorPerVertex;if(colorInd.length==indexes.length)
{hasColorInd=true;}
var positions,colors;var coordNode=this._cf.coord.node;x3dom.debug.assert(coordNode);positions=coordNode.getPoints();var numColComponents=3;var colorNode=this._cf.color.node;if(colorNode)
{hasColor=true;colors=colorNode._vf.color;if(x3dom.isa(colorNode,x3dom.nodeTypes.ColorRGBA)){numColComponents=4;}}
else{hasColor=false;}
this._mesh._indices[0]=[];this._mesh._positions[0]=[];this._mesh._colors[0]=[];var i,t,cnt,lineCnt;var p0,p1,c0,c1;if((hasColor&&hasColorInd)||positions.length>x3dom.Utils.maxIndexableCoords)
{t=0;cnt=0;lineCnt=0;for(i=0;i<indexes.length;++i)
{if(indexes[i]>positions.length-1)
{continue;}
if(indexes[i]===-1){t=0;continue;}
if(hasColorInd){x3dom.debug.assert(colorInd[i]!=-1);}
switch(t)
{case 0:p0=+indexes[i];if(hasColorInd&&colPerVert){c0=+colorInd[i];}
else{c0=p0;}
t=1;break;case 1:p1=+indexes[i];if(hasColorInd&&colPerVert){c1=+colorInd[i];}
else if(hasColorInd&&!colPerVert){c1=+colorInd[lineCnt];}
else{c1=p1;}
this._mesh._indices[0].push(cnt++,cnt++);this._mesh._positions[0].push(positions[p0].x);this._mesh._positions[0].push(positions[p0].y);this._mesh._positions[0].push(positions[p0].z);this._mesh._positions[0].push(positions[p1].x);this._mesh._positions[0].push(positions[p1].y);this._mesh._positions[0].push(positions[p1].z);if(hasColor){if(!colPerVert){c0=c1;}
this._mesh._colors[0].push(colors[c0].r);this._mesh._colors[0].push(colors[c0].g);this._mesh._colors[0].push(colors[c0].b);this._mesh._colors[0].push(colors[c1].r);this._mesh._colors[0].push(colors[c1].g);this._mesh._colors[0].push(colors[c1].b);}
t=2;lineCnt++;break;case 2:p0=p1;c0=c1;p1=+indexes[i];if(hasColorInd&&colPerVert){c1=+colorInd[i];}
else if(hasColorInd&&!colPerVert){c1=+colorInd[lineCnt];}
else{c1=p1;}
this._mesh._indices[0].push(cnt++,cnt++);this._mesh._positions[0].push(positions[p0].x);this._mesh._positions[0].push(positions[p0].y);this._mesh._positions[0].push(positions[p0].z);this._mesh._positions[0].push(positions[p1].x);this._mesh._positions[0].push(positions[p1].y);this._mesh._positions[0].push(positions[p1].z);if(hasColor){if(!colPerVert){c0=c1;}
this._mesh._colors[0].push(colors[c0].r);this._mesh._colors[0].push(colors[c0].g);this._mesh._colors[0].push(colors[c0].b);this._mesh._colors[0].push(colors[c1].r);this._mesh._colors[0].push(colors[c1].g);this._mesh._colors[0].push(colors[c1].b);}
lineCnt++;break;default:}}
if(positions.length>x3dom.Utils.maxIndexableCoords)
this._mesh.splitMesh(2);}
else
{var n=indexes.length;t=0;for(i=0;i<n;++i)
{if(indexes[i]==-1){t=0;continue;}
switch(t){case 0:p0=+indexes[i];t=1;break;case 1:p1=+indexes[i];t=2;this._mesh._indices[0].push(p0,p1);break;case 2:p0=p1;p1=+indexes[i];this._mesh._indices[0].push(p0,p1);break;}}
this._mesh._positions[0]=positions.toGL();if(hasColor){this._mesh._colors[0]=colors.toGL();this._mesh._numColComponents=numColComponents;}}
this.invalidateVolume();this._mesh._numCoords=0;for(i=0;i<this._mesh._indices.length;i++){this._mesh._numCoords+=this._mesh._positions[i].length/3;}
var time1=new Date().getTime()-time0;},nodeChanged:function()
{this._buildGeometry();},fieldChanged:function(fieldName)
{var pnts=null;if(fieldName=="coord")
{this._buildGeometry();Array.forEach(this._parentNodes,function(node){node._dirty.positions=true;node.invalidateVolume();});}
else if(fieldName=="color")
{this._buildGeometry();Array.forEach(this._parentNodes,function(node){node._dirty.colors=true;});}
else if(fieldName=="coordIndex"){this._buildGeometry();Array.forEach(this._parentNodes,function(node){node._dirty.indexes=true;node.invalidateVolume();});}
else if(fieldName=="colorIndex"){this._buildGeometry();Array.forEach(this._parentNodes,function(node){node._dirty.colors=true;node.invalidateVolume();});}}}));x3dom.registerNodeType("IndexedTriangleSet","Rendering",defineClass(x3dom.nodeTypes.X3DComposedGeometryNode,function(ctx){x3dom.nodeTypes.IndexedTriangleSet.superClass.call(this,ctx);this.addField_MFInt32(ctx,'index',[]);},{nodeChanged:function()
{var time0=new Date().getTime();this.handleAttribs();var colPerVert=this._vf.colorPerVertex;var normPerVert=this._vf.normalPerVertex;var ccw=this._vf.ccw;var indexes=this._vf.index;var hasNormal=false,hasTexCoord=false,hasColor=false;var positions,normals,texCoords,colors;var coordNode=this._cf.coord.node;x3dom.debug.assert(coordNode);positions=coordNode._vf.point;var normalNode=this._cf.normal.node;if(normalNode){hasNormal=true;normals=normalNode._vf.vector;}
else{hasNormal=false;}
var texMode="",numTexComponents=2;var texCoordNode=this._cf.texCoord.node;if(x3dom.isa(texCoordNode,x3dom.nodeTypes.MultiTextureCoordinate)){if(texCoordNode._cf.texCoord.nodes.length)
texCoordNode=texCoordNode._cf.texCoord.nodes[0];}
if(texCoordNode){if(texCoordNode._vf.point){hasTexCoord=true;texCoords=texCoordNode._vf.point;if(x3dom.isa(texCoordNode,x3dom.nodeTypes.TextureCoordinate3D)){numTexComponents=3;}}
else if(texCoordNode._vf.mode){texMode=texCoordNode._vf.mode;}}
else{hasTexCoord=false;}
var numColComponents=3;var colorNode=this._cf.color.node;if(colorNode){hasColor=true;colors=colorNode._vf.color;if(x3dom.isa(colorNode,x3dom.nodeTypes.ColorRGBA)){numColComponents=4;}}
else{hasColor=false;}
this._mesh._indices[0]=[];this._mesh._positions[0]=[];this._mesh._normals[0]=[];this._mesh._texCoords[0]=[];this._mesh._colors[0]=[];var i,t,cnt,faceCnt,posMax;var p0,p1,p2,n0,n1,n2,t0,t1,t2,c0,c1,c2;while(positions.length%3>0){positions.push(positions.length-1);}
posMax=positions.length;if(!normPerVert||!colPerVert||posMax>x3dom.Utils.maxIndexableCoords)
{t=0;cnt=0;faceCnt=0;this._mesh._multiIndIndices=[];this._mesh._posSize=positions.length;for(i=0;i<indexes.length;++i)
{if((i>0)&&(i%3===0)){t=0;faceCnt++;}
switch(t)
{case 0:p0=+indexes[i];if(normPerVert){n0=p0;}else if(!normPerVert){n0=faceCnt;}
t0=p0;if(colPerVert){c0=p0;}else if(!colPerVert){c0=faceCnt;}
t=1;break;case 1:p1=+indexes[i];if(normPerVert){n1=p1;}else if(!normPerVert){n1=faceCnt;}
t1=p1;if(colPerVert){c1=p1;}else if(!colPerVert){c1=faceCnt;}
t=2;break;case 2:p2=+indexes[i];if(normPerVert){n2=p2;}else if(!normPerVert){n2=faceCnt;}
t2=p2;if(colPerVert){c2=p2;}else if(!colPerVert){c2=faceCnt;}
t=3;this._mesh._indices[0].push(cnt++,cnt++,cnt++);this._mesh._positions[0].push(positions[p0].x);this._mesh._positions[0].push(positions[p0].y);this._mesh._positions[0].push(positions[p0].z);this._mesh._positions[0].push(positions[p1].x);this._mesh._positions[0].push(positions[p1].y);this._mesh._positions[0].push(positions[p1].z);this._mesh._positions[0].push(positions[p2].x);this._mesh._positions[0].push(positions[p2].y);this._mesh._positions[0].push(positions[p2].z);if(hasNormal){this._mesh._normals[0].push(normals[n0].x);this._mesh._normals[0].push(normals[n0].y);this._mesh._normals[0].push(normals[n0].z);this._mesh._normals[0].push(normals[n1].x);this._mesh._normals[0].push(normals[n1].y);this._mesh._normals[0].push(normals[n1].z);this._mesh._normals[0].push(normals[n2].x);this._mesh._normals[0].push(normals[n2].y);this._mesh._normals[0].push(normals[n2].z);}
else{this._mesh._multiIndIndices.push(p0,p1,p2);}
if(hasColor){this._mesh._colors[0].push(colors[c0].r);this._mesh._colors[0].push(colors[c0].g);this._mesh._colors[0].push(colors[c0].b);if(numColComponents===4){this._mesh._colors[0].push(colors[c0].a);}
this._mesh._colors[0].push(colors[c1].r);this._mesh._colors[0].push(colors[c1].g);this._mesh._colors[0].push(colors[c1].b);if(numColComponents===4){this._mesh._colors[0].push(colors[c1].a);}
this._mesh._colors[0].push(colors[c2].r);this._mesh._colors[0].push(colors[c2].g);this._mesh._colors[0].push(colors[c2].b);if(numColComponents===4){this._mesh._colors[0].push(colors[c2].a);}}
if(hasTexCoord){this._mesh._texCoords[0].push(texCoords[t0].x);this._mesh._texCoords[0].push(texCoords[t0].y);if(numTexComponents===3){this._mesh._texCoords[0].push(texCoords[t0].z);}
this._mesh._texCoords[0].push(texCoords[t1].x);this._mesh._texCoords[0].push(texCoords[t1].y);if(numTexComponents===3){this._mesh._texCoords[0].push(texCoords[t1].z);}
this._mesh._texCoords[0].push(texCoords[t2].x);this._mesh._texCoords[0].push(texCoords[t2].y);if(numTexComponents===3){this._mesh._texCoords[0].push(texCoords[t2].z);}}
break;default:}}
if(!hasNormal){this._mesh.calcNormals(normPerVert?Math.PI:0);}
if(!hasTexCoord){this._mesh.calcTexCoords(texMode);}
this._mesh.splitMesh();}
else
{faceCnt=0;for(i=0;i<indexes.length;i++)
{if((i>0)&&(i%3===0)){faceCnt++;}
this._mesh._indices[0].push(indexes[i]);if(!normPerVert&&hasNormal){this._mesh._normals[0].push(normals[faceCnt].x);this._mesh._normals[0].push(normals[faceCnt].y);this._mesh._normals[0].push(normals[faceCnt].z);}
if(!colPerVert&&hasColor){this._mesh._colors[0].push(colors[faceCnt].r);this._mesh._colors[0].push(colors[faceCnt].g);this._mesh._colors[0].push(colors[faceCnt].b);if(numColComponents===4){this._mesh._colors[0].push(colors[faceCnt].a);}}}
this._mesh._positions[0]=positions.toGL();if(hasNormal){this._mesh._normals[0]=normals.toGL();}
else{this._mesh.calcNormals(normPerVert?Math.PI:0,ccw);}
if(hasTexCoord){this._mesh._texCoords[0]=texCoords.toGL();this._mesh._numTexComponents=numTexComponents;}
else{this._mesh.calcTexCoords(texMode);}
if(hasColor&&colPerVert){this._mesh._colors[0]=colors.toGL();this._mesh._numColComponents=numColComponents;}}
this.invalidateVolume();this._mesh._numFaces=0;this._mesh._numCoords=0;for(i=0;i<this._mesh._indices.length;i++){this._mesh._numFaces+=this._mesh._indices[i].length/3;this._mesh._numCoords+=this._mesh._positions[i].length/3;}
var time1=new Date().getTime()-time0;},fieldChanged:function(fieldName)
{var pnts=this._cf.coord.node._vf.point;if(pnts.length>x3dom.Utils.maxIndexableCoords)
{x3dom.debug.logWarning("IndexedTriangleSet: fieldChanged with "+"too many coordinates not yet implemented!");return;}
if(fieldName=="coord")
{this._mesh._positions[0]=pnts.toGL();this.invalidateVolume();Array.forEach(this._parentNodes,function(node){node._dirty.positions=true;node.invalidateVolume();});}
else if(fieldName=="color")
{pnts=this._cf.color.node._vf.color;if(this._vf.colorPerVertex){this._mesh._colors[0]=pnts.toGL();}else if(!this._vf.colorPerVertex){var faceCnt=0;var numColComponents=3;if(x3dom.isa(this._cf.color.node,x3dom.nodeTypes.ColorRGBA)){numColComponents=4;}
this._mesh._colors[0]=[];var indexes=this._vf.index;for(var i=0;i<indexes.length;++i)
{if((i>0)&&(i%3===0)){faceCnt++;}
this._mesh._colors[0].push(pnts[faceCnt].r);this._mesh._colors[0].push(pnts[faceCnt].g);this._mesh._colors[0].push(pnts[faceCnt].b);if(numColComponents===4){this._mesh._colors[0].push(pnts[faceCnt].a);}}}
Array.forEach(this._parentNodes,function(node){node._dirty.colors=true;});}
else if(fieldName=="normal")
{pnts=this._cf.normal.node._vf.vector;if(this._vf.normalPerVertex){this._mesh._normals[0]=pnts.toGL();}else if(!this._vf.normalPerVertex){var indexes=this._vf.index;this._mesh._normals[0]=[];var faceCnt=0;for(var i=0;i<indexes.length;++i)
{if((i>0)&&(i%3===0)){faceCnt++;}
this._mesh._normals[0].push(pnts[faceCnt].x);this._mesh._normals[0].push(pnts[faceCnt].y);this._mesh._normals[0].push(pnts[faceCnt].z);}}
Array.forEach(this._parentNodes,function(node){node._dirty.normals=true;});}
else if(fieldName=="texCoord")
{var texCoordNode=this._cf.texCoord.node;if(x3dom.isa(texCoordNode,x3dom.nodeTypes.MultiTextureCoordinate)){if(texCoordNode._cf.texCoord.nodes.length)
texCoordNode=texCoordNode._cf.texCoord.nodes[0];}
pnts=texCoordNode._vf.point;this._mesh._texCoords[0]=pnts.toGL();Array.forEach(this._parentNodes,function(node){node._dirty.texcoords=true;});}}}));x3dom.registerNodeType("IndexedTriangleStripSet","Rendering",defineClass(x3dom.nodeTypes.X3DComposedGeometryNode,function(ctx){x3dom.nodeTypes.IndexedTriangleStripSet.superClass.call(this,ctx);this.addField_MFInt32(ctx,'index',[]);this._hasIndexOffset=false;this._indexOffset=null;},{hasIndexOffset:function(){return this._hasIndexOffset;},nodeChanged:function()
{this.handleAttribs();var hasNormal=false,hasTexCoord=false,hasColor=false;var colPerVert=this._vf.colorPerVertex;var normPerVert=this._vf.normalPerVertex;var indexes=this._vf.index;if(indexes.length&&indexes[indexes.length-1]!=-1)
{indexes.push(-1);}
var positions,normals,texCoords,colors;var coordNode=this._cf.coord.node;x3dom.debug.assert(coordNode);positions=coordNode._vf.point;var normalNode=this._cf.normal.node;if(normalNode){hasNormal=true;normals=normalNode._vf.vector;}
else{hasNormal=false;}
var texMode="",numTexComponents=2;var texCoordNode=this._cf.texCoord.node;if(x3dom.isa(texCoordNode,x3dom.nodeTypes.MultiTextureCoordinate)){if(texCoordNode._cf.texCoord.nodes.length)
texCoordNode=texCoordNode._cf.texCoord.nodes[0];}
if(texCoordNode){if(texCoordNode._vf.point){hasTexCoord=true;texCoords=texCoordNode._vf.point;if(x3dom.isa(texCoordNode,x3dom.nodeTypes.TextureCoordinate3D)){numTexComponents=3;}}
else if(texCoordNode._vf.mode){texMode=texCoordNode._vf.mode;}}
else{hasTexCoord=false;}
this._mesh._numTexComponents=numTexComponents;var numColComponents=3;var colorNode=this._cf.color.node;if(colorNode){hasColor=true;colors=colorNode._vf.color;if(x3dom.isa(colorNode,x3dom.nodeTypes.ColorRGBA)){numColComponents=4;}}
else{hasColor=false;}
this._mesh._numColComponents=numColComponents;this._mesh._indices[0]=[];this._mesh._positions[0]=[];this._mesh._normals[0]=[];this._mesh._texCoords[0]=[];this._mesh._colors[0]=[];this.invalidateVolume();this._mesh._numFaces=0;this._mesh._numCoords=0;var faceCnt=0,cnt=0;if(hasNormal&&positions.length<=x3dom.Utils.maxIndexableCoords)
{this._hasIndexOffset=true;this._indexOffset=[];this._mesh._primType='TRIANGLESTRIP';var indexOffset=[0];for(i=0;i<indexes.length;i++)
{if(indexes[i]==-1){faceCnt++;indexOffset.push(this._mesh._indices[0].length);}
else{this._mesh._indices[0].push(+indexes[i]);if(!normPerVert){this._mesh._normals[0].push(normals[faceCnt].x);this._mesh._normals[0].push(normals[faceCnt].y);this._mesh._normals[0].push(normals[faceCnt].z);}
if(!colPerVert){this._mesh._colors[0].push(colors[faceCnt].r);this._mesh._colors[0].push(colors[faceCnt].g);this._mesh._colors[0].push(colors[faceCnt].b);if(numColComponents===4){this._mesh._colors[0].push(colors[faceCnt].a);}}}}
this._mesh._positions[0]=positions.toGL();if(normPerVert){this._mesh._normals[0]=normals.toGL();}
if(hasTexCoord){this._mesh._texCoords[0]=texCoords.toGL();this._mesh._numTexComponents=numTexComponents;}
else{x3dom.debug.logWarning("IndexedTriangleStripSet: no texCoords given and won't calculate!");}
if(hasColor){if(colPerVert){this._mesh._colors[0]=colors.toGL();}
this._mesh._numColComponents=numColComponents;}
for(i=1;i<indexOffset.length;i++){var triCnt=indexOffset[i]-indexOffset[i-1];this._indexOffset.push({count:triCnt,offset:2*indexOffset[i-1]});this._mesh._numFaces+=(triCnt-2);}
this._mesh._numCoords=this._mesh._positions[0].length/3;}
else
{this._hasIndexOffset=false;var p1,p2,p3,n1,n2,n3,t1,t2,t3,c1,c2,c3;var swapOrder=false;for(var i=1;i<indexes.length-2;++i)
{if(indexes[i+1]==-1){i=i+2;faceCnt++;continue;}
if(swapOrder){p1=indexes[i];p2=indexes[i-1];p3=indexes[i+1];}
else{p1=indexes[i-1];p2=indexes[i];p3=indexes[i+1];}
swapOrder=!swapOrder;if(normPerVert){n1=p1;n2=p2;n3=p3;}else if(!normPerVert){n1=n2=n3=faceCnt;}
t1=p1;t2=p2;t3=p3;if(colPerVert){c1=p1;c2=p2;c3=p3;}else if(!colPerVert){c1=c2=c3=faceCnt;}
this._mesh._indices[0].push(cnt++,cnt++,cnt++);this._mesh._positions[0].push(positions[p1].x);this._mesh._positions[0].push(positions[p1].y);this._mesh._positions[0].push(positions[p1].z);this._mesh._positions[0].push(positions[p2].x);this._mesh._positions[0].push(positions[p2].y);this._mesh._positions[0].push(positions[p2].z);this._mesh._positions[0].push(positions[p3].x);this._mesh._positions[0].push(positions[p3].y);this._mesh._positions[0].push(positions[p3].z);if(hasNormal){this._mesh._normals[0].push(normals[n1].x);this._mesh._normals[0].push(normals[n1].y);this._mesh._normals[0].push(normals[n1].z);this._mesh._normals[0].push(normals[n2].x);this._mesh._normals[0].push(normals[n2].y);this._mesh._normals[0].push(normals[n2].z);this._mesh._normals[0].push(normals[n3].x);this._mesh._normals[0].push(normals[n3].y);this._mesh._normals[0].push(normals[n3].z);}
if(hasColor){this._mesh._colors[0].push(colors[c1].r);this._mesh._colors[0].push(colors[c1].g);this._mesh._colors[0].push(colors[c1].b);if(numColComponents===4){this._mesh._colors[0].push(colors[c1].a);}
this._mesh._colors[0].push(colors[c2].r);this._mesh._colors[0].push(colors[c2].g);this._mesh._colors[0].push(colors[c2].b);if(numColComponents===4){this._mesh._colors[0].push(colors[c2].a);}
this._mesh._colors[0].push(colors[c3].r);this._mesh._colors[0].push(colors[c3].g);this._mesh._colors[0].push(colors[c3].b);if(numColComponents===4){this._mesh._colors[0].push(colors[c3].a);}}
if(hasTexCoord){this._mesh._texCoords[0].push(texCoords[t1].x);this._mesh._texCoords[0].push(texCoords[t1].y);if(numTexComponents===3){this._mesh._texCoords[0].push(texCoords[t1].z);}
this._mesh._texCoords[0].push(texCoords[t2].x);this._mesh._texCoords[0].push(texCoords[t2].y);if(numTexComponents===3){this._mesh._texCoords[0].push(texCoords[t2].z);}
this._mesh._texCoords[0].push(texCoords[t3].x);this._mesh._texCoords[0].push(texCoords[t3].y);if(numTexComponents===3){this._mesh._texCoords[0].push(texCoords[t3].z);}}}
if(!hasNormal){this._mesh.calcNormals(Math.PI);}
if(!hasTexCoord){this._mesh.calcTexCoords(texMode);}
this._mesh.splitMesh();this.invalidateVolume();for(i=0;i<this._mesh._indices.length;i++){this._mesh._numFaces+=this._mesh._indices[i].length/3;this._mesh._numCoords+=this._mesh._positions[i].length/3;}}},fieldChanged:function(fieldName)
{if(fieldName!="coord"&&fieldName!="normal"&&fieldName!="texCoord"&&fieldName!="color")
{x3dom.debug.logWarning("IndexedTriangleStripSet: fieldChanged for "+
fieldName+" not yet implemented!");return;}
var pnts=this._cf.coord.node._vf.point;if((this._cf.normal.node===null)||(pnts.length>x3dom.Utils.maxIndexableCoords))
{if(fieldName=="coord"){this._mesh._positions[0]=[];this._mesh._indices[0]=[];this._mesh._normals[0]=[];this._mesh._texCoords[0]=[];var hasNormal=false,hasTexCoord=false,hasColor=false;var colPerVert=this._vf.colorPerVertex;var normPerVert=this._vf.normalPerVertex;var indexes=this._vf.index;var positions,normals,texCoords,colors;var coordNode=this._cf.coord.node;x3dom.debug.assert(coordNode);positions=coordNode._vf.point;var normalNode=this._cf.normal.node;if(normalNode){hasNormal=true;normals=normalNode._vf.vector;}
else{hasNormal=false;}
var texMode="",numTexComponents=2;var texCoordNode=this._cf.texCoord.node;if(x3dom.isa(texCoordNode,x3dom.nodeTypes.MultiTextureCoordinate)){if(texCoordNode._cf.texCoord.nodes.length)
texCoordNode=texCoordNode._cf.texCoord.nodes[0];}
if(texCoordNode){if(texCoordNode._vf.point){hasTexCoord=true;texCoords=texCoordNode._vf.point;if(x3dom.isa(texCoordNode,x3dom.nodeTypes.TextureCoordinate3D)){numTexComponents=3;}}
else if(texCoordNode._vf.mode){texMode=texCoordNode._vf.mode;}}
else{hasTexCoord=false;}
this._mesh._numTexComponents=numTexComponents;var numColComponents=3;var colorNode=this._cf.color.node;if(colorNode){hasColor=true;colors=colorNode._vf.color;if(x3dom.isa(colorNode,x3dom.nodeTypes.ColorRGBA)){numColComponents=4;}}
else{hasColor=false;}
this._mesh._numColComponents=numColComponents;this._mesh._indices[0]=[];this._mesh._positions[0]=[];this._mesh._normals[0]=[];this._mesh._texCoords[0]=[];this._mesh._colors[0]=[];var faceCnt=0,cnt=0;var p1,p2,p3,n1,n2,n3,t1,t2,t3,c1,c2,c3;var swapOrder=false;if(hasNormal||hasTexCoord||hasColor){for(var i=1;i<indexes.length-2;++i)
{if(indexes[i+1]==-1){i=i+2;faceCnt++;continue;}
if(swapOrder){p1=indexes[i];p2=indexes[i-1];p3=indexes[i+1];}
else{p1=indexes[i-1];p2=indexes[i];p3=indexes[i+1];}
swapOrder=!swapOrder;if(normPerVert){n1=p1;n2=p2;n3=p3;}else if(!normPerVert){n1=n2=n3=faceCnt;}
t1=p1;t2=p2;t3=p3;if(colPerVert){c1=p1;c2=p2;c3=p3;}else if(!colPerVert){c1=c2=c3=faceCnt;}
this._mesh._indices[0].push(cnt++,cnt++,cnt++);this._mesh._positions[0].push(positions[p1].x);this._mesh._positions[0].push(positions[p1].y);this._mesh._positions[0].push(positions[p1].z);this._mesh._positions[0].push(positions[p2].x);this._mesh._positions[0].push(positions[p2].y);this._mesh._positions[0].push(positions[p2].z);this._mesh._positions[0].push(positions[p3].x);this._mesh._positions[0].push(positions[p3].y);this._mesh._positions[0].push(positions[p3].z);if(hasNormal){this._mesh._normals[0].push(normals[n1].x);this._mesh._normals[0].push(normals[n1].y);this._mesh._normals[0].push(normals[n1].z);this._mesh._normals[0].push(normals[n2].x);this._mesh._normals[0].push(normals[n2].y);this._mesh._normals[0].push(normals[n2].z);this._mesh._normals[0].push(normals[n3].x);this._mesh._normals[0].push(normals[n3].y);this._mesh._normals[0].push(normals[n3].z);}
if(hasColor){this._mesh._colors[0].push(colors[c1].r);this._mesh._colors[0].push(colors[c1].g);this._mesh._colors[0].push(colors[c1].b);if(numColComponents===4){this._mesh._colors[0].push(colors[c1].a);}
this._mesh._colors[0].push(colors[c2].r);this._mesh._colors[0].push(colors[c2].g);this._mesh._colors[0].push(colors[c2].b);if(numColComponents===4){this._mesh._colors[0].push(colors[c2].a);}
this._mesh._colors[0].push(colors[c3].r);this._mesh._colors[0].push(colors[c3].g);this._mesh._colors[0].push(colors[c3].b);if(numColComponents===4){this._mesh._colors[0].push(colors[c3].a);}}
if(hasTexCoord){this._mesh._texCoords[0].push(texCoords[t1].x);this._mesh._texCoords[0].push(texCoords[t1].y);if(numTexComponents===3){this._mesh._texCoords[0].push(texCoords[t1].z);}
this._mesh._texCoords[0].push(texCoords[t2].x);this._mesh._texCoords[0].push(texCoords[t2].y);if(numTexComponents===3){this._mesh._texCoords[0].push(texCoords[t2].z);}
this._mesh._texCoords[0].push(texCoords[t3].x);this._mesh._texCoords[0].push(texCoords[t3].y);if(numTexComponents===3){this._mesh._texCoords[0].push(texCoords[t3].z);}}}
if(!hasNormal){this._mesh.calcNormals(Math.PI);}
if(!hasTexCoord){this._mesh.calcTexCoords(texMode);}
this._mesh.splitMesh();}else{var swapOrder=false;for(var i=1;i<indexes.length;++i)
{if(indexes[i+1]==-1){i=i+2;continue;}
if(swapOrder){this._mesh._indices[0].push(indexes[i]);this._mesh._indices[0].push(indexes[i-1]);this._mesh._indices[0].push(indexes[i+1]);}
else{this._mesh._indices[0].push(indexes[i-1]);this._mesh._indices[0].push(indexes[i]);this._mesh._indices[0].push(indexes[i+1]);}
swapOrder=!swapOrder;}
this._mesh._positions[0]=positions.toGL();if(hasNormal){this._mesh._normals[0]=normals.toGL();}
else{this._mesh.calcNormals(Math.PI);}
if(hasTexCoord){this._mesh._texCoords[0]=texCoords.toGL();this._mesh._numTexComponents=numTexComponents;}
else{this._mesh.calcTexCoords(texMode);}
if(hasColor){this._mesh._colors[0]=colors.toGL();this._mesh._numColComponents=numColComponents;}}
this.invalidateVolume();this._mesh._numFaces=0;this._mesh._numCoords=0;for(i=0;i<this._mesh._indices.length;i++){this._mesh._numFaces+=this._mesh._indices[i].length/3;this._mesh._numCoords+=this._mesh._positions[i].length/3;}
Array.forEach(this._parentNodes,function(node){node.setAllDirty();node.invalidateVolume();});}
else if(fieldName=="color"){var col=this._cf.color.node._vf.color;var faceCnt=0;var c1=c2=c3=0;var numColComponents=3;if(x3dom.isa(this._cf.color.node,x3dom.nodeTypes.ColorRGBA)){numColComponents=4;}
this._mesh._colors[0]=[];var indexes=this._vf.index;var swapOrder=false;for(i=1;i<indexes.length-2;++i)
{if(indexes[i+1]==-1){i=i+2;faceCnt++;continue;}
if(this._vf.colorPerVertex){if(swapOrder){c1=indexes[i];c2=indexes[i-1];c3=indexes[i+1];}
else{c1=indexes[i-1];c2=indexes[i];c3=indexes[i+1];}
swapOrder=!swapOrder;}else if(!this._vf.colorPerVertex){c1=c2=c3=faceCnt;}
this._mesh._colors[0].push(col[c1].r);this._mesh._colors[0].push(col[c1].g);this._mesh._colors[0].push(col[c1].b);if(numColComponents===4){this._mesh._colors[0].push(col[c1].a);}
this._mesh._colors[0].push(col[c2].r);this._mesh._colors[0].push(col[c2].g);this._mesh._colors[0].push(col[c2].b);if(numColComponents===4){this._mesh._colors[0].push(col[c2].a);}
this._mesh._colors[0].push(col[c3].r);this._mesh._colors[0].push(col[c3].g);this._mesh._colors[0].push(col[c3].b);if(numColComponents===4){this._mesh._colors[0].push(col[c3].a);}}
Array.forEach(this._parentNodes,function(node){node._dirty.colors=true;});}
else if(fieldName=="normal"){var nor=this._cf.normal.node._vf.vector;var faceCnt=0;var n1=n2=n3=0;this._mesh._normals[0]=[];var indexes=this._vf.index;var swapOrder=false;for(i=1;i<indexes.length-2;++i)
{if(indexes[i+1]==-1){i=i+2;faceCnt++;continue;}
if(this._vf.normalPerVertex){if(swapOrder){n1=indexes[i];n2=indexes[i-1];n3=indexes[i+1];}
else{n1=indexes[i-1];n2=indexes[i];n3=indexes[i+1];}
swapOrder=!swapOrder;}else if(!this._vf.normalPerVertex){n1=n2=n3=faceCnt;}
this._mesh._normals[0].push(nor[n1].x);this._mesh._normals[0].push(nor[n1].y);this._mesh._normals[0].push(nor[n1].z);this._mesh._normals[0].push(nor[n2].x);this._mesh._normals[0].push(nor[n2].y);this._mesh._normals[0].push(nor[n2].z);this._mesh._normals[0].push(nor[n3].x);this._mesh._normals[0].push(nor[n3].y);this._mesh._normals[0].push(nor[n3].z);}
Array.forEach(this._parentNodes,function(node){node._dirty.normals=true;});}
else if(fieldName=="texCoord"){var texCoordNode=this._cf.texCoord.node;if(x3dom.isa(texCoordNode,x3dom.nodeTypes.MultiTextureCoordinate)){if(texCoordNode._cf.texCoord.nodes.length)
texCoordNode=texCoordNode._cf.texCoord.nodes[0];}
var tex=texCoordNode._vf.point;var t1=t2=t3=0;var numTexComponents=2;if(x3dom.isa(texCoordNode,x3dom.nodeTypes.TextureCoordinate3D)){numTexComponents=3;}
this._mesh._texCoords[0]=[];var indexes=this._vf.index;var swapOrder=false;for(i=1;i<indexes.length-2;++i)
{if(indexes[i+1]==-1){i=i+2;continue;}
if(swapOrder){t1=indexes[i];t2=indexes[i-1];t3=indexes[i+1];}
else{t1=indexes[i-1];t2=indexes[i];t3=indexes[i+1];}
swapOrder=!swapOrder;this._mesh._texCoords[0].push(tex[t1].x);this._mesh._texCoords[0].push(tex[t1].y);if(numTexComponents===3){this._mesh._texCoords[0].push(tex[t1].z);}
this._mesh._texCoords[0].push(tex[t2].x);this._mesh._texCoords[0].push(tex[t2].y);if(numTexComponents===3){this._mesh._texCoords[0].tex(col[t2].z);}
this._mesh._texCoords[0].push(tex[t3].x);this._mesh._texCoords[0].push(tex[t3].y);if(numTexComponents===3){this._mesh._texCoords[0].push(tex[t3].z);}}
Array.forEach(this._parentNodes,function(node){node._dirty.texcoords=true;});}}
else
{if(fieldName=="coord")
{this._mesh._positions[0]=pnts.toGL();this.invalidateVolume();Array.forEach(this._parentNodes,function(node){node._dirty.positions=true;node.invalidateVolume();});}
else if(fieldName=="color")
{pnts=this._cf.color.node._vf.color;if(this._vf.colorPerVertex){this._mesh._colors[0]=pnts.toGL();}else if(!this._vf.colorPerVertex){var faceCnt=0;var numColComponents=3;if(x3dom.isa(this._cf.color.node,x3dom.nodeTypes.ColorRGBA)){numColComponents=4;}
this._mesh._colors[0]=[];var indexes=this._vf.index;for(i=0;i<indexes.length;++i)
{if(indexes[i]==-1){faceCnt++;continue;}
this._mesh._colors[0].push(pnts[faceCnt].r);this._mesh._colors[0].push(pnts[faceCnt].g);this._mesh._colors[0].push(pnts[faceCnt].b);if(numColComponents===4){this._mesh._colors[0].push(pnts[faceCnt].a);}}}
Array.forEach(this._parentNodes,function(node){node._dirty.colors=true;});}
else if(fieldName=="normal")
{pnts=this._cf.normal.node._vf.vector;if(this._vf.normalPerVertex){this._mesh._normals[0]=pnts.toGL();}else if(!this._vf.normalPerVertex){var indexes=this._vf.index;this._mesh._normals[0]=[];var faceCnt=0;for(i=0;i<indexes.length;++i)
{if(indexes[i]==-1){faceCnt++;continue;}
this._mesh._normals[0].push(pnts[faceCnt].x);this._mesh._normals[0].push(pnts[faceCnt].y);this._mesh._normals[0].push(pnts[faceCnt].z);}}
Array.forEach(this._parentNodes,function(node){node._dirty.normals=true;});}
else if(fieldName=="texCoord")
{var texCoordNode=this._cf.texCoord.node;if(x3dom.isa(texCoordNode,x3dom.nodeTypes.MultiTextureCoordinate)){if(texCoordNode._cf.texCoord.nodes.length)
texCoordNode=texCoordNode._cf.texCoord.nodes[0];}
pnts=texCoordNode._vf.point;this._mesh._texCoords[0]=pnts.toGL();Array.forEach(this._parentNodes,function(node){node._dirty.texcoords=true;});}}}}));x3dom.registerNodeType("TriangleSet","Rendering",defineClass(x3dom.nodeTypes.X3DComposedGeometryNode,function(ctx){x3dom.nodeTypes.TriangleSet.superClass.call(this,ctx);},{_buildGeometry:function()
{var colPerVert=this._vf.colorPerVertex;var normPerVert=this._vf.normalPerVertex;var ccw=this._vf.ccw;var hasNormal=false,hasTexCoord=false,hasColor=false;var positions,normals,texCoords,colors;var coordNode=this._cf.coord.node;x3dom.debug.assert(coordNode);if(!coordNode||coordNode._vf.point.length<3)
{this._vf.render=false;return;}
positions=coordNode._vf.point;var normalNode=this._cf.normal.node;if(normalNode){hasNormal=true;normals=normalNode._vf.vector;}
else{hasNormal=false;}
var texMode="",numTexComponents=2;var texCoordNode=this._cf.texCoord.node;if(x3dom.isa(texCoordNode,x3dom.nodeTypes.MultiTextureCoordinate)){if(texCoordNode._cf.texCoord.nodes.length)
texCoordNode=texCoordNode._cf.texCoord.nodes[0];}
if(texCoordNode){if(texCoordNode._vf.point){hasTexCoord=true;texCoords=texCoordNode._vf.point;if(x3dom.isa(texCoordNode,x3dom.nodeTypes.TextureCoordinate3D)){numTexComponents=3;}}
else if(texCoordNode._vf.mode){texMode=texCoordNode._vf.mode;}}
else{hasTexCoord=false;}
var numColComponents=3;var colorNode=this._cf.color.node;if(colorNode){hasColor=true;colors=colorNode._vf.color;if(x3dom.isa(colorNode,x3dom.nodeTypes.ColorRGBA)){numColComponents=4;}}
else{hasColor=false;}
while(positions.length%3>0){positions.pop();}
this._mesh._indices[0]=new Array(positions.length);this._mesh._positions[0]=[];this._mesh._normals[0]=[];this._mesh._texCoords[0]=[];this._mesh._colors[0]=[];var posMax=positions.length/3;var faceCnt,i=0;for(faceCnt=0;faceCnt<posMax;faceCnt++)
{this._mesh._indices[0][i]=i++;this._mesh._indices[0][i]=i++;this._mesh._indices[0][i]=i++;if(!normPerVert&&hasNormal){this._mesh._normals[0].push(normals[faceCnt].x);this._mesh._normals[0].push(normals[faceCnt].y);this._mesh._normals[0].push(normals[faceCnt].z);}
if(!colPerVert&&hasColor){this._mesh._colors[0].push(colors[faceCnt].r);this._mesh._colors[0].push(colors[faceCnt].g);this._mesh._colors[0].push(colors[faceCnt].b);if(numColComponents===4){this._mesh._colors[0].push(colors[faceCnt].a);}}}
this._mesh._positions[0]=positions.toGL();if(hasNormal){this._mesh._normals[0]=normals.toGL();}
else{this._mesh.calcNormals(normPerVert?Math.PI:0,ccw);}
if(hasTexCoord){this._mesh._texCoords[0]=texCoords.toGL();this._mesh._numTexComponents=numTexComponents;}
else{this._mesh.calcTexCoords(texMode);}
if(hasColor&&colPerVert){this._mesh._colors[0]=colors.toGL();this._mesh._numColComponents=numColComponents;}
this._mesh._numFaces=posMax;this._mesh._numCoords=positions.length;this.invalidateVolume();},nodeChanged:function()
{this._buildGeometry();},fieldChanged:function(fieldName)
{if(fieldName=="coord")
{this._buildGeometry();Array.forEach(this._parentNodes,function(node){node._dirty.positions=true;node.invalidateVolume();});}
else if(fieldName=="color")
{this._buildGeometry();Array.forEach(this._parentNodes,function(node){node._dirty.colors=true;});}
else if(fieldName=="normal")
{this._buildGeometry();Array.forEach(this._parentNodes,function(node){node._dirty.normals=true;});}
else if(fieldName=="texCoord")
{this._buildGeometry();Array.forEach(this._parentNodes,function(node){node._dirty.texcoords=true;});}}}));x3dom.registerNodeType("X3DGeometricPropertyNode","Rendering",defineClass(x3dom.nodeTypes.X3DNode,function(ctx){x3dom.nodeTypes.X3DGeometricPropertyNode.superClass.call(this,ctx);}));x3dom.registerNodeType("X3DCoordinateNode","Rendering",defineClass(x3dom.nodeTypes.X3DGeometricPropertyNode,function(ctx){x3dom.nodeTypes.X3DCoordinateNode.superClass.call(this,ctx);},{fieldChanged:function(fieldName){if(fieldName==="coord"||fieldName==="point"){Array.forEach(this._parentNodes,function(node){node.fieldChanged("coord");});}},parentAdded:function(parent){if(parent._mesh&&parent._cf.coord.node!==this){parent.fieldChanged("coord");}}}));x3dom.registerNodeType("Coordinate","Rendering",defineClass(x3dom.nodeTypes.X3DCoordinateNode,function(ctx){x3dom.nodeTypes.Coordinate.superClass.call(this,ctx);this.addField_MFVec3f(ctx,'point',[]);},{getPoints:function(){return this._vf.point;}}));x3dom.registerNodeType("Normal","Rendering",defineClass(x3dom.nodeTypes.X3DGeometricPropertyNode,function(ctx){x3dom.nodeTypes.Normal.superClass.call(this,ctx);this.addField_MFVec3f(ctx,'vector',[]);},{fieldChanged:function(fieldName){if(fieldName==="normal"||fieldName==="vector"){Array.forEach(this._parentNodes,function(node){node.fieldChanged("normal");});}},parentAdded:function(parent){if(parent._mesh&&parent._cf.normal.node!==this){parent.fieldChanged("normal");}}}));x3dom.registerNodeType("X3DColorNode","Rendering",defineClass(x3dom.nodeTypes.X3DGeometricPropertyNode,function(ctx){x3dom.nodeTypes.X3DColorNode.superClass.call(this,ctx);},{fieldChanged:function(fieldName){if(fieldName==="color"){Array.forEach(this._parentNodes,function(node){node.fieldChanged("color");});}},parentAdded:function(parent){if(parent._mesh&&parent._cf.color.node!==this){parent.fieldChanged("color");}}}));x3dom.registerNodeType("Color","Rendering",defineClass(x3dom.nodeTypes.X3DColorNode,function(ctx){x3dom.nodeTypes.Color.superClass.call(this,ctx);this.addField_MFColor(ctx,'color',[]);}));x3dom.registerNodeType("ColorRGBA","Rendering",defineClass(x3dom.nodeTypes.X3DColorNode,function(ctx){x3dom.nodeTypes.ColorRGBA.superClass.call(this,ctx);this.addField_MFColorRGBA(ctx,'color',[]);}));x3dom.registerNodeType("ParticleSet","Rendering",defineClass(x3dom.nodeTypes.PointSet,function(ctx){x3dom.nodeTypes.ParticleSet.superClass.call(this,ctx);this.addField_SFString(ctx,'mode','ViewDirQuads');this.addField_SFString(ctx,'drawOrder','Any');this.addField_SFNode('normal',x3dom.nodeTypes.Normal);this.addField_MFVec3f(ctx,'size',[]);this.addField_MFInt32(ctx,'index',[]);this.addField_MFFloat(ctx,'textureZ',[]);this._mesh._primType='POINTS';},{drawOrder:function(){return this._vf.drawOrder.toLowerCase();},nodeChanged:function()
{var coordNode=this._cf.coord.node;x3dom.debug.assert(coordNode,"ParticleSet without coord node!");var positions=coordNode.getPoints();var numColComponents=3;var colorNode=this._cf.color.node;var colors=new x3dom.fields.MFColor();if(colorNode){colors=colorNode._vf.color;x3dom.debug.assert(positions.length==colors.length,"Size of color and coord array differs!");if(x3dom.isa(colorNode,x3dom.nodeTypes.ColorRGBA)){numColComponents=4;}}
var normalNode=this._cf.normal.node;var normals=new x3dom.fields.MFVec3f();if(normalNode){normals=normalNode._vf.vector;}
var indices=[];if(this.drawOrder()!="any"){indices=this._vf.index.toGL();if(indices.length==0){var i,n=positions.length;indices=new Array(n);for(i=0;i<n;i++){indices[i]=i;}}}
this._mesh._numColComponents=numColComponents;this._mesh._lit=false;this._mesh._indices[0]=indices;this._mesh._positions[0]=positions.toGL();this._mesh._colors[0]=colors.toGL();this._mesh._normals[0]=normals.toGL();this._mesh._texCoords[0]=[];this.invalidateVolume();this._mesh._numCoords=this._mesh._positions[0].length/3;},fieldChanged:function(fieldName)
{var pnts=null;if(fieldName=="index")
{this._mesh._indices[0]=this._vf.index.toGL();Array.forEach(this._parentNodes,function(node){node._dirty.indexes=true;});}
else if(fieldName=="size")
{Array.forEach(this._parentNodes,function(node){node._dirty.specialAttribs=true;});}
else if(fieldName=="coord")
{pnts=this._cf.coord.node.getPoints();this._mesh._positions[0]=pnts.toGL();var indices=[];if(this.drawOrder()!="any"){indices=this._vf.index.toGL();if(indices.length==0){var i,n=pnts.length;indices=new Array(n);for(i=0;i<n;i++){indices[i]=i;}}}
this._mesh._indices[0]=indices;this.invalidateVolume();Array.forEach(this._parentNodes,function(node){node._dirty.positions=true;node._dirty.indexes=true;node.invalidateVolume();});}
else if(fieldName=="color")
{pnts=this._cf.color.node._vf.color;this._mesh._colors[0]=pnts.toGL();Array.forEach(this._parentNodes,function(node){node._dirty.colors=true;});}}}));x3dom.registerNodeType("ClipPlane","Rendering",defineClass(x3dom.nodeTypes.X3DChildNode,function(ctx){x3dom.nodeTypes.ClipPlane.superClass.call(this,ctx);this.addField_SFBool(ctx,'enabled',true);this.addField_SFVec4f(ctx,'plane',0,1,0,0);this.addField_SFFloat(ctx,'cappingStrength',0.0);this.addField_SFColor(ctx,'cappingColor',1.0,1.0,1.0);this.addField_SFBool(ctx,'on',true);},{fieldChanged:function(fieldName){if(fieldName=="enabled"||fieldName=="on"){}},nodeChanged:function(){x3dom.nodeTypes.ClipPlane.count++;},onRemove:function(){x3dom.nodeTypes.ClipPlane.count--;},parentAdded:function(parent){},parentRemoved:function(parent){}}));x3dom.nodeTypes.ClipPlane.count=0;x3dom.registerNodeType("X3DAppearanceNode","Shape",defineClass(x3dom.nodeTypes.X3DNode,function(ctx){x3dom.nodeTypes.X3DAppearanceNode.superClass.call(this,ctx);}));x3dom.registerNodeType("Appearance","Shape",defineClass(x3dom.nodeTypes.X3DAppearanceNode,function(ctx){x3dom.nodeTypes.Appearance.superClass.call(this,ctx);this.addField_SFNode('material',x3dom.nodeTypes.X3DMaterialNode);this.addField_SFNode('texture',x3dom.nodeTypes.X3DTextureNode);this.addField_SFNode('textureTransform',x3dom.nodeTypes.X3DTextureTransformNode);this.addField_SFNode('lineProperties',x3dom.nodeTypes.LineProperties);this.addField_SFNode('colorMaskMode',x3dom.nodeTypes.ColorMaskMode);this.addField_SFNode('blendMode',x3dom.nodeTypes.BlendMode);this.addField_SFNode('depthMode',x3dom.nodeTypes.DepthMode);this.addField_MFNode('shaders',x3dom.nodeTypes.X3DShaderNode);this.addField_SFString(ctx,'sortType','auto');this.addField_SFInt32(ctx,'sortKey',0);this.addField_SFFloat(ctx,'alphaClipThreshold',0.1);this._shader=null;},{fieldChanged:function(fieldName){if(fieldName=="alphaClipThreshold"){Array.forEach(this._parentNodes,function(shape){shape.setAppDirty();});}},nodeChanged:function(){if(!this._cf.material.node){}
if(this._cf.shaders.nodes.length){this._shader=this._cf.shaders.nodes[0];}
else if(this._shader)
this._shader=null;Array.forEach(this._parentNodes,function(shape){shape.setAppDirty();});this.checkSortType();},checkSortType:function(){if(this._vf.sortType=='auto'){if(this._cf.material.node&&(this._cf.material.node._vf.transparency>0||this._cf.material.node._vf.backTransparency&&this._cf.material.node._vf.backTransparency>0)){this._vf.sortType='transparent';}
else if(this._cf.texture.node&&this._cf.texture.node._vf.url.length){if(this._cf.texture.node._vf.url[0].toLowerCase().indexOf('.'+'png')>=0){this._vf.sortType='transparent';}
else{this._vf.sortType='opaque';}}
else{this._vf.sortType='opaque';}}},texTransformMatrix:function(){if(this._cf.textureTransform.node===null){return x3dom.fields.SFMatrix4f.identity();}
else{return this._cf.textureTransform.node.texTransformMatrix();}},parentAdded:function(parent){if(this!=x3dom.nodeTypes.Appearance._defaultNode){parent.setAppDirty();}}}));x3dom.nodeTypes.Appearance.defaultNode=function(){if(!x3dom.nodeTypes.Appearance._defaultNode){x3dom.nodeTypes.Appearance._defaultNode=new x3dom.nodeTypes.Appearance();x3dom.nodeTypes.Appearance._defaultNode.nodeChanged();}
return x3dom.nodeTypes.Appearance._defaultNode;};x3dom.registerNodeType("X3DAppearanceChildNode","Shape",defineClass(x3dom.nodeTypes.X3DNode,function(ctx){x3dom.nodeTypes.X3DAppearanceChildNode.superClass.call(this,ctx);}));x3dom.registerNodeType("BlendMode","Shape",defineClass(x3dom.nodeTypes.X3DAppearanceChildNode,function(ctx){x3dom.nodeTypes.BlendMode.superClass.call(this,ctx);this.addField_SFString(ctx,'srcFactor',"src_alpha");this.addField_SFString(ctx,'destFactor',"one_minus_src_alpha");this.addField_SFColor(ctx,'color',1,1,1);this.addField_SFFloat(ctx,'colorTransparency',0);this.addField_SFString(ctx,'alphaFunc',"none");this.addField_SFFloat(ctx,'alphaFuncValue',0);this.addField_SFString(ctx,'equation',"none");}));x3dom.registerNodeType("DepthMode","Shape",defineClass(x3dom.nodeTypes.X3DAppearanceChildNode,function(ctx){x3dom.nodeTypes.DepthMode.superClass.call(this,ctx);this.addField_SFBool(ctx,'enableDepthTest',true);this.addField_SFString(ctx,'depthFunc',"none");this.addField_SFBool(ctx,'readOnly',false);this.addField_SFFloat(ctx,'zNearRange',-1);this.addField_SFFloat(ctx,'zFarRange',-1);}));x3dom.registerNodeType("ColorMaskMode","Shape",defineClass(x3dom.nodeTypes.X3DAppearanceChildNode,function(ctx){x3dom.nodeTypes.ColorMaskMode.superClass.call(this,ctx);this.addField_SFBool(ctx,'maskR',true);this.addField_SFBool(ctx,'maskG',true);this.addField_SFBool(ctx,'maskB',true);this.addField_SFBool(ctx,'maskA',true);}));x3dom.registerNodeType("LineProperties","Shape",defineClass(x3dom.nodeTypes.X3DAppearanceChildNode,function(ctx){x3dom.nodeTypes.LineProperties.superClass.call(this,ctx);this.addField_SFBool(ctx,'applied',true);this.addField_SFInt32(ctx,'linetype',1);this.addField_SFFloat(ctx,'linewidthScaleFactor',0);}));x3dom.registerNodeType("X3DMaterialNode","Shape",defineClass(x3dom.nodeTypes.X3DAppearanceChildNode,function(ctx){x3dom.nodeTypes.X3DMaterialNode.superClass.call(this,ctx);}));x3dom.registerNodeType("Material","Shape",defineClass(x3dom.nodeTypes.X3DMaterialNode,function(ctx){x3dom.nodeTypes.Material.superClass.call(this,ctx);this.addField_SFFloat(ctx,'ambientIntensity',0.2);this.addField_SFColor(ctx,'diffuseColor',0.8,0.8,0.8);this.addField_SFColor(ctx,'emissiveColor',0,0,0);this.addField_SFFloat(ctx,'shininess',0.2);this.addField_SFColor(ctx,'specularColor',0,0,0);this.addField_SFFloat(ctx,'transparency',0);},{fieldChanged:function(fieldName){if(fieldName=="ambientIntensity"||fieldName=="diffuseColor"||fieldName=="emissiveColor"||fieldName=="shininess"||fieldName=="specularColor"||fieldName=="transparency")
{Array.forEach(this._parentNodes,function(app){Array.forEach(app._parentNodes,function(shape){shape._dirty.material=true;});app.checkSortType();});}}}));x3dom.nodeTypes.Material.defaultNode=function(){if(!x3dom.nodeTypes.Material._defaultNode){x3dom.nodeTypes.Material._defaultNode=new x3dom.nodeTypes.Material();x3dom.nodeTypes.Material._defaultNode.nodeChanged();}
return x3dom.nodeTypes.Material._defaultNode;};x3dom.registerNodeType("TwoSidedMaterial","Shape",defineClass(x3dom.nodeTypes.Material,function(ctx){x3dom.nodeTypes.TwoSidedMaterial.superClass.call(this,ctx);this.addField_SFFloat(ctx,'backAmbientIntensity',0.2);this.addField_SFColor(ctx,'backDiffuseColor',0.8,0.8,0.8);this.addField_SFColor(ctx,'backEmissiveColor',0,0,0);this.addField_SFFloat(ctx,'backShininess',0.2);this.addField_SFColor(ctx,'backSpecularColor',0,0,0);this.addField_SFFloat(ctx,'backTransparency',0);this.addField_SFBool(ctx,'separateBackColor',false);},{fieldChanged:function(fieldName){if(fieldName=="ambientIntensity"||fieldName=="diffuseColor"||fieldName=="emissiveColor"||fieldName=="shininess"||fieldName=="specularColor"||fieldName=="transparency"||fieldName=="backAmbientIntensity"||fieldName=="backDiffuseColor"||fieldName=="backEmissiveColor"||fieldName=="backShininess"||fieldName=="backSpecularColor"||fieldName=="backTransparency"||fieldName=="separateBackColor")
{Array.forEach(this._parentNodes,function(app){Array.forEach(app._parentNodes,function(shape){shape._dirty.material=true;});app.checkSortType();});}}}));x3dom.registerNodeType("X3DShapeNode","Shape",defineClass(x3dom.nodeTypes.X3DBoundedObject,function(ctx){x3dom.nodeTypes.X3DShapeNode.superClass.call(this,ctx);this.addField_SFBool(ctx,'isPickable',true);this.addField_SFInt32(ctx,'idOffset',0);this.addField_SFNode('appearance',x3dom.nodeTypes.X3DAppearanceNode);this.addField_SFNode('geometry',x3dom.nodeTypes.X3DGeometryNode);this._objectID=0;this._shaderProperties=null;this._clipPlanes=[];this._cleanupGLObjects=null;this._dirty={positions:true,normals:true,texcoords:true,colors:true,specialAttribs:true,indexes:true,texture:true,material:true,text:true,shader:true,ids:true};this._coordStrideOffset=[0,0];this._normalStrideOffset=[0,0];this._texCoordStrideOffset=[0,0];this._colorStrideOffset=[0,0];this._idStrideOffset=[0,0];this._tessellationProperties=[];},{collectDrawableObjects:function(transform,drawableCollection,singlePath,invalidateCache,planeMask,clipPlanes)
{var graphState=this.graphState();if(singlePath&&(this._parentNodes.length>1))
singlePath=false;if(singlePath&&(invalidateCache=invalidateCache||this.cacheInvalid()))
this.invalidateCache();if(!this._cf.geometry.node||drawableCollection.cull(transform,graphState,singlePath,planeMask)<=0){return false;}
if(singlePath&&!this._graph.globalMatrix)
this._graph.globalMatrix=transform;if(this._clipPlanes.length!=clipPlanes.length)
{this._dirty.shader=true;}
this._clipPlanes=clipPlanes;drawableCollection.addShape(this,transform,graphState);return true;},getVolume:function()
{var vol=this._graph.volume;if(!this.volumeValid()&&this._vf.render)
{var geo=this._cf.geometry.node;var childVol=geo?geo.getVolume():null;if(childVol&&childVol.isValid())
vol.extendBounds(childVol.min,childVol.max);}
return vol;},getCenter:function(){var geo=this._cf.geometry.node;return(geo?geo.getCenter():new x3dom.fields.SFVec3f(0,0,0));},getDiameter:function(){var geo=this._cf.geometry.node;return(geo?geo.getDiameter():0);},doIntersect:function(line){return this._cf.geometry.node.doIntersect(line);},forceUpdateCoverage:function()
{var geo=this._cf.geometry.node;return(geo?geo.forceUpdateCoverage():false);},tessellationProperties:function()
{var geo=this._cf.geometry.node;if(geo&&geo._indexOffset)
return geo._indexOffset;else
return this._tessellationProperties;},isLit:function(){return this._cf.geometry.node._vf.lit;},isSolid:function(){var twoSidedMat=(this._cf.appearance.node&&this._cf.appearance.node._cf.material.node&&x3dom.isa(this._cf.appearance.node._cf.material.node,x3dom.nodeTypes.TwoSidedMaterial));return this._cf.geometry.node._vf.solid&&!twoSidedMat;},isCCW:function(){return this._cf.geometry.node._vf.ccw;},parentRemoved:function(parent){for(var i=0,n=this._childNodes.length;i<n;i++){var child=this._childNodes[i];if(child){child.parentRemoved(this);}}
if(parent)
parent.invalidateVolume();if(this._parentNodes.length>0)
this.invalidateVolume();if(this._cleanupGLObjects){this._cleanupGLObjects();}},unsetDirty:function(){this._dirty.positions=false;this._dirty.normals=false;this._dirty.texcoords=false;this._dirty.colors=false;this._dirty.specialAttribs=false;this._dirty.indexes=false;this._dirty.texture=false;this._dirty.material=false;this._dirty.text=false;this._dirty.shader=false;},unsetGeoDirty:function(){this._dirty.positions=false;this._dirty.normals=false;this._dirty.texcoords=false;this._dirty.colors=false;this._dirty.specialAttribs=false;this._dirty.indexes=false;},setAllDirty:function(){this._dirty.positions=true;this._dirty.normals=true;this._dirty.texcoords=true;this._dirty.colors=true;this._dirty.specialAttribs=true;this._dirty.indexes=true;this._dirty.texture=true;this._dirty.material=true;this._dirty.text=true;this._dirty.shader=true;this.invalidateVolume();},setAppDirty:function(){this._dirty.texture=true;this._dirty.material=true;this._dirty.shader=true;},setGeoDirty:function(){this._dirty.positions=true;this._dirty.normals=true;this._dirty.texcoords=true;this._dirty.colors=true;this._dirty.specialAttribs=true;this._dirty.indexes=true;this.invalidateVolume();},getShaderProperties:function(viewarea)
{if(this._shaderProperties==null||this._dirty.shader==true||(this._webgl!==undefined&&this._webgl.dirtyLighting!=x3dom.Utils.checkDirtyLighting(viewarea))||x3dom.Utils.checkDirtyEnvironment(viewarea,this._shaderProperties)==true)
{this._shaderProperties=x3dom.Utils.generateProperties(viewarea,this);this._dirty.shader=false;if(this._webgl!==undefined)
{this._webgl.dirtyLighting=x3dom.Utils.checkDirtyLighting(viewarea);}}
return this._shaderProperties;},getTextures:function(){var textures=[];var appearance=this._cf.appearance.node;if(appearance){var tex=appearance._cf.texture.node;if(tex){if(x3dom.isa(tex,x3dom.nodeTypes.MultiTexture)){textures=textures.concat(tex.getTextures());}
else{textures.push(tex);}}
var shader=appearance._cf.shaders.nodes[0];if(shader){if(x3dom.isa(shader,x3dom.nodeTypes.CommonSurfaceShader)){textures=textures.concat(shader.getTextures());}}}
var geometry=this._cf.geometry.node;if(geometry){if(x3dom.isa(geometry,x3dom.nodeTypes.ImageGeometry)){textures=textures.concat(geometry.getTextures());}
else if(x3dom.isa(geometry,x3dom.nodeTypes.Text)){textures=textures.concat(geometry);}}
return textures;}}));x3dom.registerNodeType("Shape","Shape",defineClass(x3dom.nodeTypes.X3DShapeNode,function(ctx){x3dom.nodeTypes.Shape.superClass.call(this,ctx);},{nodeChanged:function(){if(!this._cf.appearance.node){}
if(!this._cf.geometry.node){if(this._DEF)
x3dom.debug.logError("No geometry given in Shape/"+this._DEF);}
else if(!this._objectID){this._objectID=++x3dom.nodeTypes.Shape.objectID;x3dom.nodeTypes.Shape.idMap.nodeID[this._objectID]=this;}
this.invalidateVolume();}}));x3dom.nodeTypes.Shape.shaderPartID=0;x3dom.nodeTypes.Shape.objectID=0;x3dom.nodeTypes.Shape.idMap={nodeID:{},remove:function(obj){for(var prop in this.nodeID){if(this.nodeID.hasOwnProperty(prop)){var val=this.nodeID[prop];if(val._objectID&&obj._objectID&&val._objectID===obj._objectID)
{delete this.nodeID[prop];x3dom.debug.logInfo("Unreg "+val._objectID);}}}}};x3dom.registerNodeType("ExternalShape","Shape",defineClass(x3dom.nodeTypes.X3DShapeNode,function(ctx){x3dom.nodeTypes.ExternalShape.superClass.call(this,ctx);this.addField_MFString(ctx,'url',[]);this._currentURLIdx=0;this._cf.geometry.node=new x3dom.nodeTypes.X3DSpatialGeometryNode(ctx);this.loaded=false;},{update:function(shape,shaderProgram,gl,viewarea,context){var that=this;if(this._vf['url'].length==0||this._currentURLIdx>=this._vf['url'].length)
{return;}
var xhr=new XMLHttpRequest();xhr.open("GET",this._nameSpace.getURL(this._vf['url'][this._currentURLIdx]),true);xhr.responseType="arraybuffer";xhr.send(null);xhr.onerror=function(){x3dom.debug.logError("Unable to load SRC data from URL \""+that._vf['url'][that._currentURLIdx]+"\"");};xhr.onload=function(){if((xhr.status==200||xhr.status==0)){var glTF=new x3dom.glTF.glTFLoader(xhr.response);if(glTF.header.sceneLength>0)
{glTF.loaded={};glTF.loaded.meshes={};glTF.loaded.meshCount=0;that.glTF=glTF;var url=that._vf['url'][that._currentURLIdx];if(url.includes('#'))
{var split=url.split('#');var meshName=split[split.length-1];glTF.getMesh(shape,shaderProgram,gl,meshName);}
else
{glTF.getScene(shape,shaderProgram,gl);}
for(var key in glTF._mesh){if(!glTF._mesh.hasOwnProperty(key))continue;that._cf.geometry.node._mesh[key]=glTF._mesh[key];}}
else
{if((that._currentURLIdx+1)<that._vf['url'].length)
{x3dom.debug.logWarning("Invalid SRC data, loaded from URL \""+
that._vf['url'][that._currentURLIdx]+"\", trying next specified URL");++that._currentURLIdx;that.update(shape,shaderProgram,gl,viewarea,context);}
else
{x3dom.debug.logError("Invalid SRC data, loaded from URL \""+
that._vf['url'][that._currentURLIdx]+"\","+" no other URLs left to try.");}}}
else
{if((that._currentURLIdx+1)<that._vf['url'].length)
{x3dom.debug.logWarning("Invalid SRC data, loaded from URL \""+
that._vf['url'][that._currentURLIdx]+"\", trying next specified URL");++that._currentURLIdx;that.update(shape,shaderProgram,gl,viewarea,context);}
else
{x3dom.debug.logError("Invalid SRC data, loaded from URL \""+
that._vf['url'][that._currentURLIdx]+"\","+" no other URLs left to try.");}}};},collectDrawableObjects:function(transform,drawableCollection,singlePath,invalidateCache,planeMask,clipPlanes)
{var graphState=this.graphState();if(singlePath&&(this._parentNodes.length>1))
singlePath=false;if(singlePath&&(invalidateCache=invalidateCache||this.cacheInvalid()))
this.invalidateCache();if(singlePath&&!this._graph.globalMatrix)
this._graph.globalMatrix=transform;if(this._clipPlanes.length!=clipPlanes.length)
{this._dirty.shader=true;}
this._clipPlanes=clipPlanes;drawableCollection.addShape(this,transform,graphState);return true;},getShaderProperties:function(viewarea)
{var properties=x3dom.Utils.generateProperties(viewarea,this);properties.CSHADER=-1;properties.LIGHTS=viewarea.getLights().length+(viewarea._scene.getNavigationInfo()._vf.headlight);properties.EMPTY_SHADER=1;return properties;},nodeChanged:function()
{if(!this._objectID){this._objectID=++x3dom.nodeTypes.Shape.objectID;x3dom.nodeTypes.Shape.idMap.nodeID[this._objectID]=this;}}}));x3dom.registerNodeType("X3DLightNode","Lighting",defineClass(x3dom.nodeTypes.X3DChildNode,function(ctx){x3dom.nodeTypes.X3DLightNode.superClass.call(this,ctx);if(ctx)
ctx.doc._nodeBag.lights.push(this);else
x3dom.debug.logWarning("X3DLightNode: No runtime context found!");this._lightID=0;this._dirty=true;this.addField_SFFloat(ctx,'ambientIntensity',0);this.addField_SFColor(ctx,'color',1,1,1);this.addField_SFFloat(ctx,'intensity',1);this.addField_SFBool(ctx,'global',false);this.addField_SFBool(ctx,'on',true);this.addField_SFFloat(ctx,'shadowIntensity',0);this.addField_SFInt32(ctx,'shadowMapSize',1024);this.addField_SFInt32(ctx,'shadowFilterSize',0);this.addField_SFFloat(ctx,'shadowOffset',0);this.addField_SFFloat(ctx,'zNear',-1);this.addField_SFFloat(ctx,'zFar',-1);},{getViewMatrix:function(vec){return x3dom.fields.SFMatrix4f.identity;},nodeChanged:function(){if(!this._lightID){this._lightID=++x3dom.nodeTypes.X3DLightNode.lightID;}},fieldChanged:function(fieldName)
{if(this._vf.hasOwnProperty(fieldName)){this._dirty=true;}},parentRemoved:function(parent)
{if(this._parentNodes.length===1&&this._parentNodes[0]==parent){var doc=this.findX3DDoc();for(var i=0,n=doc._nodeBag.lights.length;i<n;i++){if(doc._nodeBag.lights[i]===this){doc._nodeBag.lights.splice(i,1);}}}},onRemove:function()
{}}));x3dom.nodeTypes.X3DLightNode.lightID=0;x3dom.registerNodeType("DirectionalLight","Lighting",defineClass(x3dom.nodeTypes.X3DLightNode,function(ctx){x3dom.nodeTypes.DirectionalLight.superClass.call(this,ctx);this.addField_SFVec3f(ctx,'direction',0,0,-1);this.addField_SFInt32(ctx,'shadowCascades',1);this.addField_SFFloat(ctx,'shadowSplitFactor',1);this.addField_SFFloat(ctx,'shadowSplitOffset',0.1);},{getViewMatrix:function(vec){var dir=this.getCurrentTransform().multMatrixVec(this._vf.direction).normalize();var orientation=x3dom.fields.Quaternion.rotateFromTo(new x3dom.fields.SFVec3f(0,0,-1),dir);return orientation.toMatrix().transpose().mult(x3dom.fields.SFMatrix4f.translation(vec.negate()));}}));x3dom.registerNodeType("PointLight","Lighting",defineClass(x3dom.nodeTypes.X3DLightNode,function(ctx){x3dom.nodeTypes.PointLight.superClass.call(this,ctx);this.addField_SFVec3f(ctx,'attenuation',1,0,0);this.addField_SFVec3f(ctx,'location',0,0,0);this.addField_SFFloat(ctx,'radius',100);this._vf.global=true;},{getViewMatrix:function(vec){var pos=this.getCurrentTransform().multMatrixPnt(this._vf.location);var orientation=x3dom.fields.Quaternion.rotateFromTo(new x3dom.fields.SFVec3f(0,0,-1),vec);return orientation.toMatrix().transpose().mult(x3dom.fields.SFMatrix4f.translation(pos.negate()));}}));x3dom.registerNodeType("SpotLight","Lighting",defineClass(x3dom.nodeTypes.X3DLightNode,function(ctx){x3dom.nodeTypes.SpotLight.superClass.call(this,ctx);this.addField_SFVec3f(ctx,'direction',0,0,-1);this.addField_SFVec3f(ctx,'attenuation',1,0,0);this.addField_SFVec3f(ctx,'location',0,0,0);this.addField_SFFloat(ctx,'radius',100);this.addField_SFFloat(ctx,'beamWidth',1.5707963);this.addField_SFFloat(ctx,'cutOffAngle',1.5707963);this.addField_SFInt32(ctx,'shadowCascades',1);this.addField_SFFloat(ctx,'shadowSplitFactor',1);this.addField_SFFloat(ctx,'shadowSplitOffset',0.1);this._vf.global=true;},{getViewMatrix:function(vec){var pos=this.getCurrentTransform().multMatrixPnt(this._vf.location);var dir=this.getCurrentTransform().multMatrixVec(this._vf.direction).normalize();var orientation=x3dom.fields.Quaternion.rotateFromTo(new x3dom.fields.SFVec3f(0,0,-1),dir);return orientation.toMatrix().transpose().mult(x3dom.fields.SFMatrix4f.translation(pos.negate()));}}));x3dom.registerNodeType("X3DFollowerNode","Followers",defineClass(x3dom.nodeTypes.X3DChildNode,function(ctx){x3dom.nodeTypes.X3DFollowerNode.superClass.call(this,ctx);if(ctx)
ctx.doc._nodeBag.followers.push(this);else
x3dom.debug.logWarning("X3DFollowerNode: No runtime context found!");this.addField_SFBool(ctx,'isActive',false);this._eps=x3dom.fields.Eps;},{parentRemoved:function(parent)
{if(this._parentNodes.length===0){var doc=this.findX3DDoc();for(var i=0,n=doc._nodeBag.followers.length;i<n;i++){if(doc._nodeBag.followers[i]===this){doc._nodeBag.followers.splice(i,1);}}}},tick:function(t){return false;},stepResponse:function(t)
{if(t<=0){return 0;}
if(t>=this._vf.duration){return 1;}
return this.stepResponseCore(t/this._vf.duration);},stepResponseCore:function(T)
{return 0.5-0.5*Math.cos(T*Math.PI);}}));x3dom.registerNodeType("X3DChaserNode","Followers",defineClass(x3dom.nodeTypes.X3DFollowerNode,function(ctx){x3dom.nodeTypes.X3DChaserNode.superClass.call(this,ctx);this.addField_SFTime(ctx,'duration',1);this._initDone=false;this._stepTime=0;this._currTime=0;this._bufferEndTime=0;this._numSupports=60;}));x3dom.registerNodeType("X3DDamperNode","Followers",defineClass(x3dom.nodeTypes.X3DFollowerNode,function(ctx){x3dom.nodeTypes.X3DDamperNode.superClass.call(this,ctx);this.addField_SFTime(ctx,'tau',0.3);this.addField_SFFloat(ctx,'tolerance',-1);this.addField_SFInt32(ctx,'order',3);this._eps=this._vf.tolerance<0?this._eps:this._vf.tolerance;this._lastTick=0;}));x3dom.registerNodeType("ColorChaser","Followers",defineClass(x3dom.nodeTypes.X3DChaserNode,function(ctx){x3dom.nodeTypes.ColorChaser.superClass.call(this,ctx);this.addField_SFColor(ctx,'initialDestination',0.8,0.8,0.8);this.addField_SFColor(ctx,'initialValue',0.8,0.8,0.8);this.addField_SFColor(ctx,'value',0,0,0);this.addField_SFColor(ctx,'destination',0,0,0);this._buffer=new x3dom.fields.MFColor();this._previousValue=new x3dom.fields.SFColor(0,0,0);this._value=new x3dom.fields.SFColor(0,0,0);this.initialize();},{fieldChanged:function(fieldName)
{if(fieldName.indexOf("destination")>=0)
{this.initialize();this.updateBuffer(this._currTime);if(!this._vf.isActive){this.postMessage('isActive',true);}}
else if(fieldName.indexOf("value")>=0)
{this.initialize();this._previousValue.setValues(this._vf.value);for(var C=1;C<this._buffer.length;C++){this._buffer[C].setValues(this._vf.value);}
this.postMessage('value',this._vf.value);if(!this._vf.isActive){this.postMessage('isActive',true);}}},initialize:function()
{if(!this._initDone)
{this._initDone=true;this._vf.destination=this._vf.initialDestination;this._buffer.length=this._numSupports;this._buffer[0]=this._vf.initialDestination;for(var C=1;C<this._buffer.length;C++){this._buffer[C]=this._vf.initialValue;}
this._previousValue=this._vf.initialValue;this._stepTime=this._vf.duration/this._numSupports;var active=!this._buffer[0].equals(this._buffer[1],this._eps);if(this._vf.isActive!==active){this.postMessage('isActive',active);}}},tick:function(now)
{this.initialize();this._currTime=now;if(!this._bufferEndTime)
{this._bufferEndTime=now;this._value=this._vf.initialValue;this.postMessage('value',this._value);return true;}
var Frac=this.updateBuffer(now);var Output=this._previousValue;var DeltaIn=this._buffer[this._buffer.length-1].subtract(this._previousValue);var DeltaOut=DeltaIn.multiply(this.stepResponse((this._buffer.length-1+Frac)*this._stepTime));Output=Output.add(DeltaOut);for(var C=this._buffer.length-2;C>=0;C--)
{DeltaIn=this._buffer[C].subtract(this._buffer[C+1]);DeltaOut=DeltaIn.multiply(this.stepResponse((C+Frac)*this._stepTime));Output=Output.add(DeltaOut);}
if(!Output.equals(this._value,this._eps)){this._value.setValues(Output);this.postMessage('value',this._value);}
else{this.postMessage('isActive',false);}
return this._vf.isActive;},updateBuffer:function(now)
{var Frac=(now-this._bufferEndTime)/this._stepTime;var C;var NumToShift;var Alpha;if(Frac>=1)
{NumToShift=Math.floor(Frac);Frac-=NumToShift;if(NumToShift<this._buffer.length)
{this._previousValue=this._buffer[this._buffer.length-NumToShift];for(C=this._buffer.length-1;C>=NumToShift;C--){this._buffer[C]=this._buffer[C-NumToShift];}
for(C=0;C<NumToShift;C++)
{Alpha=C/NumToShift;this._buffer[C]=this._buffer[NumToShift].multiply(Alpha).add(this._vf.destination.multiply((1-Alpha)));}}
else
{this._previousValue=(NumToShift==this._buffer.length)?this._buffer[0]:this._vf.destination;for(C=0;C<this._buffer.length;C++){this._buffer[C]=this._vf.destination;}}
this._bufferEndTime+=NumToShift*this._stepTime;}
return Frac;}}));x3dom.registerNodeType("ColorDamper","Followers",defineClass(x3dom.nodeTypes.X3DDamperNode,function(ctx){x3dom.nodeTypes.ColorDamper.superClass.call(this,ctx);this.addField_SFColor(ctx,'initialDestination',0.8,0.8,0.8);this.addField_SFColor(ctx,'initialValue',0.8,0.8,0.8);this.addField_SFColor(ctx,'value',0,0,0);this.addField_SFColor(ctx,'destination',0,0,0);this._value0=new x3dom.fields.SFColor(0,0,0);this._value1=new x3dom.fields.SFColor(0,0,0);this._value2=new x3dom.fields.SFColor(0,0,0);this._value3=new x3dom.fields.SFColor(0,0,0);this._value4=new x3dom.fields.SFColor(0,0,0);this._value5=new x3dom.fields.SFColor(0,0,0);this.initialize();},{fieldChanged:function(fieldName)
{if(fieldName==="tolerance")
{this._eps=this._vf.tolerance<0?0.001:this._vf.tolerance;}
else if(fieldName.indexOf("destination")>=0)
{if(!this._value0.equals(this._vf.destination,this._eps)){this._value0=this._vf.destination;if(!this._vf.isActive){this.postMessage('isActive',true);}}}
else if(fieldName.indexOf("value")>=0)
{this._value1.setValues(this._vf.value);this._value2.setValues(this._vf.value);this._value3.setValues(this._vf.value);this._value4.setValues(this._vf.value);this._value5.setValues(this._vf.value);this._lastTick=0;this.postMessage('value',this._value5);if(!this._vf.isActive){this._lastTick=0;this.postMessage('isActive',true);}}},initialize:function()
{this._value0.setValues(this._vf.initialDestination);this._value1.setValues(this._vf.initialValue);this._value2.setValues(this._vf.initialValue);this._value3.setValues(this._vf.initialValue);this._value4.setValues(this._vf.initialValue);this._value5.setValues(this._vf.initialValue);this._lastTick=0;var active=!this._value0.equals(this._value1,this._eps);if(this._vf.isActive!==active){this.postMessage('isActive',active);}},distance:function(a,b)
{var diff=a.subtract(b);return Math.sqrt(diff.r*diff.r+diff.g*diff.g+diff.b*diff.b);},tick:function(now)
{if(!this._lastTick)
{this._lastTick=now;return false;}
var delta=now-this._lastTick;var alpha=Math.exp(-delta/this._vf.tau);this._value1=this._vf.order>0&&this._vf.tau?this._value0.add(this._value1.subtract(this._value0).multiply(alpha)):new x3dom.fields.SFColor(this._value0.r,this._value0.g,this._value0.b);this._value2=this._vf.order>1&&this._vf.tau?this._value1.add(this._value2.subtract(this._value1).multiply(alpha)):new x3dom.fields.SFColor(this._value1.r,this._value1.g,this._value1.b);this._value3=this._vf.order>2&&this._vf.tau?this._value2.add(this._value3.subtract(this._value2).multiply(alpha)):new x3dom.fields.SFColor(this._value2.r,this._value2.g,this._value2.b);this._value4=this._vf.order>3&&this._vf.tau?this._value3.add(this._value4.subtract(this._value3).multiply(alpha)):new x3dom.fields.SFColor(this._value3.r,this._value3.g,this._value3.b);this._value5=this._vf.order>4&&this._vf.tau?this._value4.add(this._value5.subtract(this._value4).multiply(alpha)):new x3dom.fields.SFColor(this._value4.r,this._value4.g,this._value4.b);var dist=this.distance(this._value1,this._value0);if(this._vf.order>1)
{var dist2=this.distance(this._value2,this._value1);if(dist2>dist){dist=dist2;}}
if(this._vf.order>2)
{var dist3=this.distance(this._value3,this._value2);if(dist3>dist){dist=dist3;}}
if(this._vf.order>3)
{var dist4=this.distance(this._value4,this._value3);if(dist4>dist){dist=dist4;}}
if(this._vf.order>4)
{var dist5=this.distance(this._value5,this._value4);if(dist5>dist){dist=dist5;}}
if(dist<=this._eps)
{this._value1.setValues(this._value0);this._value2.setValues(this._value0);this._value3.setValues(this._value0);this._value4.setValues(this._value0);this._value5.setValues(this._value0);this.postMessage('value',this._value0);this.postMessage('isActive',false);this._lastTick=0;return false;}
this.postMessage('value',this._value5);this._lastTick=now;return true;}}));x3dom.registerNodeType("OrientationChaser","Followers",defineClass(x3dom.nodeTypes.X3DChaserNode,function(ctx){x3dom.nodeTypes.OrientationChaser.superClass.call(this,ctx);this.addField_SFRotation(ctx,'initialDestination',0,1,0,0);this.addField_SFRotation(ctx,'initialValue',0,1,0,0);this.addField_SFRotation(ctx,'value',0,1,0,0);this.addField_SFRotation(ctx,'destination',0,1,0,0);this._numSupports=30;this._buffer=new x3dom.fields.MFRotation();this._previousValue=new x3dom.fields.Quaternion(0,1,0,0);this._value=new x3dom.fields.Quaternion(0,1,0,0);this.initialize();},{fieldChanged:function(fieldName)
{if(fieldName.indexOf("destination")>=0)
{this.initialize();this.updateBuffer(this._currTime);if(!this._vf.isActive){this.postMessage('isActive',true);}}
else if(fieldName.indexOf("value")>=0)
{this.initialize();this._previousValue.setValues(this._vf.value);for(var C=1;C<this._buffer.length;C++){this._buffer[C].setValues(this._vf.value);}
this.postMessage('value',this._vf.value);if(!this._vf.isActive){this.postMessage('isActive',true);}}},initialize:function()
{if(!this._initDone)
{this._initDone=true;this._vf.destination=x3dom.fields.Quaternion.copy(this._vf.initialDestination);this._buffer.length=this._numSupports;this._buffer[0]=x3dom.fields.Quaternion.copy(this._vf.initialDestination);for(var C=1;C<this._buffer.length;C++){this._buffer[C]=x3dom.fields.Quaternion.copy(this._vf.initialValue);}
this._previousValue=x3dom.fields.Quaternion.copy(this._vf.initialValue);this._stepTime=this._vf.duration/this._numSupports;var active=!this._buffer[0].equals(this._buffer[1],this._eps);if(this._vf.isActive!==active){this.postMessage('isActive',active);}}},tick:function(now)
{this.initialize();this._currTime=now;if(!this._bufferEndTime)
{this._bufferEndTime=now;this._value=x3dom.fields.Quaternion.copy(this._vf.initialValue);this.postMessage('value',this._value);return true;}
var Frac=this.updateBuffer(now);var Output=x3dom.fields.Quaternion.copy(this._previousValue);var DeltaIn=this._previousValue.inverse().multiply(this._buffer[this._buffer.length-1]);Output=Output.slerp(Output.multiply(DeltaIn),this.stepResponse((this._buffer.length-1+Frac)*this._stepTime));for(var C=this._buffer.length-2;C>=0;C--)
{DeltaIn=this._buffer[C+1].inverse().multiply(this._buffer[C]);Output=Output.slerp(Output.multiply(DeltaIn),this.stepResponse((C+Frac)*this._stepTime));}
if(!Output.equals(this._value,this._eps)){Output=Output.normalize(Output);this._value.setValues(Output);this.postMessage('value',this._value);}
else{this.postMessage('isActive',false);}
return this._vf.isActive;},updateBuffer:function(now)
{var Frac=(now-this._bufferEndTime)/this._stepTime;var C;var NumToShift;var Alpha;if(Frac>=1)
{NumToShift=Math.floor(Frac);Frac-=NumToShift;if(NumToShift<this._buffer.length)
{this._previousValue=x3dom.fields.Quaternion.copy(this._buffer[this._buffer.length-NumToShift]);for(C=this._buffer.length-1;C>=NumToShift;C--){this._buffer[C]=x3dom.fields.Quaternion.copy(this._buffer[C-NumToShift]);}
for(C=0;C<NumToShift;C++)
{Alpha=C/NumToShift;this._buffer[C]=this._vf.destination.slerp(this._buffer[NumToShift],Alpha);}}
else
{this._previousValue=x3dom.fields.Quaternion.copy((NumToShift==this._buffer.length)?this._buffer[0]:this._vf.destination);for(C=0;C<this._buffer.length;C++){this._buffer[C]=x3dom.fields.Quaternion.copy(this._vf.destination);}}
this._bufferEndTime+=NumToShift*this._stepTime;}
return Frac;}}));x3dom.registerNodeType("OrientationDamper","Followers",defineClass(x3dom.nodeTypes.X3DDamperNode,function(ctx){x3dom.nodeTypes.OrientationDamper.superClass.call(this,ctx);this.addField_SFRotation(ctx,'initialDestination',0,1,0,0);this.addField_SFRotation(ctx,'initialValue',0,1,0,0);this.addField_SFRotation(ctx,'value',0,1,0,0);this.addField_SFRotation(ctx,'destination',0,1,0,0);this._value0=new x3dom.fields.Quaternion(0,1,0,0);this._value1=new x3dom.fields.Quaternion(0,1,0,0);this._value2=new x3dom.fields.Quaternion(0,1,0,0);this._value3=new x3dom.fields.Quaternion(0,1,0,0);this._value4=new x3dom.fields.Quaternion(0,1,0,0);this._value5=new x3dom.fields.Quaternion(0,1,0,0);this.initialize();},{fieldChanged:function(fieldName)
{if(fieldName==="tolerance")
{this._eps=this._vf.tolerance<0?0.001:this._vf.tolerance;}
else if(fieldName.indexOf("destination")>=0)
{if(!this._value0.equals(this._vf.destination,this._eps)){this._value0=this._vf.destination;if(!this._vf.isActive){this.postMessage('isActive',true);}}}
else if(fieldName.indexOf("value")>=0)
{this._value1.setValues(this._vf.value);this._value2.setValues(this._vf.value);this._value3.setValues(this._vf.value);this._value4.setValues(this._vf.value);this._value5.setValues(this._vf.value);this._lastTick=0;this.postMessage('value',this._value5);if(!this._vf.isActive){this._lastTick=0;this.postMessage('isActive',true);}}},initialize:function()
{this._value0.setValues(this._vf.initialDestination);this._value1.setValues(this._vf.initialValue);this._value2.setValues(this._vf.initialValue);this._value3.setValues(this._vf.initialValue);this._value4.setValues(this._vf.initialValue);this._value5.setValues(this._vf.initialValue);this._lastTick=0;var active=!this._value0.equals(this._value1,this._eps);if(this._vf.isActive!==active){this.postMessage('isActive',active);}},tick:function(now)
{if(!this._lastTick)
{this._lastTick=now;return false;}
var delta=now-this._lastTick;var alpha=Math.exp(-delta/this._vf.tau);this._value1=this._vf.order>0&&this._vf.tau?this._value0.slerp(this._value1,alpha):new x3dom.fields.Quaternion(this._value0.x,this._value0.y,this._value0.z,this._value0.w);this._value2=this._vf.order>1&&this._vf.tau?this._value1.slerp(this._value2,alpha):new x3dom.fields.Quaternion(this._value1.x,this._value1.y,this._value1.z,this._value1.w);this._value3=this._vf.order>2&&this._vf.tau?this._value2.slerp(this._value3,alpha):new x3dom.fields.Quaternion(this._value2.x,this._value2.y,this._value2.z,this._value2.w);this._value4=this._vf.order>3&&this._vf.tau?this._value3.slerp(this._value4,alpha):new x3dom.fields.Quaternion(this._value3.x,this._value3.y,this._value3.z,this._value3.w);this._value5=this._vf.order>4&&this._vf.tau?this._value4.slerp(this._value5,alpha):new x3dom.fields.Quaternion(this._value4.x,this._value4.y,this._value4.z,this._value4.w);var dist=Math.abs(this._value1.inverse().multiply(this._value0).angle());if(this._vf.order>1)
{var dist2=Math.abs(this._value2.inverse().multiply(this._value1).angle());if(dist2>dist){dist=dist2;}}
if(this._vf.order>2)
{var dist3=Math.abs(this._value3.inverse().multiply(this._value2).angle());if(dist3>dist){dist=dist3;}}
if(this._vf.order>3)
{var dist4=Math.abs(this._value4.inverse().multiply(this._value3).angle());if(dist4>dist){dist=dist4;}}
if(this._vf.order>4)
{var dist5=Math.abs(this._value5.inverse().multiply(this._value4).angle());if(dist5>dist){dist=dist5;}}
if(dist<=this._eps)
{this._value1.setValues(this._value0);this._value2.setValues(this._value0);this._value3.setValues(this._value0);this._value4.setValues(this._value0);this._value5.setValues(this._value0);this.postMessage('value',this._value0);this.postMessage('isActive',false);this._lastTick=0;return false;}
this.postMessage('value',this._value5);this._lastTick=now;return true;}}));x3dom.registerNodeType("PositionChaser","Followers",defineClass(x3dom.nodeTypes.X3DChaserNode,function(ctx){x3dom.nodeTypes.PositionChaser.superClass.call(this,ctx);this.addField_SFVec3f(ctx,'initialDestination',0,0,0);this.addField_SFVec3f(ctx,'initialValue',0,0,0);this.addField_SFVec3f(ctx,'value',0,0,0);this.addField_SFVec3f(ctx,'destination',0,0,0);this._buffer=new x3dom.fields.MFVec3f();this._previousValue=new x3dom.fields.SFVec3f(0,0,0);this._value=new x3dom.fields.SFVec3f(0,0,0);this.initialize();},{fieldChanged:function(fieldName)
{if(fieldName.indexOf("destination")>=0)
{this.initialize();this.updateBuffer(this._currTime);if(!this._vf.isActive){this.postMessage('isActive',true);}}
else if(fieldName.indexOf("value")>=0)
{this.initialize();this._previousValue.setValues(this._vf.value);for(var C=1;C<this._buffer.length;C++){this._buffer[C].setValues(this._vf.value);}
this.postMessage('value',this._vf.value);if(!this._vf.isActive){this.postMessage('isActive',true);}}},initialize:function()
{if(!this._initDone)
{this._initDone=true;this._vf.destination=x3dom.fields.SFVec3f.copy(this._vf.initialDestination);this._buffer.length=this._numSupports;this._buffer[0]=x3dom.fields.SFVec3f.copy(this._vf.initialDestination);for(var C=1;C<this._buffer.length;C++){this._buffer[C]=x3dom.fields.SFVec3f.copy(this._vf.initialValue);}
this._previousValue=x3dom.fields.SFVec3f.copy(this._vf.initialValue);this._stepTime=this._vf.duration/this._numSupports;var active=!this._buffer[0].equals(this._buffer[1],this._eps);if(this._vf.isActive!==active){this.postMessage('isActive',active);}}},tick:function(now)
{this.initialize();this._currTime=now;if(!this._bufferEndTime)
{this._bufferEndTime=now;this._value=x3dom.fields.SFVec3f.copy(this._vf.initialValue);this.postMessage('value',this._value);return true;}
var Frac=this.updateBuffer(now);var Output=x3dom.fields.SFVec3f.copy(this._previousValue);var DeltaIn=this._buffer[this._buffer.length-1].subtract(this._previousValue);var DeltaOut=DeltaIn.multiply(this.stepResponse((this._buffer.length-1+Frac)*this._stepTime));Output=Output.add(DeltaOut);for(var C=this._buffer.length-2;C>=0;C--)
{DeltaIn=this._buffer[C].subtract(this._buffer[C+1]);DeltaOut=DeltaIn.multiply(this.stepResponse((C+Frac)*this._stepTime));Output=Output.add(DeltaOut);}
if(!Output.equals(this._value,this._eps)){this._value.setValues(Output);this.postMessage('value',this._value);}
else{this.postMessage('isActive',false);}
return this._vf.isActive;},updateBuffer:function(now)
{var Frac=(now-this._bufferEndTime)/this._stepTime;var C;var NumToShift;var Alpha;if(Frac>=1)
{NumToShift=Math.floor(Frac);Frac-=NumToShift;if(NumToShift<this._buffer.length)
{this._previousValue=x3dom.fields.SFVec3f.copy(this._buffer[this._buffer.length-NumToShift]);for(C=this._buffer.length-1;C>=NumToShift;C--){this._buffer[C]=x3dom.fields.SFVec3f.copy(this._buffer[C-NumToShift]);}
for(C=0;C<NumToShift;C++)
{Alpha=C/NumToShift;this._buffer[C]=this._buffer[NumToShift].multiply(Alpha).add(this._vf.destination.multiply((1-Alpha)));}}
else
{this._previousValue=x3dom.fields.SFVec3f.copy((NumToShift==this._buffer.length)?this._buffer[0]:this._vf.destination);for(C=0;C<this._buffer.length;C++){this._buffer[C]=x3dom.fields.SFVec3f.copy(this._vf.destination);}}
this._bufferEndTime+=NumToShift*this._stepTime;}
return Frac;}}));x3dom.registerNodeType("PositionChaser2D","Followers",defineClass(x3dom.nodeTypes.X3DChaserNode,function(ctx){x3dom.nodeTypes.PositionChaser2D.superClass.call(this,ctx);this.addField_SFVec2f(ctx,'initialDestination',0,0);this.addField_SFVec2f(ctx,'initialValue',0,0);this.addField_SFVec2f(ctx,'value',0,0);this.addField_SFVec2f(ctx,'destination',0,0);this._buffer=new x3dom.fields.MFVec2f();this._previousValue=new x3dom.fields.SFVec2f(0,0);this._value=new x3dom.fields.SFVec2f(0,0);this.initialize();},{fieldChanged:function(fieldName)
{if(fieldName.indexOf("destination")>=0)
{this.initialize();this.updateBuffer(this._currTime);if(!this._vf.isActive){this.postMessage('isActive',true);}}
else if(fieldName.indexOf("value")>=0)
{this.initialize();this._previousValue.setValues(this._vf.value);for(var C=1;C<this._buffer.length;C++){this._buffer[C].setValues(this._vf.value);}
this.postMessage('value',this._vf.value);if(!this._vf.isActive){this.postMessage('isActive',true);}}},initialize:function()
{if(!this._initDone)
{this._initDone=true;this._vf.destination=x3dom.fields.SFVec2f.copy(this._vf.initialDestination);this._buffer.length=this._numSupports;this._buffer[0]=x3dom.fields.SFVec2f.copy(this._vf.initialDestination);for(var C=1;C<this._buffer.length;C++){this._buffer[C]=x3dom.fields.SFVec2f.copy(this._vf.initialValue);}
this._previousValue=x3dom.fields.SFVec2f.copy(this._vf.initialValue);this._stepTime=this._vf.duration/this._numSupports;var active=!this._buffer[0].equals(this._buffer[1],this._eps);if(this._vf.isActive!==active){this.postMessage('isActive',active);}}},tick:function(now)
{this.initialize();this._currTime=now;if(!this._bufferEndTime)
{this._bufferEndTime=now;this._value=x3dom.fields.SFVec2f.copy(this._vf.initialValue);this.postMessage('value',this._value);return true;}
var Frac=this.updateBuffer(now);var Output=x3dom.fields.SFVec2f.copy(this._previousValue);var DeltaIn=this._buffer[this._buffer.length-1].subtract(this._previousValue);var DeltaOut=DeltaIn.multiply(this.stepResponse((this._buffer.length-1+Frac)*this._stepTime));Output=Output.add(DeltaOut);for(var C=this._buffer.length-2;C>=0;C--)
{DeltaIn=this._buffer[C].subtract(this._buffer[C+1]);DeltaOut=DeltaIn.multiply(this.stepResponse((C+Frac)*this._stepTime));Output=Output.add(DeltaOut);}
if(!Output.equals(this._value,this._eps)){this._value.setValues(Output);this.postMessage('value',this._value);}
else{this.postMessage('isActive',false);}
return this._vf.isActive;},updateBuffer:function(now)
{var Frac=(now-this._bufferEndTime)/this._stepTime;var C;var NumToShift;var Alpha;if(Frac>=1)
{NumToShift=Math.floor(Frac);Frac-=NumToShift;if(NumToShift<this._buffer.length)
{this._previousValue=x3dom.fields.SFVec2f.copy(this._buffer[this._buffer.length-NumToShift]);for(C=this._buffer.length-1;C>=NumToShift;C--){this._buffer[C]=x3dom.fields.SFVec2f.copy(this._buffer[C-NumToShift]);}
for(C=0;C<NumToShift;C++)
{Alpha=C/NumToShift;this._buffer[C]=this._buffer[NumToShift].multiply(Alpha).add(this._vf.destination.multiply((1-Alpha)));}}
else
{this._previousValue=x3dom.fields.SFVec2f.copy((NumToShift==this._buffer.length)?this._buffer[0]:this._vf.destination);for(C=0;C<this._buffer.length;C++){this._buffer[C]=x3dom.fields.SFVec2f.copy(this._vf.destination);}}
this._bufferEndTime+=NumToShift*this._stepTime;}
return Frac;}}));x3dom.registerNodeType("PositionDamper","Followers",defineClass(x3dom.nodeTypes.X3DDamperNode,function(ctx){x3dom.nodeTypes.PositionDamper.superClass.call(this,ctx);this.addField_SFVec3f(ctx,'initialDestination',0,0,0);this.addField_SFVec3f(ctx,'initialValue',0,0,0);this.addField_SFVec3f(ctx,'value',0,0,0);this.addField_SFVec3f(ctx,'destination',0,0,0);this._value0=new x3dom.fields.SFVec3f(0,0,0);this._value1=new x3dom.fields.SFVec3f(0,0,0);this._value2=new x3dom.fields.SFVec3f(0,0,0);this._value3=new x3dom.fields.SFVec3f(0,0,0);this._value4=new x3dom.fields.SFVec3f(0,0,0);this._value5=new x3dom.fields.SFVec3f(0,0,0);this.initialize();},{fieldChanged:function(fieldName)
{if(fieldName==="tolerance")
{this._eps=this._vf.tolerance<0?0.001:this._vf.tolerance;}
else if(fieldName.indexOf("destination")>=0)
{if(!this._value0.equals(this._vf.destination,this._eps)){this._value0=this._vf.destination;if(!this._vf.isActive){this.postMessage('isActive',true);}}}
else if(fieldName.indexOf("value")>=0)
{this._value1.setValues(this._vf.value);this._value2.setValues(this._vf.value);this._value3.setValues(this._vf.value);this._value4.setValues(this._vf.value);this._value5.setValues(this._vf.value);this._lastTick=0;this.postMessage('value',this._value5);if(!this._vf.isActive){this._lastTick=0;this.postMessage('isActive',true);}}},initialize:function()
{this._value0.setValues(this._vf.initialDestination);this._value1.setValues(this._vf.initialValue);this._value2.setValues(this._vf.initialValue);this._value3.setValues(this._vf.initialValue);this._value4.setValues(this._vf.initialValue);this._value5.setValues(this._vf.initialValue);this._lastTick=0;var active=!this._value0.equals(this._value1,this._eps);if(this._vf.isActive!==active){this.postMessage('isActive',active);}},tick:function(now)
{if(!this._lastTick)
{this._lastTick=now;return false;}
var delta=now-this._lastTick;var alpha=Math.exp(-delta/this._vf.tau);this._value1=this._vf.order>0&&this._vf.tau?this._value0.add(this._value1.subtract(this._value0).multiply(alpha)):new x3dom.fields.SFVec3f(this._value0.x,this._value0.y,this._value0.z);this._value2=this._vf.order>1&&this._vf.tau?this._value1.add(this._value2.subtract(this._value1).multiply(alpha)):new x3dom.fields.SFVec3f(this._value1.x,this._value1.y,this._value1.z);this._value3=this._vf.order>2&&this._vf.tau?this._value2.add(this._value3.subtract(this._value2).multiply(alpha)):new x3dom.fields.SFVec3f(this._value2.x,this._value2.y,this._value2.z);this._value4=this._vf.order>3&&this._vf.tau?this._value3.add(this._value4.subtract(this._value3).multiply(alpha)):new x3dom.fields.SFVec3f(this._value3.x,this._value3.y,this._value3.z);this._value5=this._vf.order>4&&this._vf.tau?this._value4.add(this._value5.subtract(this._value4).multiply(alpha)):new x3dom.fields.SFVec3f(this._value4.x,this._value4.y,this._value4.z);var dist=this._value1.subtract(this._value0).length();if(this._vf.order>1)
{var dist2=this._value2.subtract(this._value1).length();if(dist2>dist){dist=dist2;}}
if(this._vf.order>2)
{var dist3=this._value3.subtract(this._value2).length();if(dist3>dist){dist=dist3;}}
if(this._vf.order>3)
{var dist4=this._value4.subtract(this._value3).length();if(dist4>dist){dist=dist4;}}
if(this._vf.order>4)
{var dist5=this._value5.subtract(this._value4).length();if(dist5>dist){dist=dist5;}}
if(dist<=this._eps)
{this._value1.setValues(this._value0);this._value2.setValues(this._value0);this._value3.setValues(this._value0);this._value4.setValues(this._value0);this._value5.setValues(this._value0);this.postMessage('value',this._value0);this.postMessage('isActive',false);this._lastTick=0;return false;}
this.postMessage('value',this._value5);this._lastTick=now;return true;}}));x3dom.registerNodeType("PositionDamper2D","Followers",defineClass(x3dom.nodeTypes.X3DDamperNode,function(ctx){x3dom.nodeTypes.PositionDamper2D.superClass.call(this,ctx);this.addField_SFVec2f(ctx,'initialDestination',0,0);this.addField_SFVec2f(ctx,'initialValue',0,0);this.addField_SFVec2f(ctx,'value',0,0);this.addField_SFVec2f(ctx,'destination',0,0);this._value0=new x3dom.fields.SFVec2f(0,0);this._value1=new x3dom.fields.SFVec2f(0,0);this._value2=new x3dom.fields.SFVec2f(0,0);this._value3=new x3dom.fields.SFVec2f(0,0);this._value4=new x3dom.fields.SFVec2f(0,0);this._value5=new x3dom.fields.SFVec2f(0,0);this.initialize();},{fieldChanged:function(fieldName)
{if(fieldName==="tolerance")
{this._eps=this._vf.tolerance<0?0.001:this._vf.tolerance;}
else if(fieldName.indexOf("destination")>=0)
{if(!this._value0.equals(this._vf.destination,this._eps)){this._value0=this._vf.destination;if(!this._vf.isActive){this.postMessage('isActive',true);}}}
else if(fieldName.indexOf("value")>=0)
{this._value1.setValues(this._vf.value);this._value2.setValues(this._vf.value);this._value3.setValues(this._vf.value);this._value4.setValues(this._vf.value);this._value5.setValues(this._vf.value);this._lastTick=0;this.postMessage('value',this._value5);if(!this._vf.isActive){this._lastTick=0;this.postMessage('isActive',true);}}},initialize:function()
{this._value0.setValues(this._vf.initialDestination);this._value1.setValues(this._vf.initialValue);this._value2.setValues(this._vf.initialValue);this._value3.setValues(this._vf.initialValue);this._value4.setValues(this._vf.initialValue);this._value5.setValues(this._vf.initialValue);this._lastTick=0;var active=!this._value0.equals(this._value1,this._eps);if(this._vf.isActive!==active){this.postMessage('isActive',active);}},tick:function(now)
{if(!this._lastTick)
{this._lastTick=now;return false;}
var delta=now-this._lastTick;var alpha=Math.exp(-delta/this._vf.tau);this._value1=this._vf.order>0&&this._vf.tau?this._value0.add(this._value1.subtract(this._value0).multiply(alpha)):new x3dom.fields.SFVec2f(this._value0.x,this._value0.y,this._value0.z);this._value2=this._vf.order>1&&this._vf.tau?this._value1.add(this._value2.subtract(this._value1).multiply(alpha)):new x3dom.fields.SFVec2f(this._value1.x,this._value1.y,this._value1.z);this._value3=this._vf.order>2&&this._vf.tau?this._value2.add(this._value3.subtract(this._value2).multiply(alpha)):new x3dom.fields.SFVec2f(this._value2.x,this._value2.y,this._value2.z);this._value4=this._vf.order>3&&this._vf.tau?this._value3.add(this._value4.subtract(this._value3).multiply(alpha)):new x3dom.fields.SFVec2f(this._value3.x,this._value3.y,this._value3.z);this._value5=this._vf.order>4&&this._vf.tau?this._value4.add(this._value5.subtract(this._value4).multiply(alpha)):new x3dom.fields.SFVec2f(this._value4.x,this._value4.y,this._value4.z);var dist=this._value1.subtract(this._value0).length();if(this._vf.order>1)
{var dist2=this._value2.subtract(this._value1).length();if(dist2>dist){dist=dist2;}}
if(this._vf.order>2)
{var dist3=this._value3.subtract(this._value2).length();if(dist3>dist){dist=dist3;}}
if(this._vf.order>3)
{var dist4=this._value4.subtract(this._value3).length();if(dist4>dist){dist=dist4;}}
if(this._vf.order>4)
{var dist5=this._value5.subtract(this._value4).length();if(dist5>dist){dist=dist5;}}
if(dist<=this._eps)
{this._value1.setValues(this._value0);this._value2.setValues(this._value0);this._value3.setValues(this._value0);this._value4.setValues(this._value0);this._value5.setValues(this._value0);this.postMessage('value',this._value0);this.postMessage('isActive',false);this._lastTick=0;return false;}
this.postMessage('value',this._value5);this._lastTick=now;return true;}}));x3dom.registerNodeType("ScalarChaser","Followers",defineClass(x3dom.nodeTypes.X3DChaserNode,function(ctx){x3dom.nodeTypes.ScalarChaser.superClass.call(this,ctx);this.addField_SFFloat(ctx,'initialDestination',0);this.addField_SFFloat(ctx,'initialValue',0);this.addField_SFFloat(ctx,'value',0);this.addField_SFFloat(ctx,'destination',0);this._buffer=[];this._previousValue=0;this._value=0;this.initialize();},{fieldChanged:function(fieldName)
{if(fieldName.indexOf("destination")>=0)
{this.initialize();this.updateBuffer(this._currTime);if(!this._vf.isActive){this.postMessage('isActive',true);}}
else if(fieldName.indexOf("value")>=0)
{this.initialize();this._previousValue=this._vf.value;for(var C=1;C<this._buffer.length;C++){this._buffer[C]=this._vf.value;}
this.postMessage('value',this._vf.value);if(!this._vf.isActive){this.postMessage('isActive',true);}}},initialize:function()
{if(!this._initDone)
{this._initDone=true;this._vf.destination=this._vf.initialDestination;this._buffer.length=this._numSupports;this._buffer[0]=this._vf.initialDestination;for(var C=1;C<this._buffer.length;C++){this._buffer[C]=this._vf.initialValue;}
this._previousValue=this._vf.initialValue;this._stepTime=this._vf.duration/this._numSupports;var active=(Math.abs(this._buffer[0]-this._buffer[1])>this._eps);if(this._vf.isActive!==active){this.postMessage('isActive',active);}}},tick:function(now)
{this.initialize();this._currTime=now;if(!this._bufferEndTime)
{this._bufferEndTime=now;this._value=this._vf.initialValue;this.postMessage('value',this._value);return true;}
var Frac=this.updateBuffer(now);var Output=this._previousValue;var DeltaIn=this._buffer[this._buffer.length-1]-this._previousValue;var DeltaOut=DeltaIn*(this.stepResponse((this._buffer.length-1+Frac)*this._stepTime));Output=Output+DeltaOut;for(var C=this._buffer.length-2;C>=0;C--)
{DeltaIn=this._buffer[C]-this._buffer[C+1];DeltaOut=DeltaIn*(this.stepResponse((C+Frac)*this._stepTime));Output=Output+DeltaOut;}
if(Math.abs(Output-this._value)>this._eps){this._value=Output;this.postMessage('value',this._value);}
else{this.postMessage('isActive',false);}
return this._vf.isActive;},updateBuffer:function(now)
{var Frac=(now-this._bufferEndTime)/this._stepTime;var C;var NumToShift;var Alpha;if(Frac>=1)
{NumToShift=Math.floor(Frac);Frac-=NumToShift;if(NumToShift<this._buffer.length)
{this._previousValue=this._buffer[this._buffer.length-NumToShift];for(C=this._buffer.length-1;C>=NumToShift;C--){this._buffer[C]=this._buffer[C-NumToShift];}
for(C=0;C<NumToShift;C++)
{Alpha=C/NumToShift;this._buffer[C]=this._buffer[NumToShift]*Alpha+this._vf.destination*(1-Alpha);}}
else
{this._previousValue=(NumToShift==this._buffer.length)?this._buffer[0]:this._vf.destination;for(C=0;C<this._buffer.length;C++){this._buffer[C]=this._vf.destination;}}
this._bufferEndTime+=NumToShift*this._stepTime;}
return Frac;}}));x3dom.registerNodeType("ScalarDamper","Followers",defineClass(x3dom.nodeTypes.X3DDamperNode,function(ctx){x3dom.nodeTypes.ScalarDamper.superClass.call(this,ctx);this.addField_SFFloat(ctx,'initialDestination',0);this.addField_SFFloat(ctx,'initialValue',0);this.addField_SFFloat(ctx,'value',0);this.addField_SFFloat(ctx,'destination',0);this._value0=0;this._value1=0;this._value2=0;this._value3=0;this._value4=0;this._value5=0;this.initialize();},{fieldChanged:function(fieldName)
{if(fieldName==="tolerance")
{this._eps=this._vf.tolerance<0?0.001:this._vf.tolerance;}
else if(fieldName.indexOf("destination")>=0)
{if(Math.abs(this._value0-this._vf.destination)>this._eps){this._value0=this._vf.destination;if(!this._vf.isActive){this.postMessage('isActive',true);}}}
else if(fieldName.indexOf("value")>=0)
{this._value1=this._vf.value;this._value2=this._vf.value;this._value3=this._vf.value;this._value4=this._vf.value;this._value5=this._vf.value;this._lastTick=0;this.postMessage('value',this._value5);if(!this._vf.isActive){this._lastTick=0;this.postMessage('isActive',true);}}},initialize:function()
{this._value0=this._vf.initialDestination;this._value1=this._vf.initialValue;this._value2=this._vf.initialValue;this._value3=this._vf.initialValue;this._value4=this._vf.initialValue;this._value5=this._vf.initialValue;this._lastTick=0;var active=(Math.abs(this._value0-this._value1)>this._eps);if(this._vf.isActive!==active){this.postMessage('isActive',active);}},tick:function(now)
{if(!this._lastTick)
{this._lastTick=now;return false;}
var delta=now-this._lastTick;var alpha=Math.exp(-delta/this._vf.tau);this._value1=this._vf.order>0&&this._vf.tau?this._value0+alpha*(this._value1-this._value0):this._value0;this._value2=this._vf.order>1&&this._vf.tau?this._value1+alpha*(this._value2-this._value1):this._value1;this._value3=this._vf.order>2&&this._vf.tau?this._value2+alpha*(this._value3-this._value2):this._value2;this._value4=this._vf.order>3&&this._vf.tau?this._value3+alpha*(this._value4-this._value3):this._value3;this._value5=this._vf.order>4&&this._vf.tau?this._value4+alpha*(this._value5-this._value4):this._value4;var dist=Math.abs(this._value1-this._value0);if(this._vf.order>1)
{var dist2=Math.abs(this._value2-this._value1);if(dist2>dist){dist=dist2;}}
if(this._vf.order>2)
{var dist3=Math.abs(this._value3-this._value2);if(dist3>dist){dist=dist3;}}
if(this._vf.order>3)
{var dist4=Math.abs(this._value4-this._value3);if(dist4>dist){dist=dist4;}}
if(this._vf.order>4)
{var dist5=Math.abs(this._value5-this._value4);if(dist5>dist){dist=dist5;}}
if(dist<=this._eps)
{this._value1=this._value0;this._value2=this._value0;this._value3=this._value0;this._value4=this._value0;this._value5=this._value0;this.postMessage('value',this._value0);this.postMessage('isActive',false);this._lastTick=0;return false;}
this.postMessage('value',this._value5);this._lastTick=now;return true;}}));x3dom.registerNodeType("CoordinateDamper","Followers",defineClass(x3dom.nodeTypes.X3DDamperNode,function(ctx){x3dom.nodeTypes.CoordinateDamper.superClass.call(this,ctx);this.addField_MFVec3f(ctx,'initialDestination',[]);this.addField_MFVec3f(ctx,'initialValue',[]);this.addField_MFVec3f(ctx,'value',[]);this.addField_MFVec3f(ctx,'destination',[]);x3dom.debug.logWarning("CoordinateDamper NYI");}));x3dom.registerNodeType("TexCoordDamper2D","Followers",defineClass(x3dom.nodeTypes.X3DDamperNode,function(ctx){x3dom.nodeTypes.TexCoordDamper2D.superClass.call(this,ctx);this.addField_MFVec2f(ctx,'initialDestination',[]);this.addField_MFVec2f(ctx,'initialValue',[]);this.addField_MFVec2f(ctx,'value',[]);this.addField_MFVec2f(ctx,'destination',[]);x3dom.debug.logWarning("TexCoordDamper2D NYI");}));x3dom.registerNodeType("X3DInterpolatorNode","Interpolation",defineClass(x3dom.nodeTypes.X3DChildNode,function(ctx){x3dom.nodeTypes.X3DInterpolatorNode.superClass.call(this,ctx);this.addField_MFFloat(ctx,'key',[]);this.addField_SFFloat(ctx,'set_fraction',0);},{linearInterp:function(time,interp){if(time<=this._vf.key[0])
return this._vf.keyValue[0];else if(time>=this._vf.key[this._vf.key.length-1])
return this._vf.keyValue[this._vf.key.length-1];for(var i=0;i<this._vf.key.length-1;++i){if((this._vf.key[i]<time)&&(time<=this._vf.key[i+1]))
return interp(this._vf.keyValue[i],this._vf.keyValue[i+1],(time-this._vf.key[i])/(this._vf.key[i+1]-this._vf.key[i]));}
return this._vf.keyValue[0];}}));x3dom.registerNodeType("OrientationInterpolator","Interpolation",defineClass(x3dom.nodeTypes.X3DInterpolatorNode,function(ctx){x3dom.nodeTypes.OrientationInterpolator.superClass.call(this,ctx);this.addField_MFRotation(ctx,'keyValue',[]);},{fieldChanged:function(fieldName)
{if(fieldName==="set_fraction")
{var value=this.linearInterp(this._vf.set_fraction,function(a,b,t){return a.slerp(b,t);});this.postMessage('value_changed',value);}}}));x3dom.registerNodeType("PositionInterpolator","Interpolation",defineClass(x3dom.nodeTypes.X3DInterpolatorNode,function(ctx){x3dom.nodeTypes.PositionInterpolator.superClass.call(this,ctx);this.addField_MFVec3f(ctx,'keyValue',[]);},{fieldChanged:function(fieldName)
{if(fieldName==="set_fraction")
{var value=this.linearInterp(this._vf.set_fraction,function(a,b,t){return a.multiply(1.0-t).add(b.multiply(t));});this.postMessage('value_changed',value);}}}));x3dom.registerNodeType("NormalInterpolator","Interpolation",defineClass(x3dom.nodeTypes.X3DInterpolatorNode,function(ctx){x3dom.nodeTypes.NormalInterpolator.superClass.call(this,ctx);this.addField_MFVec3f(ctx,'keyValue',[]);},{fieldChanged:function(fieldName)
{if(fieldName==="set_fraction")
{var value=this.linearInterp(this._vf.set_fraction,function(a,b,t){return a.multiply(1.0-t).add(b.multiply(t)).normalize();});this.postMessage('value_changed',value);}}}));x3dom.registerNodeType("ColorInterpolator","Interpolation",defineClass(x3dom.nodeTypes.X3DInterpolatorNode,function(ctx){x3dom.nodeTypes.ColorInterpolator.superClass.call(this,ctx);this.addField_MFColor(ctx,'keyValue',[]);},{fieldChanged:function(fieldName)
{if(fieldName==="set_fraction")
{var value=this.linearInterp(this._vf.set_fraction,function(a,b,t){return a.multiply(1.0-t).add(b.multiply(t));});this.postMessage('value_changed',value);}}}));x3dom.registerNodeType("ScalarInterpolator","Interpolation",defineClass(x3dom.nodeTypes.X3DInterpolatorNode,function(ctx){x3dom.nodeTypes.ScalarInterpolator.superClass.call(this,ctx);this.addField_MFFloat(ctx,'keyValue',[]);},{fieldChanged:function(fieldName)
{if(fieldName==="set_fraction")
{var value=this.linearInterp(this._vf.set_fraction,function(a,b,t){return(1.0-t)*a+t*b;});this.postMessage('value_changed',value);}}}));x3dom.registerNodeType("CoordinateInterpolator","Interpolation",defineClass(x3dom.nodeTypes.X3DInterpolatorNode,function(ctx){x3dom.nodeTypes.CoordinateInterpolator.superClass.call(this,ctx);this.addField_MFVec3f(ctx,'keyValue',[]);if(ctx&&ctx.xmlNode.hasAttribute('keyValue')){this._vf.keyValue=[];var arr=x3dom.fields.MFVec3f.parse(ctx.xmlNode.getAttribute('keyValue'));var key=this._vf.key.length>0?this._vf.key.length:1;var len=arr.length/key;for(var i=0;i<key;i++){var val=new x3dom.fields.MFVec3f();for(var j=0;j<len;j++){val.push(arr[i*len+j]);}
this._vf.keyValue.push(val);}}},{fieldChanged:function(fieldName)
{if(fieldName==="set_fraction")
{var value=this.linearInterp(this._vf.set_fraction,function(a,b,t){var val=new x3dom.fields.MFVec3f();for(var i=0;i<a.length;i++)
val.push(a[i].multiply(1.0-t).add(b[i].multiply(t)));return val;});this.postMessage('value_changed',value);}}}));x3dom.registerNodeType("SplinePositionInterpolator","Interpolation",defineClass(x3dom.nodeTypes.X3DInterpolatorNode,function(ctx){x3dom.nodeTypes.SplinePositionInterpolator.superClass.call(this,ctx);this.addField_MFVec3f(ctx,'keyValue',[]);this.addField_MFVec3f(ctx,'keyVelocity',[]);this.addField_SFBool(ctx,'closed',false);this.addField_SFBool(ctx,'normalizeVelocity',false);this.dtot=0.0;this.T0=[];this.T1=[];this.checkSanity=function(){var sane=(this._vf.key.length==this._vf.keyValue.length)&&((this._vf.key.length==this._vf.keyVelocity.length)||(this._vf.keyVelocity.length==2&&this._vf.key.length>=2)||(this._vf.keyVelocity.length==0));if(!sane)
x3dom.debug.logWarning("SplinePositionInterpolator Node: 'key' , 'keyValue' and/or 'keyVelocity' fields have inappropriate sizes");};this.calcDtot=function()
{this.dtot=0.0;for(var i=0;i<this._vf.key.length-1;i++)
{this.dtot+=Math.abs(this._vf.key[i]-this._vf.key[i+1]);}};this.calcAdjustedKeyVelocity=function()
{var i,Ti,F_plus_i,F_minus_i;var N=this._vf.key.length;if(this._vf.keyVelocity.length==N)
{for(i=0;i<N;i++)
{Ti=this._vf.keyVelocity[i];if(this._vf.normalizeVelocity)
Ti=Ti.multiply(this.dtot/Ti.length());F_plus_i=(i==0||i==N-1)?1.0:2.0*(this._vf.key[i]-this._vf.key[i-1])/(this._vf.key[i+1]-this._vf.key[i-1]);F_minus_i=(i==0||i==N-1)?1.0:2.0*(this._vf.key[i+1]-this._vf.key[i])/(this._vf.key[i+1]-this._vf.key[i-1]);this.T0[i]=Ti.multiply(F_plus_i);this.T1[i]=Ti.multiply(F_minus_i);}}
else if(this._vf.keyVelocity.length==2&&N>2)
{for(i=0;i<N;i++)
{if(i==0)
Ti=this._vf.keyVelocity[0];else if(i==N-1)
Ti=this._vf.keyVelocity[1];else
Ti=this._vf.keyValue[i+1].subtract(this._vf.keyValue[i-1]).multiply(0.5);if(this._vf.normalizeVelocity)
Ti=Ti.multiply(this.dtot/Ti.length());F_plus_i=(i==0||i==N-1)?1.0:2.0*(this._vf.key[i]-this._vf.key[i-1])/(this._vf.key[i+1]-this._vf.key[i-1]);F_minus_i=(i==0||i==N-1)?1.0:2.0*(this._vf.key[i+1]-this._vf.key[i])/(this._vf.key[i+1]-this._vf.key[i-1]);this.T0[i]=Ti.multiply(F_plus_i);this.T1[i]=Ti.multiply(F_minus_i);}}
else
{var closed=this._vf.closed&&this._vf.keyValue[0].equals(this._vf.keyValue[N-1],0.00001);for(i=0;i<N;i++)
{if((i==0||i==N-1)&&!closed)
{this.T0[i]=new x3dom.fields.SFVec3f(0,0,0);this.T1[i]=new x3dom.fields.SFVec3f(0,0,0);continue;}
else if((i==0||i==N-1)&&closed)
{Ti=this._vf.keyValue[1].subtract(this._vf.keyValue[N-2]).multiply(0.5);if(i==0){F_plus_i=2.0*(this._vf.key[0]-this._vf.key[N-2])/(this._vf.key[1]-this._vf.key[N-2]);F_minus_i=2.0*(this._vf.key[1]-this._vf.key[0])/(this._vf.key[1]-this._vf.key[N-2]);}
else{F_plus_i=2.0*(this._vf.key[N-1]-this._vf.key[N-2])/(this._vf.key[1]-this._vf.key[N-2]);F_minus_i=2.0*(this._vf.key[1]-this._vf.key[N-1])/(this._vf.key[1]-this._vf.key[N-2]);}
F_plus_i=2.0*(this._vf.key[N-1]-this._vf.key[N-2])/(this._vf.key[N-2]-this._vf.key[1]);F_minus_i=2.0*(this._vf.key[1]-this._vf.key[0])/(this._vf.key[N-2]-this._vf.key[1]);}
else
{Ti=this._vf.keyValue[i+1].subtract(this._vf.keyValue[i-1]).multiply(0.5);F_plus_i=2.0*(this._vf.key[i]-this._vf.key[i-1])/(this._vf.key[i+1]-this._vf.key[i-1]);F_minus_i=2.0*(this._vf.key[i+1]-this._vf.key[i])/(this._vf.key[i+1]-this._vf.key[i-1]);}
this.T0[i]=Ti.multiply(F_plus_i);this.T1[i]=Ti.multiply(F_minus_i);}}};this.checkSanity();this.calcDtot();this.calcAdjustedKeyVelocity();},{fieldChanged:function(fieldName)
{switch(fieldName)
{case'key':case'keyValue':case'keyVelocity':{this.checkSanity();this.calcDtot();this.calcAdjustedKeyVelocity();break;}
case'closed':case'normalizeVelocity':{this.calcAdjustedKeyVelocity();break;}
case'set_fraction':{var value;if(this._vf.key.length>0.0){if(this._vf.set_fraction<=this._vf.key[0])
value=x3dom.fields.SFVec3f.copy(this._vf.keyValue[0]);else if(this._vf.set_fraction>=this._vf.key[this._vf.key.length-1])
value=x3dom.fields.SFVec3f.copy(this._vf.keyValue[this._vf.key.length-1]);}
for(var i=0;i<this._vf.key.length-1;i++){if((this._vf.key[i]<this._vf.set_fraction)&&(this._vf.set_fraction<=this._vf.key[i+1])){var s=(this._vf.set_fraction-this._vf.key[i])/(this._vf.key[i+1]-this._vf.key[i]);var S_H=new x3dom.fields.SFVec4f(2.0*s*s*s-3.0*s*s+1.0,-2.0*s*s*s+3.0*s*s,s*s*s-2.0*s*s+s,s*s*s-s*s);value=new x3dom.fields.SFVec3f(S_H.x*this._vf.keyValue[i].x+S_H.y*this._vf.keyValue[i+1].x+S_H.z*this.T0[i].x+S_H.w*this.T1[i+1].x,S_H.x*this._vf.keyValue[i].y+S_H.y*this._vf.keyValue[i+1].y+S_H.z*this.T0[i].y+S_H.w*this.T1[i+1].y,S_H.x*this._vf.keyValue[i].z+S_H.y*this._vf.keyValue[i+1].z+S_H.z*this.T0[i].z+S_H.w*this.T1[i+1].z);break;}}
if(value!==undefined)
this.postMessage('value_changed',value);else
x3dom.debug.logWarning("SplinePositionInterpolator Node: value_changed is undefined!");}}}}));x3dom.registerNodeType("TimeSensor","Time",defineClass(x3dom.nodeTypes.X3DSensorNode,function(ctx){x3dom.nodeTypes.TimeSensor.superClass.call(this,ctx);if(ctx)
ctx.doc._nodeBag.timer.push(this);else
x3dom.debug.logWarning("TimeSensor: No runtime context found!");this.addField_SFTime(ctx,'cycleInterval',1);this.addField_SFBool(ctx,'loop',false);this.addField_SFTime(ctx,'startTime',0);this.addField_SFTime(ctx,'stopTime',0);this.addField_SFTime(ctx,'pauseTime',0);this.addField_SFTime(ctx,'resumeTime',0);this.addField_SFTime(ctx,'cycleTime',0);this.addField_SFTime(ctx,'elapsedTime',0);this.addField_SFFloat(ctx,'fraction_changed',0);this.addField_SFBool(ctx,'isActive',false);this.addField_SFBool(ctx,'isPaused',false);this.addField_SFTime(ctx,'time',0);this.addField_SFBool(ctx,'first',true);this.addField_SFFloat(ctx,'firstCycle',0.0);this._prevCycle=-1;this._lastTime=0;this._cycleStopTime=0;this._activatedTime=0;if(this._vf.startTime>0){this._updateCycleStopTime();}
this._backupStartTime=this._vf.startTime;this._backupStopTime=this._vf.stopTime;this._backupCycleInterval=this._vf.cycleInterval;},{tick:function(time)
{if(!this._vf.enabled){this._lastTime=time;return false;}
var isActive=(this._vf.cycleInterval>0&&time>=this._vf.startTime&&(time<this._vf.stopTime||this._vf.stopTime<=this._vf.startTime)&&(this._vf.loop==true||(this._vf.loop==false&&time<this._cycleStopTime)));if(isActive&&!this._vf.isActive){this.postMessage('isActive',true);this._activatedTime=time;}
if(isActive||this._vf.isActive){this.postMessage('elapsedTime',time-this._activatedTime);var isPaused=(time>=this._vf.pauseTime&&this._vf.pauseTime>this._vf.resumeTime);if(isPaused&&!this._vf.isPaused){this.postMessage('isPaused',true);this.postMessage('pauseTime',time);}else if(!isPaused&&this._vf.isPaused){this.postMessage('isPaused',false);this.postMessage('resumeTime',time);}
if(!isPaused){var cycleFrac=this._getCycleAt(time);var cycle=Math.floor(cycleFrac);var cycleTime=this._vf.startTime+cycle*this._vf.cycleInterval;var adjustTime=0;if(this._vf.stopTime>this._vf.startTime&&this._lastTime<this._vf.stopTime&&time>=this._vf.stopTime)
adjustTime=this._vf.stopTime;else if(this._lastTime<cycleTime&&time>=cycleTime)
adjustTime=cycleTime;if(adjustTime>0){time=adjustTime;cycleFrac=this._getCycleAt(time);cycle=Math.floor(cycleFrac);}
var fraction=cycleFrac-cycle;if(fraction<x3dom.fields.Eps){fraction=(this._lastTime<this._vf.startTime?0.0:1.0);this.postMessage('cycleTime',time);}
this.postMessage('fraction_changed',fraction);this.postMessage('time',time);}}
if(!isActive&&this._vf.isActive)
this.postMessage('isActive',false);this._lastTime=time;return true;},fieldChanged:function(fieldName)
{if(fieldName=="enabled"){if(!this._vf.enabled&&this._vf.isActive){this.postMessage('isActive',false);}}
else if(fieldName=="startTime"){if(this._vf.isActive){this._vf.startTime=this._backupStartTime;return;}
this._backupStartTime=this._vf.startTime;this._updateCycleStopTime();}
else if(fieldName=="stopTime"){if(this._vf.isActive&&this._vf.stopTime<=this._vf.startTime){this._vf.stopTime=this._backupStopTime;return;}
this._backupStopTime=this._vf.stopTime;}
else if(fieldName=="cycleInterval"){if(this._vf.isActive){this._vf.cycleInterval=this._backupCycleInterval;return;}
this._backupCycleInterval=this._vf.cycleInterval;}
else if(fieldName=="loop"){this._updateCycleStopTime();}},parentRemoved:function(parent)
{if(this._parentNodes.length===0){var doc=this.findX3DDoc();for(var i=0,n=doc._nodeBag.timer.length;i<n;i++){if(doc._nodeBag.timer[i]===this){doc._nodeBag.timer.splice(i,1);}}}},_getCycleAt:function(time)
{return Math.max(0.0,time-this._vf.startTime)/this._vf.cycleInterval;},_updateCycleStopTime:function()
{if(this._vf.loop==false){var now=new Date().getTime()/1000;var cycleToStop=Math.floor(this._getCycleAt(now))+1;this._cycleStopTime=this._vf.startTime+cycleToStop*this._vf.cycleInterval;}
else{this._cycleStopTime=0;}}}));x3dom.registerNodeType("X3DTimeDependentNode","Time",defineClass(x3dom.nodeTypes.X3DChildNode,function(ctx){x3dom.nodeTypes.X3DTimeDependentNode.superClass.call(this,ctx);this.addField_SFBool(ctx,'loop',false);}));x3dom.registerNodeType("Anchor","Networking",defineClass(x3dom.nodeTypes.X3DGroupingNode,function(ctx){x3dom.nodeTypes.Anchor.superClass.call(this,ctx);this.addField_MFString(ctx,'url',[]);this.addField_MFString(ctx,'parameter',[]);this.addField_SFString(ctx,'description',"");},{doIntersect:function(line){var isect=false;for(var i=0;i<this._childNodes.length;i++){if(this._childNodes[i]){isect=this._childNodes[i].doIntersect(line)||isect;}}
return isect;},handleTouch:function(){var url=this._vf.url.length?this._vf.url[0]:"";var aPos=url.search("#");var anchor="";if(aPos>=0)
anchor=url.slice(aPos+1);var param=this._vf.parameter.length?this._vf.parameter[0]:"";var tPos=param.search("target=");var target="";if(tPos>=0)
target=param.slice(tPos+7);x3dom.debug.logInfo("Anchor url="+url+", target="+target+", #viewpoint="+anchor);if(target.length!=0||target!="_self"){window.open(this._nameSpace.getURL(url),target);}
else{window.location=this._nameSpace.getURL(url);}}}));x3dom.registerNodeType("Inline","Networking",defineClass(x3dom.nodeTypes.X3DGroupingNode,function(ctx){x3dom.nodeTypes.Inline.superClass.call(this,ctx);this.addField_MFString(ctx,'url',[]);this.addField_SFBool(ctx,'load',true);this.addField_MFString(ctx,'nameSpaceName',[]);this.addField_SFBool(ctx,'mapDEFToID',false);this.initDone=false;this.count=0;this.numRetries=x3dom.nodeTypes.Inline.MaximumRetries;},{fieldChanged:function(fieldName)
{if(fieldName=="url"){for(var i=0;i<this._childNodes.length;i++)
{this.removeChild(this._childNodes[i]);}
if(this._vf.nameSpaceName.length!=0){var node=this._xmlNode;if(node&&node.hasChildNodes())
{while(node.childNodes.length>=1)
{node.removeChild(node.firstChild);}}}
this.loadInline();}
else if(fieldName=="render"){this.invalidateVolume();}},nodeChanged:function()
{if(!this.initDone){this.initDone=true;this.loadInline();}},fireEvents:function(eventType)
{if(this._xmlNode&&(this._xmlNode['on'+eventType]||this._xmlNode.hasAttribute('on'+eventType)||this._listeners[eventType]))
{var event={target:this._xmlNode,type:eventType,error:(eventType=="error")?"XMLHttpRequest Error":"",cancelBubble:false,stopPropagation:function(){this.cancelBubble=true;}};try{var attrib=this._xmlNode["on"+eventType];if(typeof(attrib)==="function"){attrib.call(this._xmlNode,event);}
else{var funcStr=this._xmlNode.getAttribute("on"+eventType);var func=new Function('event',funcStr);func.call(this._xmlNode,event);}
var list=this._listeners[eventType];if(list){for(var i=0;i<list.length;i++){list[i].call(this._xmlNode,event);}}}
catch(ex){x3dom.debug.logException(ex);}}},loadInline:function()
{var that=this;var xhr=new window.XMLHttpRequest();if(xhr.overrideMimeType)
xhr.overrideMimeType('text/xml');xhr.onreadystatechange=function()
{if(xhr.readyState!=4){return xhr;}
if(xhr.status===x3dom.nodeTypes.Inline.AwaitTranscoding){if(that.count<that.numRetries)
{that.count++;var refreshTime=+xhr.getResponseHeader("Refresh")||5;x3dom.debug.logInfo('XHR status: '+xhr.status+' - Await Transcoding ('+that.count+'/'+that.numRetries+'): '+'Next request in '+refreshTime+' seconds');window.setTimeout(function(){that._nameSpace.doc.downloadCount-=1;that.loadInline();},refreshTime*1000);return xhr;}
else
{x3dom.debug.logError('XHR status: '+xhr.status+' - Await Transcoding ('+that.count+'/'+that.numRetries+'): '+'No Retries left');that._nameSpace.doc.downloadCount-=1;that.count=0;return xhr;}}
else if((xhr.status!==200)&&(xhr.status!==0)){that.fireEvents("error");x3dom.debug.logError('XHR status: '+xhr.status+' - XMLHttpRequest requires web server running!');that._nameSpace.doc.downloadCount-=1;that.count=0;return xhr;}
else if((xhr.status==200)||(xhr.status==0)){that.count=0;}
x3dom.debug.logInfo('Inline: downloading '+that._vf.url[0]+' done.');var inlScene=null,newScene=null,nameSpace=null,xml=null;if(navigator.appName!="Microsoft Internet Explorer")
xml=xhr.responseXML;else
xml=new DOMParser().parseFromString(xhr.responseText,"text/xml");if(xml!==undefined&&xml!==null)
{inlScene=xml.getElementsByTagName('Scene')[0]||xml.getElementsByTagName('scene')[0];}
else{that.fireEvents("error");}
if(inlScene)
{var nsName=(that._vf.nameSpaceName.length!=0)?that._vf.nameSpaceName.toString().replace(' ',''):"";nameSpace=new x3dom.NodeNameSpace(nsName,that._nameSpace.doc);var url=that._vf.url.length?that._vf.url[0]:"";if((url[0]==='/')||(url.indexOf(":")>=0))
nameSpace.setBaseURL(url);else
nameSpace.setBaseURL(that._nameSpace.baseURL+url);newScene=nameSpace.setupTree(inlScene);that._nameSpace.addSpace(nameSpace);if(that._vf.nameSpaceName.length!=0)
{Array.forEach(inlScene.childNodes,function(childDomNode)
{if(childDomNode instanceof Element)
{setNamespace(that._vf.nameSpaceName,childDomNode,that._vf.mapDEFToID);that._xmlNode.appendChild(childDomNode);}});}}
else{if(xml&&xml.localName)
x3dom.debug.logError('No Scene in '+xml.localName);else
x3dom.debug.logError('No Scene in resource');}
var global=x3dom.getGlobal();if(that._childNodes.length>0&&that._childNodes[0]&&that._childNodes[0]._nameSpace)
that._nameSpace.removeSpace(that._childNodes[0]._nameSpace);while(that._childNodes.length!==0)
global['_remover']=that.removeChild(that._childNodes[0]);delete global['_remover'];if(newScene)
{that.addChild(newScene);that.invalidateVolume();that._nameSpace.doc.downloadCount-=1;that._nameSpace.doc.needRender=true;x3dom.debug.logInfo('Inline: added '+that._vf.url[0]+' to scene.');var theScene=that._nameSpace.doc._scene;if(theScene){theScene.invalidateVolume();window.setTimeout(function(){that.invalidateVolume();theScene.updateVolume();that._nameSpace.doc.needRender=true;},1000);}
that.fireEvents("load");}
newScene=null;nameSpace=null;inlScene=null;xml=null;return xhr;};if(this._vf.url.length&&this._vf.url[0].length)
{var xhrURI=this._nameSpace.getURL(this._vf.url[0]);xhr.open('GET',xhrURI,true);this._nameSpace.doc.downloadCount+=1;try{x3dom.RequestManager.addRequest(xhr);}
catch(ex){this.fireEvents("error");x3dom.debug.logError(this._vf.url[0]+": "+ex);}}}}));x3dom.nodeTypes.Inline.AwaitTranscoding=202;x3dom.nodeTypes.Inline.MaximumRetries=15;function setNamespace(prefix,childDomNode,mapDEFToID)
{if(childDomNode instanceof Element&&childDomNode.__setAttribute!==undefined){if(childDomNode.hasAttribute('id')){childDomNode.__setAttribute('id',prefix.toString().replace(' ','')+'__'+childDomNode.getAttribute('id'));}else if(childDomNode.hasAttribute('DEF')&&mapDEFToID){childDomNode.__setAttribute('id',prefix.toString().replace(' ','')+'__'+childDomNode.getAttribute('DEF'));if(!childDomNode.id)
childDomNode.id=childDomNode.__getAttribute('id');}}
if(childDomNode.hasChildNodes()){Array.forEach(childDomNode.childNodes,function(children){setNamespace(prefix,children,mapDEFToID);});}}
x3dom.registerNodeType("MultiPart","Networking",defineClass(x3dom.nodeTypes.Inline,function(ctx){x3dom.nodeTypes.MultiPart.superClass.call(this,ctx);this.addField_MFString(ctx,'urlIDMap',[]);this.addField_SFBool(ctx,'isPickable',true);this.addField_SFString(ctx,'sortType','auto');this.addField_SFBool(ctx,'solid',false);this.addField_SFInt32(ctx,'sortKey',0);this.addField_SFString(ctx,'initialVisibility','auto');this._idMap=null;this._inlineNamespace=null;this._highlightedParts=[];this._minId=0;this._maxId=0;this._lastId=-1;this._lastClickedId=-1;this._lastButton=0;this._identifierToPartId=[];this._identifierToAppId=[];this._visiblePartsPerShape=[];this._partVolume=[];this._partVisibility=[];this._originalColor=[];this._materials=[];},{fieldChanged:function(fieldName)
{if(fieldName=="url"){if(this._vf.nameSpaceName.length!=0){var node=this._xmlNode;if(node&&node.hasChildNodes())
{while(node.childNodes.length>=1)
{node.removeChild(node.firstChild);}}}
this.loadInline();}
else if(fieldName=="render"){this.invalidateVolume();}},nodeChanged:function()
{if(!this.initDone){this.initDone=true;this.loadIDMap();}},getVolume:function()
{var vol=this._graph.volume;if(!this.volumeValid()&&this._vf.render)
{for(var i=0;i<this._partVisibility.length;i++)
{if(!this._partVisibility[i])
continue;var childVol=this._partVolume[i];if(childVol&&childVol.isValid())
vol.extendBounds(childVol.min,childVol.max);}}
if(!vol.equals(this._graph.lastVolume))
{this._graph.lastVolume=x3dom.fields.BoxVolume.copy(vol);var event={target:this._xmlNode,type:"volumechanged",volume:x3dom.fields.BoxVolume.copy(vol)};this.callEvtHandler("onvolumechanged",event);}
return vol;},handleEvents:function(e)
{if(this._inlineNamespace){var colorMap=this._inlineNamespace.defMap["MultiMaterial_ColorMap"];var emissiveMap=this._inlineNamespace.defMap["MultiMaterial_EmissiveMap"];var specularMap=this._inlineNamespace.defMap["MultiMaterial_SpecularMap"];var visibilityMap=this._inlineNamespace.defMap["MultiMaterial_VisibilityMap"];if(e.pickedId==-1&&e.button!=0){this._lastClickedId=-1;this._lastButton=e.button;}else if(e.pickedId==-1&&e.button==0){this._lastClickedId=-1;this._lastButton=0;}
if(e.pickedId!=-1){e.part=new x3dom.Parts(this,[e.pickedId-this._minId],colorMap,emissiveMap,specularMap,visibilityMap);e.partID=this._idMap.mapping[e.pickedId-this._minId].name;e.type="mousemove";this.callEvtHandler("onmousemove",e);e.type="mouseover";this.callEvtHandler("onmouseover",e);if(!e.mouseup&&e.button&&e.button!=this._lastButton){e.type="mousedown";this._lastButton=e.button;if(this._lastClickedId==-1){this._lastClickedId=e.pickedId;}
this.callEvtHandler("onmousedown",e);}
if(e.mouseup||(this._lastButton!=0&&e.button==0)){e.type="mouseup";this.callEvtHandler("onmouseup",e);this._lastButton=0;if(e.pickedId==this._lastClickedId){this._lastClickedId=-1;e.type="click";this.callEvtHandler("onclick",e);}
this._lastClickedId=-1;}
if(e.pickedId!=this._lastId){if(this._lastId!=-1){e.part=new x3dom.Parts(this,[this._lastId-this._minId],colorMap,emissiveMap,specularMap,visibilityMap);e.partID=this._idMap.mapping[this._lastId-this._minId].name;e.type="mouseleave";this.callEvtHandler("onmouseleave",e);}
e.part=new x3dom.Parts(this,[e.pickedId-this._minId],colorMap,emissiveMap,specularMap,visibilityMap);e.partID=this._idMap.mapping[e.pickedId-this._minId].name;e.type="mouseenter";this.callEvtHandler("onmouseenter",e);this._lastId=e.pickedId;}
this._lastId=e.pickedId;}
else if(this._lastId!=-1){e.part=new x3dom.Parts(this,[this._lastId-this._minId],colorMap,emissiveMap,specularMap,visibilityMap);e.partID=this._idMap.mapping[this._lastId-this._minId].name;e.type="mouseout";this.callEvtHandler("onmouseout",e);e.type="mouseleave";this.callEvtHandler("onmouseleave",e);this._lastId=-1;}}},loadIDMap:function()
{if(this._vf.urlIDMap.length&&this._vf.urlIDMap[0].length)
{var i;var that=this;var idMapURI=this._nameSpace.getURL(this._vf.urlIDMap[0]);var xhr=new XMLHttpRequest();xhr.open("GET",idMapURI,true);xhr.onload=function()
{that._idMap=JSON.parse(this.responseText);if(that._nameSpace.doc._scene._multiPartMap==null){that._nameSpace.doc._scene._multiPartMap={numberOfIds:0,multiParts:[]};}
that._minId=that._nameSpace.doc._scene._multiPartMap.numberOfIds;that._maxId=that._minId+that._idMap.numberOfIDs-1;that._nameSpace.doc._scene._multiPartMap.numberOfIds+=that._idMap.numberOfIDs;that._nameSpace.doc._scene._multiPartMap.multiParts.push(that);for(i=0;i<that._idMap.mapping.length;i++)
{if(!that._identifierToPartId[that._idMap.mapping[i].name]){that._identifierToPartId[that._idMap.mapping[i].name]=[];}
if(!that._identifierToPartId[that._idMap.mapping[i].appearance]){that._identifierToPartId[that._idMap.mapping[i].appearance]=[];}
that._identifierToPartId[that._idMap.mapping[i].name].push(i);that._identifierToPartId[that._idMap.mapping[i].appearance].push(i);if(!that._partVolume[i]){var min=x3dom.fields.SFVec3f.parse(that._idMap.mapping[i].min);var max=x3dom.fields.SFVec3f.parse(that._idMap.mapping[i].max);that._partVolume[i]=new x3dom.fields.BoxVolume(min,max);}}
for(i=0;i<that._idMap.appearance.length;i++)
{that._identifierToAppId[that._idMap.appearance[i].name]=i;}
that.loadInline();};x3dom.RequestManager.addRequest(xhr);}},createMaterialData:function()
{var diffuseColor,transparency,specularColor,shininess,emissiveColor,ambientIntensity;var backDiffuseColor,backTransparency,backSpecularColor,backShininess,backEmissiveColor,backAmbientIntensity;var rgba_DT="",rgba_SS="",rgba_EA="";var rgba_DT_B="",rgba_SS_B="",rgba_EA_B="";var size=Math.ceil(Math.sqrt(this._idMap.numberOfIDs));size=x3dom.Utils.nextHighestPowerOfTwo(size);var sizeTwo=size*2.0;var diffuseTransparencyData=size+" "+sizeTwo+" 4";var specularShininessData=size+" "+sizeTwo+" 4";var emissiveAmbientIntensityData=size+" "+sizeTwo+" 4";for(var i=0;i<size*size;i++)
{if(i<this._idMap.mapping.length)
{var appName=this._idMap.mapping[i].appearance;var appID=this._identifierToAppId[appName];if(this._idMap.appearance[appID].material.ambientIntensity){ambientIntensity=this._idMap.appearance[appID].material.ambientIntensity}else{ambientIntensity="0.2";}
if(this._idMap.appearance[appID].material.backAmbientIntensity){backAmbientIntensity=this._idMap.appearance[appID].material.backAmbientIntensity}else{backAmbientIntensity=ambientIntensity;}
if(this._idMap.appearance[appID].material.diffuseColor){diffuseColor=this._idMap.appearance[appID].material.diffuseColor}else{diffuseColor="0.8 0.8 0.8";}
if(this._idMap.appearance[appID].material.backDiffuseColor){backDiffuseColor=this._idMap.appearance[appID].material.backDiffuseColor}else{backDiffuseColor=diffuseColor;}
if(this._idMap.appearance[appID].material.emissiveColor){emissiveColor=this._idMap.appearance[appID].material.emissiveColor}else{emissiveColor="0.0 0.0 0.0";}
if(this._idMap.appearance[appID].material.backEmissiveColor){backEmissiveColor=this._idMap.appearance[appID].material.backEmissiveColor}else{backEmissiveColor=emissiveColor;}
if(this._idMap.appearance[appID].material.shininess){shininess=this._idMap.appearance[appID].material.shininess;}else{shininess="0.2";}
if(this._idMap.appearance[appID].material.backShininess){backShininess=this._idMap.appearance[appID].material.backShininess;}else{backShininess=shininess;}
if(this._idMap.appearance[appID].material.specularColor){specularColor=this._idMap.appearance[appID].material.specularColor;}else{specularColor="0 0 0";}
if(this._idMap.appearance[appID].material.backSpecularColor){backSpecularColor=this._idMap.appearance[appID].material.backSpecularColor;}else{backSpecularColor=specularColor;}
if(this._idMap.appearance[appID].material.transparency){transparency=this._idMap.appearance[appID].material.transparency;}else{transparency=0.0;}
if(this._idMap.appearance[appID].material.backTransparency){backTransparency=this._idMap.appearance[appID].material.backTransparency;}else{backTransparency=transparency;}
rgba_DT+=" "+x3dom.fields.SFColorRGBA.parse(diffuseColor+" "+(1.0-transparency)).toUint();rgba_SS+=" "+x3dom.fields.SFColorRGBA.parse(specularColor+" "+shininess).toUint();rgba_EA+=" "+x3dom.fields.SFColorRGBA.parse(emissiveColor+" "+ambientIntensity).toUint();rgba_DT_B+=" "+x3dom.fields.SFColorRGBA.parse(backDiffuseColor+" "+(1.0-backTransparency)).toUint();rgba_SS_B+=" "+x3dom.fields.SFColorRGBA.parse(backSpecularColor+" "+backShininess).toUint();rgba_EA_B+=" "+x3dom.fields.SFColorRGBA.parse(backEmissiveColor+" "+backAmbientIntensity).toUint();this._originalColor[i]=rgba_DT;this._materials[i]=new x3dom.MultiMaterial({"ambientIntensity":ambientIntensity,"diffuseColor":x3dom.fields.SFColor.parse(diffuseColor),"emissiveColor":x3dom.fields.SFColor.parse(emissiveColor),"shininess":shininess,"specularColor":x3dom.fields.SFColor.parse(specularColor),"transparency":transparency,"backAmbientIntensity":backAmbientIntensity,"backDiffuseColor":x3dom.fields.SFColor.parse(backDiffuseColor),"backEmissiveColor":x3dom.fields.SFColor.parse(backEmissiveColor),"backShininess":backShininess,"backSpecularColor":x3dom.fields.SFColor.parse(backSpecularColor),"backTransparency":backTransparency});}
else
{rgba_DT+=" 255";rgba_SS+=" 255";rgba_EA+=" 255";rgba_DT_B+=" 255";rgba_SS_B+=" 255";rgba_EA_B+=" 255";}}
diffuseTransparencyData+=rgba_DT+rgba_DT_B;specularShininessData+=rgba_SS+rgba_SS_B;emissiveAmbientIntensityData+=rgba_EA+rgba_EA_B;return{"diffuseTransparency":diffuseTransparencyData,"specularShininess":specularShininessData,"emissiveAmbientIntensity":emissiveAmbientIntensityData};},createVisibilityData:function()
{var i,j;var size=Math.ceil(Math.sqrt(this._idMap.numberOfIDs));size=x3dom.Utils.nextHighestPowerOfTwo(size);var visibilityData=size+" "+size+" 1";for(i=0;i<size*size;i++)
{if(i<this._idMap.mapping.length)
{if(this._vf.initialVisibility=='auto')
{visibilityData+=" 255";if(!this._partVisibility[i]){this._partVisibility[i]=true;}
for(j=0;j<this._idMap.mapping[i].usage.length;j++)
{if(!this._visiblePartsPerShape[this._idMap.mapping[i].usage[j]]){this._visiblePartsPerShape[this._idMap.mapping[i].usage[j]]={val:0,max:0};}
this._visiblePartsPerShape[this._idMap.mapping[i].usage[j]].val++;this._visiblePartsPerShape[this._idMap.mapping[i].usage[j]].max++;}}
else if(this._vf.initialVisibility=='visible')
{visibilityData+=" 255";if(!this._partVisibility[i]){this._partVisibility[i]=true;}
for(j=0;j<this._idMap.mapping[i].usage.length;j++)
{if(!this._visiblePartsPerShape[this._idMap.mapping[i].usage[j]]){this._visiblePartsPerShape[this._idMap.mapping[i].usage[j]]={val:0,max:0};}
this._visiblePartsPerShape[this._idMap.mapping[i].usage[j]].val++;this._visiblePartsPerShape[this._idMap.mapping[i].usage[j]].max++;}}
else if(this._vf.initialVisibility=='invisible')
{visibilityData+=" 0";if(!this._partVisibility[i]){this._partVisibility[i]=false;}
for(j=0;j<this._idMap.mapping[i].usage.length;j++)
{if(!this._visiblePartsPerShape[this._idMap.mapping[i].usage[j]]){this._visiblePartsPerShape[this._idMap.mapping[i].usage[j]]={val:0,max:0};}
this._visiblePartsPerShape[this._idMap.mapping[i].usage[j]].max++;}}}
else
{visibilityData+=" 0";}}
return visibilityData;},replaceMaterials:function(inlScene)
{var css,shapeDEF,materialData,visibilityData,appearance;var firstMat=true;if(inlScene&&inlScene.hasChildNodes())
{materialData=this.createMaterialData();visibilityData=this.createVisibilityData();var shapes=inlScene.getElementsByTagName("Shape");for(var s=0;s<shapes.length;s++)
{shapeDEF=shapes[s].getAttribute("DEF")||shapes[s].getAttribute("def");if(shapeDEF&&this._visiblePartsPerShape[shapeDEF]&&this._visiblePartsPerShape[shapeDEF].val==0)
{shapes[s].setAttribute("render","false");}
shapes[s].setAttribute("idOffset",this._minId);shapes[s].setAttribute("isPickable",this._vf.isPickable);var geometries=shapes[s].getElementsByTagName("BinaryGeometry");if(geometries&&geometries.length){for(var g=0;g<geometries.length;g++){geometries[g].setAttribute("solid",this._vf.solid);}}
var appearances=shapes[s].getElementsByTagName("Appearance");if(appearances.length)
{for(var a=0;a<appearances.length;a++)
{appearances[a].removeAttribute("DEF");appearances[a].removeAttribute("USE");appearances[a].setAttribute("sortType",this._vf.sortType);appearances[a].setAttribute("sortKey",this._vf.sortKey);var materials=appearances[a].getElementsByTagName("Material");if(materials.length)
{if(firstMat){firstMat=false;css=document.createElement("CommonSurfaceShader");css.setAttribute("DEF","MultiMaterial");var ptDA=document.createElement("PixelTexture");ptDA.setAttribute("containerField","multiDiffuseAlphaTexture");ptDA.setAttribute("id","MultiMaterial_ColorMap");ptDA.setAttribute("image",materialData.diffuseTransparency);var ptEA=document.createElement("PixelTexture");ptEA.setAttribute("containerField","multiEmissiveAmbientTexture");ptEA.setAttribute("id","MultiMaterial_EmissiveMap");ptEA.setAttribute("image",materialData.emissiveAmbientIntensity);var ptSS=document.createElement("PixelTexture");ptSS.setAttribute("containerField","multiSpecularShininessTexture");ptSS.setAttribute("id","MultiMaterial_SpecularMap");ptSS.setAttribute("image",materialData.specularShininess);var ptV=document.createElement("PixelTexture");ptV.setAttribute("containerField","multiVisibilityTexture");ptV.setAttribute("id","MultiMaterial_VisibilityMap");ptV.setAttribute("image",visibilityData);css.appendChild(ptDA);css.appendChild(ptEA);css.appendChild(ptSS);css.appendChild(ptV);}
else
{css=document.createElement("CommonSurfaceShader");css.setAttribute("USE","MultiMaterial");}
appearances[a].replaceChild(css,materials[0]);}
else
{if(firstMat){firstMat=false;css=document.createElement("CommonSurfaceShader");css.setAttribute("DEF","MultiMaterial");var ptDA=document.createElement("PixelTexture");ptDA.setAttribute("containerField","multiDiffuseAlphaTexture");ptDA.setAttribute("id","MultiMaterial_ColorMap");ptDA.setAttribute("image",materialData.diffuseTransparency);var ptEA=document.createElement("PixelTexture");ptEA.setAttribute("containerField","multiEmissiveAmbientTexture");ptEA.setAttribute("id","MultiMaterial_EmissiveMap");ptEA.setAttribute("image",materialData.emissiveAmbientIntensity);var ptSS=document.createElement("PixelTexture");ptSS.setAttribute("containerField","multiSpecularShininessTexture");ptSS.setAttribute("id","MultiMaterial_SpecularMap");ptSS.setAttribute("image",materialData.specularShininess);var ptV=document.createElement("PixelTexture");ptV.setAttribute("containerField","multiVisibilityTexture");ptV.setAttribute("id","MultiMaterial_VisibilityMap");ptV.setAttribute("image",visibilityData);css.appendChild(ptDA);css.appendChild(ptEA);css.appendChild(ptSS);css.appendChild(ptV);}
else
{css=document.createElement("CommonSurfaceShader");css.setAttribute("USE","MultiMaterial");}
appearances[a].appendChild(css);}}}
else
{appearance=document.createElement("Appearance");if(firstMat){firstMat=false;css=document.createElement("CommonSurfaceShader");css.setAttribute("DEF","MultiMaterial");var ptDA=document.createElement("PixelTexture");ptDA.setAttribute("containerField","multiDiffuseAlphaTexture");ptDA.setAttribute("id","MultiMaterial_ColorMap");ptDA.setAttribute("image",materialData.diffuseTransparency);var ptEA=document.createElement("PixelTexture");ptEA.setAttribute("containerField","multiEmissiveAmbientTexture");ptEA.setAttribute("id","MultiMaterial_EmissiveMap");ptEA.setAttribute("image",materialData.emissiveAmbientIntensity);var ptSS=document.createElement("PixelTexture");ptSS.setAttribute("containerField","multiSpecularShininessTexture");ptSS.setAttribute("id","MultiMaterial_SpecularMap");ptSS.setAttribute("image",materialData.specularShininess);var ptV=document.createElement("PixelTexture");ptV.setAttribute("containerField","multiVisibilityTexture");ptV.setAttribute("id","MultiMaterial_VisibilityMap");ptV.setAttribute("image",visibilityData);css.appendChild(ptDA);css.appendChild(ptEA);css.appendChild(ptSS);css.appendChild(ptV);}
else
{css=document.createElement("CommonSurfaceShader");css.setAttribute("USE","MultiMaterial");}
appearance.appendChild(css);geometries[g].appendChild(appearance);}}}},appendAPI:function()
{var multiPart=this;this._xmlNode.getIdList=function()
{var i,ids=[];for(i=0;i<multiPart._idMap.mapping.length;i++){ids.push(multiPart._idMap.mapping[i].name);}
return ids;};this._xmlNode.getAppearanceIdList=function()
{var i,ids=[];for(i=0;i<multiPart._idMap.appearance.length;i++){ids.push(multiPart._idMap.appearance[i].name);}
return ids;};this._xmlNode.getParts=function(selector)
{var i,m;var selection=[];if(selector==undefined){for(m=0;m<multiPart._idMap.mapping.length;m++){selection.push(m);}}else if(selector instanceof Array){for(i=0;i<selector.length;i++){if(multiPart._identifierToPartId[selector[i]]){selection=selection.concat(multiPart._identifierToPartId[selector[i]]);}}}else if(selector instanceof RegExp){for(var key in multiPart._identifierToPartId){if(key.match(selector)){selection=selection.concat(multiPart._identifierToPartId[key]);}}}
var colorMap=multiPart._inlineNamespace.defMap["MultiMaterial_ColorMap"];var emissiveMap=multiPart._inlineNamespace.defMap["MultiMaterial_EmissiveMap"];var specularMap=multiPart._inlineNamespace.defMap["MultiMaterial_SpecularMap"];var visibilityMap=multiPart._inlineNamespace.defMap["MultiMaterial_VisibilityMap"];if(selection.length==0){return null;}else{return new x3dom.Parts(multiPart,selection,colorMap,emissiveMap,specularMap,visibilityMap);}};this._xmlNode.getPartsByRect=function(left,right,bottom,top)
{var viewarea=multiPart._nameSpace.doc._viewarea;var viewpoint=viewarea._scene.getViewpoint();var origViewMatrix=viewarea.getViewMatrix();var origProjMatrix=viewarea.getProjectionMatrix();var upDir=new x3dom.fields.SFVec3f(origViewMatrix._01,origViewMatrix._11,origViewMatrix._21);var viewDir=new x3dom.fields.SFVec3f(origViewMatrix._02,origViewMatrix._12,origViewMatrix._22);var pos=new x3dom.fields.SFVec3f(origViewMatrix._03,origViewMatrix._13,origViewMatrix._23);var normalizedLeft=(left-viewarea._width/2)/(viewarea._width/2);var normalizedRight=(right-viewarea._width/2)/(viewarea._width/2);var normalizedTop=(top-viewarea._height/2)/(viewarea._height/2);var normalizedBottom=(bottom-viewarea._height/2)/(viewarea._height/2);var fov=viewpoint._vf.fieldOfView;var factorH=Math.tan(fov/2)*viewpoint._zNear;var factorW=Math.tan(fov/2)*viewpoint._lastAspect*viewpoint._zNear;var projMatrix=x3dom.fields.SFMatrix4f.perspectiveFrustum(normalizedLeft*factorW,normalizedRight*factorW,normalizedBottom*factorH,normalizedTop*factorH,viewpoint.getNear(),viewpoint.getFar());var viewMatrix=x3dom.fields.SFMatrix4f.lookAt(pos,pos.subtract(viewDir.multiply(5.0)),upDir);var frustum=new x3dom.fields.FrustumVolume(projMatrix.mult(viewMatrix));var selection=[];var volumes=this._x3domNode._partVolume;for(id in volumes){if(!volumes.hasOwnProperty(id))
continue;var intersect=frustum.intersect(volumes[id],0);if(intersect>0)
selection.push(id);}
var colorMap=multiPart._inlineNamespace.defMap["MultiMaterial_ColorMap"];var emissiveMap=multiPart._inlineNamespace.defMap["MultiMaterial_EmissiveMap"];var specularMap=multiPart._inlineNamespace.defMap["MultiMaterial_SpecularMap"];var visibilityMap=multiPart._inlineNamespace.defMap["MultiMaterial_VisibilityMap"];if(selection.length==0){return null;}else{return new x3dom.Parts(multiPart,selection,colorMap,emissiveMap,specularMap,visibilityMap);}};},loadInline:function()
{var that=this;var xhr=new window.XMLHttpRequest();if(xhr.overrideMimeType)
xhr.overrideMimeType('text/xml');xhr.onreadystatechange=function()
{if(xhr.readyState!=4){return xhr;}
if(xhr.status===x3dom.nodeTypes.Inline.AwaitTranscoding){if(that.count<that.numRetries)
{that.count++;var refreshTime=+xhr.getResponseHeader("Refresh")||5;x3dom.debug.logInfo('XHR status: '+xhr.status+' - Await Transcoding ('+that.count+'/'+that.numRetries+'): '+'Next request in '+refreshTime+' seconds');window.setTimeout(function(){that._nameSpace.doc.downloadCount-=1;that.loadInline();},refreshTime*1000);return xhr;}
else
{x3dom.debug.logError('XHR status: '+xhr.status+' - Await Transcoding ('+that.count+'/'+that.numRetries+'): '+'No Retries left');that._nameSpace.doc.downloadCount-=1;that.count=0;return xhr;}}
else if((xhr.status!==200)&&(xhr.status!==0)){that.fireEvents("error");x3dom.debug.logError('XHR status: '+xhr.status+' - XMLHttpRequest requires web server running!');that._nameSpace.doc.downloadCount-=1;that.count=0;return xhr;}
else if((xhr.status==200)||(xhr.status==0)){that.count=0;}
x3dom.debug.logInfo('Inline: downloading '+that._vf.url[0]+' done.');var inlScene=null,newScene=null,nameSpace=null,xml=null;if(navigator.appName!="Microsoft Internet Explorer")
xml=xhr.responseXML;else
xml=new DOMParser().parseFromString(xhr.responseText,"text/xml");if(xml!==undefined&&xml!==null)
{inlScene=xml.getElementsByTagName('Scene')[0]||xml.getElementsByTagName('scene')[0];}
else{that.fireEvents("error");}
if(inlScene)
{var nsDefault="ns"+that._nameSpace.childSpaces.length;var nsName=(that._vf.nameSpaceName.length!=0)?that._vf.nameSpaceName.toString().replace(' ',''):nsDefault;that._inlineNamespace=new x3dom.NodeNameSpace(nsName,that._nameSpace.doc);var url=that._vf.url.length?that._vf.url[0]:"";if((url[0]==='/')||(url.indexOf(":")>=0))
{that._inlineNamespace.setBaseURL(url);}
else
{that._inlineNamespace.setBaseURL(that._nameSpace.baseURL+url);}
that.replaceMaterials(inlScene);newScene=that._inlineNamespace.setupTree(inlScene);that._nameSpace.addSpace(that._inlineNamespace);if(that._vf.nameSpaceName.length!=0)
{Array.forEach(inlScene.childNodes,function(childDomNode)
{if(childDomNode instanceof Element)
{setNamespace(that._vf.nameSpaceName,childDomNode,that._vf.mapDEFToID);that._xmlNode.appendChild(childDomNode);}});}}
else{if(xml&&xml.localName){x3dom.debug.logError('No Scene in '+xml.localName);}else{x3dom.debug.logError('No Scene in resource');}}
var global=x3dom.getGlobal();if(that._childNodes.length>0&&that._childNodes[0]&&that._childNodes[0]._nameSpace){that._nameSpace.removeSpace(that._childNodes[0]._nameSpace);}
while(that._childNodes.length!==0){global['_remover']=that.removeChild(that._childNodes[0]);}
delete global['_remover'];if(newScene)
{that.addChild(newScene);that.invalidateVolume();that._nameSpace.doc.downloadCount-=1;that._nameSpace.doc.needRender=true;x3dom.debug.logInfo('Inline: added '+that._vf.url[0]+' to scene.');var theScene=that._nameSpace.doc._scene;if(theScene){theScene.invalidateVolume();window.setTimeout(function(){that.invalidateVolume();theScene.updateVolume();that._nameSpace.doc.needRender=true;},1000);}
that.appendAPI();that.fireEvents("load");}
newScene=null;inlScene=null;xml=null;return xhr;};if(this._vf.url.length&&this._vf.url[0].length)
{var xhrURI=this._nameSpace.getURL(this._vf.url[0]);xhr.open('GET',xhrURI,true);this._nameSpace.doc.downloadCount+=1;try{x3dom.RequestManager.addRequest(xhr);}
catch(ex){this.fireEvents("error");x3dom.debug.logError(this._vf.url[0]+": "+ex);}}}}));x3dom.registerNodeType("X3DBackgroundNode","EnvironmentalEffects",defineClass(x3dom.nodeTypes.X3DBindableNode,function(ctx){x3dom.nodeTypes.X3DBackgroundNode.superClass.call(this,ctx);var trans=(ctx&&ctx.autoGen)?1:0;this.addField_SFString(ctx,'crossOrigin','');this.addField_MFColor(ctx,'groundColor',[]);this.addField_MFFloat(ctx,'groundAngle',[]);this.addField_MFColor(ctx,'skyColor',[new x3dom.fields.SFColor(0,0,0)]);this.addField_MFFloat(ctx,'skyAngle',[]);this.addField_SFFloat(ctx,'transparency',trans);this._dirty=true;},{getSkyColor:function(){return new x3dom.fields.SFColor(0,0,0);},getTransparency:function(){return 0;},getTexUrl:function(){return[];}}));x3dom.registerNodeType("X3DFogNode","EnvironmentalEffects",defineClass(x3dom.nodeTypes.X3DBindableNode,function(ctx){x3dom.nodeTypes.X3DFogNode.superClass.call(this,ctx);this.addField_SFColor(ctx,'color',1,1,1);this.addField_SFString(ctx,'fogType',"LINEAR");this.addField_SFFloat(ctx,'visibilityRange',0);},{}));x3dom.registerNodeType("Fog","EnvironmentalEffects",defineClass(x3dom.nodeTypes.X3DFogNode,function(ctx){x3dom.nodeTypes.Fog.superClass.call(this,ctx);},{}));x3dom.registerNodeType("Background","EnvironmentalEffects",defineClass(x3dom.nodeTypes.X3DBackgroundNode,function(ctx){x3dom.nodeTypes.Background.superClass.call(this,ctx);this.addField_MFString(ctx,'backUrl',[]);this.addField_MFString(ctx,'bottomUrl',[]);this.addField_MFString(ctx,'frontUrl',[]);this.addField_MFString(ctx,'leftUrl',[]);this.addField_MFString(ctx,'rightUrl',[]);this.addField_MFString(ctx,'topUrl',[]);this.addField_SFBool(ctx,'scaling',false);},{fieldChanged:function(fieldName)
{if(fieldName.indexOf("Url")>0||fieldName=="transparency"||fieldName.search("sky")>=0||fieldName.search("ground")>=0){this._dirty=true;}
else if(fieldName.indexOf("bind")>=0){this.bind(this._vf.bind);}},getSkyColor:function(){return this._vf.skyColor;},getGroundColor:function(){return this._vf.groundColor;},getTransparency:function(){return this._vf.transparency;},getTexUrl:function(){return[this._nameSpace.getURL(this._vf.backUrl[0]),this._nameSpace.getURL(this._vf.frontUrl[0]),this._nameSpace.getURL(this._vf.bottomUrl[0]),this._nameSpace.getURL(this._vf.topUrl[0]),this._nameSpace.getURL(this._vf.leftUrl[0]),this._nameSpace.getURL(this._vf.rightUrl[0])];}}));x3dom.registerNodeType("X3DEnvironmentNode","EnvironmentalEffects",defineClass(x3dom.nodeTypes.X3DBindableNode,function(ctx){x3dom.nodeTypes.X3DEnvironmentNode.superClass.call(this,ctx);}));x3dom.registerNodeType("Environment","EnvironmentalEffects",defineClass(x3dom.nodeTypes.X3DEnvironmentNode,function(ctx){x3dom.nodeTypes.Environment.superClass.call(this,ctx);this.addField_SFBool(ctx,'sortTrans',true);this.addField_SFBool(ctx,'shadowExcludeTransparentObjects',false);this.addField_SFString(ctx,'gammaCorrectionDefault',"linear");this.addField_SFBool(ctx,'frustumCulling',true);this.addField_SFBool(ctx,'smallFeatureCulling',false);this.addField_SFFloat(ctx,'smallFeatureThreshold',1.0);this.addField_SFBool(ctx,'occlusionCulling',false);this.addField_SFFloat(ctx,'occlusionVisibilityThreshold',0.0);this.addField_SFBool(ctx,'lowPriorityCulling',false);this.addField_SFFloat(ctx,'lowPriorityThreshold',1.0);this.addField_SFBool(ctx,'tessellationDetailCulling',false);this.addField_SFFloat(ctx,'tessellationErrorThreshold',0.0);this.addField_SFBool(ctx,'enableARC',false);this.addField_SFFloat(ctx,'minFrameRate',1.0);this.addField_SFFloat(ctx,'maxFrameRate',62.5);this.addField_SFFloat(ctx,'userDataFactor',-1);this.addField_SFFloat(ctx,'smallFeatureFactor',-1);this.addField_SFFloat(ctx,'occlusionVisibilityFactor',-1);this.addField_SFFloat(ctx,'lowPriorityFactor',-1);this.addField_SFFloat(ctx,'tessellationErrorFactor',-1);this.addField_SFBool(ctx,'SSAO',false);this.addField_SFFloat(ctx,'SSAOradius',0.7);this.addField_SFFloat(ctx,'SSAOamount',0.3);this.addField_SFInt32(ctx,'SSAOrandomTextureSize',4);this.addField_SFInt32(ctx,'SSAOblurDepthTreshold',1);this._validGammaCorrectionTypes=["none","fastlinear","linear"];this.checkSanity();},{checkSanity:function()
{var checkParam=function(flag,value,defaultOn,defaultOff)
{if(flag&&(value==defaultOff))
return defaultOn;if(!flag&&(value!=defaultOff))
return defaultOff;return value;};this._smallFeatureThreshold=checkParam(this._vf.smallFeatureCulling,this._vf.smallFeatureThreshold,10,0);this._lowPriorityThreshold=checkParam(this._vf.lowPriorityCulling,this._vf.lowPriorityThreshold,0.5,1);this._occlusionVisibilityThreshold=checkParam(this._vf.occlusionCulling,this._vf.occlusionVisibilityThreshold,1,0);this._tessellationErrorThreshold=checkParam(this._vf.tessellationDetailCulling,this._vf.tessellationErrorThreshold,1,0);var checkGamma=function(field,that){field=field.toLowerCase();if(that._validGammaCorrectionTypes.indexOf(field)>-1){return field;}
else{x3dom.debug.logWarning(field+" gammaCorrectionDefault may only be linear (default), fastLinear, or none");return that._validGammaCorrectionTypes[0];}};this._vf.gammaCorrectionDefault=checkGamma(this._vf.gammaCorrectionDefault,this);}}));x3dom.registerNodeType("X3DViewpointNode","Navigation",defineClass(x3dom.nodeTypes.X3DBindableNode,function(ctx){x3dom.nodeTypes.X3DViewpointNode.superClass.call(this,ctx);if(ctx&&ctx.xmlNode){var domNode=ctx.xmlNode;if(!domNode.resetView&&!domNode.getFieldOfView&&!domNode.getNear&&!domNode.getFar)
{domNode.resetView=function(){var that=this._x3domNode;that.resetView();that._nameSpace.doc.needRender=true;};domNode.getFieldOfView=function(){return this._x3domNode.getFieldOfView();};domNode.getNear=function(){return this._x3domNode.getNear();};domNode.getFar=function(){return this._x3domNode.getFar();};}}},{activate:function(prev){var viewarea=this._nameSpace.doc._viewarea;if(prev&&this._bindAnimation){viewarea.animateTo(this,prev._autoGen?null:prev);}
viewarea._needNavigationMatrixUpdate=true;x3dom.nodeTypes.X3DBindableNode.prototype.activate.call(this,prev);},deactivate:function(prev){x3dom.nodeTypes.X3DBindableNode.prototype.deactivate.call(this,prev);},getTransformation:function(){return this.getCurrentTransform();},getCenterOfRotation:function(){return new x3dom.fields.SFVec3f(0,0,0);},setCenterOfRotation:function(cor){this._vf.centerOfRotation.setValues(cor);},getFieldOfView:function(){return 1.57079633;},setView:function(newView){var mat=this.getCurrentTransform();this._viewMatrix=newView.mult(mat);},setViewAbsolute:function(newView)
{this._viewMatrix=newView;},setProjectionMatrix:function(matrix)
{},resetView:function(){},getNear:function(){return 0.1;},getFar:function(){return 10000;},getImgPlaneHeightAtDistOne:function(){return 2.0;},getViewMatrix:function(){return null;},getProjectionMatrix:function(aspect){return null;},setZoom:function(value){}}));x3dom.registerNodeType("Viewpoint","Navigation",defineClass(x3dom.nodeTypes.X3DViewpointNode,function(ctx){x3dom.nodeTypes.Viewpoint.superClass.call(this,ctx);this.addField_SFFloat(ctx,'fieldOfView',0.785398);this.addField_SFVec3f(ctx,'position',0,0,10);this.addField_SFRotation(ctx,'orientation',0,0,0,1);this.addField_SFVec3f(ctx,'centerOfRotation',0,0,0);this.addField_SFFloat(ctx,'zNear',-1);this.addField_SFFloat(ctx,'zFar',-1);this._viewMatrix=x3dom.fields.SFMatrix4f.translation(this._vf.position).mult(this._vf.orientation.toMatrix()).inverse();this._projMatrix=null;this._lastAspect=1.0;this._zRatio=10000;this._zNear=this._vf.zNear;this._zFar=this._vf.zFar;this._imgPlaneHeightAtDistOne=2.0*Math.tan(this._vf.fieldOfView/2.0);},{fieldChanged:function(fieldName){if(fieldName=="position"||fieldName=="orientation"){this.resetView();}
else if(fieldName=="fieldOfView"||fieldName=="zNear"||fieldName=="zFar"){this._projMatrix=null;this._zNear=this._vf.zNear;this._zFar=this._vf.zFar;this._imgPlaneHeightAtDistOne=2.0*Math.tan(this._vf.fieldOfView/2.0);}
else if(fieldName.indexOf("bind")>=0){this.bind(this._vf.bind);}},setProjectionMatrix:function(matrix)
{this._projMatrix=matrix;},getCenterOfRotation:function(){return this.getCurrentTransform().multMatrixPnt(this._vf.centerOfRotation);},getViewMatrix:function(){return this._viewMatrix;},getFieldOfView:function(){return this._vf.fieldOfView;},resetView:function(){this._viewMatrix=x3dom.fields.SFMatrix4f.translation(this._vf.position).mult(this._vf.orientation.toMatrix()).inverse();if(this._vf.isActive&&this._nameSpace&&this._nameSpace.doc._viewarea){this._nameSpace.doc._viewarea.resetNavHelpers();}},getNear:function(){return this._zNear;},getFar:function(){return this._zFar;},getImgPlaneHeightAtDistOne:function(){return this._imgPlaneHeightAtDistOne;},getProjectionMatrix:function(aspect)
{var fovy=this._vf.fieldOfView;var zfar=this._vf.zFar;var znear=this._vf.zNear;if(znear<=0||zfar<=0)
{var nearScale=0.8,farScale=1.2;var viewarea=this._nameSpace.doc._viewarea;var scene=viewarea._scene;var min=x3dom.fields.SFVec3f.copy(scene._lastMin);var max=x3dom.fields.SFVec3f.copy(scene._lastMax);var dia=max.subtract(min);var sRad=dia.length()/2;var mat=viewarea.getViewMatrix().inverse();var vp=mat.e3();var translation=new x3dom.fields.SFVec3f(0,0,0),scaleFactor=new x3dom.fields.SFVec3f(1,1,1);var rotation=new x3dom.fields.Quaternion(0,0,1,0),scaleOrientation=new x3dom.fields.Quaternion(0,0,1,0);mat.getTransform(translation,rotation,scaleFactor,scaleOrientation);var minScal=scaleFactor.x,maxScal=scaleFactor.x;if(maxScal<scaleFactor.y)maxScal=scaleFactor.y;if(minScal>scaleFactor.y)minScal=scaleFactor.y;if(maxScal<scaleFactor.z)maxScal=scaleFactor.z;if(minScal>scaleFactor.z)minScal=scaleFactor.z;if(maxScal>1)
nearScale/=maxScal;else if(minScal>x3dom.fields.Eps&&minScal<1)
farScale/=minScal;var sCenter=min.add(dia.multiply(0.5));var vDist=(vp.subtract(sCenter)).length();if(sRad){if(vDist>sRad)
znear=(vDist-sRad)*nearScale;else
znear=0;zfar=(vDist+sRad)*farScale;}
else{znear=0.1;zfar=100000;}
var zNearLimit=zfar/this._zRatio;znear=Math.max(znear,Math.max(x3dom.fields.Eps,zNearLimit));if(zfar>this._vf.zNear&&this._vf.zNear>0)
znear=this._vf.zNear;if(this._vf.zFar>znear)
zfar=this._vf.zFar;if(zfar<=znear)
zfar=znear+1;}
if(this._projMatrix==null)
{this._projMatrix=x3dom.fields.SFMatrix4f.perspective(fovy,aspect,znear,zfar);}
else if(this._zNear!=znear||this._zFar!=zfar)
{var div=znear-zfar;this._projMatrix._22=(znear+zfar)/div;this._projMatrix._23=2*znear*zfar/div;}
else if(this._lastAspect!=aspect)
{this._projMatrix._00=(1/Math.tan(fovy/2))/aspect;this._lastAspect=aspect;}
this._zNear=znear;this._zFar=zfar;return this._projMatrix;}}));x3dom.registerNodeType("OrthoViewpoint","Navigation",defineClass(x3dom.nodeTypes.X3DViewpointNode,function(ctx){x3dom.nodeTypes.OrthoViewpoint.superClass.call(this,ctx);this.addField_MFFloat(ctx,'fieldOfView',[-1,-1,1,1]);this.addField_SFVec3f(ctx,'position',0,0,10);this.addField_SFRotation(ctx,'orientation',0,0,0,1);this.addField_SFVec3f(ctx,'centerOfRotation',0,0,0);this.addField_SFFloat(ctx,'zNear',-1);this.addField_SFFloat(ctx,'zFar',-1);this._viewMatrix=null;this._projMatrix=null;this._lastAspect=1.0;this._zRatio=10000;this._zNear=this._vf.zNear;this._zFar=this._vf.zFar;this._fieldOfView=this._vf.fieldOfView.slice(0);this.resetView();},{fieldChanged:function(fieldName){if(fieldName=="position"||fieldName=="orientation"){this.resetView();}
else if(fieldName=="fieldOfView")
{this._fieldOfView=this._vf.fieldOfView;this._projMatrix=null;}
else if(fieldName=="zNear"||fieldName=="zFar"){this._projMatrix=null;this.resetView();}
else if(fieldName.indexOf("bind")>=0){this.bind(this._vf.bind);}},getCenterOfRotation:function(){return this.getCurrentTransform().multMatrixPnt(this._vf.centerOfRotation);},getViewMatrix:function(){return this._viewMatrix;},resetView:function(){var offset=x3dom.fields.SFMatrix4f.translation(new x3dom.fields.SFVec3f((this._vf.fieldOfView[0]+this._vf.fieldOfView[2])/2,(this._vf.fieldOfView[1]+this._vf.fieldOfView[3])/2,0));this._viewMatrix=x3dom.fields.SFMatrix4f.translation(this._vf.position).mult(this._vf.orientation.toMatrix());this._viewMatrix=this._viewMatrix.inverse();this._projMatrix=null;if(this._vf.isActive&&this._nameSpace&&this._nameSpace.doc._viewarea){this._nameSpace.doc._viewarea.resetNavHelpers();}},getNear:function(){return this._vf.zNear;},getFar:function(){return this._vf.zFar;},getFieldOfView:function(){return 0.785;},setZoom:function(value){this._fieldOfView[0]=-value;this._fieldOfView[1]=-value;this._fieldOfView[2]=value;this._fieldOfView[3]=value;this._projMatrix=null;},getZoom:function(value){return this._fieldOfView;},getProjectionMatrix:function(aspect)
{var fov=this.getFieldOfView();var zfar=this._vf.zFar;var znear=this._vf.zNear;if(znear<=0||zfar<=0)
{var scene=this._nameSpace.doc._viewarea._scene;var min=x3dom.fields.SFVec3f.copy(scene._lastMin);var max=x3dom.fields.SFVec3f.copy(scene._lastMax);var dia=max.subtract(min);var tanfov2=Math.tan(fov/2.0);var dist1=(dia.y/2.0)/tanfov2+dia.z;var dist2=(dia.x/2.0)/tanfov2+dia.z;znear=0.00001;zfar=(dist1>dist2)?dist1*4:dist2*4;}
if(this._projMatrix==null||this._lastAspect!=aspect||this._zNear!=znear||this._zFar!=zfar)
{var near=this._zNear=znear;var far=this._zFar=zfar;var left=this._fieldOfView[0];var bottom=this._fieldOfView[1];var right=this._fieldOfView[2];var top=this._fieldOfView[3];this._projMatrix=x3dom.fields.SFMatrix4f.ortho(left,right,bottom,top,near,far,aspect);}
this._lastAspect=aspect;return this._projMatrix;}}));x3dom.registerNodeType("Viewfrustum","Navigation",defineClass(x3dom.nodeTypes.X3DViewpointNode,function(ctx){x3dom.nodeTypes.Viewfrustum.superClass.call(this,ctx);this.addField_SFMatrix4f(ctx,'modelview',1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1);this.addField_SFMatrix4f(ctx,'projection',1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1);this._viewMatrix=this._vf.modelview.transpose().inverse();this._projMatrix=this._vf.projection.transpose();this._centerOfRotation=new x3dom.fields.SFVec3f(0,0,0);},{fieldChanged:function(fieldName){if(fieldName=="modelview"){this.resetView();}
else if(fieldName=="projection"){this._projMatrix=this._vf.projection.transpose();}
else if(fieldName.indexOf("bind")>=0){this.bind(this._vf.bind);}},getCenterOfRotation:function(){return this.getCurrentTransform().multMatrixPnt(this._centerOfRotation);},setCenterOfRotation:function(cor){this._centerOfRotation.setValues(cor);},getViewMatrix:function(){return this._viewMatrix;},getFieldOfView:function(){return(2.0*Math.atan(1.0/this._projMatrix._11));},getImgPlaneHeightAtDistOne:function(){return 2.0/this._projMatrix._11;},resetView:function(){this._viewMatrix=this._vf.modelview.transpose().inverse();this._centerOfRotation=new x3dom.fields.SFVec3f(0,0,0);},getProjectionMatrix:function(aspect){return this._projMatrix;}}));x3dom.registerNodeType("X3DNavigationInfoNode","Navigation",defineClass(x3dom.nodeTypes.X3DBindableNode,function(ctx){x3dom.nodeTypes.X3DNavigationInfoNode.superClass.call(this,ctx);}));x3dom.registerNodeType("NavigationInfo","Navigation",defineClass(x3dom.nodeTypes.X3DNavigationInfoNode,function(ctx){x3dom.nodeTypes.NavigationInfo.superClass.call(this,ctx);this.addField_SFBool(ctx,'headlight',true);this.addField_MFString(ctx,'type',["EXAMINE","ANY"]);this.addField_MFFloat(ctx,'typeParams',[-0.4,60,0.05,2.8]);this.addField_SFString(ctx,'explorationMode','all');this.addField_MFFloat(ctx,'avatarSize',[0.25,1.6,0.75]);this.addField_SFFloat(ctx,'visibilityLimit',0.0);this.addField_SFFloat(ctx,'speed',1.0);this.addField_SFTime(ctx,'transitionTime',1.0);this.addField_MFString(ctx,'transitionType',["LINEAR"]);this._validTypes=["none","examine","turntable","fly","freefly","lookat","lookaround","walk","game","helicopter","any"];this._typeMapping={"default":x3dom.DefaultNavigation,"turntable":x3dom.TurntableNavigation};this._heliUpdated=false;var type=this.setType(this.getType());x3dom.debug.logInfo("NavType: "+type);},{fieldChanged:function(fieldName){if(fieldName=="typeParams"){this._heliUpdated=false;}
else if(fieldName=="type"){this.setType(this.getType());}},setType:function(type,viewarea){var navType=this.checkType(type.toLowerCase());var oldType=this.checkType(this.getType());if(oldType!==navType||this._impl==null){if(this._typeMapping[navType]==null)
this._impl=new this._typeMapping['default'](this);else
this._impl=new this._typeMapping[navType](this);switch(navType){case'game':if(viewarea)
viewarea.initMouseState();else
this._nameSpace.doc._viewarea.initMouseState();break;case'helicopter':this._heliUpdated=false;break;case"turntable":if(viewarea){viewarea.initMouseState();}
else if(this._nameSpace.doc._viewarea){this._nameSpace.doc._viewarea.initMouseState();}
break;default:break;}
if(this._nameSpace.doc._viewarea)
this._impl.init(this._nameSpace.doc._viewarea,false);}
this._vf.type[0]=navType;x3dom.debug.logInfo("Switch to "+navType+" mode.");},getType:function(){var type=this._vf.type[0].toLowerCase();if(type.length<=1)
type="none";else if(type=="any")
type="examine";return type;},getTypeParams:function(){var length=this._vf.typeParams.length;var theta=(length>=1)?this._vf.typeParams[0]:-0.4;var height=(length>=2)?this._vf.typeParams[1]:60.0;var minAngle=(length>=3)?this._vf.typeParams[2]:x3dom.fields.Eps;var maxAngle=(length>=4)?this._vf.typeParams[3]:Math.PI-x3dom.fields.Eps;var params=[theta,height,minAngle,maxAngle];if(length>=5)
{params=params.concat(this._vf.typeParams.slice(4));}
console.log(params);return params;},setTypeParams:function(params){for(var i=0;i<params.length;i++){this._vf.typeParams[i]=params[i];}},checkType:function(type){if(this._validTypes.indexOf(type)>-1){return type;}
else{x3dom.debug.logWarning(type+" is no valid navigation type, use one of "+
this._validTypes.toString());return"examine";}},getExplorationMode:function(){switch(this._vf.explorationMode.toLowerCase()){case"all":return 7;case"rotate":return 1;case"zoom":return 2;case"pan":return 4;case"none":return 0;default:return 7;}}}));x3dom.registerNodeType("Billboard","Navigation",defineClass(x3dom.nodeTypes.X3DGroupingNode,function(ctx){x3dom.nodeTypes.Billboard.superClass.call(this,ctx);this.addField_SFVec3f(ctx,'axisOfRotation',0,1,0);this._eye=new x3dom.fields.SFVec3f(0,0,0);this._eyeViewUp=new x3dom.fields.SFVec3f(0,0,0);this._eyeLook=new x3dom.fields.SFVec3f(0,0,0);},{collectDrawableObjects:function(transform,drawableCollection,singlePath,invalidateCache,planeMask,clipPlanes)
{if(singlePath&&(this._parentNodes.length>1))
singlePath=false;if(singlePath&&(invalidateCache=invalidateCache||this.cacheInvalid()))
this.invalidateCache();planeMask=drawableCollection.cull(transform,this.graphState(),singlePath,planeMask);if(planeMask<0){return;}
singlePath=false;var vol=this.getVolume();var min=x3dom.fields.SFVec3f.MAX();var max=x3dom.fields.SFVec3f.MIN();vol.getBounds(min,max);var mat_view=drawableCollection.viewMatrix;var center=new x3dom.fields.SFVec3f(0,0,0);center=mat_view.inverse().multMatrixPnt(center);var mat_view_model=mat_view.mult(transform);this._eye=transform.inverse().multMatrixPnt(center);this._eyeViewUp=new x3dom.fields.SFVec3f(mat_view_model._10,mat_view_model._11,mat_view_model._12);this._eyeLook=new x3dom.fields.SFVec3f(mat_view_model._20,mat_view_model._21,mat_view_model._22);var rotMat=x3dom.fields.SFMatrix4f.identity();var mid=max.add(min).multiply(0.5);var billboard_to_viewer=this._eye.subtract(mid);if(this._vf.axisOfRotation.equals(new x3dom.fields.SFVec3f(0,0,0),x3dom.fields.Eps)){var rot1=x3dom.fields.Quaternion.rotateFromTo(billboard_to_viewer,new x3dom.fields.SFVec3f(0,0,1));rotMat=rot1.toMatrix().transpose();var yAxis=rotMat.multMatrixPnt(new x3dom.fields.SFVec3f(0,1,0)).normalize();var zAxis=rotMat.multMatrixPnt(new x3dom.fields.SFVec3f(0,0,1)).normalize();if(!this._eyeViewUp.equals(new x3dom.fields.SFVec3f(0,0,0),x3dom.fields.Eps)){var rot2=x3dom.fields.Quaternion.rotateFromTo(this._eyeLook,zAxis);var rotatedyAxis=rot2.toMatrix().transpose().multMatrixVec(yAxis);var rot3=x3dom.fields.Quaternion.rotateFromTo(this._eyeViewUp,rotatedyAxis);rotMat=rot2.toMatrix().transpose().mult(rotMat);rotMat=rot3.toMatrix().transpose().mult(rotMat);}}
else{var normalPlane=this._vf.axisOfRotation.cross(billboard_to_viewer).normalize();if(this._eye.z<0){normalPlane=normalPlane.multiply(-1);}
var degreesToRotate=Math.asin(normalPlane.dot(new x3dom.fields.SFVec3f(0,0,1)));if(this._eye.z<0){degreesToRotate+=Math.PI;}
rotMat=x3dom.fields.SFMatrix4f.parseRotation(this._vf.axisOfRotation.x+", "+this._vf.axisOfRotation.y+", "+
this._vf.axisOfRotation.z+", "+degreesToRotate*(-1));}
var childTransform=this.transformMatrix(transform.mult(rotMat));for(var i=0,i_n=this._childNodes.length;i<i_n;i++)
{var cnode=this._childNodes[i];if(cnode){cnode.collectDrawableObjects(childTransform,drawableCollection,singlePath,invalidateCache,planeMask,clipPlanes);}}}}));x3dom.registerNodeType("Collision","Navigation",defineClass(x3dom.nodeTypes.X3DGroupingNode,function(ctx){x3dom.nodeTypes.Collision.superClass.call(this,ctx);this.addField_SFBool(ctx,"enabled",true);this.addField_SFNode("proxy",x3dom.nodeTypes.X3DGroupingNode);this.addField_SFTime(ctx,"collideTime",0);this.addField_SFBool(ctx,"isActive",true);},{collectDrawableObjects:function(transform,drawableCollection,singlePath,invalidateCache,planeMask,clipPlanes)
{if(singlePath&&(this._parentNodes.length>1))
singlePath=false;if(singlePath&&(invalidateCache=invalidateCache||this.cacheInvalid()))
this.invalidateCache();planeMask=drawableCollection.cull(transform,this.graphState(),singlePath,planeMask);if(planeMask<0){return;}
var cnode,childTransform;if(singlePath){if(!this._graph.globalMatrix){this._graph.globalMatrix=this.transformMatrix(transform);}
childTransform=this._graph.globalMatrix;}
else{childTransform=this.transformMatrix(transform);}
for(var i=0,n=this._childNodes.length;i<n;i++)
{if((cnode=this._childNodes[i])&&(cnode!==this._cf.proxy.node)){cnode.collectDrawableObjects(childTransform,drawableCollection,singlePath,invalidateCache,planeMask,clipPlanes);}}}}));x3dom.registerNodeType("X3DLODNode","Navigation",defineClass(x3dom.nodeTypes.X3DGroupingNode,function(ctx){x3dom.nodeTypes.X3DLODNode.superClass.call(this,ctx);this.addField_SFBool(ctx,"forceTransitions",false);this.addField_SFVec3f(ctx,"center",0,0,0);this._eye=new x3dom.fields.SFVec3f(0,0,0);},{collectDrawableObjects:function(transform,drawableCollection,singlePath,invalidateCache,planeMask,clipPlanes)
{if(singlePath&&(this._parentNodes.length>1))
singlePath=false;if(singlePath&&(invalidateCache=invalidateCache||this.cacheInvalid()))
this.invalidateCache();planeMask=drawableCollection.cull(transform,this.graphState(),singlePath,planeMask);if(planeMask<0){return;}
singlePath=false;this.visitChildren(transform,drawableCollection,singlePath,invalidateCache,planeMask,clipPlanes);},visitChildren:function(transform,drawableCollection,singlePath,invalidateCache,planeMask,clipPlanes){}}));x3dom.registerNodeType("LOD","Navigation",defineClass(x3dom.nodeTypes.X3DLODNode,function(ctx){x3dom.nodeTypes.LOD.superClass.call(this,ctx);this.addField_MFFloat(ctx,"range",[]);this._lastRangePos=-1;},{visitChildren:function(transform,drawableCollection,singlePath,invalidateCache,planeMask,clipPlanes)
{var i=0,n=this._childNodes.length;var mat_view=drawableCollection.viewMatrix;var center=new x3dom.fields.SFVec3f(0,0,0);center=mat_view.inverse().multMatrixPnt(center);this._eye=transform.inverse().multMatrixPnt(center);var len=this._vf.center.subtract(this._eye).length();while(i<this._vf.range.length&&len>this._vf.range[i]){i++;}
if(i&&i>=n){i=n-1;}
this._lastRangePos=i;var cnode=this._childNodes[i];if(n&&cnode)
{var childTransform=this.transformMatrix(transform);cnode.collectDrawableObjects(childTransform,drawableCollection,singlePath,invalidateCache,planeMask,clipPlanes);}},getVolume:function()
{var vol=this._graph.volume;if(!this.volumeValid()&&this._vf.render)
{var child,childVol;if(this._lastRangePos>=0){child=this._childNodes[this._lastRangePos];childVol=(child&&child._vf.render===true)?child.getVolume():null;if(childVol&&childVol.isValid())
vol.extendBounds(childVol.min,childVol.max);}
else{for(var i=0,n=this._childNodes.length;i<n;i++)
{if(!(child=this._childNodes[i])||child._vf.render!==true)
continue;childVol=child.getVolume();if(childVol&&childVol.isValid())
vol.extendBounds(childVol.min,childVol.max);}}}
return vol;},nodeChanged:function(){this.invalidateVolume();},fieldChanged:function(fieldName){if(fieldName=="render"||fieldName=="center"||fieldName=="range"){this.invalidateVolume();}}}));x3dom.registerNodeType("DynamicLOD","Navigation",defineClass(x3dom.nodeTypes.X3DLODNode,function(ctx){x3dom.nodeTypes.DynamicLOD.superClass.call(this,ctx);this.addField_SFFloat(ctx,'subScale',0.5);this.addField_SFVec2f(ctx,'size',2,2);this.addField_SFVec2f(ctx,'subdivision',1,1);this.addField_SFNode('root',x3dom.nodeTypes.X3DShapeNode);this.addField_SFString(ctx,'urlHead',"http://r");this.addField_SFString(ctx,'urlCenter',".ortho.tiles.virtualearth.net/tiles/h");this.addField_SFString(ctx,'urlTail',".png?g=-1");this.rootGeometry=new x3dom.nodeTypes.Plane(ctx);this.level=0;this.quadrant=4;this.cell="";},{nodeChanged:function()
{var root=this._cf.root.node;if(root==null||root._cf.geometry.node!=null)
return;this.rootGeometry._vf.size.setValues(this._vf.size);this.rootGeometry._vf.subdivision.setValues(this._vf.subdivision);this.rootGeometry._vf.center.setValues(this._vf.center);this.rootGeometry.fieldChanged("subdivision");this._cf.root.node.addChild(this.rootGeometry);this.rootGeometry.nodeChanged();this._cf.root.node.nodeChanged();this._nameSpace.doc.needRender=true;},visitChildren:function(transform,drawableCollection,singlePath,invalidateCache,planeMask,clipPlanes)
{var root=this._cf.root.node;if(root==null)
return;var mat_view=drawableCollection.viewMatrix;var center=new x3dom.fields.SFVec3f(0,0,0);center=mat_view.inverse().multMatrixPnt(center);this._eye=transform.inverse().multMatrixPnt(center);var l,len=this._vf.center.subtract(this._eye).length();if(len>x3dom.fields.Eps&&len*this._vf.subScale<=this._vf.size.length()){if(this._childNodes.length<=1){var offset=new Array(new x3dom.fields.SFVec3f(-0.25*this._vf.size.x,0.25*this._vf.size.y,0),new x3dom.fields.SFVec3f(0.25*this._vf.size.x,0.25*this._vf.size.y,0),new x3dom.fields.SFVec3f(-0.25*this._vf.size.x,-0.25*this._vf.size.y,0),new x3dom.fields.SFVec3f(0.25*this._vf.size.x,-0.25*this._vf.size.y,0));for(l=0;l<4;l++){var node=new x3dom.nodeTypes.DynamicLOD();node._nameSpace=this._nameSpace;node._eye.setValues(this._eye);node.level=this.level+1;node.quadrant=l;node.cell=this.cell+l;node._vf.urlHead=this._vf.urlHead;node._vf.urlCenter=this._vf.urlCenter;node._vf.urlTail=this._vf.urlTail;node._vf.center=this._vf.center.add(offset[l]);node._vf.size=this._vf.size.multiply(0.5);node._vf.subdivision.setValues(this._vf.subdivision);var app=new x3dom.nodeTypes.Appearance();var tex=new x3dom.nodeTypes.ImageTexture();tex._nameSpace=this._nameSpace;tex._vf.url[0]=this._vf.urlHead+node.quadrant+this._vf.urlCenter+node.cell+this._vf.urlTail;app.addChild(tex);tex.nodeChanged();var shape=new x3dom.nodeTypes.Shape();shape._nameSpace=this._nameSpace;shape.addChild(app);app.nodeChanged();node.addChild(shape,"root");shape.nodeChanged();this.addChild(node);node.nodeChanged();}}
else{for(l=1;l<this._childNodes.length;l++){this._childNodes[l].collectDrawableObjects(transform,drawableCollection,singlePath,invalidateCache,planeMask,clipPlanes);}}}
else{root.collectDrawableObjects(transform,drawableCollection,singlePath,invalidateCache,planeMask,clipPlanes);}},getVolume:function(){var vol=this._graph.volume;if(!vol.isValid()){vol.min.setValues(this._vf.center);vol.min.x-=0.5*this._vf.size.x;vol.min.y-=0.5*this._vf.size.y;vol.min.z-=x3dom.fields.Eps;vol.max.setValues(this._vf.center);vol.max.x+=0.5*this._vf.size.x;vol.max.y+=0.5*this._vf.size.y;vol.max.z+=x3dom.fields.Eps;}
return vol;}}));x3dom.DefaultNavigation=function(navigationNode)
{this.navi=navigationNode;};x3dom.DefaultNavigation.prototype.onMousePress=function(view,x,y,buttonState)
{};x3dom.DefaultNavigation.prototype.onMouseReleased=function(view,x,y,buttonState,prevButton)
{};x3dom.DefaultNavigation.prototype.init=function(view,flyTo)
{};x3dom.DefaultNavigation.prototype.zoom=function(view,zoomAmount)
{var navi=this.navi;var viewpoint=view._scene.getViewpoint();var d=(view._scene._lastMax.subtract(view._scene._lastMin)).length();d=((d<x3dom.fields.Eps)?1:d)*navi._vf.speed;vec=new x3dom.fields.SFVec3f(0,0,d*(zoomAmount)/view._height);if(x3dom.isa(viewpoint,x3dom.nodeTypes.OrthoViewpoint))
{viewpoint.setZoom(Math.abs(viewpoint._fieldOfView[0])-vec.z);}
else
{if(navi._vf.typeParams.length>=6){var min=-navi._vf.typeParams[5];var max=navi._vf.typeParams[4];view._movement.z=Math.min(Math.max(view._movement.z,min),max);}
view._movement=view._movement.add(vec);mat=view.getViewpointMatrix().mult(view._transMat);view._transMat=mat.inverse().mult(x3dom.fields.SFMatrix4f.translation(view._movement)).mult(mat);}}
x3dom.DefaultNavigation.prototype.moveForward=function(view)
{var navi=this.navi;if(navi.getType()==="game")
{var avatarRadius=0.25;var avatarHeight=1.6;if(navi._vf.avatarSize.length>2){avatarRadius=navi._vf.avatarSize[0];avatarHeight=navi._vf.avatarSize[1];}
var speed=5*view._deltaT*navi._vf.speed;var yRotRad=(view._yaw/180*Math.PI);var xRotRad=(view._pitch/180*Math.PI);var dist=0;var fMat=view._flyMat.inverse();view._scene._nameSpace.doc.ctx.pickValue(view,view._width/2,view._height/2,view._lastButton);if(view._pickingInfo.pickObj)
{dist=view._pickingInfo.pickPos.subtract(fMat.e3()).length();if(dist<=2*avatarRadius){}
else{view._eyePos.x-=Math.sin(yRotRad)*speed;view._eyePos.z+=Math.cos(yRotRad)*speed;view._eyePos.y+=Math.sin(xRotRad)*speed;}}}};x3dom.DefaultNavigation.prototype.moveBackwards=function(view)
{var navi=this.navi;if(navi.getType()==="game")
{var speed=5*view._deltaT*navi._vf.speed;var yRotRad=(view._yaw/180*Math.PI);var xRotRad=(view._pitch/180*Math.PI);view._eyePos.x+=Math.sin(yRotRad)*speed;view._eyePos.z-=Math.cos(yRotRad)*speed;view._eyePos.y-=Math.sin(xRotRad)*speed;}};x3dom.DefaultNavigation.prototype.strafeLeft=function(view)
{var navi=this.navi;if(navi.getType()==="game")
{var speed=5*view._deltaT*navi._vf.speed;var yRotRad=(view._yaw/180*Math.PI);view._eyePos.x+=Math.cos(yRotRad)*speed;view._eyePos.z+=Math.sin(yRotRad)*speed;}};x3dom.DefaultNavigation.prototype.strafeRight=function(view)
{var navi=this.navi;if(navi.getType()==="game")
{var speed=5*view._deltaT*navi._vf.speed;var yRotRad=(view._yaw/180*Math.PI);view._eyePos.x-=Math.cos(yRotRad)*speed;view._eyePos.z-=Math.sin(yRotRad)*speed;}};x3dom.DefaultNavigation.prototype.navigateTo=function(view,timeStamp)
{var navi=this.navi;var navType=navi.getType();var savedPickingInfo=null;var needNavAnim=(view._currentInputType==x3dom.InputTypes.NAVIGATION)&&(navType==="game"||(view._lastButton>0&&(navType.indexOf("fly")>=0||navType==="walk"||navType==="helicopter"||navType.substr(0,5)==="looka")));view._deltaT=timeStamp-view._lastTS;var removeZeroMargin=function(val,offset){if(val>0){if(val<=offset){return 0;}else{return val-offset;}}else if(val<=0){if(val>=-offset){return 0;}else{return val+offset;}}};var humanizeDiff=function(scale,diff){return((diff<0)?-1:1)*Math.pow(scale*Math.abs(diff),1.65);};if(needNavAnim)
{if(view._pickingInfo.pickObj!==null){savedPickingInfo={pickPos:view._pickingInfo.pickPos,pickNorm:view._pickingInfo.pickNorm,pickObj:view._pickingInfo.pickObj,firstObj:view._pickingInfo.firstObj,lastObj:view._pickingInfo.lastObj,lastClickObj:view._pickingInfo.lastClickObj,shadowObjectId:view._pickingInfo.shadowObjectId};}
var avatarRadius=0.25;var avatarHeight=1.6;var avatarKnee=0.75;if(navi._vf.avatarSize.length>2){avatarRadius=navi._vf.avatarSize[0];avatarHeight=navi._vf.avatarSize[1];avatarKnee=navi._vf.avatarSize[2];}
var currViewMat=view.getViewMatrix();var dist=0;var screenSize=Math.min(view._width,view._height);var rdeltaX=removeZeroMargin((view._pressX-view._lastX)/screenSize,0.01);var rdeltaY=removeZeroMargin((view._pressY-view._lastY)/screenSize,0.01);var userXdiff=humanizeDiff(1,rdeltaX);var userYdiff=humanizeDiff(1,rdeltaY);var step=(view._lastButton&2)?-1:1;step*=(view._deltaT*navi._vf.speed);var userXstep=view._deltaT*navi._vf.speed*userXdiff;var userYstep=view._deltaT*navi._vf.speed*userYdiff;var phi=Math.PI*view._deltaT*userXdiff;var theta=Math.PI*view._deltaT*userYdiff;if(view._needNavigationMatrixUpdate===true)
{view._needNavigationMatrixUpdate=false;view._rotMat=x3dom.fields.SFMatrix4f.identity();view._transMat=x3dom.fields.SFMatrix4f.identity();view._movement=new x3dom.fields.SFVec3f(0,0,0);var angleX=0;var angleY=Math.asin(currViewMat._02);var C=Math.cos(angleY);if(Math.abs(C)>0.0001){angleX=Math.atan2(-currViewMat._12/C,currViewMat._22/C);}
view._flyMat=currViewMat.inverse();view._from=view._flyMat.e3();view._at=view._from.subtract(view._flyMat.e2());if(navType==="helicopter")
view._at.y=view._from.y;view._up=view._flyMat.e1();view._pitch=angleX*180/Math.PI;view._yaw=angleY*180/Math.PI;view._eyePos=view._from.negate();}
var tmpAt=null,tmpUp=null,tmpMat=null;var q,temp,fin;var lv,sv,up;if(navType==="game")
{view._pitch+=view._dy;view._yaw+=view._dx;if(view._pitch>=89)view._pitch=89;if(view._pitch<=-89)view._pitch=-89;if(view._yaw>=360)view._yaw-=360;if(view._yaw<0)view._yaw=360+view._yaw;view._dx=0;view._dy=0;var xMat=x3dom.fields.SFMatrix4f.rotationX(view._pitch/180*Math.PI);var yMat=x3dom.fields.SFMatrix4f.rotationY(view._yaw/180*Math.PI);var fPos=x3dom.fields.SFMatrix4f.translation(view._eyePos);view._flyMat=xMat.mult(yMat).mult(fPos);var flyMat=view._flyMat.inverse();var tmpFrom=flyMat.e3();tmpUp=new x3dom.fields.SFVec3f(0,-1,0);tmpAt=tmpFrom.add(tmpUp);tmpUp=flyMat.e0().cross(tmpUp).normalize();tmpMat=x3dom.fields.SFMatrix4f.lookAt(tmpFrom,tmpAt,tmpUp);tmpMat=tmpMat.inverse();view._scene._nameSpace.doc.ctx.pickValue(view,view._width/2,view._height/2,view._lastButton,tmpMat,view.getProjectionMatrix().mult(tmpMat));if(view._pickingInfo.pickObj)
{dist=view._pickingInfo.pickPos.subtract(tmpFrom).length();tmpFrom.y+=(avatarHeight-dist);flyMat.setTranslate(tmpFrom);view._eyePos=flyMat.e3().negate();view._flyMat=flyMat.inverse();view._pickingInfo.pickObj=null;}
view._scene.getViewpoint().setView(view._flyMat);return needNavAnim;}
else if(navType==="helicopter"){var typeParams=navi.getTypeParams();if(view._lastButton&2)
{var stepUp=200*userYstep;typeParams[1]+=stepUp;navi.setTypeParams(typeParams);}
if(view._lastButton&1){step=300*userYstep;}
else{step=0;}
theta=typeParams[0];view._from.y=typeParams[1];view._at.y=view._from.y;q=x3dom.fields.Quaternion.axisAngle(view._up,phi);temp=q.toMatrix();fin=x3dom.fields.SFMatrix4f.translation(view._from);fin=fin.mult(temp);temp=x3dom.fields.SFMatrix4f.translation(view._from.negate());fin=fin.mult(temp);view._at=fin.multMatrixPnt(view._at);lv=view._at.subtract(view._from).normalize();sv=lv.cross(view._up).normalize();up=sv.cross(lv).normalize();lv=lv.multiply(step);view._from=view._from.add(lv);view._at=view._at.add(lv);q=x3dom.fields.Quaternion.axisAngle(sv,theta);temp=q.toMatrix();fin=x3dom.fields.SFMatrix4f.translation(view._from);fin=fin.mult(temp);temp=x3dom.fields.SFMatrix4f.translation(view._from.negate());fin=fin.mult(temp);var at=fin.multMatrixPnt(view._at);view._flyMat=x3dom.fields.SFMatrix4f.lookAt(view._from,at,up);view._scene.getViewpoint().setView(view._flyMat.inverse());return needNavAnim;}
q=x3dom.fields.Quaternion.axisAngle(view._up,phi);temp=q.toMatrix();fin=x3dom.fields.SFMatrix4f.translation(view._from);fin=fin.mult(temp);temp=x3dom.fields.SFMatrix4f.translation(view._from.negate());fin=fin.mult(temp);view._at=fin.multMatrixPnt(view._at);lv=view._at.subtract(view._from).normalize();sv=lv.cross(view._up).normalize();up=sv.cross(lv).normalize();q=x3dom.fields.Quaternion.axisAngle(sv,theta);temp=q.toMatrix();fin=x3dom.fields.SFMatrix4f.translation(view._from);fin=fin.mult(temp);temp=x3dom.fields.SFMatrix4f.translation(view._from.negate());fin=fin.mult(temp);view._at=fin.multMatrixPnt(view._at);if(navType.substr(0,5)!=="looka")
{var currProjMat=view.getProjectionMatrix();if(navType!=="freefly"){if(step<0){tmpMat=new x3dom.fields.SFMatrix4f();tmpMat.setValue(view._last_mat_view.e0(),view._last_mat_view.e1(),view._last_mat_view.e2().negate(),view._last_mat_view.e3());view._scene._nameSpace.doc.ctx.pickValue(view,view._width/2,view._height/2,view._lastButton,tmpMat,currProjMat.mult(tmpMat));}
else{view._scene._nameSpace.doc.ctx.pickValue(view,view._width/2,view._height/2,view._lastButton);}
if(view._pickingInfo.pickObj)
{dist=view._pickingInfo.pickPos.subtract(view._from).length();if(dist<=avatarRadius){step=0;}}}
lv=view._at.subtract(view._from).normalize().multiply(step);view._at=view._at.add(lv);view._from=view._from.add(lv);if(navType==="walk")
{tmpAt=view._from.addScaled(up,-1.0);tmpUp=sv.cross(up.negate()).normalize();tmpMat=x3dom.fields.SFMatrix4f.lookAt(view._from,tmpAt,tmpUp);tmpMat=tmpMat.inverse();view._scene._nameSpace.doc.ctx.pickValue(view,view._width/2,view._height/2,view._lastButton,tmpMat,currProjMat.mult(tmpMat));if(view._pickingInfo.pickObj)
{dist=view._pickingInfo.pickPos.subtract(view._from).length();view._at=view._at.add(up.multiply(avatarHeight-dist));view._from=view._from.add(up.multiply(avatarHeight-dist));}}
view._pickingInfo.pickObj=null;}
view._flyMat=x3dom.fields.SFMatrix4f.lookAt(view._from,view._at,up);view._scene.getViewpoint().setView(view._flyMat.inverse());if(savedPickingInfo!==null){view._pickingInfo=savedPickingInfo;}}
return needNavAnim;};x3dom.DefaultNavigation.prototype.animateTo=function(view,target,prev,dur)
{var navi=this.navi;if(x3dom.isa(target,x3dom.nodeTypes.X3DViewpointNode)){target=target.getViewMatrix().mult(target.getCurrentTransform().inverse());}
if(navi._vf.transitionType[0].toLowerCase()!=="teleport"&&dur!=0&&navi.getType()!=="game")
{if(prev&&x3dom.isa(prev,x3dom.nodeTypes.X3DViewpointNode)){prev=prev.getViewMatrix().mult(prev.getCurrentTransform().inverse()).mult(view._transMat).mult(view._rotMat);view._mixer.beginTime=view._lastTS;if(arguments.length>=4&&arguments[3]!=null){view._mixer.endTime=view._lastTS+dur;}
else{view._mixer.endTime=view._lastTS+navi._vf.transitionTime;}
view._mixer.setBeginMatrix(prev);view._mixer.setEndMatrix(target);view._scene.getViewpoint().setView(prev);}
else{view._scene.getViewpoint().setView(target);}}
else
{view._scene.getViewpoint().setView(target);}
view._rotMat=x3dom.fields.SFMatrix4f.identity();view._transMat=x3dom.fields.SFMatrix4f.identity();view._movement=new x3dom.fields.SFVec3f(0,0,0);view._needNavigationMatrixUpdate=true;};x3dom.DefaultNavigation.prototype.orthoAnimateTo=function(view,target,prev,duration)
{var navi=this.navi;duration=duration||navi._vf.transitionTime;view._interpolator.beginValue=prev;view._interpolator.endValue=target;view._interpolator.beginTime=view._lastTS;view._interpolator.endTime=view._lastTS+duration;};x3dom.DefaultNavigation.prototype.resetView=function(view)
{var navi=this.navi;if(navi._vf.transitionType[0].toLowerCase()!=="teleport"&&navi.getType()!=="game")
{var viewpoint=view._scene.getViewpoint();view._mixer.beginTime=view._lastTS;view._mixer.endTime=view._lastTS+navi._vf.transitionTime;view._mixer.setBeginMatrix(view.getViewMatrix());if(x3dom.isa(viewpoint,x3dom.nodeTypes.OrthoViewpoint))
{this.orthoAnimateTo(view,Math.abs(viewpoint._vf.fieldOfView[0]),Math.abs(viewpoint._fieldOfView[0]));}
var target=view._scene.getViewpoint();target.resetView();target=target.getViewMatrix().mult(target.getCurrentTransform().inverse());view._mixer.setEndMatrix(target);}else
{view._scene.getViewpoint().resetView();}
view.resetNavHelpers();navi._heliUpdated=false;};x3dom.DefaultNavigation.prototype.onDrag=function(view,x,y,buttonState)
{var navi=this.navi;var navType=navi.getType();var navRestrict=navi.getExplorationMode();if(navType==="none"||navRestrict==0){return;}
var viewpoint=view._scene.getViewpoint();var dx=x-view._lastX;var dy=y-view._lastY;var d,vec,cor,mat=null;var alpha,beta;buttonState=(!navRestrict||(navRestrict!=7&&buttonState==1))?navRestrict:buttonState;if(buttonState&1)
{alpha=(dy*2*Math.PI)/view._width;beta=(dx*2*Math.PI)/view._height;mat=view.getViewMatrix();var mx=x3dom.fields.SFMatrix4f.rotationX(alpha);var my=x3dom.fields.SFMatrix4f.rotationY(beta);var center=viewpoint.getCenterOfRotation();mat.setTranslate(new x3dom.fields.SFVec3f(0,0,0));view._rotMat=view._rotMat.mult(x3dom.fields.SFMatrix4f.translation(center)).mult(mat.inverse()).mult(mx).mult(my).mult(mat).mult(x3dom.fields.SFMatrix4f.translation(center.negate()));}
if(buttonState&4)
{d=(view._scene._lastMax.subtract(view._scene._lastMin)).length();d=((d<x3dom.fields.Eps)?1:d)*navi._vf.speed;vec=new x3dom.fields.SFVec3f(d*dx/view._width,d*(-dy)/view._height,0);view._movement=view._movement.add(vec);mat=view.getViewpointMatrix().mult(view._transMat);view._transMat=mat.inverse().mult(x3dom.fields.SFMatrix4f.translation(view._movement)).mult(mat);}
if(buttonState&2)
{d=(view._scene._lastMax.subtract(view._scene._lastMin)).length();d=((d<x3dom.fields.Eps)?1:d)*navi._vf.speed;vec=new x3dom.fields.SFVec3f(0,0,d*(dx+dy)/view._height);if(x3dom.isa(viewpoint,x3dom.nodeTypes.OrthoViewpoint))
{viewpoint.setZoom(Math.abs(viewpoint._fieldOfView[0])-vec.z);}
else
{if(navi._vf.typeParams.length>=6){var min=-navi._vf.typeParams[5];var max=navi._vf.typeParams[4];view._movement.z=Math.min(Math.max(view._movement.z,min),max);}
view._movement=view._movement.add(vec);mat=view.getViewpointMatrix().mult(view._transMat);view._transMat=mat.inverse().mult(x3dom.fields.SFMatrix4f.translation(view._movement)).mult(mat);}}
view._isMoving=true;view._dx=dx;view._dy=dy;view._lastX=x;view._lastY=y;};x3dom.DefaultNavigation.prototype.onTouchStart=function(view,evt,touches)
{};x3dom.DefaultNavigation.prototype.onTouchDrag=function(view,evt,touches,translation,rotation)
{if(view._currentInputType==x3dom.InputTypes.NAVIGATION)
{var navi=this.navi;var viewpoint=view._scene.getViewpoint();if(navi.getType()==="examine")
{if(translation)
{var distance=(view._scene._lastMax.subtract(view._scene._lastMin)).length();distance=((distance<x3dom.fields.Eps)?1:distance)*navi._vf.speed;translation=translation.multiply(distance);view._movement=view._movement.add(translation);view._transMat=viewpoint.getViewMatrix().inverse().mult(x3dom.fields.SFMatrix4f.translation(view._movement)).mult(viewpoint.getViewMatrix());}
if(rotation)
{var center=viewpoint.getCenterOfRotation();var mat=view.getViewMatrix();mat.setTranslate(new x3dom.fields.SFVec3f(0,0,0));view._rotMat=view._rotMat.mult(x3dom.fields.SFMatrix4f.translation(center)).mult(mat.inverse()).mult(rotation).mult(mat).mult(x3dom.fields.SFMatrix4f.translation(center.negate()));}
view._isMoving=true;}}};x3dom.DefaultNavigation.prototype.onTouchEnd=function(evt,touches)
{};x3dom.DefaultNavigation.prototype.onDoubleClick=function(view,x,y)
{if(view._doc._x3dElem.hasAttribute('disableDoubleClick')&&view._doc._x3dElem.getAttribute('disableDoubleClick')==='true'){return;}
var navi=view._scene.getNavigationInfo();if(navi.getType()=="none"){return;}
var pickMode=view._scene._vf.pickMode.toLowerCase();if((pickMode=="color"||pickMode=="texcoord")){return;}
var viewpoint=view._scene.getViewpoint();viewpoint.setCenterOfRotation(view._pick);x3dom.debug.logInfo("New center of Rotation:  "+view._pick);var mat=view.getViewMatrix().inverse();var from=mat.e3();var at=view._pick;var up=mat.e1();var norm=mat.e0().cross(up).normalize();var dist=norm.dot(view._pick.subtract(from));from=at.addScaled(norm,-dist);mat=x3dom.fields.SFMatrix4f.lookAt(from,at,up);x3dom.debug.logInfo("New camera position:  "+from);view.animateTo(mat.inverse(),viewpoint);};x3dom.TurntableNavigation=function(navigationNode)
{x3dom.DefaultNavigation.call(this,navigationNode);this.panAxisX=null;this.panAxisY=null;this.panEnabled=true;};x3dom.TurntableNavigation.prototype=Object.create(x3dom.DefaultNavigation.prototype);x3dom.TurntableNavigation.prototype.constructor=x3dom.TurntableNavigation;x3dom.TurntableNavigation.prototype.onDrag=function(view,x,y,buttonState)
{navi=this.navi;if(!view._flyMat)
this.initTurnTable(view,false);var navType=navi.getType();var navRestrict=navi.getExplorationMode();if(navType==="none"||navRestrict==0){return;}
var dx=x-view._lastX;var dy=y-view._lastY;var d=null;var alpha,beta;buttonState=(!navRestrict||(navRestrict!=7&&buttonState==1))?navRestrict:buttonState;if(buttonState&1)
{alpha=(dy*2*Math.PI)/view._height;beta=(dx*2*Math.PI)/view._width;this.rotate(view,alpha,beta);}
else if(buttonState&2)
{d=(view._scene._lastMax.subtract(view._scene._lastMin)).length();d=((d<x3dom.fields.Eps)?1:d)*navi._vf.speed;var zoomAmount=d*(dx+dy)/view._height;this.zoom(view,zoomAmount);}
else if((buttonState&4)&&this.panEnabled==true)
{d=(view._scene._lastMax.subtract(view._scene._lastMin)).length();d=((d<x3dom.fields.Eps)?1:d)*navi._vf.speed*0.75;var tx=-d*dx/view._width;var ty=d*dy/view._height;this.pan(view,tx,ty);}
view._isMoving=true;view._dx=dx;view._dy=dy;view._lastX=x;view._lastY=y;};x3dom.TurntableNavigation.prototype.pan=function(view,tx,ty)
{if(this.target!=null){var target=this.target;var bbox=target._x3domNode.getVolume();var viewpoint=view._scene.getViewpoint();view._up=view._flyMat.e1();view._from=view._flyMat.e3();var cor=view._at;cor=cor.addScaled(this.panAxisY,ty);var temp=cor;if(cor.y>bbox.max.y||cor.y<bbox.min.y)
temp=view._at;else
view._from=view._from.addScaled(this.panAxisY,ty);cor=temp.addScaled(this.panAxisX,tx);if(cor.x>bbox.max.x||cor.x<bbox.min.x)
cor=temp;else
view._from=view._from.addScaled(this.panAxisX,tx);view._at=cor;view._flyMat=x3dom.fields.SFMatrix4f.lookAt(view._from,cor,view._up);viewpoint.setViewAbsolute(view._flyMat.inverse());}else if(this.panAxisX!=null&&this.panAxisY!=null){var viewpoint=view._scene.getViewpoint();view._up=view._flyMat.e1();view._from=view._flyMat.e3();var cor=view._at;cor=cor.addScaled(this.panAxisY,ty);var temp=cor;view._from=view._from.addScaled(this.panAxisY,ty);cor=temp.addScaled(this.panAxisX,tx);view._from=view._from.addScaled(this.panAxisX,tx);view._at=cor;view._flyMat=x3dom.fields.SFMatrix4f.lookAt(view._from,cor,view._up);viewpoint.setViewAbsolute(view._flyMat.inverse());}else{var vec=new x3dom.fields.SFVec3f(-tx*navi._vf.speed,-ty*navi._vf.speed,0);view._movement=view._movement.add(vec);var mat=view.getViewpointMatrix().mult(view._transMat);view._transMat=mat.inverse().mult(x3dom.fields.SFMatrix4f.translation(view._movement)).mult(mat);}};x3dom.TurntableNavigation.prototype.rotate=function(view,alpha,beta)
{var viewpoint=view._scene.getViewpoint();view._flyMat=this.calcOrbit(view,alpha,beta);viewpoint.setView(view._flyMat.inverse());};x3dom.TurntableNavigation.prototype.zoom=function(view,zoomAmount)
{var navi=this.navi;var viewpoint=view._scene.getViewpoint();view._up=view._flyMat.e1();view._from=view._flyMat.e3();cor=view._at;var lastDir=cor.subtract(view._from);var lastDirL=lastDir.length();lastDir=lastDir.normalize();var newDist=Math.min(zoomAmount,lastDirL-navi._vf.typeParams[6]);newDist=Math.max(newDist,lastDirL-navi._vf.typeParams[7]);view._from=view._from.addScaled(lastDir,newDist);view._flyMat=x3dom.fields.SFMatrix4f.lookAt(view._from,cor,view._up);viewpoint.setView(view._flyMat.inverse());};x3dom.TurntableNavigation.prototype.calcOrbit=function(view,alpha,beta)
{navi=this.navi;view._up=view._flyMat.e1();view._from=view._flyMat.e3();var offset=view._from.subtract(view._at);var phi=Math.atan2(offset.x,offset.z);var theta=Math.atan2(Math.sqrt(offset.x*offset.x+offset.z*offset.z),offset.y);phi-=beta;theta-=alpha;theta=Math.max(navi._vf.typeParams[2],Math.min(navi._vf.typeParams[3],theta));if(navi._vf.typeParams[4]<=navi._vf.typeParams[5])
phi=Math.max(navi._vf.typeParams[4],Math.min(navi._vf.typeParams[5],phi));else{if(beta>0&&phi<navi._vf.typeParams[4]&&phi>navi._vf.typeParams[5])phi=navi._vf.typeParams[4];else if(beta<0&&phi>navi._vf.typeParams[5]&&phi<navi._vf.typeParams[4])phi=navi._vf.typeParams[5];}
var radius=offset.length();var rSinPhi=radius*Math.sin(theta);offset.x=rSinPhi*Math.sin(phi);offset.y=radius*Math.cos(theta);offset.z=rSinPhi*Math.cos(phi);offset=view._at.add(offset);theta-=Math.PI/2;var sinPhi=Math.sin(theta);var cosPhi=Math.cos(theta);var up=new x3dom.fields.SFVec3f(sinPhi*Math.sin(phi),cosPhi,sinPhi*Math.cos(phi));if(up.y<0)
up=up.negate();return x3dom.fields.SFMatrix4f.lookAt(offset,view._at,up);};x3dom.TurntableNavigation.prototype.initTurnTable=function(view,flyTo)
{var navi=this.navi;flyTo=(flyTo==undefined)?true:flyTo;var currViewMat=view.getViewMatrix();var viewpoint=view._scene.getViewpoint();var center=x3dom.fields.SFVec3f.copy(viewpoint.getCenterOfRotation());view._flyMat=currViewMat.inverse();view._from=viewpoint._vf.position;view._at=center;view._up=new x3dom.fields.SFVec3f(0,1,0);view._flyMat=x3dom.fields.SFMatrix4f.lookAt(view._from,view._at,view._up);view._flyMat=this.calcOrbit(view,0,0);var dur=0.0;if(flyTo){dur=0.2/navi._vf.speed;}
view.animateTo(view._flyMat.inverse(),viewpoint,dur);view.resetNavHelpers();};x3dom.TurntableNavigation.prototype.onMousePress=function(view,x,y,buttonState)
{if(!view._flyMat)
this.initTurnTable(view,false);};x3dom.TurntableNavigation.prototype.init=function(view,flyTo)
{this.initTurnTable(view,false);};x3dom.TurntableNavigation.prototype.resetView=function(view)
{view._mixer.beginTime=view._lastTS;view._mixer.endTime=view._lastTS+this.navi._vf.transitionTime;view._mixer.setBeginMatrix(view.getViewMatrix());var target=view._scene.getViewpoint();target.resetView();target=x3dom.fields.SFMatrix4f.lookAt(target._vf.position,target.getCenterOfRotation(),new x3dom.fields.SFVec3f(0,1,0));view._mixer.setEndMatrix(target.inverse());this.updateFlyMat(view);}
x3dom.TurntableNavigation.prototype.updateFlyMat=function(view,nextViewpoint)
{if(!view._flyMat)
this.initTurnTable(view,false);var currViewMat=view.getViewMatrix();var viewpoint=nextViewpoint;if(viewpoint==null||!x3dom.isa(viewpoint,x3dom.nodeTypes.X3DViewpointNode))
viewpoint=view._scene.getViewpoint();var center=x3dom.fields.SFVec3f.copy(viewpoint.getCenterOfRotation());view._flyMat=currViewMat.inverse();view._from=viewpoint._vf.position;view._at=center;view._up=new x3dom.fields.SFVec3f(0,1,0);view._flyMat=x3dom.fields.SFMatrix4f.lookAt(view._from,view._at,view._up);}
x3dom.TurntableNavigation.prototype.animateTo=function(view,target,prev,dur)
{var navi=this.navi;var targetMat;if(x3dom.isa(target,x3dom.nodeTypes.X3DViewpointNode)){targetMat=x3dom.fields.SFMatrix4f.lookAt(target._vf.position,target.getCenterOfRotation(),new x3dom.fields.SFVec3f(0,1,0));}else
targetMat=target;if(navi._vf.transitionType[0].toLowerCase()!=="teleport"&&dur!=0&&navi.getType()!=="game")
{if(prev&&x3dom.isa(prev,x3dom.nodeTypes.X3DViewpointNode)){prev=prev.getViewMatrix().mult(prev.getCurrentTransform().inverse()).mult(view._transMat).mult(view._rotMat);view._mixer.beginTime=view._lastTS;if(arguments.length>=4&&arguments[3]!=null){view._mixer.endTime=view._lastTS+dur;}
else{view._mixer.endTime=view._lastTS+navi._vf.transitionTime;}
view._mixer.setBeginMatrix(prev);view._mixer.setEndMatrix(targetMat.inverse());view._scene.getViewpoint().setView(prev);}
else{view._scene.getViewpoint().setView(targetMat.inverse());}}
else
{view._scene.getViewpoint().setView(target);}
view._rotMat=x3dom.fields.SFMatrix4f.identity();view._transMat=x3dom.fields.SFMatrix4f.identity();view._movement=new x3dom.fields.SFVec3f(0,0,0);view._needNavigationMatrixUpdate=true;this.updateFlyMat(view,target);}
x3dom.TurntableNavigation.prototype.onTouchStart=function(view,evt,touches)
{console.log("touchStart "+evt.touches.length);console.log(evt);view._numTouches=evt.touches.length;view._lastX=evt.touches[0].screenX;view._lastY=evt.touches[0].screenY;};x3dom.TurntableNavigation.prototype.onTouchDrag=function(view,evt,touches,translation,rotation)
{if(view._currentInputType==x3dom.InputTypes.NAVIGATION)
{if(evt.touches.length==1)
{var dx=(evt.touches[0].screenX-view._lastX);var dy=(evt.touches[0].screenY-view._lastY);var alpha=(dy*2*Math.PI)/view._height;var beta=(dx*2*Math.PI)/view._width;this.rotate(view,alpha,beta);view._lastX=evt.touches[0].screenX;view._lastY=evt.touches[0].screenY;}
else if(evt.touches.length>=2)
{if(this.panEnabled==true)
this.pan(view,-translation.x*4.0,-translation.y*4.0);this.zoom(view,translation.z*4.0);}}};x3dom.TurntableNavigation.prototype.onTouchEnd=function(view,evt,touches)
{console.log("touchEnd "+evt.touches.length);console.log(evt);if(view._numTouches==2&&evt.touches.length==1){view._lastX=evt.touches[0].screenX;view._lastY=evt.touches[0].screenY;}
view._numTouches=evt.touches.length;};x3dom.TurntableNavigation.prototype.onDoubleClick=function(view,x,y)
{};x3dom.TurntableNavigation.prototype.setPanTarget=function(target)
{this.target=target;};x3dom.TurntableNavigation.prototype.setPanAxis=function(a,b)
{this.panAxisX=a;this.panAxisY=b;};x3dom.TurntableNavigation.prototype.setPanEnabled=function(enabled)
{this.panEnabled=enabled;};x3dom.registerNodeType("X3DFontStyleNode","Text",defineClass(x3dom.nodeTypes.X3DNode,function(ctx){x3dom.nodeTypes.X3DFontStyleNode.superClass.call(this,ctx);}));x3dom.registerNodeType("FontStyle","Text",defineClass(x3dom.nodeTypes.X3DFontStyleNode,function(ctx){x3dom.nodeTypes.FontStyle.superClass.call(this,ctx);this.addField_MFString(ctx,'family',['SERIF']);this.addField_SFBool(ctx,'horizontal',true);this.addField_MFString(ctx,'justify',['MIDDLE','MIDDLE']);this.addField_SFString(ctx,'language',"");this.addField_SFBool(ctx,'leftToRight',true);this.addField_SFFloat(ctx,'size',1.0);this.addField_SFFloat(ctx,'spacing',1.0);this.addField_SFString(ctx,'style',"PLAIN");this.addField_SFBool(ctx,'topToBottom',true);this.addField_SFFloat(ctx,'quality',2.0);},{fieldChanged:function(fieldName){if(fieldName=='family'||fieldName=='horizontal'||fieldName=='justify'||fieldName=='language'||fieldName=='leftToRight'||fieldName=='size'||fieldName=='spacing'||fieldName=='style'||fieldName=='topToBottom'){Array.forEach(this._parentNodes,function(node){Array.forEach(node._parentNodes,function(textnode){textnode.setAllDirty();});});}}}));x3dom.nodeTypes.FontStyle.defaultNode=function(){if(!x3dom.nodeTypes.FontStyle._defaultNode){x3dom.nodeTypes.FontStyle._defaultNode=new x3dom.nodeTypes.FontStyle();x3dom.nodeTypes.FontStyle._defaultNode.nodeChanged();}
return x3dom.nodeTypes.FontStyle._defaultNode;};x3dom.registerNodeType("Text","Text",defineClass(x3dom.nodeTypes.X3DGeometryNode,function(ctx){x3dom.nodeTypes.Text.superClass.call(this,ctx);this.addField_MFString(ctx,'string',[]);this.addField_MFFloat(ctx,'length',[]);this.addField_SFFloat(ctx,'maxExtent',0.0);this.addField_SFNode('fontStyle',x3dom.nodeTypes.X3DFontStyleNode);this._mesh._positions[0]=[0,0,0,1,0,0,1,1,0,0,1,0];this._mesh._normals[0]=[0,0,1,0,0,1,0,0,1,0,0,1];this._mesh._texCoords[0]=[0,0,1,0,1,1,0,1];this._mesh._colors[0]=[];this._mesh._indices[0]=[0,1,2,2,3,0];this._mesh._invalidate=true;this._mesh._numFaces=2;this._mesh._numCoords=4;},{nodeChanged:function(){if(!this._cf.fontStyle.node){this.addChild(x3dom.nodeTypes.FontStyle.defaultNode());}
this.invalidateVolume();},fieldChanged:function(fieldName){if(fieldName=='string'||fieldName=='length'||fieldName=='maxExtent'){this.invalidateVolume();Array.forEach(this._parentNodes,function(node){node.setAllDirty();});}}}));x3dom.registerNodeType("X3DSoundNode","Sound",defineClass(x3dom.nodeTypes.X3DChildNode,function(ctx){x3dom.nodeTypes.X3DSoundNode.superClass.call(this,ctx);}));x3dom.registerNodeType("Sound","Sound",defineClass(x3dom.nodeTypes.X3DSoundNode,function(ctx){x3dom.nodeTypes.Sound.superClass.call(this,ctx);this.addField_SFNode('source',x3dom.nodeTypes.X3DSoundSourceNode);},{nodeChanged:function()
{if(this._cf.source.node||!this._xmlNode){return;}
x3dom.debug.logInfo("No AudioClip child node given, searching for &lt;audio&gt; elements...");try{Array.forEach(this._xmlNode.childNodes,function(childDomNode){if(childDomNode.nodeType===1)
{x3dom.debug.logInfo("### Found &lt;"+childDomNode.nodeName+"&gt; tag.");if(childDomNode.localName.toLowerCase()==="audio")
{var loop=childDomNode.getAttribute("loop");loop=loop?(loop.toLowerCase()==="loop"):false;var newNode=childDomNode.cloneNode(false);childDomNode.parentNode.removeChild(childDomNode);childDomNode=null;if(navigator.appName!="Microsoft Internet Explorer"){document.body.appendChild(newNode);}
var startAudio=function(){newNode.play();};var audioDone=function(){if(loop){newNode.play();}};newNode.addEventListener("canplaythrough",startAudio,true);newNode.addEventListener("ended",audioDone,true);}}});}
catch(e){x3dom.debug.logException(e);}}}));x3dom.registerNodeType("X3DSoundSourceNode","Sound",defineClass(x3dom.nodeTypes.X3DTimeDependentNode,function(ctx){x3dom.nodeTypes.X3DSoundSourceNode.superClass.call(this,ctx);}));x3dom.registerNodeType("AudioClip","Sound",defineClass(x3dom.nodeTypes.X3DSoundSourceNode,function(ctx){x3dom.nodeTypes.AudioClip.superClass.call(this,ctx);this.addField_MFString(ctx,'url',[]);this.addField_SFBool(ctx,'enabled',false);this.addField_SFBool(ctx,'loop',false);this._audio=document.createElement('audio');if(navigator.appName!="Microsoft Internet Explorer"){document.body.appendChild(this._audio);}
this._sources=[];},{nodeChanged:function()
{this._createSources=function()
{this._sources=[];for(var i=0;i<this._vf.url.length;i++)
{var audioUrl=this._nameSpace.getURL(this._vf.url[i]);x3dom.debug.logInfo('Adding sound file: '+audioUrl);var src=document.createElement('source');src.setAttribute('src',audioUrl);this._sources.push(src);this._audio.appendChild(src);}};var that=this;this._startAudio=function()
{that._audio.loop=that._vf.loop?"loop":"";if(that._vf.enabled===true)
{that._audio.play();}};this._stopAudio=function()
{that._audio.pause();};this._audioEnded=function()
{if(that._vf.enabled===true&&that._vf.loop===true)
{that._startAudio();}};var log=function(e)
{x3dom.debug.logWarning("MediaEvent error:"+e);};this._audio.addEventListener("canplaythrough",this._startAudio,true);this._audio.addEventListener("ended",this._audioEnded,true);this._audio.addEventListener("error",log,true);this._audio.addEventListener("pause",this._audioEnded,true);this._createSources();},fieldChanged:function(fieldName)
{if(fieldName==="enabled")
{if(this._vf.enabled===true)
{this._startAudio();}
else
{this._stopAudio();}}
else if(fieldName==="loop")
{}
else if(fieldName==="url")
{this._stopAudio();while(this._audio.hasChildNodes())
{this._audio.removeChild(this._audio.firstChild);}
for(var i=0;i<this._vf.url.length;i++)
{var audioUrl=this._nameSpace.getURL(this._vf.url[i]);x3dom.debug.logInfo('Adding sound file: '+audioUrl);var src=document.createElement('source');src.setAttribute('src',audioUrl);this._audio.appendChild(src);}}},shutdown:function(){if(this._audio){this._audio.pause();while(this._audio.hasChildNodes()){this._audio.removeChild(this._audio.firstChild);}
document.body.removeChild(this._audio);this._audio=null;}}}));x3dom.registerNodeType("X3DTextureTransformNode","Texturing",defineClass(x3dom.nodeTypes.X3DAppearanceChildNode,function(ctx){x3dom.nodeTypes.X3DTextureTransformNode.superClass.call(this,ctx);}));x3dom.registerNodeType("TextureTransform","Texturing",defineClass(x3dom.nodeTypes.X3DTextureTransformNode,function(ctx){x3dom.nodeTypes.TextureTransform.superClass.call(this,ctx);this.addField_SFVec2f(ctx,'center',0,0);this.addField_SFFloat(ctx,'rotation',0);this.addField_SFVec2f(ctx,'scale',1,1);this.addField_SFVec2f(ctx,'translation',0,0);var negCenter=new x3dom.fields.SFVec3f(-this._vf.center.x,-this._vf.center.y,1);var posCenter=new x3dom.fields.SFVec3f(this._vf.center.x,this._vf.center.y,0);var trans3=new x3dom.fields.SFVec3f(this._vf.translation.x,this._vf.translation.y,0);var scale3=new x3dom.fields.SFVec3f(this._vf.scale.x,this._vf.scale.y,0);this._trafo=x3dom.fields.SFMatrix4f.translation(negCenter).mult(x3dom.fields.SFMatrix4f.scale(scale3)).mult(x3dom.fields.SFMatrix4f.rotationZ(this._vf.rotation)).mult(x3dom.fields.SFMatrix4f.translation(posCenter.add(trans3)));},{fieldChanged:function(fieldName){if(fieldName=='center'||fieldName=='rotation'||fieldName=='scale'||fieldName=='translation'){var negCenter=new x3dom.fields.SFVec3f(-this._vf.center.x,-this._vf.center.y,1);var posCenter=new x3dom.fields.SFVec3f(this._vf.center.x,this._vf.center.y,0);var trans3=new x3dom.fields.SFVec3f(this._vf.translation.x,this._vf.translation.y,0);var scale3=new x3dom.fields.SFVec3f(this._vf.scale.x,this._vf.scale.y,0);this._trafo=x3dom.fields.SFMatrix4f.translation(negCenter).mult(x3dom.fields.SFMatrix4f.scale(scale3)).mult(x3dom.fields.SFMatrix4f.rotationZ(this._vf.rotation)).mult(x3dom.fields.SFMatrix4f.translation(posCenter.add(trans3)));}},texTransformMatrix:function(){return this._trafo;}}));x3dom.registerNodeType("TextureProperties","Texturing",defineClass(x3dom.nodeTypes.X3DNode,function(ctx){x3dom.nodeTypes.TextureProperties.superClass.call(this,ctx);this.addField_SFFloat(ctx,'anisotropicDegree',1.0);this.addField_SFColorRGBA(ctx,'borderColor',0,0,0,0);this.addField_SFInt32(ctx,'borderWidth',0);this.addField_SFString(ctx,'boundaryModeS',"REPEAT");this.addField_SFString(ctx,'boundaryModeT',"REPEAT");this.addField_SFString(ctx,'boundaryModeR',"REPEAT");this.addField_SFString(ctx,'magnificationFilter',"FASTEST");this.addField_SFString(ctx,'minificationFilter',"FASTEST");this.addField_SFString(ctx,'textureCompression',"FASTEST");this.addField_SFFloat(ctx,'texturePriority',0);this.addField_SFBool(ctx,'generateMipMaps',false);},{fieldChanged:function(fieldName)
{if(this._vf.hasOwnProperty(fieldName)){Array.forEach(this._parentNodes,function(texture){Array.forEach(texture._parentNodes,function(app){Array.forEach(app._parentNodes,function(shape){shape._dirty.texture=true;});});});this._nameSpace.doc.needRender=true;}}}));x3dom.registerNodeType("X3DTextureNode","Texturing",defineClass(x3dom.nodeTypes.X3DAppearanceChildNode,function(ctx){x3dom.nodeTypes.X3DTextureNode.superClass.call(this,ctx);this.addField_SFInt32(ctx,'origChannelCount',0);this.addField_MFString(ctx,'url',[]);this.addField_SFBool(ctx,'repeatS',true);this.addField_SFBool(ctx,'repeatT',true);this.addField_SFBool(ctx,'scale',true);this.addField_SFString(ctx,'crossOrigin','');this.addField_SFNode('textureProperties',x3dom.nodeTypes.TextureProperties);this._needPerFrameUpdate=false;this._isCanvas=false;this._type="diffuseMap";this._blending=(this._vf.origChannelCount==1||this._vf.origChannelCount==2);},{invalidateGLObject:function()
{Array.forEach(this._parentNodes,function(app){Array.forEach(app._parentNodes,function(shape){if(x3dom.isa(shape,x3dom.nodeTypes.X3DShapeNode)){shape._dirty.texture=true;}
else{Array.forEach(shape._parentNodes,function(realShape){if(x3dom.isa(realShape,x3dom.nodeTypes.X3DShapeNode)){realShape._dirty.texture=true;}else{Array.forEach(realShape._parentNodes,function(realShape2){if(x3dom.isa(realShape2,x3dom.nodeTypes.X3DShapeNode)){realShape2._dirty.texture=true;}});}});}});});this._nameSpace.doc.needRender=true;},parentAdded:function(parent)
{Array.forEach(parent._parentNodes,function(shape){if(x3dom.isa(shape,x3dom.nodeTypes.Shape)){shape._dirty.texture=true;}
else{Array.forEach(shape._parentNodes,function(realShape){realShape._dirty.texture=true;});}});},parentRemoved:function(parent)
{Array.forEach(parent._parentNodes,function(shape){if(x3dom.isa(shape,x3dom.nodeTypes.Shape)){shape._dirty.texture=true;}
else{Array.forEach(shape._parentNodes,function(realShape){realShape._dirty.texture=true;});}});},fieldChanged:function(fieldName)
{if(fieldName=="url"||fieldName=="origChannelCount"||fieldName=="repeatS"||fieldName=="repeatT"||fieldName=="scale"||fieldName=="crossOrigin")
{var that=this;Array.forEach(this._parentNodes,function(app){if(x3dom.isa(app,x3dom.nodeTypes.X3DAppearanceNode)){app.nodeChanged();Array.forEach(app._parentNodes,function(shape){shape._dirty.texture=true;});}
else if(x3dom.isa(app,x3dom.nodeTypes.MultiTexture)){Array.forEach(app._parentNodes,function(realApp){realApp.nodeChanged();Array.forEach(realApp._parentNodes,function(shape){shape._dirty.texture=true;});});}
else if(x3dom.isa(app,x3dom.nodeTypes.ComposedCubeMapTexture)){Array.forEach(app._parentNodes,function(realApp){realApp.nodeChanged();Array.forEach(realApp._parentNodes,function(shape){shape._dirty.texture=true;});});}
else if(x3dom.isa(app,x3dom.nodeTypes.ImageGeometry)){var cf=null;if(that._xmlNode&&that._xmlNode.hasAttribute('containerField')){cf=that._xmlNode.getAttribute('containerField');app._dirty[cf]=true;}}
else if(x3dom.nodeTypes.X3DVolumeDataNode!==undefined){if(x3dom.isa(app,x3dom.nodeTypes.X3DVolumeRenderStyleNode)){if(that._xmlNode&&that._xmlNode.hasAttribute('containerField')){if(app._volumeDataParent){app._volumeDataParent._dirty.texture=true;}else{var volumeDataParent=app._parentNodes[0];while(!x3dom.isa(volumeDataParent,x3dom.nodeTypes.X3DVolumeDataNode)&&x3dom.isa(volumeDataParent,x3dom.nodeTypes.X3DNode)){volumeDataParent=volumeDataParent._parentNodes[0];}
if(x3dom.isa(volumeDataParent,x3dom.nodeTypes.X3DNode)){volumeDataParent._dirty.texture=true;}}}}else if(x3dom.isa(app,x3dom.nodeTypes.X3DVolumeDataNode)){if(that._xmlNode&&that._xmlNode.hasAttribute('containerField')){app._dirty.texture=true;}}}});}},getTexture:function(pos){if(pos===0){return this;}
return null;},size:function(){return 1;}}));x3dom.registerNodeType("MultiTexture","Texturing",defineClass(x3dom.nodeTypes.X3DTextureNode,function(ctx){x3dom.nodeTypes.MultiTexture.superClass.call(this,ctx);this.addField_MFNode('texture',x3dom.nodeTypes.X3DTextureNode);},{getTexture:function(pos){if(pos>=0&&pos<this._cf.texture.nodes.length){return this._cf.texture.nodes[pos];}
return null;},getTextures:function(){return this._cf.texture.nodes;},size:function(){return this._cf.texture.nodes.length;}}));x3dom.registerNodeType("Texture","Texturing",defineClass(x3dom.nodeTypes.X3DTextureNode,function(ctx){x3dom.nodeTypes.Texture.superClass.call(this,ctx);this.addField_SFBool(ctx,'hideChildren',true);this._video=null;this._intervalID=0;this._canvas=null;},{nodeChanged:function()
{if(this._vf.url.length||!this._xmlNode){return;}
x3dom.debug.logInfo("No Texture URL given, searching for &lt;img&gt; elements...");var that=this;try{Array.forEach(this._xmlNode.childNodes,function(childDomNode){if(childDomNode.nodeType===1){var url=childDomNode.getAttribute("src");if(url){that._vf.url.push(url);x3dom.debug.logInfo(that._vf.url[that._vf.url.length-1]);if(childDomNode.localName.toLowerCase()==="video"){that._needPerFrameUpdate=true;that._video=document.createElement('video');that._video.setAttribute('preload','auto');that._video.setAttribute('muted','muted');var p=document.getElementsByTagName('body')[0];p.appendChild(that._video);that._video.style.display="none";that._video.style.visibility="hidden";}}
else if(childDomNode.localName.toLowerCase()==="canvas"){that._needPerFrameUpdate=true;that._isCanvas=true;that._canvas=childDomNode;}
if(childDomNode.style&&that._vf.hideChildren){childDomNode.style.display="none";childDomNode.style.visibility="hidden";}
x3dom.debug.logInfo("### Found &lt;"+childDomNode.nodeName+"&gt; tag.");}});}
catch(e){x3dom.debug.logException(e);}},shutdown:function(){if(this._video){this._video.pause();while(this._video.hasChildNodes()){this._video.removeChild(this._video.firstChild);}
document.body.removeChild(this._video);this._video=null;}}}));x3dom.registerNodeType("RenderedTexture","Texturing",defineClass(x3dom.nodeTypes.X3DTextureNode,function(ctx){x3dom.nodeTypes.RenderedTexture.superClass.call(this,ctx);if(ctx)
ctx.doc._nodeBag.renderTextures.push(this);else
x3dom.debug.logWarning("RenderedTexture: No runtime context found!");this.addField_SFNode('viewpoint',x3dom.nodeTypes.X3DViewpointNode);this.addField_SFNode('background',x3dom.nodeTypes.X3DBackgroundNode);this.addField_SFNode('fog',x3dom.nodeTypes.X3DFogNode);this.addField_SFNode('scene',x3dom.nodeTypes.X3DNode);this.addField_MFNode('excludeNodes',x3dom.nodeTypes.X3DNode);this.addField_MFInt32(ctx,'dimensions',[128,128,4]);this.addField_SFString(ctx,'update','NONE');this.addField_SFBool(ctx,'showNormals',false);this.addField_SFString(ctx,'stereoMode','NONE');this.addField_SFFloat(ctx,'interpupillaryDistance',0.064);this.addField_SFFloat(ctx,'eyeToScreenDistance',0.041);this.addField_SFFloat(ctx,'vScreenSize',0.07074);this.addField_SFVec3f(ctx,'lensCenter',0.15197,0,0);this.addField_SFBool(ctx,'depthMap',false);this.addField_SFBool(ctx,'oculusRiftVersion',1);x3dom.debug.assert(this._vf.dimensions.length>=3,"RenderedTexture.dimensions requires at least 3 entries.");this._clearParents=true;this._needRenderUpdate=true;this.checkDepthTextureSupport=function(){if(this._vf.depthMap&&x3dom.caps.DEPTH_TEXTURE===null)
x3dom.debug.logWarning("RenderedTexture Node: depth texture extension not supported");};this.checkDepthTextureSupport();},{nodeChanged:function()
{this._clearParents=true;this._needRenderUpdate=true;},fieldChanged:function(fieldName)
{switch(fieldName)
{case"excludeNodes":this._clearParents=true;break;case"update":if(this._vf.update.toUpperCase()=="NEXT_FRAME_ONLY"||this._vf.update.toUpperCase()=="ALWAYS"){this._needRenderUpdate=true;}
break;case"depthMap":this.checkDepthTextureSupport();this._x3domTexture.updateTexture();this._needRenderUpdate=true;default:break;}},getViewMatrix:function()
{if(this._clearParents&&this._cf.excludeNodes.nodes.length){var that=this;Array.forEach(this._cf.excludeNodes.nodes,function(node){for(var i=0,n=node._parentNodes.length;i<n;i++){if(node._parentNodes[i]===that){node._parentNodes.splice(i,1);node.parentRemoved(that);}}});this._clearParents=false;}
var locScene=this._cf.scene.node;var scene=this._nameSpace.doc._scene;var vbP=scene.getViewpoint();var view=this._cf.viewpoint.node;var ret_mat=null;if(view===null||view===vbP){ret_mat=this._nameSpace.doc._viewarea.getViewMatrix();}
else if(locScene&&locScene!==scene){ret_mat=view.getViewMatrix()}
else{var mat_viewpoint=view.getCurrentTransform();ret_mat=view.getViewMatrix().mult(mat_viewpoint.inverse());}
var stereoMode=this._vf.stereoMode.toUpperCase();if(stereoMode!="NONE"){var d=this._vf.interpupillaryDistance/2;if(stereoMode=="RIGHT_EYE"){d=-d;}
var modifier=new x3dom.fields.SFMatrix4f(1,0,0,d,0,1,0,0,0,0,1,0,0,0,0,1);ret_mat=modifier.mult(ret_mat);}
return ret_mat;},getProjectionMatrix:function()
{var doc=this._nameSpace.doc;var vbP=doc._scene.getViewpoint();var view=this._cf.viewpoint.node;var ret_mat=null;var f,w=this._vf.dimensions[0],h=this._vf.dimensions[1];var stereoMode=this._vf.stereoMode.toUpperCase();var stereo=(stereoMode!="NONE");if(view===null||view===vbP){ret_mat=x3dom.fields.SFMatrix4f.copy(doc._viewarea.getProjectionMatrix());if(stereo){f=2*Math.atan(this._vf.vScreenSize/(2*this._vf.eyeToScreenDistance));f=1/Math.tan(f/2);}
else{f=1/Math.tan(vbP._vf.fieldOfView/2);}
ret_mat._00=f/(w/h);ret_mat._11=f;}
else{ret_mat=view.getProjectionMatrix(w/h);}
if(stereo){var lensCenter=this._vf.lensCenter.copy();if(stereoMode=="RIGHT_EYE"){lensCenter.x=-lensCenter.x;}
var modifier=new x3dom.fields.SFMatrix4f(1,0,0,lensCenter.x,0,1,0,lensCenter.y,0,0,1,lensCenter.z,0,0,0,1);ret_mat=modifier.mult(ret_mat);}
return ret_mat;},getWCtoCCMatrix:function()
{var view=this.getViewMatrix();var proj=this.getProjectionMatrix();return proj.mult(view);},parentRemoved:function(parent)
{if(this._parentNodes.length===0){var doc=this.findX3DDoc();for(var i=0,n=doc._nodeBag.renderTextures.length;i<n;i++){if(doc._nodeBag.renderTextures[i]===this){doc._nodeBag.renderTextures.splice(i,1);}}}
if(this._cf.scene.node){this._cf.scene.node.parentRemoved(this);}},requirePingPong:function()
{return false;}}));x3dom.registerNodeType("RefinementTexture","Texturing",defineClass(x3dom.nodeTypes.RenderedTexture,function(ctx){x3dom.nodeTypes.RefinementTexture.superClass.call(this,ctx);this.addField_SFString(ctx,'stamp0',"gpuii/stamps/0.gif");this.addField_SFString(ctx,'stamp1',"gpuii/stamps/1.gif");this.addField_SFBool(ctx,'autoRefinement',true);this.addField_SFString(ctx,'format','jpg');this.addField_SFInt32(ctx,'iterations',7);this.addField_SFInt32(ctx,'maxLevel',this._vf.iterations);if(this._vf.iterations%2===0){var temp=this._vf.stamp0;this._vf.stamp0=this._vf.stamp1;this._vf.stamp1=temp;}
this._vf.iterations=(this._vf.iterations>11)?11:this._vf.iterations;this._vf.iterations=(this._vf.iterations<3)?3:this._vf.iterations;this._vf.maxLevel=(this._vf.maxLevel>11)?11:this._vf.maxLevel;this._vf.maxLevel=(this._vf.maxLevel<3)?3:this._vf.maxLevel;this._vf.maxLevel=(this._vf.maxLevel>this._vf.iterations)?this._vf.iterations:this._vf.maxLevel;var repeatConfig=[{x:4,y:8},{x:8,y:8},{x:8,y:16},{x:16,y:16},{x:16,y:32},{x:32,y:32},{x:32,y:64},{x:64,y:64},{x:64,y:128}];this._repeat=new x3dom.fields.SFVec2f(this._vf.dimensions[0]/repeatConfig[this._vf.iterations-3].x,this._vf.dimensions[1]/repeatConfig[this._vf.iterations-3].y);this._renderedImage=0;this._currLoadLevel=0;this._loadLevel=1;},{nextLevel:function(){if(this._loadLevel<this._vf.maxLevel){this._loadLevel++;this._nameSpace.doc.needRender=true;}},requirePingPong:function(){return(this._currLoadLevel<=this._vf.maxLevel&&this._renderedImage<this._loadLevel);}}));x3dom.registerNodeType("PixelTexture","Texturing",defineClass(x3dom.nodeTypes.X3DTextureNode,function(ctx){x3dom.nodeTypes.PixelTexture.superClass.call(this,ctx);this.addField_SFImage(ctx,'image',0,0,0);},{fieldChanged:function(fieldName)
{if(fieldName=="image"){this.invalidateGLObject();}},getWidth:function(){return this._vf.image.width;},getHeight:function(){return this._vf.image.height;},getComponents:function(){return this._vf.image.comp;},setPixel:function(x,y,color,update){update=(update==undefined)?true:update;if(this._x3domTexture){this._x3domTexture.setPixel(x,y,[color.r*255,color.g*255,color.b*255,color.a*255],update);this._vf.image.setPixel(x,y,color);}
else
{this._vf.image.setPixel(x,y,color);if(update){this.invalidateGLObject();}}},getPixel:function(x,y){return this._vf.image.getPixel(x,y);},setPixels:function(pixels,update){update=(update==undefined)?true:update;this._vf.image.setPixels(pixels);if(update){this.invalidateGLObject();}},getPixels:function(){return this._vf.image.getPixels();}}));x3dom.registerNodeType("ImageTexture","Texturing",defineClass(x3dom.nodeTypes.Texture,function(ctx){x3dom.nodeTypes.ImageTexture.superClass.call(this,ctx);}));x3dom.registerNodeType("MovieTexture","Texturing",defineClass(x3dom.nodeTypes.Texture,function(ctx){x3dom.nodeTypes.MovieTexture.superClass.call(this,ctx);this.addField_SFBool(ctx,'loop',false);this.addField_SFFloat(ctx,'speed',1.0);this.addField_SFTime(ctx,'pauseTime',0);this.addField_SFFloat(ctx,'pitch',1.0);this.addField_SFTime(ctx,'resumeTime',0);this.addField_SFTime(ctx,'startTime',0);this.addField_SFTime(ctx,'stopTime',0);}));x3dom.registerNodeType("X3DTextureCoordinateNode","Texturing",defineClass(x3dom.nodeTypes.X3DGeometricPropertyNode,function(ctx){x3dom.nodeTypes.X3DTextureCoordinateNode.superClass.call(this,ctx);},{fieldChanged:function(fieldName){if(fieldName==="texCoord"||fieldName==="point"||fieldName==="parameter"||fieldName==="mode")
{Array.forEach(this._parentNodes,function(node){node.fieldChanged("texCoord");});}},parentAdded:function(parent){if(parent._mesh&&parent._cf.texCoord.node!==this){parent.fieldChanged("texCoord");}}}));x3dom.registerNodeType("TextureCoordinate","Texturing",defineClass(x3dom.nodeTypes.X3DTextureCoordinateNode,function(ctx){x3dom.nodeTypes.TextureCoordinate.superClass.call(this,ctx);this.addField_MFVec2f(ctx,'point',[]);}));x3dom.registerNodeType("TextureCoordinateGenerator","Texturing",defineClass(x3dom.nodeTypes.X3DTextureCoordinateNode,function(ctx){x3dom.nodeTypes.TextureCoordinateGenerator.superClass.call(this,ctx);this.addField_SFString(ctx,'mode',"SPHERE");this.addField_MFFloat(ctx,'parameter',[]);}));x3dom.registerNodeType("MultiTextureCoordinate","Texturing",defineClass(x3dom.nodeTypes.X3DTextureCoordinateNode,function(ctx){x3dom.nodeTypes.MultiTextureCoordinate.superClass.call(this,ctx);this.addField_MFNode('texCoord',x3dom.nodeTypes.X3DTextureCoordinateNode);}));x3dom.registerNodeType("ImageTextureAtlas","Texturing",defineClass(x3dom.nodeTypes.Texture,function(ctx){x3dom.nodeTypes.ImageTextureAtlas.superClass.call(this,ctx);this.addField_SFInt32(ctx,'numberOfSlices',0);this.addField_SFInt32(ctx,'slicesOverX',0);this.addField_SFInt32(ctx,'slicesOverY',0);}));x3dom.registerNodeType("X3DEnvironmentTextureNode","CubeMapTexturing",defineClass(x3dom.nodeTypes.X3DTextureNode,function(ctx){x3dom.nodeTypes.X3DEnvironmentTextureNode.superClass.call(this,ctx);},{getTexUrl:function(){return[];},getTexSize:function(){return-1;}}));x3dom.registerNodeType("ComposedCubeMapTexture","CubeMapTexturing",defineClass(x3dom.nodeTypes.X3DEnvironmentTextureNode,function(ctx){x3dom.nodeTypes.ComposedCubeMapTexture.superClass.call(this,ctx);this.addField_SFNode('back',x3dom.nodeTypes.Texture);this.addField_SFNode('front',x3dom.nodeTypes.Texture);this.addField_SFNode('bottom',x3dom.nodeTypes.Texture);this.addField_SFNode('top',x3dom.nodeTypes.Texture);this.addField_SFNode('left',x3dom.nodeTypes.Texture);this.addField_SFNode('right',x3dom.nodeTypes.Texture);this._type="environmentMap";},{getTexUrl:function(){return[this._nameSpace.getURL(this._cf.back.node._vf.url[0]),this._nameSpace.getURL(this._cf.front.node._vf.url[0]),this._nameSpace.getURL(this._cf.bottom.node._vf.url[0]),this._nameSpace.getURL(this._cf.top.node._vf.url[0]),this._nameSpace.getURL(this._cf.left.node._vf.url[0]),this._nameSpace.getURL(this._cf.right.node._vf.url[0])];}}));x3dom.registerNodeType("GeneratedCubeMapTexture","CubeMapTexturing",defineClass(x3dom.nodeTypes.X3DEnvironmentTextureNode,function(ctx){x3dom.nodeTypes.GeneratedCubeMapTexture.superClass.call(this,ctx);this.addField_SFInt32(ctx,'size',128);this.addField_SFString(ctx,'update','NONE');this._type="cubeMap";x3dom.debug.logWarning("GeneratedCubeMapTexture NYI");},{getTexSize:function(){return this._vf.size;}}));x3dom.registerNodeType("Uniform","Shaders",defineClass(x3dom.nodeTypes.Field,function(ctx){x3dom.nodeTypes.Uniform.superClass.call(this,ctx);}));x3dom.registerNodeType("SurfaceShaderTexture","Shaders",defineClass(x3dom.nodeTypes.X3DTextureNode,function(ctx){x3dom.nodeTypes.SurfaceShaderTexture.superClass.call(this,ctx);this.addField_SFInt32(ctx,'textureCoordinatesId',0);this.addField_SFString(ctx,'channelMask',"DEFAULT");this.addField_SFBool(ctx,'isSRGB',false);this.addField_SFNode('texture',x3dom.nodeTypes.X3DTextureNode);this.addField_SFNode('textureTransform',x3dom.nodeTypes.X3DTextureTransformNode);}));x3dom.registerNodeType("X3DShaderNode","Shaders",defineClass(x3dom.nodeTypes.X3DAppearanceChildNode,function(ctx){x3dom.nodeTypes.X3DShaderNode.superClass.call(this,ctx);this.addField_SFString(ctx,'language',"");}));x3dom.registerNodeType("CommonSurfaceShader","Shaders",defineClass(x3dom.nodeTypes.X3DShaderNode,function(ctx){x3dom.nodeTypes.CommonSurfaceShader.superClass.call(this,ctx);this.addField_SFInt32(ctx,'tangentTextureCoordinatesId',-1);this.addField_SFInt32(ctx,'binormalTextureCoordinatesId',-1);this.addField_SFVec3f(ctx,'emissiveFactor',0,0,0);this.addField_SFInt32(ctx,'emissiveTextureId',-1);this.addField_SFInt32(ctx,'emissiveTextureCoordinatesId',0);this.addField_SFString(ctx,'emissiveTextureChannelMask','rgb');this.addField_SFVec3f(ctx,'ambientFactor',0.2,0.2,0.2);this.addField_SFInt32(ctx,'ambientTextureId',-1);this.addField_SFInt32(ctx,'ambientTextureCoordinatesId',0);this.addField_SFString(ctx,'ambientTextureChannelMask','rgb');this.addField_SFVec3f(ctx,'diffuseFactor',0.8,0.8,0.8);this.addField_SFInt32(ctx,'diffuseTextureId',-1);this.addField_SFInt32(ctx,'diffuseTextureCoordinatesId',0);this.addField_SFString(ctx,'diffuseTextureChannelMask','rgb');this.addField_SFVec3f(ctx,'specularFactor',0,0,0);this.addField_SFInt32(ctx,'specularTextureId',-1);this.addField_SFInt32(ctx,'specularTextureCoordinatesId',0);this.addField_SFString(ctx,'specularTextureChannelMask','rgb');this.addField_SFFloat(ctx,'shininessFactor',0.2);this.addField_SFInt32(ctx,'shininessTextureId',-1);this.addField_SFInt32(ctx,'shininessTextureCoordinatesId',0);this.addField_SFString(ctx,'shininessTextureChannelMask','a');this.addField_SFString(ctx,'normalFormat','UNORM');this.addField_SFString(ctx,'normalSpace','TANGENT');this.addField_SFInt32(ctx,'normalTextureId',-1);this.addField_SFInt32(ctx,'normalTextureCoordinatesId',0);this.addField_SFString(ctx,'normalTextureChannelMask','rgb');this.addField_SFVec3f(ctx,'reflectionFactor',0,0,0);this.addField_SFInt32(ctx,'reflectionTextureId',-1);this.addField_SFInt32(ctx,'reflectionTextureCoordinatesId',0);this.addField_SFString(ctx,'reflectionTextureChannelMask','rgb');this.addField_SFVec3f(ctx,'transmissionFactor',0,0,0);this.addField_SFInt32(ctx,'transmissionTextureId',-1);this.addField_SFInt32(ctx,'transmissionTextureCoordinatesId',0);this.addField_SFString(ctx,'transmissionTextureChannelMask','rgb');this.addField_SFVec3f(ctx,'environmentFactor',1,1,1);this.addField_SFInt32(ctx,'environmentTextureId',-1);this.addField_SFInt32(ctx,'environmentTextureCoordinatesId',0);this.addField_SFString(ctx,'environmentTextureChannelMask','rgb');this.addField_SFFloat(ctx,'relativeIndexOfRefraction',1);this.addField_SFFloat(ctx,'fresnelBlend',0);this.addField_SFString(ctx,'displacementAxis','y');this.addField_SFFloat(ctx,'displacementFactor',255.0);this.addField_SFInt32(ctx,'displacementTextureId',-1);this.addField_SFInt32(ctx,'displacementTextureCoordinatesId',0);this.addField_SFNode('emissiveTexture',x3dom.nodeTypes.X3DTextureNode);this.addField_SFNode('ambientTexture',x3dom.nodeTypes.X3DTextureNode);this.addField_SFNode('diffuseTexture',x3dom.nodeTypes.X3DTextureNode);this.addField_SFNode('specularTexture',x3dom.nodeTypes.X3DTextureNode);this.addField_SFNode('shininessTexture',x3dom.nodeTypes.X3DTextureNode);this.addField_SFNode('normalTexture',x3dom.nodeTypes.X3DTextureNode);this.addField_SFNode('reflectionTexture',x3dom.nodeTypes.X3DTextureNode);this.addField_SFNode('transmissionTexture',x3dom.nodeTypes.X3DTextureNode);this.addField_SFNode('environmentTexture',x3dom.nodeTypes.X3DTextureNode);this.addField_SFNode('displacementTexture',x3dom.nodeTypes.X3DTextureNode);this.addField_SFNode('diffuseDisplacementTexture',x3dom.nodeTypes.X3DTextureNode);this.addField_SFNode('multiDiffuseAlphaTexture',x3dom.nodeTypes.X3DTextureNode);this.addField_SFNode('multiSpecularShininessTexture',x3dom.nodeTypes.X3DTextureNode);this.addField_SFNode('multiEmissiveAmbientTexture',x3dom.nodeTypes.X3DTextureNode);this.addField_SFNode('multiVisibilityTexture',x3dom.nodeTypes.X3DTextureNode);this.addField_SFVec3f(ctx,'normalScale',2,2,2);this.addField_SFVec3f(ctx,'normalBias',-1,-1,-1);this.addField_SFFloat(ctx,'alphaFactor',1);this.addField_SFBool(ctx,'invertAlphaTexture',false);this.addField_SFInt32(ctx,'alphaTextureId',-1);this.addField_SFInt32(ctx,'alphaTextureCoordinatesId',0);this.addField_SFString(ctx,'alphaTextureChannelMask','a');this.addField_SFNode('alphaTexture',x3dom.nodeTypes.X3DTextureNode);this._dirty={};},{getDiffuseMap:function()
{if(this._cf.diffuseTexture.node){if(x3dom.isa(this._cf.diffuseTexture.node,x3dom.nodeTypes.SurfaceShaderTexture)){this._cf.diffuseTexture.node._cf.texture.node._type="diffuseMap";return this._cf.diffuseTexture.node._cf.texture.node;}else{this._cf.diffuseTexture.node._type="diffuseMap";return this._cf.diffuseTexture.node;}}else{return null;}},getEnvironmentMap:function()
{if(this._cf.environmentTexture.node){if(x3dom.isa(this._cf.environmentTexture.node,x3dom.nodeTypes.SurfaceShaderTexture)){this._cf.environmentTexture.node._cf.texture.node._type="environmentMap";return this._cf.environmentTexture.node._cf.texture.node;}else{this._cf.environmentTexture.node._type="environmentMap";return this._cf.environmentTexture.node;}}else{return null;}},getNormalMap:function()
{if(this._cf.normalTexture.node){if(x3dom.isa(this._cf.normalTexture.node,x3dom.nodeTypes.SurfaceShaderTexture)){this._cf.normalTexture.node._cf.texture.node._type="normalMap";return this._cf.normalTexture.node._cf.texture.node;}else{this._cf.normalTexture.node._type="normalMap";return this._cf.normalTexture.node;}}else{return null;}},getAmbientMap:function()
{if(this._cf.ambientTexture.node){if(x3dom.isa(this._cf.ambientTexture.node,x3dom.nodeTypes.SurfaceShaderTexture)){this._cf.ambientTexture.node._cf.texture.node._type="ambientMap";return this._cf.ambientTexture.node._cf.texture.node;}else{this._cf.ambientTexture.node._type="ambientMap";return this._cf.ambientTexture.node;}}else{return null;}},getSpecularMap:function()
{if(this._cf.specularTexture.node){if(x3dom.isa(this._cf.specularTexture.node,x3dom.nodeTypes.SurfaceShaderTexture)){this._cf.specularTexture.node._cf.texture.node._type="specularMap";return this._cf.specularTexture.node._cf.texture.node;}else{this._cf.specularTexture.node._type="specularMap";return this._cf.specularTexture.node;}}else{return null;}},getShininessMap:function()
{if(this._cf.shininessTexture.node){if(x3dom.isa(this._cf.shininessTexture.node,x3dom.nodeTypes.SurfaceShaderTexture)){this._cf.shininessTexture.node._cf.texture.node._type="shininessMap";return this._cf.shininessTexture.node._cf.texture.node;}else{this._cf.shininessTexture.node._type="shininessMap";return this._cf.shininessTexture.node;}}else{return null;}},getAlphaMap:function()
{if(this._cf.alphaTexture.node){if(x3dom.isa(this._cf.alphaTexture.node,x3dom.nodeTypes.SurfaceShaderTexture)){this._cf.alphaTexture.node._cf.texture.node._type="alphaMap";return this._cf.alphaTexture.node._cf.texture.node;}else{this._cf.alphaTexture.node._type="alphaMap";return this._cf.alphaTexture.node;}}else{return null;}},getDisplacementMap:function()
{if(this._cf.displacementTexture.node){if(x3dom.isa(this._cf.displacementTexture.node,x3dom.nodeTypes.SurfaceShaderTexture)){this._cf.displacementTexture.node._cf.texture.node._type="displacementMap";return this._cf.displacementTexture.node._cf.texture.node;}else{this._cf.displacementTexture.node._type="displacementMap";return this._cf.displacementTexture.node;}}else{return null;}},getDiffuseDisplacementMap:function()
{if(this._cf.diffuseDisplacementTexture.node){if(x3dom.isa(this._cf.diffuseDisplacementTexture.node,x3dom.nodeTypes.SurfaceShaderTexture)){this._cf.diffuseDisplacementTexture.node._cf.texture.node._type="diffuseDisplacementMap";return this._cf.diffuseDisplacementTexture.node._cf.texture.node;}else{this._cf.diffuseDisplacementTexture.node._type="diffuseDisplacementMap";return this._cf.diffuseDisplacementTexture.node;}}else{return null;}},getMultiDiffuseAlphaMap:function()
{if(this._cf.multiDiffuseAlphaTexture.node){if(x3dom.isa(this._cf.multiDiffuseAlphaTexture.node,x3dom.nodeTypes.SurfaceShaderTexture)){this._cf.multiDiffuseAlphaTexture.node._cf.texture.node._type="multiDiffuseAlphaMap";return this._cf.multiDiffuseAlphaTexture.node._cf.texture.node;}else{this._cf.multiDiffuseAlphaTexture.node._type="multiDiffuseAlphaMap";return this._cf.multiDiffuseAlphaTexture.node;}}else{return null;}},getMultiEmissiveAmbientMap:function()
{if(this._cf.multiEmissiveAmbientTexture.node){if(x3dom.isa(this._cf.multiEmissiveAmbientTexture.node,x3dom.nodeTypes.SurfaceShaderTexture)){this._cf.multiEmissiveAmbientTexture.node._cf.texture.node._type="multiEmissiveAmbientMap";return this._cf.multiEmissiveAmbientTexture.node._cf.texture.node;}else{this._cf.multiEmissiveAmbientTexture.node._type="multiEmissiveAmbientMap";return this._cf.multiEmissiveAmbientTexture.node;}}else{return null;}},getMultiSpecularShininessMap:function()
{if(this._cf.multiSpecularShininessTexture.node){if(x3dom.isa(this._cf.multiSpecularShininessTexture.node,x3dom.nodeTypes.SurfaceShaderTexture)){this._cf.multiSpecularShininessTexture.node._cf.texture.node._type="multiSpecularShininessMap";return this._cf.multiSpecularShininessTexture.node._cf.texture.node;}else{this._cf.multiSpecularShininessTexture.node._type="multiSpecularShininessMap";return this._cf.multiSpecularShininessTexture.node;}}else{return null;}},getMultiVisibilityMap:function()
{if(this._cf.multiVisibilityTexture.node){if(x3dom.isa(this._cf.multiVisibilityTexture.node,x3dom.nodeTypes.SurfaceShaderTexture)){this._cf.multiVisibilityTexture.node._cf.texture.node._type="multiVisibilityMap";return this._cf.multiVisibilityTexture.node._cf.texture.node;}else{this._cf.multiVisibilityTexture.node._type="multiVisibilityMap";return this._cf.multiVisibilityTexture.node;}}else{return null;}},getTextures:function()
{var textures=[];var diff=this.getDiffuseMap();if(diff)textures.push(diff);var norm=this.getNormalMap();if(norm)textures.push(norm);var spec=this.getSpecularMap();if(spec)textures.push(spec);var shin=this.getShininessMap();if(shin)textures.push(shin);var env=this.getEnvironmentMap();if(env)textures.push(env);var displacement=this.getDisplacementMap();if(displacement)textures.push(displacement);var diffuseDisplacement=this.getDiffuseDisplacementMap();if(diffuseDisplacement)textures.push(diffuseDisplacement);var multiDiffuseAlpha=this.getMultiDiffuseAlphaMap();if(multiDiffuseAlpha)textures.push(multiDiffuseAlpha);var multiEmissiveAmbient=this.getMultiEmissiveAmbientMap();if(multiEmissiveAmbient)textures.push(multiEmissiveAmbient);var multiSpecularShininess=this.getMultiSpecularShininessMap();if(multiSpecularShininess)textures.push(multiSpecularShininess);var multiVisibility=this.getMultiVisibilityMap();if(multiVisibility)textures.push(multiVisibility);return textures;},needTexcoords:function()
{return(this.getDiffuseMap()||this.getNormalMap()||this.getSpecularMap()||this.getShininessMap()||this.getDisplacementMap()||this.getDiffuseDisplacementMap()||this.getEnvironmentMap())?true:false;}}));x3dom.registerNodeType("ComposedShader","Shaders",defineClass(x3dom.nodeTypes.X3DShaderNode,function(ctx){x3dom.nodeTypes.ComposedShader.superClass.call(this,ctx);this.addField_MFNode('fields',x3dom.nodeTypes.Field);this.addField_MFNode('parts',x3dom.nodeTypes.ShaderPart);this._vertex=null;this._fragment=null;this._id=null;if(!x3dom.nodeTypes.ComposedShader.ShaderInfoMsgShown){x3dom.debug.logInfo("Current ComposedShader node implementation limitations:\n"+"Vertex attributes (if given in the standard X3D fields 'coord', 'color', "+"'normal', 'texCoord'), matrices and texture are provided as follows...\n"+"(see also <a href='http://x3dom.org/x3dom/doc/help/composedShader.html'>"+"http://x3dom.org/x3dom/doc/help/composedShader.html</a>)\n"+"    attribute vec3 position;\n"+"    attribute vec3 normal;\n"+"    attribute vec2 texcoord;\n"+"    attribute vec3 color;\n"+"    uniform mat4 modelViewProjectionMatrix;\n"+"    uniform mat4 modelViewMatrix;\n"+"    uniform mat4 normalMatrix;\n"+"    uniform mat4 viewMatrix;\n"+"    uniform sampler2D tex;\n");x3dom.nodeTypes.ComposedShader.ShaderInfoMsgShown=true;}},{nodeChanged:function()
{var i,n=this._cf.parts.nodes.length;for(i=0;i<n;i++)
{if(this._cf.parts.nodes[i]._vf.type.toLowerCase()=='vertex'){this._vertex=this._cf.parts.nodes[i];this._id=this._cf.parts.nodes[i]._id;}
else if(this._cf.parts.nodes[i]._vf.type.toLowerCase()=='fragment'){this._fragment=this._cf.parts.nodes[i];this._id+=" - "+this._cf.parts.nodes[i]._id;}}
var ctx={};n=this._cf.fields.nodes.length;for(i=0;i<n;i++)
{var fieldName=this._cf.fields.nodes[i]._vf.name;ctx.xmlNode=this._cf.fields.nodes[i]._xmlNode;var needNode=false;if(ctx.xmlNode===undefined||ctx.xmlNode===null){ctx.xmlNode=document.createElement("field");needNode=true;}
ctx.xmlNode.setAttribute(fieldName,this._cf.fields.nodes[i]._vf.value);var funcName="this.addField_"+this._cf.fields.nodes[i]._vf.type+"(ctx, name);";var func=new Function('ctx','name',funcName);func.call(this,ctx,fieldName);if(needNode){ctx.xmlNode=null;}}
Array.forEach(this._parentNodes,function(app){Array.forEach(app._parentNodes,function(shape){if(shape._cleanupGLObjects)
shape._cleanupGLObjects();shape.setAllDirty();});});},fieldChanged:function(fieldName)
{var i,n=this._cf.fields.nodes.length;for(i=0;i<n;i++)
{var field=this._cf.fields.nodes[i]._vf.name;if(field===fieldName)
{var msg=this._cf.fields.nodes[i]._vf.value;try{this._vf[field].setValueByStr(msg);}
catch(exc1){try{switch((typeof(this._vf[field])).toString()){case"number":this._vf[field]=+msg;break;case"boolean":this._vf[field]=(msg.toLowerCase()==="true");break;case"string":this._vf[field]=msg;break;}}
catch(exc2){x3dom.debug.logError("setValueByStr() NYI for "+typeof(this._vf[field]));}}
break;}}
if(field==='url')
{Array.forEach(this._parentNodes,function(app){Array.forEach(app._parentNodes,function(shape){shape._dirty.shader=true;});});}},parentAdded:function(parent)
{parent.nodeChanged();}}));x3dom.nodeTypes.ComposedShader.ShaderInfoMsgShown=false;x3dom.registerNodeType("ShaderPart","Shaders",defineClass(x3dom.nodeTypes.X3DNode,function(ctx){x3dom.nodeTypes.ShaderPart.superClass.call(this,ctx);this.addField_MFString(ctx,'url',[]);this.addField_SFString(ctx,'type',"VERTEX");this._id=(ctx&&ctx.xmlNode&&ctx.xmlNode.id!="")?ctx.xmlNode.id:++x3dom.nodeTypes.Shape.shaderPartID;x3dom.debug.assert(this._vf.type.toLowerCase()=='vertex'||this._vf.type.toLowerCase()=='fragment',"Unknown shader part type!");},{nodeChanged:function()
{var ctx={};ctx.xmlNode=this._xmlNode;if(ctx.xmlNode!==undefined&&ctx.xmlNode!==null)
{var that=this;if(that._vf.url.length&&that._vf.url[0].indexOf('\n')==-1)
{var xhr=new XMLHttpRequest();xhr.open("GET",that._nameSpace.getURL(that._vf.url[0]),false);xhr.onload=function(){that._vf.url=new x3dom.fields.MFString([]);that._vf.url.push(xhr.response);};xhr.onerror=function(){x3dom.debug.logError("Could not load file '"+that._vf.url[0]+"'.");};x3dom.RequestManager.addRequest(xhr);}
else
{if(that._vf.url.length){that._vf.url=new x3dom.fields.MFString([]);}
try{that._vf.url.push(ctx.xmlNode.childNodes[1].nodeValue);ctx.xmlNode.removeChild(ctx.xmlNode.childNodes[1]);}
catch(e){Array.forEach(ctx.xmlNode.childNodes,function(childDomNode){if(childDomNode.nodeType===3){that._vf.url.push(childDomNode.nodeValue);}
else if(childDomNode.nodeType===4){that._vf.url.push(childDomNode.data);}
childDomNode.parentNode.removeChild(childDomNode);});}}}
Array.forEach(this._parentNodes,function(shader){shader.nodeChanged();});},fieldChanged:function(fieldName)
{if(fieldName==="url"){Array.forEach(this._parentNodes,function(shader){shader.fieldChanged("url");});}},parentAdded:function(parent)
{parent.nodeChanged();}}));x3dom.registerNodeType("X3DVertexAttributeNode","Shaders",defineClass(x3dom.nodeTypes.X3DGeometricPropertyNode,function(ctx){x3dom.nodeTypes.X3DVertexAttributeNode.superClass.call(this,ctx);this.addField_SFString(ctx,'name',"");}));x3dom.registerNodeType("FloatVertexAttribute","Shaders",defineClass(x3dom.nodeTypes.X3DVertexAttributeNode,function(ctx){x3dom.nodeTypes.FloatVertexAttribute.superClass.call(this,ctx);this.addField_SFInt32(ctx,'numComponents',4);this.addField_MFFloat(ctx,'value',[]);}));x3dom.registerNodeType("X3DSpatialGeometryNode","Geometry3D",defineClass(x3dom.nodeTypes.X3DGeometryNode,function(ctx){x3dom.nodeTypes.X3DSpatialGeometryNode.superClass.call(this,ctx);}));x3dom.registerNodeType("Plane","Geometry3D",defineClass(x3dom.nodeTypes.X3DSpatialGeometryNode,function(ctx){x3dom.nodeTypes.Plane.superClass.call(this,ctx);this.addField_SFVec2f(ctx,'size',2,2);this.addField_SFVec2f(ctx,'subdivision',1,1);this.addField_SFVec3f(ctx,'center',0,0,0);this.addField_MFString(ctx,'primType',['TRIANGLES']);if(this._vf.primType.length)
this._mesh._primType=this._vf.primType[0];var sx=this._vf.size.x,sy=this._vf.size.y;var subx=this._vf.subdivision.x,suby=this._vf.subdivision.y;var geoCacheID='Plane_'+sx+'-'+sy+'-'+subx+'-'+suby+'-'+
this._vf.center.x+'-'+this._vf.center.y+'-'+this._vf.center.z;if(ctx&&this._vf.useGeoCache&&x3dom.geoCache[geoCacheID]!==undefined){this._mesh=x3dom.geoCache[geoCacheID];}
else{var x=0,y=0;var xstep=sx/subx;var ystep=sy/suby;sx/=2;sy/=2;for(y=0;y<=suby;y++){for(x=0;x<=subx;x++){this._mesh._positions[0].push(this._vf.center.x+x*xstep-sx);this._mesh._positions[0].push(this._vf.center.y+y*ystep-sy);this._mesh._positions[0].push(this._vf.center.z);this._mesh._normals[0].push(0);this._mesh._normals[0].push(0);this._mesh._normals[0].push(1);this._mesh._texCoords[0].push(x/subx);this._mesh._texCoords[0].push(y/suby);}}
for(y=1;y<=suby;y++){for(x=0;x<subx;x++){this._mesh._indices[0].push((y-1)*(subx+1)+x);this._mesh._indices[0].push((y-1)*(subx+1)+x+1);this._mesh._indices[0].push(y*(subx+1)+x);this._mesh._indices[0].push(y*(subx+1)+x);this._mesh._indices[0].push((y-1)*(subx+1)+x+1);this._mesh._indices[0].push(y*(subx+1)+x+1);}}
this._mesh._invalidate=true;this._mesh._numFaces=this._mesh._indices[0].length/3;this._mesh._numCoords=this._mesh._positions[0].length/3;x3dom.geoCache[geoCacheID]=this._mesh;}},{fieldChanged:function(fieldName){if(fieldName=="size"||fieldName=="center"){this._mesh._positions[0]=[];var sx=this._vf.size.x,sy=this._vf.size.y;var subx=this._vf.subdivision.x,suby=this._vf.subdivision.y;var x=0,y=0;var xstep=sx/subx;var ystep=sy/suby;sx/=2;sy/=2;for(y=0;y<=suby;y++){for(x=0;x<=subx;x++){this._mesh._positions[0].push(this._vf.center.x+x*xstep-sx);this._mesh._positions[0].push(this._vf.center.y+y*ystep-sy);this._mesh._positions[0].push(this._vf.center.z);}}
this.invalidateVolume();this._mesh._numCoords=this._mesh._positions[0].length/3;Array.forEach(this._parentNodes,function(node){node._dirty.positions=true;node.invalidateVolume();});}
else if(fieldName=="subdivision"){this._mesh._positions[0]=[];this._mesh._indices[0]=[];this._mesh._normals[0]=[];this._mesh._texCoords[0]=[];var sx=this._vf.size.x,sy=this._vf.size.y;var subx=this._vf.subdivision.x,suby=this._vf.subdivision.y;var x=0,y=0;var xstep=sx/subx;var ystep=sy/suby;sx/=2;sy/=2;for(y=0;y<=suby;y++){for(x=0;x<=subx;x++){this._mesh._positions[0].push(this._vf.center.x+x*xstep-sx);this._mesh._positions[0].push(this._vf.center.y+y*ystep-sy);this._mesh._positions[0].push(this._vf.center.z);this._mesh._normals[0].push(0);this._mesh._normals[0].push(0);this._mesh._normals[0].push(1);this._mesh._texCoords[0].push(x/subx);this._mesh._texCoords[0].push(y/suby);}}
for(y=1;y<=suby;y++){for(x=0;x<subx;x++){this._mesh._indices[0].push((y-1)*(subx+1)+x);this._mesh._indices[0].push((y-1)*(subx+1)+x+1);this._mesh._indices[0].push(y*(subx+1)+x);this._mesh._indices[0].push(y*(subx+1)+x);this._mesh._indices[0].push((y-1)*(subx+1)+x+1);this._mesh._indices[0].push(y*(subx+1)+x+1);}}
this.invalidateVolume();this._mesh._numFaces=this._mesh._indices[0].length/3;this._mesh._numCoords=this._mesh._positions[0].length/3;Array.forEach(this._parentNodes,function(node){node.setAllDirty();node.invalidateVolume();});}}}));x3dom.registerNodeType("Box","Geometry3D",defineClass(x3dom.nodeTypes.X3DSpatialGeometryNode,function(ctx){x3dom.nodeTypes.Box.superClass.call(this,ctx);this.addField_SFVec3f(ctx,'size',2,2,2);this.addField_SFBool(ctx,'hasHelperColors',false);var sx=this._vf.size.x,sy=this._vf.size.y,sz=this._vf.size.z;var geoCacheID='Box_'+sx+'-'+sy+'-'+sz;if(this._vf.useGeoCache&&x3dom.geoCache[geoCacheID]!==undefined)
{this._mesh=x3dom.geoCache[geoCacheID];}
else
{sx/=2;sy/=2;sz/=2;this._mesh._positions[0]=[-sx,-sy,-sz,-sx,sy,-sz,sx,sy,-sz,sx,-sy,-sz,-sx,-sy,sz,-sx,sy,sz,sx,sy,sz,sx,-sy,sz,-sx,-sy,-sz,-sx,-sy,sz,-sx,sy,sz,-sx,sy,-sz,sx,-sy,-sz,sx,-sy,sz,sx,sy,sz,sx,sy,-sz,-sx,sy,-sz,-sx,sy,sz,sx,sy,sz,sx,sy,-sz,-sx,-sy,-sz,-sx,-sy,sz,sx,-sy,sz,sx,-sy,-sz];this._mesh._normals[0]=[0,0,-1,0,0,-1,0,0,-1,0,0,-1,0,0,1,0,0,1,0,0,1,0,0,1,-1,0,0,-1,0,0,-1,0,0,-1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,0,1,0,0,1,0,0,1,0,0,1,0,0,-1,0,0,-1,0,0,-1,0,0,-1,0];this._mesh._texCoords[0]=[1,0,1,1,0,1,0,0,0,0,0,1,1,1,1,0,0,0,1,0,1,1,0,1,1,0,0,0,0,1,1,1,0,1,0,0,1,0,1,1,0,0,0,1,1,1,1,0];if(this._vf.hasHelperColors){this._mesh._colors[0]=[0,0,0,0,1,0,1,1,0,1,0,0,0,0,1,0,1,1,1,1,1,1,0,1,0,0,0,0,0,1,0,1,1,0,1,0,1,0,0,1,0,1,1,1,1,1,1,0,0,1,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,0,1,1,0,0];}
this._mesh._indices[0]=[0,1,2,2,3,0,4,7,5,5,7,6,8,9,10,10,11,8,12,14,13,14,12,15,16,17,18,18,19,16,20,22,21,22,20,23];this._mesh._invalidate=true;this._mesh._numFaces=12;this._mesh._numCoords=24;x3dom.geoCache[geoCacheID]=this._mesh;}},{fieldChanged:function(fieldName)
{if(fieldName==="size"){var sx=this._vf.size.x/2,sy=this._vf.size.y/2,sz=this._vf.size.z/2;this._mesh._positions[0]=[-sx,-sy,-sz,-sx,sy,-sz,sx,sy,-sz,sx,-sy,-sz,-sx,-sy,sz,-sx,sy,sz,sx,sy,sz,sx,-sy,sz,-sx,-sy,-sz,-sx,-sy,sz,-sx,sy,sz,-sx,sy,-sz,sx,-sy,-sz,sx,-sy,sz,sx,sy,sz,sx,sy,-sz,-sx,sy,-sz,-sx,sy,sz,sx,sy,sz,sx,sy,-sz,-sx,-sy,-sz,-sx,-sy,sz,sx,-sy,sz,sx,-sy,-sz];this.invalidateVolume();Array.forEach(this._parentNodes,function(node){node._dirty.positions=true;node.invalidateVolume();});}
else if(fieldName==="hasHelperColors"){if(this._vf.hasHelperColors){this._mesh._colors[0]=[0,0,0,0,1,0,1,1,0,1,0,0,0,0,1,0,1,1,1,1,1,1,0,1,0,0,0,0,0,1,0,1,1,0,1,0,1,0,0,1,0,1,1,1,1,1,1,0,0,1,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,0,1,1,0,0];}
else{this._mesh._colors[0]=[];}
Array.forEach(this._parentNodes,function(node){node._dirty.colors=true;});}}}));x3dom.registerNodeType("Sphere","Geometry3D",defineClass(x3dom.nodeTypes.X3DSpatialGeometryNode,function(ctx){x3dom.nodeTypes.Sphere.superClass.call(this,ctx);this.addField_SFFloat(ctx,'radius',ctx?1:10000);this.addField_SFVec2f(ctx,'subdivision',24,24);var qfactor=1.0;var r=this._vf.radius;var subx=this._vf.subdivision.x,suby=this._vf.subdivision.y;var geoCacheID='Sphere_'+r+'-'+subx+'-'+suby;if(this._vf.useGeoCache&&x3dom.geoCache[geoCacheID]!==undefined){this._mesh=x3dom.geoCache[geoCacheID];}
else{if(ctx){qfactor=ctx.doc.properties.getProperty("PrimitiveQuality","Medium");}
if(!x3dom.Utils.isNumber(qfactor)){switch(qfactor.toLowerCase()){case"low":qfactor=0.3;break;case"medium":qfactor=0.5;break;case"high":qfactor=1.0;break;}}else{qfactor=parseFloat(qfactor);}
this._quality=qfactor;var latNumber,longNumber;var latitudeBands=Math.floor(subx*qfactor);var longitudeBands=Math.floor(suby*qfactor);var theta,sinTheta,cosTheta;var phi,sinPhi,cosPhi;var x,y,z,u,v;for(latNumber=0;latNumber<=latitudeBands;latNumber++){theta=(latNumber*Math.PI)/latitudeBands;sinTheta=Math.sin(theta);cosTheta=Math.cos(theta);for(longNumber=0;longNumber<=longitudeBands;longNumber++){phi=(longNumber*2.0*Math.PI)/longitudeBands;sinPhi=Math.sin(phi);cosPhi=Math.cos(phi);x=-cosPhi*sinTheta;y=-cosTheta;z=-sinPhi*sinTheta;u=0.25-(longNumber/longitudeBands);v=latNumber/latitudeBands;this._mesh._positions[0].push(r*x);this._mesh._positions[0].push(r*y);this._mesh._positions[0].push(r*z);this._mesh._normals[0].push(x);this._mesh._normals[0].push(y);this._mesh._normals[0].push(z);this._mesh._texCoords[0].push(u);this._mesh._texCoords[0].push(v);}}
var first,second;for(latNumber=0;latNumber<latitudeBands;latNumber++){for(longNumber=0;longNumber<longitudeBands;longNumber++){first=(latNumber*(longitudeBands+1))+longNumber;second=first+longitudeBands+1;this._mesh._indices[0].push(first);this._mesh._indices[0].push(second);this._mesh._indices[0].push(first+1);this._mesh._indices[0].push(second);this._mesh._indices[0].push(second+1);this._mesh._indices[0].push(first+1);}}
this._mesh._invalidate=true;this._mesh._numFaces=this._mesh._indices[0].length/3;this._mesh._numCoords=this._mesh._positions[0].length/3;x3dom.geoCache[geoCacheID]=this._mesh;}},{fieldChanged:function(fieldName){if(fieldName==="radius"){this._mesh._positions[0]=[];var r=this._vf.radius;var subx=this._vf.subdivision.x,suby=this._vf.subdivision.y;var qfactor=this._quality;var latNumber,longNumber;var latitudeBands=Math.floor(subx*qfactor);var longitudeBands=Math.floor(suby*qfactor);var theta,sinTheta,cosTheta;var phi,sinPhi,cosPhi;var x,y,z;for(latNumber=0;latNumber<=latitudeBands;latNumber++){theta=(latNumber*Math.PI)/latitudeBands;sinTheta=Math.sin(theta);cosTheta=Math.cos(theta);for(longNumber=0;longNumber<=longitudeBands;longNumber++){phi=(longNumber*2.0*Math.PI)/longitudeBands;sinPhi=Math.sin(phi);cosPhi=Math.cos(phi);x=-cosPhi*sinTheta;y=-cosTheta;z=-sinPhi*sinTheta;this._mesh._positions[0].push(r*x);this._mesh._positions[0].push(r*y);this._mesh._positions[0].push(r*z);}}
this.invalidateVolume();this._mesh._numCoords=this._mesh._positions[0].length/3;Array.forEach(this._parentNodes,function(node){node._dirty.positions=true;node.invalidateVolume();});}
else if(fieldName==="subdivision"){this._mesh._positions[0]=[];this._mesh._indices[0]=[];this._mesh._normals[0]=[];this._mesh._texCoords[0]=[];var r=this._vf.radius;var subx=this._vf.subdivision.x,suby=this._vf.subdivision.y;var qfactor=this._quality;var latNumber,longNumber;var latitudeBands=Math.floor(subx*qfactor);var longitudeBands=Math.floor(suby*qfactor);var theta,sinTheta,cosTheta;var phi,sinPhi,cosPhi;var x,y,z,u,v;for(latNumber=0;latNumber<=latitudeBands;latNumber++){theta=(latNumber*Math.PI)/latitudeBands;sinTheta=Math.sin(theta);cosTheta=Math.cos(theta);for(longNumber=0;longNumber<=longitudeBands;longNumber++){phi=(longNumber*2.0*Math.PI)/longitudeBands;sinPhi=Math.sin(phi);cosPhi=Math.cos(phi);x=-cosPhi*sinTheta;y=-cosTheta;z=-sinPhi*sinTheta;u=0.25-(longNumber/longitudeBands);v=latNumber/latitudeBands;this._mesh._positions[0].push(r*x);this._mesh._positions[0].push(r*y);this._mesh._positions[0].push(r*z);this._mesh._normals[0].push(x);this._mesh._normals[0].push(y);this._mesh._normals[0].push(z);this._mesh._texCoords[0].push(u);this._mesh._texCoords[0].push(v);}}
var first,second;for(latNumber=0;latNumber<latitudeBands;latNumber++){for(longNumber=0;longNumber<longitudeBands;longNumber++){first=(latNumber*(longitudeBands+1))+longNumber;second=first+longitudeBands+1;this._mesh._indices[0].push(first);this._mesh._indices[0].push(second);this._mesh._indices[0].push(first+1);this._mesh._indices[0].push(second);this._mesh._indices[0].push(second+1);this._mesh._indices[0].push(first+1);}}
this.invalidateVolume();this._mesh._numFaces=this._mesh._indices[0].length/3;this._mesh._numCoords=this._mesh._positions[0].length/3;Array.forEach(this._parentNodes,function(node){node.setAllDirty();node.invalidateVolume();});}}}));x3dom.registerNodeType("Torus","Geometry3D",defineClass(x3dom.nodeTypes.X3DSpatialGeometryNode,function(ctx){x3dom.nodeTypes.Torus.superClass.call(this,ctx);var twoPi=2.0*Math.PI;this.addField_SFFloat(ctx,'innerRadius',0.5);this.addField_SFFloat(ctx,'outerRadius',1.0);this.addField_SFFloat(ctx,'angle',twoPi);this.addField_SFBool(ctx,'caps',true);this.addField_SFVec2f(ctx,'subdivision',24,24);this.addField_SFBool(ctx,'insideOutsideRadius',false);if(this._vf.angle<0)
this._vf.angle=0;else if(this._vf.angle>twoPi)
this._vf.angle=twoPi;this._origCCW=this._vf.ccw;var innerRadius=this._vf.innerRadius;var outerRadius=this._vf.outerRadius;if(this._vf.insideOutsideRadius==true)
{if(innerRadius>outerRadius){var tmp=innerRadius;innerRadius=outerRadius;outerRadius=tmp;}
var rad=(outerRadius-innerRadius)/2;outerRadius=innerRadius+rad;innerRadius=rad;this._vf.ccw=!this._origCCW;}
var rings=this._vf.subdivision.x,sides=this._vf.subdivision.y;rings=Math.max(3,Math.round((this._vf.angle/twoPi)*rings));var geoCacheID='Torus_'+innerRadius+'_'+outerRadius+'_'+this._vf.angle+'_'+
this._vf.subdivision+'-'+this._vf.caps;if(this._vf.useGeoCache&&x3dom.geoCache[geoCacheID]!==undefined)
{this._mesh=x3dom.geoCache[geoCacheID];}
else
{var ringDelta=this._vf.angle/rings;var sideDelta=twoPi/sides;var a,b,theta,phi;var cosTheta,sinTheta,cosPhi,sinPhi,dist;for(a=0,theta=0;a<=rings;a++,theta+=ringDelta)
{cosTheta=Math.cos(theta);sinTheta=Math.sin(theta);for(b=0,phi=0;b<=sides;b++,phi+=sideDelta)
{cosPhi=Math.cos(phi);sinPhi=Math.sin(phi);dist=outerRadius+innerRadius*cosPhi;if(this._vf.insideOutsideRadius){this._mesh._positions[0].push(cosTheta*dist,innerRadius*sinPhi,-sinTheta*dist);this._mesh._normals[0].push(cosTheta*cosPhi,sinPhi,-sinTheta*cosPhi);}
else{this._mesh._positions[0].push(cosTheta*dist,-sinTheta*dist,innerRadius*sinPhi);this._mesh._normals[0].push(cosTheta*cosPhi,-sinTheta*cosPhi,sinPhi);}
this._mesh._texCoords[0].push(-a/rings,b/sides);}}
for(a=0;a<sides;a++)
{for(b=0;b<rings;b++)
{this._mesh._indices[0].push(b*(sides+1)+a);this._mesh._indices[0].push(b*(sides+1)+a+1);this._mesh._indices[0].push((b+1)*(sides+1)+a);this._mesh._indices[0].push(b*(sides+1)+a+1);this._mesh._indices[0].push((b+1)*(sides+1)+a+1);this._mesh._indices[0].push((b+1)*(sides+1)+a);}}
if(this._vf.angle<twoPi&&this._vf.caps==true)
{var origPos=this._mesh._positions[0].length/3;if(this._vf.insideOutsideRadius){this._mesh._positions[0].push(outerRadius,0,0);this._mesh._normals[0].push(0,0,1);}
else{this._mesh._positions[0].push(outerRadius,0,0);this._mesh._normals[0].push(0,1,0);}
this._mesh._texCoords[0].push(0.5,0.5);for(b=0,phi=0;b<=sides;b++,phi+=sideDelta)
{cosPhi=Math.cos(phi);sinPhi=Math.sin(phi);dist=outerRadius+innerRadius*cosPhi;if(this._vf.insideOutsideRadius){this._mesh._positions[0].push(dist,sinPhi*innerRadius,0);this._mesh._normals[0].push(0,0,1);}
else{this._mesh._positions[0].push(dist,0,sinPhi*innerRadius);this._mesh._normals[0].push(0,1,0);}
this._mesh._texCoords[0].push((1+cosPhi)*0.5,(1-sinPhi)*0.5);if(b>0){this._mesh._indices[0].push(origPos);this._mesh._indices[0].push(origPos+b);this._mesh._indices[0].push(origPos+b-1);}
if(b==sides){this._mesh._indices[0].push(origPos);this._mesh._indices[0].push(origPos+1);this._mesh._indices[0].push(origPos+b);}}
cosTheta=Math.cos(this._vf.angle);sinTheta=Math.sin(this._vf.angle);origPos=this._mesh._positions[0].length/3;var nx=-sinTheta,ny=-cosTheta;if(this._vf.insideOutsideRadius){this._mesh._positions[0].push(cosTheta*outerRadius,0,-sinTheta*outerRadius);this._mesh._normals[0].push(nx,0,ny);}
else{this._mesh._positions[0].push(cosTheta*outerRadius,-sinTheta*outerRadius,0);this._mesh._normals[0].push(nx,ny,0);}
this._mesh._texCoords[0].push(0.5,0.5);for(b=0,phi=0;b<=sides;b++,phi+=sideDelta)
{cosPhi=Math.cos(phi);sinPhi=Math.sin(phi);dist=outerRadius+innerRadius*cosPhi;if(this._vf.insideOutsideRadius){this._mesh._positions[0].push(cosTheta*dist,sinPhi*innerRadius,-sinTheta*dist);this._mesh._normals[0].push(nx,0,ny);}
else{this._mesh._positions[0].push(cosTheta*dist,-sinTheta*dist,sinPhi*innerRadius);this._mesh._normals[0].push(nx,ny,0);}
this._mesh._texCoords[0].push(1-(1+cosPhi)*0.5,(1-sinPhi)*0.5);if(b>0){this._mesh._indices[0].push(origPos);this._mesh._indices[0].push(origPos+b-1);this._mesh._indices[0].push(origPos+b);}
if(b==sides){this._mesh._indices[0].push(origPos);this._mesh._indices[0].push(origPos+b);this._mesh._indices[0].push(origPos+1);}}}
this._mesh._invalidate=true;this._mesh._numFaces=this._mesh._indices[0].length/3;this._mesh._numCoords=this._mesh._positions[0].length/3;x3dom.geoCache[geoCacheID]=this._mesh;}},{fieldChanged:function(fieldName)
{if(fieldName=="innerRadius"||fieldName=="outerRadius"||fieldName=="subdivision"||fieldName=="angle"||fieldName=="insideOutsideRadius"||fieldName=="caps")
{var twoPi=2.0*Math.PI;if(this._vf.angle<0)
this._vf.angle=0;else if(this._vf.angle>twoPi)
this._vf.angle=twoPi;var innerRadius=this._vf.innerRadius;var outerRadius=this._vf.outerRadius;if(this._vf.insideOutsideRadius==true)
{if(innerRadius>outerRadius){var tmp=innerRadius;innerRadius=outerRadius;outerRadius=tmp;}
var rad=(outerRadius-innerRadius)/2;outerRadius=innerRadius+rad;innerRadius=rad;this._vf.ccw=!this._origCCW;}
else
this._vf.ccw=this._origCCW;var rings=this._vf.subdivision.x,sides=this._vf.subdivision.y;rings=Math.max(3,Math.round((this._vf.angle/twoPi)*rings));var ringDelta=this._vf.angle/rings;var sideDelta=twoPi/sides;var a,b,theta,phi;var cosTheta,sinTheta,cosPhi,sinPhi,dist;this._mesh._positions[0]=[];this._mesh._normals[0]=[];this._mesh._texCoords[0]=[];this._mesh._indices[0]=[];for(a=0,theta=0;a<=rings;a++,theta+=ringDelta)
{cosTheta=Math.cos(theta);sinTheta=Math.sin(theta);for(b=0,phi=0;b<=sides;b++,phi+=sideDelta)
{cosPhi=Math.cos(phi);sinPhi=Math.sin(phi);dist=outerRadius+innerRadius*cosPhi;if(this._vf.insideOutsideRadius){this._mesh._positions[0].push(cosTheta*dist,innerRadius*sinPhi,-sinTheta*dist);this._mesh._normals[0].push(cosTheta*cosPhi,sinPhi,-sinTheta*cosPhi);}
else{this._mesh._positions[0].push(cosTheta*dist,-sinTheta*dist,innerRadius*sinPhi);this._mesh._normals[0].push(cosTheta*cosPhi,-sinTheta*cosPhi,sinPhi);}
this._mesh._texCoords[0].push(-a/rings,b/sides);}}
for(a=0;a<sides;a++)
{for(b=0;b<rings;b++)
{this._mesh._indices[0].push(b*(sides+1)+a);this._mesh._indices[0].push(b*(sides+1)+a+1);this._mesh._indices[0].push((b+1)*(sides+1)+a);this._mesh._indices[0].push(b*(sides+1)+a+1);this._mesh._indices[0].push((b+1)*(sides+1)+a+1);this._mesh._indices[0].push((b+1)*(sides+1)+a);}}
if(this._vf.angle<twoPi&&this._vf.caps==true)
{var origPos=this._mesh._positions[0].length/3;if(this._vf.insideOutsideRadius){this._mesh._positions[0].push(outerRadius,0,0);this._mesh._normals[0].push(0,0,1);}
else{this._mesh._positions[0].push(outerRadius,0,0);this._mesh._normals[0].push(0,1,0);}
this._mesh._texCoords[0].push(0.5,0.5);for(b=0,phi=0;b<=sides;b++,phi+=sideDelta)
{cosPhi=Math.cos(phi);sinPhi=Math.sin(phi);dist=outerRadius+innerRadius*cosPhi;if(this._vf.insideOutsideRadius){this._mesh._positions[0].push(dist,sinPhi*innerRadius,0);this._mesh._normals[0].push(0,0,1);}
else{this._mesh._positions[0].push(dist,0,sinPhi*innerRadius);this._mesh._normals[0].push(0,1,0);}
this._mesh._texCoords[0].push((1+cosPhi)*0.5,(1-sinPhi)*0.5);if(b>0){this._mesh._indices[0].push(origPos);this._mesh._indices[0].push(origPos+b);this._mesh._indices[0].push(origPos+b-1);}
if(b==sides){this._mesh._indices[0].push(origPos);this._mesh._indices[0].push(origPos+1);this._mesh._indices[0].push(origPos+b);}}
cosTheta=Math.cos(this._vf.angle);sinTheta=Math.sin(this._vf.angle);origPos=this._mesh._positions[0].length/3;var nx=-sinTheta,ny=-cosTheta;if(this._vf.insideOutsideRadius){this._mesh._positions[0].push(cosTheta*outerRadius,0,-sinTheta*outerRadius);this._mesh._normals[0].push(nx,0,ny);}
else{this._mesh._positions[0].push(cosTheta*outerRadius,-sinTheta*outerRadius,0);this._mesh._normals[0].push(nx,ny,0);}
this._mesh._texCoords[0].push(0.5,0.5);for(b=0,phi=0;b<=sides;b++,phi+=sideDelta)
{cosPhi=Math.cos(phi);sinPhi=Math.sin(phi);dist=outerRadius+innerRadius*cosPhi;if(this._vf.insideOutsideRadius){this._mesh._positions[0].push(cosTheta*dist,sinPhi*innerRadius,-sinTheta*dist);this._mesh._normals[0].push(nx,0,ny);}
else{this._mesh._positions[0].push(cosTheta*dist,-sinTheta*dist,sinPhi*innerRadius);this._mesh._normals[0].push(nx,ny,0);}
this._mesh._texCoords[0].push(1-(1+cosPhi)*0.5,(1-sinPhi)*0.5);if(b>0){this._mesh._indices[0].push(origPos);this._mesh._indices[0].push(origPos+b-1);this._mesh._indices[0].push(origPos+b);}
if(b==sides){this._mesh._indices[0].push(origPos);this._mesh._indices[0].push(origPos+b);this._mesh._indices[0].push(origPos+1);}}}
this.invalidateVolume();this._mesh._numFaces=this._mesh._indices[0].length/3;this._mesh._numCoords=this._mesh._positions[0].length/3;Array.forEach(this._parentNodes,function(node){node.setAllDirty();node.invalidateVolume();});}}}));x3dom.registerNodeType("Cone","Geometry3D",defineClass(x3dom.nodeTypes.X3DSpatialGeometryNode,function(ctx){x3dom.nodeTypes.Cone.superClass.call(this,ctx);this.addField_SFFloat(ctx,'bottomRadius',1.0);this.addField_SFFloat(ctx,'topRadius',0);this.addField_SFFloat(ctx,'height',2.0);this.addField_SFBool(ctx,'bottom',true);this.addField_SFBool(ctx,'side',true);this.addField_SFBool(ctx,'top',true);this.addField_SFFloat(ctx,'subdivision',32);var geoCacheID='Cone_'+this._vf.bottomRadius+'_'+this._vf.height+'_'+this._vf.top+'_'+
this._vf.bottom+'_'+this._vf.side+'_'+this._vf.topRadius+'_'+this._vf.subdivision;if(this._vf.useGeoCache&&x3dom.geoCache[geoCacheID]!==undefined){this._mesh=x3dom.geoCache[geoCacheID];}
else{var bottomRadius=this._vf.bottomRadius,height=this._vf.height;var topRadius=this._vf.topRadius,sides=this._vf.subdivision;var beta,x,z;var delta=2.0*Math.PI/sides;var incl=(bottomRadius-topRadius)/height;var nlen=1.0/Math.sqrt(1.0+incl*incl);var j=0,k=0;var h,base;if(this._vf.side&&height>0){var px=0,pz=0;for(j=0,k=0;j<=sides;j++){beta=j*delta;x=Math.sin(beta);z=-Math.cos(beta);if(topRadius>x3dom.fields.Eps){px=x*topRadius;pz=z*topRadius;}
this._mesh._positions[0].push(px,height/2,pz);this._mesh._normals[0].push(x/nlen,incl/nlen,z/nlen);this._mesh._texCoords[0].push(1.0-j/sides,1);this._mesh._positions[0].push(x*bottomRadius,-height/2,z*bottomRadius);this._mesh._normals[0].push(x/nlen,incl/nlen,z/nlen);this._mesh._texCoords[0].push(1.0-j/sides,0);if(j>0){this._mesh._indices[0].push(k);this._mesh._indices[0].push(k+2);this._mesh._indices[0].push(k+1);this._mesh._indices[0].push(k+1);this._mesh._indices[0].push(k+2);this._mesh._indices[0].push(k+3);k+=2;}}}
if(this._vf.bottom&&bottomRadius>0){base=this._mesh._positions[0].length/3;for(j=sides-1;j>=0;j--){beta=j*delta;x=bottomRadius*Math.sin(beta);z=-bottomRadius*Math.cos(beta);this._mesh._positions[0].push(x,-height/2,z);this._mesh._normals[0].push(0,-1,0);this._mesh._texCoords[0].push(x/bottomRadius/2+0.5,z/bottomRadius/2+0.5);}
h=base+1;for(j=2;j<sides;j++){this._mesh._indices[0].push(h);this._mesh._indices[0].push(base);h=base+j;this._mesh._indices[0].push(h);}}
if(this._vf.top&&topRadius>x3dom.fields.Eps){base=this._mesh._positions[0].length/3;for(j=sides-1;j>=0;j--){beta=j*delta;x=topRadius*Math.sin(beta);z=-topRadius*Math.cos(beta);this._mesh._positions[0].push(x,height/2,z);this._mesh._normals[0].push(0,1,0);this._mesh._texCoords[0].push(x/topRadius/2+0.5,1.0-z/topRadius/2+0.5);}
h=base+1;for(j=2;j<sides;j++){this._mesh._indices[0].push(base);this._mesh._indices[0].push(h);h=base+j;this._mesh._indices[0].push(h);}}
this._mesh._invalidate=true;this._mesh._numFaces=this._mesh._indices[0].length/3;this._mesh._numCoords=this._mesh._positions[0].length/3;x3dom.geoCache[geoCacheID]=this._mesh;}},{fieldChanged:function(fieldName)
{if(fieldName=="bottomRadius"||fieldName=="topRadius"||fieldName=="height"||fieldName=="subdivision"||fieldName=="bottom"||fieldName=="top"||fieldName=="side")
{this._mesh._positions[0]=[];this._mesh._indices[0]=[];this._mesh._normals[0]=[];this._mesh._texCoords[0]=[];var bottomRadius=this._vf.bottomRadius,height=this._vf.height;var topRadius=this._vf.topRadius,sides=this._vf.subdivision;var beta,x,z;var delta=2.0*Math.PI/sides;var incl=(bottomRadius-topRadius)/height;var nlen=1.0/Math.sqrt(1.0+incl*incl);var j=0,k=0;var h,base;if(this._vf.side&&height>0)
{var px=0,pz=0;for(j=0,k=0;j<=sides;j++){beta=j*delta;x=Math.sin(beta);z=-Math.cos(beta);if(topRadius>x3dom.fields.Eps){px=x*topRadius;pz=z*topRadius;}
this._mesh._positions[0].push(px,height/2,pz);this._mesh._normals[0].push(x/nlen,incl/nlen,z/nlen);this._mesh._texCoords[0].push(1.0-j/sides,1);this._mesh._positions[0].push(x*bottomRadius,-height/2,z*bottomRadius);this._mesh._normals[0].push(x/nlen,incl/nlen,z/nlen);this._mesh._texCoords[0].push(1.0-j/sides,0);if(j>0){this._mesh._indices[0].push(k);this._mesh._indices[0].push(k+2);this._mesh._indices[0].push(k+1);this._mesh._indices[0].push(k+1);this._mesh._indices[0].push(k+2);this._mesh._indices[0].push(k+3);k+=2;}}}
if(this._vf.bottom&&bottomRadius>0)
{base=this._mesh._positions[0].length/3;for(j=sides-1;j>=0;j--){beta=j*delta;x=bottomRadius*Math.sin(beta);z=-bottomRadius*Math.cos(beta);this._mesh._positions[0].push(x,-height/2,z);this._mesh._normals[0].push(0,-1,0);this._mesh._texCoords[0].push(x/bottomRadius/2+0.5,z/bottomRadius/2+0.5);}
h=base+1;for(j=2;j<sides;j++){this._mesh._indices[0].push(h);this._mesh._indices[0].push(base);h=base+j;this._mesh._indices[0].push(h);}}
if(this._vf.top&&topRadius>x3dom.fields.Eps)
{base=this._mesh._positions[0].length/3;for(j=sides-1;j>=0;j--){beta=j*delta;x=topRadius*Math.sin(beta);z=-topRadius*Math.cos(beta);this._mesh._positions[0].push(x,height/2,z);this._mesh._normals[0].push(0,1,0);this._mesh._texCoords[0].push(x/topRadius/2+0.5,1.0-z/topRadius/2+0.5);}
h=base+1;for(j=2;j<sides;j++){this._mesh._indices[0].push(base);this._mesh._indices[0].push(h);h=base+j;this._mesh._indices[0].push(h);}}
this.invalidateVolume();this._mesh._numFaces=this._mesh._indices[0].length/3;this._mesh._numCoords=this._mesh._positions[0].length/3;Array.forEach(this._parentNodes,function(node){node.setAllDirty();node.invalidateVolume();});}}}));x3dom.registerNodeType("Cylinder","Geometry3D",defineClass(x3dom.nodeTypes.X3DSpatialGeometryNode,function(ctx){x3dom.nodeTypes.Cylinder.superClass.call(this,ctx);this.addField_SFFloat(ctx,'radius',1.0);this.addField_SFFloat(ctx,'height',2.0);this.addField_SFBool(ctx,'bottom',true);this.addField_SFBool(ctx,'top',true);this.addField_SFFloat(ctx,'subdivision',32);this.addField_SFBool(ctx,'side',true);var sides=this._vf.subdivision;var geoCacheID='Cylinder_'+this._vf.radius+'_'+this._vf.height+'_'+this._vf.bottom+'_'+this._vf.top+'_'+
this._vf.side+'_'+this._vf.subdivision;if(this._vf.useGeoCache&&x3dom.geoCache[geoCacheID]!==undefined)
{this._mesh=x3dom.geoCache[geoCacheID];}
else
{var radius=this._vf.radius;var height=this._vf.height/2;var beta,x,z;var delta=2.0*Math.PI/sides;var j,k;if(this._vf.side)
{for(j=0,k=0;j<=sides;j++)
{beta=j*delta;x=Math.sin(beta);z=-Math.cos(beta);this._mesh._positions[0].push(x*radius,-height,z*radius);this._mesh._normals[0].push(x,0,z);this._mesh._texCoords[0].push(1.0-j/sides,0);this._mesh._positions[0].push(x*radius,height,z*radius);this._mesh._normals[0].push(x,0,z);this._mesh._texCoords[0].push(1.0-j/sides,1);if(j>0)
{this._mesh._indices[0].push(k);this._mesh._indices[0].push(k+1);this._mesh._indices[0].push(k+2);this._mesh._indices[0].push(k+2);this._mesh._indices[0].push(k+1);this._mesh._indices[0].push(k+3);k+=2;}}}
if(radius>0)
{var h,base=this._mesh._positions[0].length/3;if(this._vf.top)
{for(j=sides-1;j>=0;j--)
{beta=j*delta;x=radius*Math.sin(beta);z=-radius*Math.cos(beta);this._mesh._positions[0].push(x,height,z);this._mesh._normals[0].push(0,1,0);this._mesh._texCoords[0].push(x/radius/2+0.5,-z/radius/2+0.5);}
h=base+1;for(j=2;j<sides;j++)
{this._mesh._indices[0].push(base);this._mesh._indices[0].push(h);h=base+j;this._mesh._indices[0].push(h);}
base=this._mesh._positions[0].length/3;}
if(this._vf.bottom)
{for(j=sides-1;j>=0;j--)
{beta=j*delta;x=radius*Math.sin(beta);z=-radius*Math.cos(beta);this._mesh._positions[0].push(x,-height,z);this._mesh._normals[0].push(0,-1,0);this._mesh._texCoords[0].push(x/radius/2+0.5,z/radius/2+0.5);}
h=base+1;for(j=2;j<sides;j++)
{this._mesh._indices[0].push(h);this._mesh._indices[0].push(base);h=base+j;this._mesh._indices[0].push(h);}}}
this._mesh._invalidate=true;this._mesh._numFaces=this._mesh._indices[0].length/3;this._mesh._numCoords=this._mesh._positions[0].length/3;x3dom.geoCache[geoCacheID]=this._mesh;}},{fieldChanged:function(fieldName){if(fieldName==="radius"||fieldName==="height")
{this._mesh._positions[0]=[];var radius=this._vf.radius,height=this._vf.height/2;var sides=this._vf.subdivision;var beta,x,z,j;var delta=2.0*Math.PI/sides;if(this._vf.side)
{for(j=0;j<=sides;j++)
{beta=j*delta;x=Math.sin(beta);z=-Math.cos(beta);this._mesh._positions[0].push(x*radius,-height,z*radius);this._mesh._positions[0].push(x*radius,height,z*radius);}}
if(radius>0)
{var h,base=this._mesh._positions[0].length/3;if(this._vf.top)
{for(j=sides-1;j>=0;j--)
{beta=j*delta;x=radius*Math.sin(beta);z=-radius*Math.cos(beta);this._mesh._positions[0].push(x,height,z);}}}
if(this._vf.bottom)
{for(j=sides-1;j>=0;j--)
{beta=j*delta;x=radius*Math.sin(beta);z=-radius*Math.cos(beta);this._mesh._positions[0].push(x,-height,z);}}
this.invalidateVolume();this._mesh._numCoords=this._mesh._positions[0].length/3;Array.forEach(this._parentNodes,function(node){node._dirty.positions=true;node.invalidateVolume();});}
else if(fieldName==="subdivision"||fieldName==="bottom"||fieldName==="top"||fieldName==="side")
{this._mesh._positions[0]=[];this._mesh._indices[0]=[];this._mesh._normals[0]=[];this._mesh._texCoords[0]=[];var radius=this._vf.radius,height=this._vf.height/2;var sides=this._vf.subdivision;var beta,x,z,j;var delta=2.0*Math.PI/sides;var k=0;if(this._vf.side)
{for(j=0,k=0;j<=sides;j++)
{beta=j*delta;x=Math.sin(beta);z=-Math.cos(beta);this._mesh._positions[0].push(x*radius,-height,z*radius);this._mesh._normals[0].push(x,0,z);this._mesh._texCoords[0].push(1.0-j/sides,0);this._mesh._positions[0].push(x*radius,height,z*radius);this._mesh._normals[0].push(x,0,z);this._mesh._texCoords[0].push(1.0-j/sides,1);if(j>0)
{this._mesh._indices[0].push(k+0);this._mesh._indices[0].push(k+1);this._mesh._indices[0].push(k+2);this._mesh._indices[0].push(k+2);this._mesh._indices[0].push(k+1);this._mesh._indices[0].push(k+3);k+=2;}}}
if(radius>0)
{var h,base=this._mesh._positions[0].length/3;if(this._vf.top)
{for(j=sides-1;j>=0;j--)
{beta=j*delta;x=radius*Math.sin(beta);z=-radius*Math.cos(beta);this._mesh._positions[0].push(x,height,z);this._mesh._normals[0].push(0,1,0);this._mesh._texCoords[0].push(x/radius/2+0.5,-z/radius/2+0.5);}
h=base+1;for(j=2;j<sides;j++)
{this._mesh._indices[0].push(base);this._mesh._indices[0].push(h);h=base+j;this._mesh._indices[0].push(h);}
base=this._mesh._positions[0].length/3;}
if(this._vf.bottom)
{for(j=sides-1;j>=0;j--)
{beta=j*delta;x=radius*Math.sin(beta);z=-radius*Math.cos(beta);this._mesh._positions[0].push(x,-height,z);this._mesh._normals[0].push(0,-1,0);this._mesh._texCoords[0].push(x/radius/2+0.5,z/radius/2+0.5);}
h=base+1;for(j=2;j<sides;j++)
{this._mesh._indices[0].push(h);this._mesh._indices[0].push(base);h=base+j;this._mesh._indices[0].push(h);}}}
this.invalidateVolume();this._mesh._numFaces=this._mesh._indices[0].length/3;this._mesh._numCoords=this._mesh._positions[0].length/3;Array.forEach(this._parentNodes,function(node){node.setAllDirty();node.invalidateVolume();});}}}));x3dom.registerNodeType("ExternalGeometry","Geometry3D",defineClass(x3dom.nodeTypes.X3DSpatialGeometryNode,function(ctx){x3dom.nodeTypes.ExternalGeometry.superClass.call(this,ctx);this.addField_MFString(ctx,'url',[]);this._mesh._invalidate=false;this._mesh._numCoords=0;this._mesh._numFaces=0;this._currentURLIdx=0;},{update:function(shape,shaderProgram,gl,viewarea,context){var that=this;var xhr;if(this._vf['url'].length==0||this._currentURLIdx>=this._vf['url'].length)
{return;}
if(x3dom.BinaryContainerLoader.outOfMemory){return;}
shape._webgl.internalDownloadCount=1;shape._nameSpace.doc.downloadCount=1;xhr=new XMLHttpRequest();xhr.open("GET",shape._nameSpace.getURL(this._vf['url'][this._currentURLIdx]),true);xhr.responseType="arraybuffer";x3dom.RequestManager.addRequest(xhr);xhr.onerror=function(){x3dom.debug.logError("Unable to load SRC data from URL \""+that._vf['url'][that._currentURLIdx]+"\"");};xhr.onload=function(){shape._webgl.internalDownloadCount=0;shape._nameSpace.doc.downloadCount=0;shape._webgl.primType=[];shape._webgl.indexOffset=[];shape._webgl.drawCount=[];if((xhr.status==200||xhr.status==0)){var glTF=new x3dom.glTF.glTFLoader(xhr.response,true);if(glTF.header.sceneLength>0)
{glTF.loaded={};glTF.loaded.meshes={};glTF.loaded.meshCount=0;var url=that._vf['url'][that._currentURLIdx];if(url.includes('#'))
{var split=url.split('#');var meshName=split[split.length-1];glTF.getMesh(shape,shaderProgram,gl,meshName);}
else
{glTF.getScene(shape,shaderProgram,gl);}
for(var key in glTF._mesh){if(!glTF._mesh.hasOwnProperty(key))continue;that._mesh[key]=glTF._mesh[key];}}
else
{if((that._currentURLIdx+1)<that._vf['url'].length)
{x3dom.debug.logWarning("Invalid SRC data, loaded from URL \""+
that._vf['url'][that._currentURLIdx]+"\", trying next specified URL");++that._currentURLIdx;that.update(shape,shaderProgram,gl,viewarea,context);}
else
{x3dom.debug.logError("Invalid SRC data, loaded from URL \""+
that._vf['url'][that._currentURLIdx]+"\","+" no other URLs left to try.");}}}
else
{if((that._currentURLIdx+1)<that._vf['url'].length)
{x3dom.debug.logWarning("Invalid SRC data, loaded from URL \""+
that._vf['url'][that._currentURLIdx]+"\", trying next specified URL");++that._currentURLIdx;that.update(shape,shaderProgram,gl,viewarea,context);}
else
{x3dom.debug.logError("Invalid SRC data, loaded from URL \""+
that._vf['url'][that._currentURLIdx]+"\","+" no other URLs left to try.");}}};},getVolume:function()
{var vol=this._mesh._vol;var shapeNode;if(!vol.isValid())
{shapeNode=this._parentNodes[0];if(typeof shapeNode._vf["bboxCenter"]!='undefined'&&typeof shapeNode._vf["bboxSize"]!='undefined')
{vol.setBoundsByCenterSize(shapeNode._vf["bboxCenter"],shapeNode._vf["bboxSize"]);}
else
{}}
return vol;}}));x3dom.registerNodeType("X3DBinaryContainerGeometryNode","Geometry3D",defineClass(x3dom.nodeTypes.X3DSpatialGeometryNode,function(ctx){x3dom.nodeTypes.X3DBinaryContainerGeometryNode.superClass.call(this,ctx);this.addField_SFVec3f(ctx,'position',0,0,0);this.addField_SFVec3f(ctx,'size',1,1,1);this.addField_MFInt32(ctx,'vertexCount',[0]);this.addField_MFString(ctx,'primType',['TRIANGLES']);this._mesh._invalidate=false;this._mesh._numCoords=0;this._mesh._numFaces=0;this._diameter=this._vf.size.length();},{getMin:function(){var vol=this._mesh._vol;if(!vol.isValid()){vol.setBoundsByCenterSize(this._vf.position,this._vf.size);}
return vol.min;},getMax:function(){var vol=this._mesh._vol;if(!vol.isValid()){vol.setBoundsByCenterSize(this._vf.position,this._vf.size);}
return vol.max;},getVolume:function(){var vol=this._mesh._vol;if(!vol.isValid()){vol.setBoundsByCenterSize(this._vf.position,this._vf.size);}
return vol;},invalidateVolume:function(){},getCenter:function(){return this._vf.position;},getDiameter:function(){return this._diameter;},needLighting:function(){var hasTris=(this._vf.primType.length&&this._vf.primType[0].indexOf("TRIANGLE")>=0);return(this._vf.lit&&hasTris);}}));x3dom.registerNodeType("BinaryGeometry","Geometry3D",defineClass(x3dom.nodeTypes.X3DBinaryContainerGeometryNode,function(ctx){x3dom.nodeTypes.BinaryGeometry.superClass.call(this,ctx);this.addField_SFString(ctx,'index',"");this.addField_SFString(ctx,'coord',"");this.addField_SFString(ctx,'normal',"");this.addField_SFString(ctx,'texCoord',"");this.addField_SFString(ctx,'color',"");this.addField_SFString(ctx,'tangent',"");this.addField_SFString(ctx,'binormal',"");this.addField_SFString(ctx,'indexType',"Uint16");this.addField_SFString(ctx,'coordType',"Float32");this.addField_SFString(ctx,'normalType',"Float32");this.addField_SFString(ctx,'texCoordType',"Float32");this.addField_SFString(ctx,'colorType',"Float32");this.addField_SFString(ctx,'tangentType',"Float32");this.addField_SFString(ctx,'binormalType',"Float32");this.addField_SFBool(ctx,'normalAsSphericalCoordinates',false);this.addField_SFBool(ctx,'rgbaColors',false);this.addField_SFInt32(ctx,'numTexCoordComponents',2);this.addField_SFBool(ctx,'normalPerVertex',true);this.addField_SFBool(ctx,'idsPerVertex',false);this.addField_SFBool(ctx,'compressed',false);this._hasStrideOffset=false;this._mesh._numPosComponents=this._vf.normalAsSphericalCoordinates?4:3;this._mesh._numTexComponents=this._vf.numTexCoordComponents;this._mesh._numColComponents=this._vf.rgbaColors?4:3;this._mesh._numNormComponents=this._vf.normalAsSphericalCoordinates?2:3;this._vertexCountSum=0;for(var i=0;i<this._vf.vertexCount.length;++i){this._vertexCountSum+=this._vf.vertexCount[i];}},{parentAdded:function(parent)
{var offsetInd,strideInd,offset,stride;offsetInd=this._vf.coord.lastIndexOf('#');strideInd=this._vf.coord.lastIndexOf('+');if(offsetInd>=0&&strideInd>=0){offset=+this._vf.coord.substring(++offsetInd,strideInd);stride=+this._vf.coord.substring(strideInd);parent._coordStrideOffset=[stride,offset];this._hasStrideOffset=true;if((offset/8)-Math.floor(offset/8)==0){this._mesh._numPosComponents=4;}}
else if(strideInd>=0){stride=+this._vf.coord.substring(strideInd);parent._coordStrideOffset=[stride,0];if((stride/8)-Math.floor(stride/8)==0){this._mesh._numPosComponents=4;}}
offsetInd=this._vf.normal.lastIndexOf('#');strideInd=this._vf.normal.lastIndexOf('+');if(offsetInd>=0&&strideInd>=0){offset=+this._vf.normal.substring(++offsetInd,strideInd);stride=+this._vf.normal.substring(strideInd);parent._normalStrideOffset=[stride,offset];}
else if(strideInd>=0){stride=+this._vf.normal.substring(strideInd);parent._normalStrideOffset=[stride,0];}
offsetInd=this._vf.texCoord.lastIndexOf('#');strideInd=this._vf.texCoord.lastIndexOf('+');if(offsetInd>=0&&strideInd>=0){offset=+this._vf.texCoord.substring(++offsetInd,strideInd);stride=+this._vf.texCoord.substring(strideInd);parent._texCoordStrideOffset=[stride,offset];}
else if(strideInd>=0){stride=+this._vf.texCoord.substring(strideInd);parent._texCoordStrideOffset=[stride,0];}
offsetInd=this._vf.color.lastIndexOf('#');strideInd=this._vf.color.lastIndexOf('+');if(offsetInd>=0&&strideInd>=0){offset=+this._vf.color.substring(++offsetInd,strideInd);stride=+this._vf.color.substring(strideInd);parent._colorStrideOffset=[stride,offset];}
else if(strideInd>=0){stride=+this._vf.color.substring(strideInd);parent._colorStrideOffset=[stride,0];}
if(this._vf.indexType!="Uint16"&&!x3dom.caps.INDEX_UINT)
x3dom.debug.logWarning("Index type "+this._vf.indexType+" problematic");},doIntersect:function(line)
{var min=this.getMin();var max=this.getMax();var isect=line.intersect(min,max);if(isect&&line.enter<line.dist){line.dist=line.enter;line.hitObject=this;line.hitPoint=line.pos.add(line.dir.multiply(line.enter));return true;}
else{return false;}},getPrecisionMax:function(type)
{switch(this._vf[type])
{case"Int8":return 127.0;case"Uint8":return 255.0;case"Int16":return 32767.0;case"Uint16":return 65535.0;case"Int32":return 2147483647.0;case"Uint32":return 4294967295.0;case"Float32":case"Float64":default:return 1.0;}}}));x3dom.registerNodeType("PopGeometryLevel","Geometry3D",defineClass(x3dom.nodeTypes.X3DGeometricPropertyNode,function(ctx){x3dom.nodeTypes.PopGeometryLevel.superClass.call(this,ctx);this.addField_SFString(ctx,'src',"");this.addField_SFInt32(ctx,'numIndices',0);this.addField_SFInt32(ctx,'vertexDataBufferOffset',0);},{getSrc:function(){return this._vf.src;},getNumIndices:function(){return this._vf.numIndices;},getVertexDataBufferOffset:function(){return this._vf.vertexDataBufferOffset;}}));x3dom.registerNodeType("PopGeometry","Geometry3D",defineClass(x3dom.nodeTypes.X3DBinaryContainerGeometryNode,function(ctx){x3dom.nodeTypes.PopGeometry.superClass.call(this,ctx);this.addField_SFVec3f(ctx,'tightSize',1,1,1);this.addField_SFVec3f(ctx,'maxBBSize',1,1,1);this.addField_SFVec3f(ctx,'bbMinModF',0,0,0);this.addField_SFVec3f(ctx,'bbMaxModF',1,1,1);this.addField_SFVec3f(ctx,'bbMin',0,0,0);this.addField_SFVec3f(ctx,'bbShiftVec',0,0,0);if(this._vf.bbMinModF.x>this._vf.bbMaxModF.x)
this._vf.bbShiftVec.x=1.0;if(this._vf.bbMinModF.y>this._vf.bbMaxModF.y)
this._vf.bbShiftVec.y=1.0;if(this._vf.bbMinModF.z>this._vf.bbMaxModF.z)
this._vf.bbShiftVec.z=1.0;this.addField_MFNode('levels',x3dom.nodeTypes.PopGeometryLevel);this.addField_SFInt32(ctx,'attributeStride',0);this.addField_SFInt32(ctx,'positionOffset',0);this.addField_SFInt32(ctx,'normalOffset',0);this.addField_SFInt32(ctx,'texcoordOffset',0);this.addField_SFInt32(ctx,'colorOffset',0);this.addField_SFInt32(ctx,'numAnchorVertices',0);this.addField_SFInt32(ctx,'positionPrecision',2);this.addField_SFInt32(ctx,'normalPrecision',1);this.addField_SFInt32(ctx,'texcoordPrecision',2);this.addField_SFInt32(ctx,'colorPrecision',1);this.addField_SFInt32(ctx,'minPrecisionLevel',-1);this.addField_SFInt32(ctx,'maxPrecisionLevel',-1);this.addField_SFFloat(ctx,'precisionFactor',1.0);this.addField_SFString(ctx,'coordType',"Uint16");this.addField_SFString(ctx,'normalType',"Uint8");this.addField_SFString(ctx,'texCoordType',"Uint16");this.addField_SFString(ctx,'colorType',"Uint8");this.addField_SFInt32(ctx,'vertexBufferSize',0);this.addField_SFBool(ctx,'indexedRendering',true);this.addField_SFBool(ctx,'sphericalNormals',false);this.addField_MFInt32(ctx,'originalVertexCount',[0]);for(var i=0;i<this._vf.vertexCount.length;++i){this._vf.originalVertexCount[i]=this._vf.vertexCount[i];}
this._vf.maxBBSize=x3dom.fields.SFVec3f.copy(this._vf.size);this._vf.size=this._vf.tightSize;this._diameter=this._vf.size.length();this._bbMinBySize=[Math.floor(this._vf.bbMin.x/this._vf.maxBBSize.x),Math.floor(this._vf.bbMin.y/this._vf.maxBBSize.y),Math.floor(this._vf.bbMin.z/this._vf.maxBBSize.z)];this._volRadius=this._vf.size.length()/2;this._volLargestRadius=this._vf.maxBBSize.length()/2;this._mesh._numPosComponents=this._vf.sphericalNormals?4:3;this._mesh._numNormComponents=this._vf.sphericalNormals?2:3;this._mesh._numTexComponents=2;this._mesh._numColComponents=3;x3dom.nodeTypes.PopGeometry.numTotalVerts+=this.getVertexCount();x3dom.nodeTypes.PopGeometry.numTotalTris+=(this.hasIndex()?this.getTotalNumberOfIndices():this.getVertexCount())/3;},{forceUpdateCoverage:function(){return true;},getBBoxShiftVec:function(){return this._vf.bbShiftVec;},getBBoxSize:function(){return this._vf.size;},hasIndex:function(){return this._vf.indexedRendering;},getTotalNumberOfIndices:function(){if(this._vf.indexedRendering){var sum=0;for(var i=0;i<this._vf.originalVertexCount.length;++i){sum+=this._vf.originalVertexCount[i];}
return sum;}
else{return 0;}},getVertexCount:function(){var sum=0;for(var i=0;i<this._vf.originalVertexCount.length;++i){sum+=this._vf.originalVertexCount[i];}
return sum;},adaptVertexCount:function(numVerts){var verts=0;for(var i=0;i<this._vf.originalVertexCount.length;++i){if((this._vf.originalVertexCount[i]+verts)<=numVerts){this._vf.vertexCount[i]=this._vf.originalVertexCount[i];verts+=this._vf.originalVertexCount[i];}
else{this._vf.vertexCount[i]=numVerts-verts;break;}}},hasNormal:function(){return(this._vf.normalOffset!=0)&&!this._vf.sphericalNormals;},hasTexCoord:function(){return(this._vf.texcoordOffset!=0);},hasColor:function(){return(this._vf.colorOffset!=0);},getPositionPrecision:function(){return this._vf.positionPrecision;},getNormalPrecision:function(){return this._vf.normalPrecision;},getTexCoordPrecision:function(){return this._vf.texcoordPrecision;},getColorPrecision:function(){return this._vf.colorPrecision;},getAttributeStride:function(){return this._vf.attributeStride;},getPositionOffset:function(){return this._vf.positionOffset;},getNormalOffset:function(){return this._vf.normalOffset;},getTexCoordOffset:function(){return this._vf.texcoordOffset;},getColorOffset:function(){return this._vf.colorOffset;},getBufferTypeStringFromByteCount:function(bytes){switch(bytes)
{case 1:return"Uint8";case 2:return"Uint16";default:return 0;}},getDataURLs:function(){var urls=[];for(var i=0;i<this._cf.levels.nodes.length;++i){urls.push(this._cf.levels.nodes[i].getSrc());}
return urls;},getNumIndicesByLevel:function(lvl){return this._cf.levels.nodes[lvl].getNumIndices();},getNumLevels:function(lvl){return this._cf.levels.nodes.length;},getVertexDataBufferOffset:function(lvl){return this._cf.levels.nodes[lvl].getVertexDataBufferOffset();},getPrecisionMax:function(type){switch(this._vf[type])
{case"Uint8":return 255.0;case"Uint16":return 65535.0;default:return 1.0;}}}));x3dom.nodeTypes.PopGeometry.ErrorToleranceFactor=1;x3dom.nodeTypes.PopGeometry.PrecisionFactorOnMove=1;x3dom.nodeTypes.PopGeometry.numRenderedVerts=0;x3dom.nodeTypes.PopGeometry.numRenderedTris=0;x3dom.nodeTypes.PopGeometry.numTotalVerts=0;x3dom.nodeTypes.PopGeometry.numTotalTris=0;x3dom.nodeTypes.PopGeometry.powLUT=[32768,16384,8192,4096,2048,1024,512,256,128,64,32,16,8,4,2,1];x3dom.registerNodeType("ImageGeometry","Geometry3D",defineClass(x3dom.nodeTypes.X3DBinaryContainerGeometryNode,function(ctx){x3dom.nodeTypes.ImageGeometry.superClass.call(this,ctx);this.addField_SFVec2f(ctx,'implicitMeshSize',256,256);this.addField_SFInt32(ctx,'numColorComponents',3);this.addField_SFInt32(ctx,'numTexCoordComponents',2);this.addField_SFNode('index',x3dom.nodeTypes.X3DTextureNode);this.addField_MFNode('coord',x3dom.nodeTypes.X3DTextureNode);this.addField_SFNode('normal',x3dom.nodeTypes.X3DTextureNode);this.addField_SFNode('texCoord',x3dom.nodeTypes.X3DTextureNode);this.addField_SFNode('color',x3dom.nodeTypes.X3DTextureNode);this._mesh._numColComponents=this._vf.numColorComponents;this._mesh._numTexComponents=this._vf.numTexCoordComponents;if(this._vf.implicitMeshSize.y==0)
this._vf.implicitMeshSize.y=this._vf.implicitMeshSize.x;if(x3dom.caps.BACKEND=='webgl'&&x3dom.caps.MAX_VERTEX_TEXTURE_IMAGE_UNITS>0){var geoCacheID='ImageGeometry_'+this._vf.implicitMeshSize.x+'_'+this._vf.implicitMeshSize.y;if(this._vf.useGeoCache&&x3dom.geoCache[geoCacheID]!==undefined)
{this._mesh=x3dom.geoCache[geoCacheID];}
else
{for(var y=0;y<this._vf.implicitMeshSize.y;y++)
{for(var x=0;x<this._vf.implicitMeshSize.x;x++)
{this._mesh._positions[0].push(x/this._vf.implicitMeshSize.x,y/this._vf.implicitMeshSize.y,0);}}
this._mesh._numFaces=this._mesh._indices[0].length/3;this._mesh._numCoords=this._mesh._positions[0].length/3;x3dom.geoCache[geoCacheID]=this._mesh;}}
this._vol=new x3dom.fields.BoxVolume();this._dirty={coord:true,normal:true,texCoord:true,color:true,index:true};},{setGeoDirty:function(){this._dirty.coord=true;this._dirty.normal=true;this._dirty.texCoords=true;this._dirty.color=true;this._dirty.index=true;},unsetGeoDirty:function(){this._dirty.coord=false;this._dirty.normal=false;this._dirty.texCoords=false;this._dirty.color=false;this._dirty.index=false;},nodeChanged:function()
{Array.forEach(this._parentNodes,function(node){node._dirty.positions=true;node._dirty.normals=true;node._dirty.texcoords=true;node._dirty.colors=true;});this._vol.invalidate();},fieldChanged:function(fieldName)
{if(fieldName=="index"||fieldName=="coord"||fieldName=="normal"||fieldName=="texCoord"||fieldName=="color"){this._dirty[fieldName]=true;this._vol.invalidate();}
else if(fieldName=="implicitMeshSize"){this._vol.invalidate();}},getMin:function(){var vol=this._vol;if(!vol.isValid()){vol.setBoundsByCenterSize(this._vf.position,this._vf.size);}
return vol.min;},getMax:function(){var vol=this._vol;if(!vol.isValid()){vol.setBoundsByCenterSize(this._vf.position,this._vf.size);}
return vol.max;},getVolume:function(){var vol=this._vol;if(!vol.isValid()){vol.setBoundsByCenterSize(this._vf.position,this._vf.size);}
return vol;},numCoordinateTextures:function()
{return this._cf.coord.nodes.length;},getIndexTexture:function()
{if(this._cf.index.node){this._cf.index.node._type="IG_index";return this._cf.index.node;}else{return null;}},getIndexTextureURL:function()
{if(this._cf.index.node){return this._cf.index.node._vf.url;}else{return null;}},getCoordinateTexture:function(pos)
{if(this._cf.coord.nodes[pos]){this._cf.coord.nodes[pos]._type="IG_coords"+pos;return this._cf.coord.nodes[pos];}else{return null;}},getCoordinateTextureURL:function(pos)
{if(this._cf.coord.nodes[pos]){return this._cf.coord.nodes[pos]._vf.url;}else{return null;}},getCoordinateTextureURLs:function()
{var urls=[];for(var i=0;i<this._cf.coord.nodes.length;i++)
{urls.push(this._cf.coord.nodes[i]._vf.url);}
return urls;},getNormalTexture:function()
{if(this._cf.normal.node){this._cf.normal.node._type="IG_normals";return this._cf.normal.node;}else{return null;}},getNormalTextureURL:function()
{if(this._cf.normal.node){return this._cf.normal.node._vf.url;}else{return null;}},getTexCoordTexture:function()
{if(this._cf.texCoord.node){this._cf.texCoord.node._type="IG_texCoords";return this._cf.texCoord.node;}else{return null;}},getTexCoordTextureURL:function()
{if(this._cf.texCoord.node){return this._cf.texCoord.node._vf.url;}else{return null;}},getColorTexture:function()
{if(this._cf.color.node){this._cf.color.node._type="IG_colors";return this._cf.color.node;}else{return null;}},getColorTextureURL:function()
{if(this._cf.color.node){return this._cf.color.node._vf.url;}else{return null;}},getTextures:function()
{var textures=[];var index=this.getIndexTexture();if(index)textures.push(index);for(i=0;i<this.numCoordinateTextures();i++){var coord=this.getCoordinateTexture(i);if(coord)textures.push(coord);}
var normal=this.getNormalTexture();if(normal)textures.push(normal);var texCoord=this.getTexCoordTexture();if(texCoord)textures.push(texCoord);var color=this.getColorTexture();if(color)textures.push(color);return textures;}}));x3dom.registerNodeType("IndexedFaceSet","Geometry3D",defineClass(x3dom.nodeTypes.X3DComposedGeometryNode,function(ctx){x3dom.nodeTypes.IndexedFaceSet.superClass.call(this,ctx);this.addField_SFFloat(ctx,'creaseAngle',0);this.addField_SFBool(ctx,'convex',true);this.addField_MFInt32(ctx,'coordIndex',[]);this.addField_MFInt32(ctx,'normalIndex',[]);this.addField_MFInt32(ctx,'colorIndex',[]);this.addField_MFInt32(ctx,'texCoordIndex',[]);},{nodeChanged:function()
{var time0=new Date().getTime();this.handleAttribs();var indexes=this._vf.coordIndex;if(indexes.length&&indexes[indexes.length-1]!=-1)
{indexes.push(-1);}
var normalInd=this._vf.normalIndex;var texCoordInd=this._vf.texCoordIndex;var colorInd=this._vf.colorIndex;var hasNormal=false,hasNormalInd=false;var hasTexCoord=false,hasTexCoordInd=false;var hasColor=false,hasColorInd=false;var colPerVert=this._vf.colorPerVertex;var normPerVert=this._vf.normalPerVertex;if(normalInd.length>0)
{hasNormalInd=true;}
if(texCoordInd.length>0)
{hasTexCoordInd=true;}
if(colorInd.length>0)
{hasColorInd=true;}
var positions,normals,texCoords,colors;var coordNode=this._cf.coord.node;x3dom.debug.assert(coordNode);positions=coordNode.getPoints();var normalNode=this._cf.normal.node;if(normalNode)
{hasNormal=true;normals=normalNode._vf.vector;}
else{hasNormal=false;}
var texMode="",numTexComponents=2;var texCoordNode=this._cf.texCoord.node;if(x3dom.isa(texCoordNode,x3dom.nodeTypes.MultiTextureCoordinate)){if(texCoordNode._cf.texCoord.nodes.length)
texCoordNode=texCoordNode._cf.texCoord.nodes[0];}
if(texCoordNode)
{if(texCoordNode._vf.point){hasTexCoord=true;texCoords=texCoordNode._vf.point;if(x3dom.isa(texCoordNode,x3dom.nodeTypes.TextureCoordinate3D)){numTexComponents=3;}}
else if(texCoordNode._vf.mode){texMode=texCoordNode._vf.mode;}}
else{hasTexCoord=false;}
this._mesh._numTexComponents=numTexComponents;var numColComponents=3;var colorNode=this._cf.color.node;if(colorNode)
{hasColor=true;colors=colorNode._vf.color;if(x3dom.isa(colorNode,x3dom.nodeTypes.ColorRGBA)){numColComponents=4;}}
else{hasColor=false;}
this._mesh._numColComponents=numColComponents;this._mesh._indices[0]=[];this._mesh._positions[0]=[];this._mesh._normals[0]=[];this._mesh._texCoords[0]=[];this._mesh._colors[0]=[];var i,j,t,cnt,faceCnt;var p0,p1,p2,n0,n1,n2,t0,t1,t2,c0,c1,c2;if((this._vf.creaseAngle<=x3dom.fields.Eps)||(positions.length>x3dom.Utils.maxIndexableCoords)||(hasNormal&&hasNormalInd)||(hasTexCoord&&hasTexCoordInd)||(hasColor&&hasColorInd))
{if(this._vf.creaseAngle<=x3dom.fields.Eps)
x3dom.debug.logWarning('Fallback to inefficient multi-index mode since creaseAngle=0.');if(this._vf.convex){t=0;cnt=0;faceCnt=0;this._mesh._multiIndIndices=[];this._mesh._posSize=positions.length;for(i=0;i<indexes.length;++i)
{if(indexes[i]==-1){t=0;faceCnt++;continue;}
if(hasNormalInd){x3dom.debug.assert(normalInd[i]!=-1);}
if(hasTexCoordInd){x3dom.debug.assert(texCoordInd[i]!=-1);}
if(hasColorInd){x3dom.debug.assert(colorInd[i]!=-1);}
switch(t)
{case 0:p0=+indexes[i];if(hasNormalInd&&normPerVert){n0=+normalInd[i];}
else if(hasNormalInd&&!normPerVert){n0=+normalInd[faceCnt];}
else if(normPerVert){n0=p0;}
else{n0=faceCnt;}
if(hasTexCoordInd){t0=+texCoordInd[i];}
else{t0=p0;}
if(hasColorInd&&colPerVert){c0=+colorInd[i];}
else if(hasColorInd&&!colPerVert){c0=+colorInd[faceCnt];}
else if(colPerVert){c0=p0;}
else{c0=faceCnt;}
t=1;break;case 1:p1=+indexes[i];if(hasNormalInd&&normPerVert){n1=+normalInd[i];}
else if(hasNormalInd&&!normPerVert){n1=+normalInd[faceCnt];}
else if(normPerVert){n1=p1;}
else{n1=faceCnt;}
if(hasTexCoordInd){t1=+texCoordInd[i];}
else{t1=p1;}
if(hasColorInd&&colPerVert){c1=+colorInd[i];}
else if(hasColorInd&&!colPerVert){c1=+colorInd[faceCnt];}
else if(colPerVert){c1=p1;}
else{c1=faceCnt;}
t=2;break;case 2:p2=+indexes[i];if(hasNormalInd&&normPerVert){n2=+normalInd[i];}
else if(hasNormalInd&&!normPerVert){n2=+normalInd[faceCnt];}
else if(normPerVert){n2=p2;}
else{n2=faceCnt;}
if(hasTexCoordInd){t2=+texCoordInd[i];}
else{t2=p2;}
if(hasColorInd&&colPerVert){c2=+colorInd[i];}
else if(hasColorInd&&!colPerVert){c2=+colorInd[faceCnt];}
else if(colPerVert){c2=p2;}
else{c2=faceCnt;}
t=3;this._mesh._positions[0].push(positions[p0].x);this._mesh._positions[0].push(positions[p0].y);this._mesh._positions[0].push(positions[p0].z);this._mesh._positions[0].push(positions[p1].x);this._mesh._positions[0].push(positions[p1].y);this._mesh._positions[0].push(positions[p1].z);this._mesh._positions[0].push(positions[p2].x);this._mesh._positions[0].push(positions[p2].y);this._mesh._positions[0].push(positions[p2].z);if(hasNormal){this._mesh._normals[0].push(normals[n0].x);this._mesh._normals[0].push(normals[n0].y);this._mesh._normals[0].push(normals[n0].z);this._mesh._normals[0].push(normals[n1].x);this._mesh._normals[0].push(normals[n1].y);this._mesh._normals[0].push(normals[n1].z);this._mesh._normals[0].push(normals[n2].x);this._mesh._normals[0].push(normals[n2].y);this._mesh._normals[0].push(normals[n2].z);}
this._mesh._multiIndIndices.push(p0,p1,p2);if(hasColor){this._mesh._colors[0].push(colors[c0].r);this._mesh._colors[0].push(colors[c0].g);this._mesh._colors[0].push(colors[c0].b);if(numColComponents===4){this._mesh._colors[0].push(colors[c0].a);}
this._mesh._colors[0].push(colors[c1].r);this._mesh._colors[0].push(colors[c1].g);this._mesh._colors[0].push(colors[c1].b);if(numColComponents===4){this._mesh._colors[0].push(colors[c1].a);}
this._mesh._colors[0].push(colors[c2].r);this._mesh._colors[0].push(colors[c2].g);this._mesh._colors[0].push(colors[c2].b);if(numColComponents===4){this._mesh._colors[0].push(colors[c2].a);}}
if(hasTexCoord){this._mesh._texCoords[0].push(texCoords[t0].x);this._mesh._texCoords[0].push(texCoords[t0].y);if(numTexComponents===3){this._mesh._texCoords[0].push(texCoords[t0].z);}
this._mesh._texCoords[0].push(texCoords[t1].x);this._mesh._texCoords[0].push(texCoords[t1].y);if(numTexComponents===3){this._mesh._texCoords[0].push(texCoords[t1].z);}
this._mesh._texCoords[0].push(texCoords[t2].x);this._mesh._texCoords[0].push(texCoords[t2].y);if(numTexComponents===3){this._mesh._texCoords[0].push(texCoords[t2].z);}}
break;case 3:p1=p2;t1=t2;if(normPerVert){n1=n2;}
if(colPerVert){c1=c2;}
p2=+indexes[i];if(hasNormalInd&&normPerVert){n2=+normalInd[i];}else if(hasNormalInd&&!normPerVert){}else if(normPerVert){n2=p2;}else{n2=faceCnt;}
if(hasTexCoordInd){t2=+texCoordInd[i];}else{t2=p2;}
if(hasColorInd&&colPerVert){c2=+colorInd[i];}else if(hasColorInd&&!colPerVert){}else if(colPerVert){c2=p2;}else{c2=faceCnt;}
this._mesh._positions[0].push(positions[p0].x);this._mesh._positions[0].push(positions[p0].y);this._mesh._positions[0].push(positions[p0].z);this._mesh._positions[0].push(positions[p1].x);this._mesh._positions[0].push(positions[p1].y);this._mesh._positions[0].push(positions[p1].z);this._mesh._positions[0].push(positions[p2].x);this._mesh._positions[0].push(positions[p2].y);this._mesh._positions[0].push(positions[p2].z);if(hasNormal){this._mesh._normals[0].push(normals[n0].x);this._mesh._normals[0].push(normals[n0].y);this._mesh._normals[0].push(normals[n0].z);this._mesh._normals[0].push(normals[n1].x);this._mesh._normals[0].push(normals[n1].y);this._mesh._normals[0].push(normals[n1].z);this._mesh._normals[0].push(normals[n2].x);this._mesh._normals[0].push(normals[n2].y);this._mesh._normals[0].push(normals[n2].z);}
this._mesh._multiIndIndices.push(p0,p1,p2);if(hasColor){this._mesh._colors[0].push(colors[c0].r);this._mesh._colors[0].push(colors[c0].g);this._mesh._colors[0].push(colors[c0].b);if(numColComponents===4){this._mesh._colors[0].push(colors[c0].a);}
this._mesh._colors[0].push(colors[c1].r);this._mesh._colors[0].push(colors[c1].g);this._mesh._colors[0].push(colors[c1].b);if(numColComponents===4){this._mesh._colors[0].push(colors[c1].a);}
this._mesh._colors[0].push(colors[c2].r);this._mesh._colors[0].push(colors[c2].g);this._mesh._colors[0].push(colors[c2].b);if(numColComponents===4){this._mesh._colors[0].push(colors[c2].a);}}
if(hasTexCoord){this._mesh._texCoords[0].push(texCoords[t0].x);this._mesh._texCoords[0].push(texCoords[t0].y);if(numTexComponents===3){this._mesh._texCoords[0].push(texCoords[t0].z);}
this._mesh._texCoords[0].push(texCoords[t1].x);this._mesh._texCoords[0].push(texCoords[t1].y);if(numTexComponents===3){this._mesh._texCoords[0].push(texCoords[t1].z);}
this._mesh._texCoords[0].push(texCoords[t2].x);this._mesh._texCoords[0].push(texCoords[t2].y);if(numTexComponents===3){this._mesh._texCoords[0].push(texCoords[t2].z);}}
break;default:}}}
else{var linklist=new x3dom.DoublyLinkedList();var data={};cnt=0;faceCnt=0;for(i=0;i<indexes.length;++i)
{if(indexes[i]==-1){var multi_index_data=x3dom.EarClipping.getMultiIndexes(linklist);for(j=0;j<multi_index_data.indices.length;j++)
{this._mesh._indices[0].push(cnt);cnt++;this._mesh._positions[0].push(multi_index_data.point[j].x,multi_index_data.point[j].y,multi_index_data.point[j].z);if(hasNormal){this._mesh._normals[0].push(multi_index_data.normals[j].x,multi_index_data.normals[j].y,multi_index_data.normals[j].z);}
if(hasColor){this._mesh._colors[0].push(multi_index_data.colors[j].r,multi_index_data.colors[j].g,multi_index_data.colors[j].b);if(numColComponents===4){this._mesh._colors[0].push(multi_index_data.colors[j].a);}}
if(hasTexCoord){this._mesh._texCoords[0].push(multi_index_data.texCoords[j].x,multi_index_data.texCoords[j].y);if(numTexComponents===3){this._mesh._texCoords[0].push(multi_index_data.texCoords[j].z);}}}
linklist=new x3dom.DoublyLinkedList();faceCnt++;continue;}
if(hasNormal){if(hasNormalInd&&normPerVert){data.normals=normals[normalInd[i]];}else if(hasNormalInd&&!normPerVert){data.normals=normals[normalInd[faceCnt]];}else{data.normals=normals[indexes[i]];}}
if(hasColor){if(hasColorInd&&colPerVert){data.colors=colors[colorInd[i]];}else if(hasColorInd&&!colPerVert){data.colors=colors[colorInd[faceCnt]];}else if(colPerVert){data.colors=colors[indexes[i]];}else{data.colors=colors[faceCnt];}}
if(hasTexCoord){if(hasTexCoordInd){data.texCoords=texCoords[texCoordInd[i]];}else{data.texCoords=texCoords[indexes[i]];}}
linklist.appendNode(new x3dom.DoublyLinkedList.ListNode(positions[indexes[i]],indexes[i],data.normals,data.colors,data.texCoords));}
this._mesh.splitMesh();}
if(!hasNormal){this._mesh.calcNormals(this._vf.creaseAngle,this._vf.ccw);}
if(!hasTexCoord){this._mesh.calcTexCoords(texMode);}}
else
{t=0;if(this._vf.convex){for(i=0;i<indexes.length;++i)
{if(indexes[i]==-1){t=0;continue;}
switch(t){case 0:n0=+indexes[i];t=1;break;case 1:n1=+indexes[i];t=2;break;case 2:n2=+indexes[i];t=3;this._mesh._indices[0].push(n0,n1,n2);break;case 3:n1=n2;n2=+indexes[i];this._mesh._indices[0].push(n0,n1,n2);break;}}}
else{linklist=new x3dom.DoublyLinkedList();for(i=0;i<indexes.length;++i)
{if(indexes[i]==-1){var linklist_indices=x3dom.EarClipping.getIndexes(linklist);for(j=0;j<linklist_indices.length;j++){this._mesh._indices[0].push(linklist_indices[j]);}
linklist=new x3dom.DoublyLinkedList();continue;}
linklist.appendNode(new x3dom.DoublyLinkedList.ListNode(positions[indexes[i]],indexes[i]));}}
this._mesh._positions[0]=positions.toGL();if(hasNormal){this._mesh._normals[0]=normals.toGL();}
else{this._mesh.calcNormals(this._vf.creaseAngle,this._vf.ccw);}
if(hasTexCoord){this._mesh._texCoords[0]=texCoords.toGL();this._mesh._numTexComponents=numTexComponents;}
else{this._mesh.calcTexCoords(texMode);}
if(hasColor){this._mesh._colors[0]=colors.toGL();this._mesh._numColComponents=numColComponents;}}
this.invalidateVolume();this._mesh._numFaces=0;this._mesh._numCoords=0;for(i=0;i<this._mesh._positions.length;i++){var indexLength=this._mesh._indices[i].length;var numCoords=this._mesh._positions[i].length/3;this._mesh._numCoords+=numCoords;if(indexLength>0)
this._mesh._numFaces+=indexLength/3;else
this._mesh._numFaces+=numCoords/3;}},fieldChanged:function(fieldName)
{if(fieldName!="coord"&&fieldName!="normal"&&fieldName!="texCoord"&&fieldName!="color"&&fieldName!="coordIndex")
{x3dom.debug.logWarning("IndexedFaceSet: fieldChanged for "+
fieldName+" not yet implemented!");return;}
var pnts=this._cf.coord.node._vf.point;var n=pnts.length;var texCoordNode=this._cf.texCoord.node;if(x3dom.isa(texCoordNode,x3dom.nodeTypes.MultiTextureCoordinate)){if(texCoordNode._cf.texCoord.nodes.length)
texCoordNode=texCoordNode._cf.texCoord.nodes[0];}
if(((this._vf.creaseAngle<=x3dom.fields.Eps)||(n>x3dom.Utils.maxIndexableCoords)||(this._vf.normalIndex.length>0&&this._cf.normal.node)||(this._vf.texCoordIndex.length>0&&texCoordNode)||(this._vf.colorIndex.length>0&&this._cf.color.node))&&this._mesh._multiIndIndices)
{var needNormals=!this._cf.normal.node&&this._vf.normalUpdateMode.toLowerCase()!='none';n=this._mesh._multiIndIndices.length;this._mesh._positions[0]=[];this._mesh._indices[0]=[];if(fieldName=="coord"&&n)
{if(needNormals){this._mesh._normals[0]=[];}
for(i=0;i<n;i+=3){var ind0=this._mesh._multiIndIndices[i];var ind1=this._mesh._multiIndIndices[i+1];var ind2=this._mesh._multiIndIndices[i+2];var pos0=pnts[ind0];var pos1=pnts[ind1];var pos2=pnts[ind2];this._mesh._positions[0].push(pos0.x,pos0.y,pos0.z);this._mesh._positions[0].push(pos1.x,pos1.y,pos1.z);this._mesh._positions[0].push(pos2.x,pos2.y,pos2.z);if(needNormals){var a=pos0.subtract(pos1);var b=pos1.subtract(pos2);var norm=a.cross(b).normalize();if(!this._vf.ccw)
norm=norm.negate();this._mesh._normals[0].push(norm.x,norm.y,norm.z);this._mesh._normals[0].push(norm.x,norm.y,norm.z);this._mesh._normals[0].push(norm.x,norm.y,norm.z);}}
this.invalidateVolume();Array.forEach(this._parentNodes,function(node){node._dirty.positions=true;if(needNormals)
node._dirty.normals=true;});return;}
this._mesh._normals[0]=[];this._mesh._texCoords[0]=[];this._mesh._colors[0]=[];var indexes=this._vf.coordIndex;var normalInd=this._vf.normalIndex;var texCoordInd=this._vf.texCoordIndex;var colorInd=this._vf.colorIndex;var hasNormal=false,hasNormalInd=false;var hasTexCoord=false,hasTexCoordInd=false;var hasColor=false,hasColorInd=false;var colPerVert=this._vf.colorPerVertex;var normPerVert=this._vf.normalPerVertex;if(normalInd.length>0)
{hasNormalInd=true;}
if(texCoordInd.length>0)
{hasTexCoordInd=true;}
if(colorInd.length>0)
{hasColorInd=true;}
var positions,normals,texCoords,colors;var coordNode=this._cf.coord.node;x3dom.debug.assert(coordNode);positions=coordNode.getPoints();var normalNode=this._cf.normal.node;if(normalNode)
{hasNormal=true;normals=normalNode._vf.vector;}
else{hasNormal=false;}
var texMode="",numTexComponents=2;texCoordNode=this._cf.texCoord.node;if(x3dom.isa(texCoordNode,x3dom.nodeTypes.MultiTextureCoordinate)){if(texCoordNode._cf.texCoord.nodes.length)
texCoordNode=texCoordNode._cf.texCoord.nodes[0];}
if(texCoordNode)
{if(texCoordNode._vf.point){hasTexCoord=true;texCoords=texCoordNode._vf.point;if(x3dom.isa(texCoordNode,x3dom.nodeTypes.TextureCoordinate3D)){numTexComponents=3;}}
else if(texCoordNode._vf.mode){texMode=texCoordNode._vf.mode;}}
else{hasTexCoord=false;}
this._mesh._numTexComponents=numTexComponents;var numColComponents=3;var colorNode=this._cf.color.node;if(colorNode)
{hasColor=true;colors=colorNode._vf.color;if(x3dom.isa(colorNode,x3dom.nodeTypes.ColorRGBA)){numColComponents=4;}}
else{hasColor=false;}
this._mesh._numColComponents=numColComponents;var i,j,t,cnt,faceCnt;var p0,p1,p2,n0,n1,n2,t0,t1,t2,c0,c1,c2;if(this._vf.convex){t=0;cnt=0;faceCnt=0;this._mesh._multiIndIndices=[];this._mesh._posSize=positions.length;for(i=0;i<indexes.length;++i)
{if(indexes[i]==-1){t=0;faceCnt++;continue;}
if(hasNormalInd){x3dom.debug.assert(normalInd[i]!=-1);}
if(hasTexCoordInd){x3dom.debug.assert(texCoordInd[i]!=-1);}
if(hasColorInd){x3dom.debug.assert(colorInd[i]!=-1);}
switch(t)
{case 0:p0=+indexes[i];if(hasNormalInd&&normPerVert){n0=+normalInd[i];}
else if(hasNormalInd&&!normPerVert){n0=+normalInd[faceCnt];}
else if(normPerVert){n0=p0;}
else{n0=faceCnt;}
if(hasTexCoordInd){t0=+texCoordInd[i];}
else{t0=p0;}
if(hasColorInd&&colPerVert){c0=+colorInd[i];}
else if(hasColorInd&&!colPerVert){c0=+colorInd[faceCnt];}
else if(colPerVert){c0=p0;}
else{c0=faceCnt;}
t=1;break;case 1:p1=+indexes[i];if(hasNormalInd&&normPerVert){n1=+normalInd[i];}
else if(hasNormalInd&&!normPerVert){n1=+normalInd[faceCnt];}
else if(normPerVert){n1=p1;}
else{n1=faceCnt;}
if(hasTexCoordInd){t1=+texCoordInd[i];}
else{t1=p1;}
if(hasColorInd&&colPerVert){c1=+colorInd[i];}
else if(hasColorInd&&!colPerVert){c1=+colorInd[faceCnt];}
else if(colPerVert){c1=p1;}
else{c1=faceCnt;}
t=2;break;case 2:p2=+indexes[i];if(hasNormalInd&&normPerVert){n2=+normalInd[i];}
else if(hasNormalInd&&!normPerVert){n2=+normalInd[faceCnt];}
else if(normPerVert){n2=p2;}
else{n2=faceCnt;}
if(hasTexCoordInd){t2=+texCoordInd[i];}
else{t2=p2;}
if(hasColorInd&&colPerVert){c2=+colorInd[i];}
else if(hasColorInd&&!colPerVert){c2=+colorInd[faceCnt];}
else if(colPerVert){c2=p2;}
else{c2=faceCnt;}
t=3;this._mesh._positions[0].push(positions[p0].x);this._mesh._positions[0].push(positions[p0].y);this._mesh._positions[0].push(positions[p0].z);this._mesh._positions[0].push(positions[p1].x);this._mesh._positions[0].push(positions[p1].y);this._mesh._positions[0].push(positions[p1].z);this._mesh._positions[0].push(positions[p2].x);this._mesh._positions[0].push(positions[p2].y);this._mesh._positions[0].push(positions[p2].z);if(hasNormal){this._mesh._normals[0].push(normals[n0].x);this._mesh._normals[0].push(normals[n0].y);this._mesh._normals[0].push(normals[n0].z);this._mesh._normals[0].push(normals[n1].x);this._mesh._normals[0].push(normals[n1].y);this._mesh._normals[0].push(normals[n1].z);this._mesh._normals[0].push(normals[n2].x);this._mesh._normals[0].push(normals[n2].y);this._mesh._normals[0].push(normals[n2].z);}
this._mesh._multiIndIndices.push(p0,p1,p2);if(hasColor){this._mesh._colors[0].push(colors[c0].r);this._mesh._colors[0].push(colors[c0].g);this._mesh._colors[0].push(colors[c0].b);if(numColComponents===4){this._mesh._colors[0].push(colors[c0].a);}
this._mesh._colors[0].push(colors[c1].r);this._mesh._colors[0].push(colors[c1].g);this._mesh._colors[0].push(colors[c1].b);if(numColComponents===4){this._mesh._colors[0].push(colors[c1].a);}
this._mesh._colors[0].push(colors[c2].r);this._mesh._colors[0].push(colors[c2].g);this._mesh._colors[0].push(colors[c2].b);if(numColComponents===4){this._mesh._colors[0].push(colors[c2].a);}}
if(hasTexCoord){this._mesh._texCoords[0].push(texCoords[t0].x);this._mesh._texCoords[0].push(texCoords[t0].y);if(numTexComponents===3){this._mesh._texCoords[0].push(texCoords[t0].z);}
this._mesh._texCoords[0].push(texCoords[t1].x);this._mesh._texCoords[0].push(texCoords[t1].y);if(numTexComponents===3){this._mesh._texCoords[0].push(texCoords[t1].z);}
this._mesh._texCoords[0].push(texCoords[t2].x);this._mesh._texCoords[0].push(texCoords[t2].y);if(numTexComponents===3){this._mesh._texCoords[0].push(texCoords[t2].z);}}
break;case 3:p1=p2;t1=t2;if(normPerVert){n1=n2;}
if(colPerVert){c1=c2;}
p2=+indexes[i];if(hasNormalInd&&normPerVert){n2=+normalInd[i];}else if(hasNormalInd&&!normPerVert){}else if(normPerVert){n2=p2;}else{n2=faceCnt;}
if(hasTexCoordInd){t2=+texCoordInd[i];}else{t2=p2;}
if(hasColorInd&&colPerVert){c2=+colorInd[i];}else if(hasColorInd&&!colPerVert){}else if(colPerVert){c2=p2;}else{c2=faceCnt;}
this._mesh._positions[0].push(positions[p0].x);this._mesh._positions[0].push(positions[p0].y);this._mesh._positions[0].push(positions[p0].z);this._mesh._positions[0].push(positions[p1].x);this._mesh._positions[0].push(positions[p1].y);this._mesh._positions[0].push(positions[p1].z);this._mesh._positions[0].push(positions[p2].x);this._mesh._positions[0].push(positions[p2].y);this._mesh._positions[0].push(positions[p2].z);if(hasNormal){this._mesh._normals[0].push(normals[n0].x);this._mesh._normals[0].push(normals[n0].y);this._mesh._normals[0].push(normals[n0].z);this._mesh._normals[0].push(normals[n1].x);this._mesh._normals[0].push(normals[n1].y);this._mesh._normals[0].push(normals[n1].z);this._mesh._normals[0].push(normals[n2].x);this._mesh._normals[0].push(normals[n2].y);this._mesh._normals[0].push(normals[n2].z);}
this._mesh._multiIndIndices.push(p0,p1,p2);if(hasColor){this._mesh._colors[0].push(colors[c0].r);this._mesh._colors[0].push(colors[c0].g);this._mesh._colors[0].push(colors[c0].b);if(numColComponents===4){this._mesh._colors[0].push(colors[c0].a);}
this._mesh._colors[0].push(colors[c1].r);this._mesh._colors[0].push(colors[c1].g);this._mesh._colors[0].push(colors[c1].b);if(numColComponents===4){this._mesh._colors[0].push(colors[c1].a);}
this._mesh._colors[0].push(colors[c2].r);this._mesh._colors[0].push(colors[c2].g);this._mesh._colors[0].push(colors[c2].b);if(numColComponents===4){this._mesh._colors[0].push(colors[c2].a);}}
if(hasTexCoord){this._mesh._texCoords[0].push(texCoords[t0].x);this._mesh._texCoords[0].push(texCoords[t0].y);if(numTexComponents===3){this._mesh._texCoords[0].push(texCoords[t0].z);}
this._mesh._texCoords[0].push(texCoords[t1].x);this._mesh._texCoords[0].push(texCoords[t1].y);if(numTexComponents===3){this._mesh._texCoords[0].push(texCoords[t1].z);}
this._mesh._texCoords[0].push(texCoords[t2].x);this._mesh._texCoords[0].push(texCoords[t2].y);if(numTexComponents===3){this._mesh._texCoords[0].push(texCoords[t2].z);}}
break;default:}}}
else{var linklist=new x3dom.DoublyLinkedList();var data={};cnt=0;faceCnt=0;for(i=0;i<indexes.length;++i)
{if(indexes[i]==-1){var multi_index_data=x3dom.EarClipping.getMultiIndexes(linklist);for(j=0;j<multi_index_data.indices.length;j++)
{this._mesh._indices[0].push(cnt);cnt++;this._mesh._positions[0].push(multi_index_data.point[j].x,multi_index_data.point[j].y,multi_index_data.point[j].z);if(hasNormal){this._mesh._normals[0].push(multi_index_data.normals[j].x,multi_index_data.normals[j].y,multi_index_data.normals[j].z);}
if(hasColor){this._mesh._colors[0].push(multi_index_data.colors[j].r,multi_index_data.colors[j].g,multi_index_data.colors[j].b);if(numColComponents===4){this._mesh._colors[0].push(multi_index_data.colors[j].a);}}
if(hasTexCoord){this._mesh._texCoords[0].push(multi_index_data.texCoords[j].x,multi_index_data.texCoords[j].y);if(numTexComponents===3){this._mesh._texCoords[0].push(multi_index_data.texCoords[j].z);}}}
linklist=new x3dom.DoublyLinkedList();faceCnt++;continue;}
if(hasNormal){if(hasNormalInd&&normPerVert){data.normals=normals[normalInd[i]];}else if(hasNormalInd&&!normPerVert){data.normals=normals[normalInd[faceCnt]];}else{data.normals=normals[indexes[i]];}}
if(hasColor){if(hasColorInd&&colPerVert){data.colors=colors[colorInd[i]];}else if(hasColorInd&&!colPerVert){data.colors=colors[colorInd[faceCnt]];}else{data.colors=colors[indexes[i]];}}
if(hasTexCoord){if(hasTexCoordInd){data.texCoords=texCoords[texCoordInd[i]];}else{data.texCoords=texCoords[indexes[i]];}}
linklist.appendNode(new x3dom.DoublyLinkedList.ListNode(positions[indexes[i]],indexes[i],data.normals,data.colors,data.texCoords));}
this._mesh.splitMesh();}
if(!hasNormal){this._mesh.calcNormals(this._vf.creaseAngle,this._vf.ccw);}
if(!hasTexCoord){this._mesh.calcTexCoords(texMode);}
this.invalidateVolume();this._mesh._numFaces=0;this._mesh._numCoords=0;for(i=0;i<this._mesh._positions.length;i++){var indexLength=this._mesh._indices[i].length;var numCoords=this._mesh._positions[i].length/3;this._mesh._numCoords+=numCoords;if(indexLength>0)
this._mesh._numFaces+=indexLength/3;else
this._mesh._numFaces+=numCoords/3;}
Array.forEach(this._parentNodes,function(node){node.setGeoDirty();});}
else{if(fieldName=="coord")
{var needNormals=!this._cf.normal.node&&this._vf.normalUpdateMode.toLowerCase()!='none';this._mesh._positions[0]=pnts.toGL();if(needNormals){this._mesh.calcNormals(this._vf.creaseAngle,this._vf.ccw);}
this.invalidateVolume();Array.forEach(this._parentNodes,function(node){node._dirty.positions=true;if(needNormals)
node._dirty.normals=true;node.invalidateVolume();});}
else if(fieldName=="color")
{pnts=this._cf.color.node._vf.color;this._mesh._colors[0]=pnts.toGL();Array.forEach(this._parentNodes,function(node){node._dirty.colors=true;});}
else if(fieldName=="normal")
{pnts=this._cf.normal.node._vf.vector;this._mesh._normals[0]=pnts.toGL();Array.forEach(this._parentNodes,function(node){node._dirty.normals=true;});}
else if(fieldName=="texCoord")
{texCoordNode=this._cf.texCoord.node;if(x3dom.isa(texCoordNode,x3dom.nodeTypes.MultiTextureCoordinate)){if(texCoordNode._cf.texCoord.nodes.length)
texCoordNode=texCoordNode._cf.texCoord.nodes[0];}
pnts=texCoordNode._vf.point;this._mesh._texCoords[0]=pnts.toGL();Array.forEach(this._parentNodes,function(node){node._dirty.texcoords=true;});}
else if(fieldName=="coordIndex")
{needNormals=!this._cf.normal.node&&this._vf.normalUpdateMode.toLowerCase()!='none';indexes=this._vf.coordIndex;t=0;n=indexes.length;this._mesh._indices[0]=[];for(i=0;i<n;++i){if(indexes[i]==-1){t=0;}
else{switch(t){case 0:p0=+indexes[i];t=1;break;case 1:p1=+indexes[i];t=2;break;case 2:p2=+indexes[i];t=3;this._mesh._indices[0].push(p0,p1,p2);break;case 3:p1=p2;p2=+indexes[i];this._mesh._indices[0].push(p0,p1,p2);break;}}}
if(needNormals){this._mesh.calcNormals(this._vf.creaseAngle,this._vf.ccw);}
Array.forEach(this._parentNodes,function(node){node._dirty.indexes=true;if(needNormals)
node._dirty.normals=true;});}}}}));x3dom.registerNodeType("X3DTexture3DNode","Texturing3D",defineClass(x3dom.nodeTypes.X3DTextureNode,function(ctx){x3dom.nodeTypes.X3DTexture3DNode.superClass.call(this,ctx);}));x3dom.registerNodeType("ComposedTexture3D","Texturing3D",defineClass(x3dom.nodeTypes.X3DTexture3DNode,function(ctx){x3dom.nodeTypes.ComposedTexture3D.superClass.call(this,ctx);this.addField_MFNode('texture',x3dom.nodeTypes.X3DTexture3DNode);}));x3dom.registerNodeType("ImageTexture3D","Texturing3D",defineClass(x3dom.nodeTypes.X3DTexture3DNode,function(ctx){x3dom.nodeTypes.ImageTexture3D.superClass.call(this,ctx);}));x3dom.registerNodeType("PixelTexture3D","Texturing3D",defineClass(x3dom.nodeTypes.X3DTexture3DNode,function(ctx){x3dom.nodeTypes.PixelTexture3D.superClass.call(this,ctx);}));x3dom.registerNodeType("TextureCoordinate3D","Texturing3D",defineClass(x3dom.nodeTypes.X3DTextureCoordinateNode,function(ctx){x3dom.nodeTypes.TextureCoordinate3D.superClass.call(this,ctx);this.addField_MFVec3f(ctx,'point',[]);}));x3dom.registerNodeType("TextureTransform3D","Texturing3D",defineClass(x3dom.nodeTypes.X3DTextureTransformNode,function(ctx){x3dom.nodeTypes.TextureTransform3D.superClass.call(this,ctx);this.addField_SFVec3f(ctx,'center',0,0,0);this.addField_SFRotation(ctx,'rotation',0,0,1,0);this.addField_SFVec3f(ctx,'scale',1,1,1);this.addField_SFVec3f(ctx,'translation',0,0,0);this.addField_SFRotation(ctx,'scaleOrientation',0,0,1,0);}));x3dom.registerNodeType("TextureTransformMatrix3D","Texturing3D",defineClass(x3dom.nodeTypes.X3DTextureTransformNode,function(ctx){x3dom.nodeTypes.TextureTransformMatrix3D.superClass.call(this,ctx);this.addField_SFMatrix4f(ctx,'matrix',1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1);}));x3dom.registerNodeType("X3DPointingDeviceSensorNode","PointingDeviceSensor",defineClass(x3dom.nodeTypes.X3DSensorNode,function(ctx)
{x3dom.nodeTypes.X3DPointingDeviceSensorNode.superClass.call(this,ctx);},{pointerPressedOverSibling:function(event)
{if(this._vf.enabled)
{this._vf.isActive=true;this.postMessage('isActive',true);}},pointerMoved:function(event)
{},pointerMovedOver:function(event)
{if(this._vf.enabled)
{this.postMessage('isOver',true);}},pointerMovedOut:function(event)
{if(this._vf.enabled)
{this.postMessage('isOver',false);}},pointerReleased:function()
{if(this._vf.enabled)
{this._vf.isActive=false;this.postMessage('isActive',false);}}}));x3dom.registerNodeType("X3DDragSensorNode","PointingDeviceSensor",defineClass(x3dom.nodeTypes.X3DPointingDeviceSensorNode,function(ctx)
{x3dom.nodeTypes.X3DDragSensorNode.superClass.call(this,ctx);this.addField_SFBool(ctx,'autoOffset',true);this._lastX=-1;this._lastY=-1;},{pointerPressedOverSibling:function(event)
{x3dom.nodeTypes.X3DPointingDeviceSensorNode.prototype.pointerPressedOverSibling.call(this,event);this._lastX=event.layerX;this._lastY=event.layerY;this._startDragging(event.viewarea,event.layerX,event.layerX,event.worldX,event.worldY,event.worldZ);},pointerMoved:function(event)
{x3dom.nodeTypes.X3DPointingDeviceSensorNode.prototype.pointerMoved.call(this,event);if(this._vf.isActive&&this._vf.enabled)
{this._process2DDrag(event.layerX,event.layerY,event.layerX-this._lastX,event.layerY-this._lastY);}},pointerReleased:function()
{x3dom.nodeTypes.X3DPointingDeviceSensorNode.prototype.pointerReleased.call(this);this._stopDragging();},_startDragging:function(viewarea,x,y,wx,wy,wz)
{},_process2DDrag:function(x,y,dx,dy)
{},_stopDragging:function()
{}}));x3dom.registerNodeType("X3DTouchSensorNode","PointingDeviceSensor",defineClass(x3dom.nodeTypes.X3DPointingDeviceSensorNode,function(ctx)
{x3dom.nodeTypes.X3DTouchSensorNode.superClass.call(this,ctx);},{}));x3dom.registerNodeType("TouchSensor","PointingDeviceSensor",defineClass(x3dom.nodeTypes.X3DTouchSensorNode,function(ctx)
{x3dom.nodeTypes.TouchSensor.superClass.call(this,ctx);},{}));x3dom.registerNodeType("PlaneSensor","PointingDeviceSensor",defineClass(x3dom.nodeTypes.X3DDragSensorNode,function(ctx)
{x3dom.nodeTypes.PlaneSensor.superClass.call(this,ctx);this.addField_SFRotation(ctx,'axisRotation',0,0,1,0);this.addField_SFVec2f(ctx,'minPosition',0,0);this.addField_SFVec2f(ctx,'maxPosition',-1,-1);this.addField_SFVec3f(ctx,'offset',0,0,0);this._rotationMatrix=this._vf.axisRotation.toMatrix();this._worldToLocalMatrix=null;this._initialPlaneIntersection=null;this._planeNormal=null;this._viewArea=null;this._currentTranslation=new x3dom.fields.SFVec3f(0.0,0.0,0.0);},{getCurrentTransform:function()
{var parentTransform=x3dom.nodeTypes.X3DDragSensorNode.prototype.getCurrentTransform.call(this);return parentTransform.mult(this._rotationMatrix);},_startDragging:function(viewarea,x,y,wx,wy,wz)
{x3dom.nodeTypes.X3DDragSensorNode.prototype._startDragging.call(this,viewarea,x,y,wx,wy,wz);this._viewArea=viewarea;this._currentTranslation=new x3dom.fields.SFVec3f(0.0,0.0,0.0).add(this._vf.offset);this._worldToLocalMatrix=this.getCurrentTransform().inverse();this._initialPlaneIntersection=this._worldToLocalMatrix.multMatrixPnt(new x3dom.fields.SFVec3f(wx,wy,wz));this._planeNormal=new x3dom.fields.SFVec3f(0.0,0.0,1.0);},_process2DDrag:function(x,y,dx,dy)
{x3dom.nodeTypes.X3DDragSensorNode.prototype._process2DDrag.call(this,x,y,dx,dy);var intersectionPoint;var minPos,maxPos;if(this._initialPlaneIntersection)
{var viewRay=this._viewArea.calcViewRay(x,y);viewRay.pos=this._worldToLocalMatrix.multMatrixPnt(viewRay.pos);viewRay.dir=this._worldToLocalMatrix.multMatrixVec(viewRay.dir.normalize());intersectionPoint=viewRay.intersectPlane(this._initialPlaneIntersection,this._planeNormal);if(!intersectionPoint)
{intersectionPoint=viewRay.intersectPlane(this._initialPlaneIntersection,this._planeNormal.negate());}
if(intersectionPoint)
{this._currentTranslation=intersectionPoint.subtract(this._initialPlaneIntersection);this._currentTranslation=this._currentTranslation.add(this._vf.offset);minPos=this._vf.minPosition;maxPos=this._vf.maxPosition;if(minPos.x<=maxPos.x)
{this._currentTranslation.x=Math.min(this._currentTranslation.x,maxPos.x);this._currentTranslation.x=Math.max(this._currentTranslation.x,minPos.x);}
if(minPos.y<=maxPos.y)
{this._currentTranslation.y=Math.min(this._currentTranslation.y,maxPos.y);this._currentTranslation.y=Math.max(this._currentTranslation.y,minPos.y);}
this.postMessage('translation_changed',x3dom.fields.SFVec3f.copy(this._currentTranslation));}}},_stopDragging:function()
{x3dom.nodeTypes.X3DDragSensorNode.prototype._stopDragging.call(this);if(this._vf.autoOffset)
{this._vf.offset=x3dom.fields.SFVec3f.copy(this._currentTranslation);this.postMessage('offset_changed',this._vf.offset);}}}));x3dom.registerNodeType("SphereSensor","PointingDeviceSensor",defineClass(x3dom.nodeTypes.X3DDragSensorNode,function(ctx)
{x3dom.nodeTypes.SphereSensor.superClass.call(this,ctx);this.addField_SFRotation(ctx,'offset',0,1,0,0);this._currentRotation=null;this._rotationMatrix=this._vf.offset.toMatrix();},{getCurrentTransform:function()
{var parentTransform=x3dom.nodeTypes.X3DDragSensorNode.prototype.getCurrentTransform.call(this);return parentTransform.mult(this._rotationMatrix);},_startDragging:function(viewarea,x,y,wx,wy,wz)
{x3dom.nodeTypes.X3DDragSensorNode.prototype._startDragging.call(this,viewarea,x,y,wx,wy,wz);this._currentRotation=new x3dom.fields.Quaternion();this._viewArea=viewarea;this._localOrigin=new x3dom.fields.SFVec3f(0.0,0.0,0.0);this._inverseToWorldMatrix=this.getCurrentTransform().inverse();var firstIntersection=this._inverseToWorldMatrix.multMatrixPnt(new x3dom.fields.SFVec3f(wx,wy,wz));this._initialSphereIntersectionVector=firstIntersection.subtract(this._localOrigin);this._sphereRadius=this._initialSphereIntersectionVector.length();this._initialSphereIntersectionVector=this._initialSphereIntersectionVector.normalize();},_process2DDrag:function(x,y,dx,dy)
{x3dom.nodeTypes.X3DDragSensorNode.prototype._process2DDrag.call(this,x,y,dx,dy);var viewRay=this._viewArea.calcViewRay(x,y);viewRay.pos=this._inverseToWorldMatrix.multMatrixPnt(viewRay.pos);viewRay.dir=this._inverseToWorldMatrix.multMatrixVec(viewRay.dir);var A=viewRay.dir.dot(viewRay.dir);var B=2.0*(viewRay.dir.dot(viewRay.pos.subtract(this._localOrigin)));var C=viewRay.pos.dot(viewRay.pos)-2.0*this._localOrigin.dot(viewRay.pos)+
this._localOrigin.dot(this._localOrigin)-this._sphereRadius*this._sphereRadius;var determinant=(B*B)-(4.0*A*C);var alpha_1;var alpha_2;if(determinant>=0.0){alpha_1=(-B+Math.sqrt(determinant))/(2.0*A);alpha_2=(-B-Math.sqrt(determinant))/(2.0*A);alpha_1=Math.min(alpha_1,alpha_2);if(alpha_1>=1.0){var hitPoint=viewRay.pos.add(viewRay.dir.multiply(alpha_1));var vecToHitPoint=hitPoint.subtract(this._localOrigin).normalize();this._currentRotation=x3dom.fields.Quaternion.rotateFromTo(this._initialSphereIntersectionVector,vecToHitPoint);this._currentRotation=this._currentRotation.multiply(this._vf.offset);this.postMessage('rotation_changed',this._currentRotation);}}
else{}},_stopDragging:function()
{x3dom.nodeTypes.X3DDragSensorNode.prototype._stopDragging.call(this);if(this._vf.autoOffset)
{this._vf.offset=this._currentRotation;this.postMessage('offset_changed',this._vf.offset);}
this._currentRotation=new x3dom.fields.Quaternion();}}));x3dom.registerNodeType("CylinderSensor","PointingDeviceSensor",defineClass(x3dom.nodeTypes.X3DDragSensorNode,function(ctx)
{x3dom.nodeTypes.CylinderSensor.superClass.call(this,ctx);this.addField_SFFloat(ctx,'offset',0);this.addField_SFRotation(ctx,'axisRotation',0,1,0,0);this.addField_SFFloat(ctx,'diskAngle',0.262);this.addField_SFFloat(ctx,'minAngle',0);this.addField_SFFloat(ctx,'maxAngle',-1);this._rotationMatrix=this._vf.axisRotation.toMatrix();this._inverseToWorldMatrix=null;this._initialCylinderIntersectionVector=null;this._viewArea=null;this._cylinderRadius=0.0;this._yAxisLine=null;this._cylinderMode=true;this._currentRotationAngle=0.0;},{getCurrentTransform:function()
{var parentTransform=x3dom.nodeTypes.X3DDragSensorNode.prototype.getCurrentTransform.call(this);return parentTransform.mult(this._rotationMatrix);},_startDragging:function(viewarea,x,y,wx,wy,wz)
{x3dom.nodeTypes.X3DDragSensorNode.prototype._startDragging.call(this,viewarea,x,y,wx,wy,wz);this._currentRotation=new x3dom.fields.Quaternion();this._viewArea=viewarea;this._yAxisLine=new x3dom.fields.Line(new x3dom.fields.SFVec3f(0.0,0.0,0.0),new x3dom.fields.SFVec3f(0.0,1.0,0.0));this._inverseToWorldMatrix=this.getCurrentTransform().inverse();var firstIntersection=this._inverseToWorldMatrix.multMatrixPnt(new x3dom.fields.SFVec3f(wx,wy,wz));var closestPointOnYAxis=this._yAxisLine.closestPoint(firstIntersection);this._initialCylinderIntersectionVector=firstIntersection.subtract(closestPointOnYAxis);this._cylinderRadius=this._initialCylinderIntersectionVector.length();this._initialCylinderIntersectionVector=this._initialCylinderIntersectionVector.normalize();},_process2DDrag:function(x,y,dx,dy)
{x3dom.nodeTypes.X3DDragSensorNode.prototype._process2DDrag.call(this,x,y,dx,dy);if(this._cylinderMode)
{var viewRay=this._viewArea.calcViewRay(x,y);viewRay.pos=this._inverseToWorldMatrix.multMatrixPnt(viewRay.pos);viewRay.dir=this._inverseToWorldMatrix.multMatrixVec(viewRay.dir);var A=viewRay.dir.subtract(this._yAxisLine.dir.multiply(viewRay.dir.dot(this._yAxisLine.dir)));var B=viewRay.pos.subtract(this._yAxisLine.pos).add(this._yAxisLine.dir.multiply(this._yAxisLine.dir.dot(this._yAxisLine.pos.subtract(viewRay.pos))));var p=2*A.dot(B)/A.dot(A);var q=(B.dot(B)-this._cylinderRadius*this._cylinderRadius)/A.dot(A);var sqrt_part=p*p*0.25-q;var alpha_1;var alpha_2;if(sqrt_part>=0)
{sqrt_part=Math.sqrt(sqrt_part);alpha_1=-p*0.5+sqrt_part;alpha_2=-p*0.5-sqrt_part;alpha_1=Math.min(alpha_1,alpha_2);if(alpha_1>0.0)
{var hitPoint=viewRay.pos.add(viewRay.dir.multiply(alpha_1));var closestPointOnYAxis=this._yAxisLine.closestPoint(hitPoint);var vecToHitPoint=hitPoint.subtract(closestPointOnYAxis).normalize();this._currentRotation=x3dom.fields.Quaternion.rotateFromTo(this._initialCylinderIntersectionVector,vecToHitPoint);var offsetQuat=x3dom.fields.Quaternion.axisAngle(this._yAxisLine.dir,this._vf.offset);this._currentRotation=this._currentRotation.multiply(offsetQuat);this.postMessage('rotation_changed',this._currentRotation);}}}
else
{}},_stopDragging:function()
{x3dom.nodeTypes.X3DDragSensorNode.prototype._stopDragging.call(this);if(this._vf.autoOffset)
{this._vf.offset=this._currentRotation.angle();this.postMessage('offset_changed',this._vf.offset);}}}));x3dom.versionInfo={version:'1.7.2',revision:'61a235203deb34329fe615cbbf21314db6ebf49f',date:'Mon Dec 19 19:17:05 2016 +0100'};