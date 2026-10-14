#macro spinrgb(sx,sy,sz,cr,cg,cb)
#declare ireal = sz;
#declare ssred=-ireal*0.5+0.1;
#declare ssblue=ireal*0.5+0.2;
#declare ssgreen=0.0;
#switch(ireal)
#range(-1.0,-0.9)
  #declare ssred=(4)/254;
  #declare ssgreen=(35)/254;
  #declare ssblue=(59)/254;
  #break
#range(-0.9,-0.7)
  #declare frac=(0.9+ireal)/0.2;
  #declare ssred=(frac*(41-4)+4)/254;
  #declare ssgreen=(frac*(62-35)+35)/254;
  #declare ssblue=(frac*(77-59)+59)/254;
  #break
#range(-0.7,-0.5)
  #declare frac=(0.7+ireal)/0.2;
  #declare ssred=(frac*(32-41)+41)/254;
  #declare ssgreen=(frac*(84-62)+62)/254;
  #declare ssblue=(frac*(124-77)+77)/254;
  #break
#range(-0.5,-0.3)
  #declare frac=(0.5+ireal)/0.2;
  #declare ssred=(frac*(119-32)+32)/254;
  #declare ssgreen=(frac*(177-84)+84)/254;
  #declare ssblue=(frac*(222-124)+124)/254;
  #break
#range(-0.3,-0.1)
  #declare frac=(0.3+ireal)/0.2;
  #declare ssred=(frac*(169-119)+119)/254;
  #declare ssgreen=(frac*(199-177)+177)/254;
  #declare ssblue=(frac*(222-222)+222)/254;
  #break
#range(-0.1,0.1)
  #declare frac=(0.1+ireal)/0.2;
  #declare ssred=(frac*(239-169)+169)/254;
  #declare ssgreen=(frac*(213-199)+199)/254;
  #declare ssblue=(frac*(178-222)+222)/254;
  #break
#range(0.1,0.3)
  #declare frac=(ireal-0.1)/0.2;
  #declare ssred=(frac*(239-239)+239)/254;
  #declare ssgreen=(frac*(189-213)+213)/254;
  #declare ssblue=(frac*(121-178)+178)/254;
  #break
#range(0.3,0.5)
  #declare frac=(ireal-0.3)/0.2;
  #declare ssred=(frac*(192-239)+239)/254;
  #declare ssgreen=(frac*(127-189)+189)/254;
  #declare ssblue=(frac*(40-121)+121)/254;
  #break
#range(0.5,0.7)
  #declare frac=(ireal-0.5)/0.2;
  #declare ssred=(frac*(120-192)+192)/254;
  #declare ssgreen=(frac*(95-127)+127)/254;
  #declare ssblue=(frac*(61-40)+40)/254;
  #break
#range(0.7,0.9)
  #declare frac=(ireal-0.7)/0.2;
  #declare ssred=(frac*(91-120)+120)/254;
  #declare ssgreen=(frac*(54-95)+95)/254;
  #declare ssblue=(frac*(5-61)+61)/254;
  #break
#range(0.9,1.0)
  #declare ssred=(91)/254;
  #declare ssgreen=(54)/254;
  #declare ssblue=(5)/254;
  #break
#else
  #declare ssred = 1.0;
  #declare ssgreen = 1.0;
  #declare ssblue = 1.0;
#end
#if(ssblue<0.0) #declare ssblue=0.0; #end
#if(ssred<0.0) #declare ssred=0.0; #end
#if(ssgreen<0.0) #declare ssgreen=0.0; #end
  pigment {color rgb <ssred ssgreen ssblue>}
#end
