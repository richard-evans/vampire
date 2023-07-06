
//0,0,144
#declare r1 = 0/255;
#declare g1 = 0/255;
#declare b1 = 144/255;

//0,15,255
#declare r2 = 0/255;
#declare g2 = 15/255;
#declare b2 = 255/255;

//0,144,255
#declare r3 = 0/255;
#declare g3 = 144/255;
#declare b3 = 255/255;

//15,255,238
#declare r4 = 15/255;
#declare g4 = 255/255;
#declare b4 = 238/255;

//144,255,112
#declare r5 = 144/255;
#declare g5 = 255/255;
#declare b5 = 112/255;

//255,238,0
#declare r6 = 255/255;
#declare g6 = 238/255;
#declare b6 = 0/255;

//255,112,0
#declare r7 = 255/255;
#declare g7 = 112/255;
#declare b7 = 0/255;

//238,0,0
#declare r8 = 238/255;
#declare g8 = 0/255;
#declare b8 = 0/255;

//127,0,0
#declare r9 = 127/255;
#declare g9 = 0/255;
#declare b9 = 0/255;

#macro spinrgb(sx,sy,sz,cr,cg,cb)
#declare ireal = sz*0.5+0.5;
#declare ssred=-ireal*0.5;
#declare ssblue=ireal*0.5;
#declare ssgreen=0.0;
#switch(ireal)
#range(0,0.1)
  #declare ssred=r9;
  #declare ssgreen=g9;
  #declare ssblue=b9;
  #break
#range(0.1,0.2)
  #declare frac=(ireal-0.1)/0.1;
  #declare ssred=(frac*(r8-r9)+r9);
  #declare ssgreen=(frac*(g8-g9)+g9);
  #declare ssblue=(frac*(b8-b9)+b9);
  #break
#range(0.2,0.3)
  #declare frac=(ireal-0.2)/0.1;
  #declare ssred=(frac*(r7-r8)+r8);
  #declare ssgreen=(frac*(g7-g8)+g8);
  #declare ssblue=(frac*(b7-b8)+b8);
  #break
#range(0.3,0.4)
  #declare frac=(ireal-0.3)/0.1;
  #declare ssred=(frac*(r6-r7)+r7);
  #declare ssgreen=(frac*(g6-g7)+g7);
  #declare ssblue=(frac*(b6-b7)+b7);
  #break
#range(0.4,0.5)
  #declare frac=(ireal-0.4)/0.1;
  #declare ssred=(frac*(r5-r6)+r6);
  #declare ssgreen=(frac*(g5-g6)+g6);
  #declare ssblue=(frac*(b5-b6)+b6);
  #break
#range(0.5,0.6)
  #declare frac=(ireal-0.5)/0.1;
  #declare ssred=(frac*(r4-r5)+r5);
  #declare ssgreen=(frac*(g4-g5)+g5);
  #declare ssblue=(frac*(b4-b5)+b5);
  #break
#range(0.6,0.7)
  #declare frac=(ireal-0.6)/0.1;
  #declare ssred=(frac*(r3-r4)+r4);
  #declare ssgreen=(frac*(g3-g4)+g4);
  #declare ssblue=(frac*(b3-b4)+b4);
  #break
#range(0.7,0.8)
  #declare frac=(ireal-0.7)/0.1;
  #declare ssred=(frac*(r2-r3)+r3);
  #declare ssgreen=(frac*(g2-g3)+g3);
  #declare ssblue=(frac*(b2-b3)+b3);
  #break
#range(0.8,0.9)
  #declare frac=(ireal-0.8)/0.1;
  #declare ssred=(frac*(r1-r2)+r2);
  #declare ssgreen=(frac*(g1-g2)+g2);
  #declare ssblue=(frac*(b1-b2)+b2);
  #break
#range(0.9,1.0)
  #declare ssred=r1;
  #declare ssgreen=g1;
  #declare ssblue=b1;
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
