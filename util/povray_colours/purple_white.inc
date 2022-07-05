
//35
#declare r1 = 50/255;
#declare g1 = 50/255;
#declare b1 = 50/255;

//50
#declare r2 = 50/255;
#declare g2 = 20/255;
#declare b2 = 75/255;

//75
#declare r3 = 74/255;
#declare g3 = 32/255;
#declare b3 = 113/255;

//100
#declare r4 = 102/255;
#declare g4 = 44/255;
#declare b4 = 146/255;

//125
#declare r5 = 154/255;
#declare g5 = 67/255;
#declare b5 = 139/255;

//150
#declare r6 = 207/255;
#declare g6 = 92/255;
#declare b6 = 133/255;

//175
#declare r7 = 222/255;
#declare g7 = 137/255;
#declare b7 = 150/255;

//200
#declare r8 = 237/255;
#declare g8 = 189/255;
#declare b8 = 169/255;

//235
#declare r9 = 255/255;
#declare g9 = 255/255;
#declare b9 = 255/255;

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
