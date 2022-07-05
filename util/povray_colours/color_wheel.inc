//----------------------------------------------------------------
//
//  POVRAY Macro to generate 3D colour space based on angles
//  from the x-axis (theta) and z (phi).
//
//  (c) R F L Evans 2014
//
//-----------------------------------------------------------------
#macro colour_wheel(theta, phi) // theta from x, phi from z
//-------------------------------
// Declare r,g,b components
//-------------------------------
#declare rr=1.0;
#declare gg=1.0;
#declare bb=1.0;
//-------------------
// Determine colour
//-------------------
#switch(theta)

#range(0,120)
  #declare sf=2.-abs(60-theta)/60.0;
  #declare rr=sf*(1-(theta/120.0));
  #declare gg=sf*((theta/120.0));
  #declare bb=0.0;
#break

#range(120,240)
  #declare sf=2.-abs(180-theta)/60.0;
  #declare rr=0.0;
  #declare gg=sf*(1-((theta-120)/120.0));
  #declare bb=sf*(theta-120)/120.0;
#break

#range(240,360)
  #declare sf=2.0-abs(300-theta)/60.0;
  #declare rr=sf*(theta-240.)/120.0;
  #declare gg=0.0;
  #declare bb=sf*(1-((theta-240)/120.0));
#end

//---------------------
// Determine grading
//---------------------
//#declare factor=pow(cos(phi*pi/180),2);
//#declare factor=pow(cos(phi*pi/180),3);
#declare factor=cos(phi*pi/180);

#switch(phi)
  #range(0,90)
    #declare factor=cos(phi*pi/180);
    #switch(rr)
    #range(0,factor)
      #declare rr=factor;
    #end

    #switch(gg)
    #range(0,factor)
      #declare gg=factor;
    #end

    #switch(bb)
    #range(0,factor)
      #declare bb=factor;
    #end
  #break

  #range(90,180)
    #declare factor=1.+cos(phi*pi/180);
    //#debug concat("Phi: ",str(phi,5,1),"\n")
    //#debug concat("factor: ",str(factor,5,1),"\n")
    #switch(rr)
    #range(factor,10)
      #declare rr=factor;
    #end

    #switch(gg)
    #range(factor,10)
      #declare gg=factor;
    #end

    #switch(bb)
    #range(factor,10)
      #declare bb=factor;
    #end
    //#debug concat("red   : ",str(rr,5,2),"\n")
    //#debug concat("green : ",str(gg,5,2),"\n")
    //#debug concat("blue  : ",str(bb,5,2),"\n")
  #end

  // Ensure colours in range 0 -> 1
  #switch(rr)
  #range(1,10)
    #declare rr=1.0;
  #end

  #switch(gg)
  #range(1,10)
    #declare gg=1.0;
  #end

  #switch(bb)
  #range(1,10)
    #declare bb=1.0;
  #end


//---------------------
// Set final pigment
//---------------------
pigment {color rgb <rr gg bb>}

#end // End macro

//-------------------------------------------
// Macro for povray spins
//-------------------------------------------
#macro spinrgb(sx,sy,sz,cr,cg,cb)
   #declare theta = acos(sx/(sqrt(sx*sx+sy*sy)));
   #declare phi = acos(sz);
   colour_wheel(theta, phi)
#end
