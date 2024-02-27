#include "create.hpp"
#include "errors.hpp"
#include "grains.hpp"
#include "material.hpp"
#include "qvoronoi.hpp"
#include "random.hpp"
#include "vio.hpp"
#include "vmpi.hpp"
#include "vmath.hpp"
#include "voro++.hh"

#include <cmath>
#include <list>
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <algorithm>

#include "internal.hpp"

namespace create_radical_voronoi{
	bool rounded=false;
	double voronoi_sd=0.15;			/// Standard Deviation of voronoi grains
}

namespace cs{

bool Drop_and_Roll(double HozMin,double HozMax,double VerMin,double VerMax,
		std::mt19937_64 VoroGen,
		  std::vector<double>*Hoz,std::vector<double>*Ver,std::vector<double>*r){

	// TODO: Implement AABB tree to improve performance

	 double dTh = ((2.0*M_PI)/1000.0);
	 double Boundary_check = (HozMax-HozMin)*0.5;
	 double dVert= ((VerMax-VerMin)/1000.0);
	 std::uniform_real_distribution<double> UNI(0.0,HozMax);
	 double HozPos=UNI(VoroGen);
	 double VerPos=VerMax+r->back();
	 // Force particle within a box and set starting height
	 if(HozPos-r->back()<HozMin){HozPos=HozMin+r->back();}
	 else if(HozPos+r->back()>HozMax){HozPos=HozMax-r->back();}
	 bool placed=false;
	 bool falling=true;
	 bool rolling=false;
	 unsigned int contact=0;

	 while(!placed){
		  double HozCont=0.0;
		  double VerCont=0.0;
		  double rCont=0.0;
		  while(falling){
				// Drop
				VerPos -= dVert;
				// Check -- NB: for improved efficiency split box into grid for coarse check before fine
				for(unsigned int neigh=0; neigh<(Ver->size()-1);++neigh){
					 double HozNeigh=Hoz->at(neigh);
					 double VerNeigh=Ver->at(neigh);
					 double RNeigh=r->at(neigh);
					 double Hdist	 = (HozPos-HozNeigh);
					 double radSumSq = (r->back()+RNeigh)*(r->back()+RNeigh);

					 //########### Periodic Boundary Checks ###########//
					 if(Hdist>Boundary_check){HozNeigh+=HozMax;}// Move neigh to the right
					 else if(Hdist<-Boundary_check){HozNeigh-=HozMax;}// Move neigh to the left
					 //################################################//

					 double distSq = (HozPos-HozNeigh)*(HozPos-HozNeigh)
										+ (VerPos-VerNeigh)*(VerPos-VerNeigh);

					 if(distSq==radSumSq || distSq<(radSumSq-r->back()*0.1)){ // Touching or intersecting
						  HozCont=HozNeigh;VerCont=VerNeigh;rCont=RNeigh;
						  falling=false;
						  rolling=true;
						  contact=neigh;
						  goto Rolling_label;
					 }
				}
				// Hit bottom
				if((VerPos-r->back())<=VerMin){ // Good improvement here would be maybe round up the comparison so it counts as touching when ~1e-12
					 placed=true;
					 falling=false;
				}
		  }
Rolling_label:
		  while(rolling){
				// Particle exactly on top
				if(HozPos==HozCont){
					 placed=true;
					 falling=true;
					 break;
				}
				// Particle grazing past
				if(VerPos<=VerCont){
					 rolling=false;
					 falling=true;
					 break;
				}
				// Roll the particle
				// Determine incident angle
				double angle = atan2((VerPos-VerCont),(HozPos-HozCont));
				if(angle>1.570796){ //Roll left
					 while(angle<=M_PI && !placed){
						  HozPos = HozCont + (r->back()+rCont)*cos(angle);
						  VerPos = VerCont + (r->back()+rCont)*sin(angle);
						  // Check contact
						  for(unsigned int neigh=0; neigh<(Ver->size()-1);++neigh){
								if(neigh==contact){continue;}
								double HozNeigh=Hoz->at(neigh);
								double VerNeigh=Ver->at(neigh);
								double RNeigh=r->at(neigh);

								double Hdist	  = (HozPos-HozNeigh);
								double radSumSq = (r->back()+RNeigh)*(r->back()+RNeigh);

								//########### Periodic Boundary Checks ###########//
								if(Hdist>Boundary_check){HozNeigh+=HozMax;}// Move neigh to the right
								else if(Hdist<-Boundary_check){HozNeigh-=HozMax;}// Move neigh to the left
								//################################################//

								double distSq	= (HozPos-HozNeigh)*(HozPos-HozNeigh)
													 + (VerPos-VerNeigh)*(VerPos-VerNeigh);

								if(distSq<radSumSq){ // We have hit another disc
									 // Is current disc's CoM in between the two contacting discs?
									 double GAP_LIMIT=1.0e-9; // Set a gap limit to prevent infinite loop where the gap is very close but not exactly equal
									 if( (fabs(HozNeigh-HozPos)+fabs(HozCont-HozPos) - fabs(HozNeigh-HozCont)) < GAP_LIMIT ){
										  placed=true;
										  rolling=false;
										  break;
									 }
									 // Disc should now Roll around Neigh NOT Cont
									 HozCont=HozNeigh;VerCont=VerNeigh;rCont=RNeigh;
									 contact=neigh;
									 goto Rolling_label; //NOLINT
								}
						  }
						  if(!placed){
								// Check for floor
								if((VerPos-r->back())<=VerMin){
									 placed=true;
									 rolling=false;
									 break;
								}
								angle += dTh;
						  }
					 }
					 if(!placed){
						  rolling=false;
						  falling=true;
					 }
				}
				else{ // Roll right
					 while(angle>=0.0 && !placed){
						  HozPos = HozCont + (r->back()+rCont)*cos(angle);
						  VerPos = VerCont + (r->back()+rCont)*sin(angle);
						  for(unsigned int neigh=0; neigh<(Ver->size()-1);++neigh){
								if(neigh==contact){continue;}
								double HozNeigh=Hoz->at(neigh);
								double VerNeigh=Ver->at(neigh);
								double RNeigh=r->at(neigh);

								double Hdist	 = (HozPos-HozNeigh);
								double radSumSq = (r->back()+RNeigh)*(r->back()+RNeigh);

								//########### Periodic Boundary Checks ###########//
								if(Hdist>Boundary_check){HozNeigh+=HozMax;}// Move neigh to the right
								else if(Hdist<-Boundary_check){HozNeigh-=HozMax;}// Move neigh to the left
								//################################################//

								double distSq	= (HozPos-HozNeigh)*(HozPos-HozNeigh)
													 + (VerPos-VerNeigh)*(VerPos-VerNeigh);

								if(distSq<radSumSq){
									 // Is current disc's CoM in between the two contacting discs?
								double GAP_LIMIT=1.0e-9; // Set a gap limit to prevent infinite loop where the gap is very close but not exactly equal
									 if( (fabs(HozNeigh-HozPos)+fabs(HozCont-HozPos) - fabs(HozNeigh-HozCont)) < GAP_LIMIT ){
										  placed=true;
										  rolling=false;
										  break;
									 }
									 // Disc should now Roll around Neigh NOT Cont
									 HozCont=HozNeigh;VerCont=VerNeigh;rCont=RNeigh;
									 contact=neigh;
									 goto Rolling_label; //NOLINT
								}
						  }
						  if(!placed){ // No collision with other particles so check for floor collision
								// Check for floor
								if((VerPos-r->back())<=VerMin){
									 placed=true;
									 rolling=false;
									 break;
								}
								angle -= dTh;
						  }
					 }
					 if(!placed){
						  rolling=false;
						  falling=true;
					 }
				}
		  }
		  // Check that we aren't full
		  if(VerPos+r->back()>VerMax){
				placed=false;
				break;
		  }
	 }

	 Hoz->back()=HozPos, Ver->back()=VerPos;
	return placed;
}

int voronoi_radical_film(std::vector<cs::catom_t> & catom_array){

	// check calling of routine if error checking is activated
	if(err::check==true){
		terminaltextcolor(RED);
		std::cerr << "cs::voronoi_radical_film has been called" << std::endl;
		terminaltextcolor(WHITE);
	}
	// We start by packing our system with circles
	// TODO: Implement AABB Tree to improve performance

	double grain_sd=create_radical_voronoi::voronoi_sd;
	std::vector<double> X, Y, r;
	std::mt19937_64 VoroGen;
	VoroGen.seed(mtrandom::voronoi_seed);
	double TinyGrainChance = create::internal::voronoi_tiny_grain_chance;
	//double Max_width_limit = 100; // Used to prevent rare occurance of very large grains
	double Grain_width = create::internal::voronoi_grain_size;
	double mulg = log(create::internal::voronoi_grain_size);
	double Grain_spacing = create::internal::voronoi_grain_spacing;
	double Packing_fraction = 1.0 - Grain_spacing/Grain_width;
	double Dim_x=cs::system_dimensions[0];
	double Dim_y=cs::system_dimensions[1];
	double Max_x=Dim_x;
	double Max_y=Dim_y;

	if(grain_sd > 0){
		bool filled=false;
		int attempt=0;
		int attempt_limit=100;
		std::uniform_real_distribution<double> Tiny(0,1.0);
		std::uniform_real_distribution<double> TinyGrainSize(0.0,Grain_width);
		std::lognormal_distribution<double> LN(mulg,grain_sd);

		unsigned int Disc=0;
		while(!filled){
			if(attempt==attempt_limit){filled=true;}

			X.push_back(0.0);
			Y.push_back(0.0);
			double radius=0.0;
			if(Tiny(VoroGen) < TinyGrainChance){
				radius = TinyGrainSize(VoroGen);
			}
			else{
				do{radius=LN(VoroGen)*0.5;}
				while(/*radius > Max_width_limit*(Grain_width*0.5) ||*/ radius < 0.0);
			}

			r.push_back(radius);
			bool Placed=false;
			std::vector<double> Hoz(X.size(),0.0);
			std::vector<double> Ver(X.size(),0.0);

			Hoz=X;
			Ver=Y;
			Placed = Drop_and_Roll(0.0,Max_x,0.0,Max_y+Grain_width*10,VoroGen,&Hoz,&Ver,&r);

			if(Placed){
				attempt=0;
				X.back()=Hoz.back();
				Y.back()=Ver.back();
			}else{
				--Disc;
				++attempt;
				X.pop_back();
				Y.pop_back();
				r.pop_back();
			}
			++Disc;
		}
		  // Cut into system
		  std::vector<double> Xtmp;
		  std::vector<double> Ytmp;
		  std::vector<double> rtmp;

		  for(size_t i=0; i<Y.size(); ++i){
				if(Y[i]>Grain_width*5.0 && Y[i]<(Max_y+Grain_width*5.0)){
					 Xtmp.push_back(X[i]);
					 Ytmp.emplace_back(Y[i]-Grain_width*5.0);
					 rtmp.push_back(r[i]);
				}
		  }
		  X=Xtmp;
		  Y=Ytmp;
		  r=rtmp;

	}else{
		double Local_Grain_width = Grain_width;
		double Local_Grain_height = (2.0/sqrt(3.0))/Local_Grain_width;
		int Bx=int(Dim_x/(Local_Grain_width));
		  int By=int(Dim_y/(1.75*Local_Grain_height));
		  Max_x = Bx*Local_Grain_width;
		  Max_y = 1.5*By*Local_Grain_height;
		  for(int j=0;j<By;++j){
				for(int i=0;i<Bx;++i){

					 double x = 0.5*Local_Grain_width  + Local_Grain_width*double(i);
					 double y = 0.5*Local_Grain_height + Local_Grain_height*(double(j)*1.5);
					 X.push_back(x);
					 Y.push_back(y);
					 r.push_back(Local_Grain_width);
					 x = 0.5*Local_Grain_width + Local_Grain_width*0.5 + Local_Grain_width*double(i);
					 y = 0.5*Local_Grain_height+Local_Grain_height*0.75 + Local_Grain_height*(double(j)*1.5);
					 X.push_back(x);
					 Y.push_back(y);
					 r.push_back(Local_Grain_width);
				}
		  }
	}
	double Block_size=0.0;
	if(Max_x<Max_y){
		Block_size = Max_x*0.1;
	}else{
		Block_size = Max_y*0.1;
	}
	const int Blocks_x = int(Dim_x/Block_size);
	const int Blocks_y = int(Dim_y/Block_size);
	const int Blocks_z = 1;

	voro::container_poly Container(0.0,Max_x,0.0,Max_y,-1.0,1.0,Blocks_x,Blocks_y,Blocks_z,true,true,false,8);
	voro::particle_order Particle_Order;

	for(uint Disc=0; Disc<r.size(); ++Disc){Container.put(Particle_Order,Disc,X[Disc],Y[Disc],0.0,r[Disc]);}
	X.clear();
	Y.clear();
	r.clear();

	voro::c_loop_order Container_loop(Container, Particle_Order);
	voro::voronoicell_neighbor Cell_with_neighbour_info;

	std::vector<unsigned int> FailedGrains; // Store the grains which cannot be computed
	double x, y, z;
	std::vector<double> tmp_vertices;
	std::vector<double> tmp_vx_hold, tmp_vy_hold;
	std::vector<int> Grain_ID;
	std::vector<double> tmp_x_hold, tmp_y_hold;
	std::vector<std::vector<double>> tmp_Vertex_X, tmp_Vertex_Y;
	std::vector<int> tmp_Grain_neighbours;
	std::vector<std::vector<int>> tmp_Neighbour;

	 if(Container_loop.start()){
		  do{
				if(Container.compute_cell(Cell_with_neighbour_info,Container_loop)){
					 Container_loop.pos(x,y,z);
					 Grain_ID.emplace_back(Container_loop.pid()-FailedGrains.size());

					 tmp_x_hold.push_back(x);tmp_y_hold.push_back(y);
					 Cell_with_neighbour_info.vertices(x,y,z,tmp_vertices);
					 // Voro++ works in 3D so here we want to set the z-coord to 1 (flatten the system to 2D)
					 for(size_t i=0; i<tmp_vertices.size(); i+=3){
						  if(tmp_vertices[i+2] == 1){
								tmp_vx_hold.push_back(tmp_vertices[i]);
								tmp_vy_hold.push_back(tmp_vertices[i+1]);
						  }
					 }
					 tmp_Vertex_X.push_back(tmp_vx_hold);
					 tmp_Vertex_Y.push_back(tmp_vy_hold);
					 tmp_vx_hold.clear();
					 tmp_vy_hold.clear();

					 std::vector<int> Neighbours;
					 Cell_with_neighbour_info.neighbors(tmp_Grain_neighbours);
					 for(int Gneigh : tmp_Grain_neighbours){
						  if(Gneigh>=0){
								Neighbours.push_back(Gneigh);
						  }
					 }
					 tmp_Neighbour.push_back(Neighbours);
				}
				else{
					 FailedGrains.emplace_back(Container_loop.pid());
				}
		  }
		  while(Container_loop.inc());
	 }

// ############################################
	size_t NumGrains = tmp_x_hold.size();
	tmp_Neighbour.resize(NumGrains);
	std::vector<std::vector<double>> Grain_Pos(NumGrains, std::vector<double>(2));
	std::vector<std::vector<std::vector<double>>> Grain_Vertices(NumGrains);
	std::vector<std::vector<int>> Grain_neigh_final(NumGrains);

//########################REORDER LISTS SO INDEX MATCHES GRAIN ID######################//
	for(size_t g=0; g<NumGrains; ++g){
		int ID = Grain_ID[g];
		Grain_Pos[g][0] = tmp_x_hold[g];
		Grain_Pos[g][1] = tmp_y_hold[g];

		Grain_Vertices[g].resize(tmp_Vertex_X[g].size(), std::vector<double>(2));
		for(size_t v=0; v<tmp_Vertex_X[g].size(); ++v){
			Grain_Vertices[g][v][0] = tmp_Vertex_X[g][v];
			Grain_Vertices[g][v][1] = tmp_Vertex_Y[g][v];
		}

		Grain_neigh_final[g].resize(tmp_Neighbour[g].size());
		for(size_t n=0; n<tmp_Neighbour[g].size(); ++n){
			Grain_neigh_final[g][n] = tmp_Neighbour[g][n];
		}
	}
    tmp_x_hold.clear();
    tmp_y_hold.clear();

//################################REORDER VERTEX LIST #################################//
	// We reorder the vertex list to enable determination of the grain area later
	std::vector<double> angle_hold;
	for(size_t g=0; g<NumGrains; ++g){
		size_t NumVert = Grain_Vertices[g].size();
		double Grain_x=Grain_Pos[g][0];
		double Grain_y=Grain_Pos[g][1];

		tmp_vx_hold.resize(NumVert);
		tmp_vy_hold.resize(NumVert);
		for(size_t v=0;v<NumVert;++v){
			tmp_x_hold.push_back(Grain_Vertices[g][v][0]);
			tmp_y_hold.push_back(Grain_Vertices[g][v][1]);
		}

		for(size_t v=0; v<NumVert; ++v){
			double Angle_temp = atan2((tmp_y_hold[v]-Grain_y),(tmp_x_hold[v]-Grain_x));
			if(Angle_temp<0){Angle_temp+=2.0*M_PI;}
			angle_hold.push_back(Angle_temp);
		}

		double min_pos=0;
		for(size_t j=0; j<angle_hold.size(); ++j){
			double min_angle=7.0;
			for(size_t k=0; k<angle_hold.size(); ++k){
				if(angle_hold[k]<min_angle){
					min_angle=angle_hold[k];
					min_pos=k;
			}    }
			tmp_vx_hold[j]=tmp_x_hold[min_pos];
			tmp_vy_hold[j]=tmp_y_hold[min_pos];
			angle_hold[min_pos]=7.0;
		}

		for(size_t j=0;j<tmp_vx_hold.size();++j){
		Grain_Vertices[g][j][0]=tmp_vx_hold[j];
		Grain_Vertices[g][j][1]=tmp_vy_hold[j];
		}
		tmp_x_hold.clear();
		tmp_y_hold.clear();
		tmp_vx_hold.clear();
		tmp_vy_hold.clear();
		angle_hold.clear();
	}

//#############################DETERMINE SYSTEM CENTRE#################################//
	double Vx=0.0;
	double Vy=0.0;
	double VxMAX=0.0;
	double VxMIN=0.0;
	double VyMAX=0.0;
	double VyMIN=0.0;
	std::vector<double> SystemCentre(2);
	for(size_t grain=0; grain<NumGrains; ++grain){
		for(size_t vertex=0;vertex<Grain_Vertices[grain].size();++vertex){
			Vx = Grain_Vertices[grain][vertex][0];
			Vy = Grain_Vertices[grain][vertex][1];
			if(Vx>VxMAX){VxMAX=Vx;}
			else if(Vx<VxMIN){VxMIN=Vx;}
			if(Vy>VyMAX){VyMAX=Vy;}
			else if(Vy<VyMIN){VyMIN=Vy;}
	}    }
	SystemCentre[0] = (VxMAX-VxMIN)/2.0+VxMIN;
	SystemCentre[1] = (VyMAX-VyMIN)/2.0+VyMIN;

//################################APPLY GRAIN SPACING##################################//
	for(size_t grain=0; grain<NumGrains; ++grain){
		double local_Px=Grain_Pos[grain][0];
		double local_Py=Grain_Pos[grain][1];
		for(size_t vert=0; vert<Grain_Vertices[grain].size(); ++vert){
			double local_Vx=0.0;
			double local_Vy=0.0;
			local_Vx = Grain_Vertices[grain][vert][0];
			local_Vy = Grain_Vertices[grain][vert][1];
			Grain_Vertices[grain][vert][0] = local_Px+(local_Vx-local_Px)*sqrt(Packing_fraction);
			Grain_Vertices[grain][vert][1] = local_Py+(local_Vy-local_Py)*sqrt(Packing_fraction);
		}
	}

//######################################### GRAIN AREA AND GEO CENTRE ##########################################//
	std::vector<double> Grain_Area(NumGrains);
	std::vector<std::vector<double>> Grain_Geo_Pos(NumGrains, std::vector<double>(2));
	for(size_t grain=0; grain<NumGrains; ++grain){
		double Area_temp=0.0;
		double CentroidX=0.0;
		double CentroidY=0.0;
		for(size_t i=0; i<Grain_Vertices[grain].size(); ++i){
			size_t k=(i+1)%Grain_Vertices[grain].size();

			Area_temp += (Grain_Vertices[grain][i][0]*Grain_Vertices[grain][k][1]);
			Area_temp -= (Grain_Vertices[grain][k][0]*Grain_Vertices[grain][i][1]);

			CentroidX += (Grain_Vertices[grain][i][0]+Grain_Vertices[grain][k][0])
                       * (Grain_Vertices[grain][i][0]*Grain_Vertices[grain][k][1]
                       -  Grain_Vertices[grain][k][0]*Grain_Vertices[grain][i][1]);

            CentroidY += (Grain_Vertices[grain][i][1]+Grain_Vertices[grain][k][1])
                       * (Grain_Vertices[grain][i][0]*Grain_Vertices[grain][k][1]
                       -  Grain_Vertices[grain][k][0]*Grain_Vertices[grain][i][1]);
		}
		Grain_Area[grain]=Area_temp*0.5;
		Grain_Geo_Pos[grain][0] = CentroidX/(6.0*Grain_Area[grain]);
		Grain_Geo_Pos[grain][1] = CentroidY/(6.0*Grain_Area[grain]);
	}
	Grain_Pos = Grain_Geo_Pos;

	// Shift all coordinates to ensure no negative values.
	double MinX=0.0;
	double MinY=0.0;
	for(size_t grain=0; grain<NumGrains; ++grain){
		for(size_t vert=0; vert<Grain_Vertices[grain].size(); ++vert){
			if(Grain_Vertices[grain][vert][0]<MinX){MinX=Grain_Vertices[grain][vert][0];}
			if(Grain_Vertices[grain][vert][1]<MinY){MinY=Grain_Vertices[grain][vert][1];}
		}
	}

	for(size_t grain=0; grain<NumGrains; ++grain){
		for(size_t vert=0; vert<Grain_Vertices[grain].size(); ++vert){
			Grain_Vertices[grain][vert][0] -= MinX;
			Grain_Vertices[grain][vert][1] -= MinY;
		}
		Grain_Pos[grain][0] -= MinX;
		Grain_Pos[grain][1] -= MinY;
	}

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	// Check for any grains which extend beyond the specified system dimensions and delete those grains
	std::vector<size_t> GrainToDelete;
	for(int g=NumGrains-1; g>=0; --g){
		// If grain pos beyond system then remove and skip next check
		if(Grain_Pos[g][0]>Dim_x || Grain_Pos[g][1]>Dim_y){
			GrainToDelete.push_back(g);
			continue;
		}
		// Check if any vertices are outside the system dimensions
		for(size_t v=0; v<Grain_Vertices[g].size(); ++v){
			if(Grain_Vertices[g][v][0]>Dim_x || Grain_Vertices[g][v][1]>Dim_y){
				GrainToDelete.push_back(g);
				break;
			}
		}
	}
	// Ensure we have no duplicates
	GrainToDelete.erase(std::unique(GrainToDelete.begin(), GrainToDelete.end()), GrainToDelete.end());

	for(size_t id : GrainToDelete){
		Grain_Pos.erase(Grain_Pos.begin()+id);
		Grain_Vertices.erase(Grain_Vertices.begin()+id);
		Grain_neigh_final.erase(Grain_neigh_final.begin()+id);
	}
	NumGrains = Grain_Pos.size();

	if(NumGrains==0){
		std::cerr << "No grains generated! Exiting..." << std::endl;
		exit(1);
	}

	// Organise data in the way that VAMPIRE expects
	// Reserve space for pointers
	std::vector<std::vector<double>> grain_coord_array(NumGrains, std::vector<double>(2));
	std::vector<std::vector<std::vector<double>>> grain_vertices_array(NumGrains);

	grain_coord_array = Grain_Pos;
	grain_vertices_array = Grain_Vertices;
	// Convert vertices to be relative to grain centre - required for rounding and supercell code
	for(size_t grain=0; grain<NumGrains; ++grain){
		for(size_t v=0;v<Grain_Vertices[grain].size();++v){
			grain_vertices_array[grain][v][0] -= grain_coord_array[grain][0];
			grain_vertices_array[grain][v][1] -= grain_coord_array[grain][1];
		}
	}
	if(create_radical_voronoi::rounded) create::internal::voronoi_grain_rounding(grain_coord_array, grain_vertices_array);
	// Create a 2D supercell array of atom numbers to improve performance for systems with many grains
	std::vector<std::vector<std::vector<int>>> supercell_array;

	int min_bounds[3];
	int max_bounds[3];

	min_bounds[0]=0;
	min_bounds[1]=0;
	min_bounds[2]=0;
	max_bounds[0]=cs::total_num_unit_cells[0];
	max_bounds[1]=cs::total_num_unit_cells[1];
	max_bounds[2]=cs::total_num_unit_cells[2];

	// allocate supercell array
	int dx = max_bounds[0]-min_bounds[0];
	int dy = max_bounds[1]-min_bounds[1];

	supercell_array.resize(dx);
	for(int i=0;i<dx;i++) supercell_array[i].resize(dy);

	// loop over atoms and populate supercell array
	for(unsigned int atom=0;atom<catom_array.size();atom++){
		int cx = static_cast<int>(catom_array[atom].x/unit_cell.dimensions[0]);
		int cy = static_cast<int>(catom_array[atom].y/unit_cell.dimensions[1]);
		supercell_array.at(cx).at(cy).push_back(atom);
	}

	// Determine order for core-shell grains
	std::list<create::internal::core_radius_t> material_order(0);
	for(int mat=0;mat<mp::num_materials;mat++){
		create::internal::core_radius_t tmp;
		tmp.mat=mat;
		tmp.radius=mp::material[mat].core_shell_size;
		material_order.push_back(tmp);
	}
	// sort by increasing radius
	material_order.sort(create::internal::compare_radius);

	// arrays to store list of grain vertices
	int max_vertices = 50;
	double tmp_grain_pointx_array[max_vertices];
	double tmp_grain_pointy_array[max_vertices];

	// loop over all grains with vertices
	for(unsigned int grain=0;grain<grain_coord_array.size();grain++){
		// Output progress indicator
		if((grain%(grain_coord_array.size()/10))==0){
		  std::cout << "." << std::flush;
		  zlog << "." << std::flush;
		}

		// Exclude grains with zero vertices
		if(grain_vertices_array[grain].size()!=0){

			// initialise minimum and max supercell coordinates for grain
			int minx=10000000;
			int maxx=0;
			int miny=10000000;
			int maxy=0;

			// Set temporary vertex coordinates (real) and compute cell ranges
			int num_vertices = grain_vertices_array[grain].size();
			for(int vertex=0;vertex<num_vertices;vertex++){
				// determine vertex coordinates
				tmp_grain_pointx_array[vertex]=grain_vertices_array[grain][vertex][0];
				tmp_grain_pointy_array[vertex]=grain_vertices_array[grain][vertex][1];
				// determine unit cell coordinates encompassed by grain
				int x = int((tmp_grain_pointx_array[vertex]+grain_coord_array[grain][0])/unit_cell.dimensions[0]);
				int y = int((tmp_grain_pointy_array[vertex]+grain_coord_array[grain][1])/unit_cell.dimensions[1]);
				if(x < minx) minx = x;
				if(x > maxx) maxx = x;
				if(y < miny) miny = y;
				if(y > maxy) maxy = y;
			}

			// determine coordinate offset for grains
			const double x0 = grain_coord_array[grain][0];
			const double y0 = grain_coord_array[grain][1];

			// loop over cells
			for(int i=minx;i<=maxx;i++){
				for(int j=miny;j<=maxy;j++){
					// loop over atoms in cells;
					for(unsigned int id=0;id<supercell_array[i][j].size();id++){
						int atom = supercell_array[i][j][id];
						// Get atomic position
						double x = catom_array[atom].x;
						double y = catom_array[atom].y;
						if(mp::material[catom_array[atom].material].core_shell_size>0.0){
							// Iterate over materials
							for(std::list<create::internal::core_radius_t>::iterator it = material_order.begin(); it !=  material_order.end(); it++){
								int mat = (it)->mat;
								double factor = mp::material[mat].core_shell_size;
								double maxz=create::internal::mp[mat].max*cs::system_dimensions[2];
								double minz=create::internal::mp[mat].min*cs::system_dimensions[2];
								double cz=catom_array[atom].z;
								const int atom_uc_cat = catom_array[atom].uc_category;
								const int mat_uc_cat = create::internal::mp[mat].unit_cell_category;
								if(vmath::point_in_polygon_factor(x-x0,y-y0,factor, tmp_grain_pointx_array,tmp_grain_pointy_array,num_vertices)==true){
									if((cz>=minz) && (cz<maxz) && (atom_uc_cat == mat_uc_cat) ){
										catom_array[atom].include=true;
										catom_array[atom].material=mat;
										catom_array[atom].grain=grain;
									}
									// if set to clear atoms then remove atoms within radius
									else if(cs::fill_core_shell==false){
										catom_array[atom].include=false;
									}
								}
							}
						}
						// Check to see if site is within polygon
						else if(vmath::point_in_polygon_factor(x-x0,y-y0,1.0,tmp_grain_pointx_array,tmp_grain_pointy_array,num_vertices)==true){
							catom_array[atom].include=true;
							catom_array[atom].grain=grain;
						}
					}
				}
			}
		}
	}
	terminaltextcolor(GREEN);
	std::cout << "done!" << std::endl;
	terminaltextcolor(WHITE);
	zlog << "done!" << std::endl;

	// add final grain for continuous layer
	grain_coord_array.push_back(std::vector <double>());
	grain_coord_array[grain_coord_array.size()-1].push_back(0.0); // x
	grain_coord_array[grain_coord_array.size()-1].push_back(0.0); // y

	// check for continuous layer
	for(unsigned int atom=0; atom < catom_array.size(); atom++){
		if(mp::material[catom_array[atom].material].continuous==true && catom_array[atom].include == false ){
			catom_array[atom].include=true;
			catom_array[atom].grain=int(grain_coord_array.size()-1);
		}
	}

	// set number of grains
	grains::num_grains = int(grain_coord_array.size());

	// sort atoms by grain number
	create::internal::sort_atoms_by_grain(catom_array);

	/*
	// Used to print out grain vertices for checking
    std::ofstream VERT_GNU_FILE("gnuplot_vert_file_3.dat");
    for(unsigned int i=0;i<NumGrains;++i){
        for(size_t j=0;j<grain_vertices_array[i].size();++j){
            VERT_GNU_FILE << grain_vertices_array[i][j][0] << " " << grain_vertices_array[i][j][1] << " "
						  << grain_coord_array[i][0] << " " << grain_coord_array[i][1] << "\n";
        }
        VERT_GNU_FILE << grain_vertices_array[i][0][0] << " " << grain_vertices_array[i][0][1] << " "
						  << grain_coord_array[i][0] << " " << grain_coord_array[i][1] << "\n";
        VERT_GNU_FILE << "\n\n";
    }
    VERT_GNU_FILE.flush();
    VERT_GNU_FILE.close();
	*/

	return EXIT_SUCCESS;
}

} // End of cs namespace
