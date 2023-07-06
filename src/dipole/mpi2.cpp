//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins and Richard F L Evans 2020. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <cmath>
#include <cstdlib>
#include <sstream>

// Vampire headers
#include "vio.hpp"
#include "vmpi.hpp"
#include "vutil.hpp"

// cells module headers
#include "internal.hpp"

namespace dipole{
namespace internal{

      //-------------------------------------------------------------------------
      // Function to send and receive data between cells on different processors
      //-------------------------------------------------------------------------
      void initialise_atomistic_cell_data(
         const int num_cells,
         const int num_local_cells,
         const double cutoff,                                         // cutoff range for dipole tensor construction (Angstroms)
         const std::vector<int>& num_atoms_in_cell,                   // number of atoms in each cell (local CPU)
         const std::vector<int>& list_of_local_cells,                 // numerical list of cells containing atoms on local processor
         const std::vector<int>& global_atoms_in_cell_count,          // number of atoms in each cell (all CPUs)
         const std::vector<double>& pos_and_mom_array,                // array of positions and cell moments
         const std::vector < std::vector <int> >& index_atoms_array,  // 2D array of [cells][atomID]
         const std::vector<double>& atoms_coords_x,                   // input arrays of atom coordinates
         const std::vector<double>& atoms_coords_y,                   //
         const std::vector<double>& atoms_coords_z,                   //
         const std::vector<double>& atoms_moments,                    // input array of atom moments (Bohr magnetons)
         std::vector<int>& list_of_cells_with_atoms,                  // list of cells to access atoms
         std::vector< std::vector<double> >& atoms_in_cells_array     // output array of positions and moments of atoms in cells
      ){

         #ifdef MPICF
         // print informative message to user
         std::cout << "Transferring atomistic data for parallel dipole initialisation " << std::flush;
         zlog << zTs() << "Transferring atomistic data for parallel dipole initialisation" << std::endl;

         // wait for all processes to reach here
         vmpi::barrier();
         //std::cerr << "Rank " << vmpi::my_rank << " is here" << std::endl;

         // instantiate timer
         vutil::vtimer_t timer;

         // start timer
         timer.start();

         //---------------------------------------------------
         // Determine cells I require on local processor
         //---------------------------------------------------
         // bool array to check whether a cell is needed on local processor
         std::vector<bool> i_need_this_cell(num_cells, false);

         const double cutoffsq = cutoff * cutoff;

         // loop over local cells where the dipole field must be calculated
         // since I have atoms in this cell
         for(int lc = 0; lc < num_local_cells; lc++){

            // extract local cell ID
            const int cell_id = list_of_local_cells[lc];

            // load coordinates of local cell
            const double cxi = cells_pos_and_mom_array[ 4*cell_id + 0 ];
            const double cyi = cells_pos_and_mom_array[ 4*cell_id + 1 ];
            const double czi = cells_pos_and_mom_array[ 4*cell_id + 2 ];

            // loop over all other cells in the system to determine if it is within
            // cutoff range for atomistic dipole tensor calculation
            for(int cell = 0; cell < num_cells; cell++){

               // check that the remote cell has any atoms
               if( global_atoms_in_cell_count[cell] > 0 ){

                  // load coordinates of local cell
                  const double cxj = cells_pos_and_mom_array[ 4*cell + 0 ];
                  const double cyj = cells_pos_and_mom_array[ 4*cell + 1 ];
                  const double czj = cells_pos_and_mom_array[ 4*cell + 2 ];

                  // compute displacement
                  const double dx = cxi - cxj;
                  const double dy = cyi - cyj;
                  const double dz = czi - czj;

                  // calculate distance between cells squared
                  const double rij2 = dx*dx + dy*dy + dz*dz;

                  // determine of cell is in range and tag as needed (previously included local cells explicitly but not needed since all cells are looped over)
                  if( rij2 <= cutoffsq ) i_need_this_cell[cell] = true;

               } // end of if statement checking num_atoms_in_cell
            } // end of loop over num_cells
         } // end of loop over local cells

         // now determine list of cells that I need (could combine this into one step)
         std::vector<int> list_of_cells_i_need(0);
         for(int cell = 0; cell < num_cells; cell++){
            if( i_need_this_cell[cell] ) list_of_cells_i_need.push_back( cell );
         }

         //vmpi::barrier();
         //std::cerr << "Rank " << vmpi::my_rank << " is starting to accumulate" << std::endl;

         //-----------------------------------------------------------
         // Accumulate atomic positions for all cells in linear order
         // (so each CPU only has complete cell lists) (distributed)
         //-----------------------------------------------------------

         // determine number of cells on local process
         int all_num_cells = num_cells / vmpi::num_processors;
         int my_num_cells = all_num_cells;
         // allocate remaining cells to last process
         if ( vmpi::my_rank == vmpi::num_processors - 1 ) my_num_cells = num_cells - ( my_num_cells * ( vmpi::num_processors - 1 ) );

         // calculate the first cell to process on my rank
         const int my_first_cell = vmpi::my_rank * all_num_cells;

         // determine parallel list of cell locations (for processing) reduced on all CPUs
         std::vector<int> cell_location_list(num_cells,0);
         for ( int cell = my_first_cell; cell < my_first_cell + my_num_cells; cell++ ) cell_location_list[ cell ] = vmpi::my_rank;
         vmpi::all_reduce_sum(cell_location_list);

         vmpi::barrier();
         //std::cerr << "Rank " << vmpi::my_rank << " first cell " <<  my_first_cell << " my last cell " << my_first_cell + my_num_cells << " total cells " << num_cells << std::endl;


         // determine number of atoms to accumulate
         int total_num_atoms = 0;
         for ( int cell = my_first_cell; cell < my_first_cell + my_num_cells; cell++ ){
            total_num_atoms += global_atoms_in_cell_count[cell];
            if(cell >= num_cells) std::cerr << "************************** " << cell << "\t" << num_cells << std::endl;
         }

         // initialise storage for distributed atom data
         std::vector<double> recv_atom_data( 4 * total_num_atoms );

         std::stringstream textssa;
         //textssa << "Rank: " << vmpi::my_rank << " cell atoms ";
         //for ( int cell = 0; cell < num_cells; cell++ ) textssa << cell << " " <<  global_atoms_in_cell_count[cell] << " ";
         //std::cerr << textssa.str() << std::endl;
         vmpi::barrier();
         //std::cerr << "Rank " << vmpi::my_rank << " total num atoms to recv " <<  total_num_atoms << std::endl;

         // initialise temporary data structures (this is memory efficient since data is only created when atoms need to be sent)
         std::vector<int> num_atoms_to_send( vmpi::num_processors, 0); // number of atoms to be sent to each process
         std::vector<int> num_cells_to_send( vmpi::num_processors, 0); // number of cells to be sent to each process
         std::vector< std::vector <int> > cell_ids_to_send( vmpi::num_processors ); // 2D list of cellID and atom count data [ cellID ][num_atoms_in_cell_from_proc]

         int atoms_send_buffer_size = 0;
         int cells_send_buffer_size = 0;

         // loop over local cells to determine number of atoms to send to each processor
         for(int lc = 0; lc < num_local_cells; lc++){

            // extract local cell ID
            const int cell_id = list_of_local_cells[lc];

            // get location of cell
            const int location = cell_location_list[cell_id];

            // increment counters with number of local atoms in the cell and number of cells
            num_atoms_to_send[location] += num_atoms_in_cell[cell_id];
            atoms_send_buffer_size += num_atoms_in_cell[cell_id];

            // save cell IDs and number of atoms to send for non-zero values only
            if( num_atoms_in_cell[cell_id] > 0){
               num_cells_to_send[location]++;
               cell_ids_to_send[location].push_back( cell_id ); // first value is cell location of atoms
               cell_ids_to_send[location].push_back( num_atoms_in_cell[cell_id] ); // second value is number of atoms to send
               cells_send_buffer_size++;
            }

         }

         //std::cerr << "Rank " << vmpi::my_rank << " Packing data " << atoms_send_buffer_size << " atoms " << std::endl;

         //std::stringstream textssa;
         //textssa << "Rank: " << vmpi::my_rank << " cell atom counts:\n";

         // pack data for sending
         std::vector<double> atom_send_buffer( 4 * atoms_send_buffer_size );
         std::vector<int> cell_send_buffer( 2 * cells_send_buffer_size,-1 );
         int counter = 0; // counter for adding data
         int cells_counter = 0; // counter for cells data
         for ( int p = 0; p < vmpi::num_processors; p++){
            //std::cerr << "Rank " << vmpi::my_rank << " cells to send " << cell_ids_to_send[p].size() << "to proc " << p << "\n";
            // loop over all cells to send to processor p
            for(int id = 0; id < cell_ids_to_send[p].size(); id+=2){

               const int cell = cell_ids_to_send[p][id]; // get cell ID

               // pack cell buffer for sending
               cell_send_buffer[2*cells_counter + 0] = cell_ids_to_send[p][id];   // cell ID where atoms are located
               cell_send_buffer[2*cells_counter + 1] = cell_ids_to_send[p][id+1]; // number of atoms in cell on local CPU

               // increment cells counter
               cells_counter++;

               //textssa << cell << " " << num_atoms_in_cell[cell] << "\n";
               // loop over all atoms in cell
               for( int atomid = 0; atomid < num_atoms_in_cell[cell]; atomid++ ){
                  const int atom = index_atoms_array[cell][atomid];
                  atom_send_buffer[4*counter + 0] = atoms_coords_x[atom];
                  atom_send_buffer[4*counter + 1] = atoms_coords_y[atom];
                  atom_send_buffer[4*counter + 2] = atoms_coords_z[atom];
                  atom_send_buffer[4*counter + 3] = atoms_moments[atom];
                  counter++;
               }
            }
         }
         //std::cerr << textssa.str() << std::endl;
         if(counter != atoms_send_buffer_size) std::cerr << "Error! mismatching buffer size "  << atoms_send_buffer_size << " and data count " << counter << " on rank " << vmpi::my_rank << std::endl;

         //std::stringstream textssd;
   //      for(int i=0; i< cell_send_buffer.size()/2; i++){
   //         textssd << "Rank " << vmpi::my_rank << " cells buffer: " << i << " " << cell_send_buffer[2*i+0] << " " << cell_send_buffer[2*i+1] << "\n";
   //      }
         //for(int i=0; i< cell_send_buffer.size(); i++){
         //   textssd << "Rank " << vmpi::my_rank << " cells buffer: " << i << " " << cell_send_buffer[i] << "\n";
         //}
         //std::cerr << textssd.str() << std::endl;

         // determine matrix for send counts
         std::vector<int> num_atoms_to_recv( vmpi::num_processors );
         std::vector<int> num_cells_to_recv( vmpi::num_processors );

         // exchange send and receive counts
         #ifdef MPICF
            MPI_Alltoall(&num_atoms_to_send[0], 1, MPI_INT, &num_atoms_to_recv[0], 1, MPI_INT, MPI_COMM_WORLD);
            MPI_Alltoall(&num_cells_to_send[0], 1, MPI_INT, &num_cells_to_recv[0], 1, MPI_INT, MPI_COMM_WORLD);
         #endif

         /*std::stringstream textss;
         for ( int p = 0; p < vmpi::num_processors; p++){
            textss << "Rank " << vmpi::my_rank << "[" << p << "]: send cells count " << num_cells_to_send[p] << " send atoms count " << num_atoms_to_send[p] << "\n";
            textss << "Rank " << vmpi::my_rank << "[" << p << "]: recv cells count " << num_cells_to_recv[p] << " recv atoms count " << num_atoms_to_recv[p] << "\n";
         }

         std::cerr << textss.str() << std::endl;*/

         // store compressed data offsets and counts for send and recieve operations
         std::vector<int> send_cell_offsets(0);
         std::vector<int> recv_cell_offsets(0);
         std::vector<int> send_cell_counts(0);
         std::vector<int> recv_cell_counts(0);

         int send_cell_offset_counter = 0;
         //int send_cells_data_counter   = 0;
         int recv_cell_offset_counter = 0;
         //int recv_cells_data_counter   = 0;

         std::vector<int> send_atom_offsets(0);
         std::vector<int> recv_atom_offsets(0);
         std::vector<int> send_atom_counts(0);
         std::vector<int> recv_atom_counts(0);

         int send_atom_offset_counter  = 0;
         //int send_atom_data_counter    = 0;
         int recv_atom_offset_counter  = 0;
         //int recv_atom_data_counter    = 0;

         // counter for total num cells to recieve (greater than actual number of
         // cells as atoms from partial cells are distributed)
         int total_num_cells = 0;

         // determine where to place data from each process (in message order)
         for ( int p = 0; p < vmpi::num_processors; p++){
            if( num_atoms_to_send[p] > 0){
               send_atom_offsets.push_back( send_atom_offset_counter );
               send_atom_counts.push_back(  4*num_atoms_to_send[p]   );
               send_atom_offset_counter  += 4*num_atoms_to_send[p];
               send_cell_offsets.push_back( send_cell_offset_counter );
               send_cell_counts.push_back(  2*num_cells_to_send[p]   );
               send_cell_offset_counter  += 2*num_cells_to_send[p];
            }
            if( num_atoms_to_recv[p] > 0 ){
               recv_atom_offsets.push_back( recv_atom_offset_counter );
               recv_atom_counts.push_back(  4*num_atoms_to_recv[p]   );
               recv_atom_offset_counter  += 4*num_atoms_to_recv[p];
               recv_cell_offsets.push_back( recv_cell_offset_counter );
               recv_cell_counts.push_back(  2*num_cells_to_recv[p]   );
               recv_cell_offset_counter  += 2*num_cells_to_recv[p];
               total_num_cells += num_cells_to_recv[p];
            }
         }

         //for(int i=0; i<recv_cell_offsets.size(); i++ ){
         //   std::cerr << "Rank " << vmpi::my_rank << " SC " << i << " " << send_cell_offsets[i] << " " << send_cell_counts[i] << "\n";
         //   std::cerr << "Rank " << vmpi::my_rank << " RC " << i << " " << recv_cell_offsets[i] << " " << recv_cell_counts[i] << "\n";
         //}

         std::vector<int> recv_cell_data( 2*total_num_cells,123 );
         //std::cerr << "Rank " << vmpi::my_rank << " num cells to recv " << total_num_cells << " offset counter " << recv_cell_offset_counter << "\n";

         vmpi::barrier();
         //std::cerr << "Rank " << vmpi::my_rank << " swapping atom data " << std::endl;

         // send receive all data asynchronously
         {

            // array of MPI requests
            std::vector<MPI_Request> requests(0);

            MPI_Request req = MPI_REQUEST_NULL; // temporary variable for push_back operations
            //MPI_Status status; // temporary variable for stati

            //std::stringstream textsse;
            //for(int i=0; i< cell_send_buffer.size()/2; i++){
            //   textsse << "Rank " << vmpi::my_rank << " cells buffer: " << i << " " << cell_send_buffer[2*i+0] << " " << cell_send_buffer[2*i+1] << "\n";
            //}
            //std::cerr << textsse.str() << std::endl;

            int recv_message_ID = 0;
            // loop over all processors and dispatch only necessary sends and recieves
            for ( int cpu = 0; cpu < vmpi::num_processors; cpu++){
               //int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request)
               // check that atoms are to be sent or receieved, and post receives for atom and cell data
               if ( num_atoms_to_recv[cpu] > 0 ){
                  int offset = recv_atom_offsets[recv_message_ID];
                  int count  = recv_atom_counts [recv_message_ID];
                  //std::cerr << "Rank " << vmpi::my_rank << " posted recv from rank " << cpu << " with count " << count/4 << " and offset " << offset/4 << "\n";
                  requests.push_back(req); // add storage for request handle
                  MPI_Irecv(&recv_atom_data[offset], count, MPI_DOUBLE, cpu, 654, MPI_COMM_WORLD, &requests.back());
                  // recieve for cells data and first message
                  int cell_offset = recv_cell_offsets[recv_message_ID];
                  int cell_count  = recv_cell_counts [recv_message_ID];
                  //std::cerr << "Rank " << vmpi::my_rank << " posted recv from rank " << cpu << " with count " << cell_count << " and offset " << cell_offset << "\n";
                  requests.push_back(req); // add storage for request handle
                  MPI_Irecv(&recv_cell_data[cell_offset], cell_count, MPI_INT, cpu, 634, MPI_COMM_WORLD, &requests.back());

                  // increment message ID counter
                  recv_message_ID++;
               }
            }

            MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
            //int error;
            int send_message_ID = 0;
            // loop over all processors and dispatch only necessary sends and recieves
            for ( int cpu = 0; cpu < vmpi::num_processors; cpu++){
               vmpi::barrier();
               //int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request)
               if ( num_atoms_to_send[cpu] > 0 ){
                  int offset = send_atom_offsets[send_message_ID];
                  int count  = send_atom_counts [send_message_ID];
                  //std::cerr << "Rank " << vmpi::my_rank << " posted send to rank " << cpu << " with count " << count/4 << " and offset " << offset/4 << "\n";
                  requests.push_back(req); // add storage for request handle
                  MPI_Isend(&atom_send_buffer[offset], count, MPI_DOUBLE, cpu, 654, MPI_COMM_WORLD, &requests.back());
                  int cell_offset = send_cell_offsets[send_message_ID];
                  int cell_count  = send_cell_counts [send_message_ID];
                  //std::cerr << "Rank " << vmpi::my_rank << " posted send to rank " << cpu << " with count " << cell_count << " and offset " << cell_offset << " num cells to send " << num_cells_to_send[cpu] << "\n";
                  //for(int i=cell_offset; i < cell_offset+cell_count; i++) std::cerr << "  -> " << i << " " << cell_send_buffer[i] << "\n";
                  requests.push_back(req); // add storage for request handle
                  MPI_Isend(&cell_send_buffer[cell_offset], cell_count, MPI_INT, cpu, 634, MPI_COMM_WORLD, &requests.back());
                  //int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
                  // if (error != MPI_SUCCESS) {
                  //    char error_string[256];
                  //    int length_of_error_string, error_class;
                  //    MPI_Error_class(error, &error_class);
                  //    MPI_Error_string(error_class, error_string, &length_of_error_string);
                  //    fprintf(stderr, "%3d: %s\n", vmpi::my_rank, error_string);
                  //    MPI_Error_string(error, error_string, &length_of_error_string);
                  //    fprintf(stderr, "%3d: %s\n", vmpi::my_rank, error_string);
                  //    //send_error = TRUE;
                  // }
                  // increment message ID counter
                  send_message_ID++;
               }
            }

            //vmpi::barrier();

            // Wait for all communications to complete
            int num_requests = requests.size();
            std::vector<MPI_Status> stati(num_requests);
            MPI_Waitall(num_requests, &requests[0], &stati[0]);
            //std::stringstream textssf;
            /*for(int i=0; i< stati.size(); i++){
               vmpi::barrier();
               int num_messages=3456;
               MPI_Get_count(&stati[i], MPI_INT, &num_messages);
               int num_d_messages=3456;
               MPI_Get_count(&stati[i], MPI_DOUBLE, &num_d_messages);
               //textssf << "Rank " << vmpi::my_rank << " num_messages " << i << " " << num_messages << " " << num_d_messages << "\n";
            }
            //std::cerr << textssf.str() << std::endl;*/
         }

         /*std::vector<int> fbuff;
         if(vmpi::my_rank==0){
            int cell_offset = send_cell_offsets[1];
            int cell_count  = send_cell_counts [1];
            std::cerr << "Rank " << vmpi::my_rank << " posted send to rank " << 1 << " with count " << cell_count << " and offset " << cell_offset << " buffer size " << cell_send_buffer.size() << "\n";
            //MPI_Send(&cell_send_buffer[cell_offset], cell_count, MPI_INT, 1, 222, MPI_COMM_WORLD);
            MPI_Send(&cell_send_buffer[0], cell_send_buffer.size(), MPI_INT, 1, 222, MPI_COMM_WORLD);
            vmpi::barrier();
            //std::vector<int> buff(1000);
            int rcell_offset = recv_cell_offsets[1];
            int rcell_count  = recv_cell_counts [1];
            MPI_Status req;
            std::cerr << "Rank " << vmpi::my_rank << " posted recv from rank " << 1 << " with count " << rcell_count << " and offset " << rcell_offset << "\n";
            MPI_Status status;
            MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            // Allocate memory to receive data
            int count;
            MPI_Get_count(&status, MPI_INT, &count);
            std::cerr << "Message contains " << count << " values " << std::endl;
            std::vector<int> buff(count);
            MPI_Recv(&buff[0], count, MPI_INT, 1, 222, MPI_COMM_WORLD, &req);
            fbuff.resize(buff.size());
            for(int i=0; i< fbuff.size(); i++) fbuff[i] = buff[i];
                        //MPI_Recv(&recv_cell_data[rcell_offset], rcell_count, MPI_INT, 1, 223, MPI_COMM_WORLD, &req);
         }
         if(vmpi::my_rank==1){
            //std::vector<int> buff(1000);
            int rcell_offset = recv_cell_offsets[0];
            int rcell_count  = recv_cell_counts [0];
            MPI_Status req;
            std::cerr << "Rank " << vmpi::my_rank << " posted recv from rank " << 0 << " with count " << rcell_count << " and offset " << rcell_offset << "\n";
            MPI_Status status;
            MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            // Allocate memory to receive data
            int count;
            MPI_Get_count(&status, MPI_INT, &count);
            std::cerr << "Message contains " << count << " values " << std::endl;
            std::vector<int> buff(count);
            //MPI_Recv(&recv_cell_data[rcell_offset], rcell_count, MPI_INT, 0, 222, MPI_COMM_WORLD, &req);
            MPI_Recv(&buff[0], count, MPI_INT, 0, 222, MPI_COMM_WORLD, &req);
            vmpi::barrier();
            int cell_offset = send_cell_offsets[0];
            int cell_count  = send_cell_counts [0];
            std::cerr << "Rank " << vmpi::my_rank << " posted send to rank " << 0 << " with count " << cell_count << " and offset " << cell_offset << "\n";
            //MPI_Send(&cell_send_buffer[cell_offset], cell_count, MPI_INT, 0, 223, MPI_COMM_WORLD);
            MPI_Send(&cell_send_buffer[0], cell_send_buffer.size(), MPI_INT, 0, 222, MPI_COMM_WORLD);
            fbuff.resize(buff.size());
            for(int i=0; i< fbuff.size(); i++) fbuff[i] = buff[i];

         }*/

         /*std::stringstream textssc;
         for ( int i = 0; i < recv_cell_data.size()/2 ; i++){
            textssc << vmpi::my_rank << " " << recv_cell_data[2*i+0] << " " << recv_cell_data[2*i+1] << "\n";
         }
         std::cerr << textssc.str() << std::endl;

         std::stringstream textssb;
         for ( int i = 0; i < recv_atom_data.size()/4 ; i++){
            textssb << vmpi::my_rank << " " << recv_atom_data[4*i+0] << " " << recv_atom_data[4*i+1] << " " << recv_atom_data[4*i+2] << " " << recv_atom_data[4*i+3] << "\n";
         }
         //std::cerr << textssb.str() << std::endl;

         std::vector<double> all_atoms;
         vmpi::collate(recv_atom_data, all_atoms);

         if(vmpi::my_rank == 0){

            std::ofstream ofile("data.txt");
            for ( int i = 0; i < all_atoms.size()/4 ; i++){
               ofile << all_atoms[4*i+0] << " " << all_atoms[4*i+1] << " " << all_atoms[4*i+2] << " " << all_atoms[4*i+3] << "\n";
            }
         }*/

         /*std::stringstream textssg;
         int count=0;
         for ( int i = 0; i < recv_cell_data.size()/2 ; i++){
            const int cellID = recv_cell_data[2*i+0];
            const int num_atoms = recv_cell_data[2*i+1];
            const int offset = count;
            textssg << vmpi::my_rank << " " << cellID << " " << num_atoms << " ";
            for (int atom = offset; atom < offset+num_atoms; atom++){
               textssg << " [ " << recv_atom_data[4*atom+0] << " " << recv_atom_data[4*atom+1] << " " << recv_atom_data[4*atom+2] << " " << recv_atom_data[4*atom+3] << "] ";
            }
            textssg << "\n";
            count+=num_atoms;
         }
         std::cerr << textssg.str() << std::endl;*/

         //---------------------------------------------------
         // now organise data linear in cells
         //---------------------------------------------------
         // [ x y z m x y z m ][ x y z m x y z m]
         //      cell 3             cell 4
         //---------------------------------------------------


         std::vector< std::vector<double> > ordered_cells_atoms_list_2D( my_num_cells ); // array to store sorted list of atom indexes
         //std::vector<int> ordered_atoms_count( my_num_cells ); // array to store sorted list of cell counts

         // loop over all cells and place into temporary 2D storage
         int count=0;
         for ( int i = 0; i < recv_cell_data.size()/2 ; i++){
            const int cellID = recv_cell_data[2*i+0]; // get global cell ID
            const int cell_index = cellID - my_first_cell; // convert to cell index on local processor
            const int num_atoms = recv_cell_data[2*i+1]; // get number of atoms associated with this cell
            const int offset = count; // set the offset in the atoms array
            // now loop over all atoms, adding list of indices
            for (int atom = offset; atom < offset+num_atoms; atom++){
               ordered_cells_atoms_list_2D[cell_index].push_back( recv_atom_data[4*atom+0] );
               ordered_cells_atoms_list_2D[cell_index].push_back( recv_atom_data[4*atom+1] );
               ordered_cells_atoms_list_2D[cell_index].push_back( recv_atom_data[4*atom+2] );
               ordered_cells_atoms_list_2D[cell_index].push_back( recv_atom_data[4*atom+3] );
            }
            count+=num_atoms;
         }

         //actually 2D is more useful at this stage...

         //std::vector<int> ordered_atoms_cellID( 4*my_num_cells ); // array to store sorted list of cellIDs

         // now squash 2D data into 1D for sending
         //std::vector<double> ordered_cells_atoms_list_1D( recv_atom_data.size() );
         //int index=0;
                  std::stringstream textssh;
         for(int cell=0; cell<my_num_cells; cell++){



            // check total atom counts for each cell agree with global number
            if( ordered_cells_atoms_list_2D[cell].size()/4 != global_atoms_in_cell_count[my_first_cell + cell]){
               std::cerr << "Programmer error on rank " << vmpi::my_rank << " in accumulation of atomic posisitions in cells. Expecting ";
               std::cerr << global_atoms_in_cell_count[my_first_cell + cell] << " but only have " << ordered_cells_atoms_list_2D[cell].size() << std::endl;
            }

            /*textssh << vmpi::my_rank << " " << cell+my_first_cell << " " << num_atoms << " ";
            // if we have the right number of atoms, go ahead and store these in cell order for sending to all other CPUs
            for (int atom = 0; atom < ordered_cells_atoms_list_2D[cell].size()/4; atom++){
               //const int atom_ID = ordered_cells_atoms_list_2D[cell][atom];
               //ordered_cells_atoms_list_1D[4*index+0] = recv_atom_data[4*atom_ID+0];
               //ordered_cells_atoms_list_1D[4*index+1] = recv_atom_data[4*atom_ID+1];
               //ordered_cells_atoms_list_1D[4*index+2] = recv_atom_data[4*atom_ID+2];
               //ordered_cells_atoms_list_1D[4*index+3] = recv_atom_data[4*atom_ID+3];
               //index++;
               textssh << " [ " << ordered_cells_atoms_list_2D[cell][4*atom+0] << " " << ordered_cells_atoms_list_2D[cell][4*atom+1] << " " << ordered_cells_atoms_list_2D[cell][4*atom+2] << " " << ordered_cells_atoms_list_2D[cell][4*atom+3] << "] ";
            }
            textssh << "\n";*/
         }
         //std::cerr << textssh.str() << std::endl;

         // output cell data to screen for posterity

         /*count=0;
         for(int cell=0; cell<my_num_cells; cell++){
            const int num_atoms = global_atoms_in_cell_count[my_first_cell + cell]; // get number of atoms associated with this cell
            const int offset = count;
            textssh << vmpi::my_rank << " " << cell+my_first_cell << " " << num_atoms << " ";
            for (int atom = offset; atom < num_atoms+offset; atom++){
               textssh << " [ " << ordered_cells_atoms_list_1D[4*atom+0] << " " << ordered_cells_atoms_list_1D[4*atom+1] << " " << ordered_cells_atoms_list_1D[4*atom+2] << " " << ordered_cells_atoms_list_1D[4*atom+3] << "] ";
            }
            textssh << "\n";
            count+=num_atoms;
         }
         std::cerr << textssh.str() << std::endl;*/

            // determine destination of cell atoms and number of atoms to send num_data_to_send[proc]

         //all_to_all cell_start, cell_count

         //int MPI_Alltoallv(const void *sendbuf, const int *sendcounts,
         //                  const int *sdispls, MPI_Datatype sendtype, void *recvbuf,
         //                  const int *recvcounts, const int *rdispls, MPI_Datatype recvtype,
         //                  MPI_Comm comm)

         //---------------------------------------------------------------------------
         // determine number of cells and atoms i need from each process
         //---------------------------------------------------------------------------

         // reset send/recv counts reusing previous memory
         std::fill (num_atoms_to_send.begin(),num_atoms_to_send.end(), 0);
         std::fill (num_cells_to_send.begin(),num_cells_to_send.end(), 0);
         std::fill (num_atoms_to_recv.begin(),num_atoms_to_recv.end(), 0);
         std::fill (num_cells_to_recv.begin(),num_cells_to_recv.end(), 0);

         // 2D array for lists of cells from each processor
         std::vector< std::vector<int> > cells_i_need_from_cpu( vmpi::num_processors );


         std::stringstream textssk;
         //textssk << "Rank " << vmpi::my_rank << " " << list_of_cells_i_need.size() << " cells i need : ";
         for( int i = 0; i < list_of_cells_i_need.size() ; i++ ){

            // get global cell number
            const int cell = list_of_cells_i_need[i];

            //textssk << cell << " ";

            // get location (processor) of cell
            const int location = cell_location_list[cell];

            // increment cell and atom counters
            num_cells_to_recv[location] ++;
            num_atoms_to_recv[location] += global_atoms_in_cell_count[cell];

            // save cell in processor ordered list
            cells_i_need_from_cpu[location].push_back(cell);

         }
         //textssk << "\n";
         //std::cerr << textssk.str() << std::endl;

         // Swap lists of atom and processor counts using MPI_alltoall magic
         #ifdef MPICF
            MPI_Alltoall(&num_atoms_to_recv[0], 1, MPI_INT, &num_atoms_to_send[0], 1, MPI_INT, MPI_COMM_WORLD);
            MPI_Alltoall(&num_cells_to_recv[0], 1, MPI_INT, &num_cells_to_send[0], 1, MPI_INT, MPI_COMM_WORLD);
         #endif

         /*std::stringstream textssi;
         for ( int p = 0; p < vmpi::num_processors; p++){
            textssi << "Rank " << vmpi::my_rank << "[" << p << "]: send cells count " << num_cells_to_send[p] << " send atoms count " << num_atoms_to_send[p] << "\n";
            textssi << "Rank " << vmpi::my_rank << "[" << p << "]: recv cells count " << num_cells_to_recv[p] << " recv atoms count " << num_atoms_to_recv[p] << "\n";
         }

         std::cerr << textssi.str() << std::endl;*/

         //---------------------------------------------------------------------------
         // send the list of cells I need to each process (asynchronously P2P)
         // note that send and recieve data structures being backwards is actually
         // correct...
         //---------------------------------------------------------------------------

         // data for storing list of cells I need to send to other processors (2D)
         std::vector< std::vector<int> > list_of_cells_to_send_2D( vmpi::num_processors );
         for(int p = 0; p < vmpi::num_processors; p++ ) list_of_cells_to_send_2D[p].resize( num_cells_to_send[p] );

         // send receive all data
         {
            // array of MPI requests
            std::vector<MPI_Request> requests(0);

            MPI_Request req = MPI_REQUEST_NULL; // temporary variable for push_back operations

            int recv_message_ID = 0;
            // loop over all processors and dispatch only necessary sends and recieves
            for ( int cpu = 0; cpu < vmpi::num_processors; cpu++){
               // check that atoms are to be sent or receieved, and post receives for atom and cell data
               if ( num_cells_to_send[cpu] > 0 ){

                  // add storage for request handle
                  requests.push_back(req);

                  // recieve list of cells I need to send back to cpu
                  MPI_Irecv(&list_of_cells_to_send_2D[cpu][0], num_cells_to_send[cpu], MPI_INT, cpu, 650, MPI_COMM_WORLD, &requests.back());

                  // increment message ID counter
                  recv_message_ID++;

               }
            }

            int send_message_ID = 0;
            // loop over all processors and dispatch only necessary sends and recieves
            for ( int cpu = 0; cpu < vmpi::num_processors; cpu++){
               if ( num_cells_to_recv[cpu] > 0 ){

                  // add storage for request handle
                  requests.push_back(req);

                  // send list of cells I need from cpu
                  MPI_Isend(&cells_i_need_from_cpu[cpu][0], num_cells_to_recv[cpu], MPI_INT, cpu, 650, MPI_COMM_WORLD, &requests.back());

                  // increment message ID counter
                  send_message_ID++;

               }
            }

            // Wait for all communications to complete
            int num_requests = requests.size();
            MPI_Waitall(num_requests, &requests[0], MPI_STATUS_IGNORE);

         }

         // wait for everyone before transferring data
         vmpi::barrier();
         //std::cerr << vmpi::my_rank << " Transferring atom data " << std::endl;

         // print out list of cells I need to send here
         /*std::stringstream textssj;
         for ( int cpu = 0; cpu < vmpi::num_processors; cpu++){
            textssj << vmpi::my_rank << " " << list_of_cells_to_send_2D[cpu].size() << " : ";
            for(int i=0; i<list_of_cells_to_send_2D[cpu].size();i++ ) textssj << list_of_cells_to_send_2D[cpu][i] << " ";
            textssj << "\n";
         }
         std::cerr << textssj.str() << std::endl;*/


         //---------------------------------------------------------------------------
         // pack atoms into 2D buffers for sending
         //---------------------------------------------------------------------------

         // allocate storage for all atoms and offsets
         std::vector < std::vector <double> > atom_send_buffers_2D( vmpi::num_processors );

         //std::stringstream textssm;
         // loop over all processors and pack data into buffer
         for ( int cpu = 0; cpu < vmpi::num_processors; cpu++){
            //textssm << " Buffer for transfer from " << vmpi::my_rank << " to " << cpu << " ";

            const int num_cells = list_of_cells_to_send_2D[cpu].size(); // same as num_cells_to_send

            // check for on-zero arrays before memory operations
            if( num_cells > 0 ){

               // pre-reserve memory for efficiency
               atom_send_buffers_2D[cpu].reserve( 4*num_atoms_to_send[cpu] );

               // loop over all cells we need to send to cpu
               for(int id = 0; id < num_cells; id++ ){

                  // get actual cell_ID relative to start
                  const int cell = list_of_cells_to_send_2D[cpu][id] - my_first_cell;

                  // get number of atoms in cell to pack
                  const int num_atoms = ordered_cells_atoms_list_2D[cell].size()/4;

                  // loop over all atoms and pack
                  for (int atom = 0; atom < num_atoms; atom++){
                     atom_send_buffers_2D[cpu].push_back( ordered_cells_atoms_list_2D[cell][4*atom+0] );
                     atom_send_buffers_2D[cpu].push_back( ordered_cells_atoms_list_2D[cell][4*atom+1] );
                     atom_send_buffers_2D[cpu].push_back( ordered_cells_atoms_list_2D[cell][4*atom+2] );
                     atom_send_buffers_2D[cpu].push_back( ordered_cells_atoms_list_2D[cell][4*atom+3] );
                  }
               }

               //for(int i=0; i<atom_send_buffers_2D[cpu].size()/4; i++){
               //   textssm << "        -> " << atom_send_buffers_2D[cpu][4*i+0] << " " << atom_send_buffers_2D[cpu][4*i+1] << " " << atom_send_buffers_2D[cpu][4*i+2] << " " << atom_send_buffers_2D[cpu][4*i+3] << "\n";
               //}
            }
         }

         //std::cerr << textssm.str() << std::endl;
         //for(int cell=0; cell<my_num_cells; cell++){
         //textssh << vmpi::my_rank << " " << cell+my_first_cell << " " << num_atoms << " ";
         //for (int atom = 0; atom < ordered_cells_atoms_list_2D[cell].size()/4; atom++){

         // allocate storage to recieve [proc][atoms]
         std::vector< std::vector< double > > recv_buffer_2D( vmpi::num_processors );
         for ( int cpu = 0; cpu < vmpi::num_processors; cpu++) recv_buffer_2D[cpu].resize( 4*num_atoms_to_recv[cpu] );

         //---------------------------------------------------------------------------
         // send by p2p comms
         //---------------------------------------------------------------------------
         {
            // array of MPI requests
            std::vector<MPI_Request> requests(0);

            MPI_Request req = MPI_REQUEST_NULL; // temporary variable for push_back operations

            // loop over all processors and dispatch only necessary sends and recieves
            for ( int cpu = 0; cpu < vmpi::num_processors; cpu++){
               // check that atoms are to be sent or receieved, and post receives for atom and cell data
               if ( num_atoms_to_recv[cpu] > 0 ){
                  //std::cerr << "Recv on rank " << vmpi::my_rank << " from rank " << cpu << " with data size " << num_atoms_to_recv[cpu] << " " << recv_buffer_2D[cpu].size() << std::endl;

                  // add storage for request handle
                  requests.push_back(req);

                  // recieve list of cells I need to send back to cpu
                  MPI_Irecv(&recv_buffer_2D[cpu][0], 4*num_atoms_to_recv[cpu], MPI_DOUBLE, cpu, 651, MPI_COMM_WORLD, &requests.back());

               }
            }

            // loop over all processors and dispatch only necessary sends and recieves
            for ( int cpu = 0; cpu < vmpi::num_processors; cpu++){
               if ( num_atoms_to_send[cpu] > 0 ){
                  //std::cerr << "Send on rank " << vmpi::my_rank << " to rank " << cpu << " with data size " << num_atoms_to_send[cpu] << " " << atom_send_buffers_2D[cpu].size() << std::endl;

                  // add storage for request handle
                  requests.push_back(req);

                  // send list of cells I need from cpu
                  MPI_Isend(&atom_send_buffers_2D[cpu][0], 4*num_atoms_to_send[cpu], MPI_DOUBLE, cpu, 651, MPI_COMM_WORLD, &requests.back());

               }
            }

            // Wait for all communications to complete
            int num_requests = requests.size();
            MPI_Waitall(num_requests, &requests[0], MPI_STATUS_IGNORE);

         }

         // check data as received
         /*std::stringstream textssn;
         for ( int cpu = 0; cpu < vmpi::num_processors; cpu++){
            textssn << " Recv buffer for transfer to " << vmpi::my_rank << " from " << cpu << "\n";
            for( int i=0; i<recv_buffer_2D[cpu].size()/4; i++){
               textssn << "        -> " << recv_buffer_2D[cpu][4*i+0] << " " << recv_buffer_2D[cpu][4*i+1] << " " << recv_buffer_2D[cpu][4*i+2] << " " << recv_buffer_2D[cpu][4*i+3] << "\n";
            }
         }
         std::cerr << textssn.str() << std::endl;*/

         // finally unpack data to make it accessible in [cells][atoms] format
         //cells_i_need_from_cpu - list of cells and in the order that they are sent to me
         //std::vector < std::vector <double > > atoms_in_cells_array ( list_of_cells_i_need.size() );
         atoms_in_cells_array.resize( list_of_cells_i_need.size() );

         // check that cells from processors match up with expected order from my list of cells I need
         int cell_index = 0;

         // loop over all cells in processor order
         for ( int cpu = 0; cpu < vmpi::num_processors; cpu++){

            // buffer index counter to keep track of which atoms have been extracted (reset for each CPU)
            int buff_index = 0;

            for( int i = 0; i < cells_i_need_from_cpu[cpu].size(); i++ ){

               // get cell ID in order to find number of atoms to unpack
               const int cell = cells_i_need_from_cpu[cpu][i];

               // check that cell IDs match up for sanity
               if( cell != list_of_cells_i_need[cell_index]){
                  std::cerr << "Programmer error - mismatch in received cell ID " << cell << " from that expected from cell list " << list_of_cells_i_need[cell_index] << std::endl;
               }

               // now get number of atoms in that cell
               const int num_atoms_in_cell = global_atoms_in_cell_count[cell];

               // reserve memory for atoms
               atoms_in_cells_array[i].reserve( 4*num_atoms_in_cell );

               // now unpack atoms into cell storage
               for( int atom = 0 ; atom < num_atoms_in_cell ; atom++ ){
                  atoms_in_cells_array[cell_index].push_back( recv_buffer_2D[cpu][4*buff_index+0] );
                  atoms_in_cells_array[cell_index].push_back( recv_buffer_2D[cpu][4*buff_index+1] );
                  atoms_in_cells_array[cell_index].push_back( recv_buffer_2D[cpu][4*buff_index+2] );
                  atoms_in_cells_array[cell_index].push_back( recv_buffer_2D[cpu][4*buff_index+3] );
                  buff_index++; // increment buffer index to point to next atom in buffer
               }

               // increment cell index counter
               cell_index++;

            }
         }

         // output final data for checking
         /*std::stringstream textssl;
         for (int i=0; i< atoms_in_cells_array.size(); i++){
            const int cell = list_of_cells_i_need[i];

            textssl << vmpi::my_rank << " " << cell << " " << atoms_in_cells_array[i].size()/4 << " : ";
            for(int atom=0; atom< atoms_in_cells_array[i].size()/4; atom++ ){
               textssl << " [ " << atoms_in_cells_array[i][4*atom+0] << " " << atoms_in_cells_array[i][4*atom+1] << " " << atoms_in_cells_array[i][4*atom+2] << " " << atoms_in_cells_array[i][4*atom+3] << " ]";
            }
            textssl << "\n";
         }
         std::cerr << textssl.str() << std::endl;*/

         // save list of atoms with cells
         list_of_cells_with_atoms = list_of_cells_i_need;

         // hold parallel calculation until all processors have completed the initialisation
         vmpi::barrier();

         // stop timer
         timer.stop();

         std::cout << "done! [ " << timer.elapsed_time() << " s ]" << std::endl;
         zlog << zTs() << "Transfer of atomistic data for dipole calculation complete. Time taken: " << timer.elapsed_time() << " s"<< std::endl;

         #else

            // loop over all cell identifying all cells with atoms
            for(int cell=0; cell<num_cells; cell++){
               if(global_atoms_in_cell_count[cell] > 0){
                  // add cell with > 0 atoms to list
                  list_of_cells_with_atoms.push_back(cell);
                  atoms_in_cells_array.push_back(std::vector<double>());
                  const int last_cell_index = atoms_in_cells_array.size()-1;
                  atoms_in_cells_array[last_cell_index].reserve(4*global_atoms_in_cell_count[cell]);
                  // loop over all atoms in cell and add to condensed list
                  for(int atom_index = 0; atom_index < global_atoms_in_cell_count[cell]; atom_index++){
                     const int atom = index_atoms_array[cell][atom_index];
                     atoms_in_cells_array[last_cell_index].push_back(atoms_coords_x[atom]);
                     atoms_in_cells_array[last_cell_index].push_back(atoms_coords_y[atom]);
                     atoms_in_cells_array[last_cell_index].push_back(atoms_coords_z[atom]);
                     atoms_in_cells_array[last_cell_index].push_back(atoms_moments [atom]);
                  }
               }
            }

         #endif

         return;

      }

   } // end of namespace internal

} // end namespace dipole
