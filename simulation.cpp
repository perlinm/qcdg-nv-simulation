#include <sys/stat.h>

#include <iostream>
#include <fstream>
using namespace std;

#include <eigen3/Eigen/Dense> // linear algebra library
using namespace Eigen;

#include <boost/algorithm/string.hpp> // string manipulation library
#include <boost/program_options.hpp> // options parsing library
namespace po = boost::program_options;

#include "constants.h"
#include "nv-math.h"
#include "gates.h"

int main(int arg_num, const char *arg_vec[]) {

  // -----------------------------------------------------------------------------------------
  // Parse and process input options
  // -----------------------------------------------------------------------------------------

  // define parameters
  int cell_radius;
  double c13_abundance = c13_natural_abundance;
  int max_cluster_size;
  int ms;
  string input_lattice;
  string output_lattice = "lattice-cr[cell_radius]-s[seed].txt";
  string output_dir = "data";
  bool fast_decomposition;
  int seed;


  double w_DD = 1e5;
  double Bz = 1000;

  // define options
  po::options_description options("allowed options");
  options.add_options()
    ("help,h", "produce help message")
    ("cell_radius", po::value<int>(&cell_radius)->default_value(6),
     "number of unit cells to simulate out from the NV center")
    ("c13_abundance", po::value<double>(&c13_abundance),
     "relative isotopic abundance of C-13")
    ("max_cluster_size", po::value<int>(&max_cluster_size)->default_value(6),
     "maximum allowable size of C-13 clusters")
    ("ms", po::value<int>(&ms)->default_value(1),
     "NV center spin used for effective two-level system (must be +/-1)")
    ("input_lattice", po::value<string>(&input_lattice), "input file definint initial system")
    ("output_lattice", po::value<string>(&output_lattice),
     "name for output file defining initial system")
    ("output_dir", po::value<string>(&output_dir), "directory for storing data")
    ("fast_decomposition,f", po::value<bool>(&fast_decomposition)->default_value(true),
     "use fast, but less accurate algorithm for matrix decomposition")
    ("seed,s", po::value<int>(&seed)->default_value(1), "seed for random number generator")


    (",w", po::value<double>(&w_DD), "set w_DD (in kHz)")
    (",B", po::value<double>(&Bz), "set Bz (in gauss)")


    ;

  // collect inputs
  po::variables_map inputs;
  po::store(parse_command_line(arg_num, arg_vec, options), inputs);
  po::notify(inputs);

  // if requested, print help text
  if(inputs.count("help")){
    cout << options;
    return 0;
  }

  // print summary of inputs
  if(arg_num > 1){
    cout << "------------------------------------------------------------------" << endl;
    cout << "running " << arg_vec[0] << " with parameters" << endl;
    for(int i = 1; i < arg_num; i++) {
      if(arg_vec[i][0] == '-') cout << endl;
      cout << arg_vec[i] << " ";
    }
    cout << endl;
    cout << "------------------------------------------------------------------" << endl;
    cout << endl;
  }

  bool set_cell_radius = !inputs["cell_radius"].defaulted();
  bool set_c13_abundance = inputs.count("c13_abundance");
  bool using_input_lattice = inputs.count("input_lattice");
  bool set_output_lattice = inputs.count("output_lattice");
  bool set_output_dir = inputs.count("output_dir");

  // run a sanity check on inputs
  assert(cell_radius > 0);
  assert(c13_abundance >= 0 && c13_abundance <= 1);
  assert(ms == 1 || ms == -1);
  assert(seed > 0);
  assert(!(using_input_lattice && set_cell_radius));
  assert(!(using_input_lattice && set_c13_abundance));
  assert(!(using_input_lattice && set_output_lattice));

  if(using_input_lattice){
    // check that the given input file exists
    if(access(input_lattice.c_str(),F_OK) == -1){
      cout << "file does not exist: " << input_lattice << endl;
      return 1;
    }
  }

  if(!using_input_lattice && !set_output_lattice){
    boost::replace_all(output_lattice, "[cell_radius]", to_string(cell_radius));
    boost::replace_all(output_lattice, "[seed]", to_string(seed));
  }

  string output_lattice_path = output_dir+"/"+output_lattice;

  if(!set_output_lattice && !set_output_dir){
    cout << "using default output lattice path: " << output_lattice_path << endl;
  } else {
    cout << "using output lattice path: " << output_lattice_path << endl;
  }

  cout << "total number of lattice sites: "
       << int(pow(2*cell_radius,3)*cell_sites.size()) << endl;

  srand(seed); // initialize random number generator
  mkdir(output_dir.c_str(),0777); // create data directory

  // -----------------------------------------------------------------------------------------
  // Construct spin lattice
  // -----------------------------------------------------------------------------------------

  // vector of C-13 nuclei
  vector<spin> nuclei;

  if(!using_input_lattice){

    // set positions of C-13 nuclei at lattice sites
    Vector3d cell_pos; // position of unit cell indexed by i,j,k
    for(int i = -cell_radius; i < cell_radius; i++){ // loop over all unit cell indices
      cell_pos(0) = i;
      for(int j = -cell_radius; j < cell_radius; j++){
        cell_pos(1) = j;
        for(int k = -cell_radius; k < cell_radius; k++){
          cell_pos(2) = k;
          for(int ls = 0; ls < cell_sites.size(); ls++){ // loop over sites in unit cell
            if(rnd() <= c13_abundance){ // if we pass a check for C-13 isotopic abundance
              nuclei.push_back(spin(cell_pos+cell_sites.at(ls),gC13));
            }
          }
        }
      }
    }

    // remove any C-13 atoms at the NV lattice sites
    nuclei.erase(remove(nuclei.begin(), nuclei.end(), spin(n.pos,gC13)), nuclei.end());
    nuclei.erase(remove(nuclei.begin(), nuclei.end(), spin(e.pos,gC13)), nuclei.end());

    // write cell radius and C-13 positions to file
    ofstream output(output_lattice_path);
    output << "cell_radius: " << cell_radius << endl;
    for(int c = 0; c < nuclei.size(); c++){
      output << nuclei.at(c).pos[0] << " "
             << nuclei.at(c).pos[1] << " "
             << nuclei.at(c).pos[2] << endl;
    }
    output.close();

  } else { // if using_input_lattice

    string line;
    ifstream input(input_lattice);

    // get cell_radius
    getline(input,line,' ');
    getline(input,line);
    cell_radius = stoi(line);

    // get C-13 positions
    double x,y,z;
    while(getline(input,line,' ')){
      x = stod(line);
      getline(input,line,' ');
      y = stod(line);
      getline(input,line);
      z = stod(line);
      nuclei.push_back(spin((Vector3d() << x,y,z).finished(),gC13));
    }
    input.close();

  }

  // assert that no C-13 spins lie at the NV lattice sites
  if(in_vector(spin(n.pos,gC13),nuclei) || in_vector(spin(e.pos,gC13),nuclei)){
    cout << "we have a C-13 nucleus at one of the NV lattice sites!" << endl;
    return 2;
  }

  // -----------------------------------------------------------------------------------------
  // Cluster C-13 nuclei
  // -----------------------------------------------------------------------------------------

  double cluster_coupling_guess = 100;
  double dcc_cutoff = 1e-5;

  // get cluster_coupling for which the largest cluster size is >= max_cluster_size
  double cluster_coupling = find_target_coupling(nuclei,max_cluster_size,
                                                 cluster_coupling_guess,dcc_cutoff);
  vector<vector<spin>> clusters = get_clusters(nuclei,cluster_coupling);

  // if largest cluster size > max_cluster_size,
  //   which can occur because it is impossible to exactly match max_cluster_size,
  //   find largest cluster_coupling for which the largest cluster size is < max_cluster_size
  if(largest_cluster_size(clusters) > max_cluster_size){
    int cluster_size_target
      = largest_cluster_size(get_clusters(nuclei,cluster_coupling+dcc_cutoff));
    cluster_coupling = find_target_coupling(nuclei,cluster_size_target,
                                            cluster_coupling,dcc_cutoff);
    clusters = get_clusters(nuclei,cluster_coupling);
  }

  cout << endl;
  cout << nuclei.size() << " nuclei grouped into " << clusters.size() << " clusters" << endl;
  cout << "largest cluster size: " << largest_cluster_size(clusters) << endl;
  cout << "cluster coupling factor: " << cluster_coupling << " Hz" << endl;
  cout << endl;

  // -----------------------------------------------------------------------------------------
  // Compute coherence
  // -----------------------------------------------------------------------------------------

  double f_DD = 0.06;
  harmonic k_DD = first;

  double coherence = compute_coherence(clusters, w_DD, k_DD, f_DD, Bz, ms);

  cout << coherence << endl;

}

