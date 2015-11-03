#include <sys/stat.h>
#include <float.h>

#include <iostream> // for standard output
#include <fstream> // for file input
#include <iomanip> // some nice printing functions
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
  double c13_abundance;
  int max_cluster_size;
  int ms;
  double Bz_in_gauss;

  string input_lattice;
  string output_lattice = "lattice-cr[cell_radius]-s[seed].txt";
  string output_dir;
  int seed;

  // define options
  po::options_description options("allowed options", 90);
  options.add_options()
    ("help,h", "produce help message")
    ("cell_radius", po::value<int>(&cell_radius)->default_value(6),
     "number of unit cells to simulate out from the NV center (along cell axes)")
    ("c13_abundance", po::value<double>(&c13_abundance)->default_value(0.0107,"0.0107"),
     "relative isotopic abundance of C-13")
    ("max_cluster_size", po::value<int>(&max_cluster_size)->default_value(6),
     "maximum allowable size of C-13 clusters")
    ("ms", po::value<int>(&ms)->default_value(1),
     "NV center spin state used with |0> for an effective two-level system (must be +/-1)")
    ("Bz", po::value<double>(&Bz_in_gauss)->default_value(140.1,"140.1"),
     "strength of static magnetic field along the NV axis (in gauss)")
    ("input_lattice", po::value<string>(&input_lattice),
     "input file defining system configuration")
    ("output_lattice",
     po::value<string>(&output_lattice),
     "output file defining system configuration")
    ("output_dir", po::value<string>(&output_dir)->default_value("./data"),
     "directory for storing data")
    ("seed", po::value<int>(&seed)->default_value(1),
     "seed for random number generator (must be >=1)")
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
  bool set_c13_abundance = !inputs["c13_abundance"].defaulted();
  bool using_input_lattice = inputs.count("input_lattice");
  bool set_output_lattice = inputs.count("output_lattice");
  bool set_output_dir = inputs.count("output_dir");

  // run a sanity check on inputs
  assert(cell_radius > 0);
  assert(c13_abundance >= 0 && c13_abundance <= 1);
  assert(ms == 1 || ms == -1);
  assert(Bz_in_gauss >= 0);
  assert(!(using_input_lattice && set_cell_radius));
  assert(!(using_input_lattice && set_c13_abundance));
  assert(!(using_input_lattice && set_output_lattice));
  assert(seed > 0);

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

  // set variables based on iputs
  double Bz = Bz_in_gauss*gauss;

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
              nuclei.push_back(spin(cell_pos+cell_sites.at(ls),gC13,s_vec));
            }
          }
        }
      }
    }

    // remove any C-13 atoms at the NV lattice sites
    nuclei.erase(remove(nuclei.begin(), nuclei.end(), spin(n.pos,gC13,s_vec)),
                 nuclei.end());
    nuclei.erase(remove(nuclei.begin(), nuclei.end(), spin(e(ms).pos,gC13,s_vec)),
                 nuclei.end());

    // write cell radius and C-13 positions to file
    ofstream output(output_lattice_path);
    output << "cell_radius: " << cell_radius << endl;
    for(int n = 0; n < nuclei.size(); n++){
      output << nuclei.at(n).pos(0) << " "
             << nuclei.at(n).pos(1) << " "
             << nuclei.at(n).pos(2) << endl;
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
      nuclei.push_back(spin((Vector3d() << x,y,z).finished(),gC13,s_vec));
    }
    input.close();

  }

  // assert that no C-13 spins lie at the NV lattice sites
  if(in_vector(spin(n.pos,gC13,s_vec),nuclei)
     || in_vector(spin(e(ms).pos,gC13,s_vec),nuclei)){
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
  double scan_time = 1e-3;

  double max_Ax = 0;
  double max_w = 0, min_w = DBL_MAX;
  for(int n = 0; n < nuclei.size(); n++){
    spin s = nuclei.at(n);
    double w = (s.g*Bz*zhat - ms/2.*A(s,ms)).norm();

    if(w < min_w) min_w = w;
    if(w > max_w) max_w = w;
  }

  int scans = 100;

  for(int i = 0; i < scans; i++){
    double w_scan = min_w + i*(max_w-min_w)/scans;
    cout << w_scan << "   "
         << coherence_scan(clusters, w_scan, k_DD, f_DD, Bz, ms, scan_time) << endl;
  }

  // double coherence = compute_coherence(clusters, w_scan, k_DD, f_DD, Bz, ms, scan_time);
  // cout << coherence << endl;

}

