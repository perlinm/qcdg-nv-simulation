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
#include "qp-math.h"
#include "nv-math.h"
#include "printing.h"
#include "gates.h"

// random double from 0 to 1
inline double rnd(){ return std::rand()/(double)RAND_MAX; }

int main(int arg_num, const char *arg_vec[]) {

  // -----------------------------------------------------------------------------------------
  // Parse and process input options
  // -----------------------------------------------------------------------------------------

  // define parameters
  int cell_radius;
  double c13_abundance;
  string lattice_file = "lattice-cr[cell_radius]-s[seed].txt";

  uint max_cluster_size;
  int ms;
  double Bz_in_gauss;

  bool perform_scan;
  uint scan_bins;
  double f_DD;
  uint k_DD_int;
  double scan_time_in_ms;
  string larmor_file = "larmor-cr[cell_radius]-s[seed]-[ms].txt";
  string scan_file = "scan-cr[cell_radius]-s[seed]-[ms].txt";

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
    ("lattice_file", po::value<string>(&lattice_file),
     "specify file defining system configuration")

    ("max_cluster_size", po::value<uint>(&max_cluster_size)->default_value(6),
     "maximum allowable size of C-13 clusters")
    ("ms", po::value<int>(&ms)->default_value(1),
     "NV center spin state used with |0> for an effective two-level system (must be +/-1)")
    ("Bz", po::value<double>(&Bz_in_gauss)->default_value(140.1,"140.1"),
     "strength of static magnetic field along the NV axis (in gauss)")

    ("scan", po::value<bool>(&perform_scan)->default_value(true),
     "perform coherence scan of effective larmor frequencies?")
    ("scan_bins", po::value<uint>(&scan_bins)->default_value(200),
     "number of bins in coherence scanning range")
    ("f_DD", po::value<double>(&f_DD)->default_value(0.06,"0.06"),
     "magnitude of fourier component used in coherence scanning")
    ("k_DD", po::value<uint>(&k_DD_int)->default_value(1),
     "resonance harmonic used in coherence scanning (must be 1 or 3)")
    ("scan_time", po::value<double>(&scan_time_in_ms)->default_value(1,"1"),
     "time for each coherence measurement (in milliseconds)")

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
    printf("------------------------------------------------------------------\n");
    printf("running %s with parameters:\n", arg_vec[0]);
    for(int i = 1; i < arg_num; i++) {
      if(arg_vec[i][0] == '-') printf("\n");
      printf("%s ", arg_vec[i]);
    }
    printf("\n------------------------------------------------------------------\n\n");
  }

  bool set_cell_radius = !inputs["cell_radius"].defaulted();
  bool set_c13_abundance = !inputs["c13_abundance"].defaulted();
  bool input_lattice = inputs.count("lattice_file");

  // run a sanity check on inputs
  assert(cell_radius > 0);
  assert(c13_abundance >= 0 && c13_abundance <= 1);

  assert(ms == 1 || ms == -1);
  assert(Bz_in_gauss >= 0);

  if(perform_scan){
    assert(scan_bins > 0);
    assert((f_DD > 0) && (f_DD < 1));
    assert((k_DD_int == 1) || (k_DD_int == 3));
    assert(scan_time_in_ms > 0);
  }

  assert(!(input_lattice && set_cell_radius));
  assert(!(input_lattice && set_c13_abundance));

  assert(seed > 0);

  string lattice_path;
  if(input_lattice){
    lattice_path = lattice_file;
    // check that the given input file exists
    if(access(lattice_path.c_str(),F_OK) == -1){
      printf("file does not exist: %s\n",lattice_file.c_str());
      return 1;
    }
  } else{ // if !input_lattice
    boost::replace_all(lattice_file, "[cell_radius]", to_string(cell_radius));
    boost::replace_all(lattice_file, "[seed]", to_string(seed));
    lattice_path = output_dir+"/"+lattice_file;
  }

  printf("total number of lattice sites: %d\n", int(pow(2*cell_radius,3)*cell_sites.size()));

  boost::replace_all(larmor_file, "[cell_radius]", to_string(cell_radius));
  boost::replace_all(larmor_file, "[seed]", to_string(seed));
  boost::replace_all(larmor_file, "[ms]", (ms > 0)?"up":"dn");
  string larmor_path = output_dir+"/"+larmor_file;

  boost::replace_all(scan_file, "[cell_radius]", to_string(cell_radius));
  boost::replace_all(scan_file, "[seed]", to_string(seed));
  boost::replace_all(scan_file, "[ms]", (ms > 0)?"up":"dn");
  string scan_path = output_dir+"/"+scan_file;


  // set variables based on iputs
  double Bz = Bz_in_gauss*gauss;
  harmonic k_DD;
  if(k_DD_int == 1) k_DD = first;
  if(k_DD_int == 3) k_DD = third;
  double scan_time = scan_time_in_ms*1e-3;

  srand(seed); // initialize random number generator
  mkdir(output_dir.c_str(),0777); // create data directory

  // -----------------------------------------------------------------------------------------
  // Construct spin lattice
  // -----------------------------------------------------------------------------------------

  // vector of C-13 nuclei
  vector<spin> nuclei;

  if(!input_lattice){

    // set positions of C-13 nuclei at lattice sites
    Vector3d cell_pos; // position of unit cell indexed by i,j,k
    for(int i = -cell_radius; i < cell_radius; i++){ // loop over all unit cell indices
      cell_pos(0) = i;
      for(int j = -cell_radius; j < cell_radius; j++){
        cell_pos(1) = j;
        for(int k = -cell_radius; k < cell_radius; k++){
          cell_pos(2) = k;
          for(uint ls = 0; ls < cell_sites.size(); ls++){ // loop over sites in unit cell
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
    ofstream lattice(lattice_path);
    lattice << "# cell_radius: " << cell_radius << endl;
    for(uint n = 0; n < nuclei.size(); n++){
      lattice << nuclei.at(n).pos(0) << ' '
             << nuclei.at(n).pos(1) << ' '
             << nuclei.at(n).pos(2) << endl;
    }
    lattice.close();

  } else { // if input_lattice

    string line;
    ifstream lattice(lattice_path);

    // get cell_radius
    getline(lattice,line,' ');
    getline(lattice,line,' ');
    getline(lattice,line);
    cell_radius = stoi(line);

    // get C-13 positions
    double x,y,z;
    while(getline(lattice,line,' ')){
      x = stod(line);
      getline(lattice,line,' ');
      y = stod(line);
      getline(lattice,line);
      z = stod(line);
      nuclei.push_back(spin((Vector3d() << x,y,z).finished(),gC13,s_vec));
    }
    lattice.close();

  }

  // assert that no C-13 spins lie at the NV lattice sites
  if(in_vector(spin(n.pos,gC13,s_vec),nuclei)
     || in_vector(spin(e(ms).pos,gC13,s_vec),nuclei)){
    printf("we have a C-13 nucleus at one of the NV lattice sites!\n");
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

  printf("\n");
  printf("%d nuclei grouped into %d clusters\n",nuclei.size(),clusters.size());
  printf("largest cluster size: %d\n",largest_cluster_size(clusters));
  printf("cluster coupling factor: %g Hz\n\n",cluster_coupling);
  printf("\n");

  // -----------------------------------------------------------------------------------------
  // Coherence scan
  // -----------------------------------------------------------------------------------------

  if(perform_scan){

    // idenfy effictive larmor frequencies and NV coupling strengths
    VectorXd w_larmor = VectorXd::Zero(nuclei.size());
    VectorXd A_perp = VectorXd::Zero(nuclei.size());

    double w_max = 0, w_min = DBL_MAX;
    for(uint n = 0; n < nuclei.size(); n++){
      Vector3d A_n = A(nuclei.at(n),ms);
      w_larmor(n) = effective_larmor(nuclei.at(n),Bz*zhat,A_n,ms);
      A_perp(n) = (A_n-dot(A_n,zhat)*zhat).norm();

      if(w_larmor(n) < w_min) w_min = w_larmor(n);
      if(w_larmor(n) > w_max) w_max = w_larmor(n);
    }

    ofstream larmor(larmor_path);
    larmor << "# f_DD: " << f_DD << endl;
    larmor << "# k_DD: " << k_DD_int << endl;
    larmor << "# w A_perp" << endl;
    for(uint n = 0; n < nuclei.size(); n++){
      larmor << w_larmor(n) << " " << A_perp(n) << endl;
    }
    larmor.close();

    // perform coherence scan
    printf("beginning coherence scan\n");
    VectorXd w_scan = VectorXd::Zero(scan_bins);
    VectorXd coherence = VectorXd::Zero(scan_bins);

    double w_range = w_max - w_min;
    double w_start = w_min - w_range/10;
    double w_end = w_max + w_range/10;
    for(uint i = 0; i < scan_bins; i++){
      w_scan(i) = w_start + i*(w_end-w_start)/scan_bins;
      coherence(i) =
        coherence_measurement(ms, clusters, w_scan(i), k_DD, f_DD, Bz, scan_time);
      printf("(%d/%d) %g %g\n",i+1,scan_bins,w_scan(i),coherence(i));
    }

    ofstream scan(scan_path);
    scan << "# f_DD: " << f_DD << endl;
    scan << "# k_DD: " << k_DD_int << endl;
    scan << "# w L" << endl;
    for(uint i = 0; i < scan_bins; i++){
      scan << w_scan(i) << " " << coherence(i) << endl;
    }
    scan.close();
  }

}

