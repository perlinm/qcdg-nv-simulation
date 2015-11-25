#include <sys/stat.h>
#include <float.h>

#include <iostream> // for standard output
#include <fstream> // for file input
#include <iomanip> // some nice printing functions
using namespace std;

#include <eigen3/Eigen/Dense> // linear algebra library
using namespace Eigen;

#include <boost/algorithm/string.hpp> // string manipulation library
#include <boost/filesystem.hpp> // filesystem path manipulation library
namespace fs = boost::filesystem;
#include <boost/program_options.hpp> // options parsing library
namespace po = boost::program_options;

#include "constants.h"
#include "qp-math.h"
#include "nv-math.h"
#include "printing.h"
#include "gates.h"

// random double from 0 to 1
inline double rnd(){ return rand()/(double)RAND_MAX; }

int main(int arg_num, const char *arg_vec[]) {

  // -----------------------------------------------------------------------------------------
  // Parse and process input options
  // -----------------------------------------------------------------------------------------

  // define inputs and parameters
  int cell_radius;
  double c13_abundance;
  string lattice_file = "lattice-r[cell_radius]-s[seed].txt";

  uint max_cluster_size;
  int ms;
  Vector3d B_static;
  double B_static_norm_in_gauss;

  bool perform_scan;
  int scan_bins;
  double f_DD;
  uint k_DD;
  double scan_time;
  double scan_time_in_ms;

  fs::path output_dir;
  int seed;
  bool no_output;

  fs::path lattice_path;
  fs::path larmor_path;
  fs::path scan_path;
  string output_suffix = "r[cell_radius]-s[seed]-c[cluster_size]-k[k_DD]-[ms].txt";
  string output_suffix_with_input_lattice = "c[cluster_size]-k[k_DD]-[ms].txt";

  // define input options
  po::options_description options("Allowed options", 95);
  options.add_options()
    ("help,h", "produce help message")
    ("cell_radius,r", po::value<int>(&cell_radius)->default_value(7),
     "number of unit cells to simulate out from the NV center (along cell axes)")
    ("c13_abundance", po::value<double>(&c13_abundance)->default_value(0.0107,"0.0107"),
     "relative isotopic abundance of C-13")
    ("lattice_file", po::value<string>(&lattice_file),
     "specify file defining system configuration")
    ("output_suffix", po::value<string>(&output_suffix), "output file suffix")

    ("max_cluster_size,c", po::value<uint>(&max_cluster_size)->default_value(6),
     "maximum allowable size of C-13 clusters")
    ("ms,m", po::value<int>(&ms)->default_value(1),
     "NV center spin state used with |0> for an effective two-level system (+/-1)")
    ("B_static,B", po::value<double>(&B_static_norm_in_gauss)->default_value(140.1,"140.1"),
     "strength of static magnetic field along the NV axis (in gauss)")

    ("scan", po::value<bool>(&perform_scan)->default_value(false)->implicit_value(true),
     "perform coherence scan of effective larmor frequencies?")
    ("scan_bins", po::value<int>(&scan_bins)->default_value(100),
     "number of bins in coherence scanning range")
    ("f_DD", po::value<double>(&f_DD)->default_value(0.06,"0.06"),
     "magnitude of fourier component used in coherence scanning")
    ("k_DD,k", po::value<uint>(&k_DD)->default_value(1),
     "resonance harmonic used in coherence scanning (1 or 3)")
    ("scan_time", po::value<double>(&scan_time_in_ms)->default_value(1),
     "time for each coherence measurement (in microseconds)")

    ("output_dir", po::value<fs::path>(&output_dir)->default_value("./data"),
     "directory for storing data")
    ("seed,s", po::value<int>(&seed)->default_value(1),
     "seed for random number generator (>=1)")
    ("no_output", po::value<bool>(&no_output)->default_value(false)->implicit_value(true),
     "don't generate output files")
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

  // determine whether certain options were used
  bool set_cell_radius = !inputs["cell_radius"].defaulted();
  bool set_c13_abundance = !inputs["c13_abundance"].defaulted();
  bool using_input_lattice = inputs.count("lattice_file");
  bool set_output_suffix = inputs.count("output_suffix");

  // run a sanity check on inputs and parameter values
  assert(cell_radius > 0);
  assert(c13_abundance >= 0 && c13_abundance <= 1);
  assert(!(using_input_lattice && set_cell_radius));
  assert(!(using_input_lattice && set_c13_abundance));

  assert(max_cluster_size > 0);
  assert(ms == 1 || ms == -1);
  assert(B_static_norm_in_gauss >= 0);

  if(perform_scan){
    assert(scan_bins > 0);
    assert((f_DD > 0) && (f_DD < 1));
    assert((k_DD == 1) || (k_DD == 3));
    assert(scan_time_in_ms > 0);
  }

  assert(seed > 0);

  // define path of input or output file defining system configuration
  if(using_input_lattice){
    if(!fs::exists(lattice_file)){
      cout << "file does not exist: " << lattice_file << endl;
      return 1;
    }
    lattice_path = lattice_file;

    if(!set_output_suffix){
      output_suffix = lattice_path.stem().string() + "-" + output_suffix_with_input_lattice;
    }

  } else{ // if !using_input_lattice
    boost::replace_all(lattice_file, "[cell_radius]", to_string(cell_radius));
    boost::replace_all(lattice_file, "[seed]", to_string(seed));
    lattice_path = output_dir/fs::path(lattice_file);
  }

  // set some variables based on iputs
  B_static = B_static_norm_in_gauss*gauss*zhat;
  scan_time = scan_time_in_ms*1e-3;

  srand(seed); // initialize random number generator
  fs::create_directory(output_dir); // create data directory

  // -----------------------------------------------------------------------------------------
  // Construct spin lattice
  // -----------------------------------------------------------------------------------------

  cout << "Total number of lattice sites: " << int(pow(2*cell_radius,3)*cell_sites.size());
  cout << endl << endl;

  // vector of C-13 nuclei
  vector<spin> nuclei;

  if(!using_input_lattice){ // place nuclei at lattice cites

    // set positions of nuclei at lattice sites
    Vector3d cell_pos; // position of unit cell indexed by i,j,k
    for(int i = -cell_radius; i < cell_radius; i++){ // loop over all unit cell indices
      cell_pos(0) = i;
      for(int j = -cell_radius; j < cell_radius; j++){
        cell_pos(1) = j;
        for(int k = -cell_radius; k < cell_radius; k++){
          cell_pos(2) = k;
          for(uint ls = 0; ls < cell_sites.size(); ls++){ // loop over lattice sites in cell
            if(rnd() <= c13_abundance){ // if we pass a check for C-13 isotopic abundance
              nuclei.push_back(spin(cell_pos+cell_sites.at(ls),gC13,s_vec));
            }
          }
        }
      }
    }

    // remove any nuclei at the NV lattice sites
    for(uint i = 0; i < nuclei.size(); i++){
      if((nuclei.at(i).pos == n.pos) || (nuclei.at(i).pos == e(ms).pos)){
        nuclei.erase(nuclei.begin()+i);
      }
    }

    // write cell radius and nucleus positions to file
    if(!no_output){
      ofstream lattice(lattice_path.string());
      lattice << "# cell radius: " << cell_radius << endl;
      for(uint i = 0; i < nuclei.size(); i++){
        lattice << nuclei.at(i).pos(0) << ' '
                << nuclei.at(i).pos(1) << ' '
                << nuclei.at(i).pos(2) << endl;
      }
      lattice.close();
    }

  } else { // if using_input_lattice, read in the lattice

    string line;
    ifstream lattice(lattice_path.string());

    // get cell_radius
    getline(lattice,line,' ');
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

    // assert that no C-13 spins lie at the NV lattice sites
    for(uint i = 0; i < nuclei.size(); i++){
      if((nuclei.at(i).pos == n.pos) || (nuclei.at(i).pos == e(ms).pos)){
        cout << "we have a C-13 nucleus at one of the NV lattice sites!" << endl;
        return 2;
      }
    }
  }

  // -----------------------------------------------------------------------------------------
  // Cluster C-13 nuclei
  // -----------------------------------------------------------------------------------------

  cout << "Clustering " << nuclei.size() << " nuclei" << endl;
  double cluster_coupling_guess = 100; // this value doesn't really matter
  double dcc_cutoff = 1e-5; // cutoff for tuning of cluster_coupling in clustering algorithm

  // get cluster_coupling for which the largest cluster size is >= max_cluster_size
  double cluster_coupling = find_target_coupling(nuclei,max_cluster_size,
                                                 cluster_coupling_guess,dcc_cutoff);
  vector<vector<spin>> clusters = get_clusters(nuclei,cluster_coupling);

  // if (largest cluster size > max_cluster_size),
  //   which can occur when (largest cluster size == max_cluster_size) is impossible,
  //   find largest cluster_coupling for which (largest cluster size < max_cluster_size)
  while(largest_cluster_size(clusters) > max_cluster_size){
    int cluster_size_target
      = largest_cluster_size(get_clusters(nuclei,cluster_coupling+dcc_cutoff));
    cluster_coupling = find_target_coupling(nuclei,cluster_size_target,
                                            cluster_coupling,dcc_cutoff);
    clusters = get_clusters(nuclei,cluster_coupling);
  }
  cout << "Nuclei grouped into " << clusters.size() << " clusters"
       << " with a coupling factor of "  << cluster_coupling << " Hz" << endl;

  // collect and print histogram of cluster sizes
  max_cluster_size = largest_cluster_size(clusters);
  vector<uint> size_hist(max_cluster_size);
  for(uint i = 0; i < clusters.size(); i++){
    size_hist.at(clusters.at(i).size()-1) += 1;
  }
  cout << "Cluster size histogram: " << endl;
  for(uint i = 0; i < size_hist.size(); i++){
    cout << "  " << i+1 << ": " << size_hist.at(i) << endl;
  }
  cout << endl;

  // now that we are done with initialization, fix up the output filename suffix
  boost::replace_all(output_suffix, "[cell_radius]", to_string(cell_radius));
  boost::replace_all(output_suffix, "[seed]", to_string(seed));
  boost::replace_all(output_suffix, "[cluster_size]", to_string(max_cluster_size));
  boost::replace_all(output_suffix, "[k_DD]", to_string(k_DD));
  boost::replace_all(output_suffix, "[ms]", (ms > 0)?"up":"dn");

  // -----------------------------------------------------------------------------------------
  // Coherence scan
  // -----------------------------------------------------------------------------------------

  if(perform_scan){
    // define paths of output files
    larmor_path = output_dir/fs::path("larmor-"+output_suffix);
    scan_path = output_dir/fs::path("scan-"+output_suffix);

    // define header to output files
    stringstream file_header;
    file_header << "# cluster coupling factor (Hz): " << cluster_coupling << endl;
    file_header << "# f_DD: " << f_DD << endl;
    file_header << "# scan time: " << scan_time << endl << endl;

    // identify effictive larmor frequencies and NV coupling strengths
    vector<double> w_larmor(nuclei.size());
    vector<double> A_perp(nuclei.size());

    double w_max = 0, w_min = DBL_MAX; // maximum and minimum effective larmor frequencies
    for(uint i = 0; i < nuclei.size(); i++){
      Vector3d A_i = A(nuclei.at(i));
      w_larmor.at(i) = effective_larmor(nuclei.at(i), B_static, A_i, ms);
      A_perp.at(i) = (A_i-dot(A_i,zhat)*zhat).norm();

      if(w_larmor.at(i) < w_min) w_min = w_larmor.at(i);
      if(w_larmor.at(i) > w_max) w_max = w_larmor.at(i);
    }

    // print effective larmor frequencies and NV couping strengths to output file
    if(!no_output){
      ofstream larmor(larmor_path.string());
      larmor << file_header.str();
      larmor << "# w_larmor A_perp\n";
      for(uint i = 0; i < nuclei.size(); i++){
        larmor << w_larmor.at(i) << " " << A_perp.at(i) << endl;
      }
      larmor.close();
    }

    // perform coherence scan
    cout << "Beginning coherence scan" << endl;
    vector<double> w_scan(scan_bins);
    vector<double> coherence(scan_bins);
    control_fields controls;

    // double w = 130.712*2*pi*1e3;
    // cout << coherence_measurement(ms, clusters, w, k_DD, f_DD, scan_time, B_static) << endl;
    // cout << coherence_measurement(ms, clusters, w, k_DD, f_DD, scan_time,
                                  // B_static, controls) << endl;
    // return 5;


    double w_range = w_max - w_min;
    double w_start = max(w_min - w_range/10, 0.);
    double w_end = w_max + w_range/10;
    for(int i = 0; i < scan_bins; i++){
      w_scan.at(i) = w_start + i*(w_end-w_start)/scan_bins;
      coherence.at(i) = coherence_measurement(ms, clusters, w_scan.at(i), k_DD, f_DD,
                                              scan_time, B_static, controls);
      cout << "(" << i+1 << "/" << scan_bins << ") "
           << w_scan.at(i)/(2*pi*1e3) << " " << coherence.at(i) << endl;
      // if(i+1 >= 10) break;
    }

    // print coherence scan results to output file
    if(!no_output){
      ofstream scan(scan_path.string());
      scan << file_header.str();
      scan << "# w_scan coherence\n";
      for(int i = 0; i < scan_bins; i++){
        scan << w_scan.at(i) << " " << coherence.at(i) << endl;
      }
      scan.close();
    }
  }

}

