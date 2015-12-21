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
#include <boost/program_options.hpp> // options parsing library
namespace fs = boost::filesystem;
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
  double hyperfine_cutoff;
  double hyperfine_cutoff_in_kHz;
  double c13_abundance;
  string lattice_file = "lattice-r[cell_radius]-s[seed].txt";

  bool pair_search;
  bool perform_scan;
  bool compute_fidelities;

  uint max_cluster_size;
  int ms;
  double static_B;
  double static_B_norm_in_gauss;
  uint k_DD;
  double scale_factor;

  uint scan_bins;
  double f_DD;
  double scan_time;
  double scan_time_in_ms;

  fs::path output_dir;
  int seed;
  bool no_output;

  fs::path lattice_path;
  fs::path larmor_path;
  fs::path scan_path;
  fs::path fidelities_path;
  string output_suffix = "r[cell_radius]-s[seed]-c[cluster_size]-k[k_DD]-[ms].txt";
  string output_suffix_with_input_lattice = "c[cluster_size]-k[k_DD]-[ms].txt";

  // define input options
  po::options_description options("Allowed options", 95);
  options.add_options()
    ("help,h", "produce help message")

    ("pair_search", po::value<bool>(&pair_search)->default_value(false)->implicit_value(true),
     "search for larmor pairs of nuclei")
    ("scan", po::value<bool>(&perform_scan)->default_value(false)->implicit_value(true),
     "perform coherence scan of effective larmor frequencies?")
    ("fidelities",
     po::value<bool>(&compute_fidelities)->default_value(false)->implicit_value(true),
     "compute expected iswap fidelities?")

    ("hyperfine_cutoff,c", po::value<double>(&hyperfine_cutoff_in_kHz)->default_value(10),
     "set cutoff scale for hyperfine field (in kHz)")
    ("c13_abundance", po::value<double>(&c13_abundance)->default_value(0.0107,"0.0107"),
     "relative isotopic abundance of C-13")
    ("lattice_file", po::value<string>(&lattice_file),
     "specify file defining system configuration")
    ("output_suffix", po::value<string>(&output_suffix), "output file suffix")

    ("max_cluster_size,m", po::value<uint>(&max_cluster_size)->default_value(6),
     "maximum allowable size of C-13 clusters")
    ("ms", po::value<int>(&ms)->default_value(1),
     "NV center spin state used with |0> for an effective two-level system (+/-1)")
    ("static_B,B", po::value<double>(&static_B_norm_in_gauss)->default_value(140.1,"140.1"),
     "strength of static magnetic field along the NV axis (in gauss)")
    ("k_DD,k", po::value<uint>(&k_DD)->default_value(1),
     "resonance harmonic used in spin addressing (1 or 3)")
    ("scale_factor", po::value<double>(&scale_factor)->default_value(100,"100"),
     "factor used to define different scales")

    ("scan_bins", po::value<uint>(&scan_bins)->default_value(100),
     "number of bins in coherence scanning range")
    ("f_DD", po::value<double>(&f_DD)->default_value(0.06,"0.06"),
     "magnitude of fourier component used in coherence scanning")
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
  bool set_hyperfine_cutoff = !inputs["hyperfine_cutoff"].defaulted();
  bool set_c13_abundance = !inputs["c13_abundance"].defaulted();
  bool using_input_lattice = inputs.count("lattice_file");
  bool set_output_suffix = inputs.count("output_suffix");

  // run a sanity check on inputs and parameter values
  assert(hyperfine_cutoff_in_kHz > 0);
  assert(c13_abundance >= 0 && c13_abundance <= 1);
  assert(!(using_input_lattice && set_hyperfine_cutoff));
  assert(!(using_input_lattice && set_c13_abundance));

  assert(max_cluster_size > 0);
  assert(ms == 1 || ms == -1);
  assert(static_B_norm_in_gauss >= 0);
  assert((k_DD == 1) || (k_DD == 3));
  assert(scale_factor > 1);

  if(perform_scan){
    assert(scan_bins > 0);
    assert(scan_time_in_ms > 0);
  }

  assert(seed > 0); // seeds of 0 and 1 give the same result, so don't allow nonpositive seeds

  // set some variables based on iputs
  hyperfine_cutoff = hyperfine_cutoff_in_kHz*kHz;
  static_B = static_B_norm_in_gauss*gauss;
  scan_time = scan_time_in_ms*1e-3;

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
    cell_radius = pow(abs(ge*gC13)/(4*pi*a0*a0*a0*hyperfine_cutoff),1.0/3);
    cout << "Setting cell radius to: " << cell_radius << endl;

    boost::replace_all(lattice_file, "[cell_radius]", to_string(cell_radius));
    boost::replace_all(lattice_file, "[seed]", to_string(seed));
    lattice_path = output_dir/fs::path(lattice_file);
  }

  srand(seed); // initialize random number generator
  fs::create_directory(output_dir); // create data directory

  // -----------------------------------------------------------------------------------------
  // Construct spin lattice
  // -----------------------------------------------------------------------------------------

  // vector of C-13 nuclei
  vector<spin> nuclei;

  if(!using_input_lattice){ // place nuclei at lattice cites

    // set positions of nuclei at lattice sites
    for(int b = 0; b <= 1; b++){
      for(int l = -2*cell_radius; l <= 2*cell_radius; l++){
        for(int m = -2*cell_radius; m <= 2*cell_radius; m++){
          for(int n = -2*cell_radius; n <= 2*cell_radius; n++){
            if(rnd() <= c13_abundance){ // if we pass a check for C-13 isotopic abundance
              nuclei.push_back(spin(b*ao+l*a1+m*a2+n*a3, gC13, s_vec/2));
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

    const uint lattice_sites = 2 * pow(2*cell_radius+1 ,3) - 2;
    cout << "C-13 nuclei occupy " << nuclei.size() << " of "
         << lattice_sites << " lattice sites" << endl << endl;

    // write cell radius and nucleus positions to file
    if(!no_output && !pair_search){
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
      nuclei.push_back(spin((Vector3d() << x,y,z).finished(),gC13,s_vec/2));
    }
    lattice.close();

    // assert that no C-13 nuclei lie at the NV lattice sites
    for(uint i = 0; i < nuclei.size(); i++){
      if((nuclei.at(i).pos == n.pos) || (nuclei.at(i).pos == e(ms).pos)){
        cout << "we have a C-13 nucleus at one of the NV lattice sites!" << endl;
        return 2;
      }
    }
  }

  // -----------------------------------------------------------------------------------------
  // Perform search for larmor pairs
  // -----------------------------------------------------------------------------------------

  if(pair_search){

    const double tolerance = 1e-5; // to account for numerical error
    for(uint i = 0; i < nuclei.size(); i++){

      const Vector3d A_i = A(nuclei.at(i));
      if((A_i).norm() < hyperfine_cutoff) continue;

      const double A_i_z = dot(A_i,zhat);
      const double A_i_xy = (A_i - dot(A_i,zhat)*zhat).norm();

      for(uint j = i+1; j < nuclei.size(); j++){

        const Vector3d A_j = A(nuclei.at(j));
        if((A_j).norm() < hyperfine_cutoff) continue;

        const double A_j_z = dot(A_j,zhat);
        const double A_j_xy = (A_j - dot(A_j,zhat)*zhat).norm();

        if(abs(A_i_z/A_j_z-1) < tolerance && abs(A_i_xy/A_j_xy-1) < tolerance){
          cout << "found larmor pair: " << i << ", " << j << endl;
          return 1;
        }
      }
    }

    cout << "no larmor pairs found" << endl;
    return 0;
  }

  // -----------------------------------------------------------------------------------------
  // Cluster C-13 nuclei
  // -----------------------------------------------------------------------------------------

  const double cluster_coupling_guess = 100; // this value doesn't really matter
  const double dcc_cutoff = 1e-5; // cutoff for tuning of cluster_coupling

  // get cluster_coupling for which the largest cluster size is >= max_cluster_size
  double cluster_coupling = find_target_coupling(nuclei, max_cluster_size,
                                                 cluster_coupling_guess, dcc_cutoff);
  vector<vector<uint>> ind_clusters = get_index_clusters(nuclei, cluster_coupling);
  // if (largest cluster size > max_cluster_size),
  //   which can occur when (largest cluster size == max_cluster_size) is impossible,
  //   find largest cluster_coupling for which (largest cluster size < max_cluster_size)
  while(largest_cluster_size(ind_clusters) > max_cluster_size){
    int cluster_size_target
      = largest_cluster_size(get_index_clusters(nuclei, cluster_coupling+dcc_cutoff));
    cluster_coupling = find_target_coupling(nuclei, cluster_size_target,
                                            cluster_coupling, dcc_cutoff);
    ind_clusters = get_index_clusters(nuclei, cluster_coupling);
  }

  if(max_cluster_size > 1){
    cout << "Nuclei grouped into " << ind_clusters.size() << " clusters"
         << " with a coupling factor of "  << cluster_coupling << " Hz" << endl;

    // collect and print histogram of cluster sizes
    max_cluster_size = largest_cluster_size(ind_clusters);
    vector<uint> size_hist(max_cluster_size);
    for(uint i = 0; i < ind_clusters.size(); i++){
      size_hist.at(ind_clusters.at(i).size()-1) += 1;
    }
    cout << "Cluster size histogram: " << endl;
    for(uint i = 0; i < size_hist.size(); i++){
      cout << "  " << i+1 << ": " << size_hist.at(i) << endl;
    }
    cout << endl;
  } else {
    cout << "Largest internuclear coupling: " << cluster_coupling << " Hz" << endl;
  }

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

    // identify effictive larmor frequencies and NV coupling strengths
    vector<double> w_larmor(nuclei.size());
    vector<double> A_perp(nuclei.size());

    double w_max = 0, w_min = DBL_MAX; // maximum and minimum effective larmor frequencies
    for(uint i = 0; i < nuclei.size(); i++){
      const Vector3d A_i = A(nuclei.at(i));
      A_perp.at(i) = (A_i-dot(A_i,zhat)*zhat).norm();
      w_larmor.at(i) = effective_larmor(nuclei.at(i), static_B, ms).norm();

      if(w_larmor.at(i) < w_min) w_min = w_larmor.at(i);
      if(w_larmor.at(i) > w_max) w_max = w_larmor.at(i);
    }

    // print effective larmor frequencies and NV couping strengths to output file
    if(!no_output){
      ofstream larmor(larmor_path.string());
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

    const double w_range = w_max - w_min;
    const double w_start = max(w_min - w_range/10, 0.);
    const double w_end = w_max + w_range/10;
    for(uint i = 0; i < scan_bins; i++){
      w_scan.at(i) = w_start + i*(w_end-w_start)/scan_bins;
      coherence.at(i) = coherence_measurement(nuclei, ind_clusters, scan_time,
                                              w_scan.at(i), k_DD, f_DD, static_B, ms);
      cout << "(" << i+1 << "/" << scan_bins << ") "
           << w_scan.at(i)/(2*pi*1e3) << " " << coherence.at(i) << endl;
    }

    // print coherence scan results to output file
    if(!no_output){
      ofstream scan(scan_path.string());
      scan << "# cluster coupling factor: " << cluster_coupling << endl;
      scan << "# f_DD: " << f_DD << endl;
      scan << "# scan time: " << scan_time << endl;
      cout << endl;
      scan << "# w_scan coherence\n";
      for(uint i = 0; i < scan_bins; i++){
        scan << w_scan.at(i) << " " << coherence.at(i) << endl;
      }
      scan.close();
    }
  }

  // -----------------------------------------------------------------------------------------
  // NV/nucleus SWAP fidelity
  // -----------------------------------------------------------------------------------------

  if(compute_fidelities){
    fidelities_path = output_dir/fs::path("fidelities-"+output_suffix);

    vector<uint> addressable_targets;
    vector<double> larmor_effs;
    vector<double> hyperfines;
    vector<double> hyperfine_perps;
    vector<double> dw_mins;
    vector<double> f_DDs;
    vector<double> operation_times;
    vector<double> fidelities;

    for(uint target = 0; target < nuclei.size(); target++){

      // larmor frequency of and hyperfine field at target nucleus
      const Vector3d larmor_eff = effective_larmor(nuclei.at(target), static_B, ms);
      const Vector3d hyperfine = A(nuclei.at(target));
      const Vector3d hyperfine_w = dot(hyperfine,hat(larmor_eff))*hat(larmor_eff);
      const Vector3d hyperfine_perp = hyperfine - hyperfine_w;

      // if our interaction strength is too weak, we can't address this nucleus
      if(hyperfine_perp.norm() < cluster_coupling) continue;

      // minimum difference in larmor frequencies between target nucleus and other nuclei
      double dw_min = DBL_MAX;
      for(uint s = 0; s < nuclei.size(); s++){
        if(s == target) continue;
        const double dw = abs(larmor_eff.norm()
                              - effective_larmor(nuclei.at(s), static_B, ms).norm());
        if(dw < dw_min) dw_min = dw;
      }

      // if this larmor frequency is too close to another, we cannot (yet) address the nucleus
      if(dw_min < cluster_coupling/scale_factor) continue;

      // AXY sequence parameters
      const double f_DD = -ms*dw_min/(hyperfine_perp.norm()*scale_factor);
      const double w_DD = larmor_eff.norm()/k_DD; // AXY protocol angular frequency
      const double t_DD = 2*pi/w_DD; // AXY protocol period
      const double operation_time = 2*pi/abs(f_DD*hyperfine_perp.norm());

      // iSWAP gate fidelity
      const double fidelity = iswap_fidelity(target, nuclei, ind_clusters, static_B, ms,
                                             k_DD, cluster_coupling, scale_factor);

      addressable_targets.push_back(target);
      larmor_effs.push_back(larmor_eff.norm());
      hyperfines.push_back(hyperfine.norm());
      hyperfine_perps.push_back(hyperfine_perp.norm());
      dw_mins.push_back(dw_min);
      f_DDs.push_back(f_DD);
      operation_times.push_back(operation_time);
      fidelities.push_back(fidelity);
    }

    if(!no_output){
      ofstream info_file(fidelities_path.string());
      info_file << "# index fidelity\n";
      info_file << "# cluster coupling factor: " << cluster_coupling << endl;
      info_file << "# scale factor: " << scale_factor << endl;
      cout << endl;
      for(uint t = 0; t < addressable_targets.size(); t++){
        info_file << addressable_targets.at(t) << " "
                  << larmor_effs.at(t) << " "
                  << hyperfines.at(t) << " "
                  << hyperfine_perps.at(t) << " "
                  << dw_mins.at(t) << " "
                  << f_DDs.at(t) << " "
                  << operation_times.at(t) << " "
                  << fidelities.at(t) << endl;
      }
      info_file.close();
    }
  }

}

