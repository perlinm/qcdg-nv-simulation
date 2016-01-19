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

  const uint help_text_length = 95;

  int seed;

  po::options_description general("General options", help_text_length);
  general.add_options()
    ("help,h", "produce help message")
    ("seed", po::value<int>(&seed)->default_value(1),
     "seed for random number generator (>=1)")
    ;

  string lattice_file = "lattice-r[cell_radius]-s[seed].txt";
  fs::path output_dir;
  string output_suffix = "r[cell_radius]-s[seed]-c[cluster_size]-k[k_DD]-[ms].txt";
  string output_suffix_with_input_lattice = "c[cluster_size]-k[k_DD]-[ms].txt";
  bool no_output;

  po::options_description file_io("File IO",help_text_length);
  file_io.add_options()
    ("lattice_file", po::value<string>(&lattice_file),
     "input file defining system configuration")
    ("output_dir", po::value<fs::path>(&output_dir)->default_value("./data"),
     "directory for storing data")
    ("output_suffix", po::value<string>(&output_suffix), "output file suffix")
    ("no_output", po::value<bool>(&no_output)->default_value(false)->implicit_value(true),
     "don't generate output files")
    ;

  bool pair_search;
  bool coherence_scan;
  bool compute_iswap_fidelities;
  bool testing;

  po::options_description simulations("Available simulations",help_text_length);
  simulations.add_options()
    ("pairs", po::value<bool>(&pair_search)->default_value(false)->implicit_value(true),
     "search for larmor pairs")
    ("scan", po::value<bool>(&coherence_scan)->default_value(false)->implicit_value(true),
     "perform coherence scan of effective larmor frequencies")
    ("iswap",
     po::value<bool>(&compute_iswap_fidelities)->default_value(false)->implicit_value(true),
     "compute expected iswap fidelities")
    ("test" ,po::value<bool>(&testing)->default_value(false)->implicit_value(true),
     "enable testing mode")
    ;

  double c13_abundance;
  uint max_cluster_size;
  double hyperfine_cutoff;
  double hyperfine_cutoff_in_kHz;
  int cell_radius; // determined by hyperfine_cutoff
  int ms;
  uint k_DD;
  double static_Bz_in_gauss;
  double scale_factor;

  po::options_description simulation_options("Simulation options",help_text_length);
  simulation_options.add_options()
    ("c13_abundance", po::value<double>(&c13_abundance)->default_value(0.0107,"0.0107"),
     "relative isotopic abundance of C-13")
    ("max_cluster_size", po::value<uint>(&max_cluster_size)->default_value(6),
     "maximum allowable size of C-13 clusters")
    ("hyperfine_cutoff", po::value<double>(&hyperfine_cutoff_in_kHz)->default_value(100),
     "set cutoff scale for hyperfine field (kHz)")
    ("ms", po::value<int>(&ms)->default_value(1),
     "NV center spin state used with |0> for an effective two-level system (+/-1)")
    ("k_DD", po::value<uint>(&k_DD)->default_value(1),
     "resonance harmonic used in spin addressing (1 or 3)")
    ("static_Bz", po::value<double>(&static_Bz_in_gauss)->default_value(140.1,"140.1"),
     "strength of static magnetic field along the NV axis (gauss)")
    ("scale_factor", po::value<double>(&scale_factor)->default_value(100,"100"),
     "factor used to define different scales (i.e. if a << b, then a = b/scale_factor)")
    ;

  uint scan_bins;
  double f_DD;
  double scan_time;
  double scan_time_in_ms;

  po::options_description scan_options("Coherence scanning options",help_text_length);
  scan_options.add_options()
    ("scan_bins", po::value<uint>(&scan_bins)->default_value(100),
     "number of bins in coherence scanning range")
    ("f_DD", po::value<double>(&f_DD)->default_value(0.06,"0.06"),
     "magnitude of fourier component used in coherence scanning")
    ("scan_time", po::value<double>(&scan_time_in_ms)->default_value(1),
     "time for each coherence measurement (microseconds)")
    ;

  po::options_description all("Allowed options");
  all.add(general).add(file_io).add(simulations).add(simulation_options).add(scan_options);

  // collect inputs
  po::variables_map inputs;
  po::store(parse_command_line(arg_num, arg_vec, all), inputs);
  po::notify(inputs);

  // if requested, print help text
  if(inputs.count("help")){
    cout << all;
    return 0;
  }

  // determine whether certain options were used
  bool using_input_lattice = inputs.count("lattice_file");
  bool set_output_suffix = inputs.count("output_suffix");
  bool set_hyperfine_cutoff = !inputs["hyperfine_cutoff_in_kHz"].defaulted();
  bool set_c13_abundance = !inputs["c13_abundance"].defaulted();

  // run a sanity check on inputs
  if(!testing){
    assert(int(pair_search)+int(coherence_scan)+int(compute_iswap_fidelities) == 1);
  }
  assert(!(using_input_lattice && set_hyperfine_cutoff));
  assert(!(using_input_lattice && set_c13_abundance));
  assert(hyperfine_cutoff_in_kHz > 0);
  assert(c13_abundance >= 0 && c13_abundance <= 1);

  assert(max_cluster_size > 0);
  assert(ms == 1 || ms == -1);
  assert((k_DD == 1) || (k_DD == 3));
  assert(scale_factor > 10);

  if(coherence_scan){
    assert(scan_bins > 0);
    assert(scan_time_in_ms > 0);
  }

  assert(seed > 0); // seeds of 0 and 1 give the same result

  // set some variables based on iputs
  hyperfine_cutoff = hyperfine_cutoff_in_kHz*kHz;
  scan_time = scan_time_in_ms*1e-3;

  // define path of lattice file defining system configuration
  fs::path lattice_path;
  if(using_input_lattice){
    if(!fs::exists(lattice_file)){
      cout << "file does not exist: " << lattice_file << endl;
      return 2;
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

  // initialize nv_system object
  nv_system nv(ms, k_DD, static_Bz_in_gauss*gauss, scale_factor);

  // -----------------------------------------------------------------------------------------
  // Construct lattice of nuclei
  // -----------------------------------------------------------------------------------------

  if(!using_input_lattice){ // place nuclei at lattice sites

    // set positions of nuclei at lattice sites
    for(int b = 0; b <= 1; b++){
      for(int l = -2*cell_radius; l <= 2*cell_radius; l++){
        for(int m = -2*cell_radius; m <= 2*cell_radius; m++){
          for(int n = -2*cell_radius; n <= 2*cell_radius; n++){
            if(rnd() < c13_abundance){ // if we pass a check for C-13 isotopic abundance
              const spin nucleus = spin(b*ao+l*a1+m*a2+n*a3, gC13, s_vec/2);
              // only place C-13 nuclei with a hyperfine field strength above the cutoff
              if(A(nv,nucleus).norm() > hyperfine_cutoff){
                nv.nuclei.push_back(nucleus);
              }
            }
          }
        }
      }
    }

    // remove any nuclei at the NV lattice sites
    for(uint i = 0; i < nv.nuclei.size(); i++){
      if((nv.nuclei.at(i).pos == nv.n.pos) || (nv.nuclei.at(i).pos == nv.e.pos)){
        nv.nuclei.erase(nv.nuclei.begin()+i);
        i--;
      }
    }

    cout << "Placed " << nv.nuclei.size() << " C-13 nuclei\n\n";

    // write cell radius and nucleus positions to file
    if(!no_output && !pair_search){
      ofstream lattice(lattice_path.string());
      lattice << "# cell radius: " << cell_radius << endl;
      for(uint i = 0; i < nv.nuclei.size(); i++){
        lattice << nv.nuclei.at(i).pos(0) << ' '
                << nv.nuclei.at(i).pos(1) << ' '
                << nv.nuclei.at(i).pos(2) << endl;
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
      nv.nuclei.push_back(spin((Vector3d() << x,y,z).finished(),gC13,s_vec/2));
    }
    lattice.close();

    // assert that no C-13 nuclei lie at the NV lattice sites
    for(uint i = 0; i < nv.nuclei.size(); i++){
      if((nv.nuclei.at(i).pos == nv.n.pos) || (nv.nuclei.at(i).pos == nv.e.pos)){
        cout << "input lattice places a C-13 nucleus at one of the NV lattice sites!\n";
        return 3;
      }
    }
  }

  // -----------------------------------------------------------------------------------------
  // Perform search for larmor pairs
  // -----------------------------------------------------------------------------------------

  if(pair_search){
    uint pairs_found = 0;
    for(uint i = 0; i < nv.nuclei.size(); i++){
      for(uint j = i+1; j < nv.nuclei.size(); j++){
        if(larmor_pair(nv.nuclei.at(i),nv.nuclei.at(j))){
          cout << "larmor pair: " << i << ", " << j << endl;
          pairs_found++;
        }
      }
    }

    if(!pairs_found){ cout << "no larmor pairs found\n"; }
    return pairs_found;
  }

  // -----------------------------------------------------------------------------------------
  // Cluster C-13 nuclei
  // -----------------------------------------------------------------------------------------

  const double cluster_coupling_guess = 100; // this value doesn't really matter
  const double dcc_cutoff = 1e-5; // cutoff for tuning of cluster_coupling

  uint cluster_size_target = max_cluster_size;
  do{
    assert(cluster_size_target > 0);
    // get cluster_coupling for which the largest cluster size is >= cluster_size_target
    nv.cluster_coupling = find_target_coupling(nv.nuclei, cluster_coupling_guess,
                                               cluster_size_target, dcc_cutoff);
    nv.clusters = get_index_clusters(nv.nuclei, nv.cluster_coupling);
    // if (largest cluster size > cluster_size_target),
    //   which can occur when (largest cluster size == cluster_size_target) is impossible,
    //   find largest cluster_coupling for which (largest cluster size < cluster_size_target)
    while(largest_cluster_size(nv.clusters) > cluster_size_target){
      int cluster_size_target
        = largest_cluster_size(get_index_clusters(nv.nuclei, nv.cluster_coupling+dcc_cutoff));
      nv.cluster_coupling = find_target_coupling(nv.nuclei, nv.cluster_coupling,
                                                 cluster_size_target, dcc_cutoff);
      nv.clusters = get_index_clusters(nv.nuclei, nv.cluster_coupling);
    }

    // if we are going to perform an actual simulation instead of just a coherence scan,
    //   we want to group together clusters sharing nuclei with similar larmor frequencies.
    if(!coherence_scan) nv.clusters = group_clusters(nv);
    cluster_size_target--; // decrease cluster_size_target every time we perform this loop

    // grouping together clusters might create a cluster greater than max_cluster_size,
    //   so we check to make sure that we do not go over the limit
  } while(largest_cluster_size(nv.clusters) > max_cluster_size);

  if(max_cluster_size > 1){
    cout << "Nuclei grouped into " << nv.clusters.size() << " clusters"
         << " with a coupling factor of "  << nv.cluster_coupling << " Hz\n";

    // collect and print histogram of cluster sizes
    max_cluster_size = largest_cluster_size(nv.clusters);
    vector<uint> size_hist(max_cluster_size);
    for(uint i = 0; i < nv.clusters.size(); i++){
      size_hist.at(nv.clusters.at(i).size()-1) += 1;
    }
    cout << "Cluster size histogram:\n";
    for(uint i = 0; i < size_hist.size(); i++){
      cout << "  " << i+1 << ": " << size_hist.at(i) << endl;
    }
    cout << endl;
  } else {
    cout << "Largest internuclear coupling: " << nv.cluster_coupling << " Hz\n";
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

  if(coherence_scan){
    // identify effictive larmor frequencies and NV coupling strengths
    vector<double> w_larmor(nv.nuclei.size());
    vector<double> A_perp(nv.nuclei.size());

    double w_max = 0, w_min = DBL_MAX; // maximum and minimum effective larmor frequencies
    for(uint i = 0; i < nv.nuclei.size(); i++){
      const Vector3d A_i = A(nv,i);
      A_perp.at(i) = (A_i-dot(A_i,zhat)*zhat).norm();
      w_larmor.at(i) = effective_larmor(nv,i).norm();

      if(w_larmor.at(i) < w_min) w_min = w_larmor.at(i);
      if(w_larmor.at(i) > w_max) w_max = w_larmor.at(i);
    }

    // print effective larmor frequencies and NV couping strengths to output file
    if(!no_output){
      fs::path larmor_path = output_dir/fs::path("larmor-"+output_suffix);
      ofstream larmor_file(larmor_path.string());
      larmor_file << "# w_larmor A_perp\n";
      for(uint i = 0; i < nv.nuclei.size(); i++){
        larmor_file << w_larmor.at(i) << " " << A_perp.at(i) << endl;
      }
      larmor_file.close();
    }

    // perform coherence scan
    cout << "Beginning coherence scan\n";
    vector<double> w_scan(scan_bins);
    vector<double> coherence(scan_bins);

    const double w_range = w_max - w_min;
    const double w_start = max(w_min - w_range/10, 0.);
    const double w_end = w_max + w_range/10;
    for(uint i = 0; i < scan_bins; i++){
      w_scan.at(i) = w_start + i*(w_end-w_start)/scan_bins;
      coherence.at(i) = coherence_measurement(nv, w_scan.at(i), f_DD, scan_time);
      cout << "(" << i+1 << "/" << scan_bins << ") "
           << w_scan.at(i)/(2*pi*1e3) << " " << coherence.at(i) << endl;
    }

    // print coherence scan results to output file
    if(!no_output){
      fs::path scan_path = output_dir/fs::path("scan-"+output_suffix);
      ofstream scan_file(scan_path.string());
      scan_file << "# cluster coupling factor (Hz): " << nv.cluster_coupling << endl;
      scan_file << "# f_DD: " << f_DD << endl;
      scan_file << "# scan time: " << scan_time << endl;
      scan_file << endl;
      scan_file << "# w_scan coherence\n";
      for(uint i = 0; i < scan_bins; i++){
        scan_file << w_scan.at(i) << " " << coherence.at(i) << endl;
      }
      scan_file.close();
    }
  }

  // -----------------------------------------------------------------------------------------
  // NV/nucleus SWAP fidelity
  // -----------------------------------------------------------------------------------------

  if(compute_iswap_fidelities){

    // collect iswap_fidelity summaries for individually addressable nuclei
    vector<uint> addressable_targets;
    vector<fidelity_info> iswap_summaries;

    for(uint target = 0; target < nv.nuclei.size(); target++){
      const fidelity_info target_info = iswap_fidelity(nv,target);
      cout << "(" << target << "/" << nv.nuclei.size() << ") ";
      if(target_info.valid){
        addressable_targets.push_back(target);
        iswap_summaries.push_back(target_info);
        cout << target_info.operation_time*1e3 << " " << target_info.fidelity << endl;
      } else{
        cout << "---\n";
      }
    }

    if(!no_output){
      fs::path fidelity_path = output_dir/fs::path("fidelities-"+output_suffix);
      ofstream fidelity_file(fidelity_path.string());
      fidelity_file << "# cluster coupling factor: " << nv.cluster_coupling << endl;
      fidelity_file << "# scale factor: " << nv.scale_factor << endl;
      fidelity_file << endl;
      fidelity_file << "# target_index larmor_eff hyperfine hyperfine_perp dw_min f_DD"
                    << " operation_time iswap_idelity\n";
      for(uint t = 0; t < addressable_targets.size(); t++){
        fidelity_file << addressable_targets.at(t) << " "
                      << iswap_summaries.at(t).larmor_eff << " "
                      << iswap_summaries.at(t).hyperfine << " "
                      << iswap_summaries.at(t).hyperfine_perp << " "
                      << iswap_summaries.at(t).dw_min << " "
                      << iswap_summaries.at(t).f_DD << " "
                      << iswap_summaries.at(t).operation_time << " "
                      << iswap_summaries.at(t).fidelity << endl;
      }
      fidelity_file.close();
    }
  }


  // -----------------------------------------------------------------------------------------
  // Testing code
  // -----------------------------------------------------------------------------------------

  uint target_cluster;
  for(uint c = 0; c < nv.clusters.size(); c++){
    if(nv.clusters.at(c).size() == 1){
      target_cluster = c;
      break;
    }
  }
  uint index = nv.clusters.at(target_cluster).at(0);

  Vector3d axis = yhat;
  MatrixXcd U_full = U_ctl(nv,index,-pi/2,axis);

  MatrixXcd U_desired = dot(s_vec,axis);
  cout << "desired propagator:\n";
  cout << U_desired << endl << endl;

  uint start = 0;
  if(U_full.rows() > 2 && abs(U_full(0,2)) > abs(U_full(0,0))){
    start = 2;
  }

  MatrixXcd U = U_full.block<2,2>(0,start);
  U /= sqrt(real(trace(U.adjoint()*U)/double(U.rows())));

  cout << "reduced propagator:\n";
  cout << clean(U) << endl << endl;

  cout << "gate fidelity: " << gate_fidelity(U,U_desired) << endl;
}
