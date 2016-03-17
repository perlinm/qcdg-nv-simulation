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
#include "nv-gates.h"
#include "nv-control.h"

int main(const int arg_num, const char *arg_vec[]) {

  // -----------------------------------------------------------------------------------------
  // Parse and process input options
  // -----------------------------------------------------------------------------------------

  const uint help_text_length = 95;

  unsigned long long int seed;

  po::options_description general("General options", help_text_length);
  general.add_options()
    ("help,h", "produce help message")
    ("seed", po::value<unsigned long long int>(&seed)->default_value(0),
     "seed for random number generator")
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
  bool single_control;
  bool single_coupling;
  bool iswap_fidelities;
  bool swap_fidelities;
  bool swap_nvst_fidelity;
  bool testing;

  po::options_description simulations("Available simulations",help_text_length);
  simulations.add_options()
    ("pairs", po::value<bool>(&pair_search)->default_value(false)->implicit_value(true),
     "search for larmor pairs")
    ("scan", po::value<bool>(&coherence_scan)->default_value(false)->implicit_value(true),
     "perform coherence scan of effective larmor frequencies")
    ("control", po::value<bool>(&single_control)->default_value(false)->implicit_value(true),
     "control individual nucleus")
    ("couple", po::value<bool>(&single_coupling)->default_value(false)->implicit_value(true),
     "couple individual nucleus to NV center")
    ("iswap",
     po::value<bool>(&iswap_fidelities)->default_value(false)->implicit_value(true),
     "compute expected iSWAP fidelities")
    ("swap",
     po::value<bool>(&swap_fidelities)->default_value(false)->implicit_value(true),
     "compute expected SWAP fidelities")
    ("swap_nvst",
     po::value<bool>(&swap_nvst_fidelity)->default_value(false)->implicit_value(true),
     "compute expected SWAP_NVST fidelity")
    ("test" ,po::value<bool>(&testing)->default_value(false)->implicit_value(true),
     "enable testing mode")
    ;

  double c13_abundance;
  double c13_percentage;
  uint max_cluster_size;
  double hyperfine_cutoff;
  double hyperfine_cutoff_in_kHz;
  int cell_radius; // determined by hyperfine_cutoff
  int ms;
  uint k_DD_int;
  axy_harmonic k_DD;
  double static_Bz_in_gauss;
  double scale_factor;
  double integration_factor;

  po::options_description simulation_options("Simulation options",help_text_length);
  simulation_options.add_options()
    ("c13_abundance", po::value<double>(&c13_percentage)->default_value(1.07,"1.07"),
     "relative isotopic abundance of C-13 (percentage)")
    ("max_cluster_size", po::value<uint>(&max_cluster_size)->default_value(6),
     "maximum allowable size of C-13 clusters")
    ("hyperfine_cutoff", po::value<double>(&hyperfine_cutoff_in_kHz)->default_value(100),
     "set cutoff scale for hyperfine field (kHz)")
    ("ms", po::value<int>(&ms)->default_value(1),
     "NV center spin state used with |0> for an effective two-level system (+/-1)")
    ("k_DD", po::value<uint>(&k_DD_int)->default_value(1),
     "resonance harmonic used in spin addressing (1 or 3)")
    ("static_Bz", po::value<double>(&static_Bz_in_gauss)->default_value(140.1,"140.1"),
     "strength of static magnetic field along the NV axis (gauss)")
    ("scale_factor", po::value<double>(&scale_factor)->default_value(10),
     "factor used to define different scales (i.e. if a << b, then a = b/scale_factor)")
    ("integration_factor", po::value<double>(&integration_factor)->default_value(10000),
     "factor used to determine size of integration step size")
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

  vector<uint> target_nuclei;
  double phase;
  double phase_over_pi;
  double target_polar;
  double target_polar_over_pi;
  double target_azimuth;
  double target_azimuth_over_pi;
  double nv_polar;
  double nv_polar_over_pi;
  double nv_azimuth;
  double nv_azimuth_over_pi;

  po::options_description addressing_options("Nucleus addressing options", help_text_length);
  addressing_options.add_options()
    ("targets", po::value<vector<uint>>(&target_nuclei)->multitoken(),
     "indices of nuclei to target (if applicable)")
    ("phase", po::value<double>(&phase_over_pi)->default_value(0.5,"0.5"),
     "operation phase in units of pi")
    ("target_polar", po::value<double>(&target_polar_over_pi)->default_value(0.5,"0.5"),
     "polar angle of target rotation axis (in units of pi radians)")
    ("target_azimuth", po::value<double>(&target_azimuth_over_pi)->default_value(0),
     "azimuthal angle of target rotation axis (in units of pi radians)")
    ("nv_polar", po::value<double>(&nv_polar_over_pi)->default_value(0),
     "polar angle of NV rotation axis (in units of pi radians)")
    ("nv_azimuth", po::value<double>(&nv_azimuth_over_pi)->default_value(0),
     "azimuthal angle of NV rotation axis (in units of pi radians)")
    ;

  po::options_description all("Allowed options");
  all.add(general);
  all.add(file_io);
  all.add(simulations);
  all.add(simulation_options);
  all.add(scan_options);
  all.add(addressing_options);

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
  bool set_hyperfine_cutoff = !inputs["hyperfine_cutoff"].defaulted();
  bool set_c13_abundance = !inputs["c13_abundance"].defaulted();
  bool set_target_nuclei = inputs.count("targets");

  // run a sanity check on inputs
  if(!testing){
    assert(int(pair_search)
           + int(coherence_scan)
           + int(single_control)
           + int(single_coupling)
           + int(iswap_fidelities)
           + int(swap_fidelities)
           + int(swap_nvst_fidelity)
           == 1);
  }
  assert(!(using_input_lattice && set_hyperfine_cutoff));
  assert(!(using_input_lattice && set_c13_abundance));
  assert(hyperfine_cutoff_in_kHz > 0);
  assert(c13_percentage >= 0 && c13_percentage <= 100);

  assert(max_cluster_size > 0);
  assert(ms == 1 || ms == -1);
  assert((k_DD_int == 1) || (k_DD_int == 3));
  assert(scale_factor > 1);
  assert(integration_factor > 1);

  if(coherence_scan){
    assert(scan_bins > 0);
    assert(scan_time_in_ms > 0);
  }

  // set some variables based on iputs
  c13_abundance = c13_percentage/100;
  hyperfine_cutoff = hyperfine_cutoff_in_kHz*kHz;
  k_DD = (k_DD_int == 1 ? first : third);
  scan_time = scan_time_in_ms*1e-3;

  phase = phase_over_pi*pi;
  target_polar = target_polar_over_pi*pi;
  target_azimuth = target_azimuth_over_pi*pi;
  nv_polar = nv_polar_over_pi*pi;
  nv_azimuth = nv_azimuth_over_pi*pi;

  // define path of lattice file defining system configuration
  fs::path lattice_path;
  if(using_input_lattice){
    if(!fs::exists(lattice_file)){
      cout << "file does not exist: " << lattice_file << endl;
      return -1;
    }
    lattice_path = lattice_file;

    if(!set_output_suffix){
      output_suffix = lattice_path.stem().string() + "-" + output_suffix_with_input_lattice;
    }

  } else{ // if !using_input_lattice
    cell_radius = round(pow(abs(ge*gC13)/(4*pi*a0*a0*a0*hyperfine_cutoff),1.0/3));
    cout << "Setting cell radius to: " << cell_radius << endl;

    boost::replace_all(lattice_file, "[cell_radius]", to_string(cell_radius));
    boost::replace_all(lattice_file, "[seed]", to_string(seed));
    lattice_path = output_dir/fs::path(lattice_file);
  }

  uniform_real_distribution<double> rnd(0.0,1.0); // uniform distribution on the range [0,1)
  mt19937_64 generator(seed); // use and seed the 64-bit Mersenne Twister 19937 generator

  fs::create_directory(output_dir); // create data directory

  // initialize nv_system object
  nv_system nv(ms, static_Bz_in_gauss*gauss, k_DD, scale_factor, integration_factor);

  // -----------------------------------------------------------------------------------------
  // Construct lattice of nuclei
  // -----------------------------------------------------------------------------------------

  if(!using_input_lattice){ // place nuclei at lattice sites
    // set positions of nuclei at lattice sites
    for(uint b: {0,1}){
      for(int l = -2*cell_radius; l <= 2*cell_radius; l++){
        for(int m = -2*cell_radius; m <= 2*cell_radius; m++){
          for(int n = -2*cell_radius; n <= 2*cell_radius; n++){
            if(rnd(generator) < c13_abundance){ // pass a check for C-13 isotopic abundance
              if(l != 0 || m != 0 || n != 0){ // don't place C-13 nucleus on NV lattice site
                const spin nucleus(b*ao+l*a1+m*a2+n*a3, gC13, s_vec/2);
                // only place C-13 nuclei with a hyperfine field strength above the cutoff
                if(hyperfine(nv,nucleus).norm() > hyperfine_cutoff){
                  nv.nuclei.push_back(nucleus);
                }
              }
            }
          }
        }
      }
    }
    cout << "Placed " << nv.nuclei.size() << " C-13 nuclei\n";
    if(nv.nuclei.size() == 0) return 0;
    if(!set_target_nuclei){
      for(uint i = 0; i < nv.nuclei.size(); i++){
        target_nuclei.push_back(i);
      }
    }

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
        return -2;
      }
    }
  }

  cout << endl;

  // -----------------------------------------------------------------------------------------
  // Perform search for larmor pairs
  // -----------------------------------------------------------------------------------------

  if(pair_search){
    uint pairs_found = 0;
    for(uint i = 0; i < nv.nuclei.size(); i++){
      for(uint j = i+1; j < nv.nuclei.size(); j++){
        if(is_larmor_pair(nv,i,j)){
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

  const double cluster_coupling_guess = 100; // this value doesn't actually matter
  const double dcc_cutoff = 1e-5; // cutoff for tuning of cluster_coupling

  // if we are going to perform an actual simulation instead of just a coherence scan,
  //   we want to group together clusters sharing nuclei with similar larmor frequencies.
  const bool group_larmor_pairs = !coherence_scan && !testing;
  if(group_larmor_pairs){
    nv.clusters = cluster_nuclei(nv.nuclei, DBL_MAX);
    nv.clusters = group_clusters(nv);
    const uint min_cluster_size_cap = largest_cluster_size(nv.clusters);
    cout << "The minimum cluster size cap is " << min_cluster_size_cap << endl;
    if(max_cluster_size < min_cluster_size_cap) return -3;
  }

  uint cluster_size_target = min(max_cluster_size,uint(nv.nuclei.size()));
  do{
    assert(cluster_size_target > 0);
    // get cluster_coupling for which the largest cluster size is >= cluster_size_target
    nv.cluster_coupling = find_target_coupling(nv.nuclei, cluster_coupling_guess,
                                               cluster_size_target, dcc_cutoff);
    nv.clusters = cluster_nuclei(nv.nuclei, nv.cluster_coupling);
    // if (largest cluster size > cluster_size_target),
    //   which can occur when (largest cluster size == cluster_size_target) is impossible,
    //   find largest cluster_coupling for which (largest cluster size < cluster_size_target)
    while(largest_cluster_size(nv.clusters) > cluster_size_target){
      const uint cluster_size_target
        = largest_cluster_size(cluster_nuclei(nv.nuclei, nv.cluster_coupling+dcc_cutoff));
      nv.cluster_coupling = find_target_coupling(nv.nuclei, nv.cluster_coupling,
                                                 cluster_size_target, dcc_cutoff);
      nv.clusters = cluster_nuclei(nv.nuclei, nv.cluster_coupling);
    }

    // if we are going to perform an actual simulation instead of just a coherence scan,
    //   we want to group together clusters sharing nuclei with similar larmor frequencies.
    if(group_larmor_pairs) nv.clusters = group_clusters(nv);
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
  boost::replace_all(output_suffix, "[k_DD]", to_string(k_DD_int));
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
      const Vector3d A_i = hyperfine(nv,i);
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
  // Individual addressing -- control
  // -----------------------------------------------------------------------------------------

  if(single_control){
    for(uint target: target_nuclei){
      const Vector3d rotation = 2*phase * axis(target_polar,target_azimuth);
      vector<MatrixXcd> U(2);
      for(bool exact : {true,false}){
        U.at(exact) = rotate_target(nv, target, rotation, exact);
      }
      const uint cluster = get_cluster_containing_index(nv, target);
      const uint target_in_cluster = get_index_in_cluster(target, nv.clusters.at(cluster));



      // rotate into the frame of the nuclei
      const uint spins = nv.clusters.at(cluster).size()+1;
      MatrixXcd R = MatrixXcd::Identity(pow(2,spins),pow(2,spins));
      for(uint index: nv.clusters.at(cluster)){
        const uint index_in_cluster = get_index_in_cluster(index, nv.clusters.at(cluster));
        const Matrix2cd R_index = rotate(natural_basis(nv,index), {xhat,yhat,zhat});
        R = (act(R_index, {index_in_cluster+1},spins) * R).eval();
      }
      const MatrixXcd U_e = R.adjoint() * U.at(true) * R;
      const MatrixXcd U_a = R.adjoint() * U.at(false) * R;
      const double err = 1e-2;
      cout << endl;
      U_print(j*log(U_e)/pi, err);
      cout << endl << endl;
      U_print(j*log(U_a)/pi, err);
      cout << endl;




      cout << target << ": " << gate_fidelity(U.at(0), U.at(1), {target_in_cluster}) << endl;
    }
  }

  // -----------------------------------------------------------------------------------------
  // Individual addressing -- NV coupling
  // -----------------------------------------------------------------------------------------

  if(single_coupling){
    const Vector3d nv_axis = axis(nv_polar, nv_azimuth);
    for(uint target: target_nuclei){
      const Vector3d target_axis = axis(target_polar, target_azimuth);
      vector<MatrixXcd> U(2);
      for(bool exact : {true,false}){
        U.at(exact) = couple_target(nv, target, phase, nv_axis, target_axis, exact);
      }
      const uint cluster = get_cluster_containing_index(nv, target);
      const uint target_in_cluster = get_index_in_cluster(target, nv.clusters.at(cluster));



      // rotate into the frame of the nuclei
      const uint spins = nv.clusters.at(cluster).size()+1;
      MatrixXcd R = MatrixXcd::Identity(pow(2,spins),pow(2,spins));
      for(uint index: nv.clusters.at(cluster)){
        const uint index_in_cluster = get_index_in_cluster(index, nv.clusters.at(cluster));
        const Matrix2cd R_index = rotate(natural_basis(nv,index), {xhat,yhat,zhat});
        R = (act(R_index, {index_in_cluster+1},spins) * R).eval();
      }
      const MatrixXcd U_e = R.adjoint() * U.at(true) * R;
      const MatrixXcd U_a = R.adjoint() * U.at(false) * R;
      const double err = 1e-2;
      cout << endl;
      U_print(j*log(U_e)/pi, err);
      cout << endl << endl;
      U_print(j*log(U_a)/pi, err);
      cout << endl;



      cout << target << ": " << gate_fidelity(U.at(0), U.at(1), {target_in_cluster}) << endl;
    }
  }

  // -----------------------------------------------------------------------------------------
  // NV/nucleus iSWAP fidelity
  // -----------------------------------------------------------------------------------------

  if(iswap_fidelities){
    for(uint target: target_nuclei){
      vector<MatrixXcd> U(2);
      for(bool exact : {true,false}){
        U.at(exact) = iSWAP(nv, target, exact);
      }
      const uint cluster = get_cluster_containing_index(nv, target);
      const uint target_in_cluster = get_index_in_cluster(target, nv.clusters.at(cluster));
      cout << target << ": " << gate_fidelity(U.at(0), U.at(1), {target_in_cluster}) << endl;
    }
  }

  // -----------------------------------------------------------------------------------------
  // NV/nucleus SWAP fidelity
  // -----------------------------------------------------------------------------------------

  if(swap_fidelities){
    for(uint target: target_nuclei){
      vector<MatrixXcd> U(2);
      for(bool exact : {true,false}){
        U.at(exact) = SWAP(nv, target, exact);
      }
      const uint cluster = get_cluster_containing_index(nv, target);
      const uint target_in_cluster = get_index_in_cluster(target, nv.clusters.at(cluster));
      cout << target << ": " << gate_fidelity(U.at(0), U.at(1), {target_in_cluster}) << endl;
    }
  }

  // -----------------------------------------------------------------------------------------
  // NV/ST SWAP operation fidelity
  // -----------------------------------------------------------------------------------------

  if(swap_nvst_fidelity){
    assert(nv.nuclei.size() >= 2);
    assert(target_nuclei.size() >= 2);
    const uint idx1 = target_nuclei.at(0);
    const uint idx2 = target_nuclei.at(1);
    vector<MatrixXcd> U(2);
    for(bool exact : {true,false}){
      U.at(exact) = SWAP_NVST(nv, idx1, idx2, exact);
    }
    const uint cluster = get_cluster_containing_index(nv,idx1);
    const uint idx1_in_cluster = get_index_in_cluster(idx1,nv.clusters.at(cluster));
    const uint idx2_in_cluster = get_index_in_cluster(idx2,nv.clusters.at(cluster));
    cout << idx1 << " " << idx2 << ": "
         << gate_fidelity(U.at(0), U.at(1), {idx1_in_cluster, idx2_in_cluster})
         << endl;
  }

}
