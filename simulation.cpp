#define EIGEN_USE_MKL_ALL

#include <iostream> // for standard output
#include <fstream> // for file input
#include <iomanip> // some nice printing functions
#include <random> // for randomness
#include <boost/algorithm/string.hpp> // string manipulation library
#include <boost/filesystem.hpp> // filesystem path manipulation library
#include <boost/program_options.hpp> // options parsing library
#include <eigen3/Eigen/Dense> // linear algebra library

#include "constants.h"
#include "qp-math.h"
#include "gates.h"
#include "nv-math.h"
#include "nv-control.h"

using namespace std;
using namespace Eigen;
namespace fs = boost::filesystem;
namespace po = boost::program_options;

int main(const int arg_num, const char *arg_vec[]) {

  // -------------------------------------------------------------------------------------
  // Parse and process input options
  // -------------------------------------------------------------------------------------

  const uint help_text_length = 90;

  unsigned long long int seed;

  po::options_description general("General options", help_text_length);
  general.add_options()
    ("help,h", "produce help message")
    ("seed", po::value<unsigned long long int>(&seed)->default_value(0),
     "seed for random number generator")
    ;

  bool print_pairs;
  bool coherence_scan;
  bool single_control;
  bool single_coupling;
  bool iswap_fidelities;
  bool swap_fidelities;
  bool swap_nvst_fidelity;
  bool testing;

  po::options_description simulations("Available simulations",help_text_length);
  simulations.add_options()
    ("pairs", po::value<bool>(&print_pairs)->default_value(false)->implicit_value(true),
     "search for larmor pairs")
    ("scan", po::value<bool>(&coherence_scan)->default_value(false)->implicit_value(true),
     "perform coherence scan for effective larmor frequencies")
    ("control",
     po::value<bool>(&single_control)->default_value(false)->implicit_value(true),
     "control individual nucleus")
    ("couple",
     po::value<bool>(&single_coupling)->default_value(false)->implicit_value(true),
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
  int cell_radius = 0; // requires a default value which will be overwritten later
  int ms;
  uint k_DD_int;
  axy_harmonic k_DD;
  double static_Bz_in_gauss;
  double scale_factor;
  double integration_factor;
  bool no_nn;

  po::options_description simulation_options("Simulation options",help_text_length);
  simulation_options.add_options()
    ("c13_percentage", po::value<double>(&c13_percentage)->default_value(1.07,"1.07"),
     "unweigheted isotopic abundance of C-13 as a percentage")
    ("max_cluster_size", po::value<uint>(&max_cluster_size)->default_value(6),
     "maximum allowable size of C-13 clusters")
    ("hyperfine_cutoff", po::value<double>(&hyperfine_cutoff_in_kHz)->default_value(10),
     "set cutoff scale for hyperfine field (kHz)")
    ("ms", po::value<int>(&ms)->default_value(1),
     "NV center spin state used with |0> for an effective two-level system (+/-1)")
    ("k_DD", po::value<uint>(&k_DD_int)->default_value(1),
     "resonance harmonic used in spin addressing (1 or 3)")
    ("static_Bz", po::value<double>(&static_Bz_in_gauss)->default_value(140.1,"140.1"),
     "strength of static magnetic field along the NV axis (gauss)")
    ("scale_factor", po::value<double>(&scale_factor)->default_value(10),
     "factor used to define different scales (i.e. if a << b, then a = b/scale_factor)")
    ("integration_factor", po::value<double>(&integration_factor)->default_value(10),
     "factor used to determine size of integration step size")
    ("no_nn" ,po::value<bool>(&no_nn)->default_value(false)->implicit_value(true),
     "turn off internuclear couplings")
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

  po::options_description addressing_options("Nuclear spin addressing options"
                                             " (all angles are in units of pi radians)",
                                             help_text_length);
  addressing_options.add_options()
    ("targets", po::value<vector<uint>>(&target_nuclei)->multitoken(),
     "indices of nuclei to target")
    ("phase", po::value<double>(&phase_over_pi)->default_value(1),
     "operation phase")
    ("target_polar", po::value<double>(&target_polar_over_pi)->default_value(0.5),
     "polar angle of target rotation axis")
    ("target_azimuth", po::value<double>(&target_azimuth_over_pi)->default_value(0),
     "azimuthal angle of target rotation axis")
    ("nv_polar", po::value<double>(&nv_polar_over_pi)->default_value(0),
     "polar angle of NV rotation axis")
    ("nv_azimuth", po::value<double>(&nv_azimuth_over_pi)->default_value(0),
     "azimuthal angle of NV rotation axis")
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

  po::options_description all("Allowed options");
  all.add(general);
  all.add(simulations);
  all.add(simulation_options);
  all.add(addressing_options);
  all.add(scan_options);
  all.add(file_io);

  // collect inputs
  po::variables_map inputs;
  po::store(parse_command_line(arg_num, arg_vec, all), inputs);
  po::notify(inputs);

  // if requested, print help text
  if (inputs.count("help")) {
    cout << all;
    return 0;
  }

  // determine whether certain options were used
  bool using_input_lattice = inputs.count("lattice_file");
  bool set_output_suffix = inputs.count("output_suffix");
  bool set_hyperfine_cutoff = !inputs["hyperfine_cutoff"].defaulted();
  bool set_c13_abundance = !inputs["c13_percentage"].defaulted();
  bool set_target_nuclei = inputs.count("targets");

  // run a sanity check on inputs
  if (!testing && (int(print_pairs)
                   + int(coherence_scan)
                   + int(single_control)
                   + int(single_coupling)
                   + int(iswap_fidelities)
                   + int(swap_fidelities)
                   + int(swap_nvst_fidelity)
                   != 1)) {
    cout << "Please choose one simulation to perform\n";
    return -1;
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

  if (coherence_scan) {
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
  if (using_input_lattice) {
    if (!fs::exists(lattice_file)) {
      cout << "File does not exist: " << lattice_file << endl;
      return -1;
    }
    lattice_path = lattice_file;

    if (!set_output_suffix) {
      output_suffix =
        lattice_path.stem().string() + "-" + output_suffix_with_input_lattice;
    }

  } else { // if !using_input_lattice
    cell_radius = round(pow(abs(g_e*g_C13)/(4*pi*a0*a0*a0*hyperfine_cutoff),1.0/3));
    cout << "Setting cell radius to: " << cell_radius << endl;

    boost::replace_all(lattice_file, "[cell_radius]", to_string(cell_radius));
    boost::replace_all(lattice_file, "[seed]", to_string(seed));
    lattice_path = output_dir/fs::path(lattice_file);
  }

  uniform_real_distribution<double> rnd(0.0,1.0); // uniform distribution on [0,1)
  mt19937_64 generator(seed); // use and seed the 64-bit Mersenne Twister 19937 generator

  fs::create_directory(output_dir); // create data directory

  // -------------------------------------------------------------------------------------
  // Construct lattice of nuclei
  // -------------------------------------------------------------------------------------

  vector<Vector3d> nuclei;

  if (!using_input_lattice) { // place nuclei at lattice sites
    // set positions of nuclei at lattice sites
    for (uint b: {0,1}) {
      for (int l = -2*cell_radius; l <= 2*cell_radius; l++) {
        for (int m = -2*cell_radius; m <= 2*cell_radius; m++) {
          for (int n = -2*cell_radius; n <= 2*cell_radius; n++) {
            if (rnd(generator) < c13_abundance) { // check for C-13 isotopic abundance
              if (l != 0 || m != 0 || n != 0) { // don't place C-13 nucleus on NV sites
                const Vector3d pos = b*ao+l*a1+m*a2+n*a3;
                // only place C-13 nuclei with a hyperfine field strength above the cutoff
                if (hyperfine(pos).norm() > hyperfine_cutoff) {
                  nuclei.push_back(pos);
                }
              }
            }
          }
        }
      }
    }
    cout << "Placed " << nuclei.size() << " C-13 nuclei\n\n";
    if (nuclei.size() == 0) return 0;

    // write cell radius and nucleus positions to file
    if (!no_output && !print_pairs) {
      ofstream lattice(lattice_path.string());
      lattice << "# cell radius: " << cell_radius << endl;
      for (uint i = 0; i < nuclei.size(); i++) {
        lattice << nuclei.at(i)(0) << ' '
                << nuclei.at(i)(1) << ' '
                << nuclei.at(i)(2) << endl;
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
    while (getline(lattice,line,' ')) {
      x = stod(line);
      getline(lattice,line,' ');
      y = stod(line);
      getline(lattice,line);
      z = stod(line);
      nuclei.push_back((Vector3d() << x,y,z).finished());
    }
    lattice.close();

    // assert that no C-13 nuclei lie at the NV lattice sites
    for (uint i = 0; i < nuclei.size(); i++) {
      if ((nuclei.at(i) == n_pos) || (nuclei.at(i) == e_pos)) {
        cout << "The input lattice places a C-13 nucleus at one of the NV lattice sites!"
             << endl;
        return -1;
      }
    }
  }

  // identify nuclei which cannot be addressed
  vector<uint> unaddressable_nuclei;
  for (uint n = 0; n < nuclei.size(); n++) {
    if (!can_address(nuclei,n)) {
      unaddressable_nuclei.push_back(n);
    }
  }
  if (unaddressable_nuclei.size() > 0) {
    cout << "The following nuclei cannot be addressed:";
    for (uint n: unaddressable_nuclei) {
      cout << " " << n;
    }
    cout << endl << endl;
  }

  // remove nuclei which cannot be addressed from the list of target nuclei
  if (!set_target_nuclei) {
    for (uint n = 0; n < nuclei.size(); n++) {
      if (!in_vector(n,unaddressable_nuclei)) {
        target_nuclei.push_back(n);
      }
    }
  } else { // if (set_target_nuclei)
    vector<uint> unaddressable_targets;
    for (uint n = 0; n < target_nuclei.size(); n++) {
      if (!can_address(nuclei,target_nuclei.at(n))) {
        unaddressable_targets.push_back(target_nuclei.at(n));
        target_nuclei.erase(target_nuclei.begin()+n);
        n--;
      }
    }
    if (unaddressable_targets.size() > 0) {
      cout << "(WARNING) Ignoring following target nuclei:";
      for (uint n: unaddressable_targets) {
        cout << " " << n;
      }
      cout << endl << endl;
    }
  }

  // -------------------------------------------------------------------------------------
  // Identify larmor pairs
  // -------------------------------------------------------------------------------------

  vector<vector<uint>> larmor_pairs;
  for (uint i = 0; i < target_nuclei.size(); i++) {
    const uint n_i = target_nuclei.at(i);
    for (uint j = i+1; j < target_nuclei.size(); j++) {
      const uint n_j = target_nuclei.at(j);
      if (is_larmor_pair(nuclei,n_i,n_j)) {
        larmor_pairs.push_back({n_i,n_j});
      }
    }
  }

  if (print_pairs) {
    if (larmor_pairs.size() > 0) {
      cout << "Larmor pairs:\n";
      for (vector<uint> pair: larmor_pairs) {
        cout << " " << pair.at(0) << " " << pair.at(1) << endl;
      }
    } else {
      cout << "No larmor pairs found\n";
    }
    return larmor_pairs.size();
  }

  // -------------------------------------------------------------------------------------
  // Cluster C-13 nuclei
  // -------------------------------------------------------------------------------------

  // unless we are performing a coherence scan, we will be grouping together clusters by
  //  the larmor frequencies of the nuclei, so first we check whether doing so is possible
  //  for the given min_cluster_size_cap
  const uint min_cluster_size_cap = smallest_possible_cluster_size(nuclei);
  const bool cluster_larmor_pairs =
    !(coherence_scan || (max_cluster_size < min_cluster_size_cap));
  if (!coherence_scan) {
    cout << "The minimum cluster size cap is " << min_cluster_size_cap << endl;
    if (!testing && (max_cluster_size < min_cluster_size_cap)) return -1;
  }

  const vector<vector<uint>> clusters = cluster_nuclei(nuclei, max_cluster_size,
                                                       cluster_larmor_pairs);
  const double cluster_coupling = get_cluster_coupling(nuclei,clusters);

  cout << "Nuclei grouped into " << clusters.size() << " clusters"
       << " with a coupling factor of " << cluster_coupling << " Hz\n";

  // collect and print histogram of cluster sizes
  vector<uint> size_hist(largest_cluster_size(clusters));
  for (uint i = 0; i < clusters.size(); i++) {
    size_hist.at(clusters.at(i).size()-1) += 1;
  }
  cout << "Cluster size histogram:\n";
  for (uint i = 0; i < size_hist.size(); i++) {
    cout << "  " << i+1 << ": " << size_hist.at(i) << endl;
  }
  cout << endl;

  // now that we are done with initialization, fix up the output filename suffix
  boost::replace_all(output_suffix, "[cell_radius]", to_string(cell_radius));
  boost::replace_all(output_suffix, "[seed]", to_string(seed));
  boost::replace_all(output_suffix, "[cluster_size]", to_string(max_cluster_size));
  boost::replace_all(output_suffix, "[k_DD]", to_string(k_DD_int));
  boost::replace_all(output_suffix, "[ms]", (ms > 0)?"up":"dn");

  // initialize nv_system object
  const nv_system nv(nuclei, clusters, ms, g_C13*static_Bz_in_gauss*gauss, k_DD,
                     scale_factor, integration_factor, no_nn);

  // -------------------------------------------------------------------------------------
  // Coherence scan
  // -------------------------------------------------------------------------------------

  if (coherence_scan) {
    // identify effictive larmor frequencies and NV coupling strengths
    vector<double> w_larmor(nv.nuclei.size());
    vector<double> A_perp(nv.nuclei.size());

    double w_max = 0, w_min = DBL_MAX; // maximum and minimum effective larmor frequencies
    for (uint i = 0; i < nv.nuclei.size(); i++) {
      A_perp.at(i) = hyperfine_perp(nv,i).norm();
      w_larmor.at(i) = effective_larmor(nv,i).norm();

      if (w_larmor.at(i) < w_min) w_min = w_larmor.at(i);
      if (w_larmor.at(i) > w_max) w_max = w_larmor.at(i);
    }

    // print effective larmor frequencies and NV couping strengths to output file
    if (!no_output) {
      fs::path larmor_path = output_dir/fs::path("larmor-"+output_suffix);
      ofstream larmor_file(larmor_path.string());
      larmor_file << "# w_larmor A_perp\n";
      for (uint i = 0; i < nv.nuclei.size(); i++) {
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
    for (uint i = 0; i < scan_bins; i++) {
      w_scan.at(i) = w_start + i*(w_end-w_start)/scan_bins;
      coherence.at(i) = coherence_measurement(nv, w_scan.at(i), f_DD, scan_time);
      cout << "(" << i+1 << "/" << scan_bins << ") "
           << w_scan.at(i)/(2*pi*1e3) << " " << coherence.at(i) << endl;
    }

    // print coherence scan results to output file
    if (!no_output) {
      fs::path scan_path = output_dir/fs::path("scan-"+output_suffix);
      ofstream scan_file(scan_path.string());
      scan_file << "# cluster coupling factor (Hz): " << cluster_coupling << endl;
      scan_file << "# f_DD: " << f_DD << endl;
      scan_file << "# scan time: " << scan_time << endl;
      scan_file << endl;
      scan_file << "# w_scan coherence\n";
      for (uint i = 0; i < scan_bins; i++) {
        scan_file << w_scan.at(i) << " " << coherence.at(i) << endl;
      }
      scan_file.close();
    }
  }

  // -------------------------------------------------------------------------------------
  // Individual addressing -- control
  // -------------------------------------------------------------------------------------

  if (single_control) {
    cout << "target fidelity time pulses\n";
    for (uint target: target_nuclei) {
      const Vector3d target_axis = axis(target_polar,target_azimuth);
      vector<protocol> P(2);
      for (bool exact : {true,false}) {
        P.at(exact) = rotate_target(nv, target, phase, target_axis, exact);
      }
      const uint target_in_cluster = get_index_in_cluster(nv, target);
      cout << target << " "
           << gate_fidelity(P.at(0), P.at(1), {0, target_in_cluster}) << " "
           << P.at(false).time << " "
           << P.at(false).pulses << endl;
    }
  }

  // -------------------------------------------------------------------------------------
  // Individual addressing -- NV coupling
  // -------------------------------------------------------------------------------------

  if (single_coupling) {
    cout << "target fidelity time pulses\n";
    const Vector3d nv_axis = axis(nv_polar, nv_azimuth);
    for (uint target: target_nuclei) {
      const Vector3d target_axis = axis(target_polar, target_azimuth);
      vector<protocol> P(2);
      for (bool exact : {true,false}) {
        P.at(exact) = couple_target(nv, target, phase, nv_axis, target_axis, exact);
      }
      const uint target_in_cluster = get_index_in_cluster(nv, target);
      cout << target << " "
           << gate_fidelity(P.at(0), P.at(1), {0, target_in_cluster}) << " "
           << P.at(false).time << " "
           << P.at(false).pulses << endl;
    }
  }

  // -------------------------------------------------------------------------------------
  // NV/nucleus iSWAP fidelity
  // -------------------------------------------------------------------------------------

  if (iswap_fidelities) {
    cout << "target fidelity time pulses\n";
    for (uint target: target_nuclei) {
      vector<protocol> P(2);
      for (bool exact : {true,false}) {
        P.at(exact) = iSWAP(nv, target, exact);
      }
      const uint target_in_cluster = get_index_in_cluster(nv, target);
      cout << target << " "
           << gate_fidelity(P.at(0), P.at(1), {0, target_in_cluster}) << " "
           << P.at(false).time << " "
           << P.at(false).pulses << endl;
    }
  }

  // -------------------------------------------------------------------------------------
  // SWAP fidelity: NV electron spin and single nuclear spin
  // -------------------------------------------------------------------------------------

  if (swap_fidelities) {
    cout << "target fidelity time pulses\n";
    for (uint target: target_nuclei) {
      vector<protocol> P(2);
      for (bool exact : {true,false}) {
        P.at(exact) = SWAP(nv, target, exact);
      }
      const uint target_in_cluster = get_index_in_cluster(nv, target);
      cout << target << " "
           << gate_fidelity(P.at(0), P.at(1), {0, target_in_cluster}) << " "
           << P.at(false).time << " "
           << P.at(false).pulses << endl;
    }
  }

  // -------------------------------------------------------------------------------------
  // SWAP fidelity: NV electron spin and singlet-triplet subspace of two nuclear spins
  // -------------------------------------------------------------------------------------

  if (swap_nvst_fidelity) {
    if (larmor_pairs.size() == 0) {
      cout << "There are no larmor pairs in this system\n";
      return -1;
    }
    cout << "idx1 idx2 fidelity time pulses\n";
    for (vector<uint> idxs: larmor_pairs) {
      const uint idx1 = idxs.at(0);
      const uint idx2 = idxs.at(1);
      vector<protocol> P(2);
      for (bool exact : {true,false}) {
        P.at(exact) = SWAP_NVST(nv, idx1, idx2, exact);
      }
      const uint idx1_in_cluster = get_index_in_cluster(nv,idx1);
      const uint idx2_in_cluster = get_index_in_cluster(nv,idx2);
      cout << idx1 << " " << idx2 << " "
           << gate_fidelity(P.at(0), P.at(1), {0, idx1_in_cluster, idx2_in_cluster}) << " "
           << P.at(false).time << " "
           << P.at(false).pulses << endl;
    }
  }

  // -------------------------------------------------------------------------------------
  // Print info about target nuclei
  // -------------------------------------------------------------------------------------

  if (testing) {

    for (uint n: target_nuclei) {
      cout << endl
           << "index: " << n << endl
           << "position (nm): "
           << in_crystal_basis(nv.nuclei.at(n)).transpose() * a0/2 / nm << endl
           << "hyperfine (kHz): " << hyperfine(nv,n).norm() / kHz << endl
           << "hyperfine_perp (kHz): " << hyperfine_perp(nv,n).norm() / kHz << endl;
    }

    for (uint i = 0; i < target_nuclei.size(); i++) {
      const Vector3d pos_i = nv.nuclei.at(target_nuclei.at(i));
      for (uint j = i + 1; j < target_nuclei.size(); j++) {
        const Vector3d pos_j = nv.nuclei.at(target_nuclei.at(j));
        cout << endl
             << "indices: " << i << " " << j << endl
             << "displacement (nm): "
             << in_crystal_basis(pos_j-pos_i).transpose() * a0/2 / nm << endl
             << "coupling (Hz): " << coupling_strength(pos_i,pos_j) / Hz << endl;
      }
    }

  }

}
