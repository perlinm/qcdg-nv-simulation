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
  // Set input options
  // -------------------------------------------------------------------------------------

  const uint help_text_length = 90;

  unsigned long long int seed;

  po::options_description general("General options", help_text_length);
  general.add_options()
    ("help,h", "produce help message")
    ("seed", po::value<unsigned long long int>(&seed)->default_value(0),
     "seed for random number generator")
    ;

  bool coherence_scan;
  bool coherence_signal;
  bool angular_coherence_signal;
  bool rotation;
  bool coupling;
  bool iswap;
  bool swap;
  bool identity;
  bool swap_nvst;
  bool larmor_identity;
  bool initialize;
  bool initialize_x;
  bool initialize_larmor;
  bool testing_mode;

  po::options_description simulations("Available simulations",help_text_length);
  simulations.add_options()
    ("scan", po::value<bool>(&coherence_scan)->default_value(false)->implicit_value(true),
     "perform NV coherence scan for effective larmor frequencies")
    ("signal",
     po::value<bool>(&coherence_signal)->default_value(false)->implicit_value(true),
     "measure NV coherence signal as a function of the coupling constant"
     " at the frequency of a target nucleus")
    ("angular_signal",
     po::value<bool>(&angular_coherence_signal)->default_value(false)->implicit_value(true),
     "measure NV coherence signal as a function of the coupling constant"
     " and angle of a control field at the frequency of a target nucleus")
    ("rotate", po::value<bool>(&rotation)->default_value(false)->implicit_value(true),
     "rotate an individual nucleus")
    ("couple", po::value<bool>(&coupling)->default_value(false)->implicit_value(true),
     "couple an individual nucleus to NV center")
    ("iswap", po::value<bool>(&iswap)->default_value(false)->implicit_value(true),
     "compute iSWAP fidelities")
    ("swap", po::value<bool>(&swap)->default_value(false)->implicit_value(true),
     "compute SWAP fidelities")
    ("identity", po::value<bool>(&identity)->default_value(false)->implicit_value(true),
     "compute fidelity of an identity operation on a single spin")
    ("swap_nvst", po::value<bool>(&swap_nvst)->default_value(false)->implicit_value(true),
     "compute SWAP_NVST fidelity")
    ("larmor_identity",
     po::value<bool>(&larmor_identity)->default_value(false)->implicit_value(true),
     "compute fidelity of an identity operation on a larmor qubit")
    ("initialize",
     po::value<bool>(&initialize)->default_value(false)->implicit_value(true),
     "compute fidelity of a deterministic initialization of a thermalized nucleus"
     " into |u> or |d>")
    ("initialize_x",
     po::value<bool>(&initialize_x)->default_value(false)->implicit_value(true),
     "compute fidelity of a probabalistic initialization of a thermalized nucleus"
     " into |u> +/- |d>")
    ("initialize_larmor",
     po::value<bool>(&initialize_larmor)->default_value(false)->implicit_value(true),
     "compute fidelity of a probabalistic initialization of a larmor pair"
     " from |dd> into |ud> +/- |du>")
    ("test" ,po::value<bool>(&testing_mode)->default_value(false)->implicit_value(true),
     "enable testing mode")
    ;

  double c13_abundance;
  double c13_factor;
  uint max_cluster_size;
  double hyperfine_cutoff;
  double hyperfine_cutoff_in_kHz;
  int ms;
  uint k_DD_int;
  axy_harmonic k_DD;
  double static_Bz_in_gauss;
  double identity_time;
  double scale_factor;
  double integration_factor;
  bool no_nn;

  po::options_description simulation_options("Simulation options",help_text_length);
  simulation_options.add_options()
    ("c13_factor", po::value<double>(&c13_factor)->default_value(1),
     "abundance of C-13 relative to its natural abundance")
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
    ("identity_time", po::value<double>(&identity_time)->default_value(1),
     "time of identity operation (s)")
    ("scale_factor", po::value<double>(&scale_factor)->default_value(10),
     "factor used to define different scales (i.e. if a << b, then a = b/scale_factor)")
    ("integration_factor", po::value<double>(&integration_factor)->default_value(10),
     "factor used to determine size of integration step size")
    ("no_nn" ,po::value<bool>(&no_nn)->default_value(false)->implicit_value(true),
     "turn off internuclear couplings")
    ;

  vector<uint> target_nuclei;
  bool target_pairs;
  double angle;
  double angle_over_pi;
  double target_polar;
  double target_pitch_over_pi;
  double target_azimuth;
  double target_azimuth_over_pi;
  double nv_polar;
  double nv_pitch_over_pi;
  double nv_azimuth;
  double nv_azimuth_over_pi;

  po::options_description addressing_options("Nuclear spin addressing options"
                                             " (all angles are in units of pi radians)",
                                             help_text_length);
  addressing_options.add_options()
    ("target", po::value<vector<uint>>(&target_nuclei)->multitoken(),
     "indices of nuclei to target")
    ("target_pairs",
     po::value<bool>(&target_pairs)->default_value(false)->implicit_value(true),
     "target only one nucleus in each larmor pair")
    ("angle", po::value<double>(&angle_over_pi)->default_value(1),
     "angle of operation to perform")
    ("target_pitch", po::value<double>(&target_pitch_over_pi)->default_value(0),
     "pitch (angle above x-y plane) of target rotation axis")
    ("target_azimuth", po::value<double>(&target_azimuth_over_pi)->default_value(0),
     "azimuthal angle of target rotation axis")
    ("nv_pitch", po::value<double>(&nv_pitch_over_pi)->default_value(0.5),
     "pitch (angle above x-y plane) of NV rotation axis")
    ("nv_azimuth", po::value<double>(&nv_azimuth_over_pi)->default_value(0),
     "azimuthal angle of NV rotation axis")
    ;

  uint coherence_bins;
  double scan_time;
  double scan_time_in_ms;
  double f_DD;
  uint angular_resolution;
  double signal_gB_factor;
  double max_f_factor;

  po::options_description scan_options("Coherence measurement options",help_text_length);
  scan_options.add_options()
    ("bins", po::value<uint>(&coherence_bins)->default_value(500),
     "number of bins in coherence measurements")
    ("scan_time", po::value<double>(&scan_time_in_ms)->default_value(1),
     "time for each measurement in the coherence scan (microseconds)")
    ("f_DD", po::value<double>(&f_DD)->default_value(0.06,"0.06"),
     "magnitude of fourier component used in coherence scanning")
    ("angular_resolution", po::value<uint>(&angular_resolution)->default_value(20),
     "angular resolution for magnetic field direction during signal measurement"
     " (1/[pi radians])")
    ("signal_gB_factor", po::value<double>(&signal_gB_factor)->default_value(100),
     "sets signal field strength to static_gBz/signal_gB_factor")
    ("max_f_factor", po::value<double>(&max_f_factor)->default_value(0.5),
     "factor scaling maximum value of f_DD in coherence signal measurement")
    ;

  string lattice_file;

  po::options_description file_io("File IO",help_text_length);
  file_io.add_options()
    ("lattice_file", po::value<string>(&lattice_file),
     "input file defining system configuration")
    ;

  bool print_lattice;
  bool print_pairs;
  bool target_info;

  po::options_description print_options("Available printing options",help_text_length);
  print_options.add_options()
    ("print_lattice",
     po::value<bool>(&print_lattice)->default_value(false)->implicit_value(true),
     "print C-13 lattice to a file")
    ("print_pairs",
     po::value<bool>(&print_pairs)->default_value(false)->implicit_value(true),
     "print larmor pairs")
    ("target_info",
     po::value<bool>(&target_info)->default_value(false)->implicit_value(true),
     "print information about target nuclei")
    ;

  po::options_description all("Allowed options");
  all.add(general);
  all.add(simulations);
  all.add(simulation_options);
  all.add(addressing_options);
  all.add(scan_options);
  all.add(file_io);
  all.add(print_options);

  // collect inputs
  po::variables_map inputs;
  po::store(parse_command_line(arg_num, arg_vec, all), inputs);
  po::notify(inputs);


  // if requested, print help text
  if (inputs.count("help")) {
    cout << all;
    return 0;
  }


  // -------------------------------------------------------------------------------------
  // Run a sanity check on inputs
  // -------------------------------------------------------------------------------------

  // determine whether certain options were used
  bool using_input_lattice = inputs.count("lattice_file");
  bool set_hyperfine_cutoff = !inputs["hyperfine_cutoff"].defaulted();
  bool set_c13_factor = !inputs["c13_factor"].defaulted();
  bool set_target_nuclei = inputs.count("target");

  // make sure we are either printing something, or performing a simulation
  if (!testing_mode) {
    const bool printing = print_lattice || print_pairs || target_info;
    if (!printing) {
      if (int(coherence_scan)
          + int(coherence_signal)
          + int(angular_coherence_signal)
          + int(rotation)
          + int(coupling)
          + int(iswap)
          + int(swap)
          + int(identity)
          + int(swap_nvst)
          + int(larmor_identity)
          + int(initialize)
          + int(initialize_x)
          + int(initialize_larmor)
          != 1) {
        cout << "Please choose one simulation to perform\n";
        return -1;
      }
    }
  }

  // check lattice options
  assert(!(print_lattice && using_input_lattice));
  assert(!(using_input_lattice && !fs::exists(lattice_file)));
  assert(!(using_input_lattice && set_hyperfine_cutoff));
  assert(!(using_input_lattice && set_c13_factor));

  // check targeting options
  assert(!(set_target_nuclei && target_pairs));

  // verify validity of other values
  assert(hyperfine_cutoff_in_kHz > 0);
  assert(c13_factor >= 0);
  assert(max_cluster_size > 0);
  assert(ms == 1 || ms == -1);
  assert((k_DD_int == 1) || (k_DD_int == 3));
  assert(scale_factor > 1);
  assert(integration_factor > 1);

  if (coherence_scan) {
    assert(coherence_bins > 0);
    assert(scan_time_in_ms > 0);
  }
  assert(angular_resolution > 2);
  assert(signal_gB_factor > 1);
  assert(max_f_factor <= 1);

  // set some variables based on iputs
  c13_abundance = c13_factor*c13_natural_abundance;
  hyperfine_cutoff = hyperfine_cutoff_in_kHz*kHz;
  k_DD = (k_DD_int == 1 ? first : third);
  scan_time = scan_time_in_ms*1e-3;

  angle = angle_over_pi*pi;
  target_polar = pi/2 - target_pitch_over_pi*pi;
  target_azimuth = target_azimuth_over_pi*pi;
  nv_polar = pi/2 - nv_pitch_over_pi*pi;
  nv_azimuth = nv_azimuth_over_pi*pi;

  uniform_real_distribution<double> rnd(0.0,1.0); // uniform distribution on [0,1)
  mt19937_64 generator(seed); // use and seed the 64-bit Mersenne Twister 19937 generator

  // -------------------------------------------------------------------------------------
  // Construct lattice of nuclei
  // -------------------------------------------------------------------------------------

  vector<Vector3d> nuclei;

  if (!using_input_lattice) { // place nuclei at lattice sites

    // number of fcc cells to simulate out from the origin
    const int cell_radius
      = round(pow(abs(g_e*g_C13)/(4*pi*a0*a0*a0*hyperfine_cutoff),1.0/3));

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

    if (print_lattice) {
      for (uint i = 0; i < nuclei.size(); i++) {
        cout << nuclei.at(i)(0) << " "
             << nuclei.at(i)(1) << " "
             << nuclei.at(i)(2) << endl;
      }
    }

  } else { // if using_input_lattice, read in the lattice

    string line;
    ifstream lattice(lattice_file);

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

  // -------------------------------------------------------------------------------------
  // Characterize nuclei and targets
  // -------------------------------------------------------------------------------------

  // identify nuclei which cannot be addressed
  vector<uint> addressable_nuclei;
  vector<uint> unaddressable_nuclei;
  for (uint n = 0; n < nuclei.size(); n++) {
    if (can_address(nuclei,n)) {
      addressable_nuclei.push_back(n);
    } else {
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

  // identify larmor pairs
  vector<uint> larmor_nuclei;
  vector<vector<uint>> larmor_pairs;
  for (uint i = 0; i < addressable_nuclei.size(); i++) {
    const uint n_i = addressable_nuclei.at(i);
    for (uint j = i+1; j < addressable_nuclei.size(); j++) {
      const uint n_j = addressable_nuclei.at(j);
      if (is_larmor_pair(nuclei,n_i,n_j)) {
        larmor_nuclei.push_back(n_i);
        larmor_nuclei.push_back(n_j);
        larmor_pairs.push_back({n_i,n_j});
        continue;
      }
    }
  }

  // if we wish to print pairs, do so
  if (print_pairs) {
    if (larmor_pairs.size() > 0) {
      cout << "Larmor pairs:\n";
      for (vector<uint> pair: larmor_pairs) {
        cout << " " << pair.at(0) << " " << pair.at(1) << endl;
      }
      cout << endl;
    } else {
      cout << "No larmor pairs found\n";
    }
  }

  // determine which nuclei to target
  if (!set_target_nuclei) {
    if (target_pairs) target_nuclei = larmor_nuclei;
    else target_nuclei = addressable_nuclei;

  } else { // if (set_target_nuclei)
    vector<uint> unaddressable_targets;
    for (uint t_i = 0; t_i < target_nuclei.size(); t_i++) {
      if (in_vector(target_nuclei.at(t_i),unaddressable_nuclei)) {
        unaddressable_targets.push_back(target_nuclei.at(t_i));
        target_nuclei.erase(target_nuclei.begin()+t_i);
        t_i--;
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
  if (target_nuclei.size() == 0) {
    cout << "There are no target nuclei." << endl;
    return 0;
  }

  // identify targeted larmor pairs
  vector<vector<uint>> targeted_larmor_pairs;
  for (uint i = 0; i < target_nuclei.size(); i++) {
    const uint n_i = target_nuclei.at(i);
    for (uint j = i+1; j < target_nuclei.size(); j++) {
      const uint n_j = target_nuclei.at(j);
      if (is_larmor_pair(nuclei,n_i,n_j)) {
        targeted_larmor_pairs.push_back({n_i,n_j});
        continue;
      }
    }
  }

  if (swap_nvst || larmor_identity) {
    if (targeted_larmor_pairs.size() == 0) {
      if (larmor_pairs.size() == 0) {
        cout << "There are no larmor pairs in this system\n";
      } else {
        cout << "No larmor pairs are targeted\n";
      }
      return -1;
    }
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
    if (!testing_mode && (max_cluster_size < min_cluster_size_cap)) return -1;
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

  // initialize nv_system object
  const nv_system nv(nuclei, clusters, ms, g_C13*static_Bz_in_gauss*gauss, k_DD,
                     scale_factor, integration_factor, no_nn);

  // -------------------------------------------------------------------------------------
  // Print info about target nuclei
  // -------------------------------------------------------------------------------------

  if (target_info) {
    cout << "Target info:\n\n";
    for (uint n: target_nuclei) {
      cout << "index: " << n << endl
           << "position (nm): "
           << in_crystal_basis(nv.nuclei.at(n)).transpose() * a0/2 / nm << endl
           << "effective_larmor (kHz): " << effective_larmor(nv,n).norm() / kHz << endl
           << "hyperfine (kHz): " << hyperfine(nv,n).norm() / kHz << endl
           << "hyperfine_perp (kHz): " << hyperfine_perp(nv,n).norm() / kHz << endl
           << endl;
    }
    for (uint i = 0; i < target_nuclei.size(); i++) {
      const Vector3d pos_i = nv.nuclei.at(target_nuclei.at(i));
      for (uint j = i + 1; j < target_nuclei.size(); j++) {
        const Vector3d pos_j = nv.nuclei.at(target_nuclei.at(j));
        cout << "indices: " << i << " " << j << endl
             << "displacement (nm): "
             << in_crystal_basis(pos_j-pos_i).transpose() * a0/2 / nm << endl
             << "coupling (Hz): " << coupling_strength(pos_i,pos_j) / Hz << endl
             << endl;
      }
    }
  }


  // -------------------------------------------------------------------------------------
  // Coherence scan
  // -------------------------------------------------------------------------------------

  if (coherence_scan) {
    cout << endl
         << "Larmor and hyperfine frequency data:" << endl
         << "# format: w_larmor A_perp" << endl;

    double w_max = 0, w_min = DBL_MAX;
    for (uint i = 0; i < nv.nuclei.size(); i++) {
      const double A_perp = hyperfine_perp(nv,i).norm();
      const double w_larmor = effective_larmor(nv,i).norm();
      cout << i << " " << w_larmor << " " << A_perp << endl;

      if (w_larmor < w_min) w_min = w_larmor;
      if (w_larmor > w_max) w_max = w_larmor;
    }

    cout << endl
         << "Coherence scan results:" << endl
         << "# format: w_scan coherence" << endl;
    const double w_range = w_max - w_min;
    const double w_start = max(w_min - w_range/10, 0.);
    const double w_end = w_max + w_range/10;
    for (uint i = 0; i < coherence_bins; i++) {
      const double w_scan = w_start + (i+0.5)*(w_end-w_start)/coherence_bins;
      const double coherence = coherence_measurement(nv, w_scan, f_DD, scan_time);
      cout << w_scan << " " << coherence << endl;
    }

  }

  // -------------------------------------------------------------------------------------
  // Coherence signal
  // -------------------------------------------------------------------------------------

  if (coherence_signal || angular_coherence_signal) {
    cout << "Coherence signal results:" << endl
         << "# k_DD: " << nv.k_DD << endl;
    cout << "# format: f_DD coherence(s)" << endl;

    for (uint target: target_nuclei) {
      // only perform this measurement once for larmor pairs
      vector<uint> targets = {target};
      for (uint nn: target_nuclei) {
        if (nn == target) continue;
        if (is_larmor_pair(nv,nn,target)) targets.push_back(nn);
      }
      if (targets.size() == 2 && targets.at(0) > targets.at(1)) continue;

      cout << "# targets:";
      for (uint nn: targets) cout << " " << nn;
      cout << endl;
      if (angular_coherence_signal) {
        cout << "# azimuths/pi:";
        for (uint nn: targets) {
          const double angular_pos = atan2(dot(nv.nuclei.at(nn),yhat),
                                           dot(nv.nuclei.at(nn),xhat));
          cout << " " << angular_pos/pi;
        }
        cout << endl;
      }

      const double w_signal = effective_larmor(nv,target).norm();
      const double measurement_time = 4*pi / (axy_f_max(nv.k_DD) * max_f_factor
                                              * hyperfine_perp(nv,target).norm() / 4);
      for (uint ii = 0; ii < coherence_bins; ii++) {
        const double f_DD = (ii+0.5)/coherence_bins * axy_f_max(nv.k_DD) * max_f_factor;

        if (!angular_coherence_signal) {
          const double coherence =
            coherence_measurement(nv, w_signal, f_DD, measurement_time);
          cout << f_DD << " " << coherence << endl;

        } else {
          cout << f_DD;
          for (uint jj = 0; jj < angular_resolution; jj++) {
            const double phi_dec = -(jj+0.5)/angular_resolution * pi;
            const control_fields controls(nv.static_gBz/signal_gB_factor*xhat,
                                          w_signal, phi_dec);
            const double coherence =
              coherence_measurement(nv, w_signal, f_DD, measurement_time, controls);
            cout << " " << coherence;
          }
          cout << endl;
        }
      }
    }
  }

  // -------------------------------------------------------------------------------------
  // Individual addressing -- control
  // -------------------------------------------------------------------------------------

  if (rotation) {
    cout << "target fidelity time pulses\n";
    for (uint target: target_nuclei) {
      const Vector3d target_axis = axis(target_polar,target_azimuth);
      vector<protocol> P(2);
      for (bool exact : {true,false}) {
        P.at(exact) = rotate_target(nv, target, angle, target_axis, exact);
      }
      const uint subsystem_target = get_index_in_subsystem(nv, target);
      cout << target << " "
           << protocol_fidelity(P, {subsystem_target}) << " "
           << P.at(false).time << " "
           << P.at(false).pulses << endl;
    }
  }

  // -------------------------------------------------------------------------------------
  // Individual addressing -- NV coupling
  // -------------------------------------------------------------------------------------

  if (coupling) {
    cout << "target fidelity time pulses\n";
    const Vector3d nv_axis = axis(nv_polar, nv_azimuth);
    for (uint target: target_nuclei) {
      const Vector3d target_axis = axis(target_polar, target_azimuth);
      vector<protocol> P(2);
      for (bool exact : {true,false}) {
        P.at(exact) = couple_target(nv, target, angle, nv_axis, target_axis, exact);
      }
      const uint subsystem_target = get_index_in_subsystem(nv, target);
      cout << target << " "
           << protocol_fidelity(P, {0, subsystem_target}) << " "
           << P.at(false).time << " "
           << P.at(false).pulses << endl;
    }
  }

  // -------------------------------------------------------------------------------------
  // NV/nucleus iSWAP fidelity
  // -------------------------------------------------------------------------------------

  if (iswap) {
    cout << "target fidelity time pulses\n";
    for (uint target: target_nuclei) {
      vector<protocol> P(2);
      for (bool exact : {true,false}) {
        P.at(exact) = iSWAP(nv, target, exact);
      }
      const uint subsystem_target = get_index_in_subsystem(nv, target);
      cout << target << " "
           << protocol_fidelity(P, {0, subsystem_target}) << " "
           << P.at(false).time << " "
           << P.at(false).pulses << endl;
    }
  }

  // -------------------------------------------------------------------------------------
  // SWAP fidelity: NV electron spin and single nuclear spin
  // -------------------------------------------------------------------------------------

  if (swap) {
    cout << "target fidelity time pulses\n";
    for (uint target: target_nuclei) {
      vector<protocol> P(2);
      for (bool exact : {true,false}) {
        P.at(exact) = SWAP(nv, target, exact);
      }
      const uint subsystem_target = get_index_in_subsystem(nv, target);
      cout << target << " "
           << protocol_fidelity(P, {0, subsystem_target}) << " "
           << P.at(false).time << " "
           << P.at(false).pulses << endl;
    }
  }

  // -------------------------------------------------------------------------------------
  // Fidelity of identity operation on a single spin
  // -------------------------------------------------------------------------------------

  if (identity) {
    cout << "target fidelity\n";
    for (uint target: target_nuclei) {
      vector<protocol> P(2);
      for (bool exact : {true,false}) {
        P.at(exact) = target_identity(nv, target, identity_time, exact);
      }
      const uint cluster_target = get_index_in_cluster(nv, target);
      cout << target << " "
           << protocol_fidelity(P, {cluster_target}) << endl;
    }
  }

  // -------------------------------------------------------------------------------------
  // SWAP fidelity: NV electron spin and singlet-triplet subspace of two nuclear spins
  // -------------------------------------------------------------------------------------

  if (swap_nvst) {
    cout << "idx1 idx2 fidelity time pulses\n";
    for (vector<uint> idxs: targeted_larmor_pairs) {
      const uint idx1 = idxs.at(0);
      const uint idx2 = idxs.at(1);
      vector<protocol> P(2);
      for (bool exact : {true,false}) {
        P.at(exact) = SWAP_NVST(nv, idx1, idx2, exact);
      }
      const uint ss_idx1 = get_index_in_subsystem(nv,idx1);
      const uint ss_idx2 = get_index_in_subsystem(nv,idx2);
      cout << idx1 << " " << idx2 << " "
           << protocol_fidelity(P, {0, ss_idx1, ss_idx2}) << " "
           << P.at(false).time << " "
           << P.at(false).pulses << endl;
    }
  }

  // -------------------------------------------------------------------------------------
  // Fidelity of identity operation on a larmor qubit
  // -------------------------------------------------------------------------------------

  if (larmor_identity) {
    cout << "idx1 idx2 fidelity\n";
    const bool targeting_pair = true;
    for (vector<uint> idxs: targeted_larmor_pairs) {
      const uint idx1 = idxs.at(0);
      const uint idx2 = idxs.at(1);
      vector<protocol> P(2);
      for (bool exact : {true,false}) {
        P.at(exact) = target_identity(nv, idx1, identity_time, exact, targeting_pair);
      }
      const uint cluster_idx1 = get_index_in_cluster(nv,idx1);
      const uint cluster_idx2 = get_index_in_cluster(nv,idx2);
      cout << idx1 << " " << idx2 << " "
           << protocol_fidelity(P, {cluster_idx1, cluster_idx2}) << endl;
    }
  }

  // -------------------------------------------------------------------------------------
  // Fidelity of deterministically initializing a thermalized nucleus into |u> or |d>
  // -------------------------------------------------------------------------------------

  if (initialize) {
    cout << "target fidelity time pulses\n";
    for (uint target: target_nuclei) {
      vector<protocol> P(2);
      for (bool exact : {true,false}) {
        P.at(exact) = initialize_spin(nv, target, exact);
      }
      const uint subsystem_target = get_index_in_subsystem(nv, target);
      cout << target << " "
           << protocol_fidelity(P, {0, subsystem_target}) << " "
           << P.at(false).time << " "
           << P.at(false).pulses << endl;
    }
  }

  // -------------------------------------------------------------------------------------
  // Fidelity of probabalistically initializing a thermalized nucleus into |u> +/- |d>
  // -------------------------------------------------------------------------------------

  if (initialize_x) {
    cout << "target fidelity time pulses\n";
    for (uint target: target_nuclei) {
      vector<protocol> P(2);
      for (bool exact : {true,false}) {
        P.at(exact) = initialize_spin_X(nv, target, exact);
      }
      const uint subsystem_target = get_index_in_subsystem(nv, target);
      cout << target << " "
           << protocol_fidelity(P, {0, subsystem_target}) << " "
           << P.at(false).time << " "
           << P.at(false).pulses << endl;
    }
  }

  // -------------------------------------------------------------------------------------
  // Fidelity of probabalistically initializing a larmor pair from |dd> into |ud> +/- |du>
  // -------------------------------------------------------------------------------------

  if (initialize_larmor) {
    cout << "idx1 idx2 fidelity time pulses\n";
    for (vector<uint> idxs: targeted_larmor_pairs) {
      const uint idx1 = idxs.at(0);
      const uint idx2 = idxs.at(1);
      vector<protocol> P(2);
      for (bool exact : {true,false}) {
        P.at(exact) = initialize_larmor_qubit(nv, idx1, idx2, exact);
      }
      const uint ss_idx1 = get_index_in_subsystem(nv,idx1);
      const uint ss_idx2 = get_index_in_subsystem(nv,idx2);
      cout << idx1 << " " << idx2 << " "
           << protocol_fidelity(P, {0, ss_idx1, ss_idx2}) << " "
           << P.at(false).time << " "
           << P.at(false).pulses << endl;
    }
  }

}
