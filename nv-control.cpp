#include <iostream> // for standard output
using namespace std;

#include <eigen3/Eigen/Dense> // linear algebra library
using namespace Eigen;

#include "constants.h"
#include "qp-math.h"
#include "nv-math.h"
#include "nv-control.h"

//--------------------------------------------------------------------------------------------
// Coordinate systems
//--------------------------------------------------------------------------------------------

// return "natural" basis of a nucleus
vector<Vector3d> natural_basis(const nv_system& nv, const uint index){
  // larmor frequency of and hyperfine field at target nucleus
  const Vector3d larmor_eff = effective_larmor(nv,index);
  const Vector3d hyperfine_perp = A_perp(nv,index);
  // "natural" basis vectors for target nucleus
  const Vector3d target_zhat = hat(larmor_eff);
  const Vector3d target_xhat = hat(hyperfine_perp);
  const Vector3d target_yhat = target_zhat.cross(target_xhat);
  return {target_xhat, target_yhat, target_zhat};
}

//--------------------------------------------------------------------------------------------
// General control methods
//--------------------------------------------------------------------------------------------

// propagator U = exp(-i * rotation_angle * sigma_{axis}^{index})
MatrixXcd U_ctl(const nv_system& nv, const uint index, const double target_axis_azimuth,
                const double rotation_angle, const bool exact, const bool adjust_AXY){
  // identify cluster of target nucleus
  const uint cluster = get_cluster_containing_index(nv,index);
  const uint index_in_cluster = get_index_in_cluster(index,nv.clusters.at(cluster));
  const uint spins = nv.clusters.at(cluster).size()+1;

  // target axis of rotation
  const Vector3d axis_ctl = natural_axis(nv, index, target_axis_azimuth);

  if(exact){
    // return exact propagator
    const MatrixXcd G = exp(-j * rotation_angle * dot(s_vec,axis_ctl));
    return act(G, {index_in_cluster+1}, spins);
  }

  for(uint i = 0; i < nv.clusters.at(cluster).size(); i++){
    if(nv.clusters.at(cluster).at(i) == index) continue;
    if(is_larmor_pair(nv, index, nv.clusters.at(cluster).at(i))){
      cout << "Cannot address nuclei with larmor pairs: "
           << index << ", " << nv.clusters.at(cluster).at(i) << endl;
      return MatrixXcd::Identity(pow(2,spins),pow(2,spins));
    }
  }

  // larmor frequency of target nucleus
  const double w_larmor = effective_larmor(nv,index).norm();
  const double t_larmor = 2*pi/w_larmor;

  // AXY protocol parameters
  const double sA = nv.scale_factor * A(nv,index).norm();
  const double w_DD = [&]() -> double {
    const double w_DD_large = (w_larmor+sA)/3.;
    if(w_larmor < sA){
      return w_DD_large;
    } else{
      const uint k_m = 2*int(0.5 * (w_larmor/sA-1) );
      const double w_DD_small = (w_larmor-sA)/k_m;
      if(w_DD_small > sA && isfinite(w_DD_small)){
        return w_DD_small;
      } else{
        return w_DD_large;
      }
    }
  }();

  const double t_DD = 2*pi/w_DD;
  const axy_harmonic k_DD = abs(w_DD - w_larmor) < abs(3*w_DD - w_larmor) ? first : third;
  const double f_DD = 0;

  // control field frequency = effective larmor frequency of target nucleus
  const double w_ctl = w_larmor; // control field frequency
  const double dw_min = larmor_resolution(nv,index);
  double g_B_ctl = dw_min/nv.scale_factor; // ctl field strength * gyromangnetic ratio

  const double control_period = 2*pi*4/g_B_ctl;
  double control_time = -4*rotation_angle/g_B_ctl; // control operation time
  while(control_time >= control_period) control_time -= control_period;
  while(control_time < 0) control_time += control_period;
  if(control_time > control_period/2){
    g_B_ctl *= -1;
    control_time = control_period-control_time;
  }

  const double B_ctl = g_B_ctl/nv.nuclei.at(index).g; // control field strength
  const control_fields controls(B_ctl*axis_ctl, w_ctl, 0.); // control field object

  MatrixXcd U_ctl = MatrixXcd::Identity(pow(2,spins),pow(2,spins));
  if(!adjust_AXY){
    U_ctl = simulate_propagator(nv, cluster, w_DD, f_DD, k_DD, control_time, controls);
  } else{ // if(adjust_AXY)
    assert(w_DD != w_larmor);
    if(w_DD < w_larmor){
      const uint freq_ratio = 2*round(0.5*w_larmor/w_DD);
      const double w_DD_adjusted = w_larmor/freq_ratio;
      const double t_DD_adjusted = 2*pi/w_DD_adjusted;
      const uint cycles = int(control_time/t_DD_adjusted);

      const double leading_time = control_time - cycles*t_DD_adjusted;
      const double trailing_time = t_DD_adjusted - leading_time;

      const MatrixXcd U_leading = simulate_propagator(nv, cluster, w_DD_adjusted, f_DD, k_DD,
                                                      leading_time, controls);
      const MatrixXcd U_trailing = simulate_propagator(nv, cluster, w_DD_adjusted, f_DD, k_DD,
                                                       trailing_time, controls, leading_time);

      U_ctl = U_leading * pow(U_trailing*U_leading,cycles);
    } else{ // if(w_DD > w_larmor)
      const uint freq_ratio = round(w_DD/w_larmor);
      const double w_DD_adjusted = w_larmor*freq_ratio;
      const double t_DD_adjusted = 2*pi/w_DD_adjusted;
      const uint cycles = int(control_time/t_larmor);

      const double leading_time = control_time - cycles*t_larmor;
      const double trailing_time = t_larmor - leading_time;

      const MatrixXcd U_leading = simulate_propagator(nv, cluster, w_DD_adjusted, f_DD, k_DD,
                                                      leading_time, controls);
      const MatrixXcd U_trailing = simulate_propagator(nv, cluster, w_DD_adjusted, f_DD, k_DD,
                                                       trailing_time, controls, leading_time);

      U_ctl = U_leading * pow(U_trailing*U_leading,cycles);
    }
  }

  const double flush_time = ceil(control_time/t_larmor)*t_larmor - control_time;
  const MatrixXcd U_flush =
    simulate_propagator(nv, cluster, w_DD, f_DD, k_DD, flush_time, control_time);

  return U_flush * U_ctl;
}

// propagator U = exp(-i * rotation_angle * sigma_{n_1}^{NV}*sigma_{n_2}^{index})
MatrixXcd U_int(const nv_system& nv, const uint index, const axy_harmonic k_DD,
                const double nv_axis_azimuth, const double nv_axis_polar,
                const double target_axis_azimuth, const double rotation_angle,
                const bool exact){
  // identify cluster of target nucleus
  const uint cluster = get_cluster_containing_index(nv,index);
  const uint index_in_cluster = get_index_in_cluster(index,nv.clusters.at(cluster));
  const uint spins = nv.clusters.at(cluster).size()+1;

  // NV spin axis
  const Vector3d nv_axis = axis(nv_axis_azimuth,nv_axis_polar);

  if(exact){
    // return exact propagator
    const Vector3d target_axis = natural_axis(nv,index,target_axis_azimuth);
    const MatrixXcd G = exp(-j * rotation_angle *
                            tp(dot(s_vec,nv_axis), dot(s_vec,target_axis)));
    return act(G, {0,index_in_cluster+1}, spins);
  }

  // verify that we can address this nucleus
  if(round(4*dot(nv.nuclei.at(index).pos,ao)) == 0.){
    cout << "Cannot address nuclei without hyperfine coupling perpendicular to the NV axis: "
         << index << endl;
    return MatrixXcd::Identity(pow(2,spins),pow(2,spins));
  }
  for(uint i = 0; i < nv.clusters.at(cluster).size(); i++){
    if(nv.clusters.at(cluster).at(i) == index) continue;
    if(is_larmor_pair(nv, index, nv.clusters.at(cluster).at(i))){
      cout << "Cannot address nuclei with larmor pairs: "
           << index << ", " << nv.clusters.at(cluster).at(i) << endl;
      return MatrixXcd::Identity(pow(2,spins),pow(2,spins));
    }
  }

  // larmor frequency of and perpendicular component of hyperfine field at target nucleus
  const double w_larmor = effective_larmor(nv,index).norm();
  const double t_larmor = 2*pi/w_larmor;
  const double dw_min = larmor_resolution(nv,index);
  const Vector3d hyperfine_perp = A_perp(nv,index);

  // AXY sequence parameters
  const double w_DD = w_larmor/k_DD; // AXY protocol angular frequency
  const double t_DD = 2*pi/w_DD; // AXY protocol period
  double f_DD = min(dw_min/(hyperfine_perp.norm()*nv.scale_factor), axy_f_max(k_DD));

  const double interaction_period = 2*pi/abs(f_DD*hyperfine_perp.norm()/8);
  double interaction_time = rotation_angle/(nv.ms*f_DD*hyperfine_perp.norm()/8);
  while(interaction_time >= interaction_period) interaction_time -= interaction_period;
  while(interaction_time < 0) interaction_time += interaction_period;
  if(interaction_time > interaction_period/2){
    f_DD *= -1;
    interaction_time = interaction_period-interaction_time;
  }

  const MatrixXcd U_int = simulate_propagator(nv, cluster, w_DD, f_DD, k_DD,
                                              interaction_time, -target_axis_azimuth/w_DD);

  const double flush_time = ceil(interaction_time/t_larmor)*t_larmor - interaction_time;
  const MatrixXcd U_flush = simulate_propagator(nv, cluster, w_DD, 0, k_DD, flush_time,
                                                interaction_time - target_axis_azimuth/w_DD);

  // rotate the NV spin between the desired axis and zhat
  const MatrixXcd nv_axis_to_zhat = act(rotate(zhat,nv_axis), {0}, spins);

  return U_flush * nv_axis_to_zhat.adjoint() * U_int * nv_axis_to_zhat;
}

//--------------------------------------------------------------------------------------------
// Specific operations
//--------------------------------------------------------------------------------------------

// iSWAP operation
MatrixXcd iSWAP(const nv_system& nv, const uint index, const axy_harmonic k_DD,
                const bool exact){
  const double iswap_angle = -pi/4;
  const double xy_polar = pi/2;
  const double xhat_azimuth = 0;
  const double yhat_azimuth = pi/2;
  return
    U_int(nv, index, k_DD, xhat_azimuth, xy_polar, xhat_azimuth, iswap_angle, exact) *
    U_int(nv, index, k_DD, yhat_azimuth, xy_polar, yhat_azimuth, iswap_angle, exact);
};

MatrixXcd SWAP_NVST(const nv_system& nv, const uint idx1, const uint idx2,
                    const axy_harmonic k_DD, const bool exact){
  // assert that both target nuclei are in the same cluster
  const vector<uint> cluster = nv.clusters.at(get_cluster_containing_index(nv,idx1));
  assert(in_vector(idx2,cluster));

  // define angles
  const double angle = pi/4;
  const double z_polar = 0;
  const double xy_polar = pi/2;
  const double xhat_azimuth = 0;
  const double yhat_azimuth = pi/2;

  // compute components of SWAP_NVST
  const MatrixXcd Rz_NV = R_NV(nv, zhat, 2*nv.ms*angle, cluster.size()+1);
  const MatrixXcd Rx_1 = U_ctl(nv, idx1, xhat_azimuth, angle, exact);
  const MatrixXcd Ry_1 = U_ctl(nv, idx1, yhat_azimuth, angle, exact);
  const MatrixXcd Rz_1 = Rx_1 * Ry_1 * Rx_1.adjoint();
  const MatrixXcd iSWAP_NV_1 =
    U_int(nv, idx1, k_DD, xhat_azimuth, xy_polar, xhat_azimuth, -angle, exact) *
    U_int(nv, idx1, k_DD, yhat_azimuth, xy_polar, yhat_azimuth, -angle, exact);
  const MatrixXcd cNOT_NV_1 =
    Rz_NV * Rx_1 * U_int(nv, idx1, k_DD, xhat_azimuth, z_polar, xhat_azimuth, -angle, exact);
  const MatrixXcd E_NV_2 =
    U_int(nv, idx2, k_DD, yhat_azimuth, xy_polar, xhat_azimuth, -angle, exact);

  // combine componenets into full SWAP_NVST operation
  const MatrixXcd M = E_NV_2.adjoint() * iSWAP_NV_1 * Rz_1.adjoint() * Rz_NV;
  return M.adjoint() * cNOT_NV_1 * M;
}
