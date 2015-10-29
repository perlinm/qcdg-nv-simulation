#include <sys/stat.h>

#include <iostream>
#include <fstream>
using namespace std;

#include <eigen3/Eigen/Dense> // linear algebra library
using namespace Eigen;

#include <boost/algorithm/string.hpp> // string manipulation library
#include <boost/program_options.hpp> // options parsing library
namespace po = boost::program_options;

#include "nv-math.h"

inline double rn(){ return (double)rand()/(double)RAND_MAX; }

int main(int arg_num, const char *arg_vec[]) {

  // -----------------------------------------------------------------------------------------
  // Parse and process input options
  // -----------------------------------------------------------------------------------------

  // define parameters
  int cell_radius;
  bool fast_decomposition;
  string input_lattice;
  string output_lattice;
  string output_dir;
  double c13_abundance = c13_natural_abundance;
  int seed;

  // define options
  po::options_description options("allowed options");
  options.add_options()
    ("help,h", "produce help message")
    ("fast_decomposition,f", po::value<bool>(&fast_decomposition)->default_value(true),
     "use fast, but less accurate algorithm for matrix decomposition")
    ("cell_radius", po::value<int>(&cell_radius)->default_value(2),
     "number of unit cells to simulate out from the NV center")
    ("c13_abundance", po::value<double>(&c13_abundance), "relative isotopic abundance of C13")
    ("input_lattice", po::value<string>(&input_lattice), "input file definint initial system")
    ("output_lattice",
     po::value<string>(&output_lattice)->default_value("lattice-cr[cell_radius]-s[seed].txt"),
     "name for output file defining initial system")
    ("output_dir", po::value<string>(&output_dir)->default_value("data"),
     "directory for storing data")
    ("seed", po::value<int>(&seed)->default_value(1), "seed for random number generator")
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
  bool set_output_lattice = !inputs["output_lattice"].defaulted();
  bool set_output_dir = !inputs["output_dir"].defaulted();

  // run a sanity check on inputs
  assert(cell_radius > 0);
  assert(c13_abundance >= 0 && c13_abundance <= 1);
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

  srand(seed); // initialize random number generator
  mkdir(output_dir.c_str(),0777); // create data directory

  // -----------------------------------------------------------------------------------------
  // Construct spin lattice
  // -----------------------------------------------------------------------------------------

  // positions of C-13 atoms with a0 = 1
  vector<Vector3d> c13_pos;

  if(!using_input_lattice){

    // set positions of C-13 atoms at lattice sites
    Vector3d cell_pos; // position of unit cell indexed by i,j,k
    for(int i = -cell_radius; i < cell_radius; i++){ // loop over all unit cell indices
      cell_pos(0) = i;
      for(int j = -cell_radius; j < cell_radius; j++){
        cell_pos(1) = j;
        for(int k = -cell_radius; k < cell_radius; k++){
          cell_pos(2) = k;
          for(int ls = 0; ls < cell_sites.size(); ls++){ // loop over sites in unit cell
            if(rn() <= c13_abundance){ // if we pass a check for C-13 isotopic abundance
              c13_pos.push_back(cell_pos+cell_sites.at(ls));
            }
          }
        }
      }
    }

    // remove any C-13 atoms at the NV lattice sites
    c13_pos.erase(remove(c13_pos.begin(), c13_pos.end(), n_pos), c13_pos.end());
    c13_pos.erase(remove(c13_pos.begin(), c13_pos.end(), e_pos), c13_pos.end());

    // write cell radius and C-13 positions to file
    ofstream output(output_lattice_path);
    output << "cell_radius: " << cell_radius << endl;
    for(int c = 0; c < c13_pos.size(); c++){
      output << c13_pos.at(c)[0] << " "
             << c13_pos.at(c)[1] << " "
             << c13_pos.at(c)[2] << endl;
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
      c13_pos.push_back((Vector3d() << x,y,z).finished());
    }
    input.close();

  }

  // assert that no C-13 atoms lie at the NV lattice sites
  if((find(c13_pos.begin(), c13_pos.end(), n_pos) != c13_pos.end()) ||
     (find(c13_pos.begin(), c13_pos.end(), e_pos) != c13_pos.end())){
    cout << "we have a C-13 atom at one of the NV lattice sites!" << endl;
    return 2;
  }

  // print summary of C-13 atom placement
  cout << "C-13 atoms occupy " << c13_pos.size() << " of "
       << pow(2*cell_radius,3)*cell_sites.size() << " lattice sites" << endl;

  int qbits = c13_pos.size()+1;

  // -----------------------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------------------



}

