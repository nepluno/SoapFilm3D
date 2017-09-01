#pragma once
/**
 * Storage class for all FMM options
 */

#include <iostream>
#include <fstream>
#include <cmath>

#include "fmmtl/config.hpp"

/** Class to define compile-time and run-time FMM options */
class FMMOptions {
  struct DefaultMAC {
    double theta_sq_;
    DefaultMAC(double theta) : theta_sq_(theta * theta) {}

    template <typename BOX>
    bool operator()(const BOX& b1, const BOX& b2) const {
      double r0_sq = norm_2_sq(b1.center() - b2.center());
      double r1_sq = b1.radius_sq();
      double r2_sq = b2.radius_sq();
      double r_sq = r1_sq + r2_sq + 2*std::sqrt(r1_sq*r2_sq);
      return theta_sq_ * r0_sq > r_sq;
    }
  };

 public:
  // Standard algorithm parameters
  unsigned ncrit;  // The maximum number of particles per box in the tree
  double theta;    // The aperture of the standard multipole acceptance criteria

  // DEBUGGING FLAGS
  bool print_tree;

  // OTHER
  //! Evaluation type
  enum EvalType {FMM, TREECODE};
  EvalType evaluator;

  FMMOptions()
      : ncrit(128),
        theta(0.5),
        print_tree(false),
        evaluator(FMM) {
  };

  // TODO: Generalize type/construction
  DefaultMAC MAC() const {
    return DefaultMAC(theta);
  }
};


#include <cstdio>
#include <cstdlib>
#include <cstring>

// Get the FMMOptions from command line arguments
FMMOptions get_options(int argc, char** argv) {
  FMMOptions opts = FMMOptions();

  // parse command line args
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i],"-theta") == 0) {
      opts.theta = (double) atof(argv[++i]);
    } else if (strcmp(argv[i],"-ncrit") == 0) {
      opts.ncrit = (unsigned) atoi(argv[++i]);
    } else if (strcmp(argv[i],"-printtree") == 0) {
      opts.print_tree = true;
    }
  }

  return opts;
}



/*
  #include <boost/program_options.hpp>
  namespace PO = boost::program_options;

  // XXX: This requires a linker... find a header-only solution.

  FMMOptions get_options(int argc, char** argv) {
  try {
  std::string config_file;

  // Declare options only allowed on command line
  PO::options_description cmdline("Command-line options");
  cmdline.add_options()
  ("help",
  "Produce help message")

  ("config,c",
  PO::value<std::string>(&config_file)->default_value("fmm_config.cfg"),
  "Name of configuration file")
  ;

  // Declare options allowed on command line and config file
  PO::options_description config("Configuration");
  config.add_options()
  ("ncrit",
  PO::value<unsigned>(),
  "Maximum points/box")

  ("theta",
  PO::value<double>(),
  "Multipole acceptance criteria parameter")

  ("verbose,v",
  PO::value<unsigned>()->implicit_value(0),
  "Verbosity level")
  ;

  // Declare options that are available, but hidden
  PO::options_description hidden("Hidden options");
  hidden.add_options()
  ;

  PO::options_description cmdline_options;
  cmdline_options.add(cmdline).add(config).add(hidden);

  PO::options_description config_options;
  config_options.add(config).add(hidden);

  PO::options_description visible("Allowed options");
  visible.add(cmdline).add(config);

  PO::variables_map vm;
  PO::store(PO::command_line_parser(argc, argv).
  options(cmdline_options).allow_unregistered().run(), vm);
  PO::notify(vm);

  if (vm.count("help")) {
  std::cout << "Usage: " << argv[0] << " [options]" << std::endl;
  std::cout << visible;
  exit(0);
  }

  std::ifstream ifs(config_file.c_str());
  if (ifs) {
  PO::store(PO::parse_config_file(ifs, config_options), vm);
  PO::notify(vm);
  } else {
  std::cout << "Can not open config file: " << config_file << std::endl;
  }

  FMMOptions opts;

  if (vm.count("ncrit")) {
  opts.set_max_per_box(vm["ncrit"].as<unsigned>());
  }

  if (vm.count("theta")) {
  opts.set_mac_theta(vm["theta"].as<double>());
  }

  return opts;

  } catch(std::exception& e) {
  std::cout << e.what() << "\n";
  exit(1);
  }
  }
*/
