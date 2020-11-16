//
//  Adapted from FOSSSim 2010
//
#ifndef __COLORMAP_H__
#define __COLORMAP_H__

#include <cassert>
#include <cmath>
#include <iostream>

class Colormap {
 public:
  // Colormaps supported by this class
  enum ColorScheme {
    MATLAB_JET,
    // MATLAB_HSV,
    MATLAB_HOT,
    MATLAB_COOL,
    MATLAB_SPRING,
    MATLAB_SUMMER,
    MATLAB_AUTUMN,
    MATLAB_WINTER,
    MATLAB_GRAY,
    MATLAB_BONE,
    MATLAB_COPPER,
    MATLAB_PINK
  };

  // Generate a color map with a given number of samples
  Colormap(const ColorScheme& color_scheme, const int& num_samples);

  // Deallocate any memory used by this object
  ~Colormap();

  // Set the number of samples in the color map
  void changeNumSamples(const int& num_samples);

  // Change the color map
  void changeColormap(const ColorScheme& color_scheme);

  // Change to the 'next' supported color map
  void incrementColormap();
  // Change to the 'previous' supported color map
  void decrementColormap();
  // TODO: refactor the above two methods so as not to repeat code.

  // Returns the number of samples of the current color map
  int getNumSamples() const;

  // Returns the current color map in use
  ColorScheme getColorScheme() const;

  // Given a scalar in [0,1.0], returns the corresponding color
  void getColorByDensity(const double& rho, double& r, double& g,
                         double& b) const;

  // Given an integer in [0,num_samples], returns the corresponding color
  void getColorByIndex(const int& i, double& r, double& g, double& b) const;

  // Given a scalar in [0,1.0], returns the red channel
  double getRByDensity(const double& rho) const;
  // Given a scalar in [0,1.0], returns the green channel
  double getGByDensity(const double& rho) const;
  // Given a scalar in [0,1.0], returns the blue channel
  double getBByDensity(const double& rho) const;

  // Given an integer in [0,num_samples], returns the red channel
  double getRByIndex(const int& i) const;
  // Given an integer in [0,num_samples], returns the green channel
  double getGByIndex(const int& i) const;
  // Given an integer in [0,num_samples], returns the blue channel
  double getBByIndex(const int& i) const;

  // Prints the values of the current color map
  friend std::ostream& operator<<(std::ostream& os, Colormap& rhs);

 private:
  void generateColormap(const ColorScheme& color_scheme);
  void generateMatlabJet();
  void generateMatlabHot();
  void generateMatlabCool();
  void generateMatlabSpring();
  void generateMatlabSummer();
  void generateMatlabAutumn();
  void generateMatlabWinter();
  void generateMatlabGray();
  void generateMatlabBone();
  void generateMatlabCopper();
  void generateMatlabPink();

  int m_num_samples;
  ColorScheme m_color_scheme;

  // TODO: Move these to eigen types, clean up generation
  //  code accordingly.
  double* m_r;
  double* m_g;
  double* m_b;
};

#endif
