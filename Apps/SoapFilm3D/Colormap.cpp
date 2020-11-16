#include "Colormap.h"

#include <algorithm>

Colormap::Colormap(const ColorScheme& color_scheme, const int& num_samples)
    : m_num_samples(num_samples),
      m_color_scheme(color_scheme),
      m_r(NULL),
      m_g(NULL),
      m_b(NULL) {
  assert(m_num_samples > 0);

  m_r = new double[m_num_samples];
  m_g = new double[m_num_samples];
  m_b = new double[m_num_samples];

  this->generateColormap(m_color_scheme);
}

Colormap::~Colormap() {
  if (m_r != NULL) {
    delete[] m_r;
    m_r = NULL;
  }
  if (m_g != NULL) {
    delete[] m_g;
    m_g = NULL;
  }
  if (m_b != NULL) {
    delete[] m_b;
    m_b = NULL;
  }
}

void Colormap::changeNumSamples(const int& num_samples) {
  assert(num_samples > 0);
  m_num_samples = num_samples;

  assert(m_r != NULL);
  delete[] m_r;
  assert(m_g != NULL);
  delete[] m_g;
  assert(m_b != NULL);
  delete[] m_b;

  m_r = new double[m_num_samples];
  m_g = new double[m_num_samples];
  m_b = new double[m_num_samples];

  this->generateColormap(m_color_scheme);
}

void Colormap::changeColormap(const ColorScheme& color_scheme) {
  m_color_scheme = color_scheme;
  this->generateColormap(m_color_scheme);
}

void Colormap::incrementColormap() {
  switch (m_color_scheme) {
    case MATLAB_JET:
      m_color_scheme = MATLAB_HOT;
      break;
    case MATLAB_HOT:
      m_color_scheme = MATLAB_COOL;
      break;
    case MATLAB_COOL:
      m_color_scheme = MATLAB_SPRING;
      break;
    case MATLAB_SPRING:
      m_color_scheme = MATLAB_SUMMER;
      break;
    case MATLAB_SUMMER:
      m_color_scheme = MATLAB_AUTUMN;
      break;
    case MATLAB_AUTUMN:
      m_color_scheme = MATLAB_WINTER;
      break;
    case MATLAB_WINTER:
      m_color_scheme = MATLAB_GRAY;
      break;
    case MATLAB_GRAY:
      m_color_scheme = MATLAB_BONE;
      break;
    case MATLAB_BONE:
      m_color_scheme = MATLAB_COPPER;
      break;
    case MATLAB_COPPER:
      m_color_scheme = MATLAB_PINK;
      break;
    case MATLAB_PINK:
      m_color_scheme = MATLAB_JET;
      break;
    default:
      break;
  }

  changeColormap(m_color_scheme);
}

void Colormap::decrementColormap() {
  switch (m_color_scheme) {
    case MATLAB_JET:
      m_color_scheme = MATLAB_PINK;
      break;
    case MATLAB_HOT:
      m_color_scheme = MATLAB_JET;
      break;
    case MATLAB_COOL:
      m_color_scheme = MATLAB_HOT;
      break;
    case MATLAB_SPRING:
      m_color_scheme = MATLAB_COOL;
      break;
    case MATLAB_SUMMER:
      m_color_scheme = MATLAB_SPRING;
      break;
    case MATLAB_AUTUMN:
      m_color_scheme = MATLAB_SUMMER;
      break;
    case MATLAB_WINTER:
      m_color_scheme = MATLAB_AUTUMN;
      break;
    case MATLAB_GRAY:
      m_color_scheme = MATLAB_WINTER;
      break;
    case MATLAB_BONE:
      m_color_scheme = MATLAB_GRAY;
      break;
    case MATLAB_COPPER:
      m_color_scheme = MATLAB_BONE;
      break;
    case MATLAB_PINK:
      m_color_scheme = MATLAB_COPPER;
      break;
    default:
      break;
  }

  changeColormap(m_color_scheme);
}

Colormap::ColorScheme Colormap::getColorScheme() const {
  return m_color_scheme;
}

int Colormap::getNumSamples() const { return m_num_samples; }

void Colormap::getColorByDensity(const double& rho, double& r, double& g,
                                 double& b) const {
  assert(rho >= 0.0);
  assert(rho <= 1.0);
  int idx = floor(((double)(m_num_samples - 1)) * rho + 0.5);
  getColorByIndex(idx, r, g, b);
}

void Colormap::getColorByIndex(const int& i, double& r, double& g,
                               double& b) const {
  assert(i >= 0);
  assert(i < m_num_samples);

  r = m_r[i];
  g = m_g[i];
  b = m_b[i];
}

double Colormap::getRByDensity(const double& rho) const {
  assert(rho >= 0.0);
  assert(rho <= 1.0);
  int idx = floor(((double)(m_num_samples - 1)) * rho + 0.5);
  return getRByIndex(idx);
}

double Colormap::getGByDensity(const double& rho) const {
  assert(rho >= 0.0);
  assert(rho <= 1.0);
  int idx = floor(((double)(m_num_samples - 1)) * rho + 0.5);
  return getGByIndex(idx);
}

double Colormap::getBByDensity(const double& rho) const {
  assert(rho >= 0.0);
  assert(rho <= 1.0);
  int idx = floor(((double)(m_num_samples - 1)) * rho + 0.5);
  return getBByIndex(idx);
}

double Colormap::getRByIndex(const int& i) const {
  assert(i >= 0);
  assert(i < m_num_samples);
  return m_r[i];
}

double Colormap::getGByIndex(const int& i) const {
  assert(i >= 0);
  assert(i < m_num_samples);
  return m_g[i];
}

double Colormap::getBByIndex(const int& i) const {
  assert(i >= 0);
  assert(i < m_num_samples);
  return m_b[i];
}

void Colormap::generateColormap(const ColorScheme& color_scheme) {
  switch (color_scheme) {
    case MATLAB_JET:
      generateMatlabJet();
      break;
    case MATLAB_HOT:
      generateMatlabHot();
      break;
    case MATLAB_COOL:
      generateMatlabCool();
      break;
    case MATLAB_SPRING:
      generateMatlabSpring();
      break;
    case MATLAB_SUMMER:
      generateMatlabSummer();
      break;
    case MATLAB_AUTUMN:
      generateMatlabAutumn();
      break;
    case MATLAB_WINTER:
      generateMatlabWinter();
      break;
    case MATLAB_GRAY:
      generateMatlabGray();
      break;
    case MATLAB_BONE:
      generateMatlabBone();
      break;
    case MATLAB_COPPER:
      generateMatlabCopper();
      break;
    case MATLAB_PINK:
      generateMatlabPink();
      break;
    default:
      break;
  }
}

void Colormap::generateMatlabJet() {
  // J = zeros(m,3);
  for (int i = 0; i < m_num_samples; ++i) {
    m_r[i] = m_g[i] = m_b[i] = 0.0;
  }

  // n = ceil(m/4);
  int n = ceil(((double)m_num_samples) / 4.0);

  // u = [(1:1:n)/n ones(1,n-1) (n:-1:1)/n]';
  double* u = new double[3 * n - 1];
  for (int i = 0; i < n; ++i) {
    u[i] = ((double)(i + 1)) / ((double)n);
  }
  for (int i = 1; i < n; ++i) {
    u[n - 1 + i] = 1.0;
  }
  for (int i = 0; i < n; ++i) {
    u[2 * n - 1 + i] = ((double)(n - i)) / ((double)n);
  }

  // g = ceil(n/2) - (mod(m,4)==1) + (1:length(u))';
  int* g = new int[3 * n - 1];
  int addedtoall = ((int)ceil(((double)n) / 2.0)) - ((m_num_samples % 4) == 1);
  for (int i = 0; i < 3 * n - 1; ++i) {
    g[i] = addedtoall + i + 1;
  }

  // r = g + n;
  int* r = new int[3 * n - 1];
  for (int i = 0; i < 3 * n - 1; ++i) {
    r[i] = g[i] + n;
  }

  // b = g - n;
  int* b = new int[3 * n - 1];
  for (int i = 0; i < 3 * n - 1; ++i) {
    b[i] = g[i] - n;
  }

  // g(g>m) = [];
  // r(r>m) = [];
  // b(b<1) = [];
  // J(r,1) = u(1:length(r));
  // J(g,2) = u(1:length(g));
  // J(b,3) = u(end-length(b)+1:end);
  int q = 0;
  int o = 0;
  int s = 0;
  int b_len = 0;
  for (int i = 0; i < 3 * n - 1; ++i) {
    if (!(b[i] < 1)) {
      ++b_len;
    }
  }
  for (int i = 0; i < 3 * n - 1; ++i) {
    if (!(r[i] > m_num_samples)) {
      m_r[r[o] - 1] = u[o];
      ++o;
    }
    if (!(g[i] > m_num_samples)) {
      m_g[g[q] - 1] = u[q];
      ++q;
    }
    if (!(b[i] < 1)) {
      m_b[b[(3 * n - 1) - b_len + s] - 1] = u[(3 * n - 1) - b_len + s];
      ++s;
    }
  }

  delete[] b;
  delete[] r;
  delete[] g;
  delete[] u;
}

void Colormap::generateMatlabHot() {
  // n = fix(3/8*m);
  int n = int((3.0 / 8.0) * m_num_samples);

  // r = [(1:n)'/n; ones(m-n,1)];
  for (int i = 0; i < n; ++i) {
    m_r[i] = ((double)i + 1) / ((double)n);
  }
  for (int i = n; i < m_num_samples; ++i) {
    m_r[i] = 1.0;
  }

  // g = [zeros(n,1); (1:n)'/n; ones(m-2*n,1)];
  for (int i = 0; i < n; ++i) {
    m_g[i] = 0.0;
  }
  for (int i = n; i < 2 * n; ++i) {
    m_g[i] = m_r[i - n];
  }
  for (int i = 2 * n; i < m_num_samples; ++i) {
    m_g[i] = 1.0;
  }

  // b = [zeros(2*n,1); (1:m-2*n)'/(m-2*n)];
  for (int i = 0; i < 2 * n; ++i) {
    m_b[i] = 0.0;
  }
  for (int i = 2 * n; i < m_num_samples; ++i) {
    m_b[i] = ((double)(i - 2 * n + 1)) / ((double)(m_num_samples - 2 * n));
  }
}

void Colormap::generateMatlabCool() {
  // r = (0:m-1)'/max(m-1,1);
  // c = [r 1-r ones(m,1)];
  double denom = std::max(m_num_samples - 1, 1);
  for (int i = 0; i < m_num_samples; ++i) {
    m_r[i] = ((double)i) / denom;
  }
  for (int i = 0; i < m_num_samples; ++i) {
    m_g[i] = 1.0 - m_r[i];
  }
  for (int i = 0; i < m_num_samples; ++i) {
    m_b[i] = 1.0;
  }
}

void Colormap::generateMatlabSpring() {
  // r = (0:m-1)'/max(m-1,1);
  // c = [ones(m,1) r 1-r];
  double denom = std::max(m_num_samples - 1, 1);
  for (int i = 0; i < m_num_samples; ++i) {
    m_g[i] = ((double)i) / denom;
  }
  for (int i = 0; i < m_num_samples; ++i) {
    m_b[i] = 1.0 - m_g[i];
  }
  for (int i = 0; i < m_num_samples; ++i) {
    m_r[i] = 1.0;
  }
}

void Colormap::generateMatlabSummer() {
  // r = (0:m-1)'/max(m-1,1);
  // c = [r .5+r/2 .4*ones(m,1)];
  double denom = std::max(m_num_samples - 1, 1);
  for (int i = 0; i < m_num_samples; ++i) {
    m_r[i] = ((double)i) / denom;
  }
  for (int i = 0; i < m_num_samples; ++i) {
    m_g[i] = 0.5 + 0.5 * m_r[i];
  }
  for (int i = 0; i < m_num_samples; ++i) {
    m_b[i] = 0.4;
  }
}

void Colormap::generateMatlabAutumn() {
  // r = (0:m-1)'/max(m-1,1);
  // c = [ones(m,1) r zeros(m,1)];
  double denom = std::max(m_num_samples - 1, 1);
  for (int i = 0; i < m_num_samples; ++i) {
    m_g[i] = ((double)i) / denom;
  }
  for (int i = 0; i < m_num_samples; ++i) {
    m_r[i] = 1.0;
  }
  for (int i = 0; i < m_num_samples; ++i) {
    m_b[i] = 0.0;
  }
}

void Colormap::generateMatlabWinter() {
  // r = (0:m-1)'/max(m-1,1);
  // c = [zeros(m,1) r .5+(1-r)/2];
  double denom = std::max(m_num_samples - 1, 1);
  for (int i = 0; i < m_num_samples; ++i) {
    m_g[i] = ((double)i) / denom;
  }
  for (int i = 0; i < m_num_samples; ++i) {
    m_b[i] = 0.5 + 0.5 * (1.0 - m_g[i]);
  }
  for (int i = 0; i < m_num_samples; ++i) {
    m_r[i] = 0.0;
  }
}

void Colormap::generateMatlabGray() {
  // g = (0:m-1)'/max(m-1,1);
  // g = [g g g];
  double denom = std::max(m_num_samples - 1, 1);
  for (int i = 0; i < m_num_samples; ++i) {
    m_r[i] = m_g[i] = m_b[i] = ((double)i) / denom;
  }
}

void Colormap::generateMatlabBone() {
  // b = (7*gray(m) + fliplr(hot(m)))/8;
  double* accm_r = new double[m_num_samples];
  double* accm_g = new double[m_num_samples];
  double* accm_b = new double[m_num_samples];
  this->generateMatlabHot();
  for (int i = 0; i < m_num_samples; ++i) {
    accm_b[i] = m_r[i];
    accm_g[i] = m_g[i];
    accm_r[i] = m_b[i];
  }
  this->generateMatlabGray();
  for (int i = 0; i < m_num_samples; ++i) {
    accm_r[i] += 7.0 * m_r[i];
    accm_g[i] += 7.0 * m_g[i];
    accm_b[i] += 7.0 * m_b[i];
  }
  for (int i = 0; i < m_num_samples; ++i) {
    m_r[i] = accm_r[i] / 8.0;
    m_g[i] = accm_g[i] / 8.0;
    m_b[i] = accm_b[i] / 8.0;
  }
  delete accm_r;
  delete accm_g;
  delete accm_b;
}

void Colormap::generateMatlabCopper() {
  // c = min(1,gray(m)*diag([1.2500 0.7812 0.4975]));
  this->generateMatlabGray();
  for (int i = 0; i < m_num_samples; ++i) {
    m_r[i] *= 1.25;
    m_g[i] *= 0.7812;
    m_b[i] *= 0.4975;
    m_r[i] = std::min(1.0, m_r[i]);
    m_g[i] = std::min(1.0, m_g[i]);
    m_b[i] = std::min(1.0, m_b[i]);
  }
}

void Colormap::generateMatlabPink() {
  // p = sqrt((2*gray(m) + hot(m))/3);
  double* accm_r = new double[m_num_samples];
  double* accm_g = new double[m_num_samples];
  double* accm_b = new double[m_num_samples];
  this->generateMatlabHot();
  for (int i = 0; i < m_num_samples; ++i) {
    accm_r[i] = m_r[i];
    accm_g[i] = m_g[i];
    accm_b[i] = m_b[i];
  }
  this->generateMatlabGray();
  for (int i = 0; i < m_num_samples; ++i) {
    accm_r[i] += 2.0 * m_r[i];
    accm_g[i] += 2.0 * m_g[i];
    accm_b[i] += 2.0 * m_b[i];
  }
  for (int i = 0; i < m_num_samples; ++i) {
    m_r[i] = accm_r[i] / 3.0;
    m_g[i] = accm_g[i] / 3.0;
    m_b[i] = accm_b[i] / 3.0;
    m_r[i] = sqrt(m_r[i]);
    m_g[i] = sqrt(m_g[i]);
    m_b[i] = sqrt(m_b[i]);
  }
  delete accm_r;
  delete accm_g;
  delete accm_b;
}

std::ostream& operator<<(std::ostream& os, Colormap& rhs) {
  // os <<  "Red\t\tGreen\t\tBlue" << std::endl;
  // for( int i = 0; i < rhs.m_num_samples; ++i )
  //{
  //	os << rhs.m_r[i] << "\t\t" << rhs.m_g[i] << "\t\t" << rhs.m_b[i] <<
  //std::endl;
  //}
  // os << dt.mo << '/' << dt.da << '/' << dt.yr;

  std::cout << "{";
  for (int i = 0; i < rhs.m_num_samples; ++i) {
    std::cout << rhs.m_r[i] << "," << rhs.m_g[i] << "," << rhs.m_b[i];
    if (i != rhs.m_num_samples - 1) std::cout << ",";
  }
  std::cout << "}";

  return os;
}

// function map = hsv(m)
//
// if nargin < 1, m = size(get(gcf,'colormap'),1); end
// h = (0:m-1)'/max(m,1);
// if isempty(h)
// map = [];
// else
// map = hsv2rgb([h ones(m,2)]);
// end
//
//
// function map = lines(n)
//
// if nargin<1, n = size(get(gcf,'Colormap'),1); end
//
// c = get(0,'defaultaxescolororder');
//
// map = c(rem(0:n-1,size(c,1))+1,:);
