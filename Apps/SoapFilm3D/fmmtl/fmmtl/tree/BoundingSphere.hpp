#include <algorithm>
#include <cmath>
#include <iostream>
#include "fmmtl/numeric/norm.hpp"
#include "fmmtl/meta/dimension.hpp"

namespace fmmtl {

template <typename POINT>
class BoundingSphere {
 public:
	typedef POINT point_type;

	// Constructor
	BoundingSphere(const point_type& center, double radius_sq)
      : center_(center), radius_sq_(radius_sq) {}

	const point_type& center() const {
		return center_;
	}

	double radius() const {
		return std::sqrt(radius_sq_);
	}

	double radius_sq() const {
		return radius_sq_;
	}

	bool constains(const point_type& p) const {
		return norm_2_sq(p - center()) < radius_sq();
	}

 private:
	point_type center_;
	double radius_sq_;
};  // end class BoundingSphere


template <typename P>
inline std::ostream& operator<<(std::ostream& s, const BoundingSphere<P>& b) {
	const unsigned dim = b.center().size();
	s << "[Center: " << b.center()[0];
	for (unsigned i = 1; i != dim; ++i)
		s << ", " << b.center()[i];
	s << " ; Radius: " << b.radius();
	return s << "]";
}

} //end namespace fmmtl
