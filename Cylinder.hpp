#ifndef CYLINDER_HPP
#define CYLINDER_HPP

#include <ostream>
#include <functional>
#include <cmath>
#include "Vector3.hpp"
#include "Shape.hpp"
#include "linear_algebra.hpp"

// Todo. Make sure cylinder is never larger than 1 cellsize or something.  
template<typename T_>
class Cylinder
{
public:
    typedef T_ value_type;
    typedef Vector3<T_> position_type;
    typedef T_ length_type;

public:
    // constructors
    Cylinder()
        : position_(), radius_(0), unit_z_(), half_length_(0) {}

    Cylinder(position_type const& position, length_type const& radius,
             position_type const& unit_z, length_type const& half_length )
        : position_(position), radius_(radius), unit_z_(unit_z),
          half_length_(half_length) {}

    
    bool operator==(const Cylinder& rhs) const
    {
        return position_ == rhs.position() && radius_ == rhs.radius() && unit_z_ == rhs.unit_z() && half_length_ == rhs.half_length();
    }

    bool operator!=(const Cylinder& rhs) const
    {
        return !operator==(rhs);
    }

    position_type const& position() const
    {
        return position_;
    }

    position_type& position()
    {
        return position_;
    }

    length_type const& radius() const
    {
        return radius_;
    }

    length_type& radius()
    {
        return radius_;
    }

    position_type const& unit_z() const
    {
        return unit_z_;
    }

    position_type& unit_z()
    {
        return unit_z_;
    }

    length_type const& half_length() const
    {
        return half_length_;
    }

    length_type& half_length()
    {
        return half_length_;
    }
    
    const int dof() const
    {   // degrees of freedom for particle movement
        return 1;
    }

    std::string show(int precision)
    {
        std::ostringstream strm;
        strm.precision(precision);
        strm << *this;
        return strm.str();
    }

/////// Member variables
private:
    position_type position_;    // centre.
    length_type radius_;
    position_type unit_z_;      // Z-unit_z. should be normalized.
    length_type half_length_;
};


//////// Inline functions
template<typename Tstrm_, typename T_>
inline std::basic_ostream<Tstrm_>& operator<<(std::basic_ostream<Tstrm_>& strm,
        const Cylinder<T_>& v)
{
    strm << "{" << v.position() <<  ", " << v.radius() << ", " << v.unit_z() << ", " << v.half_length() << "}";
    return strm;
}

template<typename T_>
inline std::pair<typename Cylinder<T_>::length_type,
                 typename Cylinder<T_>::length_type>
to_internal(Cylinder<T_> const& obj, typename Cylinder<T_>::position_type const& pos)
{
    // Return pos relative to position of cylinder. 
    typedef typename Cylinder<T_>::position_type    position_type;
    typedef typename Cylinder<T_>::length_type      length_type;

    const position_type pos_vector(subtract(pos, obj.position()));

    const length_type z(dot_product(pos_vector, obj.unit_z()));             // z can be < 0
    const length_type r(length(subtract(pos_vector, multiply(obj.unit_z(), z))));    // r is always >= 0

    return std::make_pair(r, z);
}

// The following projects 'pos' on the cylinder
// It returns a pair of which the first component is the projected position.
// The second is again a pair of which the first entry is the distance of 'pos'
// from the cylinder axis, the second a length l which indicates whether the
// projected position is in the cylinder (true if l negative; l is the negative
// distance of the projected position to the cylinder edge).
template<typename T_>
inline std::pair<typename Cylinder<T_>::position_type,
                 std::pair<typename Cylinder<T_>::length_type,
                           typename Cylinder<T_>::length_type> >
project_point(Cylinder<T_> const& obj, typename Cylinder<T_>::position_type const& pos)
{
    typedef typename Cylinder<T_>::length_type length_type;

    // The projection lies on the z-axis.
    std::pair<length_type, length_type> r_z(to_internal(obj, pos));

    return std::make_pair( add(obj.position(), multiply(obj.unit_z(), r_z.second)),
                           std::make_pair(r_z.first,
                                          subtract(abs(r_z.second), obj.half_length())) );
}

//Almost equal to projected point method, but for the substraction of the cylinder radius from the radial distance r.
//And projected point now lies on the surface, not on the central axis.
template<typename T_>
inline std::pair<typename Cylinder<T_>::position_type,
                 std::pair<typename Cylinder<T_>::length_type,
                           typename Cylinder<T_>::length_type> >
project_point_on_surface(Cylinder<T_> const& obj,
                typename Cylinder<T_>::position_type const& pos)
{
    typedef typename Cylinder<T_>::length_type length_type;
    typedef typename Cylinder<T_>::position_type position_type;

    // Here we do not call 'to_internal' for efficiency
    const position_type pos_vector(subtract(pos, obj.position()));

    const length_type   z ( dot_product(pos_vector, obj.unit_z()) );
    const position_type z_vector (multiply(obj.unit_z(), z));
    const position_type r_vector (subtract(pos_vector, z_vector));
    const length_type   r (length(r_vector));

    const position_type projected_point( add(obj.position(), z_vector) );

    return std::make_pair( add(projected_point, multiply( normalize(r_vector), obj.radius() )),
                           std::make_pair(subtract(r, obj.radius()),
                                          subtract(abs(z), obj.half_length())) );
}

template<typename T_>
inline typename Cylinder<T_>::length_type
distance(Cylinder<T_> const& obj,
                typename Cylinder<T_>::position_type const& pos)
{
    typedef typename Cylinder<T_>::length_type length_type;

    /* First compute the (r,z) components of pos in a coordinate system 
     * defined by the vectors unitR and unit_z, where unitR is
     * choosen such that unitR and unit_z define a plane in which
     * pos lies. */
    const std::pair<length_type, length_type> r_z(to_internal(obj, pos));

    /* Then compute distance to cylinder. */
    const length_type dz(std::fabs(r_z.second) - obj.half_length());
    const length_type dr(r_z.first - obj.radius());
    length_type distance;
    if (dz > 0)
    {
        // pos is (either) to the right or to the left of the cylinder.
        if (r_z.first > obj.radius())
        {
            // Compute distance to edge.
            distance = std::sqrt( dz * dz + dr * dr );
        }
        else
        {
            distance = dz;
        }
    }
    else
    {
        if (dr > obj.radius())
        {
            // pos is somewhere 'parallel' to the cylinder.
            distance = dr;
        }
        else
        {
            // Inside cylinder. 
            distance = std::max(dr, dz);
        }
    }
    return distance;
}

template<typename T_>
inline std::pair<typename Cylinder<T_>::position_type, bool>
deflect(Cylinder<T_> const& obj,
        typename Cylinder<T_>::position_type const& r0,
        typename Cylinder<T_>::position_type const& d  )
{
    // Displacements are not deflected on cylinders (yet),
    // but this function has to be defined for every shape to be used in structure.
    // For now it just returns original pos. + displacement. The changeflage = false.
    return std::make_pair( add(r0, d), false );
}
/*
template<typename T_>
inline typename Cylinder<T_>::position_type
deflect_back(Cylinder<T_> const& obj,
        typename Cylinder<T_>::position_type const& r,
        typename Cylinder<T_>::position_type const& u_z  )
{
    // Return the vector r without any changes
    return r;
}
*/
template<typename T, typename Trng>
inline typename Cylinder<T>::position_type
random_position(Cylinder<T> const& shape, Trng& rng)
{
    // -1 < rng() < 1. See for example CylindricalSurface.hpp.
    return add(shape.position(),
               multiply(shape.unit_z(), rng() * shape.half_length()));
}

template<typename T_>
inline Cylinder<T_> const& shape(Cylinder<T_> const& shape)
{
    return shape;
}

template<typename T_>
inline Cylinder<T_>& shape(Cylinder<T_>& shape)
{
    return shape;
}

template<typename T_>
struct is_shape<Cylinder<T_> >: public boost::mpl::true_ {};

template<typename T_>
struct shape_position_type<Cylinder<T_> >
{
    typedef typename Cylinder<T_>::position_type type;
};

template<typename T_>
struct shape_length_type<Cylinder<T_> > {
    typedef typename Cylinder<T_>::length_type type;
};

template<typename T>
inline typename shape_length_type<Cylinder<T> >::type const& shape_size(Cylinder<T> const& shape)
{
    return shape.radius();
} 

template<typename T>
inline typename shape_length_type<Cylinder<T> >::type& shape_size(Cylinder<T>& shape)
{
    return shape.radius();
} 

namespace std {
template<typename T_>
struct hash<Cylinder<T_> >
{
    typedef Cylinder<T_> argument_type;

    std::size_t operator()(argument_type const& val)
    {
        return hash<typename argument_type::position_type>()(val.position()) ^
            hash<typename argument_type::length_type>()(val.radius()) ^
            hash<typename argument_type::position_type>()(val.unit_z()) ^
            hash<typename argument_type::length_type>()(val.half_length());
    }
};
} // namespace std

#endif /* CYLINDER_HPP */
