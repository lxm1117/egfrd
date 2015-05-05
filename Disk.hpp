#ifndef DISK_HPP
#define DISK_HPP

#include <ostream>
#include <functional>
#include <cmath>
#include "Vector3.hpp"
#include "Shape.hpp"
#include "linear_algebra.hpp"
#include "utils/math.hpp"

template<typename T_>
class Disk
{
public:
    typedef T_          value_type;
    typedef Vector3<T_> position_type;
    typedef T_          length_type;

public:
    // constructors
    Disk() : position_(), radius_(0), unit_z_() {}

    Disk(position_type const& position, length_type const& radius, position_type const& unit_z)
        : position_(position), radius_(radius), unit_z_(unit_z) {}

    bool operator==(const Disk& rhs) const
    {
        return position_ == rhs.position() && radius_ == rhs.radius() && unit_z_ == rhs.unit_z();
    }

    bool operator!=(const Disk& rhs) const
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

    const int dof() const
    {   // degrees of freedom for particle movement
        return 0;
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
    position_type   position_;    // centre.
    length_type     radius_;
    position_type   unit_z_;      // Z-unit_z. should be normalized.
};


template<typename Tstrm_, typename T_>
inline std::basic_ostream<Tstrm_>& operator<<(std::basic_ostream<Tstrm_>& strm, const Disk<T_>& v)
{
    strm << "{" << v.position() << ", " << v.radius() << ", " << v.unit_z() << ", " << "}";
    return strm;
}

template<typename T_>
inline boost::array<typename Disk<T_>::length_type, 2>
to_internal(Disk<T_> const& obj, typename Disk<T_>::position_type const& pos)
{
    // Same as for the cylinder:
    // Return pos relative to position of the disk.
    typedef typename Disk<T_>::position_type    position_type;
    typedef typename Disk<T_>::length_type      length_type;

    const position_type pos_vector(subtract(pos, obj.position()));

    const length_type z(dot_product(pos_vector, obj.unit_z()));               // z can be < 0
    const length_type r(length(subtract(pos_vector, multiply(obj.unit_z(), z))));// r is always >= 0

    return array_gen<typename Disk<T_>::length_type>(r, z);
}

template<typename T_>
inline std::pair<typename Disk<T_>::position_type, std::pair<typename Disk<T_>::length_type, typename Disk<T_>::length_type> >
project_point(Disk<T_> const& obj, typename Disk<T_>::position_type const& pos)
// Calculates the projection of 'pos' onto the disk 'obj' and also returns the coefficient
// for the normal component (z) of 'pos' in the basis of the disk and the distance of the
// projected point to the 'edge' of the disk. Here a positive number means the projected
// point is outside the disk, and a negative numbers means it is 'inside' the disk.
// As a special case, we calculate the projection differently for a position that is in 
// the plane of the disk. Then we imagine that the particle is always 'above' the structure
// and return -1 as a standard, instead of the distance to the disk center.
{
    typedef typename Disk<T_>::length_type   length_type;
    typedef typename Disk<T_>::position_type position_type;

    // Here we do not call 'to_internal' for efficiency
    const position_type pos_vector(subtract(pos, obj.position()));

    const length_type   z(dot_product(pos_vector, obj.unit_z()));
    const position_type r_vector(subtract(pos_vector, multiply(obj.unit_z(), z)));
    const length_type   r(length(r_vector));
    assert(r >= 0.0);

    // The quantities that will be return, with default values
    // for the standard case (pos is not in plane of disk)
    position_type proj_pos(add(obj.position(), r_vector));
    length_type   normal_comp(z);
    length_type   dist_to_edge(r - obj.radius());

    // Special case: pos is in the same plane as the disk
    if (feq(z, 0.0, obj.radius())){ // third argument is typical scale

        proj_pos = obj.position(); // projected position = disk center
        normal_comp = r;
        dist_to_edge = -1.0;
    }

    return std::make_pair(proj_pos, std::make_pair(normal_comp, dist_to_edge));
}

// The same as in case of the plane: project_point_on_surface = project_point
template<typename T_>
inline std::pair<typename Disk<T_>::position_type, std::pair<typename Disk<T_>::length_type, typename Disk<T_>::length_type> >
project_point_on_surface(Disk<T_> const& obj, typename Disk<T_>::position_type const& pos)
{
    return project_point(obj, pos);
}

template<typename T_>
inline typename Disk<T_>::length_type distance(Disk<T_> const& obj, typename Disk<T_>::position_type const& pos)
{
    typedef typename Disk<T_>::length_type length_type;

    /* First compute the (r,z) components of pos in a coordinate system
     * defined by the vectors unitR and unit_z, where unitR is
     * choosen such that unitR and unit_z define a plane in which
     * pos lies. */
    const boost::array<length_type, 2> r_z(to_internal(obj, pos));

    /* Then compute distance to cylinder with zero length. */
    const length_type dz(std::fabs(r_z[1]));
    const length_type dr(r_z[0] - obj.radius());
    length_type distance;

    if (dr > 0)
    {
        // pos is not above the disk.
        // Compute distance to edge.
        distance = std::sqrt(dz * dz + dr * dr);
    }
    else
    {
        // pos is above the disk.
        distance = dz;
    }

    return distance;
}

template<typename T_>
inline std::pair<typename Disk<T_>::position_type, bool> deflect(Disk<T_> const& obj, typename Disk<T_>::position_type const& r0, typename Disk<T_>::position_type const& d)
{
    // Displacements are not deflected on disks (yet),
    // but this function has to be defined for every shape to be used in structure.
    // For now it just returns original pos. + displacement. The changeflage = false.
    return std::make_pair(add(r0, d), false);
}
/*
template<typename T_>
inline typename Disk<T_>::position_type
deflect_back(Disk<T_> const& obj,
typename Disk<T_>::position_type const& r,
typename Disk<T_>::position_type const& u_z  )
{
// Return the vector r without any changes
return r;
}
*/

template<typename T, typename Trng>
inline typename Disk<T>::position_type random_position(Disk<T> const& shape, Trng& rng)
{
    // The disk has only one "legal" position = its center
    return shape.position();
}

template<typename T_>
inline Disk<T_> const& shape(Disk<T_> const& shape)
{
    return shape;
}

template<typename T_>
inline Disk<T_>& shape(Disk<T_>& shape)
{
    return shape;
}

template<typename T_>
struct is_shape<Disk<T_> > : public boost::mpl::true_{};

template<typename T_>
struct shape_position_type < Disk<T_> >
{
    typedef typename Disk<T_>::position_type type;
};

template<typename T_>
struct shape_length_type < Disk<T_> > 
{
    typedef typename Disk<T_>::length_type type;
};

template<typename T>
inline typename shape_length_type<Disk<T> >::type const& shape_size(Disk<T> const& shape)
{
    return shape.radius();
}

template<typename T>
inline typename shape_length_type<Disk<T> >::type& shape_size(Disk<T>& shape)
{
    return shape.radius();
}

namespace std {
    template<typename T_>
    struct hash < Disk<T_> >
    {
        typedef Disk<T_> argument_type;

        std::size_t operator()(argument_type const& val)
        {
            return hash<typename argument_type::position_type>()(val.position()) ^
                hash<typename argument_type::length_type>()(val.radius()) ^
                hash<typename argument_type::position_type>()(val.unit_z());
        }
    };
} // namespace std

#endif /* DISK_HPP */
