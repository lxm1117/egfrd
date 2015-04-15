#ifndef SPHERE_HPP
#define SPHERE_HPP

#include <ostream>
#include <functional>
#include "Vector3.hpp"
#include "Shape.hpp"
#include "linear_algebra.hpp"

template<typename T_>
class Sphere
{
public:
    typedef T_ value_type;
    typedef Vector3<T_> position_type;
    typedef T_ length_type;

public:
    Sphere()
        : position_(), radius_(0) {}

    Sphere(const position_type& position, const length_type& radius)
        : position_(position), radius_(radius) {}

    bool operator==(const Sphere& rhs) const
    {
        return position_ == rhs.position_ && radius_ == rhs.radius_;
    }

    bool operator!=(const Sphere& rhs) const
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
    
    const int dof() const
    {   // degrees of freedom for particle movement
        return 2;
    }

    std::string show(int precision)
    {
        std::ostringstream strm;
        strm.precision(precision);
        strm << *this;
        return strm.str();
    }

private:
    position_type position_;
    length_type radius_;
};

template<typename Tstrm_, typename Ttraits_, typename T_>
inline std::basic_ostream<Tstrm_, Ttraits_>& operator<<(std::basic_ostream<Tstrm_, Ttraits_>& strm,
        const Sphere<T_>& v)
{
    strm << "{" << v.position() <<  ", " << v.radius() << "}";
    return strm;
}

template<typename T_>
inline typename Sphere<T_>::length_type
to_internal(Sphere<T_> const& obj, typename Sphere<T_>::position_type const& pos)
// The function calculates the coefficients to express 'pos' into the base of the sphere 'obj'
{
    typedef typename Sphere<T_>::position_type position_type;
    position_type pos_vector(subtract(pos, obj.position()));

    // Todo. If we ever need it.
    return length(pos_vector);
}

template<typename T_>
inline std::pair<typename Sphere<T_>::position_type,
                 std::pair<typename Sphere<T_>::length_type,
                           typename Sphere<T_>::length_type> >
project_point(Sphere<T_> const& obj, typename Sphere<T_>::position_type const& pos)
{
    typename Sphere<T_>::length_type r(to_internal(obj, pos));

    // The projection of a point on a sphere is always the centerpoint of the sphere.
    return std::make_pair(obj.position(),
                          std::make_pair(r, 0.0) );
}

template<typename T_>
inline std::pair<typename Sphere<T_>::position_type,
                 std::pair<typename Sphere<T_>::length_type,
                           typename Sphere<T_>::length_type> >
project_point_on_surface(Sphere<T_> const& obj,
                typename Sphere<T_>::position_type const& pos)
{
    // Todo. If we ever need it.
    // The projection of a point on a sphere.
    return std::make_pair(typename Sphere<T_>::position_type(),
                          std::make_pair(typename Sphere<T_>::length_type(),
                                         typename Sphere<T_>::length_type()) );
}

template<typename T_>
inline typename Sphere<T_>::length_type
distance(Sphere<T_> const& obj, typename Sphere<T_>::position_type const& pos)
{
    return distance(pos, obj.position()) - obj.radius();
}

template<typename T_>
inline std::pair<typename Sphere<T_>::position_type, bool>
deflect(Sphere<T_> const& obj, typename Sphere<T_>::position_type const& r0, typename Sphere<T_>::position_type const& d)
{
    // Displacements are not deflected on spheres (yet),
    // but this function has to be defined for every shape to be used in structure.
    // For now it just returns the new position. The changeflag = 0.
    return std::make_pair( add(r0, d), false );
}
/*
template<typename T_>
inline typename Sphere<T_>::position_type
deflect_back(Sphere<T_> const& obj,
        typename Sphere<T_>::position_type const& r,
        typename Sphere<T_>::position_type const& u_z  )
{
    // Return the vector r without any changes
    return r;
}
*/
template<typename T_>
inline Sphere<T_> const& shape(Sphere<T_> const& shape)
{
    return shape;
}

template<typename T_>
inline Sphere<T_>& shape(Sphere<T_>& shape)
{
    return shape;
}

template<typename T, typename Trng>
inline typename Sphere<T>::position_type
random_position(Sphere<T> const& shape, Trng& rng)
{
    return add(shape.position(),
                create_vector<typename Sphere<T>::position_type>(
                    shape.radius() * rng(),
                    shape.radius() * rng(),
                    shape.radius() * rng())); 
}

template<typename T_>
struct is_shape<Sphere<T_> >: public boost::mpl::true_ {};

template<typename T_>
struct shape_position_type<Sphere<T_> > {
    typedef typename Sphere<T_>::position_type type;
};

template<typename T_>
struct shape_length_type<Sphere<T_> > {
    typedef typename Sphere<T_>::length_type type;
};

template<typename T>
inline typename shape_length_type<Sphere<T> >::type const& shape_size(Sphere<T> const& shape)
{
    return shape.radius();
} 

template<typename T>
inline typename shape_length_type<Sphere<T> >::type& shape_size(Sphere<T>& shape)
{
    return shape.radius();
} 

namespace std {
template<typename T_>
struct hash<Sphere<T_> >
{
    typedef Sphere<T_> argument_type;

    std::size_t operator()(argument_type const& val)
    {
        return hash<typename argument_type::position_type>()(val.position()) ^
            hash<typename argument_type::length_type>()(val.radius());
    }
};
} // namespace std

#endif /* SPHERE_HPP */
