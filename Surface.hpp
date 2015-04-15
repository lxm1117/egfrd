#ifndef SURFACE_HPP
#define SURFACE_HPP

#include <ostream>
#include <functional>
#include <sstream>
#include "ParticleSimulationStructure.hpp"
#include "Disk.hpp"
#include "Cylinder.hpp"
#include "Sphere.hpp"
#include "Plane.hpp"

template<typename Ttraits_>
class Surface: public ParticleSimulationStructure<Ttraits_>
{
public:
    typedef ParticleSimulationStructure<Ttraits_>       base_type;
    typedef typename base_type::structure_name_type     structure_name_type;    // This is just a string
    typedef typename base_type::structure_id_type       structure_id_type;
    typedef typename base_type::structure_type_id_type  structure_type_id_type;
    typedef typename base_type::length_type             length_type;

public:
    virtual ~Surface() {}

    // Constructor
    Surface(structure_name_type const& name, structure_type_id_type const& sid, structure_id_type const& parent_struct_id):
         base_type(name, sid, parent_struct_id) {}

//    virtual length_type minimal_distance(length_type const& radius) const = 0;
};





template<typename Ttraits_, typename Tshape_>
class BasicSurfaceImpl: public Surface<Ttraits_>
{
public:
  typedef Surface<Ttraits_>                           base_type;
  typedef Tshape_                                     shape_type;

    typedef typename base_type::structure_name_type         structure_name_type;
    typedef typename base_type::structure_id_type           structure_id_type;
    typedef typename base_type::structure_type_id_type      structure_type_id_type;
    typedef typename base_type::length_type                 length_type;
    typedef typename base_type::position_type               position_type;
    typedef std::pair<length_type, length_type>             components_pair_type;
    typedef std::pair<position_type, components_pair_type>  projected_type;
    typedef std::pair<position_type, bool>                  position_flag_pair_type;

public:
    virtual ~BasicSurfaceImpl() {}

    shape_type& shape()
    {
        return shape_;
    }

    shape_type const& shape() const
    {
        return shape_;
    }
    
    virtual bool operator==(Structure<Ttraits_> const& rhs) const
    {
        BasicSurfaceImpl const* _rhs(dynamic_cast<BasicSurfaceImpl const*>(&rhs));
        return _rhs && base_type::id_ == rhs.id() && base_type::sid_ == rhs.sid() && shape_ == _rhs->shape();
    }
    
    virtual std::size_t hash() const
    {
        return std::hash<structure_name_type>()(base_type::name_) ^ std::hash<structure_type_id_type>()(base_type::sid_) ^ std::hash<shape_type>()(shape());
    }

    virtual std::string as_string() const
    {
        std::ostringstream out;
        out << "Surface(" << base_type::id_ << ":" << shape() << ")";
        return out.str();
    }

    virtual projected_type project_point(position_type const& pos) const
    {
        return ::project_point(shape(), pos);
    }
    
    virtual projected_type project_point_on_surface(position_type const& pos) const
    {
        return ::project_point_on_surface(shape(), pos);
    }
    
    virtual length_type distance(position_type const& pos) const
    {
        return ::distance(shape(), pos);
    }

    virtual position_type const& position() const
    {
        return ::shape_position(shape());
    }

    virtual position_flag_pair_type deflect(position_type const& pos0, position_type const& displacement) const
    {
        return ::deflect(shape(), pos0, displacement);
    }
/*    
    virtual position_type deflect_back(position_type const& pos, position_type const& u_z) const
    {
        return ::deflect_back(shape(), pos, u_z);
    }
*/
    // Constructor
    BasicSurfaceImpl(structure_name_type const& name, structure_type_id_type const& sid, structure_id_type const& parent_struct_id, shape_type const& shape)
        : base_type(name, sid, parent_struct_id), shape_(shape) {}

protected:
    shape_type shape_;
};

#endif /* SURFACE_HPP */
