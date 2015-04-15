#ifndef REGION_HPP
#define REGION_HPP

#include <ostream>
#include <functional>
#include <sstream>

#include "ParticleSimulationStructure.hpp"
#include "Box.hpp"

template<typename Ttraits_>
class Region: public ParticleSimulationStructure<Ttraits_>
{
public:
    typedef ParticleSimulationStructure<Ttraits_>       base_type;
    typedef typename base_type::structure_name_type     structure_name_type;
    typedef typename base_type::structure_id_type       structure_id_type;
    typedef typename base_type::structure_type_id_type  structure_type_id_type;

public:
    virtual ~Region() {}

    // Constructor
    Region(structure_name_type const& name, structure_type_id_type const& sid, structure_id_type const& parent_struct_id)
        : base_type(name, sid, parent_struct_id) {}
};



template<typename Ttraits_, typename Tshape_>
class BasicRegionImpl: public Region<Ttraits_>
{
public:
    typedef Region<Ttraits_>                            base_type;
    typedef Tshape_                                     shape_type;

    typedef typename base_type::structure_name_type         structure_name_type;
    typedef typename base_type::structure_id_type           structure_id_type;
    typedef typename base_type::structure_type_id_type      structure_type_id_type;
    typedef typename base_type::length_type                 length_type;
    typedef typename base_type::position_type               position_type;
    typedef std::pair<position_type, bool>                  position_flag_pair_type;
    typedef std::pair<length_type, length_type>             components_pair_type;
    typedef std::pair<position_type, components_pair_type>  projected_type;

public:
    virtual ~BasicRegionImpl() {}

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
        BasicRegionImpl const* _rhs(dynamic_cast<BasicRegionImpl const*>(&rhs));
        return _rhs &&
               base_type::id_ == rhs.id() &&
               base_type::sid_ == rhs.sid() &&
               shape_ == _rhs->shape();
    }

    virtual std::size_t hash() const
    {    // Displacements are not deflected on cylinders (yet),
    // but this function has to be defined for every shape to be used in structure
    // For now it just returns the new position
        return std::hash<structure_name_type>()(base_type::name_) ^
            std::hash<structure_type_id_type>()(base_type::sid_) ^
            std::hash<shape_type>()(shape());
    }

    virtual std::string as_string() const
    {
        std::ostringstream out;
        out << "Region(" << base_type::id_ << ":" << shape() << ")";
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
    // The constructor
    BasicRegionImpl(structure_name_type const& name, structure_type_id_type const& sid, structure_id_type const& parent_struct_id, shape_type const& shape)
        : base_type(name, sid, parent_struct_id), shape_(shape) {}

protected:
    shape_type shape_;
};

#endif /* REGION_HPP */
