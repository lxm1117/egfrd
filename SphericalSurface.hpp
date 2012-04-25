#ifndef SPHERICAL_SURFACE_HPP
#define SPHERICAL_SURFACE_HPP

#include "Surface.hpp"
#include "Sphere.hpp"
#include "freeFunctions.hpp"

template<typename Ttraits_>
class SphericalSurface
    : public BasicSurfaceImpl<Ttraits_, Sphere<typename Ttraits_::world_type::traits_type::length_type> >
{
public:
    typedef BasicSurfaceImpl<Ttraits_, Sphere<typename Ttraits_::world_type::traits_type::length_type> > base_type;
    typedef typename base_type::traits_type traits_type;
    typedef typename base_type::structure_name_type     structure_name_type;        // THis is just the name of the structure
    typedef typename base_type::structure_id_type       structure_id_type;
    typedef typename base_type::structure_type_id_type  structure_type_id_type;
    typedef typename base_type::shape_type              shape_type;
    typedef typename base_type::rng_type                rng_type;
    typedef typename base_type::position_type           position_type;
    typedef typename base_type::length_type             length_type;
    typedef typename Ttraits_::world_type::species_type species_type;
    typedef std::pair<position_type, position_type>     position_pair_type;

    virtual position_type random_position(rng_type& rng) const
    {
        return position_type(); // TODO
    }

    virtual position_type random_vector(length_type const& r, rng_type& rng) const
    {
        return position_type(); // TODO
    }

    virtual position_type bd_displacement(length_type const& mean, length_type const& r, rng_type& rng) const
    {
        return position_type(); // TODO
    }

    virtual length_type drawR_gbd(Real const& rnd, length_type const& r01, Real const& dt, Real const& D01, Real const& v) const
    {    
        return length_type(); // TODO
    }

    virtual Real p_acceptance(Real const& k_a, Real const& dt, length_type const& r01, position_type const& ipv, 
                                Real const& D0, Real const& D1, Real const& v0, Real const& v1) const
    {    
        return Real(); //TODO
    }

    virtual position_type dissociation_vector( rng_type& rng, length_type const& r01, Real const& dt, 
                                                Real const& D01, Real const& v ) const
    {
        return position_type(); //TODO
    }

    virtual Real get_1D_rate_geminate( Real const& k, length_type const& r01) const
    {
        return Real(); //TODO
    }
    
    virtual Real get_1D_rate_surface( Real const& k, length_type const& r0 ) const
    {
        return Real(); //TODO
    }

    virtual Real particle_reaction_volume( length_type const& r01, length_type const& rl ) const
    {
        return Real(); //TODO
    }
    
    virtual Real surface_reaction_volume( length_type const& r0, length_type const& rl ) const
    {
        return Real(); //TODO
    }

    virtual position_type surface_dissociation_vector( rng_type& rng, length_type const& r0, length_type const& rl ) const
    {
        return position_type(); //TODO  
    }
    
    virtual position_pair_type geminate_dissociation_positions( rng_type& rng, species_type const& s0, species_type const& s1, position_type const& op, length_type const& rl ) const
    {
        return position_pair_type(); //TODO
    }
    
    virtual position_pair_type special_geminate_dissociation_positions( rng_type& rng, species_type const& s_surf, species_type const& s_bulk, position_type const& op_surf, length_type const& rl ) const
    {
        return position_pair_type(); //TODO
    }
    
/*
    virtual bool bounced(position_type const& old_pos, position_type const& new_pos, 
        length_type const& dist_to_surface, length_type const& particle_radius) const
    {  
        return false; //TODO
    }
    
    virtual bool in_reaction_volume( length_type const& dist_to_surface, length_type const& particle_radius, length_type const& rl ) const
    {
        return false; //TODO
    }
*/

    // This should replace above two methods.
    virtual length_type newBD_distance(position_type const& new_pos, length_type const& radius, position_type const& old_pos, length_type const& sigma) const
    {
        return base_type::distance(new_pos);
    }

    virtual length_type minimal_distance(length_type const& radius) const
    {
        return 0.; // TODO
    }

    virtual void accept(ImmutativeStructureVisitor<traits_type> const& visitor) const
    {
        visitor(*this);
    }

    virtual void accept(MutativeStructureVisitor<traits_type> const& visitor)
    {
        visitor(*this);
    }

    // The Constructor
    SphericalSurface(structure_name_type const& name, structure_type_id_type const& sid, structure_id_type const& parent_struct_id, shape_type const& shape)
        : base_type(name, sid, parent_struct_id, shape) {}
};

#endif /* SPHERICAL_SURFACE_HPP */
