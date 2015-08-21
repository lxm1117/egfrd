#ifndef EGFRD_SIMULATOR_FACTORY_HPP
#define EGFRD_SIMULATOR_FACTORY_HPP

#include "ParticleSimulatorFactory.hpp"
#include "EGFRDSimulator.hpp"
#include "CuboidalRegion.hpp"
#include "linear_algebra.hpp"

template<typename Ttraits_>
class EGFRDSimulatorFactory: public ParticleSimulatorFactory<Ttraits_>
{
public:
    typedef ParticleSimulatorFactory<Ttraits_>  base_type;
    typedef Ttraits_                            traits_type;
    typedef typename traits_type::world_type    world_type;
    typedef typename world_type::traits_type    world_traits_type;

    // shorthand typedefs
    typedef typename world_traits_type::length_type     length_type;
    typedef typename world_traits_type::size_type       size_type;
    typedef typename world_traits_type::position_type   position_type;
    typedef typename world_traits_type::rng_type        rng_type;
    typedef typename traits_type::network_rules_type    network_rules_type;

    typedef CuboidalRegion<world_traits_type>           cuboidal_region_type;


public:
    EGFRDSimulatorFactory(rng_type& rng): rng_(rng) {}

    virtual ~EGFRDSimulatorFactory() {}

    virtual EGFRDSimulator<traits_type>* operator()(ParticleModel const& model) const
    {
        length_type const world_size(boost::lexical_cast<length_type>(model["size"]));
        size_type matrix_size(3);
        int dissociation_retry_moves(3);

        try
        {
            matrix_size = boost::lexical_cast<length_type>(model["matrix_size"]);
        }
        catch (not_found const&) {}

        try
        {
            dissociation_retry_moves = boost::lexical_cast<length_type>(model["dissociation_retry_moves"]);
        }
        catch (not_found const&) {}

        position_type const x(divide(position_type(world_size, world_size, world_size), 2));
        boost::shared_ptr<world_type> world(new world_type(world_size, matrix_size));
        
        world->add_structure_type(*model.get_structure_type_by_id(model.get_def_structure_type_id()));   // or else next command fails with Unknown structure_type!
        world->set_def_structure_type_id(model.get_def_structure_type_id());

        world->add_structure(boost::shared_ptr<cuboidal_region_type>(new cuboidal_region_type("world", model.get_def_structure_type_id(), 0,typename cuboidal_region_type::shape_type(x, x))));

        BOOST_FOREACH (boost::shared_ptr<StructureType> st, model.get_structure_types())
        {
            std::string const& type((*st)["type"]);
            // TODO: add surfaces to world
        }

        // Making sure that all the species have a structure_type defined?
        BOOST_FOREACH (boost::shared_ptr<SpeciesType> st, model.get_species_types())
        {
            auto structure_type_id((*st)["structure_type"]);
            world->add_species( typename world_traits_type::species_type(
                    st->id(),
                    //boost::lexical_cast<typename world_traits_type::structure_type_id_type>(structure_type_id.empty() ? model.get_def_structure_type_id() : (*st)["structure_type"]),     <- TODO MS, conversion does not compile?
                    world->get_def_structure_type_id(),
                    boost::lexical_cast<typename world_traits_type::D_type>((*st)["D"]),
                    boost::lexical_cast<length_type>((*st)["radius"])
                    ));
        }

        return new EGFRDSimulator<traits_type>( world, boost::shared_ptr<network_rules_type>( new network_rules_type(model.network_rules())), rng_, dissociation_retry_moves);
    }


///// Member variables
protected:
    rng_type& rng_;
};

#endif /* EGFRD_SIMULATION_HPP */
