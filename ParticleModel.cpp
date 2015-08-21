#include <algorithm>
#include <boost/type_traits/remove_pointer.hpp>
#include <boost/lexical_cast.hpp>
#include "Defs.hpp"
#include "utils/fun_composition.hpp"
#include "utils/fun_wrappers.hpp"
#include "ParticleModel.hpp"
#include "BasicNetworkRulesImpl.hpp"

// Constructor
GFRD_CLASS ParticleModel::ParticleModel() 
{
    // TODO add default structure_type for the bulk?
    boost::shared_ptr<ParticleModel::structure_type_type> default_structure_type(new StructureType());
    add_structure_type(default_structure_type);
    default_structure_type_id_ = default_structure_type->id();
    
}

// Add a structure type to the model
GFRD_CLASS void ParticleModel::add_structure_type(boost::shared_ptr<structure_type_type> const& structure_type)
{
    // std::pair<structure_type_map_type::iterator, bool> r(
    //     structure_type_map_.insert(std::make_pair(structure_type->id(), structure_type)));
    // if (!r.second)
    // {
    //     throw already_exists(
    //         (boost::format("structure_type id \"%s\" is already used by %s") %
    //             structure_type->id() %
    //             boost::lexical_cast<std::string>(*(*(r.first)).second)).str());
    // }
    structure_type->bind_to_model(this, species_type_id_generator_());
    structure_type_map_.insert(std::make_pair(structure_type->id(), structure_type));
}

// Get a structure type from the model
GFRD_CLASS boost::shared_ptr<ParticleModel::structure_type_type> ParticleModel::get_structure_type_by_id(structure_type_id_type const& id) const
{
    structure_type_map::const_iterator i(structure_type_map_.find(id));
    if (structure_type_map_.end() == i)
    {
        throw not_found(std::string("Unknown structure_type (id=") + boost::lexical_cast<std::string>(id)+")");
    }

    return (*i).second;
}
