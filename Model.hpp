#ifndef MODEL_HPP
#define MODEL_HPP

#include <map>
#include <unordered_map>

#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/noncopyable.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/range/iterator_range.hpp>

#include "Defs.hpp"
#include "SerialIDGenerator.hpp"
#include "SpeciesTypeID.hpp"
#include "SpeciesType.hpp"
#include "utils/range.hpp"
#include "utils/pair.hpp"

#include "NetworkRules.hpp"

//class NetworkRules;

class Model: private boost::noncopyable
{
public:
    typedef SpeciesType                         species_type_type;
    typedef species_type_type::identifier_type  species_id_type;

private:
    typedef SerialIDGenerator<species_id_type>                                  species_type_id_generator_type;
    typedef std::map<species_id_type, boost::shared_ptr<species_type_type> >    species_type_map_type;
    typedef select_second<species_type_map_type::value_type>                    species_second_selector_type;

    typedef std::unordered_map<std::string, std::string>                        string_map_type;

public:
    typedef boost::transform_iterator<species_second_selector_type,
            species_type_map_type::const_iterator>                  species_type_iterator;
    typedef boost::iterator_range<species_type_iterator>            species_type_range;
    typedef NetworkRules                                            network_rules_type;
    typedef string_map_type::const_iterator                         string_map_iterator;
    typedef boost::iterator_range<string_map_iterator>              attributes_range;

public:
    // The model class is more or less an implementation of the network rules
    // The network rules define the possible reactions between species.

    // Constructor
    GFRD_CLASS Model();

    virtual ~Model() {};

    NetworkRules& network_rules() const { return *network_rules_; };

    // Add a species to the model
    GFRD_CLASS void add_species_type(boost::shared_ptr<species_type_type> const& species);

    // Get a species by the species 'id'
    GFRD_CLASS boost::shared_ptr<species_type_type> get_species_type_by_id(species_id_type const& id) const;

    // Get all the species in the model
    GFRD_CLASS species_type_range get_species_types() const;

    GFRD_CLASS std::string const& operator[](std::string const& name) const;

    GFRD_CLASS std::string& operator[](std::string const& name) { return attrs_[name]; };
    // Not sure. Get attribute by name?

    // Get all the attributes
    attributes_range attributes() const
    {
        return attributes_range(attrs_.begin(), attrs_.end());
    }

//////// Member variables (why are they public?)
public:
    species_type_id_generator_type  species_type_id_generator_; // The id generator which makes sure that all species have a unique id
    species_type_map_type           species_type_map_;          // mapping: species_type_id->species_type
    boost::scoped_ptr<NetworkRules> network_rules_;             // The network rules in the model
    string_map_type                 attrs_;                     // All the attributes
};


#endif /* MODEL_HPP */
