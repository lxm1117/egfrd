#ifndef MULTI_HPP
#define MULTI_HPP

#include <boost/scoped_ptr.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/foreach.hpp>

#include <gsl/gsl_math.h>

#include "exceptions.hpp"
#include "Domain.hpp"
#include "ParticleContainer.hpp"
#include "ParticleContainerBase.hpp"
#include "Sphere.hpp"
#include "BDSimulator.hpp"
#include "BDPropagator.hpp"
#include "Logger.hpp"
#include "PairGreensFunction.hpp"
#include "Transaction.hpp"
#include "VolumeClearer.hpp"
#include "utils/array_helper.hpp"
#include "utils/range.hpp"


template<typename Ttraits_>
class MultiParticleContainer
    : public Ttraits_::world_type::particle_container_type
{
public:
    typedef ParticleContainerUtils<typename Ttraits_::world_type::traits_type> utils;
    typedef typename Ttraits_::world_type       world_type;
    typedef typename world_type::traits_type    traits_type;

    // shorthands types for used types
    typedef typename traits_type::length_type               length_type;
    typedef typename traits_type::size_type                 size_type;
    typedef typename traits_type::position_type             position_type;
    typedef typename traits_type::species_type              species_type;
    typedef typename traits_type::species_id_type           species_id_type;
    typedef typename traits_type::particle_type             particle_type;
    typedef typename traits_type::particle_id_type          particle_id_type;
    typedef typename traits_type::structure_type_type       structure_type_type;
    typedef typename traits_type::structure_type_id_type    structure_type_id_type;
    typedef typename traits_type::structure_type            structure_type;
    typedef typename traits_type::structure_id_type         structure_id_type;

    typedef typename particle_type::shape_type              particle_shape_type;
    typedef Transaction<traits_type>                        transaction_type;

    // types inherited from the ParticleContainerBase class
    typedef typename world_type::particle_container_type::structure_types_range             structure_types_range;
    typedef typename world_type::particle_container_type::structures_range                  structures_range;
    typedef typename world_type::particle_container_type::particle_id_pair                  particle_id_pair;
    typedef typename world_type::particle_container_type::structure_id_pair                 structure_id_pair;
    typedef typename world_type::particle_container_type::structure_id_set                  structure_id_set;
    typedef typename world_type::particle_container_type::particle_id_pair_and_distance_list    particle_id_pair_and_distance_list;
    typedef typename world_type::particle_container_type::structure_id_pair_and_distance_list   structure_id_pair_and_distance_list;
    typedef typename world_type::particle_container_type::position_structid_pair_type       position_structid_pair_type;

    typedef typename Ttraits_::network_rules_type       network_rules_type;
    typedef typename Ttraits_::reaction_rule_type       reaction_rule_type;
    typedef typename network_rules_type::reaction_rules reaction_rules;

    typedef abstract_limited_generator<particle_id_pair>                    particle_id_pair_generator;
    typedef std::pair<particle_id_pair, length_type>                        particle_id_pair_and_distance;
    //    typedef unassignable_adapter<particle_id_pair_and_distance>             particle_id_pair_and_distance_list;
    typedef std::map<particle_id_type, particle_type>                       particle_map;
    typedef sized_iterator_range<typename particle_map::const_iterator>     particle_id_pair_range;
    typedef std::pair<const Real, const Real> real_pair;

private:
    typedef std::map<species_id_type, species_type>                         species_map;
    typedef select_second<typename species_map::value_type>                 species_second_selector_type;

public:
    typedef boost::transform_iterator < species_second_selector_type,
        typename species_map::const_iterator > species_iterator;
    typedef sized_iterator_range<species_iterator>                          species_range;



    //
    virtual ~MultiParticleContainer() {}

    virtual size_type num_particles() const
    {
        return particles_.size();
    }

    virtual length_type world_size() const
    {
        return world_.world_size();
    }

    // Species stuff
    virtual species_type const& get_species(species_id_type const& id) const
    {
        return world_.get_species(id);
    }

    // Start structure stuff
    virtual boost::shared_ptr<structure_type> get_structure(structure_id_type const& id) const
    {
        return world_.get_structure(id);
    }
    virtual structures_range get_structures() const
    {
        return world_.get_structures();     // TODO now gets all structures in world, -> make structure local to Multi
    }
    virtual boost::shared_ptr<structure_type> get_some_structure_of_type(structure_type_id_type const& sid) const
    {
        return world_.get_some_structure_of_type(sid);
    }
    // virtual structure_id_type add_structure(structure_type const& structure); // TODO add structure from the world to multi
    template <typename Tstructid_pair_>
    bool update_structure(Tstructid_pair_ const& structid_pair)
    {
        return world_.update_structure(structid_pair);
    }
    virtual bool remove_structure(structure_id_type const& id)
    {
        return world_.remove_structure(id);
    }
    virtual structure_id_set get_structure_ids(structure_type_id_type const& sid) const
    {
        return world_.get_structure_ids(sid);
    }
    virtual structure_id_type get_def_structure_id() const
    {
        return world_.get_def_structure_id();
    }
    virtual structure_id_pair_and_distance_list* get_close_structures(position_type const& pos, structure_id_type const& current_struct_id,
        structure_id_type const& ignore) const
    {
        return world_.get_close_structures(pos, current_struct_id, ignore);
    }
    // End structure stuff

    // StructureType stuff
    // virtual bool add_structure_type(structure_type_type const& structure_type); // TODO add structure_type from the world to multi
    virtual structure_type_type get_structure_type(structure_type_id_type const& sid) const
    {
        return world_.get_structure_type(sid);      // TODO make local just like structures.
    }
    virtual structure_types_range get_structure_types() const
    {
        return world_.get_structure_types();        // TODO ditto
    }
    virtual structure_type_id_type get_def_structure_type_id() const
    {
        return world_.get_def_structure_type_id();  // This information stays in the world.
    }


    virtual particle_id_pair new_particle(species_id_type const& sid, structure_id_type const& structure_id,
        position_type const& pos)
    {
        particle_id_pair const retval(world_.new_particle(sid, structure_id, pos));
        particles_.insert(retval);
        add_species_to_multi(sid);
        return retval;
    }

    virtual bool update_particle(particle_id_pair const& pi_pair)
    {
        world_.update_particle(pi_pair);
        typename particle_map::iterator const i(particles_.find(pi_pair.first));
        if (i != particles_.end())
        {
            (*i).second = pi_pair.second;
            add_species_to_multi(pi_pair.second.sid());
            return false;
        }
        else
        {
            particles_.insert(i, pi_pair);
            add_species_to_multi(pi_pair.second.sid());
            return true;
        }
    }

    virtual bool remove_particle(particle_id_type const& id)
    {
        species_.erase(get_particle(id).second.sid());
        world_.remove_particle(id);
        return 0!=particles_.erase(id);     //return true if actually erased something
    }

    virtual particle_id_pair get_particle(particle_id_type const& id) const
    {
        typename particle_map::const_iterator i(particles_.find(id));
        if (particles_.end() == i)
        {
            throw not_found(std::string("No such particle: id=")
                + boost::lexical_cast<std::string>(id));
        }
        return *i;
    }

    virtual bool has_particle(particle_id_type const& id) const
    {
        return particles_.end() != particles_.find(id);
    }

    virtual particle_id_pair_and_distance_list* check_overlap(particle_shape_type const& s) const
    {
        return check_overlap(s, array_gen<particle_id_type>());
    }

    virtual particle_id_pair_and_distance_list* check_overlap(particle_shape_type const& s, particle_id_type const& ignore) const
    {
        return check_overlap(s, array_gen(ignore));
    }

    virtual particle_id_pair_and_distance_list* check_overlap(particle_shape_type const& s, particle_id_type const& ignore1, particle_id_type const& ignore2) const
    {
        return check_overlap(s, array_gen(ignore1, ignore2));
    }

    template<typename Tsph_, typename Tset_>
    particle_id_pair_and_distance_list* check_overlap(Tsph_ const& s, Tset_ const& ignore) const
    {
        typename utils::template overlap_checker<particle_id_pair_and_distance_list, Tset_> checker(ignore);
        for (typename particle_map::const_iterator i(particles_.begin()),
            e(particles_.end());
            i != e; ++i)
        {
            length_type const dist(world_.distance(shape((*i).second), s.position()));
            if (dist < s.radius())
            {
                checker(i, dist);
            }
        }
        return checker.result();
    }

    virtual structure_id_pair_and_distance_list* check_surface_overlap(particle_shape_type const& s, position_type const& old_pos, structure_id_type const& current,
        length_type const& sigma) const
    {
        return world_.check_surface_overlap(s, old_pos, current, sigma);
    }

    virtual structure_id_pair_and_distance_list* check_surface_overlap(particle_shape_type const& s, position_type const& old_pos, structure_id_type const& current,
        length_type const& sigma, structure_id_type const& ignore) const
    {
        return world_.check_surface_overlap(s, old_pos, current, sigma, ignore);
    }

    virtual structure_id_pair_and_distance_list* check_surface_overlap(particle_shape_type const& s, position_type const& old_pos, structure_id_type const& current,
        length_type const& sigma, structure_id_type const& ignore1, structure_id_type const& ignore2) const
    {
        return world_.check_surface_overlap(s, old_pos, current, sigma, ignore1, ignore2);
    }

    virtual particle_id_pair_generator* get_particles() const
    {
        return make_range_generator<particle_id_pair>(particles_);
    }

    virtual transaction_type* create_transaction()
    {
        return new TransactionImpl<MultiParticleContainer>(*this);
    }

    virtual length_type distance(position_type const& lhs,
        position_type const& rhs) const
    {
        return world_.distance(lhs, rhs);
    }

    virtual position_type apply_boundary(position_type const& v) const
    {
        return world_.apply_boundary(v);
    }

    virtual length_type apply_boundary(length_type const& v) const
    {
        return world_.apply_boundary(v);
    }

    virtual position_structid_pair_type apply_boundary(position_structid_pair_type const& pos_struct_id) const
    {
        return world_.apply_boundary(pos_struct_id);
    }

    virtual position_type cyclic_transpose(position_type const& p0, position_type const& p1) const
    {
        return world_.cyclic_transpose(p0, p1);
    }

    virtual length_type cyclic_transpose(length_type const& p0, length_type const& p1) const
    {
        return world_.cyclic_transpose(p0, p1);
    }

    virtual position_structid_pair_type cyclic_transpose(position_structid_pair_type const& pos_struct_id,
        structure_type const& structure) const
    {
        return world_.cyclic_transpose(pos_struct_id, structure);
    }

    particle_id_pair_range get_particles_range() const
    {
        return particle_id_pair_range(particles_.begin(), particles_.end(),
            particles_.size());
    }

    void add_species_to_multi(species_id_type const& sid)
    {
        typename species_map::const_iterator i(species_.find(sid));
        if (species_.end() == i)
        {
            species_[sid] = world_.get_species(sid);
        }
    }

    species_range get_species_in_multi() const
    {
        return species_range(
            species_iterator(species_.begin(), species_second_selector_type()),
            species_iterator(species_.end(), species_second_selector_type()),
            species_.size());
    }

    /* Determines the largest diffusion constant and largest drift coefficient inside the mpc.
       and whether movement is dominated by diffusion or drift for a given length scale.
       Then it returns the typical time to travel that length scale based on which type
       of movement is dominant.                                                             */
    Real get_min_tau_Dv(Real r_typical) const
    {
        Real tau_D(0.0), tau_v(0.0); // the estimated times to travel r_typical by diffusion / convection
        Real tau_dominant(0.0), tau_min(-1.0); // always keep tau_min negative initially!

        const Real LARGE_PREF(GSL_POSINF);
        assert(r_typical > 0.0);

        BOOST_FOREACH(species_type s, get_species_in_multi())
        {
            // First we calculate the typical times to travel r_typical by
            // diffusion and convection, respectively.
            if (s.D() == 0 && s.v() == 0)

                tau_dominant = GSL_POSINF; // static particle, will never move

            else{

                tau_D = s.D() > 0.0 ? 2.0 * gsl_pow_2(r_typical) / s.D() : LARGE_PREF * r_typical / s.v();
                tau_v = s.v() > 0.0 ? r_typical / s.v() : LARGE_PREF * 2.0 * gsl_pow_2(r_typical) / s.D();

                tau_dominant = tau_v < 0.1 * tau_D ? tau_v : tau_D;
                // We believe that convective movement is dominant if it takes a tenth of the
                // time to travel the typical distance by convection as compared to diffusion,
                // and then also the associated time is the relevant timescale.
                // Note that the 0.1 prefactor is to ensure that we pick the convective
                // timescale only when the movement is clearly dominated by convection,
                // i.e. when both timescales are comparable we prefer the diffusion timescale.
                if (tau_dominant == tau_D){ LOG_DEBUG(("Diffusion dominates on given reaction length")); }
                else                     { LOG_DEBUG(("Convection dominates on given reaction length")); }
            }

            // Now compare if within the current set of species this is the smallest time scale.
            assert(tau_dominant > 0.0);
            if (tau_min < 0.0 || tau_dominant < tau_min)
                tau_min = tau_dominant;
        }

        assert(tau_min > 0.0);

        return tau_min;
    }

    /* Returns smallest particle radius inside the mpc. */
    Real get_min_radius() const
    {
        Real radius_min(std::numeric_limits<Real>::max());

        BOOST_FOREACH(species_type s, get_species_in_multi())
        {
            if (radius_min > s.radius())
                radius_min = s.radius();
        }

        return radius_min;
    }

    /* Functions returns largest _1D_ intrinsic reaction rate in the multi. */
    Real get_max_rate(network_rules_type const& rules) const
    {
        Real k_max(0.);
        int i = 0, j = 0;

        // Rates for particle-surface interactions
        BOOST_FOREACH(particle_id_pair pp, get_particles_range())
        {
            // s = species of the current particle
            species_type const s(get_species(pp.second.sid()));

            // Get list of close structures with ids and distances from current particle's position           
            const boost::scoped_ptr<const structure_id_pair_and_distance_list>   close_struct_id_distance(
                get_close_structures(pp.second.position(), pp.second.structure_id(), pp.second.structure_id()));

            // Take the closest one and pass it as an id-and-distance pair
            const std::pair<boost::shared_ptr<structure_type>, length_type>   struct_and_dist(
                close_struct_id_distance ? std::make_pair(close_struct_id_distance->at(0).first.second, close_struct_id_distance->at(0).second)
                : std::make_pair(get_structure(pp.second.structure_id()), std::numeric_limits<length_type>::max()));

            //structure_id_and_distance_pair const struct_id_and_dist( 
            //    get_closest_surface( pp.second.position(), pp.second.structure_id() ) );    // only ignore structure that the particle is on.

            // If the structure is within a specified range
            // Here we assume that the user defined the reaction rules for allowed combinations of origin and
            // target structure type (if this is not the case the simulation will fail at a later stage when
            // propagation is attempted).
            if (struct_and_dist.second < 2.0 * s.radius())
            {
                // Get the reaction rule for this particle-structure interaction
                structure_type_id_type const struct_sid(struct_and_dist.first->sid());
                reaction_rules const& rrules(rules.query_reaction_rule(s.id(), struct_sid));
                if (::size(rrules) == 0)
                    continue; // no reaction rule for this structure type

                // If there is rules, determine the largest on-rate for this interaction
                for (typename boost::range_const_iterator<reaction_rules>::type
                    it(boost::begin(rrules)), e(boost::end(rrules)); it != e; ++it)
                {
                    Real const k(struct_and_dist.first->get_1D_rate_surface((*it).k(), s.radius()));

                    if (k_max < k)
                        k_max = k;
                }
            }

        }

        // since surface rates are not devided by 2 to compensate for double reaction attempts.
        k_max *= 2.0;

        // Rates for particle-particle interactions
        BOOST_FOREACH(species_type s0, get_species_in_multi())
        {
            j = 0;
            BOOST_FOREACH(species_type s1, get_species_in_multi())
            {
                if (j++ < i)
                    continue;

                reaction_rules const& rrules(rules.query_reaction_rule(s0.id(), s1.id()));
                if (::size(rrules) == 0)
                    continue;   // no rule for reactions between these two species

                for (typename boost::range_const_iterator<reaction_rules>::type
                    it(boost::begin(rrules)), e(boost::end(rrules)); it != e; ++it)
                {
                    // We have found a valid rule - get the rate
                    const length_type r01(s0.radius() + s1.radius());
                    Real k(0.0), k0, k1;

                    // To access the rate functions we first find some structure of the structure type of the 
                    // species considered and then calculate the modified rate. This is a necessary workaround
                    // because as yet the rate functions are methods to the structure, not the structure type.

                    // If one of the two particles lives in the default structure (bulk) we always take the 
                    // 3D rate, i.e. call the get_1D_rate_geminate() method of the bulk structure.
                    if (s0.structure_type_id() == get_def_structure_type_id())

                        k = get_some_structure_of_type(s0.structure_type_id())->get_1D_rate_geminate((*it).k(), r01);

                    else if (s1.structure_type_id() == get_def_structure_type_id())

                        k = get_some_structure_of_type(s1.structure_type_id())->get_1D_rate_geminate((*it).k(), r01);

                    else
                        // If both particles live on lower dimensionality structures more care has to be taken in determining
                        // what is the right modifier function for the bare rate.
                        // For now by default we take the maximal rate as determined from get_1D_rate_geminate() of both structures.
                        // However, if one of the two interacting species is static (D=0, v=0), the "dimensionality" of the motion
                        // is exclusively determined by the mobile species, and we have to take the modifier function of that one.
                        // FIXME Taking the max by default may lead to a rate which actually is higher than the one that should
                        // be used in BD propagation and therefore wastes resources.
                        // Until now this case basically only comprises the rod-particle/cap-particle interaction. Since
                        // both rates are equal in this case we do not make any approximation in that case. However, if we
                        // allow for other types of lower dimensionality particle-particle reactions we should fix this.
                    {
                        // In case that one of the two particles is static, make sure we modify the rate using the
                        // modifier function of the structure that holds the *mobile* species:
                        if (s0.D() == 0.0 && s0.v() == 0.0) // s1 is the mobile species
                            k = get_some_structure_of_type(s1.structure_type_id())->get_1D_rate_geminate((*it).k(), r01);

                        else if (s1.D() == 0.0 && s1.v() == 0.0) // s0 is the mobile species
                            k = get_some_structure_of_type(s0.structure_type_id())->get_1D_rate_geminate((*it).k(), r01);

                        else{
                            // Pick the higher of the modified rates
                            k0 = get_some_structure_of_type(s0.structure_type_id())->get_1D_rate_geminate((*it).k(), r01);
                            k1 = get_some_structure_of_type(s1.structure_type_id())->get_1D_rate_geminate((*it).k(), r01);

                            k = k0 > k1 ? k0 : k1;
                            // FIXME This is not very elegant yet and may cause that a rate that is way too high is picked 
                            // if the species' structures have different dimensionalities
                        }
                    }
                    // NOTE: If a structure of the required structure type can not be found a not_found exception
                    // is risen. This however should never happen, because whenever a particle of species s0 (s1)
                    // is in the system also at least one structure of the associated structure type should exist.

                    // Compare with the fastest rate found so far
                    k_max = k > k_max ? k : k_max;
                }
            }
            i++;
        }

        return k_max;
    }

    /* Function returns the timestep and reaction length (rl) for BD propegator.

       --dt is calulated with the constraints:
       (1) The reaction length is equal to step_size_factor (ssf) * r_min,
       where r_min is the radius of the smallest particle in the multi.
       (2) The largest acceptance probability in the multi is smaller than Pacc_max (.01).
       (3) particles escape the multi with a maximum step size in the order of the
       reaction length. (Dmax * dt ~ (ssf * r_min)**2 ).

       TODO: This function should be a method of the Multi Class,
       but I put it in the mpc (multi particle container) such that we can use it in python.

       PROBLEM: for certain parameters (large k) dt can be very small and the simulation will slow down.
       */
    real_pair determine_dt_and_reaction_length(network_rules_type const& rules, Real const& step_size_factor, Real const& dt_hardcore_min = -1.0) const
    {

        const Real k_max(get_max_rate(rules));
        const Real r_min(get_min_radius());

        // The following gives us the typical timescale to travel the length step_size_factor * r_min,
        // taking into account all species in the Multi and automatically comparing whether their movement
        // is dominated by convection or diffusion, respectively.
        const Real tau_Dv(get_min_tau_Dv(step_size_factor * r_min));

        const Real Pacc_max(0.1); // Maximum allowed value of the acceptance probability. // TESTING was 0.01
        // This should be kept very low (max. 0.01), otherwise the approximation of
        // treating the reaction as two sequential attempts (first move, then react)
        // might break down!
        Real dt, dt_temp;

        if (k_max > 0)
        {
            // step_size_factor * r_min is the reaction length
            // Here it is assumed that the RL is linear in any dimension,
            // which requires it to be very small!
            dt_temp = 2. * Pacc_max * step_size_factor * r_min / k_max;
            dt = std::min(dt_temp, tau_Dv); // tau_Dv is upper limit of dt.

            if (dt == tau_Dv){ LOG_DEBUG(("Time step set by timescale of motion")); }
            else            { LOG_DEBUG(("Time step set by largest reaction rate, k_max = %e, r_min = %e", k_max, r_min)); }
        }
        else
            dt = tau_Dv;

        if (dt < dt_hardcore_min){
            dt = dt_hardcore_min;
            LOG_WARNING(("Setting timestep to hard-coded minimal bound, dt = %e", dt));
        }

        return real_pair(dt, step_size_factor * r_min);
    }

    // The constructor
    MultiParticleContainer(world_type& world) : world_(world) {}

    ////// Member variables
private:
    world_type&         world_;     // reference to the world. The MultiParticleContainer reflects a part of the World(ParticleContainer)
    particle_map        particles_; // local copy of the particles (the particles in the Multi are a subset of the particles in the world)
    species_map         species_;   // local copy of the species
};




////////////////////////
template<typename Tsim_>
class Multi : public Domain < typename Tsim_::traits_type >
{
public:
    typedef Tsim_ simulator_type;
    typedef typename simulator_type::traits_type    traits_type;
    typedef typename traits_type::world_type        world_type;
    typedef Domain<traits_type>                     base_type;

    // shorthand typenames that we use a lot
    typedef typename world_type::length_type           length_type;
    typedef typename world_type::size_type             size_type;
    typedef typename world_type::position_type         position_type;
    typedef typename world_type::particle_type         particle_type;
    typedef typename world_type::particle_id_type      particle_id_type;
    typedef typename particle_type::shape_type         particle_shape_type;
    typedef typename world_type::species_type          species_type;
    typedef typename world_type::species_id_type       species_id_type;
    typedef typename world_type::structure_type        structure_type;
    typedef typename world_type::particle_id_pair      particle_id_pair;
    typedef typename traits_type::network_rules_type                network_rules_type;
    typedef typename traits_type::reaction_rule_type                reaction_rule_type;
    typedef typename network_rules_type::reaction_rules             reaction_rules;
    typedef typename traits_type::shell_id_type                     shell_id_type;
    typedef typename traits_type::domain_id_type                    identifier_type;
    typedef typename traits_type::template shell_generator<
        typename simulator_type::sphere_type>::type                 spherical_shell_type;

    typedef std::pair<const typename traits_type::shell_id_type, spherical_shell_type>  spherical_shell_id_pair;
    typedef std::pair<particle_id_pair, length_type>                                            particle_id_pair_and_distance;
    typedef unassignable_adapter<particle_id_pair_and_distance>  particle_id_pair_and_distance_list;
    typedef typename traits_type::reaction_record_type reaction_record_type;
    typedef std::pair<const Real, const Real>                                                   real_pair;

private:
    typedef std::map<species_id_type, species_type>                             species_map;
    typedef select_second<typename species_map::value_type>                     species_second_selector_type;

public:
    typedef boost::transform_iterator < species_second_selector_type,
        typename species_map::const_iterator > species_iterator;
    typedef sized_iterator_range<species_iterator>                              species_range;

    typedef std::map<shell_id_type, spherical_shell_type>                       spherical_shell_map;
    typedef sized_iterator_range<typename spherical_shell_map::const_iterator>  spherical_shell_id_pair_range;
    typedef MultiParticleContainer<traits_type>                                 multi_particle_container_type;

    enum class event_kind
    {
        NONE,
        ESCAPE,
        REACTION,
        NUM_MULTI_EVENT_KINDS
    };

private:
    struct last_reaction_setter : ReactionRecorder < reaction_record_type >
    {
        virtual ~last_reaction_setter() {}

        virtual void operator()(reaction_record_type const& rec)
        {
            outer_.last_reaction_.swap(const_cast<reaction_record_type&>(rec));
        }

        last_reaction_setter(Multi& outer) : outer_(outer) {}

        Multi& outer_;
    };

    struct volume_clearer : VolumeClearer < particle_shape_type, particle_id_type >
    {
        virtual ~volume_clearer() {}

        virtual bool operator()(particle_shape_type const& shape, particle_id_type const& ignore)
        {
            if (!outer_.within_shell(shape))
            {
                outer_.last_event_ = event_kind::ESCAPE;
                return outer_.clear_volume(shape, ignore);
            }
            return true;
        }

        virtual bool operator()(particle_shape_type const& shape, particle_id_type const& ignore0, particle_id_type const& ignore1)
        {
            if (!outer_.within_shell(shape))
            {
                outer_.last_event_ = event_kind::ESCAPE;
                return outer_.clear_volume(shape, ignore0, ignore1);
            }
            return true;
        }

        volume_clearer(Multi& outer) : outer_(outer) {}

        Multi& outer_;
    };

    friend struct volume_clearer;

public:
    virtual ~Multi() {}

    virtual char const* type_name() const
    {
        return "Multi";
    }

    virtual std::string as_string() const
    {
        return (boost::format(
            "%s(id=%s, event=%s, last_time=%.16g, dt=%.16g, particles=[%s])") %
            type_name() %
            boost::lexical_cast<std::string>(base_type::id_).c_str() %
            boost::lexical_cast<std::string>(base_type::event_.first).c_str() %
            base_type::last_time_ % base_type::dt_ %
            stringize_and_join(
            make_select_first_range(pc_.get_particles_range()),
            ", ")).str();
    }

    Multi(identifier_type const& id, simulator_type& main, Real dt_factor)
        : base_type(id), main_(main), pc_(*main.world()), dt_factor_(dt_factor),
        shells_(), last_event_(event_kind::NONE)
    {
        //TODO Do not base dt and rl on all particles in the world but only those in the multi.
        BOOST_ASSERT(dt_factor > 0.);
        base_type::dt_ = dt_factor_ * BDSimulator<traits_type>::determine_dt(*main_.world());
    }

    event_kind const& last_event() const
    {
        return last_event_;
    }

    reaction_record_type const& last_reaction() const
    {
        return last_reaction_;
    }

    bool has_particle(particle_id_type const& pid) const
    {
        return pc_.has_particle(pid);
    }

    bool add_particle(particle_id_pair const& pp)
    {
        return pc_.update_particle(pp);
    }

    bool add_shell(spherical_shell_id_pair const& sp)
    {
        spherical_shell_id_pair new_sp(sp);
        new_sp.second.did() = base_type::id();
        return shells_.insert(new_sp).second;
    }

    spherical_shell_id_pair_range get_shells() const
    {
        return spherical_shell_id_pair_range(shells_.begin(), shells_.end(), shells_.size());
    }

    virtual typename Domain<traits_type>::size_type num_shells() const
    {
        return shells_.size();
    }

    virtual typename Domain<traits_type>::size_type multiplicity() const
    {
        return pc_.num_particles();
    }

    virtual void accept(ImmutativeDomainVisitor<traits_type> const& visitor) const
    {
        visitor(*this);
    }

    virtual void accept(MutativeDomainVisitor<traits_type> const& visitor)
    {
        visitor(*this);
    }

    bool within_shell(particle_shape_type const& sphere) const
    {
        for (typename spherical_shell_map::const_iterator
            i(shells_.begin()), e(shells_.end()); i != e; ++i)
        {
            spherical_shell_id_pair const& sp(*i);
            position_type ppos(main_.world()->cyclic_transpose(sphere.position(), (sp).second.position()));
            //TODO for newBD scheme, add reaction length to sphere.radius().
            if (distance(ppos, (sp).second.shape().position()) < (sp).second.shape().radius() - sphere.radius())
            {
                return true;
            }
        }
        return false;
    }

    bool clear_volume(particle_shape_type const& shape, particle_id_type const& ignore) const
    {
        LOG_DEBUG(("clear_volume was called here."));
        main_.clear_volume(shape, base_type::id_);
        boost::scoped_ptr<particle_id_pair_and_distance_list> overlapped(
            main_.world()->check_overlap(shape, ignore));
        if (overlapped && ::size(*overlapped))
        {
            return false;
        }
        return true;
    }

    bool clear_volume(particle_shape_type const& shape, particle_id_type const& ignore0, particle_id_type const& ignore1) const
    {
        LOG_DEBUG(("clear_volume was called here."));
        main_.clear_volume(shape, base_type::id_);
        boost::scoped_ptr<particle_id_pair_and_distance_list> overlapped(
            main_.world()->check_overlap(shape, ignore0, ignore1));
        if (overlapped && ::size(*overlapped))
        {
            return false;
        }
        return true;
    }

    typename multi_particle_container_type::particle_id_pair_range
        get_particles_range() const
    {
        return pc_.get_particles_range();
    }

    void step()
    {
        boost::scoped_ptr<typename multi_particle_container_type::transaction_type> tx(pc_.create_transaction());
        last_reaction_setter rs(*this);
        volume_clearer vc(*this);

        BDPropagator<traits_type> ppg(
            *tx, *main_.network_rules(), main_.rng(),
            base_type::dt_,
            1 /* FIXME: dissociation_retry_moves */, &rs, &vc,
            make_select_first_range(pc_.get_particles_range()));

        last_event_ = event_kind::NONE;

        while (ppg())
        {
            if (last_reaction_)
            {
                last_event_ = event_kind::REACTION;
                break;
            }
        }
    }

protected:
    simulator_type& main_;
    multi_particle_container_type pc_;
    Real dt_factor_;
    spherical_shell_map shells_;
    event_kind last_event_;
    reaction_record_type last_reaction_;

    static Logger& log_;
};

template<typename Tsim_>
Logger& Multi<Tsim_>::log_(Logger::get_logger("ecell.Multi"));

#endif /* MULTI_HPP */
