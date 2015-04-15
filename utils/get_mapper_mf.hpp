#ifndef GET_MAPPER_MF_HPP
#define GET_MAPPER_MF_HPP

#include <unordered_map>
template<typename Tkey_, typename Tval_>
struct get_mapper_mf
{
    typedef std::unordered_map<Tkey_, Tval_> type;
};

#endif /* GET_MAPPER_MF_HPP */
