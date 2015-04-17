#ifndef ARRAY_TRAITS_HPP
#define ARRAY_TRAITS_HPP

#include <boost/array.hpp>
#include <boost/multi_array.hpp>

template< typename T_ >
struct element_type_of
{
    //typedef typename T_::value_type type;
};

template< typename T_, std::size_t N_ >
struct element_type_of< T_[N_] >
{
    typedef T_ type;
};

template< typename T_, std::size_t N_ >
struct element_type_of< boost::array< T_, N_ > >
{
    typedef T_ type;
};

template< typename T_, typename Talloc_ >
struct element_type_of< boost::multi_array< T_, 1, Talloc_ > >
{
    typedef T_ type;
};

#endif /* ARRAY_TRAITS_HPP */
