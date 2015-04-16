#ifndef BOOST_TEST_MODULE
#define BOOST_TEST_MODULE "linear_algebra_test"
#include <boost/test/included/unit_test.hpp>
#endif

#include "linear_algebra.hpp"

BOOST_AUTO_TEST_CASE(test_is_matrix)
{
    BOOST_CHECK((!is_matrix<boost::multi_array<int, 1>, 0>::value));
    BOOST_CHECK((is_matrix<boost::multi_array<int, 1>, 1>::value));
    BOOST_CHECK((!is_matrix<boost::multi_array<int, 1>, 2>::value));
    BOOST_CHECK((!is_matrix<boost::multi_array<int, 2>, 0>::value));
    BOOST_CHECK((!is_matrix<boost::multi_array<int, 2>, 1>::value));
    BOOST_CHECK((is_matrix<boost::multi_array<int, 2>, 2>::value));
    BOOST_CHECK((!is_matrix<boost::multi_array<int, 2>, 3>::value));
}
