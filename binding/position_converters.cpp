#include <boost/python.hpp>
#include "peer/utils.hpp"
#include "peer/numpy/ndarray_converters.hpp"
#include "peer/numpy/wrapped_multi_array.hpp"

#include "binding_common.hpp"

namespace binding {

    //// Converters. These convert C++ types to something that Python can handle.

    typedef peer::util::detail::array_to_ndarray_converter<WorldTraits::position_type, WorldTraits::position_type::value_type, 3> position_to_ndarray_converter;

    struct ndarray_to_position_converter
    {
        typedef Position native_type;

        static void* convertible(PyObject* ptr)
        {
            if (!PyArray_Check(ptr))
            {
                return nullptr;
            }

            PyObject* retval(PyArray_CastToType(reinterpret_cast<PyArrayObject*>(ptr), PyArray_DescrFromType(peer::util::get_numpy_typecode<native_type::value_type >::value), 0));
            if (!retval)
            {
                return nullptr;
            }

            if (PyObject_Size(retval) != 3)
            {
                boost::python::decref(retval);
                return nullptr;
            }

            return retval;
        }

        static void construct(PyObject* ptr, boost::python::converter::rvalue_from_python_storage<native_type>* data)
        {
            PyArrayObject* array_obj = static_cast<PyArrayObject*>(data->stage1.convertible);
            data->stage1.convertible = new(data->storage.bytes) native_type(reinterpret_cast<native_type::value_type*>(PyArray_DATA(array_obj)));
            boost::python::decref(reinterpret_cast<PyObject*>(array_obj));
        }
    };

    //
    struct tuple_to_position_converter
    {
        typedef Position native_type;

        static void* convertible(PyObject* ptr)
        {
            if (!PyTuple_Check(ptr))
            {
                return nullptr;
            }

            if (PyTuple_GET_SIZE(ptr) != 3)
            {
                return nullptr;
            }

            return ptr;
        }

        static void construct(PyObject* ptr, boost::python::converter::rvalue_from_python_storage<native_type>* data)
        {
            data->stage1.convertible = new(data->storage.bytes) native_type(
                PyFloat_AsDouble(PyTuple_GET_ITEM(ptr, 0)),
                PyFloat_AsDouble(PyTuple_GET_ITEM(ptr, 1)),
                PyFloat_AsDouble(PyTuple_GET_ITEM(ptr, 2)));
        }
    };


    struct list_to_position_converter
    {
        typedef Position native_type;

        static void* convertible(PyObject* ptr)
        {
            if (!PyList_Check(ptr))
            {
                return nullptr;
            }

            if (PyList_GET_SIZE(ptr) != 3)
            {
                return nullptr;
            }

            return ptr;
        }

        static void construct(PyObject* ptr, boost::python::converter::rvalue_from_python_storage<native_type>* data)
        {
            data->stage1.convertible = new(data->storage.bytes) native_type(
                PyFloat_AsDouble(PyList_GET_ITEM(ptr, 0)),
                PyFloat_AsDouble(PyList_GET_ITEM(ptr, 1)),
                PyFloat_AsDouble(PyList_GET_ITEM(ptr, 2)));
        }
    };


    ////// Registering master function
    void register_position_converters()
    {
        boost::python::to_python_converter<WorldTraits::position_type, position_to_ndarray_converter>();
        peer::util::to_native_converter<WorldTraits::position_type, ndarray_to_position_converter>();
        peer::util::to_native_converter<WorldTraits::position_type, tuple_to_position_converter>();
        peer::util::to_native_converter<WorldTraits::position_type, list_to_position_converter>();
    }

} // namespace binding
