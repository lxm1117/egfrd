#include <map>
#include <cstdio>
#include <boost/foreach.hpp>
#include "PythonAppender.hpp"
#include "../Logger.hpp"

namespace binding {

    typedef std::map<Logger::loglevel, boost::python::object> loglevel_map_type;

    static boost::python::object logging_module;
    static loglevel_map_type loglevel_map;

    static void import_logging_module()
    {
        if (logging_module.ptr() != Py_None) return;

        logging_module = boost::python::object(boost::python::borrowed(PyImport_Import(boost::python::handle<>(PyString_FromString("logging")).get())));
        if (PyErr_Occurred()) boost::python::throw_error_already_set();

        loglevel_map[Logger::loglevel::L_OFF] = getattr(logging_module, "NOTSET");
        loglevel_map[Logger::loglevel::L_DEBUG] = getattr(logging_module, "DEBUG");
        loglevel_map[Logger::loglevel::L_INFO] = getattr(logging_module, "INFO");
        loglevel_map[Logger::loglevel::L_WARNING] = getattr(logging_module, "WARNING");
        loglevel_map[Logger::loglevel::L_ERROR] = getattr(logging_module, "ERROR");
        loglevel_map[Logger::loglevel::L_FATAL] = getattr(logging_module, "CRITICAL");
    }

    class PythonAppender : public LogAppender
    {
    public:
        virtual ~PythonAppender() {}

        virtual void operator()(Logger::loglevel lv, char const* name, char const** chunks)
        {
            std::string msg;
            for (char const** p = chunks; *p; ++p)
                msg.append(*p);
            handle_(makeRecord_(name, lv, "", 0, msg.c_str(), nullptr, nullptr, nullptr));
        }

        virtual void flush()
        {
            flush_();
        }

        PythonAppender(boost::python::object handler) : handler_(handler), flush_(handler_.attr("flush")),
            handle_(handler_.attr("handle")), makeRecord_(handler_.attr("makeRecord")) {}

    private:
        boost::python::object handler_;
        boost::python::object flush_;
        boost::python::object handle_;
        boost::python::object makeRecord_;
    };

    class CppLoggerHandler
    {
    public:
        void* operator new(size_t)
        {
            PyObject* retval = PyObject_New(PyObject, &__class__);
            return retval;
        }

            void operator delete(void* ptr)
        {
            reinterpret_cast<PyObject*>(ptr)->ob_type->tp_free(reinterpret_cast<PyObject*>(ptr));
        }

        static PyObject* __class_init__(const char* name, PyObject* mod)
        {
            using namespace boost::python;

            import_logging_module();
            __name__ = static_cast<std::string>(extract<std::string>(object(borrowed(mod)).attr("__name__"))) + "." + name;
            __class__.tp_name = const_cast<char*>(__name__.c_str());
            __class__.tp_base = reinterpret_cast<PyTypeObject*>(getattr(logging_module, "Handler").ptr());
            if (PyType_Ready(&__class__) < 0) return nullptr;
            return reinterpret_cast<PyObject*>(&__class__);
        }

        static CppLoggerHandler* get_self(PyObject* _self)
        {
            return reinterpret_cast<CppLoggerHandler*>(_self);
        }

        static PyObject* __new__(PyTypeObject* klass, PyObject* arg, PyObject* kwarg)
        {
            try
            {
                if (PyTuple_Size(arg) != 1)
                {
                    PyErr_SetString(PyExc_TypeError, "the number of arguments must be 1");
                    return nullptr;
                }

                CppLoggerHandler* retval(new CppLoggerHandler(*boost::python::extract<Logger*>(PyTuple_GetItem(arg, 0))()));
                if (PyErr_Occurred())
                {
                    boost::python::decref(reinterpret_cast<PyObject*>(retval));
                    return nullptr;
                }
                return reinterpret_cast<PyObject*>(retval);
            }
            catch (boost::python::error_already_set const&) { }
            return nullptr;
        }

        static int __init__(PyObject* self, PyObject* arg, PyObject* kwarg)
        {
            using namespace boost::python;
            return handle<>(allow_null(PyObject_Call(PyObject_GetAttrString(reinterpret_cast<PyObject*>(Py_TYPE(self)->tp_base),
                const_cast<char*>("__init__")), boost::python::handle<>(PyTuple_Pack(1, self)).get(), nullptr))) ? 0 : 1;
        }

        static int __traverse__(PyObject *self, visitproc visit, void *arg)
        {
            Py_VISIT(reinterpret_cast<CppLoggerHandler*>(self)->__dict__.get());
            return 0;
        }

        ~CppLoggerHandler()
        {
            if (__weakreflist__)
            {
                PyObject_ClearWeakRefs(reinterpret_cast<PyObject*>(this));
            }
        }

        static void __dealloc__(void* ptr)
        {
            delete reinterpret_cast<CppLoggerHandler*>(ptr);
        }

        static PyObject* createLock(PyObject* _self)
        {
            CppLoggerHandler* const self(get_self(_self));
            if (!self) return nullptr;
            return boost::python::incref(Py_None);
        }

        static PyObject* acquire(PyObject* _self)
        {
            CppLoggerHandler* const self(get_self(_self));
            if (!self) return nullptr;
            return boost::python::incref(Py_None);
        }

        static PyObject* release(PyObject* _self)
        {
            CppLoggerHandler* const self(get_self(_self));
            if (!self) return nullptr;
            return boost::python::incref(Py_None);
        }

        static PyObject* setLevel(PyObject* _self, PyObject* _level)
        {
            CppLoggerHandler* const self(get_self(_self));
            if (!self) return nullptr;
            Logger::loglevel const level(translate_level(_level));
            if (level == Logger::loglevel::L_OFF)
            {
                PyErr_SetString(PyExc_ValueError, "invalid loglevel");
                return nullptr;
            }
            self->impl_.level(level);
            return boost::python::incref(Py_None);
        }

        static PyObject* getLevel(PyObject* _self)
        {
            CppLoggerHandler* const self(get_self(_self));
            if (!self) return nullptr;
            return boost::python::incref(loglevel_map[self->impl_.level()].ptr());
        }

        static PyObject* emit(PyObject* _self, PyObject* _record)
        {
            CppLoggerHandler* const self(get_self(_self));
            if (!self) return nullptr;
            boost::python::object record(boost::python::borrowed(_record));

            try
            {
                boost::python::handle<> _level(PyObject_GetAttrString(_record, "levelno"));
                boost::python::object msg(boost::python::getattr(boost::python::object(boost::python::borrowed(_self)), "format")(record));
                self->impl_.log(translate_level(_level.get()), "%s", boost::python::extract<char const*>(msg)());
            }
            catch (boost::python::error_already_set const&)
            {
                return nullptr;
            }
            return boost::python::incref(Py_None);
        }

        static PyObject* flush(PyObject* _self, PyObject* arg)
        {
            CppLoggerHandler* const self(get_self(_self));
            if (!self) return nullptr;
            self->impl_.flush();
            return boost::python::incref(Py_None);
        }

        static PyObject* close(PyObject* _self, PyObject* arg)
        {
            CppLoggerHandler* const self(get_self(_self));
            if (!self) return nullptr;
            return boost::python::incref(Py_None);
        }

        static PyObject* translateLevelValue(PyObject* _self, PyObject* arg)
        {
            try
            {
                return boost::python::incref(boost::python::object(translate_level(arg)).ptr());
            }
            catch (boost::python::error_already_set const&) { }
            return nullptr;
        }

        static PyObject* get_logger(PyObject* _self, void* context)
        {
            CppLoggerHandler* const self(get_self(_self));
            if (!self) return nullptr;
            return boost::python::reference_existing_object::apply<Logger*>::type()(self->impl_);
        }


        static PyObject* get___dict__(PyObject *_self, void* context)
        {
            CppLoggerHandler* const self(get_self(_self));
            if (!self) return nullptr;
            return boost::python::incref(self->__dict__.get());
        }

        CppLoggerHandler(Logger& impl) : __weakreflist__(0), __dict__(PyDict_New()), impl_(impl) {}

    protected:
        static Logger::loglevel translate_level(PyObject* _level)
        {
            Logger::loglevel retval(Logger::loglevel::L_OFF);
            boost::python::object level(boost::python::borrowed(_level));
            boost::python::object closest;

            BOOST_FOREACH(loglevel_map_type::value_type const& i, loglevel_map)
            {
                if (i.second <= level && closest < i.second)
                {
                    retval = i.first;
                    closest = i.second;
                }
            }
            return retval;
        }


    protected:
        PyObject_VAR_HEAD static PyTypeObject __class__;
        static std::string __name__;
        PyObject *__weakreflist__;
        boost::python::handle<> __dict__;
        static PyMethodDef __methods__[];
        static PyGetSetDef __getsets__[];
        Logger& impl_;
    };

    std::string CppLoggerHandler::__name__;

    PyMethodDef CppLoggerHandler::__methods__[] = {
        { "createLock", (PyCFunction)CppLoggerHandler::createLock, METH_NOARGS, "" },
        { "acquire", (PyCFunction)CppLoggerHandler::acquire, METH_NOARGS, "" },
        { "release", (PyCFunction)CppLoggerHandler::release, METH_NOARGS, "" },
        { "setLevel", (PyCFunction)CppLoggerHandler::setLevel, METH_O, "" },
        { "getLevel", (PyCFunction)CppLoggerHandler::getLevel, METH_NOARGS, "" },
        { "emit", (PyCFunction)CppLoggerHandler::emit, METH_O, "" },
        { "flush", (PyCFunction)CppLoggerHandler::flush, METH_NOARGS, "" },
        { "close", (PyCFunction)CppLoggerHandler::close, METH_NOARGS, "" },
        { "translateLevelValue", (PyCFunction)CppLoggerHandler::translateLevelValue, METH_O | METH_STATIC, "" },
        { nullptr, nullptr }
    };

    PyGetSetDef CppLoggerHandler::__getsets__[] = {
        { const_cast<char*>("__dict__"), (getter)&CppLoggerHandler::get___dict__, nullptr, nullptr },
        { const_cast<char*>("logger"), (getter)&CppLoggerHandler::get_logger, nullptr, nullptr },
        { nullptr, nullptr }
    };

    PyTypeObject CppLoggerHandler::__class__ = {
        PyObject_HEAD_INIT(&PyType_Type)
        0,					/* ob_size */
        0,                  /* tp_name */
        sizeof(CppLoggerHandler), /* tp_basicsize */
        0,					/* tp_itemsize */
        /* methods */
        (destructor)&CppLoggerHandler::__dealloc__, /* tp_dealloc */
        0,					/* tp_print */
        0,					/* tp_getattr */
        0,					/* tp_setattr */
        0,					/* tp_compare */
        0,					/* tp_repr */
        0,					/* tp_as_number */
        0,                  /* tp_as_sequence */
        0,					/* tp_as_mapping */
        0,					/* tp_hash */
        0,					/* tp_call */
        0,					/* tp_str */
        PyObject_GenericGetAttr,		/* tp_getattro */
        0,					/* tp_setattro */
        0,					/* tp_as_buffer */
        Py_TPFLAGS_HAVE_CLASS | Py_TPFLAGS_HAVE_WEAKREFS | Py_TPFLAGS_BASETYPE, /* tp_flags */
        0,					/* tp_doc */
        CppLoggerHandler::__traverse__,              	/* tp_traverse */
        0,					/* tp_clear */
        0,                  /* tp_richcompare */
        offsetof(CppLoggerHandler, __weakreflist__),    /* tp_weaklistoffset */
        0,                  /* tp_iter */
        0,                  /* tp_iternext */
        CppLoggerHandler::__methods__,		        	/* tp_methods */
        0,					/* tp_members */
        CppLoggerHandler::__getsets__, /* tp_getset */
        0,                  /* tp_base */
        0,                  /* tp_dict */
        0,                  /* tp_descr_get */
        0,                  /* tp_descr_set */
        offsetof(CppLoggerHandler, __dict__),           /* tp_dictoffset */
        CppLoggerHandler::__init__, /* tp_init */
        0,                  /* tp_alloc */
        CppLoggerHandler::__new__,  /*tp_new */
        0                   /* tp_free */
    };

    boost::python::object register_logger_handler_class(char const* name)
    {
        using namespace boost::python;

        import_logging_module();

        PyObject* klass(CppLoggerHandler::__class_init__(name, reinterpret_cast<PyObject*>(scope().ptr())));
        boost::python::object retval(borrowed(klass));
        scope().attr(name) = retval;
        return retval;
    }

    boost::python::objects::class_base register_python_appender_class(char const* name)
    {
        using namespace boost::python;
        typedef PythonAppender impl_type;

        import_logging_module();

        return class_<impl_type, bases<LogAppender>, boost::shared_ptr<impl_type>, boost::noncopyable>(name, init<boost::python::object>())
            ;
    }

} // namespace binding
