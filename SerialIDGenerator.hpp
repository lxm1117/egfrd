#ifndef SERIAL_ID_GENERATOR_HPP
#define SERIAL_ID_GENERATOR_HPP

template<typename Tid_>
struct SerialIDGenerator
{
    typedef Tid_ identifier_type;
    SerialIDGenerator() : next_(identifier_type(1)) {}
    SerialIDGenerator(int lot) : next_(identifier_type(1)) {}           // backward compatible constructor (for Python binding)
    identifier_type operator()() { return ++next_; }

private:
    identifier_type next_;
};


#endif /* SERIAL_ID_GENERATOR_HPP */
