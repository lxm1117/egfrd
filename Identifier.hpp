#ifndef IDENTIFIER_HPP
#define IDENTIFIER_HPP

template<typename Tbase_, typename Tserial_>
struct Identifier
{
    typedef Tserial_ value_type;

    Identifier(value_type const& value) : value_(value) { }

    value_type operator++() { return value_++; }

    Tbase_& operator=(Tbase_ const& rhs) { value_ = rhs.value_; }

    value_type const& operator()() const { return value_; }

    bool operator!() const { return value_ == 0; }

    bool operator==(Tbase_ const& rhs) const { return value_ == rhs.value_; }

    bool operator!=(Tbase_ const& rhs) const { return value_ != rhs.value_; }

    bool operator<(Tbase_ const& rhs) const { return value_ < rhs.value_; }

    bool operator>=(Tbase_ const& rhs) const { return value_ >= rhs.value_; }

    bool operator>(Tbase_ const& rhs) const { return value_ > rhs.value_; }

    bool operator<=(Tbase_ const& rhs) const  { return value_ <= rhs.value_; }

    unsigned int serial() { return (unsigned int)value_; }       // this breaks TMP, but makes Python-binding a lot easier

private:
    value_type value_;
};

#endif /* IDENTIFIER_HPP */
