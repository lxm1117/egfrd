#ifndef UNASSIGNABLE_ADAPTER_HPP
#define UNASSIGNABLE_ADAPTER_HPP

#include <algorithm>
#include <boost/utility/enable_if.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/range/size.hpp>
#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>
#include <boost/range/rbegin.hpp>
#include <boost/range/rend.hpp>
#include <boost/range/value_type.hpp>
#include <boost/range/size_type.hpp>
#include <boost/range/difference_type.hpp>
#include <boost/range/iterator.hpp>
#include <boost/range/const_iterator.hpp>
#include <boost/range/reverse_iterator.hpp>
#include <boost/range/const_reverse_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>

template<typename T_, template<typename> class TT_, bool Bas_reference_ = false>
struct unassignable_adapter
{
public:
    struct placeholder { char _[sizeof(T_)]; };

    typedef typename TT_<placeholder>::type container_type;
    typedef typename boost::mpl::if_<boost::mpl::bool_<Bas_reference_>,
        container_type&, container_type>::type container_holder_type;
    typedef T_ value_type;

    typedef value_type* pointer;
    typedef value_type const* const_pointer;
    typedef value_type& reference;
    typedef value_type const& const_reference;

public:
    typedef reinterpret_caster<reference, typename container_type::value_type&> caster;
    typedef reinterpret_caster<typename container_type::value_type&, reference> reverse_caster;
    typedef reinterpret_caster<const_reference, typename container_type::value_type const&> const_caster;
    typedef reinterpret_caster<typename container_type::value_type const&, const_reference> const_reverse_caster;

    typedef typename boost::range_size<container_type>::type size_type;
    typedef typename boost::range_difference<container_type>::type difference_type;
    typedef typename boost::range_iterator<container_type>::type placeholder_iterator;
    typedef typename boost::range_const_iterator<container_type>::type const_placeholder_iterator;
    typedef typename boost::range_reverse_iterator<container_type>::type placeholder_reverse_iterator;
    typedef typename boost::range_const_reverse_iterator<container_type>::type const_placeholder_reverse_iterator;

    typedef typename boost::transform_iterator<caster, placeholder_iterator> iterator;

    struct const_iterator: boost::transform_iterator<const_caster, const_placeholder_iterator>
    {
        typedef boost::transform_iterator<const_caster, const_placeholder_iterator> base_type;

        const_iterator(const_iterator const& that): base_type(that) {}

        const_iterator(base_type const& that): base_type(that) {}

        const_iterator(iterator const& that): base_type(that.base(), const_caster()) {}

        const_iterator(const_placeholder_iterator const& iter, const_caster const& functor): base_type(iter, functor) {}

        bool operator==(iterator const& rhs) const
        {
            return base_type::base() == rhs.base();
        }

        bool operator!=(iterator const& rhs) const
        {
            return base_type::base() != rhs.base();
        }

        bool operator<(iterator const& rhs) const
        {
            return base_type::base() < rhs.base();
        }

        bool operator>=(iterator const& rhs) const
        {
            return base_type::base() >= rhs.base();
        }

        bool operator>(iterator const& rhs) const
        {
            return base_type::base() > rhs.base();
        }

        bool operator<=(iterator const& rhs) const
        {
            return base_type::base() <= rhs.base();
        }
    };

    typedef typename boost::transform_iterator<caster, placeholder_reverse_iterator> reverse_iterator;

    struct const_reverse_iterator: boost::transform_iterator<const_caster, const_placeholder_reverse_iterator>
    {
        typedef boost::transform_iterator<const_caster, const_placeholder_reverse_iterator> base_type;

        const_reverse_iterator(const_reverse_iterator const& that): base_type(that) {}

        const_reverse_iterator(base_type const& that): base_type(that) {}

        const_reverse_iterator(reverse_iterator const& that): base_type(that.base(), const_caster()) {}

        const_reverse_iterator(const_placeholder_reverse_iterator const& iter, const_caster const& functor): base_type(iter, functor) {}

        bool operator==(iterator const& rhs) const
        {
            return base_type::base() == rhs.base();
        }

        bool operator!=(iterator const& rhs) const
        {
            return base_type::base() != rhs.base();
        }

        bool operator<(iterator const& rhs) const
        {
            return base_type::base() < rhs.base();
        }

        bool operator>=(iterator const& rhs) const
        {
            return base_type::base() >= rhs.base();
        }

        bool operator>(iterator const& rhs) const
        {
            return base_type::base() > rhs.base();
        }

        bool operator<=(iterator const& rhs) const
        {
            return base_type::base() <= rhs.base();
        }
    };


    size_type size() const
    {
        return boost::size(cntnr_);
    }

    iterator begin()
    {
        return iterator(boost::begin(cntnr_), caster());
    }

    placeholder_iterator pbegin()
    {
        return boost::begin(cntnr_);
    }

    const_iterator begin() const
    {
        return const_iterator(boost::begin(cntnr_), const_caster());
    }

    const_placeholder_iterator pbegin() const
    {
        return boost::begin(cntnr_);
    }

    iterator end()
    {
        return iterator(boost::end(cntnr_), caster());
    }

    placeholder_iterator pend()
    {
        return boost::end(cntnr_);
    }

    const_iterator end() const
    {
        return const_iterator(boost::end(cntnr_), const_caster());
    }

    const_placeholder_iterator pend() const
    {
        return boost::end(cntnr_);
    }

    reverse_iterator rbegin()
    {
        return reverse_iterator(boost::rbegin(cntnr_), caster());
    }

    placeholder_reverse_iterator prbegin()
    {
        return boost::rbegin(cntnr_);
    }

    const_reverse_iterator rbegin() const
    {
        return const_reverse_iterator(boost::rbegin(cntnr_), caster());
    }

    const_placeholder_reverse_iterator prbegin() const
    {
        return boost::rbegin(cntnr_);
    }

    reverse_iterator rend()
    {
        return reverse_iterator(boost::rend(cntnr_), caster());
    }

    placeholder_reverse_iterator prend()
    {
        return boost::rend(cntnr_);
    }

    const_reverse_iterator rend() const
    {
        return const_reverse_iterator(boost::rend(cntnr_), const_caster());
    }

    const_placeholder_reverse_iterator prend() const
    {
        return boost::rend(cntnr_);
    }

    void clear()
    {
        cntnr_.clear();
    }

    void resize(size_type n, T_ const& v)
    {
        cntnr_.resize(n, reinterpret_cast<typename container_type::value_type const&>(v));
    }

    void reserve(size_type n)
    {
        cntnr_.reserve(n);
    }

    size_type capacity()
    {
        return cntnr_.capacity();
    }

    value_type const& at(size_type const& pos) const
    {
        return reinterpret_cast<value_type const&>(cntnr_.at(pos));
    }

    value_type& at(size_type const& pos)
    {
        return reinterpret_cast<value_type&>(cntnr_.at(pos));
    }

    value_type const& operator[](size_type const& pos) const
    {
        return reinterpret_cast<value_type const&>(cntnr_[pos]);
    }

    void set(size_type const& pos, value_type const& v)
    {
        cntnr_[pos] = reinterpret_cast<typename container_type::value_type const&>(v);
    }

    void insert(iterator const& pos, value_type const& v)
    {
        cntnr_.insert(pos.base(), reinterpret_cast<typename container_type::value_type const&>(v));
    }

    void insert(iterator const& pos, size_type n, value_type const v)
    {
        cntnr_.insert(pos.base(), n, reinterpret_cast<typename container_type::value_type const&>(v));
    }

    template<typename Titer_>
    void insert(iterator const& pos, Titer_ const& b, Titer_ const& e)
    {
        typedef boost::transform_iterator<reverse_caster, Titer_> transform_iterator;
        cntnr_.insert(pos.base(), transform_iterator(b, const_reverse_caster()), transform_iterator(e, const_reverse_caster()));
    }

    void push_front(value_type const& v)
    {
        cntnr_.push_front(reinterpret_cast<typename container_type::value_type const&>(v));
    }

    void push_back(value_type const& v)
    {
        cntnr_.push_back(reinterpret_cast<typename container_type::value_type const&>(v));
    }

    void pop_front()
    {
        cntnr_.pop_front();
    }

    void pop_back()
    {
        cntnr_.pop_back();
    }

    reference front()
    {
        return reinterpret_cast<reference>(cntnr_.front());
    }

    const_reference front() const
    {
        return reinterpret_cast<const_reference>(cntnr_.front());
    }

    reference back()
    {
        return reinterpret_cast<reference>(cntnr_.back());
    }

    const_reference back() const
    {
        return reinterpret_cast<const_reference>(cntnr_.back());
    }

    iterator erase(iterator const& pos)
    {
        return iterator(cntnr_.erase(pos.base()), caster());
    }

    iterator erase(iterator const& b, iterator const& e)
    {
        return iterator(cntnr_.erase(b.base(), e.base()), caster());
    }

    unassignable_adapter() {}

    unassignable_adapter(container_holder_type cntnr): cntnr_(cntnr) {}

private:
    container_holder_type cntnr_;
};

#endif /* UNASSIGNABLE_ADAPTER_HPP */
