#include <set>
#include <queue>
#include <assert.h>

#ifndef UNIQUE_PRIORITY_QUEUE_H
#define UNIQUE_PRIORITY_QUEUE_H

template<class T>
struct do_nothing {
    void operator()(const T&){
        // intentionally left blank
    }
};

/**
 * Priority queue with unique (according to FuncCompare) elements of type T where the sorting is based on CostCompare.
 * If NDEBUG is *not* defined, there are some assertions that help catching errors in the provided comparision functions.
 */
template<class T, class CleanObsoleteElement = do_nothing<T>, class CostCompare = std::less<T>, class FuncCompare = CostCompare>
class unique_priority_queue
{
public:
    typedef typename std::priority_queue<T, std::vector<T>, CostCompare>::size_type size_type;

    /**
     * Return true if the element was inserted into the queue.
     * This happens if equivalent element is present or if the new element has a lower cost associated to it.
     * False is returned if no insertion into the queue took place.
     */
    bool push(const T& v)
    {
        const auto& insertion_pair = membership_.insert(v);
        if(insertion_pair.second)
        {
            queue_.push(v);
        }
        else if (CostCompare()(*(insertion_pair.first), v))
        {
            const auto number_erased = membership_.erase(*(insertion_pair.first));
            assert(number_erased == 1);
            CleanObsoleteElement()(*(insertion_pair.first));
            const auto inserted = membership_.insert(v);
            assert(inserted.second);
            queue_ = std::priority_queue<T, std::vector<T>, CostCompare>();
            for(const auto& element : membership_) {
                queue_.push(element);
            }
            assert(queue_.size() == membership_.size());
            return true;

        }
        assert(queue_.size() == membership_.size());
        return insertion_pair.second;
    }

    void pop()
    {
        assert(!queue_.empty() && queue_.size() == membership_.size());

        const auto& top_element = queue_.top();
        const auto number_erased = membership_.erase(top_element);

        assert(number_erased == 1);

        queue_.pop();
        assert(queue_.size() == membership_.size());
    }

    const T& top() const
    {
        assert(!queue_.empty());
        return queue_.top();
    }

    bool empty() const
    {
        assert(queue_.size() == membership_.size());
        return queue_.empty();
    }

    size_type size() const
    {
        return queue_.size();
    }

private:
    std::priority_queue<T, std::vector<T>, CostCompare> queue_;
    std::set<T, FuncCompare> membership_;
};
#endif
