#ifndef __TUPLE_H_
#define __TUPLE_H_

/*
 * Something vaguely similar to std::tuple.
 * We need something that works on both CUDA devices and host, for holding a
 * list of scattering mechanisms.
 * 
 * NOTE THAT THIS IS NOT BY ANY MEANS A FULLY COMPLIANT TUPLE IMPLEMENTATION.
 * It cannot hold references, move-only types, or any of that.
 * 
 * There is no perfect forwarding anywhere.
 */

#include "variadic.h"

namespace nbl { namespace tuple {


namespace detail
{
	template<size_t, typename T>
	struct tuple_element
	{
		tuple_element(T value) : value_(value) {}
		T value_;
	};

	template<typename sequence, typename... types>
	struct tuple_impl;
	template<size_t... indices, typename... types>
	struct tuple_impl<index_sequence<indices...>, types...>
		: tuple_element<indices, types>...
	{
		explicit tuple_impl(types... elements)
			: tuple_element<indices, types>(elements)...
		{}
	};

	#ifdef __NVCC__
	// TODO: Workaround for nvcc before CUDA 9
	template<size_t i, typename... types>
	struct ___get_helper { using type = type_at_index<i, types...>; };
	#endif
}


// The tuple class itself
template<typename... types>
struct tuple
	: detail::tuple_impl<make_index_sequence<sizeof...(types)>, types...>
{
	using base_t = detail::tuple_impl<make_index_sequence<sizeof...(types)>, types...>;
	explicit tuple(types... elements) :
		base_t(elements...)
	{}
};


// get<i>(tuple)
template<size_t i, typename... types>
#ifdef __NVCC__
PHYSICS typename detail::___get_helper<i, types...>::type& get(tuple<types...>& tup)
#else
PHYSICS type_at_index<i, types...>& get(tuple<types...>& tup)
#endif
{
	detail::tuple_element<i, type_at_index<i, types...>>& base = tup;
	return base.value_;
}

template<size_t i, typename... types>
#ifdef __NVCC__
PHYSICS typename detail::___get_helper<i, types...>::type const & get(tuple<types...> const & tup)
#else
PHYSICS type_at_index<i, types...> const & get(tuple<types...> const & tup)
#endif
{
	detail::tuple_element<i, type_at_index<i, types...>> const & base = tup;
	return base.value_;
}


namespace detail
{
	// visit()
	template<size_t i, typename T, typename F>
	inline PHYSICS typename std::enable_if<i != 0, void>::type visit_impl(T const & tup, F& fun)
	{
		fun(get<i - 1>(tup), i - 1);
		visit_impl<i - 1>(tup, fun);
	}
	template<size_t i, typename T, typename F>
	inline PHYSICS typename std::enable_if<i == 0, void>::type visit_impl(T const & tup, F& fun)
	{}


	// visit_at()
	template<size_t i, typename T, typename F>
	inline PHYSICS typename std::enable_if<i != 0, void>::type visit_at_impl(T const & tup, size_t idx, F& fun)
	{
		if (idx == i - 1) fun(get<i - 1>(tup));
		else visit_at_impl<i - 1>(tup, idx, fun);
	}
	template<size_t i, typename T, typename F>
	inline PHYSICS typename std::enable_if<i == 0, void>::type visit_at_impl(T const & tup, size_t idx, F& fun)
	{}
}

// visit() calls fun on every element of the tuple.
template<typename F, typename... Ts>
PHYSICS void visit(tuple<Ts...> const & tup, F& fun)
{
	detail::visit_impl<sizeof...(Ts)>(tup, fun);
}

// visit_at() calls fun at the element of the tuple indicated by the index.
// If the index is known at compile-time, fun(get<idx>(tuple)) is better.
template<typename F, typename... Ts>
PHYSICS void visit_at(tuple<Ts...> const & tup, size_t idx, F& fun)
{
	detail::visit_at_impl<sizeof...(Ts)>(tup, idx, fun);
}

}} // namespace nbl::tuple

#endif // __TUPLE_H_
