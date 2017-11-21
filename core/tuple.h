#ifndef __TUPLE_H_
#define __TUPLE_H_

/*
 * Something vaguely similar to std::tuple.
 * We need something that works on both CUDA devices and host, for holding a
 * list of scattering mechanisms.
 * 
 * NOTE THAT THIS IS NOT BY ANY MEANS A FULLY COMPLIANT TUPLE IMPLEMENTATION.
 * It was cannot hold references, move-only types, or any of that.
 * 
 * There is no perfect forwarding, or any of that, anywhere.
 */

namespace nbl { namespace tuple {

/*
 * index_sequence
 */
template <size_t... Indices>
struct index_sequence
{};

template <size_t Start, typename index_sequence, size_t End>
struct make_index_sequence_impl;

template <size_t Start, size_t... Indices, size_t End>
struct make_index_sequence_impl<Start, index_sequence<Indices...>, End>
{
	typedef typename make_index_sequence_impl<Start + 1,
		index_sequence<Indices..., Start>, End>::type type;
};

template <size_t End, size_t... Indices>
struct make_index_sequence_impl<End, index_sequence<Indices...>, End>
{
	typedef index_sequence<Indices...> type;
};

template <size_t End, size_t Start = 0>
struct make_index_sequence
{
	static_assert(Start <= End, "make_index_sequence: invalid parameters");
	typedef typename make_index_sequence_impl<Start, index_sequence<>, End>::type type;
};

template <size_t End, size_t Start = 0>
using make_index_sequence_t = typename make_index_sequence<End, Start>::type;


/*
 * type_at_index
 */
template <size_t I, typename Head, typename... Tail>
struct type_at_index
{
	using type = typename type_at_index<I - 1, Tail...>::type;
};
template <typename Head, typename... Tail>
struct type_at_index<0, Head, Tail...>
{
	using type = Head;
};
template <size_t I, typename... Types>
using type_at_index_t = typename type_at_index<I, Types...>::type;


/*
 * tuple implementation
 */
template <size_t, typename T>
struct tuple_element
{
	tuple_element(T value) : value_(value) {}
	T value_;
};

template <typename Sequences, typename... Types>
struct tuple_impl; // undefined base; parameters are named for clarity only
template <size_t... Indices, typename... Types>
struct tuple_impl<index_sequence<Indices...>, Types...>
	: tuple_element<Indices, Types>...
{
	explicit tuple_impl(Types... elements)
		: tuple_element<Indices, Types>(elements)...
	{}
};

#ifdef __NVCC__
// TODO: Workaround for nvcc before CUDA 8
template<size_t i, typename... Ts>
struct ___get_helper { using type = type_at_index_t<i, Ts...>; };
#endif

template <typename... Types>
struct tuple
	: tuple_impl<make_index_sequence_t<sizeof...(Types)>, Types...>
{
	using base_t = tuple_impl<make_index_sequence_t<sizeof...(Types)>, Types...>;
	explicit tuple(Types... elements) :
		base_t(elements...)
	{}

	template <size_t I, typename... Ts>
#ifdef __NVCC__
	friend PHYSICS typename ___get_helper<I, Ts...>::type& get(tuple<Ts...>& tup);
#else
	friend PHYSICS type_at_index_t<I, Ts...>& get(tuple<Ts...>& tup);
#endif
};


/*
 * get()
 */
template <size_t I, typename... Types>
#ifdef __NVCC__
PHYSICS typename ___get_helper<I, Types...>::type& get(tuple<Types...>& tup)
#else
PHYSICS type_at_index_t<I, Types...>& get(tuple<Types...>& tup)
#endif
{
	tuple_element<I, type_at_index_t<I, Types...>>& base = tup;
	return base.value_;
}

template <size_t I, typename... Types>
#ifdef __NVCC__
PHYSICS typename ___get_helper<I, Types...>::type const & get(tuple<Types...> const & tup)
#else
PHYSICS type_at_index_t<I, Types...> const & get(tuple<Types...> const & tup)
#endif
{
	tuple_element<I, type_at_index_t<I, Types...>> const & base = tup;
	return base.value_;
}


/*
 * visit()
 */
template <size_t I>
struct visit_impl
{
	template <typename T, typename F>
	static PHYSICS void visit(T& tup, F& fun)
	{
		fun(get<I - 1>(tup), I-1);
		visit_impl<I - 1>::visit(tup, fun);
	}
};
template <>
struct visit_impl<0>
{
	template <typename T, typename F>
	static PHYSICS void visit(T& tup, F& fun) {}
};
template <typename F, typename... Ts>
PHYSICS void visit(tuple<Ts...> const & tup, F& fun)
{
	visit_impl<sizeof...(Ts)>::visit(tup, fun);
}

/*
 * visit_at()
 */
template <size_t I>
struct visit_at_impl
{
	template <typename T, typename F>
	static PHYSICS void visit(T& tup, size_t idx, F& fun)
	{
		if (idx == I - 1) fun(get<I - 1>(tup));
		else visit_at_impl<I - 1>::visit(tup, idx, fun);
	}
};
template <>
struct visit_at_impl<0>
{
	template <typename T, typename F>
	static PHYSICS void visit(T& tup, size_t idx, F& fun) { /* Too large index provided! */ }
};
template <typename F, typename... Ts>
PHYSICS void visit_at(tuple<Ts...> const & tup, size_t idx, F& fun)
{
	visit_at_impl<sizeof...(Ts)>::visit(tup, idx, fun);
}

}} // namespace nbl::tuple

#endif // __TUPLE_H_