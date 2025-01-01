/*
  ==============================================================================

    TypeHelper.h
    Created: 28 Dec 2024 12:39:33pm
    Author:  John

  ==============================================================================
*/

#pragma once

#include <tuple>
#include <type_traits>
#include <juce_dsp/juce_dsp.h>

// ====== FILTER TYPE REPLACEMENT BEGIN ======

// Helper to replace a type at a specific index in a type list
template <std::size_t Index, typename NewType, typename... Types>
struct ReplaceAtIndex;

// Recursive construct our tuple type by decrementing index
template <std::size_t Index, typename NewType, typename First, typename... Rest>
struct ReplaceAtIndex<Index, NewType, First, Rest...> {
    using type = typename std::conditional_t <
        Index == 0,
        std::tuple<NewType, Rest...>, // Replace at index
        decltype(std::tuple_cat(
            std::tuple<First>{},
            typename ReplaceAtIndex<Index - 1, NewType, Rest...>::type{})) // Keep going
    > ;
};

// Specialization for Index = 0
template <typename NewType, typename First, typename... Rest>
struct ReplaceAtIndex<0, NewType, First, Rest...> {
    using type = std::tuple<NewType, Rest...>;
};

// Specialization for an empty list. This implies we called replaceFilterAtIndex with out of bounds index
template <std::size_t Index, typename NewType>
struct ReplaceAtIndex<Index, NewType> {
    static_assert(Index == 0, "Index out of range");
};

// Utility alias
template <std::size_t Index, typename NewType, typename... Types>
using ReplaceAtIndex_t = typename ReplaceAtIndex<Index, NewType, Types...>::type;

template <typename Tuple>
struct ProcessorChainFromTuple;

template <typename... Types>
struct ProcessorChainFromTuple<std::tuple<Types...>> {
    using type = juce::dsp::ProcessorChain<Types...>;
};

// Function to replace the filter at the specified index
template <std::size_t FilterIndex, typename FilterType, typename... Filters>
auto replaceFilterAtIndex(std::tuple<Filters...>) {
    using newFilterTuple = ReplaceAtIndex_t<FilterIndex, FilterType, Filters...>
    using NewProcessorChain = ProcessorChainFromTuple<newFilterTuple>;
    return NewProcessorChain{};
}

// ====== FILTER TYPE REPLACEMENT END ======