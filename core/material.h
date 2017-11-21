#ifndef __MATERIAL_H_
#define __MATERIAL_H_

#include "scatter_list.h"
#include "../legacy_thomas/material.hh"

template<typename...>
struct material;

template<typename... scatter_types>
struct material<scatter_list<scatter_types...>>
	: public scatter_list<scatter_types...>
{
public:
	HOST material(material_legacy_thomas const & mat);

	inline PHYSICS bool can_reach_vacuum(real kineticEnergy) const;

	real barrier;
};

#include "material.inl"

#endif // __MATERIAL_H_