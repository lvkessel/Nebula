#ifndef __PHYSICS_CONFIG_H_
#define __PHYSICS_CONFIG_H_

#include "physics/inelastic_thomas.h"
#include "physics/inelastic_penn.h"
#include "physics/elastic_thomas.h"
#include "physics/intersect_thomas.h"
#include "core/scatter_list.h"

/*
 * Physics definitions below
 */

// Old e-scatter inelastic models
//template<bool gpu_flag>
//using inelastic_scatter = nbl::scatter::inelastic_thomas<gpu_flag,
//	true, // Optical phonon loss
//	true, // Generate secondaries
//	true, // Random instantaneous momentum for secondary
//	true  // Momentum conservation
//>;

// Penn inelastic model
template<bool gpu_flag>
using inelastic_scatter = nbl::scatter::inelastic_penn<gpu_flag,
	true, // Optical phonon loss
	true, // Generate secondaries
	true, // Random instantaneous momentum for secondary (inner shells only)
	true  // Momentum conservation (inner shells only)
>;

// Default elastic model
template<bool gpu_flag>
using elastic_scatter = nbl::scatter::elastic_thomas<gpu_flag,
	true, // Acoustic phonon loss
	true  // Atomic recoil loss
>;

// Material boundary intersection
using intersect_t = intersect_thomas<
	true, // Quantum-mechanical (probabilistic) transmission
	true, // Interface refraction
	false // Kieft & Bosch empirical interface absorption
>;


// Putting it all together
template<bool gpu_flag>
using scatter_physics = scatter_list<
	inelastic_scatter<gpu_flag>,
	elastic_scatter<gpu_flag>
>;

#endif
