#ifndef __MATERIAL_H_
#define __MATERIAL_H_

#include "scatter_list.h"
#include "../legacy_thomas/material.hh"
#include "../material/hdf5_file.h"

// TODO namespace

template<typename...>
struct material;

/**
 * \brief Material class.
 *
 * Template parameter is a ::scatter_list, i.e. a list of physical scattering
 * mechanisms. These scattering mechansims contain the physics of what actually
 * goes on inside the material.
 */
template<typename... scatter_types>
struct material<scatter_list<scatter_types...>>
	: public scatter_list<scatter_types...>
{
public:
	/**
	 * \brief Constructor. Read from an `e-scatter`-style `.mat` file.
	 * \deprecated Old file format is deprecated and not supported by all
	 *             scattering mechanisms. This function will be removed soon.
	 */
	CPU material(material_legacy_thomas const & mat);

	/**
	 * \brief Constructor. Read from a HDF5 material file.
	 */
	CPU material(nbl::hdf5_file const & mat);

	/**
	 * \brief Return whether an electron with given kinetic energy can leave
	 * the material and reach the vacuum.
	 */
	inline PHYSICS bool can_reach_vacuum(real kineticEnergy) const;

	/**
	 * \brief Barrier energy, that is, the work function plus the Fermi energy.
	 */
	real barrier;
};

#include "material.inl"

#endif // __MATERIAL_H_
