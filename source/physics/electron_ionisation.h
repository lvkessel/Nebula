#ifndef __ELECTRON_IONISATION_H_
#define __ELECTRON_IONISATION_H_

#include <functional>
#include "../common/util/table_2D.h"
#include "../common/util/range.h"
#include "../common/util/random.h"
#include "../legacy_thomas/material.hh"
#include "../material/hdf5_file.h"

namespace nbl { namespace scatter {

/**
 * \brief Stores electron ionization probabilities.
 *
 * TODO: this class is only used by ::nbl::scatter::inelastic_thomas. It should
 * really just be part of that class.
 */
template<bool gpu_flag>
class electron_ionisation
{
public:
	/**
	 * \brief Get inner shell energy.
	 */
	PHYSICS real sample(real K, util::random_generator<gpu_flag>& rng) const
	{
		const real x = logr(K);
		const real y = rng.unit();
		return _ionisation_table.get_rounddown(x, y);
	}

	static CPU electron_ionisation create(material_legacy_thomas const & mat)
	{
		util::geomspace<double> K_range(K_min, K_max, K_cnt);
		util::linspace<double> P_range(0, 1, P_cnt);

		electron_ionisation ei;
		ei._ionisation_table = util::table_2D<real, gpu_flag>::create(logr(K_min), logr(K_max), K_cnt, 0, 1, P_cnt);
		ei._ionisation_table.mem_scope([&](real** ionisation_vector)
		{
			for (int y = 0; y < P_cnt; ++y)
			{
				const double P = P_range[y];
				for (int x = 0; x < K_cnt; ++x)
				{
					const double omega0 = K_range[x]; // in eV
					const double margin = 10; // Magic KB number, in eV

					double binding = -1;
					if (omega0 > 100)
					{
						binding = mat.ionization_energy((omega0 + margin)*constant::ec, P) / constant::ec;
						if (binding < 50)
							binding = -1;
					}
					if (binding < 0)
					{
						binding = mat.outer_shell_ionization_energy(omega0*constant::ec) / constant::ec;
						if (binding < 0)
							binding = -1;
					}

					ionisation_vector[y][x] = (real)binding;
				}
			}
		});

		return ei;
	}

	static CPU electron_ionisation create(hdf5_file const & mat)
	{
		util::geomspace<units::quantity<double>> K_range(K_min*units::eV, K_max*units::eV, K_cnt);
		util::linspace<units::quantity<double>> P_range(0*units::dimensionless, 1*units::dimensionless, P_cnt);

		electron_ionisation ei;
		ei._ionisation_table = util::table_2D<real, gpu_flag>::create(logr(K_min), logr(K_max), K_cnt, 0, 1, P_cnt);
		ei._ionisation_table.mem_scope([&](real** ionisation_vector)
		{
			const auto icdf_table = mat.get_table_axes<2>("ionization/binding_icdf");

			// Get outer_shells vector, containing outer-shell energies sorted from low to high.
			const auto outer_shell_table = mat.get_table_axes<1>("ionization/outer_shells");
			std::vector<units::quantity<double>> outer_shells;
			for (double osi : outer_shell_table)
			{
				const units::quantity<double> value = osi * outer_shell_table.unit;
				if (value < 100*units::eV)
					outer_shells.push_back(value);
			}
			std::sort(outer_shells.begin(), outer_shells.end());

			// Create the simulation table
			for (int y = 0; y < P_cnt; ++y)
			{
				const units::quantity<double> P = P_range[y];
				for (int x = 0; x < K_cnt; ++x)
				{
					const units::quantity<double> K = K_range[x];
					const units::quantity<double> margin = 10 * units::eV; // Magic KB number

					units::quantity<double> binding = -1 * units::eV;
					if (K > 100*units::eV)
					{
						binding = icdf_table.get_rounddown(K + margin, P);
						if (binding < 50*units::eV || !std::isfinite(binding.value))
							binding = -1 * units::eV;
					}
					if (binding < 0*units::eV)
					{
						// Find largest outer shell less than or equal to K
						auto outer_shell_iterator = std::lower_bound(outer_shells.rbegin(),
							outer_shells.rend(), K, std::greater<units::quantity<double>>{});
						if (outer_shell_iterator != outer_shells.rend())
							binding = (*outer_shell_iterator);
					}

					ionisation_vector[y][x] = real(binding / units::eV);
				}
			}
		});
		return ei;
	}

	template<bool source_gpu_flag>
	static CPU electron_ionisation create(electron_ionisation<source_gpu_flag> const & source)
	{
		electron_ionisation target;
		target._ionisation_table = util::table_2D<real, gpu_flag>::create(source._ionisation_table);
		return target;
	}

	static CPU void destroy(electron_ionisation & ei)
	{
		util::table_2D<real, gpu_flag>::destroy(ei._ionisation_table);
	}

private:
	/**
	 * \brief Ionization table, representing probability of ionizing a given
	 * inner or outer shell.
	 *
	 * Specifically, stores `binding energy / eV` as function of
	 *   - x axis: `log(kinetic energy / eV)`
	 *   - y axis: cumulative probability (between 0 and 1)
	 *
	 * Be sure to use .get_rounddown() for non-interpolated values
	 */
	util::table_2D<real, gpu_flag> _ionisation_table;

	template<bool>
	friend class electron_ionisation;
};

}} // namespace nbl::scatter

#endif // __ELECTRON_IONISATION_H_
