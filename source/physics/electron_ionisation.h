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

					ionisation_vector[x][y] = (real)binding;
				}
			}
		});

		return ei;
	}

	static CPU electron_ionisation create(hdf5_file const & mat)
	{
		electron_ionisation ei;

		ei._ionisation_table = mat.fill_table2D<real>("ionization/binding_icdf");
		auto K_range = mat.get_log_dimscale("ionization/binding_icdf", 0,
			ei._ionisation_table.width());
		auto P_range = mat.get_lin_dimscale("ionization/binding_icdf", 1,
			ei._ionisation_table.height());
		ei._ionisation_table.set_scale(
			std::log(K_range.front()/units::eV), std::log(K_range.back()/units::eV),
			P_range.front(), P_range.back());
		ei._ionisation_table.mem_scope([&](real** ionisation_vector)
		{
			// Get outer_shells vector, containing outer-shell energies sorted from low to high.
			std::vector<units::quantity<double>> outer_shells;
			{
				const auto outer_shell_table = mat.fill_table1D<double>("ionization/outer_shells");
				const auto unit = mat.get_unit("ionization/outer_shells");
				for (int i = 0; i < outer_shell_table.width(); ++i)
				{
					const units::quantity<double> value = outer_shell_table(i) * unit;
					if (value < 100*units::eV)
						outer_shells.push_back(value);
				}
				std::sort(outer_shells.begin(), outer_shells.end());
			}

			// Create the simulation table
			const auto unit = mat.get_unit("ionization/binding_icdf");
			for (int x = 0; x < K_range.size(); ++x)
			{
				auto K = K_range[x];
				for (int y = 0; y < P_range.size(); ++y)
				{
					const real P = P_range[y];

					// Magic KB number
					const units::quantity<double> margin = 10 * units::eV;

					units::quantity<double> binding = -1 * units::eV;
					if (K > 100*units::eV)
					{
						binding = ei._ionisation_table.get_rounddown(
							std::log((K+margin)/units::eV), P) * units::eV;
						if (binding < 50*units::eV || !std::isfinite(binding.value))
							binding = -1 * units::eV;
					}
					if (binding < 0*units::eV)
					{
						// Find largest outer shell less than or equal to K
						auto outer_shell_iterator = std::lower_bound(outer_shells.rbegin(),
							outer_shells.rend(), K, std::greater<units::quantity<double>>{});
						if (outer_shell_iterator != outer_shells.rend())
							binding = *outer_shell_iterator;
					}

					ionisation_vector[x][y] = real(binding / units::eV);
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
