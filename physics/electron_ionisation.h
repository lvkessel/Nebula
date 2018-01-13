#ifndef __ELECTRON_IONISATION_H_
#define __ELECTRON_IONISATION_H_

#include <functional>
#include "../common/util/table_2D.h"
#include "../common/util/random.h"
#include "../legacy_thomas/material.hh"
#include "../material/hdf5_file.h"

namespace nbl { namespace scatter {

template<bool gpu_flag>
class electron_ionisation
{
public:
	PHYSICS real sample(real K, util::random_generator<gpu_flag>& rng) const
	{
		const real x = logr(K);
		const real y = rng.unit();
		return _ionisation_table.get_rounddown(x, y);
	}

	static CPU electron_ionisation create(material_legacy_thomas const & mat)
	{
		/*
		 * TODO: move this somewhere else (see inelastic, elastic, binding)
		 * Translate index in a table to the relevant physical value.
		 * K_at refers to kinetic energy (log space) and P_at refers to
		 * the differential cross section (linear space).
		 */
		auto __logspace_K_at = [&](int x)
		{
			return K_min*std::exp(1.0*x / (K_cnt - 1)*std::log(K_max / K_min));
		};
		auto __linspace_P_at = [&](int y)
		{
			return 1.0*y / (P_cnt - 1);
		};

		electron_ionisation ei;
		ei._ionisation_table = util::table_2D<real, gpu_flag>::create(logr(K_min), logr(K_max), K_cnt, 0, 1, P_cnt);
		ei._ionisation_table.mem_scope([&](real** ionisation_vector)
		{
			for (int y = 0; y < P_cnt; ++y)
			{
				const double P = __linspace_P_at(y);
				for (int x = 0; x < K_cnt; ++x)
				{
					const double omega0 = __logspace_K_at(x); // in eV
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
		auto __logspace_K_at = [&](int x)
		{
			return K_min * std::exp(1.0*x / (K_cnt - 1)*std::log(K_max / K_min));
		};
		auto __linspace_P_at = [&](int y)
		{
			return 1.0*y / (P_cnt - 1);
		};

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
				const units::quantity<double> P = __linspace_P_at(y) * units::dimensionless;
				for (int x = 0; x < K_cnt; ++x)
				{
					const units::quantity<double> K = __logspace_K_at(x)*units::eV;
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

	static CPU void destroy(electron_ionisation & ei)
	{
		util::table_2D<real, gpu_flag>::destroy(ei._ionisation_table);
	}

private:
	util::table_2D<real, gpu_flag> _ionisation_table;
};

}} // namespace nbl::scatter

#endif // __ELECTRON_IONISATION_H_
