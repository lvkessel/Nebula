#ifndef __ELASTIC_THOMAS_H_
#define __ELASTIC_THOMAS_H_

#include "../common/util/table_1D.h"
#include "../common/util/table_2D.h"
#include "../material/hdf5_file.h"

namespace nbl { namespace scatter {

/**
 * \brief Elastic scattering, as described in Thomas Verduin's thesis.
 *
 * This is a combination of Mott scattering (for energy > 200 eV), acoustic
 * phonon scattering (< 100 eV) and interpolation in between. The cross sections
 * are combined by the cross section tool.
 *
 * T.V.'s thesis: doi:10.4233/uuid:f214f594-a21f-4318-9f29-9776d60ab06c
 *
 * \tparam gpu_flag                 Is the code to be run on a GPU?
 * \tparam opt_acoustic_phonon_loss Lose energy to acoustic phonons
 * \tparam opt_atomic_recoil_loss   Lose energy to atomic recoil for Mott scattering.
 */
template<bool gpu_flag,
	bool opt_acoustic_phonon_loss = true,
	bool opt_atomic_recoil_loss = true>
class elastic_thomas
{
public:
	/**
	 * \brief Indicate that this class never generates secondary electrons.
	 */
	constexpr static bool may_create_se = false;

	/**
	 * \brief Sample a random free path length
	 */
	template<typename particle_t>
	inline PHYSICS real sample_path(particle_t const & this_particle, util::random_generator<gpu_flag> & rng) const
	{
		// Get inverse mean free path for this kinetic energy
		const real imfp = expr(_log_imfp_table.get(logr(this_particle.kin_energy)));

		// Draw a distance
		return rng.exponential(1 / imfp);
	}

	/**
	 * \brief Perform a scattering event
	 */
	template<typename particle_manager>
	inline PHYSICS void execute(
		particle_manager& particle_mgr,
		typename particle_manager::particle_index_t particle_idx,
		util::random_generator<gpu_flag>& rng) const
	{
		// Retrieve current particle from global memory
		auto this_particle = particle_mgr[particle_idx];

		real cos_theta, sin_theta;
		{// draw a random elastic scatter angle by interpolating tables
			const float x = logr(particle_mgr[particle_idx].kin_energy);
			const float y = rng.unit();
			cos_theta = clampr(_icdf_table.get(x, y), -1, 1);
			sin_theta = sqrtr(1 - cos_theta*cos_theta);
		}

		// normalize current direction
		this_particle.dir = normalised(this_particle.dir);

		// find a random normal vector to the current direction of flight and normalize
		vec3 normal_dir = normalised(make_normal_vec(this_particle.dir, rng.phi()));

		// determine the scattered direction
		this_particle.dir = this_particle.dir * cos_theta + normal_dir * sin_theta;

		// special cases for phonon scattering and atom recoil energy loss.
		// the energy domain for phonon scattering is exactly the same as in the
		//   original Kieft & Bosch code.
		// the amount of energy loss for phonons can be found in the thesis of T.V. Eq. 3.116.
		// TODO: set variable domain for phonon scattering
		if (this_particle.kin_energy < 200)
		{
			if (opt_acoustic_phonon_loss)
				this_particle.kin_energy -= minr(50e-3f, _phonon_loss);
		}
		else
		{
			// account for atomic recoil (only added for compliance with Kieft & Bosch code)
			// There is no reference for this formula, can only be found in the Kieft & Bosch code.
			if (opt_atomic_recoil_loss)
			{
				this_particle.kin_energy -= _recoil_const*(1 - cos_theta) * this_particle.kin_energy;
			}
		}

		// Store the scattered particle in memory
		particle_mgr[particle_idx] = this_particle;
	}

	/**
	 * \brief Create, given a legacy material file
	 *
	 * \deprecated Old file format is deprecated and not supported by all
	 *             scattering mechanisms. This function will be removed soon.
	 */
	static CPU elastic_thomas create(material_legacy_thomas const & mat)
	{
		elastic_thomas el;

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


		/*
		 * Set IMFP table
		 */
		el._log_imfp_table = util::table_1D<real, gpu_flag>::create(logr(K_min), logr(K_max), K_cnt);
		el._log_imfp_table.mem_scope([&](real* imfp_vector)
		{
			for (int x = 0; x < K_cnt; ++x)
			{
				const double K = __logspace_K_at(x)*constant::ec; // in Joules
				imfp_vector[x] = (real)std::log(mat.density()*mat.elastic_tcs(K)*1e-9);
			}
		});


		/*
		 * Set ICDF table
		 */
		el._icdf_table = util::table_2D<real, gpu_flag>::create(logr(K_min), logr(K_max), K_cnt, 0, 1, P_cnt);
		el._icdf_table.mem_scope([&](real** icdf_vector)
		{
			for (int y = 0; y < P_cnt; ++y)
			{
				const double P = __linspace_P_at(y);
				for (int x = 0; x < K_cnt; ++x)
				{
					const double K = __logspace_K_at(x)*constant::ec; // in Joules

					icdf_vector[y][x] = (real)std::cos(std::max(0.0, std::min((double)pi, mat.elastic_icdf(K, P))));
				}
			}
		});

		el._phonon_loss = mat.phonon_loss()/constant::ec;
		el._recoil_const = real(2 * 9.109282e-28 * 6.02214086e+23 / 28.1);
		//                      2 * el mass      * Avogadro const / silicon A
		// Hard-coded A for silicon: not in the material file. This is also what e-scatter does.

		return el;
	}

	/**
	 * \brief Create, given a material file
	 */
	static CPU elastic_thomas create(hdf5_file const & mat)
	{
		auto __logspace_K_at = [&](int x)
		{
			return K_min * std::exp(1.0*x / (K_cnt - 1)*std::log(K_max / K_min));
		};
		auto __linspace_P_at = [&](int y)
		{
			return 1.0*y / (P_cnt - 1);
		};

		elastic_thomas el;

		el._log_imfp_table = util::table_1D<real, gpu_flag>::create(logr(K_min), logr(K_max), K_cnt);
		el._log_imfp_table.mem_scope([&](real* imfp_vector)
		{
			auto elastic_imfp = mat.get_table_axes<1>("elastic/imfp");
			for (int x = 0; x < K_cnt; ++x)
			{
				imfp_vector[x] = (real)std::log(elastic_imfp.get_loglog(__logspace_K_at(x) * units::eV) * units::nm);
			}
		});

		el._icdf_table = util::table_2D<real, gpu_flag>::create(logr(K_min), logr(K_max), K_cnt, 0, 1, P_cnt);
		el._icdf_table.mem_scope([&](real** icdf_vector)
		{
			auto elastic_icdf = mat.get_table_axes<2>("elastic/costheta_icdf");
			for (int y = 0; y < P_cnt; ++y)
			{
				const double P = __linspace_P_at(y);
				for (int x = 0; x < K_cnt; ++x)
				{
					// TODO: support creation of dimensionless quantities from scalars
					icdf_vector[y][x] = (real)std::max(-1.0, std::min<double>(1.0, elastic_icdf.get_linear(__logspace_K_at(x)*units::eV, P*units::dimensionless)));
				}
			}
		});

		el._phonon_loss = static_cast<real>(mat.get_property_quantity("phonon_loss") / units::eV);
		el._recoil_const = 2 * static_cast<real>(
			2 * (9.109383e-28 * units::g) // 2 * (electron mass)
			/ mat.get_property_quantity("effective_A"));

		return el;
	}

	/**
	 * \brief Dealllocate data held by an instance of this class.
	 */
	static CPU void destroy(elastic_thomas & el)
	{
		util::table_1D<real, gpu_flag>::destroy(el._log_imfp_table);
		util::table_2D<real, gpu_flag>::destroy(el._icdf_table);
	}

private:
	/**
	 * \brief Table storing the inverse mean free path as function of energy.
	 *
	 * Actually, stores `log(inverse mean free path / nm^-1)` as function of
	 * `log(kinetic energy / eV)`.
	 */
	util::table_1D<real, gpu_flag> _log_imfp_table;

	/**
	 * \brief Table storing the probability distribution for the scattering angle.
	 *
	 * "ICDF" is short for Inverse Cumulative Distribution Function, which is
	 * what this table stores.
	 *
	 * Specifically, stores `cos(theta)` as function of
	 *   - x axis: `log(kinetic energy / eV)`
	 *   - y axis: cumulative probability (between 0 and 1)
	 */
	util::table_2D<real, gpu_flag> _icdf_table;

	real _phonon_loss;  ///< Amount of energy lost in a phonon event (eV)
	real _recoil_const; ///< Amount of energy lost in a Mott event (eV)
};

}} // namespace nbl::scatter

#endif // __ELASTIC_THOMAS_H_
