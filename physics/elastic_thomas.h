#ifndef __ELASTIC_THOMAS_H_
#define __ELASTIC_THOMAS_H_

#include "../common/util/table_1D.h"
#include "../common/util/table_2D.h"

/*
 * TODO: acoustic phonon loss
 * TODO: atomic recoil loss
 */

namespace nbl { namespace scatter {

template<bool gpu_flag>
class elastic_thomas
{
public:
	constexpr static bool may_create_se = false;

	template<typename particle_t>
	inline PHYSICS real sample_path(particle_t const & this_particle, util::random_generator<gpu_flag> & rng) const
	{
		// Get inverse mean free path for this kinetic energy
		const real imfp = expr(_log_imfp_table.get(logr(this_particle.kin_energy)));

		// Draw a distance
		return rng.exponential(1 / imfp);
	}

	template<typename particle_manager>
	inline PHYSICS void execute(
		particle_manager& particle_mgr,
		typename particle_manager::particle_index_t particle_idx,
		util::random_generator<gpu_flag>& rng) const
	{
		// Retrieve current particle from global memory
		//auto this_particle = particle_mgr[particle_idx];

		real cos_theta, sin_theta;
		{// draw a random elastic scatter angle by interpolating tables
			const float x = logr(particle_mgr[particle_idx].kin_energy);
			const float y = rng.unit();
			cos_theta = clampr(_icdf_table.get(x, y), -1, 1);
			sin_theta = sqrtr(1 - cos_theta*cos_theta);
		}

		// normalize current direction
		const auto this_particle_dir = normalised(particle_mgr[particle_idx].dir);

		// find a random normal vector to the current direction of flight and normalize
		vec3 normal_dir = normalised(make_normal_vec(this_particle_dir, rng.phi()));

		// determine the scattered direction
		particle_mgr[particle_idx].dir = this_particle_dir * cos_theta + normal_dir * sin_theta;

/*		// special cases for phonon scattering and atom recoil energy loss.
		// the energy domain for phonon scattering is exactly the same as in the
		//   original Kieft & Bosch code.
		// the amount of energy loss for phonons can be found in the thesis of T.V. Eq. 3.116.
		// TODO: set variable domain for phonon scattering
		if (K < 200) {
			if (opt.acoustic_phonon_loss_flag)
				pstruct.K_energy_dev_p[particle_idx] = K - fminf(50e-3f, mstruct.phonon_loss_dev_p[material_idx]);
		}
		else {
			// account for atomic recoil (only added for compliance with Kieft & Bosch code)
			// There is no reference for this formula, can only be found in the Kieft & Bosch code.
			// Default behaviour does not include this effect.
			if (opt.atomic_recoil_loss_flag) {
				//            #warning "fixed energy loss due to atom recoil assumes silicon material"
				pstruct.K_energy_dev_p[particle_idx] = K - 2.0f*(_me*_NA)*K*(1.0f - cos_theta) / 28.1f;
			}
		}*/
	}

	static HOST elastic_thomas create(material_legacy_thomas const & mat)
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
		std::vector<real> imfp_vector(K_cnt, 0);
		for (int x = 0; x < K_cnt; ++x)
		{
			const double K = __logspace_K_at(x)*constant::ec; // in Joules
			imfp_vector[x] = (real)std::log(mat.density()*mat.elastic_tcs(K)*1e-9);
		}
		el._log_imfp_table = util::table_1D<real, gpu_flag>::create(imfp_vector.data(), logr(K_min), logr(K_max), K_cnt);


		/*
		 * Set ICDF table
		 */
		std::vector<real> icdf_vector(K_cnt*P_cnt, 0);
		for (int y = 0; y < P_cnt; ++y)
		{
			const double P = __linspace_P_at(y);
			for (int x = 0; x < K_cnt; ++x)
			{
				const double K = __logspace_K_at(x)*constant::ec; // in Joules

				icdf_vector[y*K_cnt + x] = (real)std::cos(std::max(0.0, std::min((double)pi, mat.elastic_icdf(K, P))));
			}
		}
		el._icdf_table = util::table_2D<real, gpu_flag>::create(icdf_vector.data(), logr(K_min), logr(K_max), K_cnt, 0, 1, P_cnt);


		return el;
	}

	static HOST void destroy(elastic_thomas & el)
	{
		util::table_1D<real, gpu_flag>::destroy(el._log_imfp_table);
		util::table_2D<real, gpu_flag>::destroy(el._icdf_table);
	}

private:
	// log(inverse mean free path / nm^-1) as function of log(kinetic energy / eV).
	util::table_1D<real, gpu_flag> _log_imfp_table;

	// cos(theta) as function of
	//   - x axis: log(kinetic energy / eV)
	//   - y axis: cum. probability
	util::table_2D<real, gpu_flag> _icdf_table;
};

}} // namespace nbl::scatter

#endif // __ELASTIC_THOMAS_H_