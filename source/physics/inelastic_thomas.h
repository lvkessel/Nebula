#ifndef __INELASTIC_THOMAS_H_
#define __INELASTIC_THOMAS_H_

#include "../core/particle.h"
#include "../common/util/table_1D.h"
#include "../common/util/table_2D.h"
#include "../common/util/range.h"
#include "electron_ionisation.h"

namespace nbl { namespace scatter {

/**
 * \brief Inelastic scattering, according to the Kieft & Bosch model.
 *
 *   - See doi:10.1088/0022-3727/41/21/215310 (Kieft paper)
 *   - See doi:10.4233/uuid:f214f594-a21f-4318-9f29-9776d60ab06c (Verduin thesis)
 *
 * \tparam gpu_flag                   Is the code to be run on a GPU?
 * \tparam opt_optical_phonon_loss    Assume optical phonon loss for energy loss less than band gap
 * \tparam opt_generate_secondary     Generate secondary electrons
 * \tparam opt_instantaneous_momentum Consider instantaneous momentum for SE
 * \tparam opt_momentum_conservation  Obey conservation of momentum
 */
template<bool gpu_flag,
	bool opt_optical_phonon_loss = true,
	bool opt_generate_secondary = true,
	bool opt_instantaneous_momentum = true,
	bool opt_momentum_conservation = true>
class inelastic_thomas
{
public:
	/**
	 * \brief Indicate when this class generates secondary electrons
	 */
	constexpr static bool may_create_se = opt_generate_secondary;

	/**
	 * \brief Sample a random free path length
	 */
	inline PHYSICS real sample_path(particle const & this_particle, util::random_generator<gpu_flag> & rng) const
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

		// draw a random zero-momentum energy loss of the primary electron
		real omega0;
		{// see thesis T.V. Eq. 3.82.
			const real x = logr(this_particle.kin_energy);
			const real y = rng.unit();
			omega0 = expr(_log_icdf_table.get(x, y));
		}

		// draw a random binding energy of the secondary electron.
		real binding = _binding.sample(omega0, rng);

		// draw a random total energy loss for the primary electron
		real omega;
		{// see thesis T.V. Eq. 3.85.
			real omega_max = 0.5_r*(this_particle.kin_energy + omega0 - _fermi); // upper limit of eq. 9 in Ashley, but corrected for the fermi energy
			real omega_min = omega0;
			real w0 = minr(omega0 - 1, maxr(0, binding) - _fermi);
			if (this_particle.kin_energy > 2*omega0)
			{
				// equation 10 in Ashley
				omega_min = 0.5_r*(this_particle.kin_energy + omega0
					- sqrtr(this_particle.kin_energy*(this_particle.kin_energy - 2*omega0)));
				w0 = omega0;
			}

			const real U = rng.unit();
			if ((w0 > 0) && (omega_min > w0) && (omega_min < omega_max)) {
				// For nonzero binding energy, sample omega according to equation 7 in Ashley,
				// using the lower and upper limits as defined above.
				// For inner-shell ionization (Ebind > 50 eV) we substitute the Fermi-energy corrected
				// binding energy for omegaprime (so that the differential cross section becomes inversely
				// proportional to both the total energy transfer and the kinetic energy of the secondary
				// electron).
				const real f_min = 1 / w0*logr((omega_min - w0) / omega_min);
				const real f_max = 1 / w0*logr((omega_max - w0) / omega_max);
				omega = -w0 / expm1r(w0*(f_min*(1 - U) + f_max*U));
			}
			else {
				// In some cases (typically only occuring for binding < 50 eV) we get omega_min > omega_max.
				// This is due to our Fermi energy correction in the definition of omega_max. Physically, this
				// means that momentum cannot be conserved because the primary electron cannot have a final
				// kinetic energy that is lower than the Fermi energy. In this (relatively rare) case we have
				// to ignore momentum conservation and probe omega according to a 1/(omega)^2 distribution
				// with omega0 and omega_max as lower and upper limits respectively.
				omega = omega0*omega_max / (omega0*(1 - U) + omega_max*U);
			}
		}

		// special cases if there is no binding energy:
		if (binding < 0) {
			if (_band_gap < 0) {
				// TODO for metals: excitation of a fermi sea electron
			}
			else if (omega0 > _band_gap) {
				// electron excitation across the band gap (see page 78 thesis T.V.)
				binding = _band_gap;
			}
			else {
				// sub-band gap energy loss in semiconductors and insulators (see page 78 thesis T.V.)
				// energy loss due to longitudinal optical phonon excitation is assumed
				// update energy and EXIT
				if (opt_optical_phonon_loss)
				{
					this_particle.kin_energy -= omega0;
					particle_mgr[particle_idx] = this_particle;
				}
				return;
			}
		}
		binding = maxr(0, binding);

		// Normalise direction
		normalise(this_particle.dir);

		// determine random normal vector to determine the scattering direction
		vec3 normal_dir = normalised(make_normal_vec(this_particle.dir, rng.phi()));

		// Determine the inelastic scattering angles.
		// We have strictly followed the method of Ivanchenko (and thus Kieft and Bosch)
		// See for more details thesis T.V. page 80-85.
		const real _K = this_particle.kin_energy - _fermi + 2*binding; // see thesis T.V. Eq. 3.105
		const real dK = binding + omega;        // see thesis T.V. Eq. 3.106
		const real cos_theta = sqrtr(dK / _K);  // see thesis T.V. Eq. 3.100
		const real sin_theta = sqrtr(1 - saturater(cos_theta*cos_theta));

		// Determine initial secondary direction (see thesis T.V. Eq. 3.107)
		// The initial direction is determined by assuming that the secondary electron
		// is at rest
		vec3 secondary_dir = this_particle.dir*cos_theta + normal_dir*sin_theta;

		// Add (optional) direction to account for the (intrinsic) instantaneous momentum
		//  of the secondary electron.
		// See thesis T.V. Eq. 3.108
		if (opt_instantaneous_momentum)
		{
			normalise(secondary_dir);
			secondary_dir += sqrtr(binding / dK) * rng.uniform_vector();
		}

		// ensure proper normalization of the secondary directional vector.
		normalise(secondary_dir);

		if (opt_generate_secondary)
		{
			particle secondary_particle;
			secondary_particle.kin_energy = _fermi + omega - binding; // See thesis T.V. Eq. 3.86
			secondary_particle.pos = this_particle.pos;
			secondary_particle.dir = secondary_dir;

			particle_mgr.create_secondary(particle_idx, secondary_particle);
		}

		this_particle.kin_energy -= omega;

		if (opt_momentum_conservation)
		{
			// primary direction determined by non-relativistic momentum-conservation, i.e.:
			//   sin(theta)*primary_dir_2 = primary_dir - cos(theta)*secondary_dir;
			// See thesis T.V. Eq. 3.111
			this_particle.dir -= cos_theta * secondary_dir;
		}

		// Store the scattered particle in memory
		particle_mgr[particle_idx] = this_particle;
	}

	/**
	 * \brief Create, given a legacy material file.
	 *
	 * \deprecated Old file format is deprecated and not supported by all
	 *             scattering mechanisms. This function will be removed soon.
	 */
	static CPU inelastic_thomas create(material_legacy_thomas const & mat)
	{
		inelastic_thomas inel;

		/*
		 * Set general data
		 */
		inel._fermi = static_cast<real>(mat.fermi() / constant::ec);
		inel._band_gap = static_cast<real>(mat.band_gap()() / constant::ec);
		inel._binding = electron_ionisation<gpu_flag>::create(mat);

		util::geomspace<double> K_range(K_min, K_max, K_cnt);
		util::linspace<double> P_range(0, 1, P_cnt);


		/*
		 * Set IMFP table
		 */
		inel._log_imfp_table = util::table_1D<real, gpu_flag>::create(logr(K_min), logr(K_max), K_cnt);
		inel._log_imfp_table.mem_scope([&](real* imfp_vector)
		{
			for (int x = 0; x < K_cnt; ++x)
			{
				const double K = K_range[x]*constant::ec; // in Joules
				imfp_vector[x] = (real)std::log(mat.density()*mat.inelastic_tcs(K)*1e-9);
			}
		});


		/*
		 * Set ICDF table
		 */
		inel._log_icdf_table = util::table_2D<real, gpu_flag>::create(logr(K_min), logr(K_max), K_cnt, 0, 1, P_cnt);
		inel._log_icdf_table.mem_scope([&](real** icdf_vector)
		{
			for (int y = 0; y < P_cnt; ++y)
			{
				const double P = P_range[y];
				for (int x = 0; x < K_cnt; ++x)
				{
					const double K = K_range[x]*constant::ec; // in Joules

					icdf_vector[y][x] = (real)std::log(std::max(0.0, std::min(
						(K - mat.fermi()) / constant::ec,
						mat.inelastic_icdf(K, P) / constant::ec
					)));
				}
			}
		});

		return inel;
	}

	/**
	 * \brief Create, given a material file.
	 */
	static CPU inelastic_thomas create(hdf5_file const & mat)
	{
		util::geomspace<units::quantity<double>> K_range(K_min*units::eV, K_max*units::eV, K_cnt);
		util::linspace<units::quantity<double>> P_range(0*units::dimensionless, 1*units::dimensionless, P_cnt);

		inelastic_thomas inel;

		inel._fermi = static_cast<real>(mat.get_property_quantity("fermi") / units::eV);
		inel._band_gap = static_cast<real>(mat.get_property_quantity("band_gap", -1*units::eV) / units::eV);
		inel._binding = electron_ionisation<gpu_flag>::create(mat);

		inel._log_imfp_table = util::table_1D<real, gpu_flag>::create(logr(K_min), logr(K_max), K_cnt);
		inel._log_imfp_table.mem_scope([&](real* imfp_vector)
		{
			auto inelastic_imfp = mat.get_table_axes<1>("inelastic_kieft/imfp");
			for (int x = 0; x < K_cnt; ++x)
			{
				imfp_vector[x] = (real)std::log(inelastic_imfp.get_loglog(K_range[x]) * units::nm);
			}
		});

		inel._log_icdf_table = util::table_2D<real, gpu_flag>::create(logr(K_min), logr(K_max), K_cnt, 0, 1, P_cnt);
		inel._log_icdf_table.mem_scope([&](real** icdf_vector)
		{
			auto fermi = mat.get_property_quantity("fermi");
			auto inelastic_icdf = mat.get_table_axes<2>("inelastic_kieft/w0_icdf");
			for (int y = 0; y < P_cnt; ++y)
			{
				const units::quantity<double> P = P_range[y];
				for (int x = 0; x < K_cnt; ++x)
				{
					units::quantity<double> K = K_range[x];

					// TODO: support creation of dimensionless quantities from scalars
					icdf_vector[y][x] = (real)std::log(std::max(0.0, std::min<double>(
						(K - fermi) / units::eV,
						inelastic_icdf.get_linear(K, P) / units::eV
					)));
				}
			}
		});

		return inel;
	}

	/**
	 * \brief Clone from another instance.
	 */
	template<bool source_gpu_flag>
	static CPU inelastic_thomas create(inelastic_thomas<source_gpu_flag, opt_optical_phonon_loss, opt_generate_secondary, opt_instantaneous_momentum, opt_momentum_conservation> const & source)
	{
		inelastic_thomas target;

		target._fermi = source._fermi;
		target._band_gap = source._band_gap;
		target._binding = electron_ionisation<gpu_flag>::create(source._binding);

		target._log_imfp_table = util::table_1D<real, gpu_flag>::create(source._log_imfp_table);
		target._log_icdf_table = util::table_2D<real, gpu_flag>::create(source._log_icdf_table);

		return target;
	}

	/**
	 * \brief Dealllocate data held by an instance of this class.
	 */
	static CPU void destroy(inelastic_thomas & inel)
	{
		util::table_1D<real, gpu_flag>::destroy(inel._log_imfp_table);
		util::table_2D<real, gpu_flag>::destroy(inel._log_icdf_table);
		electron_ionisation<gpu_flag>::destroy(inel._binding);
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
	 * \brief Table storing the probability distribution for "zero-momentum
	 * energy loss", omega prime.
	 *
	 * "ICDF" is short for Inverse Cumulative Distribution Function, which is
	 * what this table stores.
	 *
	 * Specifically, stores `log(omega' / eV)` as function of
	 *   - x axis: `log(kinetic energy / eV)`
	 *   - y axis: cumulative probability (between 0 and 1)
	 */
	util::table_2D<real, gpu_flag> _log_icdf_table;

	real _fermi;    ///< Fermi energy (eV)
	real _band_gap; ///< Band gap (eV)
	electron_ionisation<gpu_flag> _binding; ///< Table for inner and outer-shell binding energies

	template<bool, bool, bool, bool, bool>
	friend class inelastic_thomas;
};

}} // namespace nbl::scatter

#endif // __INELASTIC_THOMAS_H_
