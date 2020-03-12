#ifndef __INELASTIC_PENN_H_
#define __INELASTIC_PENN_H_

#include <functional>
#include "../core/particle.h"
#include "../common/util/table_1D.h"
#include "../common/util/table_2D.h"
#include "../common/util/table_3D.h"
#include "../common/util/range.h"

namespace nbl { namespace scatter {

/**
 * \brief Inelastic scattering, according to the full Penn model.
 *
 * We distinguish between inner-shell and Fermi-sea excitations. Inner-shells
 * are those with binding > 50 eV.
 *   - Fermi sea: same as Mao et al, doi:10.1063/1.3033564
 *   - Inner shells: same as Kieft et al, doi:10.1088/0022-3727/41/21/215310
 *                      (see also ::nbl::scatter::inelastic_thomas)
 *
 * \tparam gpu_flag                   Is the code to be run on a GPU?
 * \tparam opt_optical_phonon_loss    Assume optical phonon loss for energy loss less than band gap
 * \tparam opt_generate_secondary     Generate secondary electrons
 * \tparam opt_instantaneous_momentum Large losses: consider instantaneous momentum for SE
 * \tparam opt_momentum_conservation  Large losses: obey conservation of momentum
 */
template<bool gpu_flag,
	bool opt_optical_phonon_loss = true,
	bool opt_generate_secondary = true,
	bool opt_instantaneous_momentum = true,
	bool opt_momentum_conservation = true>
class inelastic_penn
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

		// "omega" is the energy lost by the primary. Represented as hbar*omega, units eV
		real omega;
		{
			const real x = logr(this_particle.kin_energy);
			const real y = rng.unit();
			omega = expr(_log_omega_icdf_table.get(x, y));
		}

		// Binding energy of the secondary electron
		real binding;
		{
			// We sample the ionisation table at (omega, P) instead of (K, P).
			// This is incorrect, it would be better to sample at (K, P) and
			// reject when binding > omega.
			// Physically, though, the shell with largest binding < omega is
			// expected to dominate the result. In addition, the ratio between
			// different shells stays fairly similar through energy. So this is
			// a good approximation and gives a significant speed advantage.
			const real x = logr(this_particle.kin_energy);
			const real y = logr(omega/this_particle.kin_energy);
			const real z = rng.unit();
			binding = _ionisation_table.get_rounddown(x, y, z);
		}

		// Normalise direction
		normalise(this_particle.dir);

		// No inner-shell binding energy: Fermi sea excitation
		if (binding <= 0)
		{
			if (omega > _band_gap)
			{
				// Electron excitation from the Fermi sea

				// "q" is represented as hbar*q / sqrt(2m), m electron mass; units eV^-1/2
				real q;
				{
					const real x = logr(this_particle.kin_energy);
					const real y = omega / this_particle.kin_energy;
					const real z = rng.unit();
					q = _q_icdf_table.get(x, y, z);
				}

				const vec3 prim_normal_dir = normalised(make_normal_vec(this_particle.dir, rng.phi()));
				vec3 q_dir, q_normal_dir;
				{
					const real costheta_prim_q = clampr((omega + q * q) / (2 * sqrtr(this_particle.kin_energy)*q), -1, 1);
					const real sintheta_prim_q = sqrtr(1 - costheta_prim_q * costheta_prim_q);

					q_dir = normalised(this_particle.dir*costheta_prim_q + prim_normal_dir * sintheta_prim_q);
					q_normal_dir = normalised(make_normal_vec(q_dir, rng.phi()));
				}


				// Deflect primary
				const real costheta_pi_pf = clampr((2 * this_particle.kin_energy - omega - q * q) / (2 * sqrtr(this_particle.kin_energy*(this_particle.kin_energy - omega))), -1, 1);
				const real sintheta_pi_pf = sqrtr(1 - costheta_pi_pf * costheta_pi_pf);
				this_particle.kin_energy -= omega;
				this_particle.dir = this_particle.dir*costheta_pi_pf + prim_normal_dir * sintheta_pi_pf;
				particle_mgr[particle_idx] = this_particle;


				// Create secondary
				if (opt_generate_secondary)
				{
					particle secondary_particle;
					secondary_particle.pos = this_particle.pos;

					const real k_fermi = sqrtr(_fermi);
					const real k_fermi_omega = sqrtr(_fermi + omega);
					if (q < k_fermi_omega - k_fermi)
					{
						// Plasmon decay.

						// Assume that the probability for the SE's energy before
						// excitation, Ei, is distributed according to the DOS
						// before multiplied the DOS after excitation; that the
						// bands are parabolic; and that the SE's direction after
						// excitation is uniformly distributed.

						// Define x = Ei/omega; f = fermi/omega.
						// We have p(Ei) ~ sqrt(Ei) * sqrt(Ei + omega), so we
						// need to sample from sqrt(x(x+1)) between x=0 and x=f.
						// Knowing sqrt(x^2+x) < sqrt(x^2+x+1/4) = x + 1/2, we
						// sample from p(x) ~ x + 1/2 and accept a fraction
						// (x+1/2 - sqrt(x^2+x)) / (x+1/2).
						// If f < 0.1, we sample from p(x) ~ sqrt(x) instead for
						// efficiency, and also to handle the case where fermi = 0.

						const real f = _fermi / omega;

						real x;
						if (f < 0.1)
						{
							const real U = cbrtr(rng.unit());
							x = f * U*U;
						}
						else
						{
							// Lower and upper bounds on x
							const real lower = maxr(0, f - 1);
							const real upper = f;

							// Utility variables
							const real d1 = upper - lower;
							const real d2 = upper * upper - lower * lower;
							const real alpha = d2 / (d1 + d2);

							do
							{
								// Sample x from p ~ x+1/2 between lower and upper
								const real U1 = rng.unit();
								const real U2 = rng.unit();
								x = (U1 < alpha)
									? sqrtr(d2*U2 + lower * lower) // x
									: d1 * U2 + lower;             // 1
							}
							// Reject
							while (rng.unit() > sqrtr(x*x + x) / (x + .5_r));
						}

						secondary_particle.kin_energy = x*omega + omega;
						secondary_particle.dir = rng.uniform_vector();
					}
					else if (q < k_fermi_omega + k_fermi)
					{
						// Single excitation.

						// We want to get the SE from the Fermi sphere. Its k
						// vector before excitation is ki. There is an annulus
						// of allowed initial k_i, centered around the q axis.
						// ki_parr is the value of ki in the direction parallel
						// to q, ki_perp is the value perpendicular to q.

						// Note that ki_parr < k_fermi
						const real ki_parr = (omega - q * q) / (2 * q);

						// Minimum ki_perp: because the SE's final energy must be above Fermi
						const real kf_parr = ki_parr + q;
						const real ki_perp_min_sq = maxr(0, _fermi - kf_parr * kf_parr);
						// Maximum ki_perp: because the SE's initial energy must be below Fermi
						const real ki_perp_max_sq = _fermi - ki_parr*ki_parr;

						// Sample a radius [probability distribution: p(r) = 2r / (max^2 - min^2)]
						const real ki_perp_sq = maxr(0, // For possible round-off errors
							ki_perp_min_sq + (ki_perp_max_sq - ki_perp_min_sq)*rng.unit());

						const vec3 ki = ki_parr*q_dir + sqrtr(ki_perp_sq)*q_normal_dir;

						secondary_particle.kin_energy = magnitude_squared(ki) + omega;
						secondary_particle.dir = ki + q*q_dir;
					}
					else
					{
						// q is too large for either plasmon decay or single
						// excitation: do not create a secondary electron.
						return;
					}

					particle_mgr.create_secondary(particle_idx, secondary_particle);
				}
				return;
			}
			else
			{
				// sub-band gap energy loss in semiconductors and insulators (see page 78 thesis T.V.)
				// energy loss due to longitudinal optical phonon excitation is assumed
				// update energy and EXIT
				if (opt_optical_phonon_loss)
				{
					this_particle.kin_energy -= omega;
					particle_mgr[particle_idx] = this_particle;
				}
				return;
			}
		}

		// Inner-shell excitation.

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
	 * \brief Create, given a legacy material file: this function immediately
	 *        throws an exception.
	 *
	 * The Penn model is not supported by legacy material files.
	 *
	 * \deprecated Old file format is deprecated and not supported by all
	 *             scattering mechanisms. This function will be removed soon.
	 */
	static CPU inelastic_penn create(material_legacy_thomas const & mat)
	{
		throw std::runtime_error("Cannot get Penn inelastic model from old-style .mat files");
	}

	/**
	 * \brief Create, given a material file.
	 */
	static CPU inelastic_penn create(hdf5_file const & mat)
	{
		util::geomspace<units::quantity<double>> K_range(K_min*units::eV, K_max*units::eV, K_cnt);
		util::linspace<units::quantity<double>> P_range(0*units::dimensionless, 1*units::dimensionless, P_cnt);

		inelastic_penn inel;

		inel._fermi = static_cast<real>(mat.get_property_quantity("fermi") / units::eV);
		inel._band_gap = static_cast<real>(mat.get_property_quantity("band_gap", -1*units::eV) / units::eV);

		inel._log_imfp_table = util::table_1D<real, gpu_flag>::create(logr(K_min), logr(K_max), K_cnt);
		inel._log_imfp_table.mem_scope([&](real* imfp_vector)
		{
			auto inelastic_imfp = mat.get_table_axes<1>("full_penn/imfp");
			for (int x = 0; x < K_cnt; ++x)
			{
				imfp_vector[x] = (real)std::log(inelastic_imfp.get_loglog(K_range[x]) * units::nm);
			}
		});

		inel._log_omega_icdf_table = util::table_2D<real, gpu_flag>::create(logr(K_min), logr(K_max), K_cnt, 0, 1, P_cnt);
		inel._log_omega_icdf_table.mem_scope([&](real** icdf_vector)
		{
			auto fermi = mat.get_property_quantity("fermi");
			auto inelastic_icdf = mat.get_table_axes<2>("full_penn/omega_icdf");
			for (int y = 0; y < P_cnt; ++y)
			{
				const units::quantity<double> P = P_range[y];
				for (int x = 0; x < K_cnt; ++x)
				{
					const units::quantity<double> K = K_range[x];

					icdf_vector[x][y] = (real)std::log(std::max(0.0, std::min<double>(
						(K - fermi) / units::eV,
						inelastic_icdf.get_linear(K, P) / units::eV
					)));
				}
			}
		});


		util::geomspace<units::quantity<double>> K512_range(K_min*units::eV, K_max*units::eV, 512);
		util::linspace<units::quantity<double>> P512_range(0*units::dimensionless, 1*units::dimensionless, 512);
		inel._q_icdf_table = util::table_3D<real, gpu_flag>::create(logr(K_min), logr(K_max), 512, 0, 1, 512, 0, 1, 512);
		inel._q_icdf_table.mem_scope([&](real*** icdf_vector)
		{
			auto inelastic_icdf = mat.get_table_axes<3>("full_penn/q_icdf");
			for (int z = 0; z < 512; ++z)
			{
				const units::quantity<double> P = P512_range[z];
				for (int y = 0; y < 512; ++y)
				{
					const units::quantity<double> Q = P512_range[y];
					for (int x = 0; x < 512; ++x)
					{
						units::quantity<double> K = K512_range[x];

						// hbar / sqrt(2*electron mass) == 0.19519 nm eV^1/2
						icdf_vector[x][y][z] = static_cast<real>(0.19519 *
							inelastic_icdf.get_linear(K, Q, P) * units::nm);
					}
				}
			}
		});


		inel._ionisation_table = util::table_3D<real, gpu_flag>::create(logr(K_min), logr(K_max), 128, logr(1e-4), logr(1), 1024, 0, 1, 1024);
		util::geomspace<units::quantity<double>> kk_range(K_min*units::eV, K_max*units::eV, 128);
		util::geomspace<units::quantity<double>> omega_range(1e-4*units::dimensionless, 1*units::dimensionless, 1024);
		util::linspace<units::quantity<double>> pp_range(0*units::dimensionless, 1*units::dimensionless, 1024);
		inel._ionisation_table.mem_scope([&](real*** ionisation_vector)
		{
			const auto icdf_table = mat.get_table_axes<3>("ionization/binding_icdf");

			// Create the simulation table
			for (int z = 0; z < 1024; ++z)
			{
				const units::quantity<double> P = pp_range[z];
				for (int y = 0; y < 1024; ++y)
				{
					const units::quantity<double> omega = omega_range[y];
					for (int x = 0; x < 128; ++x)
					{
						const units::quantity<double> K = kk_range[x];

						units::quantity<double> binding = icdf_table.get_rounddown(K, omega, P);
						if (binding < 50 * units::eV || !std::isfinite(binding.value))
							binding = -1 * units::eV;

						ionisation_vector[x][y][z] = real(binding / units::eV);
					}
				}
			}
		});

		return inel;
	}

	/**
	 * \brief Clone from another instance.
	 */
	template<bool source_gpu_flag>
	static CPU inelastic_penn create(inelastic_penn<source_gpu_flag, opt_optical_phonon_loss, opt_generate_secondary, opt_instantaneous_momentum, opt_momentum_conservation> const & source)
	{
		inelastic_penn target;

		target._fermi = source._fermi;
		target._band_gap = source._band_gap;

		target._log_imfp_table = util::table_1D<real, gpu_flag>::create(source._log_imfp_table);
		target._log_omega_icdf_table = util::table_2D<real, gpu_flag>::create(source._log_omega_icdf_table);
		target._q_icdf_table = util::table_3D<real, gpu_flag>::create(source._q_icdf_table);
		target._ionisation_table = util::table_3D<real, gpu_flag>::create(source._ionisation_table);

		return target;
	}

	/**
	 * \brief Dealllocate data held by an instance of this class.
	 */
	static CPU void destroy(inelastic_penn & inel)
	{
		util::table_1D<real, gpu_flag>::destroy(inel._log_imfp_table);
		util::table_2D<real, gpu_flag>::destroy(inel._log_omega_icdf_table);
		util::table_3D<real, gpu_flag>::destroy(inel._q_icdf_table);
		util::table_3D<real, gpu_flag>::destroy(inel._ionisation_table);
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
	 * \brief Table storing the probability distribution for energy loss.
	 *
	 * "ICDF" is short for Inverse Cumulative Distribution Function, which is
	 * what this table stores.
	 *
	 * Specifically, stores `log(omega / eV)` as function of
	 *   - x axis: `log(kinetic energy / eV)`
	 *   - y axis: cumulative probability (between 0 and 1)
	 */
	util::table_2D<real, gpu_flag> _log_omega_icdf_table;

	/**
	 * \brief Table storing the probability distribution for momentum transfer.
	 *
	 * "ICDF" is short for Inverse Cumulative Distribution Function, which is
	 * what this table stores.
	 *
	 * Specifically, stores `hbar*q/sqrt(2m) / eV^1/2` as function of
	 *   - x axis: `log(kinetic energy / eV)`
	 *   - y axis: omega / kinetic energy (between 0 and 1)
	 *   - y axis: cumulative probability (between 0 and 1)
	 */
	util::table_3D<real, gpu_flag> _q_icdf_table;

	/**
	 * \brief Inner-shell ionization table.
	 *
	 * Specifically, stores `binding energy / eV` as function of
	 *   - x axis: `log(kinetic energy / eV)`
	 *   - y axis: `log(omega / kinetic energy)`
	 *   - z axis: cumulative probability (between 0 and 1)
	 *
	 * "No binding energy" is denoted by -1 in the table.
	 *
	 * Be sure to use .get_rounddown() for non-interpolated values
	 */
	util::table_3D<real, gpu_flag> _ionisation_table;

	real _fermi;    ///< Fermi energy (eV)
	real _band_gap; ///< Band gap (eV) (no band gap is denoted by -1)

	template<bool, bool, bool, bool, bool>
	friend class inelastic_penn;
};

}} // namespace nbl::scatter

#endif // __INELASTIC_PENN_H_
