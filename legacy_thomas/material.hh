/**
 * @file src/common/material.hh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef eSCATTER__CDSEM__MATERIAL__HEADER_INCLUDED
#define eSCATTER__CDSEM__MATERIAL__HEADER_INCLUDED

#include <map>
#include <string>
#include <vector>

#include "archive.hh"
#include "optional.hh"
#include "constant.hh"

class material_legacy_thomas {
friend archive::ostream& operator<<(archive::ostream&, const material_legacy_thomas&);
friend archive::istream& operator>>(archive::istream&, material_legacy_thomas&);
public:
	material_legacy_thomas() = default;
	material_legacy_thomas(const std::string& name, double fermi, double barrier, double phononloss, double density);
    /*!
     * @param[in] name
     *  Human readable identification string for the material.
     * @param[in] fermi
     *  Fermi energy.
     * @param[in] barrier
     *  Minimum energy required for an electron to escape from the material.
     * @param[in] bandgap
     *  Energy gap between the valence band and the conduction band.
     * @param[in] phononloss
     *  Energy loss in a elastic collision caused by phonon interaction.
     * @param[in] density
     *  The number density of the material.
     */
	material_legacy_thomas(const std::string& name, double fermi, double barrier, double phononloss, double bandgap, double density);
    inline const std::string& name() const;
    inline double fermi() const;
    inline double barrier() const;
    inline const optional<double>& band_gap() const;
    inline double phonon_loss() const;
    inline double density() const;
    inline double elastic_tcs(double K) const;
    inline double elastic_icdf(double K, double P) const;
    inline double inelastic_tcs(double K) const;
    inline double inelastic_icdf(double K, double P) const;
    double ionization_energy(double K, double P) const;
    double outer_shell_ionization_energy(double omega0) const;
	material_legacy_thomas& set_elastic_data(double K, const std::map<double,double>& dcs_map);
	material_legacy_thomas& set_inelastic_data(double K, const std::map<double,double>& dcs_map);
	material_legacy_thomas& set_ionization_data(double B, const std::map<double,double>& tcs_map);
	material_legacy_thomas& set_outer_shell_ionization_data(const std::vector<double>& osi_vector);
//#warning "members of material class are public!"
public:
    std::string _name;
    double _fermi = 0;
    double _barrier = 0;
    optional<double> _band_gap;
    double _phonon_loss = 0;
    double _density = 0;
    std::map<double,double> _elastic_tcs;
    std::map<double,std::map<double,double>> _elastic_icdf;
    std::map<double,double> _inelastic_tcs;
    std::map<double,std::map<double,double>> _inelastic_icdf;
    std::map<double,std::map<double,double>> _ionization_tcs;
    std::vector<double> _osi_energies;
};

archive::ostream& operator<<(archive::ostream&, const material_legacy_thomas&);
archive::istream& operator>>(archive::istream&, material_legacy_thomas&);

#include "material.inc"

#endif
