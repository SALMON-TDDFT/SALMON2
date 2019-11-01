!
!  Copyright 2019 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!
!-----------------------------------------------------------------------------------------
module phys_constants
  implicit none
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Important physical constants used in SALMON code
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    public :: au_fs      ! Atomic unit of time in femtosecond([fs]) unit.
    public :: au_aa      ! Atomic unit of length in Angstrom([Ang.]) unit.
    public :: au_nm      ! Atomic unit of length in [nm] unit.
    public :: au_ev      ! Atomic unit of energy in [eV] unit.
    public :: au_va      ! Atomic unit of electric field in [V/Ang.] unit
    public :: cspeed_au  ! Speed of light in atomic units
    public :: kB_au      ! Boltzman constant in atomic units
    public :: Debye_au   ! Debye electric dipole moment (1[D]) in atomic units

    public :: atomic_unit_of_1st_hyperpolarizability
    public :: atomic_unit_of_2nd_hyperpolarizability


  private
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Fundamental Physical Constants --- Complete Listing
    !! 2018 CODATA adjustment
    !!
    !! Followings are generated from http://physics.nist.gov/constants
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8), parameter :: alpha_particle_mass = 6.6446573357d-27 ! [kg] (+/- 2.0e-36)
    real(8), parameter :: alpha_particle_mass_energy_equivalent = 5.9719201914d-10 ! [J] (+/- 1.8e-19)
    real(8), parameter :: alpha_particle_mass_energy_equivalent_in_mev = 3727.3794066d0 ! [MeV] (+/- 1.1e-06)
    real(8), parameter :: alpha_particle_mass_in_u = 4.001506179127d0 ! [u] (+/- 6.3e-11)
    real(8), parameter :: alpha_particle_molar_mass = 4.0015061777d-3 ! [kg mol^-1] (+/- 1.2e-12)
    real(8), parameter :: alpha_particle_proton_mass_ratio = 3.97259969009d0 ! (+/- 2.2e-10)
    real(8), parameter :: alpha_particle_relative_atomic_mass = 4.001506179127d0 ! (+/- 6.3e-11)
    real(8), parameter :: angstrom_star = 1.00001495d-10 ! [m] (+/- 9.0e-17)
    real(8), parameter :: atomic_mass_constant = 1.66053906660d-27 ! [kg] (+/- 5.0e-37)
    real(8), parameter :: atomic_mass_constant_energy_equivalent = 1.49241808560d-10 ! [J] (+/- 4.5e-20)
    real(8), parameter :: atomic_mass_constant_energy_equivalent_in_mev = 931.49410242d0 ! [MeV] (+/- 2.8e-07)
    real(8), parameter :: atomic_mass_unit_electron_volt_relationship = 9.3149410242d8 ! [eV] (+/- 2.8e-01)
    real(8), parameter :: atomic_mass_unit_hartree_relationship = 3.4231776874d7 ! [E_h] (+/- 1.0e-02)
    real(8), parameter :: atomic_mass_unit_hertz_relationship = 2.25234271871d23 ! [Hz] (+/- 6.8e+13)
    real(8), parameter :: atomic_mass_unit_inverse_meter_relationship = 7.5130066104d14 ! [m^-1] (+/- 2.3e+05)
    real(8), parameter :: atomic_mass_unit_joule_relationship = 1.49241808560d-10 ! [J] (+/- 4.5e-20)
    real(8), parameter :: atomic_mass_unit_kelvin_relationship = 1.08095401916d13 ! [K] (+/- 3.3e+03)
    real(8), parameter :: atomic_mass_unit_kilogram_relationship = 1.66053906660d-27 ! [kg] (+/- 5.0e-37)
    real(8), parameter :: atomic_unit_of_1st_hyperpolarizability = 3.2063613061d-53 ! [C^3 m^3 J^-2] (+/- 1.5e-62)
    real(8), parameter :: atomic_unit_of_2nd_hyperpolarizability = 6.2353799905d-65 ! [C^4 m^4 J^-3] (+/- 3.8e-74)
    real(8), parameter :: atomic_unit_of_action = 1.054571817d-34 ! [J s]
    real(8), parameter :: atomic_unit_of_charge = 1.602176634d-19 ! [C]
    real(8), parameter :: atomic_unit_of_charge_density = 1.08120238457d12 ! [C m^-3] (+/- 4.9e+02)
    real(8), parameter :: atomic_unit_of_current = 6.623618237510d-3 ! [A] (+/- 1.3e-14)
    real(8), parameter :: atomic_unit_of_electric_dipole_mom = 8.4783536255d-30 ! [C m] (+/- 1.3e-39)
    real(8), parameter :: atomic_unit_of_electric_field = 5.14220674763d11 ! [V m^-1] (+/- 7.8e+01)
    real(8), parameter :: atomic_unit_of_electric_field_gradient = 9.7173624292d21 ! [V m^-2] (+/- 2.9e+12)
    real(8), parameter :: atomic_unit_of_electric_polarizability = 1.64877727436d-41 ! [C^2 m^2 J^-1] (+/- 5.0e-51)
    real(8), parameter :: atomic_unit_of_electric_potential = 27.211386245988d0 ! [V] (+/- 5.3e-11)
    real(8), parameter :: atomic_unit_of_electric_quadrupole_mom = 4.4865515246d-40 ! [C m^2] (+/- 1.4e-49)
    real(8), parameter :: atomic_unit_of_energy = 4.3597447222071d-18 ! [J] (+/- 8.5e-30)
    real(8), parameter :: atomic_unit_of_force = 8.2387234983d-8 ! [N] (+/- 1.2e-17)
    real(8), parameter :: atomic_unit_of_length = 5.29177210903d-11 ! [m] (+/- 8.0e-21)
    real(8), parameter :: atomic_unit_of_mag_dipole_mom = 1.85480201566d-23 ! [J T^-1] (+/- 5.6e-33)
    real(8), parameter :: atomic_unit_of_mag_flux_density = 2.35051756758d5 ! [T] (+/- 7.1e-05)
    real(8), parameter :: atomic_unit_of_magnetizability = 7.8910366008d-29 ! [J T^-2] (+/- 4.8e-38)
    real(8), parameter :: atomic_unit_of_mass = 9.1093837015d-31 ! [kg] (+/- 2.8e-40)
    real(8), parameter :: atomic_unit_of_momentum = 1.99285191410d-24 ! [kg m s^-1] (+/- 3.0e-34)
    real(8), parameter :: atomic_unit_of_permittivity = 1.11265005545d-10 ! [F m^-1] (+/- 1.7e-20)
    real(8), parameter :: atomic_unit_of_time = 2.4188843265857d-17 ! [s] (+/- 4.7e-29)
    real(8), parameter :: atomic_unit_of_velocity = 2.18769126364d6 ! [m s^-1] (+/- 3.3e-04)
    real(8), parameter :: avogadro_constant = 6.02214076d23 ! [mol^-1]
    real(8), parameter :: bohr_magneton = 9.2740100783d-24 ! [J T^-1] (+/- 2.8e-33)
    real(8), parameter :: bohr_magneton_in_ev_t = 5.7883818060d-5 ! [eV T^-1] (+/- 1.7e-14)
    real(8), parameter :: bohr_magneton_in_hz_t = 1.39962449361d10 ! [Hz T^-1] (+/- 4.2e+00)
    real(8), parameter :: bohr_magneton_in_inverse_meter_per_tesla = 46.686447783d0 ! [m^-1 T^-1] (+/- 1.4e-08)
    real(8), parameter :: bohr_magneton_in_k_t = 0.67171381563d0 ! [K T^-1] (+/- 2.0e-10)
    real(8), parameter :: bohr_radius = 5.29177210903d-11 ! [m] (+/- 8.0e-21)
    real(8), parameter :: boltzmann_constant = 1.380649d-23 ! [J K^-1]
    real(8), parameter :: boltzmann_constant_in_ev_k = 8.617333262d-5 ! [eV K^-1]
    real(8), parameter :: boltzmann_constant_in_hz_k = 2.083661912d10 ! [Hz K^-1]
    real(8), parameter :: boltzmann_constant_in_inverse_meter_per_kelvin = 69.50348004d0 ! [m^-1 K^-1]
    real(8), parameter :: classical_electron_radius = 2.8179403262d-15 ! [m] (+/- 1.3e-24)
    real(8), parameter :: compton_wavelength = 2.42631023867d-12 ! [m] (+/- 7.3e-22)
    real(8), parameter :: conductance_quantum = 7.748091729d-5 ! [S]
    real(8), parameter :: conventional_value_of_ampere_90 = 1.00000008887d0 ! [A]
    real(8), parameter :: conventional_value_of_coulomb_90 = 1.00000008887d0 ! [C]
    real(8), parameter :: conventional_value_of_farad_90 = 0.99999998220d0 ! [F]
    real(8), parameter :: conventional_value_of_henry_90 = 1.00000001779d0 ! [H]
    real(8), parameter :: conventional_value_of_josephson_constant = 483597.9d9 ! [Hz V^-1]
    real(8), parameter :: conventional_value_of_ohm_90 = 1.00000001779d0 ! [ohm]
    real(8), parameter :: conventional_value_of_volt_90 = 1.00000010666d0 ! [V]
    real(8), parameter :: conventional_value_of_von_klitzing_constant = 25812.807d0 ! [ohm]
    real(8), parameter :: conventional_value_of_watt_90 = 1.00000019553d0 ! [W]
    real(8), parameter :: cu_x_unit = 1.00207697d-13 ! [m] (+/- 2.8e-20)
    real(8), parameter :: deuteron_electron_mag_mom_ratio = -4.664345551d-4 ! (+/- 1.2e-12)
    real(8), parameter :: deuteron_electron_mass_ratio = 3670.48296788d0 ! (+/- 1.3e-07)
    real(8), parameter :: deuteron_g_factor = 0.8574382338d0 ! (+/- 2.2e-09)
    real(8), parameter :: deuteron_mag_mom = 4.330735094d-27 ! [J T^-1] (+/- 1.1e-35)
    real(8), parameter :: deuteron_mag_mom_to_bohr_magneton_ratio = 4.669754570d-4 ! (+/- 1.2e-12)
    real(8), parameter :: deuteron_mag_mom_to_nuclear_magneton_ratio = 0.8574382338d0 ! (+/- 2.2e-09)
    real(8), parameter :: deuteron_mass = 3.3435837724d-27 ! [kg] (+/- 1.0e-36)
    real(8), parameter :: deuteron_mass_energy_equivalent = 3.00506323102d-10 ! [J] (+/- 9.1e-20)
    real(8), parameter :: deuteron_mass_energy_equivalent_in_mev = 1875.61294257d0 ! [MeV] (+/- 5.7e-07)
    real(8), parameter :: deuteron_mass_in_u = 2.013553212745d0 ! [u] (+/- 4.0e-11)
    real(8), parameter :: deuteron_molar_mass = 2.01355321205d-3 ! [kg mol^-1] (+/- 6.1e-13)
    real(8), parameter :: deuteron_neutron_mag_mom_ratio = -0.44820653d0 ! (+/- 1.1e-07)
    real(8), parameter :: deuteron_proton_mag_mom_ratio = 0.30701220939d0 ! (+/- 7.9e-10)
    real(8), parameter :: deuteron_proton_mass_ratio = 1.99900750139d0 ! (+/- 1.1e-10)
    real(8), parameter :: deuteron_relative_atomic_mass = 2.013553212745d0 ! (+/- 4.0e-11)
    real(8), parameter :: deuteron_rms_charge_radius = 2.12799d-15 ! [m] (+/- 7.4e-19)
    real(8), parameter :: electron_charge_to_mass_quotient = -1.75882001076d11 ! [C kg^-1] (+/- 5.3e+01)
    real(8), parameter :: electron_deuteron_mag_mom_ratio = -2143.9234915d0 ! (+/- 5.6e-06)
    real(8), parameter :: electron_deuteron_mass_ratio = 2.724437107462d-4 ! (+/- 9.6e-15)
    real(8), parameter :: electron_g_factor = -2.00231930436256d0 ! (+/- 3.5e-13)
    real(8), parameter :: electron_gyromag_ratio = 1.76085963023d11 ! [s^-1 T^-1] (+/- 5.3e+01)
    real(8), parameter :: electron_gyromag_ratio_in_mhz_t = 28024.9514242d0 ! [MHz T^-1] (+/- 8.5e-06)
    real(8), parameter :: electron_helion_mass_ratio = 1.819543074573d-4 ! (+/- 7.9e-15)
    real(8), parameter :: electron_mag_mom = -9.2847647043d-24 ! [J T^-1] (+/- 2.8e-33)
    real(8), parameter :: electron_mag_mom_anomaly = 1.15965218128d-3 ! (+/- 1.8e-13)
    real(8), parameter :: electron_mag_mom_to_bohr_magneton_ratio = -1.00115965218128d0 ! (+/- 1.8e-13)
    real(8), parameter :: electron_mag_mom_to_nuclear_magneton_ratio = -1838.28197188d0 ! (+/- 1.1e-07)
    real(8), parameter :: electron_mass = 9.1093837015d-31 ! [kg] (+/- 2.8e-40)
    real(8), parameter :: electron_mass_energy_equivalent = 8.1871057769d-14 ! [J] (+/- 2.5e-23)
    real(8), parameter :: electron_mass_energy_equivalent_in_mev = 0.51099895000d0 ! [MeV] (+/- 1.5e-10)
    real(8), parameter :: electron_mass_in_u = 5.48579909065d-4 ! [u] (+/- 1.6e-14)
    real(8), parameter :: electron_molar_mass = 5.4857990888d-7 ! [kg mol^-1] (+/- 1.7e-16)
    real(8), parameter :: electron_muon_mag_mom_ratio = 206.7669883d0 ! (+/- 4.6e-06)
    real(8), parameter :: electron_muon_mass_ratio = 4.83633169d-3 ! (+/- 1.1e-10)
    real(8), parameter :: electron_neutron_mag_mom_ratio = 960.92050d0 ! (+/- 2.3e-04)
    real(8), parameter :: electron_neutron_mass_ratio = 5.4386734424d-4 ! (+/- 2.6e-13)
    real(8), parameter :: electron_proton_mag_mom_ratio = -658.21068789d0 ! (+/- 2.0e-07)
    real(8), parameter :: electron_proton_mass_ratio = 5.44617021487d-4 ! (+/- 3.3e-14)
    real(8), parameter :: electron_relative_atomic_mass = 5.48579909065d-4 ! (+/- 1.6e-14)
    real(8), parameter :: electron_tau_mass_ratio = 2.87585d-4 ! (+/- 1.9e-08)
    real(8), parameter :: electron_to_alpha_particle_mass_ratio = 1.370933554787d-4 ! (+/- 4.5e-15)
    real(8), parameter :: electron_to_shielded_helion_mag_mom_ratio = 864.058257d0 ! (+/- 1.0e-05)
    real(8), parameter :: electron_to_shielded_proton_mag_mom_ratio = -658.2275971d0 ! (+/- 7.2e-06)
    real(8), parameter :: electron_triton_mass_ratio = 1.819200062251d-4 ! (+/- 9.0e-15)
    real(8), parameter :: electron_volt = 1.602176634d-19 ! [J]
    real(8), parameter :: electron_volt_atomic_mass_unit_relationship = 1.07354410233d-9 ! [u] (+/- 3.2e-19)
    real(8), parameter :: electron_volt_hartree_relationship = 3.6749322175655d-2 ! [E_h] (+/- 7.1e-14)
    real(8), parameter :: electron_volt_hertz_relationship = 2.417989242d14 ! [Hz]
    real(8), parameter :: electron_volt_inverse_meter_relationship = 8.065543937d5 ! [m^-1]
    real(8), parameter :: electron_volt_joule_relationship = 1.602176634d-19 ! [J]
    real(8), parameter :: electron_volt_kelvin_relationship = 1.160451812d4 ! [K]
    real(8), parameter :: electron_volt_kilogram_relationship = 1.782661921d-36 ! [kg]
    real(8), parameter :: elementary_charge = 1.602176634d-19 ! [C]
    real(8), parameter :: elementary_charge_over_h_bar = 1.519267447d15 ! [A J^-1]
    real(8), parameter :: faraday_constant = 96485.33212d0 ! [C mol^-1]
    real(8), parameter :: fermi_coupling_constant = 1.1663787d-5 ! [GeV^-2] (+/- 6.0e-12)
    real(8), parameter :: fine_structure_constant = 7.2973525693d-3 ! (+/- 1.1e-12)
    real(8), parameter :: first_radiation_constant = 3.741771852d-16 ! [W m^2]
    real(8), parameter :: first_radiation_constant_for_spectral_radiance = 1.191042972d-16 ! [W m^2 sr^-1]
    real(8), parameter :: hartree_atomic_mass_unit_relationship = 2.92126232205d-8 ! [u] (+/- 8.8e-18)
    real(8), parameter :: hartree_electron_volt_relationship = 27.211386245988d0 ! [eV] (+/- 5.3e-11)
    real(8), parameter :: hartree_energy = 4.3597447222071d-18 ! [J] (+/- 8.5e-30)
    real(8), parameter :: hartree_energy_in_ev = 27.211386245988d0 ! [eV] (+/- 5.3e-11)
    real(8), parameter :: hartree_hertz_relationship = 6.579683920502d15 ! [Hz] (+/- 1.3e+04)
    real(8), parameter :: hartree_inverse_meter_relationship = 2.1947463136320d7 ! [m^-1] (+/- 4.3e-05)
    real(8), parameter :: hartree_joule_relationship = 4.3597447222071d-18 ! [J] (+/- 8.5e-30)
    real(8), parameter :: hartree_kelvin_relationship = 3.1577502480407d5 ! [K] (+/- 6.1e-07)
    real(8), parameter :: hartree_kilogram_relationship = 4.8508702095432d-35 ! [kg] (+/- 9.4e-47)
    real(8), parameter :: helion_electron_mass_ratio = 5495.88528007d0 ! (+/- 2.4e-07)
    real(8), parameter :: helion_g_factor = -4.255250615d0 ! (+/- 5.0e-08)
    real(8), parameter :: helion_mag_mom = -1.074617532d-26 ! [J T^-1] (+/- 1.3e-34)
    real(8), parameter :: helion_mag_mom_to_bohr_magneton_ratio = -1.158740958d-3 ! (+/- 1.4e-11)
    real(8), parameter :: helion_mag_mom_to_nuclear_magneton_ratio = -2.127625307d0 ! (+/- 2.5e-08)
    real(8), parameter :: helion_mass = 5.0064127796d-27 ! [kg] (+/- 1.5e-36)
    real(8), parameter :: helion_mass_energy_equivalent = 4.4995394125d-10 ! [J] (+/- 1.4e-19)
    real(8), parameter :: helion_mass_energy_equivalent_in_mev = 2808.39160743d0 ! [MeV] (+/- 8.5e-07)
    real(8), parameter :: helion_mass_in_u = 3.014932247175d0 ! [u] (+/- 9.7e-11)
    real(8), parameter :: helion_molar_mass = 3.01493224613d-3 ! [kg mol^-1] (+/- 9.1e-13)
    real(8), parameter :: helion_proton_mass_ratio = 2.99315267167d0 ! (+/- 1.3e-10)
    real(8), parameter :: helion_relative_atomic_mass = 3.014932247175d0 ! (+/- 9.7e-11)
    real(8), parameter :: helion_shielding_shift = 5.996743d-5 ! (+/- 1.0e-10)
    real(8), parameter :: hertz_atomic_mass_unit_relationship = 4.4398216652d-24 ! [u] (+/- 1.3e-33)
    real(8), parameter :: hertz_electron_volt_relationship = 4.135667696d-15 ! [eV]
    real(8), parameter :: hertz_hartree_relationship = 1.5198298460570d-16 ! [E_h] (+/- 2.9e-28)
    real(8), parameter :: hertz_inverse_meter_relationship = 3.335640951d-9 ! [m^-1]
    real(8), parameter :: hertz_joule_relationship = 6.62607015d-34 ! [J]
    real(8), parameter :: hertz_kelvin_relationship = 4.799243073d-11 ! [K]
    real(8), parameter :: hertz_kilogram_relationship = 7.372497323d-51 ! [kg]
    real(8), parameter :: hyperfine_transition_frequency_of_cs_133 = 9192631770d0 ! [Hz]
    real(8), parameter :: inverse_fine_structure_constant = 137.035999084d0 ! (+/- 2.1e-08)
    real(8), parameter :: inverse_meter_atomic_mass_unit_relationship = 1.33102505010d-15 ! [u] (+/- 4.0e-25)
    real(8), parameter :: inverse_meter_electron_volt_relationship = 1.239841984d-6 ! [eV]
    real(8), parameter :: inverse_meter_hartree_relationship = 4.5563352529120d-8 ! [E_h] (+/- 8.8e-20)
    real(8), parameter :: inverse_meter_hertz_relationship = 299792458d0 ! [Hz]
    real(8), parameter :: inverse_meter_joule_relationship = 1.986445857d-25 ! [J]
    real(8), parameter :: inverse_meter_kelvin_relationship = 1.438776877d-2 ! [K]
    real(8), parameter :: inverse_meter_kilogram_relationship = 2.210219094d-42 ! [kg]
    real(8), parameter :: inverse_of_conductance_quantum = 12906.40372d0 ! [ohm]
    real(8), parameter :: josephson_constant = 483597.8484d9 ! [Hz V^-1]
    real(8), parameter :: joule_atomic_mass_unit_relationship = 6.7005352565d9 ! [u] (+/- 2.0e+00)
    real(8), parameter :: joule_electron_volt_relationship = 6.241509074d18 ! [eV]
    real(8), parameter :: joule_hartree_relationship = 2.2937122783963d17 ! [E_h] (+/- 4.5e+05)
    real(8), parameter :: joule_hertz_relationship = 1.509190179d33 ! [Hz]
    real(8), parameter :: joule_inverse_meter_relationship = 5.034116567d24 ! [m^-1]
    real(8), parameter :: joule_kelvin_relationship = 7.242970516d22 ! [K]
    real(8), parameter :: joule_kilogram_relationship = 1.112650056d-17 ! [kg]
    real(8), parameter :: kelvin_atomic_mass_unit_relationship = 9.2510873014d-14 ! [u] (+/- 2.8e-23)
    real(8), parameter :: kelvin_electron_volt_relationship = 8.617333262d-5 ! [eV]
    real(8), parameter :: kelvin_hartree_relationship = 3.1668115634556d-6 ! [E_h] (+/- 6.1e-18)
    real(8), parameter :: kelvin_hertz_relationship = 2.083661912d10 ! [Hz]
    real(8), parameter :: kelvin_inverse_meter_relationship = 69.50348004d0 ! [m^-1]
    real(8), parameter :: kelvin_joule_relationship = 1.380649d-23 ! [J]
    real(8), parameter :: kelvin_kilogram_relationship = 1.536179187d-40 ! [kg]
    real(8), parameter :: kilogram_atomic_mass_unit_relationship = 6.0221407621d26 ! [u] (+/- 1.8e+17)
    real(8), parameter :: kilogram_electron_volt_relationship = 5.609588603d35 ! [eV]
    real(8), parameter :: kilogram_hartree_relationship = 2.0614857887409d34 ! [E_h] (+/- 4.0e+22)
    real(8), parameter :: kilogram_hertz_relationship = 1.356392489d50 ! [Hz]
    real(8), parameter :: kilogram_inverse_meter_relationship = 4.524438335d41 ! [m^-1]
    real(8), parameter :: kilogram_joule_relationship = 8.987551787d16 ! [J]
    real(8), parameter :: kilogram_kelvin_relationship = 6.509657260d39 ! [K]
    real(8), parameter :: lattice_parameter_of_silicon = 5.431020511d-10 ! [m] (+/- 8.9e-18)
    real(8), parameter :: lattice_spacing_of_ideal_si_220 = 1.920155716d-10 ! [m] (+/- 3.2e-18)
    real(8), parameter :: loschmidt_constant_273_15_k_100_kpa = 2.651645804d25 ! [m^-3]
    real(8), parameter :: loschmidt_constant_273_15_k_101_325_kpa = 2.686780111d25 ! [m^-3]
    real(8), parameter :: luminous_efficacy = 683d0 ! [lm W^-1]
    real(8), parameter :: mag_flux_quantum = 2.067833848d-15 ! [Wb]
    real(8), parameter :: molar_gas_constant = 8.314462618d0 ! [J mol^-1 K^-1]
    real(8), parameter :: molar_mass_constant = 0.99999999965d-3 ! [kg mol^-1] (+/- 3.0e-13)
    real(8), parameter :: molar_mass_of_carbon_12 = 11.9999999958d-3 ! [kg mol^-1] (+/- 3.6e-12)
    real(8), parameter :: molar_planck_constant = 3.990312712d-10 ! [J Hz^-1 mol^-1]
    real(8), parameter :: molar_volume_of_ideal_gas_273_15_k_100_kpa = 22.71095464d-3 ! [m^3 mol^-1]
    real(8), parameter :: molar_volume_of_ideal_gas_273_15_k_101_325_kpa = 22.41396954d-3 ! [m^3 mol^-1]
    real(8), parameter :: molar_volume_of_silicon = 1.205883199d-5 ! [m^3 mol^-1] (+/- 6.0e-13)
    real(8), parameter :: mo_x_unit = 1.00209952d-13 ! [m] (+/- 5.3e-20)
    real(8), parameter :: muon_compton_wavelength = 1.173444110d-14 ! [m] (+/- 2.6e-22)
    real(8), parameter :: muon_electron_mass_ratio = 206.7682830d0 ! (+/- 4.6e-06)
    real(8), parameter :: muon_g_factor = -2.0023318418d0 ! (+/- 1.3e-09)
    real(8), parameter :: muon_mag_mom = -4.49044830d-26 ! [J T^-1] (+/- 1.0e-33)
    real(8), parameter :: muon_mag_mom_anomaly = 1.16592089d-3 ! (+/- 6.3e-10)
    real(8), parameter :: muon_mag_mom_to_bohr_magneton_ratio = -4.84197047d-3 ! (+/- 1.1e-10)
    real(8), parameter :: muon_mag_mom_to_nuclear_magneton_ratio = -8.89059703d0 ! (+/- 2.0e-07)
    real(8), parameter :: muon_mass = 1.883531627d-28 ! [kg] (+/- 4.2e-36)
    real(8), parameter :: muon_mass_energy_equivalent = 1.692833804d-11 ! [J] (+/- 3.8e-19)
    real(8), parameter :: muon_mass_energy_equivalent_in_mev = 105.6583755d0 ! [MeV] (+/- 2.3e-06)
    real(8), parameter :: muon_mass_in_u = 0.1134289259d0 ! [u] (+/- 2.5e-09)
    real(8), parameter :: muon_molar_mass = 1.134289259d-4 ! [kg mol^-1] (+/- 2.5e-12)
    real(8), parameter :: muon_neutron_mass_ratio = 0.1124545170d0 ! (+/- 2.5e-09)
    real(8), parameter :: muon_proton_mag_mom_ratio = -3.183345142d0 ! (+/- 7.1e-08)
    real(8), parameter :: muon_proton_mass_ratio = 0.1126095264d0 ! (+/- 2.5e-09)
    real(8), parameter :: muon_tau_mass_ratio = 5.94635d-2 ! (+/- 4.0e-06)
    real(8), parameter :: natural_unit_of_action = 1.054571817d-34 ! [J s]
    real(8), parameter :: natural_unit_of_action_in_ev_s = 6.582119569d-16 ! [eV s]
    real(8), parameter :: natural_unit_of_energy = 8.1871057769d-14 ! [J] (+/- 2.5e-23)
    real(8), parameter :: natural_unit_of_energy_in_mev = 0.51099895000d0 ! [MeV] (+/- 1.5e-10)
    real(8), parameter :: natural_unit_of_length = 3.8615926796d-13 ! [m] (+/- 1.2e-22)
    real(8), parameter :: natural_unit_of_mass = 9.1093837015d-31 ! [kg] (+/- 2.8e-40)
    real(8), parameter :: natural_unit_of_momentum = 2.73092453075d-22 ! [kg m s^-1] (+/- 8.2e-32)
    real(8), parameter :: natural_unit_of_momentum_in_mev_c = 0.51099895000d0 ! [MeV/c] (+/- 1.5e-10)
    real(8), parameter :: natural_unit_of_time = 1.28808866819d-21 ! [s] (+/- 3.9e-31)
    real(8), parameter :: natural_unit_of_velocity = 299792458d0 ! [m s^-1]
    real(8), parameter :: neutron_compton_wavelength = 1.31959090581d-15 ! [m] (+/- 7.5e-25)
    real(8), parameter :: neutron_electron_mag_mom_ratio = 1.04066882d-3 ! (+/- 2.5e-10)
    real(8), parameter :: neutron_electron_mass_ratio = 1838.68366173d0 ! (+/- 8.9e-07)
    real(8), parameter :: neutron_g_factor = -3.82608545d0 ! (+/- 9.0e-07)
    real(8), parameter :: neutron_gyromag_ratio = 1.83247171d8 ! [s^-1 T^-1] (+/- 4.3e+01)
    real(8), parameter :: neutron_gyromag_ratio_in_mhz_t = 29.1646931d0 ! [MHz T^-1] (+/- 6.9e-06)
    real(8), parameter :: neutron_mag_mom = -9.6623651d-27 ! [J T^-1] (+/- 2.3e-33)
    real(8), parameter :: neutron_mag_mom_to_bohr_magneton_ratio = -1.04187563d-3 ! (+/- 2.5e-10)
    real(8), parameter :: neutron_mag_mom_to_nuclear_magneton_ratio = -1.91304273d0 ! (+/- 4.5e-07)
    real(8), parameter :: neutron_mass = 1.67492749804d-27 ! [kg] (+/- 9.5e-37)
    real(8), parameter :: neutron_mass_energy_equivalent = 1.50534976287d-10 ! [J] (+/- 8.6e-20)
    real(8), parameter :: neutron_mass_energy_equivalent_in_mev = 939.56542052d0 ! [MeV] (+/- 5.4e-07)
    real(8), parameter :: neutron_mass_in_u = 1.00866491595d0 ! [u] (+/- 4.9e-10)
    real(8), parameter :: neutron_molar_mass = 1.00866491560d-3 ! [kg mol^-1] (+/- 5.7e-13)
    real(8), parameter :: neutron_muon_mass_ratio = 8.89248406d0 ! (+/- 2.0e-07)
    real(8), parameter :: neutron_proton_mag_mom_ratio = -0.68497934d0 ! (+/- 1.6e-07)
    real(8), parameter :: neutron_proton_mass_difference = 2.30557435d-30 ! [kg] (+/- 8.2e-37)
    real(8), parameter :: neutron_proton_mass_difference_energy_equivalent = 2.07214689d-13 ! [J] (+/- 7.4e-20)
    real(8), parameter :: neutron_proton_mass_difference_energy_equivalent_in_mev = 1.29333236d0 ! [MeV] (+/- 4.6e-07)
    real(8), parameter :: neutron_proton_mass_difference_in_u = 1.38844933d-3 ! [u] (+/- 4.9e-10)
    real(8), parameter :: neutron_proton_mass_ratio = 1.00137841931d0 ! (+/- 4.9e-10)
    real(8), parameter :: neutron_relative_atomic_mass = 1.00866491595d0 ! (+/- 4.9e-10)
    real(8), parameter :: neutron_tau_mass_ratio = 0.528779d0 ! (+/- 3.6e-05)
    real(8), parameter :: neutron_to_shielded_proton_mag_mom_ratio = -0.68499694d0 ! (+/- 1.6e-07)
    real(8), parameter :: newtonian_constant_of_gravitation = 6.67430d-11 ! [m^3 kg^-1 s^-2] (+/- 1.5e-15)
    real(8), parameter :: newtonian_constant_of_gravitation_over_h_bar_c = 6.70883d-39 ! [(GeV/c^2)^-2] (+/- 1.5e-43)
    real(8), parameter :: nuclear_magneton = 5.0507837461d-27 ! [J T^-1] (+/- 1.5e-36)
    real(8), parameter :: nuclear_magneton_in_ev_t = 3.15245125844d-8 ! [eV T^-1] (+/- 9.6e-18)
    real(8), parameter :: nuclear_magneton_in_inverse_meter_per_tesla = 2.54262341353d-2 ! [m^-1 T^-1] (+/- 7.8e-12)
    real(8), parameter :: nuclear_magneton_in_k_t = 3.6582677756d-4 ! [K T^-1] (+/- 1.1e-13)
    real(8), parameter :: nuclear_magneton_in_mhz_t = 7.6225932291d0 ! [MHz T^-1] (+/- 2.3e-09)
    real(8), parameter :: planck_constant = 6.62607015d-34 ! [J Hz^-1]
    real(8), parameter :: planck_constant_in_ev_hz = 4.135667696d-15 ! [eV Hz^-1]
    real(8), parameter :: planck_length = 1.616255d-35 ! [m] (+/- 1.8e-40)
    real(8), parameter :: planck_mass = 2.176434d-8 ! [kg] (+/- 2.4e-13)
    real(8), parameter :: planck_mass_energy_equivalent_in_gev = 1.220890d19 ! [GeV] (+/- 1.4e+14)
    real(8), parameter :: planck_temperature = 1.416784d32 ! [K] (+/- 1.6e+27)
    real(8), parameter :: planck_time = 5.391247d-44 ! [s] (+/- 6.0e-49)
    real(8), parameter :: proton_charge_to_mass_quotient = 9.5788331560d7 ! [C kg^-1] (+/- 2.9e-02)
    real(8), parameter :: proton_compton_wavelength = 1.32140985539d-15 ! [m] (+/- 4.0e-25)
    real(8), parameter :: proton_electron_mass_ratio = 1836.15267343d0 ! (+/- 1.1e-07)
    real(8), parameter :: proton_g_factor = 5.5856946893d0 ! (+/- 1.6e-09)
    real(8), parameter :: proton_gyromag_ratio = 2.6752218744d8 ! [s^-1 T^-1] (+/- 1.1e-01)
    real(8), parameter :: proton_gyromag_ratio_in_mhz_t = 42.577478518d0 ! [MHz T^-1] (+/- 1.8e-08)
    real(8), parameter :: proton_mag_mom = 1.41060679736d-26 ! [J T^-1] (+/- 6.0e-36)
    real(8), parameter :: proton_mag_mom_to_bohr_magneton_ratio = 1.52103220230d-3 ! (+/- 4.6e-13)
    real(8), parameter :: proton_mag_mom_to_nuclear_magneton_ratio = 2.79284734463d0 ! (+/- 8.2e-10)
    real(8), parameter :: proton_mag_shielding_correction = 2.5689d-5 ! (+/- 1.1e-08)
    real(8), parameter :: proton_mass = 1.67262192369d-27 ! [kg] (+/- 5.1e-37)
    real(8), parameter :: proton_mass_energy_equivalent = 1.50327761598d-10 ! [J] (+/- 4.6e-20)
    real(8), parameter :: proton_mass_energy_equivalent_in_mev = 938.27208816d0 ! [MeV] (+/- 2.9e-07)
    real(8), parameter :: proton_mass_in_u = 1.007276466621d0 ! [u] (+/- 5.3e-11)
    real(8), parameter :: proton_molar_mass = 1.00727646627d-3 ! [kg mol^-1] (+/- 3.1e-13)
    real(8), parameter :: proton_muon_mass_ratio = 8.88024337d0 ! (+/- 2.0e-07)
    real(8), parameter :: proton_neutron_mag_mom_ratio = -1.45989805d0 ! (+/- 3.4e-07)
    real(8), parameter :: proton_neutron_mass_ratio = 0.99862347812d0 ! (+/- 4.9e-10)
    real(8), parameter :: proton_relative_atomic_mass = 1.007276466621d0 ! (+/- 5.3e-11)
    real(8), parameter :: proton_rms_charge_radius = 8.414d-16 ! [m] (+/- 1.9e-18)
    real(8), parameter :: proton_tau_mass_ratio = 0.528051d0 ! (+/- 3.6e-05)
    real(8), parameter :: quantum_of_circulation = 3.6369475516d-4 ! [m^2 s^-1] (+/- 1.1e-13)
    real(8), parameter :: quantum_of_circulation_times_2 = 7.2738951032d-4 ! [m^2 s^-1] (+/- 2.2e-13)
    real(8), parameter :: reduced_compton_wavelength = 3.8615926796d-13 ! [m] (+/- 1.2e-22)
    real(8), parameter :: reduced_muon_compton_wavelength = 1.867594306d-15 ! [m] (+/- 4.2e-23)
    real(8), parameter :: reduced_neutron_compton_wavelength = 2.1001941552d-16 ! [m] (+/- 1.2e-25)
    real(8), parameter :: reduced_planck_constant = 1.054571817d-34 ! [J s]
    real(8), parameter :: reduced_planck_constant_in_ev_s = 6.582119569d-16 ! [eV s]
    real(8), parameter :: reduced_planck_constant_times_c_in_mev_fm = 197.3269804d0 ! [MeV fm]
    real(8), parameter :: reduced_proton_compton_wavelength = 2.10308910336d-16 ! [m] (+/- 6.4e-26)
    real(8), parameter :: reduced_tau_compton_wavelength = 1.110538d-16 ! [m] (+/- 7.5e-21)
    real(8), parameter :: rydberg_constant = 10973731.568160d0 ! [m^-1] (+/- 2.1e-05)
    real(8), parameter :: rydberg_constant_times_c_in_hz = 3.2898419602508d15 ! [Hz] (+/- 6.4e+03)
    real(8), parameter :: rydberg_constant_times_hc_in_ev = 13.605693122994d0 ! [eV] (+/- 2.6e-11)
    real(8), parameter :: rydberg_constant_times_hc_in_j = 2.1798723611035d-18 ! [J] (+/- 4.2e-30)
    real(8), parameter :: sackur_tetrode_constant_1_k_100_kpa = -1.15170753706d0 ! (+/- 4.5e-10)
    real(8), parameter :: sackur_tetrode_constant_1_k_101_325_kpa = -1.16487052358d0 ! (+/- 4.5e-10)
    real(8), parameter :: second_radiation_constant = 1.438776877d-2 ! [m K]
    real(8), parameter :: shielded_helion_gyromag_ratio = 2.037894569d8 ! [s^-1 T^-1] (+/- 2.4e+00)
    real(8), parameter :: shielded_helion_gyromag_ratio_in_mhz_t = 32.43409942d0 ! [MHz T^-1] (+/- 3.8e-07)
    real(8), parameter :: shielded_helion_mag_mom = -1.074553090d-26 ! [J T^-1] (+/- 1.3e-34)
    real(8), parameter :: shielded_helion_mag_mom_to_bohr_magneton_ratio = -1.158671471d-3 ! (+/- 1.4e-11)
    real(8), parameter :: shielded_helion_mag_mom_to_nuclear_magneton_ratio = -2.127497719d0 ! (+/- 2.5e-08)
    real(8), parameter :: shielded_helion_to_proton_mag_mom_ratio = -0.7617665618d0 ! (+/- 8.9e-09)
    real(8), parameter :: shielded_helion_to_shielded_proton_mag_mom_ratio = -0.7617861313d0 ! (+/- 3.3e-09)
    real(8), parameter :: shielded_proton_gyromag_ratio = 2.675153151d8 ! [s^-1 T^-1] (+/- 2.9e+00)
    real(8), parameter :: shielded_proton_gyromag_ratio_in_mhz_t = 42.57638474d0 ! [MHz T^-1] (+/- 4.6e-07)
    real(8), parameter :: shielded_proton_mag_mom = 1.410570560d-26 ! [J T^-1] (+/- 1.5e-34)
    real(8), parameter :: shielded_proton_mag_mom_to_bohr_magneton_ratio = 1.520993128d-3 ! (+/- 1.7e-11)
    real(8), parameter :: shielded_proton_mag_mom_to_nuclear_magneton_ratio = 2.792775599d0 ! (+/- 3.0e-08)
    real(8), parameter :: shielding_difference_of_d_and_p_in_hd = 2.0200d-8 ! (+/- 2.0e-11)
    real(8), parameter :: shielding_difference_of_t_and_p_in_ht = 2.4140d-8 ! (+/- 2.0e-11)
    real(8), parameter :: speed_of_light_in_vacuum = 299792458d0 ! [m s^-1]
    real(8), parameter :: standard_acceleration_of_gravity = 9.80665d0 ! [m s^-2]
    real(8), parameter :: standard_atmosphere = 101325d0 ! [Pa]
    real(8), parameter :: standard_state_pressure = 100000d0 ! [Pa]
    real(8), parameter :: stefan_boltzmann_constant = 5.670374419d-8 ! [W m^-2 K^-4]
    real(8), parameter :: tau_compton_wavelength = 6.97771d-16 ! [m] (+/- 4.7e-20)
    real(8), parameter :: tau_electron_mass_ratio = 3477.23d0 ! (+/- 2.3e-01)
    real(8), parameter :: tau_energy_equivalent = 1776.86d0 ! [MeV] (+/- 1.2e-01)
    real(8), parameter :: tau_mass = 3.16754d-27 ! [kg] (+/- 2.1e-31)
    real(8), parameter :: tau_mass_energy_equivalent = 2.84684d-10 ! [J] (+/- 1.9e-14)
    real(8), parameter :: tau_mass_in_u = 1.90754d0 ! [u] (+/- 1.3e-04)
    real(8), parameter :: tau_molar_mass = 1.90754d-3 ! [kg mol^-1] (+/- 1.3e-07)
    real(8), parameter :: tau_muon_mass_ratio = 16.8170d0 ! (+/- 1.1e-03)
    real(8), parameter :: tau_neutron_mass_ratio = 1.89115d0 ! (+/- 1.3e-04)
    real(8), parameter :: tau_proton_mass_ratio = 1.89376d0 ! (+/- 1.3e-04)
    real(8), parameter :: thomson_cross_section = 6.6524587321d-29 ! [m^2] (+/- 6.0e-38)
    real(8), parameter :: triton_electron_mass_ratio = 5496.92153573d0 ! (+/- 2.7e-07)
    real(8), parameter :: triton_g_factor = 5.957924931d0 ! (+/- 1.2e-08)
    real(8), parameter :: triton_mag_mom = 1.5046095202d-26 ! [J T^-1] (+/- 3.0e-35)
    real(8), parameter :: triton_mag_mom_to_bohr_magneton_ratio = 1.6223936651d-3 ! (+/- 3.2e-12)
    real(8), parameter :: triton_mag_mom_to_nuclear_magneton_ratio = 2.9789624656d0 ! (+/- 5.9e-09)
    real(8), parameter :: triton_mass = 5.0073567446d-27 ! [kg] (+/- 1.5e-36)
    real(8), parameter :: triton_mass_energy_equivalent = 4.5003878060d-10 ! [J] (+/- 1.4e-19)
    real(8), parameter :: triton_mass_energy_equivalent_in_mev = 2808.92113298d0 ! [MeV] (+/- 8.5e-07)
    real(8), parameter :: triton_mass_in_u = 3.01550071621d0 ! [u] (+/- 1.2e-10)
    real(8), parameter :: triton_molar_mass = 3.01550071517d-3 ! [kg mol^-1] (+/- 9.2e-13)
    real(8), parameter :: triton_proton_mass_ratio = 2.99371703414d0 ! (+/- 1.5e-10)
    real(8), parameter :: triton_relative_atomic_mass = 3.01550071621d0 ! (+/- 1.2e-10)
    real(8), parameter :: triton_to_proton_mag_mom_ratio = 1.0666399191d0 ! (+/- 2.1e-09)
    real(8), parameter :: unified_atomic_mass_unit = 1.66053906660d-27 ! [kg] (+/- 5.0e-37)
    real(8), parameter :: vacuum_electric_permittivity = 8.8541878128d-12 ! [F m^-1] (+/- 1.3e-21)
    real(8), parameter :: vacuum_mag_permeability = 1.25663706212d-6 ! [N A^-2] (+/- 1.9e-16)
    real(8), parameter :: von_klitzing_constant = 25812.80745d0 ! [ohm]
    real(8), parameter :: weak_mixing_angle = 0.22290d0 ! (+/- 3.0e-04)
    real(8), parameter :: wien_frequency_displacement_law_constant = 5.878925757d10 ! [Hz K^-1]
    real(8), parameter :: wien_wavelength_displacement_law_constant = 2.897771955d-3 ! [m K]
    real(8), parameter :: w_to_z_mass_ratio = 0.88153d0 ! (+/- 1.7e-04)

    !! Physical constants used in SALMON program:
    real(8), parameter :: au_fs = atomic_unit_of_time * 1d+15
    real(8), parameter :: au_aa = atomic_unit_of_length * 1d+10
    real(8), parameter :: au_nm = atomic_unit_of_length * 1d+9
    real(8), parameter :: au_va = atomic_unit_of_electric_field * 1d-10
    real(8), parameter :: au_ev = atomic_unit_of_energy / electron_volt
    real(8), parameter :: cspeed_au = inverse_fine_structure_constant
    real(8), parameter :: kB_au = boltzmann_constant / atomic_unit_of_energy

    real(8), parameter :: debye_au = (1d-21 / speed_of_light_in_vacuum) / (atomic_unit_of_charge * atomic_unit_of_length)
end module phys_constants

