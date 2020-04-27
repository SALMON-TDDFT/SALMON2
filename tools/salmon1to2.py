#!/usr/bin/env python
#
#   Copyright 2019-2020 SALMON developers
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
#

USAGE = '''
SALMON v.1 to v.2 input file converter:
This script helps to translates the old SALMON input files to
present input file formats.

Specify the input/output file by options:
    %prog -i old_inputfile.inp -o new_inputfile.inp
Or use stdin as below:
    %prog < old_inputfile.inp 
'''

keyword_table = [
    # Inputfile group and parameter table:
    ['calculation', [
        ['theory',               ('calculation', 'theory')],
        ['calc_mode',            ('calculation', 'calc_mode')],
        ['use_ehrenfest_md',     ('calculation', 'use_ehrenfest_md')],
        ['use_adiabatic_md',     ('calculation', 'use_adiabatic_md')],
        ['use_ms_maxwell',       ('calculation', 'use_ms_maxwell')],
        ['use_geometry_opt',     ('calculation', 'use_geometry_opt')],
        ['use_potential_model',  ('calculation', 'use_potential_model')],
        ['use_singlescale',      ('calculation', 'use_singlescale')],
    ]],
    ['control', [
        ['restart_option',       ('control', 'restart_option')],
        ['backup_frequency',     ('control', 'backup_frequency')],
        ['time_shutdown',        ('control', 'time_shutdown')],
        ['sysname',              ('control', 'sysname')],
        ['directory',            ('control', 'directory')],
        ['dump_filename',        ('control', 'dump_filename')],
        ['modify_gs_wfn_k',      ('control', 'modify_gs_wfn_k')],
        ['read_gs_wfn_k',        ('control', 'read_gs_wfn_k')],
        ['read_gs_dns_cube',     ('control', 'read_gs_dns_cube')],
        ['read_rt_wfn_k',        ('control', 'read_rt_wfn_k')],
        ['write_gs_wfn_k',       ('control', 'write_gs_wfn_k')],
        ['write_rt_wfn_k',       ('control', 'write_rt_wfn_k')],
        ['read_gs_wfn_k_ms',     ('control', 'read_gs_wfn_k_ms')],
        ['read_rt_wfn_k_ms',     ('control', 'read_rt_wfn_k_ms')],
        ['write_gs_wfn_k_ms',    ('control', 'write_gs_wfn_k_ms')],
        ['write_rt_wfn_k_ms',    ('control', 'write_rt_wfn_k_ms')],
    ]],
    ['units', [
        ['unit_system',          ('units', 'unit_system')],
    ]],
    ['parallel', [
        ['nproc_k',              ('parallel', 'nproc_k')],
        ['nproc_ob',             ('parallel', 'nproc_ob')],
        ['nproc_domain_orbital', ('parallel', 'nproc_domain')],
        ['nproc_domain_general', ('parallel', 'nproc_domain_s')],
        ['num_datafiles_in',     ('parallel', 'num_datafiles_in')],
        ['num_datafiles_out',    ('parallel', 'num_datafiles_out')],
        ['yn_ffte',              None],  # Special rule (fourier)
        ['process_allocation',   None],
    ]],
    ['system', [
        ['yn_periodic',          None],  # Special rule (iperiodic)
        ['ispin',                ('system', 'ispin')],
        ['al',                   ('system', 'al')],
        ['al_vec1',              ('system', 'al_vec1')],
        ['al_vec2',              ('system', 'al_vec2')],
        ['al_vec3',              ('system', 'al_vec3')],
        ['isym',                 ('system', 'isym')],
        ['crystal_structure',    ('system', 'crystal_structure')],
        ['nstate',               ('system', 'nstate')],
        ['nstate_spin',          ('system', 'nstate_spin')],
        ['nelec',                ('system', 'nelec')],
        ['nelec_spin',           ('system', 'nelec_spin')],
        ['temperature',          ('system', 'temperature')],
        ['temperature_k',        ('system', 'temperature_k')],
        ['nelem',                ('system', 'nelem')],
        ['natom',                ('system', 'natom')],
        ['file_atom_coor',       ('system', 'file_atom_coor')],
        ['file_atom_red_coor',   ('system', 'file_atom_red_coor')],
    ]],
    ['pseudo', [
        ['file_pseudo',          ('pseudo', 'pseudo_file')],
        ['lmax_ps',              ('pseudo', 'lmax_ps')],
        ['lloc_ps',              ('pseudo', 'lloc_ps')],
        ['izatom',               ('pseudo', 'izatom')],
        ['yn_psmask',            ('pseudo', 'psmask_option')],
        ['alpha_mask',           ('pseudo', 'alpha_mask')],
        ['gamma_mask',           ('pseudo', 'gamma_mask')],
        ['eta_mask',             ('pseudo', 'eta_mask')],
    ]],
    ['functional', [
        ['xc',                   ('functional', 'xc')],
        ['cname',                ('functional', 'cname')],
        ['xname',                ('functional', 'xname')],
        ['alibx',                ('functional', 'alibx')],
        ['alibc',                ('functional', 'alibc')],
        ['alibxc',               ('functional', 'alibxc')],
        ['cval',                 ('functional', 'cval')],
    ]],
    ['rgrid', [
        ['dl',                   ('rgrid', 'dl')],
        ['num_rgrid',            ('rgrid', 'num_rgrid')],
    ]],
    ['kgrid', [
        ['num_kgrid',            ('kgrid', 'num_kgrid')],
        ['file_kw',              ('kgrid', 'file_kw')],
    ]],
    ['tgrid', [
        ['nt',                   ('tgrid', 'nt')],
        ['dt',                   ('tgrid', 'dt')],
    ]],
    ['propagation', [
        ['n_hamil',              ('propagation', 'n_hamil')],
        ['propagator',           ('propagation', 'propagator')],
        ['yn_fix_func',          ('functional', 'no_update_func')],
    ]],
    ['scf', [
        ['method_min',           ('scf', 'amin_routine')],
        ['ncg',                  ('scf', 'ncg')],
        ['method_mixing',        ('scf', 'amixing')],
        ['mixrate',              ('scf', 'rmixrate')],
        ['nmemory_mb',           ('scf', 'nmemory_mb')],
        ['alpha_mb',             ('scf', 'alpha_mb')],
        ['fsset_option',         ('scf', 'fsset_option')],
        ['nfsset_start',         ('scf', 'nfsset_start')],
        ['nfsset_every',         ('scf', 'nfsset_every')],
        ['nscf',                 ('scf', 'nscf')],
        ['yn_subspace_diagonalization', ('scf', 'subspace_diagonalization')],
        ['convergence',          ('scf', 'convergence')],
        ['threshold',            ('scf', 'threshold')],
        ['omp_loop',             ('scf', 'omp_loop')],
        ['skip_gsortho',         ('scf', 'skip_gsortho')],
        ['iditer_notemperature', ('scf', 'iditer_notemperature')],
    ]],
    ['emfield', [
        ['trans_longi',          ('emfield', 'trans_longi')],
        ['ae_shape1',            ('emfield', 'ae_shape1')],
        ['e_impulse',            ('emfield', 'e_impulse')],
        ['e_amplitude1',         ('emfield', 'amplitude1')],
        ['i_wcm2_1',             ('emfield', 'rlaser_int_wcm2_1')],
        ['tw1',                  ('emfield', 'pulse_tw1')],
        ['omega1',               ('emfield', 'omega1')],
        ['epdir_re1',            ('emfield', 'epdir_re1')],
        ['epdir_im1',            ('emfield', 'epdir_im1')],
        ['phi_cep1',             ('emfield', 'phi_cep1')],
        ['ae_shape2',            ('emfield', 'ae_shape2')],
        ['e_amplitude2',         ('emfield', 'amplitude2')],
        ['i_wcm2_2',             ('emfield', 'rlaser_int_wcm2_2')],
        ['tw2',                  ('emfield', 'pulse_tw2')],
        ['omega2',               ('emfield', 'omega2')],
        ['epdir_re2',            ('emfield', 'epdir_re2')],
        ['epdir_im2',            ('emfield', 'epdir_im2')],
        ['phi_cep2',             ('emfield', 'phi_cep2')],
        ['t1_t2',                ('emfield', 't1_t2')],
        ['t1_start',             ('emfield', 't1_delay')],
        ['yn_local_field',       ('emfield', 'alocal_laser')],
        ['rlaserbound_sta',      ('emfield', 'rlaserbound_sta')],
        ['rlaserbound_end',      ('emfield', 'rlaserbound_end')],
        ['num_dipole_source',    ('emfield', 'nump')],
        ['vec_dipole_source',    ('emfield', 'vecp')],
        ['cood_dipole_source',   ('emfield', 'coop')],
        ['rad_dipole_source',    ('emfield', 'radp_diele')],
    ]],
    ['multiscale', [
        ['fdtddim',              ('multiscale', 'fdtddim')],
        ['twod_shape',           ('multiscale', 'twod_shape')],
        ['nx_m',                 ('multiscale', 'nx_m')],
        ['ny_m',                 ('multiscale', 'ny_m')],
        ['nz_m',                 ('multiscale', 'nz_m')],
        ['hx_m',                 ('multiscale', 'hx_m')],
        ['hy_m',                 ('multiscale', 'hy_m')],
        ['hz_m',                 ('multiscale', 'hz_m')],
        ['nksplit',              ('multiscale', 'nksplit')],
        ['nxysplit',             ('multiscale', 'nxysplit')],
        ['nxvacl_m',             ('multiscale', 'nxvacl_m')],
        ['nxvacr_m',             ('multiscale', 'nxvacr_m')],
        ['nx_origin_m',          ('multiscale', 'nx_origin_m')],
        ['ny_origin_m',          ('multiscale', 'ny_origin_m')],
        ['nz_origin_m',          ('multiscale', 'nz_origin_m')],
        ['file_macropoint',      ('multiscale', 'file_macropoint')],
        ['num_macropoint',       ('multiscale', 'num_macropoint')],
        ['set_ini_coor_vel',     ('multiscale', 'set_ini_coor_vel')],
        ['nmacro_write_group',   ('multiscale', 'nmacro_write_group')],
    ]],
    ['maxwell', [
        ['al_em',                ('maxwell', 'al_em')],
        ['dl_em',                ('maxwell', 'dl_em')],
        ['dt_em',                ('maxwell', 'dt_em')],
        ['nt_em',                ('maxwell', 'nt_em')],
        ['boundary_em',          ('maxwell', 'boundary_em')],
        ['shape_file',           ('maxwell', 'shape_file')],
        ['imedia_num',           ('maxwell', 'imedia_num')],
        ['type_media',           ('maxwell', 'type_media')],
        ['epsilon',              ('maxwell', 'epsilon')],
        ['rmu',                  ('maxwell', 'rmu')],
        ['sigma',                ('maxwell', 'sigma')],
        ['pole_num_ld',          ('maxwell', 'pole_num_ld')],
        ['omega_p_ld',           ('maxwell', 'omega_p_ld')],
        ['f_ld',                 ('maxwell', 'f_ld')],
        ['gamma_ld',             ('maxwell', 'gamma_ld')],
        ['omega_ld',             ('maxwell', 'omega_ld')],
        ['wave_input',           ('maxwell', 'wave_input')],
        ['ek_dir1',              ('maxwell', 'ek_dir1')],
        ['source_loc1',          ('maxwell', 'source_loc1')],
        ['ek_dir2',              ('maxwell', 'ek_dir2')],
        ['source_loc2',          ('maxwell', 'source_loc2')],
        ['iobs_num_em',          ('maxwell', 'iobs_num_em')],
        ['iobs_samp_em',         ('maxwell', 'iobs_samp_em')],
        ['obs_loc_em',           ('maxwell', 'obs_loc_em')],
        ['obs_plane_em',         ('maxwell', 'obs_plane_em')],
        ['wf_em',                ('maxwell', 'wf_em')],
    ]],
    ['analysis', [
        ['projection_option',    ('analysis', 'projection_option')],
        ['projection_decomp',    ('analysis', 'projection_decomp')],
        ['nenergy',              ('analysis', 'nenergy')],
        ['de',                   ('analysis', 'de')],
        ['yn_out_psi',           ('analysis', 'out_psi')],
        ['yn_out_dos',           ('analysis', 'out_dos')],
        ['yn_out_dos_set_fe_origin', ('analysis', 'out_dos_fshift')],
        ['out_dos_start',        ('analysis', 'out_dos_start')],
        ['out_dos_end',          ('analysis', 'out_dos_end')],
        ['out_dos_nenergy',      ('analysis', 'iout_dos_nenergy')],
        ['out_dos_width',        ('analysis', 'out_dos_smearing')],
        ['out_dos_function',     ('analysis', 'out_dos_method')],
        ['yn_out_pdos',          ('analysis', 'out_pdos')],
        ['yn_out_dns',           ('analysis', 'out_dns')],
        ['yn_out_dns_rt',        ('analysis', 'out_dns_rt')],
        ['out_dns_rt_step',      ('analysis', 'out_dns_rt_step')],
        ['out_old_dns',          ('analysis', 'out_old_dns')],
        ['yn_out_dns_trans',     ('analysis', 'out_dns_trans')],
        ['out_dns_trans_energy', ('analysis', 'out_dns_trans_energy')],
        ['yn_out_elf',           ('analysis', 'out_elf')],
        ['yn_out_elf_rt',        ('analysis', 'out_elf_rt')],
        ['out_elf_rt_step',      ('analysis', 'out_elf_rt_step')],
        ['yn_out_estatic_rt',    ('analysis', 'out_estatic_rt')],
        ['out_estatic_rt_step',  ('analysis', 'out_estatic_rt_step')],
        ['yn_out_rvf_rt',        ('analysis', 'out_rvf_rt')],
        ['out_rvf_rt_step',      ('analysis', 'out_rvf_rt_step')],
        ['yn_out_tm',            ('analysis', 'out_tm')],
        ['out_projection_step',  ('analysis', 'out_projection_step')],
        ['out_ms_step',          ('analysis', 'out_ms_step')],
        ['format_voxel_data',    ('analysis', 'format3d')],
        ['nsplit_voxel_data',    ('analysis', 'numfiles_out_3d')],
        ['timer_process',        ('analysis', 'timer_process')],
    ]],
    ['poisson', [
        ['layout_multipole',     ('hartree', 'meo')],
        ['num_multipole_xyz',    ('hartree', 'num_pole_xyz')],
    ]],
    ['ewald', [
        ['newald',               ('ewald', 'newald')],
        ['aewald',               ('ewald', 'aewald')],
    ]],
    ['opt', [
        ['nopt',                 ('opt', 'nopt')],
        ['cg_alpha_ini',         ('opt', 'cg_alpha_ini')],
        ['cg_alpha_up',          ('opt', 'cg_alpha_up')],
        ['cg_alpha_down',        ('opt', 'cg_alpha_down')],
        ['convrg_scf_force',     ('opt', 'convrg_scf_force')],
        ['convrg_scf_ene',       ('opt', 'convrg_scf_ene')],
        ['convrg_opt_fmax',      ('opt', 'convrg_opt_fmax')],
        ['convrg_opt_ene',       ('opt', 'convrg_opt_ene')],
    ]],
    ['md', [
        ['ensemble',             ('md', 'ensemble')],
        ['thermostat',           ('md', 'thermostat')],
        ['step_velocity_scaling', ('md', 'step_velocity_scaling')],
        ['step_update_ps',       ('md', 'step_update_ps')],
        ['step_update_ps2',      ('md', 'step_update_ps2')],
        ['temperature0_ion_k',   ('md', 'temperature0_ion')],
        ['yn_set_ini_velocity',  ('md', 'set_ini_velocity')],
        ['file_ini_velocity',    ('md', 'file_ini_velocity')],
        ['thermostat_tau',       ('md', 'thermostat_tau')],
        ['friction',             ('md', 'friction')],
        ['yn_stop_system_momt',  ('md', 'stop_system_momt')],
    ]],
    ['group_fundamental', [
        ['iditerybcg',           ('group_fundamental', 'iditerybcg')],
        ['iditer_nosubspace_diag', ('group_fundamental', 'iditer_nosubspace_diag')],
        ['ntmg',                 ('group_fundamental', 'ntmg')],
        ['idisnum',              ('group_fundamental', 'idisnum')],
        ['iwrite_projection',    ('group_fundamental', 'iwrite_projection')],
        ['itwproj',              ('group_fundamental', 'itwproj')],
        ['iwrite_projnum',       ('group_fundamental', 'iwrite_projnum')],
        ['itcalc_ene',           ('group_fundamental', 'itcalc_ene')],
    ]],
    ['group_parallel', [
        ['imesh_s_all',          ('group_parallel', 'imesh_s_all')],
        ['iflag_comm_rho',       ('group_parallel', 'iflag_comm_rho')],
    ]],
    ['group_hartree', [
        ['hconv',                ('group_hartree', 'hconv')],
        ['lmax_lmp',             ('group_hartree', 'lmax_meo')],
    ]],
    ['group_file', [
        ['ic',                   ('group_file', 'ic')],
        ['oc',                   ('group_file', 'oc')],
        ['ic_rt',                ('group_file', 'ic_rt')],
        ['oc_rt',                ('group_file', 'oc_rt')],
    ]],
    ['group_others', [
        ['iscf_order',           ('group_others', 'iscf_order')],
        ['iswitch_orbital_mesh', ('group_others', 'iswitch_orbital_mesh')],
        ['iflag_psicube',        ('group_others', 'iflag_psicube')],
        ['lambda1_diis',         ('group_others', 'lambda1_diis')],
        ['lambda2_diis',         ('group_others', 'lambda2_diis')],
        ['file_ini',             ('group_others', 'file_ini')],
        ['num_projection',       ('group_others', 'num_projection')],
        ['iwrite_projection_ob', ('group_others', 'iwrite_projection_ob')],
        ['iwrite_projection_k',  ('group_others', 'iwrite_projection_k')],
        ['filename_pot',         ('group_others', 'filename_pot')],
        ['iwrite_external',      ('group_others', 'iwrite_external')],
        ['iflag_intelectron',    ('group_others', 'iflag_intelectron')],
        ['num_dip2',             ('group_others', 'num_dip2')],
        ['dip2boundary',         ('group_others', 'dip2boundary')],
        ['dip2center',           ('group_others', 'dip2center')],
        ['itotntime2',           ('group_others', 'itotntime2')],
        ['iwdenoption',          ('group_others', 'iwdenoption')],
        ['iwdenstep',            ('group_others', 'iwdenstep')],
        ['iflag_estatic',        ('group_others', 'iflag_estatic')],
    ]],
    ['code', [
        ['want_stencil_openmp_parallelization', ('code', 'want_stencil_openmp_parallelization')],
        ['want_stencil_hand_vectorization', ('code', 'want_stencil_hand_vectorization')],
        ['force_stencil_openmp_parallelization', ('code', 'force_stencil_openmp_parallelization')],
        ['force_stencil_sequential_computation', ('code', 'force_stencil_sequential_computation')],
    ]],
]

import sys
import re
import collections
import optparse

def special_rule(group, name, old_input):
    # iperiodic -> system/yn_periodic
    if (group, name) == ('system', 'yn_periodic'):
        if int(old_input['iperiodic']['']) == 3:
            return {'': "'y'"}
        else:
            return {'': "'n'"}

    # fourier -> parallel/yn_ffte
    if (group, name) == ('parallel', 'yn_ffte'):
        if 'fourier' in old_input:
            if old_input['fourier'][''].upper() == "'FFTE'":
                return {'': "'y'"}
            else:
                return {'': "'n'"}

################################################################
################################################################
################################################################
################################################################
################################################################
################################################################
################################################################
################################################################


ptn_begin = re.compile(r'^\s*&(\w+)\s*$')
ptn_end = re.compile(r'^\s*/\s*$')
ptn_set = re.compile(r'^\s*(\w+)([\(\)\d\s,]*)\s*=\s*(.*)$')


def read_input(fh):
    buf = collections.defaultdict(dict)
    for line in fh:
        m_begin = ptn_begin.match(line)
        if m_begin:
            group = m_begin.group(1).lower()
            if group == 'atomic_coor' or group == 'atomic_red_coor':
                buf[group] = []
                for tmp in fh:
                    m_end = ptn_end.match(tmp)
                    if m_end:
                        break
                    buf[group] += [tmp.strip()]

        m_set = ptn_set.match(line)
        if m_set:
            name = m_set.group(1).lower()
            index = m_set.group(2).strip()
            var = m_set.group(3).strip()
            buf[name][index] = var
    return buf


def main():
    parser = optparse.OptionParser(usage=USAGE)
    parser.add_option('-i', '--input', dest='INPUT')
    parser.add_option('-o', '--output', dest='OUTPUT')
    opts, args = parser.parse_args()

    if opts.INPUT:
        fh_inp = open(opts.INPUT, 'r')
    else:
        fh_inp = sys.stdin
    
    old_input = read_input(fh_inp)
    fh_inp.close()

    buf = []
    for group, group_items in keyword_table:
        tmp = []
        for name, target in group_items:
            # Search values:
            if target is None:
                result = special_rule(group, name, old_input)
            else:
                old_group, old_name = target
                if old_name in old_input:
                    result = old_input[old_name]
                else:
                    result = None
            if result:
                for index, value in result.items():
                    tmp += ['%s%s = %s' % (name, index, value)]
        # Create input group
        if 0 < len(tmp):
            buf += ['&%s' % group]
            buf += tmp
            buf += ['/\n']

    if 'atomic_coor' in old_input:
        buf += ['&atomic_coor'] + old_input['atomic_coor'] + ['/', '']

    if 'atomic_red_coor' in old_input:
        buf += ['&atomic_red_coor'] + old_input['atomic_red_coor'] + ['/', '']

    if opts.OUTPUT:
        fh_out = open(opts.OUTPUT, 'w')
    else:
        fh_out = sys.stdout

    fh_out.write('! This file is automatically generated for v.2.0.0.\n')
    fh_out.write('\n'.join(buf) + '\n')
    fh_out.close()




if __name__ == '__main__':
    main()
