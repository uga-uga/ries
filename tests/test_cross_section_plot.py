import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
from scipy.integrate import cumtrapz

from ries.resonance.debye_model import effective_temperature_debye_approximation, room_temperature_T_D
from ries.resonance.voigt import Voigt
from ries.nonresonant.klein_nishina import KleinNishina
from ries.nonresonant.xrmac import xrmac

from .boron import B11, natural_boron

def test_cross_section_plot():
    cross_section_nonresonant = xrmac['B']
    cross_section_resonant = sum([
        Voigt(
                B11.ground_state,
                B11.excited_states[excited_state],
                B11.amu, effective_temperature_debye_approximation(293., room_temperature_T_D['B'])
        ) for excited_state in B11.excited_states
    ])
    cross_section = cross_section_nonresonant + cross_section_resonant
    cross_section_compton = KleinNishina(natural_boron.Z)

    ene_lim = (1e-2, 12.)
    ene = cross_section.equidistant_probability_grid(ene_lim, 1000)
    cro_sec = cross_section(ene)
    cro_sec_int = cumtrapz(cro_sec, ene)
    cro_sec_non_res = cross_section_nonresonant(ene)
    cro_sec_non_res_int = cumtrapz(cro_sec_non_res, ene)
    cro_sec_res = cross_section_resonant(ene)
    cro_sec_res_int = cumtrapz(cro_sec_res, ene)
    cro_sec_com = cross_section_compton(ene)
    cro_sec_com_int = cumtrapz(cro_sec_com, ene)

    xlabel = 'Energy (MeV)'

    fig, ax = plt.subplots(1,1, figsize=(6, 7))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis('off')
    ax_cro_sec = inset_axes(
        ax, width='100%', height='100%',
        bbox_to_anchor=[0.1, 0.75, 0.8, 0.3],
        bbox_transform=ax.transAxes
    )
    ax_cro_sec.set_xlabel(xlabel)
    ax_cro_sec.set_xlim(ene_lim)
    ax_cro_sec.set_ylabel(r'$\sigma_\mathrm{NR} + \sigma_\mathrm{R}$ (fm$^2$)')
    ax_cro_sec.semilogy(ene, cro_sec, color='black')
    ax_cro_sec.semilogy(ene, cro_sec_com, '--', color='black')

    ax_cro_sec_non_res = inset_axes(
        ax, width='100%', height='100%',
        bbox_to_anchor=[0.05, 0., 0.4, 0.3],
        bbox_transform=ax.transAxes
    )
    ax_cro_sec_non_res.set_xlabel(xlabel)
    ax_cro_sec_non_res.set_xlim(ene_lim)
    ax_cro_sec_non_res.set_ylabel(r'$\sigma_\mathrm{NR}$ (fm$^2$)')
    ax_cro_sec_non_res.semilogy(ene, cro_sec_non_res, color='black')
    ax_cro_sec_non_res.semilogy(ene, cro_sec_com, '--', color='black')

    ax_cro_sec_non_res_int = inset_axes(
        ax, width='100%', height='100%',
        bbox_to_anchor=[0.05, 0.3, 0.4, 0.3],
        bbox_transform=ax.transAxes
    )
    ax_cro_sec_non_res_int.set_xticks([])
    ax_cro_sec_non_res_int.set_xlabel('')
    ax_cro_sec_non_res_int.set_xlim(ene_lim)
    ax_cro_sec_non_res_int.set_ylabel(r'$\int \sigma_\mathrm{NR} \mathrm{d} E$ (MeV fm$^2$)')
    ax_cro_sec_non_res_int.plot(ene[:-1], cro_sec_non_res_int, color='black')
    ax_cro_sec_non_res_int.plot(ene[:-1], cro_sec_com_int, '--', color='black')

    ax_cro_sec_res = inset_axes(
        ax, width='100%', height='100%',
        bbox_to_anchor=[0.55, 0., 0.4, 0.3],
        bbox_transform=ax.transAxes
    )
    ax_cro_sec_res.set_xlabel(xlabel)
    ax_cro_sec_res.set_xlim(ene_lim)
    ax_cro_sec_res.set_yticks([])
    ax_cro_sec_res_2 = ax_cro_sec_res.twinx()
    ax_cro_sec_res_2.set_ylabel(r'$\sigma_\mathrm{R}$ (fm$^2$)')
    ax_cro_sec_res_2.semilogy(ene, cro_sec_res, color='black')

    ax_cro_sec_res_int = inset_axes(
        ax, width='100%', height='100%',
        bbox_to_anchor=[0.55, 0.3, 0.4, 0.3],
        bbox_transform=ax.transAxes
    )
    ax_cro_sec_res_int.set_xticks([])
    ax_cro_sec_res_int.set_xlim(ene_lim)
    ax_cro_sec_res_int.set_yticks([])
    ax_cro_sec_res_int_2 = ax_cro_sec_res_int.twinx()
    ax_cro_sec_res_int_2.set_ylabel(r'$\int \sigma_\mathrm{NR} \mathrm{d} E$ (MeV fm$^2$)')
    ax_cro_sec_res_int_2.plot(ene[:-1], cro_sec_res_int, color='black')

    plt.savefig('11B.pdf')