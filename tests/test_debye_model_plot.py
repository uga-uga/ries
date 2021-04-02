# This file is part of ries.

# ries is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# ries is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with ries.  If not, see <https://www.gnu.org/licenses/>.

import matplotlib.pyplot as plt
import numpy as np

from ries.resonance.debye_model import effective_temperature_debye_approximation


def test_debye_model_plot():
    T_D = 1.0

    T_over_T_D_max = 1.5
    T_over_T_D = np.linspace(1e-2, T_over_T_D_max)
    T = T_over_T_D * T_D

    _fontsize_axis_label = 16
    _fontsize_legend = 12
    _fontsize_text = 16
    _fontsize_ticks = 13

    fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    ax.set_xlabel(r"$T / \Theta_D$", fontsize=_fontsize_axis_label)
    ax.set_xlim(0.0, T_over_T_D_max)
    ax.set_ylim(0.0, T_over_T_D_max)
    ax.set_ylabel(r"$T_\mathrm{eff} / \Theta_D$", fontsize=_fontsize_axis_label)
    ax.tick_params(labelsize=_fontsize_ticks)
    ax.plot(
        T_over_T_D,
        effective_temperature_debye_approximation(T, T_D) / T_D,
        color="black",
        label="Solid",
    )
    ax.plot(T_over_T_D, T_over_T_D, "--", color="black", label="Ideal Gas")

    ax.plot([0.0, 1.0], 3.0 / 8.0 * np.array([1.0, 1.0]), ":", color="black")
    ax.annotate(
        "",
        [1.0, 0.0],
        [1.0, 3.0 / 8.0],
        arrowprops=dict(arrowstyle="<|-|>", facecolor="black"),
    )
    ax.text(
        1.05,
        3.0 / 16.0,
        r"$\frac{3}{8}$",
        verticalalignment="center",
        fontsize=_fontsize_text,
    )

    ax.legend(fontsize=_fontsize_legend)

    plt.savefig("debye_model_plot.pdf", bbox_inches="tight")
