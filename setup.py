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

import setuptools

setuptools.setup(
    name='ries',
    version='0.2.0',
    description='resonances integrated over energy and space',
    author='Udo Friman-Gayer',
    author_email='ufg@email.unc.edu',
    classifiers=[
        "Programming Language :: Python :: 3 :: Only",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
    ],
    packages=setuptools.find_packages(),
    include_package_data=True,
    python_requires='>=3',
    install_requires=['numpy', 'scipy']
)