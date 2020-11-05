import setuptools

setuptools.setup(
    name='ries',
    version='0.1.0',
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