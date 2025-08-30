from setuptools import setup, find_packages

setup(name="casparing", 
        version = "1.0",
        description='caspar rederivation and plotting tool kit',
        author = "Sarah Betti",
        author_email = "sbetti@stsci.edu",
        license='MIT',
        url = "https://github.com/sbetti22/casparing",
        classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research ',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License'
        ],
        keywords='accretion',
        packages = ['plot', 'derive', 'data'],
        install_requires=[
        'numpy', 'pandas', 'astropy', 'matplotlib', 'astroquery', 'tqdm', 'pandarallel', 'scipy', 'uncertainties', 'extinction', 'isochrones', 'banyansigma', 
        ], 
        package_data={'data': ['*']},
        zip_safe=False
        
        )