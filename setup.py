from setuptools import setup, find_packages

setup(
    name = 'lab_tools',
    version= '1.0.0',
    packages = find_packages(where='src'),
    package_dir = {'': 'src'},
    install_requires = [
        #list dependencies here
        #'re',
        #'typing',
        #'collections'
    ],
    entry_points = {
        'console_scripts':[
            # Define cmd line scripts
        ],
    },
    author= 'Dario C. Lewczyk',
    author_email='darlewczyk@gmail.com',
    description= 'A suite of tools for assisting in basic laboratory tasks like experiment design',
    long_description=open('README.md').read(),
    long_description_content_type = 'text/markdown',
    url = 'https://github.com/DarioLewczyk/lab_tools',
    classifiers = [
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',

    ],
    python_requires = '>=3.6',

)
