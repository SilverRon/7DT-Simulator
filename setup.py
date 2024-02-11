from setuptools import setup, find_packages

# requirements.txt 파일에서 의존성을 읽어옵니다.
with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name='7DT-Simulator',
    version='0.1.0',
    author='Gregory Paek',
    author_email='gregorypaek94@gmail.com',
    description='A brief description of the 7DT-Simulator package',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/yourusername/7DT-Simulator',
    packages=find_packages(),
    install_requires=required,
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.10',
)

