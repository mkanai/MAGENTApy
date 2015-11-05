# MAGENTApy
Experimental port of MAGENTA software to Python. https://www.broadinstitute.org/mpg/magenta/

**Not fully tested yet!**

Note: For maintainability, this script is **intentionally un-pythonic**, following the original code structure of MAGENTA written in MATLAB.

## Usage
```{shell}
git clone https://github.com/mkanai/MAGENTApy.git
cd MAGENTApy
python Run_MAGENTA.py
```

## Configuration
Edit `config.yml`, or specify your own configuration by `--config myconfig.yml`. All parameters are compatible with the original MAGENTA.


## Requirements
- numpy
- scipy
- pandas
- pyaml


## Copyright
MAGENTA is copyright by Ayellet Segre, Mark Daly, and David Altshuler of The Broad Institute.
