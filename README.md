[![GitHub license](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://github.com/google/eidos-audition/blob/master/LICENSE)
[![Paper](https://img.shields.io/badge/paper-SLTU%202020-blue.svg)](https://www.aclweb.org/anthology/2020.sltu-1.2/)
[![Build Tests (Linux)](https://github.com/google/eidos-audition/workflows/linux/badge.svg)](https://github.com/google/eidos-audition/actions?query=workflow%3A%22linux%22)
[![Build Tests (macOS)](https://github.com/google/eidos-audition/workflows/macos/badge.svg)](https://github.com/google/eidos-audition/actions?query=workflow%3A%22macos%22)

## About

This library represents a simple pipeline integrating various models
of auditory periphery developed over the years in academia. The pipeline
can be used to combine models representing various stages of auditory
perception for extracting features for machine learning algorithms.

Detailed documentation for the auditory library can be found
[here](eidos/audition/README.md).

## Citation

If you use the code in your research, please cite the following
[paper](https://www.aclweb.org/anthology/2020.sltu-1.2/):

```
@inproceedings{eidos2020,
  title   = {{Eidos: An Open-Source Auditory Periphery Modeling Toolkit and
             Evaluation of Cross-Lingual Phonemic Contrasts}},
  author  = {Alexander Gutkin},
  booktitle = {Proc. of 1st Joint Spoken Language Technologies for Under-Resourced
               Languages (SLTU) and Collaboration and Computing for Under-Resourced
               Languages (CCURL) Workshop (SLTU-CCURL 2020)},
  pages   = {9--20},
  year    = {2020},
  month   = may,
  address = {Marseille, France},
  publisher = {European Language Resources Association (ELRA)},
  url = {https://lrec2020.lrec-conf.org/media/proceedings/Workshops/Books/SLTUCCURLbook.pdf},
}
```

Please also cite the reference publications describing the particular algorithms
from this collection that you use. These can be found in the respective model
directories under [third_party/audition/models](third_party/audition/models).

## Mandatory Disclaimer

This is not an officially supported Google product.
