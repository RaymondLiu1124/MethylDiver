MethylDiver documentation
=========================

MethylDiver calculates methylation level, normalized pattern entropy, and
Epiallele distance from long-read modBAMs. This guide describes version 2.4.1.

Use ``compare`` for two samples, ``single`` for one sample, ``read`` to
preview or export retained results, and ``discover`` to rank high-ED windows.
Use ``cluster`` to fit and apply read clusters or compare candidate cluster counts.
Units contain three consecutive reference
CpGs by default; ``--N`` supports 1–10.

.. toctree::
   :maxdepth: 2
   :caption: User guide

   installation
   concepts
   modifications
   compare
   single
   regions
   discover
   cluster
   commands
   reading
   examples
   methods
   python-api
   development

`Software repository <https://github.com/RaymondLiu1124/MethylDiver>`_ |
`Release notes <https://github.com/RaymondLiu1124/MethylDiver/blob/main/CHANGELOG.md>`_
