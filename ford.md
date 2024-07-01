---
project: cloud_model
summary: cloud simulation.
project_github: https://github.com/Diegoolei/Fortran77-Cloud-Model
author: Diego E. Oleiarz
author_description: Computer science student.
email: diego.oleiarz@mi.unc.edu.ar
github: https://github.com/Diegoolei
src_dir: src
exclude_dir: test doc src/adiff/hyperdual.f90 src/adiff/autodiff_api/tapenade
output_dir: doc/ford_site
preprocessor: gfortran -E
display: public
         protected
         private
source: false
proc_internals: true
sort: permission-alpha
docmark_alt: !>
docmark: !
predocmark_alt: *
print_creation_date: true
creation_date: %Y-%m-%d %H:%M %z
md_extensions: markdown.extensions.toc
               markdown.extensions.smarty
graph: false
license: MPL
page_dir: doc/page
media_dir: doc/media
---

[TOC]

{!README.md!}