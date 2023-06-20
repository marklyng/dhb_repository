All genomes were assembled using Trycycler V0.5.3 and annotated with BAKTA V1.7.
The gff3-annotation files from each genome were used in Panaroo V1.3.3 as:

panaroo <input>.gff3 -c 0.4 -f 0.3 --clean-mode moderate

The roary-formatted presence-absence matrix was then imported into R V4.2.2 using Rstudio 2023.03.0 B386 with Tidyverse V2.0.0, 
and manually inspected for co-ocurring coding sequences.

Type VI secretion system predictions were performed with SecReT6 V3


References:
RStudio Team (2020). 
RStudio: Integrated Development for R. RStudio, PBC, Boston, MA 
http://www.rstudio.com/

Wick, R.R., Judd, L.M., Cerdeira, L.T. et al. 
Trycycler: consensus long-read assemblies for bacterial genomes. Genome Biol 22, 266 (2021). 
https://doi.org/10.1186/s13059-021-02483-z

Schwengers O., Jelonek L., Dieckmann M. A., Beyvers S., Blom J., Goesmann A. (2021). 
Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification. Microbial Genomics, 7(11). 
https://doi.org/10.1099/mgen.0.000685

Tonkin-Hill, G., MacAlasdair, N., Ruis, C. et al. 
Producing polished prokaryotic pangenomes with the Panaroo pipeline. Genome Biol 21, 180 (2020). 
https://doi.org/10.1186/s13059-020-02090-4

Zhang J., Guan J., Wang M, Li G, Djordjevic M, Tai C, Wang H, Deng Z, Chen Z*, Ou HY* (2022) 
SecReT6 update: a comprehensive resource of bacterial Type VI Secretion Systems. SCIENCE CHINA Life Sciences. 
doi: 10.1007/s11427-022-2172-x.