# Design-principles-growth-laws-and-competition-of-minimal-autocatalysts
Code and Data for the paper: Design principles, growth laws, and competition of minimal autocatalysts. Yann Sakref &amp; Olivier Rivoire, 
https://doi.org/10.48550/arXiv.2403.19047

Codes:

I - Molecular Dynamics (MD) Simulations:

The MD simulation codes are available in the MD_Simulations folder. These Python codes are based on the HoomD-Blue package (version 3.5.0) [JA. Anderson, J. Glaser, SC. Glotzer. "HOOMD-blue: A Python package for high-performance molecular dynamics and hard particle Monte Carlo simulations." Comput. Mater. Sci., 173:109363, 2020]. It is a very elegant package, hence the codes are straightfoward, although we use extensions of Wang-Frankel potentials not natively supported by HoomD-Blue.
- Isotropic_Particles.py: Simulates four isotropic particles as described in the main section of the paper.
- Anisotropic_Particles.py: Extends the simulations to anisotropic particles, as detailed in the Supplementary Information.

II - Figure Generation:

Codes to generate the main figures and most supplementary figures (excluding redundant ones) are provided.

Data:

Data.zip: Contains the data used to generate the figures in the paper.

Feel free to reach out if you have any questions.
